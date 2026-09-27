#include "fvm/mesh/AirfoilMesh.h"

#include "fvm/mesh/Extrude2D.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace cfd {

namespace {

// 单侧 tanh 拉伸：s(0)=0，s(1)=1，s(1/n)=d1（d1 < 1/n）
struct TanhStretch {
    scalar beta = 0;
    explicit TanhStretch(scalar d1, int n) {
        const scalar eta1 = 1.0 / n;
        if (d1 >= eta1) {
            beta = 0;
            return;
        }
        scalar lo = 1e-6, hi = 50;
        for (int it = 0; it < 200; ++it) {
            const scalar b = 0.5 * (lo + hi);
            (eval(b, eta1) > d1 ? lo : hi) = b;
        }
        beta = 0.5 * (lo + hi);
    }
    static scalar eval(scalar b, scalar eta) { return 1.0 + std::tanh(b * (eta - 1.0)) / std::tanh(b); }
    scalar operator()(scalar eta) const { return beta == 0 ? eta : eval(beta, eta); }
};

// 从 x0 出发、首段 d0、共 n 段、总长 L 的几何增长分布
std::vector<scalar> geometricDistribution(scalar d0, int n, scalar L) {
    // 求增长率 r：d0 (r^n − 1)/(r − 1) = L
    scalar lo = 1.0 + 1e-9, hi = 2.0;
    auto total = [&](scalar r) { return d0 * (std::pow(r, n) - 1) / (r - 1); };
    if (d0 * n >= L) {
        std::vector<scalar> x(n + 1);
        for (int k = 0; k <= n; ++k) x[k] = L * k / n;
        return x;
    }
    while (total(hi) < L) hi *= 1.5;
    for (int it = 0; it < 200; ++it) {
        const scalar r = 0.5 * (lo + hi);
        (total(r) < L ? lo : hi) = r;
    }
    const scalar r = 0.5 * (lo + hi);
    std::vector<scalar> x(n + 1, 0.0);
    scalar d = d0;
    for (int k = 1; k <= n; ++k) {
        x[k] = x[k - 1] + d;
        d *= r;
    }
    for (auto& v : x) v *= L / x[n];
    return x;
}

} // namespace

Vec3 nacaSurfacePoint(const AirfoilMeshSpec& s, scalar xc, bool upper) {
    if (s.naca.size() != 4) throw std::runtime_error("airfoil: only NACA 4-digit supported");
    const scalar m = (s.naca[0] - '0') / 100.0, p = (s.naca[1] - '0') / 10.0;
    const scalar t = std::stoi(s.naca.substr(2)) / 100.0;
    const scalar a4 = s.closedTE ? -0.1036 : -0.1015;
    const scalar yt = 5 * t * (0.2969 * std::sqrt(std::max(xc, 0.0)) - 0.1260 * xc - 0.3516 * xc * xc +
                               0.2843 * xc * xc * xc + a4 * xc * xc * xc * xc);
    scalar yc = 0, dyc = 0;
    if (m > 0 && p > 0) {
        if (xc < p) {
            yc = m / (p * p) * (2 * p * xc - xc * xc);
            dyc = 2 * m / (p * p) * (p - xc);
        } else {
            yc = m / ((1 - p) * (1 - p)) * ((1 - 2 * p) + 2 * p * xc - xc * xc);
            dyc = 2 * m / ((1 - p) * (1 - p)) * (p - xc);
        }
    }
    const scalar th = std::atan(dyc);
    const scalar sg = upper ? 1.0 : -1.0;
    return {s.chord * (xc - sg * yt * std::sin(th)), s.chord * (yc + sg * yt * std::cos(th)), 0};
}

RawMesh generateAirfoilCMesh(const AirfoilMeshSpec& s) {
    if (s.nAirfoil % 2) throw std::runtime_error("airfoil: nAirfoil must be even");
    const scalar c = s.chord;
    const scalar a = -s.alphaDeg * PI / 180.0;  // 绕半弦点旋转 −α（抬头为正攻角）
    const Vec3 pivot{0.5 * c, 0, 0};
    const scalar xIn = 0.5 * c - s.upstream * c, xOut = 0.5 * c + s.downstream * c;
    const scalar yLo = -s.halfHeight * c, yHi = s.halfHeight * c;

    // ------------------------------------------------ 内边界（j=0）
    const int nh = s.nAirfoil / 2, nw = s.nWake;
    std::vector<Vec3> S;  // 下尾迹 + 下表面 + 上表面 + 上尾迹
    std::vector<Vec3> surf;
    // 弦向分布：余弦（前后缘都加密）与半余弦（只加密前缘）各半，后缘不至于过密
    auto xDist = [&](int k) {
        const scalar th = PI * scalar(k) / nh;
        return 0.5 * 0.5 * (1 - std::cos(th)) + 0.5 * (1 - std::cos(0.5 * th));
    };
    for (int k = nh; k >= 0; --k) surf.push_back(nacaSurfacePoint(s, xDist(k), false));  // 下表面：后缘 → 前缘
    for (int k = 1; k <= nh; ++k) surf.push_back(nacaSurfacePoint(s, xDist(k), true));   // 上表面：前缘 → 后缘
    // 闭合尾缘：上下后缘点取同一点
    Vec3 te = 0.5 * (surf.front() + surf.back());
    surf.front() = te;
    surf.back() = te;
    if (te.x >= xOut) throw std::runtime_error("airfoil: trailing edge outside domain");
    const scalar dTE = mag(surf[1] - surf[0]);
    const auto wx = geometricDistribution(dTE, nw, xOut - te.x);
    for (int k = nw; k >= 1; --k) S.push_back({te.x + wx[k], te.y, 0});  // 下尾迹：出口 → 后缘（不含后缘）
    const int iTElo = int(S.size());
    for (auto& p : surf) S.push_back(p);
    const int iTEup = int(S.size()) - 1;
    for (int k = 1; k <= nw; ++k) S.push_back({te.x + wx[k], te.y, 0});
    const int ni = int(S.size());  // 点数（i 方向）

    // 内边界外法向（光顺）
    std::vector<Vec3> nrm(ni);
    for (int i = 0; i < ni; ++i) {
        if (i < iTElo) {
            nrm[i] = {0, -1, 0};
        } else if (i > iTEup) {
            nrm[i] = {0, 1, 0};
        } else {
            const Vec3 t = S[std::min(i + 1, iTEup)] - S[std::max(i - 1, iTElo)];
            nrm[i] = Vec3{-t.y, t.x, 0} / mag(t);  // 沿 i 走向（绕翼型顺时针）的左侧 = 外侧
        }
    }
    nrm[iTElo] = Vec3{0, -1, 0};
    nrm[iTEup] = Vec3{0, 1, 0};
    for (int pass = 0; pass < 20; ++pass) {  // 后缘附近法向过渡
        auto n0 = nrm;
        for (int i = iTElo - 6; i <= iTElo + 6; ++i)
            if (i > 0 && i < ni - 1) nrm[i] = (n0[i - 1] + 2 * n0[i] + n0[i + 1]) / mag(n0[i - 1] + 2 * n0[i] + n0[i + 1]);
        for (int i = iTEup - 6; i <= iTEup + 6; ++i)
            if (i > 0 && i < ni - 1) nrm[i] = (n0[i - 1] + 2 * n0[i] + n0[i + 1]) / mag(n0[i - 1] + 2 * n0[i] + n0[i + 1]);
    }

    // ------------------------------------------------ 外边界（j=nj）
    std::vector<Vec3> B(ni);
    {
        // 尾迹段外边界：从后缘正上/下方开始按外边界平均间距几何增长到出口，避免细网格线一直延伸到远场
        const scalar Lc = 2 * (te.x - xIn) + (yHi - yLo);
        const auto ox = geometricDistribution(Lc / (iTEup - iTElo), nw, xOut - te.x);
        for (int k = 1; k <= nw; ++k) {
            B[iTElo - k] = {te.x + ox[k], yLo, 0};
            B[iTEup + k] = {te.x + ox[k], yHi, 0};
        }
    }
    {
        // C 形路径：(te.x, yLo) → (xIn, yLo) → (xIn, yHi) → (te.x, yHi)
        const scalar L1 = te.x - xIn, L2 = yHi - yLo, L3 = te.x - xIn, Lt = L1 + L2 + L3;
        // 三段各自均匀分布，角点上一定有网格点
        const int na = iTEup - iTElo;
        const int n1 = std::max(1, int(std::lround(na * L1 / Lt))), n2 = std::max(1, int(std::lround(na * L2 / Lt)));
        const int n3 = na - n1 - n2;
        if (n3 < 1) throw std::runtime_error("airfoil: nAirfoil too small for outer boundary");
        for (int k = 0; k <= na; ++k) {
            Vec3 p;
            if (k <= n1) p = {te.x - L1 * k / n1, yLo, 0};
            else if (k <= n1 + n2) p = {xIn, yLo + L2 * (k - n1) / n2, 0};
            else p = {xIn + L3 * (k - n1 - n2) / n3, yHi, 0};
            B[iTElo + k] = p;
        }
    }

    // ------------------------------------------------ 内部点（Hermite + 法向拉伸）
    const int nj = s.nNormal;
    std::vector<Vec3> X(std::size_t(ni) * (nj + 1));
    auto at = [&](int i, int j) -> Vec3& { return X[std::size_t(j) * ni + i]; };
    for (int i = 0; i < ni; ++i) {
        const Vec3 d = B[i] - S[i];
        const scalar L = mag(d);
        constexpr scalar kT = 0.5;  // 壁面切向量长度系数：越大网格线离壁越"直"，越小越贴近直线插值
        const TanhStretch st(s.firstCell * c / (kT * L), nj);
        const Vec3 T0 = kT * L * nrm[i], T1 = d;
        for (int j = 0; j <= nj; ++j) {
            const scalar t = st(scalar(j) / nj);
            const scalar h0 = 2 * t * t * t - 3 * t * t + 1, h1 = -2 * t * t * t + 3 * t * t;
            const scalar h2 = t * t * t - 2 * t * t + t, h3 = t * t * t - t * t;
            at(i, j) = h0 * S[i] + h1 * B[i] + h2 * T0 + h3 * T1;
        }
    }
    // 拉普拉斯光顺（i 方向），保留靠壁 nFix 层、外边界与两端出口线
    const int nFix = std::min(nj / 4, 20);
    for (int it = 0; it < s.smoothIter; ++it)
        for (int j = nFix; j < nj; ++j) {
            const scalar w = std::min(1.0, scalar(j - nFix + 1) / 10.0) * 0.5;
            std::vector<Vec3> row(ni);
            for (int i = 1; i < ni - 1; ++i) row[i] = 0.5 * (at(i - 1, j) + at(i + 1, j));
            for (int i = 1; i < ni - 1; ++i) {
                // 只沿网格线切向移动，保持法向分布
                const Vec3 tg = at(i, j + 1) - at(i, j - 1);
                const Vec3 dlt = row[i] - at(i, j);
                const scalar tt = dot(tg, tg);
                const Vec3 mv = tt > 0 ? dlt - (dot(dlt, tg) / tt) * tg : dlt;
                at(i, j) += w * mv;
            }
        }

    // 攻角：网格在 α=0 下生成，再绕半弦点整体旋转内区、在环形过渡区内把旋转角平滑降到 0。
    // 近壁网格随翼型刚性转动（保持正交性），外边界不动，尾迹割线在过渡区内逐渐转回水平。
    if (a != 0.0) {
        const scalar r1 = 1.0 * c, r2 = 0.8 * std::min({s.halfHeight, s.upstream, s.downstream}) * c;
        if (r2 <= r1) throw std::runtime_error("airfoil: domain too small for rotation blending");
        for (auto& p : X) {
            const scalar r = std::hypot(p.x - pivot.x, p.y - pivot.y);
            const scalar w = r <= r1 ? 1.0 : (r >= r2 ? 0.0 : 0.5 * (1 + std::cos(PI * (r - r1) / (r2 - r1))));
            if (w == 0.0) continue;
            const scalar aw = a * w;
            const scalar dx = p.x - pivot.x, dy = p.y - pivot.y;
            p = {pivot.x + dx * std::cos(aw) - dy * std::sin(aw), pivot.y + dx * std::sin(aw) + dy * std::cos(aw), 0};
        }
    }
    if (const char* dump = std::getenv("CFD_AIRFOIL_DUMP")) {
        std::ofstream os(dump);
        os << ni << ' ' << nj << '\n';
        for (auto& p : X) os << p.x << ' ' << p.y << '\n';
    }
    // ------------------------------------------------ 单元与边界
    Mesh2D m2;
    m2.points = X;
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < ni - 1; ++i) {
            auto id = [&](int ii, int jj) { return glabel(std::size_t(jj) * ni + ii); };
            m2.cells.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)});
        }
    // 尾迹割线两侧点重合，合并
    const scalar tol = 1e-3 * std::min(s.firstCell * c, dTE);
    const std::vector<Vec3> orig = m2.points;
    mergeDuplicatePoints(m2, tol);
    // 标记翼型表面点
    std::vector<char> onWall(m2.points.size(), 0);
    {
        // 合并后的点号：重新定位表面点
        for (int i = iTElo; i <= iTEup; ++i) {
            const Vec3& p = orig[i];
            for (std::size_t q = 0; q < m2.points.size(); ++q)
                if (std::hypot(m2.points[q].x - p.x, m2.points[q].y - p.y) < tol) {
                    onWall[q] = 1;
                    break;
                }
        }
    }
    std::vector<BoundaryPatchSpec> patches{{"airfoil", PatchType::Wall},
                                           {"inlet", PatchType::Patch},
                                           {"outlet", PatchType::Patch},
                                           {"top", PatchType::Patch},
                                           {"bottom", PatchType::Patch}};
    const scalar eps = 1e-9 * (xOut - xIn);
    auto classify = [&](glabel p, glabel q) {
        if (onWall[p] && onWall[q]) return 0;
        const Vec3 mid = 0.5 * (m2.points[p] + m2.points[q]);
        if (std::abs(mid.x - xIn) < eps) return 1;
        if (std::abs(mid.x - xOut) < eps) return 2;
        if (std::abs(mid.y - yHi) < eps) return 3;
        if (std::abs(mid.y - yLo) < eps) return 4;
        throw std::runtime_error("airfoil mesh: unclassified boundary edge at (" + std::to_string(mid.x) + ", " +
                                 std::to_string(mid.y) + ")");
    };
    return extrude2D(m2, patches, classify, s.depth * c);
}

} // namespace cfd
