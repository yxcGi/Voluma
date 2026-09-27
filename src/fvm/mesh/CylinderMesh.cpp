#include "fvm/mesh/CylinderMesh.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd {

void splitPatch(RawMesh& m, const std::string& name, const std::vector<std::pair<std::string, PatchType>>& parts,
                const std::function<int(const Vec3& faceCentre)>& classify) {
    std::size_t ip = m.patches.size();
    for (std::size_t k = 0; k < m.patches.size(); ++k)
        if (m.patches[k].name == name) ip = k;
    if (ip == m.patches.size()) throw std::runtime_error("splitPatch: no patch " + name);
    const RawPatch old = m.patches[ip];
    std::vector<std::vector<glabel>> groups(parts.size());
    for (glabel f = old.start; f < old.start + old.size; ++f) {
        Vec3 c{};
        const glabel b = m.faceOffsets[f], e = m.faceOffsets[f + 1];
        for (glabel q = b; q < e; ++q) c += m.points[m.facePoints[q]];
        c = c / scalar(e - b);
        const int g = classify(c);
        if (g < 0 || g >= int(parts.size())) throw std::runtime_error("splitPatch: bad group");
        groups[g].push_back(f);
    }
    // 重排该面组内的面（点表与 owner 一起重排）
    std::vector<glabel> order;
    for (auto& g : groups) order.insert(order.end(), g.begin(), g.end());
    std::vector<glabel> pts, offs{0}, own;
    for (glabel f : order) {
        for (glabel q = m.faceOffsets[f]; q < m.faceOffsets[f + 1]; ++q) pts.push_back(m.facePoints[q]);
        offs.push_back(glabel(pts.size()));
        own.push_back(m.owner[f]);
    }
    // 拼回：[0, start) + 新顺序 + [start + size, end)
    const glabel s0 = m.faceOffsets[old.start], s1 = m.faceOffsets[old.start + old.size];
    std::vector<glabel> fp(m.facePoints.begin(), m.facePoints.begin() + s0);
    fp.insert(fp.end(), pts.begin(), pts.end());
    fp.insert(fp.end(), m.facePoints.begin() + s1, m.facePoints.end());
    std::vector<glabel> fo(m.faceOffsets.begin(), m.faceOffsets.begin() + old.start + 1);
    for (std::size_t k = 1; k < offs.size(); ++k) fo.push_back(s0 + offs[k]);
    for (glabel f = old.start + old.size + 1; f < glabel(m.faceOffsets.size()); ++f) fo.push_back(m.faceOffsets[f]);
    m.facePoints = std::move(fp);
    m.faceOffsets = std::move(fo);
    std::copy(own.begin(), own.end(), m.owner.begin() + old.start);
    std::vector<RawPatch> np;
    for (std::size_t k = 0; k < m.patches.size(); ++k) {
        if (k != ip) {
            np.push_back(m.patches[k]);
            continue;
        }
        glabel start = old.start;
        for (std::size_t g = 0; g < parts.size(); ++g) {
            if (groups[g].empty()) continue;
            RawPatch rp;
            rp.name = parts[g].first;
            rp.type = parts[g].second;
            rp.start = start;
            rp.size = glabel(groups[g].size());
            start += rp.size;
            np.push_back(rp);
        }
    }
    m.patches = std::move(np);
}

RawMesh generateCylinderOMesh(const CylinderMeshSpec& s) {
    const scalar D = s.D, R = 0.5 * D;
    const scalar xIn = -s.upstream * D, xOut = s.downstream * D, yHi = s.halfHeight * D;
    if (s.upstream <= 0.5 || s.downstream <= 0.5 || s.halfHeight <= 0.5) throw std::runtime_error("cylinder: domain too small");
    const int n = s.nRadial;
    const scalar d0 = s.firstCell * D;
    // 射线到矩形外边界的长度
    auto rayLength = [&](scalar c, scalar sn) {
        scalar t = GREAT;
        if (c > 0) t = std::min(t, xOut / c);
        if (c < 0) t = std::min(t, xIn / c);
        if (sn > 0) t = std::min(t, yHi / sn);
        if (sn < 0) t = std::min(t, -yHi / sn);
        return t - R;
    };
    // 几何增长率：d0 (r^n − 1)/(r − 1) = L
    auto growth = [&](scalar L) {
        if (d0 * n >= L) return 1.0;
        scalar lo = 1.0 + 1e-12, hi = 2.0;
        auto total = [&](scalar r) { return d0 * (std::pow(r, n) - 1) / (r - 1); };
        while (total(hi) < L) hi *= 1.5;
        for (int it = 0; it < 100; ++it) {
            const scalar r = 0.5 * (lo + hi);
            (total(r) < L ? lo : hi) = r;
        }
        return 0.5 * (lo + hi);
    };
    BoxSpec b;
    b.n[0] = s.nTheta;
    b.n[1] = n;
    b.n[2] = s.nSpan;
    b.lo = {0, 0, 0};
    b.hi = {1, 1, s.span * D};
    b.periodic[0] = true;
    b.names[0] = "seam0";
    b.names[1] = "seam1";
    b.names[2] = "cylinder";
    b.types[2] = PatchType::Wall;
    b.names[3] = "farfield";
    b.twoD = s.twoD;
    b.periodic[2] = !s.twoD;
    b.names[4] = "front";
    b.names[5] = "back";
    b.map = [&](const Vec3& p) {
        // θ 顺时针（保证单元右手定向），t = 0 与 1 重合
        const scalar th = -2.0 * M_PI * p.x;
        const scalar c = std::cos(th), sn = std::sin(th);
        const scalar L = rayLength(c, sn);
        const scalar r = growth(L);
        const scalar j = p.y * n;
        const scalar f = r == 1.0 ? p.y : (std::pow(r, j) - 1) / (std::pow(r, n) - 1);
        const scalar rho = R + f * L;
        return Vec3{rho * c, rho * sn, p.z};
    };
    RawMesh m = generateBox(b);
    const scalar eps = 1e-6 * D;
    splitPatch(m, "farfield",
               {{"inlet", PatchType::Patch}, {"outlet", PatchType::Patch}, {"top", PatchType::Patch}, {"bottom", PatchType::Patch}},
               [&](const Vec3& c) {
                   // 取离面心最近的矩形边（跨角点的面归最近边）
                   const scalar d[4] = {std::abs(c.x - xIn), std::abs(c.x - xOut), std::abs(c.y - yHi), std::abs(c.y + yHi)};
                   (void)eps;
                   return int(std::min_element(d, d + 4) - d);
               });
    return m;
}

} // namespace cfd
