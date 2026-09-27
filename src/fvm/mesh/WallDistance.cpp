#include "fvm/mesh/WallDistance.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>

namespace cfd {

namespace {

struct WallFace {
    Vec3 c, n;
    scalar r;  // 面的等效半径（面积圆半径）
};

// 简单 k-d 树（最近 K 个点）
class KdTree {
public:
    explicit KdTree(const std::vector<WallFace>& f) : f_(f), idx_(f.size()) {
        for (std::size_t i = 0; i < idx_.size(); ++i) idx_[i] = i;
        nodes_.reserve(2 * f.size() / LEAF + 2);
        if (!f.empty()) build(0, idx_.size());
    }
    static constexpr int K = 8;
    // 返回最近的至多 K 个面编号
    std::vector<std::size_t> nearest(const Vec3& x) const {
        Best b;
        if (!nodes_.empty()) query(0, x, b);
        std::vector<std::size_t> out(b.id.begin(), b.id.begin() + b.n);
        return out;
    }

private:
    static constexpr std::size_t LEAF = 8;
    struct Node {
        std::size_t b, e;
        int axis = -1;
        scalar split = 0;
        int left = -1, right = -1;
    };
    struct Best {
        std::array<std::size_t, K> id{};
        std::array<scalar, K> d2{};
        int n = 0;
        scalar worst() const { return n < K ? GREAT : d2[n - 1]; }
        void insert(std::size_t i, scalar d) {
            if (n == K && d >= d2[K - 1]) return;
            int k = n < K ? n++ : K - 1;
            while (k > 0 && d2[k - 1] > d) {
                d2[k] = d2[k - 1];
                id[k] = id[k - 1];
                --k;
            }
            d2[k] = d;
            id[k] = i;
        }
    };
    int build(std::size_t b, std::size_t e) {
        const int me = int(nodes_.size());
        nodes_.push_back({b, e});
        if (e - b <= LEAF) return me;
        Vec3 lo{GREAT, GREAT, GREAT}, hi{-GREAT, -GREAT, -GREAT};
        for (std::size_t i = b; i < e; ++i)
            for (int a = 0; a < 3; ++a) {
                lo[a] = std::min(lo[a], f_[idx_[i]].c[a]);
                hi[a] = std::max(hi[a], f_[idx_[i]].c[a]);
            }
        int axis = 0;
        for (int a = 1; a < 3; ++a)
            if (hi[a] - lo[a] > hi[axis] - lo[axis]) axis = a;
        const std::size_t mid = (b + e) / 2;
        std::nth_element(idx_.begin() + b, idx_.begin() + mid, idx_.begin() + e,
                         [&](std::size_t x, std::size_t y) { return f_[x].c[axis] < f_[y].c[axis]; });
        nodes_[me].axis = axis;
        nodes_[me].split = f_[idx_[mid]].c[axis];
        const int l = build(b, mid);
        const int r = build(mid, e);
        nodes_[me].left = l;
        nodes_[me].right = r;
        return me;
    }
    void query(int ni, const Vec3& x, Best& best) const {
        const Node& nd = nodes_[ni];
        if (nd.axis < 0) {
            for (std::size_t i = nd.b; i < nd.e; ++i) best.insert(idx_[i], magSqr(x - f_[idx_[i]].c));
            return;
        }
        const scalar dx = x[nd.axis] - nd.split;
        const int first = dx < 0 ? nd.left : nd.right, second = dx < 0 ? nd.right : nd.left;
        query(first, x, best);
        if (dx * dx < best.worst()) query(second, x, best);
    }
    const std::vector<WallFace>& f_;
    std::vector<std::size_t> idx_;
    std::vector<Node> nodes_;
};

// 点到壁面面的近似距离（面近似为以 c 为中心、半径 r 的圆盘）
scalar faceDistance(const WallFace& w, const Vec3& x) {
    const Vec3 d = x - w.c;
    const scalar dn = dot(d, w.n);
    const scalar lat = mag(d - dn * w.n);
    const scalar out = std::max(lat - w.r, scalar(0));
    return std::sqrt(dn * dn + out * out);
}

WallDistance compute(const Mesh& m) {
    WallDistance wd;
    const label nt = m.nTotalCells();
    wd.y.assign(nt, GREAT);
    wd.n.assign(nt, Vec3{});
    wd.nearest.assign(nt, -1);
    // 二维：壁面面是沿 empty 方向拉伸的条带，等效半径取边长一半 = 面积 / (2·厚度)
    scalar depth = 0;
    if (m.twoD()) {
        const int ed = m.emptyDir();
        scalar lo = GREAT, hi = -GREAT;
        for (const auto& p : m.vtk().points) {
            lo = std::min(lo, p[ed]);
            hi = std::max(hi, p[ed]);
        }
        depth = par::allMax(hi) - par::allMin(lo);
    }
    // 收集所有壁面面到每个进程（按全局顺序，与分区无关）
    std::vector<WallFace> all;
    for (const auto& p : m.patches()) {
        if (p.type != PatchType::Wall) continue;
        std::vector<WallFace> loc(p.size);
        for (label i = 0; i < p.size; ++i) {
            const label f = p.start + i;
            const scalar A = m.magSf()[f];
            // 边界面法向朝外，取反使其指向流体
            loc[i] = {m.Cf()[f], -1.0 * m.Sf()[f] / A, m.twoD() ? 0.5 * A / depth : std::sqrt(A / PI)};
        }
        const auto g = p.ordering.gather(loc.data());
        glabel n = glabel(g.size());
        par::broadcast(&n, 1);
        std::vector<double> buf(std::size_t(n) * 7);
        if (par::master())
            for (glabel i = 0; i < n; ++i) {
                const auto& w = g[i];
                double* b = buf.data() + 7 * i;
                b[0] = w.c.x, b[1] = w.c.y, b[2] = w.c.z, b[3] = w.n.x, b[4] = w.n.y, b[5] = w.n.z, b[6] = w.r;
            }
        par::broadcast(buf.data(), int(buf.size()));
        for (glabel i = 0; i < n; ++i) {
            const double* b = buf.data() + 7 * i;
            all.push_back({{b[0], b[1], b[2]}, {b[3], b[4], b[5]}, b[6]});
        }
    }
    if (all.empty()) return wd;
    wd.hasWalls = true;
    const KdTree tree(all);
    for (label c = 0; c < nt; ++c) {
        const Vec3& x = m.C()[c];
        scalar best = GREAT;
        std::size_t bi = 0;
        for (auto i : tree.nearest(x)) {
            const scalar d = faceDistance(all[i], x);
            if (d < best || (d == best && i < bi)) {
                best = d;
                bi = i;
            }
        }
        wd.y[c] = best;
        wd.n[c] = all[bi].n;
        wd.nearest[c] = glabel(bi);
    }
    return wd;
}

} // namespace

std::vector<scalar> gatherWallFaces(const Mesh& m, const std::vector<scalar>& values) {
    std::vector<scalar> out;
    const label nI = m.nInternalFaces();
    for (const auto& p : m.patches()) {
        if (p.type != PatchType::Wall) continue;
        std::vector<double> loc(p.size);
        for (label i = 0; i < p.size; ++i) loc[i] = values[p.start - nI + i];
        auto g = p.ordering.gather(loc.data());
        glabel n = glabel(g.size());
        par::broadcast(&n, 1);
        g.resize(std::size_t(n));
        par::broadcast(g.data(), int(n));
        out.insert(out.end(), g.begin(), g.end());
    }
    return out;
}

const std::vector<glabel>& wallFaceGlobalIndex(const Mesh& m) {
    return m.cached<std::vector<glabel>>("wallFaceGlobalIndex", [&] {
        // 用全局面编号做键：收集全局壁面面的全局面号，再在本进程查表
        const label nI = m.nInternalFaces();
        std::vector<scalar> gid(m.nBoundaryFaces(), -1.0);
        for (label bf = 0; bf < m.nBoundaryFaces(); ++bf) gid[bf] = double(m.faceGlobal()[nI + bf]);
        const auto all = gatherWallFaces(m, gid);
        std::unordered_map<glabel, glabel> pos;
        for (std::size_t i = 0; i < all.size(); ++i) pos[glabel(all[i])] = glabel(i);
        std::vector<glabel> idx(m.nBoundaryFaces(), -1);
        for (const auto& p : m.patches()) {
            if (p.type != PatchType::Wall) continue;
            for (label f = p.start; f < p.end(); ++f) idx[f - nI] = pos.at(m.faceGlobal()[f]);
        }
        return idx;
    });
}

const WallDistance& wallDistance(const Mesh& m) {
    return m.cached<WallDistance>("wallDistance", [&] { return compute(m); });
}

} // namespace cfd
