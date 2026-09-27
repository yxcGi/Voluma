#pragma once
// 点集 k-d 树：最近点查询（单元中心定位、采样等）。

#include "fvm/core/Types.h"

#include <algorithm>
#include <vector>

namespace cfd {

class PointTree {
public:
    explicit PointTree(std::vector<Vec3> pts) : p_(std::move(pts)), idx_(p_.size()) {
        for (std::size_t i = 0; i < idx_.size(); ++i) idx_[i] = i;
        if (!p_.empty()) build(0, idx_.size());
    }
    // 最近点编号（空树返回 −1），d2 返回距离平方
    long nearest(const Vec3& x, scalar& d2) const {
        d2 = GREAT;
        long best = -1;
        if (!nodes_.empty()) query(0, x, best, d2);
        return best;
    }

private:
    struct Node {
        std::size_t b, e;
        int axis = -1;
        scalar split = 0;
        int left = -1, right = -1;
    };
    int build(std::size_t b, std::size_t e) {
        const int me = int(nodes_.size());
        nodes_.push_back({b, e});
        if (e - b <= 8) return me;
        Vec3 lo{GREAT, GREAT, GREAT}, hi{-GREAT, -GREAT, -GREAT};
        for (std::size_t i = b; i < e; ++i)
            for (int a = 0; a < 3; ++a) {
                lo[a] = std::min(lo[a], p_[idx_[i]][a]);
                hi[a] = std::max(hi[a], p_[idx_[i]][a]);
            }
        int axis = 0;
        for (int a = 1; a < 3; ++a)
            if (hi[a] - lo[a] > hi[axis] - lo[axis]) axis = a;
        const std::size_t mid = (b + e) / 2;
        std::nth_element(idx_.begin() + b, idx_.begin() + mid, idx_.begin() + e,
                         [&](std::size_t x, std::size_t y) { return p_[x][axis] < p_[y][axis]; });
        nodes_[me].axis = axis;
        nodes_[me].split = p_[idx_[mid]][axis];
        const int l = build(b, mid);
        const int r = build(mid, e);
        nodes_[me].left = l;
        nodes_[me].right = r;
        return me;
    }
    void query(int ni, const Vec3& x, long& best, scalar& d2) const {
        const Node& nd = nodes_[ni];
        if (nd.axis < 0) {
            for (std::size_t i = nd.b; i < nd.e; ++i) {
                const scalar d = magSqr(x - p_[idx_[i]]);
                if (d < d2 || (d == d2 && long(idx_[i]) < best)) {
                    d2 = d;
                    best = long(idx_[i]);
                }
            }
            return;
        }
        const scalar dx = x[nd.axis] - nd.split;
        query(dx < 0 ? nd.left : nd.right, x, best, d2);
        if (dx * dx <= d2) query(dx < 0 ? nd.right : nd.left, x, best, d2);
    }
    std::vector<Vec3> p_;
    std::vector<std::size_t> idx_;
    std::vector<Node> nodes_;
};

} // namespace cfd
