#include "fvm/mesh/Extrude2D.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>

namespace cfd {

void mergeDuplicatePoints(Mesh2D& m, scalar tol) {
    // 按网格桶查找近邻
    const std::size_t n = m.points.size();
    std::unordered_map<std::int64_t, std::vector<glabel>> buckets;
    auto key = [&](scalar x, scalar y) {
        return (std::int64_t(std::floor(x / tol)) * 73856093) ^ (std::int64_t(std::floor(y / tol)) * 19349663);
    };
    std::vector<glabel> newId(n, -1);
    std::vector<Vec3> pts;
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3& p = m.points[i];
        glabel found = -1;
        const std::int64_t bx = std::int64_t(std::floor(p.x / tol)), by = std::int64_t(std::floor(p.y / tol));
        for (std::int64_t dx = -1; dx <= 1 && found < 0; ++dx)
            for (std::int64_t dy = -1; dy <= 1 && found < 0; ++dy) {
                auto it = buckets.find(((bx + dx) * 73856093) ^ ((by + dy) * 19349663));
                if (it == buckets.end()) continue;
                for (glabel q : it->second)
                    if (std::hypot(pts[q].x - p.x, pts[q].y - p.y) < tol) {
                        found = q;
                        break;
                    }
            }
        if (found < 0) {
            found = glabel(pts.size());
            pts.push_back(p);
            buckets[key(p.x, p.y)].push_back(found);
        }
        newId[i] = found;
    }
    for (auto& c : m.cells) {
        for (auto& p : c) p = newId[p];
        c.erase(std::unique(c.begin(), c.end()), c.end());
        while (c.size() > 1 && c.front() == c.back()) c.pop_back();
    }
    m.points = std::move(pts);
}

RawMesh extrude2D(const Mesh2D& m, const std::vector<BoundaryPatchSpec>& patches,
                  const std::function<int(glabel, glabel)>& classify, scalar depth) {
    const glabel nP = glabel(m.points.size());
    const glabel nC = glabel(m.cells.size());
    RawMesh r;
    r.points.resize(2 * nP);
    for (glabel i = 0; i < nP; ++i) {
        r.points[i] = {m.points[i].x, m.points[i].y, 0};
        r.points[nP + i] = {m.points[i].x, m.points[i].y, depth};
    }
    r.nCells = nC;
    // 检查单元为逆时针
    for (glabel c = 0; c < nC; ++c) {
        const auto& pl = m.cells[c];
        if (pl.size() < 3) throw std::runtime_error("extrude2D: degenerate cell");
        scalar a = 0;
        for (std::size_t k = 0; k < pl.size(); ++k) {
            const Vec3& p = m.points[pl[k]];
            const Vec3& q = m.points[pl[(k + 1) % pl.size()]];
            a += p.x * q.y - q.x * p.y;
        }
        if (a <= 0) throw std::runtime_error("extrude2D: cell " + std::to_string(c) + " is not counter-clockwise");
    }
    // 边 → (单元, a, b)
    struct EdgeUse {
        glabel cell, a, b;
    };
    std::map<std::pair<glabel, glabel>, std::vector<EdgeUse>> edges;
    for (glabel c = 0; c < nC; ++c) {
        const auto& pl = m.cells[c];
        for (std::size_t k = 0; k < pl.size(); ++k) {
            const glabel a = pl[k], b = pl[(k + 1) % pl.size()];
            edges[{std::min(a, b), std::max(a, b)}].push_back({c, a, b});
        }
    }
    struct Internal {
        glabel own, nei, a, b;  // (a,b) 为 owner 的逆时针走向
    };
    std::vector<Internal> internal;
    std::vector<std::vector<EdgeUse>> bnd(patches.size());
    for (auto& [k, uses] : edges) {
        if (uses.size() == 2) {
            const EdgeUse& u0 = uses[0].cell < uses[1].cell ? uses[0] : uses[1];
            const EdgeUse& u1 = uses[0].cell < uses[1].cell ? uses[1] : uses[0];
            internal.push_back({u0.cell, u1.cell, u0.a, u0.b});
        } else if (uses.size() == 1) {
            const int p = classify(uses[0].a, uses[0].b);
            if (p < 0 || p >= int(patches.size())) throw std::runtime_error("extrude2D: unclassified boundary edge");
            bnd[p].push_back(uses[0]);
        } else {
            throw std::runtime_error("extrude2D: non-manifold edge");
        }
    }
    std::sort(internal.begin(), internal.end(),
              [](const Internal& x, const Internal& y) { return std::tie(x.own, x.nei) < std::tie(y.own, y.nei); });
    for (const auto& f : internal) {
        r.addFace({f.a, f.b, nP + f.b, nP + f.a});
        r.owner.push_back(f.own);
        r.neighbour.push_back(f.nei);
    }
    for (std::size_t p = 0; p < patches.size(); ++p) {
        auto& list = bnd[p];
        std::sort(list.begin(), list.end(), [](const EdgeUse& x, const EdgeUse& y) { return x.cell < y.cell; });
        RawPatch rp;
        rp.name = patches[p].name;
        rp.type = patches[p].type;
        rp.start = r.nFaces();
        for (const auto& e : list) {
            r.addFace({e.a, e.b, nP + e.b, nP + e.a});
            r.owner.push_back(e.cell);
        }
        rp.size = r.nFaces() - rp.start;
        r.patches.push_back(rp);
    }
    RawPatch fb;
    fb.name = "frontAndBack";
    fb.type = PatchType::Empty;
    fb.start = r.nFaces();
    for (glabel c = 0; c < nC; ++c) {
        const auto& pl = m.cells[c];
        std::vector<glabel> f(pl.rbegin(), pl.rend());
        r.addFace(f);
        r.owner.push_back(c);
        std::vector<glabel> b;
        for (auto p : pl) b.push_back(nP + p);
        r.addFace(b);
        r.owner.push_back(c);
    }
    fb.size = r.nFaces() - fb.start;
    r.patches.push_back(fb);
    return r;
}

} // namespace cfd
