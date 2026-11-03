#include "fvm/mesh/PipeMesh.h"

#include "fvm/mesh/Extrude2D.h"

#include <cmath>
#include <stdexcept>

namespace cfd {

RawMesh generatePipeMesh(const PipeMeshSpec& s) {
    const int nc = s.nCore, nr = s.nRing;
    if (nc < 1 || nr < 1 || s.nAxial < 1) throw std::runtime_error("pipe: nCore, nRing and nAxial must be >= 1");
    if (!(s.coreFrac > 0 && s.coreFrac < 1)) throw std::runtime_error("pipe: coreFrac must be in (0, 1)");
    if (!(s.radius > 0 && s.length > 0)) throw std::runtime_error("pipe: radius and length must be > 0");
    const scalar R = s.radius, a = s.coreFrac * R;

    Mesh2D m;
    // 中央方块节点 (nc+1)²
    auto core = [&](int i, int j) { return glabel(j * (nc + 1) + i); };
    for (int j = 0; j <= nc; ++j)
        for (int i = 0; i <= nc; ++i) m.points.push_back({-a + 2 * a * i / nc, -a + 2 * a * j / nc, 0});
    // 方块周界节点：逆时针，从 (−a,−a) 出发，底边 → 右边 → 顶边 → 左边
    const int np = 4 * nc;
    auto perim = [&](int k) -> glabel {
        k = ((k % np) + np) % np;
        if (k < nc) return core(k, 0);
        if (k < 2 * nc) return core(nc, k - nc);
        if (k < 3 * nc) return core(3 * nc - k, nc);
        return core(0, 4 * nc - k);
    };
    // 环向节点：ring[r][k]，r = 0 为方块周界，r = nr 为圆周；角度沿周界均匀，从 −135° 起
    std::vector<glabel> ring(std::size_t(nr + 1) * np);
    for (int k = 0; k < np; ++k) ring[k] = perim(k);
    for (int r = 1; r <= nr; ++r) {
        scalar t = scalar(r) / nr;
        if (s.stretch > 0) t = 1 - std::tanh(s.stretch * (1 - t)) / std::tanh(s.stretch);
        for (int k = 0; k < np; ++k) {
            const Vec3& q = m.points[perim(k)];
            const scalar th = -0.75 * M_PI + 0.5 * M_PI * k / nc;
            ring[r * np + k] = glabel(m.points.size());
            m.points.push_back({(1 - t) * q.x + t * R * std::cos(th), (1 - t) * q.y + t * R * std::sin(th), 0});
        }
    }
    // 单元（逆时针）：方块，然后由内向外各环
    for (int j = 0; j < nc; ++j)
        for (int i = 0; i < nc; ++i) m.cells.push_back({core(i, j), core(i + 1, j), core(i + 1, j + 1), core(i, j + 1)});
    for (int r = 0; r < nr; ++r)
        for (int k = 0; k < np; ++k) {
            const int k1 = (k + 1) % np;
            m.cells.push_back({ring[r * np + k], ring[(r + 1) * np + k], ring[(r + 1) * np + k1], ring[r * np + k1]});
        }
    // 壁面边：两端点都在最外环（最外环节点是最后追加的 np 个点）
    const glabel outer0 = ring[std::size_t(nr) * np];
    auto classify = [outer0](glabel p, glabel q) { return (p >= outer0 && q >= outer0) ? 0 : -1; };
    const std::vector<BoundaryPatchSpec> patches{{"wall", PatchType::Wall}};
    return extrudeLayers(m, patches, classify, s.length, s.nAxial, s.periodic, s.periodic ? "zMin" : "inlet",
                         s.periodic ? "zMax" : "outlet");
}

} // namespace cfd
