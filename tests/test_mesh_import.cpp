// 网格导入：Gmsh（2.2 二维混合、4.1 三维四面体/六面体/金字塔混合、二阶四面体）、Fluent（ASCII / 二进制）、
// 二进制 polyMesh。检查体积、单元闭合、patch 面积，并要求 Fluent/二进制 polyMesh 与原 ASCII polyMesh 完全相同。
#include "fvm/mesh/Mesh.h"
#include "fvm/mesh/MeshReaders.h"

#include <cmath>
#include <iostream>

using namespace cfd;

static int fails = 0;
#define CHECK(c)                                                                                                       \
    do {                                                                                                               \
        if (!(c)) {                                                                                                    \
            std::cerr << "FAIL " #c " at line " << __LINE__ << std::endl;                                              \
            ++fails;                                                                                                   \
        }                                                                                                              \
    } while (0)

namespace {

const std::string dir = "meshes/";  // 测试在源码根目录运行

// 原始网格上：每个单元 ΣS_f = 0，边界面朝外（面中心相对单元中心）
scalar maxOpenness(const RawMesh& m) {
    std::vector<Vec3> sum(m.nCells);
    std::vector<scalar> area(m.nCells, 0);
    for (glabel f = 0; f < m.nFaces(); ++f) {
        const glabel b = m.faceOffsets[f], e = m.faceOffsets[f + 1];
        Vec3 c;
        for (glabel k = b; k < e; ++k) c += m.points[m.facePoints[k]];
        c /= scalar(e - b);
        Vec3 S;
        for (glabel k = b; k < e; ++k)
            S += 0.5 * cross(m.points[m.facePoints[k]] - c, m.points[m.facePoints[k + 1 < e ? k + 1 : b]] - c);
        sum[m.owner[f]] += S;
        area[m.owner[f]] += mag(S);
        if (f < m.nInternalFaces()) {
            sum[m.neighbour[f]] -= S;
            area[m.neighbour[f]] += mag(S);
        }
    }
    scalar r = 0;
    for (glabel c = 0; c < m.nCells; ++c) r = std::max(r, mag(sum[c]) / area[c]);
    return r;
}

scalar patchArea(const Mesh& m, const std::string& name) {
    const label p = m.findPatch(name);
    if (p < 0) return -1;
    scalar a = 0;
    const auto& pt = m.patches()[p];
    for (label f = pt.start; f < pt.end(); ++f) a += m.magSf()[f];
    par::allSumInPlace(&a, 1);
    return a;
}

void checkCube(const std::string& file, bool twoD, glabel nCellsExpected, const MeshImportOptions& o = {}) {
    RawMesh raw;
    if (par::master()) {
        raw = readMeshFile(dir + file, o);
        const scalar open = maxOpenness(raw);
        std::cout << file << ": " << raw.nCells << " cells, max openness " << open << '\n';
        CHECK(open < 1e-12);
        if (nCellsExpected > 0) CHECK(raw.nCells == nCellsExpected);
    }
    auto m = Mesh::build(par::master() ? &raw : nullptr);
    const scalar depth = twoD ? o.depth2D : 1.0;
    CHECK(m->twoD() == twoD);
    CHECK(std::abs(m->totalVolume() - depth) < 1e-12);
    // movingWall 为 y=1（二维）或 z=1（三维）的单位面；fixedWalls 为其余边界
    CHECK(std::abs(patchArea(*m, "movingWall") - depth) < 1e-12);
    CHECK(std::abs(patchArea(*m, "fixedWalls") - (twoD ? 3 * depth : 5.0)) < 1e-12);
    scalar vmin = GREAT;
    for (label c = 0; c < m->nCells(); ++c) vmin = std::min(vmin, m->V()[c]);
    CHECK(par::allMin(vmin) > 0);
}

bool sameMesh(const RawMesh& a, const RawMesh& b) {
    if (a.nCells != b.nCells || a.points.size() != b.points.size() || a.owner != b.owner || a.neighbour != b.neighbour)
        return false;
    if (a.faceOffsets != b.faceOffsets || a.facePoints != b.facePoints || a.patches.size() != b.patches.size()) return false;
    for (std::size_t i = 0; i < a.points.size(); ++i)
        if (mag(a.points[i] - b.points[i]) != 0) return false;
    for (std::size_t p = 0; p < a.patches.size(); ++p)
        if (a.patches[p].name != b.patches[p].name || a.patches[p].type != b.patches[p].type ||
            a.patches[p].start != b.patches[p].start || a.patches[p].size != b.patches[p].size)
            return false;
    return true;
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    checkCube("import/mixed2D_v22.msh", true, 197);
    checkCube("import/hybrid3D_v41.msh", false, 813);
    checkCube("import/tet10_v41.msh", false, 100);

    // Fluent 与二进制 polyMesh：与原 ASCII polyMesh 逐项相同（frontAndBack 在 Fluent 中是 wall，改回 empty）
    if (par::master()) {
        const RawMesh ref = readPolyMesh(dir + "cavity/polyMesh");
        MeshImportOptions o;
        o.patchTypes["frontAndBack"] = PatchType::Empty;
        CHECK(sameMesh(ref, readFluent(dir + "import/cavity_fluent.msh", o)));
        CHECK(sameMesh(ref, readFluent(dir + "import/cavity_fluent_bin.msh", o)));
        CHECK(sameMesh(ref, readPolyMesh(dir + "import/cavity_binary_polyMesh")));
    }
    par::barrier();
    if (par::master()) std::cout << (fails ? "FAILED" : "ALL PASSED") << std::endl;
    return fails ? 1 : 0;
}
