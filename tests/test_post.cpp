// 后处理与鲁棒性：(1) FieldAverage 的加权平均与二阶矩（两次样本解析可验）及续算读写；(2) CutPlane 的切割单元数与
//   线性场重构精确性（平面落在网格面上时只取一层）；(3) WallOutput 的 Couette 流壁面切应力 / Cf / Cp / y+ 与 CSV、VTK 表面输出；
//   (4) 湍流量有界化 bound()：负值 / NaN 单元取相邻有效值的平均并计数；(5) 网格质量量在正交盒子上为 0 / 0 / 1；
//   (6) 混合初始化：直通道势流投影给出均匀速度、Bernoulli 压力为 0，湍流量取入口值。
#include "fvm/core/Json.h"
#include "fvm/mesh/Mesh.h"
#include "fvm/models/TurbulenceModel.h"
#include "fvm/post/CutPlane.h"
#include "fvm/post/FieldAverage.h"
#include "fvm/post/WallOutput.h"
#include "fvm/solvers/IncompressibleFlow.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>

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

MeshPtr boxMesh(int nx, int ny, const Vec3& hi, bool periodicX, bool periodicY, bool wallsY) {
    BoxSpec s;
    s.n[0] = nx, s.n[1] = ny, s.n[2] = 1;
    s.hi = hi;
    s.periodic[0] = periodicX;
    s.periodic[1] = periodicY;
    s.twoD = true;
    if (wallsY) s.types[2] = s.types[3] = PatchType::Wall;
    RawMesh raw;
    if (par::master()) raw = generateBox(s);
    return Mesh::build(par::master() ? &raw : nullptr);
}

std::string tempDir(const std::string& name) {
    const std::string d = (std::filesystem::temp_directory_path() / ("voluma_test_post_" + name)).string();
    if (par::master()) {
        std::filesystem::remove_all(d);
        std::filesystem::create_directories(d);
    }
    par::barrier();
    return d;
}

int countLines(const std::string& file) {
    std::ifstream is(file);
    int n = 0;
    std::string line;
    while (std::getline(is, line))
        if (!line.empty()) ++n;
    return n;
}

// (1) 时间平均
void fieldAverageTest() {
    auto mesh = boxMesh(4, 4, {1, 1, 0.1}, false, false, true);
    const Mesh& m = *mesh;
    const label nT = m.nTotalCells();
    std::vector<scalar> s(nT);
    std::vector<Vec3> v(nT);
    FieldAverage fa(mesh);
    fa.add("s", s, true);
    fa.add("v", v, true);
    for (label c = 0; c < nT; ++c) {
        s[c] = m.C()[c].x;
        v[c] = {m.C()[c].x, m.C()[c].y, 0};
    }
    fa.sample(1.0);
    for (label c = 0; c < nT; ++c) {
        s[c] = 2 * m.C()[c].x + 1;
        v[c] = {2 * m.C()[c].x + 1, -m.C()[c].y, 1};
    }
    fa.sample(3.0);
    CHECK(std::abs(fa.averagedTime() - 4.0) < 1e-14);
    CHECK(fa.nSamples() == 2);
    scalar err = 0;
    const auto& is = fa.items()[0];
    const auto& iv = fa.items()[1];
    for (label c = 0; c < nT; ++c) {  // 含幽灵：平均场经过 halo 交换
        const scalar x = m.C()[c].x, y = m.C()[c].y;
        const scalar ms = (x + 3 * (2 * x + 1)) / 4, p2 = (x * x + 3 * (2 * x + 1) * (2 * x + 1)) / 4 - ms * ms;
        err = std::max(err, std::abs(is.meanS[c] - ms));
        err = std::max(err, std::abs(is.p2S[c] - p2));
        const scalar mx = ms, my = (y - 3 * y) / 4, mz = 0.75;
        err = std::max(err, std::abs(iv.meanV[c].x - mx) + std::abs(iv.meanV[c].y - my) + std::abs(iv.meanV[c].z - mz));
        const scalar vxy = (x * y + 3 * (2 * x + 1) * (-y)) / 4 - mx * my;
        err = std::max(err, std::abs(iv.p2V[c](0, 1) - vxy) + std::abs(iv.p2V[c](1, 0) - vxy));
        err = std::max(err, std::abs(iv.p2V[c](2, 2) - (0.75 - mz * mz)));
    }
    CHECK(par::allMax(err) < 1e-12);
    // 续算
    const std::string dir = tempDir("avg");
    fa.writeRestart(dir);
    FieldAverage fb(mesh);
    fb.add("s", s, true);
    fb.add("v", v, true);
    CHECK(fb.readRestart(dir));
    CHECK(std::abs(fb.averagedTime() - 4.0) < 1e-14 && fb.nSamples() == 2);
    scalar diff = 0;
    for (label c = 0; c < nT; ++c) {
        diff = std::max(diff, std::abs(fb.items()[0].meanS[c] - is.meanS[c]));
        diff = std::max(diff, std::abs(fb.items()[1].p2V[c](0, 1) - iv.p2V[c](0, 1)));
    }
    CHECK(par::allMax(diff) == 0.0);
    FieldAverage fc(mesh);
    fc.add("s", s, true);
    CHECK(!fc.readRestart(dir + "/none"));
}

// (2) 平面采样
void cutPlaneTest() {
    auto mesh = boxMesh(8, 8, {1, 1, 0.1}, false, true, false);
    ScalarField f(mesh, "phi");
    f.zeroGradient("xMin");
    f.zeroGradient("xMax");
    f.setFromFunction([](const Vec3& x) { return 2 * x.x; });
    f.correctBoundaryConditions();
    for (scalar xp : {0.5, 0.51}) {
        CutPlane pl(mesh, {xp, 0.3, 0.05}, {1, 0, 0});
        CHECK(pl.nCells() == 8);
        scalar err = 0;
        for (const auto& p : pl.points()) err = std::max(err, std::abs(p.x - xp));
        CHECK(par::allMax(err) < 1e-12);
        const std::string dir = tempDir("plane");
        const std::string file = dir + "/plane_phi.csv";
        pl.writeTable(file, f);
        par::barrier();
        if (par::master()) {
            CHECK(countLines(file) == 9);
            std::ifstream is(file);
            std::string line;
            std::getline(is, line);
            scalar e = 0;
            int n = 0;
            while (std::getline(is, line)) {
                const auto k = line.rfind(',');
                e = std::max(e, std::abs(std::stod(line.substr(k + 1)) - 2 * xp));
                ++n;
            }
            CHECK(n == 8 && e < 1e-12);
        }
        // 平面外：无单元
        CutPlane none(mesh, {1.5, 0, 0}, {1, 0, 0});
        CHECK(none.nCells() == 0);
    }
}

// (3) 壁面输出：Couette 流 u = y/2（上壁 U = 1 以 y = 2 运动），p = 2
void wallOutputTest() {
    auto mesh = boxMesh(4, 8, {1, 2, 0.1}, true, false, true);
    IncompressibleFlow flow(mesh, 0.1);
    flow.U().fixedValue("yMin", Vec3{0, 0, 0});
    flow.U().fixedValue("yMax", Vec3{1, 0, 0});
    flow.p().zeroGradient("yMin");
    flow.p().zeroGradient("yMax");
    flow.U().setFromFunction([](const Vec3& x) { return Vec3{0.5 * x.y, 0, 0}; });
    flow.p().setUniform(2.0);
    flow.initialize();
    WallOutput wo(flow, nullptr, nullptr, {});
    CHECK(wo.patches().size() == 2);
    wo.Uref = 2;
    wo.flowDir = {1, 0, 0};
    wo.update();
    const Mesh& m = *mesh;
    const label nI = m.nInternalFaces();
    scalar err = 0;
    for (label p : wo.patches()) {
        const Patch& pt = m.patches()[p];
        const bool bottom = pt.name == "yMin";
        for (label f = pt.start; f < pt.end(); ++f) {
            const label bf = f - nI;
            const scalar tauRef = bottom ? 0.05 : -0.05;  // ν du/dy，上壁相对运动反号
            err = std::max(err, std::abs(wo.tau()[bf].x - tauRef));
            err = std::max(err, std::abs(wo.Cf()[bf] - tauRef / 2.0));
            err = std::max(err, std::abs(wo.Cp()[bf] - 1.0));
            err = std::max(err, std::abs(wo.yPlus()[bf] - 0.125 * std::sqrt(0.05) / 0.1));
        }
    }
    CHECK(par::allMax(err) < 1e-12);
    const std::string dir = tempDir("wall");
    wo.writeCsv(dir);
    wo.writeVtk(dir, 0.0);
    wo.writeVtk(dir, 1.0);
    par::barrier();
    if (par::master()) {
        CHECK(countLines(dir + "/wall_yMin.csv") == 5);
        CHECK(countLines(dir + "/wall_yMax.csv") == 5);
        CHECK(std::filesystem::exists(dir + "/surface_yMin.pvd"));
        CHECK(std::filesystem::exists(dir + "/VTK/surface_yMax_00001.pvtp"));
        CHECK(std::filesystem::file_size(dir + "/VTK/surface_yMin_00000/0.vtp") > 200);
    }
}

// (4) 有界化
struct TestModel : TurbulenceModel {
    using TurbulenceModel::TurbulenceModel;
    std::string type() const override { return "test"; }
    void correctModel() override {}
    ScalarField& field(const std::string& n) { return addField(n, "zero", 1.0); }
    void doBound(ScalarField& f, scalar mn) { bound(f, mn); }
};

void boundTest() {
    auto mesh = boxMesh(4, 4, {1, 1, 0.1}, false, false, true);
    const Mesh& m = *mesh;
    IncompressibleFlow flow(mesh, 0.01);
    const Json dict = Json::parse(R"({"model": "test"})");
    const Json bnd = Json::parse(R"({"xMin": {"k": {"type": "fixedValue", "value": 1}}, "xMax": {}})");
    const Json ini = Json::parse("{}");
    TestModel tm(flow, dict, bnd, ini);
    ScalarField& k = tm.field("k");
    const std::set<glabel> bad{5, 6, 10};
    std::vector<scalar> orig(m.nTotalCells());
    for (label c = 0; c < m.nTotalCells(); ++c) {
        const glabel g = m.cellGlobal()[c];
        orig[c] = 1 + m.C()[c].x;
        k.internal()[c] = g == 5 ? -1.0 : (g == 6 ? -2.0 : (g == 10 ? std::nan("") : orig[c]));
    }
    tm.doBound(k, 1e-10);
    CHECK(tm.boundedCells().at("k") == 3);
    scalar err = 0;
    bool allBounded = true;
    for (label c = 0; c < m.nCells(); ++c) {
        if (!(k[c] >= 1e-10)) allBounded = false;
        if (!bad.count(m.cellGlobal()[c])) {
            err = std::max(err, std::abs(k[c] - orig[c]));
            continue;
        }
        // 期望：有效邻居（不在 bad 中）的平均
        scalar sum = 0;
        int n = 0;
        for (label kk = m.cellFaceOffsets()[c]; kk < m.cellFaceOffsets()[c + 1]; ++kk) {
            const label f = m.cellFaces()[kk];
            if (f >= m.nInternalFaces()) continue;
            const label o = m.owner()[f] == c ? m.neighbour()[f] : m.owner()[f];
            if (bad.count(m.cellGlobal()[o])) continue;
            sum += orig[o];
            ++n;
        }
        err = std::max(err, std::abs(k[c] - sum / n));
    }
    CHECK(allBounded);
    CHECK(par::allMax(err) < 1e-12);
    // 再次 bound 没有需要修正的单元
    tm.doBound(k, 1e-10);
    CHECK(tm.boundedCells().at("k") == 0);
}

// (5) 网格质量
void meshQualityTest() {
    auto mesh = boxMesh(4, 8, {1, 2, 0.1}, false, false, true);
    const auto no = mesh->cellNonOrthogonality(), sk = mesh->cellSkewness(), ar = mesh->cellAspectRatio();
    scalar e[3] = {0, 0, 0};
    for (label c = 0; c < mesh->nCells(); ++c) {
        e[0] = std::max(e[0], std::abs(no[c]));
        e[1] = std::max(e[1], std::abs(sk[c]));
        e[2] = std::max(e[2], std::abs(ar[c] - 1.0));
    }
    CHECK(par::allMax(e[0]) < 1e-9 && par::allMax(e[1]) < 1e-12 && par::allMax(e[2]) < 1e-12);
}

// (6) 混合初始化
void hybridInitTest() {
    auto mesh = boxMesh(8, 4, {2, 1, 0.1}, false, false, true);
    IncompressibleFlow flow(mesh, 0.01);
    flow.U().fixedValue("xMin", Vec3{1, 0, 0});
    flow.U().zeroGradient("xMax");
    flow.U().fixedValue("yMin", Vec3{0, 0, 0});
    flow.U().fixedValue("yMax", Vec3{0, 0, 0});
    flow.p().zeroGradient("xMin");
    flow.p().fixedValue("xMax", 0.0);
    flow.p().zeroGradient("yMin");
    flow.p().zeroGradient("yMax");
    flow.initialize();
    flow.projectDivergenceFree(2);
    const scalar H = flow.setBernoulliPressure();
    CHECK(std::abs(H - 0.5) < 1e-8);
    scalar eU = 0, eP = 0;
    for (label c = 0; c < mesh->nCells(); ++c) {
        eU = std::max(eU, mag(flow.U()[c] - Vec3{1, 0, 0}));
        eP = std::max(eP, std::abs(flow.p()[c]));
    }
    eU = par::allMax(eU);
    eP = par::allMax(eP);
    if (par::master()) std::cout << "hybrid init: H " << H << "  max|U - (1,0,0)| " << eU << "  max|p| " << eP << '\n';
    CHECK(eU < 1e-8);
    CHECK(eP < 1e-8);
    const Json dict = Json::parse(R"({"model": "kEpsilon"})");
    const Json bnd = Json::parse(R"({"xMin": {"turbulence": {"intensity": 0.05, "viscosityRatio": 10}}, "xMax": {}})");
    const Json ini = Json::parse("{}");
    auto turb = TurbulenceModel::create(flow, dict, bnd, ini);
    CHECK(bool(turb));
    turb->initializeFromInlets();
    scalar ek = 0, ee = 0;
    const scalar kRef = 1.5 * 0.05 * 0.05, eRef = 0.09 * kRef * kRef / (10 * 0.01);
    for (auto* f : turb->fields())
        for (label c = 0; c < mesh->nCells(); ++c) {
            if (f->name() == "k") ek = std::max(ek, std::abs((*f)[c] - kRef));
            if (f->name() == "epsilon") ee = std::max(ee, std::abs((*f)[c] - eRef) / eRef);
        }
    CHECK(par::allMax(ek) < 1e-12);
    CHECK(par::allMax(ee) < 1e-9);
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    fieldAverageTest();
    cutPlaneTest();
    wallOutputTest();
    boundTest();
    meshQualityTest();
    hybridInitTest();
    const int total = int(par::allMax(fails));
    if (par::master()) std::cout << (total ? "FAILED" : "all post-processing tests passed") << std::endl;
    return total ? 1 : 0;
}
