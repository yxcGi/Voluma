// 能量方程：(1) 静止流体、均匀热源、两壁等温的一维导热（解析抛物线，壁面半单元差分使误差为 Δ²/4），
//   壁面热流 q_w = −Q δ（热量流出）；(2) Boussinesq 静水平衡：线性稳定分层的封闭方腔速度保持为零
//   （体力进入 HbyA、壁面 fixedFluxPressure，线性体力与二次压力的离散平衡精确）；
//   (3) 热壁面律 T+(y+)：Jayatilleke 在热子层内为 Pr y+ 且在切换点连续，Kader 的近壁与对数区渐近行为。
#include "fvm/core/Json.h"
#include "fvm/mesh/Mesh.h"
#include "fvm/models/Energy.h"
#include "fvm/solvers/IncompressibleFlow.h"

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

MeshPtr boxMesh(int nx, int ny, bool periodicX) {
    BoxSpec s;
    s.n[0] = nx, s.n[1] = ny, s.n[2] = 1;
    s.hi = {1, 2, 0.1};
    s.periodic[0] = periodicX;
    s.twoD = true;
    if (!periodicX) s.types[0] = s.types[1] = PatchType::Wall;
    s.types[2] = s.types[3] = PatchType::Wall;
    RawMesh raw;
    if (par::master()) raw = generateBox(s);
    return Mesh::build(par::master() ? &raw : nullptr);
}

// (1) 导热 + 热源：T(y) = Q y (2 − y) / (2α)
void conductionWithSource() {
    auto mesh = boxMesh(4, 40, true);
    IncompressibleFlow flow(mesh, 1.0);
    flow.controls.steady = true;
    flow.controls.alphaU = 0.7;
    flow.controls.alphaP = 0.3;
    flow.controls.nCorr = 1;
    flow.ddtScheme = DdtScheme::Steady;
    flow.U().set<FixedValueBC<Vec3>>("yMin", Vec3{});
    flow.U().set<FixedValueBC<Vec3>>("yMax", Vec3{});
    flow.p().set<ZeroGradientBC<scalar>>("yMin");
    flow.p().set<ZeroGradientBC<scalar>>("yMax");
    flow.initialize();
    const Json d = Json::parse(R"({"alpha": 0.5, "source": 1.0, "relax": 1.0})");
    const Json b = Json::parse(R"({"yMin": {"T": {"type": "fixedValue", "value": 0}}, "yMax": {"T": {"type": "fixedValue", "value": 0}}})");
    const Json i = Json::parse(R"({"T": 0})");
    EnergyModel energy(flow, d, b, i);
    scalar res = 1;
    for (int it = 0; it < 400 && res > 1e-12; ++it) {
        flow.step(1.0);
        res = energy.correct();
    }
    scalar err = 0, umax = 0;
    for (label c = 0; c < mesh->nCells(); ++c) {
        const scalar y = mesh->C()[c].y;
        err = std::max(err, std::abs(energy.T()[c] - y * (2 - y) / (2 * 0.5)));
        umax = std::max(umax, mag(flow.U()[c]));
    }
    err = par::allMax(err);
    umax = par::allMax(umax);
    const auto hf = energy.wallHeatFlux();
    scalar qErr = 0;
    for (const auto& h : hf) qErr = std::max(qErr, std::abs(h.qMean + 1.0));  // 进入流体为正：q_w = −Q δ = −1
    if (par::master())
        std::cout << "conduction: max |T − exact| " << err << ", max |U| " << umax << ", wall q error " << qErr
                  << " (res " << res << ")\n";
    CHECK(err < 1e-3);  // 壁面半单元一侧差分：二阶，Δ²/4 = 6.25e-4
    CHECK(umax < 1e-10);
    CHECK(qErr < 1e-5);
}

// (2) Boussinesq 静水平衡：T = y/2（顶热底冷，稳定分层），g = (0, −1, 0)，β = 1，T_ref = 0.5 → f = (T − 0.5) ĵ
void hydrostaticBalance() {
    auto mesh = boxMesh(20, 40, false);
    IncompressibleFlow flow(mesh, 0.01);
    flow.controls.steady = false;
    flow.controls.nCorr = 2;
    flow.ddtScheme = DdtScheme::Backward;
    for (const char* p : {"xMin", "xMax", "yMin", "yMax"}) {
        flow.U().set<FixedValueBC<Vec3>>(p, Vec3{});
        flow.p().set<FixedFluxPressureBC>(p);
    }
    flow.initialize();
    const Json d = Json::parse(R"({"alpha": 0.01, "buoyancy": {"g": [0, -1, 0], "beta": 1.0, "TRef": 0.5}})");
    const Json b = Json::parse(R"({"yMin": {"T": {"type": "fixedValue", "value": 0}}, "yMax": {"T": {"type": "fixedValue", "value": 1}}})");
    const Json i = Json::parse(R"({"T": 0.5})");
    EnergyModel energy(flow, d, b, i);
    for (label c = 0; c < mesh->nTotalCells(); ++c) energy.T()[c] = 0.5 * mesh->C()[c].y;
    energy.T().correctBoundaryConditions();
    energy.updateBuoyancy();
    // 从精确静水压 p = y²/4 − y/2（∇p = f）与对应的壁面法向梯度起步：面上 φ_F 与 ∇p_f·S_f 逐面相等，
    // 预测步源项与速度修正都应精确为零（从 p = 0 起步时第一步 U* = dt f 留下的物理瞬态要靠粘性慢慢衰减）
    for (label c = 0; c < mesh->nTotalCells(); ++c) {
        const scalar y = mesh->C()[c].y;
        flow.p()[c] = y * y / 4 - y / 2;
    }
    for (label p = 0; p < label(mesh->patches().size()); ++p) {
        const Patch& pt = mesh->patches()[p];
        auto* bc = dynamic_cast<FixedGradientBC<scalar>*>(&flow.p().bc(p));
        if (!bc) continue;
        for (label i = 0; i < pt.size; ++i) {
            const label f = pt.start + i;
            const Vec3 n = mesh->Sf()[f] / mesh->magSf()[f];
            bc->gradient()[i] = dot(Vec3{0, 0.5 * mesh->Cf()[f].y - 0.5, 0}, n);
        }
    }
    flow.p().correctBoundaryConditions();
    scalar umax = 0, umax1 = 0;
    for (int n = 0; n < 20; ++n) {
        flow.step(0.1);
        energy.correct();
        umax = 0;
        for (label c = 0; c < mesh->nCells(); ++c) umax = std::max(umax, mag(flow.U()[c]));
        umax = par::allMax(umax);
        if (n == 0) umax1 = umax;
    }
    if (par::master()) std::cout << "hydrostatic: max |U| after step 1 " << umax1 << ", after 20 steps " << umax << '\n';
    CHECK(umax < 1e-8);
}

void wallLaws() {
    const scalar Pr = 0.71, Prt = 0.85;
    CHECK(std::abs(EnergyModel::TplusJayatilleke(5.0, Pr, Prt) - Pr * 5.0) < 1e-12);
    // 切换点连续：扫描找到 y+_T 两侧的值
    scalar jump = 0, yT = 0;
    for (scalar yp = 5; yp < 30; yp += 0.001) {
        const scalar a = EnergyModel::TplusJayatilleke(yp, Pr, Prt), b = EnergyModel::TplusJayatilleke(yp + 0.001, Pr, Prt);
        if (std::abs(b - a) > jump) {
            jump = std::abs(b - a);
            yT = yp;
        }
    }
    if (par::master()) std::cout << "Jayatilleke: largest step " << jump << " near y+ " << yT << '\n';
    CHECK(jump < 0.01);
    CHECK(std::abs(EnergyModel::TplusKader(0.5, Pr) / (Pr * 0.5) - 1) < 1e-3);
    const scalar beta = std::pow(3.85 * std::cbrt(Pr) - 1.3, 2) + 2.12 * std::log(Pr);
    CHECK(std::abs(EnergyModel::TplusKader(1000, Pr) / (2.12 * std::log(1001.0) + beta) - 1) < 1e-3);
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    conductionWithSource();
    hydrostaticBalance();
    wallLaws();
    if (par::master()) std::cout << (fails ? "FAILED" : "ALL PASSED") << std::endl;
    return fails ? 1 : 0;
}
