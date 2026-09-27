// RANS 槽道流（Re_b = U_b δ/ν = 6875，DNS Re_τ = 392）：k-ω SST、Spalart-Allmaras 的 Re_τ 误差 < 2%，
// 壁面距离与网格一致；并检查结果与进程数的一致性（可复现模式下由 ctest 的 np3 版本覆盖）。
#include "fvm/core/Json.h"
#include "fvm/mesh/Mesh.h"
#include "fvm/mesh/WallDistance.h"
#include "fvm/models/TurbulenceModel.h"
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

MeshPtr channelMesh() {
    BoxSpec s;
    s.n[0] = 4, s.n[1] = 120, s.n[2] = 1;
    s.hi = {1, 2, 0.1};
    s.periodic[0] = true;
    s.twoD = true;
    s.types[2] = s.types[3] = PatchType::Wall;
    s.stretch[1] = tanhStretch(3.0);
    RawMesh raw;
    if (par::master()) raw = generateBox(s);
    return Mesh::build(par::master() ? &raw : nullptr);
}

scalar reTau(MeshPtr mesh, const std::string& model) {
    const scalar nu = 1.0 / 6875;
    IncompressibleFlow flow(mesh, nu);
    flow.controls.steady = true;
    flow.controls.alphaU = 0.7;
    flow.controls.alphaP = 0.3;
    flow.controls.nCorr = 1;
    flow.controls.UControls = {"PBiCGStab", "DILU", 1e-10, 0.1, 1000, 0, 1};
    flow.controls.pControls = flow.controls.pFinalControls = {"PCG", "GAMG", 1e-10, 0.05, 2000, 0, 1};
    flow.divScheme = ConvectionScheme::parse("linearUpwind");
    flow.ddtScheme = DdtScheme::Steady;
    flow.forcing.mode = MomentumForcing::Mode::MeanVelocity;
    flow.forcing.Ubar = {1, 0, 0};
    flow.U().setUniform({1, 0, 0});
    flow.U().set<FixedValueBC<Vec3>>("yMin", Vec3{});
    flow.U().set<FixedValueBC<Vec3>>("yMax", Vec3{});
    flow.p().set<ZeroGradientBC<scalar>>("yMin");
    flow.p().set<ZeroGradientBC<scalar>>("yMax");
    flow.initialize();
    const Json d = Json::parse("{\"model\": \"" + model + "\"}");
    const Json b = Json::parse("{}");
    const Json i = Json::parse(R"({"k": 0.005, "omega": 50, "nuTilda": 5e-4})");
    auto turb = TurbulenceModel::create(flow, d, b, i);
    for (int it = 0; it < 6000; ++it) flow.step(1.0);
    const scalar rt = std::sqrt(flow.forcing.gradP) / nu;
    if (par::master()) std::cout << model << ": Re_tau " << rt << " (DNS 392.2)\n";
    return rt;
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    auto mesh = channelMesh();
    // 壁面距离：y = min(yc, 2 − yc)
    const auto& wd = wallDistance(*mesh);
    scalar err = 0;
    for (label c = 0; c < mesh->nTotalCells(); ++c) {
        const scalar yc = mesh->C()[c].y;
        err = std::max(err, std::abs(wd.y[c] - std::min(yc, 2 - yc)));
    }
    CHECK(par::allMax(err) < 1e-12);
    CHECK(std::abs(reTau(mesh, "kOmegaSST") / 392.2 - 1) < 0.02);
    CHECK(std::abs(reTau(mesh, "SpalartAllmaras") / 392.2 - 1) < 0.02);
    if (par::master()) std::cout << (fails ? "FAILED" : "ALL PASSED") << std::endl;
    return fails ? 1 : 0;
}
