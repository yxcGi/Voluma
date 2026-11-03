// 边界条件单元测试：
//   1. mixed / inletOutlet / pressureInletOutletVelocity 的系数与 evaluate 自洽：
//      φ_b = vic∘φ_P + vbc，Δ(φ_b − φ_P) = gic∘φ_P + gbc
//   2. inletOutlet 按通量符号切换；totalPressure 流入面 p0 − ½|U_b|²；flowRateInlet 面积与速度；rotatingWall 切向速度
//   3. TimeTable / ProfileTable 插值；抛物线入口面值与公式一致；ramp 缓启动
//   4. 湍流入口（强度 + 长度尺度 / 粘性比）面值与 Fluent 公式一致，ν̃ 由粘性比反解
//   5. advective 出口：稳态退化为零梯度；纯对流高斯脉冲全部流出且不反射
#include "fvm/app/CaseSetup.h"
#include "fvm/field/Functions.h"
#include "fvm/models/TurbulenceBCs.h"
#include "fvm/models/TurbulenceModel.h"
#include "fvm/solvers/IncompressibleFlow.h"
#include "fvm/solvers/ScalarTransport.h"

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

MeshPtr box2D(int nx, int ny, Vec3 hi, bool walls) {
    RawMesh raw;
    if (par::master()) {
        BoxSpec s;
        s.n[0] = nx, s.n[1] = ny, s.n[2] = 1;
        s.hi = hi;
        s.twoD = true;
        if (walls) s.types[2] = s.types[3] = PatchType::Wall;
        raw = generateBox(s);
    }
    return Mesh::build(par::master() ? &raw : nullptr);
}

template <class T> scalar maxErr(const T& a, const T& b) { return Traits<T>::magnitude(a - b); }

// 系数自洽性：对 patch p 的每个面检查 evaluate 结果与 vic/vbc、gic/gbc 一致
template <class T> scalar coeffConsistency(const VolField<T>& f, label p) {
    const auto& bc = f.bc(p);
    const Mesh& m = f.mesh();
    scalar err = 0;
    for (label i = 0; i < m.patches()[p].size; ++i) {
        const T phiP = f[bc.faceCell(i)];
        const T b = f.patchValue(p, i);
        const T bv = cmptMul(bc.valueInternalCoeff(f, i), phiP) + bc.valueBoundaryCoeff(f, i);
        const T gv = cmptMul(bc.gradientInternalCoeff(f, i), phiP) + bc.gradientBoundaryCoeff(f, i);
        err = std::max(err, maxErr(b, bv));
        err = std::max(err, maxErr(bc.deltaCoeff(i) * (b - phiP), gv));
    }
    return par::allMax(err);
}

void testCoefficientsAndSwitching() {
    auto mesh = box2D(8, 8, {1, 1, 0.1}, true);
    IncompressibleFlow flow(mesh, 0.01);
    auto& U = flow.U();
    auto& p = flow.p();
    const label pIn = mesh->findPatch("xMin"), pOut = mesh->findPatch("xMax");
    // 非均匀速度场：u = 1 + y，v = 0.3 sin(πx)
    U.setFromFunction([](const Vec3& x) { return Vec3{1 + x.y, 0.3 * std::sin(PI * x.x), 0}; });
    U.set<FixedValueBC<Vec3>>("xMin", Vec3{1, 0, 0});
    U.set<InletOutletBC<Vec3>>("xMax", &flow.context(), Vec3{0.5, 0, 0});
    U.set<FixedValueBC<Vec3>>("yMin", Vec3{});
    U.set<FixedValueBC<Vec3>>("yMax", Vec3{});
    p.set<ZeroGradientBC<scalar>>("xMin");
    p.set<FixedValueBC<scalar>>("xMax", 0.0);
    p.set<ZeroGradientBC<scalar>>("yMin");
    p.set<ZeroGradientBC<scalar>>("yMax");
    flow.initialize();  // φ = U·S，出口为流出
    U.updateBCs(0.0);
    U.correctBoundaryConditions();
    CHECK(coeffConsistency(U, pOut) < 1e-12);
    for (label i = 0; i < mesh->patches()[pOut].size; ++i) CHECK(maxErr(U.patchValue(pOut, i), U[U.bc(pOut).faceCell(i)]) < 1e-14);
    CHECK(!U.bc(pOut).fixesValue());
    // 反向流动：出口变为流入，面值取 inletValue
    for (label c = 0; c < mesh->nTotalCells(); ++c) U[c] = {-1, 0, 0};
    flow.initialize();
    U.updateBCs(1.0);
    U.correctBoundaryConditions();
    CHECK(coeffConsistency(U, pOut) < 1e-12);
    for (label i = 0; i < mesh->patches()[pOut].size; ++i) CHECK(maxErr(U.patchValue(pOut, i), Vec3{0.5, 0, 0}) < 1e-14);
    CHECK(par::allMax(U.bc(pOut).fixesValue() ? 1.0 : 0.0) > 0.5);  // 并行时某进程可能没有该 patch 的面
    CHECK(!U.bc(pOut).fixesFlux());

    // mixed：w = 0.3，refValue 2，refGrad 1
    p.set<MixedBC<scalar>>("xMax", 2.0, 1.0, 0.3);
    p.setFromFunction([](const Vec3& x) { return x.x * x.x + x.y; });
    p.correctBoundaryConditions();
    CHECK(coeffConsistency(p, pOut) < 1e-12);
    {
        const auto& bc = p.bc(pOut);
        for (label i = 0; i < mesh->patches()[pOut].size; ++i) {
            const scalar pb = 0.3 * 2.0 + 0.7 * (p[bc.faceCell(i)] + 1.0 / bc.deltaCoeff(i));
            CHECK(std::abs(p.patchValue(pOut, i) - pb) < 1e-12);
        }
    }

    // pressureInletOutletVelocity + totalPressure：入口流入，p_b = p0 − ½|U_b|²，U_b 法向外推、切向为 0
    U.setFromFunction([](const Vec3& x) { return Vec3{1 + x.y, 0.3, 0}; });
    U.set<PressureInletOutletVelocityBC>("xMin", &flow.context());
    p.set<TotalPressureBC>("xMin", &flow.context(), [](scalar) { return 2.0; });
    flow.initialize();
    U.updateBCs(0.0);
    U.correctBoundaryConditions();
    p.updateBCs(0.0);
    p.correctBoundaryConditions();
    CHECK(coeffConsistency(U, pIn) < 1e-12);
    {
        const auto& bc = U.bc(pIn);
        for (label i = 0; i < mesh->patches()[pIn].size; ++i) {
            const Vec3 ub = U.patchValue(pIn, i);
            const Vec3 uP = U[bc.faceCell(i)];
            CHECK(flow.phi()[bc.face(i)] < 0);
            CHECK(std::abs(ub.x - uP.x) < 1e-14 && std::abs(ub.y) < 1e-14);
            CHECK(std::abs(p.patchValue(pIn, i) - (2.0 - 0.5 * magSqr(ub))) < 1e-12);
        }
    }
    // 流出时 totalPressure 给 p0，速度零梯度
    for (label c = 0; c < mesh->nTotalCells(); ++c) U[c] = {-1, 0.2, 0};
    flow.initialize();
    U.updateBCs(0.0);
    U.correctBoundaryConditions();
    p.updateBCs(0.0);
    p.correctBoundaryConditions();
    CHECK(coeffConsistency(U, pIn) < 1e-12);
    for (label i = 0; i < mesh->patches()[pIn].size; ++i) {
        CHECK(std::abs(p.patchValue(pIn, i) - 2.0) < 1e-12);
        CHECK(maxErr(U.patchValue(pIn, i), Vec3{-1, 0.2, 0}) < 1e-14);
    }

    // flowRateInlet：面积 = 1 × 0.1，Q = 0.05 → |U| = 0.5，方向 +x（内法向）
    U.set<FlowRateInletBC>("xMin", [](scalar t) { return 0.05 * (1 + t); });
    U.updateBCs(1.0);
    U.correctBoundaryConditions();
    CHECK(std::abs(static_cast<FlowRateInletBC&>(U.bc(pIn)).area() - 0.1) < 1e-12);
    for (label i = 0; i < mesh->patches()[pIn].size; ++i) CHECK(maxErr(U.patchValue(pIn, i), Vec3{1.0, 0, 0}) < 1e-12);

    // freestream：来流 (1, 0, 0) 时 xMin 全流入（w=1）、xMax 全流出（w=0）；斜来流 (1, 0.5, 0) 时 yMax 上 w = ½ − ½ (U∞·n)/|U∞|
    U.set<FreestreamVelocityBC>("xMin", Vec3{1, 0, 0});
    U.set<FreestreamVelocityBC>("xMax", Vec3{1, 0, 0});
    U.set<FreestreamVelocityBC>("yMax", Vec3{1, 0.5, 0});
    p.set<FreestreamPressureBC>("xMax", &flow.context(), 3.0, nullptr);
    U.setFromFunction([](const Vec3& x) { return Vec3{1 + x.y, 0.3, 0}; });
    U.correctBoundaryConditions();
    p.updateBCs(0.0);
    p.correctBoundaryConditions();
    for (label i = 0; i < mesh->patches()[pIn].size; ++i) CHECK(maxErr(U.patchValue(pIn, i), Vec3{1, 0, 0}) < 1e-14);
    for (label i = 0; i < mesh->patches()[pOut].size; ++i) {
        CHECK(maxErr(U.patchValue(pOut, i), U[U.bc(pOut).faceCell(i)]) < 1e-14);
        CHECK(std::abs(p.patchValue(pOut, i) - 3.0) < 1e-12);  // 流出面定压
    }
    {
        const label pTop = mesh->findPatch("yMax");
        auto& bc = static_cast<MixedBC<Vec3>&>(U.bc(pTop));
        const scalar w = 0.5 - 0.5 * 0.5 / std::sqrt(1.25);
        for (label i = 0; i < mesh->patches()[pTop].size; ++i) CHECK(std::abs(bc.valueFraction()[i] - w) < 1e-12);
        CHECK(coeffConsistency(U, pTop) < 1e-12);
    }
}

void testFunctionsAndParsing() {
    auto mesh = box2D(8, 8, {1, 1, 0.1}, true);
    IncompressibleFlow flow(mesh, 0.01);
    auto& U = flow.U();
    for (const auto& pt : mesh->patches()) U.set<ZeroGradientBC<Vec3>>(pt.name);
    const label pIn = mesh->findPatch("xMin");
    const label pTop = mesh->findPatch("yMax");

    // 时间表与剖面表
    {
        TimeTable<scalar> tt = TimeTable<scalar>::parse(Json::parse("[[0, 0], [1, 2], [3, 2]]"), ".");
        CHECK(std::abs(tt(0.5) - 1.0) < 1e-14 && std::abs(tt(-1) - 0.0) < 1e-14 && std::abs(tt(10) - 2.0) < 1e-14);
        TimeTable<Vec3> tv = TimeTable<Vec3>::parse(Json::parse("[[0, [0, 0, 0]], [2, [2, 4, 0]]]"), ".");
        CHECK(maxErr(tv(1.0), Vec3{1, 2, 0}) < 1e-14);
        TimeTable<Vec3> tc = TimeTable<Vec3>::parse(Json::parse("[1, 2, 3]"), ".");
        CHECK(tc.constant() && maxErr(tc(5.0), Vec3{1, 2, 3}) < 1e-14);
        ProfileTable<scalar> pr =
            ProfileTable<scalar>::parse(Json::parse(R"({"axis": "y", "points": [[0, 0], [0.5, 1], [1, 0]], "scale": 2})"), ".");
        CHECK(std::abs(pr({0.3, 0.25, 0}) - 1.0) < 1e-14 && std::abs(pr({0, 2, 0}) - 0.0) < 1e-14);
    }
    // 抛物线入口（槽道）：u = 1.5 Ub (1 − 4 (y − 0.5)²)，方向为内法向 +x；带 ramp 时乘 smoothstep
    setBoundaryCondition(U, "xMin",
                         Json::parse(R"({"type": "fixedValue", "parabolic": {"Ubulk": 2, "halfWidth": 0.5, "axis": [0, 1, 0]},
                                         "ramp": {"duration": 2}})"),
                         &flow.context());
    U.updateBCs(1.0);  // smoothstep(0.5) = 0.5
    U.correctBoundaryConditions();
    for (label i = 0; i < mesh->patches()[pIn].size; ++i) {
        const scalar y = mesh->Cf()[U.bc(pIn).face(i)].y;
        const Vec3 ex{0.5 * 3.0 * (1 - 4 * (y - 0.5) * (y - 0.5)), 0, 0};
        CHECK(maxErr(U.patchValue(pIn, i), ex) < 1e-12);
    }
    // 圆管抛物线：以 yMax 为“入口”，半径 0.5，中心 (0.5, 1, 0.05)：u = −y 方向 2 Ub (1 − r²/R²)
    setBoundaryCondition(U, "yMax",
                         Json::parse(R"({"type": "fixedValue", "parabolic": {"Ubulk": 1, "radius": 0.5, "center": [0.5, 1, 0.05]}})"),
                         &flow.context());
    U.correctBoundaryConditions();
    for (label i = 0; i < mesh->patches()[pTop].size; ++i) {
        const Vec3 x = mesh->Cf()[U.bc(pTop).face(i)];
        const scalar r2 = (x.x - 0.5) * (x.x - 0.5) + (x.z - 0.05) * (x.z - 0.05);
        CHECK(maxErr(U.patchValue(pTop, i), Vec3{0, -2.0 * (1 - r2 / 0.25), 0}) < 1e-12);
    }
    // 旋转壁：绕 z 轴、原点 (0.5, 0.5, 0)、ω = 2 → yMax 上 U = 2 (−(y − 0.5), x − 0.5, 0) 去法向分量 → (−1, 0, 0)
    setBoundaryCondition(U, "yMax", Json::parse(R"({"type": "rotatingWall", "origin": [0.5, 0.5, 0], "axis": [0, 0, 1], "omega": 2})"),
                         &flow.context());
    U.correctBoundaryConditions();
    for (label i = 0; i < mesh->patches()[pTop].size; ++i) CHECK(maxErr(U.patchValue(pTop, i), Vec3{-1, 0, 0}) < 1e-12);
    // 时间表 fixedValue
    setBoundaryCondition(U, "xMin", Json::parse(R"({"type": "fixedValue", "timeTable": [[0, [0, 0, 0]], [1, [2, 0, 0]]]})"),
                         &flow.context());
    U.updateBCs(0.25);
    U.correctBoundaryConditions();
    for (label i = 0; i < mesh->patches()[pIn].size; ++i) CHECK(maxErr(U.patchValue(pIn, i), Vec3{0.5, 0, 0}) < 1e-12);
    // 流量入口（JSON）
    setBoundaryCondition(U, "xMin", Json::parse(R"({"type": "flowRateInlet", "massFlowRate": 0.2, "rho": 2})"), &flow.context());
    U.updateBCs(0.0);
    U.correctBoundaryConditions();
    for (label i = 0; i < mesh->patches()[pIn].size; ++i) CHECK(maxErr(U.patchValue(pIn, i), Vec3{1.0, 0, 0}) < 1e-12);
}

void testTurbulentInlet() {
    auto mesh = box2D(8, 8, {1, 1, 0.1}, true);
    const scalar nu = 1e-4;
    IncompressibleFlow flow(mesh, nu);
    auto& U = flow.U();
    auto& p = flow.p();
    U.setUniform({1, 0, 0});
    U.set<FixedValueBC<Vec3>>("xMin", Vec3{1, 0, 0});
    U.set<InletOutletBC<Vec3>>("xMax", &flow.context(), Vec3{});
    U.set<FixedValueBC<Vec3>>("yMin", Vec3{});
    U.set<FixedValueBC<Vec3>>("yMax", Vec3{});
    p.set<ZeroGradientBC<scalar>>("xMin");
    p.set<FixedValueBC<scalar>>("xMax", 0.0);
    p.set<ZeroGradientBC<scalar>>("yMin");
    p.set<ZeroGradientBC<scalar>>("yMax");
    flow.initialize();
    const label pIn = mesh->findPatch("xMin");
    const scalar I = 0.05, L = 0.07, Cmu = 0.09;
    const scalar kIn = 1.5 * I * I;
    // k-ε：强度 + 长度尺度
    {
        const Json d = Json::parse(R"({"model": "kEpsilon"})");
        const Json b = Json::parse(R"({"xMin": {"turbulence": {"intensity": 0.05, "lengthScale": 0.07}},
                                      "xMax": {"k": {"type": "inletOutlet", "inletValue": 0}, "epsilon": {"type": "inletOutlet", "inletValue": 0}}})");
        const Json i = Json::parse(R"({"k": 0.001, "epsilon": 0.001})");
        auto turb = TurbulenceModel::create(flow, d, b, i);
        const ScalarField* k = nullptr;
        const ScalarField* eps = nullptr;
        for (auto* f : turb->fields()) {
            if (f->name() == "k") k = f;
            if (f->name() == "epsilon") eps = f;
        }
        CHECK(k && eps);
        const scalar epsIn = std::pow(Cmu, 0.75) * std::pow(kIn, 1.5) / L;
        for (label n = 0; n < mesh->patches()[pIn].size; ++n) {
            CHECK(std::abs(k->patchValue(pIn, n) / kIn - 1) < 1e-12);
            CHECK(std::abs(eps->patchValue(pIn, n) / epsIn - 1) < 1e-12);
        }
    }
    // k-ω SST：强度 + 粘性比
    {
        const Json d = Json::parse(R"({"model": "kOmegaSST"})");
        const Json b = Json::parse(R"({"xMin": {"turbulence": {"intensity": 0.05, "viscosityRatio": 10}},
                                      "xMax": {"k": {"type": "inletOutlet", "inletValue": 0}, "omega": {"type": "inletOutlet", "inletValue": 0}}})");
        const Json i = Json::parse(R"({"k": 0.001, "omega": 1})");
        auto turb = TurbulenceModel::create(flow, d, b, i);
        const ScalarField* om = nullptr;
        for (auto* f : turb->fields())
            if (f->name() == "omega") om = f;
        CHECK(om);
        const scalar omIn = kIn / (10 * nu);
        for (label n = 0; n < mesh->patches()[pIn].size; ++n) CHECK(std::abs(om->patchValue(pIn, n) / omIn - 1) < 1e-12);
    }
    // Spalart-Allmaras：粘性比 5 → χ f_v1(χ) = 5
    {
        const Json d = Json::parse(R"({"model": "SpalartAllmaras"})");
        const Json b = Json::parse(R"({"xMin": {"nuTilda": {"type": "viscosityRatio", "viscosityRatio": 5}},
                                      "xMax": {"nuTilda": {"type": "inletOutlet", "inletValue": 0}}})");
        const Json i = Json::parse(R"({"nuTilda": 1e-4})");
        auto turb = TurbulenceModel::create(flow, d, b, i);
        const ScalarField* nt = nullptr;
        for (auto* f : turb->fields())
            if (f->name() == "nuTilda") nt = f;
        CHECK(nt);
        for (label n = 0; n < mesh->patches()[pIn].size; ++n) {
            const scalar chi = nt->patchValue(pIn, n) / nu;
            const scalar fv1 = chi * chi * chi / (chi * chi * chi + 7.1 * 7.1 * 7.1);
            CHECK(std::abs(chi * fv1 - 5) < 1e-9);
        }
    }
}

void testAdvective() {
    // 一维纯对流（U = 1）：高斯脉冲从出口流出，advective 出口不反射；稳态时退化为零梯度
    auto mesh = box2D(100, 1, {1, 0.01, 0.01}, false);
    ScalarTransport st(mesh, "T", 0.0);
    std::vector<scalar> phi(mesh->nFaces());
    for (label f = 0; f < mesh->nFaces(); ++f) phi[f] = dot(Vec3{1, 0, 0}, mesh->Sf()[f]);
    BCContext ctx;
    ctx.phi = &phi;
    ctx.time = &st.time();
    ctx.ddt = DdtScheme::Backward;
    st.ddtScheme = DdtScheme::Backward;
    st.divScheme = ConvectionScheme::parse("linear");
    st.T().set<FixedValueBC<scalar>>("xMin", 0.0);
    st.T().set<AdvectiveBC<scalar>>("xMax", &ctx);
    st.T().set<ZeroGradientBC<scalar>>("yMin");
    st.T().set<ZeroGradientBC<scalar>>("yMax");
    st.T().setFromFunction([](const Vec3& x) { return std::exp(-std::pow((x.x - 0.5) / 0.1, 2)); });
    st.initialize();
    const label pOut = mesh->findPatch("xMax");
    const scalar dt = 0.005;
    scalar maxDiff = 0, peakOut = 0;
    for (int n = 0; n < 200; ++n) {  // t = 1：脉冲中心已到 x = 1.5
        st.step(dt, &phi);
        // 出口面值应接近其内部单元值（纯对流、无反射）
        for (label i = 0; i < mesh->patches()[pOut].size; ++i) {
            const scalar b = st.T().patchValue(pOut, i), c = st.T()[st.T().bc(pOut).faceCell(i)];
            maxDiff = std::max(maxDiff, std::abs(b - c));
            peakOut = std::max(peakOut, b);
        }
    }
    maxDiff = par::allMax(maxDiff);
    peakOut = par::allMax(peakOut);
    scalar remain = 0;
    for (label c = 0; c < mesh->nCells(); ++c) remain = std::max(remain, std::abs(st.T()[c]));
    remain = par::allMax(remain);
    if (par::master()) std::cout << "advective: outlet peak " << peakOut << ", max |T_b - T_P| " << maxDiff << ", remaining " << remain << '\n';
    CHECK(peakOut > 0.9);        // 脉冲几乎无衰减地到达出口
    CHECK(maxDiff < 0.06);       // 面值紧跟内部（一阶空间差 ~ Δx·∂T/∂x）
    CHECK(remain < 1e-3);        // 流出后域内不残留反射
    // 稳态：w = 0 → 零梯度
    st.steady = true;
    ctx.ddt = DdtScheme::Steady;
    st.T().setFromFunction([](const Vec3& x) { return x.x; });
    st.T().updateBCs(0.0);
    st.T().correctBoundaryConditions();
    for (label i = 0; i < mesh->patches()[pOut].size; ++i)
        CHECK(std::abs(st.T().patchValue(pOut, i) - st.T()[st.T().bc(pOut).faceCell(i)]) < 1e-14);
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    try {
        testCoefficientsAndSwitching();
        testFunctionsAndParsing();
        testTurbulentInlet();
        testAdvective();
    } catch (const std::exception& e) {
        std::cerr << "exception: " << e.what() << std::endl;
        ++fails;
    }
    const int total = int(par::allMax(fails));
    if (par::master()) std::cout << (total ? "FAILED" : "all boundary condition tests passed") << std::endl;
    return total ? 1 : 0;
}
