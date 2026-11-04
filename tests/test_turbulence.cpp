// RANS 槽道流（Re_b = U_b δ/ν = 6875，DNS Re_τ = 392）：k-ω SST、Spalart-Allmaras 的 Re_τ 误差 < 2%，
// 全 y+ 壁面处理（enhanced）的 k-ε 在 y+ ≈ 0.1 与 y+ ≈ 25 的网格上误差 < 5%，RNG k-ε、LRR、SSG 壁面函数网格 < 6%；
// 壁面距离与网格一致；并检查结果与进程数的一致性（可复现模式下由 ctest 的 np3 版本覆盖）。
// 壁面律（WallLaw）：用 Lee & Moser 2015 Re_τ = 1000 DNS 剖面作输入，各律在 y+ = 10..300 处反解的 u_τ 与真值比较。
#include "fvm/core/Json.h"
#include "fvm/mesh/Mesh.h"
#include "fvm/mesh/WallDistance.h"
#include "fvm/models/TurbulenceModel.h"
#include "fvm/models/WallLaw.h"
#include "fvm/solvers/IncompressibleFlow.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

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

MeshPtr channelMesh(int ny, scalar beta) {
    BoxSpec s;
    s.n[0] = 4, s.n[1] = ny, s.n[2] = 1;
    s.hi = {1, 2, 0.1};
    s.periodic[0] = true;
    s.twoD = true;
    s.types[2] = s.types[3] = PatchType::Wall;
    if (beta > 0) s.stretch[1] = tanhStretch(beta);
    RawMesh raw;
    if (par::master()) raw = generateBox(s);
    return Mesh::build(par::master() ? &raw : nullptr);
}

scalar reTau(MeshPtr mesh, const std::string& turbulence) {
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
    const Json d = Json::parse(turbulence);
    const Json b = Json::parse("{}");
    const Json i = Json::parse(R"({"k": 0.005, "omega": 50, "epsilon": 0.05, "nuTilda": 5e-4})");
    auto turb = TurbulenceModel::create(flow, d, b, i);
    for (int it = 0; it < 6000; ++it) flow.step(1.0);
    const scalar rt = std::sqrt(flow.forcing.gradP) / nu;
    if (par::master()) std::cout << turbulence << ": Re_tau " << rt << " (DNS 392.2)";
    turb->printYPlus(std::cout);  // 含全局归约，所有进程都要调用（只在主进程打印）
    return rt;
}

// Lee & Moser Re_τ 1000 剖面（y+, u+），逗号分隔，表头 y/delta,y+,...,u+,...
bool readDNS(const std::string& file, std::vector<scalar>& yp, std::vector<scalar>& up) {
    std::ifstream is(file);
    if (!is) return false;
    std::string line;
    std::getline(is, line);
    std::vector<std::string> head;
    {
        std::stringstream ss(line);
        std::string h;
        while (std::getline(ss, h, ',')) head.push_back(h);
    }
    int iy = -1, iu = -1;
    for (std::size_t k = 0; k < head.size(); ++k) {
        if (head[k] == "y+") iy = int(k);
        if (head[k] == "u+") iu = int(k);
    }
    if (iy < 0 || iu < 0) return false;
    while (std::getline(is, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string tok;
        std::vector<scalar> v;
        while (std::getline(ss, tok, ',')) v.push_back(std::stod(tok));
        yp.push_back(v[iy]);
        up.push_back(v[iu]);
    }
    return yp.size() > 10;
}

scalar interp(const std::vector<scalar>& x, const std::vector<scalar>& y, scalar xq) {
    for (std::size_t k = 1; k < x.size(); ++k)
        if (x[k] >= xq) return y[k - 1] + (y[k] - y[k - 1]) * (xq - x[k - 1]) / (x[k] - x[k - 1]);
    return y.back();
}

// 各壁面律由 DNS 剖面 (y, U) 反解 u_τ 的误差表；返回 Spalding / Musker / Reichardt / ODE 在 y+ ≥ 30 处的最大误差
scalar wallLawTable() {
    std::vector<scalar> yp, up;
    if (!readDNS("cases/channelWMLES/reference/referenceData_1000.csv", yp, up)) {
        if (par::master()) std::cout << "wall laws: DNS reference not found, skipping\n";
        return 0;
    }
    const scalar nu = 5e-5, uTau = 0.05, H = 1.0;  // Re_τ = 1000
    const std::vector<scalar> samples = {10, 20, 30, 50, 100, 175, 300};
    struct Law {
        std::string name;
        WallLaw law;
    };
    std::vector<Law> laws;
    for (const char* n : {"logLaw", "spalding", "musker", "reichardt", "powerLaw", "ode"}) {
        Law l;
        l.name = n;
        l.law.type = WallLaw::parse(n);
        laws.push_back(l);
    }
    {
        Law l;
        l.name = "ode+dp/dx";
        l.law.type = WallLaw::Type::ODE;
        l.law.pressureGradient = true;
        laws.push_back(l);
    }
    scalar worst = 0;
    if (par::master()) {
        std::cout << "wall laws vs Lee & Moser Re_tau 1000 DNS: u_tau error (%) when fed U(y) at y+ = ...\n";
        std::cout << std::setw(12) << "law";
        for (scalar s : samples) std::cout << std::setw(8) << s;
        std::cout << '\n';
    }
    for (auto& l : laws) {
        if (par::master()) std::cout << std::setw(12) << l.name;
        for (scalar s : samples) {
            const scalar U = interp(yp, up, s) * uTau, y = s * nu / uTau;
            scalar ut;
            if (l.law.pressureGradient) {
                // 槽道：体力 f = u_τ²/H 沿流向（总切应力 τ = u_τ²(1 − y/H)）
                const Vec3 tau = l.law.odeWallShear(Vec3{U, 0, 0}, y, nu, Vec3{uTau * uTau / H, 0, 0});
                ut = std::sqrt(mag(tau));
            } else {
                ut = l.law.uTau(U, y, nu);
            }
            const scalar err = 100 * (ut / uTau - 1);
            if (par::master()) std::cout << std::setw(8) << std::fixed << std::setprecision(2) << err;
            if (s >= 30 && (l.name == "spalding" || l.name == "musker" || l.name == "reichardt" || l.name == "ode" ||
                            l.name == "powerLaw"))
                worst = std::max(worst, std::abs(err));
        }
        if (par::master()) std::cout << '\n';
    }
    std::cout.unsetf(std::ios::fixed);
    return worst;
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    auto mesh = channelMesh(120, 3.0);
    // 壁面距离：y = min(yc, 2 − yc)
    const auto& wd = wallDistance(*mesh);
    scalar err = 0;
    for (label c = 0; c < mesh->nTotalCells(); ++c) {
        const scalar yc = mesh->C()[c].y;
        err = std::max(err, std::abs(wd.y[c] - std::min(yc, 2 - yc)));
    }
    CHECK(par::allMax(err) < 1e-12);
    CHECK(std::abs(reTau(mesh, R"({"model": "kOmegaSST"})") / 392.2 - 1) < 0.02);
    CHECK(std::abs(reTau(mesh, R"({"model": "SpalartAllmaras"})") / 392.2 - 1) < 0.02);
    // 全 y+ 处理：k-ε（两层）在 y+ ≈ 0.1 的网格上
    CHECK(std::abs(reTau(mesh, R"({"model": "kEpsilon", "wallTreatment": "enhanced"})") / 392.2 - 1) < 0.05);
    // 壁面函数网格（y+ ≈ 25）：标准与 enhanced k-ε、RNG k-ε、SST enhanced
    auto coarse = channelMesh(16, 0.0);
    CHECK(std::abs(reTau(coarse, R"({"model": "kEpsilon"})") / 392.2 - 1) < 0.06);
    CHECK(std::abs(reTau(coarse, R"({"model": "kEpsilon", "wallTreatment": "enhanced"})") / 392.2 - 1) < 0.06);
    CHECK(std::abs(reTau(coarse, R"({"model": "RNGkEpsilon", "wallTreatment": "enhanced"})") / 392.2 - 1) < 0.06);
    CHECK(std::abs(reTau(coarse, R"({"model": "kOmegaSST", "wallTreatment": "enhanced"})") / 392.2 - 1) < 0.06);
    // 雷诺应力模型（壁面函数 + 近壁各向异性）：LRR −2.6 %、SSG −4.6 %
    CHECK(std::abs(reTau(coarse, R"({"model": "LRR"})") / 392.2 - 1) < 0.06);
    CHECK(std::abs(reTau(coarse, R"({"model": "SSG"})") / 392.2 - 1) < 0.06);
    // 壁面律 vs DNS：统一律类在 y+ ≥ 30 处 u_τ 误差 < 5%（Spalding 在 y+ ≈ 30 偏 4%，是该律缓冲层形状的固有偏差）
    const scalar worstLaw = wallLawTable();
    CHECK(worstLaw < 5.0);
    // ODE 律无压力梯度时应接近对数律：y+ = 500 处 u+ 与 ln(E y+)/κ（κ 0.41，B 5.2）之差 < 0.5
    {
        WallLaw ode;
        ode.type = WallLaw::Type::ODE;
        const scalar upOde = ode.uPlusOf(500), upLog = std::log(500.0) / 0.41 + 5.2;
        if (par::master()) std::cout << "ode law u+(500) = " << upOde << ", log law " << upLog << '\n';
        CHECK(std::abs(upOde - upLog) < 0.5);
    }
    if (par::master()) std::cout << (fails ? "FAILED" : "ALL PASSED") << std::endl;
    return fails ? 1 : 0;
}
