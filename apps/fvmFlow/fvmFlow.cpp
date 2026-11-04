// fvmFlow：通用不可压流动求解程序，由 JSON 算例文件驱动（稳态 SIMPLE / 瞬态 PIMPLE）。
//
//   mpirun -np 4 fvmFlow cases/pitzDaily/case.json
//
// 算例文件结构见 docs/case-format.md 与 cases/ 下的示例。
#include "fvm/app/CaseSetup.h"
#include "fvm/field/Functions.h"
#include "fvm/io/Probes.h"
#include "fvm/io/Restart.h"
#include "fvm/io/TecplotWriter.h"
#include "fvm/io/VtkWriter.h"
#include "fvm/models/Energy.h"
#include "fvm/models/TurbulenceModel.h"
#include "fvm/post/ChannelStatistics.h"
#include "fvm/post/CutPlane.h"
#include "fvm/post/FieldAverage.h"
#include "fvm/post/Forces.h"
#include "fvm/post/SpanStatistics.h"
#include "fvm/post/WallOutput.h"
#include "fvm/solvers/IncompressibleFlow.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace cfd;

// 槽道 LES 初场：1/7 次方律平均剖面（体平均 = Ubulk）+ 无散的流向涡/条带扰动 + 小幅随机扰动。
// 随机数由单元中心坐标散列得到，与进程划分无关。
static void initChannelFlow(IncompressibleFlow& flow, const Json& d) {
    const Mesh& m = flow.mesh();
    const int sx = int(d.get("streamwise", 0.0)), ny = int(d.get("normal", 1.0)), sz = 3 - sx - ny;
    const auto& w = d["walls"].array();
    if (w.size() != 2) throw std::runtime_error("initial.channel.walls must be [y0, y1]");
    const scalar y0 = w[0].number(), y1 = w[1].number(), H = 0.5 * (y1 - y0);
    const scalar Ub = d.get("Ubulk", 1.0), A = d.get("amplitude", 0.1), noise = d.get("noise", 0.02);
    scalar lo[3] = {GREAT, GREAT, GREAT}, hi[3] = {-GREAT, -GREAT, -GREAT};
    for (label c = 0; c < m.nCells(); ++c)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], m.C()[c][k]);
            hi[k] = std::max(hi[k], m.C()[c][k]);
        }
    for (int k = 0; k < 3; ++k) {
        lo[k] = par::allMin(lo[k]);
        hi[k] = par::allMax(hi[k]);
    }
    const scalar Lx = std::max(hi[sx] - lo[sx], SMALL), Lz = std::max(hi[sz] - lo[sz], SMALL);
    const scalar kx = 2 * M_PI * 2 / Lx, kz = 2 * M_PI * 4 / Lz, pi = M_PI;
    auto& U = flow.U();
    for (label c = 0; c < m.nCells(); ++c) {
        const Vec3& x = m.C()[c];
        const scalar y = x[ny], eta = (y - y0) / (2 * H);
        const scalar yw = std::max(std::min(y - y0, y1 - y), 0.0);
        Vec3 u{};
        u[sx] = Ub * 8.0 / 7.0 * std::pow(std::min(yw / H, 1.0), 1.0 / 7.0);
        // 流函数 ψ = A Ub (2H/π) sin²(πη) sin(kx x) cos(kz z)：u_x = ∂ψ/∂y，u_y = −∂ψ/∂x（无散）
        const scalar s2 = std::sin(pi * eta);
        const scalar dpsidy = A * Ub * 2 * s2 * std::cos(pi * eta) * std::sin(kx * x[sx]) * std::cos(kz * x[sz]);
        const scalar dpsidx = A * Ub * (2 * H / pi) * s2 * s2 * kx * std::cos(kx * x[sx]) * std::cos(kz * x[sz]);
        u[sx] += dpsidy + A * Ub * s2 * std::cos(kz * x[sz]);  // 加流向条带
        u[ny] -= dpsidx;
        for (int k = 0; k < 3; ++k) u[k] += noise * Ub * s2 * coordHash(x, k);
        U[c] = u;
    }
    U.correctBoundaryConditions();
}


// Taylor-Green 涡初场（三维，周期域 [0, 2πL]³）：
//   u = U0 sin(x/L) cos(y/L) cos(z/L)，v = −U0 cos(x/L) sin(y/L) cos(z/L)，w = 0，
//   p = p0 + U0²/16 (cos(2x/L) + cos(2y/L)) (cos(2z/L) + 2)
static void initTaylorGreen(IncompressibleFlow& flow, const Json& d) {
    const Mesh& m = flow.mesh();
    const scalar U0 = d.get("U0", 1.0), L = d.get("L", 1.0), p0 = d.get("p0", 0.0);
    auto& U = flow.U();
    auto& p = flow.p();
    for (label c = 0; c < m.nCells(); ++c) {
        const Vec3& x = m.C()[c];
        U[c] = {U0 * std::sin(x.x / L) * std::cos(x.y / L) * std::cos(x.z / L),
                -U0 * std::cos(x.x / L) * std::sin(x.y / L) * std::cos(x.z / L), 0.0};
        p[c] = p0 + U0 * U0 / 16 * (std::cos(2 * x.x / L) + std::cos(2 * x.y / L)) * (std::cos(2 * x.z / L) + 2);
    }
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();
}

// 全域监控量（体积平均）：动能 k = ½<|u|²>，粘性耗散 ε = 2ν<S:S>，亚格子耗散 2<ν_t S:S>，拟涡能 ½<|ω|²>
struct Monitor {
    IncompressibleFlow& flow;
    std::ofstream file;
    int interval = 1;
    // append 为 true（续算）时接着原文件写，不重写表头
    Monitor(IncompressibleFlow& f, const std::string& path, int iv, bool append) : flow(f), interval(iv) {
        if (par::master()) {
            file.open(path, append ? std::ios::app : std::ios::out);
            if (!append) file << "t,k,epsilon,epsilon_sgs,enstrophy\n";
            file << std::setprecision(10);
        }
    }
    void sample(scalar t) {
        const Mesh& m = flow.mesh();
        flow.updateGradU();
        const auto& G = flow.gradU();
        const auto& U = flow.U();
        const auto& nut = flow.nut();
        par::SumAcc s[4];
        for (label c = 0; c < m.nCells(); ++c) {
            const scalar V = m.V()[c];
            const Tensor S = symm(G[c]);
            const Vec3 w{G[c](2, 1) - G[c](1, 2), G[c](0, 2) - G[c](2, 0), G[c](1, 0) - G[c](0, 1)};
            s[0].add(0.5 * magSqr(U[c]) * V);
            s[1].add(2 * flow.nu() * magSqr(S) * V);
            s[2].add(2 * (nut.empty() ? 0.0 : nut[c]) * magSqr(S) * V);
            s[3].add(0.5 * magSqr(w) * V);
        }
        double r[4];
        par::SumAcc::allReduce(s, 4, r);
        if (par::master()) {
            file << t;
            for (double v : r) file << ',' << v / m.totalVolume();
            file << '\n' << std::flush;
        }
    }
};

namespace {

void runCase(const std::string& caseFile) {
    const Json cfg = Json::parseFile(caseFile);
    const std::string dir = caseDirectory(caseFile);

    auto mesh = buildMesh(cfg["mesh"], dir);
    mesh->printSummary();

    const Json& phys = cfg["physics"];
    IncompressibleFlow flow(mesh, phys["nu"].number());
    const bool steady = phys.get("steady", false);
    auto& ctl = flow.controls;
    ctl.steady = steady;

    // 格式
    const Json& sch = cfg["schemes"];
    flow.divScheme = ConvectionScheme::parse(sch.get("div", steady ? "linearUpwind" : "linear"));
    // 梯度格式：未指定时，偏斜网格（三角形/四面体等非结构网格）默认最小二乘，否则 Gauss 线性；
    // 设为网格默认格式后扩散项非正交修正、linearUpwind、湍流模型等都按它取梯度
    const bool skewed = mesh->meanSkewness() > 0.05;
    flow.gradScheme = parseGradScheme(sch.has("grad") ? sch["grad"].string() : (skewed ? "leastSquares" : "Gauss"));
    mesh->setGradScheme(flow.gradScheme);
    // 压力梯度（动量方程源项与速度修正）：偏斜网格默认守恒的 GaussSkew，其余跟随 grad
    flow.gradPScheme = sch.has("gradP") ? parseGradScheme(sch["gradP"].string()) : (skewed ? GradScheme::GaussSkew : flow.gradScheme);
    if (!sch.has("gradP") && skewed) std::cout << "gradP scheme: GaussSkew (conservative, default for skewed mesh)\n";
    if (!sch.has("grad") && skewed)
        std::cout << "grad scheme: leastSquares (default for skewed mesh, mean skewness " << mesh->meanSkewness()
                  << ")\n";
    flow.ddtScheme = steady ? DdtScheme::Steady : parseDdtScheme(sch.get("ddt", "backward"));
    flow.laplacianNonOrthCorr = sch.get("laplacianNonOrthCorr", true);

    // 求解控制
    const Json& sol = cfg["solution"];
    ctl.nOuter = sol.get("nOuter", 1);
    ctl.nCorr = sol.get("nCorr", steady ? 1 : 2);
    ctl.nNonOrthCorr = sol.get("nNonOrthCorr", 0);
    ctl.momentumPredictor = sol.get("momentumPredictor", true);
    ctl.alphaU = sol.get("alphaU", steady ? 0.7 : 1.0);
    ctl.alphaP = sol.get("alphaP", steady ? 0.3 : 1.0);
    ctl.ddtPhiCoeff = sol.get("ddtPhiCoeff", 0.0);
    ctl.verbose = sol.get("verbose", false);
    if (steady) {
        ctl.UControls = {"PBiCGStab", "DILU", 1e-10, 0.1, 1000, 0, 1};
        ctl.pControls = {"PCG", "GAMG", 1e-10, 0.05, 2000, 0, 1};
    }
    ctl.UControls = parseSolverControls(sol["U"], ctl.UControls);
    ctl.pControls = parseSolverControls(sol["p"], ctl.pControls);
    ctl.pFinalControls = steady ? ctl.pControls : parseSolverControls(sol["pFinal"], ctl.pFinalControls);
    ctl.pRefCell = sol.get("pRefCell", 0);
    ctl.pRefValue = sol.get("pRefValue", 0.0);

    if (cfg.has("forcing")) {
        const Json& fo = cfg["forcing"];
        const std::string t = fo.get("type", "");
        if (t == "constant") {
            flow.forcing.mode = MomentumForcing::Mode::Constant;
            flow.forcing.force = fo["force"].vec3();
        } else if (t == "meanVelocity") {
            flow.forcing.mode = MomentumForcing::Mode::MeanVelocity;
            flow.forcing.Ubar = fo["Ubar"].vec3();
            flow.forcing.relaxation = fo.get("relaxation", 1.0);
        } else {
            throw std::runtime_error("forcing.type must be constant or meanVelocity");
        }
    }

    // 边界与初值
    const Json& bnd = cfg["boundary"];
    flow.context().caseDir = dir;
    setBoundaryConditions(flow.U(), bnd, Json::parse(R"({"type":"noSlip"})"), &flow.context());
    const bool buoyant = cfg.has("energy") && cfg["energy"].has("buoyancy");
    setBoundaryConditions(flow.p(), bnd,
                          Json::parse(buoyant ? R"({"type":"fixedFluxPressure"})" : R"({"type":"zeroGradient"})"),
                          &flow.context());
    const Json& ini = cfg["initial"];
    flow.U().setUniform(ini.get("U", Vec3{0, 0, 0}));
    flow.p().setUniform(ini.get("p", 0.0));
    if (ini.has("channel")) initChannelFlow(flow, ini["channel"]);
    if (ini.has("taylorGreen")) initTaylorGreen(flow, ini["taylorGreen"]);
    if (ini.has("noise")) {
        // 叠加随机扰动（绝对幅值），用于打破对称、触发转捩
        const scalar a = ini["noise"].number();
        auto& U = flow.U();
        for (label c = 0; c < mesh->nCells(); ++c)
            for (int k = 0; k < 3; ++k) U[c][k] += a * coordHash(mesh->C()[c], k);
        U.correctBoundaryConditions();
    }
    flow.initialize();
    // 混合初始化（与 Fluent hybrid initialization 同效）：势流投影得到无散的速度场，压力按 Bernoulli 给出，
    // 湍流量取入流边界值（湍流模型创建后做）；单独的 potentialFlow 只做投影
    const bool hybrid = ini.get("hybrid", false);
    if (ini.get("potentialFlow", false) || hybrid) {
        // 初场投影成无散场（曲面域、周期域的均匀初场散度很大，稳态第一步可能发散）
        flow.projectDivergenceFree(int(ini.get("potentialFlowNonOrthCorr", 2.0)));
        if (par::master()) std::cout << "initial field projected to be divergence-free\n";
    }
    if (hybrid) {
        const scalar H = flow.setBernoulliPressure();
        if (par::master()) std::cout << "hybrid initialisation: potential flow + Bernoulli pressure (total head " << H << ")\n";
    }

    // 湍流模型（缺省为层流）
    std::unique_ptr<TurbulenceModel> turb;
    if (cfg.has("turbulence")) turb = TurbulenceModel::create(flow, cfg["turbulence"], bnd, ini);
    if (turb && par::master()) std::cout << "turbulence model: " << turb->type() << '\n';
    if (turb && hybrid) turb->initializeFromInlets();
    // 能量方程（可选）：温度输运 + Boussinesq 浮力 + 热壁面函数
    std::unique_ptr<EnergyModel> energy;
    if (cfg.has("energy")) {
        energy = std::make_unique<EnergyModel>(flow, cfg["energy"], bnd, ini);
        if (par::master())
            std::cout << "energy: alpha " << energy->alpha() << " Prt " << energy->Prt()
                      << (energy->buoyant() ? " with Boussinesq buoyancy" : "") << '\n';
    }

    // 运行控制与输出
    const Json& run = cfg["run"];
    const std::string out = resolvePath(dir, run.get("output", "output"));
    std::string restartReadDir;
    if (run.has("restartRead")) {
        const std::string rd = resolvePath(dir, run["restartRead"].string());
        restartReadDir = rd;
        if (!restart::read(rd, flow)) throw std::runtime_error("restart data not found in " + rd);
        if (turb) turb->readRestart(rd);
        if (energy) energy->readRestart(rd);
        std::cout << "restarted from " << rd << " at t = " << flow.time().time << '\n';
    }
    const std::string restartDir = run.has("restartWrite") ? resolvePath(dir, run["restartWrite"].string()) : "";

    VtkWriter vtk(mesh, out);
    vtk.add(flow.U());
    vtk.add(flow.p());
    if (turb) {
        for (auto* f : turb->fields()) vtk.add(*f);
        vtk.add("nut", flow.nut());
    }
    if (energy)
        for (auto* f : energy->fields()) vtk.add(*f);
    const bool writeVtk = run.get("vtk", true);

    // 网格质量：坏单元统计（所含面的最大非正交角、偏斜度、面积比超过阈值），有坏单元时写 badCells.csv；
    // run.meshQuality 为 true 时三个单元质量量也写进 VTK
    std::vector<scalar> qNonOrth = mesh->cellNonOrthogonality(), qSkew = mesh->cellSkewness(),
                        qAspect = mesh->cellAspectRatio();
    {
        const scalar nonOrthWarn = run.get("nonOrthWarn", 70.0), skewWarn = run.get("skewWarn", 1.0),
                     aspectWarn = run.get("aspectWarn", 1000.0);
        std::vector<double> bad;
        double cnt[3] = {0, 0, 0};
        for (label c = 0; c < mesh->nCells(); ++c) {
            const bool b0 = qNonOrth[c] > nonOrthWarn, b1 = qSkew[c] > skewWarn, b2 = qAspect[c] > aspectWarn;
            cnt[0] += b0;
            cnt[1] += b1;
            cnt[2] += b2;
            if (b0 || b1 || b2) {
                const Vec3& x = mesh->C()[c];
                bad.insert(bad.end(), {double(mesh->cellGlobal()[c]), x.x, x.y, x.z, qNonOrth[c], qSkew[c], qAspect[c]});
            }
        }
        par::allSumInPlace(cnt, 3);
        const auto all = par::gatherDoublesToMaster(bad);
        if (par::master() && cnt[0] + cnt[1] + cnt[2] > 0) {
            std::cout << "mesh quality: " << cnt[0] << " cells with non-orthogonality > " << nonOrthWarn << " deg, "
                      << cnt[1] << " with skewness > " << skewWarn << ", " << cnt[2] << " with face area ratio > "
                      << aspectWarn << " -> " << out << "/badCells.csv\n";
            std::filesystem::create_directories(out);
            std::ofstream os(out + "/badCells.csv");
            os << "cell,x,y,z,nonOrthogonality,skewness,aspectRatio\n" << std::setprecision(10);
            for (std::size_t k = 0; k + 6 < all.size(); k += 7) {
                os << glabel(all[k]);
                for (int j = 1; j < 7; ++j) os << ',' << all[k + j];
                os << '\n';
            }
        }
        if (run.get("meshQuality", false)) {
            vtk.add("nonOrthogonality", qNonOrth);
            vtk.add("skewness", qSkew);
            vtk.add("aspectRatio", qAspect);
        }
    }

    // 场的时间平均（fieldAverage）：默认 U、p 及其二阶矩；start 之后累积（瞬态按 dt 加权，稳态按迭代计数）
    std::unique_ptr<FieldAverage> favg;
    scalar favgStart = 0;
    int favgInterval = 1;
    if (cfg.has("fieldAverage")) {
        const Json& fa = cfg["fieldAverage"];
        favg = std::make_unique<FieldAverage>(mesh);
        favgStart = fa.get("start", 0.0);
        favgInterval = std::max(fa.get("interval", 1), 1);
        const bool p2 = fa.get("prime2Mean", true);
        std::vector<std::string> names{"U", "p"};
        if (fa.has("fields")) {
            names.clear();
            for (auto& n : fa["fields"].array()) names.push_back(n.string());
        }
        for (const auto& n : names) {
            if (n == "U") {
                favg->add("U", flow.U().internal(), p2);
            } else if (n == "p") {
                favg->add("p", flow.p().internal(), p2);
            } else if (n == "nut" && turb) {
                favg->add("nut", flow.nut(), false);
            } else if (n == "T" && energy) {
                favg->add("T", energy->T().internal(), p2);
            } else {
                bool found = false;
                if (turb)
                    for (auto* f : turb->fields())
                        if (f->name() == n) {
                            favg->add(n, f->internal(), p2);
                            found = true;
                        }
                if (!found) throw std::runtime_error("fieldAverage: unknown field " + n);
            }
        }
        favg->addTo(vtk);
        if (!restartReadDir.empty() && favg->readRestart(restartReadDir) && par::master())
            std::cout << "fieldAverage: resumed, averaged over " << favg->averagedTime() << '\n';
    }

    // 壁面 patch 输出：p、Cp、τ_w、Cf、y+（、q、T_w）→ wall_<patch>.csv 与 VTK 表面 surface_<patch>.pvd
    std::unique_ptr<WallOutput> wallOut;
    bool wallCsv = true, wallVtk = true, wallEveryWrite = false;
    if (cfg.has("wallOutput")) {
        const Json& wd = cfg["wallOutput"];
        std::vector<std::string> names;
        if (wd.has("patches"))
            for (auto& n : wd["patches"].array()) names.push_back(n.string());
        wallOut = std::make_unique<WallOutput>(flow, turb.get(), energy.get(), names);
        wallOut->Uref = wd.get("Uref", 1.0);
        wallOut->pRef = wd.get("pRef", 0.0);
        wallOut->flowDir = wd.get("flowDir", Vec3{0, 0, 0});
        wallCsv = wd.get("csv", true);
        wallVtk = wd.get("vtk", true);
        wallEveryWrite = wd.get("everyWrite", false);
    }
    auto writeWallFrame = [&](scalar t) {
        if (!wallOut || !wallVtk || !wallEveryWrite) return;
        wallOut->update();
        wallOut->writeVtk(out, t);
    };

    // 平面采样（planes）：被平面切过的单元在投影点上的场值，结束时写 plane_<name>_<field>.csv
    struct PlaneItem {
        std::string name;
        std::unique_ptr<CutPlane> plane;
    };
    std::vector<PlaneItem> planes;
    if (cfg.has("planes"))
        for (auto& pl : cfg["planes"].array()) {
            PlaneItem it;
            it.name = pl.get("name", "plane");
            it.plane = std::make_unique<CutPlane>(mesh, pl["point"].vec3(), pl["normal"].vec3());
            const glabel nc = it.plane->nCells();
            if (par::master()) std::cout << "plane " << it.name << ": " << nc << " cells\n";
            planes.push_back(std::move(it));
        }

    std::unique_ptr<Monitor> monitor;
    if (cfg.has("monitor")) {
        const Json& mo = cfg["monitor"];
        const bool resuming = run.has("restartRead");
        monitor = std::make_unique<Monitor>(flow, out + "/" + mo.get("file", std::string("monitor.csv")),
                                            int(mo.get("interval", 1.0)), resuming);
        if (!resuming) monitor->sample(flow.time().time);
    }
    std::unique_ptr<Probes> probes;
    int probeInterval = 1;
    if (cfg.has("probes")) {
        const Json& pr = cfg["probes"];
        std::vector<Vec3> pts;
        for (auto& p : pr["points"].array()) pts.push_back(p.vec3());
        probes = std::make_unique<Probes>(mesh, pts, out + "/probes");
        probes->add(flow.U());
        probes->add(flow.p());
        probeInterval = pr.get("interval", 1);
    }
    // 壁面力与力系数：forces.csv 时间序列（t 从 averageStart 起做时间加权平均）、结束时写表面 Cp/Cf
    struct ForceMonitor {
        std::unique_ptr<Forces> f;
        Vec3 dragDir{1, 0, 0}, liftDir{0, 1, 0}, pitchAxis{0, 0, 1};
        scalar Uref = 1, lRef = 1, span = 1, averageStart = GREAT;
        int interval = 1;
        std::ofstream csv;
        scalar sumT = 0, sumCd = 0, sumCl = 0, sumCm = 0;
    };
    std::unique_ptr<ForceMonitor> fm;
    if (cfg.has("forces")) {
        const Json& fd = cfg["forces"];
        std::vector<std::string> names;
        for (auto& n : fd["patches"].array()) names.push_back(n.string());
        fm = std::make_unique<ForceMonitor>();
        fm->f = std::make_unique<Forces>(flow, names);
        fm->f->rho = fd.get("rho", 1.0);
        fm->f->pRef = fd.get("pRef", 0.0);
        fm->f->CofR = fd.get("CofR", Vec3{0, 0, 0});
        fm->dragDir = fd.get("dragDir", fm->dragDir);
        fm->liftDir = fd.get("liftDir", fm->liftDir);
        fm->pitchAxis = fd.get("pitchAxis", fm->pitchAxis);
        fm->Uref = fd.get("Uref", 1.0);
        fm->lRef = fd.get("lRef", 1.0);
        fm->span = fd.get("span", fm->f->span());
        fm->interval = fd.get("interval", 1);
        fm->averageStart = fd.get("averageStart", steady ? GREAT : 0.0);
        if (par::master()) {
            std::filesystem::create_directories(out);
            fm->csv.open(out + "/forces.csv", flow.time().time > 0 ? std::ios::app : std::ios::trunc);
            if (flow.time().time <= 0) fm->csv << "t,Cd,Cl,Cm,Cd_p,Cd_v,Cl_p,Cl_v,Fx,Fy,Fz\n";
            fm->csv << std::setprecision(10);
        }
    }
    // 槽道湍流统计（OpenLB channel3d 同格式）：t >= start 后每步按 dt 加权累积，写出时刻与结束时写 CSV
    std::unique_ptr<ChannelStatistics> chStats;
    scalar chStart = 0;
    std::string chFile;
    if (cfg.has("channelStatistics")) {
        const Json& cd = cfg["channelStatistics"];
        // walls：两壁面在法向上的坐标 [y0, y1]
        const auto& w = cd["walls"].array();
        if (w.size() != 2) throw std::runtime_error("channelStatistics.walls must be [y0, y1]");
        chStats = std::make_unique<ChannelStatistics>(flow, int(cd.get("streamwise", 0.0)), int(cd.get("normal", 1.0)),
                                                      w[0].number(), w[1].number());
        chStart = cd.get("start", 0.0);
        chFile = out + "/" + cd.get("file", std::string("channel_stats.csv"));
        if (!restartReadDir.empty() && chStats->readRestart(restartReadDir) && par::master())
            std::cout << "channelStatistics: resumed, averaged over " << chStats->averagedTime() << '\n';
    }
    // 展向均匀流动统计（周期山等，OpenLB periodichill3d 同格式）
    std::unique_ptr<SpanStatistics> spStats;
    scalar spStart = 0;
    std::string spDir;
    std::vector<scalar> spStations;
    if (cfg.has("spanStatistics")) {
        const Json& sd = cfg["spanStatistics"];
        spStats = std::make_unique<SpanStatistics>(flow, int(sd.get("streamwise", 0.0)), int(sd.get("normal", 1.0)),
                                                   sd.get("wall", std::string("bottom")), sd.get("Uref", 1.0),
                                                   sd.get("H", 1.0));
        if (sd.has("pRef")) spStats->setPRef(sd["pRef"].number());
        spStart = sd.get("start", 0.0);
        spDir = out + "/" + sd.get("dir", std::string("spanStatistics"));
        if (sd.has("stations"))
            for (auto& v : sd["stations"].array()) spStations.push_back(v.number());
        if (!restartReadDir.empty() && spStats->readRestart(restartReadDir) && par::master())
            std::cout << "spanStatistics: resumed, averaged over " << spStats->averagedTime() << '\n';
    }
    auto writeSpanStats = [&](bool report) {
        if (!spStats) return;
        const auto sr = spStats->write(spDir, spStations);
        if (report && par::master())
            std::cout << "spanStatistics: averaged over " << spStats->averagedTime() << ", separation x/H " << sr.first
                      << ", reattachment x/H " << sr.second << " -> " << spDir << '\n';
    };
    // 返回 {Cd, Cl, Cm}
    auto sampleForces = [&](scalar t, scalar dt, bool writeRow) -> Vec3 {
        const ForceResult r = fm->f->compute();
        const scalar q = 0.5 * fm->f->rho * fm->Uref * fm->Uref * fm->lRef * fm->span;
        const Vec3 F = r.force();
        const scalar Cd = dot(F, fm->dragDir) / q, Cl = dot(F, fm->liftDir) / q;
        const scalar Cm = dot(r.moment(), fm->pitchAxis) / (q * fm->lRef);
        if (t >= fm->averageStart * (1 - 1e-12) && dt > 0) {
            fm->sumT += dt;
            fm->sumCd += Cd * dt;
            fm->sumCl += Cl * dt;
            fm->sumCm += Cm * dt;
        }
        if (writeRow && par::master())
            fm->csv << t << ',' << Cd << ',' << Cl << ',' << Cm << ',' << dot(r.pressure, fm->dragDir) / q << ','
                    << dot(r.viscous, fm->dragDir) / q << ',' << dot(r.pressure, fm->liftDir) / q << ','
                    << dot(r.viscous, fm->liftDir) / q << ',' << F.x << ',' << F.y << ',' << F.z << '\n';
        return {Cd, Cl, Cm};
    };

    // 各边界的体积通量（正为流出）；同一边界既有流入又有流出时括号内给出较小的那部分（如出口回流量）
    auto printFluxes = [&] {
        const Mesh& m = *mesh;
        const std::size_t np = m.patches().size();
        std::vector<double> q(2 * np, 0.0);  // [净通量..., 流入部分...]
        for (std::size_t p = 0; p < np; ++p) {
            const auto& pt = m.patches()[p];
            for (label f = pt.start; f < pt.end(); ++f) {
                q[p] += flow.phi()[f];
                if (flow.phi()[f] < 0) q[np + p] += flow.phi()[f];
            }
        }
        par::allSumInPlace(q.data(), int(q.size()));
        std::cout << "  boundary fluxes:";
        for (std::size_t p = 0; p < np; ++p) {
            if (m.patches()[p].type == PatchType::Cyclic) continue;
            std::cout << ' ' << m.patches()[p].name << ' ' << q[p];
            const double in = q[np + p], out = q[p] - in;
            if (in < 0 && out > 0) std::cout << (q[p] >= 0 ? " (backflow " : " (outflow ") << (q[p] >= 0 ? in : out) << ')';
        }
        std::cout << '\n';
        // 压力力一致性诊断：单元梯度的体积和（最小二乘 / Gauss）对边界面压力合力 Σ p_f S_f（精确的净压力）
        {
            const auto gLS = fvc::grad(flow.p(), GradScheme::LeastSquares);
            const auto gG = fvc::grad(flow.p(), GradScheme::GaussLinear);
            double s[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (label c = 0; c < m.nCells(); ++c)
                for (int k = 0; k < 3; ++k) {
                    s[k] += m.V()[c] * gLS[c][k];
                    s[3 + k] += m.V()[c] * gG[c][k];
                }
            for (label f = m.nInternalFaces(); f < m.nFaces(); ++f)
                for (int k = 0; k < 3; ++k) s[6 + k] += flow.p().bValue(f) * m.Sf()[f][k];
            par::allSumInPlace(s, 9);
            std::cout << "  pressure force: sum V grad(p) LS (" << s[0] << ' ' << s[1] << ' ' << s[2] << ")  Gauss (" << s[3] << ' '
                      << s[4] << ' ' << s[5] << ")  boundary sum p_f S_f (" << s[6] << ' ' << s[7] << ' ' << s[8] << ")\n";
        }
    };
    const bool fluxReport = run.get("printFluxes", false);

    auto writeLines = [&] {
        if (!cfg.has("lines")) return;
        for (auto& l : cfg["lines"].array()) {
            const Vec3 a = l["start"].vec3(), b = l["end"].vec3();
            const int n = l.get("n", 100);
            std::vector<Vec3> pts;
            for (int k = 0; k < n; ++k) pts.push_back(a + (b - a) * (n > 1 ? scalar(k) / (n - 1) : 0.0));
            Probes line(mesh, pts, out);
            if (l.get("interpolate", false)) line.setInterpolate(true);
            const std::string name = l.get("name", "line");
            line.writeTable(out + "/line_" + name + "_U.csv", flow.U());
            line.writeTable(out + "/line_" + name + "_p.csv", flow.p());
            if (turb) {
                for (auto* f : turb->fields()) line.writeTable(out + "/line_" + name + "_" + f->name() + ".csv", *f);
                line.writeTable(out + "/line_" + name + "_nut.csv", flow.nut().data(), 1, "nut");
            }
            if (energy)
                for (auto* f : energy->fields()) line.writeTable(out + "/line_" + name + "_" + f->name() + ".csv", *f);
            if (favg)
                for (const auto& v : favg->views())
                    line.writeTable(out + "/line_" + name + "_" + v.name + ".csv", v.data, v.nc, v.name);
        }
    };
    auto writePlanes = [&] {
        for (const auto& pl : planes) {
            const std::string pre = out + "/plane_" + pl.name + "_";
            pl.plane->writeTable(pre + "U.csv", flow.U());
            pl.plane->writeTable(pre + "p.csv", flow.p());
            if (turb) {
                for (auto* f : turb->fields()) pl.plane->writeTable(pre + f->name() + ".csv", *f);
                pl.plane->writeTable(pre + "nut.csv", flow.nut().data(), 1, "nut");
            }
            if (energy)
                for (auto* f : energy->fields()) pl.plane->writeTable(pre + f->name() + ".csv", *f);
            if (favg)
                for (const auto& v : favg->views()) pl.plane->writeTable(pre + v.name + ".csv", v.data, v.nc, v.name);
        }
    };
    // 发散诊断：统计非有限值的单元并给出前几个单元中心，写出当前场的 VTK 帧便于定位，然后报错退出
    auto reportDivergence = [&](const std::string& where) {
        const auto& U = flow.U();
        const auto& p = flow.p();
        std::vector<double> rec;
        label nBad = 0;
        for (label c = 0; c < mesh->nCells(); ++c) {
            if (std::isfinite(magSqr(U[c])) && std::isfinite(p[c])) continue;
            ++nBad;
            if (rec.size() < 5 * 4) {
                const Vec3& x = mesh->C()[c];
                rec.insert(rec.end(), {double(mesh->cellGlobal()[c]), x.x, x.y, x.z});
            }
        }
        const glabel nb = par::allSum(glabel(nBad));
        const auto all = par::gatherDoublesToMaster(rec);
        if (par::master()) {
            std::cerr << "diverged " << where << ": non-finite U/p in " << nb << " cells";
            for (std::size_t k = 0; k + 3 < all.size() && k < 5 * 4; k += 4)
                std::cerr << (k ? ", " : "; first at ") << "cell " << glabel(all[k]) << " (" << all[k + 1] << ' '
                          << all[k + 2] << ' ' << all[k + 3] << ')';
            std::cerr << '\n';
        }
        if (writeVtk) vtk.write(flow.time().time);
        throw std::runtime_error("solution diverged (non-finite residual) " + where);
    };

    const double t0 = par::wallTime();
    std::cout << std::setprecision(6);
    if (steady) {
        const int maxIter = run.get("maxIter", 5000);
        const scalar tol = run.get("residualTol", 1e-6);
        const int writeInterval = run.get("writeInterval", 0);
        const int printInterval = run.get("printInterval", 50);
        int it = 0;
        for (it = 1; it <= maxIter; ++it) {
            const StepInfo info = flow.step(1.0);
            const scalar TRes = energy ? energy->correct() : 0.0;
            if (!std::isfinite(info.UInitialResidual) || !std::isfinite(info.pInitialResidual) ||
                !std::isfinite(info.continuityError) || !std::isfinite(TRes))
                reportDivergence("at iteration " + std::to_string(it));
            if (favg && it >= favgStart && it % favgInterval == 0) favg->sample(1.0);
            const bool conv = it > 1 && info.UInitialResidual < tol && info.pInitialResidual < tol && TRes < tol;
            if (it == 1 || it % printInterval == 0 || conv) {
                std::cout << "iter " << it << "  U res " << info.UInitialResidual << "  p res " << info.pInitialResidual
                          << "  cont " << info.continuityError;
                if (energy) std::cout << "  T res " << TRes;
                std::cout << '\n';
            }
            if (turb && (it == 1 || it % printInterval == 0)) turb->printYPlus(std::cout);
            if (probes && it % probeInterval == 0) probes->sample(it);
            if (fm && (it % fm->interval == 0 || conv)) sampleForces(it, 0.0, true);
            if (writeVtk && writeInterval > 0 && it % writeInterval == 0) vtk.write(it);
            if (writeInterval > 0 && it % writeInterval == 0) writeWallFrame(it);
            if (conv) break;
        }
        std::cout << "finished after " << std::min(it, maxIter) << " iterations, " << par::wallTime() - t0 << " s\n";
        if (writeVtk) vtk.write(std::min(it, maxIter));
        if (spStats) spStats->sample(1.0);  // 稳态：收敛解本身就是统计量
    } else {
        const scalar endTime = run["endTime"].number();
        const scalar maxCo = run.get("maxCo", 0.0);  // >0：按库朗数自适应步长
        scalar dt = run.get("dt", 0.0);
        const scalar maxDt = run.get("maxDt", GREAT);
        const scalar writeInterval = run.get("writeInterval", 0.0);
        const int printInterval = run.get("printInterval", 10);
        if (dt <= 0 && maxCo <= 0) throw std::runtime_error("run: need dt or maxCo");
        if (dt <= 0) {
            // 初始步长：由初场库朗数估计
            const scalar co1 = flow.courantNumber(1.0);
            dt = co1 > 0 ? std::min(maxDt, 0.2 * maxCo / co1) : std::min(maxDt, 1e-3 * endTime);
        }
        scalar nextWrite = writeInterval > 0 ? flow.time().time + writeInterval : GREAT;
        if (writeVtk && writeInterval > 0) vtk.write(flow.time().time);
        int n = 0;
        while (flow.time().time < endTime * (1 - 1e-12)) {
            if (maxCo > 0) {
                const scalar co = flow.courantNumber(dt);
                // 与 OpenFOAM 相同：增幅限制 1.2 倍、降幅不限
                const scalar f = co > 0 ? maxCo / co : 1.2;
                dt = std::min({dt * std::min(f, 1.2), maxDt});
            }
            // 落在写出时刻 / 结束时刻上：提前几步把剩余时间均分，避免最后一步极小（BDF2 变步长时会产生压力尖峰）
            {
                const scalar target = std::min(endTime, nextWrite);
                const scalar rem = target - flow.time().time;
                if (rem <= dt * (1 + 1e-9)) {
                    dt = std::max(rem, 1e-12 * dt);
                } else {
                    const scalar nSteps = std::ceil(rem / dt - 1e-9);
                    if (nSteps <= 4) dt = rem / nSteps;
                }
            }
            const StepInfo info = flow.step(dt);
            if (energy) energy->correct();
            ++n;
            if (!std::isfinite(info.maxCo) || !std::isfinite(info.continuityError))
                reportDivergence("at t = " + std::to_string(flow.time().time));
            if (favg && flow.time().time >= favgStart * (1 - 1e-12) && n % favgInterval == 0)
                favg->sample(dt * favgInterval);
            if (n % printInterval == 0)
                std::cout << "t " << flow.time().time << "  dt " << dt << "  Co " << info.maxCo << "  cont "
                          << info.continuityError << '\n';
            if (probes && n % probeInterval == 0) probes->sample(flow.time().time);
            if (monitor && n % monitor->interval == 0) monitor->sample(flow.time().time);
            if (fm) {
                // 平均每步都累积（与采样间隔无关），写行按 interval
                const Vec3 c = sampleForces(flow.time().time, dt, n % fm->interval == 0);
                if (n % printInterval == 0) std::cout << "  Cd " << c.x << "  Cl " << c.y << '\n';
            }
            if (fluxReport && n % printInterval == 0) printFluxes();
            if (chStats && flow.time().time >= chStart * (1 - 1e-12)) chStats->sample(dt);
            if (spStats && flow.time().time >= spStart * (1 - 1e-12)) spStats->sample(dt);
            if (flow.time().time >= nextWrite * (1 - 1e-12)) {
                if (writeVtk) vtk.write(flow.time().time);
                writeWallFrame(flow.time().time);
                if (chStats) chStats->write(chFile);
                writeSpanStats(true);
                nextWrite += writeInterval;
            }
        }
        std::cout << "reached t = " << flow.time().time << " in " << n << " steps, " << par::wallTime() - t0 << " s\n";
        if (writeVtk && writeInterval <= 0) vtk.write(flow.time().time);
    }
    printFluxes();
    writeLines();
    writePlanes();
    writeSpanStats(true);
    if (chStats) {
        chStats->write(chFile);
        if (par::master())
            std::cout << "channelStatistics: averaged over " << chStats->averagedTime() << ", u_tau " << chStats->uTau()
                      << " -> " << chFile << '\n';
    }
    if (fm) {
        const Vec3 c = sampleForces(flow.time().time, 0.0, false);
        std::cout << "forces: final Cd " << c.x << "  Cl " << c.y << "  Cm " << c.z << '\n';
        if (fm->sumT > 0)
            std::cout << "forces: mean over t >= " << fm->averageStart << " (" << fm->sumT << "): Cd "
                      << fm->sumCd / fm->sumT << "  Cl " << fm->sumCl / fm->sumT << "  Cm " << fm->sumCm / fm->sumT
                      << '\n';
        fm->f->writeSurface(out + "/surface.csv", fm->Uref, fm->dragDir);
    }
    if (energy) {
        energy->printWallHeatFlux(std::cout);
        energy->writeWallHeatFlux(out);
    }
    if (wallOut) {
        wallOut->update();
        if (wallCsv) wallOut->writeCsv(out);
        if (wallVtk) wallOut->writeVtk(out, flow.time().time);
    }
    if (favg && par::master())
        std::cout << "fieldAverage: averaged over " << favg->averagedTime() << " (" << favg->nSamples() << " samples)\n";
    if (run.get("tecplot", false)) {
        TecplotWriter tp(mesh);
        tp.add(flow.U());
        tp.add(flow.p());
        if (turb) {
            for (auto* f : turb->fields()) tp.add(*f);
            tp.add("nut", flow.nut(), 1);
        }
        if (energy) tp.add(energy->T());
        if (favg) favg->addTo(tp);
        tp.write(out + "/fields.dat", "U_p");
    }
    if (!restartDir.empty()) {
        restart::write(restartDir, flow);
        if (turb) turb->writeRestart(restartDir);
        if (energy) energy->writeRestart(restartDir);
        if (favg) favg->writeRestart(restartDir);
        if (chStats) chStats->writeRestart(restartDir);
        if (spStats) spStats->writeRestart(restartDir);
    }
    if (turb) turb->printYPlus(std::cout);
    if (flow.forcing.mode == MomentumForcing::Mode::MeanVelocity) std::cout << "forcing gradP " << flow.forcing.gradP << '\n';
    {
        // 体平均速度（单元值）：与边界通量对照可检查单元速度与面通量的一致性（如管流 Q/A）
        par::SumAcc sx, sy, sz;
        for (label c = 0; c < mesh->nCells(); ++c) {
            const Vec3 uv = flow.U()[c] * mesh->V()[c];
            sx.add(uv.x);
            sy.add(uv.y);
            sz.add(uv.z);
        }
        const scalar vol = mesh->totalVolume();
        std::cout << "mean U (volume average) " << sx.allReduce() / vol << ' ' << sy.allReduce() / vol << ' '
                  << sz.allReduce() / vol << '\n';
    }
    std::cout << "checksum U " << checksum(flow.U()) << " p " << checksum(flow.p()) << '\n';
    if (par::master()) cfg.reportUnused(std::cout);
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    if (argc < 2) {
        std::cout << "usage: fvmFlow case.json\n";
        return 1;
    }
    try {
        runCase(argv[1]);
    } catch (const std::exception& e) {
        std::cerr << "[rank " << par::rank() << "] error: " << e.what() << '\n';
        par::abort(1);
    }
    return 0;
}
