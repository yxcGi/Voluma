// fvmFlow：通用不可压流动求解程序，由 JSON 算例文件驱动（稳态 SIMPLE / 瞬态 PIMPLE）。
//
//   mpirun -np 4 fvmFlow cases/pitzDaily/case.json
//
// 算例文件结构见 docs/case-format.md 与 cases/ 下的示例。
#include "fvm/app/CaseSetup.h"
#include "fvm/io/Probes.h"
#include "fvm/io/Restart.h"
#include "fvm/io/TecplotWriter.h"
#include "fvm/io/VtkWriter.h"
#include "fvm/models/TurbulenceModel.h"
#include "fvm/post/ChannelStatistics.h"
#include "fvm/post/Forces.h"
#include "fvm/post/SpanStatistics.h"
#include "fvm/solvers/IncompressibleFlow.h"

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
    auto hash = [](scalar a, scalar b, scalar c, int comp) {
        std::uint64_t h = 1469598103934665603ull ^ std::uint64_t(comp);
        for (scalar v : {a, b, c}) {
            std::int64_t q = std::llround(v * 1e9);
            h = (h ^ std::uint64_t(q)) * 1099511628211ull;
            h ^= h >> 29;
        }
        return double(h >> 11) / double(1ull << 53) * 2 - 1;  // [-1, 1)
    };
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
        for (int k = 0; k < 3; ++k) u[k] += noise * Ub * s2 * hash(x.x, x.y, x.z, k);
        U[c] = u;
    }
    U.correctBoundaryConditions();
}


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
    flow.gradScheme = parseGradScheme(sch.get("grad", "Gauss"));
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
    setBoundaryConditions(flow.U(), bnd, Json::parse(R"({"type":"noSlip"})"));
    setBoundaryConditions(flow.p(), bnd, Json::parse(R"({"type":"zeroGradient"})"));
    const Json& ini = cfg["initial"];
    flow.U().setUniform(ini.get("U", Vec3{0, 0, 0}));
    flow.p().setUniform(ini.get("p", 0.0));
    if (ini.has("channel")) initChannelFlow(flow, ini["channel"]);
    flow.initialize();

    // 湍流模型（缺省为层流）
    std::unique_ptr<TurbulenceModel> turb;
    if (cfg.has("turbulence")) turb = TurbulenceModel::create(flow, cfg["turbulence"], bnd, ini);
    if (turb && par::master()) std::cout << "turbulence model: " << turb->type() << '\n';

    // 运行控制与输出
    const Json& run = cfg["run"];
    const std::string out = resolvePath(dir, run.get("output", "output"));
    if (run.has("restartRead")) {
        const std::string rd = resolvePath(dir, run["restartRead"].string());
        if (!restart::read(rd, flow)) throw std::runtime_error("restart data not found in " + rd);
        if (turb) turb->readRestart(rd);
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
    const bool writeVtk = run.get("vtk", true);

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
        spStart = sd.get("start", 0.0);
        spDir = out + "/" + sd.get("dir", std::string("spanStatistics"));
        if (sd.has("stations"))
            for (auto& v : sd["stations"].array()) spStations.push_back(v.number());
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

    // 各边界的体积通量（正为流出）
    auto printFluxes = [&] {
        const Mesh& m = *mesh;
        std::vector<double> q(m.patches().size(), 0.0);
        for (std::size_t p = 0; p < q.size(); ++p) {
            const auto& pt = m.patches()[p];
            for (label f = pt.start; f < pt.end(); ++f) q[p] += flow.phi()[f];
        }
        par::allSumInPlace(q.data(), int(q.size()));
        std::cout << "  boundary fluxes:";
        for (std::size_t p = 0; p < q.size(); ++p)
            if (m.patches()[p].type != PatchType::Cyclic) std::cout << ' ' << m.patches()[p].name << ' ' << q[p];
        std::cout << '\n';
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
            const std::string name = l.get("name", "line");
            line.writeTable(out + "/line_" + name + "_U.csv", flow.U());
            line.writeTable(out + "/line_" + name + "_p.csv", flow.p());
        }
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
            const bool conv = it > 1 && info.UInitialResidual < tol && info.pInitialResidual < tol;
            if (it == 1 || it % printInterval == 0 || conv)
                std::cout << "iter " << it << "  U res " << info.UInitialResidual << "  p res " << info.pInitialResidual
                          << "  cont " << info.continuityError << '\n';
            if (turb && (it == 1 || it % printInterval == 0)) turb->printYPlus(std::cout);
            if (probes && it % probeInterval == 0) probes->sample(it);
            if (fm && (it % fm->interval == 0 || conv)) sampleForces(it, 0.0, true);
            if (writeVtk && writeInterval > 0 && it % writeInterval == 0) vtk.write(it);
            if (conv) break;
        }
        std::cout << "finished after " << std::min(it, maxIter) << " iterations, " << par::wallTime() - t0 << " s\n";
        if (writeVtk) vtk.write(std::min(it, maxIter));
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
            ++n;
            if (n % printInterval == 0)
                std::cout << "t " << flow.time().time << "  dt " << dt << "  Co " << info.maxCo << "  cont "
                          << info.continuityError << '\n';
            if (probes && n % probeInterval == 0) probes->sample(flow.time().time);
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
    if (run.get("tecplot", false)) {
        TecplotWriter tp(mesh);
        tp.add(flow.U());
        tp.add(flow.p());
        tp.write(out + "/fields.dat", "U_p");
    }
    if (!restartDir.empty()) {
        restart::write(restartDir, flow);
        if (turb) turb->writeRestart(restartDir);
    }
    if (turb) turb->printYPlus(std::cout);
    if (flow.forcing.mode == MomentumForcing::Mode::MeanVelocity) std::cout << "forcing gradP " << flow.forcing.gradP << '\n';
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
