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
#include "fvm/solvers/IncompressibleFlow.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace cfd;

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
    flow.initialize();

    // 运行控制与输出
    const Json& run = cfg["run"];
    const std::string out = resolvePath(dir, run.get("output", "output"));
    if (run.has("restartRead")) {
        const std::string rd = resolvePath(dir, run["restartRead"].string());
        if (!restart::read(rd, flow)) throw std::runtime_error("restart data not found in " + rd);
        std::cout << "restarted from " << rd << " at t = " << flow.time().time << '\n';
    }
    const std::string restartDir = run.has("restartWrite") ? resolvePath(dir, run["restartWrite"].string()) : "";

    VtkWriter vtk(mesh, out);
    vtk.add(flow.U());
    vtk.add(flow.p());
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
            if (probes && it % probeInterval == 0) probes->sample(it);
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
            dt = std::min(dt, endTime - flow.time().time);
            // 落在写出时刻上
            if (nextWrite < GREAT && flow.time().time + dt > nextWrite * (1 + 1e-12) - 1e-15)
                dt = std::max(nextWrite - flow.time().time, 1e-12 * dt);
            const StepInfo info = flow.step(dt);
            ++n;
            if (n % printInterval == 0)
                std::cout << "t " << flow.time().time << "  dt " << dt << "  Co " << info.maxCo << "  cont "
                          << info.continuityError << '\n';
            if (probes && n % probeInterval == 0) probes->sample(flow.time().time);
            if (flow.time().time >= nextWrite * (1 - 1e-12)) {
                if (writeVtk) vtk.write(flow.time().time);
                nextWrite += writeInterval;
            }
        }
        std::cout << "reached t = " << flow.time().time << " in " << n << " steps, " << par::wallTime() - t0 << " s\n";
        if (writeVtk && writeInterval <= 0) vtk.write(flow.time().time);
    }
    writeLines();
    if (run.get("tecplot", false)) {
        TecplotWriter tp(mesh);
        tp.add(flow.U());
        tp.add(flow.p());
        tp.write(out + "/fields.dat", "U_p");
    }
    if (!restartDir.empty()) restart::write(restartDir, flow);
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
