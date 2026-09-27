// fvmScalar：标量输运/导热求解程序（给定速度场或纯扩散），由 JSON 算例文件驱动。
//
//   ∂T/∂t + ∇·(U T) − ∇·(Γ∇T) = Q
//
//   mpirun -np 4 fvmScalar cases/heatSource/case.json
#include "fvm/app/CaseSetup.h"
#include "fvm/io/Probes.h"
#include "fvm/io/TecplotWriter.h"
#include "fvm/io/VtkWriter.h"
#include "fvm/solvers/ScalarTransport.h"

#include <iostream>

using namespace cfd;

namespace {

void runCase(const std::string& caseFile) {
    const Json cfg = Json::parseFile(caseFile);
    const std::string dir = caseDirectory(caseFile);
    auto mesh = buildMesh(cfg["mesh"], dir);
    mesh->printSummary();

    const Json& phys = cfg["physics"];
    const std::string name = phys.get("field", "T");
    ScalarTransport st(mesh, name, phys["diffusivity"].number());
    st.steady = phys.get("steady", true);
    if (phys.has("source")) st.setUniformSource(phys["source"].number());

    std::vector<scalar> phi;
    const bool convect = phys.has("velocity");
    if (convect) {
        const Vec3 U = phys["velocity"].vec3();
        phi.resize(mesh->nFaces());
        for (label f = 0; f < mesh->nFaces(); ++f) phi[f] = dot(U, mesh->Sf()[f]);
    }

    const Json& sch = cfg["schemes"];
    st.divScheme = ConvectionScheme::parse(sch.get("div", "linearUpwind"));
    st.ddtScheme = parseDdtScheme(sch.get("ddt", "backward"));
    st.laplacianNonOrthCorr = sch.get("laplacianNonOrthCorr", true);
    const Json& sol = cfg["solution"];
    st.nNonOrthCorr = sol.get("nNonOrthCorr", 0);
    st.relax = sol.get("relax", 1.0);
    st.verbose = sol.get("verbose", false);
    st.controls = parseSolverControls(sol[name], st.controls);

    setBoundaryConditions(st.T(), cfg["boundary"], Json());
    st.T().setUniform(cfg["initial"].get(name, 0.0));
    st.initialize();

    const Json& run = cfg["run"];
    const std::string out = resolvePath(dir, run.get("output", "output"));
    VtkWriter vtk(mesh, out);
    vtk.add(st.T());
    const double t0 = par::wallTime();
    if (st.steady) {
        const int maxIter = run.get("maxIter", 1000);
        const scalar tol = run.get("residualTol", 1e-8);
        int it = 1;
        for (; it <= maxIter; ++it) {
            const auto perf = st.step(1.0, convect ? &phi : nullptr);
            if (it % run.get("printInterval", 10) == 0 || perf.initialResidual < tol)
                std::cout << "iter " << it << "  " << name << " res " << perf.initialResidual << '\n';
            if (perf.initialResidual < tol) break;
        }
        std::cout << "finished after " << std::min(it, maxIter) << " iterations, " << par::wallTime() - t0 << " s\n";
    } else {
        const scalar endTime = run["endTime"].number();
        const scalar dt = run["dt"].number();
        const int n = int(std::ceil(endTime / dt - 1e-9));
        for (int s = 0; s < n; ++s) st.step(dt, convect ? &phi : nullptr);
        std::cout << "reached t = " << st.time().time << " in " << n << " steps, " << par::wallTime() - t0 << " s\n";
    }
    if (run.get("vtk", true)) vtk.write(st.time().time);
    if (run.get("tecplot", false)) {
        TecplotWriter tp(mesh);
        tp.add(st.T());
        tp.write(out + "/" + name + ".dat", name);
    }
    if (cfg.has("lines")) {
        for (auto& l : cfg["lines"].array()) {
            const Vec3 a = l["start"].vec3(), b = l["end"].vec3();
            const int np = l.get("n", 100);
            std::vector<Vec3> pts;
            for (int k = 0; k < np; ++k) pts.push_back(a + (b - a) * (np > 1 ? scalar(k) / (np - 1) : 0.0));
            Probes line(mesh, pts, out);
            line.writeTable(out + "/line_" + l.get("name", "line") + "_" + name + ".csv", st.T());
        }
    }
    par::SumAcc mn;
    scalar tmax = -GREAT, tmin = GREAT;
    for (label c = 0; c < mesh->nCells(); ++c) {
        mn.add(st.T()[c] * mesh->V()[c]);
        tmax = std::max(tmax, st.T()[c]);
        tmin = std::min(tmin, st.T()[c]);
    }
    std::cout << name << " mean " << mn.allReduce() / mesh->totalVolume() << " min " << par::allMin(tmin) << " max "
              << par::allMax(tmax) << '\n';
    std::cout << "checksum " << name << ' ' << checksum(st.T()) << '\n';
    if (par::master()) cfg.reportUnused(std::cout);
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    if (argc < 2) {
        std::cout << "usage: fvmScalar case.json\n";
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
