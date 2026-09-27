// 顶盖驱动方腔（稳态 SIMPLE），与 Ghia et al. (1982) 中心线速度对比。
//   --mesh <polyMesh 目录>   读 OpenFOAM 网格（方腔 [0,1]²，顶盖 patch 名 --lid），否则生成 --N×N 均匀网格
//   --Re 100|400|1000...      U_lid = 1，ν = 1/Re
#include "fvm/core/Args.h"
#include "fvm/io/Probes.h"
#include "fvm/io/VtkWriter.h"
#include "fvm/solvers/IncompressibleFlow.h"

#include <iostream>

using namespace cfd;

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    Args args(argc, argv);
    const scalar Re = args.get("Re", 100.0);
    const std::string meshDir = args.get("mesh", "");
    const std::string lid = args.get("lid", meshDir.empty() ? "yMax" : "topWalls");
    const int maxIter = args.get("maxIter", 5000);
    const scalar tol = args.get("tol", 1e-6);
    const std::string out = args.get("out", "cavity_out");

    RawMesh raw;
    if (par::master()) {
        if (!meshDir.empty()) {
            raw = readPolyMesh(meshDir);
        } else {
            const int N = args.get("N", 64);
            BoxSpec s;
            s.n[0] = s.n[1] = N;
            s.n[2] = 1;
            s.hi = {1, 1, 1.0 / N};
            s.twoD = true;
            for (auto& t : s.types) t = PatchType::Wall;
            raw = generateBox(s);
        }
    }
    auto mesh = Mesh::build(par::master() ? &raw : nullptr);
    mesh->printSummary();

    IncompressibleFlow flow(mesh, 1.0 / Re);
    flow.controls.steady = true;
    flow.controls.alphaU = args.get("alphaU", 0.7);
    flow.controls.alphaP = args.get("alphaP", 0.3);
    flow.controls.nCorr = 1;
    flow.controls.nNonOrthCorr = args.get("nNonOrth", 1);
    const std::string pc = args.get("pPrecond", "GAMG");
    flow.controls.UControls = {"PBiCGStab", "DILU", 1e-10, 0.1, 1000, 0, 1};
    flow.controls.pControls = {"PCG", pc, 1e-10, 0.05, 2000, 0, 1};
    flow.controls.pFinalControls = flow.controls.pControls;
    flow.controls.verbose = args.flag("verbose");
    flow.divScheme = ConvectionScheme::parse(args.get("div", "linearUpwind"));

    for (const auto& p : mesh->patches()) {
        if (p.name == lid)
            flow.U().fixedValue(p.name, Vec3{1, 0, 0});
        else
            flow.U().fixedValue(p.name, Vec3{0, 0, 0});
        flow.p().zeroGradient(p.name);
    }
    flow.initialize();

    int it = 0;
    StepInfo info;
    const scalar t0 = par::wallTime();
    for (it = 1; it <= maxIter; ++it) {
        info = flow.step(1.0);
        if (it % 100 == 0 || it == 1)
            std::cout << "iter " << it << "  U res " << info.UInitialResidual << "  p res " << info.pInitialResidual
                      << "  cont " << info.continuityError << '\n';
        if (info.UInitialResidual < tol && info.pInitialResidual < tol) break;
    }
    std::cout << "finished after " << std::min(it, maxIter) << " iterations, " << par::wallTime() - t0 << " s\n";

    VtkWriter vtk(mesh, out);
    vtk.add(flow.U());
    vtk.add(flow.p());
    vtk.write(1);
    // 竖直中心线 x=0.5 上的 u，水平中心线 y=0.5 上的 v
    std::vector<Vec3> vline, hline;
    const scalar zc = mesh->C().empty() ? 0.0 : mesh->C()[0].z;
    for (int k = 0; k <= 200; ++k) {
        vline.push_back({0.5, k / 200.0, zc});
        hline.push_back({k / 200.0, 0.5, zc});
    }
    Probes(mesh, vline, out).writeTable(out + "/centreline_x0.5.csv", flow.U());
    Probes(mesh, hline, out).writeTable(out + "/centreline_y0.5.csv", flow.U());
    std::cout << "checksum U " << checksum(flow.U()) << " p " << checksum(flow.p()) << '\n';
    args.warnUnused();
    return 0;
}
