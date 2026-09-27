// 二维 Taylor-Green 涡衰减：周期域 [0,2π]²，与解析解比较验证时间/空间精度
//   u =  sin x cos y F,  v = −cos x sin y F,  p = (cos2x + cos2y) F²/4,  F = exp(−2νt)
//
// 用法：tgv2d --N 64 --nu 0.01 --endTime 1 --Co 0.2 --div linear --ddt backward
//       [--vtk 1] [--restartWrite dir] [--restartRead dir] [--steps n]
#include "fvm/core/Args.h"
#include "fvm/io/Probes.h"
#include "fvm/io/Restart.h"
#include "fvm/io/VtkWriter.h"
#include "fvm/solvers/IncompressibleFlow.h"

#include <cmath>
#include <iostream>

using namespace cfd;

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    Args args(argc, argv);
    const int N = args.get("N", 32);
    const scalar nu = args.get("nu", 0.01);
    const scalar endTime = args.get("endTime", 1.0);
    const scalar Co = args.get("Co", 0.2);
    const std::string out = args.get("out", "tgv2d_out");

    RawMesh raw;
    if (par::master()) {
        BoxSpec s;
        s.n[0] = N;
        s.n[1] = N;
        s.n[2] = 1;
        s.hi = {2 * PI, 2 * PI, 2 * PI / N};
        s.periodic[0] = s.periodic[1] = true;
        s.twoD = true;
        raw = generateBox(s);
    }
    auto mesh = Mesh::build(par::master() ? &raw : nullptr);
    mesh->printSummary();

    IncompressibleFlow flow(mesh, nu);
    flow.divScheme = ConvectionScheme::parse(args.get("div", "linear"));
    flow.ddtScheme = parseDdtScheme(args.get("ddt", "backward"));
    flow.controls.nCorr = args.get("nCorr", 2);
    flow.controls.ddtPhiCoeff = args.get("ddtPhiCoeff", 0.0);
    flow.controls.pControls = {"PCG", args.get("pPrecond", "GAMG"), 1e-10, 0.0, 5000, 0, 1};
    flow.controls.pFinalControls = flow.controls.pControls;
    flow.controls.UControls = {"PBiCGStab", "DILU", 1e-10, 0.0, 1000, 0, 1};

    auto exact = [&](const Vec3& x, scalar t, Vec3& u, scalar& p) {
        const scalar F = std::exp(-2 * nu * t);
        u = {std::sin(x.x) * std::cos(x.y) * F, -std::cos(x.x) * std::sin(x.y) * F, 0};
        p = 0.25 * (std::cos(2 * x.x) + std::cos(2 * x.y)) * F * F;
    };
    for (label c = 0; c < mesh->nTotalCells(); ++c) exact(mesh->C()[c], 0, flow.U()[c], flow.p()[c]);
    flow.initialize();
    if (args.has("restartRead")) {
        if (!restart::read(args.get("restartRead", ""), flow)) throw std::runtime_error("restart not found");
        std::cout << "restarted at t = " << flow.time().time << '\n';
    }

    const scalar dx = 2 * PI / N;
    const int nTotal = int(std::ceil(endTime / (Co * dx) - 1e-9));
    const scalar dt = endTime / nTotal;
    const int nSteps = args.get("steps", nTotal - flow.time().timeIndex);

    VtkWriter vtk(mesh, out);
    vtk.add(flow.U());
    vtk.add(flow.p());
    Probes probes(mesh, {{PI / 2, 0.1, 0}, {1.0, 2.0, 0}}, out + "/probes");
    probes.add(flow.U());
    probes.add(flow.p());
    const bool writeVtk = args.flag("vtk");
    if (writeVtk) vtk.write(flow.time().time);

    StepInfo info;
    for (int s = 0; s < nSteps; ++s) {
        info = flow.step(dt);
        probes.sample(flow.time().time);
    }
    if (writeVtk) vtk.write(flow.time().time);
    if (args.has("restartWrite")) restart::write(args.get("restartWrite", ""), flow);

    par::SumAcc eu, nu2, pm, pe, ep;
    for (label c = 0; c < mesh->nCells(); ++c) {
        Vec3 u;
        scalar p;
        exact(mesh->C()[c], flow.time().time, u, p);
        pm.add(flow.p()[c] * mesh->V()[c]);
        pe.add(p * mesh->V()[c]);
    }
    const scalar pShift = (pm.allReduce() - pe.allReduce()) / mesh->totalVolume();
    for (label c = 0; c < mesh->nCells(); ++c) {
        Vec3 u;
        scalar p;
        exact(mesh->C()[c], flow.time().time, u, p);
        eu.add(magSqr(flow.U()[c] - u) * mesh->V()[c]);
        nu2.add(magSqr(u) * mesh->V()[c]);
        ep.add(std::pow(flow.p()[c] - pShift - p, 2) * mesh->V()[c]);
    }
    const scalar errU = std::sqrt(eu.allReduce() / nu2.allReduce());
    const scalar errP = std::sqrt(ep.allReduce() / mesh->totalVolume());
    std::cout.precision(10);
    std::cout << "N " << N << " steps " << flow.time().timeIndex << " t " << flow.time().time << " relL2(U) " << errU
              << " L2(p) " << errP << " contErr " << info.continuityError << " Co " << info.maxCo << '\n';
    const auto cu = checksum(flow.U()), cp = checksum(flow.p());
    std::cout << "checksum U " << cu << " p " << cp << '\n';
    args.warnUnused();
    return 0;
}
