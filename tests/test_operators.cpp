// 离散算子精度测试（有解析解）：
//   1. Poisson（制造解）在均匀盒子网格上二阶收敛
//   2. 同一问题在非结构三角形（棱柱）网格上：误差随非正交修正迭代收敛到小量
//   3. 一维对流扩散（Pe=10）：linear 二阶、upwind 一阶
//   4. Robin 边界：线性解应精确再现
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

MeshPtr box2D(int nx, int ny, Vec3 hi) {
    RawMesh raw;
    if (par::master()) {
        BoxSpec s;
        s.n[0] = nx;
        s.n[1] = ny;
        s.n[2] = 1;
        s.hi = hi;
        s.twoD = true;
        raw = generateBox(s);
    }
    return Mesh::build(par::master() ? &raw : nullptr);
}

// 体积加权 L2 误差
scalar l2(const Mesh& m, const ScalarField& T, const std::function<scalar(const Vec3&)>& ex) {
    par::SumAcc e;
    for (label c = 0; c < m.nCells(); ++c) e.add(std::pow(T[c] - ex(m.C()[c]), 2) * m.V()[c]);
    return std::sqrt(e.allReduce() / m.totalVolume());
}

scalar poisson(MeshPtr mesh, scalar L, bool nonOrthCorr) {
    const scalar k = PI / L, G = 2.0;
    auto ex = [&](const Vec3& x) { return std::sin(k * x.x) * std::sin(k * x.y); };
    ScalarTransport st(mesh, "T", G);
    st.steady = true;
    st.laplacianNonOrthCorr = nonOrthCorr;
    st.controls = {"PCG", "GAMG", 1e-13, 0.0, 2000, 0, 1};
    for (const auto& p : mesh->patches()) st.T().fixedValue(p.name, 0.0);
    st.Q().assign(mesh->nTotalCells(), 0.0);
    for (label c = 0; c < mesh->nTotalCells(); ++c) st.Q()[c] = 2 * k * k * G * ex(mesh->C()[c]);
    st.initialize();
    for (int it = 0; it < 30; ++it)
        if (st.step(1.0, nullptr).initialResidual < 1e-11) break;
    return l2(*mesh, st.T(), ex);
}

} // namespace

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    // 1. Poisson，盒子网格
    {
        const scalar e1 = poisson(box2D(16, 16, {1, 1, 0.1}), 1.0, true);
        const scalar e2 = poisson(box2D(32, 32, {1, 1, 0.1}), 1.0, true);
        const scalar order = std::log2(e1 / e2);
        std::cout << "poisson box: e16 " << e1 << " e32 " << e2 << " order " << order << '\n';
        CHECK(order > 1.9 && order < 2.2);
    }
    // 2. Poisson，三角形网格（单位方腔，棱柱，最大非正交 30°）
    {
        RawMesh raw;
        if (par::master()) raw = readPolyMesh("meshes/cavity2D_tri/polyMesh");
        auto m = Mesh::build(par::master() ? &raw : nullptr);
        const scalar e0 = poisson(m, 1.0, false);
        const scalar e3 = poisson(m, 1.0, true);
        std::cout << "poisson tri: without non-orthogonal correction " << e0 << ", with " << e3 << '\n';
        CHECK(e3 < 2e-3 && e3 < e0);
    }
    // 3. 一维对流扩散 T(0)=0, T(1)=1, u=1, Γ=0.1
    {
        const scalar Pe = 10.0;
        auto ex = [&](const Vec3& x) { return (std::exp(Pe * x.x) - 1) / (std::exp(Pe) - 1); };
        for (const char* scheme : {"linear", "upwind"}) {
            scalar err[2];
            for (int r = 0; r < 2; ++r) {
                auto mesh = box2D(r ? 80 : 40, 1, {1, 0.1, 0.1});
                ScalarTransport st(mesh, "T", 1.0 / Pe);
                st.steady = true;
                st.divScheme = ConvectionScheme::parse(scheme);
                st.T().fixedValue("xMin", 0.0);
                st.T().fixedValue("xMax", 1.0);
                st.T().zeroGradient("yMin");
                st.T().zeroGradient("yMax");
                std::vector<scalar> phi(mesh->nFaces());
                for (label f = 0; f < mesh->nFaces(); ++f) phi[f] = mesh->Sf()[f].x;
                st.initialize();
                for (int it = 0; it < 5; ++it) st.step(1.0, &phi);
                err[r] = l2(*mesh, st.T(), ex);
            }
            const scalar order = std::log2(err[0] / err[1]);
            std::cout << "convection-diffusion " << scheme << ": e40 " << err[0] << " e80 " << err[1] << " order "
                      << order << '\n';
            CHECK(std::string(scheme) == "linear" ? order > 1.8 : (order > 0.7 && order < 1.3));
        }
    }
    // 4. Robin：a T + b ∂T/∂n = c 于 x=1，T(0)=0，解 T = x c/(a+b)
    {
        auto mesh = box2D(10, 2, {1, 1, 0.1});
        ScalarTransport st(mesh, "T", 1.0);
        st.steady = true;
        const scalar a = 2, b = 0.5, c = 3;
        st.T().fixedValue("xMin", 0.0);
        st.T().set<RobinBC<scalar>>("xMax", a, b, c);
        st.T().zeroGradient("yMin");
        st.T().zeroGradient("yMax");
        st.initialize();
        for (int it = 0; it < 3; ++it) st.step(1.0, nullptr);
        const scalar e = l2(*mesh, st.T(), [&](const Vec3& x) { return x.x * c / (a + b); });
        std::cout << "robin: error " << e << '\n';
        CHECK(e < 1e-9);
    }
    if (par::master()) std::cout << (fails ? "FAILED" : "ALL PASSED") << std::endl;
    return fails ? 1 : 0;
}
