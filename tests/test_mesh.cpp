#include "fvm/mesh/Mesh.h"
#include <cassert>
#include <iostream>
#include <random>

using namespace cfd;

static int fails = 0;
#define CHECK(c) do { if (!(c)) { std::cerr << "FAIL " #c " at line " << __LINE__ << std::endl; ++fails; } } while (0)

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    // ExactSum：顺序无关
    {
        std::mt19937 g(1);
        std::uniform_real_distribution<double> u(-1, 1);
        std::vector<double> v(10000);
        for (auto& x : v) x = u(g) * std::pow(10.0, int(u(g) * 20));
        par::ExactSum a, b;
        for (auto x : v) a.add(x);
        for (auto it = v.rbegin(); it != v.rend(); ++it) b.add(*it);
        CHECK(a.value() == b.value());
        par::ExactSum c;
        c.add(1e20); c.add(1.0); c.add(-1e20);
        CHECK(c.value() == 1.0);
        par::ExactSum d; d.add(-3.5); d.add(1.25);
        CHECK(d.value() == -2.25);
    }
    // 周期盒子
    for (int twoD = 0; twoD < 2; ++twoD) {
        RawMesh raw;
        if (par::master()) {
            BoxSpec s;
            s.n[0] = 8; s.n[1] = 6; s.n[2] = twoD ? 1 : 5;
            s.hi = {2.0, 1.0, 0.5};
            s.periodic[0] = true; s.periodic[2] = !twoD;
            s.twoD = twoD;
            s.stretch[1] = tanhStretch(1.5);
            raw = generateBox(s);
        }
        auto m = Mesh::build(par::master() ? &raw : nullptr);
        m->printSummary();
        CHECK(std::abs(m->totalVolume() - 1.0) < 1e-12);
        // 每个自有单元的外法向面积之和为零（3D 闭合；2D 去掉 empty 后 z 分量不闭合）
        std::vector<Vec3> sum(m->nTotalCells());
        for (label f = 0; f < m->nFaces(); ++f) {
            sum[m->owner()[f]] += m->Sf()[f];
            if (f < m->nInternalFaces()) sum[m->neighbour()[f]] -= m->Sf()[f];
        }
        scalar err = 0;
        for (label c = 0; c < m->nCells(); ++c) {
            Vec3 s = sum[c];
            if (twoD) s[m->emptyDir()] = 0;
            err = std::max(err, mag(s));
        }
        CHECK(err < 1e-12);
        // 周期面的 d 应接近相邻格距
        scalar maxd = 0;
        for (label f = 0; f < m->nInternalFaces(); ++f) maxd = std::max(maxd, mag(m->d()[f]));
        CHECK(maxd < 0.5);
        CHECK(m->patches().size() == (twoD ? 2u : 2u));
    }
    // 读 OpenFOAM 网格
    {
        RawMesh raw;
        if (par::master()) raw = readPolyMesh("meshes/cavity2D_tri/polyMesh");
        auto m = Mesh::build(par::master() ? &raw : nullptr);
        m->printSummary();
        CHECK(m->twoD());
    }
    if (par::master()) std::cout << (fails ? "FAILED" : "ALL PASSED") << std::endl;
    return fails ? 1 : 0;
}
