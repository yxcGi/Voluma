// 翼型 C 网格：体积守恒（计算域面积 − 翼型面积）、单元有效、边界齐全，攻角 0° 与 29°
#include "fvm/mesh/AirfoilMesh.h"
#include "fvm/mesh/Mesh.h"

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

int main(int argc, char** argv) {
    par::Environment env(argc, argv);
    for (scalar alpha : {0.0, 29.0}) {
        AirfoilMeshSpec s;
        s.alphaDeg = alpha;
        s.nAirfoil = 160;
        s.nWake = 60;
        s.nNormal = 60;
        RawMesh raw;
        if (par::master()) raw = generateAirfoilCMesh(s);
        auto m = Mesh::build(par::master() ? &raw : nullptr);
        m->printSummary();
        // NACA0012（闭合尾缘）面积 ≈ 0.0822 c²；多边形近似略小
        const scalar domain = (s.upstream + s.downstream) * 2 * s.halfHeight * s.depth;
        const scalar airfoilArea = (domain - m->totalVolume()) / s.depth;
        std::cout << "alpha " << alpha << ": airfoil area " << airfoilArea << '\n';
        CHECK(std::abs(airfoilArea - 0.0822) < 1e-3);
        scalar vmin = GREAT;
        for (label c = 0; c < m->nCells(); ++c) vmin = std::min(vmin, m->V()[c]);
        CHECK(par::allMin(vmin) > 0);
        CHECK(m->findPatch("airfoil") >= 0 && m->findPatch("inlet") >= 0 && m->findPatch("outlet") >= 0 &&
              m->findPatch("top") >= 0 && m->findPatch("bottom") >= 0);
        // 翼型表面面数 = nAirfoil
        const auto& ap = m->patches()[m->findPatch("airfoil")];
        CHECK(ap.globalSize == s.nAirfoil);
    }
    if (par::master()) std::cout << (fails ? "FAILED" : "ALL PASSED") << std::endl;
    return fails ? 1 : 0;
}
