#include "fvm/mesh/PeriodicHill.h"

#include <algorithm>
#include <cmath>

namespace cfd {

namespace {

// ERCOFTAC 多项式：山高 28 mm，自变量 x ∈ [0, 54] mm（山顶 x = 0），返回 mm
double hillMM(double x) {
    double y;
    if (x < 9.0) {
        y = std::min(28.0, 2.8e1 + 6.775070969851e-03 * x * x - 2.124527775800e-03 * x * x * x);
    } else if (x < 14.0) {
        y = 2.507355893131e+01 + 9.754803562315e-01 * x - 1.016116352781e-01 * x * x + 1.889794677828e-03 * x * x * x;
    } else if (x < 20.0) {
        y = 2.579601052357e+01 + 8.206693007457e-01 * x - 9.055370274339e-02 * x * x + 1.626510569859e-03 * x * x * x;
    } else if (x < 30.0) {
        y = 4.046435022819e+01 - 1.379581654948e+00 * x + 1.945884504128e-02 * x * x - 2.070318932190e-04 * x * x * x;
    } else if (x < 40.0) {
        y = 1.792461334664e+01 + 8.743920332081e-01 * x - 5.567361123058e-02 * x * x + 6.277731764683e-04 * x * x * x;
    } else if (x <= 54.0) {
        y = std::max(0.0, 5.639011190988e+01 - 2.010520359035e+00 * x + 1.644919857549e-02 * x * x +
                              2.674976141766e-05 * x * x * x);
    } else {
        y = 0.0;
    }
    return y;
}

} // namespace

scalar periodicHillHeight(scalar x, scalar H, scalar Lx) {
    const scalar L = Lx * H;
    x = std::fmod(x, L);
    if (x < 0) x += L;
    const double s = 28.0 / H;  // 物理长度 → mm
    const double xl = x * s, xr = (L - x) * s;
    return H / 28.0 * (xl <= 54.0 ? hillMM(xl) : xr <= 54.0 ? hillMM(xr) : 0.0);
}

RawMesh generatePeriodicHillMesh(const PeriodicHillSpec& p) {
    BoxSpec s;
    for (int k = 0; k < 3; ++k) s.n[k] = p.n[k];
    const scalar Lx = p.length * p.H, Ly = p.span * p.H, top = p.height * p.H;
    s.lo = {0, 0, 0};
    s.hi = {Lx, Ly, 1};
    s.periodic[0] = s.periodic[1] = true;
    s.names[4] = "bottom";
    s.names[5] = "top";
    s.types[4] = s.types[5] = PatchType::Wall;
    if (p.stretch > 0) s.stretch[2] = tanhStretch(p.stretch);
    const scalar H = p.H, len = p.length;
    s.map = [H, len, top](const Vec3& x) {
        const scalar h = periodicHillHeight(x.x, H, len);
        return Vec3{x.x, x.y, h + x.z * (top - h)};
    };
    return generateBox(s);
}

} // namespace cfd
