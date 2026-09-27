#include "fvm/post/Forces.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace cfd {

Forces::Forces(const IncompressibleFlow& flow, const std::vector<std::string>& patches) : flow_(flow) {
    for (const auto& n : patches) {
        const label p = flow.mesh().findPatch(n);
        if (p < 0) throw std::runtime_error("forces: no patch " + n);
        patches_.push_back(p);
    }
}

scalar Forces::span() const {
    const Mesh& m = flow_.mesh();
    if (!m.twoD()) return 1.0;
    const int ed = m.emptyDir();
    scalar lo = GREAT, hi = -GREAT;
    for (const auto& p : m.vtk().points) {
        lo = std::min(lo, p[ed]);
        hi = std::max(hi, p[ed]);
    }
    return par::allMax(hi) - par::allMin(lo);
}

namespace {
// 壁面切应力（运动学，τ/ρ），指向流体对壁面的拖曳方向
Vec3 wallShear(const Mesh& m, const IncompressibleFlow& flow, const std::vector<scalar>& nuEff, label f) {
    const Vec3 n = m.Sf()[f] / m.magSf()[f];
    Vec3 du = flow.U()[m.owner()[f]] - flow.U().bValue(f);
    du -= dot(du, n) * n;
    if (m.twoD()) du[m.emptyDir()] = 0;
    return nuEff[f] * m.nonOrthDeltaCoeffs()[f] * du;
}
} // namespace

ForceResult Forces::compute() const {
    const Mesh& m = flow_.mesh();
    const auto nuEff = flow_.nuEffFaces();
    par::SumAcc acc[12];
    for (label p : patches_) {
        const Patch& pt = m.patches()[p];
        for (label i = 0; i < pt.size; ++i) {
            const label f = pt.start + i;
            const Vec3 Fp = rho * (flow_.p().bValue(f) - pRef) * m.Sf()[f];
            const Vec3 Fv = rho * m.magSf()[f] * wallShear(m, flow_, nuEff, f);
            const Vec3 r = m.Cf()[f] - CofR;
            const Vec3 Mp = cross(r, Fp), Mv = cross(r, Fv);
            for (int k = 0; k < 3; ++k) {
                acc[k].add(Fp[k]);
                acc[3 + k].add(Fv[k]);
                acc[6 + k].add(Mp[k]);
                acc[9 + k].add(Mv[k]);
            }
        }
    }
    double v[12];
    par::SumAcc::allReduce(acc, 12, v);
    ForceResult r;
    r.pressure = {v[0], v[1], v[2]};
    r.viscous = {v[3], v[4], v[5]};
    r.momentP = {v[6], v[7], v[8]};
    r.momentV = {v[9], v[10], v[11]};
    return r;
}

void Forces::writeSurface(const std::string& file, scalar Uref, const Vec3& flowDir) const {
    const Mesh& m = flow_.mesh();
    const auto nuEff = flow_.nuEffFaces();
    const scalar q = 0.5 * Uref * Uref;
    std::ofstream os;
    if (par::master()) {
        const auto parent = std::filesystem::path(file).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        os.open(file);
        os << "patch,x,y,z,nx,ny,nz,Cp,tau_x,tau_y,tau_z,Cf\n" << std::setprecision(10);
    }
    for (label p : patches_) {
        const Patch& pt = m.patches()[p];
        // 每面 10 个量
        struct Rec {
            double v[10];
        };
        std::vector<Rec> loc(pt.size);
        for (label i = 0; i < pt.size; ++i) {
            const label f = pt.start + i;
            const Vec3 n = m.Sf()[f] / m.magSf()[f];
            const Vec3 tau = wallShear(m, flow_, nuEff, f);
            const Vec3& c = m.Cf()[f];
            loc[i] = {{c.x, c.y, c.z, n.x, n.y, n.z, (flow_.p().bValue(f) - pRef) / q, tau.x, tau.y, tau.z}};
        }
        const auto all = pt.ordering.gather(loc.data());
        if (!par::master()) continue;
        for (const auto& r : all) {
            os << pt.name;
            for (double x : r.v) os << ',' << x;
            os << ',' << dot(Vec3{r.v[7], r.v[8], r.v[9]}, flowDir) / q << '\n';
        }
    }
}

} // namespace cfd
