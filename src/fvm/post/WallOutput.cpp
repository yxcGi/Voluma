#include "fvm/post/WallOutput.h"

#include "fvm/models/Energy.h"
#include "fvm/models/TurbulenceModel.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace cfd {

WallOutput::WallOutput(const IncompressibleFlow& flow, const TurbulenceModel* turb, const EnergyModel* energy,
                       const std::vector<std::string>& patches)
    : flow_(flow), turb_(turb), energy_(energy) {
    const Mesh& m = flow.mesh();
    if (patches.empty()) {
        for (label p = 0; p < label(m.patches().size()); ++p)
            if (m.patches()[p].type == PatchType::Wall) patches_.push_back(p);
    } else {
        for (const auto& n : patches) {
            const label p = m.findPatch(n);
            if (p < 0) throw std::runtime_error("wallOutput: no patch " + n);
            patches_.push_back(p);
        }
    }
    const label nB = m.nBoundaryFaces();
    p_.assign(nB, 0.0);
    Cp_.assign(nB, 0.0);
    tauMag_.assign(nB, 0.0);
    Cf_.assign(nB, 0.0);
    yPlus_.assign(nB, 0.0);
    q_.assign(nB, 0.0);
    Tw_.assign(nB, 0.0);
    tau_.assign(nB, Vec3{});
}

void WallOutput::update() {
    const Mesh& m = flow_.mesh();
    const label nI = m.nInternalFaces();
    const auto nuEff = flow_.nuEffFaces();
    const std::vector<scalar> ypTurb = turb_ ? turb_->yPlusBoundaryFaces() : std::vector<scalar>{};
    const std::vector<scalar> qE = energy_ ? energy_->heatFluxBoundaryFaces() : std::vector<scalar>{};
    const scalar qDyn = 0.5 * Uref * Uref;
    const bool proj = mag(flowDir) > 0;
    const Vec3 fd = proj ? flowDir / mag(flowDir) : Vec3{};
    const auto& U = flow_.U();
    const auto& p = flow_.p();
    for (label pi : patches_) {
        const Patch& pt = m.patches()[pi];
        for (label f = pt.start; f < pt.end(); ++f) {
            const label bf = f - nI, o = m.owner()[f];
            const Vec3 n = m.Sf()[f] / m.magSf()[f];
            Vec3 du = U[o] - U.bValue(f);
            du -= dot(du, n) * n;
            if (m.twoD()) du[m.emptyDir()] = 0;
            const Vec3 tau = nuEff[f] * m.nonOrthDeltaCoeffs()[f] * du;
            tau_[bf] = tau;
            tauMag_[bf] = mag(tau);
            p_[bf] = p.bValue(f);
            Cp_[bf] = (p_[bf] - pRef) / qDyn;
            Cf_[bf] = (proj ? dot(tau, fd) : mag(tau)) / qDyn;
            scalar yp = turb_ ? ypTurb[bf] : -1.0;
            if (yp < 0) {
                const scalar y = std::max(dot(m.Cf()[f] - m.C()[o], n), SMALL);
                yp = y * std::sqrt(mag(tau)) / flow_.nu();
            }
            yPlus_[bf] = yp;
            if (energy_) {
                q_[bf] = qE[bf];
                Tw_[bf] = energy_->T().bValue(f);
            }
        }
    }
}

void WallOutput::writeCsv(const std::string& dir) const {
    const Mesh& m = flow_.mesh();
    const label nI = m.nInternalFaces();
    if (par::master()) std::filesystem::create_directories(dir);
    struct Rec {
        double v[17];
    };
    for (label pi : patches_) {
        const Patch& pt = m.patches()[pi];
        std::vector<Rec> loc(pt.size);
        for (label i = 0; i < pt.size; ++i) {
            const label f = pt.start + i, bf = f - nI;
            const Vec3& c = m.Cf()[f];
            const Vec3 n = m.Sf()[f] / m.magSf()[f];
            loc[i] = {{c.x, c.y, c.z, n.x, n.y, n.z, m.magSf()[f], p_[bf], Cp_[bf], tau_[bf].x, tau_[bf].y, tau_[bf].z,
                       tauMag_[bf], Cf_[bf], yPlus_[bf], q_[bf], Tw_[bf]}};
        }
        const auto all = pt.ordering.gather(loc.data());
        if (!par::master()) continue;
        std::ofstream os(dir + "/wall_" + pt.name + ".csv");
        os << "x,y,z,nx,ny,nz,area,p,Cp,tau_x,tau_y,tau_z,tauMag,Cf,yPlus";
        if (energy_) os << ",q,Tw";
        os << '\n' << std::setprecision(10);
        const int nv = energy_ ? 17 : 15;
        for (const auto& r : all) {
            for (int k = 0; k < nv; ++k) os << (k ? "," : "") << r.v[k];
            os << '\n';
        }
    }
}

void WallOutput::writeVtk(const std::string& dir, scalar time) {
    if (!surf_) {
        surf_ = std::make_unique<SurfaceWriter>(flow_.meshPtr(), dir, "surface");
        for (label p : patches_) surf_->addPatch(p);
        surf_->add("p", p_);
        surf_->add("Cp", Cp_);
        surf_->add("wallShearStress", tau_);
        surf_->add("Cf", Cf_);
        surf_->add("yPlus", yPlus_);
        if (energy_) {
            surf_->add("q", q_);
            surf_->add("Tw", Tw_);
        }
    }
    surf_->write(time);
}

} // namespace cfd
