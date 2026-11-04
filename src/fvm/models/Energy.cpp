#include "fvm/models/Energy.h"

#include "fvm/app/CaseSetup.h"
#include "fvm/discretization/Fvc.h"
#include "fvm/discretization/Fvm.h"
#include "fvm/io/Restart.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace cfd {

namespace {
std::string stringOr(const Json& d, const char* key, const std::string& def) {
    if (!d.has(key)) return def;
    const Json& v = d[key];
    if (v.isString()) return v.string();
    if (v.isBool()) return v.boolean() ? def : "none";
    throw std::runtime_error(std::string("energy.") + key + " must be a string");
}
} // namespace

EnergyModel::EnergyModel(IncompressibleFlow& flow, const Json& d, const Json& boundary, const Json& initial)
    : flow_(flow), m_(flow.mesh()), nu_(flow.nu()), T_(flow.meshPtr(), d.get("field", "T")) {
    if (d.has("alpha"))
        alpha_ = d["alpha"].number();
    else
        alpha_ = nu_ / d.get("Pr", 0.71);
    if (!(alpha_ > 0)) throw std::runtime_error("energy: thermal diffusivity must be positive");
    Prt_ = d.get("Prt", 0.85);
    if (d.has("source")) Q_.assign(m_.nTotalCells(), d["source"].number());
    if (d.has("buoyancy")) {
        const Json& b = d["buoyancy"];
        buoyant_ = true;
        g_ = b["g"].vec3();
        beta_ = b["beta"].number();
        TRef_ = b.get("TRef", initial.get(T_.name(), 0.0));
    }
    const std::string wf = stringOr(d, "wallFunction", "jayatilleke");
    if (wf == "jayatilleke" || wf == "Jayatilleke") {
        wallFunction_ = true;
        useKader_ = false;
    } else if (wf == "kader" || wf == "Kader") {
        wallFunction_ = true;
        useKader_ = true;
    } else if (wf == "none") {
        wallFunction_ = false;
    } else {
        throw std::runtime_error("energy.wallFunction must be jayatilleke | kader | none");
    }
    kappa_ = d.get("kappa", 0.41);
    E_ = d.get("E", 9.8);
    if (d.has("nusselt")) {
        nusseltL_ = d["nusselt"].get("L", 1.0);
        nusseltDT_ = d["nusselt"].get("deltaT", 1.0);
    }
    divScheme_ = ConvectionScheme::parse(d.get("div", "linearUpwind"));
    relax_ = d.get("relax", flow.controls.steady ? 0.7 : 1.0);
    controls_ = parseSolverControls(d["solver"], controls_);

    // 边界条件：heatFlux 为本类处理的类型（固定梯度，梯度按当前壁面 α_eff,w 换算），其余交给通用解析；
    // 壁面缺省绝热（零梯度）
    const label nP = label(m_.patches().size());
    heatFluxQ_.assign(nP, std::numeric_limits<scalar>::quiet_NaN());
    isWallPatch_.assign(nP, 0);
    for (label p = 0; p < nP; ++p) {
        const Patch& pt = m_.patches()[p];
        if (pt.type == PatchType::Empty || pt.type == PatchType::Cyclic) continue;
        isWallPatch_[p] = pt.type == PatchType::Wall;
        const Json& pd = boundary[pt.name];
        if (pd.has(T_.name())) {
            const Json& bd = pd[T_.name()];
            if (bd.get("type", "") == "heatFlux") {
                heatFluxQ_[p] = bd["q"].number();
                T_.set<FixedGradientBC<scalar>>(pt.name, 0.0);
            } else {
                setBoundaryCondition(T_, pt.name, bd, &flow.context());
            }
        } else if (pt.type == PatchType::Wall) {
            T_.set<ZeroGradientBC<scalar>>(pt.name);
        } else if (pt.type == PatchType::Symmetry) {
            T_.set<SymmetryBC<scalar>>(pt.name);
        } else {
            throw std::runtime_error("no boundary condition for " + T_.name() + " on patch " + pt.name);
        }
    }
    T_.setUniform(initial.get(T_.name(), 0.0));
    T_.updateBCs(flow.time().time);
    T_.correctBoundaryConditions();
    wallAlpha_.assign(m_.nBoundaryFaces(), -1.0);
    updateBuoyancy();
}

scalar EnergyModel::TplusJayatilleke(scalar yPlus, scalar Pr, scalar Prt, scalar kappa, scalar E) {
    const scalar r = Pr / Prt;
    const scalar P = 9.24 * (std::pow(r, 0.75) - 1.0) * (1.0 + 0.28 * std::exp(-0.007 * r));
    // 热子层厚度 y+_T：Pr y+ = Pr_t [ln(E y+)/κ + P]（牛顿迭代）
    scalar yT = 11.0;
    for (int it = 0; it < 30; ++it) {
        const scalar f = Pr * yT - Prt * (std::log(E * yT) / kappa + P);
        const scalar df = Pr - Prt / (kappa * yT);
        const scalar dy = f / df;
        yT -= dy;
        if (yT < 1e-3) yT = 1e-3;
        if (std::abs(dy) < 1e-10 * yT) break;
    }
    if (yPlus < yT) return Pr * yPlus;
    return Prt * (std::log(E * yPlus) / kappa + P);
}

scalar EnergyModel::TplusKader(scalar yPlus, scalar Pr) {
    const scalar G = 0.01 * std::pow(Pr * yPlus, 4) / (1.0 + 5.0 * Pr * Pr * Pr * yPlus);
    const scalar beta = std::pow(3.85 * std::cbrt(Pr) - 1.3, 2) + 2.12 * std::log(Pr);
    return Pr * yPlus * std::exp(-G) + (2.12 * std::log(1.0 + yPlus) + beta) * std::exp(-1.0 / std::max(G, 1e-300));
}

void EnergyModel::updateWallAlpha() {
    wallAlpha_.assign(m_.nBoundaryFaces(), -1.0);
    const auto& nuEffW = flow_.wallNuEff();
    const auto& U = flow_.U();
    const label nI = m_.nInternalFaces();
    const scalar Pr = nu_ / alpha_;
    for (label p = 0; p < label(m_.patches().size()); ++p) {
        if (!isWallPatch_[p]) continue;
        const Patch& pt = m_.patches()[p];
        for (label f = pt.start; f < pt.end(); ++f) {
            const label bf = f - nI, o = m_.owner()[f];
            scalar aw = alpha_;
            if (wallFunction_ && !nuEffW.empty() && nuEffW[bf] > nu_ * (1 + 1e-8)) {
                // 热壁面函数：由壁面函数的切应力反算 y+，α_eff,w = ν y+ / T+(y+)
                const Vec3 n = m_.Sf()[f] / m_.magSf()[f];
                const scalar y = std::max(dot(m_.Cf()[f] - m_.C()[o], n), SMALL);
                Vec3 du = U[o] - U.bValue(f);
                du -= dot(du, n) * n;
                const scalar uTau = std::sqrt(nuEffW[bf] * mag(du) / y);
                const scalar yPlus = y * uTau / nu_;
                if (yPlus > SMALL) {
                    const scalar Tp = useKader_ ? TplusKader(yPlus, Pr) : TplusJayatilleke(yPlus, Pr, Prt_, kappa_, E_);
                    aw = std::max(nu_ * yPlus / std::max(Tp, SMALL), alpha_);
                }
            }
            wallAlpha_[bf] = aw;
        }
    }
}

std::vector<scalar> EnergyModel::gammaFaces() const {
    const label nT = m_.nTotalCells(), nI = m_.nInternalFaces();
    std::vector<scalar> g(nT, alpha_);
    const auto& nut = flow_.nut();
    if (!nut.empty())
        for (label c = 0; c < nT; ++c) g[c] = alpha_ + nut[c] / Prt_;
    auto gf = fvc::interpolate(m_, g);
    for (label bf = 0; bf < m_.nBoundaryFaces(); ++bf)
        if (wallAlpha_[bf] >= 0) gf[nI + bf] = wallAlpha_[bf];
    return gf;
}

scalar EnergyModel::correct() {
    const bool steady = flow_.controls.steady;
    const TimeState& ts = flow_.time();
    const DdtScheme ddtS = steady ? DdtScheme::Steady : flow_.ddtScheme;
    if (!steady && ts.timeIndex != lastTimeIndex_) {
        T_.storeOld();
        lastTimeIndex_ = ts.timeIndex;
    }
    updateWallAlpha();
    const auto gf = gammaFaces();
    gammaFacesLast_ = gf;
    // 给定热流的壁面：进入流体的热流 q = α_eff,w ∂T/∂n（n 为外法向，Fourier：−q⃗·n_out = α ∇T·n_out）→ 法向梯度 = q / α_eff,w
    for (label p = 0; p < label(m_.patches().size()); ++p) {
        if (std::isnan(heatFluxQ_[p])) continue;
        auto* bc = dynamic_cast<FixedGradientBC<scalar>*>(&T_.bc(p));
        if (!bc) continue;
        const Patch& pt = m_.patches()[p];
        for (label i = 0; i < pt.size; ++i) bc->gradient()[i] = heatFluxQ_[p] / gf[pt.start + i];
    }
    T_.updateBCs(ts.time);
    T_.correctBoundaryConditions();

    FvMatrix<scalar> M = fvm::ddt(T_, ddtS, ts);
    M += fvm::div(flow_.phi(), T_, divScheme_);
    if (steady) {
        // 有界形式：减去 T ∇·φ（未收敛时的连续性误差）
        const auto dv = fvc::div(m_, flow_.phi());
        std::vector<scalar> neg(dv.size());
        for (std::size_t c = 0; c < dv.size(); ++c) neg[c] = -dv[c];
        M += fvm::Sp(neg, T_);
    }
    M -= fvm::laplacian(gf, T_, flow_.laplacianNonOrthCorr);
    if (!Q_.empty()) M -= Q_;
    if (steady) M.relax(relax_);
    const auto perf = M.solve(controls_, flow_.controls.verbose);
    updateBuoyancy();
    return perf.empty() ? 0.0 : perf[0].initialResidual;
}

void EnergyModel::updateBuoyancy() {
    if (!buoyant_) return;
    auto& f = flow_.bodyForce();
    const label nT = m_.nTotalCells();
    f.resize(nT);
    for (label c = 0; c < nT; ++c) f[c] = (-beta_ * (T_[c] - TRef_)) * g_;
}

scalar EnergyModel::faceHeatFlux(label bf, const std::vector<scalar>& gf) const {
    const label nI = m_.nInternalFaces(), f = nI + bf;
    // 找到所在 patch
    for (label p = 0; p < label(m_.patches().size()); ++p) {
        const Patch& pt = m_.patches()[p];
        if (f < pt.start || f >= pt.end()) continue;
        const auto& bc = T_.bc(p);
        const label i = f - pt.start;
        const scalar gradN = bc.gradientInternalCoeff(T_, i) * T_[m_.owner()[f]] + bc.gradientBoundaryCoeff(T_, i);
        return gf[f] * gradN;  // 进入流体为正：−q⃗·n_out = α_eff ∂T/∂n_out
    }
    return 0.0;
}

std::vector<EnergyModel::PatchHeat> EnergyModel::wallHeatFlux() const {
    const auto gf = gammaFacesLast_.empty() ? gammaFaces() : gammaFacesLast_;
    const label nI = m_.nInternalFaces();
    std::vector<PatchHeat> out;
    for (label p = 0; p < label(m_.patches().size()); ++p) {
        if (!isWallPatch_[p]) continue;
        const Patch& pt = m_.patches()[p];
        par::SumAcc s[3];
        for (label f = pt.start; f < pt.end(); ++f) {
            const scalar a = m_.magSf()[f];
            s[0].add(a);
            s[1].add(faceHeatFlux(f - nI, gf) * a);
            s[2].add(T_.bValue(f) * a);
        }
        double r[3];
        par::SumAcc::allReduce(s, 3, r);
        PatchHeat h;
        h.name = pt.name;
        h.area = r[0];
        h.Q = r[1];
        h.qMean = r[0] > 0 ? r[1] / r[0] : 0.0;
        h.TwMean = r[0] > 0 ? r[2] / r[0] : 0.0;
        h.Nu = (nusseltL_ > 0 && nusseltDT_ != 0) ? h.qMean * nusseltL_ / (alpha_ * nusseltDT_) : 0.0;
        out.push_back(h);
    }
    return out;
}

void EnergyModel::printWallHeatFlux(std::ostream& os) const {
    const auto hs = wallHeatFlux();
    if (!par::master()) return;
    for (const auto& h : hs) {
        os << "wall heat flux " << h.name << ": area " << h.area << "  Q " << h.Q << "  q_mean " << h.qMean
           << "  T_w " << h.TwMean;
        if (nusseltL_ > 0) os << "  Nu " << h.Nu;
        os << '\n';
    }
}

std::vector<scalar> EnergyModel::heatFluxBoundaryFaces() const {
    const auto gf = gammaFacesLast_.empty() ? gammaFaces() : gammaFacesLast_;
    const label nI = m_.nInternalFaces();
    std::vector<scalar> q(m_.nBoundaryFaces(), 0.0);
    for (label p = 0; p < label(m_.patches().size()); ++p) {
        if (!isWallPatch_[p]) continue;
        const Patch& pt = m_.patches()[p];
        for (label f = pt.start; f < pt.end(); ++f) q[f - nI] = faceHeatFlux(f - nI, gf);
    }
    return q;
}

void EnergyModel::writeWallHeatFlux(const std::string& dir) const {
    const auto gf = gammaFacesLast_.empty() ? gammaFaces() : gammaFacesLast_;
    const label nI = m_.nInternalFaces();
    for (label p = 0; p < label(m_.patches().size()); ++p) {
        if (!isWallPatch_[p]) continue;
        const Patch& pt = m_.patches()[p];
        struct Rec {
            double v[5];
        };
        std::vector<Rec> loc(pt.size);
        for (label i = 0; i < pt.size; ++i) {
            const label f = pt.start + i;
            const Vec3& c = m_.Cf()[f];
            loc[i] = {{c.x, c.y, c.z, faceHeatFlux(f - nI, gf), T_.bValue(f)}};
        }
        const auto all = pt.ordering.gather(loc.data());
        if (!par::master()) continue;
        std::filesystem::create_directories(dir);
        std::ofstream os(dir + "/wallHeatFlux_" + pt.name + ".csv");
        os << "x,y,z,q,Tw\n" << std::setprecision(10);
        for (const auto& r : all) os << r.v[0] << ',' << r.v[1] << ',' << r.v[2] << ',' << r.v[3] << ',' << r.v[4] << '\n';
    }
}

void EnergyModel::writeRestart(const std::string& dir) {
    restart::writeCellArray(dir + "/" + T_.name() + ".bin", m_, T_.internal().data(), 1);
    if (T_.nOldTimes() >= 1) restart::writeCellArray(dir + "/" + T_.name() + "_0.bin", m_, T_.oldRef().data(), 1);
    if (T_.nOldTimes() >= 2) restart::writeCellArray(dir + "/" + T_.name() + "_00.bin", m_, T_.oldOldRef().data(), 1);
}

void EnergyModel::readRestart(const std::string& dir) {
    if (!restart::readCellArray(dir + "/" + T_.name() + ".bin", m_, T_.internal().data(), 1)) {
        if (par::master()) std::cout << "energy: no " << T_.name() << " in " << dir << ", keeping initial value\n";
        return;
    }
    int n = 0;
    T_.oldRef().resize(m_.nTotalCells());
    T_.oldOldRef().resize(m_.nTotalCells());
    if (restart::readCellArray(dir + "/" + T_.name() + "_0.bin", m_, T_.oldRef().data(), 1)) n = 1;
    if (n == 1 && restart::readCellArray(dir + "/" + T_.name() + "_00.bin", m_, T_.oldOldRef().data(), 1)) n = 2;
    T_.setNOldTimes(n);
    T_.correctBoundaryConditions();
    lastTimeIndex_ = flow_.time().timeIndex;
    updateBuoyancy();
}

} // namespace cfd
