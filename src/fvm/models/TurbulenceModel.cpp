#include "fvm/models/TurbulenceModel.h"

#include "fvm/app/CaseSetup.h"
#include "fvm/discretization/Fvc.h"
#include "fvm/io/Restart.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace cfd {

TurbulenceModel::TurbulenceModel(IncompressibleFlow& flow, const Json& dict, const Json& boundary, const Json& initial)
    : flow_(flow), m_(flow.mesh()), nu_(flow.nu()), dict_(dict), boundary_(boundary), initial_(initial),
      divScheme_(ConvectionScheme::parse(dict.get("div", "limitedLinear 1"))) {
    const bool steady = flow.controls.steady;
    controls_ = parseSolverControls(dict["solver"], {"PBiCGStab", "DILU", 1e-10, steady ? 0.1 : 0.0, 1000, 0, 1});
    relax_ = dict.get("relax", steady ? 0.7 : 1.0);
    kappa_ = coeff("kappa", 0.41);
    E_ = coeff("E", 9.8);
    Cmu_ = coeff("Cmu", 0.09);
    // 线性律与对数律交点 y+ = ln(E y+)/κ
    yPlusLam_ = 11.0;
    for (int i = 0; i < 20; ++i) yPlusLam_ = std::log(std::max(E_ * yPlusLam_, 1.0)) / kappa_;
    flow.nut().assign(m_.nTotalCells(), 0.0);
    flow.wallNuEff().assign(m_.nBoundaryFaces(), -1.0);
    flow.updateTurbulence = [this](IncompressibleFlow&) { correct(); };
}

scalar TurbulenceModel::coeff(const std::string& name, scalar def) {
    const scalar v = dict_["coeffs"].get(name, def);
    coeffs_[name] = v;
    return v;
}

ScalarField& TurbulenceModel::addField(const std::string& name, const std::string& wallDefault, scalar initDefault) {
    owned_.push_back(std::make_unique<ScalarField>(flow_.meshPtr(), name));
    ScalarField& f = *owned_.back();
    const Json zero = Json::parse(R"({"type":"fixedValue","value":0})");
    const Json zg = Json::parse(R"({"type":"zeroGradient"})");
    std::vector<std::string> defaulted;
    for (const auto& p : m_.patches()) {
        if (p.type == PatchType::Empty || p.type == PatchType::Cyclic) continue;
        const Json& pd = boundary_[p.name];
        if (pd.has(name)) {
            setBoundaryCondition(f, p.name, pd[name]);
        } else if (p.type == PatchType::Wall) {
            setBoundaryCondition(f, p.name, wallDefault == "zero" ? zero : zg);
        } else if (p.type == PatchType::Symmetry) {
            f.set<SymmetryBC<scalar>>(p.name);
        } else {
            setBoundaryCondition(f, p.name, zg);
            defaulted.push_back(p.name);
        }
    }
    if (!defaulted.empty() && par::master()) {
        std::cout << "turbulence: " << name << " has no boundary condition on";
        for (auto& d : defaulted) std::cout << ' ' << d;
        std::cout << " -> zeroGradient\n";
    }
    f.setUniform(initial_.get(name, initDefault));
    f.updateBCs(flow_.time().time);
    f.correctBoundaryConditions();
    fields_.push_back(&f);
    return f;
}

void TurbulenceModel::bound(ScalarField& f, scalar minValue) {
    for (auto& v : f.internal()) v = std::max(v, minValue);
    f.correctBoundaryConditions();
}

void TurbulenceModel::updateWallFaces() {
    wallFaces_.clear();
    const auto& U = flow_.U();
    for (label p = 0; p < label(m_.patches().size()); ++p) {
        const Patch& pt = m_.patches()[p];
        if (pt.type != PatchType::Wall) continue;
        for (label f = pt.start; f < pt.end(); ++f) {
            WallFace w;
            w.face = f;
            w.cell = m_.owner()[f];
            w.n = m_.Sf()[f] / m_.magSf()[f];
            w.y = std::max(dot(m_.Cf()[f] - m_.C()[w.cell], w.n), 1e-3 / m_.nonOrthDeltaCoeffs()[f]);
            Vec3 du = U[w.cell] - U.bValue(f);
            du -= dot(du, w.n) * w.n;
            w.magUt = mag(du);
            wallFaces_.push_back(w);
        }
    }
}

void TurbulenceModel::setWallNut(const std::vector<scalar>& nutw) {
    auto& wn = flow_.wallNuEff();
    const label nI = m_.nInternalFaces();
    for (std::size_t i = 0; i < wallFaces_.size(); ++i) wn[wallFaces_[i].face - nI] = nu_ + nutw[i];
}

void TurbulenceModel::solveTransport(ScalarField& psi, const std::vector<scalar>& gammaCell,
                                     const std::vector<scalar>& Sp, const std::vector<scalar>& Su,
                                     const std::vector<label>& fixCells, const std::vector<scalar>& fixValues) {
    const bool steady = flow_.controls.steady;
    const DdtScheme ddtS = steady ? DdtScheme::Steady : flow_.ddtScheme;
    FvMatrix<scalar> M = fvm::ddt(psi, ddtS, flow_.time());
    M += fvm::div(flow_.phi(), psi, divScheme_);
    if (steady) {
        // 有界形式：减去 ψ ∇·φ（未收敛时的连续性误差）
        const auto dv = fvc::div(m_, flow_.phi());
        std::vector<scalar> neg(dv.size());
        for (std::size_t c = 0; c < dv.size(); ++c) neg[c] = -dv[c];
        M += fvm::Sp(neg, psi);
    }
    M -= fvm::laplacian(fvc::interpolate(m_, gammaCell), psi, flow_.laplacianNonOrthCorr);
    M += fvm::Sp(Sp, psi);
    M -= Su;
    if (steady) M.relax(relax_);
    if (!fixCells.empty()) {
        // 在指定单元上固定 ψ（壁面单元的 ε、ω）：去掉该行耦合，把列耦合移到相邻行的源项
        std::vector<char> fixed(m_.nTotalCells(), 0);
        std::vector<scalar> val(m_.nTotalCells(), 0.0);
        for (std::size_t k = 0; k < fixCells.size(); ++k) {
            fixed[fixCells[k]] = 1;
            val[fixCells[k]] = fixValues[k];
            psi[fixCells[k]] = fixValues[k];
        }
        auto& up = M.upper();
        auto& lo = M.lower();
        auto& src = M.source();
        for (label f = 0; f < m_.nInternalFaces(); ++f) {
            const label o = m_.owner()[f], n = m_.neighbour()[f];
            if (fixed[o]) {
                src[n] -= lo[f] * val[o];
                lo[f] = 0;
                up[f] = 0;
            }
            if (fixed[n]) {
                src[o] -= up[f] * val[n];
                up[f] = 0;
                lo[f] = 0;
            }
        }
        const label nI = m_.nInternalFaces();
        for (label bf = 0; bf < m_.nBoundaryFaces(); ++bf)
            if (fixed[m_.owner()[nI + bf]]) {
                M.internalCoeffs()[bf] = 0;
                M.boundarySource()[bf] = 0;
            }
        for (label c : fixCells) {
            if (M.diag()[c] == 0) M.diag()[c] = m_.V()[c];
            src[c] = M.diag()[c] * val[c];
        }
    }
    M.solve(controls_, flow_.controls.verbose);
}

void TurbulenceModel::correct() {
    const TimeState& ts = flow_.time();
    const bool newStep = ts.timeIndex != lastTimeIndex_;
    if (newStep && !flow_.controls.steady)
        for (auto* f : fields_) {
            // 只保存输运场的旧时间层（nut 等派生场不需要）
            f->storeOld();
        }
    lastTimeIndex_ = ts.timeIndex;
    for (auto* f : fields_) {
        f->updateBCs(ts.time);
        f->correctBoundaryConditions();
    }
    updateWallFaces();
    correctModel();
    m_.halo().exchange(flow_.nut());
}

std::vector<scalar> TurbulenceModel::yPlusWall() const {
    std::vector<scalar> yp(wallFaces_.size());
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
        yp[i] = wallFaces_[i].y * (i < uTauWall_.size() ? uTauWall_[i] : 0.0) / nu_;
    return yp;
}

void TurbulenceModel::printYPlus(std::ostream& os) const {
    const auto yp = yPlusWall();
    scalar mn = GREAT, mx = 0, sum = 0, n = double(yp.size());
    for (scalar v : yp) {
        mn = std::min(mn, v);
        mx = std::max(mx, v);
        sum += v;
    }
    mn = par::allMin(mn);
    mx = par::allMax(mx);
    double s[2] = {sum, n};
    par::allSumInPlace(s, 2);
    if (par::master() && s[1] > 0)
        os << "  wall y+: min " << mn << "  mean " << s[0] / s[1] << "  max " << mx << '\n';
}

void TurbulenceModel::writeRestart(const std::string& dir) {
    for (auto* f : fields_) {
        restart::writeCellArray(dir + "/" + f->name() + ".bin", m_, f->internal().data(), 1);
        if (f->nOldTimes() >= 1) restart::writeCellArray(dir + "/" + f->name() + "_0.bin", m_, f->oldRef().data(), 1);
        if (f->nOldTimes() >= 2) restart::writeCellArray(dir + "/" + f->name() + "_00.bin", m_, f->oldOldRef().data(), 1);
    }
}

void TurbulenceModel::readRestart(const std::string& dir) {
    for (auto* f : fields_) {
        if (!restart::readCellArray(dir + "/" + f->name() + ".bin", m_, f->internal().data(), 1)) {
            if (par::master()) std::cout << "turbulence: no " << f->name() << " in " << dir << ", keeping initial value\n";
            continue;
        }
        int n = 0;
        f->oldRef().resize(m_.nTotalCells());
        f->oldOldRef().resize(m_.nTotalCells());
        if (restart::readCellArray(dir + "/" + f->name() + "_0.bin", m_, f->oldRef().data(), 1)) n = 1;
        if (n == 1 && restart::readCellArray(dir + "/" + f->name() + "_00.bin", m_, f->oldOldRef().data(), 1)) n = 2;
        f->setNOldTimes(n);
        f->correctBoundaryConditions();
    }
    lastTimeIndex_ = flow_.time().timeIndex;
}

std::unique_ptr<TurbulenceModel> TurbulenceModel::create(IncompressibleFlow& flow, const Json& d, const Json& boundary,
                                                         const Json& initial) {
    const std::string type = d.get("model", "laminar");
    if (type == "laminar") return nullptr;
    if (type == "Smagorinsky" || type == "WALE") return makeLES(flow, d, boundary, initial, type);
    if (type == "kEpsilon" || type == "realizableKE" || type == "kOmega" || type == "kOmegaSST" ||
        type == "SpalartAllmaras")
        return makeRAS(flow, d, boundary, initial, type);
    throw std::runtime_error("turbulence.model: unknown model '" + type +
                             "' (laminar, Smagorinsky, WALE, kEpsilon, realizableKE, kOmega, kOmegaSST, SpalartAllmaras)");
}

} // namespace cfd
