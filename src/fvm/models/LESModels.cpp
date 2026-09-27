// LES 亚格子模型（Smagorinsky、WALE）与壁面模型（WMLES）。
#include "fvm/mesh/PointSearch.h"
#include "fvm/mesh/WallDistance.h"
#include "fvm/models/TurbulenceModel.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace cfd {

namespace {

// ------------------------------------------------------------ 壁面律：由 (U, y) 求 u_τ
class WallLaw {
public:
    enum class Type { LogLaw, Spalding, Musker, PowerLaw };
    Type type = Type::Spalding;
    scalar kappa = 0.41, B = 5.2;
    scalar A = 8.3, Bpow = 1.0 / 7.0;  // Werner-Wengle 幂律 u+ = A y+^B

    static Type parse(const std::string& s) {
        if (s == "logLaw") return Type::LogLaw;
        if (s == "spalding") return Type::Spalding;
        if (s == "musker") return Type::Musker;
        if (s == "powerLaw" || s == "wernerWengle") return Type::PowerLaw;
        throw std::runtime_error("wallModel.type must be logLaw, spalding, musker or powerLaw");
    }

    // 给定 y+ 的 u+（隐式律用于牛顿迭代）
    scalar uPlusMusker(scalar yp) const {
        // Musker (1979)：κ = 0.41，s = 0.001093
        return 5.424 * std::atan((2 * yp - 8.15) / 16.7) +
               std::log10(std::pow(yp + 10.6, 9.6) / std::pow(yp * yp - 8.15 * yp + 86, 2)) - 3.52;
    }
    // Spalding：y+ = u+ + e^{−κB}(e^{κu+} − 1 − κu+ − (κu+)²/2 − (κu+)³/6)
    scalar yPlusSpalding(scalar up) const {
        const scalar k = kappa * up;
        return up + std::exp(-kappa * B) * (std::exp(k) - 1 - k - k * k / 2 - k * k * k / 6);
    }

    scalar uTau(scalar U, scalar y, scalar nu) const {
        if (U <= 0 || y <= 0) return 0;
        const scalar uLam = std::sqrt(nu * U / y);  // 线性律
        switch (type) {
        case Type::PowerLaw: {
            // Werner & Wengle (1991) 在首层内积分的显式形式
            const scalar lim = nu / (2 * y) * std::pow(A, 2 / (1 - Bpow));
            if (U <= lim) return std::sqrt(2 * nu * U / y);
            const scalar t = (1 - Bpow) / 2 * std::pow(A, (1 + Bpow) / (1 - Bpow)) * std::pow(nu / y, 1 + Bpow) +
                             (1 + Bpow) / A * std::pow(nu / y, Bpow) * U;
            return std::sqrt(std::pow(t, 2 / (1 + Bpow)));
        }
        case Type::LogLaw: {
            // U/u_τ = ln(E y u_τ/ν)/κ，y+ < 11 时取线性律
            const scalar E = std::exp(kappa * B);
            scalar ut = std::max(uLam, 1e-12);
            for (int it = 0; it < 50; ++it) {
                const scalar yp = y * ut / nu;
                if (yp < 11.0) return uLam;
                const scalar f = ut * std::log(E * yp) / kappa - U;
                const scalar df = (std::log(E * yp) + 1) / kappa;
                const scalar dut = f / df;
                ut = std::max(ut - dut, 0.5 * ut);
                if (std::abs(dut) < 1e-10 * ut) break;
            }
            return ut;
        }
        case Type::Spalding:
        case Type::Musker: {
            // 牛顿迭代（u+ 与 y+ 之间单调），失败时用二分
            auto res = [&](scalar ut) {
                const scalar up = U / ut, yp = y * ut / nu;
                // 两者都在 u_τ 很小时为正、随 u_τ 单调减小
                return type == Type::Spalding ? yPlusSpalding(up) - yp : up - uPlusMusker(yp);
            };
            scalar lo = 1e-8 * uLam + 1e-14, hi = uLam * 1.0001 + 1e-14;
            // res(lo) > 0；线性律值处 res ≥ 0，向上加倍直到变号
            while (res(hi) > 0) hi *= 2;
            for (int it = 0; it < 100; ++it) {
                const scalar mid = 0.5 * (lo + hi);
                (res(mid) > 0 ? lo : hi) = mid;
                if (hi - lo < 1e-12 * hi) break;
            }
            return 0.5 * (lo + hi);
        }
        }
        return uLam;
    }
};

// ------------------------------------------------------------ LES 基类
class LESModel : public TurbulenceModel {
public:
    LESModel(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i) : TurbulenceModel(flow, d, b, i) {
        // 几何滤波尺度：cubeRootVol（二维取 √(面积)）
        delta_.assign(m_.nTotalCells(), 0.0);
        scalar depth = 1;
        if (m_.twoD()) {
            const int ed = m_.emptyDir();
            scalar lo = GREAT, hi = -GREAT;
            for (const auto& p : m_.vtk().points) {
                lo = std::min(lo, p[ed]);
                hi = std::max(hi, p[ed]);
            }
            depth = par::allMax(hi) - par::allMin(lo);
        }
        const scalar deltaCoeff = coeff("deltaCoeff", 1.0);
        for (label c = 0; c < m_.nTotalCells(); ++c)
            delta_[c] = deltaCoeff * (m_.twoD() ? std::sqrt(m_.V()[c] / depth) : std::cbrt(m_.V()[c]));
        vanDriest_ = d.get("vanDriest", false);
        Aplus_ = coeff("Aplus", 26.0);
        Cdelta_ = coeff("Cdelta", 0.158);
        if (d.has("wallModel")) {
            const Json& w = d["wallModel"];
            law_.type = WallLaw::parse(w.get("type", "spalding"));
            law_.kappa = w.get("kappa", 0.41);
            law_.B = w.get("B", 5.2);
            law_.A = w.get("A", 8.3);
            law_.Bpow = w.get("Bpow", 1.0 / 7.0);
            samplingHeight_ = w.get("samplingHeight", 0.0);
            wallModel_ = true;
            if (samplingHeight_ > 0) setupSampling();
        }
        if (vanDriest_) (void)wallDistance(m_);
    }

protected:
    // 由速度梯度与滤波尺度求 ν_t（单元）
    virtual scalar nutCell(const Tensor& gradU, scalar delta) const = 0;

    void correctModel() override {
        const auto& G = flow_.gradU();
        const std::size_t nw = wallFaces_.size();
        uTauWall_.assign(nw, 0.0);
        std::vector<scalar> nutw(nw, 0.0);
        std::vector<Vec3> Us;
        if (wallModel_ && samplingHeight_ > 0) Us = sampleVelocity();
        for (std::size_t i = 0; i < nw; ++i) {
            const WallFace& w = wallFaces_[i];
            if (!wallModel_) {
                // 壁面解析：u_τ 由壁面单元速度梯度得到
                uTauWall_[i] = std::sqrt(nu_ * w.magUt / w.y);
                continue;
            }
            scalar Ut = w.magUt, y = w.y;
            if (samplingHeight_ > 0) {
                Vec3 du = Us[i] - flow_.U().bValue(w.face);
                du -= dot(du, w.n) * w.n;
                Ut = mag(du);
                y = sampleY_[i];
            }
            const scalar ut = law_.uTau(Ut, y, nu_);
            uTauWall_[i] = ut;
            // 使壁面扩散通量 ν_eff·Δ⁻¹·|U_P,t| 等于 τ_w = u_τ²
            const scalar ndc = m_.nonOrthDeltaCoeffs()[w.face];
            const scalar nuEff = w.magUt > 1e-12 ? ut * ut / (ndc * w.magUt) : nu_;
            nutw[i] = std::max(nuEff - nu_, 0.0);
        }
        setWallNut(nutw);

        std::vector<scalar> uTauNearest;
        const WallDistance* wd = nullptr;
        if (vanDriest_) {
            wd = &wallDistance(m_);
            std::vector<scalar> bf(m_.nBoundaryFaces(), 0.0);
            for (std::size_t i = 0; i < nw; ++i) bf[wallFaces_[i].face - m_.nInternalFaces()] = uTauWall_[i];
            uTauNearest = gatherWallFaces(m_, bf);
        }
        auto& nut = flow_.nut();
        for (label c = 0; c < m_.nCells(); ++c) {
            scalar delta = delta_[c];
            if (wd && wd->hasWalls) {
                const scalar y = wd->y[c];
                const scalar ut = uTauNearest[wd->nearest[c]];
                const scalar yPlus = y * ut / nu_;
                delta = std::min(delta, kappa_ / Cdelta_ * y * (1 - std::exp(-yPlus / Aplus_)));
            }
            nut[c] = nutCell(G[c], delta);
        }
    }

    void setupSampling() {
        // 每个壁面面的采样点 x = C_f − h n（n 为外法向），所有进程持有全局列表
        const label nI = m_.nInternalFaces();
        std::vector<scalar> px(m_.nBoundaryFaces(), 0.0), py = px, pz = px;
        for (const auto& p : m_.patches()) {
            if (p.type != PatchType::Wall) continue;
            for (label f = p.start; f < p.end(); ++f) {
                const Vec3 x = m_.Cf()[f] - samplingHeight_ * m_.Sf()[f] / m_.magSf()[f];
                px[f - nI] = x.x;
                py[f - nI] = x.y;
                pz[f - nI] = x.z;
            }
        }
        const auto gx = gatherWallFaces(m_, px), gy = gatherWallFaces(m_, py), gz = gatherWallFaces(m_, pz);
        const std::size_t n = gx.size();
        std::vector<Vec3> own(m_.C().begin(), m_.C().begin() + m_.nCells());
        const PointTree tree(own);
        std::vector<double> d2(n), id(n);
        std::vector<long> loc(n);
        for (std::size_t k = 0; k < n; ++k) {
            loc[k] = tree.nearest({gx[k], gy[k], gz[k]}, d2[k]);
            if (loc[k] < 0) d2[k] = GREAT;
        }
        std::vector<double> dmin = d2;
        par::allMinInPlace(dmin.data(), int(n));
        for (std::size_t k = 0; k < n; ++k)
            id[k] = (loc[k] >= 0 && d2[k] == dmin[k]) ? double(m_.cellGlobal()[loc[k]]) : 1e300;
        std::vector<double> idmin = id;
        par::allMinInPlace(idmin.data(), int(n));
        sampleCell_.assign(n, -1);
        std::vector<double> cx(3 * n, 0.0);
        for (std::size_t k = 0; k < n; ++k)
            if (loc[k] >= 0 && id[k] == idmin[k]) {
                sampleCell_[k] = label(loc[k]);
                const Vec3& c = m_.C()[loc[k]];
                cx[3 * k] = c.x, cx[3 * k + 1] = c.y, cx[3 * k + 2] = c.z;
            }
        par::allSumInPlace(cx.data(), int(cx.size()));
        // 本进程壁面面的采样距离（采样单元中心到壁面的法向距离）
        const auto& gidx = wallFaceGlobalIndex(m_);
        updateWallFaces();
        sampleY_.assign(wallFaces_.size(), 0.0);
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const label f = wallFaces_[i].face;
            const glabel g = gidx[f - nI];
            const Vec3 c{cx[3 * g], cx[3 * g + 1], cx[3 * g + 2]};
            sampleY_[i] = std::max(dot(m_.Cf()[f] - c, wallFaces_[i].n), wallFaces_[i].y);
        }
        scalar ymin = GREAT, ymax = 0;
        for (scalar y : sampleY_) ymin = std::min(ymin, y), ymax = std::max(ymax, y);
        ymin = par::allMin(ymin);
        ymax = par::allMax(ymax);
        if (par::master())
            std::cout << "wall model: " << n << " wall faces sampled at distance " << ymin << " .. " << ymax << '\n';
    }

    // 各本进程壁面面对应采样点的速度
    std::vector<Vec3> sampleVelocity() const {
        const std::size_t n = sampleCell_.size();
        std::vector<double> v(3 * n, 0.0);
        const auto& U = flow_.U();
        for (std::size_t k = 0; k < n; ++k)
            if (sampleCell_[k] >= 0) {
                const Vec3& u = U[sampleCell_[k]];
                v[3 * k] = u.x, v[3 * k + 1] = u.y, v[3 * k + 2] = u.z;
            }
        par::allSumInPlace(v.data(), int(v.size()));
        const auto& gidx = wallFaceGlobalIndex(m_);
        std::vector<Vec3> out(wallFaces_.size());
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const glabel g = gidx[wallFaces_[i].face - m_.nInternalFaces()];
            out[i] = {v[3 * g], v[3 * g + 1], v[3 * g + 2]};
        }
        return out;
    }

    std::vector<scalar> delta_;
    bool vanDriest_ = false;
    scalar Aplus_ = 26, Cdelta_ = 0.158;
    bool wallModel_ = false;
    WallLaw law_;
    scalar samplingHeight_ = 0;
    std::vector<label> sampleCell_;  // 全局采样点 → 本进程单元（非本进程为 −1）
    std::vector<scalar> sampleY_;    // 本进程壁面面的采样距离
};

// ν_t = (C_s Δ)² √(2 S:S)
class Smagorinsky : public LESModel {
public:
    Smagorinsky(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i) : LESModel(f, d, b, i) {
        Cs_ = coeff("Cs", 0.17);
    }
    std::string type() const override { return "Smagorinsky"; }

protected:
    scalar nutCell(const Tensor& G, scalar delta) const override {
        const Tensor S = symm(G);
        return Cs_ * Cs_ * delta * delta * std::sqrt(2 * magSqr(S));
    }
    scalar Cs_;
};

// WALE（Nicoud & Ducros 1999）：ν_t = (C_w Δ)² (Sd:Sd)^{3/2} / ((S:S)^{5/2} + (Sd:Sd)^{5/4})
class WALE : public LESModel {
public:
    WALE(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i) : LESModel(f, d, b, i) {
        Cw_ = coeff("Cw", 0.325);
    }
    std::string type() const override { return "WALE"; }

protected:
    scalar nutCell(const Tensor& G, scalar delta) const override {
        const Tensor S = symm(G);
        const Tensor Sd = dev(symm(dot(G, G)));
        const scalar sd2 = magSqr(Sd), s2 = magSqr(S);
        return Cw_ * Cw_ * delta * delta * std::pow(sd2, 1.5) / (std::pow(s2, 2.5) + std::pow(sd2, 1.25) + 1e-300);
    }
    scalar Cw_;
};

} // namespace

std::unique_ptr<TurbulenceModel> makeLES(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i,
                                         const std::string& type) {
    if (type == "Smagorinsky") return std::make_unique<Smagorinsky>(flow, d, b, i);
    return std::make_unique<WALE>(flow, d, b, i);
}

} // namespace cfd
