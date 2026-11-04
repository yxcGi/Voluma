// LES 亚格子模型（Smagorinsky、WALE、Vreman、σ、动态 Smagorinsky）与壁面模型（WMLES）。
#include "fvm/discretization/Fvc.h"
#include "fvm/mesh/PointSearch.h"
#include "fvm/mesh/WallDistance.h"
#include "fvm/models/TurbulenceModel.h"
#include "fvm/models/WallLaw.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace cfd {

namespace {

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
            law_.Aplus = w.get("Aplus", 17.0);
            law_.nPoints = int(w.get("nPoints", 40.0));
            law_.pressureGradient = w.get("pressureGradient", false);
            if (law_.pressureGradient && law_.type != WallLaw::Type::ODE)
                throw std::runtime_error("wallModel.pressureGradient is only available with type \"ode\"");
            samplingHeight_ = w.get("samplingHeight", 0.0);
            wallModel_ = true;
            if (samplingHeight_ > 0) setupSampling();
            if (par::master())
                std::cout << "wall model: " << WallLaw::name(law_.type)
                          << (law_.pressureGradient ? " (with tangential pressure gradient / body force)" : "")
                          << ", sampling height " << samplingHeight_ << '\n';
        }
        if (vanDriest_) (void)wallDistance(m_);
    }

protected:
    // 由速度梯度与滤波尺度求 ν_t（单元）；动态模型重载 correctNut
    virtual scalar nutCell(const Tensor& gradU, scalar delta) const = 0;
    virtual void correctNut(const std::vector<scalar>& delta) {
        const auto& G = flow_.gradU();
        auto& nut = flow_.nut();
        for (label c = 0; c < m_.nCells(); ++c) nut[c] = nutCell(G[c], delta[c]);
    }

    void correctModel() override {
        const std::size_t nw = wallFaces_.size();
        uTauWall_.assign(nw, 0.0);
        std::vector<scalar> nutw(nw, 0.0);
        std::vector<Vec3> Us;
        if (wallModel_ && samplingHeight_ > 0) Us = sampleVelocity();
        const bool vectorLaw = wallModel_ && law_.type == WallLaw::Type::ODE;
        std::vector<Vec3> gp;
        if (vectorLaw && law_.pressureGradient) gp = fvc::grad(flow_.p());
        Vec3 fBody{};
        if (flow_.forcing.mode == MomentumForcing::Mode::Constant) fBody = flow_.forcing.force;
        if (flow_.forcing.mode == MomentumForcing::Mode::MeanVelocity)
            fBody = flow_.forcing.Ubar * (flow_.forcing.gradP / mag(flow_.forcing.Ubar));
        auto& src = flow_.turbSource();
        if (vectorLaw) src.assign(m_.nTotalCells(), Vec3{});
        const auto& U = flow_.U();
        for (std::size_t i = 0; i < nw; ++i) {
            const WallFace& w = wallFaces_[i];
            if (!wallModel_) {
                // 壁面解析：u_τ 由壁面单元速度梯度得到
                uTauWall_[i] = std::sqrt(nu_ * w.magUt / w.y);
                continue;
            }
            // 采样点的切向相对速度
            Vec3 ut = U[w.cell] - U.bValue(w.face);
            scalar y = w.y;
            if (samplingHeight_ > 0) {
                ut = Us[i] - U.bValue(w.face);
                y = sampleY_[i];
            }
            ut -= dot(ut, w.n) * w.n;
            const scalar ndc = m_.nonOrthDeltaCoeffs()[w.face];
            if (!vectorLaw) {
                const scalar uTau = law_.uTau(mag(ut), y, nu_);
                uTauWall_[i] = uTau;
                // 使壁面扩散通量 ν_eff·Δ⁻¹·|U_P,t| 等于 τ_w = u_τ²
                const scalar nuEff = w.magUt > 1e-12 ? uTau * uTau / (ndc * w.magUt) : nu_;
                nutw[i] = std::max(nuEff - nu_, 0.0);
                continue;
            }
            // ODE 壁面模型（矢量）：切向驱动力 f_p = −∇p + 体力（壁面单元值，切向投影）
            Vec3 fp{};
            if (law_.pressureGradient) {
                fp = fBody - gp[w.cell];
                fp -= dot(fp, w.n) * w.n;
            }
            const Vec3 tau = law_.odeWallShear(ut, y, nu_, fp);
            uTauWall_[i] = std::sqrt(mag(tau));
            // 隐式部分沿壁面单元切向速度方向，其余显式加到壁面单元的动量源
            Vec3 uP = U[w.cell] - U.bValue(w.face);
            uP -= dot(uP, w.n) * w.n;
            const scalar magUP = mag(uP);
            scalar nuEff = nu_;
            if (magUP > 1e-12) {
                const scalar tauPar = dot(tau, uP) / magUP;
                nuEff = std::max(tauPar / (ndc * magUP), nu_);
            }
            nutw[i] = nuEff - nu_;
            // 壁面对流体的作用力 −τ_w|S_f|，减去隐式已包含的 −ν_eff Δ⁻¹ U_P,t |S_f|（单位质量）
            const Vec3 residual = (nuEff * ndc * uP - tau) * (m_.magSf()[w.face] / m_.V()[w.cell]);
            src[w.cell] += residual;
        }
        setWallNut(nutw);
        if (vectorLaw && m_.twoD())
            for (label c = 0; c < m_.nCells(); ++c) src[c][m_.emptyDir()] = 0.0;

        std::vector<scalar> uTauNearest;
        const WallDistance* wd = nullptr;
        std::vector<scalar> delta = delta_;
        if (vanDriest_) {
            wd = &wallDistance(m_);
            std::vector<scalar> bf(m_.nBoundaryFaces(), 0.0);
            for (std::size_t i = 0; i < nw; ++i) bf[wallFaces_[i].face - m_.nInternalFaces()] = uTauWall_[i];
            uTauNearest = gatherWallFaces(m_, bf);
            if (wd->hasWalls)
                for (label c = 0; c < m_.nCells(); ++c) {
                    const scalar y = wd->y[c];
                    const scalar ut = uTauNearest[wd->nearest[c]];
                    const scalar yPlus = y * ut / nu_;
                    delta[c] = std::min(delta[c], kappa_ / Cdelta_ * y * (1 - std::exp(-yPlus / Aplus_)));
                }
        }
        correctNut(delta);
    }

    // 检验滤波（OpenFOAM simpleFilter）：φ̂_P = Σ_f |S_f| φ_f / Σ_f |S_f|，φ_f 线性插值，边界面取单元值；
    // 有效滤波宽度约为 2Δ。结果含幽灵单元（已交换）。
    template <class T> std::vector<T> filter(const std::vector<T>& c) const {
        const label nC = m_.nCells(), nT = m_.nTotalCells(), nI = m_.nInternalFaces();
        std::vector<T> r(nT, Traits<T>::zero());
        std::vector<scalar> ws(nC, 0.0);
        const auto& w = m_.w();
        const auto& magSf = m_.magSf();
        // 进程间面的 owner 或 neighbour 可能是幽灵单元，只累加到本进程单元
        for (label f = 0; f < nI; ++f) {
            const label o = m_.owner()[f], n = m_.neighbour()[f];
            const T fv = w[f] * c[o] + (1.0 - w[f]) * c[n];
            if (o < nC) {
                r[o] += magSf[f] * fv;
                ws[o] += magSf[f];
            }
            if (n < nC) {
                r[n] += magSf[f] * fv;
                ws[n] += magSf[f];
            }
        }
        for (label f = nI; f < m_.nFaces(); ++f) {
            const label o = m_.owner()[f];
            if (o >= nC) continue;
            r[o] += magSf[f] * c[o];
            ws[o] += magSf[f];
        }
        for (label c2 = 0; c2 < nC; ++c2) r[c2] = r[c2] * (1.0 / ws[c2]);
        m_.halo().exchange(r);
        return r;
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

// Vreman (2004)：ν_t = c √(B_β/(α:α))，α_ij = ∂u_j/∂x_i，β_ij = Δ² α_mi α_mj，
// B_β = β11β22 − β12² + β11β33 − β13² + β22β33 − β23²；c ≈ 2.5 C_s²（C_s = 0.17 → 0.07）。
// 层流剪切、纯拉伸等简单流动中 ν_t = 0（与 WALE 类似，优于 Smagorinsky）。
class Vreman : public LESModel {
public:
    Vreman(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i) : LESModel(f, d, b, i) {
        c_ = coeff("c", 0.07);
    }
    std::string type() const override { return "Vreman"; }

protected:
    scalar nutCell(const Tensor& G, scalar delta) const override {
        const scalar aa = magSqr(G);
        if (aa < 1e-300) return 0.0;
        const Tensor beta = (delta * delta) * dot(transpose(G), G);
        const scalar Bb = beta.xx * beta.yy - beta.xy * beta.xy + beta.xx * beta.zz - beta.xz * beta.xz +
                          beta.yy * beta.zz - beta.yz * beta.yz;
        return c_ * std::sqrt(std::max(Bb, 0.0) / aa);
    }
    scalar c_;
};

// 对称 3×3 矩阵的特征值（降序），三角法
static void eigenSymm(const Tensor& A, scalar& l1, scalar& l2, scalar& l3) {
    const scalar p1 = A.xy * A.xy + A.xz * A.xz + A.yz * A.yz;
    if (p1 < 1e-300) {
        scalar v[3] = {A.xx, A.yy, A.zz};
        std::sort(v, v + 3);
        l1 = v[2], l2 = v[1], l3 = v[0];
        return;
    }
    const scalar q = trace(A) / 3;
    const scalar p2 = (A.xx - q) * (A.xx - q) + (A.yy - q) * (A.yy - q) + (A.zz - q) * (A.zz - q) + 2 * p1;
    const scalar p = std::sqrt(p2 / 6);
    const Tensor Bm = (A - q * identity()) * (1.0 / p);
    const scalar detB = Bm.xx * (Bm.yy * Bm.zz - Bm.yz * Bm.zy) - Bm.xy * (Bm.yx * Bm.zz - Bm.yz * Bm.zx) +
                        Bm.xz * (Bm.yx * Bm.zy - Bm.yy * Bm.zx);
    const scalar r = std::clamp(detB / 2, -1.0, 1.0);
    const scalar phi = std::acos(r) / 3;
    l1 = q + 2 * p * std::cos(phi);
    l3 = q + 2 * p * std::cos(phi + 2 * M_PI / 3);
    l2 = 3 * q - l1 - l3;
}

// σ 模型（Nicoud 等 2011）：ν_t = (C_σ Δ)² σ3(σ1 − σ2)(σ2 − σ3)/σ1²，σ_i 为 ∇u 的奇异值（降序）。
// 二维流动、轴对称、各向同性拉伸中 ν_t = 0，近壁 ν_t ~ y³。
class Sigma : public LESModel {
public:
    Sigma(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i) : LESModel(f, d, b, i) {
        Csigma_ = coeff("Csigma", 1.35);
    }
    std::string type() const override { return "sigma"; }

protected:
    scalar nutCell(const Tensor& G, scalar delta) const override {
        const Tensor GtG = dot(transpose(G), G);
        scalar l1, l2, l3;
        eigenSymm(GtG, l1, l2, l3);
        const scalar s1 = std::sqrt(std::max(l1, 0.0)), s2 = std::sqrt(std::max(l2, 0.0)),
                     s3 = std::sqrt(std::max(l3, 0.0));
        if (s1 < 1e-150) return 0.0;
        const scalar D = s3 * (s1 - s2) * (s2 - s3) / (s1 * s1);
        return Csigma_ * Csigma_ * delta * delta * std::max(D, 0.0);
    }
    scalar Csigma_;
};

// 动态 Smagorinsky（Germano 等 1991，Lilly 1992 最小二乘）：
//   L_ij = \hat{u_i u_j} − û_i û_j，M_ij = 2Δ²(\hat{|S|S_ij} − 4|Ŝ|Ŝ_ij)（检验滤波比 2），
//   C = <L:M>/<M:M>（分子分母再做一次滤波作局部平均，Zang 等 1993），ν_t = C Δ² |S|，
//   C 限制在 [0, CsMax²]（不允许负粘性/反向散射）。层流或分辨良好的区域 C → 0。
class DynamicSmagorinsky : public LESModel {
public:
    DynamicSmagorinsky(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i) : LESModel(f, d, b, i) {
        CsMax_ = coeff("CsMax", 0.3);
        if (vanDriest_ && par::master()) std::cout << "dynamicSmagorinsky: vanDriest damping is ignored\n";
    }
    std::string type() const override { return "dynamicSmagorinsky"; }
    // 最近一次的系数场（诊断）
    const std::vector<scalar>& Cs2() const { return C_; }

protected:
    scalar nutCell(const Tensor&, scalar) const override { return 0.0; }

    void correctNut(const std::vector<scalar>&) override {
        const label nT = m_.nTotalCells(), nC = m_.nCells();
        const auto& G = flow_.gradU();
        const auto& U = flow_.U();
        std::vector<Vec3> u(nT);
        std::vector<Tensor> uu(nT), S(nT), SS(nT);
        std::vector<scalar> magS(nT);
        for (label c = 0; c < nT; ++c) {
            u[c] = U[c];
            uu[c] = outer(U[c], U[c]);
            S[c] = symm(G[c]);
            magS[c] = std::sqrt(2 * magSqr(S[c]));
            SS[c] = magS[c] * S[c];
        }
        const auto uf = filter(u);
        const auto uuf = filter(uu);
        const auto Sf = filter(S);
        const auto SSf = filter(SS);
        std::vector<scalar> LM(nT, 0.0), MM(nT, 0.0);
        for (label c = 0; c < nC; ++c) {
            const Tensor L = dev(uuf[c] - outer(uf[c], uf[c]));
            const scalar magSf = std::sqrt(2 * magSqr(Sf[c]));
            const Tensor M = (2 * delta_[c] * delta_[c]) * (SSf[c] - 4 * magSf * Sf[c]);
            LM[c] = doubleDot(L, M);
            MM[c] = doubleDot(M, M);
        }
        m_.halo().exchange(LM);
        m_.halo().exchange(MM);
        const auto aLM = filter(LM), aMM = filter(MM);
        auto& nut = flow_.nut();
        C_.assign(nT, 0.0);
        const scalar Cmax = CsMax_ * CsMax_;
        for (label c = 0; c < nC; ++c) {
            C_[c] = aMM[c] > 1e-300 ? std::clamp(aLM[c] / aMM[c], 0.0, Cmax) : 0.0;
            nut[c] = C_[c] * delta_[c] * delta_[c] * magS[c];
        }
    }

    scalar CsMax_;
    std::vector<scalar> C_;
};

} // namespace

std::unique_ptr<TurbulenceModel> makeLES(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i,
                                         const std::string& type) {
    if (type == "Smagorinsky") return std::make_unique<Smagorinsky>(flow, d, b, i);
    if (type == "Vreman") return std::make_unique<Vreman>(flow, d, b, i);
    if (type == "sigma" || type == "Sigma") return std::make_unique<Sigma>(flow, d, b, i);
    if (type == "dynamicSmagorinsky") return std::make_unique<DynamicSmagorinsky>(flow, d, b, i);
    return std::make_unique<WALE>(flow, d, b, i);
}

} // namespace cfd
