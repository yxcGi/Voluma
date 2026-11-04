// RANS 模型：kEpsilon、RNGkEpsilon、realizableKE、kOmega、kOmegaSST、BSLEARSM、SpalartAllmaras（公式与系数同 OpenFOAM）。
//
// 壁面处理（wallTreatment）：
//   standard  nutkWallFunction：y+ = C_μ^¼√k y/ν，y+ > y+_lam 时按对数律给 ν_t,w、壁面单元 G 与 ε/ω（切换式）
//   enhanced  全 y+ 处理：u_τ 由 Spalding 统一律从壁面单元速度求出（任意 y+ 连续），壁面单元 G 用同一条律的
//             du+/dy+ 给 P = τ_w(1 − D)·D u_τ²/ν（粘性底层 → 0，对数层 → u_τ³/(κy)），
//             ε 用 Wolfshtein 长度尺度 k^{3/2}/ℓ_ε（两极限分别为 2νk/y² 与 C_μ^{3/4}k^{3/2}/(κy)）；
//             k-ε 系列另加两层模型：Re_y = √k y/ν < Re_y*（200）的单元 ε 取代数值，ν_t 按 λ_ε 与一方程模型混合
//             （Fluent 的 Enhanced Wall Treatment 同源）。k-ω 系列的 ω 壁面值两种处理都是粘性/对数混合值。
#include "fvm/discretization/Fvc.h"
#include "fvm/mesh/WallDistance.h"
#include "fvm/models/TurbulenceModel.h"
#include "fvm/models/WallLaw.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace cfd {

namespace {

constexpr scalar kSmall = 1e-15;

class RASModel : public TurbulenceModel {
public:
    enum class WallTreatment { Standard, Enhanced };

    RASModel(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i, bool wallFunctionsDefault)
        : TurbulenceModel(f, d, b, i) {
        wallFunctions_ = d.get("wallFunctions", wallFunctionsDefault);
        const std::string wt = d.get("wallTreatment", "standard");
        if (wt == "standard") wallTreatment_ = WallTreatment::Standard;
        else if (wt == "enhanced" || wt == "allYPlus" || wt == "blended") wallTreatment_ = WallTreatment::Enhanced;
        else throw std::runtime_error("turbulence.wallTreatment must be standard or enhanced");
        if (wallTreatment_ == WallTreatment::Enhanced) wallFunctions_ = true;
        law_.type = WallLaw::Type::Spalding;
        law_.kappa = kappa_;
        law_.B = std::log(E_) / kappa_;
    }

protected:
    // 壁面单元：每个壁面面的权重 1/（该单元壁面面数）
    std::vector<scalar> wallWeights() const {
        std::vector<int> cnt(m_.nTotalCells(), 0);
        for (auto& w : wallFaces_) ++cnt[w.cell];
        std::vector<scalar> wt(wallFaces_.size());
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) wt[i] = 1.0 / cnt[wallFaces_[i].cell];
        return wt;
    }
    // 壁面切应力：填 uTauWall_、yPlusWall_，并设置壁面 ν_eff
    //   standard：nutkWallFunction（y+ = C_μ^¼ √k y/ν，y+ > y+_lam 时 ν_t,w = ν(κ y+/ln(E y+) − 1)）
    //   enhanced：Spalding 统一律由壁面单元切向速度求 u_τ，ν_t,w = u_τ²/(Δ⁻¹|U_P,t|) − ν
    void wallShear(const ScalarField* k) {
        const scalar Cmu25 = std::pow(Cmu_, 0.25);
        const std::size_t nw = wallFaces_.size();
        std::vector<scalar> nutw(nw, 0.0);
        uTauWall_.assign(nw, 0.0);
        yPlusWall_.assign(nw, 0.0);
        for (std::size_t i = 0; i < nw; ++i) {
            const auto& w = wallFaces_[i];
            if (wallTreatment_ == WallTreatment::Enhanced || !k) {
                // 低 Re（无壁面函数）时报告解析的 u_τ = √(ν|U_P,t|/y)
                const scalar ut = wallFunctions_ ? law_.uTau(w.magUt, w.y, nu_) : std::sqrt(nu_ * w.magUt / w.y);
                uTauWall_[i] = ut;
                yPlusWall_[i] = w.y * ut / nu_;
                if (wallFunctions_ && w.magUt > 1e-12) {
                    const scalar ndc = m_.nonOrthDeltaCoeffs()[w.face];
                    nutw[i] = std::max(ut * ut / (ndc * w.magUt) - nu_, 0.0);
                }
                continue;
            }
            const scalar yPlus = Cmu25 * std::sqrt(std::max((*k)[w.cell], 0.0)) * w.y / nu_;
            yPlusWall_[i] = yPlus;
            if (wallFunctions_ && yPlus > yPlusLam_) nutw[i] = nu_ * (yPlus * kappa_ / std::log(E_ * yPlus) - 1.0);
            // 输出用 u_τ：由壁面切应力
            const scalar tau = (nu_ + nutw[i]) * w.magUt / w.y;
            uTauWall_[i] = std::sqrt(std::max(tau, 0.0));
        }
        setWallNut(nutw);
    }
    // 壁面单元的产生项（单位体积，按壁面面加权累加到 GW；useLog 标记该单元是否被替换）
    //   standard：y+ > y+_lam 时 G = ν_eff,w (|U_P,t|/y) C_μ^¼√k/(κ y)（对数律），否则保留单元值
    //   enhanced：G = τ_w(1 − D)·(τ_w/ν)D = u_τ⁴/ν · D(1 − D)，D = du+/dy+ 取自 Spalding 律（与壁面切应力同一条律；
    //             粘性底层 D → 1、P → 0，对数层 D → 1/(κy+)、P → u_τ³/(κy)，y+ ≈ 11 处 P+ ≈ 0.25 与 DNS 一致）
    void wallProduction(std::vector<scalar>& GW, std::vector<char>& useLog, const ScalarField& k) const {
        const scalar Cmu25 = std::pow(Cmu_, 0.25);
        const auto wt = wallWeights();
        const auto& wn = flow_.wallNuEff();
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const auto& w = wallFaces_[i];
            if (wallTreatment_ == WallTreatment::Enhanced) {
                const scalar ut = uTauWall_[i];
                const scalar D = std::clamp(law_.duPlusDyPlus(yPlusWall_[i]), 0.0, 1.0);
                GW[w.cell] += wt[i] * std::pow(ut, 4) / nu_ * D * (1 - D);
                useLog[w.cell] = 1;
            } else if (wallFunctions_ && yPlusWall_[i] > yPlusLam_) {
                const scalar kc = std::max(k[w.cell], 0.0);
                const scalar nuw = wn[w.face - m_.nInternalFaces()];
                GW[w.cell] += wt[i] * nuw * (w.magUt / w.y) * Cmu25 * std::sqrt(kc) / (kappa_ * w.y);
                useLog[w.cell] = 1;
            }
        }
    }
    std::vector<scalar> GbyNu() const {
        const auto& G = flow_.gradU();
        std::vector<scalar> g(m_.nTotalCells(), 0.0);
        for (label c = 0; c < m_.nCells(); ++c) g[c] = doubleDot(dev(2.0 * symm(G[c])), G[c]);
        return g;
    }
    // 显式动量源 S = −∇·dev(R) − ∇·(ν_t ∇U) − ∇·(ν_t dev2((∇U)ᵀ))：动量方程以 ν_eff = ν + ν_t 隐式求解
    // （UEqn −= laplacian(ν_eff,U) 与显式转置项），这里加上模型应力 dev(R) 的散度并减去已含的线性涡粘部分，
    // 使收敛后的总湍流应力为 dev(R)（"表观粘度" 稳定化，Lien & Leschziner 1994 / OpenFOAM ReynoldsStress）。
    // 涡粘抵消项的面梯度：内部面用单元梯度的面插值（它对 2Δ 棋盘模态不敏感，于是隐式紧致拉普拉斯对棋盘模态的
    // 阻尼得以保留——否则 R 由中心梯度生成、再经中心插值回到动量方程，U 与 R 的奇偶解耦模态无阻尼，剖面出现锯齿，
    // 收敛态还依赖松弛路径）；壁面单元相邻的面例外，改用与隐式项完全相同的紧致梯度（壁面单元的 Gauss 梯度因壁面
    // 一侧差分比真实梯度大数倍，若用其插值会给面 1–2 留下很大的虚假应力，使近壁 −uv 偏高、Re_τ 偏低；
    // 该单元的 R 已由壁面处理给定，不存在解耦问题）；非壁面边界面用边界条件的梯度系数。
    // 壁面上雷诺应力与涡粘项都不算（壁面切应力由壁面函数的 ν_eff,w 隐式给出）。
    void setStressSource(const std::vector<Tensor>& devR) {
        const label nT = m_.nTotalCells(), nI = m_.nInternalFaces();
        const auto tf = fvc::interpolate(m_, devR);
        const auto& nut = flow_.nut();
        const auto nutf = fvc::interpolate(m_, nut);
        const auto& U = flow_.U();
        const auto& gU = flow_.gradU();
        const auto gUf = fvc::interpolate(m_, gU);
        const bool corr = flow_.laplacianNonOrthCorr;
        auto& src = flow_.turbSource();
        src.assign(nT, Vec3{0, 0, 0});
        std::vector<char> wallCell(nT, 0);
        for (const auto& w : wallFaces_) wallCell[w.cell] = 1;
        // 面通量 F_f（朝 owner 外法向）：−dev(R)·S_f − ν_t[∇U·S_f] − ν_t dev2(G_fᵀ)·S_f
        auto transposed = [&](label f, scalar nt) {
            const Tensor& G = gUf[f];
            return nt * (dot(G, m_.Sf()[f]) - (2.0 / 3.0) * trace(G) * m_.Sf()[f]);
        };
        for (label f = 0; f < nI; ++f) {
            const label o = m_.owner()[f], n = m_.neighbour()[f];
            const scalar nt = nutf[f], gs = nt * m_.magSf()[f];
            Vec3 F = -1.0 * dot(tf[f], m_.Sf()[f]);
            if (wallCell[o] || wallCell[n]) {
                F -= gs * (m_.nonOrthDeltaCoeffs()[f] * (U[n] - U[o]));
                if (corr) F -= gs * gradDot(m_.corrVecs()[f], gUf[f]);
            } else {
                F -= nt * dot(m_.Sf()[f], gUf[f]);
            }
            F -= transposed(f, nt);
            src[o] += F;
            src[n] -= F;
        }
        for (label p = 0; p < label(m_.patches().size()); ++p) {
            const Patch& pt = m_.patches()[p];
            if (pt.type == PatchType::Wall || pt.type == PatchType::Empty) continue;
            const auto& bc = U.bc(p);
            for (label i = 0; i < pt.size; ++i) {
                const label f = pt.start + i, o = m_.owner()[f];
                const scalar nt = nut[o], gs = nt * m_.magSf()[f];
                Vec3 F = -1.0 * dot(tf[f], m_.Sf()[f]);
                F -= gs * (cmptMultiply(bc.gradientInternalCoeff(U, i), U[o]) + bc.gradientBoundaryCoeff(U, i));
                F -= transposed(f, nt);
                src[o] += F;
            }
        }
        for (label c = 0; c < m_.nCells(); ++c) src[c] = src[c] / m_.V()[c];
        if (m_.twoD())
            for (label c = 0; c < m_.nCells(); ++c) src[c][m_.emptyDir()] = 0.0;
    }
    std::vector<scalar> diffusivity(const std::vector<scalar>& nutCoeff) const {
        const auto& nut = flow_.nut();
        std::vector<scalar> D(m_.nTotalCells());
        for (label c = 0; c < m_.nTotalCells(); ++c) D[c] = nutCoeff[c] * nut[c] + nu_;
        return D;
    }
    std::vector<scalar> diffusivity(scalar sigmaInv) const {
        return diffusivity(std::vector<scalar>(m_.nTotalCells(), sigmaInv));
    }

    bool wallFunctions_ = true;
    WallTreatment wallTreatment_ = WallTreatment::Standard;
    WallLaw law_;                     // enhanced 处理的壁面律（Spalding，κ、E 同模型）
    std::vector<scalar> yPlusWall_;   // 按 wallFaces_，最近一次 wallShear
};

// ------------------------------------------------------------------ k-ε（标准 / RNG / realizable）
// RNG（Yakhot 等 1992，系数同 OpenFOAM RNGkEpsilon）：C_μ 0.0845，C1 1.42，C2 1.68，σ_k = σ_ε = 0.71942，
//   ε 方程多出 R = η(1 − η/η0)/(1 + βη³) 项：(C1 − R) G ε/k，η = √(2S:S) k/ε，η0 4.38，β 0.012
class KEpsilon : public RASModel {
public:
    enum class Variant { Standard, RNG, Realizable };
    KEpsilon(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i, Variant variant)
        : RASModel(f, d, b, i, true), variant_(variant) {
        if (!wallFunctions_) throw std::runtime_error("kEpsilon models need wall functions (wallFunctions: true)");
        realizable_ = variant == Variant::Realizable;
        if (variant == Variant::Realizable) {
            A0_ = coeff("A0", 4.0);
            C2_ = coeff("C2", 1.9);
            sigmak_ = coeff("sigmak", 1.0);
            sigmaEps_ = coeff("sigmaEps", 1.2);
            CmuNut_ = Cmu_;
        } else if (variant == Variant::RNG) {
            CmuNut_ = coeff("CmuRNG", 0.0845);
            C1_ = coeff("C1", 1.42);
            C2_ = coeff("C2", 1.68);
            sigmak_ = coeff("sigmak", 0.71942);
            sigmaEps_ = coeff("sigmaEps", 0.71942);
            eta0_ = coeff("eta0", 4.38);
            betaRNG_ = coeff("beta", 0.012);
        } else {
            CmuNut_ = Cmu_;
            C1_ = coeff("C1", 1.44);
            C2_ = coeff("C2", 1.92);
            sigmak_ = coeff("sigmak", 1.0);
            sigmaEps_ = coeff("sigmaEps", 1.3);
        }
        // 两层模型常数（enhanced）：C_ℓ = κ C_μ^{-3/4}，A_μ 70，A_ε 2C_ℓ，Re_y* 200，过渡宽度 ΔRe_y = 0.1 Re_y*
        Cl_ = kappa_ * std::pow(Cmu_, -0.75);
        Amu_ = coeff("Amu", 70.0);
        Aeps_ = coeff("Aeps", 2 * Cl_);
        ReyStar_ = coeff("ReyStar", 200.0);
        ReyWidth_ = coeff("ReyWidth", 0.1 * ReyStar_);
        if (wallTreatment_ == WallTreatment::Enhanced) y_ = &wallDistance(m_).y;
        k_ = &addField("k", "zeroGradient", 1e-4);
        eps_ = &addField("epsilon", "zeroGradient", 1e-5);
        bound(*k_, kMin_);
        bound(*eps_, epsMin_);
        correctNut(nullptr);
    }
    std::string type() const override {
        return variant_ == Variant::Realizable ? "realizableKE" : variant_ == Variant::RNG ? "RNGkEpsilon" : "kEpsilon";
    }

protected:
    // 两层模型的长度尺度与混合函数
    scalar lEps(scalar Rey, scalar y) const { return Cl_ * y * (1 - std::exp(-Rey / Aeps_)); }
    scalar lMu(scalar Rey, scalar y) const { return Cl_ * y * (1 - std::exp(-Rey / Amu_)); }
    scalar lambdaEps(scalar Rey) const {
        const scalar A = ReyWidth_ / std::atanh(0.98);
        return 0.5 * (1 + std::tanh((Rey - ReyStar_) / A));
    }

    void correctNut(const std::vector<scalar>* rCmu) {
        auto& nut = flow_.nut();
        const auto& k = *k_;
        const auto& e = *eps_;
        for (label c = 0; c < m_.nTotalCells(); ++c)
            nut[c] = (rCmu ? (*rCmu)[c] : CmuNut_) * k[c] * k[c] / std::max(e[c], epsMin_);
        if (wallTreatment_ == WallTreatment::Enhanced) {
            // 两层：ν_t = λ ν_t,k-ε + (1 − λ) C_μ ℓ_μ √k
            const auto& y = *y_;
            for (label c = 0; c < m_.nTotalCells(); ++c) {
                const scalar kc = std::max(k[c], 0.0);
                const scalar Rey = std::sqrt(kc) * y[c] / nu_;
                const scalar lam = lambdaEps(Rey);
                nut[c] = lam * nut[c] + (1 - lam) * Cmu_ * lMu(Rey, y[c]) * std::sqrt(kc);
            }
        }
    }

    // realizable 的 C_μ
    std::vector<scalar> rCmu(const std::vector<scalar>& magS) const {
        const auto& G = flow_.gradU();
        std::vector<scalar> r(m_.nTotalCells(), Cmu_);
        for (label c = 0; c < m_.nCells(); ++c) {
            const Tensor S = dev(symm(G[c]));
            const scalar SS = doubleDot(S, S);
            const scalar W = 2 * std::sqrt(2.0) * doubleDot(dot(S, S), S) / (magS[c] * SS + 1e-300);
            const scalar phis = std::acos(std::clamp(std::sqrt(6.0) * W, -1.0, 1.0)) / 3.0;
            const scalar As = std::sqrt(6.0) * std::cos(phis);
            const scalar Us = std::sqrt(magS[c] * magS[c] / 2.0 + magSqr(skew(G[c])));
            r[c] = 1.0 / (A0_ + As * Us * (*k_)[c] / std::max((*eps_)[c], epsMin_));
        }
        m_.halo().exchange(r);
        return r;
    }

    void correctModel() override {
        auto& k = *k_;
        auto& e = *eps_;
        const auto& nut = flow_.nut();
        const auto gbn = GbyNu();
        const label nC = m_.nCells();
        std::vector<scalar> G(m_.nTotalCells(), 0.0), magS(m_.nTotalCells(), 0.0);
        for (label c = 0; c < nC; ++c) G[c] = nut[c] * gbn[c];
        if (realizable_ || variant_ == Variant::RNG) {
            const auto& gU = flow_.gradU();
            for (label c = 0; c < nC; ++c) magS[c] = std::sqrt(2 * magSqr(dev(symm(gU[c]))));
        }

        // 壁面函数：壁面单元的 ε 与 G
        wallShear(&k);
        const scalar Cmu75 = std::pow(Cmu_, 0.75);
        const auto wt = wallWeights();
        std::vector<scalar> epsW(m_.nTotalCells(), 0.0), GW(m_.nTotalCells(), 0.0);
        std::vector<char> isWall(m_.nTotalCells(), 0), useLog(m_.nTotalCells(), 0);
        const bool enhanced = wallTreatment_ == WallTreatment::Enhanced;
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const auto& w = wallFaces_[i];
            const scalar kc = std::max(k[w.cell], 0.0);
            isWall[w.cell] = 1;
            if (enhanced) {
                // Wolfshtein 长度尺度：两极限为 2νk/y²（Re_y → 0）与 C_μ^{3/4}k^{3/2}/(κ y)（对数层）
                const scalar Rey = std::sqrt(kc) * w.y / nu_;
                epsW[w.cell] += wt[i] * std::pow(kc, 1.5) / std::max(lEps(Rey, w.y), 1e-300);
            } else if (yPlusWall_[i] < yPlusLam_) {
                epsW[w.cell] += wt[i] * 2.0 * kc * nu_ / (w.y * w.y);
            } else {
                epsW[w.cell] += wt[i] * Cmu75 * std::pow(kc, 1.5) / (kappa_ * w.y);
            }
        }
        wallProduction(GW, useLog, k);
        std::vector<label> fixCells;
        std::vector<scalar> fixVals;
        for (label c = 0; c < nC; ++c)
            if (isWall[c]) {
                if (useLog[c]) G[c] = GW[c];
                fixCells.push_back(c);
                fixVals.push_back(std::max(epsW[c], epsMin_));
            }
        if (enhanced) {
            // 两层模型：粘性影响区（Re_y < Re_y*）的非壁面单元 ε 取代数值 k^{3/2}/ℓ_ε
            const auto& y = *y_;
            for (label c = 0; c < nC; ++c) {
                if (isWall[c]) continue;
                const scalar kc = std::max(k[c], 0.0);
                const scalar Rey = std::sqrt(kc) * y[c] / nu_;
                if (Rey < ReyStar_) {
                    fixCells.push_back(c);
                    fixVals.push_back(std::max(std::pow(kc, 1.5) / std::max(lEps(Rey, y[c]), 1e-300), epsMin_));
                }
            }
        }

        // ε 方程
        std::vector<scalar> Sp(m_.nTotalCells(), 0.0), Su(m_.nTotalCells(), 0.0);
        for (label c = 0; c < nC; ++c) {
            const scalar kc = std::max(k[c], kMin_), ec = std::max(e[c], epsMin_);
            if (realizable_) {
                const scalar eta = magS[c] * kc / ec;
                const scalar C1 = std::max(eta / (5 + eta), 0.43);
                Su[c] = C1 * magS[c] * ec;
                Sp[c] = C2_ * ec / (kc + std::sqrt(nu_ * ec));
            } else if (variant_ == Variant::RNG) {
                const scalar eta = magS[c] * kc / ec;
                const scalar R = eta * (1 - eta / eta0_) / (1 + betaRNG_ * eta * eta * eta);
                const scalar C1eff = C1_ - R;
                Su[c] = std::max(C1eff, 0.0) * G[c] * ec / kc;
                Sp[c] = C2_ * ec / kc + std::max(-C1eff, 0.0) * G[c] / kc;
            } else {
                Su[c] = C1_ * G[c] * ec / kc;
                Sp[c] = C2_ * ec / kc;
            }
        }
        solveTransport(e, diffusivity(1.0 / sigmaEps_), Sp, Su, fixCells, fixVals);
        bound(e, epsMin_);

        // k 方程
        for (label c = 0; c < nC; ++c) {
            Su[c] = G[c];
            Sp[c] = e[c] / std::max(k[c], kMin_);
        }
        solveTransport(k, diffusivity(1.0 / sigmak_), Sp, Su);
        bound(k, kMin_);

        if (realizable_) {
            const auto r = rCmu(magS);
            correctNut(&r);
        } else {
            correctNut(nullptr);
        }
        wallShear(&k);
    }

    Variant variant_;
    bool realizable_ = false;
    ScalarField* k_ = nullptr;
    ScalarField* eps_ = nullptr;
    const std::vector<scalar>* y_ = nullptr;
    scalar CmuNut_ = 0.09, C1_ = 1.44, C2_ = 1.92, sigmak_ = 1, sigmaEps_ = 1.3, A0_ = 4;
    scalar eta0_ = 4.38, betaRNG_ = 0.012;
    scalar Cl_ = 2.495, Amu_ = 70, Aeps_ = 4.99, ReyStar_ = 200, ReyWidth_ = 20;
    scalar kMin_ = kSmall, epsMin_ = kSmall;
};

// ------------------------------------------------------------------ k-ω（Wilcox 1998）与 k-ω SST（Menter 2003）
// Wallin & Johansson (2000) 显式代数雷诺应力模型（EARSM）的各向异性张量 a_ij（τ_ij = k(2/3 δ_ij + a_ij)）。
// S* = τ S，Ω* = τ Ω（τ 为湍流时间尺度）；N 取二维解析解（Menter-Garbaruk-Egorov 2009 的 BSL-EARSM 做法，
// 三维流动也用它），β 系数与张量基 T1、T3、T4、T6、T9 见 WJ 2000 式 (28)–(32)。返回 β1（< 0，线性部分
// 对应涡粘 ν_t = −½ β1 τ k）。
static scalar earsmAnisotropy(const Tensor& gU, scalar tau, Tensor& a) {
    const Tensor S = tau * symm(gU);
    const Tensor W = tau * skew(gU);
    const Tensor W2 = dot(W, W);
    const scalar IIS = trace(dot(S, S));
    const scalar IIW = trace(W2);  // ≤ 0
    const scalar IV = trace(dot(S, W2));
    const scalar c1p = 1.8;  // (9/4)(C1 − 1)，C1 = 1.8
    const scalar P1 = (c1p * c1p / 27 + 0.45 * IIS - (2.0 / 3) * IIW) * c1p;
    const scalar P2 = P1 * P1 - std::pow(c1p * c1p / 9 + 0.9 * IIS + (2.0 / 3) * IIW, 3);
    scalar N;
    if (P2 >= 0) {
        const scalar sq = std::sqrt(P2);
        N = c1p / 3 + std::cbrt(P1 + sq) + std::cbrt(P1 - sq);
    } else {
        const scalar r = std::sqrt(P1 * P1 - P2);
        N = c1p / 3 + 2 * std::cbrt(r) * std::cos(std::acos(std::clamp(P1 / r, -1.0, 1.0)) / 3);
    }
    N = std::max(N, c1p);  // 无应变时 N = c1'，解析解不会小于它
    const scalar Q = (5.0 / 6) * (N * N - 2 * IIW) * (2 * N * N - IIW);
    const scalar beta1 = -N * (2 * N * N - 7 * IIW) / Q;
    const scalar beta3 = -12 * IV / (N * Q), beta4 = -2 * (N * N - 2 * IIW) / Q, beta6 = -6 * N / Q, beta9 = 6 / Q;
    const Tensor I = identity();
    const Tensor T3 = W2 - (IIW / 3) * I;
    const Tensor T4 = dot(S, W) - dot(W, S);
    const Tensor T6 = dot(S, W2) + dot(W2, S) - IIW * S - (2.0 / 3) * IV * I;
    const Tensor T9 = dot(dot(W, S), W2) - dot(dot(W2, S), W);
    a = beta1 * S + beta3 * T3 + beta4 * T4 + beta6 * T6 + beta9 * T9;
    return beta1;
}

class KOmega : public RASModel {
public:
    // sst：Menter 2003 SST；earsm：Menter BSL k-ω + Wallin-Johansson EARSM（Fluent/CFX 的 "BSL EARSM"）
    KOmega(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i, bool sst, bool earsm = false)
        : RASModel(f, d, b, i, true), sst_(sst), earsm_(earsm) {
        betaStar_ = coeff("betaStar", 0.09);
        if (earsm) {
            alphaK1_ = coeff("alphaK1", 0.5);  // BSL 常数（与 SST 只差 σ_k1）
            alphaK2_ = coeff("alphaK2", 1.0);
            alphaOmega1_ = coeff("alphaOmega1", 0.5);
            alphaOmega2_ = coeff("alphaOmega2", 0.856);
            gamma1_ = coeff("gamma1", 5.0 / 9.0);
            gamma2_ = coeff("gamma2", 0.44);
            beta1_ = coeff("beta1", 0.075);
            beta2_ = coeff("beta2", 0.0828);
            c1_ = coeff("c1", 10.0);
            Ctau_ = coeff("Ctau", 6.0);
            Cy_ = coeff("Cy", 0.026);
            nonlinear_ = coeff("nonlinear", 1.0);
            nutMax_ = coeff("nutMax", 1.0);
        } else if (sst) {
            alphaK1_ = coeff("alphaK1", 0.85);
            alphaK2_ = coeff("alphaK2", 1.0);
            alphaOmega1_ = coeff("alphaOmega1", 0.5);
            alphaOmega2_ = coeff("alphaOmega2", 0.856);
            gamma1_ = coeff("gamma1", 5.0 / 9.0);
            gamma2_ = coeff("gamma2", 0.44);
            beta1_ = coeff("beta1", 0.075);
            beta2_ = coeff("beta2", 0.0828);
            a1_ = coeff("a1", 0.31);
            b1_ = coeff("b1", 1.0);
            c1_ = coeff("c1", 10.0);
        } else {
            beta1_ = coeff("beta", 0.072);
            gamma1_ = coeff("gamma", 0.52);
            alphaK1_ = coeff("alphaK", 0.5);
            alphaOmega1_ = coeff("alphaOmega", 0.5);
        }
        k_ = &addField("k", "zeroGradient", 1e-4);
        omega_ = &addField("omega", "zeroGradient", 1.0);
        bound(*k_, kSmall);
        bound(*omega_, kSmall);
        y_ = &wallDistance(m_).y;
        F1_.assign(m_.nTotalCells(), 1.0);
        correctNut(std::vector<scalar>(m_.nTotalCells(), 0.0));
    }
    std::string type() const override { return earsm_ ? "BSLEARSM" : sst_ ? "kOmegaSST" : "kOmega"; }

protected:
    static scalar blend(scalar F1, scalar a, scalar b) { return F1 * (a - b) + b; }

    // EARSM：由当前 k、ω、∇U 计算 ν_t（线性部分）、生成项 P_k = −k a:∇U 与显式附加应力 τ_ex = k a + 2 ν_t S 的
    // 散度（写入 flow.turbSource()，作为动量方程显式源 −∇·τ_ex）。所有单元（含幽灵）逐点计算，无需交换。
    void evalEARSM(std::vector<scalar>* Pk) {
        auto& nut = flow_.nut();
        const auto& k = *k_;
        const auto& w = *omega_;
        if (flow_.gradU().empty()) flow_.updateGradU();  // 构造时还没有速度梯度
        const auto& gU = flow_.gradU();
        const label nT = m_.nTotalCells();
        std::vector<Tensor> tauEx(nT, Tensor::zero());
        if (Pk) Pk->assign(nT, 0.0);
        for (label c = 0; c < nT; ++c) {
            const scalar kc = std::max(k[c], kSmall), wc = std::max(w[c], kSmall);
            const scalar tau = std::max(1.0 / (betaStar_ * wc), Ctau_ * std::sqrt(nu_ / (betaStar_ * kc * wc)));
            Tensor a;
            const scalar beta1 = earsmAnisotropy(gU[c], tau, a);
            // 近壁阻尼（Wallin & Johansson 2000）：f1 = 1 − exp(−C_y Re_y)，Re_y = √k y/ν；对数层以上 f1 → 1
            const scalar Rey = std::sqrt(kc) * (*y_)[c] / nu_;
            const scalar f1 = Cy_ > 0 ? 1.0 - std::exp(-Cy_ * Rey) : 1.0;
            a = f1 * a;
            nut[c] = std::max(-0.5 * f1 * beta1 * tau * kc, 0.0);
            // 可选上限 ν_t ≤ nutMax·k/ω：WJ 的 C_μ^eff 在弱应变区可达 0.33（k-ω 单位下为 3.7 k/ω），
            // 分离回流区内混合偏强；nutMax = 1 把它限制到 k-ω 的水平（0 为不限制）。整个 a 同比例缩小，
            // 使总雷诺应力一起受限（只限 ν_t 会把差额转到显式项，总应力不变）
            if (nutMax_ > 0 && nut[c] > nutMax_ * kc / wc) {
                const scalar sc = nutMax_ * kc / wc / nut[c];
                a = sc * a;
                nut[c] *= sc;
            }
            // 可实现性：|a_ij| ≤ 1（分量上界 4/3，下界 −2/3）
            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q) a(p, q) = std::clamp(a(p, q), -1.0, 1.0);
            if (Pk) (*Pk)[c] = -kc * doubleDot(a, gU[c]);
            // dev(R) = k a；nonlinear = 0 时只保留线性涡粘部分（诊断用：dev(R) = −2ν_t S，源项为零）
            tauEx[c] = nonlinear_ * (kc * a) - (1 - nonlinear_) * (2 * nut[c] * symm(gU[c]));
        }
        setStressSource(tauEx);
    }

    void correctNut(const std::vector<scalar>& S2) {
        if (earsm_) {
            evalEARSM(nullptr);
            return;
        }
        auto& nut = flow_.nut();
        const auto& k = *k_;
        const auto& w = *omega_;
        const auto& y = *y_;
        for (label c = 0; c < m_.nTotalCells(); ++c) {
            if (!sst_) {
                nut[c] = k[c] / std::max(w[c], kSmall);
                continue;
            }
            const scalar yc = y[c];
            const scalar arg2 = std::min(std::max(2 * std::sqrt(k[c]) / (betaStar_ * w[c] * yc),
                                                  500 * nu_ / (yc * yc * w[c])),
                                         100.0);
            const scalar F2 = std::tanh(arg2 * arg2);
            nut[c] = a1_ * k[c] / std::max(a1_ * w[c], b1_ * F2 * std::sqrt(S2[c]));
        }
    }

    void correctModel() override {
        auto& k = *k_;
        auto& w = *omega_;
        const auto& y = *y_;
        const auto& nut = flow_.nut();
        const label nC = m_.nCells(), nT = m_.nTotalCells();
        const auto gbn0 = GbyNu();
        const auto& gU = flow_.gradU();
        std::vector<scalar> S2(nT, 0.0), G(nT, 0.0), PkE;
        if (earsm_) evalEARSM(&PkE);  // 同时更新 ν_t 与动量源
        for (label c = 0; c < nC; ++c) {
            S2[c] = 2 * magSqr(symm(gU[c]));
            G[c] = earsm_ ? PkE[c] : nut[c] * gbn0[c];
        }
        m_.halo().exchange(S2);

        // 壁面：ω 为粘性/对数层混合值（Menter）；壁面单元 G 按壁面处理方式替换（见 wallProduction）
        wallShear(&k);
        const scalar Cmu25 = std::pow(betaStar_, 0.25);
        const auto wt = wallWeights();
        std::vector<scalar> wW(nT, 0.0), GW(nT, 0.0);
        std::vector<char> isWall(nT, 0), logG(nT, 0);
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const auto& f = wallFaces_[i];
            const scalar kc = std::max(k[f.cell], 0.0);
            const scalar wVis = 6 * nu_ / (beta1_ * f.y * f.y);
            const scalar wLog = std::sqrt(kc) / (Cmu25 * kappa_ * f.y);
            isWall[f.cell] = 1;
            wW[f.cell] += wt[i] * std::sqrt(wVis * wVis + wLog * wLog);
        }
        wallProduction(GW, logG, k);
        std::vector<label> fixCells;
        std::vector<scalar> fixVals;
        for (label c = 0; c < nC; ++c)
            if (isWall[c]) {
                if (logG[c]) G[c] = GW[c];
                fixCells.push_back(c);
                fixVals.push_back(wW[c]);
            }

        std::vector<scalar> Sp(nT, 0.0), Su(nT, 0.0), Dw(nT), Dk(nT);
        if (sst_ || earsm_) {
            // 交叉扩散与混合函数 F1
            const auto gk = fvc::grad(k), gw = fvc::grad(w);
            std::vector<scalar> CD(nT, 0.0);
            for (label c = 0; c < nC; ++c) {
                CD[c] = 2 * alphaOmega2_ * dot(gk[c], gw[c]) / std::max(w[c], kSmall);
                const scalar yc = y[c];
                const scalar CDp = std::max(CD[c], 1e-10);
                const scalar arg1 = std::min(
                    std::min(std::max(std::sqrt(std::max(k[c], 0.0)) / (betaStar_ * w[c] * yc), 500 * nu_ / (yc * yc * w[c])),
                             4 * alphaOmega2_ * k[c] / (CDp * yc * yc)),
                    10.0);
                F1_[c] = std::tanh(std::pow(arg1, 4));
            }
            m_.halo().exchange(F1_);
            for (label c = 0; c < nT; ++c) {
                Dk[c] = blend(F1_[c], alphaK1_, alphaK2_) * nut[c] + nu_;
                Dw[c] = blend(F1_[c], alphaOmega1_, alphaOmega2_) * nut[c] + nu_;
            }
            for (label c = 0; c < nC; ++c) {
                const scalar gamma = blend(F1_[c], gamma1_, gamma2_), beta = blend(F1_[c], beta1_, beta2_);
                const scalar wc = std::max(w[c], kSmall);
                // 源项 γ min(G/ν_t, c1/a1 β* ω max(a1 ω, b1 F2 √S2))：F2 与 ν_t 一致
                const scalar yc = y[c];
                const scalar arg2 = std::min(std::max(2 * std::sqrt(std::max(k[c], 0.0)) / (betaStar_ * wc * yc),
                                                      500 * nu_ / (yc * yc * wc)),
                                             100.0);
                const scalar F2 = std::tanh(arg2 * arg2);
                const scalar GbyNu = std::min(gbn0[c], (c1_ / a1_) * betaStar_ * wc * std::max(a1_ * wc, b1_ * F2 * std::sqrt(S2[c])));
                // 壁面单元用修正后的 G（对数律）换算；EARSM 用 γ ω/k P_k（限 c1 β* ω²）
                const scalar gbn = earsm_ ? std::min(G[c], c1_ * betaStar_ * std::max(k[c], 0.0) * wc) * wc / std::max(k[c], kSmall)
                                 : logG[c] ? G[c] / std::max(nut[c], kSmall) : GbyNu;
                Su[c] = gamma * gbn;
                Sp[c] = beta * wc;
                // −(F1 − 1) CD/ω · ω：系数为负时作显式源项，为正时作隐式汇项
                const scalar cd = (F1_[c] - 1) * CD[c] / wc;
                if (cd > 0) Sp[c] += cd;
                else Su[c] -= cd * wc;
            }
        } else {
            for (label c = 0; c < nT; ++c) {
                Dk[c] = alphaK1_ * nut[c] + nu_;
                Dw[c] = alphaOmega1_ * nut[c] + nu_;
            }
            for (label c = 0; c < nC; ++c) {
                const scalar wc = std::max(w[c], kSmall);
                Su[c] = gamma1_ * G[c] * wc / std::max(k[c], kSmall);
                Sp[c] = beta1_ * wc;
            }
        }
        solveTransport(w, Dw, Sp, Su, fixCells, fixVals);
        bound(w, kSmall);

        for (label c = 0; c < nC; ++c) {
            Su[c] = (sst_ || earsm_) ? std::min(G[c], c1_ * betaStar_ * std::max(k[c], 0.0) * w[c]) : G[c];
            Sp[c] = betaStar_ * w[c];
        }
        solveTransport(k, Dk, Sp, Su);
        bound(k, kSmall);

        correctNut(S2);
        wallShear(&k);
    }

    bool sst_, earsm_ = false;
    scalar Ctau_ = 6.0, Cy_ = 0.026, nonlinear_ = 1.0, nutMax_ = 1.0;
    ScalarField* k_ = nullptr;
    ScalarField* omega_ = nullptr;
    const std::vector<scalar>* y_ = nullptr;
    std::vector<scalar> F1_;
    scalar betaStar_ = 0.09, alphaK1_ = 0.85, alphaK2_ = 1, alphaOmega1_ = 0.5, alphaOmega2_ = 0.856;
    scalar gamma1_ = 5.0 / 9, gamma2_ = 0.44, beta1_ = 0.075, beta2_ = 0.0828, a1_ = 0.31, b1_ = 1, c1_ = 10;
};

// ------------------------------------------------------------------ 微分雷诺应力模型（LRR / SSG，ε 方程）
// 输运 6 个雷诺应力分量 R_ij（Rxx Rxy Rxz Ryy Ryz Rzz，运动学）与 ε：
//   ∂R/∂t + ∇·(φR) − ∇·((ν + ν_t/σ_R)∇R) = P − (2/3)ε I + Φ，P_ij = −(R_ik ∂U_j/∂x_k + R_jk ∂U_i/∂x_k)
//   LRR（Launder-Reece-Rodi 1975，IP 形式，无壁面反射）：Φ = −C1 (ε/k) dev(R) − C2 dev(P)
//   SSG（Speziale-Sarkar-Gatski 1991）：b = dev(R)/(2k)，S = symm(∇U)，Ω = skew(∇U)
//       Φ = −(C1 ε + C1* P_k) b + C2 ε dev(b·b) + k[(C3 − C3*√(b:b)) dev(S) + C4 dev(b·S + S·b) + C5 (b·Ω + (b·Ω)ᵀ)]
//   ε 方程同 k-ε（P_k = ½ tr P，ν_t = C_μ k²/ε 仅用于扩散、壁面函数与动量方程的隐式稳定项）。
// 动量方程：隐式 ν_t 拉普拉斯 + 显式源 −∇·(dev(R) + 2ν_t S)（与 EARSM 同一路径，总应力为 dev(R)）。
// 壁面：R 零梯度（OpenFOAM kqRWallFunction），壁面单元 ε 与产生项按所选壁面处理（standard / enhanced），k = ½ tr R。
//   wallAnisotropy（默认开）：对数层壁面单元的雷诺应力按近壁平衡各向异性给定 R_ij = a_ij k（Fluent 壁面函数处理：
//   uu/k 1.098、nn/k 0.247、ss/k 0.655、−un/k 0.255），其中壁面单元的 k 由标准 k 方程得到（P_k 按壁面律、ε 为壁面值、
//   扩散通向内部的 ½ tr R）；该 k 场（场名 k）在内部单元即 ½ tr R。R 分量方程在这些单元里直接固定为该值，
//   使求解与边值一致、残差可收敛。
// 系数同 OpenFOAM（LRR：C1 1.8，C2 0.6，Cε1 1.44，Cε2 1.92；SSG：C1 3.4，C1* 1.8，C2 4.2，C3 0.8，C3* 1.3，C4 1.25，C5 0.4，
// Cε1 1.44，Cε2 1.83）；扩散用各向同性 σ_R 0.82（Fluent 标量选项）、σ_ε 1.0。
class ReynoldsStress : public RASModel {
public:
    ReynoldsStress(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i, bool ssg)
        : RASModel(f, d, b, i, true), ssg_(ssg) {
        if (!wallFunctions_) throw std::runtime_error("Reynolds stress models need wall functions (wallFunctions: true)");
        if (ssg) {
            C1_ = coeff("C1", 3.4);
            C1s_ = coeff("C1s", 1.8);
            C2_ = coeff("C2", 4.2);
            C3_ = coeff("C3", 0.8);
            C3s_ = coeff("C3s", 1.3);
            C4_ = coeff("C4", 1.25);
            C5_ = coeff("C5", 0.4);
            Ceps2_ = coeff("Ceps2", 1.83);
        } else {
            C1_ = coeff("C1", 1.8);
            C2_ = coeff("C2", 0.6);
            Ceps2_ = coeff("Ceps2", 1.92);
        }
        Ceps1_ = coeff("Ceps1", 1.44);
        sigmaR_ = coeff("sigmaR", 0.82);
        sigmaEps_ = coeff("sigmaEps", 1.0);
        relaxR_ = d.get("relaxR", flow_.controls.steady ? 0.5 : 1.0);
        wallAnisotropy_ = d.get("wallAnisotropy", true);
        const scalar k0 = initial_.get("k", 1e-4);
        const char* names[6] = {"Rxx", "Rxy", "Rxz", "Ryy", "Ryz", "Rzz"};
        const bool diag[6] = {true, false, false, true, false, true};
        for (int c = 0; c < 6; ++c) R_[c] = &addField(names[c], "zeroGradient", diag[c] ? 2.0 / 3 * k0 : 0.0);
        eps_ = &addField("epsilon", "zeroGradient", 1e-5);
        bound(*eps_, epsMin_);
        if (wallAnisotropy_) kAux_ = &addField("k", "zeroGradient", k0);
        realizable();
        correctNut();
    }
    std::string type() const override { return ssg_ ? "SSG" : "LRR"; }

protected:
    Tensor Rcell(label c) const {
        const scalar xx = (*R_[0])[c], xy = (*R_[1])[c], xz = (*R_[2])[c], yy = (*R_[3])[c], yz = (*R_[4])[c],
                     zz = (*R_[5])[c];
        return {xx, xy, xz, xy, yy, yz, xz, yz, zz};
    }
    scalar kCell(label c) const { return 0.5 * ((*R_[0])[c] + (*R_[3])[c] + (*R_[5])[c]); }
    // 可实现性：对角分量 ≥ kMin，|R_ij| ≤ √(R_ii R_jj)
    void realizable() {
        const int di[3] = {0, 3, 5};
        const int off[3][3] = {{1, 0, 3}, {2, 0, 5}, {4, 3, 5}};
        for (label c = 0; c < m_.nTotalCells(); ++c) {
            for (int k = 0; k < 3; ++k) (*R_[di[k]])[c] = std::max((*R_[di[k]])[c], 2.0 / 3 * kMin_);
            for (int k = 0; k < 3; ++k) {
                const scalar lim = std::sqrt((*R_[off[k][1]])[c] * (*R_[off[k][2]])[c]);
                (*R_[off[k][0]])[c] = std::clamp((*R_[off[k][0]])[c], -lim, lim);
            }
        }
        for (auto* f : R_) f->correctBoundaryConditions();
    }
    void correctNut() {
        auto& nut = flow_.nut();
        for (label c = 0; c < m_.nTotalCells(); ++c) {
            const scalar kc = std::max(kCell(c), kMin_);
            nut[c] = Cmu_ * kc * kc / std::max((*eps_)[c], epsMin_);
        }
    }

    void correctModel() override {
        auto& e = *eps_;
        const label nC = m_.nCells(), nT = m_.nTotalCells();
        if (flow_.gradU().empty()) flow_.updateGradU();
        const auto& gU = flow_.gradU();
        // k 场（供壁面处理与扩散用）
        ScalarField kField(flow_.meshPtr(), "k");
        for (label c = 0; c < nT; ++c) kField[c] = std::max(kCell(c), kMin_);
        // 产生项张量与 P_k
        std::vector<Tensor> P(nT, Tensor::zero());
        std::vector<scalar> Pk(nT, 0.0);
        for (label c = 0; c < nT; ++c) {
            const Tensor RG = dot(Rcell(c), gU[c]);
            P[c] = -1.0 * (RG + transpose(RG));
            Pk[c] = 0.5 * trace(P[c]);
        }
        // 壁面处理：壁面切应力、壁面单元 ε 与 P_k（同 k-ε）
        wallShear(&kField);
        const scalar Cmu75 = std::pow(Cmu_, 0.75);
        const auto wt = wallWeights();
        std::vector<scalar> epsW(nT, 0.0), GW(nT, 0.0);
        std::vector<char> isWall(nT, 0), useLog(nT, 0);
        const bool enhanced = wallTreatment_ == WallTreatment::Enhanced;
        const scalar Cl = kappa_ * std::pow(Cmu_, -0.75), Aeps = 2 * Cl;
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const auto& w = wallFaces_[i];
            const scalar kc = kField[w.cell];
            isWall[w.cell] = 1;
            if (enhanced) {
                const scalar Rey = std::sqrt(kc) * w.y / nu_;
                epsW[w.cell] += wt[i] * std::pow(kc, 1.5) / std::max(Cl * w.y * (1 - std::exp(-Rey / Aeps)), 1e-300);
            } else if (yPlusWall_[i] < yPlusLam_) {
                epsW[w.cell] += wt[i] * 2.0 * kc * nu_ / (w.y * w.y);
            } else {
                epsW[w.cell] += wt[i] * Cmu75 * std::pow(kc, 1.5) / (kappa_ * w.y);
            }
        }
        wallProduction(GW, useLog, kField);
        std::vector<label> fixCells;
        std::vector<scalar> fixVals;
        std::vector<scalar> PkEps = Pk;  // ε 方程用的 P_k（壁面单元按壁面律）
        for (label c = 0; c < nC; ++c)
            if (isWall[c]) {
                if (useLog[c]) {
                    PkEps[c] = GW[c];
                    // 壁面单元的解析速度梯度比对数律大数倍（壁面面值到单元中心的一侧差分），
                    // 把产生项张量整体缩放到壁面律给出的 P_k（保持张量结构）
                    const scalar fac = std::clamp(GW[c] / std::max(Pk[c], 1e-300), 0.0, 10.0);
                    P[c] = fac * P[c];
                    Pk[c] = GW[c];
                }
                fixCells.push_back(c);
                fixVals.push_back(std::max(epsW[c], epsMin_));
            }

        // ε 方程
        std::vector<scalar> Sp(nT, 0.0), Su(nT, 0.0);
        for (label c = 0; c < nC; ++c) {
            const scalar kc = kField[c], ec = std::max(e[c], epsMin_);
            Su[c] = Ceps1_ * std::max(PkEps[c], 0.0) * ec / kc;
            Sp[c] = Ceps2_ * ec / kc;
        }
        solveTransport(e, diffusivity(1.0 / sigmaEps_), Sp, Su, fixCells, fixVals);
        bound(e, epsMin_);

        // 对数层壁面单元的雷诺应力边值（Fluent 壁面函数处理）：R_ij = a_ij k，
        //   uu/k 1.098、nn/k 0.247、ss/k 0.655、−un/k 0.255（t 切向速度方向，n 指向流体，s = t × n），
        //   k 来自辅助的标准 k 方程（P_k 壁面单元按壁面律，ε 同上）
        std::vector<label> fixR;
        std::vector<Tensor> fixA;
        if (wallAnisotropy_) {
            std::vector<Tensor> A(nT, Tensor::zero());
            std::vector<scalar> wsum(nT, 0.0);
            for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
                const auto& w = wallFaces_[i];
                if (!(yPlusWall_[i] > yPlusLam_) || w.magUt < 1e-12) continue;
                const Vec3 nIn = -1.0 * w.n;
                Vec3 du = flow_.U()[w.cell] - flow_.U().bValue(w.face);
                du -= dot(du, w.n) * w.n;
                const Vec3 t = du / mag(du);
                const Vec3 sv = cross(t, nIn);
                A[w.cell] += wt[i] * (1.098 * outer(t, t) + 0.247 * outer(nIn, nIn) + 0.655 * outer(sv, sv) -
                                      0.255 * (outer(t, nIn) + outer(nIn, t)));
                wsum[w.cell] += wt[i];
            }
            // 壁面单元的 k：标准 k 方程（P_k = 壁面律 G_W，ε = 壁面值，扩散通向内部的 ½ tr R），
            // 内部单元固定为 ½ tr R，因此场 k 处处等于模型湍动能
            auto& ka = *kAux_;
            std::vector<label> fixK;
            std::vector<scalar> fixKv;
            for (label c = 0; c < nC; ++c) {
                const scalar ec = std::max(e[c], epsMin_);
                if (wsum[c] > 0.5) {
                    Su[c] = std::max(PkEps[c], 0.0);
                    Sp[c] = ec / std::max(ka[c], kMin_);
                } else {
                    Su[c] = Sp[c] = 0.0;
                    fixK.push_back(c);
                    fixKv.push_back(kField[c]);
                }
            }
            solveTransport(ka, diffusivity(1.0), Sp, Su, fixK, fixKv);
            bound(ka, kMin_);
            for (label c = 0; c < nC; ++c)
                if (wsum[c] > 0.5) {
                    fixR.push_back(c);
                    fixA.push_back((std::max(ka[c], kMin_) / wsum[c]) * A[c]);
                }
        }

        // 雷诺应力方程：各分量用同一隐式系数 Sp，显式源 Su_ij；对数层壁面单元固定为上述边值
        const auto DR = diffusivity(1.0 / sigmaR_);
        std::vector<Tensor> Sex(nT, Tensor::zero());
        std::vector<scalar> SpR(nT, 0.0);
        const Tensor I = identity();
        for (label c = 0; c < nC; ++c) {
            const scalar kc = kField[c], ec = std::max(e[c], epsMin_);
            const Tensor R = Rcell(c);
            if (ssg_) {
                const Tensor b = dev(R) / (2 * kc);
                const Tensor S = symm(gU[c]), Om = skew(gU[c]);
                const scalar magb = std::sqrt(magSqr(b));
                const Tensor bS = dot(b, S), bOm = dot(b, Om);
                // −(C1 ε + C1* P_k) b = −(C1 ε + C1* P_k)/(2k) R + (C1 ε + C1* P_k)/3 I
                const scalar cImp = (C1_ * ec + C1s_ * std::max(Pk[c], 0.0)) / (2 * kc);
                SpR[c] = cImp;
                Sex[c] = P[c] + (-(2.0 / 3) * ec + (C1_ * ec + C1s_ * std::max(Pk[c], 0.0)) / 3) * I +
                         C2_ * ec * dev(dot(b, b)) +
                         kc * ((C3_ - C3s_ * magb) * dev(S) + C4_ * dev(bS + transpose(bS)) + C5_ * (bOm + transpose(bOm)));
                if (Pk[c] < 0) Sex[c] += (C1s_ * Pk[c]) * b * (-1.0);  // P_k < 0 时的 C1* 项显式处理
            } else {
                // −C1 (ε/k) dev(R) = −C1 ε/k R + (2/3) C1 ε I；−C2 dev(P)
                SpR[c] = C1_ * ec / kc;
                Sex[c] = P[c] - C2_ * dev(P[c]) + ((2.0 / 3) * (C1_ - 1) * ec) * I;
            }
        }
        const scalar relaxSave = relax_;
        relax_ = relaxR_;
        std::vector<scalar> fixRv(fixR.size());
        for (int comp = 0; comp < 6; ++comp) {
            ScalarField& Rc = *R_[comp];
            for (label c = 0; c < nC; ++c) {
                Su[c] = Traits<Tensor>::component(Sex[c], compIndex(comp));
                Sp[c] = SpR[c];
            }
            for (std::size_t j = 0; j < fixR.size(); ++j) fixRv[j] = Traits<Tensor>::component(fixA[j], compIndex(comp));
            solveTransport(Rc, DR, Sp, Su, fixR, fixRv);
        }
        relax_ = relaxSave;
        realizable();
        correctNut();
        // 显式应力源：dev(R) 减去隐式涡粘部分（见 setStressSource）
        std::vector<Tensor> devR(nT, Tensor::zero());
        for (label c = 0; c < nT; ++c) devR[c] = dev(Rcell(c));
        setStressSource(devR);
        for (label c = 0; c < nT; ++c) kField[c] = std::max(kCell(c), kMin_);
        wallShear(&kField);
    }
    // 分量序号 → Tensor 行主序索引（xx xy xz yy yz zz → 0 1 2 4 5 8）
    static int compIndex(int comp) {
        static const int idx[6] = {0, 1, 2, 4, 5, 8};
        return idx[comp];
    }

    bool ssg_;
    ScalarField* R_[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    ScalarField* eps_ = nullptr;
    ScalarField* kAux_ = nullptr;  // 辅助 k（仅用于壁面单元 R 边值）
    scalar C1_ = 1.8, C1s_ = 1.8, C2_ = 0.6, C3_ = 0.8, C3s_ = 1.3, C4_ = 1.25, C5_ = 0.4;
    scalar Ceps1_ = 1.44, Ceps2_ = 1.92, sigmaR_ = 0.82, sigmaEps_ = 1.0, relaxR_ = 0.5;
    bool wallAnisotropy_ = true;
    scalar kMin_ = kSmall, epsMin_ = kSmall;
};

// ------------------------------------------------------------------ Spalart-Allmaras（无 ft2；可选 Spalding 壁面函数）
class SpalartAllmaras : public RASModel {
public:
    SpalartAllmaras(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i) : RASModel(f, d, b, i, false) {
        sigmaNut_ = coeff("sigmaNut", 0.66666);
        Cb1_ = coeff("Cb1", 0.1355);
        Cb2_ = coeff("Cb2", 0.622);
        Cw2_ = coeff("Cw2", 0.3);
        Cw3_ = coeff("Cw3", 2.0);
        Cv1_ = coeff("Cv1", 7.1);
        Cs_ = coeff("Cs", 0.3);
        Cw1_ = Cb1_ / (kappa_ * kappa_) + (1 + Cb2_) / sigmaNut_;
        nuTilda_ = &addField("nuTilda", "zero", 3 * nu_);
        bound(*nuTilda_, 0.0);
        y_ = &wallDistance(m_).y;
        correctNut();
    }
    std::string type() const override { return "SpalartAllmaras"; }

protected:
    scalar fv1(scalar chi) const {
        const scalar c3 = chi * chi * chi;
        return c3 / (c3 + std::pow(Cv1_, 3));
    }
    void correctNut() {
        auto& nut = flow_.nut();
        const auto& nt = *nuTilda_;
        for (label c = 0; c < m_.nTotalCells(); ++c) nut[c] = nt[c] * fv1(nt[c] / nu_);
    }
    void correctModel() override {
        auto& nt = *nuTilda_;
        const auto& y = *y_;
        const auto& gU = flow_.gradU();
        const label nC = m_.nCells(), nT = m_.nTotalCells();
        const auto gnt = fvc::grad(nt);
        std::vector<scalar> Sp(nT, 0.0), Su(nT, 0.0), D(nT);
        for (label c = 0; c < nT; ++c) D[c] = (std::max(nt[c], 0.0) + nu_) / sigmaNut_;
        for (label c = 0; c < nC; ++c) {
            const scalar ntc = std::max(nt[c], 0.0);
            const scalar chi = ntc / nu_;
            const scalar f1 = fv1(chi);
            const scalar fv2 = 1 - chi / (1 + chi * f1);
            const scalar Omega = std::sqrt(2.0) * mag(skew(gU[c]));
            const scalar ky2 = kappa_ * kappa_ * y[c] * y[c];
            const scalar Stilda = std::max(Omega + fv2 * ntc / ky2, Cs_ * Omega);
            const scalar r = std::min(ntc / (std::max(Stilda, 1e-300) * ky2), 10.0);
            const scalar g = r + Cw2_ * (std::pow(r, 6) - r);
            const scalar c6 = std::pow(Cw3_, 6);
            const scalar fw = g * std::pow((1 + c6) / (std::pow(g, 6) + c6), 1.0 / 6.0);
            Su[c] = Cb1_ * Stilda * ntc + Cb2_ / sigmaNut_ * magSqr(gnt[c]);
            Sp[c] = Cw1_ * fw * ntc / (y[c] * y[c]);
        }
        solveTransport(nt, D, Sp, Su);
        bound(nt, 0.0);
        correctNut();
        // 壁面：wallFunctions=false（默认）为低 Re 处理（壁面 ν_t = 0）；true 时 u_τ 由 Spalding 律给出壁面切应力
        wallShear(nullptr);
    }

    ScalarField* nuTilda_ = nullptr;
    const std::vector<scalar>* y_ = nullptr;
    scalar sigmaNut_, Cb1_, Cb2_, Cw1_, Cw2_, Cw3_, Cv1_, Cs_;
};

} // namespace

std::unique_ptr<TurbulenceModel> makeRAS(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i,
                                         const std::string& type) {
    if (type == "kEpsilon") return std::make_unique<KEpsilon>(flow, d, b, i, KEpsilon::Variant::Standard);
    if (type == "RNGkEpsilon") return std::make_unique<KEpsilon>(flow, d, b, i, KEpsilon::Variant::RNG);
    if (type == "realizableKE") return std::make_unique<KEpsilon>(flow, d, b, i, KEpsilon::Variant::Realizable);
    if (type == "kOmega") return std::make_unique<KOmega>(flow, d, b, i, false);
    if (type == "kOmegaSST") return std::make_unique<KOmega>(flow, d, b, i, true);
    if (type == "BSLEARSM") return std::make_unique<KOmega>(flow, d, b, i, false, true);
    if (type == "LRR") return std::make_unique<ReynoldsStress>(flow, d, b, i, false);
    if (type == "SSG") return std::make_unique<ReynoldsStress>(flow, d, b, i, true);
    return std::make_unique<SpalartAllmaras>(flow, d, b, i);
}

} // namespace cfd
