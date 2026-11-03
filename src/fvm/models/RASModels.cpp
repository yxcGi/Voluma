// RANS 模型：kEpsilon、realizableKE、kOmega、kOmegaSST、SpalartAllmaras（公式与系数同 OpenFOAM）。
#include "fvm/discretization/Fvc.h"
#include "fvm/mesh/WallDistance.h"
#include "fvm/models/TurbulenceModel.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace cfd {

namespace {

constexpr scalar kSmall = 1e-15;

class RASModel : public TurbulenceModel {
public:
    RASModel(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i, bool wallFunctionsDefault)
        : TurbulenceModel(f, d, b, i) {
        wallFunctions_ = d.get("wallFunctions", wallFunctionsDefault);
    }

protected:
    // 壁面单元：每个壁面面的权重 1/（该单元壁面面数）
    struct WallCellWeights {
        std::vector<scalar> w;  // 按 wallFaces_
    };
    std::vector<scalar> wallWeights() const {
        std::vector<int> cnt(m_.nTotalCells(), 0);
        for (auto& w : wallFaces_) ++cnt[w.cell];
        std::vector<scalar> wt(wallFaces_.size());
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) wt[i] = 1.0 / cnt[wallFaces_[i].cell];
        return wt;
    }
    // nutkWallFunction：y+ = C_μ^¼ √k y/ν，y+ > y+_lam 时 ν_t,w = ν(κ y+/ln(E y+) − 1)
    void kWallNut(const ScalarField& k) {
        const scalar Cmu25 = std::pow(Cmu_, 0.25);
        std::vector<scalar> nutw(wallFaces_.size(), 0.0);
        uTauWall_.assign(wallFaces_.size(), 0.0);
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const auto& w = wallFaces_[i];
            const scalar yPlus = Cmu25 * std::sqrt(std::max(k[w.cell], 0.0)) * w.y / nu_;
            if (wallFunctions_ && yPlus > yPlusLam_) nutw[i] = nu_ * (yPlus * kappa_ / std::log(E_ * yPlus) - 1.0);
            // 输出用 u_τ：由壁面切应力
            const scalar tau = (nu_ + nutw[i]) * w.magUt / w.y;
            uTauWall_[i] = std::sqrt(std::max(tau, 0.0));
        }
        setWallNut(nutw);
    }
    std::vector<scalar> GbyNu() const {
        const auto& G = flow_.gradU();
        std::vector<scalar> g(m_.nTotalCells(), 0.0);
        for (label c = 0; c < m_.nCells(); ++c) g[c] = doubleDot(dev(2.0 * symm(G[c])), G[c]);
        return g;
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
};

// ------------------------------------------------------------------ k-ε（标准 / realizable）
class KEpsilon : public RASModel {
public:
    KEpsilon(IncompressibleFlow& f, const Json& d, const Json& b, const Json& i, bool realizable)
        : RASModel(f, d, b, i, true), realizable_(realizable) {
        if (!wallFunctions_) throw std::runtime_error("kEpsilon models need wall functions (wallFunctions: true)");
        if (realizable) {
            A0_ = coeff("A0", 4.0);
            C2_ = coeff("C2", 1.9);
            sigmak_ = coeff("sigmak", 1.0);
            sigmaEps_ = coeff("sigmaEps", 1.2);
        } else {
            Cmu_ = coeff("Cmu", 0.09);
            C1_ = coeff("C1", 1.44);
            C2_ = coeff("C2", 1.92);
            sigmak_ = coeff("sigmak", 1.0);
            sigmaEps_ = coeff("sigmaEps", 1.3);
        }
        k_ = &addField("k", "zeroGradient", 1e-4);
        eps_ = &addField("epsilon", "zeroGradient", 1e-5);
        bound(*k_, kMin_);
        bound(*eps_, epsMin_);
        correctNut(nullptr);
    }
    std::string type() const override { return realizable_ ? "realizableKE" : "kEpsilon"; }

protected:
    void correctNut(const std::vector<scalar>* rCmu) {
        auto& nut = flow_.nut();
        const auto& k = *k_;
        const auto& e = *eps_;
        for (label c = 0; c < m_.nTotalCells(); ++c)
            nut[c] = (rCmu ? (*rCmu)[c] : Cmu_) * k[c] * k[c] / std::max(e[c], epsMin_);
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
        if (realizable_) {
            const auto& gU = flow_.gradU();
            for (label c = 0; c < nC; ++c) magS[c] = std::sqrt(2 * magSqr(dev(symm(gU[c]))));
        }

        // 壁面函数：壁面单元的 ε 与 G
        kWallNut(k);
        const scalar Cmu25 = std::pow(Cmu_, 0.25), Cmu75 = std::pow(Cmu_, 0.75);
        const auto wt = wallWeights();
        std::vector<scalar> epsW(m_.nTotalCells(), 0.0), GW(m_.nTotalCells(), 0.0);
        std::vector<char> isWall(m_.nTotalCells(), 0);
        const auto& wn = flow_.wallNuEff();
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const auto& w = wallFaces_[i];
            const scalar kc = std::max(k[w.cell], 0.0);
            const scalar yPlus = Cmu25 * std::sqrt(kc) * w.y / nu_;
            isWall[w.cell] = 1;
            if (yPlus < yPlusLam_) {
                epsW[w.cell] += wt[i] * 2.0 * kc * nu_ / (w.y * w.y);
            } else {
                epsW[w.cell] += wt[i] * Cmu75 * std::pow(kc, 1.5) / (kappa_ * w.y);
                const scalar nuw = wn[w.face - m_.nInternalFaces()];
                GW[w.cell] += wt[i] * nuw * (w.magUt / w.y) * Cmu25 * std::sqrt(kc) / (kappa_ * w.y);
            }
        }
        std::vector<label> fixCells;
        std::vector<scalar> fixVals;
        for (label c = 0; c < nC; ++c)
            if (isWall[c]) {
                G[c] = GW[c] > 0 ? GW[c] : G[c];
                fixCells.push_back(c);
                fixVals.push_back(std::max(epsW[c], epsMin_));
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
        kWallNut(k);
    }

    bool realizable_;
    ScalarField* k_ = nullptr;
    ScalarField* eps_ = nullptr;
    scalar C1_ = 1.44, C2_ = 1.92, sigmak_ = 1, sigmaEps_ = 1.3, A0_ = 4;
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
            tauEx[c] = nonlinear_ * (kc * a + 2 * nut[c] * symm(gU[c]));  // nonlinear = 0 退化为线性涡粘（诊断用）
        }
        // −∇·τ_ex：面插值后与面矢量点乘；壁面上雷诺应力为零
        const auto tf = fvc::interpolate(m_, tauEx);
        auto& src = flow_.turbSource();
        src.assign(nT, Vec3{0, 0, 0});
        const label nI = m_.nInternalFaces();
        for (label f = 0; f < nI; ++f) {
            const Vec3 F = dot(tf[f], m_.Sf()[f]);
            src[m_.owner()[f]] -= F;
            src[m_.neighbour()[f]] += F;
        }
        for (const auto& pt : m_.patches()) {
            if (pt.type == PatchType::Wall || pt.type == PatchType::Empty) continue;
            for (label f = pt.start; f < pt.end(); ++f) src[m_.owner()[f]] -= dot(tf[f], m_.Sf()[f]);
        }
        for (label c = 0; c < m_.nCells(); ++c) src[c] = src[c] / m_.V()[c];
        if (m_.twoD())
            for (label c = 0; c < m_.nCells(); ++c) src[c][m_.emptyDir()] = 0.0;
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

        // 壁面：ω 为粘性/对数层混合值，y+ > y+_lam 时 G 用对数律
        kWallNut(k);
        const scalar Cmu25 = std::pow(betaStar_, 0.25);
        const auto wt = wallWeights();
        std::vector<scalar> wW(nT, 0.0), GW(nT, 0.0);
        std::vector<char> isWall(nT, 0), logG(nT, 0);
        const auto& wn = flow_.wallNuEff();
        for (std::size_t i = 0; i < wallFaces_.size(); ++i) {
            const auto& f = wallFaces_[i];
            const scalar kc = std::max(k[f.cell], 0.0);
            const scalar yPlus = Cmu25 * std::sqrt(kc) * f.y / nu_;
            const scalar wVis = 6 * nu_ / (beta1_ * f.y * f.y);
            const scalar wLog = std::sqrt(kc) / (Cmu25 * kappa_ * f.y);
            isWall[f.cell] = 1;
            wW[f.cell] += wt[i] * std::sqrt(wVis * wVis + wLog * wLog);
            if (wallFunctions_ && yPlus > yPlusLam_) {
                const scalar nuw = wn[f.face - m_.nInternalFaces()];
                GW[f.cell] += wt[i] * nuw * (f.magUt / f.y) * Cmu25 * std::sqrt(kc) / (kappa_ * f.y);
                logG[f.cell] = 1;
            }
        }
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
        kWallNut(k);
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

// ------------------------------------------------------------------ Spalart-Allmaras（无 ft2）
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
        // 低 Re 处理：壁面 ν_t = 0
        std::vector<scalar> nutw(wallFaces_.size(), 0.0);
        setWallNut(nutw);
        uTauWall_.assign(wallFaces_.size(), 0.0);
        for (std::size_t i = 0; i < wallFaces_.size(); ++i)
            uTauWall_[i] = std::sqrt(nu_ * wallFaces_[i].magUt / wallFaces_[i].y);
    }

    ScalarField* nuTilda_ = nullptr;
    const std::vector<scalar>* y_ = nullptr;
    scalar sigmaNut_, Cb1_, Cb2_, Cw1_, Cw2_, Cw3_, Cv1_, Cs_;
};

} // namespace

std::unique_ptr<TurbulenceModel> makeRAS(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i,
                                         const std::string& type) {
    if (type == "kEpsilon") return std::make_unique<KEpsilon>(flow, d, b, i, false);
    if (type == "realizableKE") return std::make_unique<KEpsilon>(flow, d, b, i, true);
    if (type == "kOmega") return std::make_unique<KOmega>(flow, d, b, i, false);
    if (type == "kOmegaSST") return std::make_unique<KOmega>(flow, d, b, i, true);
    if (type == "BSLEARSM") return std::make_unique<KOmega>(flow, d, b, i, false, true);
    return std::make_unique<SpalartAllmaras>(flow, d, b, i);
}

} // namespace cfd
