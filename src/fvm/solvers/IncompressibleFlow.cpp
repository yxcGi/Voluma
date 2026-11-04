#include "fvm/solvers/IncompressibleFlow.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace cfd {

IncompressibleFlow::IncompressibleFlow(MeshPtr mesh, scalar nu)
    : mesh_(std::move(mesh)), nu_(nu), U_(mesh_, "U"), p_(mesh_, "p"), phi_(mesh_->nFaces(), 0.0) {
    ctx_.phi = &phi_;
    ctx_.U = &U_;
    ctx_.time = &ts_;
    ctx_.nu = nu_;
}

void IncompressibleFlow::constrainPressure(const std::vector<scalar>& phiHbyA, const std::vector<scalar>& rAUf) {
    const Mesh& m = *mesh_;
    for (label p = 0; p < label(m.patches().size()); ++p) {
        auto& bc = p_.bc(p);
        if (bc.type() != "fixedFluxPressure") continue;
        auto& g = static_cast<FixedFluxPressureBC&>(bc).gradient();
        const Patch& pt = m.patches()[p];
        for (label i = 0; i < pt.size; ++i) {
            const label f = pt.start + i;
            g[i] = (phiHbyA[f] - dot(U_.bValue(f), m.Sf()[f])) / (m.magSf()[f] * rAUf[f]);
        }
    }
}

void IncompressibleFlow::initialize() {
    U_.updateBCs(ts_.time);
    p_.updateBCs(ts_.time);
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    phi_ = mesh_->skewed() ? fvc::flux(U_, fvc::grad(U_, gradScheme)) : fvc::flux(U_);
    rAU_.assign(mesh_->nTotalCells(), 0.0);
    initialized_ = true;
}

void IncompressibleFlow::markRestarted() {
    U_.updateBCs(ts_.time);
    p_.updateBCs(ts_.time);
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    rAU_.assign(mesh_->nTotalCells(), 0.0);
    initialized_ = true;
}

void IncompressibleFlow::updateGradU() { gradU_ = fvc::grad(U_, gradScheme); }

void IncompressibleFlow::projectDivergenceFree(int nNonOrthCorr) {
    if (!initialized_) initialize();
    const Mesh& m = *mesh_;
    const label nT = m.nTotalCells();
    // 借用 p_ 作为势函数（边界条件相同），结束后恢复
    const std::vector<scalar> pSaved = p_.internal();
    p_.setUniform(0.0);
    p_.correctBoundaryConditions();
    const bool needRef = pressureNeedsReference();
    const auto divPhi = fvc::div(m, phi_);
    for (int no = 0; no <= nNonOrthCorr; ++no) {
        FvMatrix<scalar> eqn = fvm::laplacian(1.0, p_, laplacianNonOrthCorr);
        eqn -= divPhi;
        if (needRef) eqn.setReference(controls.pRefCell, 0.0);
        eqn.solve(controls.pFinalControls, controls.verbose);
        if (no == nNonOrthCorr) {
            const auto F = eqn.flux();
            for (label f = 0; f < m.nFaces(); ++f) phi_[f] -= F[f];
        }
    }
    // 单元速度由修正后的面通量重构（potentialFoam 同法）：U −= ∇φ 在零梯度边界（入口）单元上的 Gauss 梯度只有一半，
    // 会留下与面通量不一致的单元速度
    const auto Ur = fvc::reconstruct(m, phi_);
    for (label c = 0; c < nT; ++c) U_[c] = Ur[c];
    if (m.twoD())
        for (label c = 0; c < nT; ++c) U_[c][m.emptyDir()] = 0.0;
    U_.correctBoundaryConditions();
    for (label c = 0; c < nT; ++c) p_[c] = pSaved[c];
    p_.correctBoundaryConditions();
}

std::vector<scalar> IncompressibleFlow::nuEffFaces() const {
    const Mesh& m = *mesh_;
    if (nut_.empty()) return std::vector<scalar>(m.nFaces(), nu_);
    std::vector<scalar> nuEff(m.nTotalCells());
    for (label c = 0; c < m.nTotalCells(); ++c) nuEff[c] = nu_ + nut_[c];
    auto f = fvc::interpolate(m, nuEff);
    if (!wallNuEff_.empty()) {
        const label nI = m.nInternalFaces();
        for (label bf = 0; bf < m.nBoundaryFaces(); ++bf)
            if (wallNuEff_[bf] >= 0) f[nI + bf] = wallNuEff_[bf];
    }
    return f;
}

scalar IncompressibleFlow::courantNumber(scalar dt) const {
    const Mesh& m = *mesh_;
    std::vector<scalar> s(m.nTotalCells(), 0.0);
    for (label f = 0; f < m.nFaces(); ++f) {
        s[m.owner()[f]] += std::abs(phi_[f]);
        if (f < m.nInternalFaces()) s[m.neighbour()[f]] += std::abs(phi_[f]);
    }
    scalar co = 0;
    for (label c = 0; c < m.nCells(); ++c) co = std::max(co, 0.5 * s[c] / m.V()[c] * dt);
    return par::allMax(co);
}

bool IncompressibleFlow::pressureNeedsReference() const {
    bool fixes = false;
    for (label p = 0; p < label(mesh_->patches().size()); ++p) fixes = fixes || p_.bc(p).fixesValue();
    return par::allMax(fixes ? 1.0 : 0.0) < 0.5;
}

void IncompressibleFlow::applyForcingCorrection(const std::vector<scalar>& rAU, scalar ddtCoeff, scalar alphaU) {
    if (forcing.mode != MomentumForcing::Mode::MeanVelocity) return;
    const Mesh& m = *mesh_;
    const scalar magUbar = mag(forcing.Ubar);
    const Vec3 dir = forcing.Ubar / magUbar;
    par::SumAcc s[2];
    for (label c = 0; c < m.nCells(); ++c) {
        s[0].add(dot(dir, U_[c]) * m.V()[c]);
        s[1].add(rAU[c] * m.V()[c]);
    }
    double r[2];
    par::SumAcc::allReduce(s, 2, r);
    const scalar magUbarAve = r[0] / m.totalVolume();
    const scalar rAUave = r[1] / m.totalVolume();
    const scalar dU = magUbar - magUbarAve;
    // 与 OpenFOAM meanVelocityForce 相同：按 rAU 分布修正速度，使体平均速度正好等于目标值
    for (label c = 0; c < m.nTotalCells(); ++c) U_[c] += dir * (rAU[c] * dU / rAUave);
    U_.correctBoundaryConditions();
    // 驱动压力梯度的增量按“平均速度对均匀体力的响应” resp 计算，而不是 OpenFOAM 的 rAUave。
    // 均匀增量在对流、扩散项中相互抵消，只有时间项 c/Δt 与欠松弛附加项 (1−α)A 起作用，
    // 所以 resp = 1/(c/Δt + (1−α)·A)。粘性主导（rAU ≪ Δt）时 OpenFOAM 的做法增益过大而振荡发散，
    // 这里的 resp 是响应的上界，增益 ≤ 1，总是稳定。
    const scalar Aave = 1.0 / rAUave;
    const scalar resp = 1.0 / (ddtCoeff + (1.0 - alphaU) * Aave);
    forcing.gradP += forcing.relaxation * dU / resp;
    if (std::getenv("CFD_DEBUG_FORCING"))
        std::cout << "forcing: Uave " << magUbarAve << " gradP " << forcing.gradP << '\n';
}

StepInfo IncompressibleFlow::step(scalar dt) {
    if (!initialized_) initialize();
    const Mesh& m = *mesh_;
    const label nI = m.nInternalFaces(), nC = m.nCells(), nT = m.nTotalCells();
    PimpleControls& ctl = controls;
    const bool steady = ctl.steady;
    const DdtScheme ddtS = steady ? DdtScheme::Steady : ddtScheme;
    ctx_.ddt = ddtS;

    // ---------------------------------------------------------- 新时间层
    if (!steady) {
        ts_.dt0 = ts_.timeIndex > 0 ? ts_.dt : -1.0;
        ts_.dt = dt;
        ts_.time += dt;
        ++ts_.timeIndex;
        U_.storeOld();
        phiOldOld_ = phiOld_;
        phiOld_ = phi_;
    } else {
        ts_.time += 1.0;
        ++ts_.timeIndex;
    }
    U_.updateBCs(ts_.time);
    p_.updateBCs(ts_.time);
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();

    StepInfo info;
    const bool needRef = pressureNeedsReference();

    for (int outer = 0; outer < ctl.nOuter; ++outer) {
        const bool finalOuter = outer == ctl.nOuter - 1;
        if (updateTurbulence) {
            updateGradU();
            updateTurbulence(*this);
        }
        const auto nuEff = nuEffFaces();

        // ------------------------------------------------------ 动量方程
        FvMatrix<Vec3> UEqn = fvm::ddt(U_, ddtS, ts_);
        UEqn += fvm::div(phi_, U_, divScheme);
        UEqn -= fvm::laplacian(nuEff, U_, laplacianNonOrthCorr);
        if (!nut_.empty()) {
            // ∇·(ν_eff dev2((∇U)ᵀ))，显式
            if (gradU_.empty()) updateGradU();
            const auto Gf = fvc::interpolate(m, gradU_);
            std::vector<Vec3> F(m.nFaces());
            for (label f = 0; f < m.nFaces(); ++f) {
                const Tensor& G = Gf[f];
                F[f] = nuEff[f] * (dot(G, m.Sf()[f]) - (2.0 / 3.0) * trace(G) * m.Sf()[f]);
            }
            UEqn -= fvc::surfaceIntegrate(m, F);
        }
        if (!turbSource_.empty()) UEqn -= turbSource_;
        // 单元体力（浮力等）不进入 H：按 OpenFOAM 浮力求解器的做法以面通量 φ_F = f_f·S_f 进入压力方程与速度修正，
        // 预测步用 reconstruct(φ_F − ∇p_f·S_f)，使体力与压力梯度在面上严格平衡（静水平衡不产生伪速度）
        const bool hasBF = !bodyForce_.empty();
        std::vector<scalar> phiF;
        if (hasBF) {
            const auto ff = fvc::interpolate(m, bodyForce_);
            phiF.resize(m.nFaces());
            for (label f = 0; f < m.nFaces(); ++f) phiF[f] = dot(ff[f], m.Sf()[f]);
        }
        if (forcing.mode == MomentumForcing::Mode::Constant) UEqn -= forcing.force;
        if (forcing.mode == MomentumForcing::Mode::MeanVelocity)
            UEqn -= forcing.Ubar * (forcing.gradP / mag(forcing.Ubar));

        const scalar aU = steady ? ctl.alphaU : (finalOuter ? 1.0 : ctl.alphaU);
        UEqn.relax(aU);

        if (ctl.momentumPredictor) {
            FvMatrix<Vec3> M(UEqn);
            std::vector<Vec3> ngp(nT);
            if (hasBF) {
                const auto gpF = fvm::laplacian(1.0, p_, laplacianNonOrthCorr).flux();
                std::vector<scalar> d(m.nFaces());
                for (label f = 0; f < m.nFaces(); ++f) d[f] = phiF[f] - gpF[f];
                ngp = fvc::reconstruct(m, d);
            } else {
                const auto gp = fvc::grad(p_, gradScheme);
                for (label c = 0; c < nT; ++c) ngp[c] = -gp[c];
            }
            M -= ngp;
            auto perf = M.solve(ctl.UControls, ctl.verbose);
            if (outer == 0 && !perf.empty()) {
                // 监控残差：各分量的 Σ|b − A x| 都除以分量中最大的归一化因子（Fluent 式按速度模归一），
                // 否则恒为零的分量（如周期管的横向速度）会得到噪声/噪声 ≈ O(1) 的伪残差
                scalar r = 0, nf = 0;
                for (auto& pf : perf) {
                    r = std::max(r, pf.residualSum);
                    nf = std::max(nf, pf.normFactor);
                }
                info.UInitialResidual = r / (nf + 1e-30);
            }
            info.perf.insert(info.perf.end(), perf.begin(), perf.end());
        }

        // ------------------------------------------------------ 压力修正
        for (int corr = 0; corr < ctl.nCorr; ++corr) {
            const bool finalCorr = finalOuter && corr == ctl.nCorr - 1;
            const auto A = UEqn.A();
            rAU_.assign(nT, 0.0);
            for (label c = 0; c < nT; ++c) rAU_[c] = 1.0 / A[c];
            const auto H = UEqn.H();
            std::vector<Vec3> HbyA(nT);
            for (label c = 0; c < nT; ++c) HbyA[c] = rAU_[c] * H[c];
            // rAU_f：偏斜网格上同样加偏斜修正（∇rAU 按零法向梯度边界的最小二乘），与 HbyA_f 的修正一致，
            // 使线性压力场下 φ_f = U_f·S_f 精确成立（否则 G S_z (x_f − x'_f)·∇rAU 成为伪通量）
            const auto rAUf = m.skewed() ? fvc::interpolate(m, rAU_, fvc::gradLS(m, rAU_, std::vector<char>(m.nBoundaryFaces(), 0), std::vector<scalar>(m.nBoundaryFaces(), 0.0)))
                                         : fvc::interpolate(m, rAU_);

            // 偏斜网格（四面体/三角形）：HbyA 面值加偏斜修正 (x_f − x'_f)·∇HbyA_f。线性插值点不在面心时
            // 面速度有一阶误差，连续性方程会迫使单元速度偏离（四面体圆管壁面切应力偏高约 10%，见 docs/validation.md）。
            // ∇HbyA 用最小二乘：定通量边界面取 U 的边界值，其余按零法向梯度（与 phiHbyA 的边界取法一致）
            std::vector<Tensor> gH;
            if (m.skewed()) {
                // 定通量边界面上 HbyA 的面值取 U_b + rAU_P ∇p_P（HbyA = U + rAU∇p，壁面处并不为零），
                // 这样最小二乘约束 HbyA_b − HbyA_P = U_b − U_P 与速度差一致
                const auto gpc = fvc::grad(p_, gradScheme);
                std::vector<char> bFix(m.nBoundaryFaces(), 0);
                std::vector<Vec3> bVal(m.nBoundaryFaces(), Vec3{});
                for (label p = 0; p < label(m.patches().size()); ++p) {
                    const Patch& pt = m.patches()[p];
                    const bool fixes = U_.bc(p).fixesFlux();
                    for (label i = 0; i < pt.size; ++i) {
                        const label f = pt.start + i, o = m.owner()[f];
                        bFix[f - nI] = fixes;
                        if (fixes) bVal[f - nI] = U_.bValue(f) + rAU_[o] * gpc[o];
                    }
                }
                gH = fvc::gradLS(m, HbyA, bFix, bVal);
            }
            std::vector<scalar> phiHbyA(m.nFaces());
            for (label f = 0; f < nI; ++f) {
                const label o = m.owner()[f], n = m.neighbour()[f];
                const scalar wf = m.w()[f];
                Vec3 Hf = wf * HbyA[o] + (1.0 - wf) * HbyA[n];
                if (!gH.empty()) Hf += gradDot(m.skewVecs()[f], wf * gH[o] + (1.0 - wf) * gH[n]);
                phiHbyA[f] = dot(Hf, m.Sf()[f]);
            }
            for (label p = 0; p < label(m.patches().size()); ++p) {
                const Patch& pt = m.patches()[p];
                const auto& bc = U_.bc(p);
                const bool extrapolate = !bc.fixesFlux();
                for (label i = 0; i < pt.size; ++i) {
                    const label f = pt.start + i;
                    phiHbyA[f] = extrapolate ? dot(HbyA[m.owner()[f]], m.Sf()[f]) : dot(U_.bValue(f), m.Sf()[f]);
                }
            }
            // 体力的面通量（含边界面：定通量边界上经 fixedFluxPressure 给出 ∂p/∂n = f·n）
            if (hasBF)
                for (label f = 0; f < m.nFaces(); ++f) phiHbyA[f] += rAUf[f] * phiF[f];
            if (!steady && ctl.ddtPhiCoeff != 0.0) {
                const auto k = fvm::ddtCoeffs(ddtS, ts_, U_.nOldTimes() >= 2 && !phiOldOld_.empty());
                const auto& U0 = U_.old();
                const auto& U00 = U_.oldOld();
                for (label f = 0; f < nI; ++f) {
                    const scalar wf = m.w()[f];
                    const label o = m.owner()[f], n = m.neighbour()[f];
                    const scalar Uf0 = dot(wf * U0[o] + (1.0 - wf) * U0[n], m.Sf()[f]);
                    const scalar corr0 = phiOld_[f] - Uf0;
                    const scalar coupling = ctl.ddtPhiCoeff < 0
                        ? 1.0 - std::min(std::abs(corr0) / (std::abs(phiOld_[f]) + SMALL), 1.0)
                        : ctl.ddtPhiCoeff;
                    scalar dc = k.c0 * corr0;
                    if (k.c00 != 0.0) {
                        const scalar Uf00 = dot(wf * U00[o] + (1.0 - wf) * U00[n], m.Sf()[f]);
                        dc -= k.c00 * (phiOldOld_[f] - Uf00);
                    }
                    phiHbyA[f] += rAUf[f] * coupling * dc;
                }
            }
            ctx_.phiHbyA = &phiHbyA;
            ctx_.rAUf = &rAUf;
            constrainPressure(phiHbyA, rAUf);
            const auto divPhi = fvc::div(m, phiHbyA);
            const std::vector<scalar> pPrev = steady ? p_.internal() : std::vector<scalar>{};

            for (int no = 0; no <= ctl.nNonOrthCorr; ++no) {
                const bool finalNO = no == ctl.nNonOrthCorr;
                FvMatrix<scalar> pEqn = fvm::laplacian(rAUf, p_, laplacianNonOrthCorr);
                pEqn -= divPhi;
                if (needRef) pEqn.setReference(ctl.pRefCell, ctl.pRefValue);
                auto perf = pEqn.solve((finalCorr && finalNO) ? ctl.pFinalControls : ctl.pControls, ctl.verbose);
                if (outer == 0 && corr == 0 && no == 0 && !perf.empty()) {
                    // 压力恒为常数时（周期域体力驱动）归一化因子趋于零，改用前几步的因子作下限
                    if (pNormSteps_ < 5) {
                        pNormRef_ = std::max(pNormRef_, perf[0].normFactor);
                        ++pNormSteps_;
                    }
                    info.pInitialResidual = perf[0].residualSum / (std::max(perf[0].normFactor, pNormRef_) + 1e-30);
                }
                info.perf.insert(info.perf.end(), perf.begin(), perf.end());
                if (finalNO) {
                    const auto F = pEqn.flux();
                    for (label f = 0; f < m.nFaces(); ++f) phi_[f] = phiHbyA[f] - F[f];
                }
            }
            ctx_.phiHbyA = nullptr;
            ctx_.rAUf = nullptr;
            if (steady && ctl.alphaP < 1.0) {
                for (label c = 0; c < nT; ++c) p_[c] = pPrev[c] + ctl.alphaP * (p_[c] - pPrev[c]);
                p_.correctBoundaryConditions();
            }
            if (hasBF) {
                // U = HbyA + rAU reconstruct(φ_F − ∇p_f·S_f)（压力松弛后的 p）
                const auto pF = fvm::laplacian(rAUf, p_, laplacianNonOrthCorr).flux();
                std::vector<scalar> d(m.nFaces());
                for (label f = 0; f < m.nFaces(); ++f) d[f] = phiF[f] - pF[f] / rAUf[f];
                const auto rec = fvc::reconstruct(m, d);
                for (label c = 0; c < nT; ++c) U_[c] = HbyA[c] + rAU_[c] * rec[c];
            } else {
                const auto gp = fvc::grad(p_, gradScheme);
                for (label c = 0; c < nT; ++c) U_[c] = HbyA[c] - rAU_[c] * gp[c];
            }
            if (m.twoD())
                for (label c = 0; c < nT; ++c) U_[c][m.emptyDir()] = 0.0;
            U_.correctBoundaryConditions();
        }
        applyForcingCorrection(rAU_, fvm::ddtCoeffs(ddtS, ts_, U_.nOldTimes() >= 2).c, aU);
    }

    // ---------------------------------------------------------- 统计
    const auto dv = fvc::div(m, phi_);
    scalar ce = 0;
    for (label c = 0; c < nC; ++c) ce = std::max(ce, std::abs(dv[c]));
    info.continuityError = par::allMax(ce);
    info.maxCo = steady ? 0.0 : courantNumber(dt);
    return info;
}


scalar IncompressibleFlow::setBernoulliPressure() {
    const Mesh& m = *mesh_;
    par::SumAcc s[2];
    for (label pi = 0; pi < label(m.patches().size()); ++pi) {
        const Patch& pt = m.patches()[pi];
        if (!p_.hasBC(pi) || !p_.bc(pi).fixesValue()) continue;
        for (label f = pt.start; f < pt.end(); ++f) {
            s[0].add((p_.bValue(f) + 0.5 * magSqr(U_.bValue(f))) * m.magSf()[f]);
            s[1].add(m.magSf()[f]);
        }
    }
    double r[2];
    par::SumAcc::allReduce(s, 2, r);
    const scalar H = r[1] > 0 ? r[0] / r[1] : 0.0;
    for (label c = 0; c < m.nTotalCells(); ++c) p_.internal()[c] = H - 0.5 * magSqr(U_[c]);
    p_.correctBoundaryConditions();
    return H;
}

} // namespace cfd
