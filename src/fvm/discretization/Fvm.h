#pragma once
// 隐式算子（fvm）：返回 FvMatrix，表示算子作用于 ψ 的体积分。

#include "fvm/discretization/Fvc.h"
#include "fvm/linalg/FvMatrix.h"

#include <algorithm>
#include <cmath>

namespace cfd::fvm {

// ------------------------------------------------------------------ ddt
// Euler:    V(ψ − ψ⁰)/Δt
// backward: V(c ψ − c0 ψ⁰ + c00 ψ⁰⁰)/Δt（变步长 BDF2），首步退化为 Euler
struct DdtCoeffs {
    scalar c = 0, c0 = 0, c00 = 0;  // 已除以 Δt
};
inline DdtCoeffs ddtCoeffs(DdtScheme s, const TimeState& ts, bool haveOldOld) {
    DdtCoeffs k;
    if (s == DdtScheme::Steady) return k;
    const scalar rDt = 1.0 / ts.dt;
    if (s == DdtScheme::Backward && haveOldOld && ts.dt0 > 0) {
        const scalar dt = ts.dt, dt0 = ts.dt0;
        const scalar ct = 1.0 + dt / (dt + dt0);
        const scalar ct00 = dt * dt / (dt0 * (dt + dt0));
        k.c = ct * rDt;
        k.c0 = (ct + ct00) * rDt;
        k.c00 = ct00 * rDt;
    } else {
        k.c = k.c0 = rDt;
    }
    return k;
}

template <class T> FvMatrix<T> ddt(VolField<T>& psi, DdtScheme s, const TimeState& ts) {
    FvMatrix<T> M(psi);
    if (s == DdtScheme::Steady) return M;
    const DdtCoeffs k = ddtCoeffs(s, ts, psi.nOldTimes() >= 2);
    const Mesh& m = psi.mesh();
    const auto& o = psi.old();
    const auto& oo = psi.oldOld();
    for (label c = 0; c < m.nCells(); ++c) {
        const scalar V = m.V()[c];
        M.diag()[c] += k.c * V;
        M.source()[c] += V * (k.c0 * o[c] - k.c00 * oo[c]);
    }
    return M;
}

// ------------------------------------------------------------------ Sp / Su
// Sp: 隐式源项 coeff·ψ（单位体积），Su: 显式源项
template <class T> FvMatrix<T> Sp(const std::vector<scalar>& coeff, VolField<T>& psi) {
    FvMatrix<T> M(psi);
    for (label c = 0; c < psi.mesh().nCells(); ++c) M.diag()[c] += coeff[c] * psi.mesh().V()[c];
    return M;
}
template <class T> FvMatrix<T> Sp(scalar coeff, VolField<T>& psi) {
    FvMatrix<T> M(psi);
    for (label c = 0; c < psi.mesh().nCells(); ++c) M.diag()[c] += coeff * psi.mesh().V()[c];
    return M;
}

// ------------------------------------------------------------------ 对流 ∇·(φ ψ)
namespace detail {
inline scalar limiterR(scalar F, scalar phiP, scalar phiN, const Vec3& gP, const Vec3& gN, const Vec3& d) {
    const scalar gradf = phiN - phiP;
    const scalar gradcf = F > 0 ? dot(d, gP) : dot(d, gN);
    if (std::abs(gradcf) >= 1000.0 * std::abs(gradf)) return 2.0 * 1000.0 * (gradcf >= 0 ? 1 : -1) * (gradf >= 0 ? 1 : -1) - 1.0;
    return 2.0 * (gradcf / gradf) - 1.0;
}
inline scalar limiterR(scalar F, const Vec3& phiP, const Vec3& phiN, const Tensor& gP, const Tensor& gN, const Vec3& d) {
    const Vec3 gradf = phiN - phiP;
    const Vec3 gradcf = F > 0 ? dot(d, gP) : dot(d, gN);
    const scalar a = dot(gradf, gradcf), b = dot(gradf, gradf);
    if (std::abs(a) >= 1000.0 * b) return 2.0 * 1000.0 * (a >= 0 ? 1 : -1) - 1.0;
    return 2.0 * (a / b) - 1.0;
}
inline scalar limiter(const ConvectionScheme& s, scalar r) {
    using Ty = ConvectionScheme::Type;
    switch (s.type) {
    case Ty::VanLeer: return (r + std::abs(r)) / (1.0 + std::abs(r));
    case Ty::Minmod: return std::max(std::min(r, 1.0), 0.0);
    case Ty::MUSCL: return std::max(std::min(std::min(2.0 * r, 0.5 * r + 0.5), 2.0), 0.0);
    case Ty::LimitedLinear: return std::max(std::min(2.0 / std::max(s.k, SMALL) * r, 1.0), 0.0);
    default: return 1.0;
    }
}
} // namespace detail

// F: 面体积通量（owner→neighbour 为正）
template <class T> FvMatrix<T> div(const std::vector<scalar>& F, VolField<T>& psi, const ConvectionScheme& s) {
    using Ty = ConvectionScheme::Type;
    const Mesh& m = psi.mesh();
    FvMatrix<T> M(psi);
    const label nI = m.nInternalFaces();
    const auto& own = m.owner();
    const auto& nei = m.neighbour();
    std::vector<typename Traits<T>::Grad> g;
    if (s.needsGradient()) g = fvc::grad(psi);

    auto& lo = M.lower();
    auto& up = M.upper();
    for (label f = 0; f < nI; ++f) {
        const scalar Ff = F[f];
        const scalar wUp = Ff >= 0 ? 1.0 : 0.0;
        scalar wf;
        switch (s.type) {
        case Ty::Upwind:
        case Ty::LinearUpwind: wf = wUp; break;
        case Ty::Linear: wf = m.w()[f]; break;
        case Ty::LUST: wf = 0.75 * m.w()[f] + 0.25 * wUp; break;
        default: {
            const scalar r = detail::limiterR(Ff, psi[own[f]], psi[nei[f]], g[own[f]], g[nei[f]], m.d()[f]);
            const scalar lim = detail::limiter(s, r);
            wf = lim * m.w()[f] + (1.0 - lim) * wUp;
        }
        }
        lo[f] = -wf * Ff;
        up[f] = lo[f] + Ff;
        M.diag()[own[f]] -= lo[f];
        M.diag()[nei[f]] -= up[f];
        // 梯度修正（延迟修正显式部分）
        if (s.type == Ty::LinearUpwind || s.type == Ty::LUST) {
            const scalar fac = s.type == Ty::LUST ? 0.25 : 1.0;
            T corr;
            if (Ff >= 0)
                corr = gradDot(m.Cf()[f] - m.C()[own[f]], g[own[f]]);
            else
                corr = gradDot(m.Cf()[f] - (m.C()[own[f]] + m.d()[f]), g[nei[f]]);
            const T c = (fac * Ff) * corr;
            M.source()[own[f]] -= c;
            M.source()[nei[f]] += c;
        }
    }
    // 边界：F (vic ψ_P + vbc)
    for (label p = 0; p < label(m.patches().size()); ++p) {
        const Patch& pt = m.patches()[p];
        const auto& bc = psi.bc(p);
        for (label i = 0; i < pt.size; ++i) {
            const label f = pt.start + i, bf = f - nI;
            M.internalCoeffs()[bf] += F[f] * bc.valueInternalCoeff(psi, i);
            M.boundarySource()[bf] -= F[f] * bc.valueBoundaryCoeff(psi, i);
        }
    }
    return M;
}

// ------------------------------------------------------------------ 扩散 ∇·(Γ∇ψ)
// gamma: 面上扩散系数（长度 nFaces）。nonOrthCorr=true 时加显式非正交修正。
template <class T>
FvMatrix<T> laplacian(const std::vector<scalar>& gamma, VolField<T>& psi, bool nonOrthCorr = true) {
    const Mesh& m = psi.mesh();
    FvMatrix<T> M(psi);
    const label nI = m.nInternalFaces();
    const auto& own = m.owner();
    const auto& nei = m.neighbour();
    auto& up = M.upper();
    for (label f = 0; f < nI; ++f) {
        const scalar c = gamma[f] * m.magSf()[f] * m.nonOrthDeltaCoeffs()[f];
        up[f] = c;
        M.diag()[own[f]] -= c;
        M.diag()[nei[f]] -= c;
    }
    if (nonOrthCorr) {
        bool any = false;
        for (label f = 0; f < nI && !any; ++f) any = magSqr(m.corrVecs()[f]) > 1e-20;
        any = par::allMax(any ? 1.0 : 0.0) > 0.5;
        if (any) {
            const auto g = fvc::grad(psi);
            auto& ffc = M.faceFluxCorrection();
            for (label f = 0; f < nI; ++f) {
                const auto gf = m.w()[f] * g[own[f]] + (1.0 - m.w()[f]) * g[nei[f]];
                const T fc = (gamma[f] * m.magSf()[f]) * gradDot(m.corrVecs()[f], gf);
                ffc[f] = fc;
                M.source()[own[f]] -= fc;
                M.source()[nei[f]] += fc;
            }
        }
    }
    for (label p = 0; p < label(m.patches().size()); ++p) {
        const Patch& pt = m.patches()[p];
        const auto& bc = psi.bc(p);
        for (label i = 0; i < pt.size; ++i) {
            const label f = pt.start + i, bf = f - nI;
            const scalar gs = gamma[f] * m.magSf()[f];
            M.internalCoeffs()[bf] += gs * bc.gradientInternalCoeff(psi, i);
            M.boundarySource()[bf] -= gs * bc.gradientBoundaryCoeff(psi, i);
        }
    }
    return M;
}

template <class T> FvMatrix<T> laplacian(scalar gamma, VolField<T>& psi, bool nonOrthCorr = true) {
    return laplacian(std::vector<scalar>(psi.mesh().nFaces(), gamma), psi, nonOrthCorr);
}

} // namespace cfd::fvm
