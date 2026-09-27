#pragma once
// 显式算子（fvc）：插值、梯度、散度、面法向梯度等。
// 返回的单元量只在自有单元上计算，之后统一 halo 交换，幽灵单元值与其所属进程一致。

#include "fvm/discretization/Schemes.h"
#include "fvm/field/VolField.h"

#include <vector>

namespace cfd::fvc {

// 面插值（线性），边界面取边界值
template <class T> std::vector<T> interpolate(const VolField<T>& f) {
    const Mesh& m = f.mesh();
    const label nI = m.nInternalFaces();
    std::vector<T> F(m.nFaces());
    const auto& w = m.w();
    for (label i = 0; i < nI; ++i) F[i] = w[i] * f[m.owner()[i]] + (1.0 - w[i]) * f[m.neighbour()[i]];
    for (label i = nI; i < m.nFaces(); ++i) F[i] = f.bValue(i);
    return F;
}

// 单元数组的面插值（边界面取所在单元值）
template <class T> std::vector<T> interpolate(const Mesh& m, const std::vector<T>& c) {
    const label nI = m.nInternalFaces();
    std::vector<T> F(m.nFaces());
    const auto& w = m.w();
    for (label i = 0; i < nI; ++i) F[i] = w[i] * c[m.owner()[i]] + (1.0 - w[i]) * c[m.neighbour()[i]];
    for (label i = nI; i < m.nFaces(); ++i) F[i] = c[m.owner()[i]];
    return F;
}

// 体积通量 φ_f = U_f · S_f
inline std::vector<scalar> flux(const VolField<Vec3>& U) {
    const Mesh& m = U.mesh();
    auto Uf = interpolate(U);
    std::vector<scalar> phi(m.nFaces());
    for (label f = 0; f < m.nFaces(); ++f) phi[f] = dot(Uf[f], m.Sf()[f]);
    return phi;
}

namespace detail {
// 3×3 对称矩阵求逆（最小二乘梯度用）
inline Tensor inv(const Tensor& a) {
    const scalar det = a.xx * (a.yy * a.zz - a.yz * a.zy) - a.xy * (a.yx * a.zz - a.yz * a.zx) +
                       a.xz * (a.yx * a.zy - a.yy * a.zx);
    Tensor r;
    r.xx = (a.yy * a.zz - a.yz * a.zy) / det;
    r.xy = (a.xz * a.zy - a.xy * a.zz) / det;
    r.xz = (a.xy * a.yz - a.xz * a.yy) / det;
    r.yx = (a.yz * a.zx - a.yx * a.zz) / det;
    r.yy = (a.xx * a.zz - a.xz * a.zx) / det;
    r.yz = (a.xz * a.yx - a.xx * a.yz) / det;
    r.zx = (a.yx * a.zy - a.yy * a.zx) / det;
    r.zy = (a.xy * a.zx - a.xx * a.zy) / det;
    r.zz = (a.xx * a.yy - a.xy * a.yx) / det;
    return r;
}
inline Vec3 lsApply(const Tensor& G, const Vec3& d, scalar dphi) { return dot(G, d) * dphi; }
inline Tensor lsApply(const Tensor& G, const Vec3& d, const Vec3& dphi) { return outer(dot(G, d), dphi); }
} // namespace detail

// 梯度：Gauss 线性 或 加权最小二乘。结果含幽灵单元（已交换）。
template <class T>
std::vector<typename Traits<T>::Grad> grad(const VolField<T>& f, GradScheme scheme = GradScheme::GaussLinear) {
    using G = typename Traits<T>::Grad;
    const Mesh& m = f.mesh();
    const label nI = m.nInternalFaces(), nC = m.nCells();
    std::vector<G> g(m.nTotalCells(), G{});
    const auto& own = m.owner();
    const auto& nei = m.neighbour();
    if (scheme == GradScheme::GaussLinear) {
        const auto& w = m.w();
        for (label i = 0; i < nI; ++i) {
            const T phif = w[i] * f[own[i]] + (1.0 - w[i]) * f[nei[i]];
            const G c = outerSf(m.Sf()[i], phif);
            g[own[i]] += c;
            g[nei[i]] -= c;
        }
        for (label i = nI; i < m.nFaces(); ++i) g[own[i]] += outerSf(m.Sf()[i], f.bValue(i));
        for (label c = 0; c < nC; ++c) g[c] = g[c] * (1.0 / m.V()[c]);
    } else {
        // 加权最小二乘：Σ w d⊗d，权重 1/|d|²
        std::vector<Tensor> dd(m.nTotalCells(), Tensor{});
        for (label i = 0; i < m.nFaces(); ++i) {
            const Vec3& d = m.d()[i];
            const Tensor t = outer(d, d) * (1.0 / magSqr(d));
            dd[own[i]] += t;
            if (i < nI) dd[nei[i]] += t;
        }
        std::vector<Tensor> Gi(nC);
        for (label c = 0; c < nC; ++c) {
            Tensor a = dd[c];
            if (m.twoD()) a(m.emptyDir(), m.emptyDir()) += 1.0;
            Gi[c] = detail::inv(a);
        }
        for (label i = 0; i < nI; ++i) {
            const Vec3& d = m.d()[i];
            const scalar wgt = 1.0 / magSqr(d);
            const T dphi = f[nei[i]] - f[own[i]];
            if (own[i] < nC) g[own[i]] += detail::lsApply(Gi[own[i]], d * wgt, dphi);
            if (nei[i] < nC) g[nei[i]] += detail::lsApply(Gi[nei[i]], d * wgt, dphi);
        }
        for (label i = nI; i < m.nFaces(); ++i) {
            const Vec3& d = m.d()[i];
            const T dphi = f.bValue(i) - f[own[i]];
            g[own[i]] += detail::lsApply(Gi[own[i]], d * (1.0 / magSqr(d)), dphi);
        }
    }
    if (m.twoD()) {
        // 二维：empty 方向的梯度分量置零
        const int e = m.emptyDir();
        for (label c = 0; c < nC; ++c) {
            if constexpr (std::is_same_v<T, scalar>) {
                g[c][e] = 0;
            } else {
                for (int j = 0; j < 3; ++j) g[c](e, j) = 0;
            }
        }
    }
    m.halo().exchange(g);
    return g;
}

// 散度 ∇·F（F 为面通量），单位体积
inline std::vector<scalar> div(const Mesh& m, const std::vector<scalar>& F) {
    std::vector<scalar> r(m.nTotalCells(), 0.0);
    const label nI = m.nInternalFaces();
    for (label i = 0; i < nI; ++i) {
        r[m.owner()[i]] += F[i];
        r[m.neighbour()[i]] -= F[i];
    }
    for (label i = nI; i < m.nFaces(); ++i) r[m.owner()[i]] += F[i];
    for (label c = 0; c < m.nCells(); ++c) r[c] /= m.V()[c];
    return r;
}

// 面上向量通量的散度 ∇·(F_f)，F_f 为面上的 T 值（已乘面积）
template <class T> std::vector<T> surfaceIntegrate(const Mesh& m, const std::vector<T>& F) {
    std::vector<T> r(m.nTotalCells(), Traits<T>::zero());
    const label nI = m.nInternalFaces();
    for (label i = 0; i < nI; ++i) {
        r[m.owner()[i]] += F[i];
        r[m.neighbour()[i]] -= F[i];
    }
    for (label i = nI; i < m.nFaces(); ++i) r[m.owner()[i]] += F[i];
    for (label c = 0; c < m.nCells(); ++c) r[c] = r[c] * (1.0 / m.V()[c]);
    return r;
}

// 面法向梯度（内部面含非正交修正）
template <class T>
std::vector<T> snGrad(const VolField<T>& f, const std::vector<typename Traits<T>::Grad>* gradF = nullptr) {
    const Mesh& m = f.mesh();
    const label nI = m.nInternalFaces();
    std::vector<T> s(m.nFaces());
    for (label i = 0; i < nI; ++i) {
        s[i] = m.nonOrthDeltaCoeffs()[i] * (f[m.neighbour()[i]] - f[m.owner()[i]]);
        if (gradF) {
            const auto gf = m.w()[i] * (*gradF)[m.owner()[i]] + (1.0 - m.w()[i]) * (*gradF)[m.neighbour()[i]];
            s[i] += gradDot(m.corrVecs()[i], gf);
        }
    }
    for (label i = nI; i < m.nFaces(); ++i) s[i] = m.nonOrthDeltaCoeffs()[i] * (f.bValue(i) - f[m.owner()[i]]);
    return s;
}

// 全局体积平均
template <class T> T domainAverage(const VolField<T>& f) {
    const Mesh& m = f.mesh();
    constexpr int nc = Traits<T>::nComponents;
    par::SumAcc s[nc];
    for (label c = 0; c < m.nCells(); ++c)
        for (int k = 0; k < nc; ++k) s[k].add(Traits<T>::component(f[c], k) * m.V()[c]);
    double out[nc];
    par::SumAcc::allReduce(s, nc, out);
    T r = Traits<T>::zero();
    for (int k = 0; k < nc; ++k) Traits<T>::setComponent(r, k, out[k] / m.totalVolume());
    return r;
}

} // namespace cfd::fvc
