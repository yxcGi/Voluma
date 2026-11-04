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

// 由面通量 φ_f = v_f·S_f 重构单元矢量（OpenFOAM fvc::reconstruct）：v_P = [Σ_f S_f S_f/|S_f|]⁻¹ Σ_f S_f φ_f/|S_f|，
// 两项都与面的朝向无关；用于把面上平衡的体力与压力梯度一致地回到单元（静水平衡时结果精确为零）。幽灵单元经 halo 交换。
inline std::vector<Vec3> reconstruct(const Mesh& m, const std::vector<scalar>& ssf) {
    const label nT = m.nTotalCells(), nC = m.nCells(), nI = m.nInternalFaces();
    std::vector<Tensor> A(nT, Tensor::zero());
    std::vector<Vec3> b(nT, Vec3{0, 0, 0});
    auto acc = [&](label c, label f) {
        const Vec3& S = m.Sf()[f];
        const scalar ms = m.magSf()[f];
        A[c] += outer(S, S) / ms;
        b[c] += (ssf[f] / ms) * S;
    };
    for (label f = 0; f < nI; ++f) {
        acc(m.owner()[f], f);
        acc(m.neighbour()[f], f);
    }
    for (label f = nI; f < m.nFaces(); ++f) acc(m.owner()[f], f);
    std::vector<Vec3> v(nT, Vec3{0, 0, 0});
    // 二维网格没有空方向的面，A 在该方向奇异：补一个单位对角元（该方向 b = 0，结果为零）
    const int ed = m.twoD() ? m.emptyDir() : -1;
    for (label c = 0; c < nC; ++c) {
        if (ed >= 0) (&A[c].xx)[ed * 3 + ed] += 1.0;
        const Tensor& t = A[c];
        // 3×3 对称矩阵求逆（伴随矩阵）
        const scalar c00 = t.yy * t.zz - t.yz * t.zy, c01 = t.yz * t.zx - t.yx * t.zz, c02 = t.yx * t.zy - t.yy * t.zx;
        const scalar det = t.xx * c00 + t.xy * c01 + t.xz * c02;
        if (std::abs(det) < 1e-300) continue;
        const scalar c10 = t.xz * t.zy - t.xy * t.zz, c11 = t.xx * t.zz - t.xz * t.zx, c12 = t.xy * t.zx - t.xx * t.zy;
        const scalar c20 = t.xy * t.yz - t.xz * t.yy, c21 = t.xz * t.yx - t.xx * t.yz, c22 = t.xx * t.yy - t.xy * t.yx;
        const Vec3& r = b[c];
        v[c] = Vec3{c00 * r.x + c10 * r.y + c20 * r.z, c01 * r.x + c11 * r.y + c21 * r.z,
                    c02 * r.x + c12 * r.y + c22 * r.z} / det;
    }
    m.halo().exchange(v);
    return v;
}

// 偏斜修正的面插值：φ_f = w φ_P + (1−w) φ_N + (x_f − x'_f)·∇φ_f（∇φ_f 为单元梯度的线性插值）。
// 偏斜网格（四面体、三角形）上线性插值点 x'_f 不在面心，不修正时面值有一阶误差；无偏斜网格原样返回。
template <class T>
std::vector<T> interpolate(const Mesh& m, const std::vector<T>& c, const std::vector<typename Traits<T>::Grad>& g) {
    auto F = interpolate(m, c);
    if (!m.skewed()) return F;
    const auto& w = m.w();
    const auto& sk = m.skewVecs();
    for (label i = 0; i < m.nInternalFaces(); ++i)
        F[i] += gradDot(sk[i], w[i] * g[m.owner()[i]] + (1.0 - w[i]) * g[m.neighbour()[i]]);
    return F;
}
template <class T>
std::vector<T> interpolate(const VolField<T>& f, const std::vector<typename Traits<T>::Grad>& g) {
    auto F = interpolate(f);
    const Mesh& m = f.mesh();
    if (!m.skewed()) return F;
    const auto& w = m.w();
    const auto& sk = m.skewVecs();
    for (label i = 0; i < m.nInternalFaces(); ++i)
        F[i] += gradDot(sk[i], w[i] * g[m.owner()[i]] + (1.0 - w[i]) * g[m.neighbour()[i]]);
    return F;
}

// 体积通量 φ_f = U_f · S_f（给出 ∇U 时面速度带偏斜修正）
inline std::vector<scalar> flux(const VolField<Vec3>& U, const std::vector<Tensor>& gradU) {
    const Mesh& m = U.mesh();
    auto Uf = interpolate(U, gradU);
    std::vector<scalar> phi(m.nFaces());
    for (label f = 0; f < m.nFaces(); ++f) phi[f] = dot(Uf[f], m.Sf()[f]);
    return phi;
}
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
inline scalar cmptMulT(scalar a, scalar b) { return a * b; }
inline Vec3 cmptMulT(const Vec3& a, const Vec3& b) { return cmptMultiply(a, b); }
inline Vec3 lsApply(const Tensor& G, const Vec3& d, scalar dphi) { return dot(G, d) * dphi; }
inline Tensor lsApply(const Tensor& G, const Vec3& d, const Vec3& dphi) { return outer(dot(G, d), dphi); }

// 加权最小二乘梯度核心：f 为单元值（含幽灵，已交换），de / dphi 为各边界面的约束向量与对应差值
// （定值面：d 与 φ_b − φ_P；梯度面：n (n·d) 与 g_b (n·d)），权重 1/|d|²。结果未做二维置零与幽灵交换。
template <class T>
std::vector<typename Traits<T>::Grad> leastSquaresGrad(const Mesh& m, const T* f, const std::vector<Vec3>& de,
                                                       const std::vector<T>& dphi) {
    using G = typename Traits<T>::Grad;
    const label nI = m.nInternalFaces(), nC = m.nCells();
    const auto& own = m.owner();
    const auto& nei = m.neighbour();
    std::vector<G> g(m.nTotalCells(), G{});
    std::vector<Tensor> dd(m.nTotalCells(), Tensor{});
    for (label i = 0; i < m.nFaces(); ++i) {
        const Vec3& d = i < nI ? m.d()[i] : de[i - nI];
        const Tensor t = outer(d, d) * (1.0 / magSqr(d));
        dd[own[i]] += t;
        if (i < nI) dd[nei[i]] += t;
    }
    std::vector<Tensor> Gi(nC);
    for (label c = 0; c < nC; ++c) {
        Tensor a = dd[c];
        if (m.twoD()) a(m.emptyDir(), m.emptyDir()) += 1.0;
        Gi[c] = inv(a);
    }
    for (label i = 0; i < nI; ++i) {
        const Vec3& d = m.d()[i];
        const scalar wgt = 1.0 / magSqr(d);
        const T dp = f[nei[i]] - f[own[i]];
        if (own[i] < nC) g[own[i]] += lsApply(Gi[own[i]], d * wgt, dp);
        if (nei[i] < nC) g[nei[i]] += lsApply(Gi[nei[i]], d * wgt, dp);
    }
    for (label i = nI; i < m.nFaces(); ++i) {
        const Vec3& d = de[i - nI];
        g[own[i]] += lsApply(Gi[own[i]], d * (1.0 / magSqr(d)), dphi[i - nI]);
    }
    return g;
}
// 二维 empty 方向分量置零 + 幽灵交换
template <class G> void finishGrad(const Mesh& m, std::vector<G>& g) {
    if (m.twoD()) {
        const int e = m.emptyDir();
        for (label c = 0; c < m.nCells(); ++c) {
            if constexpr (std::is_same_v<G, Vec3>) {
                g[c][e] = 0;
            } else {
                for (int j = 0; j < 3; ++j) g[c](e, j) = 0;
            }
        }
    }
    m.halo().exchange(g);
}
} // namespace detail

// 单元数组的最小二乘梯度（无边界条件对象）：bFixesValue[bf] 为真的边界面用面值 bValue[bf]，否则按零法向梯度约束
template <class T>
std::vector<typename Traits<T>::Grad> gradLS(const Mesh& m, const std::vector<T>& c, const std::vector<char>& bFixesValue,
                                             const std::vector<T>& bValue) {
    const label nI = m.nInternalFaces();
    std::vector<Vec3> de(m.nFaces() - nI);
    std::vector<T> dphi(m.nFaces() - nI, Traits<T>::zero());
    for (label i = nI; i < m.nFaces(); ++i) {
        const label bf = i - nI;
        const Vec3& d = m.d()[i];
        if (bFixesValue[bf]) {
            de[bf] = d;
            dphi[bf] = bValue[bf] - c[m.owner()[i]];
        } else {
            const Vec3 n = m.Sf()[i] * (1.0 / m.magSf()[i]);
            de[bf] = n * std::max(dot(n, d), 0.05 * mag(d));
        }
    }
    auto g = detail::leastSquaresGrad(m, c.data(), de, dphi);
    detail::finishGrad(m, g);
    return g;
}

// 梯度：Gauss 线性 或 加权最小二乘。结果含幽灵单元（已交换）。不指定格式时用网格的默认格式（Mesh::gradScheme）。
template <class T>
std::vector<typename Traits<T>::Grad> grad(const VolField<T>& f, GradScheme scheme) {
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
    } else if (scheme == GradScheme::GaussSkew) {
        // 守恒的 Gauss 梯度 + 偏斜修正：面值 φ_f = 线性插值 + (x_f − x'_f)·∇φ_f（∇φ 取最小二乘），
        // 梯度类边界面取 φ_P + ∇φ_t·d_t + g_b (n·d)（切向按单元梯度外推，法向按边界条件给的梯度）。
        // 面值由相邻两单元共享，体积和严格等于边界面合力 Σ φ_f S_f；最小二乘梯度不守恒，在四面体网格上
        // 其体积和比边界合力小约 2 %，用于动量方程的压力梯度时管流压降偏高同样的比例（docs/validation.md 六）。
        const auto gls = grad(f, GradScheme::LeastSquares);
        const auto& w = m.w();
        const auto& sk = m.skewVecs();
        for (label i = 0; i < nI; ++i) {
            const T phif = w[i] * f[own[i]] + (1.0 - w[i]) * f[nei[i]] +
                           gradDot(sk[i], w[i] * gls[own[i]] + (1.0 - w[i]) * gls[nei[i]]);
            const G c = outerSf(m.Sf()[i], phif);
            g[own[i]] += c;
            g[nei[i]] -= c;
        }
        for (label p = 0; p < label(m.patches().size()); ++p) {
            const Patch& pt = m.patches()[p];
            const auto& bc = f.bc(p);
            const bool value = bc.fixesValue();
            for (label k = 0; k < pt.size; ++k) {
                const label i = pt.start + k;
                T phif;
                if (value) {
                    phif = f.bValue(i);
                } else {
                    const Vec3& d = m.d()[i];
                    const Vec3 n = m.Sf()[i] * (1.0 / m.magSf()[i]);
                    const scalar nd = dot(n, d);
                    const Vec3 dt = d - n * nd;
                    const T gb = detail::cmptMulT(bc.gradientInternalCoeff(f, k), f[own[i]]) + bc.gradientBoundaryCoeff(f, k);
                    phif = f[own[i]] + gradDot(dt, gls[own[i]]) + gb * nd;
                }
                g[own[i]] += outerSf(m.Sf()[i], phif);
            }
        }
        for (label c = 0; c < nC; ++c) g[c] = g[c] * (1.0 / m.V()[c]);
    } else {
        // 加权最小二乘：Σ w d⊗d，权重 1/|d|²。
        // 边界面：定值类条件用面心值与 d（φ_b − φ_P = ∇φ·d）；梯度类/混合类条件（面值由单元值外推，不是面心处的真实值）
        // 改用法向导数约束 ∇φ·n = g_b，写成 ∇φ·(n (n·d)) = g_b (n·d)。否则偏斜面（d 有切向分量）上
        // zeroGradient 会错误地约束沿 d 的差分为零，使壁面单元梯度有不随加密减小的误差（四面体网格尤甚）。
        std::vector<Vec3> de(m.nFaces() - nI);
        std::vector<T> dphi(m.nFaces() - nI);
        for (label p = 0; p < label(m.patches().size()); ++p) {
            const Patch& pt = m.patches()[p];
            const auto& bc = f.bc(p);
            const bool value = bc.fixesValue();
            for (label k = 0; k < pt.size; ++k) {
                const label i = pt.start + k, bf = i - nI;
                const Vec3& d = m.d()[i];
                if (value) {
                    de[bf] = d;
                    dphi[bf] = f.bValue(i) - f[own[i]];
                } else {
                    const Vec3 n = m.Sf()[i] * (1.0 / m.magSf()[i]);
                    const scalar nd = std::max(dot(n, d), 0.05 * mag(d));
                    de[bf] = n * nd;
                    dphi[bf] = (detail::cmptMulT(bc.gradientInternalCoeff(f, k), f[own[i]]) + bc.gradientBoundaryCoeff(f, k)) * nd;
                }
            }
        }
        g = detail::leastSquaresGrad(m, f.internal().data(), de, dphi);
    }
    detail::finishGrad(m, g);
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

template <class T> std::vector<typename Traits<T>::Grad> grad(const VolField<T>& f) {
    return grad(f, f.mesh().gradScheme());
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
