#pragma once
// 有限体积矩阵：表示离散方程的残量形式 R(ψ) = A ψ − b（单元体积分形式）。
//
//   diag[c]、upper[f]=A(owner,neighbour)、lower[f]=A(neighbour,owner)：各分量共用
//   source[c]                 : b（T 值）
//   internalCoeffs[bf]        : 边界面对其所在单元对角元的逐分量贡献
//   boundarySource[bf]        : 边界面对其所在单元 b 的贡献
//   faceFluxCorrection[f]     : 显式非正交修正通量（flux() 使用）
//
// 运算约定：M += N / M -= N 合并矩阵；M -= su 表示 R − su·V（su 为单位体积源项）。
// solve() 求 R(ψ)=0。

#include "fvm/field/VolField.h"
#include "fvm/linalg/LinearSolver.h"

#include <iostream>
#include <stdexcept>

namespace cfd {

template <class T> class FvMatrix {
public:
    explicit FvMatrix(VolField<T>& psi)
        : psi_(&psi), diag_(psi.mesh().nTotalCells(), 0.0), upper_(psi.mesh().nInternalFaces(), 0.0),
          source_(psi.mesh().nTotalCells(), Traits<T>::zero()),
          internalCoeffs_(psi.mesh().nBoundaryFaces(), Traits<T>::zero()),
          boundarySource_(psi.mesh().nBoundaryFaces(), Traits<T>::zero()) {}

    VolField<T>& psi() const { return *psi_; }
    const Mesh& mesh() const { return psi_->mesh(); }

    std::vector<scalar>& diag() { return diag_; }
    const std::vector<scalar>& diag() const { return diag_; }
    std::vector<scalar>& upper() { return upper_; }
    const std::vector<scalar>& upper() const { return upper_; }
    // 首次访问非对称部分时由 upper 复制出 lower
    std::vector<scalar>& lower() {
        if (!hasLower_) {
            lower_ = upper_;
            hasLower_ = true;
        }
        return lower_;
    }
    const scalar* lowerPtr() const { return hasLower_ ? lower_.data() : upper_.data(); }
    bool symmetric() const { return !hasLower_; }
    std::vector<T>& source() { return source_; }
    const std::vector<T>& source() const { return source_; }
    std::vector<T>& internalCoeffs() { return internalCoeffs_; }
    std::vector<T>& boundarySource() { return boundarySource_; }
    std::vector<T>& faceFluxCorrection() {
        if (ffc_.empty()) ffc_.assign(mesh().nFaces(), Traits<T>::zero());
        return ffc_;
    }
    bool hasFaceFluxCorrection() const { return !ffc_.empty(); }

    // ---------------------------------------------------------- 代数运算
    FvMatrix& operator+=(const FvMatrix& o) {
        check(o);
        for (std::size_t i = 0; i < diag_.size(); ++i) diag_[i] += o.diag_[i];
        addOffDiag(o, 1.0);
        for (std::size_t i = 0; i < source_.size(); ++i) source_[i] += o.source_[i];
        for (std::size_t i = 0; i < internalCoeffs_.size(); ++i) {
            internalCoeffs_[i] += o.internalCoeffs_[i];
            boundarySource_[i] += o.boundarySource_[i];
        }
        addFfc(o, 1.0);
        return *this;
    }
    FvMatrix& operator-=(const FvMatrix& o) {
        check(o);
        for (std::size_t i = 0; i < diag_.size(); ++i) diag_[i] -= o.diag_[i];
        addOffDiag(o, -1.0);
        for (std::size_t i = 0; i < source_.size(); ++i) source_[i] -= o.source_[i];
        for (std::size_t i = 0; i < internalCoeffs_.size(); ++i) {
            internalCoeffs_[i] -= o.internalCoeffs_[i];
            boundarySource_[i] -= o.boundarySource_[i];
        }
        addFfc(o, -1.0);
        return *this;
    }
    FvMatrix& operator*=(scalar s) {
        for (auto& v : diag_) v *= s;
        for (auto& v : upper_) v *= s;
        if (hasLower_)
            for (auto& v : lower_) v *= s;
        for (auto& v : source_) v *= s;
        for (auto& v : internalCoeffs_) v *= s;
        for (auto& v : boundarySource_) v *= s;
        for (auto& v : ffc_) v *= s;
        return *this;
    }
    // R − su·V（su：单位体积源项，长度 ≥ nCells）
    FvMatrix& operator-=(const std::vector<T>& su) {
        const auto& V = mesh().V();
        for (label c = 0; c < mesh().nCells(); ++c) source_[c] += su[c] * V[c];
        return *this;
    }
    FvMatrix& operator+=(const std::vector<T>& su) {
        const auto& V = mesh().V();
        for (label c = 0; c < mesh().nCells(); ++c) source_[c] -= su[c] * V[c];
        return *this;
    }
    // 均匀源项
    FvMatrix& operator-=(const T& su) {
        const auto& V = mesh().V();
        for (label c = 0; c < mesh().nCells(); ++c) source_[c] += su * V[c];
        return *this;
    }

    // ---------------------------------------------------------- 欠松弛（隐式）
    void relax(scalar alpha) {
        if (alpha >= 1.0) return;
        const Mesh& m = mesh();
        const label n = m.nCells();
        std::vector<scalar> sumOff(m.nTotalCells(), 0.0), bD(m.nTotalCells(), 0.0);
        const scalar* lo = lowerPtr();
        for (label f = 0; f < m.nInternalFaces(); ++f) {
            sumOff[m.owner()[f]] += std::abs(upper_[f]);
            sumOff[m.neighbour()[f]] += std::abs(lo[f]);
        }
        const label nI = m.nInternalFaces();
        for (label bf = 0; bf < m.nBoundaryFaces(); ++bf) bD[m.owner()[nI + bf]] += cmptMax(internalCoeffs_[bf]);
        const auto& psi = psi_->internal();
        for (label c = 0; c < n; ++c) {
            const scalar Dfull = diag_[c] + bD[c];
            const scalar D = std::max(std::abs(Dfull), sumOff[c]);
            const scalar dD = D / alpha - Dfull;
            diag_[c] += dD;
            source_[c] += dD * psi[c];
        }
    }

    // 参考值：在全局单元 globalCell 上固定 ψ=value（纯 Neumann 压力问题）
    void setReference(glabel globalCell, scalar value) {
        const Mesh& m = mesh();
        for (label c = 0; c < m.nCells(); ++c)
            if (m.cellGlobal()[c] == globalCell) {
                source_[c] += diag_[c] * uniformT<T>(value);
                diag_[c] += diag_[c];
            }
    }

    // ---------------------------------------------------------- 求解
    std::vector<SolverPerformance> solve(const SolverControls& ctrl, bool verbose = false) {
        const Mesh& m = mesh();
        const label nt = m.nTotalCells(), nI = m.nInternalFaces();
        std::vector<SolverPerformance> out;
        std::vector<scalar> d(nt), b(nt), x(nt);
        const int nc = Traits<T>::nComponents;
        for (int cmpt = 0; cmpt < nc; ++cmpt) {
            if (nc == 3 && std::find(m.solvedComponents().begin(), m.solvedComponents().end(), cmpt) ==
                               m.solvedComponents().end())
                continue;
            for (label c = 0; c < nt; ++c) {
                d[c] = diag_[c];
                b[c] = Traits<T>::component(source_[c], cmpt);
                x[c] = Traits<T>::component(psi_->internal()[c], cmpt);
            }
            for (label bf = 0; bf < m.nBoundaryFaces(); ++bf) {
                const label c = m.owner()[nI + bf];
                d[c] += Traits<T>::component(internalCoeffs_[bf], cmpt);
                b[c] += Traits<T>::component(boundarySource_[bf], cmpt);
            }
            LduSystem sys{m, d, upper_.data(), lowerPtr(), b, symmetric()};
            auto perf = solveLdu(sys, x, ctrl);
            perf.field = psi_->name() + (nc > 1 ? std::string(1, "xyz"[cmpt]) : "");
            for (label c = 0; c < nt; ++c) Traits<T>::setComponent(psi_->internal()[c], cmpt, x[c]);
            if (verbose)
                std::cout << "  " << perf.solver << ": " << perf.field << " initial " << perf.initialResidual
                          << ", final " << perf.finalResidual << ", iters " << perf.iterations << '\n';
            out.push_back(perf);
        }
        psi_->correctBoundaryConditions();
        return out;
    }

    // ---------------------------------------------------------- 压力-速度耦合用
    // A = (diag + 边界对角的分量平均)/V
    std::vector<scalar> A() const {
        const Mesh& m = mesh();
        std::vector<scalar> a(diag_.begin(), diag_.end());
        const label nI = m.nInternalFaces();
        for (label bf = 0; bf < m.nBoundaryFaces(); ++bf) a[m.owner()[nI + bf]] += cmptAv(internalCoeffs_[bf]);
        for (label c = 0; c < m.nCells(); ++c) a[c] /= m.V()[c];
        m.halo().exchange(a);
        return a;
    }
    // H = (b − Σ 非对角·ψ_nb − (边界对角−平均)·ψ)/V
    std::vector<T> H() const {
        const Mesh& m = mesh();
        const label nI = m.nInternalFaces();
        std::vector<T> h(source_.begin(), source_.end());
        const auto& psi = psi_->internal();
        const scalar* lo = lowerPtr();
        for (label f = 0; f < nI; ++f) {
            h[m.owner()[f]] -= upper_[f] * psi[m.neighbour()[f]];
            h[m.neighbour()[f]] -= lo[f] * psi[m.owner()[f]];
        }
        for (label bf = 0; bf < m.nBoundaryFaces(); ++bf) {
            const label c = m.owner()[nI + bf];
            h[c] += boundarySource_[bf];
            const T ic = internalCoeffs_[bf];
            h[c] -= cmptMul(ic - uniformT<T>(cmptAv(ic)), psi[c]);
        }
        for (label c = 0; c < m.nCells(); ++c) h[c] = h[c] / m.V()[c];
        m.halo().exchange(h);
        return h;
    }
    // 面通量：内部面 upper·ψ_N − lower·ψ_P，边界面 internalCoeffs·ψ_P − boundarySource，外加非正交修正
    std::vector<T> flux() const {
        const Mesh& m = mesh();
        const label nI = m.nInternalFaces();
        std::vector<T> F(m.nFaces());
        const auto& psi = psi_->internal();
        const scalar* lo = lowerPtr();
        for (label f = 0; f < nI; ++f) F[f] = upper_[f] * psi[m.neighbour()[f]] - lo[f] * psi[m.owner()[f]];
        for (label bf = 0; bf < m.nBoundaryFaces(); ++bf)
            F[nI + bf] = cmptMul(internalCoeffs_[bf], psi[m.owner()[nI + bf]]) - boundarySource_[bf];
        if (!ffc_.empty())
            for (label f = 0; f < m.nFaces(); ++f) F[f] += ffc_[f];
        return F;
    }

private:
    static scalar cmptMax(scalar v) { return v; }
    static scalar cmptMax(const Vec3& v) { return std::max({v.x, v.y, v.z}); }
    static scalar cmptAv(scalar v) { return v; }
    static scalar cmptAv(const Vec3& v) { return (v.x + v.y + v.z) / 3.0; }
    void check(const FvMatrix& o) const {
        if (o.psi_ != psi_) throw std::runtime_error("FvMatrix: incompatible fields");
    }
    void addOffDiag(const FvMatrix& o, scalar s) {
        if (o.hasLower_ || hasLower_) {
            auto& lo = lower();
            const scalar* olo = o.lowerPtr();
            for (std::size_t i = 0; i < lo.size(); ++i) lo[i] += s * olo[i];
        }
        for (std::size_t i = 0; i < upper_.size(); ++i) upper_[i] += s * o.upper_[i];
    }
    void addFfc(const FvMatrix& o, scalar s) {
        if (o.ffc_.empty()) return;
        auto& f = faceFluxCorrection();
        for (std::size_t i = 0; i < f.size(); ++i) f[i] += s * o.ffc_[i];
    }

    VolField<T>* psi_;
    std::vector<scalar> diag_, upper_, lower_;
    bool hasLower_ = false;
    std::vector<T> source_, internalCoeffs_, boundarySource_, ffc_;
};

template <class T> FvMatrix<T> operator+(FvMatrix<T> a, const FvMatrix<T>& b) { return a += b; }
template <class T> FvMatrix<T> operator-(FvMatrix<T> a, const FvMatrix<T>& b) { return a -= b; }
template <class T> FvMatrix<T> operator-(FvMatrix<T> a) { return a *= -1.0; }
template <class T> FvMatrix<T> operator*(scalar s, FvMatrix<T> a) { return a *= s; }

} // namespace cfd
