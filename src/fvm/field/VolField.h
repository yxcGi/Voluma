#pragma once
// 体心场：内部值（自有 + 幽灵单元）+ 边界面值 + 各 patch 的边界条件 + 旧时间层。

#include "fvm/field/BoundaryCondition.h"
#include "fvm/mesh/Mesh.h"

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfd {

template <class T> class VolField {
public:
    VolField(MeshPtr mesh, std::string name, const T& init = Traits<T>::zero())
        : mesh_(std::move(mesh)), name_(std::move(name)), internal_(mesh_->nTotalCells(), init),
          boundary_(mesh_->nBoundaryFaces(), init), bcs_(mesh_->patches().size()) {}

    VolField(const VolField&) = delete;
    VolField& operator=(const VolField&) = delete;
    VolField(VolField&&) = default;

    const std::string& name() const { return name_; }
    const Mesh& mesh() const { return *mesh_; }
    const MeshPtr& meshPtr() const { return mesh_; }

    std::vector<T>& internal() { return internal_; }
    const std::vector<T>& internal() const { return internal_; }
    T& operator[](label c) { return internal_[c]; }
    const T& operator[](label c) const { return internal_[c]; }
    std::vector<T>& boundary() { return boundary_; }
    const std::vector<T>& boundary() const { return boundary_; }
    // 边界面 f（局部面编号 ≥ nInternalFaces）上的值
    const T& bValue(label f) const { return boundary_[f - mesh_->nInternalFaces()]; }
    // patch 第 i 个面的值
    T& patchValue(label p, label i) { return boundary_[mesh_->patches()[p].start - mesh_->nInternalFaces() + i]; }
    const T& patchValue(label p, label i) const {
        return boundary_[mesh_->patches()[p].start - mesh_->nInternalFaces() + i];
    }

    // ---------------------------------------------------------- 边界条件
    void setBC(const std::string& patch, BCPtr<T> bc) {
        const label p = mesh_->findPatch(patch);
        if (p < 0) {
            if (mesh_->isEmptyPatch(patch)) return;  // 二维 empty 面：忽略
            throw std::runtime_error(name_ + ": unknown patch " + patch);
        }
        bcs_[p] = std::move(bc);
    }
    template <class BC, class... Args> BC& set(const std::string& patch, Args&&... args) {
        const label p = mesh_->findPatch(patch);
        if (p < 0) {
            if (mesh_->isEmptyPatch(patch)) throw std::runtime_error(name_ + ": cannot set BC on empty patch " + patch);
            throw std::runtime_error(name_ + ": unknown patch " + patch);
        }
        auto bc = std::make_unique<BC>(*mesh_, p, std::forward<Args>(args)...);
        BC& ref = *bc;
        bcs_[p] = std::move(bc);
        return ref;
    }
    // 便捷接口
    void fixedValue(const std::string& patch, const T& v) {
        if (!mesh_->isEmptyPatch(patch)) set<FixedValueBC<T>>(patch, v);
    }
    void zeroGradient(const std::string& patch) {
        if (!mesh_->isEmptyPatch(patch)) set<ZeroGradientBC<T>>(patch);
    }
    // 所有尚未设置的 patch 统一设为 calculated（派生场用）
    void setCalculatedRemaining() {
        for (std::size_t p = 0; p < bcs_.size(); ++p)
            if (!bcs_[p]) bcs_[p] = std::make_unique<CalculatedBC<T>>(*mesh_, label(p));
    }
    BoundaryCondition<T>& bc(label p) {
        check(p);
        return *bcs_[p];
    }
    const BoundaryCondition<T>& bc(label p) const {
        check(p);
        return *bcs_[p];
    }
    bool hasBC(label p) const { return bool(bcs_[p]); }

    void updateBCs(scalar t) {
        for (auto& b : bcs_)
            if (b) b->update(*this, t);
    }
    // halo 交换 + 按边界条件更新边界面值
    void correctBoundaryConditions() {
        mesh_->halo().exchange(internal_);
        for (std::size_t p = 0; p < bcs_.size(); ++p) {
            check(label(p));
            const Patch& pt = mesh_->patches()[p];
            bcs_[p]->evaluate(*this, boundary_.data() + (pt.start - mesh_->nInternalFaces()));
        }
    }
    void exchangeGhosts() { mesh_->halo().exchange(internal_); }

    // ---------------------------------------------------------- 赋值
    void setUniform(const T& v) {
        std::fill(internal_.begin(), internal_.end(), v);
    }
    void setFromFunction(const std::function<T(const Vec3&)>& fn) {
        for (label c = 0; c < mesh_->nTotalCells(); ++c) internal_[c] = fn(mesh_->C()[c]);
    }

    // ---------------------------------------------------------- 旧时间层
    void storeOld() {
        if (nOld_ >= 1) oldOld_ = old_;
        old_ = internal_;
        nOld_ = std::min(nOld_ + 1, 2);
    }
    const std::vector<T>& old() const { return nOld_ >= 1 ? old_ : internal_; }
    const std::vector<T>& oldOld() const { return nOld_ >= 2 ? oldOld_ : old(); }
    std::vector<T>& oldRef() { return old_; }
    std::vector<T>& oldOldRef() { return oldOld_; }
    int nOldTimes() const { return nOld_; }
    void setNOldTimes(int n) { nOld_ = n; }

private:
    void check(label p) const {
        if (!bcs_[p])
            throw std::runtime_error("field " + name_ + ": no boundary condition on patch " + mesh_->patches()[p].name);
    }

    MeshPtr mesh_;
    std::string name_;
    std::vector<T> internal_;
    std::vector<T> boundary_;
    std::vector<BCPtr<T>> bcs_;
    std::vector<T> old_, oldOld_;
    int nOld_ = 0;
};

using ScalarField = VolField<scalar>;
using VectorField = VolField<Vec3>;
using SurfaceScalar = std::vector<scalar>;  // 面场（长度 nFaces），例如体积通量 φ

// ================================================================ BC 实现
template <class T> void ZeroGradientBC<T>::evaluate(const VolField<T>& f, T* b) {
    for (label i = 0; i < this->patch().size; ++i) b[i] = f[this->faceCell(i)];
}

template <class T> void FixedGradientBC<T>::evaluate(const VolField<T>& f, T* b) {
    for (label i = 0; i < this->patch().size; ++i) b[i] = f[this->faceCell(i)] + grad_[i] / this->deltaCoeff(i);
}

template <class T> void RobinBC<T>::evaluate(const VolField<T>& f, T* b) {
    for (label i = 0; i < this->patch().size; ++i) {
        const scalar bd = b_ * this->deltaCoeff(i);
        b[i] = (c_ + bd * f[this->faceCell(i)]) / (a_ + bd);
    }
}

// 对称面：标量为零梯度
template <> inline void SymmetryBC<scalar>::evaluate(const VolField<scalar>& f, scalar* b) {
    for (label i = 0; i < patch().size; ++i) b[i] = f[faceCell(i)];
}
template <> inline scalar SymmetryBC<scalar>::valueInternalCoeff(const VolField<scalar>&, label) const { return 1.0; }
template <> inline scalar SymmetryBC<scalar>::valueBoundaryCoeff(const VolField<scalar>&, label) const { return 0.0; }
template <> inline scalar SymmetryBC<scalar>::gradientInternalCoeff(const VolField<scalar>&, label) const { return 0.0; }
template <> inline scalar SymmetryBC<scalar>::gradientBoundaryCoeff(const VolField<scalar>&, label) const { return 0.0; }

template <> inline void SymmetryBC<Vec3>::evaluate(const VolField<Vec3>& f, Vec3* b) {
    for (label i = 0; i < patch().size; ++i) {
        const Vec3 n = nf(i);
        const Vec3& u = f[faceCell(i)];
        b[i] = u - dot(n, u) * n;
    }
}
namespace detail {
inline Vec3 cmptMagN(const Vec3& n) { return {n.x * n.x, n.y * n.y, n.z * n.z}; }
} // namespace detail
template <> inline Vec3 SymmetryBC<Vec3>::valueInternalCoeff(const VolField<Vec3>&, label i) const {
    return Vec3{1, 1, 1} - detail::cmptMagN(nf(i));
}
template <> inline Vec3 SymmetryBC<Vec3>::valueBoundaryCoeff(const VolField<Vec3>& f, label i) const {
    const Vec3 n = nf(i);
    const Vec3& u = f[faceCell(i)];
    const Vec3 ub = u - dot(n, u) * n;
    return ub - cmptMultiply(valueInternalCoeff(f, i), u);
}
template <> inline Vec3 SymmetryBC<Vec3>::gradientInternalCoeff(const VolField<Vec3>&, label i) const {
    return -deltaCoeff(i) * detail::cmptMagN(nf(i));
}
template <> inline Vec3 SymmetryBC<Vec3>::gradientBoundaryCoeff(const VolField<Vec3>& f, label i) const {
    const Vec3 n = nf(i);
    const Vec3& u = f[faceCell(i)];
    const Vec3 sn = -deltaCoeff(i) * dot(n, u) * n;
    return sn - cmptMultiply(gradientInternalCoeff(f, i), u);
}

// freestreamPressure：首次更新时从同一 patch 上 U 的 freestream 条件取 U∞
inline void FreestreamPressureBC::update(const VolField<scalar>&, scalar) {
    if (haveU_) return;
    if (!ctx_ || !ctx_->U) throw std::runtime_error("freestreamPressure: give \"U\" or use flow context");
    const auto* ub = dynamic_cast<const FreestreamVelocityBC*>(&ctx_->U->bc(patchIndex()));
    if (!ub)
        throw std::runtime_error("freestreamPressure on patch " + patch().name +
                                 ": U must use freestream on the same patch, or give \"U\" explicitly");
    setDirection(ub->freestreamValue());
}

inline Vec3 PressureInletOutletVelocityBC::faceValue(const VolField<Vec3>& f, label i) const {
    const Vec3& u = f[faceCell(i)];
    if (ctx_->phiAt(face(i)) >= 0) return u;
    const Vec3 n = nf(i);
    return dot(n, u) * n + (Ut_ - dot(n, Ut_) * n);
}
inline void PressureInletOutletVelocityBC::evaluate(const VolField<Vec3>& f, Vec3* b) {
    for (label i = 0; i < patch().size; ++i) b[i] = faceValue(f, i);
}
inline Vec3 PressureInletOutletVelocityBC::valueBoundaryCoeff(const VolField<Vec3>& f, label i) const {
    return faceValue(f, i) - cmptMultiply(valueInternalCoeff(f, i), f[faceCell(i)]);
}

inline void TotalPressureBC::update(const VolField<scalar>&, scalar t) {
    const scalar p0 = p0_(t);
    for (label i = 0; i < patch().size; ++i) {
        const label f = face(i);
        values_[i] = ctx_->phiAt(f) < 0 ? p0 - 0.5 * magSqr(ctx_->U->bValue(f)) : p0;
    }
}

inline FlowRateInletBC::FlowRateInletBC(const Mesh& m, label p, std::function<scalar(scalar)> Q)
    : FixedValueBC<Vec3>(m, p, Vec3{}), Q_(std::move(Q)) {
    double a = 0;
    for (label i = 0; i < patch().size; ++i) a += m.magSf()[face(i)];
    par::allSumInPlace(&a, 1);
    area_ = a;
    if (!(area_ > 0)) throw std::runtime_error("flowRateInlet on patch " + patch().name + ": zero area");
    FlowRateInletBC::update(0.0);
}

template <class T> void AdvectiveBC<T>::update(const VolField<T>& f, scalar) {
    const TimeState& ts = *ctx_->time;
    const label n = this->patch().size;
    if (ctx_->ddt == DdtScheme::Steady || !(ts.dt > 0)) {
        std::fill(this->w_.begin(), this->w_.end(), 0.0);
        return;
    }
    if (ts.timeIndex != lastTimeIndex_) {
        // 新时间步：当前边界值成为旧时间层
        if (nOld_ >= 1) bOldOld_ = bOld_;
        for (label i = 0; i < n; ++i) bOld_[i] = f.patchValue(this->patchIndex(), i);
        nOld_ = std::min(nOld_ + 1, 2);
        lastTimeIndex_ = ts.timeIndex;
    }
    const DdtCoeffs k = ddtCoeffs(ctx_->ddt, ts, nOld_ >= 2);
    for (label i = 0; i < n; ++i) {
        const label fc = this->face(i);
        const scalar un = std::max(ctx_->phiAt(fc) / this->mesh().magSf()[fc], 0.0);
        const scalar a = un * this->deltaCoeff(i);
        const T hist = k.c0 * bOld_[i] - k.c00 * bOldOld_[i];
        if (haveInf_) {
            const scalar kk = un / lInf_;
            this->refValue_[i] = (hist + kk * fieldInf_) / (k.c + kk);
            this->w_[i] = (k.c + kk) / (k.c + a + kk);
        } else {
            this->refValue_[i] = hist / k.c;
            this->w_[i] = k.c / (k.c + a);
        }
    }
}

} // namespace cfd
