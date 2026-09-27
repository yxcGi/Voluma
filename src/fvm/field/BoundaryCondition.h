#pragma once
// 边界条件（OpenFOAM 式系数形式）。
//
// 每个边界面上：
//   面值      φ_b      = vic ∘ φ_P + vbc        （valueInternal/BoundaryCoeffs）
//   法向梯度  ∂φ/∂n_b  = gic ∘ φ_P + gbc        （gradientInternal/BoundaryCoeffs）
// ∘ 为逐分量乘法；vic/gic 对向量是逐分量系数，因此对称（滑移）面这类
// 分量相关的条件也能隐式处理。

#include "fvm/core/Types.h"
#include "fvm/mesh/Mesh.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace cfd {

template <class T> class VolField;

template <class T> inline T cmptMul(const T& a, const T& b);
template <> inline scalar cmptMul(const scalar& a, const scalar& b) { return a * b; }
template <> inline Vec3 cmptMul(const Vec3& a, const Vec3& b) { return cmptMultiply(a, b); }
template <class T> inline T uniformT(scalar s);
template <> inline scalar uniformT(scalar s) { return s; }
template <> inline Vec3 uniformT(scalar s) { return {s, s, s}; }
template <> inline Tensor uniformT(scalar s) { return {s, s, s, s, s, s, s, s, s}; }

template <class T> class BoundaryCondition {
public:
    BoundaryCondition(const Mesh& mesh, label patch) : mesh_(mesh), patch_(patch) {}
    virtual ~BoundaryCondition() = default;
    virtual std::string type() const = 0;

    // 更新依赖时间/其他场的数据（每个时间步或每次迭代调用一次）
    virtual void update(scalar /*time*/) {}
    // 由内部场计算面值，写入 bValues（长度 = patch.size）
    virtual void evaluate(const VolField<T>& f, T* bValues) = 0;

    virtual T valueInternalCoeff(const VolField<T>& f, label i) const = 0;
    virtual T valueBoundaryCoeff(const VolField<T>& f, label i) const = 0;
    virtual T gradientInternalCoeff(const VolField<T>& f, label i) const = 0;
    virtual T gradientBoundaryCoeff(const VolField<T>& f, label i) const = 0;

    // 是否直接规定了面值（压力方程判定参考值用）
    virtual bool fixesValue() const { return false; }
    // 是否为隐式可用的条件（calculated 类不可用于离散）
    virtual bool assignable() const { return false; }
    // 直接赋面值（calculated 类场，如 ν_t）
    virtual void assign(label /*i*/, const T& /*v*/) {}

    label patchIndex() const { return patch_; }
    const Patch& patch() const { return mesh_.patches()[patch_]; }
    const Mesh& mesh() const { return mesh_; }
    label faceCell(label i) const { return mesh_.owner()[patch().start + i]; }
    label face(label i) const { return patch().start + i; }
    scalar deltaCoeff(label i) const { return mesh_.nonOrthDeltaCoeffs()[patch().start + i]; }
    Vec3 nf(label i) const { const label f = face(i); return mesh_.Sf()[f] / mesh_.magSf()[f]; }

protected:
    const Mesh& mesh_;
    label patch_;
};

template <class T> using BCPtr = std::unique_ptr<BoundaryCondition<T>>;

// ------------------------------------------------------------ fixedValue
template <class T> class FixedValueBC : public BoundaryCondition<T> {
public:
    using Fn = std::function<T(const Vec3& x, scalar t)>;
    FixedValueBC(const Mesh& m, label p, const T& v) : BoundaryCondition<T>(m, p), values_(this->patch().size, v) {}
    FixedValueBC(const Mesh& m, label p, Fn fn) : BoundaryCondition<T>(m, p), values_(this->patch().size), fn_(std::move(fn)) {
        update(0.0);
    }
    std::string type() const override { return "fixedValue"; }
    void update(scalar t) override {
        if (!fn_) return;
        for (label i = 0; i < this->patch().size; ++i) values_[i] = fn_(this->mesh().Cf()[this->face(i)], t);
    }
    void evaluate(const VolField<T>&, T* b) override {
        for (label i = 0; i < label(values_.size()); ++i) b[i] = values_[i];
    }
    T valueInternalCoeff(const VolField<T>&, label) const override { return Traits<T>::zero(); }
    T valueBoundaryCoeff(const VolField<T>&, label i) const override { return values_[i]; }
    T gradientInternalCoeff(const VolField<T>&, label i) const override { return uniformT<T>(-this->deltaCoeff(i)); }
    T gradientBoundaryCoeff(const VolField<T>&, label i) const override { return this->deltaCoeff(i) * values_[i]; }
    bool fixesValue() const override { return true; }
    std::vector<T>& values() { return values_; }

protected:
    std::vector<T> values_;
    Fn fn_;
};

// ------------------------------------------------------------ zeroGradient
template <class T> class ZeroGradientBC : public BoundaryCondition<T> {
public:
    using BoundaryCondition<T>::BoundaryCondition;
    std::string type() const override { return "zeroGradient"; }
    void evaluate(const VolField<T>& f, T* b) override;
    T valueInternalCoeff(const VolField<T>&, label) const override { return uniformT<T>(1.0); }
    T valueBoundaryCoeff(const VolField<T>&, label) const override { return Traits<T>::zero(); }
    T gradientInternalCoeff(const VolField<T>&, label) const override { return Traits<T>::zero(); }
    T gradientBoundaryCoeff(const VolField<T>&, label) const override { return Traits<T>::zero(); }
};

// ------------------------------------------------------------ fixedGradient
template <class T> class FixedGradientBC : public BoundaryCondition<T> {
public:
    FixedGradientBC(const Mesh& m, label p, const T& g) : BoundaryCondition<T>(m, p), grad_(this->patch().size, g) {}
    std::string type() const override { return "fixedGradient"; }
    void evaluate(const VolField<T>& f, T* b) override;
    T valueInternalCoeff(const VolField<T>&, label) const override { return uniformT<T>(1.0); }
    T valueBoundaryCoeff(const VolField<T>&, label i) const override { return grad_[i] / this->deltaCoeff(i); }
    T gradientInternalCoeff(const VolField<T>&, label) const override { return Traits<T>::zero(); }
    T gradientBoundaryCoeff(const VolField<T>&, label i) const override { return grad_[i]; }
    std::vector<T>& gradient() { return grad_; }

private:
    std::vector<T> grad_;
};

// ------------------------------------------------------------ Robin: a·φ + b·∂φ/∂n = c
template <class T> class RobinBC : public BoundaryCondition<T> {
public:
    RobinBC(const Mesh& m, label p, scalar a, scalar b, const T& c) : BoundaryCondition<T>(m, p), a_(a), b_(b), c_(c) {}
    std::string type() const override { return "robin"; }
    void evaluate(const VolField<T>& f, T* b) override;
    T valueInternalCoeff(const VolField<T>&, label i) const override {
        const scalar bd = b_ * this->deltaCoeff(i);
        return uniformT<T>(bd / (a_ + bd));
    }
    T valueBoundaryCoeff(const VolField<T>&, label i) const override { return c_ / (a_ + b_ * this->deltaCoeff(i)); }
    T gradientInternalCoeff(const VolField<T>& f, label i) const override {
        const scalar dc = this->deltaCoeff(i);
        return uniformT<T>(dc * (b_ * dc / (a_ + b_ * dc) - 1.0));
        (void)f;
    }
    T gradientBoundaryCoeff(const VolField<T>&, label i) const override {
        const scalar dc = this->deltaCoeff(i);
        return dc * (c_ / (a_ + b_ * dc));
    }
    bool fixesValue() const override { return a_ != 0.0; }

private:
    scalar a_, b_;
    T c_;
};

// ------------------------------------------------------------ symmetry / slip
// 向量：φ_b = φ_P − (n·φ_P)n，法向分量隐式；标量：零梯度
template <class T> class SymmetryBC : public BoundaryCondition<T> {
public:
    using BoundaryCondition<T>::BoundaryCondition;
    std::string type() const override { return "symmetry"; }
    void evaluate(const VolField<T>& f, T* b) override;
    T valueInternalCoeff(const VolField<T>& f, label i) const override;
    T valueBoundaryCoeff(const VolField<T>& f, label i) const override;
    T gradientInternalCoeff(const VolField<T>& f, label i) const override;
    T gradientBoundaryCoeff(const VolField<T>& f, label i) const override;
};

// ------------------------------------------------------------ calculated（派生场，如 ν_t）
template <class T> class CalculatedBC : public BoundaryCondition<T> {
public:
    CalculatedBC(const Mesh& m, label p) : BoundaryCondition<T>(m, p), values_(this->patch().size, Traits<T>::zero()) {}
    std::string type() const override { return "calculated"; }
    void evaluate(const VolField<T>&, T* b) override {
        for (label i = 0; i < label(values_.size()); ++i) b[i] = values_[i];
    }
    T valueInternalCoeff(const VolField<T>&, label) const override { return Traits<T>::zero(); }
    T valueBoundaryCoeff(const VolField<T>&, label i) const override { return values_[i]; }
    T gradientInternalCoeff(const VolField<T>&, label i) const override { return uniformT<T>(-this->deltaCoeff(i)); }
    T gradientBoundaryCoeff(const VolField<T>&, label i) const override { return this->deltaCoeff(i) * values_[i]; }
    bool assignable() const override { return true; }
    void assign(label i, const T& v) override { values_[i] = v; }

private:
    std::vector<T> values_;
};

} // namespace cfd
