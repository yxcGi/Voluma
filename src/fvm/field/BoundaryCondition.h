#pragma once
// 边界条件（OpenFOAM 式系数形式）。
//
// 每个边界面上：
//   面值      φ_b      = vic ∘ φ_P + vbc        （valueInternal/BoundaryCoeffs）
//   法向梯度  ∂φ/∂n_b  = gic ∘ φ_P + gbc        （gradientInternal/BoundaryCoeffs）
// ∘ 为逐分量乘法；vic/gic 对向量是逐分量系数，因此对称（滑移）面这类
// 分量相关的条件也能隐式处理。
//
// 需要流动信息（通量 φ、速度 U、时间步）的条件通过 BCContext 取得，由求解器维护其中的指针。
//
// 条件一览（JSON 名）：
//   fixedValue / noSlip / movingWall / rotatingWall / flowRateInlet（定值类，fixesValue）
//   zeroGradient / fixedGradient / fixedFluxPressure（梯度类）
//   robin、symmetry / slip、calculated
//   mixed（w·定值 + (1−w)·定梯度）及其派生：inletOutlet、freestream、freestreamPressure、advective
//   totalPressure（p）、pressureInletOutletVelocity（U）

#include "fvm/core/Types.h"
#include "fvm/discretization/Schemes.h"
#include "fvm/mesh/Mesh.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
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

// 求解器提供给边界条件的运行时数据（指向求解器内部的量，生命周期由求解器保证）
struct BCContext {
    const std::vector<scalar>* phi = nullptr;      // 面体积通量（正为流出）
    const VolField<Vec3>* U = nullptr;             // 速度场
    const TimeState* time = nullptr;               // 时间步信息
    DdtScheme ddt = DdtScheme::Euler;              // 时间格式（advective 用）
    const std::vector<scalar>* phiHbyA = nullptr;  // 压力方程组装期间：HbyA 通量（fixedFluxPressure 用）
    const std::vector<scalar>* rAUf = nullptr;     //               及面上 rAU
    scalar nu = 0;                                 // 运动粘度
    std::string caseDir = ".";                     // 算例目录（表文件相对路径）
    scalar phiAt(label f) const { return phi && f < label(phi->size()) ? (*phi)[f] : 0.0; }
};

template <class T> class BoundaryCondition {
public:
    BoundaryCondition(const Mesh& mesh, label patch) : mesh_(mesh), patch_(patch) {}
    virtual ~BoundaryCondition() = default;
    virtual std::string type() const = 0;

    // 更新依赖时间/其他场的数据（每个时间步或每次迭代调用一次）；默认转给只依赖时间的版本
    virtual void update(const VolField<T>& /*f*/, scalar t) { update(t); }
    virtual void update(scalar /*time*/) {}
    // 由内部场计算面值，写入 bValues（长度 = patch.size）
    virtual void evaluate(const VolField<T>& f, T* bValues) = 0;

    virtual T valueInternalCoeff(const VolField<T>& f, label i) const = 0;
    virtual T valueBoundaryCoeff(const VolField<T>& f, label i) const = 0;
    virtual T gradientInternalCoeff(const VolField<T>& f, label i) const = 0;
    virtual T gradientBoundaryCoeff(const VolField<T>& f, label i) const = 0;

    // 是否（至少部分）规定了面值（压力方程判定参考值用）
    virtual bool fixesValue() const { return false; }
    // 面上的体积通量是否由边界值决定（定值类、对称面）；否则压力方程用内部 HbyA 外推（OpenFOAM assignable）
    virtual bool fixesFlux() const { return fixesValue(); }
    // 是否为 calculated 类（可直接赋面值的派生场，如 ν_t）
    virtual bool assignable() const { return false; }
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
// 值可为常数或函数 fn(面心 x, 单位外法向 n, 时间 t)
template <class T> class FixedValueBC : public BoundaryCondition<T> {
public:
    using Fn = std::function<T(const Vec3& x, const Vec3& n, scalar t)>;
    FixedValueBC(const Mesh& m, label p, const T& v) : BoundaryCondition<T>(m, p), values_(this->patch().size, v) {}
    FixedValueBC(const Mesh& m, label p, Fn fn) : BoundaryCondition<T>(m, p), values_(this->patch().size), fn_(std::move(fn)) {
        FixedValueBC::update(0.0);
    }
    std::string type() const override { return "fixedValue"; }
    using BoundaryCondition<T>::update;
    void update(scalar t) override {
        if (!fn_) return;
        for (label i = 0; i < this->patch().size; ++i) values_[i] = fn_(this->mesh().Cf()[this->face(i)], this->nf(i), t);
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
    const std::vector<T>& values() const { return values_; }

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

protected:
    std::vector<T> grad_;
};

// ------------------------------------------------------------ fixedFluxPressure
// 压力梯度由求解器在压力方程组装时设定，使边界通量等于边界速度给出的通量：
//   ∂p/∂n = (φ_HbyA − U_b·S_f) / (|S_f| rAU_f)
// 没有体力进入通量时等价于 zeroGradient；浮力（Boussinesq）等体力加入 φ_HbyA 后必须用它（同 OpenFOAM）。
class FixedFluxPressureBC : public FixedGradientBC<scalar> {
public:
    FixedFluxPressureBC(const Mesh& m, label p) : FixedGradientBC<scalar>(m, p, 0.0) {}
    std::string type() const override { return "fixedFluxPressure"; }
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
    bool fixesFlux() const override { return true; }  // 法向速度为零：通量由边界值（0）决定
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

// ------------------------------------------------------------ mixed
// φ_b = w·refValue + (1 − w)(φ_P + refGrad/Δ)，w ∈ [0,1] 逐面（w=1 定值，w=0 定梯度）
template <class T> class MixedBC : public BoundaryCondition<T> {
public:
    MixedBC(const Mesh& m, label p, const T& refValue = Traits<T>::zero(), const T& refGrad = Traits<T>::zero(),
            scalar valueFraction = 1.0)
        : BoundaryCondition<T>(m, p), refValue_(this->patch().size, refValue), refGrad_(this->patch().size, refGrad),
          w_(this->patch().size, valueFraction) {}
    std::string type() const override { return "mixed"; }
    void evaluate(const VolField<T>& f, T* b) override {
        for (label i = 0; i < this->patch().size; ++i)
            b[i] = w_[i] * refValue_[i] + (1.0 - w_[i]) * (f[this->faceCell(i)] + refGrad_[i] / this->deltaCoeff(i));
    }
    T valueInternalCoeff(const VolField<T>&, label i) const override { return uniformT<T>(1.0 - w_[i]); }
    T valueBoundaryCoeff(const VolField<T>&, label i) const override {
        return w_[i] * refValue_[i] + (1.0 - w_[i]) * (refGrad_[i] / this->deltaCoeff(i));
    }
    T gradientInternalCoeff(const VolField<T>&, label i) const override { return uniformT<T>(-w_[i] * this->deltaCoeff(i)); }
    T gradientBoundaryCoeff(const VolField<T>&, label i) const override {
        return w_[i] * this->deltaCoeff(i) * refValue_[i] + (1.0 - w_[i]) * refGrad_[i];
    }
    bool fixesValue() const override {
        for (scalar w : w_)
            if (w > 0) return true;
        return false;
    }
    bool fixesFlux() const override { return false; }
    std::vector<T>& refValue() { return refValue_; }
    std::vector<T>& refGrad() { return refGrad_; }
    std::vector<scalar>& valueFraction() { return w_; }

protected:
    std::vector<T> refValue_, refGrad_;
    std::vector<scalar> w_;
};

// ------------------------------------------------------------ inletOutlet
// 流出（φ ≥ 0）零梯度，回流（φ < 0）取 inletValue（常数或 fn(x, n, t)）。出口防回流用。
template <class T> class InletOutletBC : public MixedBC<T> {
public:
    using Fn = typename FixedValueBC<T>::Fn;
    InletOutletBC(const Mesh& m, label p, const BCContext* ctx, const T& inletValue)
        : MixedBC<T>(m, p, inletValue, Traits<T>::zero(), 0.0), ctx_(ctx) {
        if (!ctx_) throw std::runtime_error("inletOutlet needs flow context (flux)");
    }
    InletOutletBC(const Mesh& m, label p, const BCContext* ctx, Fn fn)
        : MixedBC<T>(m, p, Traits<T>::zero(), Traits<T>::zero(), 0.0), ctx_(ctx), fn_(std::move(fn)) {
        if (!ctx_) throw std::runtime_error("inletOutlet needs flow context (flux)");
        InletOutletBC::update(0.0);
    }
    std::string type() const override { return "inletOutlet"; }
    using MixedBC<T>::update;
    void update(scalar t) override {
        for (label i = 0; i < this->patch().size; ++i) {
            if (fn_) this->refValue_[i] = fn_(this->mesh().Cf()[this->face(i)], this->nf(i), t);
            this->w_[i] = ctx_->phiAt(this->face(i)) < 0 ? 1.0 : 0.0;
        }
    }

protected:
    const BCContext* ctx_;
    Fn fn_;
};

// ------------------------------------------------------------ freestream（速度）
// 按来流方向定进出：w = ½ − ½ (U∞·n)/|U∞|（流入面定值 U∞，流出面零梯度），不随当地通量跳变
class FreestreamVelocityBC : public MixedBC<Vec3> {
public:
    FreestreamVelocityBC(const Mesh& m, label p, const Vec3& Uinf) : MixedBC<Vec3>(m, p, Uinf), Uinf_(Uinf) {
        const scalar mu = mag(Uinf_);
        for (label i = 0; i < patch().size; ++i) w_[i] = mu > 0 ? std::clamp(0.5 - 0.5 * dot(Uinf_, nf(i)) / mu, 0.0, 1.0) : 1.0;
    }
    std::string type() const override { return "freestream"; }
    const Vec3& freestreamValue() const { return Uinf_; }

private:
    Vec3 Uinf_;
};

// ------------------------------------------------------------ freestreamPressure
// w = ½ + ½ (U∞·n)/|U∞|：流出面定值 p∞，流入面零梯度。U∞ 取同一 patch 上 U 的 freestream 条件，或显式给出
class FreestreamPressureBC : public MixedBC<scalar> {
public:
    FreestreamPressureBC(const Mesh& m, label p, const BCContext* ctx, scalar pInf, const Vec3* Uinf)
        : MixedBC<scalar>(m, p, pInf), ctx_(ctx), haveU_(Uinf != nullptr) {
        if (Uinf) setDirection(*Uinf);
    }
    std::string type() const override { return "freestreamPressure"; }
    void update(const VolField<scalar>& f, scalar t) override;

private:
    void setDirection(const Vec3& Uinf) {
        const scalar mu = mag(Uinf);
        for (label i = 0; i < patch().size; ++i) w_[i] = mu > 0 ? std::clamp(0.5 + 0.5 * dot(Uinf, nf(i)) / mu, 0.0, 1.0) : 1.0;
        haveU_ = true;
    }
    const BCContext* ctx_;
    bool haveU_ = false;
};

// ------------------------------------------------------------ totalPressure（p）
// 流入面 p = p0 − ½|U_b|²，流出面 p = p0（不可压运动学压力）；p0 可随时间变化
class TotalPressureBC : public FixedValueBC<scalar> {
public:
    TotalPressureBC(const Mesh& m, label p, const BCContext* ctx, std::function<scalar(scalar)> p0)
        : FixedValueBC<scalar>(m, p, 0.0), ctx_(ctx), p0_(std::move(p0)) {
        if (!ctx_ || !ctx_->U) throw std::runtime_error("totalPressure needs flow context (U, flux)");
    }
    std::string type() const override { return "totalPressure"; }
    void update(const VolField<scalar>& f, scalar t) override;

private:
    const BCContext* ctx_;
    std::function<scalar(scalar)> p0_;
};

// ------------------------------------------------------------ pressureInletOutletVelocity（U）
// 与 totalPressure / 定压边界配合：流出面零梯度；流入面法向分量零梯度、切向分量取 tangentialVelocity（默认 0）
class PressureInletOutletVelocityBC : public BoundaryCondition<Vec3> {
public:
    PressureInletOutletVelocityBC(const Mesh& m, label p, const BCContext* ctx, const Vec3& Ut = {})
        : BoundaryCondition<Vec3>(m, p), ctx_(ctx), Ut_(Ut) {
        if (!ctx_) throw std::runtime_error("pressureInletOutletVelocity needs flow context (flux)");
    }
    std::string type() const override { return "pressureInletOutletVelocity"; }
    void evaluate(const VolField<Vec3>& f, Vec3* b) override;
    Vec3 valueInternalCoeff(const VolField<Vec3>&, label i) const override {
        if (ctx_->phiAt(face(i)) >= 0) return {1, 1, 1};
        const Vec3 n = nf(i);
        return {n.x * n.x, n.y * n.y, n.z * n.z};
    }
    Vec3 valueBoundaryCoeff(const VolField<Vec3>& f, label i) const override;
    Vec3 gradientInternalCoeff(const VolField<Vec3>& f, label i) const override {
        return deltaCoeff(i) * (valueInternalCoeff(f, i) - Vec3{1, 1, 1});
    }
    Vec3 gradientBoundaryCoeff(const VolField<Vec3>& f, label i) const override {
        return deltaCoeff(i) * valueBoundaryCoeff(f, i);
    }
    bool fixesFlux() const override { return false; }

private:
    Vec3 faceValue(const VolField<Vec3>& f, label i) const;
    const BCContext* ctx_;
    Vec3 Ut_;
};

// ------------------------------------------------------------ flowRateInlet（U）
// 均匀法向入流速度 = Q(t)/A，A 为 patch 总面积（全局），方向为内法向
class FlowRateInletBC : public FixedValueBC<Vec3> {
public:
    FlowRateInletBC(const Mesh& m, label p, std::function<scalar(scalar)> Q);
    std::string type() const override { return "flowRateInlet"; }
    using FixedValueBC<Vec3>::update;
    void update(scalar t) override {
        const scalar un = Q_(t) / area_;
        for (label i = 0; i < patch().size; ++i) values_[i] = -un * nf(i);
    }
    scalar area() const { return area_; }

private:
    std::function<scalar(scalar)> Q_;
    scalar area_ = 0;
};

// ------------------------------------------------------------ advective（对流出流）
// ∂φ/∂t + U_n ∂φ/∂n = 0 在边界上隐式离散（Euler/backward 与场的时间格式一致）：
//   φ_b = w·refValue + (1−w) φ_P，w = c/(c + U_n Δ)，refValue = (c0 φ_b⁰ − c00 φ_b⁰⁰)/c
// 可选向远场值松弛：lInf > 0 时加 (U_n/lInf)(φ_b − φ_inf)（OpenFOAM advective 的 fieldInf/lInf）。
// 稳态或流入（U_n < 0）时退化为零梯度。
template <class T> class AdvectiveBC : public MixedBC<T> {
public:
    AdvectiveBC(const Mesh& m, label p, const BCContext* ctx, bool haveInf = false, const T& fieldInf = Traits<T>::zero(),
                scalar lInf = 0)
        : MixedBC<T>(m, p, Traits<T>::zero(), Traits<T>::zero(), 0.0), ctx_(ctx), haveInf_(haveInf && lInf > 0),
          fieldInf_(fieldInf), lInf_(lInf), bOld_(this->patch().size, Traits<T>::zero()),
          bOldOld_(this->patch().size, Traits<T>::zero()) {
        if (!ctx_ || !ctx_->time) throw std::runtime_error("advective needs flow context (flux, time)");
    }
    std::string type() const override { return "advective"; }
    void update(const VolField<T>& f, scalar) override;

private:
    const BCContext* ctx_;
    bool haveInf_;
    T fieldInf_;
    scalar lInf_;
    std::vector<T> bOld_, bOldOld_;
    int nOld_ = 0, lastTimeIndex_ = -1;
};

} // namespace cfd
