#pragma once
// 湍流量的入口边界条件（都是 inletOutlet 型：流入面给值、流出面零梯度）：
//   turbulentIntensity      k = 1.5 (I |U_b|)²                          （k；需要速度）
//   turbulentMixingLength   ε = Cμ^{3/4} k^{3/2}/L，ω = √k/(Cμ^{1/4} L)   （ε、ω；需要同一 patch 的 k）
//   viscosityRatio          ε = Cμ k²/(r ν)，ω = k/(r ν)，ν̃ 由 χ f_v1(χ) = r 反解（ν_t/ν = r）
// 与 Fluent 的 Intensity & Length Scale / Intensity & Viscosity Ratio 入口设定对应。
// patch 级简写 "turbulence": {"intensity": 0.05, "lengthScale": 0.07, "viscosityRatio": 10} 由 TurbulenceModel::addField 展开。

#include "fvm/field/VolField.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cfd {

class TurbulentIntensityKBC : public InletOutletBC<scalar> {
public:
    TurbulentIntensityKBC(const Mesh& m, label p, const BCContext* ctx, scalar intensity)
        : InletOutletBC<scalar>(m, p, ctx, 0.0), I_(intensity) {
        if (!ctx_->U) throw std::runtime_error("turbulentIntensity needs the velocity field");
        if (!(I_ > 0)) throw std::runtime_error("turbulentIntensity: intensity must be > 0");
    }
    std::string type() const override { return "turbulentIntensity"; }
    void update(const VolField<scalar>&, scalar t) override {
        InletOutletBC<scalar>::update(t);
        for (label i = 0; i < patch().size; ++i) {
            const scalar u = mag(ctx_->U->bValue(face(i)));
            refValue_[i] = 1.5 * (I_ * u) * (I_ * u);
        }
    }

private:
    scalar I_;
};

// ε / ω 由 k 与混合长度或粘性比给出；k 取同一 patch 上 k 场的边界值
class TurbulentScaleBC : public InletOutletBC<scalar> {
public:
    enum class Field { Epsilon, Omega };
    enum class Mode { MixingLength, ViscosityRatio };
    TurbulentScaleBC(const Mesh& m, label p, const BCContext* ctx, Field field, Mode mode, scalar value,
                     const ScalarField* k, scalar Cmu)
        : InletOutletBC<scalar>(m, p, ctx, 0.0), field_(field), mode_(mode), value_(value), k_(k), Cmu_(Cmu) {
        if (!k_) throw std::runtime_error("turbulent scale inlet needs the k field");
        if (!(value_ > 0)) throw std::runtime_error("turbulent scale inlet: lengthScale / viscosityRatio must be > 0");
    }
    std::string type() const override { return mode_ == Mode::MixingLength ? "turbulentMixingLength" : "viscosityRatio"; }
    void update(const VolField<scalar>&, scalar t) override {
        InletOutletBC<scalar>::update(t);
        const scalar nu = ctx_->nu;
        for (label i = 0; i < patch().size; ++i) {
            const scalar k = std::max(k_->bValue(face(i)), 0.0);
            scalar v;
            if (mode_ == Mode::MixingLength) {
                v = field_ == Field::Epsilon ? std::pow(Cmu_, 0.75) * std::pow(k, 1.5) / value_
                                             : std::sqrt(k) / (std::pow(Cmu_, 0.25) * value_);
            } else {
                v = field_ == Field::Epsilon ? Cmu_ * k * k / (value_ * nu) : k / (value_ * nu);
            }
            refValue_[i] = v;
        }
    }

private:
    Field field_;
    Mode mode_;
    scalar value_;
    const ScalarField* k_;
    scalar Cmu_;
};

// Spalart-Allmaras ν̃ 入口：由粘性比 r = ν_t/ν 反解 χ f_v1(χ) = r（f_v1 = χ³/(χ³ + Cv1³)），
// 或由强度 + 长度尺度：ν_t = Cμ^{1/4} √k L，k = 1.5 (I|U_b|)²
class NuTildaInletBC : public InletOutletBC<scalar> {
public:
    NuTildaInletBC(const Mesh& m, label p, const BCContext* ctx, scalar viscosityRatio, scalar intensity,
                   scalar lengthScale, scalar Cmu, scalar Cv1 = 7.1)
        : InletOutletBC<scalar>(m, p, ctx, 0.0), ratio_(viscosityRatio), I_(intensity), L_(lengthScale), Cmu_(Cmu),
          Cv1_(Cv1) {
        if (ratio_ <= 0 && !(I_ > 0 && L_ > 0))
            throw std::runtime_error("nuTilda inlet needs viscosityRatio, or intensity + lengthScale");
        if (ratio_ <= 0 && !ctx_->U) throw std::runtime_error("nuTilda inlet from intensity needs the velocity field");
    }
    std::string type() const override { return "viscosityRatio"; }
    // 反解 χ：χ⁴/(χ³ + Cv1³) = r，牛顿迭代（单调）
    static scalar chiFromRatio(scalar r, scalar Cv1) {
        if (r <= 0) return 0;
        const scalar c3 = Cv1 * Cv1 * Cv1;
        scalar x = std::max(r, 1.0);
        for (int it = 0; it < 50; ++it) {
            const scalar x3 = x * x * x, den = x3 + c3;
            const scalar g = x * x3 / den - r;
            const scalar dg = (4 * x3 * den - x * x3 * 3 * x * x) / (den * den);
            const scalar dx = g / dg;
            x -= dx;
            if (x <= 0) x = 1e-6;
            if (std::abs(dx) < 1e-12 * std::max(x, 1.0)) break;
        }
        return x;
    }
    void update(const VolField<scalar>&, scalar t) override {
        InletOutletBC<scalar>::update(t);
        const scalar nu = ctx_->nu;
        for (label i = 0; i < patch().size; ++i) {
            scalar r = ratio_;
            if (r <= 0) {
                const scalar u = mag(ctx_->U->bValue(face(i)));
                const scalar k = 1.5 * (I_ * u) * (I_ * u);
                r = std::pow(Cmu_, 0.25) * std::sqrt(k) * L_ / nu;
            }
            refValue_[i] = chiFromRatio(r, Cv1_) * nu;
        }
    }

private:
    scalar ratio_, I_, L_, Cmu_, Cv1_;
};

} // namespace cfd
