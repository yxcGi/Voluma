#pragma once
// 壁面律：由离壁 y 处的切向速度 U 求壁面摩擦速度 u_τ（LES 壁面模型、RANS 全 y+ 壁面处理共用）。
//
//   logLaw     u+ = ln(E y+)/κ（y+ < 11 取线性律）
//   spalding   Spalding (1961) 统一律 y+ = u+ + e^{−κB}(e^{κu+} − 1 − κu+ − (κu+)²/2 − (κu+)³/6)
//   musker     Musker (1979)
//   reichardt  Reichardt (1951) u+ = ln(1 + κy+)/κ + (B − ln κ/κ)(1 − e^{−y+/11} − (y+/11) e^{−y+/3})
//              （原文 κ = 0.4、常数 7.8，对应 B = 5.5；这里按 κ、B 给出常数，使对数区与其他壁面律一致）
//   powerLaw   Werner & Wengle (1991) 幂律 u+ = A y+^B（A 8.3，B 1/7；y+ < 11.8 线性），按采样点的点值反解。
//              （W&W 原文给的是在首层单元内积分的形式，要求 U 为单元平均且 Δz 为单元高度；这里统一按点值处理，
//              与其他壁面律一致）
//   ode        一维薄层方程数值解（Kawai & Larsson 2012 的平衡壁面模型）：
//                d/dy[(ν + ν_t) du/dy] = −f_p，ν_t = κ y u*(1 − e^{−y u*/(ν A+)})²
//              f_p 为切向驱动力（−∇_t p + 体力，可选），u* = (u_τ³ + u_p³)^{1/3}，u_p = (ν|f_p|)^{1/3}
//              （Duprat 等 2011 的速度尺度）。无压力梯度项时退化为对数律（κ、A+ = 17 给 B ≈ 5.0–5.2）。
//              含压力梯度时为"非平衡"壁面模型：壁面切应力方向可与采样速度方向不同（三维边界层偏斜）。

#include "fvm/core/Types.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfd {

class WallLaw {
public:
    enum class Type { LogLaw, Spalding, Musker, Reichardt, PowerLaw, ODE };
    Type type = Type::Spalding;
    scalar kappa = 0.41, B = 5.2;
    scalar A = 8.3, Bpow = 1.0 / 7.0;  // Werner-Wengle 幂律 u+ = A y+^B
    scalar Aplus = 17.0;               // ODE 模型的 Van Driest 阻尼常数
    int nPoints = 40;                  // ODE 模型的一维积分点数
    bool pressureGradient = false;     // ODE 模型是否包含切向压力梯度/体力

    static Type parse(const std::string& s) {
        if (s == "logLaw") return Type::LogLaw;
        if (s == "spalding") return Type::Spalding;
        if (s == "musker") return Type::Musker;
        if (s == "reichardt") return Type::Reichardt;
        if (s == "powerLaw" || s == "wernerWengle") return Type::PowerLaw;
        if (s == "ode" || s == "tble" || s == "equilibriumODE") return Type::ODE;
        throw std::runtime_error("wallModel.type must be logLaw, spalding, musker, reichardt, powerLaw or ode");
    }
    static std::string name(Type t) {
        switch (t) {
        case Type::LogLaw: return "logLaw";
        case Type::Spalding: return "spalding";
        case Type::Musker: return "musker";
        case Type::Reichardt: return "reichardt";
        case Type::PowerLaw: return "powerLaw";
        case Type::ODE: return "ode";
        }
        return "?";
    }

    // 给定 y+ 的 u+（显式律）
    scalar uPlusMusker(scalar yp) const {
        // Musker (1979)：κ = 0.41，s = 0.001093
        return 5.424 * std::atan((2 * yp - 8.15) / 16.7) +
               std::log10(std::pow(yp + 10.6, 9.6) / std::pow(yp * yp - 8.15 * yp + 86, 2)) - 3.52;
    }
    scalar uPlusReichardt(scalar yp) const {
        const scalar C = B - std::log(kappa) / kappa;
        return std::log(1 + kappa * yp) / kappa +
               C * (1 - std::exp(-yp / 11.0) - (yp / 11.0) * std::exp(-yp / 3.0));
    }
    // Spalding：y+ = u+ + e^{−κB}(e^{κu+} − 1 − κu+ − (κu+)²/2 − (κu+)³/6)
    scalar yPlusSpalding(scalar up) const {
        const scalar k = kappa * up;
        return up + std::exp(-kappa * B) * (std::exp(k) - 1 - k - k * k / 2 - k * k * k / 6);
    }

    // 标量形式：切向速度大小 U、离壁距离 y → u_τ（ODE 模型不含压力梯度）
    scalar uTau(scalar U, scalar y, scalar nu) const {
        if (U <= 0 || y <= 0) return 0;
        const scalar uLam = std::sqrt(nu * U / y);  // 线性律
        switch (type) {
        case Type::PowerLaw: {
            // u+ = A y+^B（点值形式），交点 y+ = A^{1/(1−B)}（8.3，1/7 → 11.8）以下取线性律。
            // u_τ^{1+B} = U ν^B/(A y^B)
            const scalar yPlusLim = std::pow(A, 1 / (1 - Bpow));
            const scalar ut = std::pow(U * std::pow(nu / y, Bpow) / A, 1 / (1 + Bpow));
            return y * ut / nu < yPlusLim ? uLam : ut;
        }
        case Type::LogLaw: {
            // U/u_τ = ln(E y u_τ/ν)/κ，y+ < 11 时取线性律
            const scalar E = std::exp(kappa * B);
            scalar ut = std::max(uLam, 1e-12);
            for (int it = 0; it < 50; ++it) {
                const scalar yp = y * ut / nu;
                if (yp < 11.0) return uLam;
                const scalar f = ut * std::log(E * yp) / kappa - U;
                const scalar df = (std::log(E * yp) + 1) / kappa;
                const scalar dut = f / df;
                ut = std::max(ut - dut, 0.5 * ut);
                if (std::abs(dut) < 1e-10 * ut) break;
            }
            return ut;
        }
        case Type::ODE: {
            const Vec3 tau = odeWallShear(Vec3{U, 0, 0}, y, nu, Vec3{});
            return std::sqrt(mag(tau));
        }
        case Type::Spalding:
        case Type::Musker:
        case Type::Reichardt: {
            // 二分法（u+ 与 y+ 之间单调）
            auto res = [&](scalar ut) {
                const scalar up = U / ut, yp = y * ut / nu;
                // 都在 u_τ 很小时为正、随 u_τ 单调减小
                if (type == Type::Spalding) return yPlusSpalding(up) - yp;
                if (type == Type::Musker) return up - uPlusMusker(yp);
                return up - uPlusReichardt(yp);
            };
            scalar lo = 1e-8 * uLam + 1e-14, hi = uLam * 1.0001 + 1e-14;
            // res(lo) > 0；线性律值处 res ≥ 0，向上加倍直到变号
            int guard = 0;
            while (res(hi) > 0 && guard++ < 200) hi *= 2;
            for (int it = 0; it < 100; ++it) {
                const scalar mid = 0.5 * (lo + hi);
                (res(mid) > 0 ? lo : hi) = mid;
                if (hi - lo < 1e-12 * hi) break;
            }
            return 0.5 * (lo + hi);
        }
        }
        return uLam;
    }

    // ODE 壁面模型（矢量形式）：采样点切向速度 Ut（相对壁面）、采样高度 h、切向驱动力 fp（单位质量；
    // 仅 pressureGradient 为 true 时使用）→ 壁面切应力矢量 τ_w（单位密度，流体作用于壁面，沿近壁流向）。
    //   (ν + ν_t) du/dy = τ_w − f_p y  ⇒  U_t = τ_w I1 − f_p I2，I1 = ∫dy/(ν+ν_t)，I2 = ∫y dy/(ν+ν_t)
    //   |τ_w| 由标量方程 T I1(T) = |U_t + f_p I2(T)| 二分求解（ν_t 通过 u* 依赖 T）
    Vec3 odeWallShear(const Vec3& Ut, scalar h, scalar nu, const Vec3& fpIn) const {
        const Vec3 fp = pressureGradient ? fpIn : Vec3{};
        const scalar U = mag(Ut);
        if (h <= 0 || (U <= 0 && mag(fp) <= 0)) return Vec3{};
        const int N = std::max(nPoints, 8);
        const scalar r = 1.2;
        const scalar rn = std::pow(r, N) - 1;
        const scalar up3 = nu * mag(fp);  // u_p³
        auto integrals = [&](scalar T, scalar& I1, scalar& I2) {
            const scalar us = std::cbrt(std::pow(T, 1.5) + up3);  // u* = (u_τ³ + u_p³)^{1/3}
            I1 = I2 = 0;
            scalar yPrev = 0, gPrev = 1.0 / nu;  // 壁面：ν_t = 0
            for (int j = 1; j <= N; ++j) {
                const scalar y = h * (std::pow(r, j) - 1) / rn;
                const scalar yp = y * us / nu;
                const scalar damp = 1 - std::exp(-yp / Aplus);
                const scalar nut = kappa * y * us * damp * damp;
                const scalar g = 1.0 / (nu + nut);
                const scalar dy = y - yPrev;
                I1 += 0.5 * dy * (gPrev + g);
                I2 += 0.5 * dy * (yPrev * gPrev + y * g);
                yPrev = y;
                gPrev = g;
            }
        };
        auto resid = [&](scalar T) {
            scalar I1, I2;
            integrals(T, I1, I2);
            return T * I1 - mag(Ut + fp * I2);
        };
        // g(0) ≤ 0；寻找 g > 0 的上界
        const scalar T0 = std::max(nu * U / h, 1e-30);
        scalar lo = 0, hi = 4 * T0;
        int guard = 0;
        while (resid(hi) <= 0 && guard++ < 200) hi *= 2;
        for (int it = 0; it < 100; ++it) {
            const scalar mid = 0.5 * (lo + hi);
            (resid(mid) > 0 ? hi : lo) = mid;
            if (hi - lo < 1e-12 * hi) break;
        }
        const scalar T = 0.5 * (lo + hi);
        scalar I1, I2;
        integrals(T, I1, I2);
        return (Ut + fp * I2) / I1;
    }

    // Kader (1981) 混合函数 Γ = −0.01 y+⁴/(1 + 5 y+)：u+ = e^Γ u+_lam + e^{1/Γ} u+_log
    static scalar kaderGamma(scalar yPlus) { return -0.01 * std::pow(yPlus, 4) / (1 + 5 * yPlus); }
    // Kader 混合后的 du+/dy+：e^Γ·1 + e^{1/Γ}/(κ y+)
    scalar kaderDuPlusDyPlus(scalar yPlus) const {
        const scalar G = kaderGamma(std::max(yPlus, 1e-12));
        const scalar eG = std::exp(G), eInvG = G < -1e-300 ? std::exp(1.0 / G) : 0.0;
        return eG + eInvG / (kappa * std::max(yPlus, 1e-12));
    }
    // Spalding 律给定 y+ 的 u+（二分）
    scalar uPlusSpalding(scalar yp) const {
        if (yp <= 0) return 0;
        scalar lo = 0, hi = yp;  // u+ ≤ y+
        for (int it = 0; it < 200; ++it) {
            const scalar mid = 0.5 * (lo + hi);
            (yPlusSpalding(mid) < yp ? lo : hi) = mid;
            if (hi - lo < 1e-12 * std::max(hi, 1.0)) break;
        }
        return 0.5 * (lo + hi);
    }
    // 壁面律的 du+/dy+（RANS 全 y+ 处理中壁面单元产生项 P = u_τ⁴/ν · D(1 − D) 用）：
    // Spalding 用解析导数 1/(dy+/du+)，其余用数值微分
    scalar duPlusDyPlus(scalar yPlus) const {
        const scalar yp = std::max(yPlus, 1e-12);
        if (type == Type::Spalding) {
            const scalar up = uPlusSpalding(yp);
            const scalar k = kappa * up;
            const scalar dydu = 1 + std::exp(-kappa * B) * kappa * (std::exp(k) - 1 - k - k * k / 2);
            return 1.0 / dydu;
        }
        const scalar h = 1e-4 * yp;
        return (uPlusOf(yp + h) - uPlusOf(std::max(yp - h, 0.0))) / (yp + h - std::max(yp - h, 0.0));
    }
    // 任意壁面律给定 y+ 的 u+（logLaw 按线性/对数切换；ODE 用 ν = 1、u_τ = 1 的无量纲积分）
    scalar uPlusOf(scalar yp) const {
        switch (type) {
        case Type::Musker: return uPlusMusker(yp);
        case Type::Reichardt: return uPlusReichardt(yp);
        case Type::Spalding: return uPlusSpalding(yp);
        case Type::PowerLaw: return yp < std::pow(A, 1 / (1 - Bpow)) ? yp : A * std::pow(yp, Bpow);
        case Type::LogLaw: return yp < 11.0 ? yp : std::log(std::exp(kappa * B) * yp) / kappa;
        case Type::ODE: {
            // 无量纲：ν = 1，u_τ = 1，f_p = 0 → u+ = ∫₀^{y+} dy/(1 + ν_t+)
            const int N = std::max(nPoints, 8);
            const scalar r = 1.2, rn = std::pow(r, N) - 1;
            scalar sum = 0, yPrev = 0, gPrev = 1.0;
            for (int j = 1; j <= N; ++j) {
                const scalar y = yp * (std::pow(r, j) - 1) / rn;
                const scalar damp = 1 - std::exp(-y / Aplus);
                const scalar g = 1.0 / (1 + kappa * y * damp * damp);
                sum += 0.5 * (y - yPrev) * (gPrev + g);
                yPrev = y, gPrev = g;
            }
            return sum;
        }
        }
        return yp;
    }
};

} // namespace cfd
