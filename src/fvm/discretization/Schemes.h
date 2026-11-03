#pragma once
// 离散格式选择

#include "fvm/core/Types.h"

#include <sstream>
#include <stdexcept>
#include <string>

namespace cfd {

// 对流格式
//   隐式权重类：upwind | linear | limitedLinear k | vanLeer | MUSCL | minmod
//   梯度修正类（延迟修正）：linearUpwind | LUST（0.75 linear + 0.25 linearUpwind）
struct ConvectionScheme {
    enum class Type { Upwind, Linear, LimitedLinear, VanLeer, MUSCL, Minmod, LinearUpwind, LUST };
    Type type = Type::Linear;
    scalar k = 1.0;

    static ConvectionScheme parse(const std::string& s) {
        std::istringstream is(s);
        std::string name;
        is >> name;
        ConvectionScheme c;
        if (name == "upwind") c.type = Type::Upwind;
        else if (name == "linear") c.type = Type::Linear;
        else if (name == "limitedLinear") {
            c.type = Type::LimitedLinear;
            if (!(is >> c.k)) c.k = 1.0;
        } else if (name == "vanLeer") c.type = Type::VanLeer;
        else if (name == "MUSCL") c.type = Type::MUSCL;
        else if (name == "minmod") c.type = Type::Minmod;
        else if (name == "linearUpwind") c.type = Type::LinearUpwind;
        else if (name == "LUST") c.type = Type::LUST;
        else throw std::runtime_error("unknown convection scheme: " + s);
        return c;
    }
    bool needsGradient() const { return type != Type::Upwind && type != Type::Linear; }
};

enum class GradScheme { GaussLinear, LeastSquares };
inline GradScheme parseGradScheme(const std::string& s) {
    if (s == "Gauss linear" || s == "Gauss") return GradScheme::GaussLinear;
    if (s == "leastSquares") return GradScheme::LeastSquares;
    throw std::runtime_error("unknown gradient scheme: " + s);
}

enum class DdtScheme { Steady, Euler, Backward };
inline DdtScheme parseDdtScheme(const std::string& s) {
    if (s == "steadyState") return DdtScheme::Steady;
    if (s == "Euler") return DdtScheme::Euler;
    if (s == "backward") return DdtScheme::Backward;
    throw std::runtime_error("unknown ddt scheme: " + s);
}

// 时间步信息（ddt 系数）
struct TimeState {
    scalar time = 0;
    scalar dt = 1;      // 当前步长
    scalar dt0 = -1;    // 上一步步长（<0 表示没有更早的时间层）
    int timeIndex = 0;
};

// 时间导数离散系数（已除以 Δt）：
//   Euler:    (ψ − ψ⁰)/Δt                       c = c0 = 1/Δt
//   backward: (c ψ − c0 ψ⁰ + c00 ψ⁰⁰)/Δt（变步长 BDF2），首步或没有 ψ⁰⁰ 时退化为 Euler
struct DdtCoeffs {
    scalar c = 0, c0 = 0, c00 = 0;
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

} // namespace cfd
