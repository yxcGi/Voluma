#pragma once
// 基本类型：标量、索引、三维向量、二阶张量。
// 所有场与几何量都以扁平数组（std::vector）存储，数据类型只包含 double，
// 便于按字节打包做 MPI 通信和写重启文件。

#include <array>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <vector>

namespace cfd {

using scalar = double;
using label = std::int32_t;   // 进程内局部编号
using glabel = std::int64_t;  // 全局编号

constexpr scalar SMALL = 1.0e-15;
constexpr scalar VSMALL = 1.0e-300;
constexpr scalar GREAT = 1.0e15;
constexpr scalar PI = 3.14159265358979323846;

// ------------------------------------------------------------------ Vec3
struct Vec3 {
    scalar x = 0, y = 0, z = 0;

    constexpr Vec3() = default;
    constexpr Vec3(scalar a, scalar b, scalar c) : x(a), y(b), z(c) {}

    constexpr scalar& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    constexpr scalar operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }

    constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    constexpr Vec3& operator*=(scalar s) { x *= s; y *= s; z *= s; return *this; }
    constexpr Vec3& operator/=(scalar s) { x /= s; y /= s; z /= s; return *this; }
    static constexpr Vec3 zero() { return {}; }
    static constexpr int nComponents = 3;
};

constexpr Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
constexpr Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
constexpr Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
constexpr Vec3 operator*(Vec3 a, scalar s) { return a *= s; }
constexpr Vec3 operator*(scalar s, Vec3 a) { return a *= s; }
constexpr Vec3 operator/(Vec3 a, scalar s) { return a /= s; }
constexpr scalar dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
constexpr Vec3 cmptMultiply(const Vec3& a, const Vec3& b) { return {a.x * b.x, a.y * b.y, a.z * b.z}; }
inline scalar mag(const Vec3& a) { return std::sqrt(dot(a, a)); }
constexpr scalar magSqr(const Vec3& a) { return dot(a, a); }
inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
}

// --------------------------------------------------------------- Tensor
// 行主序 T[i][j] = xx xy xz / yx yy yz / zx zy zz
struct Tensor {
    scalar xx = 0, xy = 0, xz = 0, yx = 0, yy = 0, yz = 0, zx = 0, zy = 0, zz = 0;

    constexpr scalar& operator()(int i, int j) { return (&xx)[3 * i + j]; }
    constexpr scalar operator()(int i, int j) const { return (&xx)[3 * i + j]; }

    constexpr Tensor& operator+=(const Tensor& o) {
        for (int k = 0; k < 9; ++k) (&xx)[k] += (&o.xx)[k];
        return *this;
    }
    constexpr Tensor& operator-=(const Tensor& o) {
        for (int k = 0; k < 9; ++k) (&xx)[k] -= (&o.xx)[k];
        return *this;
    }
    constexpr Tensor& operator*=(scalar s) {
        for (int k = 0; k < 9; ++k) (&xx)[k] *= s;
        return *this;
    }
    static constexpr Tensor zero() { return {}; }
    static constexpr int nComponents = 9;
};

constexpr Tensor operator+(Tensor a, const Tensor& b) { return a += b; }
constexpr Tensor operator-(Tensor a, const Tensor& b) { return a -= b; }
constexpr Tensor operator*(Tensor a, scalar s) { return a *= s; }
constexpr Tensor operator*(scalar s, Tensor a) { return a *= s; }
constexpr Tensor operator/(Tensor a, scalar s) { return a *= (1.0 / s); }

// 外积 a⊗b
constexpr Tensor outer(const Vec3& a, const Vec3& b) {
    return {a.x * b.x, a.x * b.y, a.x * b.z, a.y * b.x, a.y * b.y, a.y * b.z, a.z * b.x, a.z * b.y, a.z * b.z};
}
constexpr Tensor transpose(const Tensor& t) {
    return {t.xx, t.yx, t.zx, t.xy, t.yy, t.zy, t.xz, t.yz, t.zz};
}
constexpr scalar trace(const Tensor& t) { return t.xx + t.yy + t.zz; }
constexpr Tensor symm(const Tensor& t) { return 0.5 * (t + transpose(t)); }
constexpr Tensor skew(const Tensor& t) { return 0.5 * (t - transpose(t)); }
constexpr Tensor identity() { return {1, 0, 0, 0, 1, 0, 0, 0, 1}; }
constexpr Tensor dev(const Tensor& t) { return t - (trace(t) / 3.0) * identity(); }
constexpr scalar doubleDot(const Tensor& a, const Tensor& b) {
    scalar s = 0;
    for (int k = 0; k < 9; ++k) s += (&a.xx)[k] * (&b.xx)[k];
    return s;
}
constexpr scalar magSqr(const Tensor& t) { return doubleDot(t, t); }
inline scalar mag(const Tensor& t) { return std::sqrt(magSqr(t)); }
// 张量点乘向量 T·v
constexpr Vec3 dot(const Tensor& t, const Vec3& v) {
    return {t.xx * v.x + t.xy * v.y + t.xz * v.z, t.yx * v.x + t.yy * v.y + t.yz * v.z,
            t.zx * v.x + t.zy * v.y + t.zz * v.z};
}
// 向量点乘张量 v·T
constexpr Vec3 dot(const Vec3& v, const Tensor& t) {
    return {v.x * t.xx + v.y * t.yx + v.z * t.zx, v.x * t.xy + v.y * t.yy + v.z * t.zy,
            v.x * t.xz + v.y * t.yz + v.z * t.zz};
}
constexpr Tensor dot(const Tensor& a, const Tensor& b) {
    Tensor r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            scalar s = 0;
            for (int k = 0; k < 3; ++k) s += a(i, k) * b(k, j);
            r(i, j) = s;
        }
    return r;
}

// ------------------------------------------------------ 类型萃取（梯度类型等）
template <class T> struct Traits;
template <> struct Traits<scalar> {
    static constexpr int nComponents = 1;
    using Grad = Vec3;
    static scalar component(scalar v, int) { return v; }
    static void setComponent(scalar& v, int, scalar s) { v = s; }
    static scalar zero() { return 0.0; }
    static scalar magnitude(scalar v) { return std::abs(v); }
};
template <> struct Traits<Vec3> {
    static constexpr int nComponents = 3;
    using Grad = Tensor;
    static scalar component(const Vec3& v, int c) { return v[c]; }
    static void setComponent(Vec3& v, int c, scalar s) { v[c] = s; }
    static Vec3 zero() { return {}; }
    static scalar magnitude(const Vec3& v) { return mag(v); }
};
template <> struct Traits<Tensor> {
    static constexpr int nComponents = 9;
    static scalar component(const Tensor& v, int c) { return (&v.xx)[c]; }
    static void setComponent(Tensor& v, int c, scalar s) { (&v.xx)[c] = s; }
    static Tensor zero() { return {}; }
    static scalar magnitude(const Tensor& v) { return mag(v); }
};

// 梯度与面法向的缩并：grad(φ)·d，对标量得标量，对向量得向量（(∇u)ᵀ 约定见下）
// 约定：Vec3 的梯度 G(i,j) = ∂u_j/∂x_i（与 OpenFOAM 一致），于是 d·G 为方向导数。
inline scalar gradDot(const Vec3& d, const Vec3& g) { return dot(d, g); }
inline Vec3 gradDot(const Vec3& d, const Tensor& g) { return dot(d, g); }

// 面值×面积向量 → 梯度贡献（高斯公式）
inline Vec3 outerSf(const Vec3& Sf, scalar phi) { return Sf * phi; }
inline Tensor outerSf(const Vec3& Sf, const Vec3& phi) { return outer(Sf, phi); }

} // namespace cfd
