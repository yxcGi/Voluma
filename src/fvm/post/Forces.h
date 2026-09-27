#pragma once
// 壁面受力：压力与粘性力积分、力系数（升力 / 阻力 / 力矩），以及表面 Cp、Cf 分布。
//
// 不可压求解器的 p 为运动学压力 p/ρ，因此 F = ρ Σ (p − pRef) S_f + ρ Σ ν_eff (∂U_t/∂n) |S_f|。
// 粘性项取壁面相邻单元的切向速度差：τ_w = ν_eff (U_P − U_b)_t · Δ⁻¹。

#include "fvm/solvers/IncompressibleFlow.h"

#include <string>
#include <vector>

namespace cfd {

struct ForceResult {
    Vec3 pressure, viscous;
    Vec3 momentP, momentV;  // 绕参考点的力矩
    Vec3 force() const { return pressure + viscous; }
    Vec3 moment() const { return momentP + momentV; }
};

class Forces {
public:
    Forces(const IncompressibleFlow& flow, const std::vector<std::string>& patches);

    scalar rho = 1.0;
    scalar pRef = 0.0;
    Vec3 CofR;  // 力矩参考点

    ForceResult compute() const;

    // 表面分布（0 号进程写，按全局面序）：x y z nx ny nz Cp tau_x tau_y tau_z Cf
    // Cf 为壁面切应力沿 flowDir 投影 / (½ρU²)
    void writeSurface(const std::string& file, scalar Uref, const Vec3& flowDir) const;

    // 二维时的展向厚度（empty 方向的网格厚度），三维时返回 1
    scalar span() const;

private:
    const IncompressibleFlow& flow_;
    std::vector<label> patches_;
};

} // namespace cfd
