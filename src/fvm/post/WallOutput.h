#pragma once
// 壁面 patch 输出：每个壁面面的 p、Cp、壁面切应力 τ_w（运动学 τ/ρ）、Cf、y+，以及有能量方程时的热流 q 与壁温 T_w。
//   Cp = (p − pRef)/(½ U_ref²)；Cf = τ_w·t/(½ U_ref²)，给定 flowDir 时 t 为其单位向量，否则取 |τ_w|；
//   τ_w = ν_eff,w (U_P − U_b)_t Δ⁻¹（壁面函数网格上 ν_eff,w 由湍流模型给出，与 Forces 相同）；
//   y+ 取湍流模型的值（壁面单元中心），层流 / 无模型时 y+ = y √|τ_w| / ν。
// 写 CSV（wall_<patch>.csv，按全局面序）与 VTK 表面（SurfaceWriter，surface_<patch>.pvd）。

#include "fvm/io/SurfaceWriter.h"
#include "fvm/solvers/IncompressibleFlow.h"

#include <memory>
#include <string>
#include <vector>

namespace cfd {

class TurbulenceModel;
class EnergyModel;

class WallOutput {
public:
    // patches 为空时取全部壁面 patch；turb / energy 可为空
    WallOutput(const IncompressibleFlow& flow, const TurbulenceModel* turb, const EnergyModel* energy,
               const std::vector<std::string>& patches);

    scalar Uref = 1.0, pRef = 0.0;
    Vec3 flowDir{0, 0, 0};

    // 按当前场重算各面量
    void update();
    void writeCsv(const std::string& dir) const;
    // 写一帧 VTK 表面（首次调用时在 dir 下创建写出器）
    void writeVtk(const std::string& dir, scalar time);

    const std::vector<label>& patches() const { return patches_; }
    // 面量按边界面编号 f − nInternalFaces 存放（非所选 patch 为 0）
    const std::vector<Vec3>& tau() const { return tau_; }
    const std::vector<scalar>& Cp() const { return Cp_; }
    const std::vector<scalar>& Cf() const { return Cf_; }
    const std::vector<scalar>& yPlus() const { return yPlus_; }
    const std::vector<scalar>& heatFlux() const { return q_; }

private:
    const IncompressibleFlow& flow_;
    const TurbulenceModel* turb_;
    const EnergyModel* energy_;
    std::vector<label> patches_;
    std::vector<scalar> p_, Cp_, tauMag_, Cf_, yPlus_, q_, Tw_;
    std::vector<Vec3> tau_;
    std::unique_ptr<SurfaceWriter> surf_;
};

} // namespace cfd
