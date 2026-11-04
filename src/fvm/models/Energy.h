#pragma once
// 能量方程（不可压、常物性；温度为运动学量，热流以 q/(ρ c_p) 计，单位 K·m/s）：
//
//   ∂T/∂t + ∇·(φT) − ∇·(α_eff ∇T) = Q,   α_eff = α + ν_t/Pr_t，α = ν/Pr（或直接给 alpha）
//
// 浮力（Boussinesq）：f = −β (T − T_ref) g 作为动量体力（IncompressibleFlow::bodyForce），
//   p 为去掉参考静水压的动压；壁面 p 应用 fixedFluxPressure（fvmFlow 在有浮力时自动取它为壁面默认）。
// 壁面扩散系数：低 Re / LES（壁面 ν_eff,w = ν）时 α_eff,w = α；RANS 壁面函数（ν_eff,w > ν）时用热壁面函数
//   α_eff,w = ν y+ / T+(y+)，y+ 由壁面函数给出的切应力反算（u_τ² = ν_eff,w U_t / y，与动量方程一致）：
//   standard：Jayatilleke（T+ = Pr y+（y+ < y+_T）或 Pr_t [ln(E y+)/κ + P]，P = 9.24[(Pr/Pr_t)^{3/4} − 1][1 + 0.28 e^{−0.007 Pr/Pr_t}]）
//   kader：Kader 1981 全 y+ 混合律 T+ = Pr y+ e^{−Γ} + [2.12 ln(1 + y+) + β(Pr)] e^{−1/Γ}，
//          Γ = 0.01 (Pr y+)⁴/(1 + 5 Pr³ y+)，β = (3.85 Pr^{1/3} − 1.3)² + 2.12 ln Pr
// 边界：T 的 fixedValue / zeroGradient（绝热）/ fixedGradient / robin 等通用类型，另加
//   heatFlux { "q": 值 }：给定进入流体的热流 q/(ρ c_p)（K·m/s），按当前壁面 α_eff,w 换算成法向梯度。
// 输出：各壁面 patch 的面积、总热流、平均热流，给定 nusselt { "L", "deltaT" } 时 Nu = q̄ L/(α ΔT)；
//   writeWallHeatFlux 写逐面热流表 wallHeatFlux_<patch>.csv（x y z q T_w）。
#include "fvm/core/Json.h"
#include "fvm/solvers/IncompressibleFlow.h"

#include <iosfwd>
#include <string>
#include <vector>

namespace cfd {

class EnergyModel {
public:
    EnergyModel(IncompressibleFlow& flow, const Json& dict, const Json& boundary, const Json& initial);

    ScalarField& T() { return T_; }
    const ScalarField& T() const { return T_; }
    scalar alpha() const { return alpha_; }
    scalar Prt() const { return Prt_; }
    bool buoyant() const { return buoyant_; }

    // 推进一步（瞬态）或一次迭代（稳态）：在 flow.step 之后调用；返回 T 方程初始残差
    scalar correct();
    // 由当前 T 更新动量体力（构造时与每次 correct 后都会调用）
    void updateBuoyancy();

    struct PatchHeat {
        std::string name;
        scalar area = 0, Q = 0, qMean = 0, Nu = 0, TwMean = 0;
    };
    // 各壁面 patch 的热流汇总（全局归约；进入流体为正）
    std::vector<PatchHeat> wallHeatFlux() const;
    // 各边界面的热流（下标 f − nInternalFaces；进入流体为正；非壁面为 0）
    std::vector<scalar> heatFluxBoundaryFaces() const;
    void printWallHeatFlux(std::ostream& os) const;
    void writeWallHeatFlux(const std::string& dir) const;

    std::vector<ScalarField*> fields() { return {&T_}; }
    void writeRestart(const std::string& dir);
    void readRestart(const std::string& dir);

    // 热壁面律 T+(y+)（静态，便于测试）
    static scalar TplusJayatilleke(scalar yPlus, scalar Pr, scalar Prt, scalar kappa = 0.41, scalar E = 9.8);
    static scalar TplusKader(scalar yPlus, scalar Pr);

private:
    // 壁面面上的 α_eff（按边界面序号，非壁面为 −1）
    void updateWallAlpha();
    // 边界面 bf 的热流（进入流体为正）
    scalar faceHeatFlux(label bf, const std::vector<scalar>& gammaFaces) const;
    std::vector<scalar> gammaFaces() const;

    IncompressibleFlow& flow_;
    const Mesh& m_;
    scalar nu_, alpha_, Prt_;
    ScalarField T_;
    std::vector<scalar> Q_;         // 均匀体积源（空：无）
    std::vector<scalar> wallAlpha_; // 边界面 α_eff,w（−1：非壁面）
    std::vector<scalar> heatFluxQ_; // 各 patch 的给定热流（heatFlux 边界），无则 NaN
    std::vector<char> isWallPatch_;
    bool buoyant_ = false;
    Vec3 g_{0, 0, 0};
    scalar beta_ = 0, TRef_ = 0;
    bool useKader_ = false;
    bool wallFunction_ = true;
    scalar nusseltL_ = 0, nusseltDT_ = 0;
    ConvectionScheme divScheme_;
    scalar relax_ = 1.0;
    SolverControls controls_{"PBiCGStab", "DILU", 1e-10, 0.1, 1000, 0, 1};
    int lastTimeIndex_ = -1;
    scalar kappa_ = 0.41, E_ = 9.8;
    std::vector<scalar> gammaFacesLast_;
};

} // namespace cfd
