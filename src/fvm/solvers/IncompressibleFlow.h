#pragma once
// 不可压 Navier-Stokes 求解器：瞬态 PIMPLE（nOuter=1 即 PISO）与稳态 SIMPLE。
//
//   ∂U/∂t + ∇·(φU) − ∇·(ν_eff ∇U) − ∇·(ν_eff (∇U)ᵀ) = −∇p + f
//   ∇·U = 0
//
// 压力-速度耦合采用 Rhie-Chow 形式（HbyA、面通量 φ 修正）；瞬态时加 ddtCorr
// 使结果与时间步长无关。ν_eff = ν + ν_t，ν_t 由湍流模型通过 setEddyViscosity 提供。

#include "fvm/discretization/Fvm.h"

#include <functional>
#include <memory>
#include <string>

namespace cfd {

struct PimpleControls {
    bool steady = false;          // true：SIMPLE 稳态
    int nOuter = 1;               // 外迭代（PIMPLE）
    int nCorr = 2;                // 压力修正次数（PISO）
    int nNonOrthCorr = 0;         // 非正交修正次数
    bool momentumPredictor = true;
    scalar alphaU = 1.0;          // 动量隐式欠松弛（稳态常用 0.7）
    scalar alphaP = 1.0;          // 压力显式欠松弛（稳态常用 0.3）
    // Rhie-Chow 时间修正（ddtCorr）系数：0 关闭（默认，瞬态 LES 推荐，数值耗散最小）；
    // <0 使用 OpenFOAM 自动系数 1−min(|φc|/|φ⁰|,1)；0~1 为固定系数。
    // 开启后准稳态解与时间步长无关，但在低粘/中心格式下会累积 Rhie-Chow 修正、增大误差。
    scalar ddtPhiCoeff = 0.0;
    SolverControls UControls{"PBiCGStab", "DILU", 1e-8, 0.0, 1000, 0, 1};
    SolverControls pControls{"PCG", "GAMG", 1e-8, 0.01, 2000, 0, 1};
    SolverControls pFinalControls{"PCG", "GAMG", 1e-8, 0.0, 2000, 0, 1};
    glabel pRefCell = 0;
    scalar pRefValue = 0.0;
    bool verbose = false;
};

// 动量体力：constant（恒定 f）或 meanVelocity（调节 f 使体平均速度等于目标值，
// 速度修正与 OpenFOAM meanVelocityForce 相同，压力梯度增量用更稳健的平均响应，见 applyForcingCorrection）
struct MomentumForcing {
    enum class Mode { None, Constant, MeanVelocity };
    Mode mode = Mode::None;
    Vec3 force;          // Constant 模式：单位质量体力
    Vec3 Ubar;           // MeanVelocity 模式：目标体平均速度
    scalar gradP = 0.0;  // MeanVelocity 模式：当前驱动压力梯度（沿 Ubar 方向，状态量）
    scalar relaxation = 1.0;
};

struct StepInfo {
    scalar continuityError = 0;
    scalar maxCo = 0;
    std::vector<SolverPerformance> perf;
    scalar UInitialResidual = 0, pInitialResidual = 0;
};

class IncompressibleFlow {
public:
    IncompressibleFlow(MeshPtr mesh, scalar nu);

    MeshPtr meshPtr() const { return mesh_; }
    const Mesh& mesh() const { return *mesh_; }
    VectorField& U() { return U_; }
    const VectorField& U() const { return U_; }
    ScalarField& p() { return p_; }
    const ScalarField& p() const { return p_; }
    std::vector<scalar>& phi() { return phi_; }
    const std::vector<scalar>& phi() const { return phi_; }
    scalar nu() const { return nu_; }

    PimpleControls controls;
    MomentumForcing forcing;
    ConvectionScheme divScheme = ConvectionScheme::parse("linear");
    GradScheme gradScheme = GradScheme::GaussLinear;
    DdtScheme ddtScheme = DdtScheme::Backward;
    bool laplacianNonOrthCorr = true;

    // 湍流/亚格子模型接口：给定单元 ν_t（长度 nTotalCells，已含幽灵），
    // 以及可选的壁面面上 ν_eff 覆盖（壁函数用），在每个外迭代开始前调用
    std::function<void(IncompressibleFlow&)> updateTurbulence;
    std::vector<scalar>& nut() { return nut_; }
    // 壁面面上的有效粘度覆盖（边界面编号 → ν_eff），为空表示不覆盖
    std::vector<scalar>& wallNuEff() { return wallNuEff_; }

    // 初始化：修正边界、由 U 计算 φ
    void initialize();
    // 续算读入后调用：保留读入的 φ，只修正边界
    void markRestarted();
    // 推进一个时间步（稳态时为一次 SIMPLE 迭代）
    StepInfo step(scalar dt);

    const TimeState& time() const { return ts_; }
    TimeState& time() { return ts_; }
    scalar courantNumber(scalar dt) const;
    std::vector<scalar> nuEffFaces() const;
    // 面上 rAU（最近一次压力修正）
    const std::vector<scalar>& rAU() const { return rAU_; }
    const std::vector<Tensor>& gradU() const { return gradU_; }
    // 更新 ∇U（湍流模型使用）
    void updateGradU();

    // 重启用：通量旧时间层
    std::vector<scalar>& phiOld() { return phiOld_; }
    std::vector<scalar>& phiOldOld() { return phiOldOld_; }

private:
    bool pressureNeedsReference() const;
    void applyForcingCorrection(const std::vector<scalar>& rAU, scalar ddtCoeff, scalar alphaU);

    MeshPtr mesh_;
    scalar nu_;
    VectorField U_;
    ScalarField p_;
    std::vector<scalar> phi_, phiOld_, phiOldOld_;
    std::vector<scalar> nut_, wallNuEff_;
    std::vector<scalar> rAU_;
    std::vector<Tensor> gradU_;
    TimeState ts_;
    bool initialized_ = false;
};

} // namespace cfd
