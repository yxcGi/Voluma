#pragma once
// 湍流模型（不可压，运动学量）。所有模型给 IncompressibleFlow 提供单元 ν_t 以及壁面面上的 ν_eff，
// 在每个外迭代开始前（使用最新的 U、φ）更新一次；输运方程的旧时间层每个时间步保存一次。
//
//   LES   Smagorinsky（可选 Van Driest 阻尼）、WALE、Vreman、sigma、dynamicSmagorinsky（Germano-Lilly）；
//         可选壁面模型（WMLES）：spalding / musker / reichardt / powerLaw（Werner-Wengle）/ logLaw /
//         ode（一维薄层方程数值解，可含切向压力梯度 = 非平衡壁面模型），采样高度可设，给出壁面切应力
//   RANS  kEpsilon（标准）、RNGkEpsilon、realizableKE（Shih）、kOmega（Wilcox 1998）、kOmegaSST（Menter 2003）、
//         BSLEARSM、SpalartAllmaras（无 ft2）
//         壁面处理 wallTreatment：standard（nutkWallFunction + y+ 切换）或 enhanced（全 y+：Spalding 律壁面切应力、
//         Kader 混合的壁面单元产生项、Wolfshtein 两层 ε）；wallFunctions=false 时为低 Re 处理（壁面 ν_t = 0）
//
// 系数与 OpenFOAM 同名模型默认值一致，可在 JSON 的 "coeffs" 中覆盖。

#include "fvm/core/Json.h"
#include "fvm/solvers/IncompressibleFlow.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace cfd {

class TurbulenceModel {
public:
    TurbulenceModel(IncompressibleFlow& flow, const Json& dict, const Json& boundary, const Json& initial);
    virtual ~TurbulenceModel() = default;

    // 由 JSON 创建："model": "laminar" 时返回空指针
    static std::unique_ptr<TurbulenceModel> create(IncompressibleFlow& flow, const Json& dict, const Json& boundary,
                                                   const Json& initial);

    virtual std::string type() const = 0;
    // 更新模型并写入 flow.nut()、flow.wallNuEff()
    void correct();

    // 模型输运/输出场（k、epsilon、omega、nuTilda；以及 nut、yPlus）
    std::vector<ScalarField*> fields() { return fields_; }
    const std::vector<scalar>& nut() const { return flow_.nut(); }
    // 各壁面面的 y+（与壁面单元中心，最近一次 correct）
    std::vector<scalar> yPlusWall() const;
    // 各边界面的 y+（下标 f − nInternalFaces；非壁面为 −1）
    std::vector<scalar> yPlusBoundaryFaces() const;
    // 打印壁面 y+ 统计（以及 ν_t 限制器、有界化作用的单元数）
    void printYPlus(std::ostream& os) const;
    // 混合初始化：没有在 initial 中给出初值的模型场取入流面（φ < 0 的非壁面边界）边界值的面积加权平均
    void initializeFromInlets();
    // 最近一次 correct 中各场被 bound() 修正的单元数（全局）
    const std::map<std::string, glabel>& boundedCells() const { return nBounded_; }
    // 续算：模型场及旧时间层（与 restart::write/read 同目录）
    void writeRestart(const std::string& dir);
    void readRestart(const std::string& dir);

protected:
    virtual void correctModel() = 0;
    // 设置字段的边界条件：非壁面按 JSON，壁面按模型要求（wallDefault）
    // 创建并登记模型场；wallDefault 为 "zeroGradient" 或 "zero"（定值 0）
    ScalarField& addField(const std::string& name, const std::string& wallDefault, scalar initDefault);
    // 湍流量入口（强度 / 长度尺度 / 粘性比），见 TurbulenceBCs.h；非该类条件时返回 false
    bool setTurbulentInlet(ScalarField& f, const std::string& patch, const Json& d);
    // 输运方程：ddt(ψ) + div(φ,ψ) − ∇·(Γ∇ψ) + Sp·ψ = Su
    // gammaCell：单元扩散系数（插值到面），Sp、Su：单位体积，Sp ≥ 0
    void solveTransport(ScalarField& psi, const std::vector<scalar>& gammaCell, const std::vector<scalar>& Sp,
                        const std::vector<scalar>& Su, const std::vector<label>& fixCells = {},
                        const std::vector<scalar>& fixValues = {});
    // 有界化（OpenFOAM bound 同法）：低于下限（或非有限）的单元取相邻单元有效值的平均，不足则取下限；统计单元数
    void bound(ScalarField& f, scalar minValue);
    scalar coeff(const std::string& name, scalar def);

    // 壁面面信息（本进程）
    struct WallFace {
        label face, cell;
        scalar y;     // 单元中心到壁面距离（法向）
        Vec3 n;       // 单位外法向
        scalar magUt; // 壁面单元切向速度大小（相对壁面）
    };
    void updateWallFaces();
    // 设置壁面面上的 ν_eff = ν + ν_t,w（nutw 按 wallFaces_ 顺序）
    void setWallNut(const std::vector<scalar>& nutw);

    IncompressibleFlow& flow_;
    const Mesh& m_;
    scalar nu_;
    const Json& dict_;
    const Json& boundary_;
    const Json& initial_;
    std::vector<ScalarField*> fields_;
    std::vector<std::unique_ptr<ScalarField>> owned_;
    std::vector<WallFace> wallFaces_;
    std::vector<scalar> uTauWall_;  // 按 wallFaces_ 顺序，最近一次 correct 的 u_τ
    SolverControls controls_;
    scalar relax_ = 1.0;
    ConvectionScheme divScheme_;
    int lastTimeIndex_ = -1;
    std::map<std::string, scalar> coeffs_;
    // ν_t/ν 上限（Fluent 的 turbulent viscosity ratio limit，默认 1e5）：防止 ε/ω → 0 时 ν_t 溢出；nLimited_ 为被限制的单元数
    scalar nutMaxRatio_ = 1e5;
    label nLimited_ = 0;
    std::map<std::string, glabel> nBounded_;

    // 近壁常数
    scalar kappa_ = 0.41, E_ = 9.8, Cmu_ = 0.09;
    scalar yPlusLam_ = 11.53;  // 线性/对数律交点
};

std::unique_ptr<TurbulenceModel> makeLES(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i,
                                         const std::string& type);
std::unique_ptr<TurbulenceModel> makeRAS(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i,
                                         const std::string& type);

} // namespace cfd
