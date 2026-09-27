#pragma once
// 湍流模型（不可压，运动学量）。所有模型给 IncompressibleFlow 提供单元 ν_t 以及壁面面上的 ν_eff，
// 在每个外迭代开始前（使用最新的 U、φ）更新一次；输运方程的旧时间层每个时间步保存一次。
//
//   LES   Smagorinsky（可选 Van Driest 阻尼）、WALE；可选壁面模型（WMLES）：
//         spalding / musker / powerLaw（Werner-Wengle）/ logLaw，采样高度可设，给出壁面切应力
//   RANS  kEpsilon（标准）、realizableKE（Shih）、kOmega（Wilcox 1998）、kOmegaSST（Menter 2003）、
//         SpalartAllmaras（无 ft2）
//         壁函数：nutkWallFunction（ν_t 壁面值）、epsilon/omega 壁面单元值（omega 为粘性/对数层混合），
//         近壁产生项 G 按对数律修正；wallFunctions=false 时为低 Re 处理（壁面 ν_t = 0）
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
    // 打印壁面 y+ 统计
    void printYPlus(std::ostream& os) const;
    // 续算：模型场及旧时间层（与 restart::write/read 同目录）
    void writeRestart(const std::string& dir);
    void readRestart(const std::string& dir);

protected:
    virtual void correctModel() = 0;
    // 设置字段的边界条件：非壁面按 JSON，壁面按模型要求（wallDefault）
    // 创建并登记模型场；wallDefault 为 "zeroGradient" 或 "zero"（定值 0）
    ScalarField& addField(const std::string& name, const std::string& wallDefault, scalar initDefault);
    // 输运方程：ddt(ψ) + div(φ,ψ) − ∇·(Γ∇ψ) + Sp·ψ = Su
    // gammaCell：单元扩散系数（插值到面），Sp、Su：单位体积，Sp ≥ 0
    void solveTransport(ScalarField& psi, const std::vector<scalar>& gammaCell, const std::vector<scalar>& Sp,
                        const std::vector<scalar>& Su, const std::vector<label>& fixCells = {},
                        const std::vector<scalar>& fixValues = {});
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

    // 近壁常数
    scalar kappa_ = 0.41, E_ = 9.8, Cmu_ = 0.09;
    scalar yPlusLam_ = 11.53;  // 线性/对数律交点
};

std::unique_ptr<TurbulenceModel> makeLES(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i,
                                         const std::string& type);
std::unique_ptr<TurbulenceModel> makeRAS(IncompressibleFlow& flow, const Json& d, const Json& b, const Json& i,
                                         const std::string& type);

} // namespace cfd
