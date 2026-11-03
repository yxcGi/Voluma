#pragma once
// 由 JSON 算例文件构建网格、边界条件、求解器参数（各求解器程序共用）。
//
// 网格：
//   "mesh": { "polyMesh": "polyMesh" }                       OpenFOAM 格式（相对算例文件目录）
//   "mesh": { "box": { "n": [64,64,1], "lo": [0,0,0], "hi": [1,1,0.1],
//                      "periodic": [false,false,false], "twoD": true,
//                      "names": {"yMax": "lid"}, "walls": ["xMin", ...],
//                      "stretch": {"y": 2.0} } }              tanh 双侧加密（beta）
//
// 边界条件（每个 patch 每个场一项，详见 docs/case-format.md）：
//   定值类   {"type": "fixedValue", "value": [1,0,0]}  值也可用 "timeTable" / "profile" / "parabolic" 给出，
//            可加 "ramp"（缓启动）、"noise"（随机扰动）；{"type": "noSlip"}、"movingWall"、"rotatingWall"、"flowRateInlet"
//   梯度类   {"type": "zeroGradient"}  {"type": "fixedGradient", "gradient": 0}  {"type": "fixedFluxPressure"}
//   混合类   {"type": "inletOutlet", "inletValue": 0}  "freestream"  "freestreamPressure"  "advective"
//   压力进出口 {"type": "totalPressure", "p0": 0} 配 U 的 {"type": "pressureInletOutletVelocity"}
//   {"type": "robin", "a": 1, "b": 0.1, "c": 0}  即 a φ + b ∂φ/∂n = c；{"type": "symmetry"} / {"type": "slip"}
// wall 类型的 patch 未给出时默认 U=noSlip、p=zeroGradient；其余 patch 必须给出。
// 需要流动信息的条件（inletOutlet、totalPressure、advective 等）通过 BCContext 取通量、速度与时间步。

#include "fvm/core/Json.h"
#include "fvm/field/VolField.h"
#include "fvm/linalg/LinearSolver.h"

#include <string>

namespace cfd {

// 算例文件所在目录（用于解析相对路径）
std::string caseDirectory(const std::string& caseFile);
std::string resolvePath(const std::string& caseDir, const std::string& p);

MeshPtr buildMesh(const Json& meshDict, const std::string& caseDir);

template <class T>
void setBoundaryCondition(VolField<T>& f, const std::string& patch, const Json& d, const BCContext* ctx = nullptr);

// 为场 f 按 "boundary" 字典设置所有 patch 的边界；defaultWall 用于 wall patch 缺省
template <class T>
void setBoundaryConditions(VolField<T>& f, const Json& boundary, const Json& defaultWall, const BCContext* ctx = nullptr);

SolverControls parseSolverControls(const Json& d, SolverControls def);

} // namespace cfd
