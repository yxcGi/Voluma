#pragma once
// 外部网格格式导入（只在 0 号进程调用），统一转成 RawMesh（OpenFOAM polyMesh 拓扑）。
//
//   readGmsh    Gmsh .msh（ASCII 2.2 / 4.1）：四面体、六面体、三棱柱、金字塔及混合网格；
//               高阶单元只取角点。纯二维网格（三角形/四边形）自动沿 z 拉伸一层（前后为 empty）。
//               边界 patch 取自 Physical 组（三维为 Physical Surface，二维为 Physical Curve）。
//   readFluent  Fluent .msh / .cas 网格段（ASCII 与二进制段均可），三维任意多面体或二维（自动拉伸）；
//               zone 名称与类型取自 (39/45 ...) 段，wall → wall，symmetry → symmetry，其余 → patch。
//   readPolyMesh（RawMesh.h）同时支持 ASCII 与二进制 polyMesh。
//
// 所有格式都可用 patchTypes 覆盖 patch 类型（如把 "frontAndBack" 设为 empty）。

#include "fvm/mesh/RawMesh.h"

#include <map>
#include <string>

namespace cfd {

struct MeshImportOptions {
    scalar scale = 1.0;                             // 坐标缩放（如 mm → m 取 0.001）
    scalar depth2D = 0.1;                           // 二维网格拉伸厚度
    std::map<std::string, PatchType> patchTypes;    // 按名称覆盖 patch 类型
};

RawMesh readGmsh(const std::string& file, const MeshImportOptions& opt = {});
RawMesh readFluent(const std::string& file, const MeshImportOptions& opt = {});

// 按扩展名选择：.msh 先看内容区分 Gmsh（$MeshFormat）与 Fluent；目录则按 polyMesh 读
RawMesh readMeshFile(const std::string& path, const MeshImportOptions& opt = {});

// ------------------------------------------------------------------ 公共装配
// 由"面 + 两侧单元"装配 RawMesh：自动按几何把面法向调整为 owner → neighbour（边界面向外），
// 内部面按 (owner, neighbour) 排序，边界面按 patch 分组（组内保持输入顺序），空 patch 去掉。
struct FaceRecord {
    std::vector<glabel> nodes;
    glabel c0 = -1, c1 = -1;  // c1 < 0 为边界面
    int patch = -1;           // 边界面所属 patch（patches 下标）
};
RawMesh assembleFaceMesh(std::vector<Vec3> points, glabel nCells, std::vector<FaceRecord> faces,
                         const std::vector<RawPatch>& patches);

} // namespace cfd
