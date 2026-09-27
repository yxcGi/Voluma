#pragma once
// 二维多边形网格 → 单层拉伸的 polyMesh（前后面为 empty），供二维网格生成器与二维网格导入使用。

#include "fvm/mesh/RawMesh.h"

#include <functional>
#include <string>
#include <vector>

namespace cfd {

struct Mesh2D {
    std::vector<Vec3> points;               // 只用 x、y
    std::vector<std::vector<glabel>> cells; // 每个单元的点，逆时针
};

struct BoundaryPatchSpec {
    std::string name;
    PatchType type = PatchType::Patch;
};

// classify(a, b) 返回边界边 (a,b) 所属 patch 的编号（patches 下标）
RawMesh extrude2D(const Mesh2D& m, const std::vector<BoundaryPatchSpec>& patches,
                  const std::function<int(glabel a, glabel b)>& classify, scalar depth);

// 合并距离小于 tol 的重复点（如 C 网格尾迹割线两侧），并更新单元点号
void mergeDuplicatePoints(Mesh2D& m, scalar tol);

} // namespace cfd
