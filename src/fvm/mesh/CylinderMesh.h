#pragma once
// 圆柱绕流 O 型结构网格：圆柱壁面到矩形外边界 [−upstream, downstream] × [−halfHeight, halfHeight]，
// 沿从圆心出发的射线布点（壁面处正交），径向几何增长、首层高度 firstCell；展向 nSpan 层、周期。
// nSpan = 1 且 twoD 时为二维（前后 empty）。单位：直径 D，圆心在原点，来流沿 +x。
// 边界：cylinder（wall）、inlet（x = −upstream）、outlet（x = downstream）、top、bottom，
// 展向 front/back（cyclic）或 frontAndBack（empty）。

#include "fvm/mesh/RawMesh.h"

namespace cfd {

struct CylinderMeshSpec {
    scalar D = 1.0;
    scalar upstream = 10, downstream = 20, halfHeight = 10;  // 以 D 计
    int nTheta = 192;       // 周向单元数
    int nRadial = 96;       // 径向单元数
    int nSpan = 32;         // 展向单元数
    scalar span = 3.141592653589793;  // 展向长度（以 D 计）
    scalar firstCell = 0.004;         // 壁面首层高度（以 D 计）
    bool twoD = false;
};

RawMesh generateCylinderOMesh(const CylinderMeshSpec& s);

// 把边界面组 name 按面心分成若干组（classify 返回 parts 下标），其余面组不变
void splitPatch(RawMesh& m, const std::string& name, const std::vector<std::pair<std::string, PatchType>>& parts,
                const std::function<int(const Vec3& faceCentre)>& classify);

} // namespace cfd
