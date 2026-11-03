#pragma once
// 直圆管的蝶形 O 型六面体网格：中央方块 + 4 个环向块（二维截面），沿轴向拉伸。
//
// 坐标：轴向为 z ∈ [0, length]，截面圆心在原点、半径 radius。中央方块半边长 coreFrac·radius，
// 环向块的网格线从方块边界向圆周按角度均匀过渡（corner 处非正交约 45°，是检验非正交修正的典型网格）。
// 边界：wall（圆柱面）；periodic 为 true 时两端为 cyclic 对 zMin/zMax，否则 inlet（z=0）/outlet（z=length）。

#include "fvm/mesh/RawMesh.h"

namespace cfd {

struct PipeMeshSpec {
    scalar radius = 0.5;
    scalar length = 1.0;
    int nCore = 8;          // 中央方块每边单元数
    int nRing = 8;          // 环向块径向单元数
    int nAxial = 16;        // 轴向单元数
    scalar coreFrac = 0.5;  // 中央方块半边长 / 半径
    scalar stretch = 0;     // 环向块径向向壁 tanh 加密参数（0 为均匀）
    bool periodic = true;
};

RawMesh generatePipeMesh(const PipeMeshSpec& s);

} // namespace cfd
