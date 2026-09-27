#pragma once
// 周期山（ERCOFTAC UFR 3-30，Mellen-Fröhlich-Rodi 2000 / Breuer 等 2009）贴体结构网格。
//
// 坐标与 OpenLB periodichill3d 算例相同：x 流向（周期，Lx = 9H），y 展向（周期，Ly = 4.5H），
// z 壁面法向：下壁为山形 z = h(x)，上壁 z = 3.036H 为平壁。山顶在 x = 0 与 x = Lx。
// 网格线沿 z 竖直，z 方向在两壁按 tanh 加密。边界：bottom、top（wall），x、y 周期。

#include "fvm/mesh/RawMesh.h"

namespace cfd {

// ERCOFTAC 山形（山高 H，周期长度 Lx）：返回 x 处的下壁高度
scalar periodicHillHeight(scalar x, scalar H = 1.0, scalar Lx = 9.0);

struct PeriodicHillSpec {
    int n[3] = {128, 64, 64};  // 流向、展向、法向单元数
    scalar H = 1.0;            // 山高
    scalar length = 9.0;       // 流向周期长度（H 的倍数）
    scalar span = 4.5;         // 展向宽度（H 的倍数）
    scalar height = 3.036;     // 上壁高度（H 的倍数）
    scalar stretch = 2.0;      // 法向 tanh 加密参数（0 为均匀）
};

RawMesh generatePeriodicHillMesh(const PeriodicHillSpec& s);

} // namespace cfd
