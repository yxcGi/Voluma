#pragma once
// 单元到最近壁面（wall 类型 patch）的距离，供 SST、Spalart-Allmaras、Van Driest 阻尼等使用。
//
// 所有壁面面（中心、单位法向、尺度）汇总到每个进程，用 k-d 树找最近的若干个面心，
// 再取到面的近似距离：法向投影距离，若投影点落在面外则加上面外的切向偏移。
// 结果与进程数无关，缓存在网格上。

#include "fvm/mesh/Mesh.h"

#include <vector>

namespace cfd {

struct WallDistance {
    std::vector<scalar> y;   // 长度 nTotalCells（含幽灵单元）；无壁面时为 GREAT
    std::vector<Vec3> n;     // 最近壁面的单位法向（指向流体）
    std::vector<glabel> nearest;  // 最近壁面面在全局壁面面列表中的序号（见 gatherWallFaces）
    bool hasWalls = false;
};

const WallDistance& wallDistance(const Mesh& m);

// 全局壁面面列表：各 wall patch 依次、patch 内按全局面序。values 按本进程边界面（长度 nBoundaryFaces）
// 给出，返回全局列表上的值（所有进程相同）。
std::vector<scalar> gatherWallFaces(const Mesh& m, const std::vector<scalar>& values);
// 本进程边界面 → 全局壁面面序号（非壁面为 −1），长度 nBoundaryFaces
const std::vector<glabel>& wallFaceGlobalIndex(const Mesh& m);

} // namespace cfd
