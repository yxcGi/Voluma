#pragma once
// 原始多面体网格（OpenFOAM polyMesh 拓扑），只在 0 号进程上存在，用于构建计算网格。

#include "fvm/core/Types.h"

#include <functional>
#include <string>
#include <vector>

namespace cfd {

enum class PatchType { Patch, Wall, Symmetry, Empty, Cyclic };

PatchType patchTypeFromString(const std::string& s);
std::string toString(PatchType t);

struct RawPatch {
    std::string name;
    PatchType type = PatchType::Patch;
    glabel start = 0;
    glabel size = 0;
    std::string neighbourPatch;  // 仅 cyclic
};

struct RawMesh {
    std::vector<Vec3> points;
    std::vector<glabel> faceOffsets{0};  // CSR：第 f 个面的点为 facePoints[faceOffsets[f] .. faceOffsets[f+1])
    std::vector<glabel> facePoints;
    std::vector<glabel> owner;      // 每个面
    std::vector<glabel> neighbour;  // 仅内部面
    std::vector<RawPatch> patches;
    glabel nCells = 0;

    glabel nFaces() const { return glabel(owner.size()); }
    glabel nInternalFaces() const { return glabel(neighbour.size()); }
    void addFace(std::initializer_list<glabel> pts) {
        facePoints.insert(facePoints.end(), pts.begin(), pts.end());
        faceOffsets.push_back(glabel(facePoints.size()));
    }
    void addFace(const std::vector<glabel>& pts) {
        facePoints.insert(facePoints.end(), pts.begin(), pts.end());
        faceOffsets.push_back(glabel(facePoints.size()));
    }
};

// 读 OpenFOAM constant/polyMesh（ASCII）
RawMesh readPolyMesh(const std::string& polyMeshDir);
// 写 OpenFOAM polyMesh（ASCII），便于在 ParaView/OpenFOAM 中检查生成的网格
void writePolyMesh(const RawMesh& m, const std::string& polyMeshDir);

// ------------------------------------------------------------ 结构网格生成
// 在 [0,1]^3 参考立方体上生成 nx×ny×nz 六面体网格，再由 map 映射到物理空间。
// 各方向可设周期（生成 cyclic 面对）。nz==1 且 twoD==true 时前后面为 empty（二维算例）。
struct BoxSpec {
    int n[3] = {10, 10, 1};
    Vec3 lo{0, 0, 0}, hi{1, 1, 1};
    bool periodic[3] = {false, false, false};
    bool twoD = false;
    // 可选：各方向一维节点分布 s∈[0,1] → [0,1]（用于壁面加密），为空时均匀
    std::function<scalar(scalar)> stretch[3];
    // 可选：整体坐标映射（在 lo/hi 线性映射之后调用）
    std::function<Vec3(const Vec3&)> map;
    // 边界面名称（xMin xMax yMin yMax zMin zMax），周期方向忽略类型
    std::string names[6] = {"xMin", "xMax", "yMin", "yMax", "zMin", "zMax"};
    PatchType types[6] = {PatchType::Patch, PatchType::Patch, PatchType::Patch,
                          PatchType::Patch, PatchType::Patch, PatchType::Patch};
};

RawMesh generateBox(const BoxSpec& spec);

// 常用一维加密函数：双侧 tanh 聚集（对称），beta 越大越贴壁
std::function<scalar(scalar)> tanhStretch(scalar beta);

} // namespace cfd
