#pragma once
// 聚合型代数多重网格（GAMG）预条件器，供 PCG 使用（压力方程）。
//
// - 每个进程独立聚合（聚合不跨进程），但各层都保留进程间耦合（接口），光顺与残差计算前交换 ghost 值；
//   各进程层数相同（按全局单元数决定是否继续粗化），以便每层做集合通信；
// - 按几何面权重 |S|/|d| 两两配对聚合（与矩阵系数无关，因此聚合层级按网格缓存）；
// - 分段常数插值，Galerkin 粗化；前光顺正向 GS、后光顺反向 GS，保证预条件器对称；
// - 最粗层：所有进程的最粗聚合连同进程间耦合组成一个全局小矩阵（≤ 约 256 行），
//   每个进程冗余地做稠密 LU，提供全局粗网格校正（并行时收敛速度基本不随进程数下降）。

#include "fvm/linalg/LinearSolver.h"

#include <memory>
#include <vector>

namespace cfd {

struct GamgLevel {
    label n = 0;
    std::vector<label> l, u;             // 本层面（l<u）
    std::vector<scalar> weight;          // 聚合用几何权重
    std::vector<label> off, adjCell, adjFace;
    std::vector<char> adjIsL;            // 行 i 是该面的 l 端（A(i,j)=upper）
    std::vector<label> agg;              // 本层单元 → 下一层单元
    std::vector<label> faceMap;          // 本层面 → 下一层面（−1：聚合内部）
    std::vector<char> faceFlip;          // 下一层面方向是否与本层相反
};

struct GamgHierarchy {
    std::vector<GamgLevel> levels;
    std::vector<label> meshFace;  // 第 0 层面对应的网格面
    std::vector<char> ownerIsL;   // 网格面 owner 是否为 l 端
    // ---- 并行：进程间耦合
    // 全局最粗层：本进程最粗单元的全局编号偏移、全局总数
    label coarseOffset = 0, nGlobalCoarse = 0;
    struct ProcFace {
        label face;          // 网格面
        label localCell;     // 第 0 层本进程单元
        label haloCell;      // 第 0 层对方单元（ghost 单元编号 ≥ nCells）
        label remoteCoarse;  // 对方单元所属最粗聚合的全局编号
        char localIsOwner;
    };
    std::vector<ProcFace> procFaces;
    std::vector<label> cellToCoarse;  // 第 0 层单元 → 本进程最粗聚合
    // 每层的进程间接口（同 OpenFOAM processorGAMGInterface）：本层向量后接 nGhost 个对方单元值，
    // 条目 e 表示 A(cell[e], ghost[e]) 的耦合，系数为映射到该条目的所有第 0 层进程间面系数之和
    struct Interface {
        label nGhost = 0;
        std::vector<label> cell, ghost;
        std::vector<label> fromProcFace;  // procFaces 下标 → 本层条目
        std::vector<label> off, idx;      // 单元 → 条目（CSR）
        par::Halo halo;                        // 第 1 层起：交换本层 ghost 值（第 0 层用网格 halo）
    };
    std::vector<Interface> iface;
};

GamgHierarchy buildHierarchy(const Mesh& m);
const GamgHierarchy& gamgHierarchy(const Mesh& m);  // 缓存在网格上

class GamgPreconditioner {
public:
    explicit GamgPreconditioner(const LduSystem& s);
    void apply(const std::vector<scalar>& r, std::vector<scalar>& w) const;

private:
    struct Coeffs {
        std::vector<scalar> diag, upper, lower;
    };
    void vcycle(int lev, const std::vector<scalar>& r, std::vector<scalar>& x) const;
    void smooth(int lev, const std::vector<scalar>& r, std::vector<scalar>& x, bool forward) const;
    void exchange(int lev, std::vector<scalar>& x) const;
    // y = A x（含进程间耦合，x 的 ghost 值须已交换）
    void residual(int lev, const std::vector<scalar>& r, const std::vector<scalar>& x, std::vector<scalar>& res) const;

    const Mesh& mesh_;
    const GamgHierarchy& H_;
    std::vector<Coeffs> A_;
    std::vector<std::vector<scalar>> ifc_;  // 各层进程间接口系数
    bool global_ = false;  // 最粗层为全局矩阵
    // 最粗层 LU（行主序，无主元：对角占优/定号矩阵）
    std::vector<scalar> lu_;
    label nCoarse_ = 0;
    mutable std::vector<std::vector<scalar>> rbuf_, xbuf_, tmp_, ebuf_, abuf_;
    mutable std::vector<scalar> coarse_;
};

} // namespace cfd
