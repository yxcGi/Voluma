#pragma once
// 聚合型代数多重网格（GAMG）预条件器，供 PCG 使用（压力方程）。
//
// - 每个进程独立聚合（进程间耦合只在最细层由外层 PCG 的 amul 处理），即分块 V 循环；
// - 按几何面权重 |S|/|d| 两两配对聚合（与矩阵系数无关，因此聚合层级按网格缓存）；
// - 分段常数插值，Galerkin 粗化；前光顺正向 GS、后光顺反向 GS，保证预条件器对称；
// - 最粗层稠密 LU 直接求解。

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

    const GamgHierarchy& H_;
    std::vector<Coeffs> A_;
    // 最粗层 LU（行主序，无主元：对角占优/定号矩阵）
    std::vector<scalar> lu_;
    label nCoarse_ = 0;
    mutable std::vector<std::vector<scalar>> rbuf_, xbuf_, tmp_;
};

} // namespace cfd
