#pragma once
// LDU 稀疏矩阵的标量线性求解器（每个向量分量单独求解）。
//
// 矩阵按面存储：diag[c]，upper[f] = A(owner, neighbour)，lower[f] = A(neighbour, owner)。
// 幽灵单元的行不参与求解，其值每次迭代由 halo 交换得到。

#include "fvm/core/Types.h"
#include "fvm/mesh/Mesh.h"

#include <string>
#include <vector>

namespace cfd {

struct SolverControls {
    std::string solver = "PCG";          // PCG | PBiCGStab | GaussSeidel | symGaussSeidel | Jacobi
    std::string preconditioner = "diagonal";  // diagonal | DIC | DILU | none
    scalar tolerance = 1e-8;
    scalar relTol = 0.0;
    int maxIter = 1000;
    int minIter = 0;
    int nSweeps = 1;  // 平滑器每次残差检查之间的扫描次数
};

struct SolverPerformance {
    std::string solver, field;
    scalar initialResidual = 0, finalResidual = 0;
    int iterations = 0;
    bool converged = false;
};

struct LduSystem {
    const Mesh& mesh;
    const std::vector<scalar>& diag;   // 已含边界贡献
    const scalar* upper;
    const scalar* lower;  // 对称矩阵时与 upper 相同
    const std::vector<scalar>& b;
    bool symmetric;
};

SolverPerformance solveLdu(const LduSystem& sys, std::vector<scalar>& x, const SolverControls& ctrl);

// y = A x（会先对 x 做 halo 交换）
void amul(const LduSystem& sys, std::vector<scalar>& x, std::vector<scalar>& y);

} // namespace cfd
