#pragma once
// 计算网格（每个进程一份局部网格）。
//
// 编号约定：
//   单元 [0, nCells)              自有单元
//        [nCells, nCells+nGhost)  幽灵单元（其他进程的单元副本，由 halo 交换更新）
//   面   [0, nInternalFaces)      内部面：两侧均为单元（含跨进程面、周期面）
//        [nInternalFaces, nFaces) 边界面，按 patch 分组连续存放
//
// 周期（cyclic）面在构建时被合并为普通内部面，owner→neighbour 的距离向量 d()
// 已包含周期平移，因此所有离散算子一律用 d() 而不要用 C()[N]-C()[P]。
// 二维网格的 empty 面不参与计算（不出现在面列表中）。
//
// 局部面顺序保持全局计算面编号的相对顺序，单元上的所有累加顺序与串行相同，
// 配合可复现归约，结果与进程数无关。

#include "fvm/core/Types.h"
#include "fvm/discretization/Schemes.h"
#include "fvm/mesh/RawMesh.h"
#include "fvm/parallel/Comm.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace cfd {

struct Patch {
    std::string name;
    PatchType type = PatchType::Patch;
    label start = 0;   // 局部面编号
    label size = 0;
    glabel globalSize = 0;
    par::GlobalOrdering ordering;  // 面在全局 patch 中的顺序（用于收集表面数据）
    label end() const { return start + size; }
};

struct VtkPiece {
    std::vector<Vec3> points;
    std::vector<std::int64_t> connectivity, offsets;
    std::vector<std::uint8_t> types;
    std::vector<std::int64_t> faces, faceOffsets;  // 仅 VTK_POLYHEDRON 使用（其余单元 faceOffsets=-1）
    bool hasPolyhedra = false;
    // 边界面的点（CSR，按局部边界面顺序 f − nInternalFaces，点号为 points 的下标，外法向绕向）：patch 表面输出用
    std::vector<std::int64_t> bFaceConn, bFaceOff{0};
};

class Mesh {
public:
    // raw 仅在 0 号进程需要（其他进程传 nullptr）
    static std::shared_ptr<Mesh> build(const RawMesh* raw);

    label nCells() const { return nCells_; }
    label nGhost() const { return nGhost_; }
    label nTotalCells() const { return nCells_ + nGhost_; }
    label nInternalFaces() const { return nInternalFaces_; }
    label nFaces() const { return label(owner_.size()); }
    label nBoundaryFaces() const { return nFaces() - nInternalFaces_; }
    glabel nGlobalCells() const { return nGlobalCells_; }

    const std::vector<scalar>& V() const { return V_; }
    const std::vector<Vec3>& C() const { return C_; }
    const std::vector<label>& owner() const { return owner_; }
    const std::vector<label>& neighbour() const { return neighbour_; }
    const std::vector<Vec3>& Sf() const { return Sf_; }
    const std::vector<scalar>& magSf() const { return magSf_; }
    const std::vector<Vec3>& Cf() const { return Cf_; }
    const std::vector<Vec3>& d() const { return d_; }                  // P→N（边界面 P→面心）
    const std::vector<scalar>& w() const { return w_; }                // 线性插值 owner 权重
    const std::vector<scalar>& deltaCoeffs() const { return dc_; }     // 1/|d|
    const std::vector<scalar>& nonOrthDeltaCoeffs() const { return ndc_; }  // 1/max(n·d, 0.05|d|)
    const std::vector<Vec3>& corrVecs() const { return corr_; }        // n − d·ndc（非正交修正向量）
    // 面偏斜度 |x_f − x'_f|/|d|（x'_f 为 PN 连线上按法向距离分割的插值点）：全局平均 / 最大。
    // Gauss 线性梯度在偏斜面上有不随加密减小的零阶误差（内插点不在面心），非结构三角形/四面体网格应用最小二乘梯度。
    scalar meanSkewness() const;
    scalar maxSkewness() const;
    // 内部面偏斜向量 x_f − x'_f（边界面为 0）与"网格存在偏斜面"标志（全局一致）：
    // 面插值的偏斜修正 φ_f = w φ_P + (1−w) φ_N + (x_f − x'_f)·∇φ_f 用它，无偏斜（六面体）网格时跳过
    const std::vector<Vec3>& skewVecs() const { return skew_; }
    bool skewed() const { return skewed_; }
    // 默认梯度格式：fvc::grad(f) 不指定格式时使用；扩散项非正交修正、linearUpwind、湍流模型等都经由它，
    // 由应用按算例 "schemes.grad" 设置（未指定时偏斜网格取最小二乘，见 fvmFlow）
    GradScheme gradScheme() const { return gradScheme_; }
    void setGradScheme(GradScheme s) { gradScheme_ = s; }

    const std::vector<Patch>& patches() const { return patches_; }
    label findPatch(const std::string& name) const;  // 不存在返回 -1
    bool isEmptyPatch(const std::string& name) const;
    const std::vector<std::string>& emptyPatchNames() const { return emptyPatches_; }

    const std::vector<glabel>& cellGlobal() const { return cellGlobal_; }
    const std::vector<glabel>& faceGlobal() const { return faceGlobal_; }  // 全局计算面编号
    glabel nGlobalFaces() const { return nGlobalFaces_; }
    // 面数据收集到 0 号进程（按全局面编号）/ 从 0 号进程分发到所有局部面（含跨进程面的两份副本）
    std::vector<double> gatherFaces(const double* local, int nc) const;
    void scatterFaces(const std::vector<double>& global, double* local, int nc) const;
    const par::Halo& halo() const { return halo_; }
    const par::GlobalOrdering& cellOrdering() const { return cellOrdering_; }
    const VtkPiece& vtk() const { return vtk_; }

    bool twoD() const { return twoD_; }
    int emptyDir() const { return emptyDir_; }
    // 需要求解的向量分量（二维时去掉 empty 方向）
    const std::vector<int>& solvedComponents() const { return solved_; }
    scalar totalVolume() const { return totalVolume_; }
    // 单元所属面（CSR，只针对自有单元）
    const std::vector<label>& cellFaceOffsets() const { return cfOff_; }
    const std::vector<label>& cellFaces() const { return cfList_; }

    void printSummary() const;
    // 单元质量（自有单元）：所含内部面的最大非正交角（度）、最大面偏斜度、面积比 max|S_f|/min|S_f|（二维不含 empty 面）
    std::vector<scalar> cellNonOrthogonality() const;
    std::vector<scalar> cellSkewness() const;
    std::vector<scalar> cellAspectRatio() const;

    // 依附于网格的派生数据缓存（如预条件器面顺序、多重网格层级），随网格一起释放
    template <class T, class Build> const T& cached(const std::string& key, Build&& build) const {
        auto it = cache_.find(key);
        if (it == cache_.end()) it = cache_.emplace(key, std::make_shared<T>(build())).first;
        return *std::static_pointer_cast<T>(it->second);
    }

private:
    mutable std::map<std::string, std::shared_ptr<void>> cache_;
    void finishLocal();

    label nCells_ = 0, nGhost_ = 0, nInternalFaces_ = 0;
    glabel nGlobalCells_ = 0;
    std::vector<scalar> V_;
    std::vector<Vec3> C_;
    std::vector<label> owner_, neighbour_;
    std::vector<Vec3> Sf_, Cf_, d_, corr_;
    std::vector<scalar> magSf_, w_, dc_, ndc_;
    std::vector<Patch> patches_;
    std::vector<std::string> emptyPatches_;
    std::vector<glabel> cellGlobal_, faceGlobal_;
    glabel nGlobalFaces_ = 0;
    std::vector<std::vector<glabel>> faceListsOnMaster_;  // 0 号：各进程局部面的全局编号
    par::Halo halo_;
    par::GlobalOrdering cellOrdering_;
    VtkPiece vtk_;
    bool twoD_ = false;
    int emptyDir_ = -1;
    std::vector<int> solved_{0, 1, 2};
    scalar totalVolume_ = 0;
    std::vector<label> cfOff_, cfList_;
    scalar skewSum_ = 0, skewMax_ = 0;  // 本进程内部面偏斜度之和 / 最大
    std::vector<Vec3> skew_;
    bool skewed_ = false;
    GradScheme gradScheme_ = GradScheme::GaussLinear;
};

using MeshPtr = std::shared_ptr<Mesh>;

} // namespace cfd
