#pragma once
// 平面采样：取被平面 (x − x0)·n = 0 切过的单元（按单元顶点到平面的有符号距离的符号变化判定，
// 平面恰好落在网格面上时只取一侧的一层），在单元中心到平面的投影点上给出场值：VolField 用单元梯度线性重构
// φ(x_s) = φ_c + ∇φ_c·(x_s − x_c)，普通单元数组取单元值。汇总到 0 号进程、按全局单元编号排序写 CSV：
//   x,y,z,<name>[_x,_y,_z]
// 结果与进程数无关；ParaView 中可用 Table To Points + Delaunay2D 显示。

#include "fvm/discretization/Fvc.h"
#include "fvm/field/VolField.h"

#include <string>
#include <vector>

namespace cfd {

class CutPlane {
public:
    CutPlane(MeshPtr mesh, const Vec3& point, const Vec3& normal);

    label nLocalCells() const { return label(cells_.size()); }
    glabel nCells() const;  // 全局被切单元数
    const std::vector<label>& cells() const { return cells_; }
    const std::vector<Vec3>& points() const { return points_; }  // 采样点（单元中心在平面上的投影）

    template <class T> void writeTable(const std::string& file, const VolField<T>& f) const {
        constexpr int nc = Traits<T>::nComponents;
        const auto g = fvc::grad(f, mesh_->gradScheme());
        std::vector<double> v(cells_.size() * nc);
        for (std::size_t k = 0; k < cells_.size(); ++k) {
            const label c = cells_[k];
            const T val = f[c] + gradDot(points_[k] - mesh_->C()[c], g[c]);
            for (int i = 0; i < nc; ++i) v[k * nc + i] = Traits<T>::component(val, i);
        }
        writeTable(file, v, nc, f.name());
    }
    // 单元数组（取单元值，不重构）
    void writeTable(const std::string& file, const double* cellData, int nc, const std::string& name) const;

private:
    void writeTable(const std::string& file, const std::vector<double>& v, int nc, const std::string& name) const;
    MeshPtr mesh_;
    Vec3 x0_, n_;
    std::vector<label> cells_;
    std::vector<Vec3> points_;
};

} // namespace cfd
