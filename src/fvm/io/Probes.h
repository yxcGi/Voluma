#pragma once
// 点探针：记录若干空间点所在单元的场值时间序列（CSV，每个场一个文件）。
// 取点所在（最近体心）单元的单元值；定位在构造时完成，结果与进程数无关。
// 沿线采样表（writeTable）可选用单元梯度线性重构到采样点：φ(x) = φ_c + ∇φ_c·(x − x_c)，
// 非结构网格上采样点一般不在体心，重构后表中坐标为采样点本身。

#include "fvm/discretization/Fvc.h"
#include "fvm/field/VolField.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace cfd {

class Probes {
public:
    Probes(MeshPtr mesh, std::vector<Vec3> points, std::string dir, std::string name = "probes");

    // 采样表取值方式：false 取所在单元值（默认）；true 用单元梯度线性重构到采样点（梯度格式默认取网格的默认格式）
    void setInterpolate(bool on) {
        interpolate_ = on;
        gradScheme_ = mesh_->gradScheme();
    }
    void setInterpolate(bool on, GradScheme scheme) {
        interpolate_ = on;
        gradScheme_ = scheme;
    }
    bool interpolate() const { return interpolate_; }

    void add(const VolField<scalar>& f) { add(f.name(), f.internal(), 1); }
    void add(const VolField<Vec3>& f) { add(f.name(), f.internal(), 3); }
    template <class T> void add(const std::string& name, const std::vector<T>& v, int nc) {
        items_.push_back({name, [&v] { return reinterpret_cast<const double*>(v.data()); }, nc, nullptr});
    }
    void sample(scalar time);
    // 当前各探针的值（所有进程相同），nc 为分量数
    std::vector<double> values(const double* data, int nc) const;
    template <class T> std::vector<double> values(const VolField<T>& f) const {
        constexpr int nc = Traits<T>::nComponents;
        if (!interpolate_) return values(reinterpret_cast<const double*>(f.internal().data()), nc);
        const auto g = fvc::grad(f, gradScheme_);
        const Mesh& m = *mesh_;
        std::vector<double> v(points_.size() * nc, 0.0);
        for (std::size_t k = 0; k < points_.size(); ++k) {
            const label c = cell_[k];
            if (c < 0) continue;
            const T val = f[c] + gradDot(points_[k] - m.C()[c], g[c]);
            for (int i = 0; i < nc; ++i) v[k * nc + i] = Traits<T>::component(val, i);
        }
        par::allSumInPlace(v.data(), int(v.size()));
        return v;
    }
    // 把当前值写成一张表（每行一个点：x y z 值...），常用于沿线采样；坐标为 locations()
    template <class T> void writeTable(const std::string& file, const VolField<T>& f) const {
        writeTable(file, values(f), Traits<T>::nComponents, f.name());
    }
    void writeTable(const std::string& file, const double* data, int nc, const std::string& name) const {
        writeTable(file, values(data, nc), nc, name);
    }
    const std::vector<Vec3>& points() const { return points_; }
    // 实际采样位置（所在单元体心）
    const std::vector<Vec3>& cellCentres() const { return centres_; }
    // 采样表中的坐标：重构时为采样点本身，否则为所在单元体心
    const std::vector<Vec3>& locations() const { return interpolate_ ? points_ : centres_; }

private:
    void writeTable(const std::string& file, const std::vector<double>& v, int nc, const std::string& name) const;
    struct Item {
        std::string name;
        std::function<const double*()> data;
        int nc;
        std::shared_ptr<std::ofstream> os;
    };
    MeshPtr mesh_;
    std::vector<Vec3> points_, centres_;
    std::vector<label> cell_;  // 本进程拥有的探针：局部单元号，否则 −1
    std::string dir_, name_;
    std::vector<Item> items_;
    bool interpolate_ = false;
    GradScheme gradScheme_ = GradScheme::LeastSquares;
};

} // namespace cfd
