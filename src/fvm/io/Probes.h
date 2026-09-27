#pragma once
// 点探针：记录若干空间点所在单元的场值时间序列（CSV，每个场一个文件）。
// 取点所在（最近体心）单元的单元值；定位在构造时完成，结果与进程数无关。

#include "fvm/field/VolField.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace cfd {

class Probes {
public:
    Probes(MeshPtr mesh, std::vector<Vec3> points, std::string dir, std::string name = "probes");

    void add(const VolField<scalar>& f) { add(f.name(), f.internal(), 1); }
    void add(const VolField<Vec3>& f) { add(f.name(), f.internal(), 3); }
    template <class T> void add(const std::string& name, const std::vector<T>& v, int nc) {
        items_.push_back({name, [&v] { return reinterpret_cast<const double*>(v.data()); }, nc, nullptr});
    }
    void sample(scalar time);
    // 当前各探针的值（所有进程相同），nc 为分量数
    std::vector<double> values(const double* data, int nc) const;
    template <class T> std::vector<double> values(const VolField<T>& f) const {
        return values(reinterpret_cast<const double*>(f.internal().data()), Traits<T>::nComponents);
    }
    // 把当前值写成一张表（每行一个点：x y z 值...），常用于沿线采样
    template <class T> void writeTable(const std::string& file, const VolField<T>& f) const {
        writeTable(file, reinterpret_cast<const double*>(f.internal().data()), Traits<T>::nComponents, f.name());
    }
    void writeTable(const std::string& file, const double* data, int nc, const std::string& name) const;
    const std::vector<Vec3>& points() const { return points_; }
    // 实际采样位置（所在单元体心）
    const std::vector<Vec3>& cellCentres() const { return centres_; }

private:
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
};

} // namespace cfd
