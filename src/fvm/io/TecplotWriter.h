#pragma once
// Tecplot ASCII 输出（有限元区、单元中心数据，与旧程序的 .dat 格式相同）：
//   二维：FEQUADRILATERAL（三角形写成退化四边形），坐标 X、Y；
//   三维：FEBRICK（楔形、四面体、金字塔写成退化六面体）。
// 数据收集到 0 号进程写一个文件。多面体单元（>8 点）三维时跳过并给出警告。

#include "fvm/field/VolField.h"

#include <functional>
#include <string>
#include <vector>

namespace cfd {

class TecplotWriter {
public:
    explicit TecplotWriter(MeshPtr mesh) : mesh_(std::move(mesh)) {}

    void add(const VolField<scalar>& f) { add(f.name(), f.internal(), 1); }
    void add(const VolField<Vec3>& f) { add(f.name(), f.internal(), 3); }
    template <class T> void add(const std::string& name, const std::vector<T>& v, int nc) {
        items_.push_back({name, [&v] { return reinterpret_cast<const double*>(v.data()); }, nc});
    }
    void write(const std::string& file, const std::string& title = "fvm") const;

private:
    struct Item {
        std::string name;
        std::function<const double*()> data;
        int nc;
    };
    MeshPtr mesh_;
    std::vector<Item> items_;
};

} // namespace cfd
