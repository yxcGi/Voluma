#pragma once
// patch 表面的 VTK 输出（PolyData）：每个进程写一个 .vtp 分块，0 号进程写 .pvtp 和 .pvd 时间序列，
// 用 ParaView 打开 <dir>/<base>_<patch>.pvd。数据按边界面编号（f − nInternalFaces）给出，长度 nBoundaryFaces，
// 只写所选 patch 的面（面的顶点来自网格的边界面点表，外法向绕向）。

#include "fvm/mesh/Mesh.h"

#include <functional>
#include <map>
#include <string>
#include <vector>

namespace cfd {

class SurfaceWriter {
public:
    SurfaceWriter(MeshPtr mesh, std::string dir, std::string base = "surface");

    void addPatch(label p) { patches_.push_back(p); }
    const std::vector<label>& patches() const { return patches_; }
    // 登记边界面数据（引用需在写出时仍有效；长度 nBoundaryFaces）
    void add(const std::string& name, const std::vector<scalar>& bf) {
        items_.push_back({name, [&bf] { return bf.data(); }, 1});
    }
    void add(const std::string& name, const std::vector<Vec3>& bf) {
        items_.push_back({name, [&bf] { return reinterpret_cast<const double*>(bf.data()); }, 3});
    }
    // 写出一帧（所有登记的 patch）
    void write(scalar time);

private:
    struct Item {
        std::string name;
        std::function<const double*()> data;
        int nc;
    };
    MeshPtr mesh_;
    std::string dir_, base_;
    std::vector<label> patches_;
    std::vector<Item> items_;
    std::map<label, std::vector<std::pair<scalar, std::string>>> frames_;
};

} // namespace cfd
