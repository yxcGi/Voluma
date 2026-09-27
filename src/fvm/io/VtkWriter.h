#pragma once
// VTK 输出：每个进程写一个 .vtu 分块（XML，二进制 appended raw），0 号进程写 .pvtu 和 .pvd 时间序列。
// 用 ParaView 打开 <dir>/<base>.pvd 即可。

#include "fvm/field/VolField.h"

#include <functional>
#include <string>
#include <vector>

namespace cfd {

class VtkWriter {
public:
    VtkWriter(MeshPtr mesh, std::string dir, std::string base = "fields");

    // 注册单元数据（引用需在写出时仍有效）
    void add(const VolField<scalar>& f) { add(f.name(), f.internal()); }
    void add(const VolField<Vec3>& f) { add(f.name(), f.internal()); }
    void add(const std::string& name, const std::vector<scalar>& v) {
        items_.push_back({name, [&v] { return v.data(); }, 1});
    }
    void add(const std::string& name, const std::vector<Vec3>& v) {
        items_.push_back({name, [&v] { return reinterpret_cast<const double*>(v.data()); }, 3});
    }
    void add(const std::string& name, const std::vector<Tensor>& v) {
        items_.push_back({name, [&v] { return reinterpret_cast<const double*>(v.data()); }, 9});
    }
    void clear() { items_.clear(); }
    // 写出一帧
    void write(scalar time);

private:
    struct Item {
        std::string name;
        std::function<const double*()> data;
        int nc;
    };
    MeshPtr mesh_;
    std::string dir_, base_;
    std::vector<Item> items_;
    std::vector<std::pair<scalar, std::string>> frames_;
};

// 整场逐位校验和（收集到 0 号按全局顺序哈希），用于比较不同进程数的结果
std::string fieldChecksum(const Mesh& m, const double* data, int nc);
inline std::string checksum(const VolField<scalar>& f) { return fieldChecksum(f.mesh(), f.internal().data(), 1); }
inline std::string checksum(const VolField<Vec3>& f) {
    return fieldChecksum(f.mesh(), reinterpret_cast<const double*>(f.internal().data()), 3);
}

} // namespace cfd
