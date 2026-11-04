#pragma once
// 场的时间平均（瞬态按 dt 加权，稳态按迭代计数）：<φ> 以及二阶矩 <φ'φ'> = <φφ> − <φ><φ>
// （与 OpenFOAM fieldAverage 的 Mean / Prime2Mean 同义；向量场的二阶矩为 3×3 对称张量，含雷诺应力）。
// 平均量命名 <name>Mean、<name>Prime2Mean，长度 nTotalCells（幽灵单元由 halo 交换保持一致），
// 可登记到 VtkWriter / TecplotWriter 与采样表；累积量写入续算目录，换进程数续算后继续累积。

#include "fvm/io/TecplotWriter.h"
#include "fvm/io/VtkWriter.h"
#include "fvm/mesh/Mesh.h"

#include <functional>
#include <string>
#include <vector>

namespace cfd {

class FieldAverage {
public:
    explicit FieldAverage(MeshPtr mesh) : mesh_(std::move(mesh)) {}

    // 登记要平均的单元数组（引用需在采样时仍有效；长度 ≥ nCells）
    void add(const std::string& name, const std::vector<scalar>& src, bool prime2 = false);
    void add(const std::string& name, const std::vector<Vec3>& src, bool prime2 = false);

    // 累积一次样本（权重 weight：瞬态为 dt，稳态为 1）
    void sample(scalar weight);
    scalar averagedTime() const { return time_; }
    int nSamples() const { return n_; }
    void reset();

    struct Item {
        std::string name;
        int nc = 1;
        bool prime2 = false;
        std::function<const double*()> src;
        std::vector<double> sum, sum2;   // 自有单元：nCells·nc、nCells·nc2
        std::vector<scalar> meanS, p2S;  // nc == 1（含幽灵）
        std::vector<Vec3> meanV;         // nc == 3
        std::vector<Tensor> p2V;         // nc == 3（对称）
        int nc2() const { return nc == 1 ? 1 : 9; }
    };
    const std::vector<Item>& items() const { return items_; }

    // 登记平均场到写出器
    void addTo(VtkWriter& w) const;
    void addTo(TecplotWriter& w) const;
    // 当前平均量视图（名字、数据、分量数），供采样表使用
    struct View {
        std::string name;
        const double* data;
        int nc;
    };
    std::vector<View> views() const;

    // 续算：累积量与累计时间（dir 与 restart::write 相同）；readRestart 文件不存在时返回 false（从零开始）
    void writeRestart(const std::string& dir) const;
    bool readRestart(const std::string& dir);

private:
    void update();
    MeshPtr mesh_;
    std::vector<Item> items_;
    double time_ = 0;
    int n_ = 0;
};

} // namespace cfd
