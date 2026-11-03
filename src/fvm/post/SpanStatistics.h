#pragma once
// 展向均匀流动的统计（周期山等）：按 (流向坐标, 法向坐标) 相同的单元列做展向 + 时间平均，
// 并对一个壁面按流向坐标统计壁面切应力、压力。输出与 OpenLB periodichill3d 算例同格式（分号分隔）：
//   profile_xh_<站>.csv：y_h;u_Ub;v_Ub;uu_Ub2;vv_Ub2;ww_Ub2;uv_Ub2;k_Ub2
//     （y_h 为离下壁高度 / H；v 为法向分量，ww 为展向，uv 为流向-法向剪应力；只含解析部分）
//   wall_<patch>.csv：x_h;Cf;Cp;ut_Ub;y_h（按壁面面心的流向、法向坐标分组）
//     （Cf = 2 τ_w·t / U_b²，t 为壁面切向且指向 +流向；Cp 以流向坐标最小处的壁压为参考（或 pRef）；
//       ut 为壁面第一层单元的切向速度；y_h 为面心法向坐标）
//   fields.csv：x_h;y_h;u_Ub;v_Ub;w_Ub;uu_Ub2;vv_Ub2;ww_Ub2;uv_Ub2;Cp（全部单元列，展向 + 时间平均）
// 并按 Cf 变号给出分离点、再附点（x/H > 0.1 起扫描）。
// 要求网格在展向上是拉伸得到的（同一列单元的流向、法向中心坐标相同），如内置 periodicHill 网格。

#include "fvm/solvers/IncompressibleFlow.h"

#include <string>
#include <vector>

namespace cfd {

class SpanStatistics {
public:
    SpanStatistics(const IncompressibleFlow& flow, int streamwise, int normal, const std::string& wallPatch,
                   scalar Uref, scalar H);

    void sample(scalar dt);
    // 写出各站位剖面与壁面分布（stations 为 x/H）；返回 {分离点, 再附点}（x/H，未找到为 −1）
    std::pair<scalar, scalar> write(const std::string& dir, const std::vector<scalar>& stations) const;
    scalar averagedTime() const { return time_; }
    // Cp 的参考压力（缺省为流向坐标最小的壁面组的压力）
    void setPRef(scalar p) {
        pRef_ = p;
        hasPRef_ = true;
    }

private:
    const IncompressibleFlow& flow_;
    int sx_, nz_, sy_;
    std::string wallName_;
    label wallPatch_;
    scalar Uref_, H_;
    // 单元列（展向平均的“格子”）：全局一致编号
    std::vector<double> binX_, binZ_;
    std::vector<label> cellBin_;
    static constexpr int NQ = 11;  // w, U(3), UU(6), p
    std::vector<double> acc_;
    // 壁面按流向坐标分组
    std::vector<double> wX_, wZ_;
    std::vector<label> faceBin_;   // 壁面面（patch 内序号）→ 组
    std::vector<Vec3> faceT_;      // 切向单位向量
    static constexpr int NW = 4;   // w, τ·t, p, U_P·t
    std::vector<double> wacc_;
    double time_ = 0;
    scalar pRef_ = 0;
    bool hasPRef_ = false;
};

} // namespace cfd
