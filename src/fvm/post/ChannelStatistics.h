#pragma once
// 槽道（两个均匀方向）湍流统计：沿壁面法向逐层做平面 + 时间平均，
// 得到平均速度、雷诺应力（解析部分与亚格子部分）、压力脉动均方根和摩擦速度。
//
// 输出格式与 OpenLB channel3d 算例相同（分号分隔），便于直接对比：
//   y+;u_tau;uAv+;uu++;uv++;uw++;vv++;vw++;ww++;pRMS;uu_sgs++;uv_sgs++;uw_sgs++;vv_sgs++;vw_sgs++;ww_sgs++
// 其中 u、v、w 依次为流向、展向、壁面法向分量（OpenLB 约定），++ 表示除以 u_τ²；
// 上下两半按到最近壁面的距离合并（上半的法向分量反号镜像）。y+ = y u_τ/ν，y 为单元中心到壁面距离。

#include "fvm/solvers/IncompressibleFlow.h"

#include <string>
#include <vector>

namespace cfd {

class ChannelStatistics {
public:
    // streamwise / normal：流向与壁面法向坐标轴（0,1,2），展向为剩下的轴；
    // yWall0 / yWall1：两壁面在法向的坐标
    ChannelStatistics(const IncompressibleFlow& flow, int streamwise, int normal, scalar yWall0, scalar yWall1);

    // 累积一次样本（权重 dt）
    void sample(scalar dt);
    // 写出当前平均剖面
    void write(const std::string& file) const;
    scalar averagedTime() const { return time_; }
    scalar uTau() const;

private:
    const IncompressibleFlow& flow_;
    int sx_, ny_, sz_;  // 流向、法向、展向轴
    std::vector<scalar> layerY_;     // 各层到壁面距离
    std::vector<label> cellLayer_;   // 本进程单元 → 层
    std::vector<char> upper_;        // 单元更靠近 yWall1（统计时镜像）
    // 每层累加量（体积 × 时间加权）：vol, U(3), UU(6), p, pp, sgs(6)
    static constexpr int NQ = 18;
    std::vector<double> acc_;
    double time_ = 0;
    double tauAcc_ = 0;  // 时间加权的平均壁面切应力
};

} // namespace cfd
