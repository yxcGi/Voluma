#pragma once
// NACA 四位数翼型的二维 C 型结构网格（矩形计算域，单层拉伸，前后 empty）。
//
// 坐标：弦长 c，前缘在 (0,0)；攻角通过把翼型绕半弦点 (c/2, 0) 旋转 −α 实现，来流保持 +x 方向
// （与 OpenLB airfoil2d 算例相同）。计算域 x ∈ [c/2 − upstream, c/2 + downstream]，y ∈ [−halfHeight, halfHeight]。
//
// 网格：i 方向依次为下尾迹（出口→后缘）、下表面（后缘→前缘）、上表面（前缘→后缘）、上尾迹（后缘→出口），
// j 方向由壁面（或尾迹割线）指向外边界。表面点前缘按余弦加密、后缘适度加密；法向按 tanh 拉伸，
// 首层高度 firstCell；网格线自壁面正交出发（Hermite 插值），再做若干次拉普拉斯光顺。
// 边界：airfoil（wall）、inlet、outlet、top、bottom、frontAndBack（empty）。

#include "fvm/mesh/RawMesh.h"

#include <string>

namespace cfd {

struct AirfoilMeshSpec {
    std::string naca = "0012";
    bool closedTE = true;      // true：尾缘闭合（系数 −0.1036），false：标准 −0.1015
    scalar chord = 1.0;
    scalar alphaDeg = 0.0;
    scalar upstream = 6.5;     // 半弦点到入口
    scalar downstream = 12.5;  // 半弦点到出口
    scalar halfHeight = 6.0;
    int nAirfoil = 300;        // 绕翼型单元数（上下表面合计，偶数）
    int nWake = 120;           // 每侧尾迹方向单元数
    int nNormal = 120;         // 法向单元数
    scalar firstCell = 1e-3;   // 壁面首层高度（弦长倍数）
    scalar wakeGrowth = 0.002; // 尾迹割线处首层高度随离后缘距离线性增大：firstCell + wakeGrowth·Δx（远尾迹不必贴壁级加密）
    int smoothIter = 0;        // 内部拉普拉斯光顺次数（靠壁若干层不动），默认不光顺
    scalar depth = 0.1;        // 展向厚度（二维）
};

RawMesh generateAirfoilCMesh(const AirfoilMeshSpec& s);

// NACA 四位数翼型表面点：theta∈[0,π] 余弦参数，upper=true 为上表面；返回未旋转坐标
Vec3 nacaSurfacePoint(const AirfoilMeshSpec& s, scalar xOverC, bool upper);

} // namespace cfd
