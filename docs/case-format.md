# 算例文件格式（JSON）

求解程序 `fvmFlow`（不可压流动）和 `fvmScalar`（标量输运 / 导热）都只读一个 JSON 算例文件，不需要改代码或重新编译。
文件允许 `//`、`/* */` 注释和尾随逗号。相对路径都相对于算例文件所在目录。
运行结束时会列出没有被读取的键（多半是拼写错误）。

```bash
fvmFlow   cases/pitzDaily/case.json
mpirun -np 8 fvmFlow cases/pitzDaily/case.json
```

## mesh

```jsonc
"mesh": { "polyMesh": "../../meshes/pitzDaily/polyMesh" }      // OpenFOAM polyMesh（ASCII 或二进制）
"mesh": { "file": "wing.msh",                                   // 外部网格文件，格式按内容自动识别：
                                                                //   Gmsh .msh（ASCII 2.2 / 4.1）、Fluent .msh/.cas（ASCII 或二进制段）、polyMesh 目录
    "format": "auto",                   // 可强制 "gmsh" / "fluent"
    "scale": 0.001,                     // 坐标缩放（如 mm → m），默认 1
    "depth": 0.1,                       // 二维网格（三角形/四边形/多边形）拉伸一层的厚度，前后面为 empty
    "patchTypes": { "airfoil": "wall", "frontAndBack": "empty" },  // 覆盖边界类型
    "writePolyMesh": "constant/polyMesh"  // 可选：顺便存成 polyMesh
} }
"mesh": { "box": {                                              // 内置长方体网格
    "n": [64, 64, 1], "lo": [0, 0, 0], "hi": [1, 1, 0.1],
    "periodic": [true, false, false],   // 周期方向
    "twoD": true,                       // 二维：z 向单层，前后面为 empty
    "walls": ["yMin", "yMax"],          // 设为 wall 类型的边界
    "names": { "yMax": "lid" },         // 重命名边界（xMin xMax yMin yMax zMin zMax）
    "stretch": { "y": 2.0 }             // 双侧 tanh 加密，参数越大越贴壁
} }
"mesh": { "airfoil": {                                          // 内置 NACA 四位数翼型二维 C 型网格
    "naca": "0012", "alpha": 4,         // 攻角：翼型绕半弦点旋转，来流保持 +x
    "upstream": 6.5, "downstream": 12.5, "halfHeight": 6,       // 半弦点到入口/出口/上下边界的距离（弦长倍数）
    "nAirfoil": 300, "nWake": 120, "nNormal": 120, "firstCell": 1e-3, "depth": 0.1
} }                                     // 边界：airfoil（wall）、inlet、outlet、top、bottom
```

外部网格的边界名：Gmsh 取 Physical 组名（三维用 Physical Surface，二维用 Physical Curve；没有分组的边界面归入
`defaultFaces`），Fluent 取 zone 名。边界类型默认：Gmsh 名称含 wall 的为 wall、含 symmetry 的为 symmetry，
Fluent 按 zone 类型（wall / symmetry / 其余为 patch），都可以用 `patchTypes` 改。支持四面体、六面体、三棱柱、
金字塔、任意多面体（Fluent / polyMesh）及其混合；Gmsh 高阶单元只取角点。周期边界和非协调（悬挂节点）网格暂不支持导入。

`meshCheck <网格> [--scale s] [--depth d] [--patch-type 名=类型] [--polymesh 目录] [--vtk 目录]` 可单独检查网格
（单元闭合性、面朝向、体积、非正交角、各边界面积）并转换成 polyMesh / VTK。

## physics

| 键 | fvmFlow | fvmScalar |
| --- | --- | --- |
| `nu` | 运动粘度 | |
| `steady` | true：SIMPLE 稳态；false：PIMPLE 瞬态 | 稳态 / 瞬态 |
| `field` | | 标量名，默认 `T` |
| `diffusivity` | | 扩散系数 Γ |
| `source` | | 均匀体积源 Q |
| `velocity` | | 给定均匀速度（缺省为纯扩散） |

## boundary

每个 patch 每个场一项。`wall` 类型的 patch 不写时，速度默认无滑移、压力默认零梯度；其余 patch 必须写。

| type | 参数 | 含义 |
| --- | --- | --- |
| `fixedValue` | `value`，可选 `ramp` | 定值（速度写 `[u, v, w]`）；`"ramp": {"duration": 2}` 时在 t∈[0, duration] 内由 0 平滑增到 value（smoothstep 6y⁵−15y⁴+10y³，`"shape": "linear"` 为线性） |
| `noSlip` | | 速度为零 |
| `zeroGradient` | | 零法向梯度 |
| `fixedGradient` | `gradient` | 给定法向梯度 |
| `robin` | `a` `b` `c` | a φ + b ∂φ/∂n = c |
| `symmetry` / `slip` | | 对称面 / 滑移壁 |

周期（cyclic）和 empty 边界由网格决定，不需要写。

## initial

`"initial": { "U": [0, 0, 0], "p": 0 }`，fvmScalar 为 `{ "T": 300 }`。

## schemes

| 键 | 可选值 | 默认 |
| --- | --- | --- |
| `div` | `upwind` `linear` `linearUpwind` `LUST` `limitedLinear 1` `vanLeer` `MUSCL` `minmod` | 稳态 `linearUpwind`，瞬态 `linear` |
| `ddt` | `Euler` `backward`（二阶，支持变步长） | `backward` |
| `grad` | `Gauss` `leastSquares` | `Gauss` |
| `laplacianNonOrthCorr` | 是否做显式非正交修正 | true |

## solution

| 键 | 含义 | 默认（稳态 / 瞬态） |
| --- | --- | --- |
| `nOuter` | PIMPLE 外迭代次数（1 即 PISO） | 1 |
| `nCorr` | 压力修正次数 | 1 / 2 |
| `nNonOrthCorr` | 非正交修正次数 | 0 |
| `alphaU` `alphaP` | 欠松弛 | 0.7, 0.3 / 1, 1 |
| `ddtPhiCoeff` | Rhie-Chow 时间修正：0 关闭，<0 OpenFOAM 自动系数，0~1 固定 | 0 |
| `U` `p` `pFinal` | 线性求解器（见下） | U：PBiCGStab+DILU；p：PCG+GAMG |
| `pRefCell` `pRefValue` | 无定值压力边界时的参考点 | 0, 0 |
| `verbose` | 打印每个线性求解 | false |

线性求解器：`{"solver": "PCG", "preconditioner": "GAMG", "tolerance": 1e-8, "relTol": 0.01, "maxIter": 2000}`。
solver 可选 `PCG`（对称）、`PBiCGStab`、`GaussSeidel`、`symGaussSeidel`、`Jacobi`；
preconditioner 可选 `GAMG`、`DIC`/`DILU`、`diagonal`、`none`。

大时间步（扩散数 νΔt/Δx² 远大于 1）时 PISO 可能不稳定，请增大 `nOuter`（PIMPLE）或 `nCorr`，或限制 `maxDt`。

## forcing（仅 fvmFlow）

```jsonc
"forcing": { "type": "constant", "force": [1.2, 0, 0] }       // 恒定体力（单位质量）
"forcing": { "type": "meanVelocity", "Ubar": [1, 0, 0] }      // 自动调节体力保持体平均速度（周期槽道）
```

## run

| 键 | 稳态 | 瞬态 |
| --- | --- | --- |
| `maxIter` `residualTol` | 最大迭代数、收敛残差 | |
| `endTime` | | 结束时间 |
| `dt` | | 固定步长 |
| `maxCo` `maxDt` | | 按库朗数自适应步长（优先于 dt）及上限 |
| `writeInterval` | 每多少步写 VTK | 每隔多少时间写 VTK（步长自动对齐写出时刻） |
| `printInterval` | 打印间隔（步） | 同左 |
| `output` | 输出目录，默认 `output` | |
| `vtk` / `tecplot` | 是否写 VTK（默认是）/ Tecplot `.dat`（默认否） | |
| `restartWrite` `restartRead` | 续算文件目录（与进程数无关，可换进程数续算） | |
| `printFluxes` | 打印时同时列出各边界体积通量（结束时总会列出） | |

## probes 与 lines

```jsonc
"probes": { "points": [[0.1, 0, 0], [0.2, 0, 0]], "interval": 1 },   // 时间序列，output/probes/*.csv
"lines":  [ { "name": "xMid", "start": [0.5, 0, 0], "end": [0.5, 1, 0], "n": 101 } ]  // 结束时沿线采样
```

探针取所在单元（最近体心）的值，结果与进程数无关。

## forces（仅 fvmFlow）

```jsonc
"forces": {
  "patches": ["airfoil"],               // 积分的壁面
  "Uref": 1, "lRef": 1,                 // 参考速度、参考长度；参考面积 = lRef × span（二维 span 自动取网格厚度）
  "dragDir": [1, 0, 0], "liftDir": [0, 1, 0], "pitchAxis": [0, 0, 1],
  "CofR": [0.25, 0, 0],                 // 力矩参考点
  "rho": 1, "pRef": 0,                  // p 为运动学压力 p/ρ
  "averageStart": 12.5,                 // 从该时刻起做时间加权平均（瞬态）
  "interval": 10                        // 每多少步写一行
}
```

输出 `output/forces.csv`（t, Cd, Cl, Cm 及压力/粘性分量、力 Fx Fy Fz），结束时打印平均 Cd、Cl、Cm，
并写 `output/surface.csv`（壁面面心坐标、法向、Cp、壁面切应力、Cf）。Cp = (p − pRef)/(½U²)，Cf 为切应力沿 dragDir 的分量 /(½U²)。

## turbulence（仅 fvmFlow）

```jsonc
"turbulence": {
  "model": "kOmegaSST",        // laminar（缺省）| Smagorinsky | WALE | kEpsilon | realizableKE | kOmega | kOmegaSST | SpalartAllmaras
  "wallFunctions": true,       // RANS：壁面函数（k-ε 必须；k-ω/SST 的 ω 壁面值对粘性/对数层自动混合，y+≈1 也可用）
  "coeffs": { "betaStar": 0.09 },   // 覆盖模型系数（名称同 OpenFOAM）
  "div": "limitedLinear 1",    // 湍流量对流格式
  "relax": 0.7,                // 稳态欠松弛
  "solver": { "solver": "PBiCGStab", "preconditioner": "DILU", "tolerance": 1e-10, "relTol": 0.1 },
  // LES：
  "vanDriest": true,           // Van Driest 近壁阻尼（Δ = min(Δ, κ y/CΔ (1 − e^{−y+/A+}))）
  "wallModel": {               // 壁面模型（WMLES），由采样速度按壁面律求 u_τ，给出壁面切应力
    "type": "spalding",        // spalding | musker | logLaw | powerLaw（Werner-Wengle）
    "samplingHeight": 0.1      // 采样点离壁距离；0（缺省）为壁面第一层单元
  }
}
```

模型场（k、epsilon、omega、nuTilda）的初值写在 `initial` 里，边界条件写在 `boundary` 里（与 U、p 相同的写法）；
wall 类型边界不写时自动设置（k、ε、ω 零梯度，壁面单元的 ε、ω 由壁面函数给定；nuTilda 为 0）；
其他边界不写时为零梯度并会给出提示，入口一般需要给定值。VTK 输出中会附带这些场和 nut，结束时打印壁面 y+。

LES 系数：Smagorinsky `Cs`（默认 0.17；OpenFOAM 默认相当于 0.168）、WALE `Cw`（0.325）、滤波尺度 `deltaCoeff`（1，Δ = 体积立方根，二维取面积平方根）。

验证（`cases/channel395`，Re_τ = 392 的 DNS）：k-ω SST 394、k-ω 399、Spalart-Allmaras 390（y+≈0.1）；
Re_τ ≈ 2000、y+≈50 用壁面函数时 k-ε 1957、SST 1949、realizable k-ε 1909（DNS 2003）。
