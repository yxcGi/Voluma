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
"mesh": { "periodicHill": { "n": [128, 64, 64], "stretch": 1.5 } }   // 周期山贴体网格（ERCOFTAC 山形，Lx 9H、Ly 4.5H、高 3.036H，下壁 bottom、上壁 top）
"mesh": { "pipe": { "radius": 0.5, "length": 1, "nCore": 8, "nRing": 8, "nAxial": 16,   // 直圆管蝶形 O 型网格：中央方块 nCore²、
    "coreFrac": 0.5, "stretch": 0, "periodic": true } }           // 环向块径向 nRing、轴向 nAxial；壁面 wall，两端 cyclic zMin/zMax（periodic=false 时为 inlet/outlet）
"mesh": { "cylinder": { "nTheta": 192, "nRadial": 96, "nSpan": 32, "span": 3.14159, "firstCell": 0.0025,
                        "upstream": 10, "downstream": 20, "halfHeight": 10, "twoD": false } }   // 圆柱 O 型网格（直径 1，圆心在原点；边界 cylinder、inlet、outlet、top、bottom）
"mesh": { "airfoil": {                                          // 内置 NACA 四位数翼型二维 C 型网格
    "naca": "0012", "alpha": 4,         // 攻角：翼型绕半弦点旋转，来流保持 +x
    "upstream": 6.5, "downstream": 12.5, "halfHeight": 6,       // 半弦点到入口/出口/上下边界的距离（弦长倍数）
    "nAirfoil": 300, "nWake": 120, "nNormal": 120, "firstCell": 1e-3, "depth": 0.1,
    "wakeGrowth": 0.002                 // 尾迹割线处首层高度 = firstCell + wakeGrowth × 离后缘距离
} }                                     // 边界：airfoil（wall）、inlet、outlet、top、bottom
```

外部网格的边界名：Gmsh 取 Physical 组名（三维用 Physical Surface，二维用 Physical Curve；没有分组的边界面归入
`defaultFaces`），Fluent 取 zone 名。边界类型默认：Gmsh 名称含 wall 的为 wall、含 symmetry 的为 symmetry，
Fluent 按 zone 类型（wall / symmetry / 其余为 patch），都可以用 `patchTypes` 改。支持四面体、六面体、三棱柱、
金字塔、任意多面体（Fluent / polyMesh）及其混合；Gmsh 高阶单元只取角点。周期边界和非协调（悬挂节点）网格暂不支持导入。

`meshCheck <网格> [--scale s] [--depth d] [--patch-type 名=类型] [--polymesh 目录] [--vtk 目录]` 可单独检查网格
（输出单元数、体积、最大非正交角、平均 / 最大偏斜度、各 patch）。偏斜网格（三角形、四面体）上面插值自动加偏斜修正、
梯度默认最小二乘，见 `docs/validation.md` 六
（单元闭合性、面朝向、体积、非正交角、各边界面积）并转换成 polyMesh / VTK。

## physics

| 键 | fvmFlow | fvmScalar |
| --- | --- | --- |
| `nu` | 运动粘度 | |
| （传热见 `energy` 节） | | |
| `steady` | true：SIMPLE 稳态；false：PIMPLE 瞬态 | 稳态 / 瞬态 |
| `field` | | 标量名，默认 `T` |
| `diffusivity` | | 扩散系数 Γ |
| `source` | | 均匀体积源 Q |
| `velocity` | | 给定均匀速度（缺省为纯扩散） |

## boundary

每个 patch 每个场一项。`wall` 类型的 patch 不写时，速度默认无滑移、压力默认零梯度；其余 patch 必须写。
周期（cyclic）和 empty 边界由网格决定，不需要写。

### 定值类

| type | 参数 | 含义 |
| --- | --- | --- |
| `fixedValue` | 值的给法见下 | 定值（速度写 `[u, v, w]`） |
| `noSlip` | | 速度为零 |
| `movingWall` | `velocity`（给法同 `value`） | 平动壁面（= 速度 fixedValue，仅为可读性） |
| `rotatingWall` | `origin` `axis` `omega`（常数或时间表） | 旋转壁面 U = ω (axis × (x − origin))，自动去掉法向分量 |
| `flowRateInlet` | `volumetricFlowRate` 或 `massFlowRate` + `rho`（常数或时间表） | 均匀法向入流，速度 = Q / patch 总面积，方向为内法向（Fluent 质量流量入口） |

定值类的值恰好用下面一种方式给出，`ramp`、`noise` 可叠加：

```jsonc
"U": { "type": "fixedValue", "value": [1, 0, 0] }                                   // 常数
"U": { "type": "fixedValue", "timeTable": [[0, [0,0,0]], [1, [1,0,0]], [5, [1,0,0]]] }  // 时间表（分段线性），也可 {"file": "u_t.csv"}
"U": { "type": "fixedValue", "profile": { "file": "inlet.csv", "axis": "y" } }      // 剖面表：文件列 s, u, v, w（标量场为 s, φ），
                                                                                    //   s 为沿 axis（x/y/z）或 direction（相对 origin）的坐标；可 "scale"
"U": { "type": "fixedValue", "parabolic": { "Ubulk": 1, "radius": 0.5 } }           // 圆管抛物线 2 Ub (1 − r²/R²)，中心默认为 patch 面心
"U": { "type": "fixedValue", "parabolic": { "Ubulk": 1, "halfWidth": 0.5, "axis": [0, 1, 0] } }  // 平面槽道 1.5 Ub (1 − s²/h²)；可用 Umax 代替 Ubulk，
                                                                                    //   "center"、"direction"（默认内法向）可选
"U": { "type": "fixedValue", "value": [1, 0, 0], "ramp": { "duration": 2 } }        // 缓启动：t∈[0, duration] 内由 0 平滑增到值
                                                                                    //   （smoothstep 6y⁵−15y⁴+10y³，"shape": "linear" 为线性）
"U": { "type": "fixedValue", "value": [1, 0, 0], "noise": 0.05 }                    // LES 入口随机扰动：每面每步 ±noise·|U|（坐标 + 时间散列，与进程数无关）
```

### 梯度类与其他

| type | 参数 | 含义 |
| --- | --- | --- |
| `zeroGradient` | | 零法向梯度 |
| `fixedGradient` | `gradient` | 给定法向梯度 |
| `fixedFluxPressure` | | 压力梯度由求解器设定，使边界通量等于边界速度的通量（无体力时同 zeroGradient；浮力等体力进入通量后须用它） |
| `robin` | `a` `b` `c` | a φ + b ∂φ/∂n = c |
| `symmetry` / `slip` | | 对称面 / 滑移壁 |
| `mixed` | （程序内部） | w·定值 + (1−w)·定梯度，下列条件的基类 |

### 进出口类（需要流动信息，由求解器提供通量 / 速度 / 时间步）

| type | 场 | 参数 | 含义 |
| --- | --- | --- | --- |
| `inletOutlet` | 任意 | `inletValue`（给法同 `value`，缺省 0） | 流出面零梯度、回流面取 inletValue：**出口防回流**的标准写法 |
| `totalPressure` | p | `p0`（常数或时间表） | 总压入口：流入面 p = p0 − ½\|U_b\|²，流出面 p = p0（Fluent 压力入口） |
| `pressureInletOutletVelocity` | U | 可选 `tangentialVelocity` | 配 totalPressure / 定压边界：流出面零梯度，流入面法向分量零梯度、切向取给定值（默认 0） |
| `freestream` | U | `freestreamValue` | 自由来流：按来流方向定进出 w = ½ − ½ (U∞·n)/\|U∞\|（流入面定值、流出面零梯度），不随当地通量跳变 |
| `freestream` | 标量 | `freestreamValue` | 同 inletOutlet |
| `freestreamPressure` | p | `freestreamValue`，可选 `U` | 流出面定值 p∞、流入面零梯度；U∞ 取同一 patch 上 U 的 freestream 条件（或显式给 `U`） |
| `advective` | 任意 | 可选 `fieldInf` + `lInf` | 对流出流 ∂φ/∂t + U_n ∂φ/∂n = 0（与场的时间格式一致地隐式离散；可向远场值 fieldInf 以长度 lInf 松弛）。稳态或流入时退化为零梯度 |

湍流量（k、epsilon、omega、nuTilda）的入口可直接按 Fluent 的方式给强度与尺度，写在 patch 下的 `turbulence` 里，
各模型自动换成对应场的条件（流出面零梯度）：

```jsonc
"inlet": { "U": { "type": "fixedValue", "value": [1, 0, 0] }, "p": { "type": "zeroGradient" },
           "turbulence": { "intensity": 0.05, "lengthScale": 0.07 } }       // 或 "viscosityRatio": 10
```

k = 1.5 (I\|U_b\|)²；ε = Cμ^{3/4} k^{3/2}/L 或 Cμ k²/(r ν)；ω = √k/(Cμ^{1/4} L) 或 k/(r ν)；ν̃ 由 χ f_v1(χ) = r 反解。
也可逐场写 `{"type": "turbulentIntensity", "intensity": 0.05}`（k）、`{"type": "turbulentMixingLength", "lengthScale": 0.07}`、
`{"type": "viscosityRatio", "viscosityRatio": 10}`（ε、ω、nuTilda）。

**常用组合**（对应 Fluent 边界类型）：

| Fluent | Voluma |
| --- | --- |
| velocity inlet | U `fixedValue`（常数 / 剖面 / 时间表 / 抛物线），p `zeroGradient`，`turbulence` 强度 + 尺度 |
| mass flow inlet | U `flowRateInlet`，p `zeroGradient` |
| pressure inlet | p `totalPressure`，U `pressureInletOutletVelocity`（稳态时对流格式用 linearUpwind / limitedLinear，中心格式会发散） |
| pressure outlet | p `fixedValue`，U `inletOutlet`（标量也用 inletOutlet 防回流） |
| outflow（瞬态） | U / 标量 `advective`，p `fixedValue` 或 `zeroGradient` |
| pressure far-field（不可压） | U `freestream`，p `freestreamPressure` |
| wall / moving wall | 不写（无滑移）/ `movingWall` / `rotatingWall` |
| symmetry | `symmetry`（或网格 patch 类型 symmetry） |

验证算例（`cases/channelInlet`、`cases/taylorCouette`、`cases/advectiveOutlet`、`cases/jetBox`，都在 ctest 中）见 `docs/validation.md`。

## initial

`"initial": { "U": [0, 0, 0], "p": 0 }`，fvmScalar 为 `{ "T": 300 }`，可叠加高斯脉冲
`"gaussian": { "center": [0.5, 0, 0], "sigma": 0.1, "amplitude": 1 }`（A exp(−|x − c|²/σ²)）。

`"noise": 0.02` 在速度初场上叠加该幅值的随机扰动（与进程数无关），用于打破对称、触发脱落或转捩。

`"potentialFlow": true` 在求解前把初场投影成无散场（解一次泊松方程修正通量，单元速度由修正后的面通量重构，p 不变，
同 OpenFOAM potentialFoam）；曲面或周期域里给均匀初速度时建议打开，否则稳态 SIMPLE 的第一步可能发散。

`"hybrid": true` 为混合初始化（与 Fluent hybrid initialization 同效）：先做上述势流投影（入口 / 出口通量按边界条件给出，
壁面无穿透），再按 Bernoulli 给压力 p = H − ½|U|²（H 取定值压力边界上 p_b + ½|U_b|² 的面积平均，没有定值压力边界时为 0），
湍流量中没有在 `initial` 里给初值的场取入流边界值（湍流入口条件给出的 k、ε、ω 等）的面积加权平均。
直通道里给出精确的均匀速度与零压力（`tests/test_post`）；对圆柱 Re 40 这类均匀初场已经很好的算例收敛步数不变（1200 对 1134）。

槽道 LES 可用 `"channel"` 生成发展湍流用的初场（1/7 次方律平均剖面，体平均为 Ubulk，叠加无散的流向涡与条带扰动
和小幅随机扰动；随机数由单元坐标散列得到，与进程数无关）：

```jsonc
"initial": { "channel": { "Ubulk": 1, "streamwise": 0, "normal": 2, "walls": [0, 2], "amplitude": 0.1, "noise": 0.05 } }
```

三维 Taylor-Green 涡初场（周期域 [0, 2πL]³，u = U0 sin(x/L) cos(y/L) cos(z/L)，v = −U0 cos(x/L) sin(y/L) cos(z/L)，
w = 0，p 为对应的解析压力）：`"initial": { "taylorGreen": { "U0": 1, "L": 1 } }`。

## schemes

| 键 | 可选值 | 默认 |
| --- | --- | --- |
| `div` | `upwind` `linear` `linearUpwind` `LUST` `limitedLinear 1` `vanLeer` `MUSCL` `minmod` | 稳态 `linearUpwind`，瞬态 `linear` |
| `ddt` | `Euler` `backward`（二阶，支持变步长） | `backward` |
| `grad` | `Gauss` `leastSquares` | `Gauss`；平均偏斜度 > 0.05 的网格（三角形 / 四面体）默认 `leastSquares` |
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

打印与 `residualTol` 使用的 `U res` / `p res` 是每步第一次线性求解的初始残差 Σ|b − Ax|，归一化方式：
U 的三个分量共用分量中最大的 OpenFOAM 式归一化因子（相当于按速度模归一，恒为零的分量不会给出噪声残差）；
p 的归一化因子取当前值与前 5 步最大值中的较大者（压力恒为常数的周期体力驱动流动也能正常收敛判定）。
`verbose` 打印的各线性求解残差仍是 OpenFOAM 式逐分量归一化值。

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
| `restartWrite` `restartRead` | 续算文件目录（与进程数无关，可换进程数续算；时间平均、槽道 / 展向统计的累积量一并写入 / 读回） | |
| `printFluxes` | 打印时同时列出各边界体积通量（结束时总会列出） | |
| `meshQuality` | 把单元质量量（`nonOrthogonality` 所含内部面的最大非正交角、`skewness` 最大面偏斜度、`aspectRatio` 面积比 max\|S_f\|/min\|S_f\|）写进 VTK | |
| `nonOrthWarn` `skewWarn` `aspectWarn` | 坏单元阈值（默认 70°、1、1000）：启动时统计超过阈值的单元数并写 `output/badCells.csv`（全局单元号、中心、三个量） | |

残差 / 库朗数 / 连续性误差出现非有限值时求解器不会静默退出：打印非有限 U、p 的单元数与前几个单元的全局号和中心坐标，
写出当前场的 VTK 帧（便于在 ParaView 里定位发散起点），再以错误退出。湍流量的有界化（k、ε、ω 低于下限或非有限时取相邻有效单元
的平均，与 OpenFOAM `bound()` 同法）会在打印 y+ 的同时报告被修正的单元数（`turbulence: bounding k in N cells`）。

## monitor（仅 fvmFlow，瞬态）

```jsonc
"monitor": { "interval": 1, "file": "monitor.csv" }   // 每隔多少步写一行 output/monitor.csv
```

列：`t,k,epsilon,epsilon_sgs,enstrophy`，均为体积平均：k = ½⟨|u|²⟩，epsilon = 2ν⟨S:S⟩（解析尺度粘性耗散），
epsilon_sgs = 2⟨ν_t S:S⟩，enstrophy = ½⟨|ω|²⟩。总耗散 −dk/dt 可由 k(t) 差分得到（`tools/compare_tgv.py`）。

## probes 与 lines

```jsonc
"probes": { "points": [[0.1, 0, 0], [0.2, 0, 0]], "interval": 1 },   // 时间序列，output/probes/*.csv
"lines":  [ { "name": "xMid", "start": [0.5, 0, 0], "end": [0.5, 1, 0], "n": 101 } ]  // 结束时沿线采样
```

探针取所在单元（最近体心）的值，结果与进程数无关。沿线采样默认也取单元值（表中坐标为单元体心）；
加 `"interpolate": true` 则用单元梯度线性重构到采样点 φ_c + ∇φ_c·(x − x_c)，表中坐标为采样点本身，非结构网格上建议开启。
沿线采样写出 U、p、湍流场、ν_t、T 以及时间平均场（见 fieldAverage）。

## planes（仅 fvmFlow）

```jsonc
"planes": [ { "name": "midZ", "point": [0, 0, 0.5], "normal": [0, 0, 1] } ]
```

结束时对每个平面写 `output/plane_<name>_<field>.csv`（`x,y,z,值`）：取被平面切过的单元（按单元顶点到平面的有符号距离判定，
平面恰落在网格面上时只取一侧的一层），在单元中心到平面的投影点上用单元梯度线性重构给出场值（时间平均场取单元值），
按全局单元编号排序，结果与进程数无关。ParaView 中可用 Table To Points + Delaunay2D 显示。

## fieldAverage（仅 fvmFlow）

```jsonc
"fieldAverage": { "start": 100, "fields": ["U", "p", "T", "k", "nut"], "prime2Mean": true, "interval": 1 }
```

从 `start`（瞬态为时间，稳态为迭代步）起累积时间平均（瞬态按 dt 加权，稳态按迭代计数）：`<name>Mean` 与二阶矩
`<name>Prime2Mean` = ⟨φφ⟩ − ⟨φ⟩⟨φ⟩（向量场为 3×3 对称张量，即解析雷诺应力；`prime2Mean: false` 关闭），与 OpenFOAM fieldAverage
的 Mean / Prime2Mean 同义。默认场为 U、p；可选湍流模型场、`nut`、`T`。平均场随 VTK 各帧一起写出，也进 Tecplot、沿线 / 平面采样；
续算时累积量写入 / 读回续算目录（可换进程数），分段运行的平均与一次跑完逐位一致。

## wallOutput（仅 fvmFlow）

```jsonc
"wallOutput": { "patches": ["cylinder"], "Uref": 1, "pRef": 0, "flowDir": [1, 0, 0], "csv": true, "vtk": true, "everyWrite": false }
```

对所选壁面 patch（缺省为全部壁面）输出每个面的 p、Cp = (p − pRef)/(½U_ref²)、壁面切应力 τ_w（运动学 τ/ρ，壁面函数网格上含模型给出的
ν_eff,w）、Cf = τ_w·t/(½U_ref²)（给 `flowDir` 时取沿它的投影，否则取 \|τ_w\|）、y+（湍流模型的壁面单元 y+，层流为 y√\|τ_w\|/ν），
有能量方程时还有热流 q 与壁温 T_w：结束时写 `output/wall_<patch>.csv`（按全局面序）与 VTK 表面 `output/surface_<patch>.pvd`
（PolyData，每个进程一个分块），`everyWrite: true` 时每个 VTK 写出时刻也写一帧表面。圆柱 Re 40 上由 Cf 变号得到的分离角与
`forces` 的 surface.csv 逐位一致（53.65°，文献 53.6–53.7°）。

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

## channelStatistics（仅 fvmFlow）

累积量随续算写入 / 读回（`restartWrite` / `restartRead`），分段运行可以继续累积；spanStatistics 同。

槽道湍流统计，输出格式与 OpenLB `channel3d` 算例相同，可直接与其结果及 Lee & Moser DNS 对比：

```jsonc
"channelStatistics": {
  "streamwise": 0,                      // 流向坐标轴（0/1/2）
  "normal": 2,                          // 壁面法向坐标轴；展向为剩下的轴
  "walls": [0, 2],                      // 两壁面在法向上的坐标
  "start": 2500,                        // 从该时刻起做时间 + 平面平均
  "file": "channel_stats.csv"
}
```

每个写出时刻和结束时写 `output/channel_stats.csv`（分号分隔）：
`y+;u_tau;uAv+;uu++;uv++;uw++;vv++;vw++;ww++;pRMS;uu_sgs++;…;ww_sgs++`。
u、v、w 依次为流向、展向、法向分量；上下两半按到最近壁面的距离合并；u_τ 由壁面平均切应力（含壁模型的 ν_eff）
的时间平均得到；亚格子应力为 −2ν_t S_ij。

## spanStatistics（仅 fvmFlow）

展向均匀流动（周期山等）的展向 + 时间平均统计，输出与 OpenLB `periodichill3d` 同格式：

```jsonc
"spanStatistics": {
  "streamwise": 0, "normal": 2,         // 展向为剩下的轴；同一列单元的流向、法向中心坐标须相同（拉伸网格）
  "wall": "bottom",                     // 统计壁面切应力、壁压的壁面
  "Uref": 1, "H": 1,                    // 归一化速度、长度
  "start": 90,                          // 开始统计的时刻
  "stations": [0.05, 0.5, 1, 2],        // 剖面站位 x/H
  "dir": "hill"
}
```

写 `output/hill/profile_xh_<站>.csv`（`y_h;u_Ub;v_Ub;uu_Ub2;vv_Ub2;ww_Ub2;uv_Ub2;k_Ub2`，v 为法向、ww 为展向）
和 `wall_<wall>.csv`（`x_h;Cf;Cp;ut_Ub`），并打印按 Cf 变号求得的分离点、再附点。
与实验对比：`tools/compare_hill.py output/hill --ref <UFR3-30 实验数据目录> --out hill.png`。

## energy（仅 fvmFlow）

能量方程（常物性不可压；温度 T 为运动学量，热流以 q/(ρ c_p) 计，单位 K·m/s）：

```jsonc
"energy": {
  "Pr": 0.71,                     // 分子普朗特数（α = ν/Pr），或直接给 "alpha"
  "Prt": 0.85,                    // 湍流普朗特数：α_eff = α + ν_t/Pr_t
  "source": 0,                    // 均匀体积热源 Q/(ρ c_p)（K/s）
  "buoyancy": { "g": [0, -9.81, 0], "beta": 3.4e-3, "TRef": 300 },  // Boussinesq：f = −β (T − T_ref) g
  "wallFunction": "jayatilleke",  // RANS 壁面函数网格上的热壁面律：jayatilleke（默认）| kader | none
  "nusselt": { "L": 1, "deltaT": 1 },  // 结束时各壁面 patch 的 Nu = q̄ L / (α ΔT)
  "div": "linearUpwind", "relax": 0.7, "solver": { "solver": "PBiCGStab", "preconditioner": "DILU" }
}
```

边界：`boundary.<patch>.T` 可用 `fixedValue`（等温）、`zeroGradient`（绝热，壁面缺省）、`fixedGradient`、`robin` 等通用类型，
以及 `{"type": "heatFlux", "q": 0.1}`（给定进入流体的热流 q/(ρ c_p)，按当前壁面 α_eff,w 换算为法向梯度）。初值 `initial.T`。

- 浮力：p 为去掉参考静水压的动压，壁面 p 自动取 `fixedFluxPressure`（体力进入 φ_HbyA，线性体力 / 二次压力的离散平衡精确，
  稳定分层的静止方腔不产生伪速度）。
- 壁面扩散系数：壁面 ν_eff,w = ν（低 Re、LES）时 α_eff,w = α；RANS 壁面函数（ν_eff,w > ν）时 α_eff,w = ν y+/T+(y+)，
  y+ 由壁面函数给出的切应力反算（u_τ² = ν_eff,w U_t/y，与动量方程一致）。Jayatilleke：T+ = Pr y+（y+ < y+_T）
  或 Pr_t[ln(E y+)/κ + P]，P = 9.24[(Pr/Pr_t)^{3/4} − 1][1 + 0.28 e^{−0.007 Pr/Pr_t}]（Fluent / OpenFOAM
  alphatJayatillekeWallFunction）；Kader 1981 全 y+ 混合律 T+ = Pr y+ e^{−Γ} + [2.12 ln(1 + y+) + β(Pr)] e^{−1/Γ}。
- 输出：T 写入 VTK / Tecplot / 采样线（`line_<name>_T.csv`）与续算；结束时打印各壁面 patch 的面积、总热流、平均热流、
  平均壁温与 Nu，并写逐面表 `output/wallHeatFlux_<patch>.csv`（x y z q Tw）。稳态收敛判据同时要求 T 方程残差 < residualTol。
- 求解顺序：每个时间步 / SIMPLE 迭代先解流动再解 T（T 对浮力滞后一步）。

验证见 `docs/validation.md`：de Vahl Davis 差分加热方腔 Ra 1e3–1e6（`cases/naturalConvection`）、Graetz 层流圆管热入口段
（`cases/graetz`）、湍流槽道传热（`cases/channelHeat`：SST 低 Re、k-ε 壁面函数 + Jayatilleke / Kader）。

## turbulence（仅 fvmFlow）

```jsonc
"turbulence": {
  "model": "kOmegaSST",        // laminar（缺省）| Smagorinsky | WALE | Vreman | sigma | dynamicSmagorinsky
                               // | kEpsilon | RNGkEpsilon | realizableKE | kOmega | kOmegaSST | BSLEARSM | SpalartAllmaras
                               // | LRR | SSG（雷诺应力输运模型，需壁面函数网格 y+ ≥ 15）
  "wallFunctions": true,       // RANS：壁面函数（k-ε 系列必须；k-ω/SST 默认开，false 为低 Re 处理；SA 默认关，true 时用 Spalding 律）
  "wallTreatment": "standard", // RANS：standard（nutkWallFunction，y+ 切换）| enhanced（全 y+ 处理，见下）
  "coeffs": { "betaStar": 0.09 },   // 覆盖模型系数（名称同 OpenFOAM）
  "nutMaxRatio": 1e5,          // ν_t/ν 上限（Fluent 的 turbulent viscosity ratio limit），超限单元数在 y+ 统计行后打印
  "div": "limitedLinear 1",    // 湍流量对流格式
  "relax": 0.7,                // 稳态欠松弛
  "solver": { "solver": "PBiCGStab", "preconditioner": "DILU", "tolerance": 1e-10, "relTol": 0.1 },
  // LES：
  "vanDriest": true,           // Van Driest 近壁阻尼（Δ = min(Δ, κ y/CΔ (1 − e^{−y+/A+}))；动态模型忽略）
  "wallModel": {               // 壁面模型（WMLES），由采样速度按壁面律求 u_τ，给出壁面切应力
    "type": "spalding",        // spalding | musker | reichardt | logLaw | powerLaw（Werner-Wengle）| ode
    "samplingHeight": 0.1,     // 采样点离壁距离；0（缺省）为壁面第一层单元
    "pressureGradient": false, // 仅 ode：薄层方程含切向压力梯度与体力（非平衡壁面模型）
    "kappa": 0.41, "B": 5.2,   // 对数律常数（spalding / logLaw / reichardt / ode）
    "A": 8.3, "Bpow": 0.142857,// 幂律 u+ = A y+^B
    "Aplus": 17, "nPoints": 40 // ode：混合长度阻尼常数、一维积分点数
  }
}
```

模型场（k、epsilon、omega、nuTilda）的初值写在 `initial` 里，边界条件写在 `boundary` 里（与 U、p 相同的写法）；
wall 类型边界不写时自动设置（k、ε、ω 零梯度，壁面单元的 ε、ω 由壁面函数给定；nuTilda 为 0）；
其他边界不写时为零梯度并会给出提示，入口一般需要给定值。VTK 输出中会附带这些场和 nut，结束时打印壁面 y+。

### RANS 壁面处理

| `wallTreatment` | 壁面切应力 | 壁面单元产生项 G | 壁面单元 ε / ω | 适用 y+ |
| --- | --- | --- | --- | --- |
| `standard` | nutkWallFunction：y+ = C_μ^¼√k y/ν，y+ > 11.5 时 ν_t,w = ν(κy+/ln(Ey+) − 1)，否则 0 | y+ > 11.5 时对数律 u_τ³/(κy)，否则单元值 | ε：y+ 切换 2νk/y² / C_μ^¾k^{3/2}/(κy)；ω：Menter 粘性/对数混合 | k-ε 要求 y+ > 15；k-ω/SST 任意（缓冲层偏差较大） |
| `enhanced` | Spalding 统一律由壁面单元切向速度反解 u_τ（任意 y+ 连续） | P = u_τ⁴/ν·D(1−D)，D = du+/dy+（Spalding 律导数；粘性底层 → 0，对数层 → u_τ³/(κy)） | ε：Wolfshtein 长度尺度 k^{3/2}/ℓ_ε（两极限同 standard），并对 Re_y = √k y/ν < 200 的单元两层处理（ε 代数值、ν_t 按 λ_ε 与一方程模型混合，Fluent Enhanced Wall Treatment 同源，常数 `Amu` 70、`Aeps` 2C_ℓ、`ReyStar` 200、`ReyWidth` 20）；ω 同 standard | 全部 y+（k-ε 系列也可用于 y+ ≈ 1 网格） |

`enhanced` 自动把 `wallFunctions` 置为 true。两种处理在 y+ ≈ 5–10（壁面单元落在缓冲层）都会偏高 5–9 %，
这是壁面单元离散的固有误差（低 Re 解析的 SA 在同样网格上偏差相同），网格设计仍应避开缓冲层；
数值见 `docs/validation.md` 的 y+ 扫描表（`tools/channel_yplus_sweep.py` 可重跑）。

### LES 系数与壁面模型

LES 系数：Smagorinsky `Cs`（默认 0.17；OpenFOAM 默认相当于 0.168）、WALE `Cw`（0.325）、Vreman `c`（0.07 ≈ 2.5 C_s²）、
σ 模型 `Csigma`（1.35）、动态 Smagorinsky `CsMax`（0.3，系数上限 C_s² ≤ CsMax²，下限 0）、滤波尺度 `deltaCoeff`（1，Δ = 体积立方根，二维取面积平方根）。
动态 Smagorinsky 用 Germano-Lilly 过程：检验滤波为面邻接平均（有效宽度 2Δ），C = ⟨L:M⟩/⟨M:M⟩，分子分母再滤波一次作局部平均；
层流 / 分辨良好区域 C → 0，不需要 Van Driest 阻尼。Vreman 与 σ 模型在层流剪切、二维流动中 ν_t = 0，近壁 ν_t 自然衰减。
Taylor-Green 涡 Re 1600、64³ 的总耗散峰值（DNS 0.0129 @ t 9.0）：无模型 0.0104 @ 8.4、WALE 0.0105 @ 8.5、Smagorinsky 0.0106 @ 7.4、
Vreman 0.0103 @ 7.7、σ 0.0102 @ 8.2、动态 Smagorinsky 0.0109 @ 8.4（图 `docs/validation/tgv3d_re1600.png`）：
动态模型在转捩早期系数自动趋零、峰值最接近 DNS；Vreman / σ 在这类无壁面转捩流中与 Smagorinsky 一样早期偏耗散，
它们的优势在近壁与层流区 ν_t 自动为零（槽道、分离流）。

壁面模型（WMLES）各律由采样点（第一层单元或 `samplingHeight` 处）的切向相对速度求 u_τ：

| `type` | 公式 | 说明 |
| --- | --- | --- |
| `spalding` | Spalding (1961) 统一律 | 全 y+，默认 |
| `musker` | Musker (1979) | 全 y+ |
| `reichardt` | Reichardt (1951) u+ = ln(1+κy+)/κ + (B − ln κ/κ)(1 − e^{−y+/11} − (y+/11)e^{−y+/3}) | 全 y+；原文常数 7.8 对应 κ 0.4、B 5.5 |
| `logLaw` | u+ = ln(E y+)/κ，y+ < 11 线性 | 采样点须在对数层 |
| `powerLaw` | u+ = A y+^B（8.3，1/7），y+ < 11.8 线性 | 按采样点点值反解（Werner-Wengle 原文的单元积分形式要求 U 为单元平均，这里不用） |
| `ode` | 一维薄层方程 d/dy[(ν+ν_t)du/dy] = −f_p 数值解，ν_t = κ y u*(1 − e^{−y u*/(νA+)})² | 等价于 Van Driest 剖面；`pressureGradient: true` 时 f_p = −∇_t p + 体力（壁面单元值），u* = (u_τ³ + u_p³)^{1/3}，u_p = (ν\|f_p\|)^{1/3}（Duprat 等 2011），壁面切应力为矢量，可与采样速度方向不同；隐式部分沿壁面单元速度方向，其余作为显式动量源 |

`ode` 不含压力梯度时与 Spalding / Reichardt 一样是平衡律（对 Lee & Moser Re_τ 1000 DNS 剖面在 y+ = 10–300 反解 u_τ 误差 < 0.8 %）；
含压力梯度的版本用于分离流（周期山等），在槽道中因总切应力随高度线性减小而使 u_τ 偏高 1–4 %（采样越高偏差越大），槽道请用平衡律。

`BSLEARSM`：Menter BSL k-ω + Wallin & Johansson (2000) 显式代数雷诺应力（Fluent/CFX 的 BSL EARSM 同源）。
各向异性张量用 WJ 的二维解析 N，线性部分给 ν_t = −½β1 τ k，非线性部分作为显式动量源 −∇·τ_ex；
时间尺度 τ = max(1/(β*ω), `Ctau`√(ν/(β*kω)))（Ctau 6），近壁阻尼 f1 = 1 − exp(−`Cy` √k y/ν)（Cy 0.026）。
`nutMax`（默认 1）把 ν_t 连同整个各向异性张量限制到 ν_t ≤ nutMax·k/ω：WJ 原模型在弱应变区 C_μ^eff 可达 0.33，
分离回流区混合偏强（周期山再附 x/H 3.2），限制后为 5.3（LES 4.8，SST 7.6）；设 0 关闭。`nonlinear` 0 可关掉非线性应力（诊断）。
用于分离流（周期山、后台阶等）时比 SST 明显改善回流区长度；槽道 Re_τ 392 给 385。

`LRR` / `SSG`：雷诺应力输运模型（6 个分量 `Rxx Rxy Rxz Ryy Ryz Rzz` + `epsilon`，均可输出与采样）。压力-应变项 LRR 为
Launder-Reece-Rodi 1975 的 IP 形式（无壁面反射项；`C1` 1.8、`C2` 0.6、`Ceps2` 1.92），SSG 为 Speziale-Sarkar-Gatski 1991 二次模型
（`C1` 3.4、`C1s` 1.8、`C2` 4.2、`C3` 0.8、`C3s` 1.3、`C4` 1.25、`C5` 0.4、`Ceps2` 1.83）；`Ceps1` 1.44，各向同性扩散 `sigmaR` 0.82、
`sigmaEps` 1.0，系数同 OpenFOAM。动量方程按"表观粘度"耦合：ν_t = C_μ k²/ε 隐式，显式源 −∇·dev(R) 再减去已含的 ν_t 部分，
收敛后总应力即 dev(R)；壁面单元相邻的面用与隐式项一致的紧致梯度，其余面用单元梯度插值以抑制 U–R 奇偶解耦。
壁面处理同 k-ε（standard / enhanced），另有 `wallAnisotropy`（默认 true）：对数层壁面单元的应力按近壁平衡各向异性给定
R_ij = a_ij k（uu/k 1.098、nn/k 0.247、ss/k 0.655、−un/k 0.255，Fluent 的 RSM 壁面处理），k 由标准 k 方程在壁面单元上求得
（场 `k` 在内部单元即 ½ tr R）；关掉时壁面单元 R 取零梯度（OpenFOAM kqRWallFunction）。`relaxR`（稳态默认 0.5）为 R 分量方程的欠松弛。
适用壁面函数网格（首层 y+ ≥ 15；y+ < 11.5 时壁面单元退化为层流关系、无低 Re 修正，误差 > 10 %）。
验证（槽道，均匀 32 层、y+ ≈ 60，DNS Re_τ 1995）：LRR 1985（−0.5 %）、SSG 1927（−3.4 %，对数区斜率略陡，是该模型 ε 方程常数
给出的有效 κ ≈ 0.38 所致）；Re_τ 392、y+ 16–24 的网格上 LRR −2 ~ −3 %、SSG −4.6 %。

`RNGkEpsilon`：Yakhot 等 1992，系数同 OpenFOAM（C_μ 0.0845、C1 1.42、C2 1.68、σ_k = σ_ε 0.71942、`eta0` 4.38、`beta` 0.012），
ε 方程的 (C1 − R) G ε/k 项中 R = η(1 − η/η0)/(1 + βη³)，η = √(2S:S) k/ε；壁面函数与两层模型的常数仍用 C_μ 0.09（`Cmu`）。

验证（`cases/channel395`，Re_τ = 392 的 DNS）：k-ω SST 394、k-ω 399、Spalart-Allmaras 390、k-ε enhanced 396（y+≈0.1）；
Re_τ ≈ 2000、y+≈50 用壁面函数时 k-ε 1957、SST 1949、realizable k-ε 1909（DNS 2003），y+≈60 时 LRR 1985、SSG 1927（DNS 1995）；
完整 y+ 扫描见 `docs/validation.md`。
