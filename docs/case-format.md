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
"mesh": { "polyMesh": "../../meshes/pitzDaily/polyMesh" }      // OpenFOAM polyMesh（ASCII）
"mesh": { "box": {                                              // 内置长方体网格
    "n": [64, 64, 1], "lo": [0, 0, 0], "hi": [1, 1, 0.1],
    "periodic": [true, false, false],   // 周期方向
    "twoD": true,                       // 二维：z 向单层，前后面为 empty
    "walls": ["yMin", "yMax"],          // 设为 wall 类型的边界
    "names": { "yMax": "lid" },         // 重命名边界（xMin xMax yMin yMax zMin zMax）
    "stretch": { "y": 2.0 }             // 双侧 tanh 加密，参数越大越贴壁
} }
```

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
| `fixedValue` | `value` | 定值（速度写 `[u, v, w]`） |
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

## probes 与 lines

```jsonc
"probes": { "points": [[0.1, 0, 0], [0.2, 0, 0]], "interval": 1 },   // 时间序列，output/probes/*.csv
"lines":  [ { "name": "xMid", "start": [0.5, 0, 0], "end": [0.5, 1, 0], "n": 101 } ]  // 结束时沿线采样
```

探针取所在单元（最近体心）的值，结果与进程数无关。
