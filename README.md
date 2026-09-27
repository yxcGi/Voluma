# Voluma

**Voluma** 是一个用 C++20 编写的通用有限体积法（FVM）CFD 平台，目标是逐步发展成类似 Fluent 的通用流动求解软件：
任意多面体非结构网格、MPI 并行、算例文件驱动，并提供湍流（RANS / LES / 壁面模型）等工程与科研所需的模型。

当前版本（第一阶段"基础平台"）已具备：

- **不可压流动**：稳态 SIMPLE、瞬态 PIMPLE/PISO（Rhie-Chow 插值），恒定体力与"平均流速保持"驱动（周期槽道）；
- **标量输运 / 导热**：对流-扩散-源项，定值、定梯度、Robin、对称等边界；
- **网格**：OpenFOAM polyMesh（任意多面体，含 cyclic、empty），内置带加密的长方体网格生成；
- **离散**：Euler / 二阶 backward（变步长）；upwind、linear、linearUpwind、LUST、TVD（vanLeer、MUSCL、minmod、limitedLinear）；Gauss / 最小二乘梯度；非正交修正；
- **线性求解器**：PCG、PBiCGStab、Gauss-Seidel、Jacobi；预条件 GAMG（代数多重网格）、DIC/DILU、对角；
- **并行**：MPI 区域分解（RCB），可选"可复现模式"使任意进程数结果逐位相同；
- **输入输出**：JSON 算例文件；VTK（ParaView）与 Tecplot；探针、沿线采样；与进程数无关的续算文件。

原先的求解器保留在 [CFD_FVM_Solver_MPI](https://github.com/yxcGi/CFD_FVM_Solver_MPI)（已冻结）。本仓库继承其提交历史，并把原有全部算例迁移到了新框架（见[算例](#算例)）。

---

## 构建

依赖：C++20 编译器、CMake ≥ 3.16；MPI 可选（找不到时自动构建串行版）。

```bash
cmake -S . -B build            # 加 -DFVM_USE_MPI=OFF 构建串行版
cmake --build build -j
ctest --test-dir build         # 单元与精度测试
```

可执行文件在 `build/bin/`：

| 程序 | 用途 |
| --- | --- |
| `fvmFlow` | 通用不可压流动求解（JSON 算例驱动） |
| `fvmScalar` | 标量输运 / 导热（JSON 算例驱动） |
| `tgv2d` | 二维 Taylor-Green 涡，与解析解比较（时间、空间精度验证） |
| `cavity` | 方腔驱动（命令行参数版，Ghia 对比） |

## 运行

```bash
./build/bin/fvmFlow cases/pitzDaily/case.json
mpirun -np 8 ./build/bin/fvmFlow cases/pitzDaily/case.json
```

结果写到算例目录下的 `output/`：`fields.pvd` 用 ParaView 打开；`line_*.csv` 为沿线采样；`probes/` 为探针时间序列；
`"tecplot": true` 时另写 Tecplot 文件。算例文件的全部选项见 [docs/case-format.md](docs/case-format.md)。

一个最小算例：

```jsonc
{
  "mesh": { "polyMesh": "../../meshes/pitzDaily/polyMesh" },
  "physics": { "nu": 0.01, "steady": true },
  "boundary": {
    "inlet":  { "U": {"type": "fixedValue", "value": [20, 0, 0]}, "p": {"type": "zeroGradient"} },
    "outlet": { "U": {"type": "zeroGradient"}, "p": {"type": "fixedValue", "value": 0} }
  },
  "schemes": { "div": "MUSCL" },
  "run": { "maxIter": 5000, "residualTol": 1e-6 }
}
```

## 算例

| 目录 | 程序 | 内容 |
| --- | --- | --- |
| `cases/pitzDaily` | fvmFlow | 后台阶层流，Re≈50（OpenFOAM pitzDaily 网格） |
| `cases/cavity` | fvmFlow | 方腔 Re=100，内置网格，与 Ghia 对比 |
| `cases/cavity2D_quad` | fvmFlow | 方腔 Re=5000，Gmsh 四边形网格 |
| `cases/cavity2D_tri` | fvmFlow | 方腔 Re=5000，三角形非结构网格 |
| `cases/cavity3D` | fvmFlow | 三维方腔 Re=1600，四面体网格 |
| `cases/poiseuille` | fvmFlow | 周期平面 Poiseuille 流（瞬态启动、平均流速驱动），与解析解对比 |
| `cases/heatSource` | fvmScalar | 带内热源的稳态导热 |

网格在 `meshes/`。

## 验证

| 项目 | 结果 |
| --- | --- |
| Poisson 制造解（`tests/test_operators`） | 盒子网格二阶收敛（阶数 2.00）；30° 非正交三角形网格上，非正交修正使误差降为约 1/8 |
| 一维对流扩散 Pe=10 | linear 二阶（2.00），upwind 一阶 |
| Robin 边界 | 线性解精确再现（误差 1e-13） |
| 二维 Taylor-Green 涡（`tgv2d`） | backward 时间二阶；速度空间误差约三阶（超收敛），压力二阶 |
| 方腔 Re=100，64² | 中心线 u 与 Ghia (1982) 最大偏差 0.0016 |
| 平面 Poiseuille（平均流速驱动） | 中心速度 1.4971（解析 1.4985），驱动压力梯度 1.1977（解析 1.2） |
| pitzDaily 与旧程序 | 同一网格、同一格式下速度场差异 rms 0.15%（最大 2%，在入口相邻单元），用时 405 s → 14 s |
| 并行 | 可复现模式下 np = 1, 2, 3, 4 结果逐位相同（TGV、pitzDaily） |

## 并行与可复现模式

默认追求速度：压力用 GAMG 预条件，全局求和用普通 `MPI_Allreduce`。与 Fluent、OpenFOAM 一样，换进程数后结果在求解容差量级内一致，但不逐位相同。

设置环境变量 `CFD_REPRODUCIBLE=1` 打开可复现模式：全局求和改为与顺序无关的精确求和，线性求解器自动把随分区变化的预条件 / 光顺（GAMG、DIC、DILU、Gauss-Seidel）换成与分区无关的对角 / Jacobi。
这时任意进程数的结果都与串行逐位相同，适合核对并行正确性，但大算例会慢几倍。

```bash
CFD_REPRODUCIBLE=1 mpirun -np 4 ./build/bin/fvmFlow cases/pitzDaily/case.json   # checksum 与 np=1 相同
```

## 代码结构

```
src/fvm/
  core/            基本类型（Vec3、Tensor）、命令行参数、JSON
  parallel/        MPI 封装、精确求和、halo 交换、全局编号
  mesh/            网格读取 / 生成（RawMesh）、分解与几何（Mesh）
  field/           体场 VolField、边界条件
  linalg/          LDU 矩阵 FvMatrix、Krylov 求解器、GAMG
  discretization/  格式，fvm（隐式）与 fvc（显式）算子
  solvers/         IncompressibleFlow（SIMPLE/PIMPLE）、ScalarTransport
  io/              VTK、Tecplot、探针、续算
  app/             算例文件 → 网格、边界、求解参数
apps/              fvmFlow、fvmScalar
examples/          tgv2d、cavity（代码驱动的验证算例）
tests/             单元与精度测试（ctest）
cases/  meshes/    算例与网格
```

矩阵语义与 OpenFOAM 的 fvMatrix 一致（LDU 存储、internalCoeffs / boundaryCoeffs、faceFluxCorrection），便于对照。
湍流模型通过 `IncompressibleFlow::updateTurbulence` 回调提供 ν_t 与壁面有效粘度，后续的 RANS、LES、壁面模型都从这里接入。

## 路线图

1. ✅ 基础平台：瞬态、周期、体力、GAMG、算例文件、VTK / Tecplot / 探针 / 续算
2. NACA0012 Re=1000 等二维层流外流，与 OpenLB 对比
3. LES：Smagorinsky（Van Driest 阻尼）、WALE；统计量输出与 OpenLB 的 CSV 格式一致
4. 壁面模型 WMLES：幂律、Spalding、Musker，可调采样高度；Re_τ=1000 槽道
5. RANS：k-ω SST、Spalart-Allmaras；周期山、圆柱绕流等
6. 通用化：传热耦合、可压缩、更多边界（对流出口、入口剖面）、GPU

已知限制：TVD 类格式（MUSCL 等）在个别算例上会让 SIMPLE 残差停在 1e-4 左右（限制器反复切换），此时可改用 linearUpwind；
瞬态大步长（扩散数远大于 1）时 PISO 可能不稳定，需要增大 `nOuter` 或限制步长。

## 许可证

MIT，见 [LICENSE](LICENSE)。
