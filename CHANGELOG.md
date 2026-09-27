# Changelog

## v1.0.0-alpha - 新框架（第一阶段：基础平台）

整体重写为通用有限体积平台，旧 API（`src/core`、`example/`）移除，全部旧算例迁移到 `cases/`。

### 新增

- 不可压流动 `IncompressibleFlow`：稳态 SIMPLE 与瞬态 PIMPLE/PISO；Euler / backward（变步长）时间格式；
  恒定体力与平均流速保持驱动；湍流模型接口（ν_t、壁面有效粘度）。
- 标量输运 `ScalarTransport`：对流、扩散、源项，稳态 / 瞬态。
- 周期（cyclic）边界、对称 / 滑移边界；Robin 边界保留。
- 对流格式 linear、linearUpwind、LUST、vanLeer、MUSCL、minmod、limitedLinear；最小二乘梯度；非正交修正。
- 线性求解器 PCG、PBiCGStab、Gauss-Seidel；GAMG、DIC/DILU 预条件（压力方程默认 GAMG）。
- JSON 算例文件驱动的 `fvmFlow`、`fvmScalar`（见 docs/case-format.md）。
- VTK（ParaView）输出、Tecplot 输出、探针、沿线采样、与进程数无关的续算文件。
- 可复现模式 `CFD_REPRODUCIBLE=1`：任意进程数结果逐位相同。
- ctest 精度测试（Poisson、对流扩散、Robin、网格）。


## v0.2.0 - MPI parallel solver

### Added

- MPI domain-decomposition parallelism: every example runs with `mpirun -np N`.
- Built-in recursive coordinate bisection (RCB) decomposition; no external partitioner needed.
- `par::` communication layer: MPI environment, halo exchange (non-blocking,
  overlapped with computation), ordered global reductions, gather to rank 0.
- Results and output files are bitwise identical to the serial program for any
  number of processes.
- `tools/compare_parallel.sh` to check outputs across process counts.
- CMake option `CFD_USE_MPI` (ON by default; falls back to a single-process build
  when MPI is not found).

### Changed

- Example `main()` functions construct `par::Environment` first.
- `Mesh::getCellNumber()` returns the number of cells owned by the process;
  cell field storage uses `Mesh::getLocalCellNumber()` (owned + ghost cells).
- `SIMPLE::Options::pressureReferenceCell` is a global cell index.
- `example/div/scalar_fud.cpp` uses the relative mesh path like the other examples.
- `example/matrix/sparse_matrix_access_test.cpp` maps global cell indices to local ones.
- Compile with `-ffp-contract=off` so floating-point results do not depend on FMA contraction.

## v0.1.0 - First runnable FVM solver release

### Added

- OpenFOAM-style `polyMesh` reading.
- Mesh topology and geometry calculation.
- Scalar, vector, and tensor field abstractions.
- Finite volume discretization operators:
  - `fvm::Laplacian`
  - `fvm::Div`
  - `fvm::Source`
- CSR sparse matrix assembly.
- Jacobi linear solver with scalar and vector unknown support.
- SIMPLE algorithm for steady incompressible flow.
- Rhie-Chow interpolation for collocated grids.
- Non-orthogonal correction for diffusion terms.
- Example cases:
  - 2D lid-driven cavity
  - 2D backward-facing step
  - 3D lid-driven cavity
  - Laplacian tests
  - convection tests
  - source term tests
  - sparse matrix tests

### Known limitations

- Mainly supports steady problems.
- Linear solver is still primarily Jacobi-based.
- Turbulence models are not implemented.
- Boundary conditions are still being extended.
- Parallelism is currently shared-memory based only.
