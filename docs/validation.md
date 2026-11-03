# Voluma 验证记录

每一行是一个算例：参考数据来源、网格、指标、误差、状态。全部结果都在本仓库的算例目录下可以重跑，
`ctest -L case` 会自动重跑其中的快速算例并按容差判定（`tools/check_case.py`）。
"通过"的含义是真正运行并与参考数据比过，不是只编译通过。

## 一、代码验证（解析解 / 制造解）

| 算例 | 参考 | 网格 | 指标 | 结果 | 状态 |
| --- | --- | --- | --- | --- | --- |
| Poisson 制造解（`tests/test_operators`） | 解析 | 盒子、30° 非正交三角形 | 收敛阶 | 2.00；非正交修正使误差降为约 1/8 | ✅ ctest |
| 一维对流扩散 Pe = 10 | 解析 | — | 收敛阶 | linear 2.00，upwind 一阶 | ✅ ctest |
| Robin 边界 | 解析 | — | 误差 | 1e-13 | ✅ ctest |
| 二维 Taylor-Green 涡（`tgv2d`） | 解析 | 盒子 | 时间/空间阶 | backward 时间二阶；速度约三阶、压力二阶 | ✅ |
| 平面 Poiseuille（`cases/poiseuille`） | 解析 | 盒子、周期 | 中心速度、压力梯度 | 1.4971（1.4985）、1.1977（1.2） | ✅ |
| 圆管 Hagen-Poiseuille（`cases/pipe`） | 解析 u = 2U_b(1 − r²/R²)，dp/dx = 32νU_b/D² | 蝶形 O 型网格，最大非正交 17°–41° | 剖面最大误差（以中心速度归一）、dp/dx | nCore=nRing=4/8/16/32：直径线 3.2 / 0.82 / 0.21 / 0.05 %（二阶）；对角线 3.4 / 1.35 / 0.80 / 0.55 %（块角偏斜单元处饱和）；dp/dx 0.3278 / 0.3233 / 0.3216 / 0.3208（圆 0.32，内接多边形面积略小） | ✅ ctest（8³，容差 1.5 %） |

## 二、层流验证（文献）

| 算例 | 参考 | 网格 | 指标 | 结果 | 状态 |
| --- | --- | --- | --- | --- | --- |
| 方腔 Re = 100（`cases/cavity`） | Ghia 1982 | 64² | 中心线 u 最大偏差 | 0.0016 | ✅ |
| 后台阶 pitzDaily Re ≈ 50 | 旧代码（一致性，非验证） | OpenFOAM pitzDaily | 速度场 rms 差 | 0.15 % | ⚠️ 仅一致性 |
| NACA0012 Re = 1000（`cases/naca0012`） | Di Ilio 2018 | C 型网格 | Cd、Cl | α=0：Cd 0.1202（0.119）；α=4：0.1251 / 0.2052（0.1252 / 0.2022） | ✅ |
| 圆柱 Re = 40 定常（`cases/cylinderRe40`） | Dennis & Chang 1970 Cd 1.522、L_w 2.35；Fornberg 1980 1.498 / 2.24；Linnick & Fasel 2005 1.54 / 2.28 / θ 53.6°；Taira & Colonius 2007 1.54 / 2.30 / 53.7°；Coutanceau & Bouard 1977 实验 L_w 2.13 | O 型 192×110，域 40D/60D/±40D | Cd（压力 + 粘性）、回流区长 L_w/D、分离角 | Cd 1.517（0.994 + 0.523）、L_w 2.22、θ_sep 53.7° | ✅ ctest |
| 同上，域大小敏感性 | — | ±15D（阻塞 3.3 %） | Cd | 1.564（偏高 3 %，阻塞效应） | 记录 |
| 圆柱 Re = 100 涡街（`cases/cylinderRe100`） | Williamson 1996 实验 St 0.164；Linnick & Fasel 2005 Cd 1.34、Cl ±0.333、St 0.166；Taira & Colonius 2007 Cd 1.35、Cl ±0.30、St 0.165 | O 型 128×110，域 40D/60D/±40D | t ≥ 100 的 Cd 均值、Cl 幅值、St | Cd 1.335、Cl ±0.334、St 0.1638 | ✅（8 分钟，2 进程） |

## 三、湍流验证

| 算例 | 参考 | 网格 | 指标 | 结果 | 状态 |
| --- | --- | --- | --- | --- | --- |
| 三维 Taylor-Green 涡 Re = 1600（`cases/tgv3d`，图 `docs/validation/tgv3d_re1600.png`） | Dairay 等 2017 DNS 512³（Zenodo 2577239） | 64³ 周期盒，dt 0.01，backward + linear | 总耗散 −dk/dt 峰值与时刻（DNS 0.0129 @ 9.0）、[0,20] 上 −dk/dt 与 k 的相对 L2 误差 | 无模型：0.0104 @ 8.4（−19 %），L2 17 % / 4.8 %；WALE（Cw 0.325）：0.0105 @ 8.5（−18 %），L2 25 % / 6.0 %；Smagorinsky（Cs 0.17）：0.0106 @ 7.4（−18 %，峰值提前 1.6），L2 33 % / 7.7 %；32³ 无模型：0.0100 @ 8.9（−22 %）。与文献中二阶有限体积在 64³ 的表现一致：无模型/隐式 LES 最接近，Smagorinsky 早期过耗散使峰值提前，WALE 介于两者之间 | ✅ ctest（32³ 至 t=12，容差 30 %） |
| RANS 槽道 Re_τ = 392（`cases/channel395`） | Moser 等 DNS | 4×120，y+ ≈ 0.1 | Re_τ、U+ 剖面 | SST 394、k-ω 399、SA 390、BSL-EARSM 385（对数区 U+ 与对数律差 < 0.3） | ✅ |
| RANS 槽道壁面函数 Re_τ ≈ 2000 | DNS 2003 | y+ ≈ 50 | Re_τ | k-ε 1957、SST 1949、realizable 1909 | ✅ |
| WMLES 槽道 Re_τ = 1000（`cases/channelWMLES`） | Lee & Moser 2015 DNS | 126×63×40 | Re_τ、U+ 误差 | 998.2（1000.5）；30 ≤ y+ ≤ 800 内平均 −0.15 %，首层 −9 % | ✅（5.7 小时，4 进程） |
| 周期山 Re_b = 10595 RANS k-ω SST（`cases/periodicHillRANS`，图 `docs/validation/periodic_hill_sst.png`） | Breuer 等 2009 LES（分离 0.19、再附 4.69）/ Rapp 2009 实验（`cases/periodicHill/reference`） | 160×2×80 贴体网格，最大非正交 41°，y+ ≤ 0.45 | 分离/再附点、10 个站位 u 剖面 rms 误差 | 分离 0.25、再附 7.6；u 剖面 rms 误差 0.06–0.15 U_b（平均 0.11）。SST 回流区过长是该模型的已知缺陷：Saïdi 等 2022（J. Comput. Phys. 457, arXiv 2104.09187，表 2）用 k-ω SST 算同一算例得分离 0.270、再附 7.642（其 LES 参考 0.215 / 4.760），与本结果一致 | ✅（与文献 SST 结果定量一致；SIMPLE 2769 步收敛到 1e-6，3.5 分钟） |
| 周期山 Re_b = 10595 RANS BSL-EARSM（`cases/periodicHillEARSM`，图 `docs/validation/periodic_hill_earsm.png`） | 同上；Saïdi 等 2022 的 BSL-EARSM：分离 0.282、再附 4.561 | 同上 | 同上 | 分离 0.24、再附 5.26；u 剖面 rms 误差平均 0.065 U_b（SST 0.111）。模型为 Wallin-Johansson EARSM + Menter BSL，加了 ν_t ≤ k/ω 的上限（整个各向异性张量同比例缩小；不加时弱应变区 C_μ^eff 达 0.33，再附只有 3.2） | ✅（分离流 RANS 的改进方案；再附仍比 LES 长 10 %） |

## 四、并行与可复现

| 项目 | 结果 |
| --- | --- |
| np = 1, 2, 3, 4 可复现模式逐位相同 | TGV、pitzDaily ✅ |
| ctest 并行（np = 3）单元测试 | 5 项 ✅ |

## 已知限制

- 周期方向至少要 2 个单元（1 个单元会通过自己的 cyclic 面耦合到自身，装配出错；现在生成器直接报错）。
- 周期山这类曲面周期域从均匀初场启动时 SST 第一步会发散（均匀场在曲面域内散度很大）；`initial.potentialFlow` 先做无散投影即可，或从静止起动。
- BSL-EARSM 的 ν_t 上限（`nutMax`）是本仓库的处理，不是 Wallin-Johansson 原文；原模型按文献实现时周期山再附 3.2、槽道 Re_τ 408，加上限后 5.3 / 385。

- 蝶形圆管网格的块角单元偏斜，该处误差随加密饱和在约 0.5 %，二阶精度只在正交块内成立（可用更大的 `coreFrac` 或后续加入的偏斜修正改善）。
- 圆柱外流域边界至少取 ±40D，±15D 时 Cd 偏高 3 %。
- 64³ 不加模型的 TGV 耗散峰值偏低 19 %，与其他二阶有限体积代码在同等分辨率下的表现一致（数值耗散不足以替代亚格子耗散）。
