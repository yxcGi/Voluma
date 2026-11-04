# 需要在大机器上跑的算例（交接清单）

云端容器只有 4 核，下面三组算例每组都要几小时到十几小时，留给本地 / 集群机器。全部用仓库里现成的算例文件与脚本，
跑完把列出的结果文件回传即可（结果判读与对表脚本已经写好）。先在仓库根目录构建：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
```

所有瞬态算例都能用 `tools/run_chunked.py` 分段续算：每段结束写续算文件，中断后重新执行同一命令即可从最后完成的段接着跑，
统计量（channelStatistics、spanStatistics、fieldAverage）随续算累积，分段结果与一次跑完一致。

## 1. 周期山 WMLES（Re_b = 10595，Breuer 等 2009 / Rapp 2009 实验）

| 项 | 内容 |
| --- | --- |
| 算例 | `cases/periodicHill/case.json`：128×64×64（52 万单元），WALE + 幂律壁模型（采样高度 0.11H），U_b 由体力保持 |
| 时长 | t = 360（T_f = 9H/U_b ≈ 12.5，约 29 个流通时间），t ≥ 90 做展向 + 时间统计 |
| 命令 | `tools/run_chunked.py cases/periodicHill/case.json --workdir runs/hill_wmles --np 16 --chunk 45` |
| 估算 | 4 核云端每步约 0.4 s，36000 步；16 进程估计 1–2 小时，8 进程 2–4 小时（粗估，看实际步数） |
| 对比 | `python3 tools/compare_hill.py runs/hill_wmles/output/hill --ref cases/periodicHill/reference --out hill_wmles.png` |
| 回传 | `runs/hill_wmles/output/hill/`（profile_xh_*.csv、wall_bottom.csv、fields.csv）、`log.txt`、`status.txt` |

关注点：分离 / 再附点（实验 0.19 / 4.69，RANS SST 0.25 / 7.6，EARSM 0.24 / 5.26，见 `docs/validation.md` 三），
10 个站位的 u、u'u'、u'v' 剖面。可以用 `--set turbulence.wallModel.type=ode --set turbulence.wallModel.pressureGradient=true`
派生非平衡壁模型的变体（`runs/hill_wmles_odepg`），这是 `ode` 含压力梯度版本的真正用途（槽道里它偏高 1–4 %，分离流才需要）。

## 2. 槽道 WMLES 采样高度扫描（Re_τ = 1000，Lee & Moser 2015 DNS）

| 项 | 内容 |
| --- | --- |
| 算例 | `cases/channelWMLES/sweep/*.json`：126×63×40 网格，Smagorinsky / WALE × 采样高度 1.5 / 2.5 / 3.5 层，t = 1000（600 起统计） |
| 命令 | `tools/run_channel_sweep.sh runs/sweep 4 3`（每例 4 进程、同时 3 例；`SWEEP_SMOKE=1` 先走一遍流程） |
| 估算 | 每例 4 进程约 5.7 小时（云端实测 Spalding 基准例），6 例并发 3 个约 12 小时 |
| 汇总 | `tools/collect_channel_sweep.sh runs/sweep > sweep_report.txt`；逐例 `python3 tools/compare_channel.py runs/sweep/<例>/output/channel_stats.csv --ref cases/channelWMLES/reference/referenceData_1000.csv` |
| 回传 | 各例 `output/channel_stats.csv`、`log.txt`，`runs/sweep/status.txt` |

关注点：Re_τ（DNS 1000.5）与 30 ≤ y+ ≤ 800 内 U+ 的平均误差；采样高度对首层 U+ 偏低（基准例 −9 %）的影响。

## 3. 槽道 WMLES 壁面律对比（完整分辩率）

粗网格（48×24×32，首层 y+ ≈ 31，`cases/channelWMLES/wallLaws/`）的结果已在云端做完，见 `docs/validation.md` 三：
spalding / reichardt / ode / odepg 的 Re_τ 与 U+ 偏差。完整分辩率的对比从 `cases/channelWMLES/case.json`（Spalding，已有结果
Re_τ 998.2）派生：

```bash
for wm in musker reichardt logLaw powerLaw ode; do
  tools/run_chunked.py cases/channelWMLES/case.json --workdir runs/wl_$wm --np 4 --chunk 100 \
      --set turbulence.wallModel.type=$wm &
done
tools/run_chunked.py cases/channelWMLES/case.json --workdir runs/wl_odepg --np 4 --chunk 100 \
    --set turbulence.wallModel.type=ode --set turbulence.wallModel.pressureGradient=true &
wait
```

每例 4 进程约 5.7 小时；对比与回传同第 2 组（`compare_channel.py` 可一次给多个 csv 画在一张图上）。

## 结果回传后要做的事

1. 把三组结果的 Re_τ / U+ 误差 / 分离再附点填进 `docs/validation.md`（三、湍流验证）对应行，WMLES 周期山新增一行；
2. 若周期山 WMLES 的再附点仍偏长，优先检查采样高度（0.11H 落在哪一层）与 `odepg` 变体的差别；
3. `docs/local_runs.md` 本文件随之更新或删除。
