#!/usr/bin/env bash
#
# WMLES 槽道采样高度扫描：批量运行 cases/channelWMLES/sweep/*.json
#
# 用法（项目根目录，已构建 build/bin/fvmFlow）：
#   tools/run_channel_sweep.sh <输出目录> [每算例进程数=4] [并发算例数=3] [算例名过滤...]
#
# 每个算例在 <输出目录>/<算例名>/ 下运行（case.json 复制过去，output/ 与 log.txt 在其中），
# 同时最多跑"并发算例数"个；全部结束后写 <输出目录>/DONE。进度看 <输出目录>/status.txt。
# 环境变量 SWEEP_SMOKE=1：改成 16×8×40 网格、算到 t=1（只用于检查流程）。
# 可用 MPIRUN 指定启动命令（默认 mpirun）。

set -u
if [ $# -lt 1 ]; then
    echo "usage: $0 <outdir> [np_per_case] [concurrent] [case names...]" >&2
    exit 2
fi
outdir=$(mkdir -p "$1" && cd "$1" && pwd); shift
np=${1:-4}; [ $# -gt 0 ] && shift
conc=${1:-3}; [ $# -gt 0 ] && shift
filter=("$@")

root=$(cd "$(dirname "$0")/.." && pwd)
exe="$root/build/bin/fvmFlow"
mpirun_cmd=${MPIRUN:-mpirun}
[ -x "$exe" ] || { echo "error: $exe not found, build first" >&2; exit 2; }

cases=()
for f in "$root"/cases/channelWMLES/sweep/*.json; do
    name=$(basename "$f" .json)
    if [ ${#filter[@]} -gt 0 ]; then
        keep=0; for p in "${filter[@]}"; do [ "$name" = "$p" ] && keep=1; done
        [ $keep = 1 ] || continue
    fi
    cases+=("$name")
done
[ ${#cases[@]} -gt 0 ] || { echo "error: no cases selected" >&2; exit 2; }

status="$outdir/status.txt"
rm -f "$outdir/DONE"
log_status() { echo "$(date '+%F %T')  $*" | tee -a "$status"; }
log_status "sweep start: ${cases[*]}  (np=$np, concurrent=$conc)"

run_case() {
    local name=$1 dir="$outdir/$1"
    mkdir -p "$dir"
    cp "$root/cases/channelWMLES/sweep/$name.json" "$dir/case.json"
    if [ "${SWEEP_SMOKE:-0}" = 1 ]; then
        sed -i -e 's/"n": \[126, 63, 40\]/"n": [16, 8, 40]/' -e 's/"endTime": 1000/"endTime": 1/' \
               -e 's/"writeInterval": 100/"writeInterval": 0.5/' -e 's/"printInterval": 200/"printInterval": 5/' \
               -e 's/"start": 600/"start": 0.5/' "$dir/case.json"
    fi
    local t0=$(date +%s)
    (cd "$dir" && rm -rf output && $mpirun_cmd -np "$np" "$exe" case.json </dev/null > log.txt 2>&1)
    local rc=$?
    local dt=$(( $(date +%s) - t0 ))
    if [ $rc = 0 ] && [ -f "$dir/output/channel_stats.csv" ]; then
        log_status "$name: done in ${dt}s, $(grep -o 'u_tau [0-9.e+-]*' "$dir/log.txt" | tail -1)"
    else
        log_status "$name: FAILED (exit $rc) after ${dt}s, see $dir/log.txt"
    fi
}

running=0
for name in "${cases[@]}"; do
    if [ $running -ge $conc ]; then
        wait -n
        running=$((running - 1))
    fi
    log_status "$name: start"
    run_case "$name" &
    running=$((running + 1))
done
wait
log_status "sweep finished"
touch "$outdir/DONE"
