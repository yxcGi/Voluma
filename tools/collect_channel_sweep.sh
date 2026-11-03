#!/usr/bin/env bash
#
# 汇总 run_channel_sweep.sh 的结果为一份纯文本报告（无 python 依赖）：
#   tools/collect_channel_sweep.sh <输出目录> [--full]
# 每个算例：状态、u_τ、Re_τ（= u_τ H/ν，H=1、ν=5e-5）、日志末尾；--full 时附完整 channel_stats.csv。
# 完整对比（与 Lee & Moser、OpenLB 结果）用 tools/compare_channel.py。

set -u
outdir=${1:?usage: $0 <outdir> [--full]}
full=${2:-}
echo "# channel WMLES sweep report  $(date -Is)"
echo "## status.txt"; cat "$outdir/status.txt" 2>/dev/null
echo
echo "## summary"
printf '%-12s %-8s %-12s %-10s %s\n' case state u_tau Re_tau averaged
for d in "$outdir"/*/; do
    name=$(basename "$d")
    csv="$d/output/channel_stats.csv"
    if [ -f "$csv" ]; then
        ut=$(sed -n 2p "$csv" | cut -d';' -f2)
        retau=$(awk "BEGIN { printf \"%.1f\", $ut / 5e-5 }")
        avg=$(grep -o 'averaged over [0-9.e+-]*' "$d/log.txt" | tail -1 | cut -d' ' -f3)
        grep -q 'reached t' "$d/log.txt" && st=done || st=partial
        printf '%-12s %-8s %-12s %-10s %s\n' "$name" "$st" "$ut" "$retau" "$avg"
    elif [ -f "$d/log.txt" ]; then
        printf '%-12s %-8s %s\n' "$name" running "$(grep '^t ' "$d/log.txt" | tail -1)"
    fi
done
for d in "$outdir"/*/; do
    name=$(basename "$d")
    [ -f "$d/log.txt" ] || continue
    echo; echo "## $name  (log tail)"; tail -4 "$d/log.txt"
    if [ "$full" = --full ] && [ -f "$d/output/channel_stats.csv" ]; then
        echo "--- channel_stats.csv"; cat "$d/output/channel_stats.csv"
    fi
done
