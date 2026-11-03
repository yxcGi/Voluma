#!/usr/bin/env python3
"""周期直圆管层流（Hagen-Poiseuille）结果与解析解对比。

用法: compare_pipe.py <output目录> [--R 0.5] [--Ub 1] [--axis 2]
读取 output/line_*_U.csv（沿线采样，给出最近单元中心坐标与速度），按半径 r 与解析解
u(r) = 2 Ub (1 − r²/R²) 比较，打印每条线的最大/均方根相对误差（以中心速度 2 Ub 归一）。
"""
import argparse
import csv
import glob
import math
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--R', type=float, default=0.5)
    ap.add_argument('--Ub', type=float, default=1.0)
    ap.add_argument('--axis', type=int, default=2, help='轴向坐标分量（0/1/2）')
    ap.add_argument('--tol', type=float, help='最大误差容差（%%），超过时退出码 1（回归测试用）')
    a = ap.parse_args()
    umax = 2 * a.Ub
    worst = 0.0
    for path in sorted(glob.glob(os.path.join(a.out, 'line_*_U.csv'))):
        with open(path) as f:
            rows = list(csv.DictReader(f))
        seen = set()
        errs = []
        for r in rows:
            x = [float(r['x']), float(r['y']), float(r['z'])]
            key = tuple(round(v, 12) for v in x)
            if key in seen:
                continue
            seen.add(key)
            rad = math.sqrt(sum(x[k] ** 2 for k in range(3) if k != a.axis))
            u = float(r['U_' + 'xyz'[a.axis]])
            ue = umax * max(0.0, 1 - (rad / a.R) ** 2)
            errs.append((rad, u, ue, (u - ue) / umax))
        rms = math.sqrt(sum(e[3] ** 2 for e in errs) / len(errs))
        mx = max(abs(e[3]) for e in errs)
        worst = max(worst, mx)
        centre = max(errs, key=lambda e: -e[0])
        print(f'{os.path.basename(path)}: {len(errs)} cells, max |err| {100 * mx:.3f}%  rms {100 * rms:.3f}%  '
              f'(r={centre[0]:.3f}: u={centre[1]:.4f}, exact {centre[2]:.4f})')
    print(f'worst max error {100 * worst:.3f}% of centreline velocity')
    if a.tol is not None and 100 * worst > a.tol:
        print(f'FAIL: exceeds tolerance {a.tol}%')
        raise SystemExit(1)
    if a.tol is not None:
        print('PASS')


if __name__ == '__main__':
    main()
