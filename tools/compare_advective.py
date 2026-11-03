#!/usr/bin/env python3
"""对流出流算例（cases/advectiveOutlet）与一维对流扩散解析解对比。

用法: compare_advective.py <output目录> [--t 1.8 --x0 0.5 --U 1 --sigma 0.1 --gamma 0.01] [--tol 0.005] [--check]
解析解 T = σ/σ(t) exp(−(x − x0 − U t)²/σ(t)²)，σ(t)² = σ² + 4Γt。打印域内最大误差（以峰值归一）与出口附近误差。
"""
import argparse
import csv
import math
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--t', type=float, default=1.8)
    ap.add_argument('--x0', type=float, default=0.5)
    ap.add_argument('--U', type=float, default=1.0)
    ap.add_argument('--sigma', type=float, default=0.1)
    ap.add_argument('--gamma', type=float, default=0.01)
    ap.add_argument('--tol', type=float, default=0.005)
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    with open(os.path.join(a.out, 'line_centre_T.csv')) as f:
        rows = list(csv.DictReader(f))
    s2 = a.sigma ** 2 + 4 * a.gamma * a.t
    peak = a.sigma / math.sqrt(s2)
    xc = a.x0 + a.U * a.t
    seen, err, errNear = set(), 0.0, 0.0
    for r in rows:
        x = float(r['x'])
        if round(x, 12) in seen:
            continue
        seen.add(round(x, 12))
        ex = peak * math.exp(-(x - xc) ** 2 / s2)
        e = abs(float(r['T']) - ex) / peak
        err = max(err, e)
        if x > xc - 2 * math.sqrt(s2):
            errNear = max(errNear, e)
    print(f't = {a.t}: pulse centre {xc}, analytic peak {peak:.4f}; max |T - exact|/peak = {err:.3e} (within 2σ of outlet: {errNear:.3e})')
    print('PASS' if err < a.tol else 'FAIL')
    if a.check and err >= a.tol:
        sys.exit(1)


if __name__ == '__main__':
    main()
