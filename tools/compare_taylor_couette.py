#!/usr/bin/env python3
"""Taylor-Couette 结果与解析解 u_θ = A r + B/r 对比（cases/taylorCouette）。

用法: compare_taylor_couette.py <output目录> [--R1 0.5 --R2 1 --omega 1] [--tol 0.01] [--check]
误差以 Ω R1 归一；--check 超过 --tol 时退出码 1。
"""
import argparse
import csv
import glob
import math
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--R1', type=float, default=0.5)
    ap.add_argument('--R2', type=float, default=1.0)
    ap.add_argument('--omega', type=float, default=1.0)
    ap.add_argument('--tol', type=float, default=0.01)
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    A = -a.omega * a.R1 ** 2 / (a.R2 ** 2 - a.R1 ** 2)
    B = a.omega * a.R1 ** 2 * a.R2 ** 2 / (a.R2 ** 2 - a.R1 ** 2)
    uref = a.omega * a.R1
    worst = 0.0
    for path in sorted(glob.glob(os.path.join(a.out, 'line_*_U.csv'))):
        with open(path) as f:
            rows = list(csv.DictReader(f))
        seen, errs = set(), []
        for r in rows:
            x, y = float(r['x']), float(r['y'])
            key = (round(x, 12), round(y, 12))
            if key in seen:
                continue
            seen.add(key)
            rad = math.hypot(x, y)
            ut = (-float(r['U_x']) * y + float(r['U_y']) * x) / rad
            ur = (float(r['U_x']) * x + float(r['U_y']) * y) / rad
            ex = A * rad + B / rad
            errs.append((rad, ut, ex, abs(ut - ex) / uref, abs(ur) / uref))
        mx = max(e[3] for e in errs)
        rms = math.sqrt(sum(e[3] ** 2 for e in errs) / len(errs))
        mr = max(e[4] for e in errs)
        print(f'{os.path.basename(path)}: {len(errs)} points, max |u_t - exact|/ΩR1 = {mx:.3e}, rms {rms:.3e}, max |u_r|/ΩR1 = {mr:.3e}')
        worst = max(worst, mx, mr)
    print('PASS' if worst < a.tol else 'FAIL', f'(worst {worst:.3e}, tol {a.tol})')
    if a.check and worst >= a.tol:
        sys.exit(1)


if __name__ == '__main__':
    main()
