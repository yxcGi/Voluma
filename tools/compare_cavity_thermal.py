#!/usr/bin/env python3
"""差分加热方腔（de Vahl Davis 1983 基准，Pr = 0.71）结果对比。

用法: compare_cavity_thermal.py <输出目录> --ra 1e5 [--log log.txt] [--check] [--tol 3]
读取 log 中热壁的 "wall heat flux hot: ... Nu X" 行，以及 output/line_xMid_U.csv（x = 0.5 竖直中线）、
output/line_yMid_U.csv（y = 0.5 水平中线），速度以 α/L 无量纲（α = 1/sqrt(Ra Pr)）。
基准（de Vahl Davis 1983）：Nu_mean、u_max(x=0.5) 及其 y、v_max(y=0.5) 及其 x。
"""
import argparse
import csv
import math
import os
import re
import sys

REF = {  # Ra: (Nu_mean, u_max, y_at_umax, v_max, x_at_vmax)
    1e3: (1.118, 3.649, 0.813, 3.697, 0.178),
    1e4: (2.243, 16.178, 0.823, 19.617, 0.119),
    1e5: (4.519, 34.73, 0.855, 68.59, 0.066),
    1e6: (8.800, 64.63, 0.850, 219.36, 0.0379),
}


def read_line(path):
    rows = list(csv.reader(open(path)))
    return [[float(x) for x in r] for r in rows[1:] if r]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--ra', type=float, required=True)
    ap.add_argument('--log', default=None)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--tol', type=float, default=3.0, help='Nu 与 u_max/v_max 的相对误差容差（%%）')
    a = ap.parse_args()
    Pr = 0.71
    alpha = 1 / math.sqrt(a.ra * Pr)
    log = a.log or os.path.join(os.path.dirname(os.path.abspath(a.out.rstrip('/'))), 'log.txt')
    text = open(log).read()
    m = re.findall(r'wall heat flux hot: .*?Nu ([-\d.eE+]+)', text)
    mc = re.findall(r'wall heat flux cold: .*?Nu ([-\d.eE+]+)', text)
    if not m:
        print('FAIL: no Nu line for patch hot in', log)
        sys.exit(1)
    nu_hot, nu_cold = float(m[-1]), (float(mc[-1]) if mc else float('nan'))
    xm = read_line(os.path.join(a.out, 'line_xMid_U.csv'))
    ym = read_line(os.path.join(a.out, 'line_yMid_U.csv'))
    u = [(r[1], r[3] / alpha) for r in xm]   # (y, u*)
    v = [(r[0], r[4] / alpha) for r in ym]   # (x, v*)
    y_umax, umax = max(u, key=lambda t: t[1])[0], max(t[1] for t in u)
    x_vmax, vmax = max(v, key=lambda t: t[1])[0], max(t[1] for t in v)
    key = min(REF, key=lambda r: abs(math.log10(r) - math.log10(a.ra)))
    ref = REF[key]
    err_nu = 100 * (nu_hot / ref[0] - 1)
    err_u = 100 * (umax / ref[1] - 1)
    err_v = 100 * (vmax / ref[3] - 1)
    print(f'Ra {a.ra:.0e}: Nu_hot {nu_hot:.4f} (ref {ref[0]}, {err_nu:+.2f} %), Nu_cold {-nu_cold:.4f}')
    print(f'  u_max(x=0.5) {umax:.3f} at y {y_umax:.3f} (ref {ref[1]} at {ref[2]}, {err_u:+.2f} %)')
    print(f'  v_max(y=0.5) {vmax:.3f} at x {x_vmax:.4f} (ref {ref[3]} at {ref[4]}, {err_v:+.2f} %)')
    if a.check:
        bad = [n for n, e in (('Nu', err_nu), ('u_max', err_u), ('v_max', err_v)) if abs(e) > a.tol]
        if bad:
            print('FAIL:', ', '.join(bad), 'outside', a.tol, '%')
            sys.exit(1)
        print('PASS')


if __name__ == '__main__':
    main()
