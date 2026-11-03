#!/usr/bin/env python3
"""槽道进出口边界条件算例（cases/channelInlet）与解析解对比。

用法: compare_channel_inlet.py <output目录> --mode poiseuille|developed|bernoulli|ramp [--nu 0.05] [--tol ...] [--check]
  poiseuille  出口剖面 u = 6 U_b y(1−y)（最大误差，以 1.5 U_b 归一）、中心线 dp/dx = −12 ν U_b/h²
  developed   出口中心速度 → 1.5 U_b、剖面与 Poiseuille 的最大误差，进出口流量守恒（读 ../log.txt 的 boundary fluxes）
  bernoulli   总压入口：U = (1,0,0)、p = 0（滑移壁塞流的精确解）
  ramp        时间表入口：t = 0.5 时 U = 0.5、p = 4 − x
--check 时按容差判定，超出退出码 1。
"""
import argparse
import csv
import math
import os
import re
import sys


def read_line(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    seen, out = set(), []
    for r in rows:
        key = (round(float(r['x']), 12), round(float(r['y']), 12))
        if key in seen:
            continue
        seen.add(key)
        out.append({k: float(v) for k, v in r.items()})
    return out


def boundary_fluxes(out):
    log = os.path.join(os.path.dirname(out.rstrip('/')), 'log.txt')
    if not os.path.isfile(log):
        return {}
    last = None
    with open(log) as f:
        for line in f:
            if 'boundary fluxes:' in line:
                last = line
    if not last:
        return {}
    toks = re.sub(r'\(.*?\)', '', last.split('boundary fluxes:')[1]).split()
    return {toks[i]: float(toks[i + 1]) for i in range(0, len(toks) - 1, 2)}


def slope(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--mode', required=True)
    ap.add_argument('--nu', type=float, default=0.05)
    ap.add_argument('--Ub', type=float, default=1.0)
    ap.add_argument('--tol', type=float, default=None, help='最大误差容差（相对）')
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    ok = True
    fl = boundary_fluxes(a.out)
    if fl:
        print('boundary fluxes:', ' '.join(f'{k} {v:.6g}' for k, v in fl.items()))
    if a.mode in ('poiseuille', 'developed'):
        prof = read_line(os.path.join(a.out, 'line_outlet_U.csv'))
        umax = 1.5 * a.Ub
        err = max(abs(r['U_x'] - 6 * a.Ub * r['y'] * (1 - r['y'])) for r in prof) / umax
        uc = max(r['U_x'] for r in prof)
        print(f'outlet profile: max |u - 6Ub y(1-y)| / 1.5Ub = {err:.3e}, max u = {uc:.6f} (1.5 Ub = {umax})')
        cen = read_line(os.path.join(a.out, 'line_centre_p.csv'))
        xs = [r['x'] for r in cen]
        ps = [r['p'] for r in cen]
        # 发展段之后的压力梯度：取后半段
        half = [(x, p) for x, p in zip(xs, ps) if x > 0.5 * max(xs)]
        dpdx = slope([h[0] for h in half], [h[1] for h in half])
        ex = -12 * a.nu * a.Ub
        print(f'dp/dx (second half) = {dpdx:.6f}, analytic {ex:.6f}, rel err {abs(dpdx / ex - 1):.3e}')
        # 壁面单元的半格梯度为一阶近似，剖面与 dp/dx 的误差随网格二阶收敛（80×20：0.21% / 0.37%；160×40：0.06% / 0.09%）
        tol = a.tol if a.tol is not None else (0.005 if a.mode == 'poiseuille' else 0.02)
        ok = err < tol and abs(dpdx / ex - 1) < (0.01 if a.mode == 'poiseuille' else 0.03)
        if fl and 'inlet' in fl and 'outlet' in fl:
            bal = abs(fl['inlet'] + fl['outlet']) / abs(fl['inlet'])
            print(f'mass balance |Q_in + Q_out| / |Q_in| = {bal:.3e}')
            ok = ok and bal < 1e-6
            if a.mode == 'developed':
                ok = ok and abs(abs(fl['inlet']) - 0.1) < 1e-8
    elif a.mode == 'bernoulli':
        cen = read_line(os.path.join(a.out, 'line_centre_U.csv'))
        cp = read_line(os.path.join(a.out, 'line_centre_p.csv'))
        eu = max(max(abs(r['U_x'] - 1), abs(r['U_y'])) for r in cen)
        ep = max(abs(r['p']) for r in cp)
        print(f'Bernoulli plug flow: max |U - (1,0,0)| = {eu:.3e}, max |p| = {ep:.3e}')
        tol = a.tol if a.tol is not None else 1e-5
        ok = eu < tol and ep < tol
    elif a.mode == 'ramp':
        cen = read_line(os.path.join(a.out, 'line_centre_U.csv'))
        cp = read_line(os.path.join(a.out, 'line_centre_p.csv'))
        L = max(r['x'] for r in cp) + 0.05
        eu = max(max(abs(r['U_x'] - 0.5), abs(r['U_y'])) for r in cen)
        ep = max(abs(r['p'] - (L - r['x'])) for r in cp)
        # 入口附近单元值有奇偶振荡（见算例说明），精确性在下游半段与边界通量上检验
        eu2 = max(max(abs(r['U_x'] - 0.5), abs(r['U_y'])) for r in cen if r['x'] > 0.5 * L)
        ep2 = max(abs(r['p'] - (L - r['x'])) for r in cp if r['x'] > 0.5 * L)
        print(f'ramp t=0.5: whole domain max |U - 0.5| = {eu:.3e}, max |p - (L - x)| = {ep:.3e}; '
              f'downstream half {eu2:.3e} / {ep2:.3e}')
        tol = a.tol if a.tol is not None else 1e-5
        ok = eu2 < tol and ep2 < 1e-4
        if fl and 'inlet' in fl and 'outlet' in fl:
            print(f'fluxes: inlet {fl["inlet"]:.8f} (exact -0.05), outlet {fl["outlet"]:.8f}')
            ok = ok and abs(fl['inlet'] + 0.05) < 1e-8 and abs(fl['outlet'] - 0.05) < 1e-8
    else:
        sys.exit('unknown mode ' + a.mode)
    print('PASS' if ok else 'FAIL')
    if a.check and not ok:
        sys.exit(1)


if __name__ == '__main__':
    main()
