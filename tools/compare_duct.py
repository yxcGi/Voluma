#!/usr/bin/env python3
"""方形直管充分发展层流与解析级数解对比（cases/squareDuct）。

用法: compare_duct.py <output目录> [--a 1 --Ub 1 --axis 2 --origin 0 0] [--tol 0.03] [--check]
解析解（边长 a 的方管，坐标以管中心为原点）：
  u(y,z) ∝ Σ_{n odd} (−1)^{(n−1)/2} [1 − cosh(nπz/a)/cosh(nπ/2)] cos(nπy/a) / n³
按数值积分得到的体平均速度归一（u_max/U_b = 2.0962，Shah & London）。读取 output/line_*_U.csv 的各条采样线，
打印以 u_max 归一的最大 / 均方根误差；--check 时超过 --tol 退出码 1。
给出 --nu 时另读 line_axis_p.csv，对充分发展段线性拟合轴向压力梯度，与 Shah & London 的 f_F Re = 14.227
（Fanning 摩擦系数）即 dp/dx = −2·14.227 ν U_b/Dh² 比较（检验扩散项离散的一致性），超过 --dpdx-tol（%）判失败。
"""
import argparse
import csv
import glob
import math
import os
import sys


def series(y, z, a, nterms=60):
    s = 0.0
    for k in range(nterms):
        n = 2 * k + 1
        s += (-1) ** k * (1 - math.cosh(n * math.pi * z / a) / math.cosh(n * math.pi / 2)) * math.cos(n * math.pi * y / a) / n ** 3
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--a', type=float, default=1.0)
    ap.add_argument('--Ub', type=float, default=1.0)
    ap.add_argument('--axis', type=int, default=2, help='轴向坐标分量')
    ap.add_argument('--origin', type=float, nargs=2, default=[0.5, 0.5], help='截面中心在另两个坐标上的位置')
    ap.add_argument('--tol', type=float, default=0.03)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--nu', type=float, help='运动粘度：给出时核对轴线压力梯度')
    ap.add_argument('--dpdx-tol', type=float, default=3.0, help='压力梯度相对误差容差（%%）')
    a = ap.parse_args()
    # 体平均（数值积分，256² 中点）
    n = 256
    tot = 0.0
    for i in range(n):
        y = -a.a / 2 + (i + 0.5) * a.a / n
        for j in range(n):
            z = -a.a / 2 + (j + 0.5) * a.a / n
            tot += series(y, z, a.a)
    mean = tot / (n * n)
    umax = series(0, 0, a.a) / mean * a.Ub
    print(f'analytic u_max/U_b = {umax / a.Ub:.4f} (Shah & London 2.0962)')
    others = [k for k in range(3) if k != a.axis]
    worst = 0.0
    for path in sorted(glob.glob(os.path.join(a.out, 'line_*_U.csv'))):
        if 'axis' in os.path.basename(path):
            continue
        with open(path) as f:
            rows = list(csv.DictReader(f))
        seen, errs = set(), []
        for r in rows:
            x = [float(r['x']), float(r['y']), float(r['z'])]
            key = tuple(round(v, 12) for v in x)
            if key in seen:
                continue
            seen.add(key)
            yy = x[others[0]] - a.origin[0]
            zz = x[others[1]] - a.origin[1]
            ue = series(yy, zz, a.a) / mean * a.Ub
            u = float(r['U_' + 'xyz'[a.axis]])
            errs.append((u - ue) / umax)
        mx = max(abs(e) for e in errs)
        rms = math.sqrt(sum(e * e for e in errs) / len(errs))
        print(f'{os.path.basename(path)}: {len(errs)} points, max |u - exact|/u_max = {mx:.3e}, rms {rms:.3e}')
        worst = max(worst, mx)
    fail = worst >= a.tol
    if a.nu is not None:
        exact = -2 * 14.227 * a.nu * a.Ub / a.a ** 2
        dpdx = axial_pressure_gradient(os.path.join(a.out, 'line_axis_p.csv'), a.axis)
        err = 100 * (dpdx - exact) / abs(exact)
        print(f'axial pressure gradient {dpdx:.4f}, exact {exact:.4f} (f_F Re = 14.227), error {err:+.2f}%')
        if abs(err) > a.dpdx_tol:
            print(f'pressure gradient error exceeds {a.dpdx_tol}%')
            fail = True
    print('FAIL' if fail else 'PASS', f'(worst {worst:.3e}, tol {a.tol})')
    if a.check and fail:
        sys.exit(1)


def axial_pressure_gradient(path, axis):
    """充分发展段的轴向压力梯度：取轴向坐标在后 60% 且离末端 10% 以上的采样点线性拟合"""
    with open(path) as f:
        rows = list(csv.DictReader(f))
    pts = sorted({(round(float(r['xyz'[axis]]), 10), float(r['p'])) for r in rows})
    lo, hi = pts[0][0], pts[-1][0]
    sel = [q for q in pts if lo + 0.4 * (hi - lo) <= q[0] <= hi - 0.1 * (hi - lo)]
    n = len(sel)
    sx = sum(q[0] for q in sel)
    sy = sum(q[1] for q in sel)
    sxx = sum(q[0] ** 2 for q in sel)
    sxy = sum(q[0] * q[1] for q in sel)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


if __name__ == '__main__':
    main()
