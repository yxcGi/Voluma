#!/usr/bin/env python3
"""周期直圆管层流（Hagen-Poiseuille）结果与解析解对比。

用法: compare_pipe.py <output目录> [--R 0.5] [--Ub 1] [--axis 2] [--nu 0.05 --dpdx-tol 3]
读取 output/line_*_U.csv（沿线采样，给出最近单元中心坐标与速度），按半径 r 与解析解
u(r) = 2 Ub (1 − r²/R²) 比较，打印每条线的最大/均方根相对误差（以中心速度 2 Ub 归一）。
给出 --nu 时另读 line_axis_p.csv，对充分发展段（轴向坐标后 60%，去掉出口附近）线性拟合压力梯度，
与 Hagen-Poiseuille 值 −8 ν Ub/R² 比较（检验扩散项离散的一致性：速度剖面由质量守恒约束，压降才反映粘性通量误差）。
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
    ap.add_argument('--nu', type=float, help='运动粘度：给出时核对轴线压力梯度')
    ap.add_argument('--dpdx-tol', type=float, default=3.0, help='压力梯度相对误差容差（%%）')
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
    fail = a.tol is not None and 100 * worst > a.tol
    if a.nu is not None:
        exact = -8 * a.nu * a.Ub / a.R ** 2
        dpdx = axial_pressure_gradient(os.path.join(a.out, 'line_axis_p.csv'), a.axis)
        err = 100 * (dpdx - exact) / abs(exact)
        print(f'axial pressure gradient {dpdx:.4f}, exact {exact:.4f}, error {err:+.2f}%')
        if a.tol is not None and abs(err) > a.dpdx_tol:
            print(f'FAIL: pressure gradient error exceeds {a.dpdx_tol}%')
            fail = True
    if fail:
        print(f'FAIL: exceeds tolerance {a.tol}%')
        raise SystemExit(1)
    if a.tol is not None:
        print('PASS')


def axial_pressure_gradient(path, axis):
    """充分发展段的轴向压力梯度：取轴向坐标在后 60% 且离末端 10% 以上的采样点线性拟合"""
    with open(path) as f:
        rows = list(csv.DictReader(f))
    pts = sorted({(round(float(r['xyz'[axis]]), 10), float(r['p'])) for r in rows})
    lo, hi = pts[0][0], pts[-1][0]
    sel = [q for q in pts if q[0] >= lo + 0.4 * (hi - lo) and q[0] <= hi - 0.1 * (hi - lo)]
    n = len(sel)
    sx = sum(q[0] for q in sel)
    sy = sum(q[1] for q in sel)
    sxx = sum(q[0] ** 2 for q in sel)
    sxy = sum(q[0] * q[1] for q in sel)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


if __name__ == '__main__':
    main()
