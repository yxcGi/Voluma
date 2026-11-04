#!/usr/bin/env python3
"""斜向顶盖（45°）三维方腔 Re = 1000 中线结果与文献对比。

用法: compare_cavity3d.py <output目录> [--line vert] [--tol 0.03] [--check]
读取 output/line_<line>_U.csv 与 line_<line>_p.csv（竖直中线 (0.5, y, 0.5)），线性插值到文献表格的 y 位置，
与 Feldman & Gelfgat 2010（152³ 网格）及 Feldman 2013（OpenFOAM icoFoam 100³）的 v_x(=v_z)、v_y、p 对比
（Feldman, arXiv:1306.3267, Table 1；p 为相对腔中心 (0.5, 0.5, 0.5) 的值，用 ρU² 归一）。
顶盖速度 U = 1 沿 (1, 0, 1)/√2，Re = U L/ν = 1000。
打印各量的最大偏差（以该量在表中的最大绝对值归一）；--check 时任一量超过 --tol 退出码 1。
"""
import argparse
import csv
import os

# y, vx(=vz) ref, vx pres, vy ref, vy pres, p ref, p pres   （表中 v×1e3，p×1e4，这里已换回原值）
TABLE = [
    (0.9766, 417.7, 417.8, 5.378, 5.458, 51.59, 51.42),
    (0.9531, 226.6, 226.5, 16.07, 16.14, 46.67, 46.60),
    (0.8516, 76.74, 76.39, 30.36, 30.36, 35.18, 35.00),
    (0.7344, 62.50, 62.00, 22.59, 22.50, 21.85, 21.73),
    (0.6172, 41.78, 41.22, 5.790, 5.561, 9.711, 9.664),
    (0.5000, -1.398, -1.395, -33.95, -34.08, 0.000, 0.000),
    (0.4531, -31.54, -31.33, -64.70, -64.35, -2.517, -2.4525),
    (0.2813, -130.7, -129.0, -160.2, -158.0, 28.24, 28.03),
    (0.1719, -134.7, -133.6, -137.9, -135.9, 107.0, 106.0),
    (0.1016, -143.1, -142.5, -86.78, -85.41, 186.5, 184.5),
    (0.0547, -162.2, -161.1, -35.52, -34.49, 232.9, 230.6),
]


def read_line(path, cols):
    """读取采样表，按 y 去重排序，返回 [(y, 值...)]"""
    rows = list(csv.DictReader(open(path)))
    seen = {}
    for r in rows:
        y = float(r['y'])
        seen[round(y, 10)] = (y,) + tuple(float(r[c]) for c in cols)
    return sorted(seen.values())


def interp(data, y, k):
    """线性插值第 k 列到 y"""
    for (y0, *v0), (y1, *v1) in zip(data, data[1:]):
        if y0 <= y <= y1:
            w = 0 if y1 == y0 else (y - y0) / (y1 - y0)
            return v0[k] * (1 - w) + v1[k] * w
    raise ValueError(f'y={y} outside sampled range [{data[0][0]}, {data[-1][0]}]')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--line', default='vert')
    ap.add_argument('--tol', type=float, default=0.03)
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    U = read_line(os.path.join(a.out, f'line_{a.line}_U.csv'), ['U_x', 'U_y', 'U_z'])
    P = read_line(os.path.join(a.out, f'line_{a.line}_p.csv'), ['p'])
    p0 = interp(P, 0.5, 0)
    print(f'{"y":>7} {"vx":>9} {"vx ref":>9} {"vx FOAM":>9} | {"vy":>9} {"vy ref":>9} {"vy FOAM":>9} | '
          f'{"p":>9} {"p ref":>9} {"p FOAM":>9}   (v×1e3, p×1e4; ref = Feldman & Gelfgat 2010 152³)')
    dev = {'vx': 0.0, 'vy': 0.0, 'p': 0.0}
    scale = {'vx': max(abs(r[1]) for r in TABLE), 'vy': max(abs(r[3]) for r in TABLE), 'p': max(abs(r[5]) for r in TABLE)}
    vz_asym = 0.0
    for y, vxr, vxp, vyr, vyp, pr, pp in TABLE:
        vx = 1e3 * interp(U, y, 0)
        vz = 1e3 * interp(U, y, 2)
        vy = 1e3 * interp(U, y, 1)
        p = 1e4 * (interp(P, y, 0) - p0)
        vz_asym = max(vz_asym, abs(vx - vz))
        dev['vx'] = max(dev['vx'], abs(vx - vxr) / scale['vx'])
        dev['vy'] = max(dev['vy'], abs(vy - vyr) / scale['vy'])
        dev['p'] = max(dev['p'], abs(p - pr) / scale['p'])
        print(f'{y:7.4f} {vx:9.2f} {vxr:9.2f} {vxp:9.2f} | {vy:9.2f} {vyr:9.2f} {vyp:9.2f} | {p:9.2f} {pr:9.2f} {pp:9.2f}')
    print(f'max deviation from Feldman & Gelfgat 2010 (normalised by column max): '
          f'vx {100 * dev["vx"]:.2f}%  vy {100 * dev["vy"]:.2f}%  p {100 * dev["p"]:.2f}%   '
          f'(icoFoam 100³ vs ref: vx {100 * max(abs(r[2] - r[1]) for r in TABLE) / scale["vx"]:.2f}%  '
          f'vy {100 * max(abs(r[4] - r[3]) for r in TABLE) / scale["vy"]:.2f}%  '
          f'p {100 * max(abs(r[6] - r[5]) for r in TABLE) / scale["p"]:.2f}%)')
    print(f'max |vx - vz| along the line (diagonal symmetry): {vz_asym:.3f} ×1e-3')
    worst = max(dev.values())
    if a.check:
        if worst > a.tol:
            print(f'FAIL: worst deviation {100 * worst:.2f}% > tol {100 * a.tol:.2f}%')
            raise SystemExit(1)
        print(f'PASS (worst {100 * worst:.2f}%, tol {100 * a.tol:.2f}%)')


if __name__ == '__main__':
    main()
