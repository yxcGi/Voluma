#!/usr/bin/env python3
"""二维层流圆柱绕流结果与文献对比（Re 40 定常、Re 100 卡门涡街）。

用法: compare_cylinder_laminar.py <output目录> --re 40|100 [--tstart 100]
Re 40：forces.csv 末行 Cd 及其压力/粘性分量，line_wake_U.csv 得回流区长度 L_w/D（U_x 过零点 − 0.5），
       surface.csv 由壁面切应力的切向分量过零得分离角 θ_sep（自后驻点量起）。
Re 100：t ≥ tstart 的 Cd 均值、Cl 幅值（max−min 的一半）、由 Cl 过零得到 Strouhal 数。
"""
import argparse
import csv
import math
import os

import numpy as np


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) if k != 'patch' else 0 for r in rows]) for k in rows[0].keys()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--re', type=int, required=True)
    ap.add_argument('--tstart', type=float, default=100.0)
    ap.add_argument('--check', action='store_true', help='按文献范围（放宽）判定，失败退出码 1（回归测试用）')
    a = ap.parse_args()
    F = read_csv(os.path.join(a.out, 'forces.csv'))
    if a.re == 40:
        print(f"Cd {F['Cd'][-1]:.4f} (pressure {F['Cd_p'][-1]:.4f}, viscous {F['Cd_v'][-1]:.4f}), Cl {F['Cl'][-1]:.1e}")
        print('  reference: Dennis&Chang 1970 Cd 1.522, Fornberg 1980 1.498, Linnick&Fasel 2005 1.54, Taira&Colonius 2007 1.54')
        W = read_csv(os.path.join(a.out, 'line_wake_U.csv'))
        x, u = W['x'], W['U_x']
        # 去重（最近单元采样）并按 x 排序
        xs, idx = np.unique(x, return_index=True)
        us = u[idx]
        i = np.where((us[:-1] < 0) & (us[1:] >= 0))[0]
        if len(i):
            k = i[0]
            xz = xs[k] + (0 - us[k]) * (xs[k + 1] - xs[k]) / (us[k + 1] - us[k])
            print(f"wake length L_w/D {xz - 0.5:.3f}   (reference 2.24-2.35)")
        else:
            print('wake length: no reattachment found on the centreline sample')
        S = read_csv(os.path.join(a.out, 'surface.csv'))
        th = np.degrees(np.arctan2(S['y'], S['x']))  # 自前驻点（−x）… 用 x 轴角：前驻点 180°，后驻点 0°
        # 切向单位向量（逆时针）与切应力的内积
        tx, ty = -np.sin(np.radians(th)), np.cos(np.radians(th))
        taut = S['tau_x'] * tx + S['tau_y'] * ty
        upper = th > 0
        o = np.argsort(th[upper])
        thu, tu = th[upper][o], taut[upper][o]
        # 上表面从后驻点(0°)向前驻点(180°)：分离点处切应力符号改变
        j = np.where(np.sign(tu[:-1]) != np.sign(tu[1:]))[0]
        if len(j):
            k = j[0]
            ths = thu[k] + (0 - tu[k]) * (thu[k + 1] - thu[k]) / (tu[k + 1] - tu[k])
            print(f"separation angle from rear stagnation point {ths:.1f} deg   (reference 53.6-53.7)")
        if a.check:
            ok = 1.49 <= F['Cd'][-1] <= 1.56 and len(i) and 2.1 <= xz - 0.5 <= 2.4 and len(j) and 51 <= ths <= 56
            print('PASS' if ok else 'FAIL: outside Cd 1.49-1.56, L_w 2.1-2.4, theta 51-56')
            if not ok:
                raise SystemExit(1)
    else:
        sel = F['t'] >= a.tstart
        t, cd, cl = F['t'][sel], F['Cd'][sel], F['Cl'][sel]
        cdm = np.mean(cd)
        clamp = 0.5 * (np.max(cl) - np.min(cl))
        # Cl 向上过零时刻 → 周期
        z = np.where((cl[:-1] < 0) & (cl[1:] >= 0))[0]
        tz = t[z] + (0 - cl[z]) * (t[z + 1] - t[z]) / (cl[z + 1] - cl[z])
        st = 1.0 / np.mean(np.diff(tz)) if len(tz) > 2 else float('nan')
        print(f"t >= {a.tstart}: Cd mean {cdm:.4f} (±{0.5 * (np.max(cd) - np.min(cd)):.4f}), Cl amplitude ±{clamp:.4f}, "
              f"St {st:.4f} ({len(tz) - 1} periods)")
        print('  reference: Williamson 1996 exp St 0.164; Linnick&Fasel 2005 Cd 1.34, Cl ±0.333, St 0.166; '
              'Taira&Colonius 2007 Cd 1.35, Cl ±0.30, St 0.165')
        if a.check:
            ok = 1.30 <= cdm <= 1.40 and 0.25 <= clamp <= 0.38 and 0.158 <= st <= 0.172
            print('PASS' if ok else 'FAIL: outside Cd 1.30-1.40, Cl 0.25-0.38, St 0.158-0.172')
            if not ok:
                raise SystemExit(1)


if __name__ == '__main__':
    main()
