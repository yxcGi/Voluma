#!/usr/bin/env python3
"""圆柱绕流（Re_D = 3900）结果与 Kravchenko & Moin (2000) Table II 对比。

用法: compare_cylinder.py <output 目录> [--start 50] [--out fig.png]

读取 fvmFlow 输出：forces.csv（Cd、Cl 时间序列）、stats/wall_cylinder.csv、stats/fields.csv
（spanStatistics，D = 1、U = 1）。计算：
  Cd（t ≥ start 平均）、St（Cl 上穿零点的平均周期）、−Cpb（θ = 180° 处，Cp 以驻点 Cp = 1 归一，
  即 p∞ = p_stag − ½U²）、Lr/D（尾迹中心线平均 u 由负变正处到圆柱表面 x = 0.5 的距离）、
  θsep（上表面 Cf 由正变负处，从前驻点量起）。
"""
import argparse
import csv
import os

import numpy as np

KM_EXP = {'Cd': (0.99, 0.05), '-Cpb': (0.88, 0.05), 'St': (0.215, 0.005), 'Lr/D': (1.4, 0.1), 'theta_sep': (86, None)}
KM_LES = {'Cd': 1.04, '-Cpb': 0.94, 'St': 0.210, 'Lr/D': 1.35, 'theta_sep': 88}


def read(path, delim):
    with open(path) as f:
        r = list(csv.reader(f, delimiter=delim))
    head = [h.strip() for h in r[0]]
    d = np.array([[float(x) for x in row] for row in r[1:] if row])
    return {h: d[:, i] for i, h in enumerate(head)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--start', type=float, default=50)
    ap.add_argument('--fig')
    a = ap.parse_args()
    res = {}
    f = read(os.path.join(a.out, 'forces.csv'), ',')
    m = f['t'] >= a.start
    t, cd, cl = f['t'][m], f['Cd'][m], f['Cl'][m]
    # 时间加权平均（步长可变）
    w = np.gradient(t)
    res['Cd'] = np.sum(cd * w) / np.sum(w)
    up = [t[i - 1] + (t[i] - t[i - 1]) * (-cl[i - 1]) / (cl[i] - cl[i - 1]) for i in range(1, len(cl)) if cl[i - 1] < 0 <= cl[i]]
    res['St'] = 1.0 / np.mean(np.diff(up)) if len(up) > 2 else float('nan')
    res['Cl_rms'] = np.sqrt(np.sum((cl - np.sum(cl * w) / np.sum(w)) ** 2 * w) / np.sum(w))
    wall = read(os.path.join(a.out, 'stats', 'wall_cylinder.csv'), ';')
    x, y = wall['x_h'], wall['y_h']
    th = np.degrees(np.arctan2(y, -x))  # 从前驻点量起：前驻点 0°，上表面 0..180°
    cp = wall['Cp']
    i0 = np.argmin(np.abs(th))
    cp = cp - cp[i0] + 1.0  # 驻点 Cp = 1
    base = np.abs(np.abs(th) - 180) < 360 / len(th) * 1.01
    res['-Cpb'] = -np.mean(cp[base])
    top = th > 0
    o = np.argsort(th[top])
    tt, cf = th[top][o], wall['Cf'][top][o]
    sep = [tt[i - 1] + (tt[i] - tt[i - 1]) * cf[i - 1] / (cf[i - 1] - cf[i]) for i in range(1, len(cf))
           if cf[i - 1] > 0 >= cf[i] and tt[i] > 30]
    res['theta_sep'] = sep[0] if sep else float('nan')
    fl = read(os.path.join(a.out, 'stats', 'fields.csv'), ';')
    # 中心线：|y| 最小的两排单元插值到 y = 0（O 网格中心线两侧对称）
    xs, ys, us = fl['x_h'], fl['y_h'], fl['u_Ub']
    sel = (xs > 0.5) & (xs < 8) & (np.abs(ys) < 0.1)
    xc, uc = xs[sel], us[sel]
    o = np.argsort(xc)
    xc, uc = xc[o], uc[o]
    # 按 x 分箱平均
    bins = np.linspace(0.5, 8, 301)
    idx = np.digitize(xc, bins)
    xb = np.array([xc[idx == k].mean() for k in np.unique(idx)])
    ub = np.array([uc[idx == k].mean() for k in np.unique(idx)])
    cross = [xb[i - 1] + (xb[i] - xb[i - 1]) * (-ub[i - 1]) / (ub[i] - ub[i - 1]) for i in range(1, len(ub))
             if ub[i - 1] < 0 <= ub[i]]
    res['Lr/D'] = cross[0] - 0.5 if cross else float('nan')
    res['u_min'] = ub.min() if len(ub) else float('nan')
    print(f'{"":10s} {"本代码":>8s} {"K&M 实验":>14s} {"K&M LES":>8s}')
    for k in ('Cd', '-Cpb', 'St', 'Lr/D', 'theta_sep'):
        e, de = KM_EXP[k]
        print(f'{k:10s} {res[k]:8.3f} {e:8.3f}' + (f' ±{de:<5.3f}' if de else ' ' * 7) + f' {KM_LES[k]:8.3f}')
    print(f'Cl rms {res["Cl_rms"]:.3f}, centreline min <u> {res["u_min"]:.3f}, averaging time {t[-1] - t[0]:.1f}')
    if a.fig:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
        o = np.argsort(th[top])
        ax[0].plot(th[top][o], cp[top][o], '-')
        ax[0].set_xlabel('θ (deg)'), ax[0].set_ylabel('Cp'), ax[0].set_title('mean surface pressure')
        ax[1].plot(xb, ub, '-'), ax[1].axhline(0, color='k', lw=0.5)
        ax[1].set_xlabel('x/D'), ax[1].set_ylabel('<u>/U'), ax[1].set_title('wake centreline')
        ax[2].plot(t, cd, label='Cd'), ax[2].plot(t, cl, label='Cl')
        ax[2].set_xlabel('t U/D'), ax[2].legend()
        fig.tight_layout()
        fig.savefig(a.fig, dpi=130)


if __name__ == '__main__':
    main()
