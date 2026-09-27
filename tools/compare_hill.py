#!/usr/bin/env python3
"""周期山结果与 ERCOFTAC UFR 3-30 实验（Rapp 2009，PIV，Re = 10600）对比。

用法: compare_hill.py output/hill --ref <实验数据目录> [--out fig.png] [--label name] [more dirs ...]

output/hill 为 fvmFlow spanStatistics 输出（或 OpenLB periodichill3d 的同名 CSV）：
  profile_xh_<站>.csv  y_h;u_Ub;v_Ub;uu_Ub2;vv_Ub2;ww_Ub2;uv_Ub2;k_Ub2（y_h 为离下壁高度）
  wall_bottom.csv      x_h;Cf;Cp;ut_Ub
实验数据为 kbwiki.ercoftac.org 上 UFR3-30_X_10600_data_CR-0NN.dat（NN = 1..10 对应 x/h = 0.05, 0.5, 1..8），
列 y/h; u/u_b; v/u_b; u'u'; v'v'; u'v'（y/h 为绝对高度，山体内为 0）。
"""
import argparse
import glob
import os

import numpy as np

STATIONS = [0.05, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]


def hill_mm(x):
    if x < 9:
        return min(28.0, 28 + 6.775070969851e-03 * x * x - 2.124527775800e-03 * x ** 3)
    if x < 14:
        return 2.507355893131e+01 + 9.754803562315e-01 * x - 1.016116352781e-01 * x * x + 1.889794677828e-03 * x ** 3
    if x < 20:
        return 2.579601052357e+01 + 8.206693007457e-01 * x - 9.055370274339e-02 * x * x + 1.626510569859e-03 * x ** 3
    if x < 30:
        return 4.046435022819e+01 - 1.379581654948e+00 * x + 1.945884504128e-02 * x * x - 2.070318932190e-04 * x ** 3
    if x < 40:
        return 1.792461334664e+01 + 8.743920332081e-01 * x - 5.567361123058e-02 * x * x + 6.277731764683e-04 * x ** 3
    if x <= 54:
        return max(0.0, 5.639011190988e+01 - 2.010520359035e+00 * x + 1.644919857549e-02 * x * x
                   + 2.674976141766e-05 * x ** 3)
    return 0.0


def hill(x):
    x = x % 9.0
    xl, xr = x * 28, (9 - x) * 28
    return (hill_mm(xl) if xl <= 54 else hill_mm(xr) if xr <= 54 else 0.0) / 28


def read_semicolon(path):
    with open(path) as f:
        head = f.readline().strip().split(';')
        data = np.array([[float(v) for v in l.split(';')] for l in f if l.strip()])
    return {h: data[:, i] for i, h in enumerate(head)}


def read_exp(path):
    rows = [l for l in open(path) if l.strip() and not l.startswith('#')]
    d = np.array([[float(v) for v in l.split(',')] for l in rows])
    keep = np.any(d[:, 1:] != 0, axis=1)  # 山体内的点全为 0
    d = d[keep]
    return {'y': d[:, 0], 'u': d[:, 1], 'v': d[:, 2], 'uu': d[:, 3], 'vv': d[:, 4], 'uv': d[:, 5]}


def sep_reat(w):
    x, cf = w['x_h'], w['Cf']
    s = r = None
    for i in range(1, len(x)):
        if x[i] < 0.1:
            continue
        xz = x[i - 1] + (x[i] - x[i - 1]) * cf[i - 1] / (cf[i - 1] - cf[i]) if cf[i] != cf[i - 1] else x[i]
        if s is None and cf[i] < 0 <= cf[i - 1]:
            s = xz
        elif s is not None and r is None and cf[i] > 0 >= cf[i - 1]:
            r = xz
    return s, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out')
    a = ap.parse_args()
    exp = {}
    for k, st in enumerate(STATIONS):
        f = glob.glob(os.path.join(a.ref, f'*_CR-{k + 1:03d}.dat'))
        if f:
            exp[st] = read_exp(f[0])
    runs = []
    for d in a.dirs:
        prof = {}
        for st in STATIONS:
            p = os.path.join(d, f'profile_xh_{st:g}.csv')
            if os.path.exists(p):
                prof[st] = read_semicolon(p)
        w = read_semicolon(os.path.join(d, 'wall_bottom.csv'))
        s, r = sep_reat(w)
        print(f'{d}: separation x/h {s}, reattachment x/h {r}  (LES literature: ~0.2 / 4.56-4.72)')
        errs = []
        for st, p in prof.items():
            if st not in exp:
                continue
            e = exp[st]
            y = p['y_h'] + hill(st)
            um = np.interp(e['y'], y, p['u_Ub'], left=np.nan, right=np.nan)
            ok = ~np.isnan(um)
            rms = np.sqrt(np.mean((um[ok] - e['u'][ok]) ** 2))
            errs.append(rms)
            print(f'  x/h {st:4g}: rms(u - u_exp)/U_b = {rms:.3f}')
        if errs:
            print(f'  mean rms over stations: {np.mean(errs):.3f}')
        runs.append((d, prof, w))
    if a.out:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        xs = np.linspace(0, 9, 400)
        fig, ax = plt.subplots(3, 1, figsize=(11, 11))
        scale = {'u': 1.0, 'uu': 10.0, 'uv': -20.0}
        for i, (key, col, lab) in enumerate([('u', 'u_Ub', 'x/h + u/U_b'), ('uu', 'uu_Ub2', "x/h + 10 u'u'/U_b²"),
                                             ('uv', 'uv_Ub2', "x/h − 20 u'v'/U_b²")]):
            axx = ax[i]
            axx.fill_between(xs, 0, [hill(x) for x in xs], color='0.8')
            axx.plot([0, 9], [3.036, 3.036], 'k-', lw=1)
            for st in STATIONS:
                if st in exp:
                    e = exp[st]
                    axx.plot(st + scale[key] * e[key], e['y'], 'k.', ms=2, label='Rapp 2009 exp.' if st == 0.05 else None)
                for j, (d, prof, w) in enumerate(runs):
                    if st in prof:
                        p = prof[st]
                        axx.plot(st + scale[key] * p[col], p['y_h'] + hill(st), f'C{j}-', lw=1.2,
                                 label=d if st == 0.05 else None)
            axx.set_xlim(-0.2, 9.5), axx.set_ylim(0, 3.1), axx.set_ylabel('y/h'), axx.set_title(lab)
            axx.legend(fontsize=8, loc='upper right')
        ax[2].set_xlabel('x/h')
        fig.tight_layout()
        fig.savefig(a.out, dpi=130)
        fig2, ax2 = plt.subplots(figsize=(9, 3.5))
        for j, (d, prof, w) in enumerate(runs):
            ax2.plot(w['x_h'], w['Cf'], f'C{j}-', label=d)
        ax2.axhline(0, color='k', lw=0.5)
        ax2.set_xlabel('x/h'), ax2.set_ylabel('Cf (bottom wall)'), ax2.legend(fontsize=8)
        fig2.tight_layout()
        fig2.savefig(os.path.splitext(a.out)[0] + '_Cf.png', dpi=130)


if __name__ == '__main__':
    main()
