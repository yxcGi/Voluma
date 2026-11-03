#!/usr/bin/env python3
"""三维 Taylor-Green 涡（Re = 1600）与 DNS 参考（Dairay 等 2017，cases/tgv3d/reference/TGV_Re1600.dat）对比。

用法: compare_tgv.py monitor.csv [more.csv ...] [--ref TGV_Re1600.dat] [--out fig.png] [--labels a,b,...]

monitor.csv 为 fvmFlow 的 "monitor" 输出（t,k,epsilon,epsilon_sgs,enstrophy）。
对每个文件打印：总耗散 −dk/dt 的峰值与峰值时刻（DNS 约 0.0126 @ t≈9.1）、解析粘性耗散峰值、
以及 t ∈ [0, 20] 上 −dk/dt 与 DNS 的相对 L2 误差和 k(t) 的相对 L2 误差。
"""
import argparse
import csv
import os

import numpy as np


def read_monitor(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    t = np.array([float(r['t']) for r in rows])
    d = {k: np.array([float(r[k]) for r in rows]) for k in rows[0].keys()}
    # 去掉重复时刻（初始采样）
    _, idx = np.unique(t, return_index=True)
    return {k: v[idx] for k, v in d.items()}


def read_ref(path):
    data = np.array([[float(x) for x in l.split()] for l in open(path) if l.strip() and not l.startswith('#')])
    return {'t': data[:, 0], 'k': data[:, 1], 'eps_t': data[:, 2], 'eps': data[:, 3], 'enstrophy': data[:, 4]}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--ref', default=os.path.join(here, '..', 'cases', 'tgv3d', 'reference', 'TGV_Re1600.dat'))
    ap.add_argument('--out')
    ap.add_argument('--labels')
    ap.add_argument('--tmax', type=float, default=20.0)
    ap.add_argument('--check', type=float, help='−dk/dt 峰值与 DNS 的相对误差容差（%%）且峰值时刻在 7.5~10.5；失败退出码 1')
    a = ap.parse_args()
    ref = read_ref(a.ref)
    labels = a.labels.split(',') if a.labels else [os.path.relpath(f) for f in a.files]
    ipk = np.argmax(ref['eps_t'])
    print(f"DNS: peak -dk/dt {ref['eps_t'][ipk]:.5f} at t = {ref['t'][ipk]:.2f}")
    series = []
    for path, label in zip(a.files, labels):
        m = read_monitor(path)
        t, k = m['t'], m['k']
        dkdt = -np.gradient(k, t)
        sel = (t >= 0) & (t <= a.tmax)
        er = np.interp(t[sel], ref['t'], ref['eps_t'])
        kr = np.interp(t[sel], ref['t'], ref['k'])
        l2e = np.sqrt(np.sum((dkdt[sel] - er) ** 2) / np.sum(er ** 2))
        l2k = np.sqrt(np.sum((k[sel] - kr) ** 2) / np.sum(kr ** 2))
        ip = np.argmax(dkdt[sel])
        tot = m['epsilon'] + m['epsilon_sgs']
        ipr = np.argmax(tot[sel])
        print(f"{label}: peak -dk/dt {dkdt[sel][ip]:.5f} at t = {t[sel][ip]:.2f} "
              f"({100 * (dkdt[sel][ip] / ref['eps_t'][ipk] - 1):+.1f}% vs DNS); "
              f"resolved+sgs eps peak {tot[sel][ipr]:.5f} at t = {t[sel][ipr]:.2f}; "
              f"L2 err -dk/dt {100 * l2e:.1f}%, k {100 * l2k:.1f}%")
        series.append((label, t, k, dkdt, m['epsilon'], m['epsilon_sgs']))
        if a.check is not None:
            dev = 100 * abs(dkdt[sel][ip] / ref['eps_t'][ipk] - 1)
            ok = dev <= a.check and 7.5 <= t[sel][ip] <= 10.5
            print('PASS' if ok else f'FAIL: peak deviation {dev:.1f}% > {a.check}% or peak time outside 7.5-10.5')
            if not ok:
                raise SystemExit(1)
    if a.out:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(ref['t'], ref['k'], 'k-', label='DNS 512³ (Dairay 2017)')
        ax[1].plot(ref['t'], ref['eps_t'], 'k-', label='DNS −dk/dt')
        for label, t, k, dkdt, eps, esgs in series:
            ax[0].plot(t, k, '--', label=label)
            (ln,) = ax[1].plot(t, dkdt, '--', label=f'{label} −dk/dt')
            ax[1].plot(t, eps + esgs, ':', color=ln.get_color(), label=f'{label} 2(ν+ν_t)⟨S:S⟩')
        ax[0].set_xlabel('t'); ax[0].set_ylabel('k'); ax[0].legend(fontsize=8)
        ax[1].set_xlabel('t'); ax[1].set_ylabel('dissipation'); ax[1].legend(fontsize=7)
        for x in ax: x.set_xlim(0, a.tmax); x.grid(alpha=0.3)
        fig.tight_layout(); fig.savefig(a.out, dpi=130)
        print('figure:', a.out)


if __name__ == '__main__':
    main()
