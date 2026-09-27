#!/usr/bin/env python3
"""槽道统计与 Lee & Moser DNS 对比。

用法: compare_channel.py channel_stats.csv [more.csv ...] --ref referenceData_1000.csv [--out fig.png] [--nu 5e-5]

channel_stats.csv 为 fvmFlow channelStatistics（或 OpenLB channel3d）输出，分号分隔：
  y+;u_tau;uAv+;uu++;uv++;uw++;vv++;vw++;ww++;pRMS;...  （u 流向、v 展向、w 法向）
LM_reference.csv 为 OpenLB reference_datas/referenceData_<Re>.csv（逗号分隔，同样 u 流向、v 展向、w 法向）。
打印 u_τ、U+ 在 y+ = 30..0.8 Re_τ 上的误差，以及解析 + 亚格子雷诺应力；可选输出对比图。
"""
import argparse
import csv

import numpy as np


def read_stats(path):
    with open(path) as f:
        rows = list(csv.reader(f, delimiter=';'))
    head = [h.strip() for h in rows[0]]
    data = np.array([[float(x) for x in r] for r in rows[1:] if r], dtype=float)
    return {h: data[:, i] for i, h in enumerate(head)}


def read_lm(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    head = [h.strip() for h in rows[0]]
    data = np.array([[float(x) for x in r] for r in rows[1:] if r], dtype=float)
    return {h: data[:, i] for i, h in enumerate(head)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stats', nargs='+')
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out')
    ap.add_argument('--retau', type=float, default=1000.512)
    ap.add_argument('--nu', type=float, default=5e-5)
    ap.add_argument('--H', type=float, default=1.0)
    a = ap.parse_args()
    lm = read_lm(a.ref)
    yl, ul = lm['y+'], lm['u+']
    series = []
    for path in a.stats:
        s = read_stats(path)
        yp, up = s['y+'], s['uAv+']
        ut = s['u_tau'][0]
        sel = (yp >= 30) & (yp <= 0.8 * a.retau)
        ref = np.interp(yp, yl, ul)
        err = (up - ref) / ref
        tot = {k: s[k + '++'] + s.get(k + '_sgs++', 0 * yp) for k in ('uu', 'vv', 'ww', 'uw')}
        print(f'{path}: u_tau {ut:.5f}  Re_tau {ut * a.H / a.nu:.1f}  (DNS {a.retau})')
        print(f'  U+ error (30 <= y+ <= {0.8 * a.retau:.0f}): mean {100 * np.mean(err[sel]):+.2f}%  '
              f'max |.| {100 * np.max(np.abs(err[sel])):.2f}%')
        print('  y+      U+     U+_DNS   uu+    uu+_DNS  -uw+   -uw+_DNS')
        lmuu, lmuv = lm["u'u'"], lm["u'w'"]
        for i in range(len(yp)):
            if i % max(1, len(yp) // 10):
                continue
            print(f'  {yp[i]:6.1f} {up[i]:7.3f} {ref[i]:7.3f} {tot["uu"][i]:6.3f} {np.interp(yp[i], yl, lmuu):6.3f}'
                  f'  {-tot["uw"][i]:6.3f} {-np.interp(yp[i], yl, lmuv):6.3f}')
        series.append((path, s, tot))
    if a.out:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        ax[0].semilogx(yl[1:], ul[1:], 'k-', label='DNS Lee & Moser')
        ax[1].plot(yl, lm["u'u'"], 'k-', label="u'u' DNS")
        ax[1].plot(yl, lm["w'w'"], 'k--', label="w'w' (normal) DNS")
        ax[1].plot(yl, lm["v'v'"], 'k:', label="v'v' (span) DNS")
        ax[1].plot(yl, -lm["u'w'"], 'k-.', label="-u'w' DNS")
        for path, s, tot in series:
            ax[0].semilogx(s['y+'], s['uAv+'], 'o-', ms=3, label=path)
            ax[1].plot(s['y+'], tot['uu'], 'o-', ms=3, label='uu ' + path)
            ax[1].plot(s['y+'], tot['ww'], 's--', ms=3, label='normal')
            ax[1].plot(s['y+'], tot['vv'], '^:', ms=3, label='span')
            ax[1].plot(s['y+'], -tot['uw'], 'd-.', ms=3, label='-uw')
        ax[0].set_xlabel('y+'), ax[0].set_ylabel('U+'), ax[0].legend(fontsize=8)
        ax[1].set_xlabel('y+'), ax[1].set_xlim(0, a.retau), ax[1].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(a.out, dpi=130)


if __name__ == '__main__':
    main()
