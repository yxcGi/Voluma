#!/usr/bin/env python3
"""Graetz 问题（层流圆管、壁温阶跃）局部 Nusselt 数 Nu_x(x*) 与级数解对比。

用法: compare_graetz.py <输出目录> --pe 100 [--check] [--tol 5]
读取 output/wallHeatFlux_wall.csv（x y z q Tw，q 为进入流体的运动学热流 q/(ρ c_p)），按轴向 z 分箱得 q̄(z)，
体温由能量平衡 T_b(z) = T_in + 4/(U_b D) ∫ q̄ dz 得到，Nu_x = q̄ D / (α (T_w − T_b))。
级数解（Shah & London 1978）：x* = x/(D Pe)，Nu_x = ½ Σ G_n e^{−2λ_n² x*} / Σ (G_n/λ_n²) e^{−2λ_n² x*}，
前 5 项 λ_n、G_n 取表值，其余用渐近式 λ_n = 4n + 8/3，G_n = 1.01276 λ_n^{−1/3}；x* → ∞ 时 Nu → λ_1²/2 = 3.657。
"""
import argparse
import csv
import math
import sys

LAM = [2.70436, 6.67903, 10.67338, 14.67108, 18.66987]
G = [0.74877, 0.54383, 0.46286, 0.41542, 0.38292]


def nu_series(xs, nterms=60):
    num = den = 0.0
    for n in range(1, nterms + 1):
        lam = LAM[n - 1] if n <= len(LAM) else 4 * n + 8.0 / 3
        g = G[n - 1] if n <= len(G) else 1.01276 * lam ** (-1.0 / 3)
        e = math.exp(-2 * lam * lam * xs)
        num += g * e
        den += g / (lam * lam) * e
    return 0.5 * num / den if den > 0 else 0.5 * LAM[0] ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--pe', type=float, required=True, help='Pe = U_b D/α')
    ap.add_argument('--D', type=float, default=1.0)
    ap.add_argument('--Ub', type=float, default=1.0)
    ap.add_argument('--alpha', type=float, default=None)
    ap.add_argument('--Tw', type=float, default=1.0)
    ap.add_argument('--Tin', type=float, default=0.0)
    ap.add_argument('--axis', type=int, default=2)
    ap.add_argument('--nbins', type=int, default=100)
    ap.add_argument('--xmin', type=float, default=0.004, help='比较用的最小 x*（入口奇点附近不比）')
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--tol', type=float, default=5.0)
    a = ap.parse_args()
    alpha = a.alpha if a.alpha else a.Ub * a.D / a.pe
    rows = [[float(v) for v in r] for r in list(csv.reader(open(a.out + '/wallHeatFlux_wall.csv')))[1:] if r]
    zs0 = [r[a.axis] for r in rows]
    z0, z1 = min(zs0), max(zs0)
    # 站位模式：有 line_st<z>_U.csv / _T.csv 时，体温取截面积分 T_b = ∫ u T r dr / ∫ u r dr，q̄ 取该 z 处壁面面的平均
    import glob
    import re
    stations = sorted((float(re.search(r'line_st([\d.]+)_U\.csv$', f).group(1)), f)
                      for f in glob.glob(a.out + '/line_st*_U.csv'))
    if stations:
        worst = 0.0
        nz = len(set(round(r[a.axis], 6) for r in rows))
        halfBin = 0.5 * (z1 - z0) / max(nz - 1, 1)  # 壁面面心沿轴向的站距之半：站位落在两排面心之间时取两排平均
        print(f'{"z":>6s} {"x*":>8s} {"T_b":>8s} {"Nu_x":>8s} {"series":>8s} {"err%":>7s}')
        for z, fU in stations:
            U = [[float(v) for v in r] for r in list(csv.reader(open(fU)))[1:] if r]
            T = [[float(v) for v in r] for r in list(csv.reader(open(fU.replace('_U.csv', '_T.csv'))))[1:] if r]
            r_ = [u[0] for u in U]
            uz = [u[3 + a.axis] for u in U]
            t = [x[3] for x in T]
            num = den = 0.0
            for i in range(len(r_) - 1):
                dr = r_[i + 1] - r_[i]
                num += 0.5 * (uz[i] * t[i] * r_[i] + uz[i + 1] * t[i + 1] * r_[i + 1]) * dr
                den += 0.5 * (uz[i] * r_[i] + uz[i + 1] * r_[i + 1]) * dr
            Tb = num / den
            near = [r for r in rows if abs(r[a.axis] - z) <= halfBin + 1e-9]
            if not near:
                continue
            q = sum(r[3] for r in near) / len(near)
            xs = z / (a.D * a.pe)
            nu = q * a.D / (alpha * (a.Tw - Tb))
            ref = nu_series(xs)
            err = 100 * (nu / ref - 1)
            if xs >= a.xmin:
                worst = max(worst, abs(err))
            print(f'{z:6.2f} {xs:8.4f} {Tb:8.4f} {nu:8.3f} {ref:8.3f} {err:+7.2f}')
        print(f'max |err| for x* >= {a.xmin}: {worst:.2f} %  (Nu_inf series 3.657)')
        if a.check:
            if worst > a.tol:
                print('FAIL')
                sys.exit(1)
            print('PASS')
        return
    # 站位模式：有 line_st<z>_U.csv / _T.csv 时，体温取截面积分 T_b = ∫ u T r dr / ∫ u r dr，q̄ 取该 z 处壁面面的平均
    import glob
    import re
    stations = sorted((float(re.search(r'line_st([\d.]+)_U\.csv$', f).group(1)), f)
                      for f in glob.glob(a.out + '/line_st*_U.csv'))
    if stations:
        worst = 0.0
        nz = len(set(round(r[a.axis], 6) for r in rows))
        halfBin = 0.5 * (z1 - z0) / max(nz - 1, 1)  # 壁面面心沿轴向的站距之半：站位落在两排面心之间时取两排平均
        print(f'{"z":>6s} {"x*":>8s} {"T_b":>8s} {"Nu_x":>8s} {"series":>8s} {"err%":>7s}')
        for z, fU in stations:
            U = [[float(v) for v in r] for r in list(csv.reader(open(fU)))[1:] if r]
            T = [[float(v) for v in r] for r in list(csv.reader(open(fU.replace('_U.csv', '_T.csv'))))[1:] if r]
            r_ = [u[0] for u in U]
            uz = [u[3 + a.axis] for u in U]
            t = [x[3] for x in T]
            num = den = 0.0
            for i in range(len(r_) - 1):
                dr = r_[i + 1] - r_[i]
                num += 0.5 * (uz[i] * t[i] * r_[i] + uz[i + 1] * t[i + 1] * r_[i + 1]) * dr
                den += 0.5 * (uz[i] * r_[i] + uz[i + 1] * r_[i + 1]) * dr
            Tb = num / den
            near = [r for r in rows if abs(r[a.axis] - z) <= halfBin + 1e-9]
            if not near:
                continue
            q = sum(r[3] for r in near) / len(near)
            xs = z / (a.D * a.pe)
            nu = q * a.D / (alpha * (a.Tw - Tb))
            ref = nu_series(xs)
            err = 100 * (nu / ref - 1)
            if xs >= a.xmin:
                worst = max(worst, abs(err))
            print(f'{z:6.2f} {xs:8.4f} {Tb:8.4f} {nu:8.3f} {ref:8.3f} {err:+7.2f}')
        print(f'max |err| for x* >= {a.xmin}: {worst:.2f} %  (Nu_inf series 3.657)')
        if a.check:
            if worst > a.tol:
                print('FAIL')
                sys.exit(1)
            print('PASS')
        return
    zs = [r[a.axis] for r in rows]
    z0, z1 = min(zs), max(zs)
    nb = a.nbins
    sq = [0.0] * nb
    cnt = [0] * nb
    for r in rows:
        b = min(int((r[a.axis] - z0) / (z1 - z0 + 1e-12) * nb), nb - 1)
        sq[b] += r[3]
        cnt[b] += 1
    dz = (z1 - z0) / nb
    Tb = a.Tin
    worst = 0.0
    print(f'{"x*":>9s} {"Nu_x":>8s} {"series":>8s} {"err%":>7s}')
    nshown = 0
    for b in range(nb):
        if cnt[b] == 0:
            continue
        q = sq[b] / cnt[b]
        zc = z0 + (b + 0.5) * dz
        # 体温：到该箱中心的能量平衡（半箱用当前 q）
        TbMid = Tb + 4 / (a.Ub * a.D) * q * dz / 2
        Tb += 4 / (a.Ub * a.D) * q * dz
        xs = zc / (a.D * a.pe)
        nu = q * a.D / (alpha * (a.Tw - TbMid))
        ref = nu_series(xs)
        err = 100 * (nu / ref - 1)
        if xs >= a.xmin:
            worst = max(worst, abs(err))
        if b % max(nb // 12, 1) == 0 or b == nb - 1:
            print(f'{xs:9.4f} {nu:8.3f} {ref:8.3f} {err:+7.2f}')
            nshown += 1
    print(f'max |err| for x* >= {a.xmin}: {worst:.2f} %  (Nu_inf series 3.657)')
    if a.check:
        if worst > a.tol:
            print('FAIL')
            sys.exit(1)
        print('PASS')


if __name__ == '__main__':
    main()
