#!/usr/bin/env python3
"""湍流槽道传热（均匀体积热源、等温壁、流向周期）结果对比。

用法: compare_channel_heat.py <输出目录> --log log.txt [--nu 1.4545e-4 --Pr 0.71 --Q 1 --delta 1] [--check] [--tol 10]
读取 output/line_profile_U.csv、line_profile_T.csv（y 从壁面到中心），log 中的 forcing gradP（u_τ² / δ）。
  T+ = (T_w − T) u_τ / q_w，q_w = Q δ（稳态能量平衡）；与 Kader 1981 律 T+ = Pr y+ e^{−Γ} + [2.12 ln(1+y+) + β] e^{−1/Γ} 比较
  （y+ ≥ 30 内 T+ 的平均相对偏差）；Nu = q_w D_h / (α (T_w − T_b))，D_h = 4δ，T_b = ∫U T dy / ∫U dy，
  与 Gnielinski 关联式（Re_Dh = 4 Re_b）比较。
"""
import argparse
import csv
import math
import re
import sys


def kader(yp, Pr):
    G = 0.01 * (Pr * yp) ** 4 / (1 + 5 * Pr ** 3 * yp)
    beta = (3.85 * Pr ** (1 / 3) - 1.3) ** 2 + 2.12 * math.log(Pr)
    return Pr * yp * math.exp(-G) + (2.12 * math.log(1 + yp) + beta) * math.exp(-1 / max(G, 1e-300))


def gnielinski(Re, Pr):
    f = (0.79 * math.log(Re) - 1.64) ** -2
    return (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2 / 3) - 1))


def read_line(path):
    return [[float(x) for x in r] for r in list(csv.reader(open(path)))[1:] if r]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--log', required=True)
    ap.add_argument('--nu', type=float, default=1 / 6875)
    ap.add_argument('--Pr', type=float, default=0.71)
    ap.add_argument('--Q', type=float, default=1.0)
    ap.add_argument('--delta', type=float, default=1.0)
    ap.add_argument('--Ub', type=float, default=1.0)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--tol', type=float, default=10.0, help='T+ 对 Kader 律平均偏差容差（%%）')
    ap.add_argument('--tolNu', type=float, default=25.0, help='Nu 对 Gnielinski 关联式偏差容差（%%，关联式本身对槽道约 ±15 %%）')
    a = ap.parse_args()
    text = open(a.log).read()
    gp = float(re.findall(r'forcing gradP ([-\d.eE+]+)', text)[-1])
    uTau = math.sqrt(gp * a.delta)
    alpha = a.nu / a.Pr
    qw = a.Q * a.delta
    U = read_line(a.out + '/line_profile_U.csv')
    T = read_line(a.out + '/line_profile_T.csv')
    y = [r[1] for r in U]
    u = [r[3] for r in U]
    t = [r[3] for r in T]
    # 体温（梯形积分，线采样等距）
    num = sum(0.5 * (u[i] * t[i] + u[i + 1] * t[i + 1]) * (y[i + 1] - y[i]) for i in range(len(y) - 1))
    den = sum(0.5 * (u[i] + u[i + 1]) * (y[i + 1] - y[i]) for i in range(len(y) - 1))
    Tb = num / den
    Tw = 0.0
    Nu = qw * 4 * a.delta / (alpha * (Tb - Tw))
    ReDh = 4 * a.Ub * a.delta / a.nu
    NuG = gnielinski(ReDh, a.Pr)
    errs = []
    print(f'u_tau {uTau:.5f}  Re_tau {uTau * a.delta / a.nu:.1f}  T_b {Tb:.4f}')
    print(f'{"y+":>8s} {"T+":>8s} {"Kader":>8s} {"err%":>7s}')
    for i in range(0, len(y), max(len(y) // 15, 1)):
        yp = y[i] * uTau / a.nu
        if yp <= 0:
            continue
        tp = (t[i] - Tw) * uTau / qw
        ref = kader(yp, a.Pr)
        e = 100 * (tp / ref - 1)
        if yp >= 30:
            errs.append(e)
        print(f'{yp:8.2f} {tp:8.3f} {ref:8.3f} {e:+7.2f}')
    meanErr = sum(errs) / len(errs) if errs else float('nan')
    errNu = 100 * (Nu / NuG - 1)
    print(f'T+ vs Kader (y+ >= 30): mean {meanErr:+.2f} %')
    print(f'Nu_Dh {Nu:.2f}  Gnielinski {NuG:.2f} ({errNu:+.2f} %)')
    if a.check:
        if abs(meanErr) > a.tol or abs(errNu) > a.tolNu:
            print('FAIL')
            sys.exit(1)
        print('PASS')


if __name__ == '__main__':
    main()
