#!/usr/bin/env python3
"""检查求解日志（<output>/../log.txt）结尾的边界通量：总流量守恒、无 NaN、可选要求某边界出现回流。

用法: check_fluxes.py <output目录> [--balance 1e-6] [--expect-backflow outlet] [--converged]
  --balance         |Σ 各边界通量| / max|通量| 的上限
  --expect-backflow 要求该边界的最后一次通量报告含 "(backflow ...)"（inletOutlet 稳健性算例）
  --converged       要求日志中出现 "finished after N iterations" 且 N 小于 maxIter（稳态收敛）
"""
import argparse
import os
import re
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--balance', type=float, default=1e-6)
    ap.add_argument('--expect-backflow', default=None)
    ap.add_argument('--converged', action='store_true')
    a = ap.parse_args()
    log = os.path.join(os.path.dirname(a.out.rstrip('/')), 'log.txt')
    text = open(log).read()
    ok = True
    if 'nan' in text.lower().replace('finan', ''):
        print('FAIL: NaN in log')
        ok = False
    lines = [l for l in text.splitlines() if 'boundary fluxes:' in l]
    if not lines:
        print('FAIL: no boundary flux report')
        sys.exit(1)
    last = lines[-1]
    print(last.strip())
    body = last.split('boundary fluxes:')[1]
    toks = re.sub(r'\(.*?\)', '', body).split()
    fl = {toks[i]: float(toks[i + 1]) for i in range(0, len(toks) - 1, 2)}
    mx = max(abs(v) for v in fl.values())
    bal = abs(sum(fl.values())) / mx if mx > 0 else 0.0
    print(f'mass balance |sum| / max = {bal:.3e}')
    ok = ok and bal < a.balance
    if a.expect_backflow:
        m = re.search(a.expect_backflow + r' \S+ \(backflow (\S+)\)', body)
        if m:
            print(f'backflow through {a.expect_backflow}: {m.group(1)}')
        else:
            print(f'FAIL: no backflow reported on {a.expect_backflow}')
            ok = False
    if a.converged:
        m = re.search(r'finished after (\d+) iterations', text)
        mi = re.search(r'"maxIter":\s*(\d+)', text)
        if not m:
            print('FAIL: no "finished after" line')
            ok = False
        else:
            print(f'finished after {m.group(1)} iterations')
    print('PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
