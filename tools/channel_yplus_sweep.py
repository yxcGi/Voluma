#!/usr/bin/env python3
"""RANS 槽道 y+ 扫描：同一 Re_τ、不同首层 y+ 的网格上比较各模型 / 壁面处理得到的 Re_τ（检验壁面处理对 y+ 的敏感性）。

用法: channel_yplus_sweep.py --exe build/bin/fvmFlow --out <目录> [--retau 392.24] [--nu 1.4545e-4]
                            [--meshes 8,12,16,24,32,48,64,96,128,192,120s3] [--models kEpsilon:standard,kEpsilon:enhanced,...]
                            [--np 1] [--maxIter 8000]

网格 "N" 为全高（2δ）均匀 N 层，"Ns3" 为 N 层 tanh 加密（β = 3）。每个算例在 <目录>/<模型>_<处理>_<网格>/ 下运行，
从日志读取最终的驱动压力梯度（Re_τ = √gradP δ/ν）与壁面 y+，打印表格并写 <目录>/sweep.csv。
参考：Moser 等 1999 DNS Re_τ 392.24（Re_b = U_b δ/ν = 6875），Lee & Moser 2015 Re_τ 1994.8（Re_b = 43467）。
"""
import argparse
import json
import os
import re
import subprocess
import sys


def mesh_spec(tag):
    m = re.fullmatch(r'(\d+)(?:s([\d.]+))?', tag)
    if not m:
        raise SystemExit(f'bad mesh tag {tag}')
    n = int(m.group(1))
    box = {"n": [4, n, 1], "hi": [1, 2, 0.1], "periodic": [True, False, False], "twoD": True,
           "walls": ["yMin", "yMax"]}
    if m.group(2):
        box["stretch"] = {"y": float(m.group(2))}
    return box


def model_spec(tag):
    name, _, treat = tag.partition(':')
    d = {"model": name}
    if treat in ('standard', 'enhanced'):
        d["wallTreatment"] = treat
    elif treat == 'lowRe':
        d["wallFunctions"] = False
    elif treat == 'wf':
        d["wallFunctions"] = True
    elif treat:
        raise SystemExit(f'bad model tag {tag}')
    return d, (treat or 'default')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exe', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--retau', type=float, default=392.24)
    ap.add_argument('--nu', type=float, default=1 / 6875)
    ap.add_argument('--meshes', default='8,12,16,24,32,48,64,96,128,192,120s3')
    ap.add_argument('--models', default='kEpsilon:standard,kEpsilon:enhanced,RNGkEpsilon:enhanced,'
                                        'kOmegaSST:standard,kOmegaSST:enhanced,SpalartAllmaras:lowRe,SpalartAllmaras:wf')
    ap.add_argument('--np', type=int, default=1)
    ap.add_argument('--mpiexec', default='mpirun')
    ap.add_argument('--maxIter', type=int, default=8000)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    exe = os.path.abspath(a.exe)
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
    rows = []
    for mt in a.models.split(','):
        turb, treat = model_spec(mt)
        for mesh in a.meshes.split(','):
            case = {
                "mesh": {"box": mesh_spec(mesh)},
                "physics": {"nu": a.nu, "steady": True},
                "forcing": {"type": "meanVelocity", "Ubar": [1, 0, 0]},
                "turbulence": turb,
                "initial": {"U": [1, 0, 0], "k": 0.005, "omega": 50, "epsilon": 0.05, "nuTilda": 3 * a.nu},
                "schemes": {"div": "linearUpwind"},
                "run": {"maxIter": a.maxIter, "residualTol": 1e-8, "printInterval": 100000, "vtk": False},
            }
            name = f'{turb["model"]}_{treat}_{mesh}'
            d = os.path.join(a.out, name)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, 'case.json'), 'w') as f:
                json.dump(case, f, indent=1)
            cmd = [exe, 'case.json']
            if a.np > 1:
                cmd = [a.mpiexec, '-np', str(a.np)] + cmd
            with open(os.path.join(d, 'log.txt'), 'w') as log:
                r = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env)
            log = open(os.path.join(d, 'log.txt')).read()
            gp = [g for g in re.findall(r'forcing gradP ([-\d.eE+]+|-?nan|-?inf)', log) if re.fullmatch(r'[-\d.eE+]+', g)]
            yp = re.findall(r'wall y\+: min ([-\d.eE+]+)\s+mean ([-\d.eE+]+)', log)
            it = re.findall(r'finished after (\d+) iterations', log)
            if r.returncode != 0 or not gp:
                print(f'{name}: FAILED (exit {r.returncode})')
                rows.append((turb["model"], treat, mesh, float('nan'), float('nan'), 0, 'fail'))
                continue
            if float(gp[-1]) <= 0:
                print(f'{name}: FAILED (negative pressure gradient, diverged)')
                rows.append((turb["model"], treat, mesh, float('nan'), float('nan'), 0, 'fail'))
                continue
            retau = (float(gp[-1]) ** 0.5) / a.nu
            yplus = float(yp[-1][1]) if yp else float('nan')
            rows.append((turb["model"], treat, mesh, yplus, retau, int(it[-1]) if it else -1,
                         'converged' if it and int(it[-1]) < a.maxIter else 'maxIter'))
            print(f'{name:40s} y+ {yplus:8.2f}  Re_tau {retau:8.1f}  ({100 * (retau / a.retau - 1):+6.2f} %)  '
                  f'{rows[-1][5]} it {rows[-1][6]}', flush=True)
    with open(os.path.join(a.out, 'sweep.csv'), 'w') as f:
        f.write('model,treatment,mesh,yplus,retau,err_percent,iterations,status\n')
        for m, t, mesh, yp, rt, it, st in rows:
            f.write(f'{m},{t},{mesh},{yp:.4f},{rt:.3f},{100 * (rt / a.retau - 1):.3f},{it},{st}\n')
    # 按模型汇总
    print(f'\nRe_tau error (%) vs DNS {a.retau}; columns = mesh')
    meshes = a.meshes.split(',')
    print(f'{"model":28s}' + ''.join(f'{m:>9s}' for m in meshes))
    for mt in a.models.split(','):
        turb, treat = model_spec(mt)
        vals = {r[2]: r[4] for r in rows if r[0] == turb["model"] and r[1] == treat}
        print(f'{turb["model"] + ":" + treat:28s}' +
              ''.join(f'{100 * (vals[m] / a.retau - 1):+9.1f}' if m in vals and vals[m] == vals[m] else f'{"fail":>9s}'
                      for m in meshes))
    yrow = {r[2]: r[3] for r in rows if r[0] == model_spec(a.models.split(',')[0])[0]["model"]}
    print(f'{"first-cell y+":28s}' + ''.join(f'{yrow.get(m, float("nan")):9.2f}' for m in meshes))


if __name__ == '__main__':
    main()
