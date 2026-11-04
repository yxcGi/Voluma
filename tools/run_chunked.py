#!/usr/bin/env python3
"""分段续算驱动：把一个瞬态算例按时间段反复运行（每段结束写续算文件，下一段从续算文件起算），
中断后重新执行同一命令即可从最后完成的段继续。适合本地长算例（周期山 WMLES、槽道壁模型对比等）。

    tools/run_chunked.py cases/periodicHill/case.json --workdir runs/hill --np 8 --chunk 45 --end 360
    tools/run_chunked.py cases/channelWMLES/case.json --workdir runs/spalding --np 4 --chunk 100 \
        --set turbulence.wallModel.type=reichardt

- 算例文件复制到 workdir（网格 / 表文件等相对路径改为绝对路径），输出在 workdir/output，续算在 workdir/restart，
  日志追加到 workdir/log.txt，每段一行进度写 workdir/status.txt，全部完成写 workdir/DONE。
- --set 用点号路径覆盖任意键（数值 / true / false / 字符串自动识别），便于从同一算例派生变体。
- 统计量（channelStatistics、spanStatistics、fieldAverage）随续算累积，分段结果与一次跑完一致。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time


def load_json(path):
    with open(path) as f:
        text = '\n'.join(l for l in f.read().splitlines() if not l.lstrip().startswith('//'))
    return json.loads(text)


def absolutize(d, base):
    """把字典里名为 file / mesh.file 等指向存在文件的相对路径改成绝对路径"""
    if isinstance(d, dict):
        for k, v in d.items():
            if isinstance(v, str) and k in ('file',) and not os.path.isabs(v):
                cand = os.path.join(base, v)
                if os.path.exists(cand):
                    d[k] = os.path.abspath(cand)
            else:
                absolutize(v, base)
    elif isinstance(d, list):
        for v in d:
            absolutize(v, base)


def parse_value(s):
    for conv in (int, float):
        try:
            return conv(s)
        except ValueError:
            pass
    if s in ('true', 'false'):
        return s == 'true'
    if s.startswith('[') or s.startswith('{'):
        return json.loads(s)
    return s


def set_path(d, path, value):
    keys = path.split('.')
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


def restart_time(rdir):
    st = os.path.join(rdir, 'state.txt')
    if not os.path.exists(st):
        return None
    for line in open(st):
        parts = line.split()
        if len(parts) == 2 and parts[0] == 'time':
            return float(parts[1])
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case')
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--np', type=int, default=4)
    ap.add_argument('--chunk', type=float, default=100.0, help='每段时长')
    ap.add_argument('--end', type=float, help='结束时间（默认取算例的 run.endTime）')
    ap.add_argument('--exe', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'build', 'bin', 'fvmFlow'))
    ap.add_argument('--mpirun', default='mpirun')
    ap.add_argument('--set', action='append', default=[], metavar='KEY=VALUE', help='覆盖键，如 turbulence.wallModel.type=ode')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    case = load_json(a.case)
    absolutize(case, os.path.dirname(os.path.abspath(a.case)))
    for kv in a.set:
        k, v = kv.split('=', 1)
        set_path(case, k, parse_value(v))
    run = case.setdefault('run', {})
    if run.get('steady') or case.get('physics', {}).get('steady', False):
        sys.exit('run_chunked: only transient cases can be chunked (steady runs converge in one go)')
    end = a.end if a.end is not None else run.get('endTime')
    if end is None:
        sys.exit('run_chunked: need --end or run.endTime')
    wd = os.path.abspath(a.workdir)
    os.makedirs(wd, exist_ok=True)
    rdir = os.path.join(wd, 'restart')
    run['output'] = 'output'
    run['restartWrite'] = 'restart'
    exe = os.path.abspath(a.exe)
    status = os.path.join(wd, 'status.txt')
    done = os.path.join(wd, 'DONE')
    if os.path.exists(done):
        os.remove(done)

    t_start = restart_time(rdir)
    t = t_start if t_start is not None else 0.0
    if t_start is not None:
        print(f'resuming from restart at t = {t_start}')
    while t < end * (1 - 1e-12):
        t_next = min(end, (int(round(t / a.chunk + 1e-9)) + 1) * a.chunk)
        if t_next <= t * (1 + 1e-12):
            t_next = min(end, t + a.chunk)
        run['endTime'] = t_next
        if restart_time(rdir) is not None:
            run['restartRead'] = 'restart'
        else:
            run.pop('restartRead', None)
        cfile = os.path.join(wd, 'case.json')
        with open(cfile, 'w') as f:
            json.dump(case, f, indent=2, ensure_ascii=False)
        cmd = [a.mpirun, '-np', str(a.np), exe, 'case.json'] if a.np > 1 else [exe, 'case.json']
        print(f'[{time.strftime("%H:%M:%S")}] chunk to t = {t_next}: {" ".join(cmd)}', flush=True)
        if a.dry_run:
            t = t_next
            continue
        t0 = time.time()
        with open(os.path.join(wd, 'log.txt'), 'a') as log:
            rc = subprocess.call(cmd, cwd=wd, stdout=log, stderr=subprocess.STDOUT)
        with open(status, 'a') as st:
            st.write(f'{time.strftime("%Y-%m-%d %H:%M:%S")} to t={t_next} rc={rc} {time.time() - t0:.0f}s\n')
        if rc != 0:
            sys.exit(f'run_chunked: chunk to t = {t_next} failed with exit code {rc} (see {wd}/log.txt)')
        t_done = restart_time(rdir)
        if t_done is None or t_done < t_next * (1 - 1e-9):
            sys.exit(f'run_chunked: restart state not advanced to {t_next} (got {t_done})')
        t = t_done
    open(done, 'w').write(time.strftime('%Y-%m-%d %H:%M:%S\n'))
    print('all chunks done ->', wd)


if __name__ == '__main__':
    main()
