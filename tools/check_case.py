#!/usr/bin/env python3
"""回归测试驱动（ctest 调用）：把算例复制到工作目录，运行 fvmFlow，再用对应的对比脚本按容差判定。

用法: check_case.py --exe <fvmFlow> --case cases/pipe/case.json --workdir build/cases/pipe [--np N --mpiexec mpirun]
                    -- <对比脚本> <对比脚本参数...>（其中以 OUTPUT 开头的参数会把 OUTPUT 替换为输出目录）
"""
import os
import shutil
import subprocess
import sys


def main():
    args = sys.argv[1:]
    sep = args.index('--')
    opts, cmp_cmd = args[:sep], args[sep + 1:]
    o = {}
    i = 0
    while i < len(opts):
        o[opts[i].lstrip('-')] = opts[i + 1]
        i += 2
    case = os.path.abspath(o['case'])
    work = os.path.abspath(o['workdir'])
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    src_dir = os.path.dirname(case)
    for f in os.listdir(src_dir):
        p = os.path.join(src_dir, f)
        if os.path.isfile(p) and os.path.getsize(p) < 50_000_000:
            shutil.copy(p, work)
    dst = os.path.join(work, os.path.basename(case))
    run = [os.path.abspath(o['exe']), dst]
    if 'np' in o and int(o['np']) > 1:
        run = [o.get('mpiexec', 'mpirun'), '-np', o['np']] + run
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
    with open(os.path.join(work, 'log.txt'), 'w') as log:
        r = subprocess.run(run, cwd=work, stdout=log, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        print(open(os.path.join(work, 'log.txt')).read()[-3000:])
        print('FAIL: solver exit code', r.returncode)
        sys.exit(1)
    out = os.path.join(work, 'output')
    cmd = [sys.executable] + [c.replace('OUTPUT', out, 1) if c.startswith('OUTPUT') else c for c in cmp_cmd]
    sys.exit(subprocess.run(cmd).returncode)


if __name__ == '__main__':
    main()
