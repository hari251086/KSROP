#!/usr/bin/env python3
"""
benchmark.py  -  Performance profiler for driver_KS.F

Runs the compiled propagator under several configurations and
reports wall-clock time, steps/second, and relative cost of
each force model component.

Usage: python benchmark.py [driver_KS.exe]
"""

import subprocess, os, shutil, time, sys

EXE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else 'driver_KS.exe')

CONST_FILE = 'const_new.dat'
INPUT_FILE = 'input.dat'
BAK_CONST  = '_bm_const.bak'
BAK_INPUT  = '_bm_input.bak'

# Initial state from production input.DAT
X0  = '0.0 -5888.9727 -3400.0'
XD0 = '9.5 0.0 0.0'
CAL = '2016 09 20 0 0 0'
DRAG_LINE = '50.0 0 7.2921150d-5 3.35281066d-3 1.0'   # IDRAG=0

CONST_TEMPLATE = (
    '3.986004415D5 6378.1363D0 1.495978707d08 '
    '1.32712440018d11 4.902801076d3\n'
    '{ngeo} {nsun} {nmoon}\n'
)

INPUT_TEMPLATE = (
    '{x0}\n{xd0}\n{nrev} {istep} 1d-15\n{cal}\n'
    '{nf1} {nf2} {nf3}\n{drag}\n'
)

N_REPEAT = 3   # runs per config (take median)

# ----------------------------------------------------------------
def write_inputs(ngeo=0, nsun=0, nmoon=0,
                 nrev=1, istep=360,
                 nf1=0, nf2=0, nf3=0):
    with open(CONST_FILE, 'w') as f:
        f.write(CONST_TEMPLATE.format(ngeo=ngeo, nsun=nsun, nmoon=nmoon))
    with open(INPUT_FILE, 'w') as f:
        f.write(INPUT_TEMPLATE.format(
            x0=X0, xd0=XD0, nrev=nrev, istep=istep,
            cal=CAL, nf1=nf1, nf2=nf2, nf3=nf3, drag=DRAG_LINE))

def run_once():
    t0 = time.perf_counter()
    r = subprocess.run([EXE], capture_output=True, timeout=300)
    t1 = time.perf_counter()
    if r.returncode != 0:
        print('  ERROR:', r.stderr.decode(errors='replace')[:200])
        return None
    return (t1 - t0) * 1000   # ms

def run_median(label, **kwargs):
    write_inputs(**kwargs)
    times = []
    for _ in range(N_REPEAT):
        t = run_once()
        if t is None:
            return None
        times.append(t)
    times.sort()
    med = times[N_REPEAT // 2]
    nsteps = kwargs.get('nrev', 1) * kwargs.get('istep', 360)
    rate   = nsteps / (med / 1000)
    print(f'  {label:<42}  {med:>8.1f} ms   {rate:>9.0f} steps/s')
    return med

# ----------------------------------------------------------------
def main():
    if not os.path.isfile(EXE):
        print(f'ERROR: {EXE} not found'); sys.exit(1)

    # backup
    for src, bak in [(CONST_FILE, BAK_CONST), (INPUT_FILE, BAK_INPUT)]:
        real = next((f for f in os.listdir('.') if f.lower()==src.lower()), None)
        if real:
            shutil.copy(real, bak)

    print('=' * 65)
    print(f' KSROP Performance Benchmark   (median of {N_REPEAT} runs)')
    print('=' * 65)
    print(f'  {"Configuration":<42}  {"Time":>10}   {"Throughput":>12}')
    print(f'  {"-"*42}  {"-"*10}   {"-"*12}')

    # ---- 1. Integrator baseline (two-body, all off) ----
    print('\n[1] Integrator cost  (two-body, no perturbations)')
    t1  = run_median('nrev=1,   istep=360   (360 steps)',
                     nrev=1,   istep=360,  ngeo=0)
    t10 = run_median('nrev=10,  istep=360   (3 600 steps)',
                     nrev=10,  istep=360,  ngeo=0)
    t50 = run_median('nrev=50,  istep=360   (18 000 steps)',
                     nrev=50,  istep=360,  ngeo=0)
    t100= run_median('nrev=100, istep=360   (36 000 steps)',
                     nrev=100, istep=360,  ngeo=0)

    # ---- 2. Step-size sensitivity ----
    print('\n[2] Steps/revolution  (nrev=10, two-body)')
    run_median('istep=100  (1 000 steps)',  nrev=10, istep=100,  ngeo=0)
    run_median('istep=360  (3 600 steps)',  nrev=10, istep=360,  ngeo=0)
    run_median('istep=1000 (10 000 steps)', nrev=10, istep=1000, ngeo=0)
    run_median('istep=3600 (36 000 steps)', nrev=10, istep=3600, ngeo=0)

    # ---- 3. Geopotential (J2 only, ngeo=2) ----
    print('\n[3] J2 geopotential  (ngeo=2, nrev=10, istep=360)')
    tgeo = run_median('J2 on  (ngeo=2, nf1=1)',
                      nrev=10, istep=360, ngeo=2, nf1=1)

    # ---- 4. Drag ----
    print('\n[4] Atmospheric drag  (nrev=10, istep=360, ngeo=0)')
    DRAG_ON = '50.0 1 7.2921150d-5 3.35281066d-3 1.0'
    # temporarily patch DRAG_LINE
    import builtins
    orig_drag = DRAG_LINE
    globals()['DRAG_LINE'] = DRAG_ON
    run_median('Drag on  (IDRAG=1)', nrev=10, istep=360, ngeo=0)
    globals()['DRAG_LINE'] = orig_drag

    # ---- 5. Summary ratios ----
    print('\n[5] Scaling analysis')
    if t1 and t100:
        ratio = t100 / t1
        print(f'  100x more steps -> {ratio:.1f}x more time  '
              f'({"linear" if abs(ratio-100)<10 else "non-linear"})')
    if t10 and tgeo:
        geo_overhead = (tgeo - t10) / t10 * 100
        print(f'  J2 overhead vs two-body (nrev=10): '
              f'{geo_overhead:+.1f}%')

    print('=' * 65)

    # restore
    for src, bak in [(CONST_FILE, BAK_CONST), (INPUT_FILE, BAK_INPUT)]:
        if os.path.isfile(bak):
            shutil.copy(bak, src)
            os.remove(bak)

if __name__ == '__main__':
    main()
