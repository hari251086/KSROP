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

CONST_FILE   = 'input/const_new.dat'
INPUT_FILE   = 'input/input.dat'
OPM_FILE     = 'input/input.opm'
BAK_CONST    = '_bm_const.bak'
BAK_INPUT    = '_bm_input.bak'
BAK_OPM      = '_bm_opm.bak'

DRAG_OFF = '50.0 0 7.2921150d-5 3.35281066d-3 1.0'
DRAG_ON  = '50.0 1 7.2921150d-5 3.35281066d-3 1.0'

HAS_EGM = os.path.isfile('input/EGM2008_to2190_TideFree')

# Fixed initial state used across all benchmark configurations
OPM_CONTENT = """\
CCSDS_OPM_VERS = 2.0
CREATION_DATE  = 2016-09-20T00:00:00.000
ORIGINATOR     = KSROP

META_START
OBJECT_NAME    = SATELLITE
CENTER_NAME    = EARTH
REF_FRAME      = EME2000
TIME_SYSTEM    = UTC
META_STOP

STATE_VECTOR
EPOCH          = 2016-09-20T00:00:00.000
X              =        0.000000 [km]
Y              =    -5888.972700 [km]
Z              =    -3400.000000 [km]
X_DOT          =        9.500000 [km/s]
Y_DOT          =        0.000000 [km/s]
Z_DOT          =        0.000000 [km/s]
"""

N_REPEAT = 3     # runs per config (median taken)
N_REPEAT_GEO = 2

# ----------------------------------------------------------------
def write_inputs(ngeo=0, nsun=0, nmoon=0,
                 nrev=1, istep=360,
                 nf1=0, nf2=0, nf3=0,
                 drag_line=DRAG_OFF):
    with open(CONST_FILE, 'w') as f:
        f.write(
            '3.986004415D5 6378.1363D0 1.495978707d08 '
            '1.32712440018d11 4.902801076d3\n'
            f'{ngeo} {nsun} {nmoon}\n'
        )
    with open(OPM_FILE, 'w') as f:
        f.write(OPM_CONTENT)
    with open(INPUT_FILE, 'w') as f:
        f.write(f'{nrev} {istep} 1d-15\n{nf1} {nf2} {nf3}\n{drag_line}\n')

def run_once():
    t0 = time.perf_counter()
    r = subprocess.run([EXE], capture_output=True, timeout=600)
    t1 = time.perf_counter()
    if r.returncode != 0:
        err = r.stderr.decode(errors='replace')
        print(f'\n    ERROR: {err[:120].strip()}')
        return None
    return (t1 - t0) * 1000   # ms

def run_median(label, n_rep=N_REPEAT, **kwargs):
    write_inputs(**kwargs)
    times = []
    for _ in range(n_rep):
        t = run_once()
        if t is None:
            return None
        times.append(t)
    times.sort()
    med = times[n_rep // 2]
    nsteps = kwargs.get('nrev', 1) * kwargs.get('istep', 360)
    rate   = nsteps / (med / 1000) if med > 0 else 0
    col_t  = f'{med:>9.1f} ms'
    col_r  = f'{rate:>10.0f} steps/s' if nsteps > 0 else ' ' * 18
    print(f'  {label:<48}  {col_t}  {col_r}')
    return med

def sep():
    print(f'  {"-"*48}  {"-"*9}    {"-"*16}')

# ----------------------------------------------------------------
def main():
    if not os.path.isfile(EXE):
        print(f'ERROR: {EXE} not found')
        print('  Compile: ifort driver_KS.F Subrouts.F Legendre.F -o driver_KS.exe')
        sys.exit(1)

    # backup originals
    for src, bak in [(CONST_FILE, BAK_CONST),
                     (INPUT_FILE, BAK_INPUT),
                     (OPM_FILE,   BAK_OPM)]:
        if os.path.isfile(src):
            shutil.copy(src, bak)

    print('=' * 82)
    print(f'  KSROP Performance Benchmark        (median of {N_REPEAT} runs)')
    print(f'  EGM2008 file present: {HAS_EGM}')
    print('=' * 82)
    print(f'  {"Configuration":<48}  {"Time":>9}    {"Throughput":>16}')
    sep()

    # ==============================================================
    # 1. Integrator baseline — two-body, no perturbations
    # ==============================================================
    print('\n[1]  Integrator baseline  (two-body, all forces off, ngeo=0)')
    t_1_360    = run_median('nrev=1,    istep=360   ->    360 steps',
                            nrev=1,   istep=360)
    t_10_360   = run_median('nrev=10,   istep=360   ->  3,600 steps',
                            nrev=10,  istep=360)
    t_50_360   = run_median('nrev=50,   istep=360   -> 18,000 steps',
                            nrev=50,  istep=360)
    t_100_360  = run_median('nrev=100,  istep=360   -> 36,000 steps',
                            nrev=100, istep=360)

    # ==============================================================
    # 2. Step-size sensitivity
    # ==============================================================
    print('\n[2]  Steps/revolution sensitivity  (nrev=10, two-body)')
    run_median('istep=100   ->  1,000 steps',  nrev=10, istep=100)
    run_median('istep=360   ->  3,600 steps',  nrev=10, istep=360)
    run_median('istep=1000  -> 10,000 steps',  nrev=10, istep=1000)
    run_median('istep=3600  -> 36,000 steps',  nrev=10, istep=3600)

    # ==============================================================
    # 3. Atmospheric drag
    # ==============================================================
    print('\n[3]  Atmospheric drag  (nrev=10, istep=360, ngeo=0)')
    t_nodrag = run_median('Drag OFF  (IDRAG=0)',
                          nrev=10, istep=360, drag_line=DRAG_OFF)
    t_drag   = run_median('Drag ON   (IDRAG=1)',
                          nrev=10, istep=360, drag_line=DRAG_ON)
    if t_nodrag and t_drag:
        ovhd = (t_drag - t_nodrag) / t_nodrag * 100
        print(f'  -> Drag overhead per run: {t_drag-t_nodrag:.1f} ms  '
              f'({ovhd:+.1f}% vs no-drag)')

    # ==============================================================
    # 4. Geopotential — EGM2008 required
    # ==============================================================
    if HAS_EGM:
        print('\n[4]  Geopotential  (EGM2008 — includes full file read)')

        # Startup-only: nrev=0 measures ONLY file I/O + initialisation
        print('    [4a] EGM2008 file-read overhead  (nrev=0, istep=1)')
        t_startup_0   = run_median('ngeo=0   (no file read)',
                                   nrev=0, istep=1, ngeo=0, nf1=0,
                                   n_rep=N_REPEAT_GEO)
        t_startup_2   = run_median('ngeo=2   (reads 2 401 333 lines)',
                                   nrev=0, istep=1, ngeo=2, nf1=1,
                                   n_rep=N_REPEAT_GEO)
        t_startup_50  = run_median('ngeo=50  (reads 2 401 333 lines)',
                                   nrev=0, istep=1, ngeo=50, nf1=1,
                                   n_rep=N_REPEAT_GEO)
        t_startup_100 = run_median('ngeo=100 (reads 2 401 333 lines)',
                                   nrev=0, istep=1, ngeo=100, nf1=1,
                                   n_rep=N_REPEAT_GEO)

        if t_startup_0 and t_startup_2 and t_startup_50 and t_startup_100:
            egm_times = [t_startup_2, t_startup_50, t_startup_100]
            spread = max(egm_times) - min(egm_times)
            print(f'  -> EGM8 file-read cost: '
                  f'{t_startup_2 - t_startup_0:.0f} ms '
                  f'(degree-independent, spread across degrees: {spread:.0f} ms)')

        # Integration cost at various degrees (nrev=10, istep=360)
        print('\n    [4b] Integration cost at different geopotential degrees')
        print('         (total time = file-read overhead + integration)')
        t_geo_2   = run_median('ngeo=2    J2 only',
                               nrev=10, istep=360, ngeo=2,   nf1=1,
                               n_rep=N_REPEAT_GEO)
        t_geo_10  = run_median('ngeo=10   J2–J10',
                               nrev=10, istep=360, ngeo=10,  nf1=1,
                               n_rep=N_REPEAT_GEO)
        t_geo_50  = run_median('ngeo=50   J2–J50  (aLegP limit)',
                               nrev=10, istep=360, ngeo=50,  nf1=1,
                               n_rep=N_REPEAT_GEO)

        if (t_startup_0 and t_startup_2 and
                t_geo_2 and t_geo_10 and t_geo_50 and t_10_360):
            egm_read  = t_startup_2  - t_startup_0
            intg_2    = t_geo_2  - t_startup_2  + t_startup_0
            intg_10   = t_geo_10 - t_startup_2  + t_startup_0
            intg_50   = t_geo_50 - t_startup_2  + t_startup_0
            intg_base = t_10_360
            print(f'\n  -> EGM8 file read:            {egm_read:>7.0f} ms '
                  f'(one-time per run, same for all degrees)')
            print(f'  -> Integration (two-body):    {intg_base:>7.0f} ms  '
                  f'for 3 600 steps')
            print(f'  -> Integration (ngeo=2,  J2): {intg_2:>7.0f} ms  '
                  f'for 3 600 steps  '
                  f'({(intg_2-intg_base)/intg_base*100:+.1f}%)')
            print(f'  -> Integration (ngeo=10):     {intg_10:>7.0f} ms  '
                  f'for 3 600 steps  '
                  f'({(intg_10-intg_base)/intg_base*100:+.1f}%)')
            print(f'  -> Integration (ngeo=50):     {intg_50:>7.0f} ms  '
                  f'for 3 600 steps  '
                  f'({(intg_50-intg_base)/intg_base*100:+.1f}%)')
    else:
        print('\n[4]  Geopotential — SKIPPED (input/EGM2008_to2190_TideFree not found)')

    # ==============================================================
    # 5. Scaling analysis
    # ==============================================================
    print('\n[5]  Scaling analysis  (two-body baseline)')
    if t_1_360 and t_50_360 and t_100_360:
        # Linear model fitted on 3 points: T = T_overhead + T_per_step * N
        n1, n2, n3 = 360, 18000, 36000
        t1, t2, t3 = t_1_360, t_50_360, t_100_360
        T_per_step  = (t3 - t1) / (n3 - n1)           # ms per step
        T_overhead  = t1 - T_per_step * n1             # ms fixed cost
        T_max_rate  = 1000.0 / T_per_step              # steps/s
        # linearity check: predicted t2 vs measured
        t2_pred = T_overhead + T_per_step * n2
        linearity = abs(t2_pred - t2) / t2 * 100

        print(f'  Fixed startup overhead:   {T_overhead:>6.1f} ms  '
              f'(file I/O + KS init)')
        print(f'  Variable cost per step:   {T_per_step*1000:>6.2f} µs/step')
        print(f'  Asymptotic throughput:    {T_max_rate:>6.0f} steps/s')
        print(f'  Linearity check (50 rev): predicted {t2_pred:.0f} ms, '
              f'measured {t2:.0f} ms  (error {linearity:.1f}%)')
        ratio = t_100_360 / t_1_360
        print(f'  100x more steps -> {ratio:.1f}x wall time  '
              f'(~linear once overhead excluded)')

    print('=' * 82)

    # restore
    for src, bak in [(CONST_FILE, BAK_CONST),
                     (INPUT_FILE, BAK_INPUT),
                     (OPM_FILE,   BAK_OPM)]:
        if os.path.isfile(bak):
            shutil.copy(bak, src)
            os.remove(bak)

if __name__ == '__main__':
    main()
