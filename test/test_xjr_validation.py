#!/usr/bin/env python3
"""
test_xjr_validation.py - Zonal geopotential cross-validation against
Xavier James Raj's PhD thesis (KSROP issue #29's "second reference",
E:\\Research\\References\\0\\29. XJR_Thesis.pdf -- not shipped in this
repo; only the already-published J_n / GEM-T2 C(n,0) numbers its
tables report are reproduced below as fixtures).

Two independent chapters of that thesis, both built on the same KS
uniformly-regular canonical formalism KSROP itself uses, published
their own RK4-Gill numerical results for the pure zonal (m=0) Earth
geopotential:

  Chapter 2 "Long Term Orbit Predictions with Earth's Oblateness":
    Table 2.1 (J2..J36), Table 2.2 (4 test cases' initial conditions),
    Tables 2.3/2.4 (semi-major axis / eccentricity after 100 revs).
  Chapter 3 "Long Term Orbit Predictions with Earth's Flattening":
    Table 3.1A (GEM-T2 normalized C(n,0)), Table 3.2 (3 cases'
    initial conditions), Table 3.4 (Delta a/e/i after 22h, J6/J10/
    J19 zonal-only columns).

This test drives KSROP's own driver_KS.exe with those initial
conditions and coefficients (via the fixture files below, temporarily
substituted for input/EGM2008_to2190_TideFree -- large, untracked,
frequently absent outside a local dev checkout) and compares the
resulting osculating elements to the thesis's own published values.
Case B of Chapter 3 is intentionally excluded: its initial state
vector in the source PDF has an OCR-corrupted velocity component and
could not be reliably transcribed (case A and C both reproduce the
thesis's own stated a/e/i to <0.01%, confirming the transcription
method itself, so this is a data-quality exclusion, not a code issue).

Tolerances are set several times wider than what was actually observed
in the original investigation (session of 2026-08-10) -- loose enough
to tolerate cross-platform/cross-compiler RK4-Gill roundoff (this runs
under gfortran in CI, ifx in local dev) without being so loose that a
genuine force-law regression (wrong sign, missing factor, wrong
normalization -- the actual bug classes issue #29/#30 hit historically)
would slip through.

Usage
-----
  python test_xjr_validation.py [executable]
"""

import subprocess, os, glob, math, shutil, sys
from datetime import datetime, timedelta

MU_CH2 = 398600.8      # km3/s2  -- thesis's own mu (Ch.2 text, p.46-47)
RE_CH2 = 6378.135      # km      -- thesis's own Re
FIXTURE_CH2 = 'test/fixture_xjr_ch2_zonal.dat'
FIXTURE_CH3 = 'test/fixture_xjr_ch3_zonal.dat'
GEO_TARGET  = 'input/EGM2008_to2190_TideFree'

# -------------------------------------------------------------------
# Chapter 2 (Oblateness): Table 2.2 initial conditions, Tables 2.3/2.4
# reference semi-major axis (km) / eccentricity after 100 revs at the
# thesis's own finest (360 steps/rev, converged for all 4 cases)
# resolution.
# -------------------------------------------------------------------
CH2_CASES = {
    'A': ((0.0, -5888.9727, -3400.0), (7.6, 0.0, 0.0)),
    'B': ((0.0, -5888.9727, -3400.0), (7.8, 0.0, 0.0)),
    'C': ((0.0, -5888.9727, -3400.0), (8.3, 0.0, 0.0)),
    'D': ((0.0, -5888.9727, -3400.0), (10.691338, 0.0, 0.0)),
}
CH2_REF = {
    'A': {2: (6705.9513, 0.0141452), 36: (6705.9665, 0.0137346)},
    'B': {2: (7071.7492, 0.0383313), 36: (7071.7653, 0.0386274)},
    'C': {2: (8245.8723, 0.1741792), 36: (8245.8644, 0.1743017)},
    'D': {2: (136127.431, 0.9500468), 36: (136159.326, 0.9500587)},
}
CH2_TRUNCS = [2, 36]
CH2_NREV = 100
CH2_ISTEP = 180  # converged for all 4 cases per the thesis's own Table 2.3
CH2_TOL_A_REL = 2.0e-3   # observed worst case 4.7e-4 -> ~4x margin
CH2_TOL_E_REL = 3.0e-2   # observed worst case 1.02e-2 -> ~3x margin

# -------------------------------------------------------------------
# Chapter 3 (Flattening): Table 3.2 cases A & C, Table 3.4 zonal-only
# (J_n,0) Delta(a,e,i) after 22h. Case B excluded (see module docstring).
# -------------------------------------------------------------------
CH3_CASES = {
    'A': ((1.0, -573.0, -6553.0), (7.9, 0.0, 0.0), 92.62),
    'C': ((1.0, -5888.972, -3400.0), (8.9, 0.0, 0.0), 178.01),
}
CH3_REF = {
    # case: {trunc_deg: (da_km, de, di_deg)}
    'A': {6: (19.3817, 23.279e-4, 67.82e-4),
          10: (19.3832, 23.283e-4, 67.82e-4),
          19: (19.3826, 23.281e-4, 67.82e-4)},
    'C': {6: (-3.3511, -6.4081e-4, 95.88e-4),
          10: (-3.3490, -6.4073e-4, 95.95e-4),
          19: (-3.3496, -6.4078e-4, 95.94e-4)},
}
CH3_TRUNCS = [6, 10, 19]
CH3_ISTEP = 240
CH3_DUR_H = 22.0
CH3_TOL_DA_KM = 3.0
CH3_TOL_DE = 2.0e-4
CH3_TOL_DI_DEG = 1.0e-2


class Counter:
    def __init__(self):
        self.npass = 0
        self.nfail = 0
    @property
    def total(self):
        return self.npass + self.nfail

def check(ctr, name, ok, detail=''):
    tag = 'PASS' if ok else 'FAIL'
    msg = f'      [{tag}]  {name}'
    if not ok and detail:
        msg += f'  ({detail})'
    print(msg)
    if ok:
        ctr.npass += 1
    else:
        ctr.nfail += 1
    return ok


def write_const(mu, re, ngeo_deg):
    with open('input/const_new.dat', 'w') as f:
        f.write(f'{mu}D0 {re}D0 1.495978707d08 '
                '1.32712440018d11 4.902801076d3\n')
        f.write(f'{ngeo_deg} 0 0\n')
        f.write('0.0d0\n')


def write_input(nrev, istep):
    with open('input/input.dat', 'w') as f:
        f.write(f'{nrev} {istep} 1d-15\n')
        f.write('1 0 0\n')
        f.write('50.0 0 7.2921150d-5 3.35281066d-3 1.0\n')
        f.write('1.2 0.01 0 0\n')


def write_opm(x, xd, epoch='1998-07-13T00:00:00.000'):
    with open('input/input.opm', 'w') as f:
        f.write('CCSDS_OPM_VERS = 2.0\n')
        f.write(f'CREATION_DATE  = {epoch}\n')
        f.write('ORIGINATOR     = KSROP\n\n')
        f.write('META_START\n')
        f.write('OBJECT_NAME    = XJR_VALIDATION\n')
        f.write('CENTER_NAME    = EARTH\n')
        f.write('REF_FRAME      = EME2000\n')
        f.write('TIME_SYSTEM    = UTC\n')
        f.write('META_STOP\n\n')
        f.write('STATE_VECTOR\n')
        f.write(f'EPOCH          = {epoch}\n')
        f.write(f'X              = {x[0]:20.9f} [km]\n')
        f.write(f'Y              = {x[1]:20.9f} [km]\n')
        f.write(f'Z              = {x[2]:20.9f} [km]\n')
        f.write(f'X_DOT          = {xd[0]:20.12f} [km/s]\n')
        f.write(f'Y_DOT          = {xd[1]:20.12f} [km/s]\n')
        f.write(f'Z_DOT          = {xd[2]:20.12f} [km/s]\n')


def read_oem_rows(path):
    rows = []
    in_data = False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s == 'DATA_START':
                in_data = True; continue
            if s == 'DATA_STOP':
                in_data = False; continue
            if not in_data:
                continue
            vals = s.split()
            if len(vals) >= 7:
                t = datetime.fromisoformat(vals[0])
                x = [float(v) for v in vals[1:4]]
                xd = [float(v) for v in vals[4:7]]
                rows.append((t, x, xd))
    return rows


def osc_elements(x, xd, mu):
    r = math.sqrt(sum(c*c for c in x))
    v2 = sum(c*c for c in xd)
    energy = v2/2.0 - mu/r
    a = -mu/(2.0*energy)
    hx = x[1]*xd[2]-x[2]*xd[1]
    hy = x[2]*xd[0]-x[0]*xd[2]
    hz = x[0]*xd[1]-x[1]*xd[0]
    h = math.sqrt(hx*hx+hy*hy+hz*hz)
    e = math.sqrt(max(0.0, 1.0 - h*h/(mu*a)))
    i = math.degrees(math.acos(max(-1.0, min(1.0, hz/h))))
    return a, e, i


def interp_at(rows, target_t):
    for k in range(len(rows) - 1):
        t0, x0, xd0 = rows[k]
        t1, x1, xd1 = rows[k+1]
        if t0 <= target_t <= t1:
            span = (t1 - t0).total_seconds()
            frac = 0.0 if span == 0 else (target_t-t0).total_seconds()/span
            x = [x0[j] + frac*(x1[j]-x0[j]) for j in range(3)]
            xd = [xd0[j] + frac*(xd1[j]-xd0[j]) for j in range(3)]
            return x, xd
    raise ValueError('target time not bracketed by trajectory buffer')


def run_driver(exe):
    for f in glob.glob('output/KSROP_*.oem'):
        os.remove(f)
    try:
        result = subprocess.run([exe], capture_output=True, text=True,
                                 timeout=180)
    except subprocess.TimeoutExpired:
        return None
    if result.returncode != 0:
        print(f'      RUNTIME ERROR: {result.stderr[:200]}')
        return None
    oem_files = sorted(glob.glob('output/KSROP_*.oem'))
    if not oem_files:
        return None
    return read_oem_rows(oem_files[-1])


def run_ch2(ctr, exe):
    print('=' * 64)
    print(' Chapter 2 (Oblateness): 100-rev zonal-only, J2 & J36')
    print('=' * 64)
    for case, (x0, xd0) in CH2_CASES.items():
        print(f'\n  --- Case {case} ---')
        for jt in CH2_TRUNCS:
            write_const(MU_CH2, RE_CH2, jt)
            write_input(CH2_NREV, CH2_ISTEP)
            write_opm(x0, xd0)
            rows = run_driver(exe)
            if not rows:
                check(ctr, f'J<={jt} case {case}: driver ran', False)
                continue
            xf, xdf = rows[-1][1], rows[-1][2]
            a, e, _ = osc_elements(xf, xdf, MU_CH2)
            a_ref, e_ref = CH2_REF[case][jt]
            rel_a = abs(a - a_ref) / a_ref
            rel_e = abs(e - e_ref) / e_ref
            check(ctr, f'J<={jt:2d} a rel err {rel_a:.2e} < {CH2_TOL_A_REL:.0e}',
                  rel_a < CH2_TOL_A_REL,
                  f'a={a:.4f} vs thesis {a_ref:.4f} km')
            check(ctr, f'J<={jt:2d} e rel err {rel_e:.2e} < {CH2_TOL_E_REL:.0e}',
                  rel_e < CH2_TOL_E_REL,
                  f'e={e:.7f} vs thesis {e_ref:.7f}')


def run_ch3(ctr, exe):
    print('\n' + '=' * 64)
    print(' Chapter 3 (Flattening): 22h zonal-only, GEM-T2, J6/J10/J19')
    print('=' * 64)
    for case, (x0, xd0, period_min) in CH3_CASES.items():
        a0, e0, i0 = osc_elements(x0, xd0, MU_CH2)
        nrev = int(math.ceil(CH3_DUR_H*60.0/period_min)) + 1
        print(f'\n  --- Case {case} (nrev={nrev}) ---')
        for jt in CH3_TRUNCS:
            write_const(MU_CH2, RE_CH2, jt)
            write_input(nrev, CH3_ISTEP)
            write_opm(x0, xd0)
            rows = run_driver(exe)
            if not rows:
                check(ctr, f'J<={jt} case {case}: driver ran', False)
                continue
            t_target = rows[0][0] + timedelta(hours=CH3_DUR_H)
            x, xd = interp_at(rows, t_target)
            a, e, i = osc_elements(x, xd, MU_CH2)
            da, de, di = a-a0, e-e0, i-i0
            da_ref, de_ref, di_ref = CH3_REF[case][jt]
            check(ctr,
                  f'J<={jt:2d} da err {abs(da-da_ref):.3f} km '
                  f'< {CH3_TOL_DA_KM} km',
                  abs(da - da_ref) < CH3_TOL_DA_KM,
                  f'da={da:+.4f} vs thesis {da_ref:+.4f} km')
            check(ctr,
                  f'J<={jt:2d} de err {abs(de-de_ref):.2e} < {CH3_TOL_DE:.0e}',
                  abs(de - de_ref) < CH3_TOL_DE,
                  f'de={de:+.3e} vs thesis {de_ref:+.3e}')
            check(ctr,
                  f'J<={jt:2d} di err {abs(di-di_ref):.2e} deg '
                  f'< {CH3_TOL_DI_DEG:.0e} deg',
                  abs(di - di_ref) < CH3_TOL_DI_DEG,
                  f'di={di:+.3e} vs thesis {di_ref:+.3e} deg')


def main():
    exe = sys.argv[1] if len(sys.argv) > 1 else 'driver_KS.exe'
    exe = os.path.abspath(exe)
    if not os.path.isfile(exe):
        print(f'ERROR: executable not found: {exe}')
        sys.exit(1)

    backup = {}
    for fname in ['input/const_new.dat', 'input/input.dat', 'input/input.opm']:
        bak = fname + '.bak_xjr'
        if os.path.isfile(fname):
            shutil.copy(fname, bak)
            backup[fname] = bak

    geo_existed = os.path.isfile(GEO_TARGET)
    geo_bak = GEO_TARGET + '.bak_xjr'
    if geo_existed:
        shutil.copy(GEO_TARGET, geo_bak)

    ctr = Counter()
    try:
        shutil.copy(FIXTURE_CH2, GEO_TARGET)
        run_ch2(ctr, exe)
        shutil.copy(FIXTURE_CH3, GEO_TARGET)
        run_ch3(ctr, exe)
    finally:
        if geo_existed:
            shutil.copy(geo_bak, GEO_TARGET)
            os.remove(geo_bak)
        elif os.path.isfile(GEO_TARGET):
            os.remove(GEO_TARGET)
        for fname, bak in backup.items():
            shutil.copy(bak, fname)
            os.remove(bak)

    print('\n' + '=' * 64)
    print(f' Total checks: {ctr.total}   Passed: {ctr.npass}'
          f'   Failed: {ctr.nfail}')
    print('=' * 64)
    sys.exit(0 if ctr.nfail == 0 else 1)


if __name__ == '__main__':
    main()
