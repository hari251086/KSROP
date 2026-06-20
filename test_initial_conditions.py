#!/usr/bin/env python3
"""
test_initial_conditions.py - Multi-case integration tests

Runs the propagator under two-body AND full-dynamics configurations
for a range of orbital regimes.

Phase 1 — Two-body (all perturbations off):
  Verifies energy, angular momentum, orbit closure, and semi-major
  axis conservation to machine precision.

Phase 2 — Full dynamics (J10 geopotential, luni-solar, drag, SRP):
  Verifies the propagator completes without divergence, all states
  remain finite, altitude stays physical, and osculating energy /
  angular momentum / semi-major axis drift stay bounded (not
  conserved, but within physically expected envelopes).

Usage
-----
  python test_initial_conditions.py [executable]
"""

import subprocess, os, shutil, math, sys, glob

MU      = 3.986004415e5   # km3/s2
R_EARTH = 6378.1363       # km

# -------------------------------------------------------------------
# Vector/orbital helpers
# -------------------------------------------------------------------
def mag(v):
    return math.sqrt(sum(c*c for c in v))

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def cross3(a, b):
    return [a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]]

def orbital_energy(x, xd):
    return dot(xd, xd)/2.0 - MU/mag(x)

def semimajor(x, xd):
    return -MU / (2.0*orbital_energy(x, xd))

def ang_momentum_mag(x, xd):
    return mag(cross3(x, xd))

def altitude(x):
    return mag(x) - R_EARTH

def is_finite(v):
    return all(math.isfinite(c) for c in v)

def oe2cart(a, e, i_deg, raan_deg, aop_deg, nu_deg):
    """Keplerian elements to Cartesian state (km, km/s)."""
    i    = math.radians(i_deg)
    raan = math.radians(raan_deg)
    aop  = math.radians(aop_deg)
    nu   = math.radians(nu_deg)

    p = a * (1.0 - e*e)
    r = p / (1.0 + e*math.cos(nu))

    rx_pf = r * math.cos(nu)
    ry_pf = r * math.sin(nu)
    vx_pf = -math.sqrt(MU/p) * math.sin(nu)
    vy_pf =  math.sqrt(MU/p) * (e + math.cos(nu))

    cosO, sinO = math.cos(raan), math.sin(raan)
    cosw, sinw = math.cos(aop),  math.sin(aop)
    cosi, sini = math.cos(i),    math.sin(i)

    l1 = cosO*cosw - sinO*sinw*cosi
    l2 = sinO*cosw + cosO*sinw*cosi
    l3 = sinw*sini
    m1 = -cosO*sinw - sinO*cosw*cosi
    m2 = -sinO*sinw + cosO*cosw*cosi
    m3 = cosw*sini

    x  = [l1*rx_pf + m1*ry_pf, l2*rx_pf + m2*ry_pf, l3*rx_pf + m3*ry_pf]
    xd = [l1*vx_pf + m1*vy_pf, l2*vx_pf + m2*vy_pf, l3*vx_pf + m3*vy_pf]
    return x, xd

# -------------------------------------------------------------------
# OEM parser
# -------------------------------------------------------------------
def read_oem(path):
    rows = []
    in_data = False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s == 'DATA_START':
                in_data = True
                continue
            if s == 'DATA_STOP':
                in_data = False
                continue
            if not in_data:
                continue
            vals = s.split()
            if len(vals) >= 7:
                x  = [float(v) for v in vals[1:4]]
                xd = [float(v) for v in vals[4:7]]
                rows.append((x, xd))
    return rows

# -------------------------------------------------------------------
# Assertion helper
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Test cases: (name, a_km, e, i_deg, raan_deg, aop_deg, nu_deg, nrev)
# -------------------------------------------------------------------
ORBITS = [
    ("LEO circular i=28.5",        6878.0,  0.001,  28.5,   0.0,   0.0,  45.0, 1),
    ("GTO high-e",                 24500.0, 0.73,    7.0,  45.0, 178.0,  45.0, 1),
    ("Polar circular",              7000.0, 0.001,  90.0,   0.0,   0.0,  45.0, 1),
    ("Retrograde i=135",            8000.0, 0.05,  135.0,  90.0,  45.0,  30.0, 1),
    ("Near-equatorial elliptical",  7500.0, 0.10,    5.0, 180.0,  90.0,  60.0, 1),
    ("Molniya-like e=0.74",        26600.0, 0.74,   63.4, 270.0, 270.0,  45.0, 1),
    ("GEO-like near-circular",     42164.0, 0.001,   0.1,   0.0,   0.0,  45.0, 1),
    ("SSO dawn-dusk",               7078.0, 0.001,  97.8,  90.0,   0.0, 120.0, 1),
    ("HEO e=0.6",                  20000.0, 0.60,   30.0,  60.0,  90.0,  90.0, 1),
    ("MEO GPS-like",               26560.0, 0.01,   55.0,   0.0,   0.0, 180.0, 1),
]

# -------------------------------------------------------------------
# Input-file writers
# -------------------------------------------------------------------
def write_opm(x, xd):
    with open('input/input.opm', 'w') as f:
        f.write('CCSDS_OPM_VERS = 2.0\n')
        f.write('CREATION_DATE  = 2016-09-20T00:00:00.000\n')
        f.write('ORIGINATOR     = KSROP\n\n')
        f.write('META_START\n')
        f.write('OBJECT_NAME    = SATELLITE\n')
        f.write('CENTER_NAME    = EARTH\n')
        f.write('REF_FRAME      = EME2000\n')
        f.write('TIME_SYSTEM    = UTC\n')
        f.write('META_STOP\n\n')
        f.write('STATE_VECTOR\n')
        f.write('EPOCH          = 2016-09-20T00:00:00.000\n')
        f.write(f'X              = {x[0]:20.9f} [km]\n')
        f.write(f'Y              = {x[1]:20.9f} [km]\n')
        f.write(f'Z              = {x[2]:20.9f} [km]\n')
        f.write(f'X_DOT          = {xd[0]:20.12f} [km/s]\n')
        f.write(f'Y_DOT          = {xd[1]:20.12f} [km/s]\n')
        f.write(f'Z_DOT          = {xd[2]:20.12f} [km/s]\n')

def write_two_body(nrev, istep=360):
    with open('input/const_new.dat', 'w') as f:
        f.write('3.986004415D5 6378.1363D0 1.495978707d08'
                ' 1.32712440018d11 4.902801076d3\n')
        f.write('0 0 0\n')
        f.write('4.56d-6\n')
    with open('input/input.dat', 'w') as f:
        f.write(f'{nrev} {istep} 1d-15\n')
        f.write('0 0 0\n')
        f.write('50.0 0 7.2921150d-5 3.35281066d-3 1.0\n')
        f.write('1.2 0.01 0 1\n')

def write_full_dynamics(nrev, istep=360, drag=False):
    """J10 geopotential + degree-2 sun/moon + optional drag + SRP."""
    with open('input/const_new.dat', 'w') as f:
        f.write('3.986004415D5 6378.1363D0 1.495978707d08'
                ' 1.32712440018d11 4.902801076d3\n')
        f.write('10 2 2\n')
        f.write('4.56d-6\n')
    idrag = 1 if drag else 0
    with open('input/input.dat', 'w') as f:
        f.write(f'{nrev} {istep} 1d-15\n')
        f.write('10 1 1\n')
        f.write(f'50.0 {idrag} 7.2921150d-5 3.35281066d-3 1.0\n')
        f.write('1.2 0.01 1 1\n')

# -------------------------------------------------------------------
# Run the driver and return OEM rows (or None on failure)
# -------------------------------------------------------------------
def run_driver(exe, label):
    for f in glob.glob('output/KSROP_*.oem'):
        os.remove(f)
    try:
        result = subprocess.run(
            [exe], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        print(f'      TIMEOUT ({label})')
        return None
    if result.returncode != 0:
        print(f'      RUNTIME ERROR ({label}): {result.stderr[:200]}')
        return None
    oem_files = sorted(glob.glob('output/KSROP_*.oem'))
    if not oem_files:
        print(f'      NO OEM OUTPUT ({label})')
        return None
    return read_oem(oem_files[-1])

# -------------------------------------------------------------------
# Phase 1: Two-body checks (exact conservation)
# -------------------------------------------------------------------
def run_two_body_checks(ctr, rows, expected_rows, failed):
    x0, xd0 = rows[0]
    xf, xdf = rows[-1]
    ok = True

    ok &= check(ctr, f'Row count = {expected_rows}',
                len(rows) == expected_rows,
                f'got {len(rows)}')

    E0 = orbital_energy(x0, xd0)
    Ef = orbital_energy(xf, xdf)
    dE = abs(Ef - E0) / abs(E0)
    ok &= check(ctr, f'Energy conserved        (rel err {dE:.2e})',
                dE < 1e-8)

    h0 = ang_momentum_mag(x0, xd0)
    hf = ang_momentum_mag(xf, xdf)
    dh = abs(hf - h0) / abs(h0)
    ok &= check(ctr, f'Ang. momentum conserved (rel err {dh:.2e})',
                dh < 1e-8)

    dr = mag([xf[j]-x0[j] for j in range(3)])
    dv = mag([xdf[j]-xd0[j] for j in range(3)])
    ok &= check(ctr, f'Orbit closure            (dr={dr:.4e}, dv={dv:.4e})',
                dr < 0.1 and dv < 1e-5)

    a0 = semimajor(x0, xd0)
    af = semimajor(xf, xdf)
    da = abs(af - a0) / a0
    ok &= check(ctr, f'SMA conserved            (rel err {da:.2e})',
                da < 1e-6)
    return ok

# -------------------------------------------------------------------
# Phase 2: Full-dynamics checks (bounded drift)
# -------------------------------------------------------------------
def run_full_dynamics_checks(ctr, rows, expected_rows, a0_nom,
                             allow_reentry=False):
    x0, xd0 = rows[0]
    xf, xdf = rows[-1]
    ok = True

    if allow_reentry:
        ok &= check(ctr, f'Propagation produced output ({len(rows)} rows)',
                    len(rows) >= 2)
    else:
        ok &= check(ctr, f'Row count = {expected_rows}',
                    len(rows) == expected_rows,
                    f'got {len(rows)}')

    all_finite = all(is_finite(x) and is_finite(xd) for x, xd in rows)
    ok &= check(ctr, 'All states finite (no NaN/Inf)', all_finite)

    min_alt = min(altitude(x) for x, xd in rows)
    max_alt = max(altitude(x) for x, xd in rows)
    if allow_reentry:
        ok &= check(ctr, f'Min altitude physical     (h={min_alt:.1f} km)',
                    min_alt > -R_EARTH)
    else:
        ok &= check(ctr, f'Min altitude > 0 km       (h={min_alt:.1f} km)',
                    min_alt > 0.0)

    # Osculating energy: should drift but stay within ~1% for 1 rev
    E0 = orbital_energy(x0, xd0)
    Ef = orbital_energy(xf, xdf)
    dE = abs(Ef - E0) / abs(E0)
    ok &= check(ctr, f'Energy drift bounded     (rel {dE:.4e})',
                dE < 1e-2,
                f'E0={E0:.6f}, Ef={Ef:.6f}')

    # Angular momentum: J2/luni-solar cause secular drift, but
    # should stay within a few percent over 1 revolution
    h0 = ang_momentum_mag(x0, xd0)
    hf = ang_momentum_mag(xf, xdf)
    dh = abs(hf - h0) / abs(h0)
    ok &= check(ctr, f'Ang. mom. drift bounded  (rel {dh:.4e})',
                dh < 1e-2,
                f'h0={h0:.4f}, hf={hf:.4f}')

    # SMA: osculating a should stay within ~5% of nominal for 1 rev
    a0 = semimajor(x0, xd0)
    af = semimajor(xf, xdf)
    da = abs(af - a0_nom) / a0_nom
    ok &= check(ctr, f'SMA physically plausible (da/a={da:.4e})',
                da < 5e-2,
                f'a_nom={a0_nom:.1f}, af={af:.1f}')

    return ok

# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
def main():
    exe = sys.argv[1] if len(sys.argv) > 1 else 'driver_KS.exe'
    exe = os.path.abspath(exe)
    if not os.path.isfile(exe):
        print(f'ERROR: executable not found: {exe}')
        sys.exit(1)

    backup = {}
    for fname in ['input/const_new.dat', 'input/input.dat', 'input/input.opm']:
        bak = fname + '.bak_ic'
        if os.path.isfile(fname):
            shutil.copy(fname, bak)
            backup[fname] = bak

    ctr = Counter()
    failed_cases = []
    istep = 360

    # ==============================================================
    # PHASE 1: Two-body
    # ==============================================================
    print('=' * 64)
    print(' Phase 1 — Two-Body (exact conservation)')
    print('=' * 64)

    for name, a, e, i, raan, aop, nu, nrev in ORBITS:
        x, xd = oe2cart(a, e, i, raan, aop, nu)
        write_opm(x, xd)
        write_two_body(nrev, istep)

        print(f'\n  --- {name} ---')
        print(f'      a={a:.1f} km  e={e:.4f}  i={i:.1f} deg')

        rows = run_driver(exe, name)
        if rows is None or len(rows) < 2:
            ctr.nfail += 5
            failed_cases.append(f'2B:{name}')
            if rows is not None:
                print(f'      TOO FEW ROWS ({len(rows)})')
            continue

        if not run_two_body_checks(ctr, rows, nrev*istep+1, failed_cases):
            failed_cases.append(f'2B:{name}')

    # ==============================================================
    # PHASE 2: Full dynamics
    # ==============================================================
    print('\n' + '=' * 64)
    print(' Phase 2 — Full Dynamics (J10 + Sun/Moon + SRP + drag)')
    print('=' * 64)

    for name, a, e, i, raan, aop, nu, nrev in ORBITS:
        x, xd = oe2cart(a, e, i, raan, aop, nu)
        perigee_alt = a*(1.0 - e) - R_EARTH
        drag_on = perigee_alt < 500.0
        write_opm(x, xd)
        write_full_dynamics(nrev, istep, drag=drag_on)

        drag_tag = '+drag' if drag_on else ''
        print(f'\n  --- {name} {drag_tag} ---')
        print(f'      a={a:.1f} km  e={e:.4f}  i={i:.1f} deg'
              f'  hp={perigee_alt:.0f} km')

        rows = run_driver(exe, name)
        if rows is None or len(rows) < 2:
            ctr.nfail += 6
            failed_cases.append(f'FD:{name}')
            if rows is not None:
                print(f'      TOO FEW ROWS ({len(rows)})')
            continue

        allow_reentry = drag_on and perigee_alt < 200.0
        if not run_full_dynamics_checks(ctr, rows, nrev*istep+1, a,
                                        allow_reentry=allow_reentry):
            failed_cases.append(f'FD:{name}')

    # ==============================================================
    # Summary
    # ==============================================================
    # Restore originals
    for fname, bak in backup.items():
        if os.path.isfile(bak):
            shutil.copy(bak, fname)
            os.remove(bak)

    print('\n' + '=' * 64)
    print(f' Total checks: {ctr.total}   Passed: {ctr.npass}'
          f'   Failed: {ctr.nfail}')
    if failed_cases:
        print(f' Failed: {", ".join(failed_cases)}')
    print('=' * 64)

    sys.exit(0 if ctr.nfail == 0 else 1)


if __name__ == '__main__':
    main()
