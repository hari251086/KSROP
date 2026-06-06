#!/usr/bin/env python3
"""
test_driver.py  -  Integration test for driver_KS.F

Runs the compiled propagator with a pure two-body setup (all
perturbations disabled) and verifies physical invariants in
the output files.

Checks
------
  1. state.out has the expected number of rows
  2. Orbital energy conserved  (rel err < 1e-8)
  3. Angular momentum conserved (rel err < 1e-8)
  4. Orbit closure  |dr| after 1 revolution  < 0.1 km
  5. KS bilinear relation Br ~ 0 at initial and final epoch
  6. Semi-major axis consistent with vis-viva  (rel err < 1e-6)

Usage
-----
  python test_driver.py [executable]
  executable: compiled driver binary (default: driver_KS.exe)
"""

import subprocess, os, shutil, math, sys

MU = 3.986004415e5      # km3/s2  (matches const_new.dat)

# ----------------------------------------------------------------
# Vector helpers
# ----------------------------------------------------------------
def mag(v):
    return math.sqrt(sum(x*x for x in v))

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]]

def orbital_energy(x, xd):
    return dot(xd, xd)/2.0 - MU/mag(x)

def semimajor(x, xd):
    E = orbital_energy(x, xd)
    return -MU / (2.0*E)

def ang_momentum_mag(x, xd):
    return mag(cross(x, xd))

# ----------------------------------------------------------------
# Parsers
# ----------------------------------------------------------------
def read_state_out(path):
    """
    state.out format (format 100):
    t  x1 x2 x3  xd1 xd2 xd3  Br
    """
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.split()
            if len(vals) >= 8:
                t  = float(vals[0])
                x  = [float(v) for v in vals[1:4]]
                xd = [float(v) for v in vals[4:7]]
                Br = float(vals[7])
                rows.append((t, x, xd, Br))
    return rows

def read_kepler_out(path):
    """
    kepler.out format (format 102):
    iE  a  e  i  RAAN  AOP  M  EA
    """
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.split()
            if len(vals) >= 3:
                a = float(vals[1])
                e = float(vals[2])
                rows.append((a, e))
    return rows

# ----------------------------------------------------------------
# Assertion helper
# ----------------------------------------------------------------
NPASS = 0
NFAIL = 0

def check(name, ok, detail=''):
    global NPASS, NFAIL
    tag = 'PASS' if ok else 'FAIL'
    print(f'  [{tag}]  {name}')
    if not ok:
        print(f'           {detail}')
        NFAIL += 1
    else:
        NPASS += 1

# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------
def main():
    exe = sys.argv[1] if len(sys.argv) > 1 else 'driver_KS.exe'
    # Resolve to absolute path so subprocess finds it on Windows
    exe = os.path.abspath(exe)
    if not os.path.isfile(exe):
        print(f'ERROR: executable not found: {exe}')
        print('  Compile with:  ifort driver_KS.F Subrouts.F Legendre.F'
              ' -o driver_KS.exe')
        sys.exit(1)

    # --------------------------------------------------------
    # Backup originals and write two-body test inputs
    # --------------------------------------------------------
    files = {'const_new.dat': '_const_new.bak',
             'input.dat':     '_input.bak'}

    for fname, bak in files.items():
        real = next((f for f in os.listdir('.')
                     if f.lower() == fname.lower()), None)
        if real:
            shutil.copy(real, bak)

    # Point-mass gravity only, no drag, 1 revolution, 360 steps
    # Same initial state as production input.DAT
    with open('const_new.dat', 'w') as f:
        f.write('3.986004415D5 6378.1363D0 1.495978707d08'
                ' 1.32712440018d11 4.902801076d3\n')
        f.write('0 0 0\n')      # ngeo_deg=0 => point mass, no EGM2008 file needed

    with open('input.dat', 'w') as f:
        f.write('0.0 -5888.9727 -3400.0\n')
        f.write('9.5 0.0 0.0\n')
        f.write('1 360 1d-15\n')
        f.write('2016 09 20 0 0 0\n')
        f.write('0 0 0\n')                                # all forces off
        f.write('50.0 0 7.2921150d-5 3.35281066d-3 1.0\n')  # IDRAG=0

    # --------------------------------------------------------
    # Run the driver
    # --------------------------------------------------------
    print('=' * 50)
    print(' driver_KS.F  Integration Test  (2-body)')
    print('=' * 50)

    try:
        result = subprocess.run(
            [exe], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        print('TIMEOUT: driver ran for more than 120 s')
        sys.exit(1)

    if result.returncode != 0:
        print('Driver returned non-zero exit code:')
        print(result.stderr or result.stdout)
        sys.exit(1)

    # --------------------------------------------------------
    # Parse outputs
    # --------------------------------------------------------
    states  = read_state_out('state.out')
    keplels = read_kepler_out('kepler.out')

    # --------------------------------------------------------
    # Test 1: correct number of output rows
    # --------------------------------------------------------
    check('state.out has 2 rows (initial + 1 revolution)',
          len(states) == 2,
          f'got {len(states)} rows')

    if len(states) < 2:
        print('Cannot continue without 2 state rows.')
        sys.exit(1)

    t0, x0, xd0, Br0 = states[0]
    t1, x1, xd1, Br1 = states[1]

    # --------------------------------------------------------
    # Test 2: orbital energy conservation
    # --------------------------------------------------------
    E0 = orbital_energy(x0, xd0)
    E1 = orbital_energy(x1, xd1)
    dE_rel = abs(E1 - E0) / abs(E0)
    check(f'Orbital energy conserved   (rel err < 1e-8, got {dE_rel:.2e})',
          dE_rel < 1e-8,
          f'E0={E0:.8f}  E1={E1:.8f}')

    # --------------------------------------------------------
    # Test 3: angular momentum conservation
    # --------------------------------------------------------
    h0 = ang_momentum_mag(x0, xd0)
    h1 = ang_momentum_mag(x1, xd1)
    dh_rel = abs(h1 - h0) / abs(h0)
    check(f'Angular momentum conserved (rel err < 1e-8, got {dh_rel:.2e})',
          dh_rel < 1e-8,
          f'h0={h0:.8f}  h1={h1:.8f}')

    # --------------------------------------------------------
    # Test 4: orbit closure after 1 revolution
    # --------------------------------------------------------
    dx    = [x1[i] - x0[i] for i in range(3)]
    dxd   = [xd1[i] - xd0[i] for i in range(3)]
    dr    = mag(dx)
    dv    = mag(dxd)
    check(f'Orbit closure |dr| < 0.1 km        (got {dr:.4e} km)',
          dr < 0.1,
          f'dr={dx}')
    check(f'Orbit closure |dv| < 1e-5 km/s     (got {dv:.4e} km/s)',
          dv < 1e-5,
          f'dv={dxd}')

    # --------------------------------------------------------
    # Test 5: KS bilinear relation Br ~ 0
    # --------------------------------------------------------
    check(f'KS bilinear Br ~ 0 at t0           (|Br|={abs(Br0):.2e})',
          abs(Br0) < 1e-6,
          f'Br0={Br0}')
    check(f'KS bilinear Br ~ 0 at t1           (|Br|={abs(Br1):.2e})',
          abs(Br1) < 1e-6,
          f'Br1={Br1}')

    # --------------------------------------------------------
    # Test 6: semi-major axis from kepler.out vs vis-viva
    # --------------------------------------------------------
    a_visviva = semimajor(x0, xd0)
    if keplels:
        a_kep_0 = keplels[0][0]
        a_kep_1 = keplels[-1][0]
        da0_rel = abs(a_kep_0 - a_visviva) / a_visviva
        da1_rel = abs(a_kep_1 - a_kep_0)  / a_kep_0
        check(f'a (kepler.out) matches vis-viva    (rel err < 1e-6, got {da0_rel:.2e})',
              da0_rel < 1e-6,
              f'a_visviva={a_visviva:.4f}  a_kep={a_kep_0:.4f}')
        check(f'a conserved over 1 revolution      (rel err < 1e-6, got {da1_rel:.2e})',
              da1_rel < 1e-6,
              f'a0={a_kep_0:.4f}  a1={a_kep_1:.4f}')
    else:
        check('kepler.out readable', False, 'no rows found')

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------
    print('=' * 50)
    print(f' Total: {NPASS+NFAIL}   Passed: {NPASS}   Failed: {NFAIL}')
    print('=' * 50)

    # Restore original input files
    for fname, bak in files.items():
        if os.path.isfile(bak):
            shutil.copy(bak, fname)
            os.remove(bak)

    sys.exit(0 if NFAIL == 0 else 1)


if __name__ == '__main__':
    main()
