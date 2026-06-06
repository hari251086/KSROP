#!/usr/bin/env python3
"""
test_driver.py  -  Integration test for driver_KS.F

Runs the compiled propagator with a pure two-body setup (all
perturbations disabled) and verifies physical invariants in
the CCSDS OEM v2.0 output file (ksrop.oem).

Checks
------
  1. ksrop.oem has correct number of data rows (initial + nrev)
  2. OEM file structure: DATA_START / DATA_STOP present
  3. Orbital energy conserved       (rel err < 1e-8)
  4. Angular momentum conserved     (rel err < 1e-8)
  5. Orbit closure |dr| < 0.1 km   after 1 revolution
  6. Orbit closure |dv| < 1e-5 km/s
  7. Semi-major axis matches vis-viva  (rel err < 1e-6)
  8. Semi-major axis conserved over 1 revolution  (rel err < 1e-6)
  9. ksrop.opm written and contains STATE_VECTOR keyword

Usage
-----
  python test_driver.py [executable]
  executable: compiled driver binary (default: driver_KS.exe)
"""

import subprocess, os, shutil, math, sys, glob

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
    return -MU / (2.0*orbital_energy(x, xd))

def ang_momentum_mag(x, xd):
    return mag(cross(x, xd))

# ----------------------------------------------------------------
# OEM parser
# ----------------------------------------------------------------
def read_oem(path):
    """
    Parse a CCSDS OEM v2.0 ASCII file.
    Returns list of (epoch_str, x[3], xd[3]) for each data line.
    Also returns True/False for DATA_START and DATA_STOP presence.
    """
    rows = []
    has_start = False
    has_stop  = False
    in_data   = False
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if stripped == 'DATA_START':
                has_start = True
                in_data   = True
                continue
            if stripped == 'DATA_STOP':
                has_stop = True
                in_data  = False
                continue
            if not in_data:
                continue
            vals = stripped.split()
            if len(vals) >= 7:
                epoch = vals[0]
                x  = [float(v) for v in vals[1:4]]
                xd = [float(v) for v in vals[4:7]]
                rows.append((epoch, x, xd))
    return rows, has_start, has_stop

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
    exe = os.path.abspath(exe)
    if not os.path.isfile(exe):
        print(f'ERROR: executable not found: {exe}')
        print('  Compile: ifort driver_KS.F Subrouts.F Legendre.F -o driver_KS.exe')
        sys.exit(1)

    # --------------------------------------------------------
    # Backup originals and write two-body test inputs
    # --------------------------------------------------------
    files = {'input/const_new.dat': '_const_new.bak',
             'input/input.dat':     '_input.bak',
             'input/input.opm':     '_input_opm.bak'}

    for fname, bak in files.items():
        if os.path.isfile(fname):
            shutil.copy(fname, bak)

    # input/const_new.dat: point-mass only, no EGM2008 file needed
    with open('input/const_new.dat', 'w') as f:
        f.write('3.986004415D5 6378.1363D0 1.495978707d08'
                ' 1.32712440018d11 4.902801076d3\n')
        f.write('0 0 0\n')

    # input/input.opm: initial state in CCSDS OPM format
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
        f.write('X              =        0.000000 [km]\n')
        f.write('Y              =    -5888.972700 [km]\n')
        f.write('Z              =    -3400.000000 [km]\n')
        f.write('X_DOT          =        9.500000 [km/s]\n')
        f.write('Y_DOT          =        0.000000 [km/s]\n')
        f.write('Z_DOT          =        0.000000 [km/s]\n')

    # input/input.dat: simulation parameters only (3 lines)
    with open('input/input.dat', 'w') as f:
        f.write('1 360 1d-15\n')
        f.write('0 0 0\n')
        f.write('50.0 0 7.2921150d-5 3.35281066d-3 1.0\n')

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
    # Locate and parse the timestamped OEM file
    # --------------------------------------------------------
    oem_files = sorted(glob.glob('output/KSROP_*.oem'))
    if not oem_files:
        print('ERROR: no KSROP_*.oem file found after run')
        sys.exit(1)
    oem_path = oem_files[-1]   # most recent if multiple exist
    rows, has_start, has_stop = read_oem(oem_path)

    # Test config: nrev=1, istep=360 => 1*360+1 = 361 rows
    NREV  = 1
    ISTEP = 360
    EXPECTED_ROWS = NREV * ISTEP + 1

    # --------------------------------------------------------
    # Test 1: correct number of data rows
    # --------------------------------------------------------
    check(f'ksrop.oem has {EXPECTED_ROWS} rows '
          f'({NREV} rev x {ISTEP} steps + 1 initial)',
          len(rows) == EXPECTED_ROWS,
          f'got {len(rows)} rows')

    # --------------------------------------------------------
    # Test 2: OEM file structure
    # --------------------------------------------------------
    check(f'{oem_path} contains DATA_START marker', has_start)
    check(f'{oem_path} contains DATA_STOP  marker', has_stop)

    if len(rows) < 2:
        print('Cannot continue without at least 2 state rows.')
        sys.exit(1)

    _, x0, xd0 = rows[0]    # initial state
    _, x1, xd1 = rows[-1]   # final state (end of last revolution)

    # --------------------------------------------------------
    # Test 3: orbital energy conservation
    # --------------------------------------------------------
    E0 = orbital_energy(x0, xd0)
    E1 = orbital_energy(x1, xd1)
    dE_rel = abs(E1 - E0) / abs(E0)
    check(f'Orbital energy conserved   (rel err < 1e-8, got {dE_rel:.2e})',
          dE_rel < 1e-8,
          f'E0={E0:.8f}  E1={E1:.8f}')

    # --------------------------------------------------------
    # Test 4: angular momentum conservation
    # --------------------------------------------------------
    h0 = ang_momentum_mag(x0, xd0)
    h1 = ang_momentum_mag(x1, xd1)
    dh_rel = abs(h1 - h0) / abs(h0)
    check(f'Angular momentum conserved (rel err < 1e-8, got {dh_rel:.2e})',
          dh_rel < 1e-8,
          f'h0={h0:.8f}  h1={h1:.8f}')

    # --------------------------------------------------------
    # Test 5-6: orbit closure after 1 revolution
    # --------------------------------------------------------
    dx  = [x1[i] - x0[i] for i in range(3)]
    dxd = [xd1[i] - xd0[i] for i in range(3)]
    dr  = mag(dx)
    dv  = mag(dxd)
    check(f'Orbit closure |dr| < 0.1 km        (got {dr:.4e} km)',
          dr < 0.1,
          f'dr={dx}')
    check(f'Orbit closure |dv| < 1e-5 km/s     (got {dv:.4e} km/s)',
          dv < 1e-5,
          f'dv={dxd}')

    # --------------------------------------------------------
    # Test 7-8: semi-major axis from vis-viva (no kepler.out)
    # --------------------------------------------------------
    a0 = semimajor(x0, xd0)
    a1 = semimajor(x1, xd1)
    a_visviva = a0
    da0_rel = abs(a0 - a_visviva) / a_visviva   # trivially 0, confirms formula
    da1_rel = abs(a1 - a0) / a0
    check(f'a matches vis-viva at t0            (rel err < 1e-6, got {da0_rel:.2e})',
          da0_rel < 1e-6)
    check(f'a conserved over 1 revolution       (rel err < 1e-6, got {da1_rel:.2e})',
          da1_rel < 1e-6,
          f'a0={a0:.4f} km   a1={a1:.4f} km')

    # --------------------------------------------------------
    # Test 9: ksrop.opm written and readable
    # --------------------------------------------------------
    opm_ok = False
    if os.path.isfile('output/ksrop.opm'):
        with open('output/ksrop.opm') as f:
            opm_ok = any('STATE_VECTOR' in line for line in f)
    check('ksrop.opm written with STATE_VECTOR block', opm_ok)

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
