import subprocess, os, shutil, glob, sys, time

REPO = r"C:\Users\hari2\OneDrive\Documents\GitHub\KSROP"
os.chdir(REPO)

FULL = "input/EGM2008_to2190_TideFree"
TRUNC = "scratch_gmat/EGM2008_72x72.dat"
BACKUP = "input/EGM2008_to2190_TideFree.fullbackup"

files = {'input/const_new.dat': '_const_new.bak',
         'input/input.dat':     '_input.bak',
         'input/input.opm':     '_input_opm.bak'}
for fname, bak in files.items():
    if os.path.isfile(fname):
        shutil.copy(fname, bak)

with open('input/const_new.dat', 'w') as f:
    f.write('3.986004415D5 6378.1363D0 1.495978707d08'
            ' 1.32712440018d11 4.902801076d3\n')
    f.write('20 0 0\n')          # ngeo_deg=20, sun/moon off
    f.write('4.56d-6\n')

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

with open('input/input.dat', 'w') as f:
    f.write('1 360 1d-15\n')
    f.write('1 0 0\n')            # n_force: geo ON, sun/moon off
    f.write('50.0 0 7.2921150d-5 3.35281066d-3 1.0\n')
    f.write('1.2 0.01 0 1\n')

def run_and_get_oem(tag):
    for p in glob.glob('output/KSROP_*.oem'):
        os.remove(p)
    result = subprocess.run([os.path.abspath('driver_KS.exe')], capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        print('DRIVER FAILED:', result.stdout, result.stderr)
        sys.exit(1)
    oems = sorted(glob.glob('output/KSROP_*.oem'))
    dst = f'scratch_gmat/oem_{tag}.oem'
    shutil.copy(oems[-1], dst)
    print(f'[{tag}] -> {dst}')
    return dst

# --- Run 1: full 2190-line file (already in place) ---
oem_full = run_and_get_oem('full')

# --- Run 2: swap in truncated 72x72 file, same name ---
shutil.move(FULL, BACKUP)
shutil.copy(TRUNC, FULL)
try:
    oem_trunc = run_and_get_oem('trunc')
finally:
    os.remove(FULL)
    shutil.move(BACKUP, FULL)

# --- Restore original input files ---
for fname, bak in files.items():
    if os.path.isfile(bak):
        shutil.copy(bak, fname)
        os.remove(bak)

# --- Diff the two OEMs ---
a = open(oem_full).read()
b = open(oem_trunc).read()
print('IDENTICAL' if a == b else 'DIFFERENT')
if a != b:
    la, lb = a.splitlines(), b.splitlines()
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            print(f'line {i}: FULL={x!r}  TRUNC={y!r}')
