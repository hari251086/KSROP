import subprocess, os, shutil, glob, math

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
    f.write('20 0 0\n')
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

shutil.move(FULL, BACKUP)
shutil.copy(TRUNC, FULL)
try:
    for istep in [360, 3600, 36000]:
        with open('input/input.dat', 'w') as f:
            f.write(f'1 {istep} 1d-15\n')
            f.write('1 0 0\n')
            f.write('50.0 0 7.2921150d-5 3.35281066d-3 1.0\n')
            f.write('1.2 0.01 0 1\n')
        for p in glob.glob('output/KSROP_*.oem'):
            os.remove(p)
        r = subprocess.run([os.path.abspath('driver_KS.exe')], capture_output=True, text=True, timeout=300)
        oems = sorted(glob.glob('output/KSROP_*.oem'))
        rows = []
        with open(oems[-1]) as f:
            in_data = False
            for line in f:
                s = line.strip()
                if s == 'DATA_START': in_data = True; continue
                if s == 'DATA_STOP': in_data = False; continue
                if not in_data: continue
                v = s.split()
                if len(v) >= 7:
                    rows.append(([float(x) for x in v[1:4]], [float(x) for x in v[4:7]]))
        x0, xd0 = rows[0]
        x1, xd1 = rows[-1]
        dr = math.sqrt(sum((a-b)**2 for a,b in zip(x0,x1)))
        dv = math.sqrt(sum((a-b)**2 for a,b in zip(xd0,xd1)))
        print(f"istep={istep:6d}  closure dr={dr:.6f} km   dv={dv:.6e} km/s")
finally:
    os.remove(FULL)
    shutil.move(BACKUP, FULL)

for fname, bak in files.items():
    if os.path.isfile(bak):
        shutil.copy(bak, fname)
        os.remove(bak)
