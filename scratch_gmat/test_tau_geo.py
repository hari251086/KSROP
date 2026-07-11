import subprocess, os, shutil, glob, math, datetime

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
    f.write('2 0 0\n')
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

istep = 3600
with open('input/input.dat', 'w') as f:
    f.write(f'1 {istep} 1d-15\n')
    f.write('1 0 0\n')
    f.write('50.0 0 7.2921150d-5 3.35281066d-3 1.0\n')
    f.write('1.2 0.01 0 1\n')

shutil.move(FULL, BACKUP)
shutil.copy(TRUNC, FULL)
try:
    for p in glob.glob('output/KSROP_*.oem'):
        os.remove(p)
    r = subprocess.run([os.path.abspath('driver_KS.exe')], capture_output=True, text=True, timeout=300)
    oems = sorted(glob.glob('output/KSROP_*.oem'))
    shutil.copy(oems[-1], 'scratch_gmat/oem_j2_fine.oem')
finally:
    os.remove(FULL)
    shutil.move(BACKUP, FULL)

for fname, bak in files.items():
    if os.path.isfile(bak):
        shutil.copy(bak, fname)
        os.remove(bak)

# --- parse ---
rows = []
with open('scratch_gmat/oem_j2_fine.oem') as f:
    in_data = False
    for line in f:
        s = line.strip()
        if s == 'DATA_START': in_data = True; continue
        if s == 'DATA_STOP': in_data = False; continue
        if not in_data: continue
        v = s.split()
        if len(v) >= 7:
            rows.append((v[0], [float(x) for x in v[1:4]], [float(x) for x in v[4:7]]))

mu = 3.986004415e5
t0 = datetime.datetime.strptime(rows[0][0], "%Y-%m-%dT%H:%M:%S.%f")

# w_kep at t=0 (two-body reference frequency)
x0, xd0 = rows[0][1], rows[0][2]
r0 = math.sqrt(sum(c*c for c in x0))
v0_2 = sum(c*c for c in xd0)
w_kep = math.sqrt(0.5*(mu/r0 - v0_2/2))
dE0 = 2*math.pi/istep
ds = dE0/w_kep   # constant fictitious-time step (Gam cancels as derived)
print("w_kep=", w_kep, " ds=", ds)

# ground truth cumulative physical time via Simpson's rule on r(s)
r_vals = [math.sqrt(sum(c*c for c in x)) for _, x, _ in rows]
T_true = [0.0]*len(rows)
acc = 0.0
for k in range(1, len(rows)):
    # trapezoid (simple, sufficient at fine istep)
    acc += 0.5*(r_vals[k-1]+r_vals[k])*ds
    T_true[k] = acc

# code-reported elapsed time
T_code = []
for e,_,_ in rows:
    t = datetime.datetime.strptime(e, "%Y-%m-%dT%H:%M:%S.%f")
    T_code.append((t-t0).total_seconds())

for i in [0, 360, 900, 1800, 2700, 3240, 3600]:
    print(f"step={i:5d}  T_code={T_code[i]:14.6f}  T_true={T_true[i]:14.6f}  diff={T_code[i]-T_true[i]:12.6e}")
