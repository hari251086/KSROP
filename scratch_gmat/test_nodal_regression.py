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

# Match KSJLSDN.F Verification/27845 case as closely as practical:
# SA=815.94725 km (perigee height), EC=0.0020592, AI=98.714729 deg,
# OMGA=210.16913, SOMGA=53.735756, AM=0 (start at perigee)
# Convert Kepler elements -> Cartesian via a quick two-body formula (J2 off)
mu = 3.986004415e5
Re = 6378.1363
hp_km = 815.94725
ec = 0.0020592
inc = math.radians(98.714729)
raan = math.radians(210.16913)
argp = math.radians(53.735756)
ta = 0.0  # true anomaly ~ mean anomaly at perigee (AM=0)

rp = Re + hp_km
a = rp/(1-ec)
p = a*(1-ec**2)
r = p/(1+ec*math.cos(ta))

# perifocal position/velocity
x_pf = r*math.cos(ta)
y_pf = r*math.sin(ta)
vx_pf = -math.sqrt(mu/p)*math.sin(ta)
vy_pf = math.sqrt(mu/p)*(ec+math.cos(ta))

cO, sO = math.cos(raan), math.sin(raan)
ci, si = math.cos(inc), math.sin(inc)
cw, sw = math.cos(argp), math.sin(argp)

R11 = cO*cw - sO*sw*ci
R12 = -cO*sw - sO*cw*ci
R21 = sO*cw + cO*sw*ci
R22 = -sO*sw + cO*cw*ci
R31 = sw*si
R32 = cw*si

x = R11*x_pf + R12*y_pf
y = R21*x_pf + R22*y_pf
z = R31*x_pf + R32*y_pf
vx = R11*vx_pf + R12*vy_pf
vy = R21*vx_pf + R22*vy_pf
vz = R31*vx_pf + R32*vy_pf

print(f"IC: x=({x:.6f},{y:.6f},{z:.6f})  v=({vx:.6f},{vy:.6f},{vz:.6f})  r={math.sqrt(x*x+y*y+z*z):.4f}")

with open('input/const_new.dat', 'w') as f:
    f.write(f'{mu}D0 {Re}D0 1.495978707d08 1.32712440018d11 4.902801076d3\n')
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
    f.write(f'X              = {x:.9f} [km]\n')
    f.write(f'Y              = {y:.9f} [km]\n')
    f.write(f'Z              = {z:.9f} [km]\n')
    f.write(f'X_DOT          = {vx:.9f} [km/s]\n')
    f.write(f'Y_DOT          = {vy:.9f} [km/s]\n')
    f.write(f'Z_DOT          = {vz:.9f} [km/s]\n')

with open('input/input.dat', 'w') as f:
    f.write('1 360 1d-15\n')
    f.write('1 0 0\n')
    f.write('50.0 0 7.2921150d-5 3.35281066d-3 1.0\n')
    f.write('1.2 0.01 0 1\n')

shutil.move(FULL, BACKUP)
shutil.copy(TRUNC, FULL)
try:
    for p_ in glob.glob('output/KSROP_*.oem'):
        os.remove(p_)
    r_ = subprocess.run([os.path.abspath('driver_KS.exe')], capture_output=True, text=True, timeout=120)
    oems = sorted(glob.glob('output/KSROP_*.oem'))
    shutil.copy(oems[-1], 'scratch_gmat/oem_nodal.oem')
finally:
    os.remove(FULL)
    shutil.move(BACKUP, FULL)

for fname, bak in files.items():
    if os.path.isfile(bak):
        shutil.copy(bak, fname)
        os.remove(bak)

# --- extract RAAN at start and end via car2oe-equivalent (vector method) ---
def state_to_raan(x,y,z,vx,vy,vz):
    rvec = [x,y,z]; vvec=[vx,vy,vz]
    h = [rvec[1]*vvec[2]-rvec[2]*vvec[1],
         rvec[2]*vvec[0]-rvec[0]*vvec[2],
         rvec[0]*vvec[1]-rvec[1]*vvec[0]]
    nvec = [-h[1], h[0], 0.0]
    nmag = math.sqrt(nvec[0]**2+nvec[1]**2)
    raan = math.acos(nvec[0]/nmag)
    if nvec[1] < 0:
        raan = 2*math.pi - raan
    return math.degrees(raan)

rows = []
with open('scratch_gmat/oem_nodal.oem') as f:
    in_data = False
    for line in f:
        s = line.strip()
        if s == 'DATA_START': in_data = True; continue
        if s == 'DATA_STOP': in_data = False; continue
        if not in_data: continue
        v = s.split()
        if len(v) >= 7:
            rows.append(([float(a) for a in v[1:4]], [float(a) for a in v[4:7]]))

x0 = rows[0][0] + rows[0][1]
x1 = rows[-1][0] + rows[-1][1]
raan0 = state_to_raan(*x0)
raan1 = state_to_raan(*x1)
print(f"RAAN start={raan0:.6f}  RAAN end={raan1:.6f}  dOmega={raan1-raan0:.6f} deg")

r0 = math.sqrt(sum(c*c for c in rows[0][0]))
v0 = math.sqrt(sum(c*c for c in rows[0][1]))
a_ = 1.0/(2.0/r0 - v0*v0/mu)
n_ = math.sqrt(mu/a_**3)
T_ = 2*math.pi/n_
J2 = 1.082626173852223e-3
p_ = a_*(1-ec**2)
dOmega_dt = -1.5*n_*J2*(Re/p_)**2*math.cos(inc)
dOmega_pred = math.degrees(dOmega_dt*T_)
print(f"a={a_:.4f} km  T={T_:.4f} s  predicted dOmega/orbit={dOmega_pred:.6f} deg")
