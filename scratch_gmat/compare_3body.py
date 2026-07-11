gmat_path = r"E:\Softwares\gmat-win-R2026a\output\crosscheck_3body_gmat.txt"
blocks = []
cur = []
with open(gmat_path) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith("DefaultSC.UTCGregorian"):
            if cur:
                blocks.append(cur)
            cur = []
            continue
        cur.append(line)
if cur:
    blocks.append(cur)

gmat_rows = []
for b in blocks:
    parts = b[-1].split()
    dt = float(parts[4])
    x, y, z, vx, vy, vz = (float(v) for v in parts[5:11])
    gmat_rows.append((dt, x, y, z, vx, vy, vz))

ksrop_rows = [
    (0.000000, -6555.799000, 0.000000, 0.000000, -0.000000000, -9.058150532, -4.816304071),
    (2495.634000, 5629.700996, -12716.382483, -6761.420503, 5.519900104, -1.920128049, -1.020950193),
    (8806.221000, 30000.700988, -12716.382483, -6761.420503, 2.564904149, 0.892216219, 0.474399779),
    (18931.761000, 42186.200985, -0.000000, -0.000000, 0.000000000, 1.407650199, 0.748460887),
    (29057.300000, 30000.700988, 12716.382483, 6761.420503, -2.564904149, 0.892216219, 0.474399779),
    (35367.887000, 5629.700996, 12716.382483, 6761.420503, -5.519900104, -1.920128049, -1.020950193),
    (37863.522000, -6555.799000, -0.000000, -0.000000, 0.000000000, -9.058150532, -4.816304071),
]

# also compare against pure two-body (no sun) to isolate the SIZE of the
# solar perturbation itself, for context
mu = 3.986004415e5
def twobody_energy(x,y,z,vx,vy,vz):
    r = (x*x+y*y+z*z)**0.5
    v2 = vx*vx+vy*vy+vz*vz
    return v2/2 - mu/r

print(f"{'dt(s)':>10} {'|dr|(km)':>12} {'|dv|(km/s)':>14}")
for (dt_k, xk, yk, zk, vxk, vyk, vzk), (dt_g, xg, yg, zg, vxg, vyg, vzg) in zip(ksrop_rows, gmat_rows):
    dr = ((xk-xg)**2 + (yk-yg)**2 + (zk-zg)**2) ** 0.5
    dv = ((vxk-vxg)**2 + (vyk-vyg)**2 + (vzk-vzg)**2) ** 0.5
    print(f"{dt_k:10.3f} {dr:12.4f} {dv:14.6e}   (gmat dt={dt_g:.3f})")
