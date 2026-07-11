import re

# --- GMAT report: last row of each block (block = between header lines) ---
gmat_path = r"E:\Softwares\gmat-win-R2026a\output\crosscheck_gmat.txt"
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
    last = b[-1]
    parts = last.split()
    # UTCGregorian is 4 tokens: "20 Sep 2016 00:00:00.000"
    dt = float(parts[4])
    x, y, z, vx, vy, vz = (float(v) for v in parts[5:11])
    gmat_rows.append((dt, x, y, z, vx, vy, vz))

# --- KSROP checkpoints (from earlier extraction) ---
ksrop_rows = [
    (0.000000, -0.000000, -5888.972700, -3400.000000, 9.500000000, -0.000000000, -0.000000000),
    (884.673000, 7309.481520, -3445.914161, -1989.499484, 6.279801693, 4.693299104, 2.709677522),
    (2113.846000, 11826.989539, 2950.096131, 1703.238809, 1.621939656, 5.134875477, 2.964621762),
    (3900.431000, 11826.989539, 10855.999637, 6267.714361, -1.158313030, 3.667086597, 2.117193450),
    (6244.429000, 7309.481520, 17252.009928, 9960.452654, -2.462936299, 1.840710470, 1.062734694),
    (8932.927000, 0.000000, 19695.068468, 11370.953170, -2.840571016, 0.000000000, 0.000000000),
    (11621.424000, -7309.481520, 17252.009928, 9960.452654, -2.462936299, -1.840710470, -1.062734694),
    (13965.422000, -11826.989539, 10855.999637, 6267.714361, -1.158313030, -3.667086597, -2.117193450),
    (15752.007000, -11826.989539, 2950.096131, 1703.238809, 1.621939656, -5.134875477, -2.964621762),
    (16981.180000, -7309.481520, -3445.914161, -1989.499484, 6.279801693, -4.693299104, -2.709677522),
    (17865.853000, 0.000000, -5888.972700, -3400.000000, 9.500000000, 0.000000000, 0.000000000),
]

print(f"{'dt(s)':>10} {'|dr|(km)':>14} {'|dv|(km/s)':>14}")
for (dt_k, xk, yk, zk, vxk, vyk, vzk), (dt_g, xg, yg, zg, vxg, vyg, vzg) in zip(ksrop_rows, gmat_rows):
    dr = ((xk-xg)**2 + (yk-yg)**2 + (zk-zg)**2) ** 0.5
    dv = ((vxk-vxg)**2 + (vyk-vyg)**2 + (vzk-vzg)**2) ** 0.5
    print(f"{dt_k:10.3f} {dr:14.6e} {dv:14.6e}   (gmat dt={dt_g:.3f})")
