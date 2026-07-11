"""
Build a shared 72x72 EGM2008 coefficient set from the full
input/EGM2008_to2190_TideFree file, usable by both:
  - KSROP driver_KS.F (geo_coeff reads raw n m Cnm Snm sigC sigS,
    zonal-only, whitespace/list-directed format, D exponents ok)
  - GMAT (needs a .cof file: header + RECOEF n m Cnm Snm, E exponents)
"""

SRC = r"C:\Users\hari2\OneDrive\Documents\GitHub\KSROP\input\EGM2008_to2190_TideFree"
RAW_OUT = r"C:\Users\hari2\OneDrive\Documents\GitHub\KSROP\scratch_gmat\EGM2008_72x72.dat"
COF_OUT = r"C:\Users\hari2\OneDrive\Documents\GitHub\KSROP\scratch_gmat\EGM2008_72x72.cof"

MU = "3.98600441500000E+14"
RE = "6.37813630000000E+06"
MAXDEG = 72

rows = []
with open(SRC) as f:
    for line in f:
        parts = line.split()
        n = int(parts[0])
        if n > MAXDEG:
            break
        rows.append(parts)

print(f"Collected {len(rows)} coefficient rows (degree 2..{MAXDEG})")

# --- 1) Raw truncated file, same format as source (KSROP-compatible) ---
with open(RAW_OUT, "w") as f:
    for n, m, cnm, snm, sc, ss in rows:
        f.write(f"{int(n):5d}{int(m):5d}   {cnm:>25}   {snm:>25}   {sc:>18}   {ss:>18}\n")

# --- 2) GMAT .cof file (header + RECOEF, fixed-width E21.14 fields
#     matching GMAT's own EGM96.cof layout exactly; m=0 rows carry
#     only Cnm (no Snm field), same as EGM96.cof) ---
def fval(s):
    x = float(s.replace("D", "E").replace("d", "e"))
    return format(x, ".14E").rjust(21)

with open(COF_OUT, "w") as f:
    f.write("COMMENT   5\n")
    f.write("C" * 80 + "\n")
    f.write("CCCCC  ------------------------------------------------------------------  CCCCC\n")
    f.write("CCCCC  EGM2008_72x72.cof : [72x72] truncated from EGM2008_to2190_TideFree   CCCCC\n")
    f.write("CCCCC  ------------------------------------------------------------------  CCCCC\n")
    f.write("C" * 80 + "\n")
    f.write(f"POTFIELD{MAXDEG}{MAXDEG:>4}  1 {MU} {RE} 1.00000000000000E+00\n")
    for n, m, cnm, snm, sc, ss in rows:
        n, m = int(n), int(m)
        if m == 0:
            f.write(f"RECOEF{n:5d}{m:3d}   {fval(cnm)}\n")
        else:
            f.write(f"RECOEF{n:5d}{m:3d}   {fval(cnm)}{fval(snm)}\n")

print(f"Wrote {RAW_OUT}")
print(f"Wrote {COF_OUT}")
