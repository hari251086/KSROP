"""
Numeric verification (no sympy needed) of the Pnm recursion from
EarthGravityPotential_KS (4).wl, at several generic phi values.
Polynomial identities that hold at more points than their degree
are identically true, so agreement at 5 generic phi values up to
degree ~8 is fully decisive.

Ground truth: scipy.special.lpmv computes the associated Legendre
function P_n^m(x) INCLUDING the Condon-Shortley phase (-1)^m. Geodesy
convention (and this .wl file, matching Rodrigues without the extra
sign) omits that phase, so ground_truth(n,m,x) = (-1)^m * lpmv(m,n,x).
"""
import numpy as np
from scipy.special import lpmv

memo_file = {}
memo_fix = {}


def file_Pnm(n, m, sp_, cp_):
    key = (n, m)
    if key in memo_file:
        f = memo_file[key]
        return f(sp_, cp_)
    if m > n:
        val = 0.0
    elif n == 0 and m == 0:
        val = 1.0
    elif n == 1 and m == 0:
        val = sp_
    elif n == 1 and m == 1:
        val = cp_
    elif m == n:
        val = (2*n - 1) * cp_ * file_Pnm(n-1, n-1, sp_, cp_)
    elif m == 0:
        val = (2*n-1)/n * sp_ * file_Pnm(n-1, 0, sp_, cp_) - (n-1)/n * file_Pnm(n-2, 0, sp_, cp_)
    else:
        val = (2*n - 1) * cp_ * file_Pnm(n-1, m-1, sp_, cp_)   # suspect branch
    return val


def fixed_Pnm(n, m, sp_, cp_):
    key = (n, m)
    if m > n:
        return 0.0
    if n == 0 and m == 0:
        return 1.0
    if n == 1 and m == 0:
        return sp_
    if n == 1 and m == 1:
        return cp_
    if m == n:
        return (2*n - 1) * cp_ * fixed_Pnm(n-1, n-1, sp_, cp_)
    # general 3-term recursion, SAME m at n-1, n-2 (also correctly covers m=0)
    return ((2*n-1) * sp_ * fixed_Pnm(n-1, m, sp_, cp_) - (n+m-1) * fixed_Pnm(n-2, m, sp_, cp_)) / (n - m)


def ground_truth(n, m, x):
    # scipy lpmv includes Condon-Shortley (-1)^m; geodesy/Rodrigues-without-CS
    # convention (matches this .wl file's P(1,1)=cos(phi), positive) needs
    # that phase removed.
    return ((-1) ** m) * lpmv(m, n, x)


print("Testing at 5 generic phi values, n=0..8, m=0..n")
print(f"{'n':>3} {'m':>3} {'max|file-truth|':>16} {'max|fixed-truth|':>17}  verdict")

phis = np.deg2rad([13.7, 41.3, 67.9, -22.4, 55.1])
sps = np.sin(phis)
cps = np.cos(phis)

any_file_mismatch = False
for n in range(0, 9):
    for m in range(0, n + 1):
        file_vals = np.array([file_Pnm(n, m, sp_, cp_) for sp_, cp_ in zip(sps, cps)])
        fixed_vals = np.array([fixed_Pnm(n, m, sp_, cp_) for sp_, cp_ in zip(sps, cps)])
        truth_vals = np.array([ground_truth(n, m, s_) for s_ in sps])
        err_file = np.max(np.abs(file_vals - truth_vals))
        err_fixed = np.max(np.abs(fixed_vals - truth_vals))
        verdict = "OK" if err_file < 1e-9 else "*** MISMATCH ***"
        if err_file >= 1e-9:
            any_file_mismatch = True
        print(f"{n:>3} {m:>3} {err_file:>16.3e} {err_fixed:>17.3e}  {verdict}")

print(f"\nAny file-recursion mismatch found: {any_file_mismatch}")
