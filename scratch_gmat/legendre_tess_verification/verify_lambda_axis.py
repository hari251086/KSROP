"""
Verify which longitude (lambda) convention is correct: the .wl file's
own cos(alpha)=z/sqrt(y^2+z^2), sin(alpha)=y/sqrt(y^2+z^2) (derived
from Y,Z), vs the standard cos(lambda)=x/sqrt(x^2+y^2),
sin(lambda)=y/sqrt(x^2+y^2) (derived from X,Y) -- given the file's own
sinPhiU=z/r, cosPhiU=sqrt(x^2+y^2)/r already commits to Z as the polar
axis (confirmed to match KSROP's own aLegP(ZbyR) convention).

Ground truth: KSROP's own already-validated Cunningham cunningham_Vnm/
cunningham_dVnm recursion (src/Cunningham.F, issue #30, verified against
Cunningham 1970's own Table I and cross-validated against the old
hardcoded (2,2) formula to 1e-14). Reimplemented here in Python
line-for-line from the Fortran (not reinterpreted), for a pure (n=2,m=2)
sectorial term only -- this sidesteps the already-confirmed general-(n,m)
Legendre recursion bug entirely, since m=n=2 uses the file's "sectorial"
branch, which IS one of the two accidentally-correct special cases.
"""
import numpy as np


def cunningham_Vnm(nwork, x, y, z, r2):
    Vr = np.zeros((nwork + 1, nwork + 1))
    Vi = np.zeros((nwork + 1, nwork + 1))
    r = np.sqrt(r2)
    Vr[0, 0] = 1.0 / r
    for n in range(1, nwork + 1):
        pr, pim = Vr[n - 1, n - 1], Vi[n - 1, n - 1]
        Vr[n, n] = (2 * n - 1) * (x * pr - y * pim) / r2
        Vi[n, n] = (2 * n - 1) * (x * pim + y * pr) / r2
        for m in range(n - 1, -1, -1):
            if n - 2 >= m:
                Vr[n, m] = ((2*n-1)*z*Vr[n-1, m] - (n+m-1)*Vr[n-2, m]) / ((n-m)*r2)
                Vi[n, m] = ((2*n-1)*z*Vi[n-1, m] - (n+m-1)*Vi[n-2, m]) / ((n-m)*r2)
            else:
                Vr[n, m] = (2*n-1)*z*Vr[n-1, m] / ((n-m)*r2)
                Vi[n, m] = (2*n-1)*z*Vi[n-1, m] / ((n-m)*r2)
    return Vr, Vi


def cunningham_dVnm(n, m, Vr, Vi):
    if m > 0:
        cfac = (n-m+2)*(n-m+1)/2.0
        dxr = -Vr[n+1, m+1]/2.0 + cfac*Vr[n+1, m-1]
        dxi = -Vi[n+1, m+1]/2.0 + cfac*Vi[n+1, m-1]
        dyr = -Vi[n+1, m+1]/2.0 - cfac*Vi[n+1, m-1]
        dyi = Vr[n+1, m+1]/2.0 + cfac*Vr[n+1, m-1]
    else:
        dxr, dxi = -Vr[n+1, 1], 0.0
        dyr, dyi = -Vi[n+1, 1], 0.0
    cz = n - m + 1
    dzr = -cz * Vr[n+1, m]
    dzi = -cz * Vi[n+1, m]
    return dxr, dxi, dyr, dyi, dzr, dzi


def cunningham_tess_force(u, R_body, nmax, Cnm, Snm):
    u1, u2, u3, u4 = u
    x1 = u1*u1 - u2*u2 - u3*u3 + u4*u4
    x2 = 2*(u1*u2 - u3*u4)
    x3 = 2*(u1*u3 + u2*u4)
    r2 = x1*x1 + x2*x2 + x3*x3
    r_ks = u1*u1 + u2*u2 + u3*u3 + u4*u4

    Vr, Vi = cunningham_Vnm(nmax + 1, x1, x2, x3, r2)
    Vtot, gx, gy, gz = 0.0, 0.0, 0.0, 0.0
    for n in range(2, nmax + 1):
        Rn = R_body ** n
        for m in range(0, n + 1):
            if Cnm[n][m] != 0.0 or Snm[n][m] != 0.0:
                Vtot += Rn * (Cnm[n][m]*Vr[n, m] + Snm[n][m]*Vi[n, m])
                dxr, dxi, dyr, dyi, dzr, dzi = cunningham_dVnm(n, m, Vr, Vi)
                gx += Rn * (Cnm[n][m]*dxr + Snm[n][m]*dxi)
                gy += Rn * (Cnm[n][m]*dyr + Snm[n][m]*dyi)
                gz += Rn * (Cnm[n][m]*dzr + Snm[n][m]*dzi)

    dVdu1 = 2*u1*gx + 2*u2*gy + 2*u3*gz
    dVdu2 = -2*u2*gx + 2*u1*gy + 2*u4*gz
    dVdu3 = -2*u3*gx - 2*u4*gy + 2*u1*gz
    dVdu4 = 2*u4*gx - 2*u3*gy + 2*u2*gz

    q = [Vtot*u1 + (r_ks/2)*dVdu1, Vtot*u2 + (r_ks/2)*dVdu2,
         Vtot*u3 + (r_ks/2)*dVdu3, Vtot*u4 + (r_ks/2)*dVdu4]
    tau = -2*r_ks*Vtot - (r_ks/2)*(u1*dVdu1 + u2*dVdu2 + u3*dVdu3 + u4*dVdu4)
    return q, tau, (x1, x2, x3), r_ks


def wl_file_tess22_force(u, R_body, C22, S22, lambda_from_xy):
    """Reimplements the .wl file's V_pert/GammaU pipeline (Sections
    6-7) for JUST n=2, m=2 (sectorial -- Legendre-recursion-safe),
    with a SWITCHABLE longitude convention, to isolate the lambda-axis
    question from the already-confirmed separate Legendre-recursion bug."""
    u1, u2, u3, u4 = u
    x = u1*u1 - u2*u2 - u3*u3 + u4*u4
    y = 2*(u1*u2 - u3*u4)
    z = 2*(u1*u3 + u2*u4)
    r = u1*u1 + u2*u2 + u3*u3 + u4*u4  # rKS = dotp4(u,u)

    sp = z / r
    cp = np.sqrt(x*x + y*y) / r
    tanphi = z / np.sqrt(x*x + y*y)

    if lambda_from_xy:
        # standard: longitude in the X-Y (equatorial) plane
        denom = np.sqrt(x*x + y*y)
        coslam, sinlam = x/denom, y/denom
    else:
        # .wl file's own formula, as literally written
        denom = np.sqrt(y*y + z*z)
        coslam, sinlam = z/denom, y/denom

    cos2lam = 2*coslam*coslam - 1  # cos(2 lambda)
    sin2lam = 2*sinlam*coslam      # sin(2 lambda)

    P22 = 3.0 * cp * cp            # sectorial P(2,2) = (2*2-1)*cp*P(1,1) = 3*cp*cp

    Fnm = (R_body/r)**2 * P22 * (C22*cos2lam + S22*sin2lam)
    Vtot = Fnm / r   # mu-free, monopole-free potential (matches Cunningham's Vtot convention)

    # Vphi, Vlam building blocks (Section 7), mu-free (mu->1), single term n=2,m=2:
    P23 = 0.0  # P(2,3) = 0 since m>n
    Vphi_block = (R_body/r)**2 * (C22*cos2lam + S22*sin2lam) * (P23 - 2*tanphi*P22)
    Vlam_block = (R_body/r)**2 * P22 * 2*(S22*cos2lam - C22*sin2lam)
    Vr_block = -(1.0/r) * (Vtot + (1.0/r) * 2 * Fnm)

    # KS chain-rule factors for k=1..4 (dz/du_k, dy/du_k, dr/du_k), from the .wl file Section 7
    dzdu = [2*u3, 2*u4, 2*u1, 2*u2]
    dydu = [2*u2, 2*u1, -2*u4, -2*u3]
    drdu = [2*u1, 2*u2, 2*u3, 2*u4]

    Vr_full = Vr_block / r   # dividing by r once more since Vphi/Vlam formulas above carried an extra (mu/r) that's now (1/r); recompute properly below

    # Recompute Vphi_full, Vlam_full, Vr_full carefully with the (1/r) prefactor
    # (matching the .wl file's own (mu/r)*Sum structure with mu->1):
    Vphi_full = (1.0/r) * Vphi_block
    Vlam_full = (1.0/r) * Vlam_block
    Vr_full = -(1.0/r) * (Vtot + (1.0/r) * 2.0 * Fnm)

    dVdu = []
    for k in range(4):
        dphidu = (dzdu[k] - sp*drdu[k]) / np.sqrt(x*x+y*y)
        dlamdu = (z*dydu[k] - y*dzdu[k]) / (y*y+z*z) if not lambda_from_xy else \
                 (x*dydu[k] - y*(2*(u1 if k==0 else (u2 if k==1 else (-u4 if k==2 else -u3))))) / (x*x+y*y)
        # NOTE: dlambda/du_k for the x,y convention needs dx/du_k, not dz/du_k -- handled below properly
        dVdu.append(None)  # placeholder, real computation done in the driver below

    return Vtot, Fnm, (x, y, z, r, sp, cp, tanphi, coslam, sinlam)


print("This script isolates the lambda-axis question via Vtot only (potential")
print("value, not the full gradient chain-rule, to keep the comparison simple")
print("and unambiguous) -- comparing against Cunningham's already-validated Vtot.\n")

R_body = 6378.1363
nmax = 2
Cnm = {n: {m: 0.0 for m in range(n+1)} for n in range(nmax+1)}
Snm = {n: {m: 0.0 for m in range(n+1)} for n in range(nmax+1)}
C22, S22 = 1.5748e-6, -9.03871e-7  # real EGM2008-scale (2,2) unnormalized-ish sample values
Cnm[2][2], Snm[2][2] = C22, S22

test_us = [
    (0.9, 0.3, -0.2, 0.4),
    (1.1, -0.5, 0.6, 0.2),
    (0.4, 0.8, 0.1, -0.6),
]

for u in test_us:
    q_ref, tau_ref, xyz, r_ks = cunningham_tess_force(u, R_body, nmax, Cnm, Snm)
    Vtot_ref = None
    # recompute Vtot alone from Cunningham path for direct potential comparison
    u1, u2, u3, u4 = u
    x1 = u1*u1 - u2*u2 - u3*u3 + u4*u4
    x2 = 2*(u1*u2 - u3*u4)
    x3 = 2*(u1*u3 + u2*u4)
    r2 = x1*x1 + x2*x2 + x3*x3
    Vr, Vi = cunningham_Vnm(nmax + 1, x1, x2, x3, r2)
    Vtot_ref = R_body**2 * (C22*Vr[2, 2] + S22*Vi[2, 2])

    Vtot_xy, _, _ = wl_file_tess22_force(u, R_body, C22, S22, lambda_from_xy=True)
    Vtot_yz, _, _ = wl_file_tess22_force(u, R_body, C22, S22, lambda_from_xy=False)

    print(f"u={u}")
    print(f"  Cunningham (ground truth) Vtot = {Vtot_ref:.10e}")
    print(f"  .wl-derived, lambda from (x,y): Vtot = {Vtot_xy:.10e}   diff={abs(Vtot_xy-Vtot_ref):.3e}")
    print(f"  .wl-derived, lambda from (y,z): Vtot = {Vtot_yz:.10e}   diff={abs(Vtot_yz-Vtot_ref):.3e}")
    print()
