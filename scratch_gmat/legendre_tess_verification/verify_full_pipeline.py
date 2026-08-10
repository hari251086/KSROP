"""
Full end-to-end verification of the CORRECTED .wl-derived classical-
Legendre-in-KS-u geopotential force (V_pert, q(j), tau), fixing both
confirmed bugs:
  (1) general (n,m) associated Legendre recursion (0<m<n-1 branch)
  (2) longitude computed from (x,y), not (y,z)

Cross-validated against KSROP's own already-shipped, already-validated
Cunningham (1970) implementation (src/Cunningham.F, issue #30) across
nmax=2..6, multiple random full-triangle Cnm/Snm coefficient sets, and
multiple random u vectors -- mirroring exactly how Cunningham.F itself
was accepted (agreement with a trusted independent formulation to near
machine precision, issue #30's own acceptance gate).
"""
import numpy as np

rng = np.random.default_rng(42)


# ---------------- Cunningham reference (ground truth) ----------------

def cunningham_Vnm(nwork, x, y, z, r2):
    Vr = np.zeros((nwork + 1, nwork + 1))
    Vi = np.zeros((nwork + 1, nwork + 1))
    r = np.sqrt(r2)
    Vr[0, 0] = 1.0 / r
    for n in range(1, nwork + 1):
        pr, pim = Vr[n - 1, n - 1], Vi[n - 1, n - 1]
        Vr[n, n] = (2*n-1) * (x*pr - y*pim) / r2
        Vi[n, n] = (2*n-1) * (x*pim + y*pr) / r2
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
    return dxr, dxi, dyr, dyi, -cz*Vr[n+1, m], -cz*Vi[n+1, m]


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
            if Cnm[n, m] != 0.0 or Snm[n, m] != 0.0:
                Vtot += Rn * (Cnm[n, m]*Vr[n, m] + Snm[n, m]*Vi[n, m])
                dxr, dxi, dyr, dyi, dzr, dzi = cunningham_dVnm(n, m, Vr, Vi)
                gx += Rn * (Cnm[n, m]*dxr + Snm[n, m]*dxi)
                gy += Rn * (Cnm[n, m]*dyr + Snm[n, m]*dyi)
                gz += Rn * (Cnm[n, m]*dzr + Snm[n, m]*dzi)

    dVdu1 = 2*u1*gx + 2*u2*gy + 2*u3*gz
    dVdu2 = -2*u2*gx + 2*u1*gy + 2*u4*gz
    dVdu3 = -2*u3*gx - 2*u4*gy + 2*u1*gz
    dVdu4 = 2*u4*gx - 2*u3*gy + 2*u2*gz

    q = np.array([Vtot*u1 + (r_ks/2)*dVdu1, Vtot*u2 + (r_ks/2)*dVdu2,
                  Vtot*u3 + (r_ks/2)*dVdu3, Vtot*u4 + (r_ks/2)*dVdu4])
    tau = -2*r_ks*Vtot - (r_ks/2)*(u1*dVdu1 + u2*dVdu2 + u3*dVdu3 + u4*dVdu4)
    return q, tau


# ---------------- Corrected .wl-derived classical-Legendre path ----------------

def alfP_general(nmax, sp, cp):
    """Corrected general (n,m) unnormalized associated Legendre P(n,m)(sin phi)."""
    P = np.zeros((nmax + 1, nmax + 1))
    P[0, 0] = 1.0
    for n in range(1, nmax + 1):
        P[n, n] = (2*n - 1) * cp * P[n-1, n-1]
        for m in range(n - 1, -1, -1):
            pnm2 = P[n-2, m] if n - 2 >= m else 0.0
            P[n, m] = ((2*n-1)*sp*P[n-1, m] - (n+m-1)*pnm2) / (n - m)
    return P


def cos_sin_mlambda(nmax, coslam, sinlam):
    cosm = np.zeros(nmax + 2)  # need one extra for m=n+1 boundary access safety
    sinm = np.zeros(nmax + 2)
    cosm[0], sinm[0] = 1.0, 0.0
    if nmax >= 1:
        cosm[1], sinm[1] = coslam, sinlam
    for m in range(2, nmax + 2):
        cosm[m] = 2*coslam*cosm[m-1] - cosm[m-2]
        sinm[m] = 2*coslam*sinm[m-1] - sinm[m-2]
    return cosm, sinm


def wl_corrected_tess_force(u, R_body, nmax, Cnm, Snm):
    """CORRECTED classical-Legendre-in-KS-u general (n,m) force, matching
    Cunningham's exact mu-free q(j)/tau convention. Both confirmed bugs
    fixed: (1) general Legendre recursion, (2) longitude from (x,y)."""
    u1, u2, u3, u4 = u
    x = u1*u1 - u2*u2 - u3*u3 + u4*u4
    y = 2*(u1*u2 - u3*u4)
    z = 2*(u1*u3 + u2*u4)
    r = u1*u1 + u2*u2 + u3*u3 + u4*u4  # rKS

    sp = z / r
    cp = np.sqrt(x*x + y*y) / r
    tanphi = z / np.sqrt(x*x + y*y)
    # CORRECTED: longitude from (x,y), z is the polar axis
    denom_xy = np.sqrt(x*x + y*y)
    coslam, sinlam = x/denom_xy, y/denom_xy

    P = alfP_general(nmax + 1, sp, cp)  # +1 for the m+1 access in Vphi
    cosm, sinm = cos_sin_mlambda(nmax + 1, coslam, sinlam)

    Vtot = 0.0
    Vphi_sum = 0.0
    Vlam_sum = 0.0
    Vr_extra = 0.0  # sum n*Fnm piece
    for n in range(2, nmax + 1):
        Rn_r = (R_body / r) ** n
        for m in range(0, n + 1):
            c, s = Cnm[n, m], Snm[n, m]
            if c == 0.0 and s == 0.0:
                continue
            trig = c*cosm[m] + s*sinm[m]
            Fnm = Rn_r * P[n, m] * trig
            Vtot += Fnm
            Vr_extra += n * Fnm
            Pnp1 = P[n, m+1] if m+1 <= n else 0.0
            Vphi_sum += Rn_r * trig * (Pnp1 - m*tanphi*P[n, m])
            Vlam_sum += Rn_r * P[n, m] * m*(s*cosm[m] - c*sinm[m])

    v_pert = Vtot / r  # matches Cunningham's mu-free Vtot exactly
    Vr_block = -(1.0/r) * (v_pert + (1.0/r) * Vr_extra)
    Vphi_block = (1.0/r) * Vphi_sum
    Vlam_block = (1.0/r) * Vlam_sum

    dzdu = [2*u3, 2*u4, 2*u1, 2*u2]
    dydu = [2*u2, 2*u1, -2*u4, -2*u3]
    dxdu = [2*u1, -2*u2, -2*u3, 2*u4]
    drdu = [2*u1, 2*u2, 2*u3, 2*u4]

    dVdu = np.zeros(4)
    for k in range(4):
        dphidu = (dzdu[k] - sp*drdu[k]) / np.sqrt(x*x + y*y)
        dlamdu = (x*dydu[k] - y*dxdu[k]) / (x*x + y*y)  # CORRECTED: x,y not y,z
        dVdu[k] = Vr_block*drdu[k] + Vphi_block*dphidu + Vlam_block*dlamdu

    q = np.array([v_pert*u1 + (r/2)*dVdu[0], v_pert*u2 + (r/2)*dVdu[1],
                  v_pert*u3 + (r/2)*dVdu[2], v_pert*u4 + (r/2)*dVdu[3]])
    tau = -2*r*v_pert - (r/2)*(u1*dVdu[0] + u2*dVdu[1] + u3*dVdu[2] + u4*dVdu[3])
    return q, tau


# ---------------- Cross-validation sweep ----------------

print(f"{'nmax':>5} {'trial':>6} {'max|dq|':>12} {'|dtau|':>12}  verdict")
worst = 0.0
for nmax in [2, 3, 4, 5, 6]:
    for trial in range(5):
        u = rng.uniform(-1.2, 1.2, size=4)
        Cnm = np.zeros((nmax+1, nmax+1))
        Snm = np.zeros((nmax+1, nmax+1))
        for n in range(2, nmax+1):
            for m in range(0, n+1):
                Cnm[n, m] = rng.uniform(-1e-5, 1e-5)
                if m > 0:
                    Snm[n, m] = rng.uniform(-1e-5, 1e-5)

        q_ref, tau_ref = cunningham_tess_force(u, 6378.1363, nmax, Cnm, Snm)
        q_new, tau_new = wl_corrected_tess_force(u, 6378.1363, nmax, Cnm, Snm)

        dq = np.max(np.abs(q_ref - q_new)) / max(np.max(np.abs(q_ref)), 1e-30)
        dtau = abs(tau_ref - tau_new) / max(abs(tau_ref), 1e-30)
        worst = max(worst, dq, dtau)
        verdict = "OK" if max(dq, dtau) < 1e-9 else "*** MISMATCH ***"
        print(f"{nmax:>5} {trial:>6} {dq:>12.3e} {dtau:>12.3e}  {verdict}")

print(f"\nWorst relative error across all trials: {worst:.3e}")
print("CORRECTED PIPELINE VALIDATED" if worst < 1e-9 else "STILL MISMATCHED -- more bugs remain")
