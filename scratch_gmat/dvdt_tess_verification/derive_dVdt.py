"""
PLAN-MODE RESEARCH ONLY (read-only, no repo files touched).

Symbolic derivation + verification of dV/dt for KSROP's general (n,m)
tesseral geopotential, needed for the wdot (energy-element rate)
correction requested by the user.

Physical setup (confirmed against actual driver_KS.F / LegendreTess.F
code, not assumed):
  - The KS position u (hence r, phi, INERTIAL longitude lambda) is
    held FIXED while asking d/dt.
  - The coefficients Cnm/Snm are defined in the body-fixed (rotating)
    frame; driver_KS.F evaluates the potential in the INERTIAL frame
    by pre-rotating them every step (rotate_tess_coeffs, called fresh
    each step at driver_KS.F:786-788 using theta(t)=gmst_deg(djulian1)):
        Cnm_rot(t) = Cnm*cos(m*theta(t)) - Snm*sin(m*theta(t))
        Snm_rot(t) = Cnm*sin(m*theta(t)) + Snm*cos(m*theta(t))
    then V is evaluated with Cnm_rot(t)/Snm_rot(t) at the INERTIAL
    lambda. This is algebraically identical to evaluating
    Cnm*cos(m*psi)+Snm*sin(m*psi) with psi = lambda - theta(t), the
    TRUE body-fixed longitude.
  - theta_dot = d(theta)/dt: confirmed gmst_deg's own dominant
    (linear) term, 360.98564736629 deg/day, converts to
    7.29211e-5 rad/s -- matching WE_rot (7.2921150d-5 rad/s) already
    used elsewhere in driver_KS.F (drag co-rotation). So WE_rot is the
    right, already-available theta_dot -- no new subroutine needed.

Goal: derive d/dt[ Cnm*cos(m*psi) + Snm*sin(m*psi) ] holding u fixed
(hence lambda fixed, psi=lambda-theta(t) time-varying only through
theta), and re-express the result in terms of the ALREADY-COMPUTED
LegendreTess.F quantities (Cnm_rot, Snm_rot -- i.e. the Cnm/Snm
ARGUMENTS tess_legendre_force actually receives -- and cos(m*lambda),
sin(m*lambda), i.e. cosm(m)/sinm(m) in the Fortran) to confirm whether
Vlam_sum (already computed inside tess_legendre_force) is directly
reusable, or a new expression is needed.
"""
import sympy as sp

m, theta, lam, t = sp.symbols('m theta lam t', real=True)
Cnm, Snm = sp.symbols('Cnm Snm', real=True)
theta_t = sp.Function('theta')(t)

psi = lam - theta_t  # true body-fixed longitude, lambda fixed (u fixed), theta varies

term = Cnm*sp.cos(m*psi) + Snm*sp.sin(m*psi)

dterm_dt = sp.diff(term, t)
dterm_dt = sp.simplify(dterm_dt)
print("d/dt[Cnm*cos(m*psi)+Snm*sin(m*psi)] wrt t, u (hence lambda) fixed:")
print(" ", dterm_dt)

# Express theta_dot = d(theta)/dt as a plain symbol for readability
theta_dot = sp.symbols('theta_dot', real=True)
dterm_dt_named = dterm_dt.subs(sp.Derivative(theta_t, t), theta_dot)
print("\nWith theta_dot substituted in:")
print(" ", sp.simplify(dterm_dt_named))

# Now: is this equal to -theta_dot * (Snm_rot*cos(m*lam) - Cnm_rot*sin(m*lam))
# where Cnm_rot/Snm_rot are defined via the ACTUAL rotate_tess_coeffs formula?
Cnm_rot = Cnm*sp.cos(m*theta) - Snm*sp.sin(m*theta)
Snm_rot = Cnm*sp.sin(m*theta) + Snm*sp.cos(m*theta)

candidate = -theta_dot * (Snm_rot*sp.cos(m*lam) - Cnm_rot*sp.sin(m*lam))
candidate = sp.expand_trig(sp.expand(candidate))

# dterm_dt_named is in terms of psi = lam - theta; rewrite theta_t -> theta (plain symbol)
# so both sides are directly comparable as functions of (lam, theta, m, Cnm, Snm, theta_dot)
lhs = dterm_dt_named.subs(theta_t, theta)
lhs = sp.expand_trig(sp.expand(lhs.rewrite(sp.cos)))

diff_check = sp.simplify(sp.expand_trig(sp.trigsimp(lhs - candidate)))
print("\nLHS (true d/dt, psi->lambda-theta explicit):")
print(" ", sp.simplify(lhs))
print("\nCandidate (-theta_dot * (Snm_rot*cos(m*lam) - Cnm_rot*sin(m*lam))):")
print(" ", sp.simplify(candidate))
print("\nDifference (should be 0 if candidate is correct):")
print(" ", diff_check)

# Numeric cross-check at random values, extra safety beyond symbolic simplify
import random
random.seed(0)
subs_vals = {m: 3, Cnm: 0.0123, Snm: -0.00456, theta: 0.7, lam: 1.9, theta_dot: 1.0}
lhs_num = float(lhs.subs(subs_vals))
cand_num = float(candidate.subs(subs_vals))
print(f"\nNumeric check: lhs={lhs_num:.12f}  candidate={cand_num:.12f}  diff={lhs_num-cand_num:.3e}")
