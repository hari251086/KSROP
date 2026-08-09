c     ============================================================
      subroutine alfP_general(nmax, sp, cp, P)
c     KSROP: general (n,m) unnormalized associated Legendre function
c     P_n^m(sin phi), via the standard three-term recurrence relation
c     (e.g. Montenbruck & Gill 2000 eq. 3.24, or the equivalent
c     "belt" recursion used throughout gravity-field synthesis):
c       P(n,n) = (2n-1) cos(phi) P(n-1,n-1)                (sectorial)
c       P(n,m) = [(2n-1) sin(phi) P(n-1,m) - (n+m-1) P(n-2,m)]
c                / (n-m),                          0 <= m < n
c       P(n,m) = 0,                                 m > n
c     (m=0/zonal is just the m=0 instance of the SAME general
c     formula -- no separate special case needed, unlike an earlier
c     derivation this replaces, see below).
c
c     Derivation source: EarthGravityPotential_KS (4).wl
c     (E:\Research\1. R&D\Mathematica) -- a from-scratch KS-
c     regularized re-derivation of the full (zonal+tesseral)
c     geopotential and its gradient, requested for evaluation as an
c     independent alternative to this repo's existing general (n,m)
c     tesseral path (src/Cunningham.F, issue #30).
c
c     BUG FOUND AND FIXED (2026-08-09): the source .wl file's own
c     "0 < m < n" recursion branch was
c       Pnm[n,m] = (2n-1) cos(phi) Pnm[n-1,m-1]
c     -- structurally different from, and NOT equivalent to, the
c     correct formula above. Verified numerically wrong against
c     scipy.special.lpmv ground truth for every interior (n,m) with
c     0 < m < n-1 (only accidentally correct at the two edges m=n
c     [sectorial] and m=n-1 [sub-sectorial], where it happens to
c     coincide with the correct recursion); error grows from O(1) at
c     (n,m)=(3,1) to O(1e5) by (n,m)=(8,1). This subroutine implements
c     the CORRECTED recursion above, re-verified to machine precision
c     (worst case 5.6e-15) against the already-validated Cunningham
c     path (cunningham_Vnm/tess_general_force) across nmax=2..6,
c     random full-triangle Cnm/Snm, random u -- see
c     test/test_legendre_tess.F.
c
c     nmax must be >= (largest degree n actually needed)+1, matching
c     aLegP's own "ngeo_deg+1" and cunningham_Vnm's "nwork>=n+1"
c     conventions, since tess_legendre_force's phi-derivative
c     identity needs P(n,m+1).
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension P(0:nmax,0:nmax)

      P(0,0) = 1.d0
      do n = 1, nmax
         P(n,n) = dble(2*n-1)*cp*P(n-1,n-1)
         do m = n-1, 0, -1
            if (n-2 .ge. m) then
               pnm2 = P(n-2,m)
            else
               pnm2 = 0.d0
            end if
            P(n,m) = (dble(2*n-1)*sp*P(n-1,m) - dble(n+m-1)*pnm2)
     &               / dble(n-m)
         end do
      end do

      return
      end

c     ============================================================
      subroutine cos_sin_mlambda(nmax, coslam, sinlam, cosm, sinm)
c     cos(m*lambda), sin(m*lambda) for m=0..nmax via the standard
c     Chebyshev multiple-angle recursion (angle-addition identity):
c       cos(m theta) = 2 cos(theta) cos((m-1)theta) - cos((m-2)theta)
c       sin(m theta) = 2 cos(theta) sin((m-1)theta) - sin((m-2)theta)
c     Purely algebraic in cos(lambda)/sin(lambda) -- no ArcTan/ArcSin,
c     no branch-cut concerns. Standard identity, unchanged from
c     EarthGravityPotential_KS (4).wl's own derivation (this piece
c     had no bug -- only the associated-Legendre recursion above and
c     the longitude convention in tess_legendre_force below did).
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension cosm(0:nmax), sinm(0:nmax)

      cosm(0) = 1.d0
      sinm(0) = 0.d0
      if (nmax .ge. 1) then
         cosm(1) = coslam
         sinm(1) = sinlam
      end if
      do m = 2, nmax
         cosm(m) = 2.d0*coslam*cosm(m-1) - cosm(m-2)
         sinm(m) = 2.d0*coslam*sinm(m-1) - sinm(m-2)
      end do

      return
      end

c     ============================================================
      subroutine tess_legendre_force(u, R_body, theta_dot, nmax,
     &                                Cnm, Snm,
     &                                q_leg, tau_leg, dvdt_leg,
     &                                v_pert_out)
c     KSROP: general (n,m) tesseral/mascon KS-element force (q_leg)
c     and time-element (tau_leg) contribution -- an INDEPENDENT
c     alternative derivation to tess_general_force (src/Cunningham.F,
c     issue #30), via the CLASSICAL spherical-coordinate (latitude/
c     longitude) associated-Legendre expansion substituted directly
c     into KS u-variables, rather than Cunningham's Cartesian solid-
c     harmonic recursion. Same q(j)/tau output convention as
c     tess_general_force (amue-free, mu-free, monopole-free) -- drop-
c     in comparable/interchangeable with it.
c
c     Derivation source: EarthGravityPotential_KS (4).wl
c     (E:\Research\1. R&D\Mathematica). TWO bugs found in the source
c     file and fixed here (2026-08-09):
c       1. General associated Legendre recursion -- see
c          alfP_general's header. Wrong for every interior (n,m),
c          0<m<n-1. Independently corroborated (2026-08-09, same day)
c          against a second reference ("Analytical and Numerical
c          Predictions for Near-Earth's Satellite Orbits with KS
c          Uniform Regular Canonical Equations", eq 3.2.2 -- the same
c          section already cited in driver_KS.F's header for the
c          geopotential physics), whose "otherwise" branch states the
c          identical incorrect formula -- exact-rational-arithmetic
c          verified (3-4-5 triple) against the Rodrigues formula, same
c          failure pattern. The bug is inherited from this literature
c          lineage, not a one-off transcription slip.
c       2. Longitude (lambda) was computed from Y,Z
c          (cos(alpha)=z/sqrt(y^2+z^2)) -- inconsistent with the same
c          file's own latitude convention (sin(phi)=z/r, i.e. Z is
c          the polar axis, confirmed to match this repo's own
c          aLegP(ZbyR) usage in driver_KS.F). Longitude must instead
c          come from X,Y: cos(lambda)=x/sqrt(x^2+y^2),
c          sin(lambda)=y/sqrt(x^2+y^2). Verified numerically: the
c          Y,Z formula disagreed with the already-validated
c          Cunningham potential by 50-100+ units (even wrong sign) on
c          a pure (2,2) sectorial test case (which sidesteps bug #1
c          entirely, isolating this as a second, independent defect);
c          the X,Y formula matches to 1e-14.
c
c     With both fixed, this subroutine's full (n,m) sum, gradient
c     chain rule (via r,phi,lambda as intermediate spherical
c     coordinates -- the standard geopotential-gradient decomposition,
c     Montenbruck & Gill sec. 3.2 or Vallado) and mu-free q(j)/tau
c     packaging are verified to reproduce tess_general_force to
c     machine precision (worst case 5.6e-15) across nmax=2..6, random
c     full-triangle Cnm/Snm, random u -- see
c     test/test_legendre_tess.F, the acceptance gate for trusting
c     this as a genuine independent cross-check (same discipline as
c     issue #30's own tess22_force acceptance gate).
c
c     Cnm/Snm are the SAME already-converted, correctly-signed,
c     fully-normalized-to-unnormalized coefficient arrays
c     geo_coeff_tess_general produces -- no separate loader needed;
c     this can be fed the exact arrays already loaded for
c     tess_general_force.
c
c     dvdt_leg (2026-08-09): the mu-free rate of change of this
c     potential w.r.t. PHYSICAL time t, holding u (hence r,phi,lambda)
c     fixed -- needed because the tesseral coefficients Cnm/Snm are
c     defined body-fixed but driver_KS.F evaluates this potential in
c     the inertial frame via a fresh per-step rotation by theta(t)
c     (rotate_tess_coeffs), making the potential explicitly time-
c     dependent: U(u,t), not U(u). Per Sellamuthu 2018 PhD thesis eq.
c     2.54, the KS energy-element rate w* has a
c     "-(r/(8w^2))*dU/dt" term that vanishes only when U does not
c     explicitly depend on t -- previously omitted here since every
c     OTHER conservative force in this codebase (zonal, third-body)
c     was treated as static; the general tesseral term breaks that
c     assumption.
c
c     Derivation: writing Cnm_rot(t)=Cnm*cos(m*theta)-Snm*sin(m*theta),
c     Snm_rot(t)=Cnm*sin(m*theta)+Snm*cos(m*theta) (the exact
c     rotate_tess_coeffs convention) is algebraically identical to
c     evaluating the potential at body-fixed longitude
c     psi=lambda-theta(t) with the RAW (un-rotated) Cnm/Snm. Since
c     lambda is a pure function of u (atan2(y,x)), holding u fixed
c     holds lambda fixed, so d/dt enters solely through theta(t):
c       d/dt[Cnm*cos(m*psi)+Snm*sin(m*psi)]
c         = -m*theta_dot*(Snm_rot*cos(m*lambda)-Cnm_rot*sin(m*lambda))
c     -- verified in sympy to zero symbolic residual, and independently
c     via finite difference (central diff at h=1e-7 in theta, u held
c     exactly fixed, matched to 1.26e-9 relative -- consistent with
c     pure O(h^2) truncation, i.e. effectively exact). The RHS is
c     exactly m times the per-term summand already accumulated in
c     vlam_sum below -- no new recursion needed:
c       d(v_pert)/dt |_u = -theta_dot * Vlam_block
c     theta_dot is Earth's rotation rate (WE_rot at the call site,
c     confirmed to match gmst_deg's own dominant rate,
c     360.98564736629 deg/day = 7.29211e-5 rad/s).
c
c     Zonal (m=0) terms contribute exactly zero to this identity (the
c     m factor vanishes) -- physically correct, since an axisymmetric
c     term cannot exchange energy with a rotating frame.
c
c     v_pert_out: the mu-free potential value itself (same v_pert
c     already used internally for q_leg/tau_leg) -- exposed so a
c     caller can build a self-consistent initial energy element w
c     (driver_KS.F's V_pot at t=t0), which needs the tesseral
c     potential's VALUE, not just its u-gradient or time-rate.
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension u(4), q_leg(4)
      dimension Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax)
      dimension P(0:nmax+1,0:nmax+1)
      dimension cosm(0:nmax), sinm(0:nmax)
      dimension dzdu(4), dydu(4), dxdu(4), drdu(4), dVdu(4)

      u1 = u(1)
      u2 = u(2)
      u3 = u(3)
      u4 = u(4)

      x = u1*u1 - u2*u2 - u3*u3 + u4*u4
      y = 2.d0*(u1*u2 - u3*u4)
      z = 2.d0*(u1*u3 + u2*u4)
      r_ks = u1*u1 + u2*u2 + u3*u3 + u4*u4

      sxy = dsqrt(x*x + y*y)
      sp = z / r_ks
      cp = sxy / r_ks
      tanphi = z / sxy
c     Longitude from X,Y (Z is the polar axis) -- see bug #2 above.
      coslam = x / sxy
      sinlam = y / sxy

      call alfP_general(nmax+1, sp, cp, P)
      call cos_sin_mlambda(nmax, coslam, sinlam, cosm, sinm)

      vtot = 0.d0
      vr_extra = 0.d0
      vphi_sum = 0.d0
      vlam_sum = 0.d0

      do n = 2, nmax
         Rn_r = (R_body/r_ks)**n
         do m = 0, n
            cc = Cnm(n,m)
            ss = Snm(n,m)
            if (cc .ne. 0.d0 .or. ss .ne. 0.d0) then
               trig = cc*cosm(m) + ss*sinm(m)
               Fnm = Rn_r*P(n,m)*trig
               vtot = vtot + Fnm
               vr_extra = vr_extra + dble(n)*Fnm
               if (m+1 .le. n) then
                  Pnp1 = P(n,m+1)
               else
                  Pnp1 = 0.d0
               end if
               vphi_sum = vphi_sum
     &            + Rn_r*trig*(Pnp1 - dble(m)*tanphi*P(n,m))
               vlam_sum = vlam_sum
     &            + Rn_r*P(n,m)*dble(m)*(ss*cosm(m) - cc*sinm(m))
            end if
         end do
      end do

c     mu-free, monopole-free potential (matches Cunningham's Vtot
c     convention exactly -- see header note).
      v_pert = vtot / r_ks
      Vr_block = -(v_pert + vr_extra/r_ks) / r_ks
      Vphi_block = vphi_sum / r_ks
      Vlam_block = vlam_sum / r_ks

      dzdu(1) =  2.d0*u3
      dzdu(2) =  2.d0*u4
      dzdu(3) =  2.d0*u1
      dzdu(4) =  2.d0*u2

      dydu(1) =  2.d0*u2
      dydu(2) =  2.d0*u1
      dydu(3) = -2.d0*u4
      dydu(4) = -2.d0*u3

      dxdu(1) =  2.d0*u1
      dxdu(2) = -2.d0*u2
      dxdu(3) = -2.d0*u3
      dxdu(4) =  2.d0*u4

      drdu(1) =  2.d0*u1
      drdu(2) =  2.d0*u2
      drdu(3) =  2.d0*u3
      drdu(4) =  2.d0*u4

      do k = 1, 4
         dphidu = (dzdu(k) - sp*drdu(k)) / sxy
         dlamdu = (x*dydu(k) - y*dxdu(k)) / (sxy*sxy)
         dVdu(k) = Vr_block*drdu(k) + Vphi_block*dphidu
     &           + Vlam_block*dlamdu
      end do

      q_leg(1) = v_pert*u1 + (r_ks/2.d0)*dVdu(1)
      q_leg(2) = v_pert*u2 + (r_ks/2.d0)*dVdu(2)
      q_leg(3) = v_pert*u3 + (r_ks/2.d0)*dVdu(3)
      q_leg(4) = v_pert*u4 + (r_ks/2.d0)*dVdu(4)

      tau_leg = -2.d0*r_ks*v_pert - (r_ks/2.d0)*
     &          (u1*dVdu(1) + u2*dVdu(2) + u3*dVdu(3) + u4*dVdu(4))

c     mu-free d(v_pert)/dt at fixed u -- see header derivation.
      dvdt_leg = -theta_dot * Vlam_block
      v_pert_out = v_pert

      return
      end
