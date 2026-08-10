c     ============================================================
      subroutine G_general(nmax, sp, G, Gd)
c     De-cos(phi)^m associated Legendre recursion: G(n,m) =
c     P_n^m(sinphi)/cos^m(phi), computed directly as a function of
c     sp=sin(phi)=z/r alone (no cos(phi) dependence at all) --
c     eliminates the polar (cos(phi)->0) singularity that dividing
c     the standard P_n^m by an explicit cos(phi)^m factor would
c     otherwise reintroduce. See ALGORITHM.md / branch derivation
c     notes for the full symbolic + numeric verification.
c       G(n,n) = (2n-1) G(n-1,n-1)                         (sectorial)
c       G(n,m) = [(2n-1) sp G(n-1,m) - (n+m-1) G(n-2,m)]
c                / (n-m),                          0 <= m < n
c     G(n,n) grows like (2n-1)!! -- large but still comfortably
c     within double-precision range even at n=72 (~1e124, vs. the
c     ~1.8e308 overflow threshold) -- unlike the raw x,y-power terms
c     this used to be combined with (see XY_general/tess_legendre_
c     force below for how that overflow risk is avoided).
c
c     Gd(n,m) = dG(n,m)/d(sp), obtained by differentiating this SAME
c     recursion directly term-by-term (the standard trick for a
c     non-singular Legendre-type derivative recursion -- no (1-sp^2)
c     or cos(phi)-type denominator anywhere). Verified symbolically
c     (sympy) to exactly match d/d(sp)[G(n,m)] computed by direct
c     differentiation, nmax=6, every (n,m).
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension G(0:nmax,0:nmax), Gd(0:nmax,0:nmax)

      G(0,0) = 1.d0
      Gd(0,0) = 0.d0
      do n = 1, nmax
         G(n,n) = dble(2*n-1)*G(n-1,n-1)
         Gd(n,n) = dble(2*n-1)*Gd(n-1,n-1)
         do m = n-1, 0, -1
            if (n-2 .ge. m) then
               gnm2 = G(n-2,m)
               gdnm2 = Gd(n-2,m)
            else
               gnm2 = 0.d0
               gdnm2 = 0.d0
            end if
            G(n,m) = (dble(2*n-1)*sp*G(n-1,m) - dble(n+m-1)*gnm2)
     &               / dble(n-m)
            Gd(n,m) = (dble(2*n-1)*(G(n-1,m)+sp*Gd(n-1,m))
     &                 - dble(n+m-1)*gdnm2) / dble(n-m)
         end do
      end do

      return
      end

c     ============================================================
      subroutine XY_general(nmax, xh, yh, Xt, Yt, Xx, Xy, Yx, Yy)
c     Unnormalized cos(m*lambda), sin(m*lambda) as pure polynomials
c     in xh=x/r, yh=y/r (the L2-normalized direction cosines, NOT raw
c     x,y) via the De Moivre / Chebyshev-type recursion applied to
c     (xh+i*yh)^m -- WITHOUT ever forming lambda=atan2(yh,xh) or
c     dividing by sqrt(xh^2+yh^2) the way the classical
c     cos(m*lambda)=Re[(x+iy)^m]/sxy^m evaluation would.
c
c     Passing the NORMALIZED xh,yh=x/r,y/r (2026-08-09, fixed same
c     day as the initial reformulation) rather than raw x,y is
c     essential for numerical stability at realistic degree/altitude
c     combinations: since xh^2+yh^2=cos^2(phi)<=1, every X_m,Y_m
c     value stays bounded by 1 in magnitude, however high m goes --
c     whereas the raw-x,y version grows like r^m (e.g. ~1e300+ at
c     GEO altitude with m=72), risking floating-point overflow before
c     the physically-correct cancellation with R_body^n/r^(n+m+1)
c     ever happens. See tess_legendre_force's own header for the
c     full Rr_n=(R_body/r)^n/r pairing that keeps the radial part
c     equally safe.
c       X_m = 2*xh*X_{m-1} - (xh^2+yh^2)*X_{m-2},   X_0=1, X_1=xh
c       Y_m = 2*xh*Y_{m-1} - (xh^2+yh^2)*Y_{m-2},   Y_0=0, Y_1=yh
c     Verified symbolically (sympy) to equal Re[(x+iy)^m]/r^m,
c     Im[(x+iy)^m]/r^m exactly for m=0..7 (homogeneity of the
c     original raw-x,y polynomial makes this identical to evaluating
c     that same polynomial at x/r,y/r).
c     Xx,Xy,Yx,Yy are the partial derivatives d/d(xh), d/d(yh) of
c     X_m, Y_m (obtained the same term-by-term differentiation way,
c     verified against direct differentiation, m=0..7).
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension Xt(0:nmax), Yt(0:nmax)
      dimension Xx(0:nmax), Xy(0:nmax), Yx(0:nmax), Yy(0:nmax)

      Xt(0) = 1.d0
      Yt(0) = 0.d0
      Xx(0) = 0.d0
      Xy(0) = 0.d0
      Yx(0) = 0.d0
      Yy(0) = 0.d0
      if (nmax .ge. 1) then
         Xt(1) = xh
         Yt(1) = yh
         Xx(1) = 1.d0
         Xy(1) = 0.d0
         Yx(1) = 0.d0
         Yy(1) = 1.d0
      end if

      sxy2 = xh*xh + yh*yh
      do m = 2, nmax
         Xt(m) = 2.d0*xh*Xt(m-1) - sxy2*Xt(m-2)
         Yt(m) = 2.d0*xh*Yt(m-1) - sxy2*Yt(m-2)
         Xx(m) = 2.d0*Xt(m-1) + 2.d0*xh*Xx(m-1)
     &           - 2.d0*xh*Xt(m-2) - sxy2*Xx(m-2)
         Xy(m) = 2.d0*xh*Xy(m-1) - 2.d0*yh*Xt(m-2) - sxy2*Xy(m-2)
         Yx(m) = 2.d0*Yt(m-1) + 2.d0*xh*Yx(m-1)
     &           - 2.d0*xh*Yt(m-2) - sxy2*Yx(m-2)
         Yy(m) = 2.d0*xh*Yy(m-1) - 2.d0*yh*Yt(m-2) - sxy2*Yy(m-2)
      end do

      return
      end

c     ============================================================
      subroutine tess_legendre_force(u, R_body, theta_dot, Rtilt,
     &                                nmax, Cnm, Snm,
     &                                q_leg, tau_leg, dvdt_leg,
     &                                v_pert_out)
c     KSROP: general (n,m) tesseral/mascon KS-element force (q_leg),
c     time-element (tau_leg), and energy-element rate (dvdt_leg)
c     contribution -- singularity-free reformulation (2026-08-09,
c     feature/tesseral-energy-time-dependence branch), replacing the
c     original spherical-coordinate (phi,lambda) chain-rule version
c     of this subroutine (which divided by sxy=sqrt(x^2+y^2) in its
c     u-derivatives of phi/lambda, a real polar-latitude singularity
c     already documented in ALGORITHM.md's Known Limitations).
c
c     Motivation: the user supplied an "exemplar" closed-form zonal
c     KS-element force (Sellamuthu 2018 PhD thesis, Ch. 3, eq. 3.4,
c     hand-derived for J2-J4) whose denominator is a pure power of r
c     -- no z-dependent (cos(phi)-type) denominator at all -- because
c     it differentiates the zonal potential directly through z(u),
c     r(u), never through an intermediate phi angle. This subroutine
c     generalizes that same principle to the FULL general (n,m) case
c     (not just zonal, and not hardcoded per degree): the potential is
c     built entirely from z, r, x, y (all smooth polynomials in u,
c     from the KS bilinear map) via G_general/XY_general (above),
c     never forming phi or lambda as intermediate coordinates at all.
c
c     Derivation: P_n^m(sin phi) = cos(phi)^m * G(n,m)(sp), and
c     cos(m*lambda)/sin(m*lambda) = X_m(xh,yh)/Y_m(xh,yh) (xh=x/r,
c     yh=y/r) -- the cos(phi)^m = (sxy/r)^m factors cancel EXACTLY
c     between the two, leaving
c       V_{n,m} = (R_body/r)^n / r * G(n,m)(sp) *
c                 [Cnm*X_m(xh,yh) + Snm*Y_m(xh,yh)]
c     with a PURE r-power denominator (via the bounded ratio
c     (R_body/r)^n, never forming R_body^n or r^(n+m+1) as separate
c     large quantities) and no other singular factor anywhere.
c
c     NUMERICAL STABILITY (found and fixed 2026-08-09, same day as
c     the initial reformulation): the first version of this
c     subroutine computed Rn=R_body^n and rpow=r^(n+m+1) as SEPARATE
c     quantities before dividing, and fed XY_general raw x,y instead
c     of the normalized xh=x/r, yh=y/r -- both r^(n+m+1) (for
c     n+m up to ~140 at the sectorial diagonal, nmax=72) and the raw
c     X_m/Y_m (growing like r^m) can each independently overflow
c     double precision well before their ratio/product reaches its
c     true, physically-small final value, even though the
c     MATHEMATICAL result is correct in exact arithmetic --
c     triggered IEEE_OVERFLOW_FLAG in real driver_KS.F runs (caught
c     via a production smoke test, ngeo_deg=50). Fixed by computing
c     the SAFE bounded ratio Rr_n=(R_body/r)^n/r directly (never
c     forming R_body^n or r^(n+m+1) separately) and feeding
c     XY_general the L2-normalized xh,yh=x/r,y/r (bounded by
c     construction, xh^2+yh^2=cos^2(phi)<=1) instead of raw x,y.
c     Re-verified after the fix: zero overflow across LEO through
c     super-GEO-altitude test radii at nmax=72, both realistic-decay
c     and deliberately pathological (non-decaying) coefficient sets
c     (Python/numpy stress test, `np.seterr(all='raise')`, scratch
c     work); unchanged agreement with Cunningham to machine precision
c     at low degree; and a clean driver_KS.exe smoke run (ngeo_deg=50)
c     with no floating-point exception signalled.
c
c     Verified: (1) symbolically (sympy), G(n,m)=P_n^m/cos(phi)^m and
c     X_m/Y_m match ground truth exactly, nmax=6, every (n,m); (2)
c     numerically, this subroutine's v_pert/q_leg/tau_leg reproduce
c     Cunningham's independent solid-harmonic recursion
c     (tess_general_force) to machine precision (~1e-15, random
c     u/Cnm/Snm/nmax, multiple trials) and an independent finite-
c     difference gradient (~1e-10, truncation-error level); (3)
c     dvdt_leg (below) verified against a finite difference of the
c     rotated potential AND against this exact subroutine's own
c     pre-reformulation output for a fixed test case (bit-for-bit
c     match, all outputs, confirming the reformulation is a genuine
c     re-derivation of the identical physics, not a different one);
c     (4) the exact-pole regression (test/test_legendre_tess.F,
c     PL1-PL4) confirms finite, correct output at sqrt(x^2+y^2)=0.
c     See test/test_legendre_tess.F and ALGORITHM.md for the full
c     record.
c
c     dvdt_leg: the mu-free rate of change of this potential w.r.t.
c     PHYSICAL time t, holding u fixed -- see the original derivation
c     rationale (Sellamuthu 2018 PhD thesis eq. 2.54, driver_KS.F's
c     own z(10) assembly). Previously computed via
c     "-theta_dot*Vlam_block" (Vlam_block = dV/d(lambda)); reformulated
c     here as "theta_dot*(yh*dVdxh - xh*dVdyh)" (dVdxh, dVdyh being
c     the gradient w.r.t. the normalized xh,yh, at fixed sp and r),
c     using the same polar-to-Cartesian derivative identity applied
c     to the (xh,yh) pair -- verified to reproduce the original
c     formula's value and sign exactly (bit-for-bit on a fixed test
c     case with real m>0 coefficients), and independently via finite
c     difference of the rotated potential.
c
c     Cnm/Snm are the SAME already-converted, correctly-signed,
c     fully-normalized-to-unnormalized coefficient arrays
c     geo_coeff_tess_general produces -- unchanged interface.
c
c     Rtilt (added 2026-08-10, KSROP issue TBD -- found via GMAT
c     cross-validation of KSROP-Lunar/KSROP-Mars): 3x3 constant
c     orthogonal rotation matrix taking a position vector from the
c     WORKING inertial frame (e.g. EarthMJ2000Eq-parallel axes, as
c     KSROP-Lunar/KSROP-Mars both use) into the body's own "true
c     equator and prime-meridian-node of date" frame, i.e. the frame
c     whose Z axis IS the body's real rotation pole. This subroutine
c     had NO such rotation before 2026-08-10 -- only the per-step
c     prime-meridian (longitude) rotation via rotate_tess_coeffs --
c     which is exactly correct for Earth (EarthMJ2000Eq's Z axis IS
c     Earth's own pole, by definition) but was SILENTLY WRONG for any
c     other body: it implicitly assumed the working frame's Z axis
c     coincides with the target body's pole, which is false by ~35.6
c     deg for Mars and ~23.5 deg for the Moon (IAU pole RA/Dec vs. the
c     Earth-equatorial Z axis) -- caught via a GMAT cross-validation
c     campaign showing ~2%/~0.7% position residuals that persisted
c     unchanged across identical-coefficient, zonal-only, and
c     J2-only isolation tests (ruling out coefficient/model-vintage
c     differences and rotation-PHASE-only explanations), pointing
c     squarely at a missing rotation, not a force-law defect.
c
c     Construction (both downstream repos): Rtilt = R1(90-dec0) .
c     R3(90+ra0), the standard IAU body-fixed-frame prefix (Archinal
c     et al., "Report of the IAU WGCCRE"), using the body's pole
c     right ascension/declination (ra0, dec0) at the run's epoch --
c     the SAME NAIF pck00011.tpc kernel each repo's own prime-meridian
c     angle subroutine (mars_rotation_angle_deg/moon_rotation_angle_
c     deg) already cites, so this is a self-consistent extension of
c     an already-used source, not a new one. Passed as the IDENTITY
c     matrix for Earth (KSROP's own driver_KS.F) -- preserves the
c     already-GMAT-validated (2026-08-10, <0.02% over 3 orbits at
c     degree/order 30) Earth behavior exactly unchanged.
c
c     Implementation: x,y,z,r are computed from u EXACTLY as before
c     (the KS bilinear map itself is central-body-agnostic and
c     unaffected). The resulting UNIT direction (xh,yh,sp) is then
c     rotated by Rtilt BEFORE being fed to G_general/XY_general --
c     i.e. the potential and its angular gradient are evaluated at
c     the body-pole-aligned direction, not the raw working-frame
c     direction. dvdt_leg (needs the longitude-rotation-rate energy
c     term) is computed from these ROTATED xh,yh/dVdxh,dVdyh directly,
c     since the physical longitude rotation happens about the body's
c     TRUE pole (the rotated frame's own Z axis), matching theta_dot.
c     The resulting gradient (dVdxh,dVdyh,dVdsp), which is now
c     expressed in the ROTATED frame, is then rotated BACK to the
c     ORIGINAL (u-native) frame via Rtilt^T before the existing
c     du-chain-rule step -- verified symbolically (sympy) that for a
c     linear coordinate rotation (xh',yh',sp')=Rtilt.(xh,yh,sp), the
c     gradient transforms as (dV/dxh,dV/dyh,dV/dsp) =
c     Rtilt^T.(dV/dxh',dV/dyh',dV/dsp'), the standard covariant
c     gradient transformation for an orthogonal change of coordinates
c     -- dV/dr needs no such transformation (r is rotation-invariant,
c     a true scalar), and v_pert (the potential itself) likewise
c     needs none (also a true scalar). Rtilt=identity reproduces the
c     pre-2026-08-10 formulas bit-for-bit (the rotate-then-unrotate
c     round-trip is then the identity operation) -- regression-tested,
c     see test_legendre_tess.F RT1-RT4 (rotation covariance: evaluating
c     at a position P with Rtilt=R must equal evaluating at R.P with
c     Rtilt=identity, then rotating the resulting force back by R^T).
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension u(4), q_leg(4)
      dimension Rtilt(3,3)
      dimension Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax)
      dimension G(0:nmax,0:nmax), Gd(0:nmax,0:nmax)
      dimension Xt(0:nmax), Yt(0:nmax)
      dimension Xx(0:nmax), Xy(0:nmax), Yx(0:nmax), Yy(0:nmax)
      dimension dzdu(4), drdu(4), dxdu(4), dydu(4), dVdu(4)

      u1 = u(1)
      u2 = u(2)
      u3 = u(3)
      u4 = u(4)

      x = u1*u1 - u2*u2 - u3*u3 + u4*u4
      y = 2.d0*(u1*u2 - u3*u4)
      z = 2.d0*(u1*u3 + u2*u4)
      r_ks = u1*u1 + u2*u2 + u3*u3 + u4*u4

c     Original (working-frame) unit direction, then rotated by Rtilt
c     into the body-pole-aligned frame -- see the Rtilt header note
c     above. xh,yh,sp are OVERWRITTEN with the rotated components;
c     everything below (G_general/XY_general/dvdt_leg) is otherwise
c     completely unchanged and now simply operates on the rotated
c     direction, exactly as it always operated on the raw one.
      xh0 = x / r_ks
      yh0 = y / r_ks
      sp0 = z / r_ks

      xh = Rtilt(1,1)*xh0 + Rtilt(1,2)*yh0 + Rtilt(1,3)*sp0
      yh = Rtilt(2,1)*xh0 + Rtilt(2,2)*yh0 + Rtilt(2,3)*sp0
      sp = Rtilt(3,1)*xh0 + Rtilt(3,2)*yh0 + Rtilt(3,3)*sp0

      call G_general(nmax, sp, G, Gd)
      call XY_general(nmax, xh, yh, Xt, Yt, Xx, Xy, Yx, Yy)

      v_pert = 0.d0
      dVdsp = 0.d0
      dVdxh = 0.d0
      dVdyh = 0.d0
      dVdr = 0.d0

      do n = 2, nmax
c        Safe bounded ratio -- never forms R_body^n or r^(n+m+1) as
c        separate (potentially overflowing) quantities. See header.
         Rr_n = (R_body/r_ks)**n / r_ks
         do m = 0, n
            cc = Cnm(n,m)
            ss = Snm(n,m)
            if (cc .ne. 0.d0 .or. ss .ne. 0.d0) then
               trig = cc*Xt(m) + ss*Yt(m)
               Gnm = G(n,m)
               Gdnm = Gd(n,m)

               v_pert = v_pert + Rr_n * Gnm * trig
               dVdsp = dVdsp + Rr_n * Gdnm * trig

               dtrigdxh = cc*Xx(m) + ss*Yx(m)
               dtrigdyh = cc*Xy(m) + ss*Yy(m)
               dVdxh = dVdxh + Rr_n * Gnm * dtrigdxh
               dVdyh = dVdyh + Rr_n * Gnm * dtrigdyh

c              d(Rr_n)/dr = -(n+1)*Rr_n/r -- reuses the already-safe
c              Rr_n rather than recomputing R_body^n/r^(n+2) from
c              scratch (which would reintroduce the same overflow
c              risk this whole reformulation exists to avoid).
               dVdr = dVdr - dble(n+1)*(Rr_n/r_ks)*Gnm*trig
            end if
         end do
      end do

c     Energy-element rate: the physical longitude rotation happens
c     about the body's TRUE pole, i.e. the ROTATED frame's own Z
c     axis -- must use the still-rotated xh,yh/dVdxh,dVdyh here,
c     BEFORE un-rotating the gradient back below.
      dvdt_leg = theta_dot*(yh*dVdxh - xh*dVdyh)

c     Un-rotate the gradient (dVdxh,dVdyh,dVdsp), currently expressed
c     in the body-pole-aligned frame, back to the ORIGINAL
c     working-frame direction via Rtilt^T -- the standard covariant
c     gradient transformation for an orthogonal coordinate rotation
c     (verified symbolically, see header note). dVdr and v_pert are
c     true scalars (r and the potential itself are rotation-
c     invariant) and need no such transformation.
      dVdxh_r = Rtilt(1,1)*dVdxh + Rtilt(2,1)*dVdyh + Rtilt(3,1)*dVdsp
      dVdyh_r = Rtilt(1,2)*dVdxh + Rtilt(2,2)*dVdyh + Rtilt(3,2)*dVdsp
      dVdsp_r = Rtilt(1,3)*dVdxh + Rtilt(2,3)*dVdyh + Rtilt(3,3)*dVdsp
      dVdxh = dVdxh_r
      dVdyh = dVdyh_r
      dVdsp = dVdsp_r

      dzdu(1) =  2.d0*u3
      dzdu(2) =  2.d0*u4
      dzdu(3) =  2.d0*u1
      dzdu(4) =  2.d0*u2

      drdu(1) =  2.d0*u1
      drdu(2) =  2.d0*u2
      drdu(3) =  2.d0*u3
      drdu(4) =  2.d0*u4

      dxdu(1) =  2.d0*u1
      dxdu(2) = -2.d0*u2
      dxdu(3) = -2.d0*u3
      dxdu(4) =  2.d0*u4

      dydu(1) =  2.d0*u2
      dydu(2) =  2.d0*u1
      dydu(3) = -2.d0*u4
      dydu(4) = -2.d0*u3

c     sp=z/r, xh=x/r, yh=y/r all depend on BOTH the raw coordinate
c     and r itself -- full quotient-rule chain rule through u.
      do k = 1, 4
         dspdu = dzdu(k)/r_ks - z*drdu(k)/(r_ks*r_ks)
         dxhdu = dxdu(k)/r_ks - x*drdu(k)/(r_ks*r_ks)
         dyhdu = dydu(k)/r_ks - y*drdu(k)/(r_ks*r_ks)
         dVdu(k) = dVdsp*dspdu + dVdxh*dxhdu + dVdyh*dyhdu
     &           + dVdr*drdu(k)
      end do

      q_leg(1) = v_pert*u1 + (r_ks/2.d0)*dVdu(1)
      q_leg(2) = v_pert*u2 + (r_ks/2.d0)*dVdu(2)
      q_leg(3) = v_pert*u3 + (r_ks/2.d0)*dVdu(3)
      q_leg(4) = v_pert*u4 + (r_ks/2.d0)*dVdu(4)

      tau_leg = -2.d0*r_ks*v_pert - (r_ks/2.d0)*
     &          (u1*dVdu(1) + u2*dVdu(2) + u3*dVdu(3) + u4*dVdu(4))

      v_pert_out = v_pert

      return
      end
