c     ============================================================
      subroutine G_general(nmax, sp, G, Gd)
c     KSROP branch feature/tesseral-energy-time-dependence (2026-08-09):
c     "de-cos(phi)^m" associated Legendre recursion --
c       G(n,m) = P_n^m(sin phi) / cos(phi)^m
c     computed directly as a function of sp=sin(phi)=z/r ALONE, with
c     no cos(phi) dependence anywhere in the recursion itself:
c       G(n,n) = (2n-1) G(n-1,n-1)                         (sectorial)
c       G(n,m) = [(2n-1) sp G(n-1,m) - (n+m-1) G(n-2,m)]
c                / (n-m),                          0 <= m < n
c     -- the SAME three-term recurrence shape as the classical
c     associated-Legendre recursion (formerly alfP_general, removed),
c     but with the sectorial step's cos(phi) factor dropped, since
c     P_n^m(sin phi) = cos(phi)^m * G(n,m) exactly (verified
c     symbolically, sympy, zero residual against the true P_n^m
c     recursion for every (n,m), nmax=6).
c
c     Gd(n,m) = dG(n,m)/d(sp), obtained by differentiating this SAME
c     recursion directly term-by-term (the standard trick for a
c     non-singular Legendre-type derivative recursion -- no (1-sp^2)
c     or cos(phi)-type denominator anywhere):
c       Gd(n,n) = (2n-1) Gd(n-1,n-1)
c       Gd(n,m) = [(2n-1)(G(n-1,m)+sp*Gd(n-1,m)) - (n+m-1) Gd(n-2,m)]
c                 / (n-m)
c     Verified symbolically (sympy) to exactly match d/d(sp)[G(n,m)]
c     computed by direct differentiation, nmax=6, every (n,m).
c
c     Why this eliminates the polar singularity: the classical
c     approach (formerly alfP_general + a phi-derivative identity
c     dividing by cos(phi), i.e. by sqrt(x^2+y^2)/r) blows up as the
c     trajectory approaches either pole. G(n,m) and Gd(n,m) are pure
c     functions of sp=z/r -- smooth for all u except r=0 (a genuine
c     physical singularity every formulation shares), with NO
c     coordinate-artifact singularity at the poles.
c
c     Both G and Gd feed tess_legendre_force's potential/gradient sum
c     alongside XY_general's own singularity-free X_m/Y_m (below) --
c     see that subroutine's header, and the branch's derivation notes
c     (ALGORITHM.md), for the full assembly and its independent
c     verification (machine-precision cross-check against Cunningham's
c     solid-harmonic recursion, plus finite-difference).
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
      subroutine XY_general(nmax, x, y, Xt, Yt, Xx, Xy, Yx, Yy)
c     KSROP branch feature/tesseral-energy-time-dependence (2026-08-09):
c     unnormalized cos(m*lambda)*sxy^m and sin(m*lambda)*sxy^m
c     (sxy=sqrt(x^2+y^2)), computed as pure polynomials in x,y via the
c     De Moivre / Chebyshev-type recursion applied to (x+iy)^m --
c     WITHOUT ever forming lambda=atan2(y,x) or dividing by sxy the
c     way the classical cos(m*lambda)=Re[(x+iy)^m]/sxy^m evaluation
c     would (the division is exactly what re-introduces a polar
c     singularity as sxy->0). Verified symbolically (sympy) to equal
c     Re[(x+iy)^m]/Im[(x+iy)^m] exactly for m=0..7.
c       X_m = 2x X_{m-1} - (x^2+y^2) X_{m-2},   X_0=1, X_1=x
c       Y_m = 2x Y_{m-1} - (x^2+y^2) Y_{m-2},   Y_0=0, Y_1=y
c     Xx,Xy,Yx,Yy are the partial derivatives d/dx, d/dy of X_m, Y_m,
c     obtained by differentiating this same recursion term-by-term
c     (verified symbolically against direct differentiation, m=0..7)
c     -- again no singular denominator anywhere.
c
c     Together with G_general (above), these replace tess_legendre_
c     force's former phi/lambda spherical-coordinate chain rule
c     (alfP_general + cos_sin_mlambda + an explicit division by sxy
c     in the u-derivative of phi and lambda) with a fully Cartesian
c     (x,y,z,r) formulation that has no coordinate-artifact
c     singularity at the poles -- see ALGORITHM.md's known-limitations
c     entry this branch resolves.
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
         Xt(1) = x
         Yt(1) = y
         Xx(1) = 1.d0
         Xy(1) = 0.d0
         Yx(1) = 0.d0
         Yy(1) = 1.d0
      end if

      sxy2 = x*x + y*y
      do m = 2, nmax
         Xt(m) = 2.d0*x*Xt(m-1) - sxy2*Xt(m-2)
         Yt(m) = 2.d0*x*Yt(m-1) - sxy2*Yt(m-2)
         Xx(m) = 2.d0*Xt(m-1) + 2.d0*x*Xx(m-1)
     &           - 2.d0*x*Xt(m-2) - sxy2*Xx(m-2)
         Xy(m) = 2.d0*x*Xy(m-1) - 2.d0*y*Xt(m-2) - sxy2*Xy(m-2)
         Yx(m) = 2.d0*Yt(m-1) + 2.d0*x*Yx(m-1)
     &           - 2.d0*x*Yt(m-2) - sxy2*Yx(m-2)
         Yy(m) = 2.d0*x*Yy(m-1) - 2.d0*y*Yt(m-2) - sxy2*Yy(m-2)
      end do

      return
      end

c     ============================================================
      subroutine tess_legendre_force(u, R_body, theta_dot, nmax,
     &                                Cnm, Snm,
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
c     cos(m*lambda)/sin(m*lambda) = X_m(x,y)/Y_m(x,y) / sxy^m -- the
c     cos(phi)^m = (sxy/r)^m factors cancel EXACTLY between the two,
c     leaving
c       V_{n,m} = (R_body^n / r^(n+m+1)) * G(n,m)(sp) *
c                 [Cnm*X_m(x,y) + Snm*Y_m(x,y)]
c     with a PURE r-power denominator (r^(n+m+1)) and no other
c     singular factor anywhere. Verified: (1) symbolically (sympy),
c     G(n,m)=P_n^m/cos(phi)^m and X_m/Y_m match ground truth exactly,
c     nmax=6, every (n,m); (2) numerically, this subroutine's v_pert/
c     q_leg/tau_leg reproduce Cunningham's independent solid-harmonic
c     recursion (tess_general_force) to machine precision (~1e-15,
c     random u/Cnm/Snm/nmax, multiple trials) and an independent
c     finite-difference gradient (~1e-10, truncation-error level);
c     (3) dvdt_leg (below) verified against a finite difference of the
c     rotated potential AND against this exact subroutine's own
c     pre-reformulation output for a fixed test case (bit-for-bit
c     match, all outputs, confirming the reformulation is a genuine
c     re-derivation of the identical physics, not a different one).
c     See test/test_legendre_tess.F and ALGORITHM.md for the full
c     record.
c
c     dvdt_leg: the mu-free rate of change of this potential w.r.t.
c     PHYSICAL time t, holding u fixed -- see the original derivation
c     rationale (Sellamuthu 2018 PhD thesis eq. 2.54, driver_KS.F's
c     own z(10) assembly). Previously computed via
c     "-theta_dot*Vlam_block" (Vlam_block = dV/d(lambda)); reformulated
c     here as "theta_dot*(y*dVdx - x*dVdy)", using the identity
c     dV/d(lambda) = x*dVdy - y*dVdx (standard polar-to-Cartesian
c     derivative chain rule, y=r_xy*sin(lambda), x=r_xy*cos(lambda))
c     -- verified to reproduce the original formula's value and sign
c     exactly (bit-for-bit on a fixed test case with real m>0
c     coefficients), and independently via finite difference of the
c     rotated potential.
c
c     Cnm/Snm are the SAME already-converted, correctly-signed,
c     fully-normalized-to-unnormalized coefficient arrays
c     geo_coeff_tess_general produces -- unchanged interface.
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension u(4), q_leg(4)
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

      sp = z / r_ks

      call G_general(nmax, sp, G, Gd)
      call XY_general(nmax, x, y, Xt, Yt, Xx, Xy, Yx, Yy)

      v_pert = 0.d0
      dVdz = 0.d0
      dVdr = 0.d0
      dVdx = 0.d0
      dVdy = 0.d0

      do n = 2, nmax
         Rn = R_body**n
         do m = 0, n
            cc = Cnm(n,m)
            ss = Snm(n,m)
            if (cc .ne. 0.d0 .or. ss .ne. 0.d0) then
               rpow = r_ks**(n+m+1)
               trig = cc*Xt(m) + ss*Yt(m)
               Gnm = G(n,m)
               Gdnm = Gd(n,m)

               v_pert = v_pert + Rn/rpow * Gnm * trig

               dVdz = dVdz + Rn/rpow * Gdnm/r_ks * trig
               dVdr = dVdr + Rn*trig*
     &               ( Gdnm*(-z/(r_ks*r_ks))/rpow
     &                 - dble(n+m+1)*Gnm/(rpow*r_ks) )

               dtrigdx = cc*Xx(m) + ss*Yx(m)
               dtrigdy = cc*Xy(m) + ss*Yy(m)
               dVdx = dVdx + Rn/rpow * Gnm * dtrigdx
               dVdy = dVdy + Rn/rpow * Gnm * dtrigdy
            end if
         end do
      end do

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

      do k = 1, 4
         dVdu(k) = dVdz*dzdu(k) + dVdr*drdu(k)
     &           + dVdx*dxdu(k) + dVdy*dydu(k)
      end do

      q_leg(1) = v_pert*u1 + (r_ks/2.d0)*dVdu(1)
      q_leg(2) = v_pert*u2 + (r_ks/2.d0)*dVdu(2)
      q_leg(3) = v_pert*u3 + (r_ks/2.d0)*dVdu(3)
      q_leg(4) = v_pert*u4 + (r_ks/2.d0)*dVdu(4)

      tau_leg = -2.d0*r_ks*v_pert - (r_ks/2.d0)*
     &          (u1*dVdu(1) + u2*dVdu(2) + u3*dVdu(3) + u4*dVdu(4))

      dvdt_leg = theta_dot*(y*dVdx - x*dVdy)
      v_pert_out = v_pert

      return
      end
