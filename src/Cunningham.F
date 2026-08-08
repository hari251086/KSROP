c     ============================================================
      subroutine cunningham_Vnm(nwork, x, y, z, r2, Vr, Vi)
c     KSROP issue #30: general (n,m) tesseral/mascon harmonic support.
c
c     Cunningham (1970, Celestial Mechanics 2:207-216) recursive
c     complex solid spherical harmonics, computed in body-fixed
c     Cartesian coordinates (no polar/equatorial singularity, unlike
c     a direct associated-Legendre-in-latitude evaluation). Returns
c     the real/imaginary parts of
c       V_{n,m} = P_n^m(sin phi) * exp(i*m*lambda) / r^(n+1)
c     for all 0<=m<=n<=nwork, via the two real recurrences (eqs 14,
c     16, 17 of the paper):
c       V_{n,n}   = (2n-1)(x+iy)/r^2 * V_{n-1,n-1},   V_{0,0}=1/r
c       (n-m)V_{n,m} = (2n-1)(z/r^2)V_{n-1,m}
c                      - (n+m-1)/r^2 * V_{n-2,m}        (n-2>=m)
c                    = (2n-1)(z/r^2)V_{n-1,m}            (n=m+1)
c     Tracked as separate real Vr/Vi arrays throughout (both
c     recurrences reduce to pure real arithmetic component-wise) --
c     no COMPLEX type, consistent with this codebase's style.
c
c     Verified (2026-08-08, see KSROP#30) against the paper's own
c     Table I closed forms for 9 (n,m) pairs up to (4,4), machine
c     precision; and against a central finite-difference gradient
c     check on the companion cunningham_dVnm derivatives below.
c
c     nwork must be >= (largest n actually needed)+1, since the
c     force-law derivatives (cunningham_dVnm) need one degree beyond
c     the potential's own degree, same convention as aLegP's own
c     "ngeo_deg+1" callers in driver_KS.F.
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension Vr(0:nwork,0:nwork), Vi(0:nwork,0:nwork)

      r = dsqrt(r2)
      Vr(0,0) = 1.d0/r
      Vi(0,0) = 0.d0

      do n = 1, nwork
         pr = Vr(n-1,n-1)
         pim = Vi(n-1,n-1)
         Vr(n,n) = dble(2*n-1)*(x*pr - y*pim)/r2
         Vi(n,n) = dble(2*n-1)*(x*pim + y*pr)/r2

         do m = n-1, 0, -1
            if (n-2 .ge. m) then
               Vr(n,m) = (dble(2*n-1)*z*Vr(n-1,m)
     &                  - dble(n+m-1)*Vr(n-2,m)) / (dble(n-m)*r2)
               Vi(n,m) = (dble(2*n-1)*z*Vi(n-1,m)
     &                  - dble(n+m-1)*Vi(n-2,m)) / (dble(n-m)*r2)
            else
               Vr(n,m) = dble(2*n-1)*z*Vr(n-1,m) / (dble(n-m)*r2)
               Vi(n,m) = dble(2*n-1)*z*Vi(n-1,m) / (dble(n-m)*r2)
            end if
         end do
      end do

      return
      end

c     ============================================================
      subroutine cunningham_dVnm(n, m, nwork, Vr, Vi,
     &                           dxr, dxi, dyr, dyi, dzr, dzi)
c     Derivatives of V_{n,m} w.r.t. body-fixed x,y,z (Cunningham
c     1970 sec. 2, simplified to elementary integer-factor form --
c     the paper's own (n-m+2)!/(2(n-m)!) and (n-m+1)!/(n-m)! ratios
c     reduce to (n-m+2)(n-m+1)/2 and (n-m+1) respectively, avoiding
c     factorial overflow at high degree):
c       dV/dx = -V_{n+1,m+1}/2 + (n-m+2)(n-m+1)/2*V_{n+1,m-1}  (m>0)
c             = -V_{n+1,1}                                     (m=0)
c       dV/dy = +i*V_{n+1,m+1}/2 + i*(n-m+2)(n-m+1)/2*V_{n+1,m-1} (m>0)
c             = -V_{n+1,1}*i -> real part -Vi(n+1,1)            (m=0)
c       dV/dz = -(n-m+1)*V_{n+1,m}                              (m>=0)
c     (m=0 cases obtained by substituting V*_{n+1,1} = conjugate of
c     V_{n+1,1} into the paper's own m=0 formulas and simplifying;
c     verified by finite difference, see cunningham_Vnm header.)
c     Requires Vr/Vi already populated up to degree n+1 (nwork>=n+1).
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension Vr(0:nwork,0:nwork), Vi(0:nwork,0:nwork)

      if (m .gt. 0) then
         cfac = dble((n-m+2)*(n-m+1))/2.d0
         dxr = -Vr(n+1,m+1)/2.d0 + cfac*Vr(n+1,m-1)
         dxi = -Vi(n+1,m+1)/2.d0 + cfac*Vi(n+1,m-1)
         dyr = -Vi(n+1,m+1)/2.d0 - cfac*Vi(n+1,m-1)
         dyi =  Vr(n+1,m+1)/2.d0 + cfac*Vr(n+1,m-1)
      else
         dxr = -Vr(n+1,1)
         dxi = 0.d0
         dyr = -Vi(n+1,1)
         dyi = 0.d0
      end if

      cz = dble(n-m+1)
      dzr = -cz*Vr(n+1,m)
      dzi = -cz*Vi(n+1,m)

      return
      end

c     ============================================================
      subroutine geo_coeff_tess_general(fname, nmax, Cnm, Snm)
c     KSROP issue #30: general (n,m) coefficient loader, replacing
c     geo_coeff_tess22's single-row (2,2) special case. Reads every
c     (n,m) row with 2<=n<=nmax, 0<=m<=n from an EGM2008-row-format
c     gravity file, converts normalized C-bar/S-bar to the
c     unnormalized, correctly-signed values this repo's force-model
c     convention expects.
c
c     Normalization: N_{n,m} = sqrt[(n-m)!(2n+1)(2-delta_0m)/(n+m)!]
c     (XJR thesis eq. 3.3.1-3.3.2, same as geo_coeff_tess22), computed
c     as a running product over the 2m-term ratio (n+m)!/(n-m)!
c     rather than separate factorials, to stay numerically stable at
c     higher degree.
c
c     Scope is deliberately m>0 (tesseral/sectorial) only -- m=0
c     (zonal) rows are always skipped, even if present in the file.
c     KSROP's existing geo_coeff/aleg pathway already handles zonal
c     terms; loading them here too would double-count them once this
c     loader's output is summed into the same force alongside that
c     existing path (a real bug caught 2026-08-08 via the
c     driver_KS.F integration test -- an unexpectedly large
c     perturbation traced directly to a first version of this loader
c     that did not skip m=0). No sign flip is needed for m>0 rows,
c     matching geo_coeff_tess22's own convention directly -- already
c     independently validated ((2,2) via #29 Phase 5 and the
c     Cunningham cross-check in #30).
c
c     Missing (n,m) rows are left at zero (auto-detect: a zonal-only
c     or lower-degree file leaves the corresponding entries inert,
c     same convention as geo_coeff_tess22).
c     ============================================================
      implicit double precision (a-h,o-z)
      character*(*) fname
      integer ideg, iord
      dimension Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax)

      do n = 0, nmax
         do m = 0, nmax
            Cnm(n,m) = 0.d0
            Snm(n,m) = 0.d0
         end do
      end do

      open (unit = 7, file = fname)

c     m=0 (zonal) rows are deliberately skipped -- KSROP's existing
c     geo_coeff/aleg pathway already handles zonal terms; loading
c     them here too would double-count them when this loader's
c     output is summed into the same force alongside that existing
c     path (caught 2026-08-08 via the driver_KS.F integration test:
c     an unexpectedly large perturbation traced directly to this).
c     This loader's scope is genuinely m>0 tesseral/sectorial terms
c     only, matching geo_coeff_tess22's own original (2,2) scope.
 20   continue
         read(7,*,end=10) ideg, iord, cnmv, snmv, sg1, sg2
         if (ideg .ge. 2 .and. ideg .le. nmax
     &       .and. iord .ge. 1 .and. iord .le. ideg) then
            ratio_fact = 1.d0
            do k = ideg-iord+1, ideg+iord
               ratio_fact = ratio_fact * dble(k)
            end do
            xnorm = dsqrt(dble(2*ideg+1)*2.d0/ratio_fact)
            Cnm(ideg,iord) = xnorm*cnmv
            Snm(ideg,iord) = xnorm*snmv
         end if
      goto 20
 10   continue

      close(7)
      return
      end

c     ============================================================
      subroutine tess_general_force(u, R_body, nmax, Cnm, Snm,
     &                               q_gen, tau_gen)
c     KSROP issue #30: general (n,m) tesseral/mascon KS-element force
c     (q_gen) and time-element (tau_gen) contribution, replacing
c     tess22_force's hardcoded (2,2)-only formula with a genuine
c     runtime recurrence over every (n,m) loaded by
c     geo_coeff_tess_general.
c
c     Method (derived and verified 2026-08-08, see KSROP#30):
c     1. Map u (KS position) to body-fixed Cartesian x via KSROP's
c        own bilinear KS map (car2ks/ks2car).
c     2. Evaluate the potential V and its Cartesian gradient at that
c        point via cunningham_Vnm/cunningham_dVnm, summed over every
c        loaded (n,m) term.
c     3. Chain-rule dV/du_j = sum_k (dV/dx_k)(dx_k/du_j), using the
c        closed-form KS position Jacobian dx_k/du_j = 2*L(u)_{kj}
c        (explicit from the same bilinear map, degree-1 in u).
c     4. Apply the same general recipe already used for every other
c        force term in this repo: q(j) = V*u(j) + (r/2)*dV/du_j,
c        tau = -2*r*V - (r/2)*sum_j u(j)*dV/du_j (r = dotp4(u,u), the
c        KS identity for |x|).
c
c     Verified to reproduce tess22_force exactly (diff ~1e-14/1e-15,
c     floating-point noise) when Cnm/Snm contain only the (2,2) term
c     -- the acceptance gate for trusting this general formula for
c     every other (n,m).
c
c     Cnm/Snm are the already-converted, correctly-signed coefficient
c     arrays from geo_coeff_tess_general -- no further sign/
c     normalization handling here.
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension u(4), q_gen(4)
      dimension Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax)
      dimension Vr(0:nmax+1,0:nmax+1), Vi(0:nmax+1,0:nmax+1)

      u1 = u(1)
      u2 = u(2)
      u3 = u(3)
      u4 = u(4)

      x1 = u1*u1 - u2*u2 - u3*u3 + u4*u4
      x2 = 2.d0*(u1*u2 - u3*u4)
      x3 = 2.d0*(u1*u3 + u2*u4)
      r2 = x1*x1 + x2*x2 + x3*x3
      r_ks = u1*u1 + u2*u2 + u3*u3 + u4*u4

      call cunningham_Vnm(nmax+1, x1, x2, x3, r2, Vr, Vi)

      Vtot = 0.d0
      gx = 0.d0
      gy = 0.d0
      gz = 0.d0

      do n = 2, nmax
         Rn = R_body**n
         do m = 0, n
            if (Cnm(n,m) .ne. 0.d0 .or. Snm(n,m) .ne. 0.d0) then
               Vtot = Vtot
     &              + Rn*(Cnm(n,m)*Vr(n,m) + Snm(n,m)*Vi(n,m))
               call cunningham_dVnm(n, m, nmax+1, Vr, Vi,
     &                              dxr, dxi, dyr, dyi, dzr, dzi)
               gx = gx + Rn*(Cnm(n,m)*dxr + Snm(n,m)*dxi)
               gy = gy + Rn*(Cnm(n,m)*dyr + Snm(n,m)*dyi)
               gz = gz + Rn*(Cnm(n,m)*dzr + Snm(n,m)*dzi)
            end if
         end do
      end do

      dVdu1 =  2.d0*u1*gx + 2.d0*u2*gy + 2.d0*u3*gz
      dVdu2 = -2.d0*u2*gx + 2.d0*u1*gy + 2.d0*u4*gz
      dVdu3 = -2.d0*u3*gx - 2.d0*u4*gy + 2.d0*u1*gz
      dVdu4 =  2.d0*u4*gx - 2.d0*u3*gy + 2.d0*u2*gz

      q_gen(1) = Vtot*u1 + (r_ks/2.d0)*dVdu1
      q_gen(2) = Vtot*u2 + (r_ks/2.d0)*dVdu2
      q_gen(3) = Vtot*u3 + (r_ks/2.d0)*dVdu3
      q_gen(4) = Vtot*u4 + (r_ks/2.d0)*dVdu4

      tau_gen = -2.d0*r_ks*Vtot - (r_ks/2.d0)*
     &          (u1*dVdu1 + u2*dVdu2 + u3*dVdu3 + u4*dVdu4)

      return
      end

c     ============================================================
      subroutine rotate_tess_coeffs(nmax, Cnm, Snm, theta_deg,
     &                               Cnm_rot, Snm_rot)
c     KSROP issue #30: generalizes the (2,2)-specific GMST-rotation
c     trick already used for A_tess/B_tess (issue #29) to every
c     loaded (n,m) term. tess_general_force operates on inertial
c     x,y,z (from the KS bilinear map, same as tess22_force always
c     has), not body-fixed coordinates, so each order-m coefficient
c     pair needs its own m*theta rotation -- not the fixed 2*theta
c     the (2,2)-only case used -- to compensate for using inertial
c     longitude instead of body-fixed longitude lambda_body =
c     lambda_inertial - theta:
c       C'_{n,m} = C_{n,m}*cos(m*theta) - S_{n,m}*sin(m*theta)
c       S'_{n,m} = C_{n,m}*sin(m*theta) + S_{n,m}*cos(m*theta)
c     Verified (2026-08-08, KSROP#30) by direct trigonometric
c     expansion that this exactly reproduces
c     C_{n,m}*cos(m*lambda_body)+S_{n,m}*sin(m*lambda_body) when
c     evaluated at inertial longitude -- the same identity the
c     (2,2)-only A_tess/B_tess formula already relies on, here
c     generalized to arbitrary m (m=0 zonal terms are correctly
c     left unrotated: cos(0)=1, sin(0)=0, identity transform).
c     ============================================================
      implicit double precision (a-h,o-z)
      double precision mth
      dimension Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax)
      dimension Cnm_rot(0:nmax,0:nmax), Snm_rot(0:nmax,0:nmax)

      pi = 4.d0*datan(1.d0)
      d2r = pi/180.d0

      do n = 0, nmax
         do m = 0, n
            mth = dble(m)*theta_deg*d2r
            cm = dcos(mth)
            sm = dsin(mth)
            Cnm_rot(n,m) = Cnm(n,m)*cm - Snm(n,m)*sm
            Snm_rot(n,m) = Cnm(n,m)*sm + Snm(n,m)*cm
         end do
      end do

      return
      end
