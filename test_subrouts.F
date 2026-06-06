c     test_subrouts.F  -  Unit tests for KSROP subroutines
c     -------------------------------------------------------
c     Compile:  ifort test_subrouts.F Subrouts.F Legendre.F
c               gfortran test_subrouts.F Subrouts.F Legendre.F
c     Run:      ./a.out   (exit 0 = all pass, exit 1 = any fail)
c     -------------------------------------------------------
c     Tests covered:
c       dotp3, dotp4, vmn, cross          (vector utilities)
c       INTPOL                            (interpolation)
c       aLegP  P2/P3/P4                  (Legendre polynomials)
c       car2ks -> ks2car  roundtrip       (KS transform)
c       car2oe            elements        (orbital elements)
c     -------------------------------------------------------
      program test_subrouts

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

c     --- declarations ---
      dimension a3(3), b3(3)
      dimension a4(4), b4(4)
      dimension v3(3)
      dimension cx(3), cy(3), cz(3)
      dimension xtab(5), ytab(5)
      dimension vleg(2500)
      dimension xr(3), xdr(3), ur(4), usr(4), xout(3), xdout(3)
      dimension xoe(3), xdoe(3), pek(7)
      integer npass, nfail

c     --- initialise common block (must match const_new.dat) ---
      pi      = 4.d0 * datan(1.d0)
      d2r     = pi / 180.d0
      r2d     = 180.d0 / pi
      amue    = 3.986004415d5    ! km3/s2
      AU      = 1.495978707d8    ! km
      R_Earth = 6378.1363d0      ! km

      npass = 0
      nfail = 0
      eps10 = 1.d-10             ! machine-precision tolerance
      eps8  = 1.d-8              ! KS roundtrip tolerance

      write(*,*) '=============================================='
      write(*,*) ' KSROP Unit Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     1. dotp3:  [1,2,3] . [4,5,6] = 4+10+18 = 32
c     ==============================================================
      a3(1) = 1.d0
      a3(2) = 2.d0
      a3(3) = 3.d0
      b3(1) = 4.d0
      b3(2) = 5.d0
      b3(3) = 6.d0
      res = dotp3(a3, b3)
      call assert('dotp3  [1,2,3].[4,5,6] = 32',
     &            res, 32.d0, eps10, npass, nfail)

c     ==============================================================
c     2. dotp4:  [1,2,3,4] . [1,2,3,4] = 1+4+9+16 = 30
c     ==============================================================
      a4(1) = 1.d0
      a4(2) = 2.d0
      a4(3) = 3.d0
      a4(4) = 4.d0
      res = dotp4(a4, a4)
      call assert('dotp4  [1,2,3,4].[1,2,3,4] = 30',
     &            res, 30.d0, eps10, npass, nfail)

c     ==============================================================
c     3. vmn:  |[3,4,0]| = 5  (3-4-5 Pythagorean triple)
c     ==============================================================
      v3(1) = 3.d0
      v3(2) = 4.d0
      v3(3) = 0.d0
      res = vmn(v3)
      call assert('vmn    |[3,4,0]| = 5',
     &            res, 5.d0, eps10, npass, nfail)

c     ==============================================================
c     4. cross:  [1,0,0] x [0,1,0] = [0,0,1]
c     ==============================================================
      cx(1) = 1.d0
      cx(2) = 0.d0
      cx(3) = 0.d0
      cy(1) = 0.d0
      cy(2) = 1.d0
      cy(3) = 0.d0
      call cross(cx, cy, cz)
      call assert('cross  x[1,0,0] result',
     &            cz(1), 0.d0, eps10, npass, nfail)
      call assert('cross  y[1,0,0] result',
     &            cz(2), 0.d0, eps10, npass, nfail)
      call assert('cross  z[1,0,0] result = 1',
     &            cz(3), 1.d0, eps10, npass, nfail)

c     ==============================================================
c     5-6. INTPOL:  table y = 2x,  x in {1,2,3,4,5}
c     ==============================================================
      do 10 i = 1, 5
         xtab(i) = dble(i)
         ytab(i) = 2.d0 * dble(i)
 10   continue

      call INTPOL(xtab, ytab, 5, 1.5d0, res)
      call assert('INTPOL midpoint  y=2x at x=1.5 -> 3.0',
     &            res, 3.d0, eps10, npass, nfail)

      call INTPOL(xtab, ytab, 5, 3.d0, res)
      call assert('INTPOL exact     y=2x at x=3.0 -> 6.0',
     &            res, 6.d0, eps10, npass, nfail)

c     ==============================================================
c     7-9. aLegP:  Legendre polynomials at x = 0.5
c     aLegP fills v such that v(i+1) = P_i(x)  for i = 0..49
c
c     Analytic values:
c       P2(0.5) = (3*0.25 - 1)/2          = -0.125
c       P3(0.5) = (5*0.125 - 1.5)/2       = -0.4375
c       P4(0.5) = (35*0.0625-30*0.25+3)/8 = -0.2890625
c     ==============================================================
      call aLegP(50, 0.5d0, vleg)

      call assert('aLegP  P2(0.5) = -0.125',
     &            vleg(3), -0.125d0, eps10, npass, nfail)
      call assert('aLegP  P3(0.5) = -0.4375',
     &            vleg(4), -0.4375d0, eps10, npass, nfail)
      call assert('aLegP  P4(0.5) = -0.2890625',
     &            vleg(5), -0.2890625d0, eps10, npass, nfail)

c     ==============================================================
c     10-15. car2ks -> ks2car roundtrip
c     Orbit at r=7000 km along x-axis, v=7.5 km/s along y-axis
c     (slightly sub-circular elliptic orbit in the equatorial plane)
c     After transform and inverse transform, x and xd must be recovered
c     to near machine precision.
c     ==============================================================
      xr(1)  = 7000.d0
      xr(2)  = 0.d0
      xr(3)  = 0.d0
      xdr(1) = 0.d0
      xdr(2) = 7.5d0
      xdr(3) = 0.d0

      rmag  = 7000.d0
      vmag2 = 7.5d0 * 7.5d0
c     KS frequency (two-body, no perturbation potential)
      ww = dsqrt(0.5d0 * (amue/rmag - vmag2/2.d0))

      call car2ks(xr, xdr, ur, usr, ww)
      call ks2car(ur, usr, xout, xdout, ww)

      call assert('car2ks/ks2car  x(1) = 7000 km',
     &            xout(1), xr(1), eps8, npass, nfail)
      call assert('car2ks/ks2car  x(2) = 0',
     &            xout(2), xr(2), eps8, npass, nfail)
      call assert('car2ks/ks2car  x(3) = 0',
     &            xout(3), xr(3), eps8, npass, nfail)
      call assert('car2ks/ks2car  xd(1) = 0',
     &            xdout(1), xdr(1), eps8, npass, nfail)
      call assert('car2ks/ks2car  xd(2) = 7.5 km/s',
     &            xdout(2), xdr(2), eps8, npass, nfail)
      call assert('car2ks/ks2car  xd(3) = 0',
     &            xdout(3), xdr(3), eps8, npass, nfail)

c     ==============================================================
c     16-18. car2oe:  known elliptic inclined orbit
c
c     Setup:  a=7000 km, e=0.1, i=45 deg, RAAN=0, AOP=0
c     At perigee (true anomaly = 0):
c       r_peri = a*(1-e) = 6300 km  along x-axis
c       v_peri = sqrt(mu*(1+e)/(a*(1-e)))  perpendicular to r
c       inclination encoded as:  xd = v_peri * [0, cos45, sin45]
c
c     Expected outputs from car2oe:
c       pek(1) = a = 7000 km
c       pek(2) = e = 0.1
c       pek(3) = i = 45 deg
c     ==============================================================
      a_test  = 7000.d0
      e_test  = 0.1d0
      ang_i   = 45.d0          ! 'i' prefix is integer; use 'a'-range name

      r_p   = a_test * (1.d0 - e_test)
      v_p   = dsqrt(amue * (1.d0 + e_test) / (a_test * (1.d0 - e_test)))

      xoe(1) = r_p
      xoe(2) = 0.d0
      xoe(3) = 0.d0
      xdoe(1) = 0.d0
      xdoe(2) = v_p * dcos(ang_i * d2r)
      xdoe(3) = v_p * dsin(ang_i * d2r)

      call car2oe(xoe, xdoe, pek)

c     pek(3..7) are returned in degrees by car2oe
      call assert('car2oe  a = 7000 km',
     &            pek(1), a_test, 1.d0, npass, nfail)
      call assert('car2oe  e = 0.1',
     &            pek(2), e_test, 1.d-4, npass, nfail)
      call assert('car2oe  i = 45 deg',
     &            pek(3), ang_i, 1.d-5, npass, nfail)

c     ==============================================================
c     Summary
c     ==============================================================
      write(*,*) '=============================================='
      write(*,200) npass + nfail, npass, nfail
 200  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='

      if (nfail .gt. 0) stop 1
      stop
      end

c     ============================================================
      subroutine assert(label, got, expected, tol, npass, nfail)
c     Absolute tolerance check: |got - expected| <= tol
c     Reports PASS/FAIL and updates counters.
c     ============================================================
      implicit double precision (a-h, o-z)
      character*(*) label
      integer npass, nfail

      diff = dabs(got - expected)
      if (diff .le. tol) then
         write(*,10) 'PASS', label
         npass = npass + 1
      else
         write(*,20) 'FAIL', label, got, expected, diff
         nfail = nfail + 1
      end if
 10   format('  [',a4,']  ',a)
 20   format('  [',a4,']  ',a,/,
     &       '           got=',e20.13,'  exp=',e20.13,
     &       '  diff=',e10.3)
      return
      end
