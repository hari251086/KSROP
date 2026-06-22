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
c       read_cdm          public sample CDM (CCSDS 508.0-B-1 Annex)
c       write_cdm/read_cdm  roundtrip
c       shadfncyl         cylindrical shadow (sunlit + eclipse)
c       shadfncone        conical shadow     (sunlit + umbra)
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
      dimension tca_cal(6), x1(3), xd1(3), cov1(6,6)
      dimension x2(3), xd2(3), cov2(6,6)
      dimension tca_cal2(6), x1b(3), xd1b(3), cov1b(6,6)
      dimension x2b(3), xd2b(3), cov2b(6,6)
      character*23 creation_str, tca_str
      integer npass, nfail, iall_ok, iq1, iq2, iq3, iq4
c     'm' is implicitly INTEGER (outside a-h, o-z) -- declare explicitly
      double precision miss_dist, miss_dist2

c     --- initialise common block from const_new.dat ---
      call init_constants()

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
c     19-30. read_cdm  on a public sample CDM
c     Sample is the CCSDS 508.0-B-1 (Annex/Section 4) example message,
c     a JSPOC conjunction notice for SATELLITE A vs FENGYUN 1C DEB,
c     stored at input/cdm_sample.kvn. Expected values below are taken
c     directly from that public message. read_cdm stores raw numeric
c     values without unit conversion (mirroring read_opm): state
c     vectors are in km/km-s as published, covariances in m**2 family.
c     ==============================================================
      open(unit=21, file='input/cdm_sample.kvn', status='old')
      call read_cdm(21, tca_cal, miss_dist, pc,
     &              x1, xd1, cov1, x2, xd2, cov2)
      close(21)

      call assert('read_cdm  TCA year  = 2010',
     &            tca_cal(1), 2010.d0, eps10, npass, nfail)
      call assert('read_cdm  TCA day   = 13',
     &            tca_cal(3), 13.d0, eps10, npass, nfail)
      call assert('read_cdm  MISS_DISTANCE = 715 m',
     &            miss_dist, 715.d0, eps10, npass, nfail)
      call assert('read_cdm  COLLISION_PROBABILITY = 4.835E-05',
     &            pc, 4.835d-05, 1.d-9, npass, nfail)

      call assert('read_cdm  OBJECT1 X     = 2570.097065 km',
     &            x1(1), 2570.097065d0, eps8, npass, nfail)
      call assert('read_cdm  OBJECT1 Z_DOT = -3.526774282 km/s',
     &            xd1(3), -3.526774282d0, eps8, npass, nfail)
      call assert('read_cdm  OBJECT1 CT_T  = 2.533E+03 m**2',
     &            cov1(2,2), 2.533d+03, 1.d-6, npass, nfail)
      call assert('read_cdm  OBJECT1 CRDOT_RDOT = 5.744E-03',
     &            cov1(4,4), 5.744d-03, 1.d-9, npass, nfail)
      call assert('read_cdm  OBJECT1 covariance symmetric (CT_R)',
     &            cov1(2,1), cov1(1,2), eps10, npass, nfail)

      call assert('read_cdm  OBJECT2 X      = 2569.540800 km',
     &            x2(1), 2569.540800d0, eps8, npass, nfail)
      call assert('read_cdm  OBJECT2 Y_DOT  = -6.007247516 km/s',
     &            xd2(2), -6.007247516d0, eps8, npass, nfail)
      call assert('read_cdm  OBJECT2 CNDOT_NDOT = 5.178E-05',
     &            cov2(6,6), 5.178d-05, 1.d-10, npass, nfail)
      call assert('read_cdm  OBJECT2 covariance symmetric (CN_T)',
     &            cov2(3,2), cov2(2,3), eps10, npass, nfail)

c     ==============================================================
c     31-40. write_cdm -> read_cdm  roundtrip
c     Re-emit the parsed sample as a new CDM (output/cdm_roundtrip.kvn)
c     and re-parse it; the recovered fields must match the originally
c     parsed values. TCA/CREATION_DATE strings are supplied directly
c     (bypassing cal2jd/jd2epoch) so this test stays isolated from the
c     known cal2jd epoch issue tracked separately for driver_KS. Like
c     read_opm/write_opm, neither routine converts units -- write_cdm
c     labels fields with its fixed default brackets ([km],[km**2],...)
c     regardless of the units the input numbers actually carry, so the
c     roundtrip checks numeric fidelity only, not unit semantics.
c     ==============================================================
      creation_str = '2026-06-07T00:00:00.000'
      tca_str      = '2010-03-13T22:37:52.618'

      open(unit=22, file='output/cdm_roundtrip.kvn', status='replace')
      call write_cdm(22, creation_str, tca_str, miss_dist,
     &               14762.d0, pc,
     &               'SATELLITE A',     x1, xd1, cov1,
     &               'FENGYUN 1C DEB',  x2, xd2, cov2)
      close(22)

      open(unit=22, file='output/cdm_roundtrip.kvn', status='old')
      call read_cdm(22, tca_cal2, miss_dist2, pc2,
     &              x1b, xd1b, cov1b, x2b, xd2b, cov2b)
      close(22)

      call assert('roundtrip  TCA day matches',
     &            tca_cal2(3), tca_cal(3), eps10, npass, nfail)
      call assert('roundtrip  TCA seconds match',
     &            tca_cal2(6), tca_cal(6), 1.d-3, npass, nfail)
      call assert('roundtrip  MISS_DISTANCE matches',
     &            miss_dist2, miss_dist, 1.d-6, npass, nfail)
      call assert('roundtrip  COLLISION_PROBABILITY matches',
     &            pc2, pc, 1.d-9, npass, nfail)
      call assert('roundtrip  OBJECT1 X matches',
     &            x1b(1), x1(1), eps8, npass, nfail)
      call assert('roundtrip  OBJECT1 Z_DOT matches',
     &            xd1b(3), xd1(3), eps8, npass, nfail)
      call assert('roundtrip  OBJECT1 CT_T matches',
     &            cov1b(2,2), cov1(2,2), 1.d-6, npass, nfail)
      call assert('roundtrip  OBJECT1 CRDOT_RDOT matches',
     &            cov1b(4,4), cov1(4,4), 1.d-9, npass, nfail)
      call assert('roundtrip  OBJECT2 X matches',
     &            x2b(1), x2(1), eps8, npass, nfail)
      call assert('roundtrip  OBJECT2 CNDOT_NDOT matches',
     &            cov2b(6,6), cov2(6,6), 1.d-10, npass, nfail)

c     ==============================================================
c     21. shadfncyl: satellite clearly on sunlit side -> ROHP = 1.0
c         sat=[8000,0,0] km, sun=[1.5e8,0,0] km  (same hemisphere as sun)
c     ==============================================================
      rohp_cyl = shadfncyl(8.d3, 0.d0, 0.d0,
     &                     1.5d8, 0.d0, 0.d0)
      call assert('shadfncyl  sunlit side -> 1.0',
     &            rohp_cyl, 1.d0, eps10, npass, nfail)

c     ==============================================================
c     22. shadfncyl: satellite inside cylindrical shadow -> ROHP = 0.0
c         sat=[-8000,0,0] km (anti-sun), sun=[1.5e8,0,0] km
c     ==============================================================
      rohp_cyl = shadfncyl(-8.d3, 0.d0, 0.d0,
     &                      1.5d8, 0.d0, 0.d0)
      call assert('shadfncyl  eclipse (anti-sun) -> 0.0',
     &            rohp_cyl, 0.d0, eps10, npass, nfail)

c     ==============================================================
c     23. shadfncone: satellite clearly on sunlit side -> ROHP = 1.0
c         sat=[8000,0,0] km, sun=[1.5e8,0,0] km
c     ==============================================================
      rohp_cone = shadfncone(8.d3, 0.d0, 0.d0,
     &                       1.5d8, 0.d0, 0.d0)
      call assert('shadfncone sunlit side -> 1.0',
     &            rohp_cone, 1.d0, eps10, npass, nfail)

c     ==============================================================
c     24. shadfncone: satellite in umbra -> ROHP = 0.0
c         sat=[-8000,0,0] km (anti-sun), sun=[1.5e8,0,0] km
c     ==============================================================
      rohp_cone = shadfncone(-8.d3, 0.d0, 0.d0,
     &                        1.5d8, 0.d0, 0.d0)
      call assert('shadfncone umbra (anti-sun) -> 0.0',
     &            rohp_cone, 0.d0, eps10, npass, nfail)

c     ==============================================================
c     Sun azimuth angle tests
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- sun_azimuth ---'

c     SA1: Sun on apsidal line → Λ_S ≈ pi/2
c     Equatorial orbit (i=0), perigee along x (ω=0, Ω=0)
c     Sun along x-axis: αS=0, δS=0
c     l3 = cos(0)cos(0)*1 + ... = 1, l4 = 0 → atan2(1,0) = pi/2
      call sun_azimuth(0.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     &                 res_atan)
      call assert('SA1 Sun on apsidal line = pi/2',
     &            res_atan, pi/2.d0, 1.d-10, npass, nfail)

c     SA2: Sun at 90° from apsidal line → Λ_S ≈ 0 or pi
c     Equatorial orbit, perigee along x, Sun along y: αS=pi/2
      call sun_azimuth(0.d0, 0.d0, 0.d0, pi/2.d0, 0.d0,
     &                 res_atan)
      call assert('SA2 Sun 90deg from apsides = 0',
     &            res_atan, 0.d0, 1.d-10, npass, nfail)

c     SA3: Sun opposite apsidal line → Λ_S ≈ 3*pi/2
c     Equatorial orbit, Sun along -x: αS=pi
      call sun_azimuth(0.d0, 0.d0, 0.d0, pi, 0.d0,
     &                 res_atan)
      call assert('SA3 Sun opposite apsides=3pi/2',
     &            res_atan, 1.5d0*pi, 1.d-10, npass, nfail)

c     SA4: Rotating RAAN by 90° shifts azimuth
c     Equatorial orbit, ω=0, Ω=pi/2, Sun at αS=pi/2
c     Perigee now along y → Sun on apsidal line → pi/2
      call sun_azimuth(0.d0, 0.d0, pi/2.d0, pi/2.d0, 0.d0,
     &                 res_atan)
      call assert('SA4 RAAN=90 shift = pi/2',
     &            res_atan, pi/2.d0, 1.d-10, npass, nfail)

c     SA5: Non-zero inclination, Sun in equatorial plane
c     i=45°, ω=0, Ω=0, Sun along x: αS=0, δS=0
c     l11=1, l12=0, l13=0 → l3=1
c     l21=0, l22=cos(45), l23=0 → l4=0
c     → atan2(1,0) = pi/2
      call sun_azimuth(45.d0*d2r, 0.d0, 0.d0, 0.d0, 0.d0,
     &                 res_atan)
      call assert('SA5 i=45 Sun on apsides=pi/2',
     &            res_atan, pi/2.d0, 1.d-10, npass, nfail)

c     SA6: Polar orbit (i=90°), ω=0, Ω=0, Sun along z
c     l13=sin(0)*sin(90)=0, l23=cos(0)*sin(90)=1
c     Sun along z: αS=0, δS=pi/2 → els=(0,0,1)
c     l3=l13*1=0, l4=l23*1=1 → atan2(0,1)=0
      call sun_azimuth(pi/2.d0, 0.d0, 0.d0, 0.d0, pi/2.d0,
     &                 res_atan)
      call assert('SA6 polar i=90 Sun@pole = 0',
     &            res_atan, 0.d0, 1.d-10, npass, nfail)

c     SA7: Output always in [0, 2*pi)
c     Various inputs should all give non-negative result
      call sun_azimuth(30.d0*d2r, 45.d0*d2r, 60.d0*d2r,
     &                 3.5d0, -0.3d0, res_atan)
      if (res_atan .ge. 0.d0 .and. res_atan .lt. 2.d0*pi) then
         write(*,*) '  [PASS]  SA7 azimuth in [0, 2pi)'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA7 azimuth in [0, 2pi)'
         nfail = nfail + 1
      end if

c     SA8: Symmetry — opposite Sun gives opposite quadrant
      call sun_azimuth(0.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     &                 res_atan)
      call sun_azimuth(0.d0, 0.d0, 0.d0, pi, 0.d0,
     &                 hmag)
      call assert('SA8 opposite Sun = pi apart',
     &            dabs(hmag - res_atan), pi, 1.d-10,
     &            npass, nfail)

c     SA9: ω rotation shifts azimuth by same amount
c     Equatorial orbit, Sun along x, ω=0 → pi/2
c     ω=pi/4 should shift by pi/4
      call sun_azimuth(0.d0, 0.d0, 0.d0, 0.d0, 0.d0,
     &                 res_atan)
      call sun_azimuth(0.d0, pi/4.d0, 0.d0, 0.d0, 0.d0,
     &                 hmag)
      call assert('SA9 omega shift = pi/4',
     &            dabs(hmag - res_atan), pi/4.d0, 1.d-10,
     &            npass, nfail)

c     SA10: Full azimuth sweep — all 4 quadrants reachable
c     Equatorial orbit, vary αS from 0 to 2pi
      iall_ok = 1
      iq1 = 0
      iq2 = 0
      iq3 = 0
      iq4 = 0
      do i = 0, 35
         call sun_azimuth(0.d0, 0.d0, 0.d0,
     &        dble(i)*10.d0*d2r, 0.d0, res_atan)
         if (res_atan .lt. 0.d0 .or. res_atan .ge. 2.d0*pi)
     &      iall_ok = 0
         if (res_atan .ge. 0.d0 .and.
     &       res_atan .lt. pi/2.d0) iq1 = 1
         if (res_atan .ge. pi/2.d0 .and.
     &       res_atan .lt. pi) iq2 = 1
         if (res_atan .ge. pi .and.
     &       res_atan .lt. 1.5d0*pi) iq3 = 1
         if (res_atan .ge. 1.5d0*pi .and.
     &       res_atan .lt. 2.d0*pi) iq4 = 1
      end do
      if (iall_ok .eq. 1 .and. iq1+iq2+iq3+iq4 .eq. 4) then
         write(*,*) '  [PASS]  SA10 all 4 quadrants reached'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA10 all 4 quadrants reached'
         nfail = nfail + 1
      end if

c     SA11: Realistic GEO orbit (i=0.05°, e≈0, Ω=75°, ω=0°)
c     Sun at vernal equinox: αS≈0, δS≈0
c     Result should be finite and in [0, 2pi)
      call sun_azimuth(0.05d0*d2r, 0.d0, 75.d0*d2r,
     &                 0.d0, 0.d0, res_atan)
      if (res_atan .ge. 0.d0 .and. res_atan .lt. 2.d0*pi) then
         write(*,*) '  [PASS]  SA11 GEO orbit valid range'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA11 GEO orbit valid range'
         nfail = nfail + 1
      end if

c     SA12: GEO — Sun at summer solstice (αS≈90°, δS≈23.4°)
c     Result should differ from equinox (SA11)
      call sun_azimuth(0.05d0*d2r, 0.d0, 75.d0*d2r,
     &                 90.d0*d2r, 23.4d0*d2r, hmag)
      if (dabs(hmag - res_atan) .gt. 0.1d0) then
         write(*,*) '  [PASS]  SA12 GEO solstice differs'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA12 GEO solstice differs'
         nfail = nfail + 1
      end if

c     SA13: Realistic GTO (i=7°, e=0.73, ω=178°, Ω=240°)
c     Ariane 5 ESC-A (case A from manuscript Table 4)
c     Sun near winter solstice: αS≈270°, δS≈-23.4°
      call sun_azimuth(7.d0*d2r, 178.d0*d2r, 240.d0*d2r,
     &                 270.d0*d2r, -23.4d0*d2r, res_atan)
      if (res_atan .ge. 0.d0 .and. res_atan .lt. 2.d0*pi) then
         write(*,*) '  [PASS]  SA13 GTO Ariane valid range'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA13 GTO Ariane valid range'
         nfail = nfail + 1
      end if

c     SA14: GTO — varying Sun RA sweeps all quadrants
c     i=7°, ω=320°, Ω=240° (Ariane GTO config)
      iq1 = 0
      iq2 = 0
      iq3 = 0
      iq4 = 0
      do i = 0, 35
         call sun_azimuth(7.d0*d2r, 320.d0*d2r, 240.d0*d2r,
     &        dble(i)*10.d0*d2r, 0.d0, res_atan)
         if (res_atan .ge. 0.d0 .and.
     &       res_atan .lt. pi/2.d0) iq1 = 1
         if (res_atan .ge. pi/2.d0 .and.
     &       res_atan .lt. pi) iq2 = 1
         if (res_atan .ge. pi .and.
     &       res_atan .lt. 1.5d0*pi) iq3 = 1
         if (res_atan .ge. 1.5d0*pi .and.
     &       res_atan .lt. 2.d0*pi) iq4 = 1
      end do
      if (iq1+iq2+iq3+iq4 .eq. 4) then
         write(*,*) '  [PASS]  SA14 GTO all quadrants'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA14 GTO all quadrants'
         nfail = nfail + 1
      end if

c     SA15: Molniya orbit (i=63.4°, e=0.74, ω=270°, Ω=0°)
c     Sun at αS=180°, δS=0° (autumn equinox)
      call sun_azimuth(63.4d0*d2r, 270.d0*d2r, 0.d0,
     &                 180.d0*d2r, 0.d0, res_atan)
      if (res_atan .ge. 0.d0 .and. res_atan .lt. 2.d0*pi) then
         write(*,*) '  [PASS]  SA15 Molniya valid range'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA15 Molniya valid range'
         nfail = nfail + 1
      end if

c     SA16: GTO with high inclination (PSLV, i=19.2°)
c     ω=128°, Ω=166°, Sun at αS=160°, δS=12°
      call sun_azimuth(19.2d0*d2r, 128.d0*d2r, 166.d0*d2r,
     &                 160.d0*d2r, 12.d0*d2r, res_atan)
      if (res_atan .ge. 0.d0 .and. res_atan .lt. 2.d0*pi) then
         write(*,*) '  [PASS]  SA16 PSLV GTO valid range'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA16 PSLV GTO valid range'
         nfail = nfail + 1
      end if

c     SA17: Proton-M Briz-M (i=48.7°, ω=167°, Ω=1.4°)
c     Sun at αS=95°, δS=23°
      call sun_azimuth(48.7d0*d2r, 167.d0*d2r, 1.4d0*d2r,
     &                 95.d0*d2r, 23.d0*d2r, res_atan)
      if (res_atan .ge. 0.d0 .and. res_atan .lt. 2.d0*pi) then
         write(*,*) '  [PASS]  SA17 Briz-M HEO valid range'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA17 Briz-M HEO valid range'
         nfail = nfail + 1
      end if

c     SA18: Seasonal variation requires meaningful inclination
c     i=30°, ω=0, Ω=75°: winter vs summer δS
      call sun_azimuth(30.d0*d2r, 0.d0, 75.d0*d2r,
     &                 90.d0*d2r, -23.4d0*d2r, res_atan)
      call sun_azimuth(30.d0*d2r, 0.d0, 75.d0*d2r,
     &                 90.d0*d2r, 23.4d0*d2r, hmag)
      if (dabs(hmag - res_atan) .gt. 0.01d0) then
         write(*,*) '  [PASS]  SA18 seasonal shift i=30'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA18 seasonal shift i=30'
         nfail = nfail + 1
      end if

c     SA19: GTO RAAN sweep — azimuth varies smoothly
c     i=7°, ω=320°, Sun fixed at equinox
      iall_ok = 1
      do i = 0, 35
         call sun_azimuth(7.d0*d2r, 320.d0*d2r,
     &        dble(i)*10.d0*d2r, 0.d0, 0.d0, res_atan)
         if (res_atan .lt. 0.d0 .or. res_atan .ge. 2.d0*pi)
     &      iall_ok = 0
      end do
      if (iall_ok .eq. 1) then
         write(*,*) '  [PASS]  SA19 GTO RAAN sweep valid'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA19 GTO RAAN sweep valid'
         nfail = nfail + 1
      end if

c     SA20: Manuscript Case A (Ariane ESC-A)
c     i=6.1°, ω=320°, Ω=240°, αS≈280°, δS≈-23°
c     Just verify finite result in [0, 2pi)
      call sun_azimuth(6.1d0*d2r, 320.d0*d2r, 240.d0*d2r,
     &                 280.d0*d2r, -23.d0*d2r, res_atan)
      if (res_atan .ge. 0.d0 .and. res_atan .lt. 2.d0*pi
     &    .and. res_atan .eq. res_atan) then
         write(*,*) '  [PASS]  SA20 Case A finite azimuth'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  SA20 Case A finite azimuth'
         nfail = nfail + 1
      end if

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
