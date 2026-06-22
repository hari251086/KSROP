c     ============================================================
c     test_tle2sv.F — Tests for TLE-to-state-vector conversion
c
c     Coverage:
c       tle_epoch2jd : known epoch dates, leap year, year boundaries
c       tle_sgp4_init: near-Earth init, deep-space rejection,
c                      perigee regimes, constant recovery
c       tle_sgp4_prop: zero-time propagation, short propagation,
c                      Kepler convergence
c       tle_teme2j2k : frame rotation identity check, J2000 epoch,
c                      precession/nutation magnitude
c       tle2sv       : ISS TLE full pipeline, sun-synchronous orbit,
c                      high-eccentricity orbit, Vanguard 1,
c                      near-circular LEO, consistency checks,
c                      roundtrip orbital elements
c     ============================================================

      program test_tle2sv
      implicit double precision (a-h, o-z)

      integer npass, nfail, ierr, iyr4, iall_ok, imax
      integer iresfl
      double precision r_j2k(3), v_j2k(3)
      double precision r_teme(3), v_teme(3)
      double precision r_out(3), v_out(3)
      double precision n0dd, a0dd, cosio, sinio
      double precision eta, c1_coeff, c4_coeff, t0, perigee_km
      double precision omgdot, xnodot, xmdot
      double precision djd, djd2, tol
      double precision rmag, vmag, hmag, energy
      double precision pi, d2r, twopi
      double precision dsstate(50)

c     Batch test arrays
      integer i_norad(3), i_yr(3), i_rev(3), i_ntle
      double precision i_eday(3), i_ndot(3), i_nddot(3), i_bstar(3)
      double precision i_inc(3), i_raan(3), i_ecc(3)
      double precision i_aop(3), i_ma(3), i_mm(3)
      character*24 i_names(3)
      character*5  i_catid(3)

      pi = 3.14159265358979324d0
      twopi = 2.d0 * pi
      d2r = pi / 180.d0

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' TLE-to-State-Vector Test Suite'
      write(*,*) '=============================================='

c     ============================================================
c     tle_epoch2jd: epoch-to-Julian-Date conversion
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle_epoch2jd ---'

c     T1: J2000.0 epoch = JD 2451545.0 (2000 Jan 1.5)
c     Year 2000, day 1.5 -> JD 2451545.0
      call tle_epoch2jd(2000, 1.5d0, djd)
      call assert_d('J2000 epoch JD',
     &              djd, 2451545.d0, 0.5d0, npass, nfail)

c     T2: 2024 Jan 1.0 = JD 2460310.5
      call tle_epoch2jd(2024, 1.0d0, djd)
      call assert_d('2024 Jan 1.0 JD',
     &              djd, 2460310.5d0, 0.5d0, npass, nfail)

c     T3: 2021 day 83.05205046 (47944 first TLE epoch)
      call tle_epoch2jd(2021, 83.05205046d0, djd)
c     2021 Jan 0.0 = JD 2459215.5, + 83.052 = 2459298.552
      call assert_d('2021 day 83 JD',
     &              djd, 2459298.55205046d0, 1.d0, npass, nfail)

c     T4: Leap year 2024, day 366
      call tle_epoch2jd(2024, 366.d0, djd)
      call tle_epoch2jd(2024, 365.d0, djd2)
      call assert_d('leap year day diff',
     &              djd - djd2, 1.d0, 1.d-10, npass, nfail)

c     T5: Year boundary: consecutive years have JD diff ~365
      call tle_epoch2jd(2025, 1.d0, djd)
      call tle_epoch2jd(2024, 1.d0, djd2)
      diff = djd - djd2
      if (dabs(diff - 366.d0) .lt. 1.d0 .or.
     &    dabs(diff - 365.d0) .lt. 1.d0) then
         call assert_i('year boundary JD diff ~365-366',
     &                 1, 1, npass, nfail)
      else
         call assert_i('year boundary JD diff ~365-366',
     &                 0, 1, npass, nfail)
      end if

c     T6: 1957 (Sputnik era)
      call tle_epoch2jd(1957, 1.d0, djd)
      if (djd .gt. 2435000.d0 .and. djd .lt. 2436500.d0) then
         call assert_i('1957 JD range ok', 1, 1, npass, nfail)
      else
         call assert_i('1957 JD range ok', 0, 1, npass, nfail)
      end if

c     T7: Day 1.0 and day 2.0 differ by exactly 1
      call tle_epoch2jd(2026, 1.d0, djd)
      call tle_epoch2jd(2026, 2.d0, djd2)
      call assert_d('day 1 vs day 2 = 1.0',
     &              djd2 - djd, 1.d0, 1.d-12, npass, nfail)

c     ============================================================
c     tle_sgp4_init: initialization tests
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle_sgp4_init ---'

c     ISS-like TLE: n=15.54 rev/day, e=0.00067, i=51.64
      call tle_epoch2jd(2026, 170.5d0, djd)
      call tle_sgp4_init(15.54d0, 0.00067d0, 51.64d0,
     &                   123.4567d0, 45.1234d0, 315.d0,
     &                   1.027d-4, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot, ierr)

c     T8: Near-Earth should succeed
      call assert_i('ISS init ierr=0', ierr, 0, npass, nfail)

c     T9: Recovered mean motion should be close to input
c     Input: 15.54 rev/day = 15.54*2pi/1440 rad/min ~ 0.0678 rad/min
      tol = 1.d-4
      call assert_d('ISS n0dd ~ 0.0678 rad/min',
     &              n0dd, 15.54d0*twopi/1440.d0, tol,
     &              npass, nfail)

c     T10: Semi-major axis ~ 6780 km / 6378.135 ~ 1.063 Earth radii
      if (a0dd .gt. 1.01d0 .and. a0dd .lt. 1.1d0) then
         call assert_i('ISS a0dd in [1.01, 1.1] Re',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS a0dd in [1.01, 1.1] Re',
     &                 0, 1, npass, nfail)
      end if

c     T11: cosio = cos(51.64 deg)
      call assert_d('ISS cosio',
     &              cosio, dcos(51.64d0*d2r), 1.d-10,
     &              npass, nfail)

c     T12: sinio = sin(51.64 deg)
      call assert_d('ISS sinio',
     &              sinio, dsin(51.64d0*d2r), 1.d-10,
     &              npass, nfail)

c     T13: Perigee altitude > 0 for ISS
      if (perigee_km .gt. 100.d0 .and.
     &    perigee_km .lt. 500.d0) then
         call assert_i('ISS perigee in [100,500] km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS perigee in [100,500] km',
     &                 0, 1, npass, nfail)
      end if

c     T14: Secular rates should be nonzero
      if (dabs(omgdot) .gt. 0.d0 .and.
     &    dabs(xnodot) .gt. 0.d0 .and.
     &    dabs(xmdot) .gt. 0.d0) then
         call assert_i('ISS secular rates nonzero',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS secular rates nonzero',
     &                 0, 1, npass, nfail)
      end if

c     T15: xmdot should be close to n0dd (mean motion dominates)
      call assert_d('ISS xmdot ~ n0dd',
     &              xmdot, n0dd, 1.d-4, npass, nfail)

c     T16: Deep-space orbit acceptance via SDP4
c     Geostationary: n ~ 1.0 rev/day -> period ~ 1440 min > 225
      call tle_sdp4_init(1.0d0, 0.001d0, 0.1d0,
     &                   0.d0, 0.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call assert_i('GEO SDP4 init ierr=0',
     &              ierr, 0, npass, nfail)
      if (n0dd .gt. 0.d0) then
         call assert_i('GEO n0dd > 0', 1, 1, npass, nfail)
      else
         call assert_i('GEO n0dd > 0', 0, 1, npass, nfail)
      end if
      if (a0dd .gt. 6.d0) then
         call assert_i('GEO a0dd > 6 Re', 1, 1, npass, nfail)
      else
         call assert_i('GEO a0dd > 6 Re', 0, 1, npass, nfail)
      end if
      if (perigee_km .gt. 30000.d0) then
         call assert_i('GEO perigee > 30000 km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('GEO perigee > 30000 km',
     &                 0, 1, npass, nfail)
      end if

c     T17: Molniya-like deep space via SDP4 (n ~ 2.006 rev/day)
      call tle_sdp4_init(2.006d0, 0.74d0, 63.4d0,
     &                   0.d0, 270.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call assert_i('Molniya SDP4 init ierr=0',
     &              ierr, 0, npass, nfail)
      if (n0dd .gt. 0.d0 .and. a0dd .gt. 1.d0) then
         call assert_i('Molniya elements physical',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Molniya elements physical',
     &                 0, 1, npass, nfail)
      end if

c     T18: Low perigee regime (very low orbit)
      call tle_sgp4_init(16.4d0, 0.001d0, 51.6d0,
     &                   0.d0, 0.d0, 0.d0,
     &                   1.d-3, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot, ierr)
      call assert_i('low-orbit init ierr=0',
     &              ierr, 0, npass, nfail)

c     T19: Zero eccentricity
      call tle_sgp4_init(15.5d0, 0.d0, 51.6d0,
     &                   0.d0, 0.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot, ierr)
      call assert_i('zero-ecc init ierr=0',
     &              ierr, 0, npass, nfail)

c     T20: Zero inclination (equatorial)
      call tle_sgp4_init(15.5d0, 0.001d0, 0.d0,
     &                   0.d0, 0.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot, ierr)
      call assert_i('equatorial init ierr=0',
     &              ierr, 0, npass, nfail)
      call assert_d('equatorial cosio = 1',
     &              cosio, 1.d0, 1.d-12, npass, nfail)
      call assert_d('equatorial sinio = 0',
     &              sinio, 0.d0, 1.d-12, npass, nfail)

c     T21: 90 deg inclination (polar)
      call tle_sgp4_init(15.5d0, 0.001d0, 90.d0,
     &                   0.d0, 0.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot, ierr)
      call assert_i('polar init ierr=0', ierr, 0, npass, nfail)
      call assert_d('polar cosio = 0',
     &              cosio, 0.d0, 1.d-12, npass, nfail)

c     ============================================================
c     tle_sgp4_prop: propagation tests
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle_sgp4_prop ---'

c     Re-initialize ISS for propagation tests
      call tle_sgp4_init(15.54d0, 0.00067d0, 51.64d0,
     &                   123.4567d0, 45.1234d0, 315.d0,
     &                   1.027d-4, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot, ierr)

c     T22: Propagate at t=0 -> should produce valid state
      call tle_sgp4_prop(0.d0,
     &              0.00067d0, 51.64d0, 123.4567d0,
     &              45.1234d0, 315.d0,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1_coeff, c4_coeff, perigee_km,
     &              omgdot, xnodot, xmdot, 1.027d-4,
     &              r_teme, v_teme, ierr)
      call assert_i('ISS prop t=0 ierr=0',
     &              ierr, 0, npass, nfail)

c     T23: Position magnitude ~ 6378+400 = 6778 km
      rmag = dsqrt(r_teme(1)**2 + r_teme(2)**2 + r_teme(3)**2)
      if (rmag .gt. 6400.d0 .and. rmag .lt. 7200.d0) then
         call assert_i('ISS rmag in [6400,7200] km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS rmag in [6400,7200] km',
     &                 0, 1, npass, nfail)
      end if

c     T24: Velocity magnitude ~ 7.5 km/s for LEO
      vmag = dsqrt(v_teme(1)**2 + v_teme(2)**2 + v_teme(3)**2)
      if (vmag .gt. 7.0d0 .and. vmag .lt. 8.5d0) then
         call assert_i('ISS vmag in [7.0,8.5] km/s',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS vmag in [7.0,8.5] km/s',
     &                 0, 1, npass, nfail)
      end if

c     T25: z-component bounded by inclination
c     |z| <= rmag * sin(inc)
      zmax = rmag * dsin(51.64d0*d2r)
      if (dabs(r_teme(3)) .le. zmax * 1.01d0) then
         call assert_i('ISS z bounded by inc',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS z bounded by inc',
     &                 0, 1, npass, nfail)
      end if

c     T26: Orbital energy should be negative (bound orbit)
      energy = 0.5d0*vmag*vmag - 398600.8d0/rmag
      if (energy .lt. 0.d0) then
         call assert_i('ISS energy < 0 (bound)',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS energy < 0 (bound)',
     &                 0, 1, npass, nfail)
      end if

c     T27: Short propagation (1 minute) should give close position
      call tle_sgp4_prop(1.d0,
     &              0.00067d0, 51.64d0, 123.4567d0,
     &              45.1234d0, 315.d0,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1_coeff, c4_coeff, perigee_km,
     &              omgdot, xnodot, xmdot, 1.027d-4,
     &              r_out, v_out, ierr)
      call assert_i('ISS prop t=1min ierr=0',
     &              ierr, 0, npass, nfail)

c     T28: Position change in 1 minute ~ v*60 ~ 450 km
      dr = dsqrt((r_out(1)-r_teme(1))**2
     &         + (r_out(2)-r_teme(2))**2
     &         + (r_out(3)-r_teme(3))**2)
      call assert_d('ISS 1min dr ~ 460 km',
     &              dr, 460.d0, 360.d0, npass, nfail)

c     T29: One full orbit ~ 92 min, should return near start
      period_min = twopi / n0dd
      call tle_sgp4_prop(period_min,
     &              0.00067d0, 51.64d0, 123.4567d0,
     &              45.1234d0, 315.d0,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1_coeff, c4_coeff, perigee_km,
     &              omgdot, xnodot, xmdot, 1.027d-4,
     &              r_out, v_out, ierr)
      call assert_i('ISS full-orbit ierr=0',
     &              ierr, 0, npass, nfail)

c     T30: After 1 orbit, rmag should be similar
      rmag2 = dsqrt(r_out(1)**2 + r_out(2)**2 + r_out(3)**2)
      call assert_d('ISS rmag after 1 orbit',
     &              rmag2, rmag, 100.d0, npass, nfail)

c     T31: Near-circular orbit velocity magnitude nearly constant
      vmag2 = dsqrt(v_out(1)**2 + v_out(2)**2 + v_out(3)**2)
      call assert_d('ISS vmag after 1 orbit',
     &              vmag2, vmag, 0.5d0, npass, nfail)

c     ============================================================
c     tle_teme2j2k: frame transformation tests
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle_teme2j2k ---'

c     T32: At J2000.0 epoch, TEME ~ J2000 (small rotation)
      r_teme(1) = 6778.d0
      r_teme(2) = 0.d0
      r_teme(3) = 0.d0
      v_teme(1) = 0.d0
      v_teme(2) = 7.5d0
      v_teme(3) = 0.d0
      call tle_teme2j2k(r_teme, v_teme, 2451545.d0,
     &                  r_j2k, v_j2k)
c     At J2000 epoch, precession is zero, nutation is small
      rmag = dsqrt(r_j2k(1)**2 + r_j2k(2)**2 + r_j2k(3)**2)
      call assert_d('J2000 epoch rmag preserved',
     &              rmag, 6778.d0, 1.d0, npass, nfail)

c     T33: Velocity magnitude preserved by rotation
      vmag = dsqrt(v_j2k(1)**2 + v_j2k(2)**2 + v_j2k(3)**2)
      call assert_d('J2000 epoch vmag preserved',
     &              vmag, 7.5d0, 0.01d0, npass, nfail)

c     T34: At J2000.0, rotation should be nearly identity
c     (only nutation, which is small ~17 arcsec max)
      call assert_d('J2000 r_x ~ 6778',
     &              r_j2k(1), 6778.d0, 5.d0, npass, nfail)
      call assert_d('J2000 r_y ~ 0',
     &              dabs(r_j2k(2)), 0.d0, 5.d0, npass, nfail)

c     T35: 2026 epoch - rotation still preserves magnitudes
      r_teme(1) = 4000.d0
      r_teme(2) = 5000.d0
      r_teme(3) = 3000.d0
      v_teme(1) = -3.d0
      v_teme(2) = 4.d0
      v_teme(3) = 5.d0
      call tle_epoch2jd(2026, 170.5d0, djd)
      call tle_teme2j2k(r_teme, v_teme, djd, r_j2k, v_j2k)

      rmag = dsqrt(r_teme(1)**2+r_teme(2)**2+r_teme(3)**2)
      rmag2 = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      call assert_d('2026 rmag preserved',
     &              rmag2, rmag, 0.01d0, npass, nfail)

      vmag = dsqrt(v_teme(1)**2+v_teme(2)**2+v_teme(3)**2)
      vmag2 = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
      call assert_d('2026 vmag preserved',
     &              vmag2, vmag, 1.d-6, npass, nfail)

c     T36: Rotation angle should be small (< 1 deg ~ 0.017 rad)
c     Dot product of unit vectors should be > cos(1 deg)
      dp = (r_teme(1)*r_j2k(1)+r_teme(2)*r_j2k(2)
     &     +r_teme(3)*r_j2k(3)) / (rmag*rmag2)
      if (dp .gt. dcos(1.d0*d2r)) then
         call assert_i('2026 rotation < 1 deg',
     &                 1, 1, npass, nfail)
      else
         call assert_i('2026 rotation < 1 deg',
     &                 0, 1, npass, nfail)
      end if

c     T37: Orthogonality: r dot v unchanged by rotation
      rv1 = r_teme(1)*v_teme(1)+r_teme(2)*v_teme(2)
     &     +r_teme(3)*v_teme(3)
      rv2 = r_j2k(1)*v_j2k(1)+r_j2k(2)*v_j2k(2)
     &     +r_j2k(3)*v_j2k(3)
      call assert_d('r.v preserved by rotation',
     &              rv2, rv1, 1.d-6, npass, nfail)

c     ============================================================
c     tle2sv: full pipeline tests
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle2sv: full pipeline ---'

c     T38-41: ISS TLE
c     1 25544U 98067A   26170.50000000  .00016717  00000-0  10270-3 0  9993
c     2 25544  51.6400 123.4567 0006700  45.1234 315.0000 15.54000000456789
      call tle2sv(2026, 170.5d0, 51.64d0, 123.4567d0, 0.00067d0,
     &            45.1234d0, 315.d0, 15.54d0, 1.027d-4,
     &            r_j2k, v_j2k, ierr)
      call assert_i('ISS tle2sv ierr=0', ierr, 0, npass, nfail)

      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      call assert_d('ISS rmag ~ 6778 km',
     &              rmag, 6778.d0, 200.d0, npass, nfail)

      vmag = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
      if (vmag .gt. 7.0d0 .and. vmag .lt. 8.5d0) then
         call assert_i('ISS vmag in [7.0,8.5] km/s',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS vmag in [7.0,8.5] km/s',
     &                 0, 1, npass, nfail)
      end if

c     T41: |z| bounded by inclination
      zmax = rmag * dsin(51.64d0*d2r)
      if (dabs(r_j2k(3)) .le. zmax * 1.01d0) then
         call assert_i('ISS J2K z bounded by inc',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS J2K z bounded by inc',
     &                 0, 1, npass, nfail)
      end if

c     T42-44: Sun-synchronous orbit (NOAA-18 like)
c     i=99.0357, e=0.0014319, n=14.12818068
      call tle2sv(2026, 170.52919894d0, 99.0357d0, 195.5267d0,
     &            0.0014319d0, 94.5046d0, 265.7602d0,
     &            14.12818068d0, 3.2076d-5,
     &            r_j2k, v_j2k, ierr)
      call assert_i('NOAA18 tle2sv ierr=0', ierr, 0, npass, nfail)

      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
c     NOAA-18 altitude ~850 km -> rmag ~ 7228 km
      if (rmag .gt. 6800.d0 .and. rmag .lt. 7600.d0) then
         call assert_i('NOAA18 rmag in [6800,7600] km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('NOAA18 rmag in [6800,7600] km',
     &                 0, 1, npass, nfail)
      end if

      vmag = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
      if (vmag .gt. 7.0d0 .and. vmag .lt. 8.0d0) then
         call assert_i('NOAA18 vmag in [7.0,8.0] km/s',
     &                 1, 1, npass, nfail)
      else
         call assert_i('NOAA18 vmag in [7.0,8.0] km/s',
     &                 0, 1, npass, nfail)
      end if

c     T45-48: Vanguard 1 (high eccentricity LEO, e=0.1835576)
      call tle2sv(2026, 170.40233228d0, 34.2451d0, 307.4935d0,
     &            0.1835576d0, 148.4938d0, 223.9981d0,
     &            10.86000951d0, 3.7765d-4,
     &            r_j2k, v_j2k, ierr)
      call assert_i('Vanguard tle2sv ierr=0',
     &              ierr, 0, npass, nfail)

      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
c     Vanguard 1: perigee ~650 km, apogee ~3800 km
c     rmag should be somewhere in [7000, 10200]
      if (rmag .gt. 6500.d0 .and. rmag .lt. 11000.d0) then
         call assert_i('Vanguard rmag in [6500,11000] km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Vanguard rmag in [6500,11000] km',
     &                 0, 1, npass, nfail)
      end if

      vmag = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
      if (vmag .gt. 4.d0 .and. vmag .lt. 10.d0) then
         call assert_i('Vanguard vmag in [4,10] km/s',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Vanguard vmag in [4,10] km/s',
     &                 0, 1, npass, nfail)
      end if

c     T49: Negative energy (bound orbit)
      energy = 0.5d0*vmag*vmag - 398600.8d0/rmag
      if (energy .lt. 0.d0) then
         call assert_i('Vanguard energy < 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Vanguard energy < 0',
     &                 0, 1, npass, nfail)
      end if

c     T50: Deep-space TLE (GEO) via SDP4 - now succeeds
      call tle2sv(2026, 170.d0, 0.1d0, 0.d0,
     &            0.001d0, 0.d0, 0.d0, 1.0d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('GEO tle2sv ierr=0', ierr, 0, npass, nfail)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      vmag = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
      if (rmag .gt. 41000.d0 .and. rmag .lt. 43000.d0) then
         call assert_i('GEO rmag in [41k,43k] km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('GEO rmag in [41k,43k] km',
     &                 0, 1, npass, nfail)
      end if
      if (vmag .gt. 2.8d0 .and. vmag .lt. 3.3d0) then
         call assert_i('GEO vmag in [2.8,3.3] km/s',
     &                 1, 1, npass, nfail)
      else
         call assert_i('GEO vmag in [2.8,3.3] km/s',
     &                 0, 1, npass, nfail)
      end if

c     ============================================================
c     Consistency and regression tests
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- consistency tests ---'

c     T51: Angular momentum conservation check
c     h = r x v, |h| should match expected for orbit
      call tle2sv(2026, 170.5d0, 51.64d0, 123.4567d0, 0.00067d0,
     &            45.1234d0, 315.d0, 15.54d0, 1.027d-4,
     &            r_j2k, v_j2k, ierr)

      hx = r_j2k(2)*v_j2k(3) - r_j2k(3)*v_j2k(2)
      hy = r_j2k(3)*v_j2k(1) - r_j2k(1)*v_j2k(3)
      hz = r_j2k(1)*v_j2k(2) - r_j2k(2)*v_j2k(1)
      hmag = dsqrt(hx*hx + hy*hy + hz*hz)

c     For circular orbit: h = r*v ~ 6778*7.5 ~ 50835
      if (hmag .gt. 40000.d0 .and. hmag .lt. 60000.d0) then
         call assert_i('ISS |h| in [40k,60k]',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ISS |h| in [40k,60k]',
     &                 0, 1, npass, nfail)
      end if

c     T52: Semi-major axis from energy
c     a = -mu/(2*E)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      vmag = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
      energy = 0.5d0*vmag*vmag - 398600.8d0/rmag
      sma = -398600.8d0 / (2.d0*energy)
c     ISS SMA ~ 6778 km
      call assert_d('ISS SMA from energy',
     &              sma, 6778.d0, 200.d0, npass, nfail)

c     T53: Orbital period from SMA
c     T = 2*pi*sqrt(a^3/mu)
      per_s = twopi * dsqrt(sma**3 / 398600.8d0)
      per_min = per_s / 60.d0
c     ISS period ~ 92 min (from 15.54 rev/day -> 1440/15.54 = 92.7)
      call assert_d('ISS period from SMA (min)',
     &              per_min, 1440.d0/15.54d0, 5.d0,
     &              npass, nfail)

c     T54: Same TLE called twice gives identical results
      call tle2sv(2026, 170.5d0, 51.64d0, 123.4567d0, 0.00067d0,
     &            45.1234d0, 315.d0, 15.54d0, 1.027d-4,
     &            r_out, v_out, ierr)
      dr = dsqrt((r_j2k(1)-r_out(1))**2 + (r_j2k(2)-r_out(2))**2
     &         + (r_j2k(3)-r_out(3))**2)
      call assert_d('repeatability dr=0',
     &              dr, 0.d0, 1.d-10, npass, nfail)

c     T55: Different mean anomaly -> different position
      call tle2sv(2026, 170.5d0, 51.64d0, 123.4567d0, 0.00067d0,
     &            45.1234d0, 180.d0, 15.54d0, 1.027d-4,
     &            r_out, v_out, ierr)
      dr = dsqrt((r_j2k(1)-r_out(1))**2 + (r_j2k(2)-r_out(2))**2
     &         + (r_j2k(3)-r_out(3))**2)
      if (dr .gt. 100.d0) then
         call assert_i('diff MA -> diff pos (dr>100km)',
     &                 1, 1, npass, nfail)
      else
         call assert_i('diff MA -> diff pos (dr>100km)',
     &                 0, 1, npass, nfail)
      end if

c     T56: MA=0 -> near perigee
      call tle2sv(2026, 170.5d0, 51.64d0, 123.4567d0, 0.00067d0,
     &            45.1234d0, 0.d0, 15.54d0, 1.027d-4,
     &            r_out, v_out, ierr)
      rmag_peri = dsqrt(r_out(1)**2+r_out(2)**2+r_out(3)**2)

c     T57: MA=180 -> near apogee
      call tle2sv(2026, 170.5d0, 51.64d0, 123.4567d0, 0.00067d0,
     &            45.1234d0, 180.d0, 15.54d0, 1.027d-4,
     &            r_out, v_out, ierr)
      rmag_apo = dsqrt(r_out(1)**2+r_out(2)**2+r_out(3)**2)

c     For e>0, apogee rmag > perigee rmag
      if (rmag_apo .ge. rmag_peri - 1.d0) then
         call assert_i('MA=180 rmag >= MA=0 rmag',
     &                 1, 1, npass, nfail)
      else
         call assert_i('MA=180 rmag >= MA=0 rmag',
     &                 0, 1, npass, nfail)
      end if

c     T58: High eccentricity: perigee/apogee ratio check
      call tle2sv(2026, 170.d0, 34.2451d0, 307.d0,
     &            0.1835576d0, 148.d0, 0.d0,
     &            10.86d0, 3.7765d-4,
     &            r_out, v_out, ierr)
      rmag_peri = dsqrt(r_out(1)**2+r_out(2)**2+r_out(3)**2)

      call tle2sv(2026, 170.d0, 34.2451d0, 307.d0,
     &            0.1835576d0, 148.d0, 180.d0,
     &            10.86d0, 3.7765d-4,
     &            r_out, v_out, ierr)
      rmag_apo = dsqrt(r_out(1)**2+r_out(2)**2+r_out(3)**2)

c     Ratio should be ~ (1+e)/(1-e) ~ 1.449
      ratio = rmag_apo / rmag_peri
      call assert_d('high-ecc apo/peri ratio',
     &              ratio, (1.d0+0.1835576d0)/(1.d0-0.1835576d0),
     &              0.1d0, npass, nfail)

c     ============================================================
c     Edge cases
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- edge cases ---'

c     T59: Zero B* (no drag)
      call tle2sv(2026, 170.5d0, 51.64d0, 123.4567d0, 0.00067d0,
     &            45.1234d0, 315.d0, 15.54d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('zero-bstar ierr=0', ierr, 0, npass, nfail)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      if (rmag .gt. 6400.d0 .and. rmag .lt. 7200.d0) then
         call assert_i('zero-bstar rmag ok',
     &                 1, 1, npass, nfail)
      else
         call assert_i('zero-bstar rmag ok',
     &                 0, 1, npass, nfail)
      end if

c     T60: Very small eccentricity
      call tle2sv(2026, 170.5d0, 51.64d0, 123.d0, 1.d-6,
     &            0.d0, 90.d0, 15.54d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('tiny-ecc ierr=0', ierr, 0, npass, nfail)

c     T61: Near-zero inclination
      call tle2sv(2026, 170.5d0, 0.01d0, 0.d0, 0.001d0,
     &            0.d0, 0.d0, 15.54d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('near-zero-inc ierr=0', ierr, 0, npass, nfail)
c     z-component should be very small
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      if (dabs(r_j2k(3))/rmag .lt. 0.01d0) then
         call assert_i('near-zero-inc z/r < 0.01',
     &                 1, 1, npass, nfail)
      else
         call assert_i('near-zero-inc z/r < 0.01',
     &                 0, 1, npass, nfail)
      end if

c     T62: Retrograde orbit (i=180 deg)
      call tle2sv(2026, 170.5d0, 179.9d0, 0.d0, 0.001d0,
     &            0.d0, 0.d0, 15.54d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('retrograde ierr=0', ierr, 0, npass, nfail)

c     T63: MA=0, AOP=0, RAAN=0 (all angles zero)
      call tle2sv(2026, 170.5d0, 51.64d0, 0.d0, 0.001d0,
     &            0.d0, 0.d0, 15.54d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('all-angles-zero ierr=0',
     &              ierr, 0, npass, nfail)

c     T64: MA=360 should equal MA=0
      call tle2sv(2026, 170.5d0, 51.64d0, 0.d0, 0.001d0,
     &            0.d0, 0.d0, 15.54d0, 0.d0,
     &            r_out, v_out, ierr)
      call tle2sv(2026, 170.5d0, 51.64d0, 0.d0, 0.001d0,
     &            0.d0, 360.d0, 15.54d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      dr = dsqrt((r_j2k(1)-r_out(1))**2 + (r_j2k(2)-r_out(2))**2
     &         + (r_j2k(3)-r_out(3))**2)
      call assert_d('MA=360 == MA=0',
     &              dr, 0.d0, 1.d0, npass, nfail)

c     T65: Large B* (decaying orbit)
      call tle2sv(2025, 208.d0, 97.5d0, 346.d0, 0.0017d0,
     &            232.d0, 185.d0, 15.066d0, 3.1149d-4,
     &            r_j2k, v_j2k, ierr)
      call assert_i('large-bstar ierr=0', ierr, 0, npass, nfail)

c     T66: 1957 epoch (earliest possible TLE)
      call tle2sv(1957, 1.d0, 65.d0, 0.d0, 0.01d0,
     &            0.d0, 0.d0, 14.d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('1957 epoch ierr=0', ierr, 0, npass, nfail)

c     T67: Critical inclination (63.4 deg) - singularity check
      call tle2sv(2026, 170.5d0, 63.4d0, 100.d0, 0.01d0,
     &            270.d0, 0.d0, 14.d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('critical-inc ierr=0', ierr, 0, npass, nfail)

c     T68: Moderate eccentricity (0.5)
      call tle2sv(2026, 170.5d0, 28.5d0, 0.d0, 0.5d0,
     &            0.d0, 0.d0, 13.d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('ecc=0.5 ierr=0', ierr, 0, npass, nfail)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
c     At perigee: r ~ a*(1-e)
      if (rmag .gt. 3000.d0) then
         call assert_i('ecc=0.5 rmag > 3000 km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('ecc=0.5 rmag > 3000 km',
     &                 0, 1, npass, nfail)
      end if

c     ============================================================
c     Regression: multiple TLEs via read_tle + batch tle2sv
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- batch TLE conversion ---'

c     T69-72: Read small inline file, convert each to state vector
      imax = 3

      open(unit=90, status='scratch')
c     ISS
      write(90,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(90,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
c     NOAA-18
      write(90,'(A)') '1 28654U 05018A   26170.52919894  .000000'
     &     // '35  00000-0  32076-4 0  9992'
      write(90,'(A)') '2 28654  99.0357 195.5267 0014319  94.504'
     &     // '6 265.7602 14.12818068 77721'
c     Vanguard 1
      write(90,'(A)') '1 00005U 58002B   26170.40233228  .000002'
     &     // '89  00000-0  37765-3 0  9990'
      write(90,'(A)') '2 00005  34.2451 307.4935 1835576 148.493'
     &     // '8 223.9981 10.86000951443521'
      rewind(90)
      call read_tle(90, imax, i_names, i_catid, i_norad,
     &              i_yr, i_eday,
     &              i_ndot, i_nddot, i_bstar,
     &              i_inc, i_raan, i_ecc, i_aop, i_ma, i_mm,
     &              i_rev, i_ntle)
      close(90)

      call assert_i('batch read 3 TLEs', i_ntle, 3, npass, nfail)

c     Convert each and check valid
      iall_ok = 1
      do i = 1, i_ntle
         call tle2sv(i_yr(i), i_eday(i),
     &               i_inc(i), i_raan(i), i_ecc(i),
     &               i_aop(i), i_ma(i), i_mm(i), i_bstar(i),
     &               r_j2k, v_j2k, ierr)
         if (ierr .ne. 0) then
            iall_ok = 0
         else
            rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
            vmag = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
            if (rmag .lt. 6000.d0 .or. rmag .gt. 15000.d0)
     &         iall_ok = 0
            if (vmag .lt. 3.d0 .or. vmag .gt. 11.d0)
     &         iall_ok = 0
            energy = 0.5d0*vmag*vmag - 398600.8d0/rmag
            if (energy .ge. 0.d0) iall_ok = 0
         end if
      end do
      call assert_i('batch all 3 valid states',
     &              iall_ok, 1, npass, nfail)

c     T71: Each TLE gives a different state
      call tle2sv(i_yr(1), i_eday(1),
     &            i_inc(1), i_raan(1), i_ecc(1),
     &            i_aop(1), i_ma(1), i_mm(1), i_bstar(1),
     &            r_j2k, v_j2k, ierr)
      call tle2sv(i_yr(2), i_eday(2),
     &            i_inc(2), i_raan(2), i_ecc(2),
     &            i_aop(2), i_ma(2), i_mm(2), i_bstar(2),
     &            r_out, v_out, ierr)
      dr = dsqrt((r_j2k(1)-r_out(1))**2 + (r_j2k(2)-r_out(2))**2
     &         + (r_j2k(3)-r_out(3))**2)
      if (dr .gt. 100.d0) then
         call assert_i('batch ISS != NOAA18 (dr>100)',
     &                 1, 1, npass, nfail)
      else
         call assert_i('batch ISS != NOAA18 (dr>100)',
     &                 0, 1, npass, nfail)
      end if

c     T72: Vanguard (high ecc) has valid energy ratio
      call tle2sv(i_yr(3), i_eday(3),
     &            i_inc(3), i_raan(3), i_ecc(3),
     &            i_aop(3), i_ma(3), i_mm(3), i_bstar(3),
     &            r_j2k, v_j2k, ierr)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      vmag = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
      energy = 0.5d0*vmag*vmag - 398600.8d0/rmag
      sma = -398600.8d0 / (2.d0*energy)
c     Vanguard SMA from n=10.86: a = (mu/(n*2pi/86400)^2)^(1/3)
c     ~ 8614 km
      if (sma .gt. 7000.d0 .and. sma .lt. 10000.d0) then
         call assert_i('Vanguard SMA in [7k,10k] km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('Vanguard SMA in [7k,10k] km',
     &                 0, 1, npass, nfail)
      end if

c     ============================================================
c     SDP4 deep-space tests
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- SDP4 deep-space ---'

c     DS1: SDP4 init GEO (typical geostationary)
      call tle_epoch2jd(2026, 170.5d0, djd)
      call tle_sdp4_init(1.00274d0, 0.0003d0, 0.5d0,
     &                   75.d0, 0.d0, 180.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call assert_i('DS1 GEO init ierr=0', ierr, 0, npass, nfail)
      call assert_i('DS1 GEO iresfl=1 (sync)',
     &              iresfl, 1, npass, nfail)
      if (a0dd .gt. 6.0d0 .and. a0dd .lt. 7.5d0) then
         call assert_i('DS1 GEO a0dd ~6.6 Re',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS1 GEO a0dd ~6.6 Re',
     &                 0, 1, npass, nfail)
      end if

c     DS2: SDP4 init Molniya
      call tle_sdp4_init(2.006d0, 0.74d0, 63.4d0,
     &                   0.d0, 270.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call assert_i('DS2 Molniya init ierr=0',
     &              ierr, 0, npass, nfail)
      call assert_i('DS2 Molniya iresfl=2 (half-day)',
     &              iresfl, 2, npass, nfail)

c     DS3: SDP4 init GPS
      call tle_sdp4_init(2.005d0, 0.005d0, 55.d0,
     &                   100.d0, 0.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call assert_i('DS3 GPS init ierr=0', ierr, 0, npass, nfail)
      call assert_i('DS3 GPS iresfl=2 (half-day)',
     &              iresfl, 2, npass, nfail)

c     DS4: SDP4 init non-resonant (n=1.5, period ~960 min)
      call tle_sdp4_init(1.5d0, 0.3d0, 28.d0,
     &                   0.d0, 0.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call assert_i('DS4 non-res init ierr=0',
     &              ierr, 0, npass, nfail)
      call assert_i('DS4 non-res iresfl=0',
     &              iresfl, 0, npass, nfail)

c     DS5: SDP4 prop GEO t=0
      call tle_sdp4_init(1.00274d0, 0.0003d0, 0.5d0,
     &                   75.d0, 0.d0, 180.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call tle_sdp4_prop(0.d0,
     &              0.0003d0, 0.5d0, 75.d0, 0.d0, 180.d0,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1_coeff, c4_coeff, perigee_km,
     &              omgdot, xnodot, xmdot, 0.d0,
     &              dsstate, iresfl,
     &              r_teme, v_teme, ierr)
      call assert_i('DS5 GEO prop t=0 ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_teme(1)**2+r_teme(2)**2+r_teme(3)**2)
      vmag = dsqrt(v_teme(1)**2+v_teme(2)**2+v_teme(3)**2)
      if (rmag .gt. 41000.d0 .and. rmag .lt. 43500.d0) then
         call assert_i('DS5 GEO rmag [41k,43.5k] km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS5 GEO rmag [41k,43.5k] km',
     &                 0, 1, npass, nfail)
      end if
      if (vmag .gt. 2.8d0 .and. vmag .lt. 3.3d0) then
         call assert_i('DS5 GEO vmag [2.8,3.3] km/s',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS5 GEO vmag [2.8,3.3] km/s',
     &                 0, 1, npass, nfail)
      end if

c     DS6: SDP4 prop Molniya t=0
      call tle_sdp4_init(2.006d0, 0.74d0, 63.4d0,
     &                   0.d0, 270.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call tle_sdp4_prop(0.d0,
     &              0.74d0, 63.4d0, 0.d0, 270.d0, 0.d0,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1_coeff, c4_coeff, perigee_km,
     &              omgdot, xnodot, xmdot, 0.d0,
     &              dsstate, iresfl,
     &              r_teme, v_teme, ierr)
      call assert_i('DS6 Molniya prop t=0 ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_teme(1)**2+r_teme(2)**2+r_teme(3)**2)
      vmag = dsqrt(v_teme(1)**2+v_teme(2)**2+v_teme(3)**2)
      if (rmag .gt. 6000.d0) then
         call assert_i('DS6 Molniya rmag > 6000 km',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS6 Molniya rmag > 6000 km',
     &                 0, 1, npass, nfail)
      end if
      energy = 0.5d0*vmag*vmag - 398600.8d0/rmag
      if (energy .lt. 0.d0) then
         call assert_i('DS6 Molniya energy < 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS6 Molniya energy < 0',
     &                 0, 1, npass, nfail)
      end if

c     DS7: SDP4 prop GEO t=1440 (1 day)
      call tle_sdp4_init(1.00274d0, 0.0003d0, 0.5d0,
     &                   75.d0, 0.d0, 180.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call tle_sdp4_prop(1440.d0,
     &              0.0003d0, 0.5d0, 75.d0, 0.d0, 180.d0,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1_coeff, c4_coeff, perigee_km,
     &              omgdot, xnodot, xmdot, 0.d0,
     &              dsstate, iresfl,
     &              r_teme, v_teme, ierr)
      call assert_i('DS7 GEO prop t=1440 ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_teme(1)**2+r_teme(2)**2+r_teme(3)**2)
      if (rmag .gt. 41000.d0 .and. rmag .lt. 43500.d0) then
         call assert_i('DS7 GEO 1day rmag [41k,43.5k]',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS7 GEO 1day rmag [41k,43.5k]',
     &                 0, 1, npass, nfail)
      end if

c     DS8: SDP4 prop Molniya t=1440
      call tle_sdp4_init(2.006d0, 0.74d0, 63.4d0,
     &                   0.d0, 270.d0, 0.d0,
     &                   0.d0, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0,
     &                   perigee_km,
     &                   omgdot, xnodot, xmdot,
     &                   dsstate, iresfl, ierr)
      call tle_sdp4_prop(1440.d0,
     &              0.74d0, 63.4d0, 0.d0, 270.d0, 0.d0,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1_coeff, c4_coeff, perigee_km,
     &              omgdot, xnodot, xmdot, 0.d0,
     &              dsstate, iresfl,
     &              r_teme, v_teme, ierr)
      call assert_i('DS8 Molniya prop t=1440 ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_teme(1)**2+r_teme(2)**2+r_teme(3)**2)
      if (rmag .gt. 6000.d0) then
         call assert_i('DS8 Molniya 1day rmag>6000',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS8 Molniya 1day rmag>6000',
     &                 0, 1, npass, nfail)
      end if
c     Check state is finite
      if (dabs(r_teme(1)).lt.1.d10 .and.
     &    dabs(r_teme(2)).lt.1.d10 .and.
     &    dabs(r_teme(3)).lt.1.d10) then
         call assert_i('DS8 Molniya state finite',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS8 Molniya state finite',
     &                 0, 1, npass, nfail)
      end if

c     DS9: tle2sv GEO full pipeline
      call tle2sv(2026, 170.5d0, 0.5d0, 75.d0, 0.0003d0,
     &            0.d0, 180.d0, 1.00274d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('DS9 tle2sv GEO ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      if (rmag .gt. 41000.d0 .and. rmag .lt. 43500.d0) then
         call assert_i('DS9 GEO rmag [41k,43.5k]',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS9 GEO rmag [41k,43.5k]',
     &                 0, 1, npass, nfail)
      end if
      energy = 0.5d0*(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
     &       - 398600.8d0/rmag
      if (energy .lt. 0.d0) then
         call assert_i('DS9 GEO energy < 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS9 GEO energy < 0',
     &                 0, 1, npass, nfail)
      end if

c     DS10: tle2sv Molniya full pipeline
      call tle2sv(2026, 170.5d0, 63.4d0, 0.d0, 0.74d0,
     &            270.d0, 0.d0, 2.006d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('DS10 tle2sv Molniya ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      vmag = dsqrt(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
      energy = 0.5d0*vmag*vmag - 398600.8d0/rmag
      if (energy .lt. 0.d0) then
         call assert_i('DS10 Molniya energy < 0',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS10 Molniya energy < 0',
     &                 0, 1, npass, nfail)
      end if
      if (dabs(r_j2k(1)).lt.1.d10 .and.
     &    dabs(r_j2k(2)).lt.1.d10 .and.
     &    dabs(r_j2k(3)).lt.1.d10) then
         call assert_i('DS10 Molniya state finite',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS10 Molniya state finite',
     &                 0, 1, npass, nfail)
      end if

c     DS11: tle2sv GPS full pipeline
      call tle2sv(2026, 170.5d0, 55.d0, 100.d0, 0.005d0,
     &            0.d0, 0.d0, 2.005d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('DS11 tle2sv GPS ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      if (rmag .gt. 26000.d0 .and. rmag .lt. 27500.d0) then
         call assert_i('DS11 GPS rmag [26k,27.5k]',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS11 GPS rmag [26k,27.5k]',
     &                 0, 1, npass, nfail)
      end if

c     DS12: Near-Earth (ISS) unchanged, still uses SGP4
      call tle2sv(2026, 170.5d0, 51.64d0, 123.4567d0, 0.00067d0,
     &            45.1234d0, 315.d0, 15.54d0, 1.027d-4,
     &            r_j2k, v_j2k, ierr)
      call assert_i('DS12 ISS still SGP4 ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      if (rmag .gt. 6400.d0 .and. rmag .lt. 7200.d0) then
         call assert_i('DS12 ISS rmag [6400,7200]',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS12 ISS rmag [6400,7200]',
     &                 0, 1, npass, nfail)
      end if

c     DS13: Period boundary (n=6.4, period ~ 225 min exactly)
      call tle2sv(2026, 170.5d0, 28.5d0, 0.d0, 0.01d0,
     &            0.d0, 0.d0, 6.4d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call assert_i('DS13 boundary n=6.4 ierr=0',
     &              ierr, 0, npass, nfail)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      if (rmag .gt. 5000.d0) then
         call assert_i('DS13 boundary rmag > 5000',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS13 boundary rmag > 5000',
     &                 0, 1, npass, nfail)
      end if

c     DS14: Deep-space z bounded by inclination
      call tle2sv(2026, 170.5d0, 0.5d0, 75.d0, 0.0003d0,
     &            0.d0, 180.d0, 1.00274d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      rmag = dsqrt(r_j2k(1)**2+r_j2k(2)**2+r_j2k(3)**2)
      zmax = rmag * dsin(0.5d0*d2r + 2.d0*d2r)
      if (dabs(r_j2k(3)) .le. zmax) then
         call assert_i('DS14 GEO z bounded by inc',
     &                 1, 1, npass, nfail)
      else
         call assert_i('DS14 GEO z bounded by inc',
     &                 0, 1, npass, nfail)
      end if

c     DS15: Repeatability - two SDP4 calls identical
      call tle2sv(2026, 170.5d0, 0.5d0, 75.d0, 0.0003d0,
     &            0.d0, 180.d0, 1.00274d0, 0.d0,
     &            r_j2k, v_j2k, ierr)
      call tle2sv(2026, 170.5d0, 0.5d0, 75.d0, 0.0003d0,
     &            0.d0, 180.d0, 1.00274d0, 0.d0,
     &            r_out, v_out, ierr)
      dr = dsqrt((r_j2k(1)-r_out(1))**2 + (r_j2k(2)-r_out(2))**2
     &         + (r_j2k(3)-r_out(3))**2)
      call assert_d('DS15 SDP4 repeatability dr=0',
     &              dr, 0.d0, 1.d-10, npass, nfail)

c     ============================================================
c     Summary
c     ============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,100) npass+nfail, npass, nfail
 100  format(' Total: ',I4,'   Passed: ',I4,'   Failed: ',I4)
      write(*,*) '=============================================='

      if (nfail .gt. 0) then
         write(*,*) 'SOME TESTS FAILED'
      else
         write(*,*) 'ALL TESTS PASSED'
      end if

      stop
      end

c     ============================================================
c     Assertion helpers (same as test_tle.F)
c     ============================================================
      subroutine assert_d(label, got, expect, tol, npass, nfail)
      implicit double precision (a-h, o-z)
      character*(*) label
      double precision got, expect, tol, aerr
      integer npass, nfail

      aerr = dabs(got - expect)
      if (aerr .le. tol) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,E20.12,A,E20.12)') '           got=',
     &        got, '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

      subroutine assert_i(label, got, expect, npass, nfail)
      character*(*) label
      integer got, expect, npass, nfail
      if (got .eq. expect) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,I10,A,I10)') '           got=',
     &        got, '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

      subroutine assert_s(label, got, expect, npass, nfail)
      character*(*) label, got, expect
      integer npass, nfail
      if (got .eq. expect) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,A,A,A,A)') '           got=[',
     &        got, ']  expect=[', expect, ']'
         nfail = nfail + 1
      end if
      return
      end
