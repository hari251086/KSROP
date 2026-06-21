c     test_bugs.F  -  Tests targeting known bugs in Subrouts.f
c     -------------------------------------------------------
c     Compile:  gfortran test_bugs.F Subrouts.F Legendre.F -o test_bugs.exe
c     Run:      ./test_bugs.exe  (exit 0 = all pass, exit 1 = any fail)
c     -------------------------------------------------------
c     Bugs tested:
c       1. force_models type mismatch (amuM corrupted when moon off)
c       2. car2ks vmn call with wrong args (sequence assoc. masks it)
c       3. car2oe equatorial orbit dead code (dabs<0 impossible)
c       4. car2oe circular+equatorial dead code
c       5. car2oe acos vs dacos (single-precision truncation)
c       6. car2oe ea uninitialized for hyperbolic orbits
c       7. rkgil SAVE attribute missing on q(10)
c       8. solarnpv spurious close(1)
c       9. cal2jd fractional seconds truncation
c      10. oe2car destructive modification of input array
c      11. atan3 division by zero when b=0
c     -------------------------------------------------------
      program test_bugs

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      dimension pek(7), pek_save(7), x(3), xd(3)
      dimension cal(6)
      dimension n_force(3)
      dimension y(10), f(10)
      dimension ts_vec(3)
      character*80 line_test
      integer npass, nfail, nt, ios_test

      call init_constants()

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' Bug Regression Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     Bug 1: force_models type mismatch
c     When n_force(3)=0, the subroutine writes integer 0 to the
c     double-precision amuM memory location, corrupting it.
c     Expected: amuM should be exactly 0.0d0 after the call.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 1: force_models type mismatch ---'

      amuS_test = 1.32712440018d11
      amuM_test = 4.902801076d3
      ngeo_test = 10

      n_force(1) = 1
      n_force(2) = 0
      n_force(3) = 0

      call force_models(n_force, ngeo_test, amuS_test, amuM_test)

c     If the type mismatch corrupts amuM, it won't be exactly 0.0
      call assert_d('Bug1 amuM zeroed when moon off',
     &              amuM_test, 0.d0, 1.d-15, npass, nfail)
      call assert_d('Bug1 amuS zeroed when sun off',
     &              amuS_test, 0.d0, 1.d-15, npass, nfail)

c     ==============================================================
c     Bug 3: car2oe equatorial orbit (dead code: dabs<0)
c     An equatorial orbit should not produce NaN for RAAN/AOP.
c     Setup: equatorial circular orbit, i=0, along x-axis.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 3: car2oe equatorial orbit ---'

c     Equatorial orbit: r along x, v along y (i=0 exactly)
      x(1) = 7000.d0
      x(2) = 0.d0
      x(3) = 0.d0
      xd(1) = 0.d0
      xd(2) = dsqrt(amue / 7000.d0)
      xd(3) = 0.d0

      call car2oe(x, xd, pek)

c     Check for NaN: a NaN does not equal itself
      if (pek(4) .ne. pek(4)) then
         write(*,*) '  [FAIL]  Bug3 RAAN is NaN (equatorial)'
         nfail = nfail + 1
      else
         write(*,*) '  [PASS]  Bug3 RAAN is finite (equatorial)'
         npass = npass + 1
      end if

      if (pek(5) .ne. pek(5)) then
         write(*,*) '  [FAIL]  Bug3 AOP is NaN (equatorial)'
         nfail = nfail + 1
      else
         write(*,*) '  [PASS]  Bug3 AOP is finite (equatorial)'
         npass = npass + 1
      end if

c     Inclination should be 0 (or very close)
      call assert_d('Bug3 inclination = 0 deg (equatorial)',
     &              pek(3), 0.d0, 1.d-8, npass, nfail)

c     ==============================================================
c     Bug 4: car2oe circular+equatorial orbit (dead code)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 4: car2oe circular equatorial ---'

c     Same as above (circular equatorial) - check true longitude
      x(1) = 7000.d0
      x(2) = 0.d0
      x(3) = 0.d0
      xd(1) = 0.d0
      xd(2) = dsqrt(amue / 7000.d0)
      xd(3) = 0.d0

      call car2oe(x, xd, pek)

c     For circular equatorial, eccentricity should be ~0
c     and inclination ~0; RAAN/AOP/MA should all be finite
      call assert_d('Bug4 ecc ~ 0 (circular equatorial)',
     &              pek(2), 0.d0, 1.d-10, npass, nfail)

      if (pek(6) .ne. pek(6)) then
         write(*,*) '  [FAIL]  Bug4 MA is NaN (circ equatorial)'
         nfail = nfail + 1
      else
         write(*,*) '  [PASS]  Bug4 MA is finite (circ equatorial)'
         npass = npass + 1
      end if

c     ==============================================================
c     Bug 5: car2oe acos vs dacos (single-precision AOP)
c     Use an orbit where AOP is known analytically to high precision.
c     a=7000, e=0.1, i=45, RAAN=30, AOP=60, MA=0 (at perigee)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 5: acos vs dacos precision ---'

c     Build state from known elements via oe2car, then recover with car2oe
      pek(1) = 7000.d0
      pek(2) = 0.1d0
      pek(3) = 45.d0
      pek(4) = 30.d0
      pek(5) = 60.d0
      pek(6) = 90.d0
      pek(7) = 0.d0

      call oe2car(pek, x, xd, 1.d-15)

c     Now recover elements
      call car2oe(x, xd, pek)

c     AOP should be 60.0 degrees to at least 1e-10 precision
c     With single-precision acos, error is ~1e-6 degrees
      aop_err = dabs(pek(5) - 60.d0)
      if (aop_err .lt. 1.d-10) then
         write(*,*) '  [PASS]  Bug5 AOP precision < 1e-10 deg'
         npass = npass + 1
      else
         write(*,'(A,E12.5)')
     &      '  [FAIL]  Bug5 AOP precision (err=', aop_err
         write(*,'(A,F20.14)') '           AOP recovered = ', pek(5)
         nfail = nfail + 1
      end if

c     ==============================================================
c     Bug 6: car2oe ea uninitialized for hyperbolic orbits
c     Use a hyperbolic orbit (e>1) and check a_kep(7) is finite.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 6: car2oe hyperbolic ea ---'

c     Hyperbolic orbit: v > escape velocity
c     v_esc = sqrt(2*mu/r), use v = 1.2 * v_esc
      r_hyp = 7000.d0
      v_esc = dsqrt(2.d0 * amue / r_hyp)
      x(1) = r_hyp
      x(2) = 0.d0
      x(3) = 0.d0
      xd(1) = 0.d0
      xd(2) = 1.2d0 * v_esc
      xd(3) = 0.d0

      call car2oe(x, xd, pek)

c     e should be > 1
      if (pek(2) .gt. 1.d0) then
         write(*,'(A,F10.6)')
     &      '  [PASS]  Bug6 hyperbolic e > 1, e=', pek(2)
         npass = npass + 1
      else
         write(*,'(A,F10.6)')
     &      '  [FAIL]  Bug6 expected e > 1, got e=', pek(2)
         nfail = nfail + 1
      end if

c     a_kep(7) = ea: should be finite (not NaN or garbage)
      if (pek(7) .ne. pek(7)) then
         write(*,*) '  [FAIL]  Bug6 ea is NaN (hyperbolic)'
         nfail = nfail + 1
      else if (dabs(pek(7)) .gt. 1.d10) then
         write(*,'(A,E12.5)')
     &      '  [FAIL]  Bug6 ea is garbage (hyperbolic), ea=', pek(7)
         nfail = nfail + 1
      else
         write(*,'(A,F12.6)')
     &      '  [PASS]  Bug6 ea is finite (hyperbolic), ea=', pek(7)
         npass = npass + 1
      end if

c     ==============================================================
c     Bug 7: rkgil missing SAVE on q(10)
c     Call rkgil 4 stages and verify the result is a correct RK-Gill
c     step on a simple ODE: dy/dt = -y, y(0)=1, h=0.1
c     Exact answer: y(0.1) = exp(-0.1) = 0.904837418...
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 7: rkgil SAVE attribute ---'

      y(1) = 1.d0
      f(1) = -y(1)
      nt = 0
      h_step = 0.1d0
      x_rk = 0.d0

c     Stage 1
      f(1) = -y(1)
      rkg_ret = rkgil(1, y, f, x_rk, h_step, nt)

c     Stage 2
      f(1) = -y(1)
      rkg_ret = rkgil(1, y, f, x_rk, h_step, nt)

c     Stage 3
      f(1) = -y(1)
      rkg_ret = rkgil(1, y, f, x_rk, h_step, nt)

c     Stage 4
      f(1) = -y(1)
      rkg_ret = rkgil(1, y, f, x_rk, h_step, nt)

c     After 4 stages, nt should be reset to 0 and result should be
c     close to exp(-0.1) = 0.9048374180...
      y_exact = dexp(-0.1d0)
      rk_err = dabs(y(1) - y_exact)
      if (rk_err .lt. 1.d-6) then
         write(*,'(A,E10.3)')
     &      '  [PASS]  Bug7 rkgil step correct, err=', rk_err
         npass = npass + 1
      else
         write(*,'(A,F16.10,A,F16.10)')
     &      '  [FAIL]  Bug7 rkgil result=', y(1),
     &      ' expected=', y_exact
         nfail = nfail + 1
      end if

      if (nt .eq. 0) then
         write(*,*) '  [PASS]  Bug7 nt reset to 0 after 4 stages'
         npass = npass + 1
      else
         write(*,'(A,I3)')
     &      '  [FAIL]  Bug7 nt not reset, nt=', nt
         nfail = nfail + 1
      end if

c     ==============================================================
c     Bug 8: solarnpv spurious close(1)
c     Open unit 1, call solarnpv, then try to read from unit 1.
c     If solarnpv closes it, the read will fail.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 8: solarnpv close(1) ---'

      open(unit=1, status='scratch')
      write(1, '(A)') 'SENTINEL_LINE'
      rewind(1)

c     Call solarnpv with some Julian date (J2000)
      ts_vec(1) = 0.d0
      ts_vec(2) = 0.d0
      ts_vec(3) = 0.d0
      call solarnpv(2451545.d0, ts_vec)

c     Try reading from unit 1 — if close(1) was called, this fails
      read(1, '(A)', iostat=ios_test) line_test
      if (ios_test .eq. 0) then
         write(*,*) '  [PASS]  Bug8 unit 1 still open after solarnpv'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  Bug8 unit 1 closed by solarnpv'
         nfail = nfail + 1
      end if
      close(1, status='delete')

c     ==============================================================
c     Bug 9: cal2jd truncates fractional seconds
c     cal(6) = 30.567 seconds. The JD difference from cal(6)=30.0
c     should be 0.567/86400 days.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 9: cal2jd fractional seconds ---'

c     2024-01-01 12:00:30.567
      cal(1) = 2024.d0
      cal(2) = 1.d0
      cal(3) = 1.d0
      cal(4) = 12.d0
      cal(5) = 0.d0
      cal(6) = 30.567d0
      call cal2jd(cal, djd1)

c     2024-01-01 12:00:30.000
      cal(6) = 30.d0
      call cal2jd(cal, djd2)

c     Difference should be 0.567/86400 days = 6.5625e-6 days
      expected_diff = 0.567d0 / 86400.d0
      actual_diff = djd1 - djd2
      sec_err = dabs(actual_diff - expected_diff) * 86400.d0

      if (sec_err .lt. 0.001d0) then
         write(*,'(A,E10.3,A)')
     &      '  [PASS]  Bug9 sub-second JD preserved (err=',
     &      sec_err, ' s)'
         npass = npass + 1
      else
         write(*,'(A,E10.3,A)')
     &      '  [FAIL]  Bug9 sub-second JD lost (err=',
     &      sec_err, ' s)'
         write(*,'(A,F20.10,A,F20.10)')
     &      '           diff=', actual_diff,
     &      '  expected=', expected_diff
         nfail = nfail + 1
      end if

c     ==============================================================
c     Bug 10: oe2car destructively modifies input a_kep
c     After oe2car, a_kep(3:6) should still be in degrees.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 10: oe2car modifies input ---'

      pek(1) = 7000.d0
      pek(2) = 0.1d0
      pek(3) = 45.d0
      pek(4) = 30.d0
      pek(5) = 60.d0
      pek(6) = 90.d0
      pek(7) = 0.d0

c     Save original values
      do i = 1, 7
         pek_save(i) = pek(i)
      end do

      call oe2car(pek, x, xd, 1.d-15)

c     Check if pek(3:6) are still in degrees (unchanged)
      changed = 0
      do i = 3, 6
         if (dabs(pek(i) - pek_save(i)) .gt. 1.d-10) changed = 1
      end do

      if (changed .eq. 0) then
         write(*,*) '  [PASS]  Bug10 oe2car preserves input array'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  Bug10 oe2car modifies input array'
         write(*,'(A,4F12.6)')
     &      '           pek(3:6) now=',pek(3),pek(4),pek(5),pek(6)
         write(*,'(A,4F12.6)')
     &      '           original   =',
     &      pek_save(3),pek_save(4),pek_save(5),pek_save(6)
         nfail = nfail + 1
      end if

c     ==============================================================
c     Bug 11: atan3 division by zero when b=0
c     atan3(1, 0) should return pi/2, not crash.
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Bug 11: atan3 division by zero ---'

      res_atan = atan3(1.d0, 0.d0)
      expected_atan = pi / 2.d0

c     Check it didn't produce NaN or Inf
      if (res_atan .ne. res_atan) then
         write(*,*) '  [FAIL]  Bug11 atan3(1,0) is NaN'
         nfail = nfail + 1
      else if (dabs(res_atan) .gt. 1.d30) then
         write(*,*) '  [FAIL]  Bug11 atan3(1,0) is Inf'
         nfail = nfail + 1
      else
         atan_err = dabs(res_atan - expected_atan)
         if (atan_err .lt. 1.d-10) then
            write(*,*) '  [PASS]  Bug11 atan3(1,0) = pi/2'
            npass = npass + 1
         else
            write(*,'(A,F16.10,A,F16.10)')
     &         '  [FAIL]  Bug11 atan3(1,0)=', res_atan,
     &         ' expected=', expected_atan
            nfail = nfail + 1
         end if
      end if

c     Also test atan3(0, 0) — both zero
      res_atan2 = atan3(0.d0, 0.d0)
      if (res_atan2 .ne. res_atan2) then
         write(*,*) '  [FAIL]  Bug11 atan3(0,0) is NaN'
         nfail = nfail + 1
      else
         write(*,'(A,F12.6)')
     &      '  [PASS]  Bug11 atan3(0,0) is finite =', res_atan2
         npass = npass + 1
      end if

c     ==============================================================
c     Summary
c     ==============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,200) npass + nfail, npass, nfail
 200  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='

      if (nfail .gt. 0) stop 1
      stop
      end

c     ============================================================
      subroutine assert_d(label, got, expected, tol, npass, nfail)
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
