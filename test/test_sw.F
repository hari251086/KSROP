c     test_sw.F — Space weather + 2-D atmosphere tests
c     -------------------------------------------------------
c     Ported from OREM's test/test_sw.F (issue #26 there) to cover
c     the loader/interpolation logic ported into src/swx.F +
c     app/driver_KS.F's sw_tinf/atm2d_interp. Needs input/SW-All.csv
c     and input/ATM2D.DAT (both committed).
c     -------------------------------------------------------
      program test_sw

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      integer npass, nfail, ierr
      double precision cal(6), jd, tinf
      double precision rho, hs, rho72, hs72, rho_hi
      integer ist
c     Issue #32 Phase 5 follow-up: diurnal density bulge
      double precision tbmax, tbmin, ts_bulge(3), xb(3), cphi

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' Space weather / ATM2D test suite'
      write(*,*) '=============================================='

c     ---------------------------------------------------------
c     W1-W2: defaults — nothing loaded, legacy contract holds
c     ---------------------------------------------------------
      call sw_tinf(2458000.d0, tinf)
      call assert_d('W1 no data: sw_tinf=0', tinf, 0.d0, 1.d-12,
     &              npass, nfail)
      call atm2d_interp(170.d3, 700.d0, rho, hs, ist)
      call assert_i('W2 no table: istat=2', ist, 2, npass, nfail)

c     ---------------------------------------------------------
c     W3-W6: sw_load + known-row lookups
c     ---------------------------------------------------------
      call sw_load('input/SW-All.csv', ierr)
      call assert_i('W3 sw_load ierr=0', ierr, 0, npass, nfail)

c     2024-05-11 (G5 storm): F10.7_ADJ=218.0, ADJ_C81=180.5,
c     KP_SUM=670 -> mean Kp=8.375
c     Tc   = 379 + 3.24*180.5 + 1.3*(218.0-180.5) = 1012.57
c     Tinf = Tc + 14*8.375 + 0.02*exp(8.375)
c          = 1012.57 + 117.25 + 86.74 = 1216.56
      cal(1) = 2024.d0
      cal(2) = 5.d0
      cal(3) = 11.d0
      cal(4) = 12.d0
      cal(5) = 0.d0
      cal(6) = 0.d0
      call cal2jd(cal, jd)
      call sw_tinf(jd, tinf)
      call assert_d('W4 2024-05-11 storm Tinf', tinf, 1216.56d0,
     &              2.d0, npass, nfail)

c     quiet solar minimum 2019-12-15: Tinf must land in [590,650] K
      cal(1) = 2019.d0
      cal(2) = 12.d0
      cal(3) = 15.d0
      call cal2jd(cal, jd)
      call sw_tinf(jd, tinf)
      if (tinf .gt. 590.d0 .and. tinf .lt. 650.d0) then
         call assert_i('W5 2019 solar-min Tinf in [590,650]',
     &                 1, 1, npass, nfail)
      else
         write(*,'(A,F8.1)') '   got Tinf=', tinf
         call assert_i('W5 2019 solar-min Tinf in [590,650]',
     &                 0, 1, npass, nfail)
      end if

c     predicted era (monthly rows): 2030-06-15 must be physical
      cal(1) = 2030.d0
      cal(2) = 6.d0
      cal(3) = 15.d0
      call cal2jd(cal, jd)
      call sw_tinf(jd, tinf)
      if (tinf .gt. 500.d0 .and. tinf .lt. 1500.d0) then
         call assert_i('W6 predicted era Tinf physical',
     &                 1, 1, npass, nfail)
      else
         write(*,'(A,F8.1)') '   got Tinf=', tinf
         call assert_i('W6 predicted era Tinf physical',
     &                 0, 1, npass, nfail)
      end if

c     ---------------------------------------------------------
c     W7-W11: atm2d_load + interpolation consistency
c     ---------------------------------------------------------
      call atm2d_load('input/ATM2D.DAT', ierr)
      call assert_i('W7 atm2d_load ierr=0', ierr, 0, npass, nfail)

c     at Tinf=626.33 (the static table's generation temperature),
c     rho(170 km) must match ATM.DAT within ~2% interpolation error
      call atm2d_interp(170.d3, 626.33d0, rho, hs, ist)
      call assert_i('W8 interp istat=0', ist, 0, npass, nfail)
      rho72 = 4.866d-10 * 1.d10
      if (dabs(rho - rho72)/rho72 .lt. 0.02d0) then
         call assert_i('W9 rho(170km,626K) ~ static table',
     &                 1, 1, npass, nfail)
      else
         write(*,'(A,E12.4,A,E12.4)') '   got ', rho,
     &      ' want ~', rho72
         call assert_i('W9 rho(170km,626K) ~ static table',
     &                 0, 1, npass, nfail)
      end if

c     monotonicity in Tinf at 300 km (scaled units)
      call atm2d_interp(300.d3, 600.d0, rho, hs, ist)
      call atm2d_interp(300.d3, 1200.d0, rho_hi, hs72, ist)
      if (rho_hi .gt. 3.d0*rho) then
         call assert_i('W10 rho(300km) grows strongly with Tinf',
     &                 1, 1, npass, nfail)
      else
         call assert_i('W10 rho(300km) grows strongly with Tinf',
     &                 0, 1, npass, nfail)
      end if

c     altitude out of range -> istat=1 (drag-inactive contract)
      call atm2d_interp(550.d3, 700.d0, rho, hs, ist)
      call assert_i('W11 h>500km istat=1', ist, 1, npass, nfail)

c     ---------------------------------------------------------
c     W12-W17: diurnal density bulge (issue #32 Phase 5 follow-up,
c     ported from OREM's test/test_sw.F). At the equinox
c     (decl_sun=0), phi_max=phi_min=0 exactly (atan(0)=0), so the
c     extrema formula collapses to an exact closed form with no
c     transcendental-function rounding risk:
c       Tmax = Tmean*(1 + 0.24*Kmax) = Tmean*1.130992
c       Tmin = Tmean*(1 + 0.24*Kmin) = Tmean*0.862672
c     ---------------------------------------------------------
      call jacchia_bulge_extrema(1000.d0, 0.d0, tbmax, tbmin)
      call assert_d('W12 bulge extrema: equinox Tmax',
     &              tbmax, 1130.992d0, 1.d-3, npass, nfail)
      call assert_d('W13 bulge extrema: equinox Tmin',
     &              tbmin, 862.672d0, 1.d-3, npass, nfail)

c     At any nonzero declination the extrema must still bracket the
c     mean (qualitative sanity, not a specific transcendental value)
      call jacchia_bulge_extrema(1000.d0, 23.439d0, tbmax, tbmin)
      if (tbmax .gt. 1000.d0 .and. tbmin .lt. 1000.d0) then
         call assert_i('W14 bulge extrema: solstice brackets mean',
     &                 1, 1, npass, nfail)
      else
         call assert_i('W14 bulge extrema: solstice brackets mean',
     &                 0, 1, npass, nfail)
      end if

c     Sun fixed on the equatorial x-axis (dec_sun=0, ra_sun=0) ->
c     bulge direction is dec=0, ra=37deg exactly. Three exact
c     geometric cases: satellite AT the bulge direction, at its
c     antipode, and 90deg away.
      ts_bulge(1) = 1.d0
      ts_bulge(2) = 0.d0
      ts_bulge(3) = 0.d0

      xb(1) = dcos(37.d0*d2r)
      xb(2) = dsin(37.d0*d2r)
      xb(3) = 0.d0
      call bulge_cosphi(xb, ts_bulge, 37.d0, cphi)
      call assert_d('W15 bulge_cosphi: at bulge direction',
     &              cphi, 1.d0, 1.d-10, npass, nfail)

      xb(1) = -dcos(37.d0*d2r)
      xb(2) = -dsin(37.d0*d2r)
      xb(3) = 0.d0
      call bulge_cosphi(xb, ts_bulge, 37.d0, cphi)
      call assert_d('W16 bulge_cosphi: antipodal',
     &              cphi, -1.d0, 1.d-10, npass, nfail)

      xb(1) = dcos(127.d0*d2r)
      xb(2) = dsin(127.d0*d2r)
      xb(3) = 0.d0
      call bulge_cosphi(xb, ts_bulge, 37.d0, cphi)
      call assert_d('W17 bulge_cosphi: 90deg away',
     &              cphi, 0.d0, 1.d-10, npass, nfail)

      write(*,*) ''
      write(*,*) '=============================================='
      write(*,100) npass + nfail, npass, nfail
 100  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='
      if (nfail .gt. 0) then
         write(*,*) 'SOME TESTS FAILED'
         stop 1
      else
         write(*,*) 'ALL TESTS PASSED'
      end if
      stop
      end

c     ============================================================
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

c     ============================================================
      subroutine assert_d(label, got, expect, tol, npass, nfail)
      implicit double precision (a-h, o-z)
      character*(*) label
      double precision got, expect, tol
      integer npass, nfail
      if (dabs(got - expect) .le. tol) then
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
