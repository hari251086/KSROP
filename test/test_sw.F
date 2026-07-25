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
