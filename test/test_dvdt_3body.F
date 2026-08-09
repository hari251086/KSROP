c     test_dvdt_3body.F -- time-dependent third-body geopotential in w*
c     -------------------------------------------------------
c     KSROP branch feature/tesseral-energy-time-dependence. Verifies
c     third_body_dvdt (src/Subrouts.F) against a genuine finite
c     difference of the REAL third-body potential, computed via the
c     actual production solarnpv/lunarpv ephemeris at t0+/-H (a
c     LARGER, independent step than the internal h_fd=60s velocity
c     estimate third_body_dvdt itself uses) -- the satellite position
c     x is held exactly fixed throughout, matching the derivation's
c     own partial-derivative convention.
c
c     DS1-DS4: Sun, four (x, epoch) cases.
c     DM1-DM4: Moon, same four cases.
c     ZD1: degree<2 guard -- dvdt_out must be exactly zero (no
c     n=2..1 terms exist), matching how the caller in driver_KS.F
c     only calls this when nsun_deg/nmoon_deg >= 2.
c     -------------------------------------------------------
      program test_dvdt_3body

      implicit double precision (a-h, o-z)
      integer npass, nfail
      integer maxdeg, deg
      double precision leg
      parameter (maxdeg = 8)
      dimension x(3), tb0(3), tb_dot(3), tb_p(3), tb_m(3)
      dimension RbyRtb(maxdeg+2), leg(maxdeg+2)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth

      npass = 0
      nfail = 0
      call init_constants()

      write(*,*) '=============================================='
      write(*,*) ' Time-dependent third-body geopotential in w*'
      write(*,*) '=============================================='

      h_fd = 60.d0/86400.d0
      H_test = 300.d0/86400.d0

      do it = 1, 4
         if (it .eq. 1) then
            x(1) = 26000.d0
            x(2) = 8000.d0
            x(3) = 3000.d0
            dj0 = 2461000.5d0
         else if (it .eq. 2) then
            x(1) = -15000.d0
            x(2) = 22000.d0
            x(3) = -5000.d0
            dj0 = 2461050.25d0
         else if (it .eq. 3) then
            x(1) = 10000.d0
            x(2) = -30000.d0
            x(3) = 6000.d0
            dj0 = 2461120.75d0
         else
            x(1) = -20000.d0
            x(2) = -12000.d0
            x(3) = 15000.d0
            dj0 = 2461200.0d0
         end if
         R1 = dsqrt(x(1)**2+x(2)**2+x(3)**2)

c        --- Sun ---
         deg = 6
         call solarnpv(dj0, tb0)
         call solarnpv(dj0+h_fd, tb_p)
         call solarnpv(dj0-h_fd, tb_m)
         do i = 1,3
            tb_dot(i) = (tb_p(i)-tb_m(i))/(2.d0*h_fd*86400.d0)
         end do
         call third_body_aux(deg, x, tb0, R1, amuS,
     &                        Rtb1, RbyRtb, first_tb, cphi, cphi2)
         call aLegP(deg+1, cphi, leg)
         call third_body_dvdt(deg, x, tb0, tb_dot, R1, amuS,
     &      Rtb1, RbyRtb, cphi, cphi2, leg, dvdt_closed)

         call third_body_U_calc(deg, x, dj0+H_test, R1, amuS,
     &                           .true., u_p)
         call third_body_U_calc(deg, x, dj0-H_test, R1, amuS,
     &                           .true., u_m)
         dvdt_fd = (u_p-u_m)/(2.d0*H_test*86400.d0)

         tol = 1.d-4*dabs(dvdt_closed) + 1.d-14
         write(*,'(A,I2)') '  --- case ', it
         call assert_d('  DS: Sun dvdt vs finite difference',
     &                 dvdt_closed, dvdt_fd, tol, npass, nfail)

c        --- Moon ---
         call lunarpv(dj0, tb0)
         call lunarpv(dj0+h_fd, tb_p)
         call lunarpv(dj0-h_fd, tb_m)
         do i = 1,3
            tb_dot(i) = (tb_p(i)-tb_m(i))/(2.d0*h_fd*86400.d0)
         end do
         call third_body_aux(deg, x, tb0, R1, amuM,
     &                        Rtb1, RbyRtb, first_tb, cphi, cphi2)
         call aLegP(deg+1, cphi, leg)
         call third_body_dvdt(deg, x, tb0, tb_dot, R1, amuM,
     &      Rtb1, RbyRtb, cphi, cphi2, leg, dvdt_closed)

         call third_body_U_calc(deg, x, dj0+H_test, R1, amuM,
     &                           .false., u_p)
         call third_body_U_calc(deg, x, dj0-H_test, R1, amuM,
     &                           .false., u_m)
         dvdt_fd = (u_p-u_m)/(2.d0*H_test*86400.d0)

         tol = 1.d-4*dabs(dvdt_closed) + 1.d-14
         call assert_d('  DM: Moon dvdt vs finite difference',
     &                 dvdt_closed, dvdt_fd, tol, npass, nfail)
      end do

c     ---------------------------------------------------------
c     ZD1: degree<2 guard
c     ---------------------------------------------------------
      x(1) = 26000.d0
      x(2) = 8000.d0
      x(3) = 3000.d0
      R1 = dsqrt(x(1)**2+x(2)**2+x(3)**2)
      call solarnpv(2461000.5d0, tb0)
      tb_dot(1) = 0.d0
      tb_dot(2) = 0.d0
      tb_dot(3) = 0.d0
      call third_body_aux(1, x, tb0, R1, amuS,
     &                     Rtb1, RbyRtb, first_tb, cphi, cphi2)
      call aLegP(2, cphi, leg)
      call third_body_dvdt(1, x, tb0, tb_dot, R1, amuS,
     &   Rtb1, RbyRtb, cphi, cphi2, leg, dvdt_closed)
      call assert_d('ZD1 deg=1 (no n=2..1 terms): dvdt exactly zero',
     &              dvdt_closed, 0.d0, 1.d-30, npass, nfail)

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

      end

c     ============================================================
      subroutine third_body_U_calc(deg, x, dj, R1, amu_tb, is_sun,
     &                              u_out)
c     Real third-body potential U(x,t) = -amu_tb/Rtb * Sum (R1/Rtb)^n
c     P_n(cphi), evaluated at epoch dj via the actual production
c     ephemeris (solarnpv/lunarpv), x held fixed -- independent
c     ground truth for the finite-difference check above (deliberately
c     does NOT reuse third_body_dvdt or its own internal machinery).
c     A subroutine (output argument), not a function, matching this
c     repo's own established style throughout -- avoids relying on an
c     external function call with no explicit interface.
c     ============================================================
      implicit double precision (a-h, o-z)
      integer deg
      logical is_sun
      double precision leg
      dimension x(3), tb(3), RbyRtb(deg+2), leg(deg+2)

      if (is_sun) then
         call solarnpv(dj, tb)
      else
         call lunarpv(dj, tb)
      end if

      call third_body_aux(deg, x, tb, R1, amu_tb,
     &                     Rtb1, RbyRtb, first_tb, cphi, cphi2)
      call aLegP(deg+1, cphi, leg)

      usum = 0.d0
      do i = 2, deg
         usum = usum + first_tb*RbyRtb(i)*leg(i+1)
      end do
      u_out = usum

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
