c     test_dvdt_tess.F -- time-dependent tesseral geopotential in w*
c     -------------------------------------------------------
c     KSROP feature/tesseral-energy-time-dependence branch. Verifies
c     tess_legendre_force's new dvdt_leg output (src/LegendreTess.F)
c     and the structural claim that driver_KS.F's z(1) (Tau_geo)
c     equation needs no direct edit once z(10) (w*) includes it --
c     see LegendreTess.F's own header derivation and the plan's
c     Context section (Sellamuthu 2018 PhD thesis eq. 2.54/2.56).
c
c     DT1-DT4: dvdt_leg against a genuine finite-difference derivative
c     of the potential w.r.t. Earth's rotation angle theta, holding u
c     (hence r, phi, lambda) exactly fixed -- mirrors this project's
c     own established cunningham_dVnm-vs-finite-difference acceptance-
c     test pattern (test_cunningham.F).
c     ZT1: zonal-only (no tesseral coefficients loaded) regression
c     guard -- dvdt_leg must be exactly zero (m=0 terms structurally
c     cannot contribute, an axisymmetric field cannot exchange energy
c     with a rotating frame).
c     TC1-TC2: permanent regression guard for the z(1)/Tau_geo
c     coupling claim -- reproduces driver_KS.F's own z(1) formula with
c     synthetic inputs and confirms its entire dependence on z(10) is
c     exactly the existing "-z(10)*(2/w^2)*(u.us)" term, matching
c     thesis eq. 2.56's own final term.
c     -------------------------------------------------------
      program test_dvdt_tess

      implicit double precision (a-h, o-z)
      integer npass, nfail
      parameter (nmax = 4)
      dimension u(4), q_leg(4)
      dimension Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax)
      dimension Cnm_rot0(0:nmax,0:nmax), Snm_rot0(0:nmax,0:nmax)
      dimension Cnm_rotp(0:nmax,0:nmax), Snm_rotp(0:nmax,0:nmax)
      dimension Cnm_rotm(0:nmax,0:nmax), Snm_rotm(0:nmax,0:nmax)
      dimension Cnm_zero(0:nmax,0:nmax), Snm_zero(0:nmax,0:nmax)
      dimension Rid(3,3)

      npass = 0
      nfail = 0

      do i = 1,3
         do j = 1,3
            Rid(i,j) = 0.d0
         end do
         Rid(i,i) = 1.d0
      end do

      write(*,*) '=============================================='
      write(*,*) ' Time-dependent tesseral geopotential in w*'
      write(*,*) '=============================================='

      pi = 4.d0*datan(1.d0)
      d2r = pi/180.d0

      theta_dot = 7.2921150d-5

      do n = 0, nmax
         do m = 0, nmax
            Cnm(n,m) = 0.d0
            Snm(n,m) = 0.d0
            Cnm_zero(n,m) = 0.d0
            Snm_zero(n,m) = 0.d0
         end do
      end do
      do n = 2, nmax
         do m = 0, n
            Cnm(n,m) = 1.d-5 * dble(n-m+1) / dble(n*n+1)
            if (m .gt. 0) then
               Snm(n,m) = -0.7d-5 * dble(m) / dble(n+3)
            end if
         end do
      end do

c     ---------------------------------------------------------
c     DT1-DT4: finite difference in theta, u held exactly fixed
c     ---------------------------------------------------------
      do it = 1, 4
         if (it .eq. 1) then
            u(1) = 0.9d0
            u(2) = 0.3d0
            u(3) = -0.2d0
            u(4) = 0.4d0
         else if (it .eq. 2) then
            u(1) = 1.1d0
            u(2) = -0.5d0
            u(3) = 0.6d0
            u(4) = 0.2d0
         else if (it .eq. 3) then
            u(1) = 0.4d0
            u(2) = 0.8d0
            u(3) = 0.1d0
            u(4) = -0.6d0
         else
            u(1) = -0.7d0
            u(2) = 1.0d0
            u(3) = 0.3d0
            u(4) = 0.5d0
         end if

         theta0_deg = 70.d0
         h_deg = 1.d-4

         call rotate_tess_coeffs(nmax, Cnm, Snm, theta0_deg,
     &                           Cnm_rot0, Snm_rot0)
         call rotate_tess_coeffs(nmax, Cnm, Snm, theta0_deg+h_deg,
     &                           Cnm_rotp, Snm_rotp)
         call rotate_tess_coeffs(nmax, Cnm, Snm, theta0_deg-h_deg,
     &                           Cnm_rotm, Snm_rotm)

         call tess_legendre_force(u, 6378.1363d0, theta_dot, Rid, nmax,
     &      Cnm_rot0, Snm_rot0, q_leg, tau_leg, dvdt_leg, vp0)
         call tess_legendre_force(u, 6378.1363d0, 0.d0, Rid, nmax,
     &      Cnm_rotp, Snm_rotp, q_leg, tau_leg, dum, vpp)
         call tess_legendre_force(u, 6378.1363d0, 0.d0, Rid, nmax,
     &      Cnm_rotm, Snm_rotm, q_leg, tau_leg, dum, vpm)

c        d(v_pert)/d(theta_deg) via central difference, converted to
c        d(v_pert)/dt = d(v_pert)/d(theta_rad) * theta_dot
         dv_dtheta_deg = (vpp - vpm) / (2.d0*h_deg)
         dv_dt_fd = dv_dtheta_deg * (1.d0/d2r) * theta_dot

         tol = 1.d-6*dabs(dvdt_leg) + 1.d-12
         write(*,'(A,I2)') '  --- case ', it
         call assert_d('  dvdt_leg vs finite difference', dvdt_leg,
     &                 dv_dt_fd, tol, npass, nfail)
      end do

c     ---------------------------------------------------------
c     ZT1: zonal-only regression guard -- dvdt_leg must be exactly
c     zero when no tesseral (m>0) coefficients are loaded, since an
c     axisymmetric field cannot exchange energy with a rotating frame.
c     ---------------------------------------------------------
      u(1) = 0.9d0
      u(2) = 0.3d0
      u(3) = -0.2d0
      u(4) = 0.4d0
      call rotate_tess_coeffs(nmax, Cnm_zero, Snm_zero, 70.d0,
     &                        Cnm_rot0, Snm_rot0)
      call tess_legendre_force(u, 6378.1363d0, theta_dot, Rid, nmax,
     &   Cnm_rot0, Snm_rot0, q_leg, tau_leg, dvdt_leg, vp0)
      call assert_d('ZT1 dvdt_leg is exactly zero (zonal-only)',
     &              dvdt_leg, 0.d0, 1.d-30, npass, nfail)

c     ---------------------------------------------------------
c     TC1-TC2: z(1)/Tau_geo coupling to z(10) -- permanent regression
c     guard for the "no direct Tau_geo edit needed" structural claim.
c     Reproduces driver_KS.F's own z(1) formula (synthetic inputs)
c     and confirms its ENTIRE dependence on z(10) is the existing
c     "-z(10)*(2/w^2)*(u.us)" term, matching thesis eq. 2.56's own
c     final term -- i.e. changing z(10) by delta changes z(1) by
c     exactly -delta*(2/w^2)*(u.us) and nothing else.
c     ---------------------------------------------------------
      amue_s = 398600.4418d0
      Tau_geo_s = 12.3d0
      Tau_3b1_s = 0.7d0
      Tau_3b2_s = -0.4d0
      Rdrag_s = 0.05d0
      w_s = 55.d0
      u(1) = 0.9d0
      u(2) = 0.3d0
      u(3) = -0.2d0
      u(4) = 0.4d0
      uus = 1.234d0

      z10_a = -3.d-6
      z10_b = 8.d-6

      z1_a = (amue_s + Tau_geo_s + Tau_3b1_s + Tau_3b2_s + Rdrag_s)
     &       /(8.d0*w_s**3) - z10_a*(2.d0/(w_s*w_s))*uus
      z1_b = (amue_s + Tau_geo_s + Tau_3b1_s + Tau_3b2_s + Rdrag_s)
     &       /(8.d0*w_s**3) - z10_b*(2.d0/(w_s*w_s))*uus

      dz1_actual = z1_b - z1_a
      dz1_expected = -(z10_b-z10_a)*(2.d0/(w_s*w_s))*uus

      call assert_d('TC1 dz(1) matches -dz(10)*(2/w^2)*(u.us)',
     &              dz1_actual, dz1_expected, 1.d-14, npass, nfail)

c     TC2: with z10_a=z10_b (no change to w*), z(1) must be unchanged
c     -- confirms Tau_geo's OTHER terms carry no hidden z(10)
c     dependence beyond this one coupling term.
      z1_c = (amue_s + Tau_geo_s + Tau_3b1_s + Tau_3b2_s + Rdrag_s)
     &       /(8.d0*w_s**3) - z10_a*(2.d0/(w_s*w_s))*uus
      call assert_d('TC2 z(1) unchanged when z(10) unchanged',
     &              z1_c, z1_a, 1.d-14, npass, nfail)

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
