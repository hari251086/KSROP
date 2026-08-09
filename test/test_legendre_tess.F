c     test_legendre_tess.F -- classical-Legendre general (n,m)
c     tesseral force: closed-form + cross-check against Cunningham
c     -------------------------------------------------------
c     Evaluates the derivation in EarthGravityPotential_KS (4).wl
c     (E:\Research\1. R&D\Mathematica), ported here as
c     src/LegendreTess.F, with two confirmed bugs fixed (see that
c     file's header comments): (1) the general associated Legendre
c     recursion, (2) the longitude convention.
c
c     AP1-AP6: alfP_general against closed-form associated Legendre
c     formulas (Rodrigues), covering m=0 (zonal), m=n (sectorial),
c     m=n-1 (sub-sectorial), and genuinely interior 0<m<n-1 cases
c     (AP5/AP6) -- the interior case is exactly where the source
c     .wl file's own recursion was wrong.
c
c     TL1-TL8: tess_legendre_force against tess_general_force
c     (already-validated, issue #30) on identical (u, R_body, Cnm,
c     Snm) inputs across nmax=2..6 -- the acceptance gate for
c     trusting this as a genuine independent cross-check, same
c     discipline as issue #30's own tess22_force acceptance gate.
c     -------------------------------------------------------
      program test_legendre_tess

      implicit double precision (a-h, o-z)
      integer npass, nfail
      dimension P(0:8,0:8)
      dimension u(4), q_leg(4), q_ref(4)
      dimension Cnm(0:6,0:6), Snm(0:6,0:6)

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' Classical-Legendre general (n,m) tesseral force'
      write(*,*) ' -- EarthGravityPotential_KS (4).wl port'
      write(*,*) '=============================================='

c     ---------------------------------------------------------
c     AP1-AP6: alfP_general vs closed-form Legendre formulas
c     ---------------------------------------------------------
      phi = 25.d0 * (4.d0*datan(1.d0)/180.d0)
      sp = dsin(phi)
      cp = dcos(phi)

      call alfP_general(8, sp, cp, P)

      call assert_d('AP1 P(2,0) zonal', P(2,0),
     &              1.5d0*sp*sp - 0.5d0, 1.d-12, npass, nfail)
      call assert_d('AP2 P(2,2) sectorial', P(2,2),
     &              3.d0*cp*cp, 1.d-12, npass, nfail)
      call assert_d('AP3 P(2,1) sub-sectorial', P(2,1),
     &              3.d0*sp*cp, 1.d-12, npass, nfail)
      call assert_d('AP4 P(3,2) sub-sectorial', P(3,2),
     &              15.d0*sp*cp*cp, 1.d-12, npass, nfail)
c     Interior (0<m<n-1) cases -- exactly where the source .wl
c     file's own recursion was wrong (verified: gave cp*(7.5sp^2
c     -2.5) instead of the correct cp*(7.5sp^2-1.5) at (3,1)).
      call assert_d('AP5 P(3,1) interior', P(3,1),
     &              cp*(7.5d0*sp*sp - 1.5d0), 1.d-12, npass, nfail)
      call assert_d('AP6 P(4,1) interior', P(4,1),
     &              cp*sp*(17.5d0*sp*sp - 7.5d0), 1.d-12,
     &              npass, nfail)

c     ---------------------------------------------------------
c     TL1-TL8: tess_legendre_force vs tess_general_force
c     ---------------------------------------------------------
      do it = 1, 4
         if (it .eq. 1) then
            nmax = 2
            u(1) = 0.9d0
            u(2) = 0.3d0
            u(3) = -0.2d0
            u(4) = 0.4d0
         else if (it .eq. 2) then
            nmax = 3
            u(1) = 1.1d0
            u(2) = -0.5d0
            u(3) = 0.6d0
            u(4) = 0.2d0
         else if (it .eq. 3) then
            nmax = 4
            u(1) = 0.4d0
            u(2) = 0.8d0
            u(3) = 0.1d0
            u(4) = -0.6d0
         else
            nmax = 6
            u(1) = -0.7d0
            u(2) = 1.0d0
            u(3) = 0.3d0
            u(4) = 0.5d0
         end if

         do n = 0, nmax
            do m = 0, nmax
               Cnm(n,m) = 0.d0
               Snm(n,m) = 0.d0
            end do
         end do
c        A representative spread of coefficients across the
c        triangle, including deliberately asymmetric magnitudes.
         do n = 2, nmax
            do m = 0, n
               Cnm(n,m) = 1.d-5 * dble(n-m+1) / dble(n*n+1)
               if (m .gt. 0) then
                  Snm(n,m) = -0.7d-5 * dble(m) / dble(n+3)
               end if
            end do
         end do

         call tess_general_force(u, 6378.1363d0, nmax, Cnm, Snm,
     &                            q_ref, tau_ref)
         call tess_legendre_force(u, 6378.1363d0, 7.2921150d-5,
     &                             nmax, Cnm, Snm,
     &                             q_leg, tau_leg, dvdt_leg, vpert_leg)

         write(*,*) ''
         write(*,'(A,I2,A,I1)') '  --- case ', it, ', nmax=', nmax
         call assert_d('  q(1) matches Cunningham', q_leg(1),
     &                 q_ref(1), 1.d-9*dabs(q_ref(1))+1.d-14,
     &                 npass, nfail)
         call assert_d('  q(2) matches Cunningham', q_leg(2),
     &                 q_ref(2), 1.d-9*dabs(q_ref(2))+1.d-14,
     &                 npass, nfail)
         call assert_d('  q(3) matches Cunningham', q_leg(3),
     &                 q_ref(3), 1.d-9*dabs(q_ref(3))+1.d-14,
     &                 npass, nfail)
         call assert_d('  q(4) matches Cunningham', q_leg(4),
     &                 q_ref(4), 1.d-9*dabs(q_ref(4))+1.d-14,
     &                 npass, nfail)
         call assert_d('  tau matches Cunningham', tau_leg,
     &                 tau_ref, 1.d-9*dabs(tau_ref)+1.d-14,
     &                 npass, nfail)
      end do

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
