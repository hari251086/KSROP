c     test_cunningham.F -- General (n,m) tesseral harmonic tests
c     -------------------------------------------------------
c     KSROP issue #30. Covers cunningham_Vnm (against the paper's
c     own Table I closed forms), cunningham_dVnm (against central
c     finite differences), geo_coeff_tess_general (against a
c     synthetic fixture matching geo_coeff_tess22's own known-good
c     EGM2008 (2,2) row), and tess_general_force (against the
c     already-validated tess22_force, at n=m=2 -- the acceptance
c     gate for trusting this general formula at every other (n,m)).
c     -------------------------------------------------------
      program test_cunningham

      implicit double precision (a-h, o-z)
      integer npass, nfail
      character*40 label
      parameter (nw = 8)
      dimension Vr(0:nw,0:nw), Vi(0:nw,0:nw)
      dimension u(4), q_gen(4), q_ref(4)
      dimension Cnm(0:4,0:4), Snm(0:4,0:4)
      dimension Cnm2(0:4,0:4), Snm2(0:4,0:4)

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' General (n,m) tesseral harmonic test suite'
      write(*,*) '=============================================='

c     ---------------------------------------------------------
c     CV1-CV9: cunningham_Vnm vs Table I closed forms (Cunningham
c     1970), test point x,y,z = 1.3,-0.7,2.1
c     ---------------------------------------------------------
      x = 1.3d0
      y = -0.7d0
      z = 2.1d0
      r2 = x*x + y*y + z*z
      r = dsqrt(r2)
      sigma2 = x*x + y*y

      call cunningham_Vnm(nw, x, y, z, r2, Vr, Vi)

      call assert_d('CV1 (0,0) Re', Vr(0,0)*r, 1.d0, 1.d-10,
     &              npass, nfail)
      call assert_d('CV2 (1,0) Re', Vr(1,0)*r**3, z, 1.d-10,
     &              npass, nfail)
      call assert_d('CV3 (1,1) Re', Vr(1,1)*r**3, x, 1.d-10,
     &              npass, nfail)
      call assert_d('CV4 (1,1) Im', Vi(1,1)*r**3, y, 1.d-10,
     &              npass, nfail)
      call assert_d('CV5 (2,0) Re', Vr(2,0)*r**5,
     &              (2.d0*z*z-sigma2)/2.d0, 1.d-9, npass, nfail)
      call assert_d('CV6 (2,1) Re', Vr(2,1)*r**5, 3.d0*z*x,
     &              1.d-9, npass, nfail)
      call assert_d('CV7 (2,2) Re', Vr(2,2)*r**5,
     &              3.d0*(x*x-y*y), 1.d-9, npass, nfail)
      call assert_d('CV8 (2,2) Im', Vi(2,2)*r**5,
     &              6.d0*x*y, 1.d-9, npass, nfail)
      call assert_d('CV9 (4,4) Re', Vr(4,4)*r**9,
     &              105.d0*((x*x-y*y)**2 - (2.d0*x*y)**2),
     &              1.d-6, npass, nfail)

c     ---------------------------------------------------------
c     CD1-CD3: cunningham_dVnm vs central finite difference
c     ---------------------------------------------------------
      h = 1.d-6

      do it = 1, 3
         if (it .eq. 1) then
            nn = 2
            mm = 2
         else if (it .eq. 2) then
            nn = 3
            mm = 1
         else
            nn = 4
            mm = 3
         end if

         call cunningham_Vnm(nw, x, y, z, r2, Vr, Vi)
         call cunningham_dVnm(nn, mm, nw, Vr, Vi,
     &                        dxr, dxi, dyr, dyi, dzr, dzi)

         r2p = (x+h)**2 + y*y + z*z
         r2m = (x-h)**2 + y*y + z*z
         call cunningham_Vnm(nw, x+h, y, z, r2p, Vr, Vi)
         vrp = Vr(nn,mm)
         call cunningham_Vnm(nw, x-h, y, z, r2m, Vr, Vi)
         vrm = Vr(nn,mm)
         fdx = (vrp - vrm)/(2.d0*h)

         write(label,'(A,I1,A,I1,A)') 'CDx n=',nn,' m=',mm,' Re'
         call assert_d(label, dxr, fdx, 1.d-6, npass, nfail)
      end do

c     ---------------------------------------------------------
c     CG1: geo_coeff_tess_general reproduces geo_coeff_tess22's
c     own already-validated (2,2) fixture exactly
c     ---------------------------------------------------------
      call geo_coeff_tess_general('test/fixture_geo_coeff.dat', 4,
     &                            Cnm, Snm)
      call geo_coeff_tess22('test/fixture_geo_coeff.dat', c22_ref,
     &                      s22_ref)
      call assert_d('CG1 C22 matches geo_coeff_tess22', Cnm(2,2),
     &              c22_ref, 1.d-14, npass, nfail)
      call assert_d('CG2 S22 matches geo_coeff_tess22', Snm(2,2),
     &              s22_ref, 1.d-14, npass, nfail)

c     ---------------------------------------------------------
c     CF1-CF5: tess_general_force at n=m=2 reproduces the
c     already-validated tess22_force exactly -- the acceptance
c     gate for the whole general method (KSROP#30)
c     ---------------------------------------------------------
      u(1) = 1.5d0
      u(2) = -1.d0/3.d0
      u(3) = 5.d0/7.d0
      u(4) = 0.4d0
      R_body = 5.5d0
      A_coef = -7.d0/3.d0
      B_coef = 4.d0/9.d0

      do n = 0, 4
         do m = 0, 4
            Cnm(n,m) = 0.d0
            Snm(n,m) = 0.d0
         end do
      end do
      Cnm(2,2) = A_coef
      Snm(2,2) = B_coef

      call tess_general_force(u, R_body, 4, Cnm, Snm, q_gen,
     &                         tau_gen)
      call tess22_force(u, R_body, A_coef, B_coef, q_ref, tau_ref)

      call assert_d('CF1 q(1) matches tess22_force', q_gen(1),
     &              q_ref(1), 1.d-9, npass, nfail)
      call assert_d('CF2 q(2) matches tess22_force', q_gen(2),
     &              q_ref(2), 1.d-9, npass, nfail)
      call assert_d('CF3 q(3) matches tess22_force', q_gen(3),
     &              q_ref(3), 1.d-9, npass, nfail)
      call assert_d('CF4 q(4) matches tess22_force', q_gen(4),
     &              q_ref(4), 1.d-9, npass, nfail)
      call assert_d('CF5 tau matches tess22_force', tau_gen,
     &              tau_ref, 1.d-9, npass, nfail)

c     ---------------------------------------------------------
c     CR1-CR4: rotate_tess_coeffs reduces exactly to the (2,2)-
c     specific A_tess/B_tess formula (issue #29) at m=2, and
c     correctly leaves m=0 (zonal) terms unrotated
c     ---------------------------------------------------------
      theta_test = 37.284d0
      c22raw = 1.7d0
      s22raw = -0.9d0
      c20raw = 3.3d0

      do n = 0, 4
         do m = 0, 4
            Cnm(n,m) = 0.d0
            Snm(n,m) = 0.d0
         end do
      end do
      Cnm(2,2) = c22raw
      Snm(2,2) = s22raw
      Cnm(2,0) = c20raw

      call rotate_tess_coeffs(4, Cnm, Snm, theta_test, Cnm2, Snm2)

      pi_t = 4.d0*datan(1.d0)
      d2r_t = pi_t/180.d0
      A_ref = c22raw*dcos(2.d0*theta_test*d2r_t)
     &      - s22raw*dsin(2.d0*theta_test*d2r_t)
      B_ref = c22raw*dsin(2.d0*theta_test*d2r_t)
     &      + s22raw*dcos(2.d0*theta_test*d2r_t)

      call assert_d('CR1 C22 rot matches A_tess formula', Cnm2(2,2),
     &              A_ref, 1.d-12, npass, nfail)
      call assert_d('CR2 S22 rot matches B_tess formula', Snm2(2,2),
     &              B_ref, 1.d-12, npass, nfail)
      call assert_d('CR3 C20 unrotated (m=0)', Cnm2(2,0), c20raw,
     &              1.d-14, npass, nfail)
      call assert_d('CR4 S20 unrotated (m=0)', Snm2(2,0), 0.d0,
     &              1.d-14, npass, nfail)

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
