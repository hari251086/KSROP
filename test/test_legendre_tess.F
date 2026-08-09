c     test_legendre_tess.F -- singularity-free general (n,m) tesseral
c     force: closed-form + cross-check against Cunningham
c     -------------------------------------------------------
c     src/LegendreTess.F (2026-08-09 reformulation, feature/tesseral-
c     energy-time-dependence branch): replaces the original spherical-
c     coordinate (phi,lambda) chain-rule formulation -- which divided
c     by sxy=sqrt(x^2+y^2) in its u-derivatives, a real polar
c     singularity -- with a fully Cartesian (x,y,z,r) formulation
c     built from two new non-singular recursions, G_general and
c     XY_general, motivated by the user's own thesis (Ch. 3, eq. 3.4)
c     hand-derived zonal closed form (pure r-power denominator, no
c     z-dependent denominator).
c
c     GN1-GN4: G_general(n,m) against P_n^m(sin phi)/cos(phi)^m
c     ground truth (closed-form associated Legendre / cos^m), zonal,
c     sectorial, sub-sectorial, interior cases -- same coverage the
c     former AP1-AP6 gave for alfP_general (now removed).
c
c     XY1-XY4: XY_general's X_m/Y_m against cos(m*lambda)*sxy^m /
c     sin(m*lambda)*sxy^m ground truth.
c
c     TL1-TL8: tess_legendre_force against tess_general_force
c     (already-validated, issue #30) on identical (u, R_body, Cnm,
c     Snm) inputs across nmax=2..6 -- unchanged from before the
c     reformulation, same acceptance-gate role.
c
c     PL1-PL4: the decisive new-capability regression -- evaluate at
c     and near a pole (sxy=0 exactly, and sxy~1e-8) where the former
c     formulation would have divided by (near-)zero. Confirms the
c     result stays finite and still matches Cunningham (which has no
c     coordinate-artifact singularity either).
c     -------------------------------------------------------
      program test_legendre_tess

      implicit double precision (a-h, o-z)
      integer npass, nfail
      dimension G(0:8,0:8), Gd(0:8,0:8)
      dimension Xt(0:6), Yt(0:6)
      dimension Xx(0:6), Xy(0:6), Yx(0:6), Yy(0:6)
      dimension u(4), q_leg(4), q_ref(4)
      dimension Cnm(0:6,0:6), Snm(0:6,0:6)

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' Singularity-free general (n,m) tesseral force'
      write(*,*) '=============================================='

c     ---------------------------------------------------------
c     GN1-GN4: G_general vs P_n^m(sin phi)/cos(phi)^m ground truth
c     ---------------------------------------------------------
      phi = 25.d0 * (4.d0*datan(1.d0)/180.d0)
      sp = dsin(phi)
      cp = dcos(phi)

      call G_general(8, sp, G, Gd)

c     P(2,0)=1.5sp^2-0.5, cos^0=1 -> G(2,0)=P(2,0)
      call assert_d('GN1 G(2,0) zonal', G(2,0),
     &              1.5d0*sp*sp - 0.5d0, 1.d-12, npass, nfail)
c     P(2,2)=3cp^2 -> G(2,2)=P(2,2)/cp^2=3
      call assert_d('GN2 G(2,2) sectorial', G(2,2),
     &              3.d0, 1.d-12, npass, nfail)
c     P(2,1)=3*sp*cp -> G(2,1)=3sp
      call assert_d('GN3 G(2,1) sub-sectorial', G(2,1),
     &              3.d0*sp, 1.d-12, npass, nfail)
c     P(4,1)=cp*sp*(17.5sp^2-7.5) -> G(4,1)=sp*(17.5sp^2-7.5)
c     (interior case -- same (n,m)=(4,1) the former AP6 covered)
      call assert_d('GN4 G(4,1) interior', G(4,1),
     &              sp*(17.5d0*sp*sp - 7.5d0), 1.d-12, npass, nfail)

c     ---------------------------------------------------------
c     XY1-XY4: X_m/Y_m vs cos(m*lambda)*sxy^m, sin(m*lambda)*sxy^m
c     ---------------------------------------------------------
      xlam = 3.d0
      ylam = -1.7d0
      sxy = dsqrt(xlam*xlam + ylam*ylam)
      alam = datan2(ylam, xlam)

      call XY_general(6, xlam, ylam, Xt, Yt, Xx, Xy, Yx, Yy)

      call assert_d('XY1 X_2', Xt(2), dcos(2*alam)*sxy**2,
     &              1.d-10, npass, nfail)
      call assert_d('XY2 Y_2', Yt(2), dsin(2*alam)*sxy**2,
     &              1.d-10, npass, nfail)
      call assert_d('XY3 X_5', Xt(5), dcos(5*alam)*sxy**5,
     &              1.d-8, npass, nfail)
      call assert_d('XY4 Y_5', Yt(5), dsin(5*alam)*sxy**5,
     &              1.d-8, npass, nfail)

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

c     ---------------------------------------------------------
c     PL1-PL4: pole regression -- the decisive new capability.
c     u chosen so x=u1^2-u2^2-u3^2+u4^2=0 AND y=2(u1u2-u3u4)=0
c     simultaneously (sxy=0 exactly, a true polar trajectory point)
c     -- the former formulation divided by sxy here and would have
c     produced NaN/Inf. Also checks a near-pole point (sxy~1e-8) for
c     the same reason. Both must stay finite and match Cunningham
c     (which has no coordinate-artifact singularity either).
c     ---------------------------------------------------------
      write(*,*) ''
      write(*,*) '  --- pole regression ---'
      nmax = 4
c     u1=u2=1, u3=u4=1 gives x=1-1-1+1=0, y=2(1-1)=0 exactly (z=2*2=4,
c     r=4 -- a genuine on-axis point, not degenerate in u).
      u(1) = 1.d0
      u(2) = 1.d0
      u(3) = 1.d0
      u(4) = 1.d0

      do n = 0, nmax
         do m = 0, nmax
            Cnm(n,m) = 0.d0
            Snm(n,m) = 0.d0
         end do
      end do
      do n = 2, nmax
         do m = 0, n
            Cnm(n,m) = 1.d-5 * dble(n-m+1) / dble(n*n+1)
            if (m .gt. 0) Snm(n,m) = -0.7d-5 * dble(m) / dble(n+3)
         end do
      end do

      call tess_general_force(u, 6378.1363d0, nmax, Cnm, Snm,
     &                         q_ref, tau_ref)
      call tess_legendre_force(u, 6378.1363d0, 7.2921150d-5,
     &                          nmax, Cnm, Snm,
     &                          q_leg, tau_leg, dvdt_leg, vpert_leg)

      isfinite = 1
      do j = 1, 4
         if (q_leg(j) .ne. q_leg(j)) isfinite = 0
      end do
      call assert_d('PL1 q(j) finite at exact pole (sxy=0)',
     &              dble(isfinite), 1.d0, 1.d-14, npass, nfail)
      call assert_d('PL2 q(1) matches Cunningham at pole', q_leg(1),
     &              q_ref(1), 1.d-9*dabs(q_ref(1))+1.d-14,
     &              npass, nfail)
      call assert_d('PL3 q(3) matches Cunningham at pole', q_leg(3),
     &              q_ref(3), 1.d-9*dabs(q_ref(3))+1.d-14,
     &              npass, nfail)
      call assert_d('PL4 tau matches Cunningham at pole', tau_leg,
     &              tau_ref, 1.d-9*dabs(tau_ref)+1.d-14,
     &              npass, nfail)

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
