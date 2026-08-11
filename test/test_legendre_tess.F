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
      double precision, allocatable :: Cnm(:,:), Snm(:,:)
c     Cnm/Snm are ALLOCATABLE, re-allocated to exactly (0:nmax,0:nmax)
c     at every site below where nmax changes -- NOT a fixed-size
c     array reused across different nmax values. tess_general_force/
c     tess_legendre_force declare their own Cnm/Snm dummy arguments as
c     explicit-shape Cnm(0:nmax,0:nmax); Fortran resolves that via
c     column-major sequence association against the ACTUAL argument's
c     linear memory, so if the actual array's leading dimension ever
c     differed from the callee's nmax+1 (e.g. a fixed (0:6,0:6) array
c     reused for nmax=2,3,4 calls), (n,m) indices would silently
c     alias the WRONG elements -- caught 2026-08-11 via intermittent
c     gfortran-only CI failures (ifx and fpm's own build happened not
c     to expose it): case nmax=3 read uninitialized memory that
c     resolved to NaN. driver_KS.F itself was never affected -- it
c     allocates Cnm_tess(0:ntess_use,0:ntess_use) fresh at the exact
c     size used every run.
      dimension Rid(3,3), Rrt(3,3)
      dimension u1v(4), up(4), um(4), usdum(4), qdum(4)
      dimension x0v(3), x1v(3), xd0v(3)
      character*48 label_rt

      npass = 0
      nfail = 0

      do i = 1,3
         do j = 1,3
            Rid(i,j) = 0.d0
         end do
         Rid(i,i) = 1.d0
      end do

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

         if (allocated(Cnm)) deallocate(Cnm, Snm)
         allocate(Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax))
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
         call tess_legendre_force(u, 6378.1363d0, 7.2921150d-5, Rid,
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

      if (allocated(Cnm)) deallocate(Cnm, Snm)
      allocate(Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax))
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
      call tess_legendre_force(u, 6378.1363d0, 7.2921150d-5, Rid,
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

c     ---------------------------------------------------------
c     RT1-RT4: Rtilt rotation covariance (2026-08-10, added when
c     Rtilt was introduced -- see LegendreTess.F header note). Two
c     independent checks that don't require deriving a KS-space
c     (u-space) rotation formula:
c
c     RT1 (forward rotation): v_pert is a true scalar (a potential
c     value at a physical point), so evaluating at u0 with Rtilt=R
c     must give EXACTLY the same v_pert as evaluating (with
c     Rtilt=identity) at u1, the KS pre-image of the ACTUALLY
c     rotated Cartesian point x1=R.x0 (obtained via car2ks, which
c     depends only on position, not velocity -- dummy xd/w are
c     safe). This directly tests that the forward xh,yh,sp rotation
c     is correct.
c
c     RT2-RT4 (backward/gradient rotation): finite-difference check
c     of q_leg's own implied dV/du_j (backed out from
c     q_leg(j)=v_pert*u_j+(r/2)*dVdu(j)) against a direct central
c     difference of v_pert itself as u_j is perturbed, at u0 with a
c     NON-trivial Rtilt=R. This tests that the gradient correctly
c     round-trips back through Rtilt^T before the du chain rule --
c     the same finite-difference acceptance gate this subroutine's
c     original 2026-08-09 derivation used (~1e-10, truncation-error
c     level).
c     ---------------------------------------------------------
      write(*,*) ''
      write(*,*) '  --- Rtilt rotation covariance ---'

c     Arbitrary well-conditioned rotation: 40deg about X after 30deg
c     about Z (R = R1(40) . R3(30), same convention as the Rtilt
c     header note's R1(90-dec0).R3(90+ra0) construction).
      ang1 = 30.d0 * (4.d0*datan(1.d0)/180.d0)
      ang2 = 40.d0 * (4.d0*datan(1.d0)/180.d0)
      c1 = dcos(ang1)
      s1 = dsin(ang1)
      c2 = dcos(ang2)
      s2 = dsin(ang2)
c     R3(ang1)
      r3_11=c1
      r3_12=s1
      r3_13=0.d0
      r3_21=-s1
      r3_22=c1
      r3_23=0.d0
      r3_31=0.d0
      r3_32=0.d0
      r3_33=1.d0
c     R1(ang2) . R3(ang1)
      rr_11 = r3_11
      rr_12 = r3_12
      rr_13 = r3_13
      rr_21 = c2*r3_21 + s2*r3_31
      rr_22 = c2*r3_22 + s2*r3_32
      rr_23 = c2*r3_23 + s2*r3_33
      rr_31 = -s2*r3_21 + c2*r3_31
      rr_32 = -s2*r3_22 + c2*r3_32
      rr_33 = -s2*r3_23 + c2*r3_33

      x0v(1) = 4200.d0
      x0v(2) = -1800.d0
      x0v(3) = 3100.d0
      xd0v(1) = 0.d0
      xd0v(2) = 0.d0
      xd0v(3) = 0.d0
      w0 = 1.d0

c     u depends only on position (car2ks), not velocity/w -- dummy
c     xd0v/w0 are safe (see LegendreTess.F Rtilt header note / this
c     block's own comment).
      call car2ks(x0v, xd0v, u, usdum, w0)

      x1v(1) = rr_11*x0v(1) + rr_12*x0v(2) + rr_13*x0v(3)
      x1v(2) = rr_21*x0v(1) + rr_22*x0v(2) + rr_23*x0v(3)
      x1v(3) = rr_31*x0v(1) + rr_32*x0v(2) + rr_33*x0v(3)

      call car2ks(x1v, xd0v, u1v, usdum, w0)

      nmax = 4
      if (allocated(Cnm)) deallocate(Cnm, Snm)
      allocate(Cnm(0:nmax,0:nmax), Snm(0:nmax,0:nmax))
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

      Rrt(1,1)=rr_11
      Rrt(1,2)=rr_12
      Rrt(1,3)=rr_13
      Rrt(2,1)=rr_21
      Rrt(2,2)=rr_22
      Rrt(2,3)=rr_23
      Rrt(3,1)=rr_31
      Rrt(3,2)=rr_32
      Rrt(3,3)=rr_33

      call tess_legendre_force(u, 6378.1363d0, 0.d0, Rrt, nmax,
     &                          Cnm, Snm, q_leg, tau_leg, dvdt_leg,
     &                          vpert_A)
      call tess_legendre_force(u1v, 6378.1363d0, 0.d0, Rid, nmax,
     &                          Cnm, Snm, q_ref, tau_ref, dvdt_dum,
     &                          vpert_B)

      call assert_d('RT1 v_pert(u0,Rtilt=R) = v_pert(R.x0,identity)',
     &              vpert_A, vpert_B, 1.d-9*dabs(vpert_B)+1.d-16,
     &              npass, nfail)

c     RT2-RT4: finite-difference dV/du_j vs q_leg's own implied value
      r_ks0 = u(1)*u(1)+u(2)*u(2)+u(3)*u(3)+u(4)*u(4)
      eps = 1.d-6
      do j = 1, 3
         up(1)=u(1)
         up(2)=u(2)
         up(3)=u(3)
         up(4)=u(4)
         um(1)=u(1)
         um(2)=u(2)
         um(3)=u(3)
         um(4)=u(4)
         up(j) = up(j) + eps
         um(j) = um(j) - eps
         call tess_legendre_force(up, 6378.1363d0, 0.d0, Rrt, nmax,
     &                             Cnm, Snm, qdum, taudum, dvdum,
     &                             vpert_p)
         call tess_legendre_force(um, 6378.1363d0, 0.d0, Rrt, nmax,
     &                             Cnm, Snm, qdum, taudum, dvdum,
     &                             vpert_m)
         dVdu_fd = (vpert_p - vpert_m) / (2.d0*eps)
         dVdu_analytic = (q_leg(j) - vpert_A*u(j)) * 2.d0 / r_ks0
         write(label_rt,'(A,I1,A)') 'RT', j+1,
     &        ' dV/du finite-diff vs analytic (post-Rtilt)'
         call assert_d(label_rt, dVdu_analytic, dVdu_fd,
     &                 1.d-6*dabs(dVdu_fd)+1.d-10, npass, nfail)
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
