c     test_drag.F  -  Unit tests for src/DragOblateCorotating.F
c     -------------------------------------------------------
c     Compile:  fpm test --compiler gfortran  (or ifx)
c     Run:      exit 0 = all pass, exit 1 = any fail
c     -------------------------------------------------------
c     2026-08-10: DragOblateCorotating.F was extracted from
c     driver_KS.F (previously inline-only) so KSROP-Mars can reuse
c     the identical oblate/co-rotating King-Hele drag force law
c     rather than forking it. These checks regression-test the three
c     shared subroutines against values independently computed in
c     Python from the SAME formulas that were previously inline in
c     driver_KS.F -- i.e. this proves the extraction is behavior-
c     preserving, on top of the (already-passing) end-to-end
c     Python integration tests (test_initial_conditions.py's
c     "+drag" regimes) that exercise the refactored driver_KS.F
c     itself.
c
c     Test case: a synthetic near-LEO Earth orbit (a=6800 km,
c     e=0.01, i=28.5 deg, argument of perigee=0). RHOPI_dg/H_dg are
c     deliberately round synthetic values (1e-11 kg/m3, 50 km) --
c     this test verifies the MATH, not a real atmosphere table.
c     -------------------------------------------------------
      program test_drag

      implicit double precision (a-h, o-z)
      dimension u(4), xd(3), qdrag(4)
      logical active, active2
      integer npass, nfail

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' DragOblateCorotating unit test suite'
      write(*,*) '=============================================='

      pi = 4.d0*datan(1.d0)
      d2r = pi/180.d0

c     ==============================================================
c     DA1-DA3: drag_perigee_active
c     ==============================================================
      A_dg     = 6.8d6
      EC_dg    = 0.01d0
      AI_dg    = 28.5d0*d2r
      SOMGA_dg = 0.d0*d2r
      R_body   = 6378.1363d0
      EPS_f    = 3.35281066d-3

      call drag_perigee_active(A_dg, EC_dg, AI_dg, SOMGA_dg, R_body,
     &                          EPS_f, 0.d0, 1.d6,
     &                          HSQ_dg, SI_dg, CSI_dg, active)

      call assert_d('DA1 HSQ_dg (perigee altitude, m)',
     &              HSQ_dg, 353863.7000000002d0, 1.d-4,
     &              npass, nfail)
      call assert_d('DA2 SI_dg = sin(28.5 deg)',
     &              SI_dg, 0.4771587602596084d0, 1.d-12,
     &              npass, nfail)
      call assert_d('DA3 CSI_dg = cos(28.5 deg)',
     &              CSI_dg, 0.8788171126619654d0, 1.d-12,
     &              npass, nfail)
      call assertl('DA4 active = .true. (in [0,1e6] m range)',
     &             active, npass, nfail)

      call drag_perigee_active(A_dg, EC_dg, AI_dg, SOMGA_dg, R_body,
     &                          EPS_f, 0.d0, 1.d4,
     &                          HSQ_dg2, SI_dg2, CSI_dg2, active2)
      call assertl('DA5 active = .false. (perigee above 10 km cap)',
     &             .not. active2, npass, nfail)

c     ==============================================================
c     DR1-DR3: drag_atm_reference
c     ==============================================================
      RHOPI_dg = 1.0d-11
      H_dg     = 50000.d0
      WE_rot   = 7.2921150d-5
      FR_rot   = 1.d0
      amue     = 398600.4415d0

      call drag_atm_reference(A_dg, EC_dg, SI_dg, CSI_dg, RHOPI_dg,
     &                         H_dg, WE_rot, FR_rot, amue, EPS_f,
     &                         F_dg, XXX_dg, C_dg)

      call assert_d('DR1 F_dg (co-rotation factor)',
     &              F_dg, 0.8915368326132731d0, 1.d-10,
     &              npass, nfail)
      call assert_d('DR2 XXX_dg (eccentricity exponent term)',
     &              XXX_dg, 1.36d0, 1.d-12,
     &              npass, nfail)
      call assert_d('DR3 C_dg (oblateness exponent term)',
     &              C_dg, 0.051390038023516306d0, 1.d-12,
     &              npass, nfail)

c     ==============================================================
c     DF1-DF4: drag_force_ks
c     ==============================================================
      u(1) = 2600.d0
      u(2) = 0.d0
      u(3) = 0.d0
      u(4) = 0.d0
      xd(1) = 0.d0
      xd(2) = 7.65d0
      xd(3) = 0.1d0
      V = 7.65d0
      EA_dg = 10.d0*d2r
      BN = 50.d0

      call drag_force_ks(u, xd, V, EC_dg, EA_dg, SOMGA_dg,
     &                    F_dg, C_dg, XXX_dg, RHOPI_dg, BN,
     &                    .true., 1.d0, qdrag)

      call assert_d('DF1 qdrag(1)',
     &              qdrag(1), 0.d0, 1.d-25, npass, nfail)
      call assert_d('DF2 qdrag(2)',
     &              qdrag(2), -1.3246148437312714d-15, 1.d-22,
     &              npass, nfail)
      call assert_d('DF3 qdrag(3)',
     &              qdrag(3), -1.7315226715441456d-17, 1.d-24,
     &              npass, nfail)
      call assert_d('DF4 qdrag(4)',
     &              qdrag(4), 0.d0, 1.d-25, npass, nfail)

c     ==============================================================
c     DF5-DF6: active=.false. must force qdrag exactly zero,
c     matching driver_KS.F's own IDRAG/range gate
c     ==============================================================
      call drag_force_ks(u, xd, V, EC_dg, EA_dg, SOMGA_dg,
     &                    F_dg, C_dg, XXX_dg, RHOPI_dg, BN,
     &                    .false., 1.d0, qdrag)
      call assert_d('DF5 qdrag(2) exactly zero when inactive',
     &              qdrag(2), 0.d0, 1.d-30, npass, nfail)
      call assert_d('DF6 qdrag(3) exactly zero when inactive',
     &              qdrag(3), 0.d0, 1.d-30, npass, nfail)

c     ==============================================================
c     Summary
c     ==============================================================
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

c     ============================================================
      subroutine assertl(label, cond, npass, nfail)
      character*(*) label
      logical cond
      integer npass, nfail

      if (cond) then
         write(*,10) 'PASS', label
         npass = npass + 1
      else
         write(*,10) 'FAIL', label
         nfail = nfail + 1
      end if
 10   format('  [',a4,']  ',a)
      return
      end
