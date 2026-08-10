c     ============================================================
c     DragOblateCorotating.F
c
c     Shared, body-agnostic implementation of KSROP's oblate,
c     co-rotating-atmosphere exponential drag model (King-Hele-style,
c     thesis eq. 5.2 / Sharma 1999):
c
c       rho = rho_p * exp{-ae/H*(1-cosE)
c             - eps*r_p*sin^2(I)/H*[sin^2(w+f)-sin^2(w)]}
c
c     where rho_p/H are the reference density/scale height at
c     perigee, eps is the CENTRAL BODY's flattening, f/E are true/
c     eccentric anomaly, and the atmosphere co-rotates with the body
c     at WE_rot*FR_rot.
c
c     Extracted 2026-08-10 from driver_KS.F (KSROP's own Earth
c     driver), where this was previously inline-only, so KSROP-Mars
c     can reuse the exact same force law rather than forking it --
c     mirrors how the geopotential force law was unified into
c     tess_legendre_force (2026-08-09) and reused unchanged by
c     KSROP-Lunar/KSROP-Mars. Nothing in this file assumes Earth:
c     every body constant (equatorial radius, flattening, rotation
c     rate) is an argument, not a literal.
c
c     Three subroutines, matching the three phases driver_KS.F's own
c     per-revolution/per-step loop already had:
c
c       drag_perigee_active : perigee altitude above the oblate body,
c                              and whether it falls inside the
c                              caller's atmosphere table's altitude
c                              range (evaluate once per revolution,
c                              BEFORE resolving a reference density).
c       drag_atm_reference   : given an ALREADY-RESOLVED reference
c                              density/scale-height at perigee (from
c                              INTPOL against a static table, or, for
c                              KSROP's own Earth driver only, an
c                              epoch-resolved space-weather lookup --
c                              this subroutine doesn't care which),
c                              the co-rotation and oblateness/
c                              eccentricity terms used every step.
c       drag_force_ks        : per-step local density (King-Hele
c                              exponential away from perigee) and the
c                              KS-space drag perturbation qdrag(4).
c
c     The diurnal density-bulge refinement (driven by KSROP's own
c     epoch-resolved space-weather subsystem, swx.F) has no Mars
c     analog in this codebase and stays driver-local: drag_force_ks
c     takes a plain bulge_mult multiplier (1.d0 = no bulge, the
c     Mars case; KSROP's own driver computes 1+f_bulge*cosphi_dg
c     itself and passes it in) rather than knowing about bulges at
c     all.
c     ============================================================

c     ============================================================
      subroutine drag_perigee_active(A_dg, EC_dg, AI_dg, SOMGA_dg,
     &                                R_body, EPS_f, ALT_lo, ALT_hi,
     &                                HSQ_dg, SI_dg, CSI_dg, active)
c     Perigee altitude above the oblate central body's reference
c     ellipsoid, and whether it falls inside [ALT_lo,ALT_hi] (the
c     caller's own atmosphere table bounds, in the same units as
c     HSQ_dg -- meters, matching driver_KS.F's own ALT_atm convention).
c
c     A_dg: semi-major axis (m). EC_dg: eccentricity (-). AI_dg,
c     SOMGA_dg: inclination, argument of perigee (rad). R_body:
c     central body equatorial radius (km). EPS_f: central body
c     flattening (-).
c     ============================================================
      implicit double precision (a-h,o-z)
      logical active

      SI_dg = dsin(AI_dg)
      CSI_dg = dcos(AI_dg)
      SIOMGA_dg = dsin(SOMGA_dg)

      RREQ_dg = R_body*1.d3*(1.d0 - SI_dg*SI_dg*SIOMGA_dg*SIOMGA_dg
     &          *EPS_f)
      SQ_dg = A_dg*(1.d0 - EC_dg)
      HSQ_dg = SQ_dg - RREQ_dg

      active = HSQ_dg .ge. ALT_lo .and. HSQ_dg .le. ALT_hi

      return
      end

c     ============================================================
      subroutine drag_atm_reference(A_dg, EC_dg, SI_dg, CSI_dg,
     &                               RHOPI_dg, H_dg, WE_rot, FR_rot,
     &                               amue, EPS_f, F_dg, XXX_dg, C_dg)
c     Per-revolution reference terms held fixed through every
c     integration step of the revolution: the co-rotating-atmosphere
c     velocity-correction factor F_dg, and the oblateness/
c     eccentricity exponent terms XXX_dg/C_dg the King-Hele exponential
c     needs every step.
c
c     RHOPI_dg (kg/m3) and H_dg (m) are the ALREADY-RESOLVED reference
c     density and scale height at this revolution's perigee altitude
c     -- resolve them however is appropriate for the body (a plain
c     INTPOL lookup against a static ALT/DEN/SCH table for both
c     KSROP-Mars and KSROP's own default path; KSROP's own Earth
c     driver additionally supports an epoch-resolved space-weather
c     override, resolved before calling this subroutine).
c     WE_rot: central body rotation rate (rad/s). FR_rot: atmosphere
c     co-rotation fraction (1.0 = full co-rotation). amue: central
c     body GM (km3/s2). EPS_f: central body flattening.
c     ============================================================
      implicit double precision (a-h,o-z)

      if (H_dg .le. 0.d0) H_dg = 1.d0
      BETA_dg = 1.d0/H_dg
      RPO_dg  = A_dg*(1.d0 - EC_dg)
      VPO_dg  = dsqrt(amue*1.d9*(2.d0/RPO_dg - 1.d0/A_dg))
      F_dg    = (1.d0 - RPO_dg*WE_rot*FR_rot*CSI_dg/VPO_dg)**2
      XXX_dg  = BETA_dg*A_dg*EC_dg
      C_dg    = EPS_f*BETA_dg*RPO_dg*SI_dg*SI_dg*0.5d0

      return
      end

c     ============================================================
      subroutine drag_force_ks(u, xd, V, EC_dg, EA_dg, SOMGA_dg,
     &                          F_dg, C_dg, XXX_dg, RHOPI_dg, BN,
     &                          active, bulge_mult, qdrag)
c     Per-step local density (King-Hele exponential away from
c     perigee) and the resulting KS-space drag perturbation
c     qdrag(4) = L(u)*p_drag.
c
c     u(4): KS elements. xd(3): inertial velocity (km/s). V: inertial
c     speed (km/s). EC_dg/EA_dg: eccentricity, eccentric anomaly
c     (rad) of the CURRENT stage state. SOMGA_dg/F_dg/C_dg/XXX_dg:
c     this revolution's reference terms (SOMGA_dg from the same
c     osculating state drag_perigee_active used; F_dg/C_dg/XXX_dg
c     from drag_atm_reference). RHOPI_dg: reference density at
c     perigee (kg/m3). BN: ballistic coefficient (kg/m2). active:
c     whether this revolution's perigee fell inside the atmosphere
c     table (from drag_perigee_active) -- qdrag is exactly zero when
c     false, matching driver_KS.F's own IDRAG/range gate. bulge_mult:
c     an external density multiplier (1.d0 = no correction; KSROP's
c     own Earth driver passes 1+f_bulge*cosphi_dg for its diurnal-
c     bulge refinement, which has no analog here and is not computed
c     by this subroutine).
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension u(4), xd(3), p_drag(3), qdrag(4)
      logical active

      if (active) then
         SEA_dg  = dsin(EA_dg)
         CEA_dg  = dcos(EA_dg)
         TER_dg  = 1.d0 - EC_dg*CEA_dg
         CSEF_dg = (CEA_dg - EC_dg)/TER_dg
         SIEF_dg = SEA_dg*dsqrt(1.d0 - EC_dg*EC_dg)/TER_dg
         FA_dg   = atan3(SIEF_dg,CSEF_dg)
         TER2_dg = dcos(SOMGA_dg)**2 - dcos(SOMGA_dg+FA_dg)**2
         exp_arg = -2.d0*C_dg*TER2_dg - XXX_dg*(1.d0-CEA_dg)
         if (exp_arg .gt. 500.d0 .or. exp_arg .lt. -500.d0) then
            TER_dg = 0.d0
         else
            TER_dg = dexp(exp_arg)
         end if

         DENS = -F_dg*RHOPI_dg*bulge_mult*TER_dg*V/(BN*2.d0)*1.d-7
         do i = 1,3
            p_drag(i) = DENS * xd(i)
         end do
      else
         do i = 1,3
            p_drag(i) = 0.d0
         end do
      end if

c     L(u) transformation: qdrag = L(u) * p_drag  (KS perturbation)
      qdrag(1) =  u(1)*p_drag(1)+u(2)*p_drag(2)+u(3)*p_drag(3)
      qdrag(2) = -u(2)*p_drag(1)+u(1)*p_drag(2)+u(4)*p_drag(3)
      qdrag(3) = -u(3)*p_drag(1)-u(4)*p_drag(2)+u(1)*p_drag(3)
      qdrag(4) =  u(4)*p_drag(1)-u(3)*p_drag(2)+u(2)*p_drag(3)

      return
      end
