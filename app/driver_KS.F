c     COPYRIGHT 2018, HARISHKUMAR SELLAMUTHU, ALL RIGHTS RESERVED
C     -------------------------------------------------------------
c     KS - Numerical driver
c       Description:
c     Orbit propagation using Kustaanheimo-Stiefel regular elements
c     with numerical integrator RKG4. (Earth)
c     Perturbations included are:
c     - Earth's oblateness ( up to 2190th degree )
c     - Luni-solar perturbations
c     - Atmospheric drag ( Oblate, co-rotating exponential
c       atmosphere referenced to perigee conditions, tabulated
c       density model via ATM.DAT -- same physical model and
c       logic as KSJLSDNP2.F )
c     - Solar radiation pressure (cannon-ball, KSJLSDNP2)
c     -------------------------------------------------------------
c     Author: Harishkumar Sellamuthu
c     hari251086@gmail.com
c     -------------------------------------------------------------
c     Input file parameteres <input.dat> by line
c     --------
c     Initial position in geocentric inertial reference frame in km
c     Initial velocity in geocentric inertial reference frame in km/s
c     Initial date in yyyy,mm,dd,hh,mm,ss (to be included)
c     Revolutions, steps, tolerance (1e-15 default)
c     Force model flags: n_force(1)=geo, n_force(2)=sun, n_force(3)=moon
c     Drag params: BN(kg/m2), IDRAG(0/1), WE_rot(rad/s), EPS_f, FR_rot
c     SRP: CR, AM(m2/kg), IPSR(0/1), ISHAD(0/1/2)
c     -------------------------------------------------------------c
c     Compiling with Intel Fortran compiler in command line:
c     ifort <driver_File>.F <Subroutines_File>.F
c     Execute with <driver_File> executable
c     -------------------------------------------------------------
c     ****Updates & Revisions****
c     --------
c     06-15-2018: New numerical program
c     Only J2
c     --------
c     06-16-2018: Trajectory dump at each revolution
c     Changing cn_j input to array,
c     Legendre Polynomial function from Numerical Recipes in F77
c     --------
c     09-13-2018: Added nth degree geopotential (Up to 2190 x 0)
c     EGM2008 model (Jn)
c     ngeo_deg value can be up to 2190 in the const_new.dat file
c     --------
c     Added atmospheric drag via oblate exponential model (ATM.DAT)
c     KS element equations updated for drag and energy correction
c     Fixed Tau_term -> Tau_geo in time element
c     -------------------------------------------------------------
c     ****Issues****
c     --------
c     Conversion from State to Keplerian elements needs to be
c     revisited for NaN outputs
c     -------------------------------------------------------------
c     Subroutines used
c     geo_coeff,car2ks,ks2car,aLegP,car2oe,ks2ksr,INTPOL
c     -------------------------------------------------------------c

      implicit double precision (a-h, o-z)
      integer cn0
      double precision iE
      parameter(n = 10, cn0 = 2190)
      parameter(ndim_atm = 291)
c     Value of pi, degree to radians and radians to degree
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
      dimension pek(7),x(3),xd(3),q(4),u(4),us(4),y(10),z(10)
      dimension qu(3),q_til(4),u_til(4),uu(4)
      dimension ZbyR(cn0*2),ObyR(cn0*2),Re(cn0)
      dimension cn_j(cn0),c_j(cn0),V_p(cn0),aleg(cn0)
      dimension sleg(cn0),oleg(cn0),R(n*2)
      dimension cal(6),ts(3),tm(3),opt(3),n_force(3)
      dimension RbyRs(cn0),RbyRm(cn0),cn_m(cn0),cn_s(cn0),Rm(n*2)
      dimension qsun(4),qmoon(4),slambda(cn0),olambda(cn0),Rs(n*2)
c     General (n,m) tesseral/mascon harmonics (issue #30, replacing
c     issue #29's (2,2)-only q22/tess22_force): q_gen is the
c     KS-element force contribution, amue-free like q(j) itself;
c     tau_gen likewise amue-free (amue applied explicitly when added
c     to Tau_geo, matching Tau_geo's own convention). Cnm_tess/
c     Snm_tess/Cnm_rot/Snm_rot are allocated once degree ntess_use is
c     known (below); a missing (n,m) row in the gravity file, or
c     ngeo_deg<=1, means that term's contribution is exactly zero --
c     the same auto-detect convention #29 already used.
c     ntess_cap bounds the general tesseral degree independently of
c     ngeo_deg (which can legitimately run to 2190 for zonal-only
c     precision): the Cunningham recursion is O(n^2) per step, not
c     O(n) like the zonal aLegP path, so reusing ngeo_deg directly
c     would make already-expensive high-degree zonal-only runs far
c     more expensive for no benefit -- real tesseral/mascon tables
c     (this issue's own Mars/Lunar motivation) don't need high
c     degree. Revisit this constant if a future consumer genuinely
c     needs more.
      integer ntess_cap, ntess_use
      parameter (ntess_cap = 10)
      dimension q_gen(4)
      double precision, allocatable :: Cnm_tess(:,:), Snm_tess(:,:)
      double precision, allocatable :: Cnm_rot(:,:), Snm_rot(:,:)
c     Drag: atmosphere tables and KS drag perturbation vectors
      dimension ALT_atm(ndim_atm),DEN_atm(ndim_atm),SCH_atm(ndim_atm)
      dimension qdrag(4),p_drag(3)
c     IDRAG: integer flag (0=off, 1=on); implicit i-n
      integer IDRAG
c     Drag model state (oblate, co-rotating exponential atmosphere
c     referenced to perigee conditions -- same physical model and
c     logic as KSJLSDNP2.F). The "_dg" quantities are the per-revolution
c     reference elements/parameters (held fixed through the revolution);
c     EA_dg (eccentric anomaly) and EC_dg (eccentricity) are refreshed
c     once per converged integration step.
      double precision A_dg,EC_dg,AI_dg,SOMGA_dg,EA_dg
      double precision SI_dg,CSI_dg,SIOMGA_dg,VIPP_dg,DE_dg
      double precision RREQ_dg,SQ_dg,HSQ_dg
      double precision RHOPI_dg,H_dg,BETA_dg,RPO_dg,VPO_dg
      double precision F_dg,XX_dg,XXX_dg,C_dg
      double precision SEA_dg,CEA_dg,TER_dg,TER2_dg
      double precision CSEF_dg,SIEF_dg,FA_dg,exp_arg
c     SRP: cannonball SRP perturbation (ported from KSJLSDNP2.F)
      dimension p_srp(3), qsrp(4), r_srp(3)
      integer IPSR, ISHAD
c     OEM/OPM epoch strings (YYYY-MM-DDTHH:MM:SS.sss = 23 chars)
      character*23 epochstr, epoch_stop, oem_start
c     Output filenames: output/KSROP_YYYYMMDDTHHMMSS.*
      character*40 oem_fname, reg_fname
      character*15 compact_epoch
c     Trajectory buffer: allocated once nrev is known; written after run
      double precision, allocatable :: traj_jd(:)
      double precision, allocatable :: traj_x(:,:)
      double precision, allocatable :: traj_xd(:,:)
c     Epoch-resolved space weather (ported from OREM issue #26):
c     auto-detected at startup, no config flag -- absent input files
c     mean the legacy static-table path runs unchanged.
      logical sw_exists, atm2d_exists
      integer ierr_sw, ierr_a2
      double precision tinf_wx
      integer ist_wx
c     Issue #32 Phase 5 follow-up (OREM): diurnal density bulge (only
c     computable when weather mode is active, since it needs the
c     atm2d T-to-density mapping). f_bulge refreshed once per
c     revolution; cosphi_dg refreshed every step since the
c     satellite's own position relative to the bulge direction
c     changes fast.
      double precision tmax_wx, tmin_wx, decl_sun_dg
      double precision rho_max_wx, rho_min_wx, h_dum_wx
      double precision f_bulge, cosphi_dg
      integer ist_bulge

      pi=4.d0*datan(1.d0)
      d2r=pi/180.d0
      r2d=180.d0/pi

c     Get UTC timestamp for output filenames
      call utc_now_epoch(epochstr, compact_epoch)
      reg_fname = 'output/KSROP_'//compact_epoch//'_Regular.out'

c     Input Files
      open (unit = 1, file = 'input/const_new.dat') !Constants
      open (unit = 2, file = 'input/input.dat')    !Sim params
      open (unit = 6, file = 'input/ATM.DAT')        !Atmosphere density
      open (unit = 8, file = 'input/input.opm')    !OPM input

c     Output files (unit 3 opened later with timestamped filename)
      open (unit = 4, file = reg_fname)            !KS elements (debug)
      open (unit = 5, file = 'output/ksrop.opm')  !OPM output

      read(1,*) amue,R_Earth,AU,amuS,amuM
      read(1,*) ngeo_deg,nsun_deg,nmoon_deg
      read(1,*) PSR_srp

c     Read initial state (position, velocity, epoch) from OPM file
      call read_opm(8, x, xd, cal)
      close(8)

c     Read simulation parameters from input.dat
      read(2,*,end=99) nrev,istep,tole !Revolutions, steps, tolerance
      read(2,*) n_force(1), n_force(2), n_force(3)
c     Drag parameters: ballistic coeff (kg/m2), drag on/off,
c     Earth rotation rate (rad/s), Earth flattening, rotation factor
      read(2,*) BN, IDRAG, WE_rot, EPS_f, FR_rot
c     SRP parameters: reflectivity coeff, area-to-mass (m^2/kg),
c     SRP on/off, shadow model (0=sunlit,1=cyl,2=cone)
      read(2,*) CR_srp, AM_srp, IPSR, ISHAD

      allocate(traj_jd(nrev*istep+1),
     &         traj_x(3,nrev*istep+1),
     &         traj_xd(3,nrev*istep+1))
      idump = 0

c     -----------------------------------------------------------
c     Atmosphere table: ALT in km (built here), SCH/DEN from ATM.DAT
c     Table covers 60 km to ~630 km. Converted to SI (m, kg/m3).
c     -----------------------------------------------------------
      ALT_atm(1) = 60.d0
      do i = 2,141
         ALT_atm(i) = ALT_atm(i-1) + 1.d0
      end do
      do i = 142,ndim_atm
         ALT_atm(i) = ALT_atm(i-1) + 2.d0
      end do

      read(6,110) (SCH_atm(i), i=1,61)
      read(6,110) (SCH_atm(i), i=62,ndim_atm)
      read(6,111) (DEN_atm(i), i=1,61)
      read(6,111) (DEN_atm(i), i=62,ndim_atm)
      close(6)

      do i = 1,ndim_atm
         ALT_atm(i) = ALT_atm(i) * 1000.d0   ! km to m
         SCH_atm(i) = SCH_atm(i) * 1000.d0   ! km to m
         DEN_atm(i) = DEN_atm(i) * 1.0d10    ! file units to kg/m3
      end do

c     -----------------------------------------------------------
c     Epoch-resolved space weather (ported from OREM issue #26):
c     auto-detect both input files, load if present, run the
c     legacy static ATM.DAT path unchanged if either is absent.
c     -----------------------------------------------------------
      inquire(file='input/SW-All.csv', exist=sw_exists)
      inquire(file='input/ATM2D.DAT', exist=atm2d_exists)
      if (sw_exists .and. atm2d_exists) then
         call sw_load('input/SW-All.csv', ierr_sw)
         call atm2d_load('input/ATM2D.DAT', ierr_a2)
         if (ierr_sw .eq. 0 .and. ierr_a2 .eq. 0) then
            write(*,'(A)') ' [SW] epoch-resolved density: ENABLED'
         else
            write(*,'(A)')
     &         ' [SW] epoch-resolved density: load failed, '//
     &         'using static ATM.DAT'
         end if
      else
         write(*,'(A)') ' [SW] epoch-resolved density: DISABLED '//
     &      '(input/SW-All.csv or input/ATM2D.DAT not found)'
      end if

      call car2oe(x,xd,pek)
      call oe2car(pek,x,xd,tole)
      call cal2jd(cal,djulian)
      call force_models(n_force,ngeo_deg,amuS,amuM)
      if (IPSR .eq. 0) PSR_srp = 0.d0
      SRPC_srp = CR_srp * AM_srp

      print *, x, xd

      write(*,*) '=============================================='

c     General (n,m) tesseral/mascon coefficients (issue #30):
c     ntess_use<2 (from ngeo_deg<=1, or a large ngeo_deg clipped by
c     ntess_cap only on its own -- ntess_use is always >=1 whenever
c     ngeo_deg>1) means geo_coeff_tess_general's own load condition
c     (ideg>=2) can never match, leaving these arrays all-zero --
c     the same genuine auto-detect #29 already relied on, not a
c     hard dependency on the row existing.
      ntess_use = min(ngeo_deg, ntess_cap)
      allocate(Cnm_tess(0:ntess_use,0:ntess_use),
     &         Snm_tess(0:ntess_use,0:ntess_use),
     &         Cnm_rot(0:ntess_use,0:ntess_use),
     &         Snm_rot(0:ntess_use,0:ntess_use))

      if(ngeo_deg.eq.0.or.ngeo_deg.eq.1) then
         do i = 1,cn0
            c_j(i) = 0.d0
         end do
c        No EGM2008 file guaranteed present in this regime (same
c        reasoning as the zonal branch above) -- zero directly
c        rather than calling geo_coeff_tess_general, which would
c        open the file unconditionally.
         do i = 0,ntess_use
            do j = 0,ntess_use
               Cnm_tess(i,j) = 0.d0
               Snm_tess(i,j) = 0.d0
            end do
         end do
      else
!     compute j_n values
         call geo_coeff(ngeo_deg,c_j)
         call geo_coeff_tess_general('input/EGM2008_to2190_TideFree',
     &        ntess_use, Cnm_tess, Snm_tess)
      end if

!     Initialization
      E_0 = 0.d0                !Gen. Eccentric Anomaly in radians
      E = E_0
      T = 0.d0                  !Physical time in seconds
      nt = 0                    !Counter in the integrator

      dE_0 = 2.d0*pi/istep      !Step-size in radians

      num_steps = nrev*istep

      call solarnpv(djulian,ts)
      call lunarpv(djulian,tm)

c     General (n,m) tesseral (issue #30): Earth's rotation angle at
c     this epoch, and the resulting time-dependent (not position-
c     dependent) rotated coefficient table tess_general_force needs
c     -- see gmst_deg/rotate_tess_coeffs/tess_general_force header
c     comments.
      call gmst_deg(djulian, theta_tess)
      call rotate_tess_coeffs(ntess_use, Cnm_tess, Snm_tess,
     &                        theta_tess, Cnm_rot, Snm_rot)

      R(1) = vmn(x)
      V = vmn(xd)

      R(2) = R(1)*R(1)
      V2 = V*V

      do i = 2, 20
         R(i) = R(i-1) * R(1)
      end do

!     Perturbing potential, V_pot

!     Oblateness terms
      Re(1) = R_Earth
      ZbyR(1) = x(3)/R(1)
      ObyR(1) = 1.d0/R(1)

!     Third-body terms
      call third_body_aux(nsun_deg,  x, ts, R(1), amuS,
     &                    Rs(1), RbyRs, first_s, cphi_s, cphi_s2)
      call third_body_aux(nmoon_deg, x, tm, R(1), amuM,
     &                    Rm(1), RbyRm, first_m, cphi_m, cphi_m2)

      do i = 1,ngeo_deg
         ObyR(i+1) = ObyR(i)*ObyR(1)
         Re(i+1) = Re(i)*Re(1)
      end do

      V_pot = 0.d0

c     +1: the force/time-element formulas below need aleg(i+2)
c     (degree i+1) for i up to ngeo_deg, i.e. one degree beyond
c     the requested ngeo_deg itself (same for sleg/oleg).
      call aLegP(ngeo_deg+1,ZbyR(1),aleg)
      call aLegP(nsun_deg+1,cphi_s,sleg)
      call aLegP(nmoon_deg+1,cphi_m,oleg)

      do i = 2,ngeo_deg
         cn_j(i) = amue*Re(i)*c_j(i)*ObyR(i+1)
         V_p(i) = cn_j(i)*aleg(i+1)
         V_pot = V_pot+V_p(i)
      end do

      do i = 2,nsun_deg
         cn_s(i) = first_s * RbyRs(i)
         V_p(i) = cn_s(i) * sleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

      do i = 2,nmoon_deg
         cn_m(i) = first_m * RbyRm(i)
         V_p(i) = cn_m(i) * oleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

!     Angular frequency
      w = dsqrt(0.5d0*(amue/R(1)-V2/2.d0-V_pot))
      w_kep = dsqrt(0.5d0*(amue/R(1)-V2/2.d0))
      Gam = w/w_kep

      dE = dE_0*Gam

      call car2ks(x,xd,u,us,w)
      call ks2car(u,us,x,xd,w)

      Tow = T+dotp4(u,us)/w     !Time-element

!     Initial KS regular elements
      y(1) = Tow
      y(10) = w

      do i = 1,4
         y(i+1) = u(i)
         y(i+5) = 2.d0*us(i)
      end do

      sie = dsin(E/2.d0)
      cse = dcos(E/2.d0)
      iE = dint(E)

      call car2oe(x,xd,pek)
      call jd2epoch(djulian, epochstr)

c     --- CCSDS OPM v2.0 (initial epoch) — written immediately ---
      call write_opm(5, epochstr, x, xd, pek)

c     --- Buffer initial state point (OEM written after run) ---
      idump = idump + 1
      traj_jd(idump) = djulian + T/86400.d0
      do i = 1,3
         traj_x(i,idump)  = x(i)
         traj_xd(i,idump) = xd(i)
      end do
      write(4,101) iE,y

c     Integration
      do ik = 1,nrev            !Revolutions

c     -----------------------------------------------------------
c     Drag: per-revolution reference setup
c     Oblate, co-rotating exponential atmosphere referenced to
c     perigee conditions (same physical model as KSJLSDNP2.F).
c
c     MODEL: thesis eq. 5.2 (Sharma 1999):
c       rho = rho_p * exp{-ae/H*(1-cosE)
c             - eps*r_p*sin^2(I)/H*[sin^2(w+f)-sin^2(w)]}
c     where H = scale height from ATM.DAT, eps = Earth
c     flattening, f = true anomaly, E = eccentric anomaly.
c
c     NOTE: the advanced model (thesis eq. 5.3) with density
c     scale height variation (R_i coefficients i=1..8, from
c     Sharma & Raj 2009/2010) is NOT implemented. The current
c     model uses constant H per revolution. For the R_i model
c     see thesis pp. 111-113 and Nair & Sharma (2003).
c     The basic model is sufficient for re-entry prediction
c     with < 5% RPE as demonstrated in the manuscript.
c     The reference orbital elements and atmosphere/rotation
c     parameters below are evaluated once per revolution from the
c     osculating elements at the start of the revolution, and held
c     fixed through all of its integration steps.
c     -----------------------------------------------------------
         call car2oe(x,xd,pek)
         A_dg     = pek(1)*1.d3       ! semi-major axis: km -> m
         EC_dg    = pek(2)
         AI_dg    = pek(3)*d2r        ! inclination: deg -> rad
         SOMGA_dg = pek(5)*d2r        ! argument of perigee: deg -> rad

c     Issue #24 (ported from OREM v1.18): the eccentric anomaly used
c     by the drag density is read from the true state (pek(7),
c     refreshed by the per-stage car2oe below), not advanced
c     analytically. The old sweep DE_dg = (VIPP*pi - EA0)/istep
c     half-rates the density phase whenever a revolution starts past
c     EA=pi (VIPP=4 targets 4pi over one rev's steps), intermittently
c     dephasing the density peak from the true perigee passage along
c     long decay arcs -- a geometry-dependent error a fitted BN
c     cannot absorb. Bit-identical on perigee-anchored windows;
c     on OREM's re-entry arcs it moved RPE from -70..-97% to
c     bracketing zero. Validated ~1% vs exact RK4 integration of the
c     same model (OREM/scratch_gmat/drag_ref.py).

         SI_dg     = dsin(AI_dg)
         CSI_dg    = dcos(AI_dg)
         SIOMGA_dg = dsin(SOMGA_dg)

c     Perigee altitude above the oblate Earth (drag is active for
c     the revolution only when this is below 500 km)
         RREQ_dg = R_Earth*1.d3*(1.d0 - SI_dg*SI_dg*SIOMGA_dg*SIOMGA_dg
     &             *EPS_f)
         SQ_dg  = A_dg*(1.d0 - EC_dg)
         HSQ_dg = SQ_dg - RREQ_dg

         if (HSQ_dg .ge. ALT_atm(1) .and.
     &       HSQ_dg .le. ALT_atm(ndim_atm)) then
c     Reference density/scale-height at perigee, the rotating-
c     atmosphere correction factor F_dg, and the oblateness/
c     eccentricity exponent terms C_dg, XXX_dg used each step.
c     Ported from OREM issue #26: when a space-weather history and
c     2-D atmosphere table are loaded (auto-detected at startup),
c     the reference is evaluated at the exospheric temperature for
c     THIS revolution's epoch instead of the static table. With
c     nothing loaded (the default), sw_tinf returns 0 and the
c     legacy path below runs unchanged.
            call sw_tinf(djulian + T/86400.d0, tinf_wx)
            ist_wx = 2
            if (tinf_wx .gt. 0.d0)
     &         call atm2d_interp(HSQ_dg, tinf_wx,
     &                           RHOPI_dg, H_dg, ist_wx)
            if (ist_wx .ne. 0) then
               call INTPOL(ALT_atm, DEN_atm, ndim_atm,
     &                     HSQ_dg, RHOPI_dg)
               call INTPOL(ALT_atm, SCH_atm, ndim_atm,
     &                     HSQ_dg, H_dg)
            end if
            if (H_dg .le. 0.d0) H_dg = 1.d0

c           Issue #32 Phase 5 follow-up (OREM): diurnal density
c           bulge amplitude, only derivable when weather mode
c           succeeded (needs the atm2d T-to-density mapping at this
c           same perigee altitude, from the same table already
c           loaded). Inactive (f_bulge=0, a pure no-op) whenever
c           weather mode is off or either extremum lookup falls
c           outside the loaded table -- same auto-detect philosophy
c           as issue #26 itself. Refreshed once per revolution (Sun
c           declination barely changes within one orbit); the
c           geometric factor that multiplies it (cosphi_dg) is
c           refreshed every step.
            f_bulge = 0.d0
            if (tinf_wx .gt. 0.d0 .and. ist_wx .eq. 0) then
               decl_sun_dg = dasin(ts(3)/dsqrt(ts(1)**2+ts(2)**2
     &                       +ts(3)**2)) * r2d
               call jacchia_bulge_extrema(tinf_wx, decl_sun_dg,
     &                                    tmax_wx, tmin_wx)
               call atm2d_interp(HSQ_dg, tmax_wx,
     &                           rho_max_wx, h_dum_wx, ist_bulge)
               if (ist_bulge .eq. 0) then
                  call atm2d_interp(HSQ_dg, tmin_wx,
     &                              rho_min_wx, h_dum_wx, ist_bulge)
                  if (ist_bulge .eq. 0 .and. rho_min_wx .gt. 0.d0)
     &               f_bulge = (rho_max_wx/rho_min_wx - 1.d0)
     &                       / (rho_max_wx/rho_min_wx + 1.d0)
               end if
            end if

            BETA_dg = 1.d0/H_dg
            RPO_dg  = SQ_dg
            VPO_dg  = dsqrt(amue*1.d9*(2.d0/RPO_dg - 1.d0/A_dg))
            F_dg    = (1.d0 - RPO_dg*WE_rot*FR_rot*CSI_dg/VPO_dg)**2
            XX_dg   = A_dg*EC_dg
            XXX_dg  = BETA_dg*XX_dg
            C_dg    = EPS_f*BETA_dg*RPO_dg*SI_dg*SI_dg*0.5d0
         end if

         do ki = 1,istep        !Steps
 60         continue
            call u2uu(u,uu)

c     -----------------------------------------------------------
c     Atmospheric drag perturbation
c     Local density = perigee-reference density (RHOPI_dg) times
c     an exponential correction (TER_dg) for the satellite's true
c     anomaly, orbital oblateness/eccentricity (C_dg, XXX_dg), and
c     a co-rotating-atmosphere correction (F_dg) -- same physical
c     model and logic as KSJLSDNP2.F.
c     Compute drag acceleration p_drag(3) in km/s^2, then
c     transform to KS space via the L(u) matrix: qdrag(4) = L(u)*p_drag
c     -----------------------------------------------------------
            if (IDRAG .eq. 1 .and.
     &          HSQ_dg .ge. ALT_atm(1) .and.
     &          HSQ_dg .le. ALT_atm(ndim_atm)) then
c              Issue #24: density at the TRUE eccentric anomaly of
c              the current (stage) state -- stage-accurate, covers
c              both wings of the perigee density bump.
               EA_dg   = pek(7)*d2r
               EC_dg   = pek(2)
               SEA_dg  = dsin(EA_dg)
               CEA_dg  = dcos(EA_dg)
               TER_dg  = 1.d0 - EC_dg*CEA_dg
               CSEF_dg = (CEA_dg - EC_dg)/TER_dg
               SIEF_dg = SEA_dg*dsqrt(1.d0 - EC_dg*EC_dg)/TER_dg
               FA_dg   = atan3(SIEF_dg,CSEF_dg)
               TER2_dg = dcos(SOMGA_dg)**2
     &                 - dcos(SOMGA_dg+FA_dg)**2
               exp_arg = -2.d0*C_dg*TER2_dg
     &                  - XXX_dg*(1.d0-CEA_dg)
               if (exp_arg .gt. 500.d0 .or.
     &             exp_arg .lt. -500.d0) then
                  TER_dg = 0.d0
               else
                  TER_dg = dexp(exp_arg)
               end if
c              Issue #32 Phase 5 follow-up (OREM): diurnal density
c              bulge -- geometric factor evaluated at the
c              satellite's CURRENT stage position (unlike f_bulge,
c              this genuinely varies fast within one revolution,
c              since that's the whole point of a day/night bulge).
c              f_bulge=0 (weather mode off, the common case) makes
c              this an exact no-op regardless of cosphi_dg's value.
               cosphi_dg = 0.d0
               if (f_bulge .ne. 0.d0)
     &            call bulge_cosphi(x, ts, 37.d0, cosphi_dg)

               DENS = -F_dg*RHOPI_dg*(1.d0 + f_bulge*cosphi_dg)
     &              *TER_dg*V/(BN*2.d0)*1.d-7
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

c     -----------------------------------------------------------
c     SRP: cannonball solar radiation pressure perturbation.
c     Shadow factor ROHP: 1=sunlit, 0=eclipse.  ISHAD=0 -> no shadow
c     model (always sunlit).  Acceleration in km/s^2 via 1e-3 scaling.
c     -----------------------------------------------------------
            if (PSR_srp .eq. 0.d0) then
               do i = 1,4
                  qsrp(i) = 0.d0
               end do
            else
               if (ISHAD .eq. 1) then
                  ROHP = shadfncyl(x(1),x(2),x(3),ts(1),ts(2),ts(3))
               else if (ISHAD .eq. 2) then
                  ROHP = shadfncone(x(1),x(2),x(3),ts(1),ts(2),ts(3))
               else
                  ROHP = 1.d0
               end if
               r_srp(1) = ts(1) - x(1)
               r_srp(2) = ts(2) - x(2)
               r_srp(3) = ts(3) - x(3)
               xsrr = vmn(r_srp)
               do i = 1,3
                  p_srp(i) = -ROHP*PSR_srp*SRPC_srp*r_srp(i)/xsrr*1.d-3
               end do
               qsrp(1) =  u(1)*p_srp(1)+u(2)*p_srp(2)+u(3)*p_srp(3)
               qsrp(2) = -u(2)*p_srp(1)+u(1)*p_srp(2)+u(4)*p_srp(3)
               qsrp(3) = -u(3)*p_srp(1)-u(4)*p_srp(2)+u(1)*p_srp(3)
               qsrp(4) =  u(4)*p_srp(1)-u(3)*p_srp(2)+u(2)*p_srp(3)
            end if

c     Time-element rate (KS regular elements eq., Sellamuthu 2018 thesis
c     eq. 2.56): tau* = [mu - 2rU]/(8w^3) - r/(16w^3)*[u.dU/du] + ...
c     which reduces, per zonal term, to +amue*(n-1)*c_j(n)*Re(n)*
c     aleg(n+1)*ObyR(n) (matches this loop's structure once amue is
c     included and the sign is flipped from the pre-2026-07-11 formula,
c     which omitted amue entirely -- verified against the thesis closed
c     form to 9 significant digits).
            Tau_geo = 0.d0
            do i = 2,ngeo_deg
               Tau_geo = Tau_geo+amue*(i-1)*c_j(i)*Re(i)*aleg(i+1)
     &                   *ObyR(i)
            end do

c     General (n,m) tesseral contribution (issue #30): q_gen/tau_gen
c     are amue-free, like q(j)/Tau_geo's own zonal accumulation, so
c     amue is applied here exactly as it already is for Tau_geo
c     above. All-zero Cnm_rot/Snm_rot (no tesseral data loaded) makes
c     q_gen and tau_gen identically zero -- see tess_legendre_force/
c     gmst_deg. The classical-Legendre derivation (src/LegendreTess.F)
c     is the sole active geopotential here -- verified to reproduce
c     Cunningham's independent Cartesian solid-harmonic recursion
c     (src/Cunningham.F, kept in the library as a standing cross-check,
c     not called from this loop) to machine precision at the force-
c     vector level and byte-for-byte in propagated output.
            call tess_legendre_force(u, R_Earth, ntess_use,
     &                                Cnm_rot, Snm_rot,
     &                                q_gen, tau_gen)
            Tau_geo = Tau_geo + amue*tau_gen

            Tau_3body1 = 0.d0
            do i = 2,nsun_deg
               Tau_3body1=Tau_3body1
     &              +(i+2)*sleg(i+1)*RbyRs(i)*RbyRs(1)*amuS
            end do
            Tau_3body2 = 0.d0
            do i = 2,nmoon_deg
               Tau_3body2=Tau_3body2
     &              +(i+2)*oleg(i+1)*RbyRm(i)*RbyRm(1)*amuM
            end do

c     -----------------------------------------------------------
c     Energy element: non-conservative drag + SRP contributions
c     dw/ds = -0.5 * (us . (qdrag + qsrp)) / w
c     -----------------------------------------------------------
            z(10) = -0.5d0*(us(1)*(qdrag(1)+qsrp(1))
     &                    + us(2)*(qdrag(2)+qsrp(2))
     &                    + us(3)*(qdrag(3)+qsrp(3))
     &                    + us(4)*(qdrag(4)+qsrp(4)))/w

c     -----------------------------------------------------------
c     Time element: conservative + drag + SRP contributions
c     Tau_term replaced by Tau_geo (fixes uninitialized variable)
c     Non-conservative: drag + SRP energy correction
c     -----------------------------------------------------------
            z(1) = (amue + Tau_geo + Tau_3body1 + Tau_3body2
     &             + R(1)*(dotp4(u,qdrag)+dotp4(u,qsrp)))/(8.d0*w*w*w)
     &             - z(10)*(2.d0/(w*w))*dotp4(u,us)

            if (ngeo_deg  .gt. 1) call aLegP(ngeo_deg+1, ZbyR(1), aleg)
            if (nsun_deg  .gt. 1) call aLegP(nsun_deg+1, cphi_s,  sleg)
            if (nmoon_deg .gt. 1) call aLegP(nmoon_deg+1,cphi_m,  oleg)

            do j = 1,4
               q(j) = 0.d0
               qsun(j) = 0.d0
               qmoon(j) = 0.d0
            end do

            do j = 1,4

!     oblateness terms
               z2mR2=x(3)**2-R(1)**2
               qj=0.d0
               do i = 2,ngeo_deg
                  qj = dble(ObyR(i+1)*Re(i)*c_j(i)*
     &                 (R(1)*(i+1)*aleg(i+2)*(R(1)*uu(j)-x(3)*u(j))+
     &                 (aleg(i+1)*(u(j)*(i*R(2)+x(3)**2)
     &                 -uu(j)*(i+1)*R(1)*x(3))))/z2mR2)
                  q(j)=qj+q(j)
               end do

c     General (n,m) tesseral contribution, computed once (all 4
c     components together) above -- amue-free like the zonal q(j)
c     sum, so it is picked up automatically by the same
c     "q(j)*amue" in z(j+1)/z(j+5) below.
               q(j) = q(j) + q_gen(j)

!     third-body terms
               call u2qu(u,qu,j)

               den_s = cphi_s2 - 1.d0
               den_m = cphi_m2 - 1.d0

               sq = ts(1)*qu(1)+ts(2)*qu(2)+ts(3)*qu(3)

c     Third-body q(j) in the KS-elements EOM convention (thesis eq.
c     2.55/2.59): the paren consumed by z(j+1)/z(j+5) equals
c     2*[U/2*u + (r/4)*dU/du], so per degree i (shape=(r/RK)^i*P_i):
c       qsun(j) = RbyRs(i)*[ P_i*u(j) - ( i*(P_i - cphi*P_(i-1))*u(j)
c                 + (i+1)*(cphi*P_i - P_(i+1))*sq/Rs )/den_s ]
c     Verified 2026-07-12: matches the independently validated
c     KSJLSDNP.F n=2 sun term (CNTS1/CNTS2 form) to machine precision
c     and a finite-difference gradient of shape*u + r*L(u)^T(grad)
c     to ~1e-8 at n=2,3. The 2026-07-11 version (commit 963968f) used
c     the wrong convention, L(u)^T(-grad shape) alone -- missing the
c     factor r and the U*u term, ~1e4 too small at GTO scale.
               do i = 2,nsun_deg
                  qs = RbyRs(i)*(sleg(i+1)*u(j)
     &                 -(dble(i)*(sleg(i+1)-cphi_s*sleg(i))*u(j)
     &                 +dble(i+1)*(cphi_s*sleg(i+1)-sleg(i+2))*sq
     &                 /Rs(1))/den_s)
                  qsun(j) = qs + qsun(j)
               end do

               oq = tm(1)*qu(1)+tm(2)*qu(2)+tm(3)*qu(3)

               do i = 2,nmoon_deg
                  qm = RbyRm(i)*(oleg(i+1)*u(j)
     &                 -(dble(i)*(oleg(i+1)-cphi_m*oleg(i))*u(j)
     &                 +dble(i+1)*(cphi_m*oleg(i+1)-oleg(i+2))*oq
     &                 /Rm(1))/den_m)
                  qmoon(j) = qm + qmoon(j)
               end do

c     -----------------------------------------------------------
c     State: conservative + drag + SRP, energy corr
c     z(j+1) = [geo + 3body - R*(drag+srp)]*sin(E/2)*0.25/w^2
c              + 2*z(10)*us(j)*sin(E/2)/w
c     z(j+5) = -[geo + 3body - R*(drag+srp)]*cos(E/2)*0.25/w^2
c              - 2*z(10)*us(j)*cos(E/2)/w
c     -----------------------------------------------------------
               aqt = sie * 0.25d0
               bqt = cse * 0.25d0

               z(j+1) = (q(j)*amue + first_s*qsun(j)
     &              + first_m*qmoon(j)
     &              - R(1)*(qdrag(j)+qsrp(j)))
     &              * aqt/(w*w)
     &              + 2.d0*z(10)*us(j)*sie/w

               z(j+5) = -(q(j)*amue + first_s*qsun(j)
     &              + first_m*qmoon(j)
     &              - R(1)*(qdrag(j)+qsrp(j)))
     &              * bqt/(w*w)
     &              - 2.d0*z(10)*us(j)*cse/w

            end do

            I_rkg = rkgil(n,y,z,E,dE,nt)
            call ks2ksr(y,u,us,E,cse,sie)
            call ks2car(u,us,x,xd,w)
            call car2oe(x,xd,pek)

            w = y(10)
            T = y(1) - dotp4(u,us)/w

            R(1) = dotp4(u,u)

            do i = 2, 20
               R(i) = R(i-1) * R(1)
            end do

            V = vmn(xd)
            ZbyR(1) = x(3)/R(1)
            ObyR(1) = 1.d0/R(1)
            do i = 1,ngeo_deg
               ZbyR(i+1) = ZbyR(i)*ZbyR(1)
               ObyR(i+1) = ObyR(i)*ObyR(1)
            end do
            FracJd = T/86400.d0
            djulian1 = djulian + FracJd
            call solarnpv(djulian1, ts)
            call lunarpv(djulian1, tm)

c     General (n,m) tesseral: refresh at the same cadence as the
c     third-body ephemerides above, since GMST (unlike Cnm_tess/
c     Snm_tess's own C/S inputs) changes with time, not position.
            call gmst_deg(djulian1, theta_tess)
            call rotate_tess_coeffs(ntess_use, Cnm_tess, Snm_tess,
     &                              theta_tess, Cnm_rot, Snm_rot)

            call third_body_aux(nsun_deg,  x, ts, R(1), amuS,
     &                          Rs(1), RbyRs, first_s, cphi_s, cphi_s2)
            call third_body_aux(nmoon_deg, x, tm, R(1), amuM,
     &                          Rm(1), RbyRm, first_m, cphi_m, cphi_m2)

            if (I_rkg .le. 1) goto 60

            E_0 = E_0+dE_0
            iE = dint(E_0*r2d)

!     Buffer state at every completed integration step
            idump = idump + 1
            traj_jd(idump) = djulian + T/86400.d0
            do i = 1,3
               traj_x(i,idump)  = x(i)
               traj_xd(i,idump) = xd(i)
            end do

c     -----------------------------------------------------------
c     Re-entry / divergence check: stop the run early, truncate
c     the trajectory buffer at the most useful final point, and
c     report why the propagation ended.
c
c     - Non-finite state (NaN): the integrator has diverged --
c       typically an unresolvable, near-singular low-perigee pass
c       during steep orbital decay (fixed per-revolution angular
c       stepping cannot resolve it). The state just buffered
c       above is garbage, so it is discarded (idump = idump - 1)
c       and the last valid epoch is reported instead.
c     - Altitude < 80 km: genuine re-entry. The state just
c       buffered above is finite and becomes the final OEM entry.
c
c     Either branch jumps (goto 90) past the remaining steps and
c     revolutions straight to the OEM write/cleanup section.
c     -----------------------------------------------------------
            h_alt = R(1) - R_Earth
            if (h_alt .ne. h_alt) then
               idump = idump - 1
               call jd2epoch(traj_jd(idump), epochstr)
               write(*,'(A,A,A)')
     &            'Integration diverged (non-finite state); ',
     &            'last valid epoch was ', epochstr
               write(*,'(A)') 'Propagation stopped; last valid '//
     &            'state saved as the last OEM data point.'
               goto 90
            else if (h_alt .lt. 80.d0) then
               call jd2epoch(traj_jd(idump), epochstr)
               write(*,'(A,F10.3,A,A)')
     &            'Re-entry occurred: altitude = ', h_alt,
     &            ' km (< 80 km) at epoch ', epochstr
               write(*,'(A)') 'Propagation stopped; final state '//
     &            'saved as the last OEM data point.'
               goto 90
            end if

         end do

         write(4,101) iE,y
      end do

  90  continue

c     --- Write CCSDS OEM v2.0 after run completes ---
c     CREATION_DATE and filename use the run-start UTC timestamp.
c     START_TIME / STOP_TIME use the actual propagated orbital epochs.
      oem_fname = 'output/KSROP_'//compact_epoch//'.oem'
      open (unit = 3, file = oem_fname)

      call jd2epoch(traj_jd(1),     oem_start)
      call jd2epoch(traj_jd(idump), epoch_stop)

      call write_oem(3, epochstr, oem_start, epoch_stop,
     &               traj_jd, traj_x, traj_xd, idump)
      close(3)
      write(*,'(A,A)') 'OEM written: ', oem_fname
      deallocate(traj_jd, traj_x, traj_xd)
      write(*,*) x

c     format specifiers
 101  format(10f20.9)           !KS-regular elements (regular.out)
 110  format(8f10.3)            !ATM.DAT scale height
 111  format(8e10.3)            !ATM.DAT density

 99   do i = 1,6
         close(i)
      end do
      stop
      end
