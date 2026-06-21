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
c     - Solar radiation pressure (cannon-ball model, ported from KSJLSDNP2.F)
c     -------------------------------------------------------------
c     Author: Harishkumar Sellamuthu
c     hari251086@gmail.com
c     -------------------------------------------------------------
c     Input file parameteres <input.dat> by line
c     --------
c     Initial position in geocentric inertial reference frame in km
c     Initial velocity in geocentric inertial reference frame in km/s
c     Initial date in yyyy,mm,dd,hh,mm,ss (to be included)
c     Number of revolutions, number of steps, tolerance (10^-15 is default)
c     Force model flags: n_force(1)=geo, n_force(2)=sun, n_force(3)=moon
c     Drag params: BN(kg/m2), IDRAG(0/1), WE_rot(rad/s), EPS_f, FR_rot
c     SRP params: CR(reflectivity), AM(m2/kg), IPSR(0/1), ISHAD(0=none/1=cyl/2=cone)
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
c     Drag: atmosphere tables and KS drag perturbation vectors
      dimension ALT_atm(ndim_atm),DEN_atm(ndim_atm),SCH_atm(ndim_atm)
      dimension qdrag(4),p_drag(3)
c     IDRAG: integer flag (0=off, 1=on); default integer via implicit i-n
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
      double precision CSEF_dg,SIEF_dg,FA_dg
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

      pi=4.d0*datan(1.d0)
      d2r=pi/180.d0
      r2d=180.d0/pi

c     Get UTC timestamp for output filenames
      call utc_now_epoch(epochstr, compact_epoch)
      reg_fname = 'output/KSROP_'//compact_epoch//'_Regular.out'

c     Input Files
      open (unit = 1, file = 'input/const_new.dat') !Constants
      open (unit = 2, file = 'input/input.dat')      !Simulation parameters
      open (unit = 6, file = 'input/ATM.DAT')        !Atmosphere density
      open (unit = 8, file = 'input/input.opm')      !CCSDS OPM initial state

c     Output files (unit 3 opened later with timestamped filename)
      open (unit = 4, file = reg_fname)            !KS elements (debug)
      open (unit = 5, file = 'output/ksrop.opm')   !CCSDS OPM initial Keplerian

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
c     SRP on/off, shadow model (0=always sunlit, 1=cylindrical, 2=conical)
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

      call car2oe(x,xd,pek)
      call oe2car(pek,x,xd,tole)
      call cal2jd(cal,djulian)
      call force_models(n_force,ngeo_deg,amuS,amuM)
      if (IPSR .eq. 0) PSR_srp = 0.d0
      SRPC_srp = CR_srp * AM_srp

      print *, x, xd

      write(*,*) '=============================================='

      if(ngeo_deg.eq.0.or.ngeo_deg.eq.1) then
         do i = 1,cn0
            c_j(i) = 0.d0
         end do
      else
!     compute j_n values
         call geo_coeff(ngeo_deg,c_j)
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

      call aLegP(ngeo_deg,ZbyR(1),aleg)
      call aLegP(nsun_deg,cphi_s,sleg)
      call aLegP(nmoon_deg,cphi_m,oleg)

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
         EA_dg    = pek(7)*d2r        ! eccentric anomaly: deg -> rad

c     Eccentric-anomaly step that advances EA_dg to the next
c     perigee/apogee crossing over the revolution's istep steps
         VIPP_dg = 4.d0
         if (EA_dg .lt. pi) VIPP_dg = 2.d0
         DE_dg = (VIPP_dg*pi - EA_dg)/istep

         SI_dg     = dsin(AI_dg)
         CSI_dg    = dcos(AI_dg)
         SIOMGA_dg = dsin(SOMGA_dg)

c     Perigee altitude above the oblate Earth (drag is active for
c     the revolution only when this is below 500 km)
         RREQ_dg = R_Earth*1.d3*(1.d0 - SI_dg*SI_dg*SIOMGA_dg*SIOMGA_dg
     &             *EPS_f)
         SQ_dg  = A_dg*(1.d0 - EC_dg)
         HSQ_dg = SQ_dg - RREQ_dg

         if (HSQ_dg .le. 500.d3) then
c     Reference density/scale-height at perigee, the rotating-
c     atmosphere correction factor F_dg, and the oblateness/
c     eccentricity exponent terms C_dg, XXX_dg used each step
            call INTPOL(ALT_atm, DEN_atm, ndim_atm, HSQ_dg, RHOPI_dg)
            call INTPOL(ALT_atm, SCH_atm, ndim_atm, HSQ_dg, H_dg)
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
            if (IDRAG .eq. 1 .and. HSQ_dg .le. 500.d3) then
               SEA_dg  = dsin(EA_dg)
               CEA_dg  = dcos(EA_dg)
               TER_dg  = 1.d0 - EC_dg*CEA_dg
               CSEF_dg = (CEA_dg - EC_dg)/TER_dg
               SIEF_dg = SEA_dg*dsqrt(1.d0 - EC_dg*EC_dg)/TER_dg
               FA_dg   = atan3(SIEF_dg,CSEF_dg)
               TER2_dg = dcos(SOMGA_dg)**2
     &                 - dcos(SOMGA_dg+FA_dg)**2
               TER_dg  = dexp(-2.d0*C_dg*TER2_dg
     &                  - XXX_dg*(1.d0-CEA_dg))
               DENS = -F_dg*RHOPI_dg*TER_dg*V/(BN*2.d0)*1.d-7
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

            Tau_geo = 0.d0
            do i = 2,ngeo_deg
               Tau_geo = Tau_geo-(i-1)*c_j(i)*Re(i)*aleg(i+1)*ObyR(i)
            end do

            Tau_3body1 = 0.d0
            do i = 2,nsun_deg
               Tau_3body1=Tau_3body1+(i+2)*sleg(i+1)*RbyRs(i)*RbyRs(1)*amuS
            end do
            Tau_3body2 = 0.d0
            do i = 2,nmoon_deg
               Tau_3body2=Tau_3body2+(i+2)*oleg(i+1)*RbyRm(i)*RbyRm(1)*amuM
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
c     Non-conservative adds: R*(u.(qdrag+qsrp))/(8w^3) - z(10)*(2/w^2)*(u.us)
c     -----------------------------------------------------------
            z(1) = (amue + Tau_geo + Tau_3body1 + Tau_3body2
     &             + R(1)*(dotp4(u,qdrag)+dotp4(u,qsrp)))/(8.d0*w*w*w)
     &             - z(10)*(2.d0/(w*w))*dotp4(u,us)

            if (ngeo_deg  .gt. 1) call aLegP(ngeo_deg, ZbyR(1), aleg)
            if (nsun_deg  .gt. 1) call aLegP(nsun_deg, cphi_s,  sleg)
            if (nmoon_deg .gt. 1) call aLegP(nmoon_deg,cphi_m,  oleg)

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

!     third-body terms
               call u2qu(u,qu,j)

               den_s = cphi_s2 - 1.d0
               den_m = cphi_m2 - 1.d0

               sq = ts(1)*qu(1)+ts(2)*qu(2)+ts(3)*qu(3)

               do i = 2,nsun_deg
                  slambda(i)=2.d0*R(i-1)*(i+1)*(sleg(i+2)
     &                 -cphi_s*sleg(i+1))/den_s

                  qs = (RbyRs(i)*(sq*slambda(i)/Rs(1)+
     &                 u(j)*(1-i)*sleg(i+1)-slambda(i)*cphi_s))
                  qsun(j) = qs + qsun(j)
               end do

               oq = tm(1)*qu(1)+tm(2)*qu(2)+tm(3)*qu(3)

               do i = 2,nmoon_deg
                  olambda(i)=2.d0*R(i-1)*(i+1)*(oleg(i+2)
     &                 -cphi_m*oleg(i+1))/den_m
                  qm = (RbyRm(i)*(oq*olambda(i)/Rm(1)+
     &                 u(j)*(1-i)*oleg(i+1)-olambda(i)*cphi_m))
                  qmoon(j) = qm + qmoon(j)
               end do

c     -----------------------------------------------------------
c     State elements: conservative + drag + SRP forces, energy correction
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

            call third_body_aux(nsun_deg,  x, ts, R(1), amuS,
     &                          Rs(1), RbyRs, first_s, cphi_s, cphi_s2)
            call third_body_aux(nmoon_deg, x, tm, R(1), amuM,
     &                          Rm(1), RbyRm, first_m, cphi_m, cphi_m2)

            if (I_rkg .le. 1) goto 60
c     Refresh the drag-model eccentricity from the latest osculating
c     elements (pek, just recomputed above for the converged state)
c     and advance the eccentric anomaly by the fixed per-revolution
c     step -- mirrors KSJLSDNP2.F's "EAR=EA; ...; EA=EAR+DE"
            EC_dg = pek(2)
            EA_dg = EA_dg + DE_dg

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
