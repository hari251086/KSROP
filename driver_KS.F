c     COPYRIGHT 2018, HARISHKUMAR SELLAMUTHU, ALL RIGHTS RESERVED
C     -------------------------------------------------------------
c     KS - Numerical driver
c       Description:
c     Orbit propagation using Kustaanheimo-Stiefel regular elements
c     with numerical integrator RKG4. (Earth)
c     Perturbations included are:
c     - Earth's oblateness ( up to 2190th degree )
c     - Luni-solar perturbations
c     - Atmospheric drag ( Oblate exponential atmosphere,
c       tabulated density model via ATM.DAT )
c     To be added
c     - Solar radiation pressure (cannon ball)
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
c     OEM/OPM epoch strings (YYYY-MM-DDTHH:MM:SS.sss = 23 chars)
      character*23 epochstr, epoch_stop
c     Trajectory buffer: allocated once nrev is known; written after run
      double precision, allocatable :: traj_jd(:), traj_x(:,:), traj_xd(:,:)

      pi=4.d0*datan(1.d0)
      d2r=pi/180.d0
      r2d=180.d0/pi

c     Input Files
      open (unit = 1, file = 'const_new.dat') !Constants
      open (unit = 2, file = 'input.dat')      !Simulation parameters
      open (unit = 6, file = 'ATM.DAT')        !Atmosphere density table
      open (unit = 8, file = 'input.opm')      !CCSDS OPM - Initial state

c     Output files
      open (unit = 3, file = 'ksrop.oem')   !CCSDS OEM - State trajectory
      open (unit = 4, file = 'regular.out') !Internal  - KS elements
      open (unit = 5, file = 'ksrop.opm')   !CCSDS OPM - Initial Keplerian

      read(1,*) amue,R_Earth,AU,amuS,amuM
      read(1,*) ngeo_deg,nsun_deg,nmoon_deg

c     Read initial state (position, velocity, epoch) from OPM file
      call read_opm(8, x, xd, cal)
      close(8)

c     Read simulation parameters from input.dat
      read(2,*,end=99) nrev,istep,tole !Revolutions, steps, tolerance
      read(2,*) n_force(1), n_force(2), n_force(3)
c     Drag parameters: ballistic coeff (kg/m2), drag on/off,
c     Earth rotation rate (rad/s), Earth flattening, rotation factor
      read(2,*) BN, IDRAG, WE_rot, EPS_f, FR_rot

      allocate(traj_jd(nrev+1), traj_x(3,nrev+1), traj_xd(3,nrev+1))
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
      Rs(1) = vmn(ts)
      Rm(1) = vmn(tm)
      RbyRs(1) = R(1)/Rs(1)
      RbyRm(1) = R(1)/Rm(1)
      first_s = amuS/Rs(1)
      first_m = amuM/Rm(1)
      do i = 1, nsun_deg
         RbyRs(i+1) = RbyRs(i)*RbyRs(1)
      end do
      do i = 1, nmoon_deg
         RbyRm(i+1) = RbyRm(i)*RbyRm(1)
      end do
      ph_s = 0.d0
      ph_m = 0.d0
      do i = 1,3
         ph_s = ph_s + x(i) * ts(i)
         ph_m = ph_m + x(i) * tm(i)
      end do
      cphi_s = ph_s/(R(1) * Rs(1))
      cphi_m = ph_m/(R(1) * Rm(1))
      cphi_s2 = cphi_s  * cphi_s
      cphi_m2 = cphi_m  * cphi_m

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
         cn_s(i) = -first_s * RbyRs(i)
         V_p(i) = cn_s(i) * sleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

      do i = 2,nmoon_deg
         cn_m(i) = -first_m * RbyRm(i)
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

      Br = y(5)*y(6)-y(4)*y(7)+y(3)*y(8)-y(2)*y(9)

      sie = dsin(E/2.d0)
      cse = dcos(E/2.d0)
      iE = dint(E)

      call car2oe(x,xd,pek)
      call jd2epoch(djulian, epochstr)

c     --- CCSDS OPM v2.0 (initial epoch) — written immediately ---
      write(5,'(A)')     'CCSDS_OPM_VERS = 2.0'
      write(5,'(A,A)')   'CREATION_DATE  = ', epochstr
      write(5,'(A)')     'ORIGINATOR     = KSROP'
      write(5,'(A)')     ''
      write(5,'(A)')     'META_START'
      write(5,'(A)')     'OBJECT_NAME    = SATELLITE'
      write(5,'(A)')     'OBJECT_ID      = UNKNOWN'
      write(5,'(A)')     'CENTER_NAME    = EARTH'
      write(5,'(A)')     'REF_FRAME      = EME2000'
      write(5,'(A)')     'TIME_SYSTEM    = UTC'
      write(5,'(A)')     'META_STOP'
      write(5,'(A)')     ''
      write(5,'(A)')     'STATE_VECTOR'
      write(5,'(A,A)')   'EPOCH          = ', epochstr
      write(5,'(A,F16.6,A)')  'X              = ', x(1),  ' [km]'
      write(5,'(A,F16.6,A)')  'Y              = ', x(2),  ' [km]'
      write(5,'(A,F16.6,A)')  'Z              = ', x(3),  ' [km]'
      write(5,'(A,F16.9,A)')  'X_DOT          = ', xd(1), ' [km/s]'
      write(5,'(A,F16.9,A)')  'Y_DOT          = ', xd(2), ' [km/s]'
      write(5,'(A,F16.9,A)')  'Z_DOT          = ', xd(3), ' [km/s]'
      write(5,'(A)')     ''
      write(5,'(A)')     'KEPLERIAN_ELEMENTS'
      write(5,'(A,F16.6,A)')  'SEMI_MAJOR_AXIS   = ', pek(1),' [km]'
      write(5,'(A,F16.9)')    'ECCENTRICITY      = ', pek(2)
      write(5,'(A,F16.6,A)')  'INCLINATION       = ', pek(3),' [deg]'
      write(5,'(A,F16.6,A)')  'RA_OF_ASC_NODE    = ', pek(4),' [deg]'
      write(5,'(A,F16.6,A)')  'ARG_OF_PERICENTER = ', pek(5),' [deg]'
      write(5,'(A,F16.6,A)')  'MEAN_ANOMALY      = ', pek(6),' [deg]'

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
         do ki = 1,istep        !Steps
 60         continue
            call u2uu(u,uu)

c     -----------------------------------------------------------
c     Atmospheric drag perturbation
c     Compute drag acceleration p_drag(3) in km/s^2, then
c     transform to KS space via the L(u) matrix: qdrag(4) = L(u)*p_drag
c     -----------------------------------------------------------
            h_drag = (R(1) - R_Earth) * 1000.d0  ! altitude in m
            if (IDRAG .eq. 1 .and. h_drag .lt. 500.d3) then
               call INTPOL(ALT_atm, DEN_atm, ndim_atm, h_drag, rho_drag)
c     DENS = -(rho/(2*BN))*|v|*1e-7 gives drag scale factor in km/s^2
c     (unit factor matches KSJLSDNP2.F convention: rho[kg/m3],v[km/s])
               DENS = -rho_drag * V / (BN * 2.d0) * 1.d-7
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
c     Energy element: dissipative drag contribution
c     dw/ds = -0.5 * (us . qdrag) / w
c     -----------------------------------------------------------
            z(10) = -0.5d0*(us(1)*qdrag(1)+us(2)*qdrag(2)
     &                    + us(3)*qdrag(3)+us(4)*qdrag(4))/w

c     -----------------------------------------------------------
c     Time element: conservative + drag contributions
c     Tau_term replaced by Tau_geo (fixes uninitialized variable)
c     Drag adds: R*(u.qdrag)/(8w^3) - z(10)*(2/w^2)*(u.us)
c     -----------------------------------------------------------
            z(1) = (amue + Tau_geo + Tau_3body1 + Tau_3body2
     &             + R(1)*dotp4(u,qdrag))/(8.d0*w*w*w)
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
c     State elements: add drag force (-R*qdrag) and energy correction
c     z(j+1) = [geo + 3body - R*drag]*sin(E/2)*0.25/w^2
c              + 2*z(10)*us(j)*sin(E/2)/w
c     z(j+5) = -[geo + 3body - R*drag]*cos(E/2)*0.25/w^2
c              - 2*z(10)*us(j)*cos(E/2)/w
c     -----------------------------------------------------------
               aqt = sie * 0.25d0
               bqt = cse * 0.25d0

               z(j+1) = (q(j)*amue + first_s*qsun(j)
     &              + first_m*qmoon(j) - R(1)*qdrag(j))
     &              * aqt/(w*w)
     &              + 2.d0*z(10)*us(j)*sie/w

               z(j+5) = -(q(j)*amue + first_s*qsun(j)
     &              + first_m*qmoon(j) - R(1)*qdrag(j))
     &              * bqt/(w*w)
     &              - 2.d0*z(10)*us(j)*cse/w

            end do

            I_rkg = rkgil(n,y,z,E,dE*Gam,nt)
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
            do i = 2,ngeo_deg
               cn_j(i) = amue*Re(i)*c_j(i)*ObyR(i+1)
            end do
            if (ngeo_deg .gt. 1) then
               call aLegP(ngeo_deg,ZbyR(1),aleg)
               V_pot = 0.d0
               do i = 2,ngeo_deg
                  V_p(i) = cn_j(i)*aleg(i+1)
                  V_pot = V_pot+V_p(i)
               end do
            end if

            FracJd = T/86400.d0
            djulian1 = djulian + FracJd
            call solarnpv(djulian1, ts)
            call lunarpv(djulian1, tm)

            Rs(1) = vmn(ts)
            Rm(1) = vmn(tm)
            RbyRs(1) = R(1)/Rs(1)
            RbyRm(1) = R(1)/Rm(1)
            first_s = -amuS/Rs(1)
            first_m = -amuM/Rm(1)

            do i = 1, nsun_deg
               RbyRs(i+1) = RbyRs(i)*RbyRs(1)
            end do
            do i = 1, nmoon_deg
               RbyRm(i+1) = RbyRm(i)*RbyRm(1)
            end do

            ph_s = 0.d0
            ph_m = 0.d0
            do i = 1,3
               ph_s = ph_s + x(i) * ts(i)
               ph_m = ph_m + x(i) * tm(i)
            end do

            cphi_s = ph_s/(R(1) * Rs(1))
            cphi_m = ph_m/(R(1) * Rm(1))

            cphi_s2 = cphi_s  * cphi_s
            cphi_m2 = cphi_m  * cphi_m

            if (nsun_deg  .gt. 1) then
               call aLegP(nsun_deg,cphi_s,sleg)
               do i = 2,nsun_deg
                  cn_s(i) = first_s * RbyRs(i)
                  V_p(i)  = cn_s(i) * sleg(i+1)
                  V_pot   = V_pot + V_p(i)
               end do
            end if
            if (nmoon_deg .gt. 1) then
               call aLegP(nmoon_deg,cphi_m,oleg)
               do i = 2,nmoon_deg
                  cn_m(i) = first_m * RbyRm(i)
                  V_p(i)  = cn_m(i) * oleg(i+1)
                  V_pot   = V_pot + V_p(i)
               end do
            end if

            Br = y(5)*y(6)-y(4)*y(7)+y(3)*y(8)-y(2)*y(9)
            if (I_rkg .le. 1) goto 60
 70         continue
            E_0 = E_0+dE_0
            print *,x
            print *,Br
            iE = dint(E_0*r2d)

         end do

!     Buffer trajectory at end of each revolution
         idump = idump + 1
         traj_jd(idump) = djulian + T/86400.d0
         do i = 1,3
            traj_x(i,idump)  = x(i)
            traj_xd(i,idump) = xd(i)
         end do
         write(4,101) iE,y
      end do

c     --- Write CCSDS OEM v2.0 after run completes ---
c     STOP_TIME is now the exact final epoch, not a Keplerian estimate.
      call jd2epoch(traj_jd(1),     epochstr)
      call jd2epoch(traj_jd(idump), epoch_stop)
      write(3,'(A)')     'CCSDS_OEM_VERS = 2.0'
      write(3,'(A,A)')   'CREATION_DATE  = ', epochstr
      write(3,'(A)')     'ORIGINATOR     = KSROP'
      write(3,'(A)')     ''
      write(3,'(A)')     'META_START'
      write(3,'(A)')     'OBJECT_NAME    = SATELLITE'
      write(3,'(A)')     'OBJECT_ID      = UNKNOWN'
      write(3,'(A)')     'CENTER_NAME    = EARTH'
      write(3,'(A)')     'REF_FRAME      = EME2000'
      write(3,'(A)')     'TIME_SYSTEM    = UTC'
      write(3,'(A,A)')   'START_TIME     = ', epochstr
      write(3,'(A,A)')   'STOP_TIME      = ', epoch_stop
      write(3,'(A)')     'META_STOP'
      write(3,'(A)')     ''
      write(3,'(A)')     'DATA_START'
      do i = 1, idump
         call jd2epoch(traj_jd(i), epochstr)
         write(3,103) epochstr,
     &                traj_x(1,i),  traj_x(2,i),  traj_x(3,i),
     &                traj_xd(1,i), traj_xd(2,i), traj_xd(3,i)
      end do
      write(3,'(A)') 'DATA_STOP'
      deallocate(traj_jd, traj_x, traj_xd)
      write(*,*) x

c     format specifiers
 101  format(10f20.9)           !KS-regular elements (regular.out)
 103  format(A23,3(1x,F16.6),3(1x,F14.9))  !OEM data line
 110  format(8f10.3)            !ATM.DAT scale height
 111  format(8e10.3)            !ATM.DAT density

 99   do i = 1,6
         close(i)
      end do
      stop
      end
