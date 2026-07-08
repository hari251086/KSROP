C KSBENCH.F -- Cartesian ECI benchmark propagator for KSROP validation
C
C Integrates EOM in Cartesian ECI using DOP853 (Dormand-Prince 8th-order
C adaptive RK, Hairer).  Force model:
C   - Zonal geopotential J2..Jn  (EGM2008, same geo_coeff as KSROP)
C   - Atmospheric drag            (ATM.DAT direct interpolation)
C   - Solar radiation pressure    (cannonball, same shadow models)
C   - Sun + Moon third-body       (same solarnpv / lunarpv ephemerides)
C
C Input/output files: identical paths and formats as driver_KS.F
C   input/const_new.dat, input/input.dat, input/input.opm
C   input/ATM.DAT, input/EGM2008_to2190_TideFree
C   output/ksbench.opm, output/ksbench.oem
C
C Build:  ifx KSBENCH.F Subrouts.F dopri8.f /exe:ksbench.exe
C
C Author: Harishkumar Sellamuthu
C hari251086@gmail.com
C KSROP issue #17
C -------------------------------------------------------------

C ==============================================================
C COMMON blocks shared between main, EOM_CART, and BENCH_SOLOUT
C ==============================================================
C /XY/    -- already used by Subrouts.F (pi,d2r,r2d,amue,AU,R_Earth)
C /BPARAM/ -- force-model and epoch parameters
C /BGEO/  -- geopotential coefficients
C /BATM/  -- atmosphere table
C /BOUT/  -- output unit and epoch for OEM writing

C ==============================================================
      program KSBENCH
C ==============================================================
      implicit double precision (a-h, o-z)
      parameter (cn0 = 2190, ndim_atm = 291)

C     Shared with Subrouts.F
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth

C     Force model parameters
      common /BPARAM/ amuS,amuM,PSR_srp,SRPC_srp,BN,WE_rot,
     &                EPS_f,FR_rot,djulian0,
     &                ngeo_deg,nsun_deg,nmoon_deg,
     &                IDRAG,IPSR,ISHAD
      integer IDRAG,IPSR,ISHAD,ngeo_deg,nsun_deg,nmoon_deg

C     Geopotential zonal coefficients  c_j(n) = Jn
      common /BGEO/ c_j(cn0)

C     Atmosphere table (ALT in m, SCH in m, DEN in internal units)
      common /BATM/ ALT_atm(ndim_atm),DEN_atm(ndim_atm),
     &              SCH_atm(ndim_atm)

C     OEM output control
      common /BOUT/ oem_unit, djulian_b0
      integer oem_unit

      dimension y(6), pek(7), cal(6)
      character*23 epochstr, epoch_stop, oem_start
      character*15 compact_epoch
      character*40 oem_fname
      double precision n_force(3)
      integer i, idid
      double precision TOL, tspan
      external EOM_CART, BENCH_SOLOUT

C     ----------------------------------------------------------
C     Constants and setup
C     ----------------------------------------------------------
      pi    = 4.d0*datan(1.d0)
      d2r   = pi/180.d0
      r2d   = 180.d0/pi

C     ----------------------------------------------------------
C     Read inputs (same files as driver_KS)
C     ----------------------------------------------------------
      open(unit=1, file='input/const_new.dat')
      open(unit=2, file='input/input.dat')
      open(unit=6, file='input/ATM.DAT')
      open(unit=8, file='input/input.opm')
      open(unit=5, file='output/ksbench.opm')

      read(1,*) amue,R_Earth,AU,amuS,amuM
      read(1,*) ngeo_deg,nsun_deg,nmoon_deg
      read(1,*) PSR_srp

      call read_opm(8, y(1), y(4), cal)
      close(8)

      read(2,*) nrev, istep, tole
      read(2,*) n_force(1), n_force(2), n_force(3)
      read(2,*) BN, IDRAG, WE_rot, EPS_f, FR_rot
      read(2,*) CR_srp, AM_srp, IPSR, ISHAD
      close(2)

      if (n_force(1) .eq. 0.d0) ngeo_deg = 0
      if (n_force(2) .eq. 0.d0) amuS     = 0.d0
      if (n_force(3) .eq. 0.d0) amuM     = 0.d0
      if (IPSR .eq. 0)           PSR_srp  = 0.d0
      SRPC_srp = CR_srp * AM_srp

C     ----------------------------------------------------------
C     Atmosphere table (same scaling as driver_KS)
C     ----------------------------------------------------------
      ALT_atm(1) = 60.d0
      do i = 2,141
        ALT_atm(i) = ALT_atm(i-1) + 1.d0
      end do
      do i = 142,ndim_atm
        ALT_atm(i) = ALT_atm(i-1) + 2.d0
      end do

      read(6,100) (SCH_atm(i), i=1,61)
      read(6,100) (SCH_atm(i), i=62,ndim_atm)
      read(6,101) (DEN_atm(i), i=1,61)
      read(6,101) (DEN_atm(i), i=62,ndim_atm)
      close(6)

      do i = 1,ndim_atm
        ALT_atm(i) = ALT_atm(i) * 1000.d0    ! km -> m
        SCH_atm(i) = SCH_atm(i) * 1000.d0    ! km -> m
        DEN_atm(i) = DEN_atm(i) * 1.0d10     ! file -> internal
      end do

C     ----------------------------------------------------------
C     Geopotential coefficients
C     ----------------------------------------------------------
      do i = 1, cn0
        c_j(i) = 0.d0
      end do
      if (ngeo_deg .ge. 2) call geo_coeff(ngeo_deg, c_j)

C     ----------------------------------------------------------
C     Epoch: convert cal -> Julian date
C     ----------------------------------------------------------
      call cal2jd(cal, djulian0)
      djulian_b0 = djulian0
      call jd2epoch(djulian0, epochstr)
      call utc_now_epoch(epochstr, compact_epoch)

C     Initial OPM
      call car2oe(y(1), y(4), pek)
      call write_opm(5, epochstr, y(1), y(4), pek)
      close(5)

C     ----------------------------------------------------------
C     Integration time span: nrev revolutions
C     Period T = 2*pi*sqrt(a^3/mu);  a from osculating elements
C     ----------------------------------------------------------
      a_sma = pek(1)
      tspan = dble(nrev) * 2.d0*pi*dsqrt((a_sma**3)/amue)

      write(*,'(A,F12.4,A)') ' Integration span: ', tspan, ' seconds'
      write(*,'(A,I6,A)')    ' Revolutions: ', nrev, ''

C     ----------------------------------------------------------
C     OEM output file: header + DATA_START written now;
C     data lines written by BENCH_SOLOUT; DATA_STOP written after.
C     ----------------------------------------------------------
      oem_fname = 'output/ksbench.oem'
      open(unit=3, file=oem_fname)
      oem_unit = 3

      oem_start = epochstr
      call jd2epoch(djulian0 + tspan/86400.d0, epoch_stop)

      write(3,'(A)')       'CCSDS_OEM_VERS = 2.0'
      write(3,'(A,A)')     'CREATION_DATE  = ', epochstr
      write(3,'(A)')       'ORIGINATOR     = KSBENCH'
      write(3,'(A)')       ''
      write(3,'(A)')       'META_START'
      write(3,'(A)')       'OBJECT_NAME    = SATELLITE'
      write(3,'(A)')       'OBJECT_ID      = UNKNOWN'
      write(3,'(A)')       'CENTER_NAME    = EARTH'
      write(3,'(A)')       'REF_FRAME      = EME2000'
      write(3,'(A)')       'TIME_SYSTEM    = UTC'
      write(3,'(A,A)')     'START_TIME     = ', oem_start
      write(3,'(A,A)')     'STOP_TIME      = ', epoch_stop
      write(3,'(A)')       'META_STOP'
      write(3,'(A)')       ''
      write(3,'(A)')       'DATA_START'

C     Write the initial state as the first OEM data point
      call jd2epoch(djulian0, epochstr)
      write(3,102) epochstr, y(1),y(2),y(3),y(4),y(5),y(6)

C     ----------------------------------------------------------
C     Integrate with DOP853
C     tolerance matched to driver_KS's tole (typically 1e-15)
C     ----------------------------------------------------------
      TOL  = tole
      tstart = 0.d0
      call DOP853(6, EOM_CART, tstart, y, tspan,
     &            TOL, BENCH_SOLOUT, 1, idid)

      write(3,'(A)') 'DATA_STOP'
      close(3)

      if (idid .eq. 1) then
        write(*,'(A)') 'KSBENCH: integration successful'
      else
        write(*,'(A,I4)') 'KSBENCH: DOP853 returned idid=', idid
      end if

      write(*,'(A,A)') 'OEM written: ', oem_fname
      write(*,'(A)') 'Final state (x,y,z km, vx,vy,vz km/s):'
      write(*,'(6F16.6)') y(1),y(2),y(3),y(4),y(5),y(6)

  100 format(8f10.3)
  101 format(8e10.3)
  102 format(A23,6F20.6)

      do i = 1,9
        close(i)
      end do
      stop
      end

C ==============================================================
C     DOP853 output callback: writes one OEM data line per step
C ==============================================================
      subroutine BENCH_SOLOUT(nr, xold, t, y, n, irtrn)
      implicit double precision (a-h, o-z)
      dimension y(n)
      integer nr, n, irtrn

      common /BOUT/ oem_unit, djulian_b0
      integer oem_unit

      character*23 epochstr

      djulian_c = djulian_b0 + t/86400.d0
      call jd2epoch(djulian_c, epochstr)
      write(oem_unit,100) epochstr, y(1),y(2),y(3),y(4),y(5),y(6)

  100 format(A23,6F20.6)
      return
      end

C ==============================================================
C     EOM_CART: Equations of motion in Cartesian ECI
C     y(1:3) = position (km), y(4:6) = velocity (km/s)
C     f(1:3) = velocity,      f(4:6) = total acceleration (km/s2)
C ==============================================================
      subroutine EOM_CART(neq, t, y, f)
      implicit double precision (a-h, o-z)
      integer neq
      dimension y(neq), f(neq)

      parameter (cn0 = 2190, ndim_atm = 291)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
      common /BPARAM/ amuS,amuM,PSR_srp,SRPC_srp,BN,WE_rot,
     &                EPS_f,FR_rot,djulian0,
     &                ngeo_deg,nsun_deg,nmoon_deg,
     &                IDRAG,IPSR,ISHAD
      integer IDRAG,IPSR,ISHAD,ngeo_deg,nsun_deg,nmoon_deg
      common /BGEO/ c_j(cn0)
      common /BATM/ ALT_atm(ndim_atm),DEN_atm(ndim_atm),
     &              SCH_atm(ndim_atm)

      double precision PLeg(0:cn0+1), PLegd(0:cn0+1)
      dimension ts(3), tm(3), r_vec(3), vrel(3)
      dimension r_tb(3)

      x1 = y(1); x2 = y(2); x3 = y(3)
      v1 = y(4); v2 = y(5); v3 = y(6)

      r2 = x1*x1 + x2*x2 + x3*x3
      r  = dsqrt(r2)
      r3 = r2*r

C     Two-body (km/s2)
      ax = -amue*x1/r3
      ay = -amue*x2/r3
      az = -amue*x3/r3

C     ----------------------------------------------------------
C     Zonal geopotential J2..ngeo_deg
C     PLeg(n) = Pn(t), PLegd(n) = Pn'(t), t = x3/r = sin(lat)
C     Recurrences:
C       Pn = ((2n-1)*t*Pn-1 - (n-1)*Pn-2) / n
C       Pn'= (2n-1)*Pn-1 + Pn-2'
C     Acceleration from n-th zonal term (Jn = c_j(n)):
C       ax_n = mu*Jn*(Re/r)^n/r^2 * (x1/r) * [(n+1)*Pn + t*Pnd]
C       ay_n = same * (x2/r)
C       az_n = mu*Jn*(Re/r)^n/r^2 * [(n+1)*t*Pn - (1-t2)*r/r * Pnd]
C     ----------------------------------------------------------
      if (ngeo_deg .ge. 2) then
        tlat = x3/r
        t2   = tlat*tlat

        PLeg(0)  = 1.d0
        PLeg(1)  = tlat
        PLegd(0) = 0.d0
        PLegd(1) = 1.d0

        do i = 2, ngeo_deg
          dn = dble(i)
          PLeg(i) = ((2.d0*dn-1.d0)*tlat*PLeg(i-1)
     &               - (dn-1.d0)*PLeg(i-2)) / dn
          PLegd(i) = (2.d0*dn-1.d0)*PLeg(i-1) + PLegd(i-2)
        end do

        ReN = R_Earth
        do i = 2, ngeo_deg
          dn  = dble(i)
          ReN = ReN * (R_Earth/r)
          Jn  = c_j(i)
          fac = amue*Jn*ReN/r2

          Pn  = PLeg(i)
          Pnd = PLegd(i)

C         Coefficient for x,y components
          Cxy = fac * ((dn+1.d0)*Pn + tlat*Pnd)

C         z component (non-singular form)
          Cz  = fac * ((dn+1.d0)*tlat*Pn - (1.d0-t2)*Pnd)

          ax = ax + Cxy*(x1/r)
          ay = ay + Cxy*(x2/r)
          az = az + Cz
        end do
      end if

C     ----------------------------------------------------------
C     Atmospheric drag (ATM.DAT direct interpolation at altitude)
C     Relative velocity accounts for co-rotating atmosphere.
C     Same 1e-7 scaling as driver_KS (DEN internal units).
C     ----------------------------------------------------------
      if (IDRAG .eq. 1) then
        h_m = (r - R_Earth) * 1000.d0     ! altitude in m
        if (h_m .ge. ALT_atm(1) .and.
     &      h_m .le. ALT_atm(ndim_atm)) then
          call INTPOL(ALT_atm, DEN_atm, ndim_atm, h_m, rho_loc)

C         Relative velocity: satellite minus rotating atmosphere
          vrel(1) = v1 + WE_rot*FR_rot*x2    ! +WE*FR*y (omega x r)_x
          vrel(2) = v2 - WE_rot*FR_rot*x1    ! -WE*FR*x
          vrel(3) = v3

          Vrel = dsqrt(vrel(1)**2+vrel(2)**2+vrel(3)**2)

          if (Vrel .gt. 0.d0) then
            DENS = -rho_loc * Vrel / (BN*2.d0) * 1.d-7
            ax = ax + DENS*vrel(1)
            ay = ay + DENS*vrel(2)
            az = az + DENS*vrel(3)
          end if
        end if
      end if

C     ----------------------------------------------------------
C     Current Julian date for ephemerides
C     ----------------------------------------------------------
      djulian_c = djulian0 + t/86400.d0

C     ----------------------------------------------------------
C     Sun third-body + SRP (one call to solarnpv for both)
C     ----------------------------------------------------------
      if (amuS .gt. 0.d0 .or. PSR_srp .gt. 0.d0) then
        call solarnpv(djulian_c, ts)

        if (amuS .gt. 0.d0) then
          rs = dsqrt(ts(1)**2+ts(2)**2+ts(3)**2)
          rs3 = rs*rs*rs

          drx = ts(1)-x1; dry = ts(2)-x2; drz = ts(3)-x3
          dr  = dsqrt(drx*drx + dry*dry + drz*drz)
          dr3 = dr*dr*dr

          ax = ax + amuS*(drx/dr3 - ts(1)/rs3)
          ay = ay + amuS*(dry/dr3 - ts(2)/rs3)
          az = az + amuS*(drz/dr3 - ts(3)/rs3)
        end if

        if (PSR_srp .gt. 0.d0) then
          if (ISHAD .eq. 1) then
            ROHP = shadfncyl(x1,x2,x3,ts(1),ts(2),ts(3))
          else if (ISHAD .eq. 2) then
            ROHP = shadfncone(x1,x2,x3,ts(1),ts(2),ts(3))
          else
            ROHP = 1.d0
          end if
          sr1 = ts(1)-x1; sr2 = ts(2)-x2; sr3 = ts(3)-x3
          xsrr = dsqrt(sr1*sr1+sr2*sr2+sr3*sr3)
          fac  = -ROHP*PSR_srp*SRPC_srp/xsrr * 1.d-3
          ax = ax + fac*sr1
          ay = ay + fac*sr2
          az = az + fac*sr3
        end if
      end if

C     ----------------------------------------------------------
C     Moon third-body
C     ----------------------------------------------------------
      if (amuM .gt. 0.d0) then
        call lunarpv(djulian_c, tm)

        rm = dsqrt(tm(1)**2+tm(2)**2+tm(3)**2)
        rm3 = rm*rm*rm

        drx = tm(1)-x1; dry = tm(2)-x2; drz = tm(3)-x3
        dr  = dsqrt(drx*drx + dry*dry + drz*drz)
        dr3 = dr*dr*dr

        ax = ax + amuM*(drx/dr3 - tm(1)/rm3)
        ay = ay + amuM*(dry/dr3 - tm(2)/rm3)
        az = az + amuM*(drz/dr3 - tm(3)/rm3)
      end if

C     Derivatives
      f(1) = v1; f(2) = v2; f(3) = v3
      f(4) = ax; f(5) = ay; f(6) = az

      return
      end
