c     ============================================================
c     swx.F — Space weather + 2-D atmosphere loaders
c
c     Ported from OREM's src/swx.F (issue #26 there) back into KSROP
c     itself: KSROP already had the ATM2D.DAT generator
c     (app/gen_atm2d_jr71.F) but nothing in driver_KS.F ever consumed
c     it at runtime -- only OREM's downstream fork did. This closes
c     that gap so KSROP's own propagation benefits from epoch-
c     resolved solar/geomagnetic activity, not just a static table.
c
c     Fills two common blocks consumed by driver_KS.F's per-
c     revolution drag setup via sw_tinf/atm2d_interp (also defined in
c     this file, below). Unlike OREM's own swx.F (where propagate_ks
c     is already a shared src/ library file, so the runtime consumers
c     live there instead), KSROP's driver_KS.F is a standalone
c     `program` in app/ -- putting the consumers there would make
c     them unreachable by test_sw's separate executable, since fpm
c     only auto-links src/ into every target, not other apps' own
c     sources. Keeping the whole loader+consumer set together in this
c     one src/ file sidesteps that. Neither loader is required: when
c     nsw=0 or n2d_t=0 (the compile-time default), driver_KS.F uses
c     the legacy single-table ALT/DEN/SCH arrays unchanged, so every
c     existing run is unaffected unless both input files are present
c     and get auto-loaded at startup.
c
c       call sw_load('input/SW-All.csv', ierr)      ! CelesTrak CSV
c       call atm2d_load('input/ATM2D.DAT', ierr)    ! gen_atm2d_jr71
c
c     /swdat/: per-row (JD, F10.7_ADJ, F10.7_ADJ_CENTER81, mean Kp)
c       from the CelesTrak file. Observed rows are daily; predicted
c       rows (PRM/PRD types) are monthly — lookup is by JD binary
c       search (latest row with jd <= t), never by index arithmetic.
c       Missing Kp fields (predicted rows) default to Kp=2.0.
c     /atm2dc/: rho(alt,Tinf) and H(alt,Tinf) tables, stored in the
c       SAME internal scaled units as the legacy path (alt in m,
c       sch in m, den x 1e10) so the drag formulas are untouched.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================
      subroutine sw_load(fname, ierr)
      implicit none
      character*(*) fname
      integer ierr

      integer NSWMAX
      parameter (NSWMAX = 31000)
      double precision sw_jd(NSWMAX), sw_f107(NSWMAX)
      double precision sw_f81(NSWMAX), sw_kp(NSWMAX)
      integer nsw
      common /swdat/ sw_jd, sw_f107, sw_f81, sw_kp, nsw

      character*400 line
      character*40 tok(31)
      double precision cal(6), jd, f107, f81, akp, ksum
      integer i, j, ntok, ip, lt, ios

      ierr = 0
      nsw = 0
      open(unit=91, file=fname, status='old', err=900)
      read(91,'(A)',iostat=ios) line          ! header
      if (ios .ne. 0) goto 900

 10   read(91,'(A)',end=50,err=50) line
      if (len_trim(line) .lt. 10) goto 10

c     --- split on commas into up to 31 tokens ---
      ntok = 1
      tok(1) = ' '
      ip = 0
      do i = 1, len_trim(line)
         if (line(i:i) .eq. ',') then
            if (ntok .ge. 31) goto 12
            ntok = ntok + 1
            tok(ntok) = ' '
            ip = 0
         else
            ip = ip + 1
            tok(ntok)(ip:ip) = line(i:i)
         end if
      end do
 12   continue
      if (ntok .lt. 30) goto 10

c     --- DATE yyyy-mm-dd -> JD (0h UTC) ---
      read(tok(1)(1:4),'(F4.0)',err=10) cal(1)
      read(tok(1)(6:7),'(F2.0)',err=10) cal(2)
      read(tok(1)(9:10),'(F2.0)',err=10) cal(3)
      cal(4) = 0.d0
      cal(5) = 0.d0
      cal(6) = 0.d0
      call cal2jd(cal, jd)

c     --- KP_SUM (token 12, x10 sum of 8) -> daily mean Kp ---
      akp = 2.0d0
      lt = len_trim(tok(12))
      if (lt .gt. 0) then
         read(tok(12)(1:lt),*,err=14) ksum
         akp = ksum / 80.d0
      end if
 14   continue

c     --- F10.7_ADJ (26) and F10.7_ADJ_CENTER81 (30) ---
      lt = len_trim(tok(26))
      if (lt .le. 0) goto 10
      read(tok(26)(1:lt),*,err=10) f107
      lt = len_trim(tok(30))
      if (lt .le. 0) goto 10
      read(tok(30)(1:lt),*,err=10) f81

      if (nsw .lt. NSWMAX) then
         nsw = nsw + 1
         sw_jd(nsw)   = jd
         sw_f107(nsw) = f107
         sw_f81(nsw)  = f81
         sw_kp(nsw)   = akp
      end if
      goto 10

 50   close(91)
      if (nsw .lt. 100) then
         nsw = 0
         ierr = 2
      end if
      return

 900  ierr = 1
      nsw = 0
      return
      end

c     ============================================================
c     atm2d_load — read ATM2D.DAT (gen_atm2d_jr71.F) and store in
c     the legacy path's internal scaled units.
c     ============================================================
      subroutine atm2d_load(fname, ierr)
      implicit none
      character*(*) fname
      integer ierr

      integer NA2D, NT2D
      parameter (NA2D = 291, NT2D = 39)
      double precision t2d(NT2D), alt2d(NA2D)
      double precision den2d(NA2D,NT2D), sch2d(NA2D,NT2D)
      integer n2d_alt, n2d_t
      common /atm2dc/ t2d, alt2d, den2d, sch2d, n2d_alt, n2d_t

      character*10 tag
      integer i, j, na, nt

      ierr = 0
      n2d_t = 0
      open(unit=92, file=fname, status='old', err=900)
      read(92,'(A10,2I6)',err=910) tag, na, nt
      if (tag .ne. 'ATM2D_JR71' .or. na .ne. NA2D .or.
     &    nt .ne. NT2D) goto 910
      read(92,*,err=910) (t2d(j), j=1,nt)
      do j = 1, nt
         read(92,*,err=910) (sch2d(i,j), i=1,na)
         read(92,*,err=910) (den2d(i,j), i=1,na)
      end do
      close(92)

      alt2d(1) = 60.d0
      do i = 2, 141
         alt2d(i) = alt2d(i-1) + 1.d0
      end do
      do i = 142, NA2D
         alt2d(i) = alt2d(i-1) + 2.d0
      end do

c     scale to the legacy internal units (m, m, x1e10)
      do i = 1, NA2D
         alt2d(i) = alt2d(i) * 1.d3
         do j = 1, nt
            sch2d(i,j) = sch2d(i,j) * 1.d3
            den2d(i,j) = den2d(i,j) * 1.d10
         end do
      end do

      n2d_alt = na
      n2d_t   = nt
      return

 900  ierr = 1
      return
 910  ierr = 2
      n2d_t = 0
      close(92)
      return
      end
c     ============================================================
c     sw_tinf — exospheric temperature at JD t from the loaded
c     space-weather history. Ported from OREM's propagate_ks.F
c     (issue #26 there). J71 nighttime-minimum plus the h<200 km
c     geomagnetic term (same convention as the ATM.DAT generator).
c     Returns tinf=0 when no history is loaded or t precedes the
c     record -- callers (driver_KS.F) then use the legacy path.
c     The common block defaults to nsw=0 regardless of whether
c     sw_load was ever called, so this is always safe to call.
c     ============================================================
      subroutine sw_tinf(t, tinf)
      implicit none
      double precision t, tinf

      integer NSWMAX
      parameter (NSWMAX = 31000)
      double precision sw_jd(NSWMAX), sw_f107(NSWMAX)
      double precision sw_f81(NSWMAX), sw_kp(NSWMAX)
      integer nsw
      common /swdat/ sw_jd, sw_f107, sw_f81, sw_kp, nsw

      double precision tc
      integer lo, hi, mid

      tinf = 0.d0
      if (nsw .le. 0) return
      if (t .lt. sw_jd(1)) return

c     binary search: largest i with sw_jd(i) <= t (rows are daily
c     in the observed era but monthly in the predicted tail, so
c     index arithmetic is never valid -- always search)
      lo = 1
      hi = nsw
      do while (hi - lo .gt. 1)
         mid = (lo + hi)/2
         if (sw_jd(mid) .le. t) then
            lo = mid
         else
            hi = mid
         end if
      end do
      if (sw_jd(hi) .le. t) lo = hi

      tc = 379.d0 + 3.24d0*sw_f81(lo)
     &   + 1.3d0*(sw_f107(lo) - sw_f81(lo))
      tinf = tc + 14.d0*sw_kp(lo) + 0.02d0*dexp(sw_kp(lo))

      return
      end

c     ============================================================
c     atm2d_interp — bilinear rho/H at (altitude h [m], tinf [K])
c     from the 2-D table loaded by atm2d_load (swx.F), in the same
c     internal scaled units as the legacy path (den x1e10, H in m).
c     Ported from OREM's propagate_ks.F (issue #26 there). Tinf
c     clamps to the grid edges; altitude outside the table returns
c     istat=1 (drag treated inactive, same contract as the legacy
c     range gate); istat=2 = table not loaded.
c     ============================================================
      subroutine atm2d_interp(h, tinf, rho_s, h_s, istat)
      implicit none
      double precision h, tinf, rho_s, h_s
      integer istat

      integer NA2D, NT2D
      parameter (NA2D = 291, NT2D = 39)
      double precision t2d(NT2D), alt2d(NA2D)
      double precision den2d(NA2D,NT2D), sch2d(NA2D,NT2D)
      integer n2d_alt, n2d_t
      common /atm2dc/ t2d, alt2d, den2d, sch2d, n2d_alt, n2d_t

      double precision tv, fa, ft, d1, d2, s1, s2
      integer ia, it, i

      istat = 0
      rho_s = 0.d0
      h_s   = 1.d0
      if (n2d_t .le. 0) then
         istat = 2
         return
      end if
      if (h .lt. alt2d(1) .or. h .gt. alt2d(n2d_alt)) then
         istat = 1
         return
      end if

      tv = tinf
      if (tv .lt. t2d(1))      tv = t2d(1)
      if (tv .gt. t2d(n2d_t))  tv = t2d(n2d_t)

      ia = 1
      do i = 2, n2d_alt
         if (alt2d(i) .ge. h) then
            ia = i - 1
            goto 20
         end if
      end do
      ia = n2d_alt - 1
 20   continue
      it = 1
      do i = 2, n2d_t
         if (t2d(i) .ge. tv) then
            it = i - 1
            goto 30
         end if
      end do
      it = n2d_t - 1
 30   continue

      fa = (h - alt2d(ia)) / (alt2d(ia+1) - alt2d(ia))
      ft = (tv - t2d(it)) / (t2d(it+1) - t2d(it))

      d1 = den2d(ia,it)   + fa*(den2d(ia+1,it)   - den2d(ia,it))
      d2 = den2d(ia,it+1) + fa*(den2d(ia+1,it+1) - den2d(ia,it+1))
      rho_s = d1 + ft*(d2 - d1)

      s1 = sch2d(ia,it)   + fa*(sch2d(ia+1,it)   - sch2d(ia,it))
      s2 = sch2d(ia,it+1) + fa*(sch2d(ia+1,it+1) - sch2d(ia,it+1))
      h_s = s1 + ft*(s2 - s1)
      if (h_s .le. 0.d0) h_s = 1.d0

      return
      end

c     ============================================================
c     jacchia_bulge_extrema — exospheric temperature extrema.
c     Ported from OREM's src/propagate_ks.F (issue #32 Phase 5
c     diurnal-bulge follow-up there). Implements Swinerd & Boulton
c     1983 eq. (28)-(29) (their own real-tracking best fit against
c     satellite 1963-27A found k_F=1.02 vs this J77-lineage formula,
c     i.e. the real diurnal amplitude is within 2% of what this
c     predicts): given the mean exospheric temperature (from
c     sw_tinf) and the Sun's current declination, returns the
c     GLOBAL max/min of the diurnal temperature distribution — the
c     inputs to the diurnal density bulge amplitude (mapped through
c     atm2d_interp by the caller; see bulge_cosphi for the separate
c     geometric factor). K_max/K_min and the 0.15/0.24 coefficients
c     are the paper's own literature-cited constants, not fitted
c     here.
c
c     tmean       : mean exospheric temperature, K (sw_tinf output)
c     decl_sun_deg: Sun declination, degrees
c     tmax, tmin  : output temperature extrema, K
c     ============================================================
      subroutine jacchia_bulge_extrema(tmean, decl_sun_deg,
     &                                  tmax, tmin)
      implicit none
      double precision tmean, decl_sun_deg, tmax, tmin

      double precision eps_ecl, ratio, kmax, kmin
      double precision phi_max, phi_min
      parameter (eps_ecl = 23.439d0)
      parameter (kmax = 0.5458d0, kmin = -0.5722d0)

      ratio = decl_sun_deg / eps_ecl

      phi_max = datan((5.d0/8.d0)*ratio/kmax)
      tmax = tmean*(1.d0 + 0.15d0*ratio*dsin(phi_max)
     &     + 0.24d0*kmax*dcos(phi_max))

      phi_min = datan((5.d0/8.d0)*ratio/kmin)
      tmin = tmean*(1.d0 + 0.15d0*ratio*dsin(phi_min)
     &     + 0.24d0*kmin*dcos(phi_min))

      return
      end

c     ============================================================
c     bulge_cosphi — geocentric angular factor for the diurnal
c     density bulge. Ported from OREM's src/propagate_ks.F (issue
c     #32 Phase 5 diurnal-bulge follow-up there). Returns cos(phi),
c     phi = angular distance from the satellite's current position
c     to the bulge-maximum direction: same declination as the Sun,
c     right ascension lagging the Sun by lam_deg (Swinerd & Boulton
c     1983's real-tracking fit, lam=37deg, matching CIRA 1972's
c     independent 14.1h local-solar-time estimate almost exactly).
c
c     x(3)    : satellite geocentric position (any consistent units)
c     ts(3)   : Sun geocentric position (same frame as x — already
c               computed via solarnpv for third-body/SRP physics)
c     lam_deg : bulge lag angle, degrees
c     ============================================================
      subroutine bulge_cosphi(x, ts, lam_deg, cosphi)
      implicit none
      double precision x(3), ts(3), lam_deg, cosphi

      double precision d2r_loc
      parameter (d2r_loc = 0.017453292519943295d0)

      double precision r_sat, r_sun, dec_sat, ra_sat
      double precision dec_sun, ra_sun, ra_bulge

      r_sat = dsqrt(x(1)**2 + x(2)**2 + x(3)**2)
      r_sun = dsqrt(ts(1)**2 + ts(2)**2 + ts(3)**2)

      dec_sat = dasin(x(3)/r_sat)
      ra_sat  = datan2(x(2), x(1))

      dec_sun = dasin(ts(3)/r_sun)
      ra_sun  = datan2(ts(2), ts(1))
      ra_bulge = ra_sun + lam_deg*d2r_loc

      cosphi = dsin(dec_sat)*dsin(dec_sun)
     &       + dcos(dec_sat)*dcos(dec_sun)*dcos(ra_sat - ra_bulge)

      return
      end
