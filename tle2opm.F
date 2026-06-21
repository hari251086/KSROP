c     ============================================================
c     tle2opm.F — Convert a TLE entry to CCSDS OPM v2.0 file
c
c     Reads a TLE file, selects the entry closest to a target epoch,
c     converts the TLE to a J2000 state vector via SGP4, computes
c     osculating orbital elements, and writes a CCSDS OPM file
c     suitable as input to driver_KS.exe.
c
c     Usage:
c       Reads configuration from input/tle2opm.cfg:
c         Line 1: TLE filename (relative path)
c         Line 2: target NORAD ID (integer, or -1 for first entry)
c         Line 3: target epoch YYYY MM DD HH MM SS.sss
c
c     Output:
c       input/KSROP_<EPOCH>_<NORAD>.opm
c       Also copies to input/input.opm for direct use by driver_KS
c
c     Compile:
c       ifort tle2opm.F Subrouts.F TLEread.F Legendre.F -o tle2opm.exe
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

      program tle2opm

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      parameter(maxtle = 100000)

      integer norad(maxtle), epoch_yr(maxtle), revnum(maxtle)
      double precision epoch_day(maxtle)
      double precision ndot_a(maxtle), nddot_a(maxtle), bstar_a(maxtle)
      double precision ainc_a(maxtle), raan_a(maxtle), ecc_a(maxtle)
      double precision aop_a(maxtle), ama_a(maxtle), amm_a(maxtle)
      character*24 names(maxtle)
      character*5  catid(maxtle)

      integer ntle, target_norad, ibest
      double precision r_j2k(3), v_j2k(3), pek(7)
      double precision cal_target(6), djd_target, djd_tle
      double precision diff_best, diff_cur

      character*120 tle_file
      character*23  epochstr
      character*60  opm_fname
      character*5   norad_str
      integer       ierr

c     --- Initialise constants ---
      pi      = 4.d0 * datan(1.d0)
      d2r     = pi / 180.d0
      r2d     = 180.d0 / pi
      amue    = 3.986004415d5
      AU      = 1.495978707d8
      R_Earth = 6378.1363d0

c     --- Read configuration ---
      open(unit=10, file='input/tle2opm.cfg', status='old', err=900)
      read(10, '(A)') tle_file
      read(10, *) target_norad
      read(10, *) cal_target(1), cal_target(2), cal_target(3),
     &            cal_target(4), cal_target(5), cal_target(6)
      close(10)

c     --- Compute target Julian Date ---
      call cal2jd(cal_target, djd_target)

      write(*,*) '=============================================='
      write(*,*) ' TLE to OPM Converter'
      write(*,*) '=============================================='
      write(*,'(A,A)') '  TLE file : ', trim(tle_file)
      write(*,'(A,I6)') '  Target NORAD : ', target_norad
      write(*,'(A,F12.4)') '  Target JD    : ', djd_target

c     --- Read all TLEs from file ---
      open(unit=11, file=tle_file, status='old', err=910)
      call read_tle(11, maxtle, names, catid, norad,
     &             epoch_yr, epoch_day,
     &             ndot_a, nddot_a, bstar_a,
     &             ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &             revnum, ntle)
      close(11)

      write(*,'(A,I6,A)') '  Read ', ntle, ' TLE entries'

      if (ntle .eq. 0) then
         write(*,*) 'ERROR: no TLE entries found'
         stop 1
      end if

c     --- Find entry closest to target epoch with matching NORAD ---
      ibest = 0
      diff_best = 1.d30

      do i = 1, ntle
         if (target_norad .gt. 0 .and. norad(i) .ne. target_norad)
     &      goto 20

         call tle_epoch2jd(epoch_yr(i), epoch_day(i), djd_tle)
         diff_cur = dabs(djd_tle - djd_target)

         if (diff_cur .lt. diff_best) then
            diff_best = diff_cur
            ibest = i
         end if

 20      continue
      end do

      if (ibest .eq. 0) then
         write(*,'(A,I6)') ' ERROR: no TLE found for NORAD ',
     &                     target_norad
         stop 1
      end if

      write(*,'(A,I6,A,I4,A,F12.6)')
     &   '  Selected entry ', ibest,
     &   '  epoch yr=', epoch_yr(ibest),
     &   '  day=', epoch_day(ibest)
      write(*,'(A,F8.4,A)') '  Epoch offset: ',
     &   diff_best, ' days from target'

c     --- Convert TLE to J2000 state vector via SGP4 ---
      call tle2sv(epoch_yr(ibest), epoch_day(ibest),
     &            ainc_a(ibest), raan_a(ibest), ecc_a(ibest),
     &            aop_a(ibest), ama_a(ibest), amm_a(ibest),
     &            bstar_a(ibest),
     &            r_j2k, v_j2k, ierr)

      if (ierr .ne. 0) then
         write(*,'(A,I2)') ' ERROR: tle2sv failed, ierr=', ierr
         stop 1
      end if

      write(*,'(A,3F14.6)') '  r_J2K (km)  : ', r_j2k
      write(*,'(A,3F14.9)') '  v_J2K (km/s): ', v_j2k

c     --- Compute osculating orbital elements ---
      call car2oe(r_j2k, v_j2k, pek)

      write(*,'(A,F12.4,A)') '  SMA  = ', pek(1), ' km'
      write(*,'(A,F12.9)')   '  ecc  = ', pek(2)
      write(*,'(A,F12.6,A)') '  inc  = ', pek(3), ' deg'
      write(*,'(A,F12.6,A)') '  RAAN = ', pek(4), ' deg'
      write(*,'(A,F12.6,A)') '  AOP  = ', pek(5), ' deg'
      write(*,'(A,F12.6,A)') '  MA   = ', pek(6), ' deg'

c     --- Build epoch string from the TLE epoch ---
      call tle_epoch2jd(epoch_yr(ibest), epoch_day(ibest), djd_tle)
      call jd2epoch(djd_tle, epochstr)

c     --- Build output filename: input/KSROP_<EPOCH>_<NORAD>.opm ---
      write(norad_str, '(I5.5)') norad(ibest)
      opm_fname = 'input/KSROP_'
     &   // epochstr(1:4) // epochstr(6:7) // epochstr(9:10)
     &   // 'T'
     &   // epochstr(12:13) // epochstr(15:16) // epochstr(18:19)
     &   // '_' // norad_str // '.opm'

      write(*,'(A,A)') '  Output: ', trim(opm_fname)

c     --- Write the OPM file ---
      open(unit=20, file=opm_fname, status='replace')
      call write_opm(20, epochstr, r_j2k, v_j2k, pek)
      close(20)

c     --- Also copy to input/input.opm for direct driver_KS use ---
      open(unit=21, file='input/input.opm', status='replace')
      call write_opm(21, epochstr, r_j2k, v_j2k, pek)
      close(21)

      write(*,*) ''
      write(*,*) '  Wrote OPM and copied to input/input.opm'
      write(*,*) '=============================================='

      stop
c     --- Error handlers ---
 900  write(*,*) 'ERROR: cannot open input/tle2opm.cfg'
      stop 1
 910  write(*,*) 'ERROR: cannot open TLE file: ', trim(tle_file)
      stop 1

      end
