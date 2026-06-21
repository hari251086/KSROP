c     test_tle2opm.F — Unit and regression tests for tle2opm pipeline
c     -------------------------------------------------------
c     Compile:  ifort test_tle2opm.F Subrouts.F TLEread.F Legendre.F
c                    -o test_tle2opm.exe
c     Run:      ./test_tle2opm.exe
c     -------------------------------------------------------
c     Tests:
c       1. TLE read from example_47944 finds entries
c       2. Epoch selection picks entry closest to 2025 May 1
c       3. tle2sv produces valid J2000 state (finite, physical)
c       4. car2oe produces valid elements from that state
c       5. write_opm / read_opm roundtrip recovers the state
c       6. Physical plausibility checks (altitude, velocity, period)
c       7. Epoch string format validation
c       8. jd2epoch roundtrip from tle_epoch2jd
c     -------------------------------------------------------
      program test_tle2opm

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

      integer ntle, npass, nfail, ibest, ierr
      double precision r_j2k(3), v_j2k(3), pek(7)
      double precision cal_target(6), djd_target, djd_tle
      double precision diff_best, diff_cur
      double precision r_read(3), v_read(3), cal_read(6)
      double precision rmag, vmag, alt, period

      character*23 epochstr

c     --- Initialise ---
      pi      = 4.d0 * datan(1.d0)
      d2r     = pi / 180.d0
      r2d     = 180.d0 / pi
      amue    = 3.986004415d5
      AU      = 1.495978707d8
      R_Earth = 6378.1363d0

      npass = 0
      nfail = 0

      write(*,*) '=============================================='
      write(*,*) ' TLE-to-OPM Pipeline Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     Test 1: Read TLE file
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- TLE file read ---'

      open(unit=11, file='input/example_47944.tle.txt', status='old',
     &     err=900)
      call read_tle(11, maxtle, names, catid, norad,
     &             epoch_yr, epoch_day,
     &             ndot_a, nddot_a, bstar_a,
     &             ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &             revnum, ntle)
      close(11)

      if (ntle .gt. 0) then
         write(*,'(A,I6)') '  [PASS]  Read TLE entries: ', ntle
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  No TLE entries read'
         nfail = nfail + 1
         goto 999
      end if

c     ==============================================================
c     Test 2: Epoch selection — find closest to 2025-05-01 12:00
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Epoch selection ---'

      cal_target(1) = 2025.d0
      cal_target(2) = 5.d0
      cal_target(3) = 1.d0
      cal_target(4) = 12.d0
      cal_target(5) = 0.d0
      cal_target(6) = 0.d0
      call cal2jd(cal_target, djd_target)

      ibest = 0
      diff_best = 1.d30
      do i = 1, ntle
         if (norad(i) .ne. 47944) goto 20
         call tle_epoch2jd(epoch_yr(i), epoch_day(i), djd_tle)
         diff_cur = dabs(djd_tle - djd_target)
         if (diff_cur .lt. diff_best) then
            diff_best = diff_cur
            ibest = i
         end if
 20      continue
      end do

      if (ibest .gt. 0) then
         write(*,'(A,I6,A,F8.4,A)')
     &      '  [PASS]  Best entry idx=', ibest,
     &      '  offset=', diff_best, ' days'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  No matching NORAD 47944 entry found'
         nfail = nfail + 1
         goto 999
      end if

c     Offset should be < 1 day (there are multiple entries per day)
      if (diff_best .lt. 1.d0) then
         write(*,*) '  [PASS]  Offset < 1 day from target'
         npass = npass + 1
      else
         write(*,'(A,F10.4)')
     &      '  [FAIL]  Offset too large: ', diff_best
         nfail = nfail + 1
      end if

c     Epoch year should be 2025
      if (epoch_yr(ibest) .eq. 2025) then
         write(*,*) '  [PASS]  Selected epoch year = 2025'
         npass = npass + 1
      else
         write(*,'(A,I5)')
     &      '  [FAIL]  Epoch year = ', epoch_yr(ibest)
         nfail = nfail + 1
      end if

c     ==============================================================
c     Test 3: tle2sv produces valid J2000 state
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- tle2sv conversion ---'

      call tle2sv(epoch_yr(ibest), epoch_day(ibest),
     &            ainc_a(ibest), raan_a(ibest), ecc_a(ibest),
     &            aop_a(ibest), ama_a(ibest), amm_a(ibest),
     &            bstar_a(ibest),
     &            r_j2k, v_j2k, ierr)

      if (ierr .eq. 0) then
         write(*,*) '  [PASS]  tle2sv ierr = 0'
         npass = npass + 1
      else
         write(*,'(A,I3)') '  [FAIL]  tle2sv ierr = ', ierr
         nfail = nfail + 1
         goto 999
      end if

c     Check finite (no NaN)
      ifinite = 1
      do i = 1, 3
         if (r_j2k(i) .ne. r_j2k(i)) ifinite = 0
         if (v_j2k(i) .ne. v_j2k(i)) ifinite = 0
      end do
      if (ifinite .eq. 1) then
         write(*,*) '  [PASS]  State vector is finite'
         npass = npass + 1
      else
         write(*,*) '  [FAIL]  State vector contains NaN'
         nfail = nfail + 1
      end if

c     ==============================================================
c     Test 4: Physical plausibility of state
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Physical plausibility ---'

      rmag = dsqrt(r_j2k(1)**2 + r_j2k(2)**2 + r_j2k(3)**2)
      vmag = dsqrt(v_j2k(1)**2 + v_j2k(2)**2 + v_j2k(3)**2)
      alt = rmag - R_Earth

c     Altitude in LEO range [100, 800] km for 47944
      if (alt .gt. 100.d0 .and. alt .lt. 800.d0) then
         write(*,'(A,F10.2,A)')
     &      '  [PASS]  Altitude = ', alt, ' km (LEO range)'
         npass = npass + 1
      else
         write(*,'(A,F10.2,A)')
     &      '  [FAIL]  Altitude = ', alt, ' km (outside LEO)'
         nfail = nfail + 1
      end if

c     Velocity in [7.0, 8.5] km/s for LEO
      if (vmag .gt. 7.d0 .and. vmag .lt. 8.5d0) then
         write(*,'(A,F10.6,A)')
     &      '  [PASS]  Velocity = ', vmag, ' km/s (LEO range)'
         npass = npass + 1
      else
         write(*,'(A,F10.6,A)')
     &      '  [FAIL]  Velocity = ', vmag, ' km/s (outside LEO)'
         nfail = nfail + 1
      end if

c     Radius in [6400, 7200] km
      if (rmag .gt. 6400.d0 .and. rmag .lt. 7200.d0) then
         write(*,'(A,F10.2,A)')
     &      '  [PASS]  Radius = ', rmag, ' km'
         npass = npass + 1
      else
         write(*,'(A,F10.2,A)')
     &      '  [FAIL]  Radius = ', rmag, ' km'
         nfail = nfail + 1
      end if

c     ==============================================================
c     Test 5: car2oe produces valid elements
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Orbital elements ---'

      call car2oe(r_j2k, v_j2k, pek)

c     SMA should be in [6500, 7200] km for LEO
      if (pek(1) .gt. 6500.d0 .and. pek(1) .lt. 7200.d0) then
         write(*,'(A,F10.4,A)')
     &      '  [PASS]  SMA = ', pek(1), ' km'
         npass = npass + 1
      else
         write(*,'(A,F10.4,A)')
     &      '  [FAIL]  SMA = ', pek(1), ' km'
         nfail = nfail + 1
      end if

c     Eccentricity should be < 0.1 (near-circular)
      if (pek(2) .ge. 0.d0 .and. pek(2) .lt. 0.1d0) then
         write(*,'(A,F12.9)')
     &      '  [PASS]  ecc = ', pek(2)
         npass = npass + 1
      else
         write(*,'(A,F12.9)')
     &      '  [FAIL]  ecc = ', pek(2)
         nfail = nfail + 1
      end if

c     Inclination should be ~97 deg (SSO) for 47944
      if (pek(3) .gt. 95.d0 .and. pek(3) .lt. 100.d0) then
         write(*,'(A,F10.4,A)')
     &      '  [PASS]  inc = ', pek(3), ' deg (SSO)'
         npass = npass + 1
      else
         write(*,'(A,F10.4,A)')
     &      '  [FAIL]  inc = ', pek(3), ' deg'
         nfail = nfail + 1
      end if

c     RAAN in [0, 360]
      if (pek(4) .ge. 0.d0 .and. pek(4) .le. 360.d0) then
         write(*,'(A,F10.4,A)')
     &      '  [PASS]  RAAN = ', pek(4), ' deg'
         npass = npass + 1
      else
         write(*,'(A,F10.4)')
     &      '  [FAIL]  RAAN = ', pek(4)
         nfail = nfail + 1
      end if

c     Period check from mean motion (should be ~90-100 min for LEO)
      period = 1440.d0 / amm_a(ibest)
      if (period .gt. 85.d0 .and. period .lt. 110.d0) then
         write(*,'(A,F8.2,A)')
     &      '  [PASS]  Period = ', period, ' min (LEO)'
         npass = npass + 1
      else
         write(*,'(A,F8.2,A)')
     &      '  [FAIL]  Period = ', period, ' min'
         nfail = nfail + 1
      end if

c     ==============================================================
c     Test 6: jd2epoch roundtrip
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Epoch string conversion ---'

      call tle_epoch2jd(epoch_yr(ibest), epoch_day(ibest), djd_tle)
      call jd2epoch(djd_tle, epochstr)

c     Year should be 2025
      if (epochstr(1:4) .eq. '2025') then
         write(*,'(A,A)')
     &      '  [PASS]  Epoch string year = ', epochstr(1:4)
         npass = npass + 1
      else
         write(*,'(A,A)')
     &      '  [FAIL]  Epoch string year = ', epochstr(1:4)
         nfail = nfail + 1
      end if

c     Month should be 05 (May) — day 121 = May 1
      if (epochstr(6:7) .eq. '05') then
         write(*,'(A,A)')
     &      '  [PASS]  Epoch string month = ', epochstr(6:7)
         npass = npass + 1
      else
         write(*,'(A,A)')
     &      '  [FAIL]  Epoch string month = ', epochstr(6:7)
         nfail = nfail + 1
      end if

c     Full string should be well-formed (YYYY-MM-DDTHH:MM:SS.mmm)
      if (epochstr(5:5) .eq. '-' .and. epochstr(8:8) .eq. '-'
     &    .and. epochstr(11:11) .eq. 'T'
     &    .and. epochstr(14:14) .eq. ':'
     &    .and. epochstr(17:17) .eq. ':') then
         write(*,'(A,A)')
     &      '  [PASS]  Epoch format valid: ', epochstr
         npass = npass + 1
      else
         write(*,'(A,A)')
     &      '  [FAIL]  Epoch format invalid: ', epochstr
         nfail = nfail + 1
      end if

c     ==============================================================
c     Test 7: write_opm / read_opm roundtrip
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- OPM write/read roundtrip ---'

      open(unit=30, file='output/test_tle2opm_rt.opm',
     &     status='replace')
      call write_opm(30, epochstr, r_j2k, v_j2k, pek)
      close(30)

      open(unit=30, file='output/test_tle2opm_rt.opm', status='old')
      call read_opm(30, r_read, v_read, cal_read)
      close(30)

c     Position roundtrip (tolerance 1e-4 km due to format width)
      dr_rt = dsqrt((r_read(1)-r_j2k(1))**2
     &            + (r_read(2)-r_j2k(2))**2
     &            + (r_read(3)-r_j2k(3))**2)
      if (dr_rt .lt. 1.d-4) then
         write(*,'(A,E10.3,A)')
     &      '  [PASS]  Position roundtrip dr=', dr_rt, ' km'
         npass = npass + 1
      else
         write(*,'(A,E10.3,A)')
     &      '  [FAIL]  Position roundtrip dr=', dr_rt, ' km'
         nfail = nfail + 1
      end if

c     Velocity roundtrip (tolerance 1e-7 km/s)
      dv_rt = dsqrt((v_read(1)-v_j2k(1))**2
     &            + (v_read(2)-v_j2k(2))**2
     &            + (v_read(3)-v_j2k(3))**2)
      if (dv_rt .lt. 1.d-7) then
         write(*,'(A,E10.3,A)')
     &      '  [PASS]  Velocity roundtrip dv=', dv_rt, ' km/s'
         npass = npass + 1
      else
         write(*,'(A,E10.3,A)')
     &      '  [FAIL]  Velocity roundtrip dv=', dv_rt, ' km/s'
         nfail = nfail + 1
      end if

c     Epoch year roundtrip
      if (nint(cal_read(1)) .eq. 2025) then
         write(*,*) '  [PASS]  Epoch year roundtrip = 2025'
         npass = npass + 1
      else
         write(*,'(A,F8.0)')
     &      '  [FAIL]  Epoch year roundtrip = ', cal_read(1)
         nfail = nfail + 1
      end if

c     ==============================================================
c     Test 8: Orbital energy is negative (bound orbit)
c     ==============================================================
      write(*,*) ''
      write(*,*) '  --- Energy check ---'

      energy = 0.5d0*(v_j2k(1)**2+v_j2k(2)**2+v_j2k(3)**2)
     &       - amue/rmag
      if (energy .lt. 0.d0) then
         write(*,'(A,F12.4,A)')
     &      '  [PASS]  Orbital energy = ', energy, ' km2/s2 (bound)'
         npass = npass + 1
      else
         write(*,'(A,F12.4)')
     &      '  [FAIL]  Orbital energy = ', energy
         nfail = nfail + 1
      end if

c     ==============================================================
c     Summary
c     ==============================================================
 999  continue
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,200) npass + nfail, npass, nfail
 200  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='

      if (nfail .gt. 0) stop 1
      stop

 900  write(*,*) '[FAIL] Cannot open input/example_47944.tle.txt'
      stop 1
      end
