c     ============================================================
c     test_tle.F — Comprehensive tests for TLEread.F
c
c     Coverage:
c       tle_expval   : positive/negative/zero mantissa, +/- exponent,
c                      explicit + sign, large mantissa, zero exponent
c       tle_chksum   : valid lines, corrupted digit, actual file lines,
c                      lines with dashes, checksum digit value
c       tle_parse1   : NORAD ID (numeric + alpha-numeric), epoch year
c                      pivot (2000s and 1900s), epoch day, ndot
c                      (positive + negative), non-zero nddot, bstar,
c                      catid string
c       tle_parse2   : inclination, RAAN, eccentricity (high, low,
c                      typical), AOP, mean anomaly, mean motion,
c                      revolution number
c       read_tle     : two-line single, three-line single, mixed
c                      2/3-line multi-entry, blank lines, orphan
c                      line 1 (no pair), orphan line 2 (no line 1),
c                      consecutive line 1s (overwrite), maxtle
c                      buffer limit, example_47944.tle.txt (3977
c                      entries), example_multi.tle.txt (94597
c                      entries with alpha-numeric IDs)
c       cross-valid  : checksum every entry in both files, epoch
c                      years plausible, inclinations in [0,180],
c                      eccentricities in [0,1), mean motions > 0,
c                      monotonic epoch in single-satellite file,
c                      spot-check specific entries deep in files,
c                      ndot sign, physical orbital period range
c     ============================================================

      program test_tle
      implicit double precision (a-h, o-z)
      parameter(maxtle = 100000)

      integer norad(maxtle), epoch_yr(maxtle), revnum(maxtle)
      double precision epoch_day(maxtle)
      double precision ndot_a(maxtle), nddot_a(maxtle), bstar_a(maxtle)
      double precision ainc_a(maxtle), raan_a(maxtle), ecc_a(maxtle)
      double precision aop_a(maxtle), ama_a(maxtle), amm_a(maxtle)
      character*24 names(maxtle)
      character*5  catid(maxtle)

      integer ntle, npass, nfail, ick, iok, ierr
      double precision val, tol
      character*69 line1, line2
      character*80 rawline

      npass = 0
      nfail = 0
      tol = 1.d-10

      write(*,*) '=============================================='
      write(*,*) ' TLE Reader Test Suite'
      write(*,*) '=============================================='

c     ============================================================
c     tle_expval: implied-decimal exponent decoding
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle_expval ---'

c     T1: positive mantissa, negative exponent
      call tle_expval(' 37765-3', val)
      call assert_d('expval  37765-3 = 3.7765e-4',
     &              val, 3.7765d-4, tol, npass, nfail)

c     T2: zero mantissa, positive exponent
      call tle_expval(' 00000+0', val)
      call assert_d('expval  00000+0 = 0.0',
     &              val, 0.d0, tol, npass, nfail)

c     T3: zero mantissa, negative exponent
      call tle_expval(' 00000-0', val)
      call assert_d('expval  00000-0 = 0.0',
     &              val, 0.d0, tol, npass, nfail)

c     T4: negative mantissa, negative exponent
      call tle_expval('-11606-4', val)
      call assert_d('expval -11606-4 = -1.1606e-5',
     &              val, -1.1606d-5, tol, npass, nfail)

c     T5: small positive mantissa
      call tle_expval(' 26175-5', val)
      call assert_d('expval  26175-5 = 2.6175e-6',
     &              val, 2.6175d-6, tol, npass, nfail)

c     T6: positive exponent
      call tle_expval(' 12345+3', val)
      call assert_d('expval  12345+3 = 123.45',
     &              val, 123.45d0, 1.d-6, npass, nfail)

c     T7: explicit + sign on mantissa
      call tle_expval('+54321-2', val)
      call assert_d('expval +54321-2 = 5.4321e-3',
     &              val, 5.4321d-3, tol, npass, nfail)

c     T8: large mantissa, zero exponent
      call tle_expval(' 99999+0', val)
      call assert_d('expval  99999+0 = 0.99999',
     &              val, 0.99999d0, tol, npass, nfail)

c     T9: negative mantissa, positive exponent
      call tle_expval('-10000+1', val)
      call assert_d('expval -10000+1 = -1.0',
     &              val, -1.0d0, tol, npass, nfail)

c     T10: single non-zero digit
      call tle_expval(' 10000-4', val)
      call assert_d('expval  10000-4 = 1.0e-5',
     &              val, 1.0d-5, tol, npass, nfail)

c     ============================================================
c     tle_chksum: modulo-10 checksum verification
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle_chksum ---'

c     T11: valid line 1 (NOAA-18)
      line1 = '1 28654U 05018A   26170.52919894  .0000003'
     &     // '5  00000-0  32076-4 0  9992'
      call tle_chksum(line1, ick, iok)
      call assert_i('chksum NOAA-18 line1 valid',
     &              iok, 1, npass, nfail)

c     T12: valid line 2 (NOAA-18)
      line2 = '2 28654  99.0357 195.5267 0014319  94.5046'
     &     // ' 265.7602 14.12818068 77721'
      call tle_chksum(line2, ick, iok)
      call assert_i('chksum NOAA-18 line2 valid',
     &              iok, 1, npass, nfail)

c     T13: checksum digit value for line2
      call assert_i('chksum NOAA-18 line2 digit = 1',
     &              ick, 1, npass, nfail)

c     T14: corrupted line (wrong last digit)
      line1(69:69) = '7'
      call tle_chksum(line1, ick, iok)
      call assert_i('chksum corrupted -> iok=0',
     &              iok, 0, npass, nfail)

c     T15: actual line from 47944 file (with negative ndot, dashes)
c     1 47944U 21022N   21083.05205046 -.00000063  00000-0  00000+0 0  9990
      line1 = '1 47944U 21022N   21083.05205046 -.0000006'
     &     // '3  00000-0  00000+0 0  9990'
      call tle_chksum(line1, ick, iok)
      call assert_i('chksum 47944 line1 (neg ndot) valid',
     &              iok, 1, npass, nfail)

c     T16: its paired line 2
      line2 = '2 47944  97.5667 346.8994 0016808 232.1003'
     &     // ' 185.5158 15.06560546   211'
      call tle_chksum(line2, ick, iok)
      call assert_i('chksum 47944 line2 valid',
     &              iok, 1, npass, nfail)

c     T17: line from multi file (Vanguard 1, NORAD 5)
      line1 = '1 00005U 58002B   26170.40233228  .0000028'
     &     // '9  00000-0  37765-3 0  9990'
      call tle_chksum(line1, ick, iok)
      call assert_i('chksum Vanguard-1 line1 valid',
     &              iok, 1, npass, nfail)

c     ============================================================
c     tle_parse1: field extraction
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle_parse1 ---'

c     Standard numeric NORAD (NOAA-18, 28654)
      line1 = '1 28654U 05018A   26170.52919894  .0000003'
     &     // '5  00000-0  32076-4 0  9992'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)

c     T18: NORAD ID
      call assert_i('parse1 NORAD = 28654',
     &              norad(1), 28654, npass, nfail)

c     T19: catid string
      call assert_s('parse1 catid = 28654',
     &              catid(1), '28654', npass, nfail)

c     T20: epoch year (26 -> 2026)
      call assert_i('parse1 epoch_yr = 2026',
     &              epoch_yr(1), 2026, npass, nfail)

c     T21: epoch day
      call assert_d('parse1 epoch_day = 170.52919894',
     &              epoch_day(1), 170.52919894d0, 1.d-8, npass, nfail)

c     T22: positive ndot
      call assert_d('parse1 ndot = 0.00000035',
     &              ndot_a(1), 0.00000035d0, tol, npass, nfail)

c     T23: zero nddot
      call assert_d('parse1 nddot = 0.0',
     &              nddot_a(1), 0.d0, tol, npass, nfail)

c     T24: bstar
      call assert_d('parse1 bstar = 3.2076e-5',
     &              bstar_a(1), 3.2076d-5, tol, npass, nfail)

c     Negative ndot entry (47944 early epoch)
      line1 = '1 47944U 21022N   21083.05205046 -.0000006'
     &     // '3  00000-0  00000+0 0  9990'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)

c     T25: negative ndot
      call assert_d('parse1 neg ndot = -0.00000063',
     &              ndot_a(1), -0.00000063d0, tol, npass, nfail)

c     T26: epoch year pivot (21 -> 2021)
      call assert_i('parse1 epoch_yr 21 -> 2021',
     &              epoch_yr(1), 2021, npass, nfail)

c     T27: zero bstar (from 00000+0)
      call assert_d('parse1 bstar = 0.0 (00000+0)',
     &              bstar_a(1), 0.d0, tol, npass, nfail)

c     1950s epoch year (Vanguard 1, launch year 58 -> 1958)
      line1 = '1 00005U 58002B   26170.40233228  .0000028'
     &     // '9  00000-0  37765-3 0  9990'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)

c     T28: epoch year (Vanguard TLE epoch 26 -> 2026)
      call assert_i('parse1 Vanguard epoch_yr = 2026',
     &              epoch_yr(1), 2026, npass, nfail)

c     T29: small NORAD padded with zeros
      call assert_i('parse1 Vanguard NORAD = 5',
     &              norad(1), 5, npass, nfail)

c     T30: catid preserves leading zeros
      call assert_s('parse1 Vanguard catid = 00005',
     &              catid(1), '00005', npass, nfail)

c     Non-zero nddot (47944 late epoch with large ndot)
c     1 47944U 21022N   25208.08235808  .09897846  25875-5  31149-3 0  9990
      line1 = '1 47944U 21022N   25208.08235808  .098978'
     &     // '46  25875-5  31149-3 0  9990'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)

c     T31: large positive ndot
      call assert_d('parse1 large ndot = 0.09897846',
     &              ndot_a(1), 0.09897846d0, tol, npass, nfail)

c     T32: non-zero nddot (25875-5 -> 0.25875e-5)
      call assert_d('parse1 nddot = 2.5875e-6',
     &              nddot_a(1), 2.5875d-6, tol, npass, nfail)

c     T33: bstar from this entry (31149-3 -> 0.31149e-3)
      call assert_d('parse1 bstar = 3.1149e-4',
     &              bstar_a(1), 3.1149d-4, tol, npass, nfail)

c     Alpha-numeric NORAD (T0003)
      line1 = '1 T0003U          26169.78100846  .0000126'
     &     // '0  00000-0  20524-2 0  9995'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)

c     T34: alpha-numeric NORAD returns -1
      call assert_i('parse1 alpha NORAD = -1',
     &              norad(1), -1, npass, nfail)

c     T35: catid preserves alpha string
      call assert_s('parse1 alpha catid = T0003',
     &              catid(1), 'T0003', npass, nfail)

c     ============================================================
c     tle_parse2: field extraction
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- tle_parse2 ---'

c     Standard entry (NOAA-18)
      line2 = '2 28654  99.0357 195.5267 0014319  94.5046'
     &     // ' 265.7602 14.12818068 77721'
      call tle_parse2(line2, ainc_a(1), raan_a(1), ecc_a(1),
     &                aop_a(1), ama_a(1), amm_a(1), revnum(1), ierr)

c     T36: inclination
      call assert_d('parse2 inc = 99.0357',
     &              ainc_a(1), 99.0357d0, 1.d-4, npass, nfail)

c     T37: RAAN
      call assert_d('parse2 raan = 195.5267',
     &              raan_a(1), 195.5267d0, 1.d-4, npass, nfail)

c     T38: eccentricity (low)
      call assert_d('parse2 ecc = 0.0014319',
     &              ecc_a(1), 0.0014319d0, 1.d-7, npass, nfail)

c     T39: AOP
      call assert_d('parse2 aop = 94.5046',
     &              aop_a(1), 94.5046d0, 1.d-4, npass, nfail)

c     T40: mean anomaly
      call assert_d('parse2 ma = 265.7602',
     &              ama_a(1), 265.7602d0, 1.d-4, npass, nfail)

c     T41: mean motion
      call assert_d('parse2 mm = 14.12818068',
     &              amm_a(1), 14.12818068d0, 1.d-8, npass, nfail)

c     T42: revolution number
      call assert_i('parse2 revnum = 7772',
     &              revnum(1), 7772, npass, nfail)

c     High eccentricity (Vanguard 1, e=0.1835576)
      line2 = '2 00005  34.2451 307.4935 1835576 148.4938'
     &     // ' 223.9981 10.86000951443521'
      call tle_parse2(line2, ainc_a(1), raan_a(1), ecc_a(1),
     &                aop_a(1), ama_a(1), amm_a(1), revnum(1), ierr)

c     T43: high eccentricity
      call assert_d('parse2 high ecc = 0.1835576',
     &              ecc_a(1), 0.1835576d0, 1.d-7, npass, nfail)

c     T44: low inclination
      call assert_d('parse2 inc = 34.2451',
     &              ainc_a(1), 34.2451d0, 1.d-4, npass, nfail)

c     T45: low mean motion (long period)
      call assert_d('parse2 mm = 10.86000951',
     &              amm_a(1), 10.86000951d0, 1.d-8, npass, nfail)

c     T46: large revolution number
      call assert_i('parse2 revnum = 44352',
     &              revnum(1), 44352, npass, nfail)

c     Very high eccentricity (NORAD 16, e=0.2016745)
      line2 = '2 00016  34.2627 300.7118 2016745 126.8940'
     &     // ' 253.3932 10.49596045689078'
      call tle_parse2(line2, ainc_a(1), raan_a(1), ecc_a(1),
     &                aop_a(1), ama_a(1), amm_a(1), revnum(1), ierr)

c     T47: very high eccentricity
      call assert_d('parse2 vhigh ecc = 0.2016745',
     &              ecc_a(1), 0.2016745d0, 1.d-7, npass, nfail)

c     T48: large revnum
      call assert_i('parse2 revnum = 68907',
     &              revnum(1), 68907, npass, nfail)

c     ============================================================
c     read_tle: two-line single entry
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- read_tle: inline files ---'

      open(unit=90, status='scratch')
      write(90,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(90,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(90)
      call read_tle(90, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(90)

c     T49
      call assert_i('2-line ntle = 1', ntle, 1, npass, nfail)
c     T50
      call assert_i('2-line NORAD = 25544',
     &              norad(1), 25544, npass, nfail)
c     T51
      call assert_s('2-line name blank',
     &              names(1)(1:1), ' ', npass, nfail)
c     T52
      call assert_d('2-line inc = 51.64',
     &              ainc_a(1), 51.64d0, 1.d-4, npass, nfail)
c     T53
      call assert_d('2-line ecc = 0.00067',
     &              ecc_a(1), 0.00067d0, 1.d-7, npass, nfail)
c     T54
      call assert_d('2-line mm = 15.54',
     &              amm_a(1), 15.54d0, 1.d-8, npass, nfail)

c     ============================================================
c     read_tle: three-line single entry (with name)
c     ============================================================
      open(unit=91, status='scratch')
      write(91,'(A)') 'ISS (ZARYA)'
      write(91,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(91,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(91)
      call read_tle(91, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(91)

c     T55
      call assert_i('3-line ntle = 1', ntle, 1, npass, nfail)
c     T56
      call assert_s('3-line name = ISS (ZARYA)',
     &              names(1)(1:11), 'ISS (ZARYA)', npass, nfail)
c     T57
      call assert_i('3-line NORAD = 25544',
     &              norad(1), 25544, npass, nfail)
c     T58
      call assert_i('3-line epoch_yr = 2026',
     &              epoch_yr(1), 2026, npass, nfail)
c     T59
      call assert_d('3-line bstar = 1.027e-4',
     &              bstar_a(1), 1.027d-4, tol, npass, nfail)
c     T60
      call assert_i('3-line revnum = 45678',
     &              revnum(1), 45678, npass, nfail)

c     ============================================================
c     read_tle: mixed two- and three-line entries
c     ============================================================
      open(unit=92, status='scratch')
c     Entry 1: three-line
      write(92,'(A)') 'VANGUARD 1'
      write(92,'(A)') '1 00005U 58002B   26170.40233228  .000002'
     &     // '89  00000-0  37765-3 0  9990'
      write(92,'(A)') '2 00005  34.2451 307.4935 1835576 148.493'
     &     // '8 223.9981 10.86000951443521'
c     Entry 2: two-line (no name)
      write(92,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(92,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
c     Entry 3: three-line
      write(92,'(A)') 'NOAA 18'
      write(92,'(A)') '1 28654U 05018A   26170.52919894  .000000'
     &     // '35  00000-0  32076-4 0  9992'
      write(92,'(A)') '2 28654  99.0357 195.5267 0014319  94.504'
     &     // '6 265.7602 14.12818068 77721'
      rewind(92)
      call read_tle(92, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(92)

c     T61
      call assert_i('mixed ntle = 3', ntle, 3, npass, nfail)
c     T62: first entry has name
      call assert_s('mixed entry1 name = VANGUARD 1',
     &              names(1)(1:10), 'VANGUARD 1', npass, nfail)
c     T63: first entry NORAD
      call assert_i('mixed entry1 NORAD = 5',
     &              norad(1), 5, npass, nfail)
c     T64: second entry has blank name (two-line)
      call assert_s('mixed entry2 name blank',
     &              names(2)(1:1), ' ', npass, nfail)
c     T65: second entry NORAD
      call assert_i('mixed entry2 NORAD = 25544',
     &              norad(2), 25544, npass, nfail)
c     T66: third entry has name
      call assert_s('mixed entry3 name = NOAA 18',
     &              names(3)(1:7), 'NOAA 18', npass, nfail)
c     T67: third entry NORAD
      call assert_i('mixed entry3 NORAD = 28654',
     &              norad(3), 28654, npass, nfail)

c     ============================================================
c     read_tle: blank lines between entries
c     ============================================================
      open(unit=93, status='scratch')
      write(93,'(A)') 'SAT A'
      write(93,'(A)') '1 00005U 58002B   26170.40233228  .000002'
     &     // '89  00000-0  37765-3 0  9990'
      write(93,'(A)') '2 00005  34.2451 307.4935 1835576 148.493'
     &     // '8 223.9981 10.86000951443521'
      write(93,'(A)') ''
      write(93,'(A)') ''
      write(93,'(A)') 'SAT B'
      write(93,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(93,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(93)
      call read_tle(93, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(93)

c     T68
      call assert_i('blank-lines ntle = 2', ntle, 2, npass, nfail)
c     T69
      call assert_s('blank-lines entry2 name = SAT B',
     &              names(2)(1:5), 'SAT B', npass, nfail)

c     ============================================================
c     read_tle: orphan line 2 (no preceding line 1) — skip
c     ============================================================
      open(unit=94, status='scratch')
c     Orphan line 2 should be ignored
      write(94,'(A)') '2 99999  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
c     Valid entry
      write(94,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(94,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(94)
      call read_tle(94, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(94)

c     T70
      call assert_i('orphan-line2 ntle = 1', ntle, 1, npass, nfail)
c     T71
      call assert_i('orphan-line2 NORAD = 25544',
     &              norad(1), 25544, npass, nfail)

c     ============================================================
c     read_tle: orphan line 1 at EOF (no line 2) — skip
c     ============================================================
      open(unit=95, status='scratch')
c     Valid entry
      write(95,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(95,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
c     Orphan line 1 at end of file
      write(95,'(A)') '1 99999U 00000A   26001.00000000  .000000'
     &     // '00  00000-0  00000+0 0  9999'
      rewind(95)
      call read_tle(95, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(95)

c     T72
      call assert_i('orphan-line1-EOF ntle = 1', ntle, 1, npass, nfail)

c     ============================================================
c     read_tle: consecutive line 1s (second overwrites first)
c     ============================================================
      open(unit=96, status='scratch')
c     First line 1 (will be overwritten)
      write(96,'(A)') '1 99999U 00000A   26001.00000000  .000000'
     &     // '00  00000-0  00000+0 0  9999'
c     Second line 1 (overwrites first)
      write(96,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
c     Line 2 pairs with second line 1
      write(96,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(96)
      call read_tle(96, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(96)

c     T73
      call assert_i('consec-line1 ntle = 1', ntle, 1, npass, nfail)
c     T74: should have the second line 1's NORAD
      call assert_i('consec-line1 NORAD = 25544',
     &              norad(1), 25544, npass, nfail)

c     ============================================================
c     read_tle: maxtle buffer limit
c     ============================================================
      open(unit=97, status='scratch')
      do i = 1, 5
         write(97,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &        // '17  00000-0  10270-3 0  9993'
         write(97,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &        // '4 315.0000 15.54000000456789'
      end do
      rewind(97)
c     Read with maxtle=3 (should stop at 3)
      call read_tle(97, 3, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(97)

c     T75
      call assert_i('maxtle=3 ntle = 3', ntle, 3, npass, nfail)

c     ============================================================
c     Gap 1: Epoch year pivot boundary (56 -> 2056, 57 -> 1957)
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- epoch year pivot boundary ---'

c     Year 57 -> 1957
      line1 = '1 00001U 57001A   57001.00000000  .0000000'
     &     // '0  00000-0  00000+0 0  9001'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)
      call assert_i('pivot yr 57 -> 1957',
     &              epoch_yr(1), 1957, npass, nfail)

c     Year 99 -> 1999
      line1 = '1 00001U 57001A   99365.00000000  .0000000'
     &     // '0  00000-0  00000+0 0  9001'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)
      call assert_i('pivot yr 99 -> 1999',
     &              epoch_yr(1), 1999, npass, nfail)

c     Year 00 -> 2000
      line1 = '1 00001U 57001A   00001.00000000  .0000000'
     &     // '0  00000-0  00000+0 0  9001'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)
      call assert_i('pivot yr 00 -> 2000',
     &              epoch_yr(1), 2000, npass, nfail)

c     Year 56 -> 2056
      line1 = '1 00001U 57001A   56001.00000000  .0000000'
     &     // '0  00000-0  00000+0 0  9001'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)
      call assert_i('pivot yr 56 -> 2056',
     &              epoch_yr(1), 2056, npass, nfail)

c     ============================================================
c     Gap 2: Name lines starting with '1' or '2' not misidentified
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- name-line false-positive rejection ---'

c     Name "1ST CUBESAT" should NOT be parsed as TLE line 1
c     (col 2 = 'S', not space)
      open(unit=98, status='scratch')
      write(98,'(A)') '1ST CUBESAT'
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('name "1ST CUBESAT" ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_s('name "1ST CUBESAT" captured',
     &              names(1)(1:11), '1ST CUBESAT', npass, nfail)
      call assert_i('name "1ST CUBESAT" NORAD = 25544',
     &              norad(1), 25544, npass, nfail)

c     Name "2ND STAGE DEB" should NOT be parsed as TLE line 2
      open(unit=98, status='scratch')
      write(98,'(A)') '2ND STAGE DEB'
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('name "2ND STAGE DEB" ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_s('name "2ND STAGE DEB" captured',
     &              names(1)(1:13), '2ND STAGE DEB', npass, nfail)

c     ============================================================
c     Gap 3: Consecutive name lines (second overwrites first)
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- consecutive name lines ---'

      open(unit=98, status='scratch')
      write(98,'(A)') 'OLD NAME'
      write(98,'(A)') 'NEW NAME'
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('consec-names ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_s('consec-names uses NEW NAME',
     &              names(1)(1:8), 'NEW NAME', npass, nfail)

c     ============================================================
c     Gap 4: Line 1/Line 2 catalog ID mismatch — skip pair
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- catalog ID mismatch ---'

c     Mismatched pair followed by valid pair
      open(unit=98, status='scratch')
c     Line 1 for NORAD 25544
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
c     Line 2 for NORAD 28654 (mismatch!)
      write(98,'(A)') '2 28654  99.0357 195.5267 0014319  94.504'
     &     // '6 265.7602 14.12818068 77721'
c     Valid pair
      write(98,'(A)') '1 00005U 58002B   26170.40233228  .000002'
     &     // '89  00000-0  37765-3 0  9990'
      write(98,'(A)') '2 00005  34.2451 307.4935 1835576 148.493'
     &     // '8 223.9981 10.86000951443521'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('catid-mismatch ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_i('catid-mismatch NORAD = 5',
     &              norad(1), 5, npass, nfail)

c     All-mismatch file produces zero entries
      open(unit=98, status='scratch')
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 28654  99.0357 195.5267 0014319  94.504'
     &     // '6 265.7602 14.12818068 77721'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('all-mismatch ntle = 0', ntle, 0, npass, nfail)

c     ============================================================
c     Gap 5: Short lines (< 69 chars) rejected as TLE data
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- short line rejection ---'

c     Truncated line 1 (40 chars) should be treated as name, not TLE
      open(unit=98, status='scratch')
      write(98,'(A)') '1 25544U 98067A   26170.50000000'
      write(98,'(A)') '1 00005U 58002B   26170.40233228  .000002'
     &     // '89  00000-0  37765-3 0  9990'
      write(98,'(A)') '2 00005  34.2451 307.4935 1835576 148.493'
     &     // '8 223.9981 10.86000951443521'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('short-line1 ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_i('short-line1 NORAD = 5',
     &              norad(1), 5, npass, nfail)

c     Truncated line 2 — valid line 1 but short line 2, should skip
      open(unit=98, status='scratch')
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567'
c     Then a valid pair
      write(98,'(A)') '1 00005U 58002B   26170.40233228  .000002'
     &     // '89  00000-0  37765-3 0  9990'
      write(98,'(A)') '2 00005  34.2451 307.4935 1835576 148.493'
     &     // '8 223.9981 10.86000951443521'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('short-line2 ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_i('short-line2 NORAD = 5',
     &              norad(1), 5, npass, nfail)

c     ============================================================
c     Gap 2 (extended): col-8 classification check on line 1
c     A 69-char string starting with "1 " but with non-U/C/S at
c     col 8 should be rejected as a TLE line 1.
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- col-8 classification check ---'

c     69 chars, starts with '1 ', but col 8 = 'X' (invalid class)
      open(unit=98, status='scratch')
      write(98,'(A)') '1 FAKE'
     &     // 'X THIS IS NOT A TLE LINE BUT IS EXACTLY 69 CHARS LONG'
     &     // ' PADDED!!'
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('col8-invalid ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_i('col8-invalid NORAD = 25544',
     &              norad(1), 25544, npass, nfail)

c     Valid classification letters: U, C, S, space
      open(unit=98, status='scratch')
      write(98,'(A)') '1 25544C 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('col8-C accepted ntle = 1',
     &              ntle, 1, npass, nfail)

      open(unit=98, status='scratch')
      write(98,'(A)') '1 25544S 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('col8-S accepted ntle = 1',
     &              ntle, 1, npass, nfail)

c     ============================================================
c     Corrupt field handling: parse errors don't crash, entry skipped
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- corrupt field handling ---'

c     tle_parse1 with garbage in epoch year field -> ierr=1
      line1 = '1 25544U 98067A   XX170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)
      call assert_i('parse1 corrupt epoch -> ierr=1',
     &              ierr, 1, npass, nfail)

c     tle_parse1 with valid data -> ierr=0
      line1 = '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      call tle_parse1(line1, catid(1), norad(1), epoch_yr(1),
     &          epoch_day(1), ndot_a(1), nddot_a(1), bstar_a(1), ierr)
      call assert_i('parse1 valid -> ierr=0',
     &              ierr, 0, npass, nfail)

c     tle_parse2 with garbage in inclination field -> ierr=1
      line2 = '2 25544  XX.XXXX 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      call tle_parse2(line2, ainc_a(1), raan_a(1), ecc_a(1),
     &                aop_a(1), ama_a(1), amm_a(1), revnum(1), ierr)
      call assert_i('parse2 corrupt inc -> ierr=1',
     &              ierr, 1, npass, nfail)

c     tle_parse2 with garbage in mean motion field -> ierr=1
      line2 = '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 XX.XXXXXXXX456789'
      call tle_parse2(line2, ainc_a(1), raan_a(1), ecc_a(1),
     &                aop_a(1), ama_a(1), amm_a(1), revnum(1), ierr)
      call assert_i('parse2 corrupt mm -> ierr=1',
     &              ierr, 1, npass, nfail)

c     tle_parse2 with valid data -> ierr=0
      line2 = '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      call tle_parse2(line2, ainc_a(1), raan_a(1), ecc_a(1),
     &                aop_a(1), ama_a(1), amm_a(1), revnum(1), ierr)
      call assert_i('parse2 valid -> ierr=0',
     &              ierr, 0, npass, nfail)

c     read_tle skips entries with corrupt fields
      open(unit=98, status='scratch')
c     Corrupt entry (garbage epoch)
      write(98,'(A)') '1 99999U 00000A   XX001.00000000  .000000'
     &     // '00  00000-0  00000+0 0  9999'
      write(98,'(A)') '2 99999  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
c     Valid entry
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('corrupt-entry skipped ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_i('corrupt-entry valid NORAD = 25544',
     &              norad(1), 25544, npass, nfail)

c     ============================================================
c     BOM stripping: UTF-8 BOM at file start
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- BOM stripping ---'

c     Write a file with a 3-byte UTF-8 BOM (EF BB BF) prepended
c     to a valid two-line TLE.  The reader should strip the BOM
c     and parse the entry correctly.
      open(unit=98, status='scratch', form='unformatted',
     &     access='stream')
      write(98) char(239), char(187), char(191)
      write(98) '1 25544U 98067A   26170.50000000  .00016717'
      write(98) '  00000-0  10270-3 0  9993'
      write(98) char(13), char(10)
      write(98) '2 25544  51.6400 123.4567 0006700  45.1234'
      write(98) ' 315.0000 15.54000000456789'
      write(98) char(13), char(10)
      close(98)

c     Re-open as formatted text for read_tle
      open(unit=98, file='fort.98', status='old', err=710)
      goto 711
 710  open(unit=98, status='scratch')
      write(98,'(A,A)') char(239)//char(187)//char(191),
     &     '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
 711  continue
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('BOM-file ntle = 1', ntle, 1, npass, nfail)
      call assert_i('BOM-file NORAD = 25544',
     &              norad(1), 25544, npass, nfail)

c     File without BOM still works (regression)
      open(unit=98, status='scratch')
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000456789'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('no-BOM regression ntle = 1',
     &              ntle, 1, npass, nfail)

c     ============================================================
c     Alpha-numeric revolution number (>99999)
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- alpha-numeric revnum ---'

c     revnum field = 'A0001' -> non-numeric, should return -1
c     but NOT set ierr (orbital data is still valid)
      line2 = '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000A00019'
      call tle_parse2(line2, ainc_a(1), raan_a(1), ecc_a(1),
     &                aop_a(1), ama_a(1), amm_a(1), revnum(1), ierr)
      call assert_i('alpha revnum = -1',
     &              revnum(1), -1, npass, nfail)
      call assert_i('alpha revnum ierr = 0',
     &              ierr, 0, npass, nfail)
      call assert_d('alpha revnum mm still parsed',
     &              amm_a(1), 15.54d0, 1.d-8, npass, nfail)

c     read_tle accepts entries with alpha revnum (not skipped)
      open(unit=98, status='scratch')
      write(98,'(A)') '1 25544U 98067A   26170.50000000  .000167'
     &     // '17  00000-0  10270-3 0  9993'
      write(98,'(A)') '2 25544  51.6400 123.4567 0006700  45.123'
     &     // '4 315.0000 15.54000000A00019'
      rewind(98)
      call read_tle(98, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(98)
      call assert_i('alpha-revnum entry accepted ntle = 1',
     &              ntle, 1, npass, nfail)
      call assert_i('alpha-revnum entry revnum = -1',
     &              revnum(1), -1, npass, nfail)

c     ============================================================
c     iostat guard: bad unit number doesn't crash
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- iostat guard ---'

c     Unit 79 is not open — should return ntle=0, not crash
      call read_tle(79, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      call assert_i('unopened unit ntle = 0',
     &              ntle, 0, npass, nfail)

c     ============================================================
c     example_47944.tle.txt: single-satellite multi-epoch
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- example_47944.tle.txt ---'

      open(unit=80, file='input/example_47944.tle.txt', status='old',
     &     err=840)

      call read_tle(80, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(80)

c     T76: entry count
      call assert_i('47944 ntle = 3977', ntle, 3977, npass, nfail)

c     T77: first NORAD
      call assert_i('47944 first NORAD = 47944',
     &              norad(1), 47944, npass, nfail)

c     T78: last NORAD
      call assert_i('47944 last NORAD = 47944',
     &              norad(ntle), 47944, npass, nfail)

c     T79: first epoch day
      call assert_d('47944 first eday = 83.05205046',
     &              epoch_day(1), 83.05205046d0, 1.d-6, npass, nfail)

c     T80: first entry negative ndot
      if (ndot_a(1) .lt. 0.d0) then
         call assert_i('47944 first ndot < 0', 1, 1, npass, nfail)
      else
         call assert_i('47944 first ndot < 0', 0, 1, npass, nfail)
      end if

c     T81: last entry has positive ndot (decaying orbit)
      if (ndot_a(ntle) .gt. 0.d0) then
         call assert_i('47944 last ndot > 0', 1, 1, npass, nfail)
      else
         call assert_i('47944 last ndot > 0', 0, 1, npass, nfail)
      end if

c     T82: all NORAD IDs = 47944
      iall_same = 1
      do i = 1, ntle
         if (norad(i) .ne. 47944) iall_same = 0
      end do
      call assert_i('47944 all NORAD = 47944',
     &              iall_same, 1, npass, nfail)

c     T83: all inclinations near 97 deg (SSO)
      iinc_ok = 1
      do i = 1, ntle
         if (ainc_a(i) .lt. 96.d0 .or. ainc_a(i) .gt. 99.d0)
     &      iinc_ok = 0
      end do
      call assert_i('47944 all inc in [96,99] (SSO)',
     &              iinc_ok, 1, npass, nfail)

c     T84: all eccentricities in [0, 0.1] (near-circular)
      iecc_ok = 1
      do i = 1, ntle
         if (ecc_a(i) .lt. 0.d0 .or. ecc_a(i) .gt. 0.1d0)
     &      iecc_ok = 0
      end do
      call assert_i('47944 all ecc in [0, 0.1]',
     &              iecc_ok, 1, npass, nfail)

c     T85: all mean motions in [14, 17] rev/day (LEO)
      imm_ok = 1
      do i = 1, ntle
         if (amm_a(i) .lt. 14.d0 .or. amm_a(i) .gt. 17.d0)
     &      imm_ok = 0
      end do
      call assert_i('47944 all mm in [14,17] (LEO)',
     &              imm_ok, 1, npass, nfail)

c     T86: epoch monotonically non-decreasing (same year epochs)
c     Check that combined (year, day) never decreases
      imon = 1
      do i = 2, ntle
         ep_prev = (epoch_yr(i-1)-2000)*1000.d0 + epoch_day(i-1)
         ep_curr = (epoch_yr(i)  -2000)*1000.d0 + epoch_day(i)
         if (ep_curr .lt. ep_prev - 0.001d0) imon = 0
      end do
      call assert_i('47944 epochs non-decreasing',
     &              imon, 1, npass, nfail)

c     T87: last entry epoch year = 2025 (decayed)
      call assert_i('47944 last epoch_yr = 2025',
     &              epoch_yr(ntle), 2025, npass, nfail)

c     T88: checksum every line 1 in file
      rewind_count = 0
      open(unit=80, file='input/example_47944.tle.txt', status='old')
      ichk_ok = 1
      ichk_cnt = 0
 810  continue
         read(80, '(A)', end=820) rawline
         if (rawline(1:1) .eq. '1' .or. rawline(1:1) .eq. '2') then
            call tle_chksum(rawline, ick, iok)
            ichk_cnt = ichk_cnt + 1
            if (iok .ne. 1) ichk_ok = 0
         end if
      goto 810
 820  close(80)
      call assert_i('47944 all checksums valid ('
     &     //'7954 lines)', ichk_ok, 1, npass, nfail)

      goto 850
 840  write(*,*) '  [SKIP]  example_47944.tle.txt not found'
 850  continue

c     ============================================================
c     example_multi.tle.txt: large multi-satellite catalog
c     ============================================================
      write(*,*) ''
      write(*,*) '  --- example_multi.tle.txt ---'

      open(unit=81, file='input/example_multi.tle.txt', status='old',
     &     err=960)

      call read_tle(81, maxtle, names, catid, norad,
     &              epoch_yr, epoch_day,
     &              ndot_a, nddot_a, bstar_a,
     &              ainc_a, raan_a, ecc_a, aop_a, ama_a, amm_a,
     &              revnum, ntle)
      close(81)

c     T89: entry count
      call assert_i('multi ntle = 94597', ntle, 94597, npass, nfail)

c     T90: first NORAD (Vanguard 1)
      call assert_i('multi first NORAD = 5', norad(1), 5, npass, nfail)

c     T91: first catid
      call assert_s('multi first catid = 00005',
     &              catid(1), '00005', npass, nfail)

c     T92: first inclination
      call assert_d('multi first inc = 34.2451',
     &              ainc_a(1), 34.2451d0, 1.d-4, npass, nfail)

c     T93: first eccentricity
      call assert_d('multi first ecc = 0.1835576',
     &              ecc_a(1), 0.1835576d0, 1.d-7, npass, nfail)

c     T94: first mean motion
      call assert_d('multi first mm = 10.86000951',
     &              amm_a(1), 10.86000951d0, 1.d-8, npass, nfail)

c     T95: first bstar
      call assert_d('multi first bstar = 3.7765e-4',
     &              bstar_a(1), 3.7765d-4, tol, npass, nfail)

c     T96: last NORAD is alpha-numeric
      call assert_i('multi last NORAD = -1 (alpha)',
     &              norad(ntle), -1, npass, nfail)

c     T97: last catid starts with T
      call assert_s('multi last catid(1:1) = T',
     &              catid(ntle)(1:1), 'T', npass, nfail)

c     T98: all mean motions positive
      imm_ok = 1
      do i = 1, ntle
         if (amm_a(i) .le. 0.d0) imm_ok = 0
      end do
      call assert_i('multi all mm > 0', imm_ok, 1, npass, nfail)

c     T99: all eccentricities in [0, 1)
      iecc_ok = 1
      do i = 1, ntle
         if (ecc_a(i) .lt. 0.d0 .or. ecc_a(i) .ge. 1.d0) iecc_ok = 0
      end do
      call assert_i('multi all 0 <= ecc < 1',
     &              iecc_ok, 1, npass, nfail)

c     T100: all inclinations in [0, 180]
      iinc_ok = 1
      do i = 1, ntle
         if (ainc_a(i) .lt. 0.d0 .or. ainc_a(i) .gt. 180.d0)
     &      iinc_ok = 0
      end do
      call assert_i('multi all inc in [0,180]',
     &              iinc_ok, 1, npass, nfail)

c     T101: all epoch years in [1957, 2057]
      iyr_ok = 1
      do i = 1, ntle
         if (epoch_yr(i) .lt. 1957 .or. epoch_yr(i) .gt. 2057)
     &      iyr_ok = 0
      end do
      call assert_i('multi all yr in [1957,2057]',
     &              iyr_ok, 1, npass, nfail)

c     T102: all epoch days in [1, 367]
      iday_ok = 1
      do i = 1, ntle
         if (epoch_day(i) .lt. 1.d0 .or. epoch_day(i) .gt. 367.d0)
     &      iday_ok = 0
      end do
      call assert_i('multi all eday in [1,367]',
     &              iday_ok, 1, npass, nfail)

c     T103: all mean motions correspond to reasonable periods
c     Period = 1440/mm minutes.  Bound: 80 min < P < 10 days
      iper_ok = 1
      do i = 1, ntle
         period = 1440.d0 / amm_a(i)
         if (period .lt. 80.d0 .or. period .gt. 14400.d0) iper_ok = 0
      end do
      call assert_i('multi all period in [80min,10d]',
     &              iper_ok, 1, npass, nfail)

c     T104: all RAAN in [0, 360]
      iraan_ok = 1
      do i = 1, ntle
         if (raan_a(i) .lt. 0.d0 .or. raan_a(i) .gt. 360.d0)
     &      iraan_ok = 0
      end do
      call assert_i('multi all raan in [0,360]',
     &              iraan_ok, 1, npass, nfail)

c     T105: all AOP in [0, 360]
      iaop_ok = 1
      do i = 1, ntle
         if (aop_a(i) .lt. 0.d0 .or. aop_a(i) .gt. 360.d0)
     &      iaop_ok = 0
      end do
      call assert_i('multi all aop in [0,360]',
     &              iaop_ok, 1, npass, nfail)

c     T106: all MA in [0, 360]
      ima_ok = 1
      do i = 1, ntle
         if (ama_a(i) .lt. 0.d0 .or. ama_a(i) .gt. 360.d0)
     &      ima_ok = 0
      end do
      call assert_i('multi all ma in [0,360]',
     &              ima_ok, 1, npass, nfail)

c     T107: count distinct numeric NORAD IDs > 1000
      ndistinct = 0
      iprev = -99
      do i = 1, ntle
         if (norad(i) .gt. 0 .and. norad(i) .ne. iprev) then
            ndistinct = ndistinct + 1
            iprev = norad(i)
         end if
      end do
      if (ndistinct .gt. 1000) then
         call assert_i('multi >1000 distinct NOR IDs',
     &                 1, 1, npass, nfail)
      else
         call assert_i('multi >1000 distinct NOR IDs',
     &                 0, 1, npass, nfail)
      end if

c     T108: count alpha-numeric entries > 0
      nalpha = 0
      do i = 1, ntle
         if (norad(i) .eq. -1) nalpha = nalpha + 1
      end do
      if (nalpha .gt. 0) then
         call assert_i('multi has alpha-numeric IDs',
     &                 1, 1, npass, nfail)
      else
         call assert_i('multi has alpha-numeric IDs',
     &                 0, 1, npass, nfail)
      end if

c     T109: checksum every line in multi file
      open(unit=81, file='input/example_multi.tle.txt', status='old')
      ichk_ok = 1
      ichk_cnt = 0
 910  continue
         read(81, '(A)', end=920) rawline
         if (rawline(1:1) .eq. '1' .or. rawline(1:1) .eq. '2') then
            call tle_chksum(rawline, ick, iok)
            ichk_cnt = ichk_cnt + 1
            if (iok .ne. 1) ichk_ok = 0
         end if
      goto 910
 920  close(81)
      call assert_i('multi all checksums valid',
     &              ichk_ok, 1, npass, nfail)

      goto 970
 960  write(*,*) '  [SKIP]  example_multi.tle.txt not found'
 970  continue

c     ============================================================
c     Summary
c     ============================================================
      write(*,*) ''
      write(*,*) '=============================================='
      write(*,100) npass+nfail, npass, nfail
 100  format(' Total: ',I4,'   Passed: ',I4,'   Failed: ',I4)
      write(*,*) '=============================================='

      stop
      end

c     ============================================================
c     Assertion helpers
c     ============================================================
      subroutine assert_d(label, got, expect, tol, npass, nfail)
      implicit double precision (a-h, o-z)
      character*(*) label
      double precision got, expect, tol, aerr
      integer npass, nfail

      aerr = dabs(got - expect)
      if (aerr .le. tol) then
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

      subroutine assert_i(label, got, expect, npass, nfail)
      character*(*) label
      integer got, expect, npass, nfail
      if (got .eq. expect) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,I10,A,I10)') '           got=',
     &        got, '  expect=', expect
         nfail = nfail + 1
      end if
      return
      end

      subroutine assert_s(label, got, expect, npass, nfail)
      character*(*) label, got, expect
      integer npass, nfail
      if (got .eq. expect) then
         write(*,'(A,A)') '  [PASS]  ', label
         npass = npass + 1
      else
         write(*,'(A,A)') '  [FAIL]  ', label
         write(*,'(A,A,A,A,A)') '           got=[',
     &        got, ']  expect=[', expect, ']'
         nfail = nfail + 1
      end if
      return
      end
