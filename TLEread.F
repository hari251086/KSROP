c     ============================================================
c     TLEread.F — Two-Line Element set reader
c
c     Reads standard NORAD/USSPACECOM TLE files in either two-line
c     (lines 1+2 only) or three-line (line 0 name + lines 1+2)
c     format.  Single and multi-entry files are handled.
c
c     Line-type detection validates column 1 digit, column 2 space,
c     minimum line length (69 chars), and — for line 1 — column 8
c     classification letter (U/C/S/space), so that satellite names
c     starting with '1' or '2' are not misidentified as TLE data.
c     Paired line 1/line 2 catalog IDs (cols 3-7) are cross-checked;
c     mismatched pairs are skipped.
c
c     All field reads in tle_parse1/tle_parse2 are guarded with err=
c     branches: a corrupt field zeros the output and sets an error
c     flag (ierr=1) rather than crashing the program.
c
c     Public subroutines
c     ------------------
c       read_tle     Main entry point
c
c     Internal helpers
c     ----------------
c       tle_parse1   Parse a TLE line 1 string
c       tle_parse2   Parse a TLE line 2 string
c       tle_expval   Decode implied-decimal exponent format
c       tle_chksum   Verify modulo-10 checksum
c       tle_isline   Validate whether a raw string is a TLE line
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine tle_isline(raw, llen, lnum, istle)
c
c     Determine whether a raw input string is a valid TLE line.
c
c     Input
c       raw   : input string (any length)
c       llen  : len_trim(raw) — actual content length
c       lnum  : expected line number (1 or 2)
c
c     Output
c       istle : 1 if the string matches TLE line lnum, 0 otherwise
c
c     Checks:
c       - content length >= 69
c       - column 1 matches the expected digit ('1' or '2')
c       - column 2 is a space
c       - for line 1: column 8 is a classification letter (U/C/S)
c         or a space (unclassified)
c     ============================================================
      character*(*) raw
      integer llen, lnum, istle
      character*1 expected, c8

      istle = 0
      if (llen .lt. 69) return
      if (lnum .eq. 1) then
         expected = '1'
      else
         expected = '2'
      end if
      if (raw(1:1) .ne. expected) return
      if (raw(2:2) .ne. ' ') return

      if (lnum .eq. 1) then
         c8 = raw(8:8)
         if (c8 .ne. 'U' .and. c8 .ne. 'C' .and.
     &       c8 .ne. 'S' .and. c8 .ne. ' ') return
      end if

      istle = 1
      return
      end

c     ============================================================
      subroutine read_tle(iunit, maxtle, names, catid, norad,
     &                    epoch_yr, epoch_day,
     &                    ndot, nddot, bstar,
     &                    ainc, raan, ecc, aop, ama, amm,
     &                    revnum, ntle)
c
c     Read all TLE entries from an already-opened file unit.
c
c     Input
c       iunit    : Fortran unit number (file must be open)
c       maxtle   : dimension of output arrays
c
c     Output — parallel arrays, one element per TLE entry:
c       names(i)    : satellite name (24 chars; blank if 2-line)
c       catid(i)    : catalog ID string (5 chars, preserves alpha-
c                     numeric IDs like 'T0003')
c       norad(i)    : NORAD catalog number (-1 if alpha-numeric)
c       epoch_yr(i) : epoch year (4-digit, expanded from 2-digit)
c       epoch_day(i): epoch day of year + fractional day
c       ndot(i)     : 1st derivative of mean motion (rev/day^2)
c       nddot(i)    : 2nd derivative of mean motion (rev/day^3)
c       bstar(i)    : B* drag term (1/Earth radii)
c       ainc(i)     : inclination (deg)
c       raan(i)     : right ascension of asc. node (deg)
c       ecc(i)      : eccentricity (dimensionless)
c       aop(i)      : argument of perigee (deg)
c       ama(i)      : mean anomaly (deg)
c       amm(i)      : mean motion (rev/day)
c       revnum(i)   : revolution number at epoch
c       ntle        : number of entries actually read
c     ============================================================
      implicit double precision (a-h, o-z)
      integer iunit, maxtle, ntle
      integer norad(maxtle), epoch_yr(maxtle), revnum(maxtle)
      double precision epoch_day(maxtle)
      double precision ndot(maxtle), nddot(maxtle), bstar(maxtle)
      double precision ainc(maxtle), raan(maxtle), ecc(maxtle)
      double precision aop(maxtle), ama(maxtle), amm(maxtle)
      character*24 names(maxtle)
      character*5  catid(maxtle)

      character*80 rawbuf
      character*69 line1
      character*24 curname
      character*5  tmp_catid
      integer      tmp_norad, tmp_yr, tmp_rev
      double precision tmp_eday, tmp_ndot, tmp_nddot, tmp_bstar
      double precision tmp_inc, tmp_raan, tmp_ecc
      double precision tmp_aop, tmp_ma, tmp_mm
      integer      llen, is1, is2, ierr1, ierr2, ios
      logical      have1, first_line

      ntle       = 0
      have1      = .false.
      first_line = .true.
      curname    = ' '
      line1      = ' '

 10   continue
         read(iunit, '(A)', iostat=ios) rawbuf
         if (ios .ne. 0) goto 99
         llen = len_trim(rawbuf)

c        Strip UTF-8 BOM (EF BB BF) from the first line if present.
c        In Fortran's default encoding the three BOM bytes appear as
c        non-printable characters before the first real column.  Detect
c        by checking whether column 1 is NOT an ASCII printable char
c        (i.e. ichar < 32 or ichar > 126) while the line is long
c        enough to still contain a TLE after skipping up to 3 bytes.
         if (first_line) then
            first_line = .false.
            if (llen .ge. 72) then
               if (ichar(rawbuf(1:1)) .gt. 126 .or.
     &             ichar(rawbuf(1:1)) .lt. 32) then
c                 Shift left by 3 to skip the 3-byte BOM
                  rawbuf = rawbuf(4:)
                  llen   = len_trim(rawbuf)
               end if
            end if
         end if

         if (llen .eq. 0) goto 10

c        --- Detect line type ---
         call tle_isline(rawbuf, llen, 1, is1)
         call tle_isline(rawbuf, llen, 2, is2)

         if (is1 .eq. 1) then
c           Valid TLE line 1.  Store it and wait for line 2.
            line1 = rawbuf(1:69)
            have1 = .true.

         else if (is2 .eq. 1 .and. have1) then
c           Valid TLE line 2 paired with a stored line 1.
c           Cross-check catalog IDs (cols 3-7) before accepting.
            if (line1(3:7) .ne. rawbuf(3:7)) then
c              Catalog ID mismatch — skip this pair.
               have1 = .false.
               goto 10
            end if

c           Parse both lines into scratch variables; skip if corrupt.
            call tle_parse1(line1, tmp_catid, tmp_norad,
     &                      tmp_yr, tmp_eday,
     &                      tmp_ndot, tmp_nddot, tmp_bstar, ierr1)
            call tle_parse2(rawbuf(1:69), tmp_inc, tmp_raan,
     &                      tmp_ecc, tmp_aop, tmp_ma,
     &                      tmp_mm, tmp_rev, ierr2)

            if (ierr1 .ne. 0 .or. ierr2 .ne. 0) then
               have1   = .false.
               curname = ' '
               goto 10
            end if

            ntle = ntle + 1
            if (ntle .gt. maxtle) then
               ntle = maxtle
               goto 99
            end if

            names(ntle)     = curname
            catid(ntle)     = tmp_catid
            norad(ntle)     = tmp_norad
            epoch_yr(ntle)  = tmp_yr
            epoch_day(ntle) = tmp_eday
            ndot(ntle)      = tmp_ndot
            nddot(ntle)     = tmp_nddot
            bstar(ntle)     = tmp_bstar
            ainc(ntle)      = tmp_inc
            raan(ntle)      = tmp_raan
            ecc(ntle)       = tmp_ecc
            aop(ntle)       = tmp_aop
            ama(ntle)       = tmp_ma
            amm(ntle)       = tmp_mm
            revnum(ntle)    = tmp_rev

            have1   = .false.
            curname = ' '

         else
c           Not a valid TLE line — treat as line 0 (satellite name).
c           If we had an unpaired line 1, discard it.
            have1   = .false.
            curname = rawbuf(1:min(llen,24))
         end if

      goto 10

 99   continue
      return
      end

c     ============================================================
      subroutine tle_parse1(line, catid, norad, iyr4, eday, dndot,
     &                      dnddot, dbstar, ierr)
c
c     Parse a TLE line 1 string (69 characters).
c
c     Output
c       catid  : catalog ID string (5 chars, as-is from columns 3-7)
c       norad  : NORAD catalog number (-1 if alpha-numeric)
c       iyr4   : epoch year (4-digit)
c       eday   : epoch day of year + fractional day
c       dndot  : 1st deriv of mean motion / 2  (rev/day^2)
c       dnddot : 2nd deriv of mean motion / 6  (rev/day^3)
c       dbstar : B* drag term (1/Earth radii)
c       ierr   : 0 = success, 1 = parse error in a required field
c     ============================================================
      implicit double precision (a-h, o-z)
      character*69 line
      character*5  catid
      integer norad, iyr4, iyr2, ierr

      ierr = 0

c     NORAD catalog number (columns 3-7): may be purely numeric
c     (e.g. '47944') or alpha-numeric (e.g. 'T0003').  Return -1
c     for non-numeric IDs.
      catid = line(3:7)
      read(catid, '(I5)', err=15) norad
      goto 16
 15   norad = -1
 16   continue

c     Epoch year (columns 19-20)
      read(line(19:20), '(I2)', err=80) iyr2
c     Epoch day (columns 21-32)
      read(line(21:32), '(F12.8)', err=80) eday

c     Expand 2-digit year: 57-99 -> 1957-1999, 00-56 -> 2000-2056
      if (iyr2 .ge. 57) then
         iyr4 = 1900 + iyr2
      else
         iyr4 = 2000 + iyr2
      end if

c     1st derivative of mean motion (columns 34-43)
      read(line(34:43), '(F10.8)', err=80) dndot

c     2nd derivative of mean motion (columns 45-52, implied-decimal)
      call tle_expval(line(45:52), dnddot)

c     B* drag term (columns 54-61, implied-decimal)
      call tle_expval(line(54:61), dbstar)

      return

c     Parse error: zero remaining fields and flag
 80   iyr4   = 0
      eday   = 0.d0
      dndot  = 0.d0
      dnddot = 0.d0
      dbstar = 0.d0
      ierr   = 1
      return
      end

c     ============================================================
      subroutine tle_parse2(line, ainc, araan, aecc, aaop, aama,
     &                      aamm, irevnum, ierr)
c
c     Parse a TLE line 2 string (69 characters).
c
c     Output
c       ainc    : inclination (deg)
c       araan   : RAAN (deg)
c       aecc    : eccentricity (dimensionless, decimal point assumed)
c       aaop    : argument of perigee (deg)
c       aama    : mean anomaly (deg)
c       aamm    : mean motion (rev/day)
c       irevnum : revolution number at epoch
c       ierr    : 0 = success, 1 = parse error in a field
c     ============================================================
      implicit double precision (a-h, o-z)
      character*69 line
      integer irevnum, ierr

      ierr = 0

      read(line(9:16),  '(F8.4)', err=80) ainc
      read(line(18:25), '(F8.4)', err=80) araan

c     Eccentricity: columns 27-33, no decimal point in file
c     (e.g. '6543210' means 0.6543210)
      read(line(27:33), '(F7.7)', err=80) aecc

      read(line(35:42), '(F8.4)', err=80) aaop
      read(line(44:51), '(F8.4)', err=80) aama
      read(line(53:63), '(F11.8)', err=80) aamm

c     Revolution number (columns 64-68): may contain alpha characters
c     for revnum > 99999.  Return -1 for non-numeric values.
      read(line(64:68), '(I5)', err=75) irevnum
      goto 76
 75   irevnum = -1
 76   continue

      return

c     Parse error: zero remaining fields and flag
 80   ainc    = 0.d0
      araan   = 0.d0
      aecc    = 0.d0
      aaop    = 0.d0
      aama    = 0.d0
      aamm    = 0.d0
      irevnum = 0
      ierr    = 1
      return
      end

c     ============================================================
      subroutine tle_expval(str, val)
c
c     Decode the TLE implied-decimal exponent format.
c
c     Format (8 characters, right-justified):
c       positions 1-6: mantissa with implied leading '0.'
c                       (first char may be sign or space/digit)
c       position  7:   exponent sign (+/-)
c       position  8:   exponent digit (single)
c
c     Examples:
c       ' 12345-6'  ->  0.12345e-6
c       '-12345-6'  -> -0.12345e-6
c       ' 00000+0'  ->  0.0
c       '+12345+0'  ->  0.12345
c     ============================================================
      implicit double precision (a-h, o-z)
      character*8 str
      double precision val, dmant
      integer imant, iexp, isign

c     Handle all-blank or zero
      val = 0.d0

c     Sign from position 1
      isign = 1
      if (str(1:1) .eq. '-') isign = -1

c     Mantissa: always positions 2-6 (5 digits), implied '0.'
      read(str(2:6), '(I5)', err=90) imant
      dmant = dble(imant) * 1.d-5

c     Exponent: position 7 = sign, position 8 = digit
      iexp = 0
      if (str(8:8) .ge. '0' .and. str(8:8) .le. '9') then
         read(str(8:8), '(I1)') iexp
      end if
      if (str(7:7) .eq. '-') iexp = -iexp

      val = dble(isign) * dmant * (10.d0 ** iexp)
      return

 90   val = 0.d0
      return
      end

c     ============================================================
      subroutine tle_chksum(line, ick, iok)
c
c     Verify the modulo-10 checksum of a TLE line.
c
c     Input
c       line  : 69-character TLE line (line 1 or line 2)
c
c     Output
c       ick   : computed checksum (mod 10)
c       iok   : 1 if checksum matches column 69, 0 otherwise
c     ============================================================
      implicit double precision (a-h, o-z)
      character*69 line
      character*1  ch
      integer      ick, iok, isum, j, id

      isum = 0
      do j = 1, 68
         ch = line(j:j)
         if (ch .ge. '0' .and. ch .le. '9') then
            read(ch, '(I1)') id
            isum = isum + id
         else if (ch .eq. '-') then
            isum = isum + 1
         end if
      end do

      ick = mod(isum, 10)

      read(line(69:69), '(I1)', err=10) id
      if (ick .eq. id) then
         iok = 1
      else
         iok = 0
      end if
      return

 10   iok = 0
      return
      end
