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
      subroutine tle_epoch2jd(iyr4, eday, djd)
c
c     Convert TLE epoch (4-digit year + fractional day-of-year)
c     to Julian Date (UTC).
c
c     Input
c       iyr4  : 4-digit year
c       eday  : day of year + fractional day (1-based)
c
c     Output
c       djd   : Julian Date (UTC)
c     ============================================================
      implicit double precision (a-h, o-z)
      integer iyr4, iy

c     JD of Jan 0.0 UT of the given year, then add day-of-year.
c     Formula from Vallado, valid for years 1901-2099.
      iy = iyr4
      djd = 2415020.5d0 + 365.d0*(iy - 1900)
     &    + int(dble(iy-1901)/4.d0) + eday - 1.d0
      return
      end

c     ============================================================
      subroutine tle_sgp4_init(amm_revday, ecc0, ainc_deg,
     &                         raan_deg, aop_deg, ama_deg,
     &                         bstar, epoch_jd,
     &                         n0dd, a0dd, cosio, sinio,
     &                         eta, c1, c4, t0, perigee,
     &                         omgdot, xnodot, xmdot, ierr)
c
c     SGP4 initialization: recover original mean motion/SMA,
c     compute secular rates and drag coefficients.
c
c     Input (TLE-parsed values)
c       amm_revday : mean motion (rev/day)
c       ecc0       : eccentricity
c       ainc_deg   : inclination (deg)
c       raan_deg   : RAAN (deg)
c       aop_deg    : argument of perigee (deg)
c       ama_deg    : mean anomaly (deg)
c       bstar      : B* drag term (1/Earth radii)
c       epoch_jd   : epoch Julian Date
c
c     Output
c       n0dd   : recovered mean motion (rad/min)
c       a0dd   : recovered semi-major axis (Earth radii)
c       cosio, sinio : trig of inclination
c       eta    : eccentricity parameter
c       c1     : drag coefficient
c       c4     : drag coefficient for perigee
c       t0     : epoch time (JD, stored for reference)
c       perigee: perigee altitude (km)
c       omgdot : secular rate of arg of perigee (rad/min)
c       xnodot : secular rate of RAAN (rad/min)
c       xmdot  : secular rate of mean anomaly (rad/min)
c       ierr   : 0=ok, 1=deep-space (unsupported), 2=bad elements
c     ============================================================
      implicit double precision (a-h, o-z)
      integer ierr

c     Explicit declarations for i-n range variables
      double precision n0dd, n0

c     WGS-72 constants (canonical for SGP4)
      double precision ke, j2, j3, j4, re_km
      parameter (re_km = 6378.135d0)
      parameter (ke    = 0.0743669161d0)
      parameter (j2    = 1.0826158d-3)
      parameter (j3    = -2.53881d-6)
      parameter (j4    = -1.65597d-6)

      double precision k2, k4
      double precision del1, del0
      double precision theta2, theta4, beta0, beta02
      double precision qoms24
      double precision coef, coef1, eeta
      double precision psisq, tsi

      ierr = 0
      pi    = 3.14159265358979324d0
      twopi = 2.d0 * pi
      d2r   = pi / 180.d0

c     Derived gravity constants
      k2  = 0.5d0 * j2 * re_km**2
      k4  = -0.375d0 * j4 * re_km**4


c     Convert inputs to SGP4 internal units
      ainc = ainc_deg * d2r
      cosio = dcos(ainc)
      sinio = dsin(ainc)
      theta2 = cosio * cosio
      theta4 = theta2 * theta2

c     Mean motion: rev/day -> rad/min
      n0 = amm_revday * twopi / 1440.d0

c     Recover original mean motion and semi-major axis
c     a is in units of Earth radii
      a1 = (ke / n0) ** (2.d0/3.d0)
      del1 = 1.5d0*(k2/(re_km**2))/(a1*a1)
     &      * (3.d0*theta2-1.d0)
     &      / (1.d0-ecc0*ecc0)**1.5d0
      a0 = a1 * (1.d0 - del1/3.d0 - del1**2
     &         - 134.d0/81.d0 * del1**3)
      del0 = 1.5d0*(k2/(re_km**2))/(a0*a0)
     &      * (3.d0*theta2-1.d0)
     &      / (1.d0-ecc0*ecc0)**1.5d0

      n0dd = n0 / (1.d0 + del0)
      a0dd = a0 / (1.d0 - del0)

c     Orbital period check: near-Earth vs deep-space
      period_min = twopi / n0dd
      if (period_min .ge. 225.d0) then
         ierr = 1
         return
      end if

c     Perigee altitude (km)
      perigee = (a0dd * (1.d0 - ecc0) - 1.d0) * re_km

c     Atmospheric parameter s
      s_param = 78.d0 / re_km + 1.d0
      if (perigee .lt. 156.d0) then
         if (perigee .lt. 98.d0) then
            s_param = 20.d0 / re_km + 1.d0
         else
            s_param = (perigee - 78.d0) / re_km + 1.d0
         end if
      end if

c     qoms2t = (q0 - s)^4 with q0 = 120 km
      qoms24 = ((120.d0 - (s_param - 1.d0)*re_km) / re_km) ** 4
      if (qoms24 .lt. 1.d-24) qoms24 = 1.d-24

      tsi = 1.d0 / (a0dd - s_param)
      eta = a0dd * ecc0 * tsi
      eeta = ecc0 * eta

      beta0 = dsqrt(1.d0 - ecc0*ecc0)
      beta02 = beta0 * beta0

      psisq = dabs(1.d0 - eta*eta)
      if (psisq .lt. 1.d-20) psisq = 1.d-20

      coef  = qoms24 * tsi**4
      coef1 = coef / psisq**3.5d0

c     C2 coefficient
      c2 = coef1 * n0dd * (a0dd * (1.d0 + 1.5d0*eta*eta + eeta*
     &     (4.d0 + eta*eta)) + 0.75d0 * (k2/(re_km**2)) * tsi /
     &     psisq * (-0.5d0 + 1.5d0*theta2) *
     &     (8.d0 + 3.d0*eta*eta*(8.d0 + eta*eta)))

c     C1 coefficient
      c1 = bstar * c2

c     C3 coefficient (for drag perigee term)
      if (ecc0 .gt. 1.d-4) then
         c3 = -2.d0 * coef * tsi * (j3/j2) * n0dd *
     &        sinio / ecc0
      else
         c3 = 0.d0
      end if

c     C4 coefficient
      c4 = 2.d0*n0dd*coef1*a0dd*beta02*
     &     (eta*(2.d0+0.5d0*eta*eta) +
     &      ecc0*(0.5d0+2.d0*eta*eta) -
     &      2.d0*(k2/(re_km**2))*tsi/(a0dd*psisq)*
     &      (-3.d0*(1.d0-3.d0*theta2)*(1.d0+1.5d0*eta*eta
     &       -2.d0*eeta*(1.d0+eta*eta)) +
     &       0.75d0*(1.d0-theta2)*(2.d0*eta*eta
     &       -eeta*(1.d0+eta*eta))))

c     Secular rates due to J2, J4 (Vallado convention: j2, j4 directly)
      pinvsq = 1.d0 / (a0dd*a0dd * beta02*beta02)
      temp1 = 1.5d0 * j2 * pinvsq * n0dd
      temp2 = 0.5d0 * temp1 * j2 * pinvsq
      temp3 = -0.46875d0 * j4 * pinvsq * pinvsq * n0dd

      xmdot = n0dd + 0.5d0*temp1*beta0*(3.d0*theta2-1.d0)
     &      + 0.0625d0*temp2*beta0*
     &        (13.d0-78.d0*theta2+137.d0*theta4)

      omgdot = -0.5d0*temp1*(1.d0-5.d0*theta2)
     &       + 0.0625d0*temp2*(7.d0-114.d0*theta2+395.d0*theta4)
     &       + temp3*(3.d0-36.d0*theta2+49.d0*theta4)

      xnodot = -(temp1*cosio)
     &       + (0.5d0*temp2*(4.d0-19.d0*theta2)
     &       + 2.d0*temp3*(3.d0-7.d0*theta2)) * cosio

      t0 = epoch_jd

      return
      end

c     ============================================================
      subroutine tle_sgp4_prop(tsince,
     &              ecc0, ainc_deg, raan_deg, aop_deg, ama_deg,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1, c4, perigee,
     &              omgdot, xnodot, xmdot, bstar,
     &              r_teme, v_teme, ierr)
c
c     SGP4 propagation at tsince minutes from epoch.
c     For TLE-to-state-vector at epoch, call with tsince = 0.
c
c     Output
c       r_teme(3) : position in TEME frame (km)
c       v_teme(3) : velocity in TEME frame (km/s)
c       ierr      : 0=ok, nonzero=error
c     ============================================================
      implicit double precision (a-h, o-z)
      integer ierr, i
      double precision r_teme(3), v_teme(3)

c     Explicit declarations for i-n range arguments/locals
      double precision n0dd, k2, mrt

      double precision j2, j3, ke, re_km
      parameter (re_km = 6378.135d0)
      parameter (ke    = 0.0743669161d0)
      parameter (j2    = 1.0826158d-3)
      parameter (j3    = -2.53881d-6)

      double precision mp, omega, xnode, ecc_t, a_t
      double precision axn, ayn, aynl
      double precision xl, u, eo1, tem5, sineo1, coseo1
      double precision ecose, esine, el2, pl, rl, rdotl, rvdotl
      double precision betal, sinu, cosu, su, sin2u, cos2u
      double precision rk, uk, xnodek, xinck
      double precision sinuk, cosuk, sinik, cosik, sinnok, cosnok
      double precision xmx, xmy, ux, uy, uz, vx, vy, vz

      ierr = 0
      pi    = 3.14159265358979324d0
      twopi = 2.d0 * pi
      d2r   = pi / 180.d0

      k2 = 0.5d0 * j2 * re_km**2


c     Convert initial angles to radians
      omega0 = aop_deg * d2r
      xnode0 = raan_deg * d2r
      xm0    = ama_deg * d2r

c     Secular updates
      mp    = xm0 + xmdot * tsince
      omega = omega0 + omgdot * tsince
      xnode = xnode0 + xnodot * tsince

c     Drag secular effects on a and e
      tempa = 1.d0 - c1 * tsince
      tempe = bstar * c4 * tsince
      a_t   = a0dd * tempa * tempa
      ecc_t = ecc0 - tempe
      if (ecc_t .lt. 1.d-6) ecc_t = 1.d-6
      if (ecc_t .gt. 0.999d0) ecc_t = 0.999d0

c     Long-period periodics
      axn  = ecc_t * dcos(omega)
      pl   = a_t * (1.d0 - ecc_t*ecc_t)
      if (pl .lt. 1.d-12) pl = 1.d-12
      aynl = -0.5d0 * (j3/j2) * sinio / pl
      ayn  = ecc_t * dsin(omega) + aynl

c     Solve Kepler's equation for (E + omega)
      xl = mp + omega + xnode
      u  = dmod(xl - xnode, twopi)
      if (u .lt. 0.d0) u = u + twopi


      eo1 = u
      do i = 1, 10
         sineo1 = dsin(eo1)
         coseo1 = dcos(eo1)
         tem5 = 1.d0 - coseo1*axn - sineo1*ayn
         if (dabs(tem5) .lt. 1.d-20) tem5 = 1.d-20
         tem5 = (u - ayn*coseo1 + axn*sineo1 - eo1) / tem5
         if (dabs(tem5) .gt. 0.95d0) then
            if (tem5 .gt. 0.d0) then
               tem5 = 0.95d0
            else
               tem5 = -0.95d0
            end if
         end if
         eo1 = eo1 + tem5
      end do

c     Short-period preliminary quantities
      sineo1 = dsin(eo1)
      coseo1 = dcos(eo1)
      ecose = axn*coseo1 + ayn*sineo1
      esine = axn*sineo1 - ayn*coseo1
      el2   = axn*axn + ayn*ayn
      pl    = a_t * (1.d0 - el2)
      if (pl .lt. 0.d0) then
         ierr = 2
         return
      end if
      rl    = a_t * (1.d0 - ecose)
      rdotl = ke * dsqrt(a_t) * esine / rl
      rvdotl = ke * dsqrt(pl) / rl

      betal = dsqrt(1.d0 - el2)
      sinu  = a_t / rl * (sineo1 - ayn - axn*esine/(1.d0+betal))
      cosu  = a_t / rl * (coseo1 - axn + ayn*esine/(1.d0+betal))
      su    = datan2(sinu, cosu)
      sin2u = 2.d0 * sinu * cosu
      cos2u = 2.d0 * cosu * cosu - 1.d0

      theta2 = cosio * cosio

c     Short-period J2 corrections
      rk = rl*(1.d0 - 1.5d0*(k2/(re_km**2))*dsqrt(1.d0-el2)/
     &     (pl)*(3.d0*theta2-1.d0)) +
     &     0.5d0*(k2/(re_km**2))*(1.d0-theta2)/pl * cos2u
      uk = su - 0.25d0*(k2/(re_km**2))/(pl*pl)*
     &     (7.d0*theta2-1.d0) * sin2u
      xnodek = xnode + 1.5d0*(k2/(re_km**2))*cosio/(pl*pl)*sin2u
      xinck = ainc_deg*d2r + 1.5d0*(k2/(re_km**2))*cosio*sinio/
     &     (pl*pl)*cos2u

c     Unit orientation vectors
      sinuk = dsin(uk)
      cosuk = dcos(uk)
      sinik = dsin(xinck)
      cosik = dcos(xinck)
      sinnok = dsin(xnodek)
      cosnok = dcos(xnodek)

      xmx = -sinnok * cosik
      xmy =  cosnok * cosik

      ux = xmx*sinuk + cosnok*cosuk
      uy = xmy*sinuk + sinnok*cosuk
      uz = sinik*sinuk

      vx = xmx*cosuk - cosnok*sinuk
      vy = xmy*cosuk - sinnok*sinuk
      vz = sinik*cosuk

c     Position (km) and velocity (km/s)
      mrt = rk * re_km
      if (mrt .lt. 1.d0) then
         ierr = 3
         return
      end if

      r_teme(1) = mrt * ux
      r_teme(2) = mrt * uy
      r_teme(3) = mrt * uz

      rdotk  = rdotl * re_km / 60.d0
      rfdotk = rvdotl * re_km / 60.d0

      v_teme(1) = rdotk*ux + rfdotk*vx
      v_teme(2) = rdotk*uy + rfdotk*vy
      v_teme(3) = rdotk*uz + rfdotk*vz

      return
      end

c     ============================================================
      subroutine tle_teme2j2k(r_teme, v_teme, djd,
     &                        r_j2k, v_j2k)
c
c     Transform TEME state vector to J2000 (EME2000) frame.
c     Uses IAU-76/FK5 precession + simplified IAU-80 nutation.
c
c     Input
c       r_teme(3) : position in TEME (km)
c       v_teme(3) : velocity in TEME (km/s)
c       djd       : Julian Date (UTC, approximated as TT)
c
c     Output
c       r_j2k(3)  : position in J2000 (km)
c       v_j2k(3)  : velocity in J2000 (km/s)
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision r_teme(3), v_teme(3), r_j2k(3), v_j2k(3)
      double precision r_mod(3), v_mod(3)

      double precision pi, d2r, as2r
      double precision ttt, omega_moon
      double precision dpsi, deps, eps0, eps_true, eqe
      double precision zeta, theta, z_prec
      double precision pmat(3,3), nmat(3,3), eqmat(3,3)
      double precision tmp(3,3), rot(3,3)
      integer i, j, kk

      pi  = 3.14159265358979324d0
      d2r = pi / 180.d0
      as2r = d2r / 3600.d0

c     Julian centuries of TT since J2000.0
c     (Approximate: ignore UTC-TT offset ~69s, negligible for SGP4)
      ttt = (djd - 2451545.d0) / 36525.d0

c     IAU-1980 nutation (4-term approximation)
c     Lunar node longitude (deg)
      omega_moon = 125.04452d0 - 1934.136261d0*ttt
      omega_moon = dmod(omega_moon, 360.d0) * d2r

c     Nutation in longitude (arcsec -> rad)
      dpsi = (-17.20d0*dsin(omega_moon)
     &        - 1.32d0*dsin(2.d0*omega_moon)) * as2r

c     Nutation in obliquity (arcsec -> rad)
      deps = (9.20d0*dcos(omega_moon)
     &       + 0.57d0*dcos(2.d0*omega_moon)) * as2r

c     Mean obliquity of date (arcsec -> rad)
      eps0 = (23.439291d0 - 0.0130042d0*ttt
     &       - 1.64d-7*ttt*ttt + 5.04d-7*ttt*ttt*ttt) * d2r

      eps_true = eps0 + deps

c     Equation of the equinoxes
      eqe = dpsi * dcos(eps0)
     &    + 0.00264d0*as2r*dsin(omega_moon)
     &    + 0.000063d0*as2r*dsin(2.d0*omega_moon)

c     Equation of equinoxes rotation: R_z(-eqe)
      eqmat(1,1) =  dcos(eqe)
      eqmat(1,2) =  dsin(eqe)
      eqmat(1,3) =  0.d0
      eqmat(2,1) = -dsin(eqe)
      eqmat(2,2) =  dcos(eqe)
      eqmat(2,3) =  0.d0
      eqmat(3,1) =  0.d0
      eqmat(3,2) =  0.d0
      eqmat(3,3) =  1.d0

c     Nutation matrix: NUT = Rx(-eps0) Rz(dpsi) Rx(eps_true)
c     Build the full nutation matrix directly:
c     NUT = Rx(-eps0) * Rz(dpsi) * Rx(eps_true)
c     NUT^T applied to TEME -> MOD
c     For simplicity, build column by column:
      se0 = dsin(eps0)
      ce0 = dcos(eps0)
      set = dsin(eps_true)
      cet = dcos(eps_true)
      sdp = dsin(dpsi)
      cdp = dcos(dpsi)

c     NUT = Rx(-eps0) Rz(dpsi) Rx(eps_true)
c     Row 1: [cdp, sdp*ce0, sdp*se0] -- wait, let me do this properly
c     Rx(-eps0) = [[1, 0, 0], [0, ce0, se0], [0, -se0, ce0]]
c     Rz(dpsi)  = [[cdp, -sdp, 0], [sdp, cdp, 0], [0, 0, 1]]
c     Rx(eps_t) = [[1, 0, 0], [0, cet, -set], [0, set, cet]]
c
c     Product A = Rz(dpsi)*Rx(eps_t):
c       A(1,1)=cdp  A(1,2)=-sdp*cet  A(1,3)=sdp*set
c       A(2,1)=sdp  A(2,2)=cdp*cet   A(2,3)=-cdp*set
c       A(3,1)=0    A(3,2)=set        A(3,3)=cet
c
c     NUT = Rx(-eps0)*A:
c       NUT(1,:) = A(1,:)
c       NUT(2,:) = ce0*A(2,:) + se0*A(3,:)
c       NUT(3,:) = -se0*A(2,:) + ce0*A(3,:)

      nmat(1,1) = cdp
      nmat(1,2) = -sdp*cet
      nmat(1,3) = sdp*set
      nmat(2,1) = ce0*sdp
      nmat(2,2) = ce0*cdp*cet + se0*set
      nmat(2,3) = -ce0*cdp*set + se0*cet
      nmat(3,1) = -se0*sdp
      nmat(3,2) = -se0*cdp*cet + ce0*set
      nmat(3,3) = se0*cdp*set + ce0*cet

c     TEME -> MOD: r_mod = NUT^T * Rz(-eqe) * r_teme
c     First: r1 = eqmat * r_teme
      do i = 1, 3
         r_mod(i) = 0.d0
         v_mod(i) = 0.d0
         do j = 1, 3
            r_mod(i) = r_mod(i) + eqmat(i,j)*r_teme(j)
            v_mod(i) = v_mod(i) + eqmat(i,j)*v_teme(j)
         end do
      end do

c     Then: r_mod2 = NUT^T * r_mod
      do i = 1, 3
         r_j2k(i) = 0.d0
         v_j2k(i) = 0.d0
         do j = 1, 3
            r_j2k(i) = r_j2k(i) + nmat(j,i)*r_mod(j)
            v_j2k(i) = v_j2k(i) + nmat(j,i)*v_mod(j)
         end do
      end do

c     Save MOD results temporarily
      do i = 1, 3
         r_mod(i) = r_j2k(i)
         v_mod(i) = v_j2k(i)
      end do

c     Precession: MOD -> J2000
c     IAU-76 precession angles (arcsec)
      zeta   = (2306.2181d0*ttt + 0.30188d0*ttt*ttt
     &        + 0.017998d0*ttt*ttt*ttt) * as2r
      theta  = (2004.3109d0*ttt - 0.42665d0*ttt*ttt
     &        - 0.041833d0*ttt*ttt*ttt) * as2r
      z_prec = (2306.2181d0*ttt + 1.09468d0*ttt*ttt
     &        + 0.018203d0*ttt*ttt*ttt) * as2r

c     PREC = Rz(-z) Ry(theta) Rz(-zeta)
c     r_j2k = PREC^T * r_mod
c     PREC^T = Rz(zeta) Ry(-theta) Rz(z)

      sz = dsin(z_prec)
      cz = dcos(z_prec)
      st = dsin(theta)
      ct = dcos(theta)
      szt = dsin(zeta)
      czt = dcos(zeta)

c     Build PREC matrix
      pmat(1,1) =  czt*ct*cz - szt*sz
      pmat(1,2) = -szt*ct*cz - czt*sz
      pmat(1,3) = -st*cz
      pmat(2,1) =  czt*ct*sz + szt*cz
      pmat(2,2) = -szt*ct*sz + czt*cz
      pmat(2,3) = -st*sz
      pmat(3,1) =  czt*st
      pmat(3,2) = -szt*st
      pmat(3,3) =  ct

c     r_j2k = PREC^T * r_mod
      do i = 1, 3
         r_j2k(i) = 0.d0
         v_j2k(i) = 0.d0
         do j = 1, 3
            r_j2k(i) = r_j2k(i) + pmat(j,i)*r_mod(j)
            v_j2k(i) = v_j2k(i) + pmat(j,i)*v_mod(j)
         end do
      end do

      return
      end

c     ============================================================
      subroutine tle_ds_common(epoch_jd, ecc0, ainc, raan, aop,
     &                         n0dd, dsstate)
c
c     SDP4 deep-space common: compute solar/lunar geometry at epoch
c     and the 24 periodic perturbation coefficients.
c
c     Based on Vallado's dscom() from sgp4unit.cpp
c
c     Input
c       epoch_jd : Julian Date of TLE epoch
c       ecc0     : eccentricity
c       ainc     : inclination (rad)
c       raan     : RAAN (rad)
c       aop      : argument of perigee (rad)
c       n0dd     : recovered mean motion (rad/min)
c
c     Output
c       dsstate(1:26) : zmol, zmos, solar coefficients (3-14),
c                       lunar coefficients (15-26)
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision dsstate(50)
      double precision n0dd

      double precision pi, twopi, d2r
      double precision tut1, day
      double precision zes, zel, zns, znl
      double precision ses, sis, sls, sghs, shs
      double precision c1ss, c1l
      double precision a1s,a2s,a3s,a4s,a5s,a6s,a7s,a8s,a9s,a10s
      double precision betasq, cc, ctem, stem, xnoi

      pi    = 3.14159265358979324d0
      twopi = 2.d0 * pi
      d2r   = pi / 180.d0

c     Solar/lunar perturbation constants from Spacetrack Report #3
      zes = 0.01675d0
      zel = 0.05490d0
      zns = 1.19459d-5
      znl = 1.5835218d-4
      c1ss = 2.9864797d-6
      c1l  = 4.7968065d-7

c     Days from epoch (Jan 0.0, 1900) for the reference
c     Convert JD to days from 1 Jan 2000 12h (J2000)
      tut1 = (epoch_jd - 2451545.0d0) / 36525.d0

c     Solar mean elements
c     Mean anomaly of the Sun
      zmos = (6898.9785d0*tut1 + 6.2565837d0)
      zmos = dmod(zmos, twopi)
      if (zmos .lt. 0.d0) zmos = zmos + twopi

c     Mean anomaly of the Moon
      zmol = (8399.709d0*tut1 + 5.8351514d0)
      zmol = dmod(zmol, twopi)
      if (zmol .lt. 0.d0) zmol = zmol + twopi

      dsstate(1) = zmol
      dsstate(2) = zmos

c     Compute intermediate trig quantities from orbital elements
      siniq = dsin(ainc)
      cosiq = dcos(ainc)
      eqsq  = ecc0 * ecc0
      betasq = 1.d0 - eqsq
      rteqsq = dsqrt(betasq)

c     Auxiliary quantities
      xnoi = 1.d0 / n0dd
      sinomo = dsin(aop)
      cosomo = dcos(aop)

c     ---- Solar terms ----
c     Solar inclination/node: fixed at ~23.4 deg, node at 0
c     From Vallado dscom: zincs = 0.39785416, zcosg = 0.91744867, etc.
      sinim = 0.39785416d0
      cosim = 0.91744867d0
      sinomm = -0.98088458d0
      cosomm = 0.19349076d0

c     Compute Z-coefficients for solar perturbations
      a1s = sinomm*cosim + cosomm*sinim
      a3s = -cosomm*cosim + sinomm*sinim
      a7s = -cosomm*sinim - sinomm*cosim
      a8s = sinomm*(-sinim)
      a9s = sinomm*cosim*cosiq
      a10s = cosomm*cosiq
      a2s = cosiq*a3s - siniq*cosomm
      a4s = cosiq*a1s + siniq*sinomm
      a5s = siniq*a3s + cosiq*cosomm
      a6s = siniq*a1s - cosiq*sinomm

      x1 = a1s*cosomo + a2s*sinomo
      x2 = a3s*cosomo + a4s*sinomo
      x3 = -a1s*sinomo + a2s*cosomo
      x4 = -a3s*sinomo + a4s*cosomo
      x5 = a5s*sinomo
      x6 = a6s*sinomo
      x7 = a5s*cosomo
      x8 = a6s*cosomo

      z31 = 12.d0*x1*x1 - 3.d0*x3*x3
      z32 = 24.d0*x1*x2 - 6.d0*x3*x4
      z33 = 12.d0*x2*x2 - 3.d0*x4*x4
      z1  = 3.d0*(a1s*a1s + a2s*a2s) + z31*eqsq
      z3  = 3.d0*(a3s*a3s + a4s*a4s) + z33*eqsq
      z11 = -6.d0*a1s*a5s + eqsq*
     &      (-24.d0*x1*x7 - 6.d0*x3*x5)
      z13 = -6.d0*a3s*a6s + eqsq*
     &      (-24.d0*x2*x8 - 6.d0*x4*x6)
      z21 = 6.d0*a2s*a5s + eqsq*
     &      (24.d0*x1*x5 - 6.d0*x3*x7)
      z23 = 6.d0*a4s*a6s + eqsq*
     &      (24.d0*x2*x6 - 6.d0*x4*x8)
      z1  = z1 + z1 + betasq*z31
      z3  = z3 + z3 + betasq*z33
      z11 = z11 + z11 + betasq*(-24.d0*x1*x7 - 6.d0*x3*x5)
      z13 = z13 + z13 + betasq*(-24.d0*x2*x8 - 6.d0*x4*x6)
      z21 = z21 + z21 + betasq*(24.d0*x1*x5 - 6.d0*x3*x7)
      z23 = z23 + z23 + betasq*(24.d0*x2*x6 - 6.d0*x4*x8)

      ses = c1ss*zns*a5s*(6.d0*(siniq*a6s - cosiq*a5s)
     &    + eqsq*(-24.d0*x2*x6+24.d0*x1*x5+6.d0*x4*x8-6.d0*x3*x7))
      sis = c1ss*zns*(z11 + z13)
      sls = -zns*c1ss*(z1 + z3 - 14.d0 - 6.d0*eqsq)
      sghs = c1ss*zns*(z31 + z33 - 6.d0)
      shs  = -zns*c1ss*(z21 + z23)
      if (dabs(siniq) .lt. 1.d-12) shs = 0.d0

c     Apply solar perturbation scale
      se2  = 2.d0*c1ss*zns*ses/(n0dd*zes*zes)
      si2  = 2.d0*c1ss*zns*sis/(n0dd*zes*zes)
      sl2  = -2.d0*zns*sls/(n0dd*zes*zes)
      sgh2 = -2.d0*zns*sghs/(n0dd*zes*zes)
      sh2  = -2.d0*zns*shs/(n0dd*zes*zes)

c     Recompute Z-coefficients using zes (solar eccentricity)
c     Adjust per Vallado formulas
      se2 = ses*xnoi*2.d0
      si2 = sis*xnoi*2.d0
      sl2 = sls*xnoi*2.d0
      sgh2 = sghs*xnoi*2.d0
      sh2  = shs*xnoi*2.d0

c     Additional solar terms (se3, si3, sl3, sl4, sgh3, sgh4, sh3)
      se3  = 0.d0
      si3  = 0.d0
      sl3  = 0.d0
      sl4  = 0.d0
      sgh3 = 0.d0
      sgh4 = -18.d0*c1ss*zns*siniq*xnoi
      sh3  = 0.d0

c     Store solar coefficients: dsstate(3:14)
c     se2, se3, si2, si3, sl2, sl3, sl4, sgh2, sgh3, sgh4, sh2, sh3
      dsstate(3)  = se2
      dsstate(4)  = se3
      dsstate(5)  = si2
      dsstate(6)  = si3
      dsstate(7)  = sl2
      dsstate(8)  = sl3
      dsstate(9)  = sl4
      dsstate(10) = sgh2
      dsstate(11) = sgh3
      dsstate(12) = sgh4
      dsstate(13) = sh2
      dsstate(14) = sh3

c     ---- Lunar terms ----
c     Lunar inclination/node from Vallado dscom
      sinim  = 0.08726646d0
      cosim  = 0.99619470d0
      sinomm = -0.27255540d0
      cosomm = 0.96215590d0

c     Recompute Z-coefficients for lunar perturbations
      a1s = sinomm*cosim + cosomm*sinim
      a3s = -cosomm*cosim + sinomm*sinim
      a7s = -cosomm*sinim - sinomm*cosim
      a8s = sinomm*(-sinim)
      a9s = sinomm*cosim*cosiq
      a10s = cosomm*cosiq
      a2s = cosiq*a3s - siniq*cosomm
      a4s = cosiq*a1s + siniq*sinomm
      a5s = siniq*a3s + cosiq*cosomm
      a6s = siniq*a1s - cosiq*sinomm

      x1 = a1s*cosomo + a2s*sinomo
      x2 = a3s*cosomo + a4s*sinomo
      x3 = -a1s*sinomo + a2s*cosomo
      x4 = -a3s*sinomo + a4s*cosomo
      x5 = a5s*sinomo
      x6 = a6s*sinomo
      x7 = a5s*cosomo
      x8 = a6s*cosomo

      z31 = 12.d0*x1*x1 - 3.d0*x3*x3
      z32 = 24.d0*x1*x2 - 6.d0*x3*x4
      z33 = 12.d0*x2*x2 - 3.d0*x4*x4
      z1  = 3.d0*(a1s*a1s + a2s*a2s) + z31*eqsq
      z3  = 3.d0*(a3s*a3s + a4s*a4s) + z33*eqsq
      z11 = -6.d0*a1s*a5s + eqsq*
     &      (-24.d0*x1*x7 - 6.d0*x3*x5)
      z13 = -6.d0*a3s*a6s + eqsq*
     &      (-24.d0*x2*x8 - 6.d0*x4*x6)
      z21 = 6.d0*a2s*a5s + eqsq*
     &      (24.d0*x1*x5 - 6.d0*x3*x7)
      z23 = 6.d0*a4s*a6s + eqsq*
     &      (24.d0*x2*x6 - 6.d0*x4*x8)
      z1  = z1 + z1 + betasq*z31
      z3  = z3 + z3 + betasq*z33
      z11 = z11 + z11 + betasq*(-24.d0*x1*x7 - 6.d0*x3*x5)
      z13 = z13 + z13 + betasq*(-24.d0*x2*x8 - 6.d0*x4*x6)
      z21 = z21 + z21 + betasq*(24.d0*x1*x5 - 6.d0*x3*x7)
      z23 = z23 + z23 + betasq*(24.d0*x2*x6 - 6.d0*x4*x8)

      sel = c1l*znl*a5s*(6.d0*(siniq*a6s - cosiq*a5s)
     &    + eqsq*(-24.d0*x2*x6+24.d0*x1*x5+6.d0*x4*x8-6.d0*x3*x7))
      sil = c1l*znl*(z11 + z13)
      sll = -znl*c1l*(z1 + z3 - 14.d0 - 6.d0*eqsq)
      sghl = c1l*znl*(z31 + z33 - 6.d0)
      shl  = -znl*c1l*(z21 + z23)
      if (dabs(siniq) .lt. 1.d-12) shl = 0.d0

c     Lunar coefficients
      ee2  = sel*xnoi*2.d0
      xi2  = sil*xnoi*2.d0
      xl2  = sll*xnoi*2.d0
      xgh2 = sghl*xnoi*2.d0
      xh2  = shl*xnoi*2.d0
      ee3  = 0.d0
      xi3  = 0.d0
      xl3  = 0.d0
      xl4  = 0.d0
      xgh3 = 0.d0
      xgh4 = -18.d0*c1l*znl*siniq*xnoi
      xh3  = 0.d0

c     Store lunar coefficients: dsstate(15:26)
      dsstate(15) = ee2
      dsstate(16) = ee3
      dsstate(17) = xi2
      dsstate(18) = xi3
      dsstate(19) = xl2
      dsstate(20) = xl3
      dsstate(21) = xl4
      dsstate(22) = xgh2
      dsstate(23) = xgh3
      dsstate(24) = xgh4
      dsstate(25) = xh2
      dsstate(26) = xh3

      return
      end

c     ============================================================
      subroutine tle_ds_init(epoch_jd, ecc0, ainc, cosio, sinio,
     &                       raan, aop, n0dd, dsstate, iresfl)
c
c     SDP4 deep-space initialization: compute secular gravity
c     perturbation rates and detect resonance.
c
c     Based on Vallado's dsinit() from sgp4unit.cpp
c
c     Input
c       epoch_jd : Julian Date of TLE epoch
c       ecc0     : eccentricity
c       ainc     : inclination (rad)
c       cosio    : cos(inclination)
c       sinio    : sin(inclination)
c       raan     : RAAN (rad)
c       aop      : argument of perigee (rad)
c       n0dd     : recovered mean motion (rad/min)
c
c     Input/Output
c       dsstate(1:26) : from tle_ds_common (read)
c       dsstate(27:50): secular rates & resonance (written)
c
c     Output
c       iresfl : resonance flag (0=none, 1=sync, 2=half-day)
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision dsstate(50)
      double precision n0dd
      integer iresfl

      double precision pi, twopi, d2r
      double precision ne, ratio
      double precision zmol, zmos
      double precision se2,se3,si2,si3,sl2,sl3,sl4
      double precision sgh2,sgh3,sgh4,sh2,sh3
      double precision ee2,ee3,xi2,xi3,xl2,xl3,xl4
      double precision xgh2,xgh3,xgh4,xh2,xh3
      double precision ses, sis, sls, sghs, shs
      double precision sel, sil, sll, sghl, shl
      double precision dedt, didt, dmdt, dnodt, domdt
      double precision zns, znl, zes, zel

      pi    = 3.14159265358979324d0
      twopi = 2.d0 * pi
      d2r   = pi / 180.d0

      zns = 1.19459d-5
      znl = 1.5835218d-4
      zes = 0.01675d0
      zel = 0.05490d0

c     Retrieve stored values
      zmol = dsstate(1)
      zmos = dsstate(2)
      se2  = dsstate(3)
      se3  = dsstate(4)
      si2  = dsstate(5)
      si3  = dsstate(6)
      sl2  = dsstate(7)
      sl3  = dsstate(8)
      sl4  = dsstate(9)
      sgh2 = dsstate(10)
      sgh3 = dsstate(11)
      sgh4 = dsstate(12)
      sh2  = dsstate(13)
      sh3  = dsstate(14)
      ee2  = dsstate(15)
      ee3  = dsstate(16)
      xi2  = dsstate(17)
      xi3  = dsstate(18)
      xl2  = dsstate(19)
      xl3  = dsstate(20)
      xl4  = dsstate(21)
      xgh2 = dsstate(22)
      xgh3 = dsstate(23)
      xgh4 = dsstate(24)
      xh2  = dsstate(25)
      xh3  = dsstate(26)

c     Compute secular rates from 3rd-body perturbations
c     These are the time-averaged effects of Sun and Moon
c     Solar secular contributions evaluated at zmos (solar mean anomaly)
      szmos = dsin(zmos)
      czmos = dcos(zmos)
      s2mos = dsin(2.d0*zmos)
      c2mos = dcos(2.d0*zmos)

      ses = se2*szmos + se3*s2mos
      sis = si2*szmos + si3*s2mos
      sls = sl2*szmos + sl3*s2mos + sl4*szmos*czmos
      sghs = sgh2*szmos + sgh3*s2mos + sgh4*szmos*czmos
      shs = sh2*szmos + sh3*s2mos

c     Lunar secular contributions evaluated at zmol (lunar mean anomaly)
      szmol = dsin(zmol)
      czmol = dcos(zmol)
      s2mol = dsin(2.d0*zmol)
      c2mol = dcos(2.d0*zmol)

      sel = ee2*szmol + ee3*s2mol
      sil = xi2*szmol + xi3*s2mol
      sll = xl2*szmol + xl3*s2mol + xl4*szmol*czmol
      sghl = xgh2*szmol + xgh3*s2mol + xgh4*szmol*czmol
      shl = xh2*szmol + xh3*s2mol

c     Total secular rates (solar + lunar)
      dedt  = ses + sel
      didt  = sis + sil
      dmdt  = sls + sll
      dnodt = sghs + sghl
      domdt = shs + shl

c     Store secular rates: dsstate(27:31)
      dsstate(27) = dedt
      dsstate(28) = didt
      dsstate(29) = dmdt
      dsstate(30) = dnodt
      dsstate(31) = domdt

c     Store epoch JD for reference
      dsstate(32) = epoch_jd

c     Resonance detection
c     Earth rotation rate (rad/min)
      ne = twopi / 1440.d0
      ratio = n0dd / ne

      if (ratio .ge. 0.9d0 .and. ratio .le. 1.1d0) then
c        Synchronous resonance (GEO)
         iresfl = 1
      else if (ratio .ge. 1.9d0 .and. ratio .le. 2.1d0) then
c        Half-day resonance (GPS, Molniya)
         iresfl = 2
      else
c        No resonance
         iresfl = 0
      end if

c     Compute resonance terms
      if (iresfl .eq. 1) then
c        Synchronous resonance terms
c        del1, del2, del3 from gravity field harmonics
c        Simplified from Vallado dsinit
         cosisq = cosio*cosio
         g200 = 1.d0 + ecc0*ecc0*(-2.5d0 + 0.8125d0*ecc0*ecc0)
         g310 = 1.d0 + 2.d0*ecc0*ecc0
         g300 = 1.d0 + ecc0*ecc0*(-6.d0 + 6.60937d0*ecc0*ecc0)

c        Synchronous resonance coefficients
         f220 = 0.75d0*(1.d0 + cosio)*(1.d0 + cosio)
         f311 = 0.9375d0*sinio*sinio*(1.d0 + 3.d0*cosio)
     &        - 0.75d0*(1.d0 + cosio)
         f330 = 1.d0 + cosio
         f330 = 1.875d0*f330*f330*f330

         dsstate(33) = 3.d0*n0dd*n0dd*1.7891679d-6*f220*g200
         dsstate(34) = 3.d0*n0dd*n0dd*2.2123015d-7*f330*g300
         dsstate(35) = 3.d0*n0dd*n0dd*3.1088075d-7*f311*g310

c        xlamo: mean longitude at epoch
         dsstate(36) = dmod(raan + aop +
     &                 dmod(n0dd*(-epoch_jd + 2451545.d0)*1440.d0,
     &                 twopi) - zmos, twopi)
         if (dsstate(36) .lt. 0.d0)
     &      dsstate(36) = dsstate(36) + twopi

c        xfact for synchronous
         dsstate(47) = dmdt + 0.5d0*n0dd*(dsstate(33) + dsstate(34)
     &               + dsstate(35)) - ne + dnodt
c        xli, xni
         dsstate(48) = dsstate(36)
         dsstate(49) = n0dd
         dsstate(50) = 0.d0

      else if (iresfl .eq. 2) then
c        Half-day resonance terms
c        From Vallado dsinit: d2201..d5433
         cosisq = cosio*cosio
         eqsq = ecc0*ecc0
         eoc = ecc0*eqsq
         g201 = -0.306d0 - (ecc0 - 0.64d0)*0.44d0

         if (ecc0 .le. 0.65d0) then
            g211 = 3.616d0 - 13.2470d0*ecc0 + 16.2900d0*eqsq
            g310 = -19.302d0 + 117.390d0*ecc0 - 228.419d0*eqsq
     &           + 156.591d0*eoc
            g322 = -18.9068d0 + 109.7927d0*ecc0 - 214.6334d0*eqsq
     &           + 146.5816d0*eoc
            g410 = -41.122d0 + 242.694d0*ecc0 - 471.094d0*eqsq
     &           + 313.953d0*eoc
            g422 = -146.407d0 + 841.880d0*ecc0 - 1629.014d0*eqsq
     &           + 1083.435d0*eoc
            g520 = -532.114d0 + 3017.977d0*ecc0 - 5740.032d0*eqsq
     &           + 3708.276d0*eoc
         else
            g211 = -72.099d0 + 331.819d0*ecc0 - 508.738d0*eqsq
     &           + 266.724d0*eoc
            g310 = -346.844d0 + 1582.851d0*ecc0 - 2415.925d0*eqsq
     &           + 1246.113d0*eoc
            g322 = -342.585d0 + 1554.908d0*ecc0 - 2366.899d0*eqsq
     &           + 1215.972d0*eoc
            g410 = -1052.797d0 + 4758.686d0*ecc0 - 7193.992d0*eqsq
     &           + 3651.957d0*eoc
            g422 = -3581.690d0 + 16178.110d0*ecc0
     &           - 24462.770d0*eqsq + 12422.520d0*eoc
            if (ecc0 .gt. 0.715d0) then
               g520 = -5149.66d0 + 29936.92d0*ecc0
     &              - 54087.36d0*eqsq + 31324.56d0*eoc
            else
               g520 = 1464.74d0 - 4664.75d0*ecc0
     &              + 3763.64d0*eqsq
            end if
         end if

         if (ecc0 .lt. 0.7d0) then
            g533 = -919.228d0 + 4988.610d0*ecc0
     &           - 9064.770d0*eqsq + 5542.210d0*eoc
            g521 = -822.710d0 + 4568.600d0*ecc0
     &           - 8491.400d0*eqsq + 5337.520d0*eoc
            g532 = -853.666d0 + 4690.250d0*ecc0
     &           - 8624.770d0*eqsq + 5341.400d0*eoc
         else
            g533 = -37995.780d0 + 161616.520d0*ecc0
     &           - 229838.200d0*eqsq + 109377.940d0*eoc
            g521 = -51752.104d0 + 218913.950d0*ecc0
     &           - 309468.160d0*eqsq + 146349.420d0*eoc
            g532 = -40023.880d0 + 170470.890d0*ecc0
     &           - 242699.480d0*eqsq + 115605.820d0*eoc
         end if

c        Inclination functions for half-day resonance
         f220 = 0.75d0*(1.d0+2.d0*cosio+cosisq)
         f221 = 1.5d0*sinio*sinio
         f321 = 1.875d0*sinio*(1.d0-2.d0*cosio-3.d0*cosisq)
         f322 = -1.875d0*sinio*(1.d0+2.d0*cosio-3.d0*cosisq)
         f441 = 35.d0*sinio*sinio*f220
         f442 = 39.3750d0*sinio*sinio*sinio*sinio
         f522 = 9.84375d0*sinio*(sinio*sinio*
     &        (1.d0-2.d0*cosio-5.d0*cosisq)+
     &        0.33333333d0*(-2.d0+4.d0*cosio+6.d0*cosisq))
         f523 = sinio*(4.92187512d0*sinio*sinio*
     &        (-2.d0-4.d0*cosio+10.d0*cosisq)+
     &        6.56250012d0*(1.d0+2.d0*cosio-3.d0*cosisq))
         f542 = 29.53125d0*sinio*(2.d0-8.d0*cosio+
     &        cosisq*(-12.d0+8.d0*cosio+10.d0*cosisq))
         f543 = 29.53125d0*sinio*(-2.d0-8.d0*cosio+
     &        cosisq*(12.d0+8.d0*cosio-10.d0*cosisq))

c        d2201..d5433 resonance coefficients
         temp1 = 3.d0*n0dd*n0dd*1.7891679d-6
         dsstate(37) = temp1*f220*g201
         dsstate(38) = temp1*f221*g211
         temp1 = 3.d0*n0dd*n0dd*2.2123015d-7
         dsstate(39) = temp1*f321*g310
         dsstate(40) = temp1*f322*g322
         temp1 = 3.d0*n0dd*n0dd*3.1088075d-7
         dsstate(41) = 2.d0*temp1*f441*g410
         dsstate(42) = temp1*f442*g422
         temp1 = 3.d0*n0dd*n0dd*1.7891679d-6
         dsstate(43) = temp1*f522*g520
         dsstate(44) = temp1*f523*g520
         temp1 = 3.d0*n0dd*n0dd*2.2123015d-7
         dsstate(45) = 2.d0*temp1*f542*g521
         dsstate(46) = 2.d0*temp1*f543*g533

c        xfact for half-day resonance
         dsstate(47) = dmdt + dnodt + domdt
     &               + n0dd*(2.d0*dsstate(37) + dsstate(38)
     &               + 2.d0*dsstate(39) + dsstate(40))
     &               - ne*2.d0
c        xli, xni
         theta_arg = raan + aop + aop
         dsstate(48) = dmod(theta_arg +
     &                 dmod(n0dd*2.d0*(-epoch_jd + 2451545.d0)
     &                 *1440.d0, twopi) - 2.d0*zmos, twopi)
         if (dsstate(48) .lt. 0.d0)
     &      dsstate(48) = dsstate(48) + twopi
         dsstate(49) = n0dd
         dsstate(50) = 0.d0
      else
c        No resonance: zero out resonance slots
         do i = 33, 50
            dsstate(i) = 0.d0
         end do
      end if

      return
      end

c     ============================================================
      subroutine tle_ds_lpper(tsince, ecc0, ainc_r, dsstate,
     &                        ep, xincp, xnodep, omegap, mp)
c
c     SDP4 long-period periodics from lunar and solar perturbations.
c     Applied at each propagation step.
c
c     Based on Vallado's dpper() from sgp4unit.cpp
c
c     Input
c       tsince : minutes from epoch
c       ecc0   : original eccentricity (not used, ep is modified)
c       ainc_r : original inclination in radians
c       dsstate: state array
c
c     Input/Output (modified in place)
c       ep     : eccentricity
c       xincp  : inclination (rad)
c       xnodep : RAAN (rad)
c       omegap : argument of perigee (rad)
c       mp     : mean anomaly (rad)
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision dsstate(50)

      double precision pi, twopi
      double precision zns, znl, zes, zel
      double precision zmol, zmos
      double precision se2,se3,si2,si3,sl2,sl3,sl4
      double precision sgh2,sgh3,sgh4,sh2,sh3
      double precision ee2,ee3,xi2,xi3,xl2,xl3,xl4
      double precision xgh2,xgh3,xgh4,xh2,xh3
      double precision zm, zf, sinzf, coszf, f2, f3

      pi    = 3.14159265358979324d0
      twopi = 2.d0 * pi

      zns = 1.19459d-5
      znl = 1.5835218d-4
      zes = 0.01675d0
      zel = 0.05490d0

c     Retrieve coefficients
      zmol = dsstate(1)
      zmos = dsstate(2)
      se2  = dsstate(3)
      se3  = dsstate(4)
      si2  = dsstate(5)
      si3  = dsstate(6)
      sl2  = dsstate(7)
      sl3  = dsstate(8)
      sl4  = dsstate(9)
      sgh2 = dsstate(10)
      sgh3 = dsstate(11)
      sgh4 = dsstate(12)
      sh2  = dsstate(13)
      sh3  = dsstate(14)
      ee2  = dsstate(15)
      ee3  = dsstate(16)
      xi2  = dsstate(17)
      xi3  = dsstate(18)
      xl2  = dsstate(19)
      xl3  = dsstate(20)
      xl4  = dsstate(21)
      xgh2 = dsstate(22)
      xgh3 = dsstate(23)
      xgh4 = dsstate(24)
      xh2  = dsstate(25)
      xh3  = dsstate(26)

c     Update solar/lunar mean anomalies for current time
c     tsince is in minutes
      zm = zmos + zns * tsince
      zf = zm + 2.d0*zes*dsin(zm)
      sinzf = dsin(zf)
      coszf = dcos(zf)
      f2 = 0.5d0*sinzf*sinzf - 0.25d0
      f3 = -0.5d0*sinzf*coszf

c     Solar periodics
      pe_s  = se2*f2 + se3*f3
      pinc_s = si2*f2 + si3*f3
      pl_s  = sl2*f2 + sl3*f3 + sl4*sinzf
      pgh_s = sgh2*f2 + sgh3*f3 + sgh4*sinzf
      ph_s  = sh2*f2 + sh3*f3

c     Lunar periodics
      zm = zmol + znl * tsince
      zf = zm + 2.d0*zel*dsin(zm)
      sinzf = dsin(zf)
      coszf = dcos(zf)
      f2 = 0.5d0*sinzf*sinzf - 0.25d0
      f3 = -0.5d0*sinzf*coszf

      pe_l  = ee2*f2 + ee3*f3
      pinc_l = xi2*f2 + xi3*f3
      pl_l  = xl2*f2 + xl3*f3 + xl4*sinzf
      pgh_l = xgh2*f2 + xgh3*f3 + xgh4*sinzf
      ph_l  = xh2*f2 + xh3*f3

c     Total periodics
      pe  = pe_s + pe_l
      pinc = pinc_s + pinc_l
      pl  = pl_s + pl_l
      pgh = pgh_s + pgh_l
      ph  = ph_s + ph_l

c     Apply corrections
      ep     = ep + pe
      xincp  = xincp + pinc
      omegap = omegap + pgh
      mp     = mp + pl

c     Apply node correction accounting for inclination
      sinip = dsin(xincp)
      cosip = dcos(xincp)
      if (dabs(sinip) .gt. 1.d-12) then
         xnodep = xnodep + ph / sinip
         omegap = omegap - ph * cosip / sinip
      end if

c     Clamp eccentricity
      if (ep .lt. 1.d-6) ep = 1.d-6
      if (ep .gt. 0.999d0) ep = 0.999d0

      return
      end

c     ============================================================
      subroutine tle_ds_resint(tsince, n0dd, ecc0, dsstate,
     &                         iresfl, ep, nm)
c
c     SDP4 resonance ODE integration using 720-second Euler steps.
c
c     Based on Vallado's dspace() from sgp4unit.cpp
c
c     Input
c       tsince : minutes from epoch
c       n0dd   : recovered mean motion (rad/min)
c       ecc0   : original eccentricity
c       dsstate: state array
c       iresfl : resonance flag (1=sync, 2=half-day)
c
c     Input/Output
c       ep     : eccentricity (modified for half-day resonance)
c       nm     : mean motion (modified by resonance)
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision dsstate(50)
      double precision n0dd, nm
      integer iresfl

      double precision pi, twopi
      double precision del1, del2, del3, xlamo
      double precision xfact, xli, xni, atime
      double precision d2201,d2211,d3210,d3222,d4410,d4422
      double precision d5220,d5232,d5421,d5433
      double precision delt, ft, xndt, xldot, xnddt

      pi    = 3.14159265358979324d0
      twopi = 2.d0 * pi

c     Integration step size (720 seconds = 12 minutes)
      delt = 720.d0 / 60.d0

c     Retrieve resonance terms
      xfact = dsstate(47)
      xli   = dsstate(48)
      xni   = dsstate(49)
      atime = dsstate(50)

      if (iresfl .eq. 1) then
c        Synchronous resonance
         del1  = dsstate(33)
         del2  = dsstate(34)
         del3  = dsstate(35)
         xlamo = dsstate(36)

c        Compute time derivatives
         ft = dmod(xlamo + xfact*tsince, twopi)
         xndt = del1*dsin(ft - xfact*atime)
     &        + del2*dsin(2.d0*(ft - xfact*atime))
     &        + del3*dsin(3.d0*(ft - xfact*atime))
         xnddt = del1*dcos(ft - xfact*atime)
     &         + 2.d0*del2*dcos(2.d0*(ft - xfact*atime))
     &         + 3.d0*del3*dcos(3.d0*(ft - xfact*atime))
         xnddt = xnddt * xfact

c        Simple integration: update mean motion
         nm = n0dd + xndt

      else if (iresfl .eq. 2) then
c        Half-day resonance
         d2201 = dsstate(37)
         d2211 = dsstate(38)
         d3210 = dsstate(39)
         d3222 = dsstate(40)
         d4410 = dsstate(41)
         d4422 = dsstate(42)
         d5220 = dsstate(43)
         d5232 = dsstate(44)
         d5421 = dsstate(45)
         d5433 = dsstate(46)

c        Euler integration from atime to tsince
c        Step in direction of tsince
         stepp = delt
         if (tsince .lt. atime) stepp = -delt

c        Integrate step by step
 100     continue
         if (dabs(tsince - atime) .lt. delt) goto 200

c        Compute resonance acceleration at current xli
         xndt = d2201*dsin(2.d0*xli - 2.d0*
     &          dsstate(2)) +
     &          d2211*dsin(xli - dsstate(2)) +
     &          d3210*dsin(xli + dsstate(2)) +
     &          d3222*dsin(-xli + dsstate(2)) +
     &          d5220*dsin(xli - dsstate(1)) +
     &          d5232*dsin(-xli + dsstate(1)) +
     &          d4410*dsin(2.d0*xli - 2.d0*dsstate(1)) +
     &          d4422*dsin(2.d0*xli - dsstate(1)
     &          - dsstate(2)) +
     &          d5421*dsin(xli + dsstate(2)
     &          - dsstate(1)) +
     &          d5433*dsin(-xli + dsstate(2)
     &          - dsstate(1))

c        Euler step
         xli = xli + stepp*(xni + xndt)
         xni = xni + stepp*xndt
         atime = atime + stepp
         goto 100

 200     continue

c        Fractional last step
         ft = tsince - atime
         xndt = d2201*dsin(2.d0*xli - 2.d0*
     &          dsstate(2)) +
     &          d2211*dsin(xli - dsstate(2)) +
     &          d3210*dsin(xli + dsstate(2)) +
     &          d3222*dsin(-xli + dsstate(2)) +
     &          d5220*dsin(xli - dsstate(1)) +
     &          d5232*dsin(-xli + dsstate(1)) +
     &          d4410*dsin(2.d0*xli - 2.d0*dsstate(1)) +
     &          d4422*dsin(2.d0*xli - dsstate(1)
     &          - dsstate(2)) +
     &          d5421*dsin(xli + dsstate(2)
     &          - dsstate(1)) +
     &          d5433*dsin(-xli + dsstate(2)
     &          - dsstate(1))

         nm = xni + xndt
         xli = xli + ft*(xni + xndt)
         atime = tsince
      end if

c     Store updated resonance state
      dsstate(48) = xli
      dsstate(49) = xni
      dsstate(50) = atime

c     Keep mean motion positive
      if (nm .le. 0.d0) nm = 1.d-6

      return
      end

c     ============================================================
      subroutine tle_sdp4_init(amm_revday, ecc0, ainc_deg,
     &                    raan_deg, aop_deg, ama_deg,
     &                    bstar, epoch_jd,
     &                    n0dd, a0dd, cosio, sinio,
     &                    eta, c1, c4, t0, perigee,
     &                    omgdot, xnodot, xmdot,
     &                    dsstate, iresfl, ierr)
c
c     SDP4 initialization for deep-space orbits (period >= 225 min).
c     Recovers mean elements, computes drag/secular coefficients,
c     then initializes deep-space perturbation state.
c
c     Input (TLE-parsed values)
c       amm_revday : mean motion (rev/day)
c       ecc0       : eccentricity
c       ainc_deg   : inclination (deg)
c       raan_deg   : RAAN (deg)
c       aop_deg    : argument of perigee (deg)
c       ama_deg    : mean anomaly (deg)
c       bstar      : B* drag term (1/Earth radii)
c       epoch_jd   : epoch Julian Date
c
c     Output
c       n0dd   : recovered mean motion (rad/min)
c       a0dd   : recovered semi-major axis (Earth radii)
c       cosio, sinio : trig of inclination
c       eta    : eccentricity parameter
c       c1     : drag coefficient
c       c4     : drag coefficient for perigee
c       t0     : epoch time (JD)
c       perigee: perigee altitude (km)
c       omgdot : secular rate of arg of perigee (rad/min)
c       xnodot : secular rate of RAAN (rad/min)
c       xmdot  : secular rate of mean anomaly (rad/min)
c       dsstate(50) : deep-space state array
c       iresfl : resonance flag (0/1/2)
c       ierr   : 0=ok, 2=bad elements
c     ============================================================
      implicit double precision (a-h, o-z)
      integer ierr, iresfl
      double precision dsstate(50)

c     Explicit declarations for i-n range variables
      double precision n0dd, n0

c     WGS-72 constants (canonical for SGP4/SDP4)
      double precision ke, j2, j3, j4, re_km
      parameter (re_km = 6378.135d0)
      parameter (ke    = 0.0743669161d0)
      parameter (j2    = 1.0826158d-3)
      parameter (j3    = -2.53881d-6)
      parameter (j4    = -1.65597d-6)

      double precision k2, k4
      double precision del1, del0
      double precision theta2, theta4, beta0, beta02
      double precision qoms24
      double precision coef, coef1, eeta
      double precision psisq, tsi

      ierr = 0
      pi    = 3.14159265358979324d0
      twopi = 2.d0 * pi
      d2r   = pi / 180.d0

c     Derived gravity constants
      k2  = 0.5d0 * j2 * re_km**2
      k4  = -0.375d0 * j4 * re_km**4

c     Convert inputs to internal units
      ainc = ainc_deg * d2r
      cosio = dcos(ainc)
      sinio = dsin(ainc)
      theta2 = cosio * cosio
      theta4 = theta2 * theta2

c     Mean motion: rev/day -> rad/min
      n0 = amm_revday * twopi / 1440.d0

c     Recover original mean motion and semi-major axis
      a1 = (ke / n0) ** (2.d0/3.d0)
      del1 = 1.5d0*(k2/(re_km**2))/(a1*a1)
     &      * (3.d0*theta2-1.d0)
     &      / (1.d0-ecc0*ecc0)**1.5d0
      a0 = a1 * (1.d0 - del1/3.d0 - del1**2
     &         - 134.d0/81.d0 * del1**3)
      del0 = 1.5d0*(k2/(re_km**2))/(a0*a0)
     &      * (3.d0*theta2-1.d0)
     &      / (1.d0-ecc0*ecc0)**1.5d0

      n0dd = n0 / (1.d0 + del0)
      a0dd = a0 / (1.d0 - del0)

c     No period rejection for deep-space (that's the whole point)

c     Perigee altitude (km)
      perigee = (a0dd * (1.d0 - ecc0) - 1.d0) * re_km

c     Atmospheric parameter s
      s_param = 78.d0 / re_km + 1.d0
      if (perigee .lt. 156.d0) then
         if (perigee .lt. 98.d0) then
            s_param = 20.d0 / re_km + 1.d0
         else
            s_param = (perigee - 78.d0) / re_km + 1.d0
         end if
      end if

c     qoms2t = (q0 - s)^4 with q0 = 120 km
      qoms24 = ((120.d0 - (s_param - 1.d0)*re_km) / re_km) ** 4
      if (qoms24 .lt. 1.d-24) qoms24 = 1.d-24

      tsi = 1.d0 / (a0dd - s_param)
      eta = a0dd * ecc0 * tsi
      eeta = ecc0 * eta

      beta0 = dsqrt(1.d0 - ecc0*ecc0)
      beta02 = beta0 * beta0

      psisq = dabs(1.d0 - eta*eta)
      if (psisq .lt. 1.d-20) psisq = 1.d-20

      coef  = qoms24 * tsi**4
      coef1 = coef / psisq**3.5d0

c     C2 coefficient
      c2 = coef1 * n0dd * (a0dd * (1.d0 + 1.5d0*eta*eta + eeta*
     &     (4.d0 + eta*eta)) + 0.75d0 * (k2/(re_km**2)) * tsi /
     &     psisq * (-0.5d0 + 1.5d0*theta2) *
     &     (8.d0 + 3.d0*eta*eta*(8.d0 + eta*eta)))

c     C1 coefficient
      c1 = bstar * c2

c     C4 coefficient
      c4 = 2.d0*n0dd*coef1*a0dd*beta02*
     &     (eta*(2.d0+0.5d0*eta*eta) +
     &      ecc0*(0.5d0+2.d0*eta*eta) -
     &      2.d0*(k2/(re_km**2))*tsi/(a0dd*psisq)*
     &      (-3.d0*(1.d0-3.d0*theta2)*(1.d0+1.5d0*eta*eta
     &       -2.d0*eeta*(1.d0+eta*eta)) +
     &       0.75d0*(1.d0-theta2)*(2.d0*eta*eta
     &       -eeta*(1.d0+eta*eta))))

c     Secular rates due to J2, J4
      pinvsq = 1.d0 / (a0dd*a0dd * beta02*beta02)
      temp1 = 1.5d0 * j2 * pinvsq * n0dd
      temp2 = 0.5d0 * temp1 * j2 * pinvsq
      temp3 = -0.46875d0 * j4 * pinvsq * pinvsq * n0dd

      xmdot = n0dd + 0.5d0*temp1*beta0*(3.d0*theta2-1.d0)
     &      + 0.0625d0*temp2*beta0*
     &        (13.d0-78.d0*theta2+137.d0*theta4)

      omgdot = -0.5d0*temp1*(1.d0-5.d0*theta2)
     &       + 0.0625d0*temp2*(7.d0-114.d0*theta2+395.d0*theta4)
     &       + temp3*(3.d0-36.d0*theta2+49.d0*theta4)

      xnodot = -(temp1*cosio)
     &       + (0.5d0*temp2*(4.d0-19.d0*theta2)
     &       + 2.d0*temp3*(3.d0-7.d0*theta2)) * cosio

      t0 = epoch_jd

c     Initialize deep-space state
      call tle_ds_common(epoch_jd, ecc0, ainc,
     &                   raan_deg*d2r, aop_deg*d2r,
     &                   n0dd, dsstate)
      call tle_ds_init(epoch_jd, ecc0, ainc, cosio, sinio,
     &                 raan_deg*d2r, aop_deg*d2r,
     &                 n0dd, dsstate, iresfl)

      return
      end

c     ============================================================
      subroutine tle_sdp4_prop(tsince,
     &              ecc0, ainc_deg, raan_deg, aop_deg, ama_deg,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1, c4, perigee,
     &              omgdot, xnodot, xmdot, bstar,
     &              dsstate, iresfl,
     &              r_teme, v_teme, ierr)
c
c     SDP4 propagation at tsince minutes from epoch.
c
c     Output
c       r_teme(3) : position in TEME frame (km)
c       v_teme(3) : velocity in TEME frame (km/s)
c       ierr      : 0=ok, nonzero=error
c     ============================================================
      implicit double precision (a-h, o-z)
      integer ierr, i, iresfl
      double precision r_teme(3), v_teme(3)
      double precision dsstate(50)

c     Explicit declarations for i-n range arguments/locals
      double precision n0dd, k2, mrt, nm

      double precision j2, j3, ke, re_km
      parameter (re_km = 6378.135d0)
      parameter (ke    = 0.0743669161d0)
      parameter (j2    = 1.0826158d-3)
      parameter (j3    = -2.53881d-6)

      double precision mp, omega, xnode, ecc_t, a_t
      double precision axn, ayn, aynl
      double precision xl, u, eo1, tem5, sineo1, coseo1
      double precision ecose, esine, el2, pl, rl, rdotl, rvdotl
      double precision betal, sinu, cosu, su, sin2u, cos2u
      double precision rk, uk, xnodek, xinck
      double precision sinuk, cosuk, sinik, cosik, sinnok, cosnok
      double precision xmx, xmy, ux, uy, uz, vx, vy, vz

      ierr = 0
      pi    = 3.14159265358979324d0
      twopi = 2.d0 * pi
      d2r   = pi / 180.d0

      k2 = 0.5d0 * j2 * re_km**2

c     Convert initial angles to radians
      omega0 = aop_deg * d2r
      xnode0 = raan_deg * d2r
      xm0    = ama_deg * d2r
      ainc   = ainc_deg * d2r

c     SGP4 secular updates (J2/J4 rates)
      mp    = xm0 + xmdot * tsince
      omega = omega0 + omgdot * tsince
      xnode = xnode0 + xnodot * tsince

c     Add deep-space secular perturbations (3rd body)
      dedt  = dsstate(27)
      didt  = dsstate(28)
      dmdt  = dsstate(29)
      dnodt = dsstate(30)
      domdt = dsstate(31)

      mp    = mp + dmdt * tsince
      omega = omega + domdt * tsince
      xnode = xnode + dnodt * tsince
      xincp = ainc + didt * tsince

c     Drag secular effects on a and e
      tempa = 1.d0 - c1 * tsince
      tempe = bstar * c4 * tsince
      a_t   = a0dd * tempa * tempa
      ecc_t = ecc0 - tempe + dedt * tsince

c     Clamp eccentricity
      if (ecc_t .lt. 1.d-6) ecc_t = 1.d-6
      if (ecc_t .gt. 0.999d0) ecc_t = 0.999d0

c     Mean motion from semi-major axis
      nm = ke / (a_t * dsqrt(a_t))

c     Resonance integration
      if (iresfl .ne. 0) then
         call tle_ds_resint(tsince, n0dd, ecc0, dsstate,
     &                      iresfl, ecc_t, nm)
c        Re-clamp after resonance
         if (ecc_t .lt. 1.d-6) ecc_t = 1.d-6
         if (ecc_t .gt. 0.999d0) ecc_t = 0.999d0
c        Recalculate a_t from nm
         a_t = (ke / nm) ** (2.d0/3.d0)
      end if

c     Long-period lunar-solar periodics
      call tle_ds_lpper(tsince, ecc0, ainc, dsstate,
     &                  ecc_t, xincp, xnode, omega, mp)

c     Re-clamp eccentricity after periodics
      if (ecc_t .lt. 1.d-6) ecc_t = 1.d-6
      if (ecc_t .gt. 0.999d0) ecc_t = 0.999d0

c     Update cosio/sinio from perturbed inclination
      cosip = dcos(xincp)
      sinip = dsin(xincp)

c     Long-period periodics (J3 term)
      axn  = ecc_t * dcos(omega)
      pl   = a_t * (1.d0 - ecc_t*ecc_t)
      if (pl .lt. 1.d-12) pl = 1.d-12
      aynl = -0.5d0 * (j3/j2) * sinip / pl
      ayn  = ecc_t * dsin(omega) + aynl

c     Solve Kepler's equation for (E + omega)
      xl = mp + omega + xnode
      u  = dmod(xl - xnode, twopi)
      if (u .lt. 0.d0) u = u + twopi

      eo1 = u
      do i = 1, 10
         sineo1 = dsin(eo1)
         coseo1 = dcos(eo1)
         tem5 = 1.d0 - coseo1*axn - sineo1*ayn
         if (dabs(tem5) .lt. 1.d-20) tem5 = 1.d-20
         tem5 = (u - ayn*coseo1 + axn*sineo1 - eo1) / tem5
         if (dabs(tem5) .gt. 0.95d0) then
            if (tem5 .gt. 0.d0) then
               tem5 = 0.95d0
            else
               tem5 = -0.95d0
            end if
         end if
         eo1 = eo1 + tem5
      end do

c     Short-period preliminary quantities
      sineo1 = dsin(eo1)
      coseo1 = dcos(eo1)
      ecose = axn*coseo1 + ayn*sineo1
      esine = axn*sineo1 - ayn*coseo1
      el2   = axn*axn + ayn*ayn
      pl    = a_t * (1.d0 - el2)
      if (pl .lt. 0.d0) then
         ierr = 2
         return
      end if
      rl    = a_t * (1.d0 - ecose)
      rdotl = ke * dsqrt(a_t) * esine / rl
      rvdotl = ke * dsqrt(pl) / rl

      betal = dsqrt(1.d0 - el2)
      sinu  = a_t / rl * (sineo1 - ayn - axn*esine/(1.d0+betal))
      cosu  = a_t / rl * (coseo1 - axn + ayn*esine/(1.d0+betal))
      su    = datan2(sinu, cosu)
      sin2u = 2.d0 * sinu * cosu
      cos2u = 2.d0 * cosu * cosu - 1.d0

      theta2 = cosip * cosip

c     Short-period J2 corrections
      rk = rl*(1.d0 - 1.5d0*(k2/(re_km**2))*dsqrt(1.d0-el2)/
     &     (pl)*(3.d0*theta2-1.d0)) +
     &     0.5d0*(k2/(re_km**2))*(1.d0-theta2)/pl * cos2u
      uk = su - 0.25d0*(k2/(re_km**2))/(pl*pl)*
     &     (7.d0*theta2-1.d0) * sin2u
      xnodek = xnode + 1.5d0*(k2/(re_km**2))*cosip/(pl*pl)*sin2u
      xinck_f = xincp + 1.5d0*(k2/(re_km**2))*cosip*sinip/
     &     (pl*pl)*cos2u

c     Unit orientation vectors
      sinuk = dsin(uk)
      cosuk = dcos(uk)
      sinik = dsin(xinck_f)
      cosik = dcos(xinck_f)
      sinnok = dsin(xnodek)
      cosnok = dcos(xnodek)

      xmx = -sinnok * cosik
      xmy =  cosnok * cosik

      ux = xmx*sinuk + cosnok*cosuk
      uy = xmy*sinuk + sinnok*cosuk
      uz = sinik*sinuk

      vx = xmx*cosuk - cosnok*sinuk
      vy = xmy*cosuk - sinnok*sinuk
      vz = sinik*cosuk

c     Position (km) and velocity (km/s)
      mrt = rk * re_km
      if (mrt .lt. 1.d0) then
         ierr = 3
         return
      end if

      r_teme(1) = mrt * ux
      r_teme(2) = mrt * uy
      r_teme(3) = mrt * uz

      rdotk  = rdotl * re_km / 60.d0
      rfdotk = rvdotl * re_km / 60.d0

      v_teme(1) = rdotk*ux + rfdotk*vx
      v_teme(2) = rdotk*uy + rfdotk*vy
      v_teme(3) = rdotk*uz + rfdotk*vz

      return
      end

c     ============================================================
      subroutine tle2sv(iyr4, eday, ainc_deg, raan_deg, ecc0,
     &                  aop_deg, ama_deg, amm_revday, bstar,
     &                  r_j2k, v_j2k, ierr)
c
c     Complete TLE-to-state-vector conversion.
c     Dispatches between SGP4 (near-Earth) and SDP4 (deep-space)
c     based on orbital period, then rotates TEME -> J2000.
c
c     Input (directly from TLE parse)
c       iyr4       : 4-digit epoch year
c       eday       : epoch day of year + fractional day
c       ainc_deg   : inclination (deg)
c       raan_deg   : RAAN (deg)
c       ecc0       : eccentricity
c       aop_deg    : argument of perigee (deg)
c       ama_deg    : mean anomaly (deg)
c       amm_revday : mean motion (rev/day)
c       bstar      : B* drag term (1/Earth radii)
c
c     Output
c       r_j2k(3)   : position in J2000 (km)
c       v_j2k(3)   : velocity in J2000 (km/s)
c       ierr       : 0=ok, 2=bad elements
c     ============================================================
      implicit double precision (a-h, o-z)
      integer iyr4, ierr, iresfl
      double precision r_j2k(3), v_j2k(3)
      double precision r_teme(3), v_teme(3)

      double precision n0dd, a0dd, cosio, sinio
      double precision eta, c1_coeff, c4_coeff, t0, perigee_km
      double precision omgdot, xnodot, xmdot
      double precision djd
      double precision dsstate(50)

c     Convert epoch to JD
      call tle_epoch2jd(iyr4, eday, djd)

c     Determine orbit type from period
      period_min = 1440.d0 / amm_revday

      if (period_min .lt. 225.d0) then
c        Near-Earth: use SGP4
         call tle_sgp4_init(amm_revday, ecc0, ainc_deg,
     &                      raan_deg, aop_deg, ama_deg,
     &                      bstar, djd,
     &                      n0dd, a0dd, cosio, sinio,
     &                      eta, c1_coeff, c4_coeff, t0,
     &                      perigee_km,
     &                      omgdot, xnodot, xmdot, ierr)
         if (ierr .ne. 0) return

         call tle_sgp4_prop(0.d0,
     &                 ecc0, ainc_deg, raan_deg, aop_deg, ama_deg,
     &                 n0dd, a0dd, cosio, sinio,
     &                 eta, c1_coeff, c4_coeff, perigee_km,
     &                 omgdot, xnodot, xmdot, bstar,
     &                 r_teme, v_teme, ierr)
         if (ierr .ne. 0) return

      else
c        Deep-space: use SDP4
         call tle_sdp4_init(amm_revday, ecc0, ainc_deg,
     &                      raan_deg, aop_deg, ama_deg,
     &                      bstar, djd,
     &                      n0dd, a0dd, cosio, sinio,
     &                      eta, c1_coeff, c4_coeff, t0,
     &                      perigee_km,
     &                      omgdot, xnodot, xmdot,
     &                      dsstate, iresfl, ierr)
         if (ierr .ne. 0) return

         call tle_sdp4_prop(0.d0,
     &                 ecc0, ainc_deg, raan_deg, aop_deg, ama_deg,
     &                 n0dd, a0dd, cosio, sinio,
     &                 eta, c1_coeff, c4_coeff, perigee_km,
     &                 omgdot, xnodot, xmdot, bstar,
     &                 dsstate, iresfl,
     &                 r_teme, v_teme, ierr)
         if (ierr .ne. 0) return
      end if

c     Transform TEME -> J2000
      call tle_teme2j2k(r_teme, v_teme, djd, r_j2k, v_j2k)

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
