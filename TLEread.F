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
      subroutine tle2sv(iyr4, eday, ainc_deg, raan_deg, ecc0,
     &                  aop_deg, ama_deg, amm_revday, bstar,
     &                  r_j2k, v_j2k, ierr)
c
c     Complete TLE-to-state-vector conversion.
c     Evaluates SGP4 at t=0 (epoch), then rotates TEME -> J2000.
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
c       ierr       : 0=ok, 1=deep-space, 2=bad elements
c     ============================================================
      implicit double precision (a-h, o-z)
      integer iyr4, ierr
      double precision r_j2k(3), v_j2k(3)
      double precision r_teme(3), v_teme(3)

      double precision n0dd, a0dd, cosio, sinio
      double precision eta, c1_coeff, c4_coeff, t0, perigee_km
      double precision omgdot, xnodot, xmdot
      double precision djd

c     Convert epoch to JD
      call tle_epoch2jd(iyr4, eday, djd)

c     Initialize SGP4
      call tle_sgp4_init(amm_revday, ecc0, ainc_deg,
     &                   raan_deg, aop_deg, ama_deg,
     &                   bstar, djd,
     &                   n0dd, a0dd, cosio, sinio,
     &                   eta, c1_coeff, c4_coeff, t0, perigee_km,
     &                   omgdot, xnodot, xmdot, ierr)
      if (ierr .ne. 0) return

c     Propagate at t=0
      call tle_sgp4_prop(0.d0,
     &              ecc0, ainc_deg, raan_deg, aop_deg, ama_deg,
     &              n0dd, a0dd, cosio, sinio,
     &              eta, c1_coeff, c4_coeff, perigee_km,
     &              omgdot, xnodot, xmdot, bstar,
     &              r_teme, v_teme, ierr)
      if (ierr .ne. 0) return

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
