      program test_ephem
c     Print KSROP analytic Sun/Moon positions at 13 epochs
c     (5-day steps from 2016-09-20), for comparison against
c     GMAT DE405 (ephemeris_check_gmat.script).
      implicit double precision (a-h,o-z)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
      dimension ts(3), tm(3)
      pi=4.d0*datan(1.d0)
      d2r=pi/180.d0
      r2d=180.d0/pi
      amue=3.986004415d5
      AU=1.495978707d08
      R_Earth=6378.1363d0
c     2016-09-20T00:00:00.000 UTC julian date
      dj0 = 2457651.5d0
      do k = 0, 12
         djulian = dj0 + 5.d0*dble(k)
         call solarnpv(djulian, ts)
         call lunarpv(djulian, tm)
         write(*,100) k, ts(1), ts(2), ts(3), tm(1), tm(2), tm(3)
      end do
 100  format(I3, 6F20.6)
      end
