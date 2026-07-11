      program test_ephem
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
      djulian = 2457651.5d0
      call solarnpv(djulian, ts)
      call lunarpv(djulian, tm)
      write(*,*) 'SUN: ', ts(1), ts(2), ts(3)
      write(*,*) 'MOON:', tm(1), tm(2), tm(3)
      end
