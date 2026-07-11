      program test_alegp
      implicit double precision (a-h,o-z)
      dimension v(10)
      dimension xs(9)
      xs(1) = -0.99d0
      xs(2) = -0.8d0
      xs(3) = -0.5d0
      xs(4) = -0.1515d0
      xs(5) = 0.0d0
      xs(6) = 0.1515d0
      xs(7) = 0.5d0
      xs(8) = 0.8d0
      xs(9) = 0.99d0
      do k = 1, 9
         x = xs(k)
         call aLegP(5, x, v)
         p0 = 1.d0
         p1 = x
         p2 = 0.5d0*(3.d0*x*x-1.d0)
         p3 = 0.5d0*(5.d0*x**3-3.d0*x)
         p4 = 0.125d0*(35.d0*x**4-30.d0*x**2+3.d0)
         p5 = 0.125d0*(63.d0*x**5-70.d0*x**3+15.d0*x)
         write(*,100) x, v(1)-p0, v(2)-p1, v(3)-p2, v(4)-p3,
     &                v(5)-p4, v(6)-p5
      end do
 100  format(F8.4, 6E14.4)
      end
