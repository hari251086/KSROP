      subroutine p_polynomial_value ( m, n, x, v )

c*********************************************************************72
c
cc P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
c
c  Discussion:
c
c    P(n,1) = 1.
c    P(n,-1) = (-1)^N.
c    | P(n,x) | .le. 1 in [-1,1].
c
c    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
c    quadrature of the integral of a function F(X) with weight function 1
c    over the interval [-1,1].
c
c    The Legendre polynomials are orthogonal under the inner product defined
c    as integration from -1 to 1:
c
c      Integral ( -1 .le. X .le. 1 ) P(I,X) * P(J,X) dX 
c        = 0 if I =/= J
c        = 2 / ( 2*I+1 ) if I = J.
c
c    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
c
c    A function F(X) defined on [-1,1] may be approximated by the series
c      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
c    where
c      C(I) = (2*I+1)/(2) * Integral ( -1 .le. X .le. 1 ) F(X) P(I,x) dx.
c
c    The formula is:
c
c      P(n,x) = (1/2^N) * sum ( 0 .le. M .le. N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
c
c  Differential equation:
c
c    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
c
c  First terms:
c
c    P( 0,x) =      1
c    P( 1,x) =      1 X
c    P( 2,x) = (    3 X^2 -       1)/2
c    P( 3,x) = (    5 X^3 -     3 X)/2
c    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
c    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
c    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
c    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
c    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
c    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
c    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
c
c  Recursion:
c
c    P(0,x) = 1
c    P(1,x) = x
c    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
c
c    P'(0,x) = 0
c    P'(1,x) = 1
c    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996.
c
c  Parameters:
c
c    Input, integer M, the number of evaluation points.
c
c    Input, integer N, the highest order polynomial to evaluate.
c    Note that polynomials 0 through N will be evaluated.
c
c    Input, double precision X(M), the evaluation points.
c
c    Output, double precision V(M,0:N), the values of the Legendre polynomials 
c    of order 0 through N at the points X.
c
      implicit none
      
      integer m
      integer n

      integer i
      integer j
      double precision v(m,0:n)
      double precision x(m)

      if ( n .lt. 0 ) then
        return
      end if

      do i = 1, m
        v(i,0) = 1.0D+00
      end do

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, m
        v(i,1) = x(i)
      end do
   
      do j = 2, n
        do i = 1, m

          v(i,j) = ( dble ( 2 * j - 1 ) * x(i) * v(i,j-1)   
     &             - dble (     j - 1 ) *        v(i,j-2) ) 
     &             / dble (     j     )

         end do
      end do
     
      return
      end

c     ============================================================      
      Subroutine aLegP(n,x,v)
c     Computes associated Legendre Polynomial P_n(x)
      implicit double precision (a-h,o-z)
c     v(3) = P_2(x)
c     Author: Harishkumar Sellamuthu, 2021
c     ============================================================
      dimension v(50),plg(36,36)
      call p_polynomial_value ( 1, 50, x, plg)
      k = 1
      do i = 1,50
         do j = 1,50
            v(k) = plg(j,i)
            k = k+1
         end do
      end do
      return
      end

      function r8_factorial ( n )

c*********************************************************************72
c
cc R8_FACTORIAL computes the factorial of N.
c
c  Discussion:
c
c    factorial ( N ) = product ( 1 .le. I .le. N ) I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c
c    Output, double precision R8_FACTORIAL, the factorial of N.
c
c     CHECKED - 16072021
      implicit none

      integer i
      integer n
      double precision r8_factorial

      r8_factorial = 1.0D+00

      do i = 1, n
        r8_factorial = r8_factorial * dble ( i )
      end do

      return
      end
