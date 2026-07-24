c     -------------------------------------------------------------
c     Shared constants initialisation
c     -------------------------------------------------------------
c     ============================================================
      subroutine init_constants()
c     Read physical constants from input/const_new.dat and populate
c     the /xy/ common block. Must be called once at program start.
c     ============================================================
      implicit double precision (a-h,o-z)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth

      pi  = 4.d0*datan(1.d0)
      d2r = pi/180.d0
      r2d = 180.d0/pi

      open(unit=99, file='input/const_new.dat', status='old')
      read(99,*) amue, R_Earth, AU
      close(99)

      return
      end

c     -------------------------------------------------------------
c     integrators
c     -------------------------------------------------------------
c     ks-orbital integrator
      function rkgil(n,y,f,x,h,nt)
      implicit double precision (a-h,o-z)
      dimension y(10),f(10),q(10)
      save q
      nt=nt+1 
      go to (1,2,3,4),nt
    1 do 11 j=1,n 
   11 q(j)=0.0d0
         a=0.5d0
      x=x+h/2.0d0 
      go to 5 
    2 a=0.29289321881d0 
      go to 5 
    3 a=1.7071067812d0
      x=x+h/2.0d0 
      go to 5 
    4 do 41 i=1,n 
  41  y(i)=y(i)+h*f(i)/6.0d0-q(i)/3.0d0 
      nt=0
      rkgil=2.0d0 
      go to 6 
    5 do 51 l=1,n 
      y(l)=y(l)+a*(h*f(l)-q(l)) 
  51  q(l)=2.0d0*a*h*f(l)+(1.0d0-3.0d0*a)*q(l)
      rkgil=1.0d0 
    6 continue
      return
      end

c     -------------------------------------------------------------
c     Coordinate Transformations
c     -------------------------------------------------------------

c     ============================================================
      subroutine car2oe(x,xd,a_kep)
c     cartesian state vector to keplerian orbital elements
c     Verified on 06-16-2018
c     ============================================================
      implicit double precision (a-h,o-z)
!     a_kep(3:7) in degrees;
!     x in km; xd in km/s;
      dimension a_k(3),a_n(3),ecc(3),h(3),x(3),xd(3),a_kep(7)

      common /xy/ pi,d2r,r2d,amue,AU,R_Earth

      rmag = vmn(x)
      vmag = vmn(xd)
      rv = dotp3(x,xd)

      call cross(x,xd,h) !Angular momentum
      a_k(1) = 0.d0
      a_k(2) = 0.d0
      a_k(3) = 1.d0
      call cross(a_k,h,a_n) !Node

      do i = 1,3
         ecc(i) = ((dotp3(xd,xd)-amue/rmag)*x(i)-rv*xd(i))/amue
      end do

      e = vmn(ecc) !Eccentricity
      a_kep(2) = e
      
      anmag = vmn(a_n)

      En = 0.5d0*dotp3(xd,xd)-amue/rmag

!     Semi-Major axis
      if(dabs(En).gt.0.d0) then
         a = -amue/(2.d0*En)
         p = a*(1.d0-e**2)
      else ! Exactly parabolic orbit
         a = -1.d0
         p = dotp3(h,h)
      end if   
      
      hmag = vmn(h)
      ai = dacos(dmax1(-1.d0,dmin1(1.d0,h(3)/hmag)))

      raan = dacos(dmax1(-1.d0,dmin1(1.d0,a_n(1)/anmag)))
      if(a_n(2).lt.0.d0) raan = 2.d0*pi-raan

      aop = dacos(dmax1(-1.d0,dmin1(1.d0,
     &      dotp3(a_n,ecc)/(anmag*e))))
      if(ecc(3).lt.0.d0) aop = 2.d0*pi-aop

      anu = dacos(dmax1(-1.d0,dmin1(1.d0,
     &      dotp3(ecc,x)/(e*rmag))))
      if(dotp3(x,xd).lt.0.d0) then
         anu = 2.d0*pi - anu
      end if
      
!     Special Cases
!     Elliptical, equatorial. aop = true longitude of periapsis
      if(dabs(dabs(h(3)/hmag)-1.d0).lt.1.d-10) then
         raan = 0.d0
         if(e.gt.1.d-10) then
            aop = dacos(ecc(1)/e)
            if(ecc(2).lt.0.d0) aop = 2.d0*pi - aop
         end if
      end if
!     Circular, inclined. anu = argument of latitude
      if(e.lt.1.d-10) then
         aop = 0.d0
         anu = dacos(dmax1(-1.d0,dmin1(1.d0,
     &         dotp3(a_n,x)/(anmag*rmag))))
         if(x(3).lt.0.d0) anu = 2.d0*pi - anu
         am = anu
      end if
!     Circular,equatorial. anu = true longitude
      if(dabs(dabs(h(3)/hmag)-1.d0).lt.1.d-10.and.e.lt.1.d-10) then
         raan = 0.d0
         aop = 0.d0
         anu = dacos(dmax1(-1.d0,dmin1(1.d0,x(1)/rmag)))
         if(x(2).lt.0.d0) anu = 2.d0*pi - anu
      end if

!     Mean anomaly
      ea = 0.d0
      sinn = dsin(anu)
      cosn = dcos(anu)
      if(e.lt.1.d0) then
         sinE = (sinn*dsqrt(1.d0-e**2))/(1.d0+e*cosn)
         cosE = (e+cosn)/(1.d0+e*cosn)
         ea = atan2(sinE,cosE)
         am = ea - e*sinE
      else if(e.gt.1.d0) then
         ga = 2.d0*datan(dsqrt((e-1.d0)/(e+1.d0))*dtan(0.5d0*anu))
         am = e*dtan(ga)-dlog(dtan(0.5d0*ga+0.25d0*pi))
      end if

      if(am.gt.2.d0*pi) am = am - 2.d0*pi
      if(e.lt.1.d0) then
         if(ea.gt.2.d0*pi) ea = ea - 2.d0*pi
         if(ea.lt.0.d0) ea = ea+2.d0*pi
      end if
      if(am.lt.0.d0) am = am+2.d0*pi
      
      a_kep(1) = a
      a_kep(2) = e
      a_kep(3) = ai
      a_kep(4) = raan
      a_kep(5) = aop
      a_kep(6) = am
      a_kep(7) = ea

      do i = 3,7
         a_kep(i)=a_kep(i)*r2d
      end do
      
      return
      end

c     ============================================================
      subroutine oe2car(a_kep,x,xd,tol)
c     keplerian orbital elements to cartesian state vector
c     Verified on 06-16-2018      
c     ============================================================
      implicit double precision (a-h,o-z)
!     a_kep(3:7) in degrees;
!     x in km; xd in km/s;
      dimension x(3),xd(3),a_kep(7),rot(3,3),V_pf(3),R_pf(3)

      common /xy/ pi,d2r,r2d,amue,AU,R_Earth

!     degree to radians (use local copies to preserve input)
      aa = a_kep(1)
      ae = a_kep(2)
      ai = a_kep(3)*d2r
      aO = a_kep(4)*d2r
      aw = a_kep(5)*d2r
      am = a_kep(6)*d2r
      p = aa*(1.d0-ae**2)

      sini = dsin(ai)
      cosi = dcos(ai)
      sinO = dsin(aO)
      cosO = dcos(aO)
      sinw = dsin(aw)
      cosw = dcos(aw)

      am = dmod(am,2.d0*pi)
      
      aea = akesolve(ae,am,tol)
      sinE = dsin(aea)
      cosE = dcos(aea)
      sinf = (dsqrt(1.d0-ae**2)*sinE)/(1.d0-ae*cosE)
      cosf = (cosE-ae)/(1.d0-ae*cosE)

      R_pf(1) = cosf/(1.d0+ae*cosf)
      R_pf(2) = sinf/(1.d0+ae*cosf)
      R_pf(3) = 0.d0
      do i = 1,3
         R_pf(i)=R_pf(i)*p
      end do

      V_pf(1) = -sinf
      V_pf(2) = ae + cosf
      V_pf(3) = 0.d0
      do i = 1,3
         V_pf(i)=V_pf(i)*dsqrt(amue/p)
      end do

      rot(1,1) = cosO*cosw - sinO*sinw*cosi
      rot(2,1) = sinO*cosw + cosO*sinw*cosi
      rot(3,1) = sinw*sini

      rot(1,2) = -cosO*sinw - sinO*cosw*cosi
      rot(2,2) = -sinO*sinw + cosO*cosw*cosi
      rot(3,2) = cosw*sini

      rot(1,3) = sinO*sini
      rot(2,3) = -cosO*sini
      rot(3,3) = cosi

      x(1) = rot(1,1)*R_pf(1)+rot(1,2)*R_pf(2)+rot(1,3)*R_pf(3)
      x(2) = rot(2,1)*R_pf(1)+rot(2,2)*R_pf(2)+rot(2,3)*R_pf(3)
      x(3) = rot(3,1)*R_pf(1)+rot(3,2)*R_pf(2)+rot(3,3)*R_pf(3)

      xd(1) = rot(1,1)*V_pf(1)+rot(1,2)*V_pf(2)+rot(1,3)*V_pf(3)
      xd(2) = rot(2,1)*V_pf(1)+rot(2,2)*V_pf(2)+rot(2,3)*V_pf(3)
      xd(3) = rot(3,1)*V_pf(1)+rot(3,2)*V_pf(2)+rot(3,3)*V_pf(3)
      
      return
      end
c     ============================================================
      function akesolve(ecc,am,tol)
c     Verified on 06-16-2018      
c     ============================================================      
      implicit double precision(a-h,o-z)
      sinM = dsin(am)
      if(dabs(sinM).le.0.d0.or.ecc.le.0.d0) then
         akesolve = am
         return
      end if

      ak0 = akepeq(ecc,am,am)
      am1 = am + sinM/dabs(sinM)*ecc
      ak1 = akepeq(ecc,am1,am)

      akprev = ak0
      eprev = am
      akcurr = ak1
      ecurr = am1

      it = 0

      do
         enext = (ecurr*akprev-eprev*akcurr)/(akprev-akcurr)
         aknext = akepeq(ecc,enext,am)

         if(aknext*akprev.gt.0.d0) then
            eprev = ecurr
            akprev = akcurr
         end if

         ecurr = enext
         akcurr = aknext

         it = it+1

         if(dabs(akcurr).le.tol) then
            akesolve = ecurr
            return
         end if

         if(it.gt.100) then
            akesolve = ecurr
            write(*,*) 'Kesolve did not converge'
            return
         end if
         
      end do

      end function

c     ============================================================      
      function akepeq(ecc,ea,am)
c     Verified on 06-16-2018            
c     ============================================================      
      implicit double precision(a-h,o-z)
      akepeq=ea-ecc*dsin(ea)-am
      end function

c     ============================================================      
      subroutine ks2car(u,us,x,xd,w)
c     Conversion from KS to Cartesian coordinates
!     w: (un)perturbed angular frequency of harmonic oscillator
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension u(4),us(4),x(3),xd(3)

      r=dotp4(u,u)
      
      x(1)=u(1)*u(1)-u(2)*u(2)-u(3)*u(3)+u(4)*u(4) 
      x(2)=2.0*(u(1)*u(2)-u(3)*u(4)) 
      x(3)=2.0*(u(1)*u(3)+u(2)*u(4)) 

      cst1=4.0*w/r 

      xd(1)=cst1*(u(1)*us(1)-u(2)*us(2)-u(3)*us(3)+u(4)*us(4))
      xd(2)=cst1*(u(2)*us(1)+u(1)*us(2)-u(4)*us(3)-u(3)*us(4))
      xd(3)=cst1*(u(3)*us(1)+u(4)*us(2)+u(1)*us(3)+u(2)*us(4))

      return
      end

c     ============================================================      
      subroutine car2ks(x,xd,u,us,w)
c     Conversion from Cartesian to KS coordinates
!     w: (un)perturbed angular frequency of harmonic oscillator
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension u(4),us(4),x(3),xd(3)

      r=vmn(x)
      if (x(1)) 12,11,11
 11   u(4)=0.0d0
      rpx1=r+x(1) 
      u(1)=.5d0*rpx1-u(4)*u(4)
      u(1)=dsqrt(u(1)) 
      u(2)=(x(2)*u(1)+x(3)*u(4))/rpx1 
      u(3)=(x(3)*u(1)-x(2)*u(4))/rpx1 

      go to 13

 12   u(3)=0.0d0
      rmx1=r-x(1)
      u(2)=dsqrt(.5d0*rmx1)
      u(1)=(x(2)*u(2)+x(3)*u(3))/rmx1 
      u(4)=(x(3)*u(2)-x(2)*u(3))/rmx1 

 13   cnst=.25d0/w 

      us(1)=(u(1)*xd(1)+u(2)*xd(2)+u(3)*xd(3))*cnst
      us(2)=(-u(2)*xd(1)+u(1)*xd(2)+u(4)*xd(3))*cnst
      us(3)=(-u(3)*xd(1)-u(4)*xd(2)+u(1)*xd(3))*cnst
      us(4)=(u(4)*xd(1)-u(3)*xd(2)+u(2)*xd(3))*cnst
      
      return
      end

c     ============================================================
      subroutine third_body_aux(deg, x, tb, R1, amu_tb,
     &                          Rtb1, RbyRtb, first_tb, cphi, cphi2)
c     Distance, power-series ratios, leading coefficient, and phase
c     cosine for one third body (sun or moon) at the satellite's
c     current position. Shared by driver_KS's initial setup and
c     per-step refresh so both use one consistent sign convention
c     for first_tb.
c     deg MUST be declared integer: callers pass nsun_deg/nmoon_deg
c     (integers). Under implicit typing 'd' made it double precision,
c     so the integer bit pattern read as a denormal ~1e-323 and the
c     power-series do-loop below ran ZERO times -- RbyRtb(2..) stayed
c     0.0 and the entire third-body force silently vanished from the
c     dynamics (2026-06-07 refactor until found 2026-07-12 via GMAT
c     ON/OFF sensitivity test, issue #23).
      implicit double precision (a-h,o-z)
      integer deg
      dimension x(3), tb(3), RbyRtb(*)

      Rtb1 = vmn(tb)
      RbyRtb(1) = R1/Rtb1
      first_tb = -amu_tb/Rtb1
      cphi  = dotp3(x,tb)/(R1*Rtb1)
      cphi2 = cphi*cphi

      do i = 1,deg
         RbyRtb(i+1) = RbyRtb(i)*RbyRtb(1)
      end do

      return
      end

c     ============================================================
      subroutine ks2ksr(y,u,us,E,cosE,sinE)
c     Conversion from Cartesian to KS coordinates
!     w: (un)perturbed angular frequency of harmonic oscillator
c     ============================================================
      implicit double precision(a - h, o - z)
      dimension u(4),us(4),y(10)

      cosE = dcos(E/2.d0)
      sinE = dsin(E/2.d0)

      do i = 1,4
         u(i) = y(i+1)*cosE+y(i+5)*sinE
         us(i) = (-y(i+1)*sinE+y(i+5)*cosE)*0.5d0
      end do
      
      return
      end

c     ----------------------------------------
c     Variable operations
c     ----------------------------------------
c     ============================================================
      subroutine u2uu(u,uu)
c     u(i): KS variables
!     uu(i): KS variables interchanged as u(i+2)
c     ============================================================
      implicit double precision(a - h, o - z)
      dimension u(4),uu(4)

      uu(1) = u(3)
      uu(2) = u(4)
      uu(3) = u(1)
      uu(4) = u(2)
      
      return
      end      

c     ============================================================      
      subroutine u2qu(u,qu,i)
c     u(i): KS variables
!     qu(i): KS variables interchanged for lunisolar potential
c     ============================================================
      implicit double precision(a - h, o - z)
      dimension u(4),qu(3)

      if(i.eq.1) then
      
         qu(1) = u(1)
         qu(2) = u(2)
         qu(3) = u(3)

      else if(i.eq.2) then
         
         qu(1) = -u(2)
         qu(2) = u(1)
         qu(3) = u(4)

      else if(i.eq.3) then
         
         qu(1) = -u(3)
         qu(2) = -u(4)
         qu(3) = u(1)
         
      else if(i.eq.4) then
                  
         qu(1) = u(4)
         qu(2) = -u(3)
         qu(3) = u(2)

      end if
         
      return
      end      
      
      
c     ----------------------------------------
c     Matrix and Vector operations
c     ----------------------------------------

c     ======================================
      function dotp3(x,y)
!     Dot product of 3-dimensional vector
c     ======================================
      implicit double precision (a-h,o-z)
      dimension x(3),y(3)
      dotp3=0.d0
      do i = 1,3
         dotp3=x(i)*y(i)+dotp3
      end do
      return
      end

c     ======================================
      function dotp4(x,y)
!     Dot product of 4-dimensional vector      
c     ======================================
      implicit double precision (a-h,o-z)
      dimension x(4),y(4)
      dotp4=0.d0
      do i = 1,4
         dotp4=x(i)*y(i)+dotp4
      end do
      return
      end

c     ===============================================
      function vmn(x)
c     magnitude of a vector - array,dot product based
!     3-dimension
c     ===============================================
      implicit double precision (a-h,o-z)
      dimension x(3)
      vmn=dsqrt(dotp3(x,x))
      return
      end

c     ===============================================      
      subroutine cross(x,y,z)
!     Cross product (3-dimensional)
!     Input: x,y; Output:z
!     z = x X y
c     ===============================================      
      implicit double precision (a-h,o-z)
      dimension x(3),y(3),z(3)
      z(1)=x(2)*y(3)-y(2)*x(3)
      z(2)=y(1)*x(3)-x(1)*y(3)
      z(3)=x(1)*y(2)-y(1)*x(2)
      return
      end
      
c     ============================================================
      subroutine cal2jd(cal,ddj)
      implicit double precision(a-h,o-z)
      dimension cal(6)
      iyr = cal(1)
      imon = cal(2)
      iday = cal(3)
      intmon1=(imon+9.0d0)/12.0d0
      intmon2=(275.0d0*imon)/9.0d0
      jyr=(7.0d0*(iyr+intmon1))/4.0d0
      tyr1=367.0d0*iyr
      ftime=(cal(4)*3600.d0 + cal(5)*60.d0 + cal(6)) / 86400.d0
      ccday1=iday+1721013.5d0
      ddj=tyr1-jyr+intmon2+ccday1+ftime
      return
      end


      function atan3(a,b)
      implicit double precision (a-h,o-z)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
      atan3 = datan2(a, b)
      if(atan3.lt.0.d0) atan3 = atan3 + 2.d0*pi
      return
      end
c     ============================================================
      
c     ============================================================
      subroutine INTPOL(XT,YT,M1,X1,Y1)
c     Linear interpolation in a sorted table XT(M1)/YT(M1)
c     Returns Y1 at X1; stops if X1 is out of table range.
c     ============================================================
      implicit double precision(a-h,o-z)
      dimension XT(M1),YT(M1)
      if(X1.lt.XT(1).or.X1.gt.XT(M1)) then
         write(*,*) 'INTPOL: altitude out of table range',X1
         Y1 = 0.d0
         return
      end if
      do 11 I=1,M1
         if(X1-XT(I)) 11,13,10
 10      if(X1-XT(I+1)) 15,11,11
 15      Y1=YT(I)+(YT(I+1)-YT(I))*(X1-XT(I))/(XT(I+1)-XT(I))
         go to 12
 13      Y1=YT(I)
 11   continue
 12   return
      end

      Function aLegendreP(l,m,x)
c     Computes associated Legendre Polynomial Plm(x)
c     0<=m<=l while x lies in the range of -1 to 1
      implicit double precision (a-h,o-z)
c     ============================================================
      pmm = 1.d0
      if(m.gt.0) then
         somx2 = dsqrt((1.d0-x)*(1.d0+x))
         fact = 1.d0
         do i = 1,m
            pmm = pmm*fact*somx2
            fact = fact+2
         end do
      end if
      if(l.eq.m) then
         aLegendreP = pmm
      else
         pmmp1 = x*(2*m+1)*pmm
         if(l.eq.m+1) then
            aLegendreP = pmmp1
         else
            do j = m+2,l
               pll = (x*(2.d0*j-1)*pmmp1-(j+m-1)*pmm)/(j-m)
               pmm = pmmp1
               pmmp1 = pll
            end do
            aLegendreP = pll
         end if
      end if
      return
      end

c     ----------------------------------------
c     Force models, Constant, Coefficients
c     ----------------------------------------
c     ============================================================
      subroutine geo_coeff(n,c_j)
c     Reads zonal harmonic coefficients Jn from EGM2008_to2190_TideFree.
c     Streams row-by-row and stops as soon as the n-1 required zonal
c     (m=0) terms for degrees 2..n have been collected.
c
c     Performance: reads only O(n^2) lines instead of all 2 401 333.
c     For n=2 (J2 only): 3 lines.  For n=10: 63 lines.  For n=100: 5148.
c
c     EGM2008 file format per line:
c       degree  order  Cnm  Snm  sigma_Cnm  sigma_Snm
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension c_j(n)
      integer ideg, iord

      open (unit = 7, file = 'input/EGM2008_to2190_TideFree')

      j = 2
 20   continue
         if (j .gt. n) goto 10
         read(7,*,end=10) ideg, iord, cnm, sval, sg1, sg2
         if (iord .eq. 0) then
            c_j(j) = -dsqrt(2.d0*dble(ideg)+1.d0)*cnm
            j = j + 1
         end if
      goto 20
 10   continue

      close(7)
      return
      end

c     ============================================================      
      subroutine solarnpv(dj,s)
c     Compute sun position vector in inertial geo-centric frame
c     (EME2000). Low-precision analytic series per Montenbruck &
c     Gill, "Satellite Orbits" (2000) sect. 3.3.2, precessed from
c     equinox-of-date to J2000. Includes the elliptic distance
c     variation -- the pre-2026-07-12 version used a FIXED 1 AU
c     distance, a +/-1.67 percent seasonal error that dominated
c     the ephemeris gap vs DE405.
      implicit double precision(a-h,o-z,m)
      dimension s(3)
c     dj = Julian day
c     s = sun position vector (x_s, y_s, z_s) in km
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
c     ============================================================

c     julian centuries since J2000
      tc=(dj-2451545.0d0)/36525.0d0

c     mean anomaly of the sun (deg)
      am=dmod(357.5256d0+35999.049d0*tc,360.0d0)
      amr=am*d2r

c     ecliptic longitude (deg). NOTE: unlike the lunar series, this
c     M&G solar series is already J2000-referenced -- adding the
c     -1.3972*T precession term here produces a constant ~0.4
c     percent offset vs DE405 (verified empirically 2026-07-12).
      tlon=282.9400d0+am+(6892.0d0*dsin(amr)
     &     +72.0d0*dsin(2.0d0*amr))/3600.0d0
      tlonr=dmod(tlon,360.0d0)*d2r

c     geocentric distance (km)
      rsun=(149.619d0-2.499d0*dcos(amr)
     &     -0.021d0*dcos(2.0d0*amr))*1.0d6

c     mean obliquity of the ecliptic (J2000)
      tobl=23.43929111d0*d2r

      s(1) = rsun*dcos(tlonr)
      s(2) = rsun*dsin(tlonr)*dcos(tobl)
      s(3) = rsun*dsin(tlonr)*dsin(tobl)

      return
      end

c     ============================================================
      subroutine sun_azimuth(ai, omega, raan, alpha_s, delta_s,
     &                       lambda_s)
c     Compute Sun azimuth angle with respect to the spacecraft
c     orbital plane. Cook (1962), thesis Appendix 2.
c
c     Input (all in radians)
c       ai      : inclination
c       omega   : argument of perigee
c       raan    : right ascension of the ascending node
c       alpha_s : Sun right ascension
c       delta_s : Sun declination
c
c     Output (radians, mapped to [0, 2*pi))
c       lambda_s: Sun azimuth angle
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision lambda_s

      pi = 4.d0*datan(1.d0)

      csw = dcos(omega)
      snw = dsin(omega)
      csi = dcos(ai)
      sni = dsin(ai)
      csO = dcos(raan)
      snO = dsin(raan)

c     Direction cosines of apsidal line (perigee) in ECI
      el11 = csw*csO - snO*snw*csi
      el12 = csw*snO + csO*snw*csi
      el13 = snw*sni

c     Direction cosines of normal to apsidal line in orbital plane
      el21 = -snw*csO - snO*csw*csi
      el22 = -snw*snO + csO*csw*csi
      el23 = csw*sni

c     Sun unit vector in ECI
      csA = dcos(alpha_s)
      snA = dsin(alpha_s)
      csD = dcos(delta_s)
      snD = dsin(delta_s)
      els1 = csA*csD
      els2 = snA*csD
      els3 = snD

c     Projections onto orbital plane
      el3 = el11*els1 + el12*els2 + el13*els3
      el4 = el21*els1 + el22*els2 + el23*els3

c     Sun azimuth angle (0 to 2*pi)
      lambda_s = datan2(el3, el4)
      if (lambda_s .lt. 0.d0) lambda_s = lambda_s + 2.d0*pi

      return
      end

c     ============================================================
      subroutine force_models(n_for,ngeo,s,amoon)
c     Selection of force models
      implicit double precision(a-h,o-z)
      double precision amoon
      dimension n_for(3)
c     n_for(1): oblateness, n_for(2): sun, n_for(3): moon
c     ============================================================
      if(n_for(1).eq.0) ngeo = 0.d0
      if(n_for(2).eq.0) s = 0.d0
      if(n_for(3).eq.0) amoon = 0.d0

      return
      end

c     ============================================================
      function shadfncyl(x1,x2,x3,xs,ys,zs)
c     Cylindrical shadow function.
c     Returns ROHP: 1 = full sunlight, 0 = eclipse (smooth tanh).
c     Ref: Hubaux et al. (2012).  All positions in km.
c     ============================================================
      implicit double precision(a-h,o-z)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
      gam    = 1.d9
      rs     = dsqrt(xs*xs+ys*ys+zs*zs)
      r      = dsqrt(x1*x1+x2*x2+x3*x3)
      cosphi = (x1*xs+x2*ys+x3*zs)/rs
      scr    = cosphi + dsqrt(r*r - R_Earth*R_Earth)
      shadfncyl = 0.5d0*(1.d0 + dtanh(gam*scr))
      return
      end

c     ============================================================
      function shadfncone(x1,x2,x3,xs,ys,zs)
c     Conical shadow function.
c     Returns ROHP: 1 = full sunlight, 0 = umbra, (0,1) = penumbra.
c     Ref: D.G. Cook (2001) Masters Thesis.  All positions in km.
c     ============================================================
      implicit double precision(a-h,o-z)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
      rds  = 696342.d0          ! solar radius (km)
      xp   = xs - x1
      yp   = ys - x2
      zp   = zs - x3
      rs   = dsqrt(xs*xs+ys*ys+zs*zs)
      r    = dsqrt(x1*x1+x2*x2+x3*x3)
      p    = dsqrt(xp*xp+yp*yp+zp*zp)
      et   = dacos(-(x1*xp+x2*yp+x3*zp)/(r*p))
      et2  = et*et
      rhoe  = dasin(R_Earth/r)
      rhos  = dasin(rds/p)
      rhos2 = rhos*rhos
      rhoe2 = rhoe*rhoe
      rh1   = dabs(rhos2 - rhoe2)

      if (et .ge. (rhos+rhoe)) then
         shadfncone = 1.d0
         return
      end if
      if (et .le. (rhoe-rhos)) then
         shadfncone = 0.d0
         return
      end if
c     penumbra
      al = 0.5d0*(rhos+rhoe+et)
      q  = 2.d0*dsqrt(al*(al-et)*(al-rhos)*(al-rhoe))/et
      if (rh1 .le. et2) then
         t = rhos2*dasin(q/rhos)+rhoe2*dasin(q/rhoe)-q*et
      else
         t = rhos2*(pi-dasin(q/rhos))+rhoe2*dasin(q/rhoe)-q*et
      end if
      shadfncone = 1.d0 - t/(pi*rhos2)
      return
      end

      subroutine lunarpv(dj,tm)
c     Compute moon position vector in inertial geo-centric frame
c     (EME2000). Truncated Brown lunar theory per Montenbruck &
c     Gill, "Satellite Orbits" (2000) sect. 3.3.2, precessed from
c     equinox-of-date to J2000. Replaces the 7-term-per-coordinate
c     Cartesian series (1999 GSFC ISSFD paper) which differed from
c     DE405 by up to ~3.6 percent -- the dominant full-force-model
c     error source vs GMAT (issue #23).
      implicit double precision(a-h,o-z)
      dimension tm(3)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth

c     time in julian centuries since J2000
      tc=(dj-2451545.0d0)/36525.0d0

c     mean elements of the lunar orbit (deg)
      al0=dmod(218.31617d0+481267.88088d0*tc,360.0d0)
      al =dmod(134.96292d0+477198.86753d0*tc,360.0d0)
      alp=dmod(357.52543d0+ 35999.04944d0*tc,360.0d0)
      af =dmod( 93.27283d0+483202.01873d0*tc,360.0d0)
      ad =dmod(297.85027d0+445267.11135d0*tc,360.0d0)
      alr=al*d2r
      alpr=alp*d2r
      afr=af*d2r
      adr=ad*d2r

c     ecliptic longitude perturbation (arcsec)
      dlam=22640.0d0*dsin(alr)+769.0d0*dsin(2.0d0*alr)
     &     -4586.0d0*dsin(alr-2.0d0*adr)+2370.0d0*dsin(2.0d0*adr)
     &     -668.0d0*dsin(alpr)-412.0d0*dsin(2.0d0*afr)
     &     -212.0d0*dsin(2.0d0*alr-2.0d0*adr)
     &     -206.0d0*dsin(alr+alpr-2.0d0*adr)
     &     +192.0d0*dsin(alr+2.0d0*adr)
     &     -165.0d0*dsin(alpr-2.0d0*adr)
     &     +148.0d0*dsin(alr-alpr)-125.0d0*dsin(adr)
     &     -110.0d0*dsin(alr+alpr)
     &     -55.0d0*dsin(2.0d0*afr-2.0d0*adr)

c     ecliptic longitude of date (deg), precessed to J2000
      tlam=al0+dlam/3600.0d0
      tlam=tlam-1.3972d0*tc
      tlamr=dmod(tlam,360.0d0)*d2r

c     ecliptic latitude (arcsec -> rad)
      ah=(af+(dlam+412.0d0*dsin(2.0d0*afr)
     &   +541.0d0*dsin(alpr))/3600.0d0)*d2r
      beta=18520.0d0*dsin(ah)-526.0d0*dsin(afr-2.0d0*adr)
     &     +44.0d0*dsin(alr+afr-2.0d0*adr)
     &     -31.0d0*dsin(-alr+afr-2.0d0*adr)
     &     -25.0d0*dsin(-2.0d0*alr+afr)
     &     -23.0d0*dsin(alpr+afr-2.0d0*adr)
     &     +21.0d0*dsin(-alr+afr)
     &     +11.0d0*dsin(-alpr+afr-2.0d0*adr)
      betar=beta/3600.0d0*d2r

c     geocentric distance (km)
      rmoon=385000.0d0-20905.0d0*dcos(alr)
     &     -3699.0d0*dcos(2.0d0*adr-alr)-2956.0d0*dcos(2.0d0*adr)
     &     -570.0d0*dcos(2.0d0*alr)
     &     +246.0d0*dcos(2.0d0*alr-2.0d0*adr)
     &     -205.0d0*dcos(alpr-2.0d0*adr)
     &     -171.0d0*dcos(alr+2.0d0*adr)
     &     -152.0d0*dcos(alr+alpr-2.0d0*adr)

c     ecliptic -> equatorial (mean obliquity J2000)
      tobl=23.43929111d0*d2r
      xe=rmoon*dcos(tlamr)*dcos(betar)
      ye=rmoon*dsin(tlamr)*dcos(betar)
      ze=rmoon*dsin(betar)
      tm(1)=xe
      tm(2)=ye*dcos(tobl)-ze*dsin(tobl)
      tm(3)=ye*dsin(tobl)+ze*dcos(tobl)

      return
      end

c     ============================================================
      subroutine parse_epoch(estr, cal)
c     Parse CCSDS epoch string YYYY-MM-DDTHH:MM:SS.sss into cal(6)
c     cal = [year, month, day, hour, minute, second]
c     ============================================================
      implicit double precision (a-h, o-z)
      character*(*) estr
      dimension cal(6)
      integer iyr, imon, iday, ih, imin

      read(estr(1:4),   '(I4)') iyr
      read(estr(6:7),   '(I2)') imon
      read(estr(9:10),  '(I2)') iday
      read(estr(12:13), '(I2)') ih
      read(estr(15:16), '(I2)') imin
      read(estr(18:),   *     ) sec_val

      cal(1) = dble(iyr)
      cal(2) = dble(imon)
      cal(3) = dble(iday)
      cal(4) = dble(ih)
      cal(5) = dble(imin)
      cal(6) = sec_val

      return
      end

c     ============================================================
      subroutine read_opm(iunit, x, xd, cal)
c     Parse a CCSDS OPM v2.0 ASCII file.
c     Scans every line for the seven STATE_VECTOR keywords and
c     extracts: EPOCH, X, Y, Z, X_DOT, Y_DOT, Z_DOT.
c     Lines without '=' and all other keywords are silently skipped.
c     Unit conversion brackets [km] / [km/s] are stripped before read.
c     ============================================================
      implicit double precision (a-h, o-z)
      dimension x(3), xd(3), cal(6)
      character*200 line
      character*80  valstr
      integer       nfound, ieq, ibr

      nfound = 0

 10   continue
         if (nfound .ge. 7) goto 99
         read(iunit, '(A)', end=99) line
         line   = adjustl(line)
         ieq    = index(line, '=')
         if (ieq .eq. 0) goto 10

c        Extract value substring (strip leading blanks and unit bracket)
         valstr = adjustl(line(ieq+1:))
         ibr    = index(valstr, '[')
         if (ibr .gt. 1) valstr = valstr(1:ibr-1)

c        Match keyword (check _DOT variants before single-letter X/Y/Z)
         if     (line(1:5) .eq. 'EPOCH') then
            call parse_epoch(adjustl(valstr), cal)
            nfound = nfound + 1
         else if (line(1:5) .eq. 'X_DOT') then
            read(valstr, *) xd(1)
            nfound = nfound + 1
         else if (line(1:5) .eq. 'Y_DOT') then
            read(valstr, *) xd(2)
            nfound = nfound + 1
         else if (line(1:5) .eq. 'Z_DOT') then
            read(valstr, *) xd(3)
            nfound = nfound + 1
         else if (line(1:1).eq.'X' .and. line(2:2).ne.'_') then
            read(valstr, *) x(1)
            nfound = nfound + 1
         else if (line(1:1).eq.'Y' .and. line(2:2).ne.'_') then
            read(valstr, *) x(2)
            nfound = nfound + 1
         else if (line(1:1).eq.'Z' .and. line(2:2).ne.'_') then
            read(valstr, *) x(3)
            nfound = nfound + 1
         end if
      goto 10

 99   return
      end

c     ============================================================
      subroutine write_opm(iunit, epochstr, x, xd, pek)
c     Write a CCSDS OPM v2.0 ASCII file.
c     Emits the header/metadata, the STATE_VECTOR block (EPOCH, X, Y,
c     Z, X_DOT, Y_DOT, Z_DOT) and a KEPLERIAN_ELEMENTS block derived
c     from the osculating elements pek (from car2oe).
c     ============================================================
      implicit double precision (a-h, o-z)
      character*23 epochstr
      dimension x(3), xd(3), pek(7)

      write(iunit,'(A)')     'CCSDS_OPM_VERS = 2.0'
      write(iunit,'(A,A)')   'CREATION_DATE  = ', epochstr
      write(iunit,'(A)')     'ORIGINATOR     = KSROP'
      write(iunit,'(A)')     ''
      write(iunit,'(A)')     'META_START'
      write(iunit,'(A)')     'OBJECT_NAME    = SATELLITE'
      write(iunit,'(A)')     'OBJECT_ID      = UNKNOWN'
      write(iunit,'(A)')     'CENTER_NAME    = EARTH'
      write(iunit,'(A)')     'REF_FRAME      = EME2000'
      write(iunit,'(A)')     'TIME_SYSTEM    = UTC'
      write(iunit,'(A)')     'META_STOP'
      write(iunit,'(A)')     ''
      write(iunit,'(A)')     'STATE_VECTOR'
      write(iunit,'(A,A)')   'EPOCH          = ', epochstr
      write(iunit,'(A,F16.6,A)')  'X              = ', x(1),  ' [km]'
      write(iunit,'(A,F16.6,A)')  'Y              = ', x(2),  ' [km]'
      write(iunit,'(A,F16.6,A)')  'Z              = ', x(3),  ' [km]'
      write(iunit,'(A,F16.9,A)')  'X_DOT          = ', xd(1), ' [km/s]'
      write(iunit,'(A,F16.9,A)')  'Y_DOT          = ', xd(2), ' [km/s]'
      write(iunit,'(A,F16.9,A)')  'Z_DOT          = ', xd(3), ' [km/s]'
      write(iunit,'(A)')     ''
      write(iunit,'(A)')     'KEPLERIAN_ELEMENTS'
      write(iunit,'(A,F16.6,A)')  'SEMI_MAJOR_AXIS   = ', pek(1),
     &                            ' [km]'
      write(iunit,'(A,F16.9)')    'ECCENTRICITY      = ', pek(2)
      write(iunit,'(A,F16.6,A)')  'INCLINATION       = ', pek(3),
     &                            ' [deg]'
      write(iunit,'(A,F16.6,A)')  'RA_OF_ASC_NODE    = ', pek(4),
     &                            ' [deg]'
      write(iunit,'(A,F16.6,A)')  'ARG_OF_PERICENTER = ', pek(5),
     &                            ' [deg]'
      write(iunit,'(A,F16.6,A)')  'MEAN_ANOMALY      = ', pek(6),
     &                            ' [deg]'

      return
      end

c     ============================================================
      subroutine utc_now_epoch(epochstr, compact)
c     Return current UTC wall-clock time via Fortran DATE_AND_TIME.
c     epochstr (23) = YYYY-MM-DDTHH:MM:SS.mmm  (CCSDS epoch string)
c     compact  (15) = YYYYMMDDTHHMMSS           (filename-safe)
c     ============================================================
      implicit double precision (a-h, o-z)
      character*23 epochstr
      character*15 compact
      character*8  datstr
      character*10 timstr
      character*5  zonstr
      integer vals(8), utc_min, uyr, umon, uday, uh, umi, usc, ums

      call date_and_time(datstr, timstr, zonstr, vals)
c     vals: (1)year (2)month (3)day (4)UTC_offset_minutes
c           (5)hour (6)minute (7)second (8)millisecond

c     Convert local time to UTC by subtracting the offset
      utc_min = vals(5)*60 + vals(6) - vals(4)
      uday = vals(3)
      umon = vals(2)
      uyr  = vals(1)

      if (utc_min .lt. 0) then
         utc_min = utc_min + 1440
         uday = uday - 1
         if (uday .eq. 0) then
            umon = umon - 1
            if (umon .eq. 0) then
               umon = 12
               uyr  = uyr - 1
            end if
            uday = 31
         end if
      else if (utc_min .ge. 1440) then
         utc_min = utc_min - 1440
         uday = uday + 1
      end if

      uh  = utc_min / 60
      umi = utc_min - uh*60
      usc = vals(7)
      ums = vals(8)

      write(epochstr,
     &      '(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)')
     &      uyr,'-',umon,'-',uday,'T',uh,':',umi,':',usc,'.',ums

      write(compact,'(I4.4,2I2.2,A1,3I2.2)')
     &      uyr,umon,uday,'T',uh,umi,usc

      return
      end

c     ============================================================
      subroutine jd2epoch(djd, epochstr)
c     Convert Julian date to CCSDS epoch string: YYYY-MM-DDTHH:MM:SS.sss
c     Meeus algorithm for JD -> calendar date.
c     ============================================================
      implicit double precision (a-h, o-z)
      character*23 epochstr
      integer iyr, imon, iday, ih, imin, isec, ims
      integer jd_int, ja, jb, jc, jd2, je, alpha, total_ms

      djd5   = djd + 0.5d0
      jd_int = int(djd5)
      frac   = djd5 - dble(jd_int)

      if (jd_int .lt. 2299161) then
         ja = jd_int
      else
         alpha = int((dble(jd_int) - 1867216.25d0) / 36524.25d0)
         ja    = jd_int + 1 + alpha - alpha/4
      end if

      jb   = ja + 1524
      jc   = int((dble(jb) - 122.1d0) / 365.25d0)
      jd2  = int(365.25d0 * dble(jc))
      je   = int(dble(jb - jd2) / 30.6001d0)

      iday = jb - jd2 - int(30.6001d0 * dble(je))
      if (je .lt. 14) then
         imon = je - 1
      else
         imon = je - 13
      end if
      if (imon .gt. 2) then
         iyr = jc - 4716
      else
         iyr = jc - 4715
      end if

c     Decompose fractional day to H:M:S.mmm using integer milliseconds
      total_ms = min(nint(frac * 86400000.d0), 86399999)
      ih       = total_ms / 3600000
      total_ms = total_ms - ih * 3600000
      imin     = total_ms / 60000
      total_ms = total_ms - imin * 60000
      isec     = total_ms / 1000
      ims      = total_ms - isec * 1000

      write(epochstr,
     &      '(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)')
     &      iyr,'-',imon,'-',iday,'T',ih,':',imin,':',isec,'.',ims

      return
      end

c     ============================================================
      subroutine read_oem(iunit, maxpts, traj_jd, traj_x, traj_xd,
     &                    npts)
c     Parse a CCSDS OEM v2.0 ASCII file.
c     Reads the ephemeris block between DATA_START and DATA_STOP:
c     each data line holds an epoch string (cols 1-23) followed by
c     the position and velocity components X,Y,Z,X_DOT,Y_DOT,Z_DOT.
c     Epochs are converted to Julian dates via parse_epoch/cal2jd.
c     Returns at most maxpts points; npts is the number actually read.
c     ============================================================
      implicit double precision (a-h, o-z)
      dimension traj_jd(maxpts), traj_x(3,maxpts), traj_xd(3,maxpts)
      dimension cal(6)
      character*200 line
      logical       in_data

      npts = 0
      in_data = .false.

 10   continue
         read(iunit, '(A)', end=99) line
         line = adjustl(line)

         if (.not. in_data) then
            if (line(1:10) .eq. 'DATA_START') in_data = .true.
            goto 10
         end if

         if (line(1:9) .eq. 'DATA_STOP') goto 99
         if (len_trim(line) .eq. 0 .or. npts .ge. maxpts) goto 10

         npts = npts + 1
         call parse_epoch(line(1:23), cal)
         call cal2jd(cal, traj_jd(npts))
         read(line(24:), *) traj_x(1,npts),  traj_x(2,npts),
     &                      traj_x(3,npts),
     &                      traj_xd(1,npts), traj_xd(2,npts),
     &                      traj_xd(3,npts)
      goto 10

 99   return
      end

c     ============================================================
      subroutine write_oem(iunit, creation_str, start_str, stop_str,
     &                     traj_jd, traj_x, traj_xd, npts)
c     Write a CCSDS OEM v2.0 ASCII file.
c     Emits the header/metadata (CREATION_DATE, START_TIME, STOP_TIME)
c     and the DATA_START/DATA_STOP ephemeris block: one line per
c     buffered point with epoch (from traj_jd via jd2epoch) and the
c     position/velocity components traj_x, traj_xd.
c     ============================================================
      implicit double precision (a-h, o-z)
      character*23 creation_str, start_str, stop_str, epochstr
      dimension traj_jd(npts), traj_x(3,npts), traj_xd(3,npts)

      write(iunit,'(A)')     'CCSDS_OEM_VERS = 2.0'
      write(iunit,'(A,A)')   'CREATION_DATE  = ', creation_str
      write(iunit,'(A)')     'ORIGINATOR     = KSROP'
      write(iunit,'(A)')     ''
      write(iunit,'(A)')     'META_START'
      write(iunit,'(A)')     'OBJECT_NAME    = SATELLITE'
      write(iunit,'(A)')     'OBJECT_ID      = UNKNOWN'
      write(iunit,'(A)')     'CENTER_NAME    = EARTH'
      write(iunit,'(A)')     'REF_FRAME      = EME2000'
      write(iunit,'(A)')     'TIME_SYSTEM    = UTC'
      write(iunit,'(A,A)')   'START_TIME     = ', start_str
      write(iunit,'(A,A)')   'STOP_TIME      = ', stop_str
      write(iunit,'(A)')     'META_STOP'
      write(iunit,'(A)')     ''
      write(iunit,'(A)')     'DATA_START'
      do i = 1, npts
         call jd2epoch(traj_jd(i), epochstr)
         write(iunit,103) epochstr,
     &                traj_x(1,i),  traj_x(2,i),  traj_x(3,i),
     &                traj_xd(1,i), traj_xd(2,i), traj_xd(3,i)
      end do
      write(iunit,'(A)') 'DATA_STOP'

 103  format(A23,3(1x,F16.6),3(1x,F14.9))  !OEM data line

      return
      end

c     ============================================================
      subroutine read_cdm(iunit, tca_cal, miss_dist, pc,
     &                    x1, xd1, cov1, x2, xd2, cov2)
c     Parse a CCSDS CDM (Conjunction Data Message) KVN v1.0 file.
c     Extracts the relative-geometry summary -- TCA epoch (-> cal(6)
c     via parse_epoch), MISS_DISTANCE and COLLISION_PROBABILITY (read
c     as raw numeric values, no unit conversion -- same convention as
c     read_opm; check the message's own [unit] brackets) -- and, for
c     each OBJECT1 / OBJECT2 block, the state vector
c     (X,Y,Z,X_DOT,Y_DOT,Z_DOT, in the block's REF_FRAME) and the
c     lower-triangular 6x6 RTN-frame covariance (CR_R ... CNDOT_NDOT),
c     expanded into a symmetric matrix ordered (R,T,N,RDOT,TDOT,NDOT).
c     Lines without '=' and unrecognized keywords are silently skipped.
c     miss_dist/pc are returned as -1 if absent from the message.
c     ============================================================
      implicit double precision (a-h, o-z)
      dimension tca_cal(6), x1(3), xd1(3), cov1(6,6)
      dimension x2(3), xd2(3), cov2(6,6)
      character*200 line
      character*80  key, valstr
      integer       ieq, ibr, iobj, i, j
c     'm' falls outside the (a-h, o-z) implicit-double ranges above
c     (it is implicitly INTEGER by default), so miss_dist must be
c     declared explicitly to avoid a type mismatch with callers/assert.
      double precision miss_dist

      pc        = -1.0d0
      miss_dist = -1.0d0
      iobj      = 0
      do i = 1, 6
         do j = 1, 6
            cov1(i,j) = 0.0d0
            cov2(i,j) = 0.0d0
         end do
      end do

 10   continue
         read(iunit, '(A)', end=99) line
         line = adjustl(line)
         ieq  = index(line, '=')
         if (ieq .eq. 0) goto 10

         key    = adjustl(line(1:ieq-1))
         valstr = adjustl(line(ieq+1:))
         ibr    = index(valstr, '[')
         if (ibr .gt. 1) valstr = adjustl(valstr(1:ibr-1))

c        Track which OBJECTn block we are in. 'OBJECT' is matched
c        exactly so OBJECT_DESIGNATOR / OBJECT_NAME / OBJECT_ID etc
c        (which carry free-text values, not OBJECT1/OBJECT2) are left
c        for the default skip below.
         if (key .eq. 'OBJECT') then
            if (valstr(1:7) .eq. 'OBJECT1') then
               iobj = 1
            else if (valstr(1:7) .eq. 'OBJECT2') then
               iobj = 2
            end if
         else if (key .eq. 'TCA') then
            call parse_epoch(valstr, tca_cal)
         else if (key .eq. 'MISS_DISTANCE') then
            read(valstr, *) miss_dist
         else if (key .eq. 'COLLISION_PROBABILITY') then
            read(valstr, *) pc
         else if (iobj .eq. 1) then
            call cdm_object_field(key, valstr, x1, xd1, cov1)
         else if (iobj .eq. 2) then
            call cdm_object_field(key, valstr, x2, xd2, cov2)
         end if
      goto 10

 99   return
      end

c     ============================================================
      subroutine cdm_object_field(key, valstr, x, xd, cov)
c     Helper for read_cdm: match a single OBJECTn keyword/value pair
c     against the state-vector and RTN-covariance keywords and store
c     the parsed value into x, xd or the symmetric 6x6 covariance
c     matrix cov (row/column order R, T, N, RDOT, TDOT, NDOT).
c     Unrecognized keywords (object metadata strings, etc.) are
c     silently ignored.
c     ============================================================
      implicit double precision (a-h, o-z)
      character*(*) key, valstr
      dimension x(3), xd(3), cov(6,6)
      integer ir, ic

      if (key .eq. 'X') then
         read(valstr, *) x(1)
      else if (key .eq. 'Y') then
         read(valstr, *) x(2)
      else if (key .eq. 'Z') then
         read(valstr, *) x(3)
      else if (key .eq. 'X_DOT') then
         read(valstr, *) xd(1)
      else if (key .eq. 'Y_DOT') then
         read(valstr, *) xd(2)
      else if (key .eq. 'Z_DOT') then
         read(valstr, *) xd(3)
      else
         call cdm_cov_index(key, ir, ic)
         if (ir .gt. 0) then
            read(valstr, *) covval
            cov(ir,ic) = covval
            cov(ic,ir) = covval
         end if
      end if

      return
      end

c     ============================================================
      subroutine cdm_cov_index(key, ir, ic)
c     Map a CDM lower-triangular RTN-covariance keyword (e.g. 'CR_R',
c     'CNDOT_TDOT') to its (row,col) position ir,ic in a symmetric
c     6x6 matrix ordered (R, T, N, RDOT, TDOT, NDOT). Returns ir = 0
c     if the keyword is not a recognized covariance element (e.g.
c     'CATALOG_NAME', 'COVARIANCE_METHOD').
c     ============================================================
      implicit double precision (a-h, o-z)
      character*(*) key
      integer ir, ic, iu, irow, icol
      integer cdm_rtn_index
      external cdm_rtn_index

      ir = 0
      ic = 0
      if (key(1:1) .ne. 'C') return
      iu = index(key, '_')
      if (iu .le. 1) return

      irow = cdm_rtn_index(key(2:iu-1))
      icol = cdm_rtn_index(key(iu+1:len_trim(key)))
      if (irow .eq. 0 .or. icol .eq. 0) return

      ir = irow
      ic = icol

      return
      end

c     ============================================================
      integer function cdm_rtn_index(lbl)
c     Map an RTN-frame axis label (R, T, N, RDOT, TDOT, NDOT) to its
c     position 1..6 in the (R,T,N,RDOT,TDOT,NDOT) ordering used for
c     the CDM covariance matrix. Returns 0 if unrecognized.
c     ============================================================
      character*(*) lbl

      if (lbl .eq. 'R') then
         cdm_rtn_index = 1
      else if (lbl .eq. 'T') then
         cdm_rtn_index = 2
      else if (lbl .eq. 'N') then
         cdm_rtn_index = 3
      else if (lbl .eq. 'RDOT') then
         cdm_rtn_index = 4
      else if (lbl .eq. 'TDOT') then
         cdm_rtn_index = 5
      else if (lbl .eq. 'NDOT') then
         cdm_rtn_index = 6
      else
         cdm_rtn_index = 0
      end if

      return
      end

c     ============================================================
      subroutine write_cdm(iunit, creation_str, tca_str,
     &                     miss_dist, rel_speed, pc,
     &                     name1, x1, xd1, cov1,
     &                     name2, x2, xd2, cov2)
c     Write a CCSDS CDM (Conjunction Data Message) v1.0 KVN file.
c     Emits the header / relative-geometry block (CREATION_DATE, TCA,
c     MISS_DISTANCE, RELATIVE_SPEED, COLLISION_PROBABILITY -- omitted
c     if pc < 0) followed by one OBJECTn block per object via
c     write_cdm_object: state vector (X,Y,Z,X_DOT,Y_DOT,Z_DOT, EME2000)
c     and the lower-triangular 6x6 RTN covariance (CR_R ... CNDOT_NDOT)
c     taken from the symmetric matrix cov, row/col order (R,T,N,
c     RDOT,TDOT,NDOT).
c     ============================================================
      implicit double precision (a-h, o-z)
      character*23  creation_str, tca_str
      character*(*) name1, name2
      dimension x1(3), xd1(3), cov1(6,6)
      dimension x2(3), xd2(3), cov2(6,6)
c     'm' is implicitly INTEGER (outside a-h, o-z) -- declare explicitly
      double precision miss_dist

      write(iunit,'(A)')   'CCSDS_CDM_VERS = 1.0'
      write(iunit,'(A,A)') 'CREATION_DATE  = ', creation_str
      write(iunit,'(A)')   'ORIGINATOR     = KSROP'
      write(iunit,'(A)')   ''
      write(iunit,'(A,A)') 'TCA                   = ', tca_str
      write(iunit,'(A,F14.6,A)') 'MISS_DISTANCE         = ',
     &                            miss_dist, ' [km]'
      write(iunit,'(A,F14.6,A)') 'RELATIVE_SPEED        = ',
     &                            rel_speed, ' [km/s]'
      if (pc .ge. 0.0d0) then
         write(iunit,'(A,1PE14.6)') 'COLLISION_PROBABILITY = ', pc
      end if
      write(iunit,'(A)')   ''

      call write_cdm_object(iunit, 1, name1, x1, xd1, cov1)
      write(iunit,'(A)')   ''
      call write_cdm_object(iunit, 2, name2, x2, xd2, cov2)

      return
      end

c     ============================================================
      subroutine write_cdm_object(iunit, iobj, name, x, xd, cov)
c     Helper for write_cdm: emit one OBJECTn metadata, state-vector
c     (X,Y,Z,X_DOT,Y_DOT,Z_DOT, EME2000) and lower-triangular RTN
c     covariance (CR_R ... CNDOT_NDOT, row/col order R,T,N,RDOT,TDOT,
c     NDOT) block, with units selected from the number of velocity
c     ([km**2], [km**2/s] or [km**2/s**2]) axes in each element.
c     ============================================================
      implicit double precision (a-h, o-z)
      character*(*) name
      character*7   objstr
      character*10  cdm_lbl(21)
      character*16  unitstr
      dimension x(3), xd(3), cov(6,6)
      dimension irow(21), icol(21)
      integer       iobj, irow, icol, k, nvel

      data cdm_lbl /
     &  'CR_R      ', 'CT_R      ', 'CT_T      ',
     &  'CN_R      ', 'CN_T      ', 'CN_N      ',
     &  'CRDOT_R   ', 'CRDOT_T   ', 'CRDOT_N   ', 'CRDOT_RDOT',
     &  'CTDOT_R   ', 'CTDOT_T   ', 'CTDOT_N   ',
     &  'CTDOT_RDOT', 'CTDOT_TDOT',
     &  'CNDOT_R   ', 'CNDOT_T   ', 'CNDOT_N   ',
     &  'CNDOT_RDOT', 'CNDOT_TDOT', 'CNDOT_NDOT' /
      data irow / 1,2,2, 3,3,3, 4,4,4,4, 5,5,5,5,5, 6,6,6,6,6,6 /
      data icol / 1,1,2, 1,2,3, 1,2,3,4, 1,2,3,4,5, 1,2,3,4,5,6 /

      if (iobj .eq. 1) then
         objstr = 'OBJECT1'
      else
         objstr = 'OBJECT2'
      end if

      write(iunit,'(A,A)')        'OBJECT               = ', objstr
      write(iunit,'(A,A)')        'OBJECT_DESIGNATOR    = ', name
      write(iunit,'(A)')          'CATALOG_NAME         = UNKNOWN'
      write(iunit,'(A,A)')        'OBJECT_NAME          = ', name
      write(iunit,'(A)')          'EPHEMERIS_NAME       = NONE'
      write(iunit,'(A)')          'COVARIANCE_METHOD    = CALCULATED'
      write(iunit,'(A)')          'MANEUVERABLE         = N/A'
      write(iunit,'(A)')          'REF_FRAME            = EME2000'
      write(iunit,'(A,F16.6,A)')  'X                    = ', x(1),
     &                            ' [km]'
      write(iunit,'(A,F16.6,A)')  'Y                    = ', x(2),
     &                            ' [km]'
      write(iunit,'(A,F16.6,A)')  'Z                    = ', x(3),
     &                            ' [km]'
      write(iunit,'(A,F16.9,A)')  'X_DOT                = ', xd(1),
     &                            ' [km/s]'
      write(iunit,'(A,F16.9,A)')  'Y_DOT                = ', xd(2),
     &                            ' [km/s]'
      write(iunit,'(A,F16.9,A)')  'Z_DOT                = ', xd(3),
     &                            ' [km/s]'

      do k = 1, 21
         nvel = 0
         if (irow(k) .gt. 3) nvel = nvel + 1
         if (icol(k) .gt. 3) nvel = nvel + 1
         if (nvel .eq. 0) then
            unitstr = ' [km**2]'
         else if (nvel .eq. 1) then
            unitstr = ' [km**2/s]'
         else
            unitstr = ' [km**2/s**2]'
         end if
         write(iunit,'(A,A,1PE14.6,A)') cdm_lbl(k), ' = ',
     &                cov(irow(k),icol(k)), unitstr
      end do

      return
      end
