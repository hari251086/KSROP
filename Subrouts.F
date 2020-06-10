c     -------------------------------------------------------------
c     integrators
c     -------------------------------------------------------------
c     ks-orbital integrator
      function rkgil(n,y,f,x,h,nt)
      implicit double precision (a-h,o-z)
      dimension y(10),f(10),q(10) 
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
      
      common /xy/ pi,d2r,r2d,amue

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
      ai = dacos(h(3)/hmag)

      raan = dacos(a_n(1)/anmag)
      if(a_n(2).lt.0.d0) raan = 2.d0*pi-raan

      aop = acos(dotp3(a_n,ecc)/(anmag*e))
      if(ecc(3).lt.0.d0) aop = 2.d0*pi-aop

      anu = dacos(dotp3(ecc,x)/(e*rmag))
      if(dotp3(x,xd).lt.0.d0) then
         anu = 2.d0*pi - anu
      end if
      
!     Special Cases
!     Elliptical, equatorial. aop = true longitude of periapsis
      if(dabs(dabs(h(3)/hmag)-1.d0).lt.0.d0) then
         raan = 0.d0
         aop = dacos(ecc(1)/e)
         if(ecc(2).lt.0.d0) aop = 2.d0*pi - aop
      end if
!     Circular, inclined. anu = argument of latitude
      if(e.le.0.d0) then
         aop = 0.d0
         anu = dacos(dotp3(a_n,x)/(anmag*rmag))
         if(x(3).lt.0.d0) anu = 2.d0*pi - anu
         am = anu
      end if
!     Circular,equatorial. anu = true longitude
      if(dabs(dabs(h(3)/hmag)-1.d0).lt.0.d0.and.e.lt.0.d0) then
         raan = 0.d0
         aop = 0.d0
         anu = dacos(x(1)/rmag)
         if(x(2).lt.0.d0) anu = 2.d0*pi - anu
      end if

!     Mean anomaly
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
      if(ea.gt.2.d0*pi) ea = ea - 2.d0*pi
      if(am.lt.0.d0) am = am+2.d0*pi
      if(ea.lt.0.d0) ea = ea+2.d0*pi
      
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
      
      common /xy/ pi,d2r,r2d,amue

!     degree to radians
      do i = 3,6
         a_kep(i)=a_kep(i)*d2r
      end do
      
      aa = a_kep(1)
      ae = a_kep(2)
      ai = a_kep(3)
      aO = a_kep(4)
      aw = a_kep(5)
      am = a_kep(6)
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

      r=vmn(x(1),x(2),x(3))
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

c     ============================================================      
      subroutine car2sph(x,sph)
c     Conversion from Cartesian to Spherical coordinates
!     sph(1): right ascension in radians
!     sph(2): declination in radians
!     sph(3): position magnitude
c     ============================================================
      implicit double precision(a - h, o - z)
      dimension x(3),sph(3)

      sph(3) = vmn(x)
      sph(1) = atan3(x(2),x(1))
      sph(2) = dasin(x(3)/sph(3))

      return
      end

c     ============================================================      
      subroutine sph2ks(sph,u)
c     Conversion from Spherical to KS coordinates
c     ============================================================
      implicit double precision(a - h, o - z)
      dimension sph(3),u(4)

      pos_mag_root = dsqrt(sph(3))
      
      alby2 = sph(1)/2.d0
      decby2 = sph(2)/2.d0

!     Normalized
      u(1) = dcos(decby2)*dcos(alby2)
      u(2) = dcos(decby2)*dsin(alby2)
      u(3) = dsin(decby2)*dcos(alby2)
      u(4) = dsin(decby2)*dsin(alby2)

      do i = 1,4
         u(i) = u(i)*pos_mag_root
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
      subroutine car2ksnew(x,xd,u,us,w)
c     Conversion from Cartesian to KS coordinates
!     w: (un)perturbed angular frequency of harmonic oscillator
c     ============================================================
      implicit double precision (a-h,o-z)
      dimension u(4),us(4),x(3),xd(3)

c      print *,x,xd
      
      pos_mag = vmn(x)

      ks_mag = dsqrt(pos_mag)

      ang_alpha = atan3(x(2),x(1))
c      print *,ang_alpha
      ang_delta = dasin(x(3)/pos_mag)
c      print *,ang_delta
      u(1) = dcos(ang_delta/2.d0)*dcos(ang_alpha/2.d0)
      u(2) = dcos(ang_delta/2.d0)*dsin(ang_alpha/2.d0)
      u(3) = dsin(ang_delta/2.d0)*dcos(ang_alpha/2.d0)
      u(4) = dsin(ang_delta/2.d0)*dsin(ang_alpha/2.d0)

      do i = 1,4

         u(i) = u(i)*ks_mag
         
      end do

      cnst=.25d0/w 

      us(1)=(u(1)*xd(1)+u(2)*xd(2)+u(3)*xd(3))*cnst
      us(2)=(-u(2)*xd(1)+u(1)*xd(2)+u(4)*xd(3))*cnst
      us(3)=(-u(3)*xd(1)-u(4)*xd(2)+u(1)*xd(3))*cnst
      us(4)=(u(4)*xd(1)-u(3)*xd(2)+u(2)*xd(3))*cnst

      
      return
      end

c     ============================================================      
      Subroutine aLegP(n,x,v)
c     Computes associated Legendre Polynomial P_n(x)
      implicit double precision (a-h,o-z)
c     v(3) = P_2(x)
c     ============================================================
      dimension v(2190)
      if(n.lt.0) return
      v(1) = 1.d0               !	P_0(x) = 1
      if(n.lt.1) return
      v(2) = x                  !	P_1(x) = x
      
      do i = 3,n+2
         v(i) = 0.d0            !P_n(x) where n > 1
      end do
      
      do i = 3,n+2
         v(i) = (2.d0*dble(i)-1.d0)/dble(i)*x*v(i-1)
     &        -dble(i-1)/dble(i)*v(i-2)
      end do
      return
      end

c     ============================================================      
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
c     Computes Jn values from Data file (Cn0 values)
      implicit double precision (a-h,o-z)
      dimension c_j(n),c_raw(2401333,6),cn_0(n,6)
c     n = degree of the geopotential
c     c_j = j_n value
c     Number of lines in EGM2008:  2401333 
c     ============================================================

c     Data Files
      open (unit = 1, file = 'Constants\EGM2008_to2190_TideFree') !EGM2008 - cnm snm coefficients
      if(n.eq.0) go to 1000
      do i = 1,2401333

         read(1,*) c_raw(i,:)

      end do

      j = 2

      do i = 1,2401333

         if(c_raw(i,2).eq.0) then
            cn_0(j,3)=c_raw(i,3)
            j = j + 1
         end if

         if(j.gt.n) go to 10

      end do

 10   continue

      do i = 2,n

         c_j(i)=-dsqrt(2.d0*i+1.d0)*cn_0(i,3)

      end do

      close(1)
 1000 continue
 100  format(f30.20)
      return
      end

c     ============================================================      
      subroutine solarnpv(dj,s)
c     Compute sun position vector in inertial geo-centric frame
      implicit double precision(a-h,o-z,m)
      dimension s(3)
c     dj = Julian day
c     s = sun position vector (x_s, y_s, z_s) in km
      common /xy/ pi,d2r,r2d,amue,AU      
c     ============================================================     

      
c     julian day
      djn=(dj-2451545.0d0)

c     omega
      omega=2.1429d0-0.0010394594d0*djn

c     mean longitude of the sun
      tlam=(4.8950630d0+0.017202791698d0*djn)

c     mean anomaly of the sun
      tmmm=(6.24006d0+0.0172019699d0*djn)

c     ecliptic longitude
      tecl=tlam+0.03341607d0*dsin(tmmm)+0.00034894d0*dsin(2*tmmm)
     &     -0.0001134d0-0.0000203d0*dsin(omega)

c     obliquity of the ecliptic
      tobl=0.4090928d0-(6.214d-09)*djn+0.0000396d0*dcos(omega)

c     right ascension and declination
      deo=dsin(tecl)*dcos(tobl)
      dec=dcos(tecl)
      ralfa = atan2(deo,dec)
      if(ralfa.lt.0.0d0) ralfa=ralfa+2.0d0*pi
      rdecl=dasin(dsin(tobl)*dsin(tecl))

      sralfa = dsin(ralfa)
      cralfa = dcos(ralfa)
      srdecl = dsin(rdecl)
      crdecl = dcos(rdecl)
!      AU = 1.495978707d08
      s(1) = AU*cralfa*crdecl
      s(2) = AU*sralfa*crdecl
      s(3) = AU*srdecl

      close(1)
      return
      end

c     ============================================================
      subroutine lunarpv(dj,o)  
c     Compute the lunar position vector in geocentric reference frame
      implicit double precision(a-h,o-z)
      dimension o(3)
c     dj = Julian day
c     o = moon position vector (x_m, y_m, z_m) in km
c     ============================================================   
c     time in julian centuries
      tlun=(dj-2451545.0d0)/36525.0d0

c     x coordinate
      xm = 383.0d0*dsin(8399.685d0*tlun+5.381d0)
      xm = xm + 31.5d0*dsin(70.99d0*tlun+6.619d0)
      xm = xm + 10.6d0*dsin(16728.377d0*tlun+1.453d0)
      xm = xm + 6.2d0*dsin(1185.622d0*tlun+0.481d0)
      xm = xm + 3.2d0*dsin(7143.070d0*tlun+5.017d0)
      xm = xm + 2.3d0*dsin(15613.745d0*tlun+0.857d0)
      xm = xm + 0.8d0*dsin(8467.263d0*tlun+1.01d0)
      xm = xm * 1.d3

c     y coordinate
      ym = 351.0d0*dsin(8399.687d0*tlun+3.811d0)
      ym = ym + 28.9d0*dsin(70.997d0*tlun+4.596d0)
      ym = ym + 13.7d0*dsin(8433.466d0*tlun+4.766d0)
      ym = ym + 9.7d0*dsin(16728.380d0*tlun+6.165d0)
      ym = ym + 5.7d0*dsin(1185.667d0*tlun+5.164d0)
      ym = ym + 2.9d0*dsin(7143.058d0*tlun+0.3d0)
      ym = ym + 2.1d0*dsin(15613.755d0*tlun+5.565d0)
      ym = ym * 1.d3

c     z coordinate
      zm = 153.2d0*dsin(8399.672d0*tlun+3.807d0)
      zm = zm + 31.5d0*dsin(8433.464d0*tlun+1.629d0)
      zm = zm + 12.5d0*dsin(70.996d0*tlun+4.595d0)
      zm = zm + 4.2d0*dsin(16728.364d0*tlun+6.162d0)
      zm = zm + 2.5d0*dsin(1185.645d0*tlun+5.167d0)
      zm = zm + 3.0d0*dsin(104.881d0*tlun+2.555d0)
      zm = zm + 1.8d0*dsin(8399.116d0*tlun+6.248d0)
      zm = zm * 1.d3

      o(1) = xm
      o(2) = ym
      o(3) = zm
      
      return
      end

c     ============================================================            
      subroutine cal2jd(cal,dj)
c     Returns julian date when calendar date is given
      implicit double precision(a-h,o-z)
      dimension cal(6)
c     cal(1): yyyy
c     cal(2): mm
c     cal(3): dd
c     cal(4): hh
c     cal(5): min
c     cal(6): seconds
c     ============================================================
      i_mon1 = (cal(2) + 9.d0)/12.d0
      i_mon2 = (275.d0 * cal(2))/9.d0
      iyr=(7.0d0*(cal(1)+i_mon1))/4.0d0
      tyr1=367.0d0*cal(1)
      ft1=(cal(6)/60.0d0)
      ft2=ft1+cal(5)
      ft3=ft2/60.0d0
      ft4=ft3+cal(4)
      ftime=ft4/24.0d0
      ccday1=cal(3)+1721013.5d0
      ddj=tyr1-jyr+i_mon2+ccday1+ftime
      return
      end
      
      
c     ============================================================      
      subroutine force_models(n_for,ngeo,s,m)
c     Selection of force models
      implicit double precision(a-h,o-z)
      dimension n_for(3)
c     n_for(1): oblateness, n_for(2): sun, n_for(3): moon
c     1: On, 0: Off
c     ============================================================
      if(n_for(1).eq.0) ngeo = 0.d0
      if(n_for(2).eq.0) s = 0.d0
      if(n_for(3).eq.0) m = 0.d0
      
      return
      end
c$$$c     ----------------------------------------
c$$$c     Matrix and Vector operations
c$$$c     ----------------------------------------
c$$$
c$$$c     ======================================
c$$$      function dotp3(x,y)
c$$$!     Dot product of 3-dimensional vector
c$$$c     ======================================
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension x(3),y(3)
c$$$      dotp3=0.d0
c$$$      do i = 1,3
c$$$         dotp3=x(i)*y(i)+dotp3
c$$$      end do
c$$$      return
c$$$      end
      
c$$$c     -----------------------------------------------------
c$$$c     ks-space to ric direct transformation
c$$$      subroutine kstoric(u,us,rr,ri,rc,w)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension u(4),us(4),rr(4),ri(4),rc(4)
c$$$      dimension rx(4),vx(4),rxv(4),cxr(4)
c$$$      call ks2car(u,us,rx,vx,w)
c$$$c     r
c$$$      do i = 1,4
c$$$         rr(i)=rx(i)/dotp4(u,u)
c$$$      end do
c$$$c     c
c$$$      call crossp4(rr,vx,rxv)
c$$$      call vunit4(rxv,rc)
c$$$c     i
c$$$      call crossp4(rc,rr,ri)
c$$$      return
c$$$      end
c$$$c     -----------------------------------------------------
c$$$c     ks-space to orbit frame transformation through eci
c$$$      subroutine kstoric(u,us,xo,yo,zo,w)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension u(4),us(4),xo(4),yo(4),zo(4)
c$$$      dimension rx(4),vx(4),vxr(4)
c$$$      call ks2car(u,us,rx,vx,w)
c$$$c     z_o
c$$$      call vunit4(-rx,zo)
c$$$c     y_o
c$$$      call crossp4(vx,rx,vxr)
c$$$      call vunit4(rxv,yo)
c$$$c     z
c$$$      call crossp4(yo,zo,xo)
c$$$      return
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$c     cartesian to ric transformation
c$$$      subroutine ecitoric(x,xd,rr,ri,rc)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension x(3),xd(3),rr(3),ri(3),rc(3),rxv(3)
c$$$c     r
c$$$      call vunit(x,rr)
c$$$c     c
c$$$      call crossp(x,xd,rxv)
c$$$      call vunit(rc,rxv)
c$$$c     i
c$$$      call crossp(rc,rr,ri)
c$$$      return
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$c     dcm from eci
c$$$c     x-axis is velocity vector
c$$$c     z-axis is anti-nadir vector
c$$$c     y-axis completes the triad
c$$$      subroutine eci2orb(r,v,dcm)
c$$$      implicit double precision(a-h,o-z)
c$$$      dimension r(3),v(3),ro(3),vo(3),ru(3),vu(3)
c$$$      dimension dcm(3,3)
c$$$      ro(1)=0
c$$$      ro(2)=0
c$$$      ro(3)=1
c$$$      vo(1)=1
c$$$      vo(2)=0
c$$$      vo(3)=0
c$$$
c$$$      call vunit(r,ru)
c$$$      call vunit(v,vu)
c$$$      
c$$$      call getattitude(ru,ro,vu,vo,dcm,quat)
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$c     direction cosine matrix and quaternion relating the body reference frame to
c$$$c     the inertial reference frame
c$$$c     dcm times the reference frame vector will give the vector in the body frame
c$$$c     a,b,c,d are unit vectors
c$$$      subroutine getattitude(a,b,c,d,dcm,quat)
c$$$      implicit double precision(a-h,o-z)
c$$$      dimension a(3),b(3),c(3),d(3),dcm(3,3)
c$$$      dimension rvi(3),rvo(3),rhi(3),rho(3),quat(4)
c$$$      dimension refi(3,3),refo(3,3)
c$$$
c$$$      call crossp(a,c,rvi)
c$$$      call crossp(b,d,rvo)
c$$$      call crossp(a,rvi,rhi)
c$$$      call crossp(b,rvo,rho)
c$$$      
c$$$      do i = 1,3
c$$$         refi(1,i)=a(i)
c$$$         refi(2,i)=rvi(i)
c$$$         refi(3,i)=rhi(i)
c$$$         refo(i,1)=b(i)
c$$$         refo(i,2)=rvo(i)
c$$$         refo(i,3)=rho(i)         
c$$$      end do
c$$$
c$$$      call mat_mul3(dcm,refo,refi)
c$$$
c$$$      call dcm2quat(quat,dcm)
c$$$      
c$$$      end      
c$$$c     -------------------------------------------------------------
c$$$c     astronomical algorithms
c$$$c     -------------------------------------------------------------
c$$$      subroutine jdate(iyr,imon,iday,ih,imin,isec,ddj)
c$$$c     conversion from calendar date and time to julian date
c$$$      implicit double precision(a-h,o-z)
c$$$      intmon1=(imon+9.0d0)/12.0d0
c$$$      intmon2=(275.0d0*imon)/9.0d0
c$$$      jyr=(7.0d0*(iyr+intmon1))/4.0d0
c$$$      tyr1=367.0d0*iyr
c$$$      ft1=(isec/60.0d0)
c$$$      ft2=ft1+imin
c$$$      ft3=ft2/60.0d0
c$$$      ft4=ft3+ih
c$$$      ftime=ft4/24.0d0
c$$$      ccday1=iday+1721013.5d0
c$$$      ddj=tyr1-jyr+intmon2+ccday1+ftime
c$$$      return
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$      subroutine jul2d(ddj,iyr,imon,iday,ih,imin,isec)
c$$$c     conversion from  julian date to calendar day
c$$$      implicit double precision(a-h,o-z)
c$$$      ddj5=ddj+0.5d0
c$$$      fint = int(ddj5)
c$$$      fdec=ddj5-fint
c$$$      if(fint.lt.2299161) then
c$$$         acc=fint
c$$$      else
c$$$         alfa=int((fint-1867216.25)/36524.25)
c$$$         acc=fint+1+alfa-int(alfa/4)
c$$$      end if
c$$$      bcc=acc+1524
c$$$      iccc=int((bcc-122.1)/365.25)
c$$$      idcc=int(365.25*iccc)
c$$$      iecc=int((bcc-idcc)/30.6001)
c$$$      day=bcc-idcc-int(30.6001*iecc)+fdec
c$$$      if(iecc.lt.14) then
c$$$         imon=iecc-1
c$$$      else
c$$$         imon=iecc-13
c$$$      end if
c$$$      if(imon.gt.2) then
c$$$         iyr=iccc-4716
c$$$      else
c$$$         iyr=iccc-4715
c$$$      end if
c$$$      iday=int(day)
c$$$      call day2hms(day,ih,imin,isec)
c$$$      return
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$      subroutine jul3d(ddj,iyr,imon,iday,ih,imin,isec)
c$$$c     conversion from  julian date to calendar day      
c$$$      implicit double precision(a-h,o-z)
c$$$      fint=int(ddj)
c$$$      fdec=ddj-fint
c$$$      if(fint.lt.2299161) then
c$$$         acc=fint
c$$$      else
c$$$         alfa=int((fint-1867216.25)/36524.25)
c$$$         acc=fint+1+alfa-int(alfa/4)
c$$$      end if
c$$$      bcc=acc+1524
c$$$      iccc=int((bcc-122.1)/365.25)
c$$$      idcc=int(365.25*iccc)
c$$$      iecc=int((bcc-idcc)/30.6001)
c$$$      day=bcc-idcc-int(30.6001*iecc)+fdec
c$$$      if(iecc.lt.14) then
c$$$         imon=iecc-1
c$$$      else
c$$$         imon=iecc-13
c$$$      end if
c$$$      if(imon.gt.2) then
c$$$         iyr=iccc-4716
c$$$      else
c$$$         iyr=iccc-4715
c$$$      end if
c$$$      iday=int(day)
c$$$      call day2hms(day,ih,imin,isec)
c$$$      
c$$$      return
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$      subroutine day2hms(d,ih,imn,isc)
c$$$c     days to hours, minutes and seconds
c$$$      implicit double precision(a-h,o-z)
c$$$      iday=int(d)
c$$$      frday=d-iday
c$$$      ih=int(frday*24)
c$$$      amn=frday*24-ih
c$$$      imn=int(amn*60)
c$$$      asc=amn*60-imn
c$$$      isc=int(asc*60)
c$$$      if(isc.eq.60) then
c$$$         isc=0
c$$$         imn=imn+1
c$$$      end if
c$$$      if(imn.eq.60) then
c$$$         imn=0
c$$$         ih=ih+1
c$$$      end if
c$$$      return
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$      subroutine cthms(time,dhour,dmin,secs)
c$$$c     seconds to hours, minutes and seconds      
c$$$      implicit double precision(a-h,o-z)
c$$$      dhour=0.0d0
c$$$      time=time-60.0d0*dhour
c$$$      dmin=dint(time/60.0d0)
c$$$      secs=time-60.0d0*dmin
c$$$      return
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$      subroutine csmhjd(dhour,dminute,seconds,djt)
c$$$c     seconds,minutes,hours to julian date      
c$$$      implicit double precision(a-h,o-z)
c$$$      dhour=0.0d0
c$$$      s1=seconds/60.0d0
c$$$      s2=(s1+dminute)/60.0d0
c$$$      djt=(s2+dhour)/24.0d0
c$$$      return
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$      function atan2(a,b) 
c$$$      implicit double precision (a-h,o-z)
c$$$      pi=4.0d0*datan(1.0d0)
c$$$      pib4=pi/4.0d0
c$$$      pi3b4=pib4*3.0d0
c$$$      absa=dabs(a)
c$$$      if(b.ge.0.0d0) atan2=pib4-pib4*((b-absa)/(b+absa))
c$$$      if(b.lt.0.0d0) atan2=pi3b4-pib4*((b+absa)/(absa-b))
c$$$      if(a.gt.0.0d0) atan2=-atan2
c$$$      return
c$$$      end 
c$$$c     -----------------------------------------------------------
c$$$      subroutine lunarpv(dj)  
c$$$c     lunar vector computation
c$$$c     1999 gsfc issfd paper
c$$$      implicit double precision(a-h,o-z)
c$$$      common /uv/ x1,x2,x3,x1dot,x2dot,x3dot,xs,ys,zs,xm,ym,zm
c$$$      common /xy/ pi,cd,cr,re
c$$$   
c$$$c     time in julian centuries
c$$$      tlun=(dj-2451545.0d0)/36525.0d0
c$$$
c$$$c     x coordinate
c$$$      a11=383.0d0*dsin(8399.685d0*tlun+5.381d0)
c$$$      a12=31.5d0*dsin(70.99d0*tlun+6.619d0)
c$$$      a13=10.6d0*dsin(16728.377d0*tlun+1.453d0)
c$$$      a14=6.2d0*dsin(1185.622d0*tlun+0.481d0)
c$$$      a15=3.2d0*dsin(7143.070d0*tlun+5.017d0)
c$$$      a16=2.3d0*dsin(15613.745d0*tlun+0.857d0)
c$$$      a17=0.8d0*dsin(8467.263d0*tlun+1.01d0)
c$$$
c$$$      xm=(a11+a12+a13+a14+a15+a16+a17)*1000.0d0
c$$$
c$$$c     y coordinate
c$$$      a21=351.0d0*dsin(8399.687d0*tlun+3.811d0)
c$$$      a22=28.9d0*dsin(70.997d0*tlun+4.596d0)
c$$$      a23=13.7d0*dsin(8433.466d0*tlun+4.766d0)
c$$$      a24=9.7d0*dsin(16728.380d0*tlun+6.165d0)
c$$$      a25=5.7d0*dsin(1185.667d0*tlun+5.164d0)
c$$$      a26=2.9d0*dsin(7143.058d0*tlun+0.3d0)
c$$$      a27=2.1d0*dsin(15613.755d0*tlun+5.565d0)
c$$$
c$$$      ym=(a21+a22+a23+a24+a25+a26+a27)*1000.0d0
c$$$
c$$$c     z coordinate
c$$$      a31=153.2d0*dsin(8399.672d0*tlun+3.807d0)
c$$$      a32=31.5d0*dsin(8433.464d0*tlun+1.629d0)
c$$$      a33=12.5d0*dsin(70.996d0*tlun+4.595d0)
c$$$      a34=4.2d0*dsin(16728.364d0*tlun+6.162d0)
c$$$      a35=2.5d0*dsin(1185.645d0*tlun+5.167d0)
c$$$      a36=3.0d0*dsin(104.881d0*tlun+2.555d0)
c$$$      a37=1.8d0*dsin(8399.116d0*tlun+6.248d0)
c$$$
c$$$      zm=(a31+a32+a33+a34+a35+a36+a37)*1000.0d0
c$$$
c$$$      return
c$$$      end
c$$$c     ----------------------------------------
c$$$      subroutine solarnpv(dj)
c$$$c     solar vector computation
c$$$c     computing solar vector paper (psa algorithm)
c$$$      implicit double precision(a-h,o-z,m)
c$$$      common /uv/ x1,x2,x3,x1dot,x2dot,x3dot,xs,ys,zs,xm,ym,zm
c$$$      common /xy/ pi,cd,cr,re
c$$$      common /sun/ ang,dip,gam
c$$$      
c$$$c     julian day
c$$$      djn=(dj-2451545.0d0)
c$$$
c$$$c     omega
c$$$      omega=2.1429d0-0.0010394594d0*djn
c$$$
c$$$c     mean longitude of the sun
c$$$      tlam=(4.8950630d0+0.017202791698d0*djn)
c$$$
c$$$c     mean anomaly of the sun
c$$$      tmmm=(6.24006d0+0.0172019699d0*djn)
c$$$
c$$$c     ecliptic longitude
c$$$
c$$$      tecl=tlam+0.03341607d0*dsin(tmmm)   
c$$$      tecl=tecl+0.00034894d0*dsin(2*tmmm)
c$$$      tecl=tecl-0.0001134d0
c$$$      tecl=tecl-0.0000203d0*dsin(omega)
c$$$
c$$$c     obliquity of the ecliptic
c$$$      tobl=0.4090928d0
c$$$      tobl=tobl-(6.214d-09)*djn
c$$$      tobl=tobl+0.0000396d0*dcos(omega)
c$$$
c$$$c     right ascension and declination
c$$$      deo=dsin(tecl)*dcos(tobl)
c$$$      dec=dcos(tecl)
c$$$      ralfa = atan2(deo,dec)
c$$$      if(ralfa.lt.0.0d0) ralfa=ralfa+2.0d0*pi
c$$$      rdecl=dasin(dsin(tobl)*dsin(tecl))
c$$$
c$$$      sralfa=dsin(ralfa)
c$$$      cralfa=dcos(ralfa)
c$$$      srdecl=dsin(rdecl)
c$$$      crdecl=dcos(rdecl)
c$$$      trs=1.495978707d08
c$$$      xs=trs*cralfa*crdecl
c$$$      ys=trs*sralfa*crdecl
c$$$      zs=trs*srdecl
c$$$
c$$$      call cpvoe(a,ec,ai,omgad,somgad,ea,am)
c$$$      ai=ai*cr
c$$$      omg=omgad*cr
c$$$      somg=somgad*cr
c$$$      
c$$$      csw=dcos(somg)
c$$$      ssw=dsin(somg)
c$$$      cbw=dcos(omg)
c$$$      sbw=dsin(omg)
c$$$      si=dsin(ai)
c$$$      ci=dcos(ai)
c$$$
c$$$      al1=csw*cbw-ci*ssw*sbw
c$$$      am1=csw*sbw+ci*ssw*cbw
c$$$      an1=ssw*si
c$$$      
c$$$      al2=-ssw*cbw-ci*csw*sbw
c$$$      am2=ci*csw*cbw-ssw*sbw
c$$$      an2=csw*si
c$$$
c$$$      al3=si*sbw
c$$$      am3=-si*cbw
c$$$      an3=ci
c$$$
c$$$      als=cralfa*crdecl
c$$$      ams=sralfa*crdecl
c$$$      ans=srdecl
c$$$
c$$$      aa=al1*als+am1*ams+an1*ans
c$$$      bb=al2*als+am2*ams+an2*ans
c$$$      cc=al3*als+am3*ams+an3*ans
c$$$
c$$$      ang=dacos(aa)*cd
c$$$      dip=dasin(cc)*cd
c$$$      gam=datan2(bb,aa)
c$$$      if(gam.lt.0.d0) gam=gam+2.d0*pi
c$$$      gam=gam*cd
c$$$      
c$$$
c$$$      return
c$$$      end
c$$$c     ----------------------------------------
c$$$      function shadfncyl(x1,x2,x3,xs,ys,zs)
c$$$      implicit double precision(a-h,o-z,m)
c$$$      common /xy/ pi,cd,cr,re
c$$$c     ----------------------------------------      
c$$$c     cylindrical shadown function obtained from
c$$$c     hubaux et al., 2012
c$$$c     ----------------------------------------      
c$$$      x=x1
c$$$      y=x2
c$$$      z=x3
c$$$      gam=1.d9
c$$$      rs=vmn(xs,ys,zs)
c$$$      r=vmn(x,y,z)
c$$$      cosphi=(x*xs+y*ys+z*zs)/rs
c$$$      scr=cosphi+dsqrt(r*r-re*re)
c$$$      shadfncyl=0.5d0*(1.d0+dtanh(gam*scr))
c$$$      return
c$$$      end
c$$$c     ----------------------------------------
c$$$      function shadfncone(x1,x2,x3,xs,ys,zs)
c$$$      implicit double precision(a-h,o-z,m)
c$$$      common /xy/ pi,cd,cr,re
c$$$c     ----------------------------------------
c$$$c     conical shadown function obtained from
c$$$c     dg cook, 2001 masters thesis
c$$$c     ----------------------------------------
c$$$      rds=696342.d0
c$$$      x=x1
c$$$      y=x2
c$$$      z=x3
c$$$      xp=xs-x
c$$$      yp=ys-y
c$$$      zp=zs-z
c$$$      rs=vmn(xs,ys,zs)
c$$$      r=vmn(x,y,z)
c$$$      p=vmn(xp,yp,zp)
c$$$      et=dacos(-(x*xp+y*yp+z*zp)/(r*p))
c$$$      et2=et*et
c$$$      rhoe=dasin(re/r)
c$$$      rhos=dasin(rds/p)
c$$$      rhoe2=rhoe*rhoe
c$$$      rhos2=rhos*rhos
c$$$      rh1=dabs(rhos2-rhoe2)
c$$$      rh2=dabs(rhos2+rhoe2)
c$$$
c$$$      if(et.ge.(rhos+rhoe)) go to 1
c$$$      if(et.le.(rhoe-rhos)) go to 2
c$$$      if(et.gt.(rhoe-rhos).and.et.lt.(rhoe+rhos)) go to 3
c$$$
c$$$      go to 4
c$$$ 1    shadfncone=1.d0
c$$$      go to 4
c$$$ 2    shadfncone=0.d0
c$$$      go to 4
c$$$ 3    al=0.5d0*(rhos+rhoe+et)
c$$$      q=2.d0*dsqrt(al*(al-et)*(al-rhos)*(al-rhoe))/et
c$$$      if(rh1.le.et2) t=rhos2*dasin(q/rhos)+rhoe2*dasin(q/rhoe)-q*et
c$$$      if(rh1.gt.et2) t=rhos2*(pi-dasin(q/rhos))+rhoe2*dasin(q/rhoe)-q*et
c$$$      shadfncone=1.d0-t/(pi*rhos2)
c$$$
c$$$ 4    continue
c$$$      return
c$$$      end
c$$$c     -----------------------------------------------------
c$$$      function shadfnnew(x1,x2,x3,xs,ys,zs)
c$$$      implicit double precision(a-h,o-z,m)
c$$$      common /xy/ pi,cd,cr,re
c$$$      dimension rei(5),wei(5),fpi(5)
c$$$c     ----------------------------------------
c$$$c     conical shadow function obtained from
c$$$c     frueh et al., 2013
c$$$c     ----------------------------------------
c$$$      rds=696342.d0
c$$$      rem=6371.0088d0
c$$$      x=x1
c$$$      y=x2
c$$$      z=x3
c$$$      xp=xs-x
c$$$      yp=ys-y
c$$$      zp=zs-z
c$$$      p=vmn(xp,yp,zp)
c$$$      rs=vmn(xs,ys,zs)
c$$$      r=vmn(x,y,z)
c$$$
c$$$      upx=xp/p
c$$$      upy=yp/p
c$$$      upz=zp/p
c$$$
c$$$      ux=x/r
c$$$      uy=y/r
c$$$      uz=z/r
c$$$
c$$$      dep=upx*ux+upy*uy+upz*uz
c$$$
c$$$      atau=dasin(rds/p)
c$$$      epp=dacos(dep)
c$$$
c$$$ 20   do i = 1,4
c$$$         rei(i)=rem+(i-1)*16.d0
c$$$         wei(i)=2.d0/9.d0         
c$$$      end do
c$$$      rei(5)=rem-40.d0
c$$$      wei(5)=1.d0/9.d0
c$$$
c$$$      do i = 1,5
c$$$         agam=dasin(rei(i)/r)
c$$$         sg=(atau+epp+agam)/2.d0
c$$$         kg1=(sg-atau)
c$$$         kg2=(sg-agam)*kg1
c$$$         kg3=(sg-epp)*kg2
c$$$         kg=dsqrt(kg3*sg)
c$$$         sdel=2.d0*kg/(epp*agam)
c$$$         cdel=(epp**2+agam**2-atau**2)/(2.d0*epp*agam)
c$$$         sbet=2.d0*kg/(epp*atau)
c$$$         cbet=(epp**2+atau**2-agam**2)/(2.d0*epp*atau)
c$$$         bet=datan2(sbet,cbet)
c$$$         del=datan2(sdel,cdel)
c$$$         s2bt=sbet*cbet
c$$$         s2dt=sdel*cdel
c$$$         bri=(del-s2dt)+(bet-s2bt)
c$$$         fpi(i)=1-(bri*(agam/atau)**2)/pi
c$$$      end do
c$$$
c$$$      do i = 1,5
c$$$         shadfnnew=shadfnnew+wei(i)*fpi(i)
c$$$      end do
c$$$      
c$$$
c$$$ 19   continue
c$$$      return
c$$$      end      
c$$$

      
c$$$c     ----------------------------------------
c$$$      subroutine vunit(x,xu)
c$$$c     unit vector - 3 dimension
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension xu(3),x(3)
c$$$      ru=vmn(x(1),x(2),x(3))
c$$$      do i = 1,3
c$$$         xu(i)=x(i)/ru
c$$$      end do
c$$$      return
c$$$      end
c$$$c     ----------------------------------------
c$$$      subroutine vunit4(x,xu)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension x(4),xu(4)
c$$$      ru=vmn4(x)
c$$$      do i = 1,4
c$$$         xu(i)=x(i)/ru
c$$$      end do
c$$$      return
c$$$      end
c$$$c     ----------------------------------------
c$$$c     magnitude of a vector - 4 dimension
c$$$      function vmn4(x)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension x(4)
c$$$      vmn4=dsqrt(dotp4(x,x))
c$$$      return
c$$$      end
c$$$c     dot product of two vectors (n = 4)
c$$$      function dotp4(x,y)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension x(4),y(4)
c$$$      dotp4=0.d0
c$$$      do i = 1,4
c$$$         dotp4=x(i)*y(i)+dotp4
c$$$      end do
c$$$      return
c$$$      end
c     ----------------------------------------
c     cross product of two vectors
c     ----------------------------------------
c$$$c     cross product of two vectors - 4 dimensions (quaternion based)
c$$$      subroutine crossp4(x,y,z)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension x(4),y(4),z(4)
c$$$      z(1)=x(1)*y(4)+x(2)*y(3)-x(3)*y(2)-x(4)*y(1)
c$$$      z(2)=x(2)*y(4)-x(1)*y(3)-x(4)*y(2)+x(3)*y(1)
c$$$      z(3)=x(3)*y(4)-x(4)*y(3)+x(1)*y(2)-x(2)*y(1)
c$$$      z(4)=0.d0
c$$$      return
c$$$      end      
c$$$c     ----------------------------------------
c$$$c     sum of the diagonal elements in a matrix
c$$$      function trace(a)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension a(3,3)
c$$$      trace=0.d0
c$$$      do i = 1,3
c$$$         trace=trace+a(i,i)
c$$$      end do
c$$$      end function
c$$$c     ----------------------------------------
c$$$c     conversion from dcm to quaternion
c$$$      subroutine dcm2quat(q,dcm)
c$$$      implicit double precision (a-h,o-z)
c$$$      dimension q(4),dcm(3,3)
c$$$      trdcm=trace(dcm)
c$$$      q(1)=dsqrt(trdcm*0.25d0)
c$$$      q(2)=dcm(3,2)-dcm(2,3)
c$$$      q(3)=dcm(1,3)-dcm(3,1)
c$$$      q(4)=dcm(2,1)-dcm(1,2)
c$$$      do i = 2,4
c$$$         q(i)=q(i)*0.25d0/q(1)
c$$$      end do
c$$$      end
c$$$c     ----------------------------------------
c$$$c     matrix multiplication (3 x 3)
c$$$c     a: output
c$$$c     b,c: inputs
c$$$      subroutine mat_mul3(a,b,c)
c$$$      implicit double precision(a-h,o-z)
c$$$      dimension a(3,3),b(3,3),c(3,3)
c$$$
c$$$      do ij = 1,3
c$$$         do jk = 1,3
c$$$            a(ij,jk)=0.d0
c$$$         end do
c$$$      end do
c$$$      
c$$$      do k = 1, 3
c$$$         do i = 1, 3
c$$$            do j = 1, 3
c$$$               a(i,j) = a(i,j) + b(i,k) * c(k,j)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      
c$$$      end
c$$$c     ----------------------------------------
c$$$
c$$$      subroutine quad(a,b,c,r1,r2)
c$$$c     quadratic equation roots
c$$$      implicit double precision(a-h,o-z,m)
c$$$      d=b*b-4.d0*a*c
c$$$      if(d.gt.0.d0) then
c$$$         r1=-b+dsqrt(d)/(2.d0*a)
c$$$         r2=-b-dsqrt(d)/(2.d0*a)
c$$$      elseif(d.eq.0.d0) then
c$$$         r1=-b/(2.d0*a)
c$$$         r2=r1
c$$$      endif         
c$$$      return
c$$$      end
c$$$c     ----------------------------------------
c$$$c     co-ordinate transformations
c$$$c     ----------------------------------------
c$$$
c$$$c     -----------------------------------------------------      
c$$$      subroutine car2arr(x,xd)
c$$$c     assigning cartesian elements into arrays      
c$$$      implicit double precision (a-h,o-z)
c$$$      common /uv/ x1,x2,x3,x1dot,x2dot,x3dot,xs,ys,zs,xm,ym,zm
c$$$      dimension x(3),xd(3)
c$$$      x(1)=x1
c$$$      x(2)=x2
c$$$      x(3)=x3
c$$$      xd(1)=x1dot
c$$$      xd(2)=x2dot
c$$$      xd(3)=x3dot
c$$$      return
c$$$      end
c$$$c     ----------------------------------------------------
c$$$      subroutine arr2car(x,xd)
c$$$c     assigning cartesian elements into arrays      
c$$$      implicit double precision (a-h,o-z)
c$$$      common /uv/ x1,x2,x3,x1dot,x2dot,x3dot,xs,ys,zs,xm,ym,zm
c$$$      dimension x(3),xd(3)
c$$$      x1=x(1)
c$$$      x2=x(2)
c$$$      x3=x(3)
c$$$      x1dot=xd(1)
c$$$      x2dot=xd(2)
c$$$      x3dot=xd(3)
c$$$      return
c$$$      end
c$$$c     --------------------------------------------------
c$$$
c$$$c     -------------------------------------------------------------
      function atan3(a,b) 
      implicit double precision (a-h,o-z)
      pi=4.0d0*datan(1.0d0)
      atan3=datan(dabs(a/b))
      if((b.lt.0.).and.(a.gt.0.))atan3=pi-atan3 
      if((b.lt.0.).and.(a.lt.0.))atan3=pi+atan3 
      if((b.gt.0.).and.(a.lt.0.))atan3=2.0d0*pi-atan3 
      return
      end 
c$$$
c$$$c     -------------------------------------------------------------
c$$$      subroutine intpol(x,y,m1,x1,y1)
c$$$      implicit double precision(a-h,o-z)
c$$$      dimension x(m1),y(m1,6),y1(1,6)
c$$$      err1=1.d-16
c$$$      cash=x1-x(1)
c$$$      if((abs(cash).lt.err1).or.x1.gt.x(m1)) go to 50
c$$$
c$$$      do 11 i=1,m1
c$$$         if(x1-x(i))11,13,10
c$$$ 10      if(x1-x(i+1))15,11,11
c$$$ 15      do ill=1,6
c$$$            y1(1,ill)=y(i,ill)+(y(i+1,ill)-y(i,ill))
c$$$     &           *(x1-x(i))/(x(i+1)-x(i))
c$$$         end do
c$$$         
c$$$         go to 12
c$$$ 13      do lil = 1,6
c$$$            y1(1,lil)=y(i,lil)            
c$$$         end do
c$$$ 11   continue
c$$$ 12   return
c$$$ 50   continue
c$$$!      write(*,600)x
c$$$ 600  format(/3x,"x1 is out of range",e15.6)
c$$$      stop
c$$$      end
c$$$c     -------------------------------------------------------------
c$$$      subroutine intpol2(x,y,m1,x1,y1,in)
c$$$      implicit double precision(a-h,o-z)
c$$$      dimension x(m1),y(m1,6),y1(1,6)
c$$$      min=in-3
c$$$      max=in+3
c$$$      err1=1.d-16
c$$$      cash=x1-x(min)
c$$$      if((abs(cash).lt.err1).or.x1.gt.x(max)) go to 50
c$$$
c$$$      do 11 i=min,max
c$$$         if(x1-x(i))11,13,10
c$$$ 10      if(x1-x(i+1))15,11,11
c$$$ 15      do ill=1,6
c$$$            y1(1,ill)=y(i,ill)+(y(i+1,ill)-y(i,ill))
c$$$     &           *(x1-x(i))/(x(i+1)-x(i))
c$$$         end do
c$$$         
c$$$         go to 12
c$$$ 13      do lil = 1,6
c$$$            y1(1,lil)=y(i,lil)            
c$$$         end do
c$$$ 11   continue
c$$$ 12   return
c$$$ 50   continue
c$$$!      write(*,600)x
c$$$ 600  format(/3x,"x1 is out of range",e15.6)
c$$$      stop
c$$$      end
c$$$      
