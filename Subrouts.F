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
      subroutine cal2jd(cal,ddj)
      implicit double precision(a-h,o-z)
      dimension cal(6)
      iyr = cal(1)
      imon = cal(2)
      iday = cal(3)
      ih = cal(1)
      imin = cal(2)
      isec = cal(3)
      intmon1=(imon+9.0d0)/12.0d0
      intmon2=(275.0d0*imon)/9.0d0
      jyr=(7.0d0*(iyr+intmon1))/4.0d0
      tyr1=367.0d0*iyr
      ft1=(isec/60.0d0)
      ft2=ft1+imin
      ft3=ft2/60.0d0
      ft4=ft3+ih
      ftime=ft4/24.0d0
      ccday1=iday+1721013.5d0
      ddj=tyr1-jyr+intmon2+ccday1+ftime
      return
      end


      function atan3(a,b) 
      implicit double precision (a-h,o-z)
      pi=4.0d0*datan(1.0d0)
      atan3=datan(dabs(a/b))
      if((b.lt.0.).and.(a.gt.0.))atan3=pi-atan3 
      if((b.lt.0.).and.(a.lt.0.))atan3=pi+atan3 
      if((b.gt.0.).and.(a.lt.0.))atan3=2.0d0*pi-atan3 
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
      open (unit = 1, file = 'EGM2008_to2190_TideFree') !EGM2008 - cnm snm coefficients

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
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth     
c     ============================================================     

c     Data files
c      open (unit = 1, file = 'const_new.dat')
c      read(1,*) amue,R_Earth,AU
      
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
      subroutine force_models(n_for,ngeo,s,m)
c     Selection of force models
      implicit double precision(a-h,o-z)
      dimension n_for(3)
c     n_for(1): oblateness, n_for(2): sun, n_for(3): moon
c     ============================================================
      if(n_for(1).eq.0) ngeo = 0.d0
      if(n_for(2).eq.0) s = 0.d0
      if(n_for(3).eq.0) m = 0.d0
      
      return
      end

      subroutine lunarpv(dj,tm)  
c     lunar vector computation
c     1999 gsfc issfd paper
      implicit double precision(a-h,o-z)
      dimension tm(3)
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
   
c     time in julian centuries
      tlun=(dj-2451545.0d0)/36525.0d0

c     x coordinate
      a11=383.0d0*dsin(8399.685d0*tlun+5.381d0)
      a12=31.5d0*dsin(70.99d0*tlun+6.619d0)
      a13=10.6d0*dsin(16728.377d0*tlun+1.453d0)
      a14=6.2d0*dsin(1185.622d0*tlun+0.481d0)
      a15=3.2d0*dsin(7143.070d0*tlun+5.017d0)
      a16=2.3d0*dsin(15613.745d0*tlun+0.857d0)
      a17=0.8d0*dsin(8467.263d0*tlun+1.01d0)

      tm(1)=(a11+a12+a13+a14+a15+a16+a17)*1000.0d0

c     y coordinate
      a21=351.0d0*dsin(8399.687d0*tlun+3.811d0)
      a22=28.9d0*dsin(70.997d0*tlun+4.596d0)
      a23=13.7d0*dsin(8433.466d0*tlun+4.766d0)
      a24=9.7d0*dsin(16728.380d0*tlun+6.165d0)
      a25=5.7d0*dsin(1185.667d0*tlun+5.164d0)
      a26=2.9d0*dsin(7143.058d0*tlun+0.3d0)
      a27=2.1d0*dsin(15613.755d0*tlun+5.565d0)

      tm(2)=(a21+a22+a23+a24+a25+a26+a27)*1000.0d0

c     z coordinate
      a31=153.2d0*dsin(8399.672d0*tlun+3.807d0)
      a32=31.5d0*dsin(8433.464d0*tlun+1.629d0)
      a33=12.5d0*dsin(70.996d0*tlun+4.595d0)
      a34=4.2d0*dsin(16728.364d0*tlun+6.162d0)
      a35=2.5d0*dsin(1185.645d0*tlun+5.167d0)
      a36=3.0d0*dsin(104.881d0*tlun+2.555d0)
      a37=1.8d0*dsin(8399.116d0*tlun+6.248d0)

      tm(3)=(a31+a32+a33+a34+a35+a36+a37)*1000.0d0

      return
      end
