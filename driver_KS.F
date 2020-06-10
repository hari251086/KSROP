c     COPYRIGHT 2018, HARISHKUMAR SELLAMUTHU, ALL RIGHTS RESERVED
C     -------------------------------------------------------------
c     KS - Numerical driver
c       Description:
c     Orbit propagation using Kustaanheimo-Stiefel regular elements
c     with numerical integrator RKG4. (Earth)
c     Perturbations included are:
c     - Earth's oblateness ( up to 2190th degree )
c     To be added
c     - Luni-solar perturbations
c     - Atmospheric drag ( Oblate exponential
c     atmosphere with varying density scale height )
c     - Solar radiation pressure (cannon ball)
c     -------------------------------------------------------------
c     Author: Harishkumar Sellamuthu
c     Ph.D candidate, Karunya Institute of Technology and Sciences
c     hari251086@gmail.com
c     -------------------------------------------------------------
c     Input file parameteres <input.dat> by line
c     --------      
c     Initial position in geocentric inertial reference frame in km
c     Initial velocity in geocentric inertial reference frame in km/s
c     Initial date in yyyy,mm,dd,hh,mm,ss (to be included)
c     Number of revolutions, number of steps, tolerance (10^-15 is default)
c     -------------------------------------------------------------c
c     Compiling with Intel Fortran compiler in command line:
c     ifort <driver_File>.F <Subroutines_File>.F
c     Execute with <driver_File> executable
c     -------------------------------------------------------------
c     ****Updates & Revisions****
c     --------
c     06-15-2018: New numerical program
c     Only J2
c     --------
c     06-16-2018: Trajectory dump at each revolution
c     Changing cn_j input to array,
c     Legendre Polynomial function from Numerical Recipes in F77
c     --------
c     09-13-2018: Added nth degree geopotential (Up to 2190 x 0)
c     EGM2008 model (Jn)
c     ngeo_deg value can be up to 2190 in the const_new.dat file
c     -------------------------------------------------------------      
c     ****Issues****
c     --------
c     Conversion from State to Keplerian elements needs to be
c     revisited for NaN outputs
c     -------------------------------------------------------------
c     Subroutines used
c     geo_coeff,car2ks,ks2car,aLegP,car2oe,ks2ksr
c     -------------------------------------------------------------c
      
      implicit double precision (a-h, o-z)
      parameter(n = 10,jn = 2401333,cn0 = 2190)
c     Value of pi, degree to radians and radians to degree
      common /xy/ pi,d2r,r2d,amue,AU
      dimension pek(7),x(3),xd(3),q(4),u(4),us(4),y(10),z(10)
      dimension qu(3),q_til(4),u_til(4),uu(4)
      dimension ZbyR(cn0*2),ObyR(cn0*2),c_raw(jn,6),Re(cn0)
      dimension cn_j(cn0),c_j(cn0),cn_0(cn0,6),V_p(cn0),aleg(cn0)
      dimension tau_sum(cn0),sleg(cn0),oleg(cn0),R(n*2)
      dimension cal(6),ts(3),tm(3),opt(3),n_force(3)
      dimension RbyRs(cn0),RbyRm(cn0),cn_m(cn0),cn_s(cn0),Rm(n*2)
      dimension qsun(4),qmoon(4),slambda(cn0),olambda(cn0),Rs(n*2)
      
      pi=4.d0*datan(1.d0)
      d2r=pi/180.d0
      r2d=180.d0/pi
      
c     Input Files
      open (unit = 1, file = 'const_new.dat') !Constants
      open (unit = 2, file = 'input.dat') !Input state vector

c     Output files
      open (unit = 3, file = 'state.out') !Output - State Vector
      open (unit = 4, file = 'regular.out') !Output - KS
      open (unit = 5, file = 'kepler.out') !Output - Kepler

      read(1,*) amue,R_Earth,AU,amuS,amuM
      read(1,*) ngeo_deg,nsun_deg,nmoon_deg
      
      read(2,*,end=99) x(1),x(2),x(3) ! Position
      read(2,*) xd(1),xd(2),xd(3) ! Velocity
      read(2,*) nrev,istep,tole !Number of revolutions, number of steps, tolerance
      read(2,*) cal(1),cal(2),cal(3),cal(4),cal(5),cal(6) !year,month,day,hour,minute,second
      read(2,*) n_force(1), n_force(2), n_force(3)

      call car2oe(x,xd,pek)
      call oe2car(pek,x,xd)
      call cal2jd(cal,djulian)
      call force_models(n_force,ngeo_deg,amuS,amuM)

      print *, x, xd
      print *, pek(1),pek(2),pek(3)*r2d,pek(4)*r2d,pek(5)*r2d,pek(6)*r2d


      
      write(*,*) '=============================================='

      if(ngeo_deg.eq.0.or.ngeo_deg.eq.1) then
         do i = 1,cn0
            c_j(i) = 0.d0
         end do
      else
!     compute j_n values
         call geo_coeff(ngeo_deg,c_j)
      end if

!     Initialization
      E_0 = 0.d0                !Gen. Eccentric Anomaly in radians (Used in dumps)
      E = E_0                   !Gen. Eccentric Anomaly in radians for integration
      T = 0.d0                  !Physical time in seconds
      nt = 0                    !Counter in the integrator
      
      dE_0 = 2.d0*pi/istep      !Step-size in radians (Used in dumps)
      
      num_steps = nrev*istep    !/(dE*r2d)          !Number of integration steps

      call solarnpv(djulian,ts) !Sun position vector in km/s (inertial geocentric frame)
      call lunarpv(djulian,tm)  !Moon position vector in km/s (inertial geocentric frame)
      
      R(1) = vmn(x)
      V = vmn(xd)

      R(2) = R(1)*R(1)
      V2 = V*V

      do i = 2, 20
         R(i) = R(i-1) * R(1)
      end do
      
      
!     Perturbing potential, V_pot

!     Oblateness terms
      Re(1) = R_Earth
      ZbyR(1) = x(3)/R(1)
      ObyR(1) = 1.d0/R(1)      

!     Third-body terms
      Rs(1) = vmn(ts)
      Rm(1) = vmn(tm)
      RbyRs(1) = R(1)/Rs(1)
      RbyRm(1) = R(1)/Rm(1)
      first_s = amuS/Rs(1)
      first_m = amuM/Rm(1)     
      do i = 1, nsun_deg
         RbyRs(i+1) = RbyRs(i)*RbyRs(1)
      end do
      do i = 1, nmoon_deg
         RbyRm(i+1) = RbyRm(i)*RbyRm(1)         
      end do
      ph_s = 0.d0
      ph_m = 0.d0
      do i = 1,3
         ph_s = ph_s + x(i) * ts(i)
         ph_m = ph_m + x(i) * tm(i)
      end do
      cphi_s = ph_s/(R(1) * Rs(1))
      cphi_m = ph_m/(R(1) * Rm(1))
      cphi_s2 = cphi_s  * cphi_s
      cphi_m2 = cphi_m  * cphi_m            
      
      do i = 1,ngeo_deg
         ObyR(i+1) = ObyR(i)*ObyR(1)
         Re(i+1) = Re(i)*Re(1)
      end do
          
      V_pot = 0.d0

      call aLegP(ngeo_deg,ZbyR(1),aleg)
      call aLegP(nsun_deg,cphi_s,sleg)
      call aLegP(nmoon_deg,cphi_m,oleg)      

      do i = 2,ngeo_deg
         cn_j(i) = amue*Re(i)*c_j(i)*ObyR(i+1)
         V_p(i) = cn_j(i)*aleg(i+1) !aLegendreP(i,0,ZbyR(1))
         V_pot = V_pot+V_p(i)
      end do

      do i = 2,nsun_deg
         cn_s(i) = -first_s * RbyRs(i)
         V_p(i) = cn_s(i) * sleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

      do i = 2,nmoon_deg
         cn_m(i) = -first_m * RbyRm(i)
         V_p(i) = cn_m(i) * oleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

!     Angular frequency

      w = dsqrt(0.5d0*(amue/R(1)-V2/2.d0-V_pot))
      w_kep = dsqrt(0.5d0*(amue/R(1)-V2/2.d0))
c     print *,w
!     Frequency Ratio
      Gam = w/w_kep

      dE = dE_0*Gam

      call car2ks(x,xd,u,us,w)
      call ks2car(u,us,x,xd,w)

c      print *,x,xd
c      print *,u,us
      
      Tow = T+dotp4(u,us)/w     !Time-element

!     Initial KS regular elements
      y(1) = Tow
      y(10) = w

      do i = 1,4
         y(i+1) = u(i)
         y(i+5) = 2.d0*us(i)
      end do

      Br = y(5)*y(6)-y(4)*y(7)+y(3)*y(8)-y(2)*y(9)

      sie = dsin(E/2.d0)
      cse = dcos(E/2.d0)
      iE = dint(E)

c     Write first line output
      write(3,100) t,x,xd,Br
      call car2oe(x,xd,pek)
      write(4,101) iE,y
      write(5,102) iE,pek

      
c     Integration
      do ik = 1,nrev            !Revolutions
         do ki = 1,istep        !Steps
 60         continue
            call u2uu(u,uu)

c$$$            do i = 1,4
c$$$               if(i.eq.1.or.i.eq.3) q_til(i)=2.d0*u(1)*u(3)
c$$$               if(i.eq.2.or.i.eq.4) q_til(i)=2.d0*u(2)*u(4)
c$$$               u_til(i) = 2.d0*u(i)**2
c$$$               qu(i) = u(i)*(q_til(i)-x(3))+uu(i)*(R-u_til(i))
c$$$            end do

            Tau_geo = 0.d0
!     Adding geo-potential to time element
            do i = 2,ngeo_deg
               tau_sum(i) = -(i-1)*c_j(i)*Re(i)*aleg(i+1)*ObyR(i)
               Tau_geo = Tau_geo+tau_sum(i)
            end do

!     Adding third-body potential to time element
            Tau_3body1 = 0.d0
            do i = 2,nsun_deg
               tau_sum(i) = (i+2)*sleg(i+1)*RbyRs(i)*RbyRs(1)*amuS
               Tau_3body1 = Tau_3body1+tau_sum(i)
            end do
            Tau_3body2 = 0.d0            
            do i = 2,nmoon_deg
               tau_sum(i) = (i+2)*oleg(i+1)*RbyRm(i)*RbyRm(1)*amuM
               Tau_3body2 = Tau_3body2+tau_sum(i)
            end do            

!     Time element
           z(1) = (amue*(1+Tau_term)+Tau_3body1+Tau_3body2)/(8.d0*w*w*w)
!     Energy element
            z(10) = 0.d0

            do i = 1,ngeo_deg
               aleg(i)=0.d0
            end do
            call aLegP(ngeo_deg,ZbyR(1),aleg)
            do i = 1,nsun_deg
               sleg(i)=0.d0
            end do
            call aLegP(nsun_deg,cphi_s,sleg)
!            print *,sleg
            do i = 1,nmoon_deg
               oleg(i)=0.d0
            end do
            call aLegP(nmoon_deg,cphi_m,oleg)
c            print *,oleg
c            stop
            do j = 1,4
               q(j) = 0.d0
               qsun(j) = 0.d0
               qmoon(j) = 0.d0
            end do
            
            do j = 1,4

!     oblateness terms
               z2mR2=x(3)**2-R(1)**2
               qj=0.d0
               do i = 2,ngeo_deg
                  qj = dble(ObyR(i+1)*Re(i)*c_j(i)*
     &                 (R(1)*(i+1)*aleg(i+2)*(R(1)*uu(j)-x(3)*u(j))+
     &                 (aleg(i+1)*(u(j)*(i*R(2)+x(3)**2)
     &                 -uu(j)*(i+1)*R(1)*x(3))))/z2mR2)
                  q(j)=qj+q(j)
               end do

!     third-body terms
               call u2qu(u,qu,j)

               den_s = cphi_s2 - 1.d0
               den_m = cphi_m2 - 1.d0

               sq = ts(1)*qu(1)+ts(2)*qu(2)+ts(3)*qu(3)
               
               do i = 2,nsun_deg
                  slambda(i)=2.d0*R(i-1)*(i+1)*(sleg(i+2)
     &                 -cphi_s*sleg(i+1))/den_s
                  
                  qs = (RbyRs(i)*(sq*slambda(i)/Rs(1)+
     &                 u(j)*(1-i)*sleg(i+1)-slambda(i)*cphi_s))
                  qsun(j) = qs + qsun(j)
               end do

               oq = tm(1)*qu(1)+tm(2)*qu(2)+tm(3)*qu(3)
               
               do i = 2,nmoon_deg
                  olambda(i)=2.d0*R(i-1)*(i+1)*(oleg(i+2)
     &                 -cphi_m*oleg(i+1))/den_m
                  qm = (RbyRm(i)*(oq*olambda(i)/Rm(1)+
     &                 u(j)*(1-i)*oleg(i+1)-olambda(i)*cphi_m))
                  qmoon(j) = qm + qmoon(j)
               end do
               
!.................to be filled............

!     drag terms
          
!.................to be filled............
               
!     State elements

               aqt = sie * 0.25d0
               bqt = cse * 0.25d0               
               
               z(j+1) = (q(j)*amue+first_s*qsun(j)+
     &              first_m*qmoon(j))*aqt/(w*w)
               z(j+5) = -(q(j)*amue+first_s*qsun(j)+
     &              first_m*qmoon(j))*bqt/(w*w)
               
            end do
            
            I_rkg = rkgil(n,y,z,E,dE*Gam,nt)
            call ks2ksr(y,u,us,E,cse,sie)
            call ks2car(u,us,x,xd,w)
            call car2oe(x,xd,pek)

            w = y(10)
            T = y(1) - dotp4(u,us)/w
           
            R(1) = dotp4(u,u)

            do i = 2, 20
               R(i) = R(i-1) * R(1)
            end do
            
            V = vmn(xd)
            ZbyR(1) = x(3)/R(1)
            ObyR(1) = 1.d0/R(1)
            do i = 1,ngeo_deg
               ZbyR(i+1) = ZbyR(i)*ZbyR(1)
               ObyR(i+1) = ObyR(i)*ObyR(1)
            end do
            do i = 2,ngeo_deg
               cn_j(i) = amue*Re(i)*c_j(i)*ObyR(i+1)
            end do
            V_pot = 0.d0
            call aLegP(ngeo_deg,ZbyR(1),aleg)
            do i = 2,ngeo_deg
               V_p(i) = cn_j(i)*aleg(i+1)
               V_pot = V_pot+V_p(i)
            end do

            FracJd = T/86400.d0
            djulian1 = djulian + FracJd
            call solarnpv(djulian1, ts)
            call lunarpv(djulian1, tm)
            
            Rs(1) = vmn(ts)
            Rm(1) = vmn(tm)
            RbyRs(1) = R(1)/Rs(1)
            RbyRm(1) = R(1)/Rm(1)
            first_s = -amuS/Rs(1)
            first_m = -amuM/Rm(1)     

            do i = 1, nsun_deg
               RbyRs(i+1) = RbyRs(i)*RbyRs(1)
            end do
            do i = 1, nmoon_deg
               RbyRm(i+1) = RbyRm(i)*RbyRm(1)         
            end do
            
            ph_s = 0.d0
            ph_m = 0.d0
            do i = 1,3
               ph_s = ph_s + x(i) * ts(i)
               ph_m = ph_m + x(i) * tm(i)
            end do

            cphi_s = ph_s/(R(1) * Rs(1))
            cphi_m = ph_m/(R(1) * Rm(1))

            cphi_s2 = cphi_s  * cphi_s
            cphi_m2 = cphi_m  * cphi_m            
            
            call aLegP(nsun_deg,cphi_s,sleg)
            call aLegP(nmoon_deg,cphi_m,oleg)      
            
            do i = 2,nsun_deg
               cn_s(i) = first_s * RbyRs(i)
               V_p(i) = cn_s(i) * sleg(i+1)
               V_pot = V_pot + V_p(i)
            end do

            do i = 2,nmoon_deg
               cn_m(i) = first_m * RbyRm(i)
               V_p(i) = cn_m(i) * oleg(i+1)
               V_pot = V_pot + V_p(i)
            end do
            
            Br = y(5)*y(6)-y(4)*y(7)+y(3)*y(8)-y(2)*y(9)         
            if(I_rkg-1) 60,60,70
 70         continue
            E_0 = E_0+dE_0
            print *,x,xd
c            print *,u,us
c            print *,Br
c            print *,z
c            stop
            iE = dint(E_0*r2d)

         end do

!     Dumping Trajectory by revolution      
         write(3,100) t,x,xd,Br
         write(4,101) iE,y
         write(5,102) iE,pek
      end do
      write(*,*) x
      write(*,*) xd
      write(*,*) pek

c     format specifiers (100 - 200)
 100  format(f20.8,1x,6f18.10,e20.5) !State Vector
 101  format(10f20.9)           !KS-regular elements
 102  format(i10.3,f18.10,f14.10,3f20.8,2f18.10) !Keplerian Elements

 99   do i = 1,5
         close(i)
      end do
      stop
      end
