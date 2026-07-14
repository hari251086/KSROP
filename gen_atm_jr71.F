C gen_atm_jr71.F - Generate ATM.DAT using the Jacchia-71 static profile
C 291-point table: 60-200 km at 1 km, 202-500 km at 2 km (KSROP grid)
C
C Replaces gen_atm_j70.F, whose hand-rolled single-exponential
C temperature profile (12 K/km forced gradient at 90 km) ran ~127 K too
C warm through 90-125 km and inflated the hydrostatic column above it:
C the resulting table was ~3.3-3.5x denser than Jacchia-Roberts in the
C 140-200 km perigee band (OREM issue #12, v1.17 investigation), which
C drove OREM's re-entry predictions ~4-5x early.
C
C This version implements the real Jacchia-71 structure, using the same
C polynomial anchors as the Jacchia-Roberts 1971 analytic formulation
C (Roberts 1971; coefficient set cross-checked against the
C SatelliteToolboxAtmosphericModels.jl reference implementation and
C validated to 0.6-6% against GMAT R2026a's JacchiaRoberts over
C 102-300 km):
C   T(z), 90-125 km : quartic with T(90)=183 K, dT/dz(90)=0,
C                     inflection at 125 km, T(125)=Tx(Tinf)
C   T(z), >125 km   : exponential-asymptotic toward Tinf
C   rho, >=125 km   : closed-form species densities anchored at the
C                     delta_ij(Tinf) polynomials (log10 n_i at 125 km)
C   rho, 100-125 km : per-species diffusive equilibrium integrated up
C                     from the zeta(Tinf) total-density anchor at 100 km
C                     (He gets thermal-diffusion exponent alpha=-0.38)
C   rho, 90-100 km  : barometric with the Aa mean-molecular-mass
C                     polynomial, anchored at rho(90)=3.46e-9 g/cm3
C   rho, <90 km     : isothermal continuation at 183 K (drag negligible)
C
C Conditions: F10.7 = F10.7B = 72, Kp = 1.0 (same convention as before).
C Tinf = 379 + 3.24*F10.7B + 1.3*(F10.7-F10.7B) + 14*Kp + 0.02*exp(Kp)
C      = 626.3 K -- nighttime-minimum exospheric temperature plus the
C J71 geomagnetic term for h<200 km (the drag-critical band). A static
C table cannot carry J71's diurnal factor (1.0-1.3 on Tc, up to ~30%
C density at these altitudes); that scale is absorbed by OREM's fitted
C BN, and only the profile *shape* drives the re-entry arc.
C
C SCH column = local density scale height H = -dz/d(ln rho), evaluated
C numerically from the generated profile (NOT R*T/(M*g)) -- this is the
C quantity propagate_ks's King-Hele drag model actually needs at
C perigee.
C
      program gen_atm_jr71
      implicit none

      double precision Rstar, Av, Ra, g0, M0, T1, z1, M1, rho1
      double precision z2, zx
      parameter (Rstar = 8.31432d0)
      parameter (Av    = 6.022045d23)
      parameter (Ra    = 6356.766d0)
      parameter (g0    = 9.80665d0)
      parameter (M0    = 28.960d0)
      parameter (T1    = 183.d0)
      parameter (z1    = 90.d0)
      parameter (M1    = 28.82678d0)
      parameter (rho1  = 3.46d-9)
      parameter (z2    = 100.d0)
      parameter (zx    = 125.d0)

      double precision F10, F10B, Kp, Tinf, Tx
      parameter (F10 = 72.d0, F10B = 72.d0, Kp = 1.0d0)

      integer ndim
      parameter (ndim = 291)
      double precision ALT(ndim), SCH(ndim), DEN(ndim)

      double precision rhoj71, tempz
      external rhoj71, tempz

      double precision rm, rp
      integer i

C     Nighttime-minimum exospheric temperature + J71 h<200 storm term
      Tinf = 379.d0 + 3.24d0*F10B + 1.3d0*(F10 - F10B)
     +     + 14.d0*Kp + 0.02d0*dexp(Kp)

      Tx = 371.6678d0 + 0.0518806d0*Tinf
     +   - 294.3505d0*dexp(-0.00216222d0*Tinf)

      write(*,'(a,f8.2,a)') '  Tinf =', Tinf, ' K  (F10.7=72, Kp=1.0)'
      write(*,'(a,f8.2,a)') '  Tx   =', Tx,   ' K  (inflection, 125 km)'
      write(*,'(a,f8.2,a)') '  T(100)=', tempz(100.d0, Tinf), ' K'

C     Altitude grid (km)
      ALT(1) = 60.d0
      do i = 2, 141
        ALT(i) = ALT(i-1) + 1.d0
      end do
      do i = 142, ndim
        ALT(i) = ALT(i-1) + 2.d0
      end do

      do i = 1, ndim
        DEN(i) = rhoj71(ALT(i), Tinf)
C       local density scale height from centered ln-density difference
        rm     = rhoj71(ALT(i) - 0.5d0, Tinf)
        rp     = rhoj71(ALT(i) + 0.5d0, Tinf)
        SCH(i) = 1.d0 / (dlog(rm) - dlog(rp))
      end do

C     Summary
      write(*,'(/,a8,2a14)') 'Alt(km)', 'DEN(kg/m3)', 'SCH(km)'
      write(*,'(i8,2e14.4)') 60,  DEN(1),   SCH(1)
      write(*,'(i8,2e14.4)') 90,  DEN(31),  SCH(31)
      write(*,'(i8,2e14.4)') 120, DEN(61),  SCH(61)
      write(*,'(i8,2e14.4)') 150, DEN(91),  SCH(91)
      write(*,'(i8,2e14.4)') 170, DEN(111), SCH(111)
      write(*,'(i8,2e14.4)') 200, DEN(141), SCH(141)
      write(*,'(i8,2e14.4)') 300, DEN(191), SCH(191)
      write(*,'(i8,2e14.4)') 500, DEN(291), SCH(291)

C     Write ATM.DAT
      open(unit=10, file='input/ATM.DAT', status='replace')
      write(10, 100) (SCH(i), i=1, 61)
      write(10, 100) (SCH(i), i=62, ndim)
      write(10, 101) (DEN(i), i=1, 61)
      write(10, 101) (DEN(i), i=62, ndim)
      close(10)

      write(*,'(/,a)') 'ATM.DAT written to input/ATM.DAT'

  100 format(8f10.3)
  101 format(8e10.3)

      stop
      end

C     ==================================================================
C     T(z): J71 quartic (90-125 km) / exponential-asymptotic (>125 km)
C     Below 90 km: held at T1 (isothermal continuation)
C     ==================================================================
      double precision function tempz(z, Tinf)
      implicit none
      double precision z, Tinf
      double precision T1, z1, zx, Ra
      parameter (T1 = 183.d0, z1 = 90.d0, zx = 125.d0)
      parameter (Ra = 6356.766d0)
      double precision Ca(5), la(5)
      data Ca / -89284375.0d0, 3542400.0d0, -52687.5d0,
     +          340.5d0, -0.8d0 /
      data la / 0.1031445d5, 0.2341230d1, 0.1579202d-2,
     +          -0.1252487d-5, 0.2462708d-9 /
      double precision Tx, p, al
      integer k

      Tx = 371.6678d0 + 0.0518806d0*Tinf
     +   - 294.3505d0*dexp(-0.00216222d0*Tinf)

      if (z .le. z1) then
        tempz = T1
      else if (z .le. zx) then
        p = 0.d0
        do k = 5, 1, -1
          p = p*z + Ca(k)
        end do
        tempz = Tx + (Tx - T1)/35.d0**4 * p
      else
        al = 0.d0
        do k = 5, 1, -1
          al = al*Tinf + la(k)
        end do
        tempz = Tinf - (Tinf - Tx)*dexp(
     +      -al * ((Tx - T1)/(Tinf - Tx))
     +          * ((z - zx)/(zx - z1)) / (Ra + z) )
      end if

      return
      end

C     ==================================================================
C     rho(z): total mass density (kg/m3) for the static J71 profile
C     ==================================================================
      double precision function rhoj71(z, Tinf)
      implicit none
      double precision z, Tinf
      double precision tempz
      external tempz

      double precision Rstar, Av, Ra, g0, M0, T1, z1, M1, rho1
      double precision z2, zx
      parameter (Rstar = 8.31432d0)
      parameter (Av    = 6.022045d23)
      parameter (Ra    = 6356.766d0)
      parameter (g0    = 9.80665d0)
      parameter (M0    = 28.960d0)
      parameter (T1    = 183.d0)
      parameter (z1    = 90.d0)
      parameter (M1    = 28.82678d0)
      parameter (rho1  = 3.46d-9)
      parameter (z2    = 100.d0)
      parameter (zx    = 125.d0)

C     Species: 1=N2, 2=O2, 3=O, 4=Ar, 5=He
      double precision Mi(5), ai(5), mui(5)
      data Mi  / 28.0134d0, 31.9988d0, 15.9994d0, 39.9480d0, 4.0026d0 /
      data ai  / 0.d0, 0.d0, 0.d0, 0.d0, -0.38d0 /
      data mui / 0.78110d0, 0.161778d0, 0.095544d0, 0.0093432d0,
     +           0.61471d-5 /

      double precision Aa(7), la(5), zet(7), dij(7,5)
      data Aa / -435093.363387d0, 28275.5646391d0, -765.33466108d0,
     +          11.043387545d0, -0.08958790995d0, 0.00038737586d0,
     +          -0.000000697444d0 /
      data la / 0.1031445d5, 0.2341230d1, 0.1579202d-2,
     +          -0.1252487d-5, 0.2462708d-9 /
      data zet / 0.1985549d-10, -0.1833490d-14, 0.1711735d-17,
     +           -0.1021474d-20, 0.3727894d-24, -0.7734110d-28,
     +           0.7026942d-32 /
C     log10 n_i(125 km) polynomials in Tinf (Roberts 1971)
      data dij /
     +  0.1093155d2,  0.1186783d-2, -0.1677341d-5,  0.1420228d-8,
     + -0.7139785d-12, 0.1969715d-15, -0.2296182d-19,
     +  0.9924237d1,  0.1600311d-2, -0.2274761d-5,  0.1938454d-8,
     + -0.9782183d-12, 0.2698450d-15, -0.3131808d-19,
     +  0.1097083d2,  0.6118742d-4, -0.1165003d-6,  0.9239354d-10,
     + -0.3490739d-13, 0.5116298d-17, 0.d0,
     +  0.8049405d1,  0.2382822d-2, -0.3391366d-5,  0.2909714d-8,
     + -0.1481702d-11, 0.4127600d-15, -0.4837461d-19,
     +  0.7646886d1, -0.4383486d-3,  0.4694319d-6, -0.2894886d-9,
     +  0.9451989d-13, -0.1270838d-16, 0.d0 /

      double precision Tx, Tz, Txz, T100, al, gam, gi
      double precision rho, r125, rho100, xint, dz, f0, f1, zj
      double precision gz, p, Mz
      integer k, isp, nstep, j

      Tx = 371.6678d0 + 0.0518806d0*Tinf
     +   - 294.3505d0*dexp(-0.00216222d0*Tinf)
      Tz = tempz(z, Tinf)

      if (z .ge. zx) then
C       ---- closed-form species densities anchored at 125 km ----
        al = 0.d0
        do k = 5, 1, -1
          al = al*Tinf + la(k)
        end do
        gam = g0*Ra*Ra*1.d3 / (Rstar*al*Tinf)
     +      * (Tinf - Tx)/(Tx - T1) * (zx - z1)/(Ra + zx)
        Txz = tempz(zx, Tinf)
        rho = 0.d0
        do isp = 1, 5
          p = 0.d0
          do k = 7, 1, -1
            p = p*Tinf + dij(k,isp)
          end do
          r125 = Mi(isp) * 10.d0**p / Av
          gi   = gam * Mi(isp) * 1.d-3
          rho  = rho + r125 * (Tx/Tz)**(1.d0 + ai(isp) + gi)
     +               * ((Tinf - Tz)/(Tinf - Txz))**gi
        end do
        rhoj71 = rho * 1.d3

      else if (z .ge. z2) then
C       ---- diffusive equilibrium up from the 100 km anchor ----
        T100 = tempz(z2, Tinf)
        p = 0.d0
        do k = 7, 1, -1
          p = p*Tinf + zet(k)
        end do
        rho100 = p * M0
C       hydrostatic integral of g/(Rstar*T) dz' (SI) from 100 to z
        nstep = max(1, nint((z - z2)/0.05d0))
        dz    = (z - z2)/dble(nstep)
        xint  = 0.d0
        f0 = g0*(Ra/(Ra + z2))**2 / (Rstar*T100)
        do j = 1, nstep
          zj = z2 + dble(j)*dz
          f1 = g0*(Ra/(Ra + zj))**2 / (Rstar*tempz(zj, Tinf))
          xint = xint + 0.5d0*(f0 + f1)*dz*1.d3
          f0 = f1
        end do
        rho = 0.d0
        do isp = 1, 5
          rho = rho + rho100*mui(isp)*Mi(isp)/M0*(Av/Mi(isp))
     +              * Mi(isp)/Av
     +              * (T100/Tz)**(1.d0 + ai(isp))
     +              * dexp(-Mi(isp)*1.d-3*xint)
        end do
        rhoj71 = rho * 1.d3

      else if (z .ge. z1) then
C       ---- barometric with Aa mean-molecular-mass poly, 90-100 ----
        Mz = 0.d0
        do k = 7, 1, -1
          Mz = Mz*z + Aa(k)
        end do
        nstep = max(1, nint((z - z1)/0.05d0))
        dz    = (z - z1)/dble(nstep)
        xint  = 0.d0
        f0 = M1*1.d-3 * g0*(Ra/(Ra + z1))**2 / (Rstar*T1)
        do j = 1, nstep
          zj = z1 + dble(j)*dz
          p = 0.d0
          do k = 7, 1, -1
            p = p*zj + Aa(k)
          end do
          f1 = p*1.d-3 * g0*(Ra/(Ra + zj))**2
     +       / (Rstar*tempz(zj, Tinf))
          xint = xint + 0.5d0*(f0 + f1)*dz*1.d3
          f0 = f1
        end do
        rhoj71 = rho1 * (Mz*T1)/(M1*Tz) * dexp(-xint) * 1.d3

      else
C       ---- below 90 km: isothermal continuation at T1 ----
        nstep = max(1, nint((z1 - z)/0.05d0))
        dz    = (z1 - z)/dble(nstep)
        xint  = 0.d0
        f0 = M1*1.d-3 * g0*(Ra/(Ra + z1))**2 / (Rstar*T1)
        do j = 1, nstep
          zj = z1 - dble(j)*dz
          gz = g0*(Ra/(Ra + zj))**2
          f1 = M1*1.d-3 * gz / (Rstar*T1)
          xint = xint + 0.5d0*(f0 + f1)*dz*1.d3
          f0 = f1
        end do
        rhoj71 = rho1 * dexp(xint) * 1.d3
      end if

      return
      end
