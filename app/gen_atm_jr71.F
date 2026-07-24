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
