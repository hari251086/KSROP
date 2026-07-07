C gen_atm_j70.F - Generate ATM.DAT using Jacchia-70 atmosphere model
C 291-point table: 60-200 km at 1 km, 202-500 km at 2 km (KSROP grid)
C Reference: Jacchia (1970, 1971), Roberts (1971)
C Conditions: F10.7=72, F10.7B=72, Kp=1.0  (matches NPOE Zone-0 test)
C
C Species: N2, O2, O, Ar, He -- diffusive equilibrium above 90 km
C Temperature: T(z) = Tinf + (T0-Tinf)*exp(-(z-90)/lam), z in km
C   lam chosen so dT/dz = 12 K/km at z=90 km
C Below 90 km: isothermal at T0=183 K (mesopause); drag negligible there
C He thermal diffusion: extra exponent 0.4 on temperature ratio
C
      program gen_atm_j70
      implicit none

C     Physical constants
      double precision Rgas, NA, g0, RE
      parameter (Rgas = 8.314462d0)
      parameter (NA   = 6.022141d23)
      parameter (g0   = 9.80665d0)
      parameter (RE   = 6.356766d6)

C     Molecular masses (kg/mol)
      double precision mN2, mO2, mO, mAr, mHe
      parameter (mN2 = 0.028014d0)
      parameter (mO2 = 0.031999d0)
      parameter (mO  = 0.015999d0)
      parameter (mAr = 0.039948d0)
      parameter (mHe = 0.004003d0)

C     Reference number densities at z0=90 km (m^-3)
C     From US Standard Atmosphere 1976: P=0.1685 Pa, T=186.87 K
C     n_total=6.53e19 m^-3, well-mixed: 78.09% N2, 20.95% O2, 0.93% Ar
      double precision nN2_0, nO2_0, nO_0, nAr_0, nHe_0
      parameter (nN2_0 = 5.100d19)
      parameter (nO2_0 = 1.368d19)
      parameter (nO_0  = 1.000d15)
      parameter (nAr_0 = 6.073d17)
      parameter (nHe_0 = 3.420d14)

C     Solar/geomagnetic inputs
      double precision F10, F10B, Kp
      parameter (F10  = 72.d0)
      parameter (F10B = 72.d0)
      parameter (Kp   = 1.0d0)

      integer ndim
      parameter (ndim = 291)

      double precision ALT(ndim), SCH(ndim), DEN(ndim)
      double precision Tinf, T0, lam, z0
      double precision h, T, zh, Th, gh, dz
      double precision intN2, intO2, intO, intAr, intHe
      double precision f0N2, f0O2, f0O, f0Ar, f0He
      double precision f1N2, f1O2, f1O, f1Ar, f1He
      double precision nN2, nO2, nO, nAr, nHe
      double precision ntot, mbar, rho, zm, g
      integer i, j, nstep

C     Exospheric temperature
      Tinf = 379.d0 + 3.24d0*F10B + 1.3d0*(F10-F10B)
     +     + 28.d0*Kp + 0.03d0*dexp(Kp)

C     Mesopause boundary
      T0 = 183.d0
      z0 = 90.d0

C     Temperature scale (km): from 12 K/km gradient at z0
      lam = (Tinf - T0) / 12.d0

      write(*,'(a,f8.2,a)') '  Tinf =', Tinf, ' K  (F10.7=72, Kp=1.0)'
      write(*,'(a,f8.2,a)') '  lam  =', lam,  ' km'

C     Altitude grid (km)
      ALT(1) = 60.d0
      do i = 2, 141
        ALT(i) = ALT(i-1) + 1.d0
      end do
      do i = 142, ndim
        ALT(i) = ALT(i-1) + 2.d0
      end do

C     Loop over all grid points
      do i = 1, ndim
        h = ALT(i)

C       Temperature at h
        if (h .ge. z0) then
          T = Tinf + (T0 - Tinf) * dexp(-(h - z0) / lam)
        else
          T = T0
        end if

C       Hydrostatic integrals from z0 to h (positive above, negative below)
        nstep = max(1, nint(abs(h - z0) / 0.1d0))
        dz    = (h - z0) / dble(nstep)

        intN2 = 0.d0
        intO2 = 0.d0
        intO  = 0.d0
        intAr = 0.d0
        intHe = 0.d0

C       Integrand at z0 (lower bound)
        zh   = z0
        Th   = T0
        gh   = g0 * (RE / (RE + zh*1.d3))**2
        f0N2 = mN2 * gh / (Rgas * Th)
        f0O2 = mO2 * gh / (Rgas * Th)
        f0O  = mO  * gh / (Rgas * Th)
        f0Ar = mAr * gh / (Rgas * Th)
        f0He = mHe * gh / (Rgas * Th)

        do j = 1, nstep
          zh = z0 + dble(j) * dz

C         Temperature at zh
          if (zh .ge. z0) then
            Th = Tinf + (T0 - Tinf) * dexp(-(zh - z0) / lam)
          else
            Th = T0
          end if

          gh   = g0 * (RE / (RE + zh*1.d3))**2
          f1N2 = mN2 * gh / (Rgas * Th)
          f1O2 = mO2 * gh / (Rgas * Th)
          f1O  = mO  * gh / (Rgas * Th)
          f1Ar = mAr * gh / (Rgas * Th)
          f1He = mHe * gh / (Rgas * Th)

C         Trapezoidal step (dz*1000: km->m for correct units)
          intN2 = intN2 + 0.5d0*(f0N2+f1N2) * dz * 1.d3
          intO2 = intO2 + 0.5d0*(f0O2+f1O2) * dz * 1.d3
          intO  = intO  + 0.5d0*(f0O +f1O ) * dz * 1.d3
          intAr = intAr + 0.5d0*(f0Ar+f1Ar) * dz * 1.d3
          intHe = intHe + 0.5d0*(f0He+f1He) * dz * 1.d3

          f0N2 = f1N2
          f0O2 = f1O2
          f0O  = f1O
          f0Ar = f1Ar
          f0He = f1He
        end do

C       Number densities (m^-3): n_i = n_i0 * (T0/T) * exp(-integral)
        nN2 = nN2_0 * (T0/T)        * dexp(-intN2)
        nO2 = nO2_0 * (T0/T)        * dexp(-intO2)
        nO  = nO_0  * (T0/T)        * dexp(-intO )
        nAr = nAr_0 * (T0/T)        * dexp(-intAr)
        nHe = nHe_0 * (T0/T)**1.4d0 * dexp(-intHe)

C       Total mass density (kg/m^3)
        rho = (nN2*mN2 + nO2*mO2 + nO*mO + nAr*mAr + nHe*mHe) / NA

C       Mean molecular mass (kg/mol)
        ntot = nN2 + nO2 + nO + nAr + nHe
        mbar = (nN2*mN2 + nO2*mO2 + nO*mO + nAr*mAr + nHe*mHe) / ntot

C       Local gravity and scale height (km)
        zm     = (RE / (RE + h*1.d3))**2
        g      = g0 * zm
        SCH(i) = Rgas * T / (mbar * g) * 1.d-3
        DEN(i) = rho
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
