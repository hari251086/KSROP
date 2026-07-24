C jr71_profile.F - Jacchia-71 static profile functions (shared)
C   tempz(z,Tinf)  : temperature at altitude z (km) for exospheric Tinf
C   rhoj71(z,Tinf) : total mass density (kg/m3)
C Extracted from gen_atm_jr71.F (2026-07-17) so both the 1-D table
C generator (gen_atm_jr71) and the 2-D T-inf-gridded generator
C (gen_atm2d_jr71, OREM issue #26) link one implementation.
C Coefficients: Roberts (1971) anchors, validated 0.80-0.95x GMAT
C JacchiaRoberts over 102-300 km (OREM v1.17).

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
      interface
        double precision function tempz(z, Tinf)
        double precision z, Tinf
        end function tempz
      end interface

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
