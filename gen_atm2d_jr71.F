C gen_atm2d_jr71.F - Generate ATM2D.DAT: Jacchia-71 density over a
C Tinf grid (OREM issue #26, epoch-resolved space weather)
C
C The J71 static profile is driven by a single parameter, the
C exospheric temperature Tinf. This generator evaluates the same
C profile functions the 1-D generator uses (jr71_profile.F) over a
C grid of Tinf values, producing a 2-D table rho(alt, Tinf) and
C H(alt, Tinf) that the propagator can interpolate per revolution
C with Tinf(t) computed from daily F10.7/Kp history:
C   Tc   = 379 + 3.24*F10.7B(81d) + 1.3*(F10.7 - F10.7B)
C   Tinf = Tc + 14*Kp + 0.02*exp(Kp)          (J71, h<200 km form)
C
C Grid: Tinf = 550(25)1500 K (39 columns) -- covers quiet nighttime
C minimum (~590 K at F10.7=65) through storm-time solar maximum
C (~1400+ K at F10.7~250, Kp~9). Altitudes: the standard KSROP
C 291-point grid (60-200 km at 1 km, 202-500 km at 2 km).
C SCH column = local density scale height -dz/dln(rho) per column.
C
C File format (input/ATM2D.DAT):
C   line 1 : 'ATM2D_JR71' nalt nT
C   then   : nT Tinf values            (8F10.2)
C   then, for each Tinf column j=1..nT:
C            291 SCH values            (8F10.3)  [km]
C            291 DEN values            (8E12.5)  [kg/m3]
C
      program gen_atm2d_jr71
      implicit none

      integer nalt, nT
      parameter (nalt = 291, nT = 39)

      double precision ALT(nalt), SCH(nalt), DEN(nalt), Tg(nT)
      double precision rhoj71
      external rhoj71
      double precision rm, rp
      integer i, j

      do j = 1, nT
        Tg(j) = 550.d0 + 25.d0*dble(j-1)
      end do

      ALT(1) = 60.d0
      do i = 2, 141
        ALT(i) = ALT(i-1) + 1.d0
      end do
      do i = 142, nalt
        ALT(i) = ALT(i-1) + 2.d0
      end do

      open(unit=10, file='input/ATM2D.DAT', status='replace')
      write(10,'(A,2I6)') 'ATM2D_JR71', nalt, nT
      write(10,'(8F10.2)') (Tg(j), j=1,nT)

      do j = 1, nT
        do i = 1, nalt
          DEN(i) = rhoj71(ALT(i), Tg(j))
          rm     = rhoj71(ALT(i) - 0.5d0, Tg(j))
          rp     = rhoj71(ALT(i) + 0.5d0, Tg(j))
          SCH(i) = 1.d0 / (dlog(rm) - dlog(rp))
        end do
        write(10,'(8F10.3)')  (SCH(i), i=1,nalt)
        write(10,'(8E12.5)') (DEN(i), i=1,nalt)
      end do
      close(10)

      write(*,'(A,I4,A,I4,A)') ' ATM2D.DAT written (', nalt,
     +   ' altitudes x', nT, ' Tinf columns)'
      write(*,'(A,F8.1,A,F8.1,A)') ' Tinf grid: ', Tg(1), ' ..',
     +   Tg(nT), ' K'
      write(*,'(A,E12.4,A,E12.4)')
     +   ' rho(170 km): T=550: ', rhoj71(170.d0, 550.d0),
     +   '   T=1500: ', rhoj71(170.d0, 1500.d0)

      stop
      end
