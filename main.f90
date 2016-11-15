!combustion looping reactor
!by: Murilo Varela
include 'kinetic_parameters.f90'
include 'solid_convertion.f90'
include 'functions.f90'

program main
      implicit none
      !step 0
      real*8 :: k0, e, r, t, rg, b, vm, c, n, temp, xs, ti
      !step 1
      real*8 :: bg, nc, phi, xs1
      !step 2
      real*8 :: rtd, tmr
      real*8 :: xsin, xsout, tm

!step 0
!determination of kinetic parameters
      k0 = 480.0              !mol1-n m3n-2 s-1
      e = 95116.4             !J/mol
      r = 8.314462            !J K−1 mol−1
      t = 5.0                !s
      rg = 1.4e-6             !m
      b = 4.0                 !mol/mol
      vm = 12.4e-6            !m3/mol
      c = 10.0 / 100.0         !vol%
      n = 0.5                 ! -
      temp = 1073.0           !k

      call kinetic_parameters (k0, e, r, t, rg, b, vm, c, n, temp, xs, ti)

!steps 1
!     solid convertion assumption (eq 30) Xs1
      bg =  4.0               ! -
      nc =  1.0               ! -
      phi = 2.2               ! - <0.7 - 2.2>

      call convertion_assumption (bg, nc, phi, xs1)

!step 2
!     average solids convertion
!           RTD curve, E(t) (eq 27)
      tmr = ti

      call rtd_curve (tmr, t, rtd)

!           solid convertion distribution (eq 28,29)
      xsin = 0.0
      xsout = xs1 - xsin
      call mean_react_time (xsout, tm, tmr, xsin)
      print*, tm
!step 3
!     reactor model
!           fluid dynamics
!                 dense bed (qe 8-17)
!                 freeboard (eq 18-19)
!           mass balance (eq 23-26, 31)
!step 4
!     solid convertion calculation (eq 30) Xs2
!step 5
!     if Xs1 = Xs2: go to step 6; else: go to step 1
!step 6
!     output data

end program main
