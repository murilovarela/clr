!combustion looping reactor
!by: Murilo Varela

!https://drive.google.com/file/d/0B6-bRnYhK1rBdWNwd0NGNU5VcWs/view?usp=sharing

include 'kinetic_parameters.f90'
include 'solid_convertion.f90'
include 'fluid_dynamics.f90'
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
      !step 3
      real*8 :: c1, c2, ar, rhog, mdp, mug, umf
      real*8 :: u, g, h, a0, ubinf
      real*8 :: dp, u0, sigb
      real*8 :: cs, z, csb, ug

!step 0
!determination of kinetic parameters
      k0 = 480.0              !mol1-n m3n-2 s-1
      e = 95116.4             !J/mol
      r = 8.314462            !J K−1 mol−1
      t = 5.0                 !s
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
      tmr = ti                !s

      call rtd_curve (tmr, t, rtd)

!           solid convertion distribution (eq 28,29)
      xsin = 0.0              !%
      xsout = xs1 - xsin      !%
      call mean_react_time (xsout, tm, tmr, xsin)

!step 3
!     reactor model
!           fluid dynamics
!                 dense bed (qe 8-17)
      c1 = 27.2
      c2 = 0.0408
      ar = 1.0                !must find a way to find this value
      rhog = 1.82              !must find a way to find this value
      mdp = 200.0e-6          !must be checked
      mug = 4.46e-6               !must find a way to find this value

      u = 0.075                                 !0.075–0.15 m/s
      a0 = (3.14/4.0) * (0.1**2.0)        !must be checked
      g = 9.8                                   !m/s2
      h = 1.2                             !value to vary

      dp = 2.0e-7                            !must be checked

      call velocity_min_fluidization (c1, c2, ar, rhog, mdp, mug, umf)
      call visible_velocity_bubble (u, g, dp, umf, h, a0, u0, sigb, ubinf)

!                 freeboard (eq 18-19)
      csb = 1.0                           !must find out what it mean
      z = 0.9
      ug = u

      call fluid_dynamics_in_freeboard (cs, z, csb, ug, dp)
      print*, cs, z
!           mass balance (eq 23-26, 31)

!step 4
!     solid convertion calculation (eq 30) Xs2
!step 5
!     if Xs1 = Xs2: go to step 6; else: go to step 1
!step 6
!     output data

end program main
