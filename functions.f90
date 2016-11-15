!equation list file
!by Murilo Varela
real*8 function f1_1 (t, ti)
      !t: time
      !ti: residence time

      implicit none

      real*8 :: t, ti

      f1_1 = t / ti

end function f1_1

real*8 function f1_2 (rg, b, vm, k, c, n)
      !rg: grains radius of metalic copper
      !b: stoichiometric coefficient
      !vm: molar volume of CuO
      !k: constant in arrhenius form (function f2)
      !c: concentration (vol%)
      !n: reaction order

      implicit none

      real*8 :: rg, b, vm, k, c, n

      f1_2 = rg / (b * vm * k * c**n)
end function f1_2

real*8 function f2 (k0, e, r, temp)
      !k0: pre exponential factor
      !e: activation energy
      !r: constant for ideal gases
      !t: temperature

      implicit none

      real*8 :: k0, e, r, temp

      f2 = k0 * exp(-e / (r * temp))

end function f2

real*8 function f30 (bg, nc, phi)
      !bg: stoichiometric coefficient to fully convert CH4 into CO2
      !nc: efficiency of combustion
      !phi: ratio of oxygen-carrier to fuel

      implicit none

      real*8 :: bg, nc, phi

      f30 = bg * nc * (1.0 / (4.0 * phi))

end function f30

real*8 function f27 (tmr, t)
      !tmr: mean residence time of particle in the whole fluidized bed reactor
      !t: time

      implicit none

      real*8 :: tmr, t

      f27 = (1.0 / tmr) * exp(-t / tmr)

end function f27

real*8 function f28 (xsin, tm, t)
      !xsin: mean convertion of particles in
      !tm: mean reacting time
      !t: time

      implicit none

      real*8 :: xsin, tm, t

      f28 = xsin + t / tm

end function f28

real*8 function f29 (xsout, tm, tmr, t, xsin)
      !xsout: mean convertion of particles going out
      !tm: mean reacting time
      !tmr: mean residence time of particle in the whole fluidized bed reactor
      !t: time
      !xsin: mean convertion of particles going in

      implicit none

      real*8 :: xsout, tm, tmr, t, xsin, f27, f28

      f29 = (1.0 - f28 (xsin, tm, t)) * f27 (tmr, t)
end function f29

real*8 function f29_i (xsout, tm, tmr, xsin)
      implicit none

      real*8 :: xsout, tm, tmr, xsin
      real*8 :: t1, t2, t3, t4, h, sol, a, b, f29
      integer :: i

      a = 0.0
      b = tm
      h = (b - a) / dble(3 * 10)
      t1 = a
      sol = 0.0

      do i=1,10,1
            t2 = t1 + h
            t3 = t1 + 2.0*h
            t4 = t1 + 3.0*h
            sol = sol   + (f29 (xsout, tm, tmr, t1, xsin) &
                        + 3.0 * f29 (xsout, tm, tmr, t2, xsin) &
                        + 3.0 * f29 (xsout, tm, tmr, t3, xsin) &
                        + f29 (xsout, tm, tmr, t4, xsin)) * 3.0 * h / 8.0
            t1 = t4
      end do

      f29_i = sol


end function f29_i
