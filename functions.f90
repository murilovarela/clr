!equation list file
!by Murilo Varela

real*8 function f1_1 (t, ti)
      !t: time (s)
      !ti: residence time (s)

      implicit none

      real*8 :: t, ti

      f1_1 = t / ti

end function f1_1

real*8 function f1_2 (rg, b, vm, k, c, n)
      !rg: grains radius of metalic copper (mu m)
      !b: stoichiometric coefficient (mol CuO / mol gas)
      !vm: molar volume of CuO (cm3//mol)
      !k: constant in arrhenius form (function f2)
      !c: concentration (vol%)
      !n: reaction order (-)

      implicit none

      real*8 :: rg, b, vm, k, c, n

      f1_2 = rg / (b * vm * k * c**n)
end function f1_2

real*8 function f2 (k0, e, r, temp)
      !k0: pre exponential factor (...)
      !e: activation energy (j/mol)
      !r: constant for ideal gases (j k−1 mol−1)
      !t: temperature (k)

      implicit none

      real*8 :: k0, e, r, temp

      f2 = k0 * exp(-e / (r * temp))

end function f2

real*8 function f30 (bg, nc, phi)
      !bg: stoichiometric coefficient to fully convert CH4 into CO2 (-)
      !nc: efficiency of combustion (%)
      !phi: ratio of oxygen-carrier to fuel (-)

      implicit none

      real*8 :: bg, nc, phi

      f30 = bg * nc * (1.0 / (4.0 * phi))

end function f30

real*8 function f27 (tmr, t)
      !tmr: mean residence time of particle in the whole fluidized bed reactor (s)
      !t: time (s)

      implicit none

      real*8 :: tmr, t

      f27 = (1.0 / tmr) * exp(-t / tmr)

end function f27

real*8 function f28 (xsin, tm, t)
      !xsin: mean convertion of particles in (%)
      !tm: mean reacting time (s)
      !t: time (s)

      implicit none

      real*8 :: xsin, tm, t

      f28 = xsin + t / tm

end function f28

real*8 function f29 (xsout, tm, tmr, t, xsin)
      !xsout: mean convertion of particles going out (%)
      !tm: mean reacting time (s)
      !tmr: mean residence time of particle in the whole fluidized bed reactor (s)
      !t: time (s)
      !xsin: mean convertion of particles going in (%)

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

real*8 function f8 (umf, uvis, utf, sigb)
      !umf: velocity of minimum fluidization (m/s)
      !uvis: visible velocity of bubble (m/s)
      !utf: velocity of the gas throughflow (m/s)
      !sigb: fraction of bubble in the dense bed (-)

      implicit none

      real*8 :: umf, uvis, utf, sigb

      f8 = (1.0 - sigb) * umf + uvis + utf

end function f8

real*8 function f9_1 (c1, c2, ar)
      !c1, c2, c3: constans for reynolds calculation (-)
      !ar: arrhenius number (-)

      implicit none

      real*8 :: c1, c2, ar

      f9_1 = sqrt(c1**2.0 + c2 * ar) - c1

end function f9_1

real*8 function f9_2 (rhog, mdp, mug, re)
      !umf: velocity of minimum fluidization (m/s)
      !rhog: gas density (kg/m3)
      !mdp: mean particle diameter (m)
      !mug: fluidizing gas viscosity (kg/m s)

      implicit none

      real*8 :: rhog, mdp, mug, re

      f9_2 = re * mug / (rhog * mdp)

end function f9_2

real*8 function f10 (psi, u0, umf, sigb)
      !umf: velocity of minimum fluidization (m/s)
      !psi:ratio of the visible bubble flow to the total flow through the bubbles
      !u0: velocity of the total gas flow
      !sigb: fraction of bubble in the dense bed (-)

      implicit none

      real*8 :: psi, u0, umf, sigb

      f10 = psi * (u0 - umf * (1.0 - sigb))

end function f10

real*8 function f11 (psi, u0, umf, sigb)
      !umf: velocity of minimum fluidization (m/s)
      !psi:ratio of the visible bubble flow to the total flow through the bubbles
      !u0: velocity of the total gas flow
      !sigb: fraction of bubble in the dense bed (-)

      implicit none

      real*8 :: psi, u0, umf, sigb

      f11 = (1.0 - psi) * (u0 - umf * (1.0 - sigb))

end function f11

real*8 function f12 (fb, h, a0)
      !fb: empirical function (eq 13)
      !h: height
      !a0: area of the gas-distributor per nozzle, m2 per nozzle

      implicit none

      real*8 :: fb, h, a0

      f12 = fb * (h + 4.0 * sqrt(a0))**0.4

end function f12

real*8 function f13 (dp, u0, umf)
      !dp: diameter of particle (m/s)
      !umf: velocity of minimum fluidization (m/s)
      !u0: velocity of the total gas flow

      implicit none

      real*8 :: dp, u0, umf

      f13 = (0.26 + 0.7 * exp(-3300.0 * dp)) / ((0.15 + u0 - umf)**(1.0/3.0))

end function f13

real*8 function f14 (uvis, ubinf)
      !uvis: visible velocity of bubble (m/s)
      !ubinf: velocity of a single bubble (m/s)

      implicit none

      real*8 :: uvis, ubinf

      f14 = uvis / (uvis + ubinf)

end function f14

real*8 function f15 (g, db)
      !g: gravity acceleration (m/s2)
      !db: diameter of bubble (m)

      implicit none

      real*8 :: g, db

      f15 = 0.71 * sqrt(g * db)

end function f15

real*8 function f16 (sigb, emf)
      !sigb: fraction of bubble in the dense bed (-)
      !emf: porosity at minimum fluidization (-)

      implicit none

      real*8 :: sigb, emf

      f16 = (1.0 - sigb) * emf + sigb

end function f16

real*8 function f17 (ar, rhog, rhos)
      !ar: arrhenius number (-)
      !rhog: gas density (kg/m3)
      !rhos:density of solid (kg/m3)

      implicit none

      real*8 :: ar, rhog, rhos

      f17 = 0.58 * ar**(-0.029) * (rhog / rhos)**0.021

end function f17

real*8 function f36 (u, umf, h, a0, g)
      !u: gas velocity
      !umf: velocity of minimum fluidization (m/s)
      !h: height
      !a0: area of the gas-distributor per nozzle, m2 per nozzle
      !g: gravity acceleration (m/s2)

      implicit none

      real*8 :: u, umf, h, a0, g

      f36 = 0.54 * ((u - umf)**0.4) * ((h + 4.0 * sqrt(a0))**0.8) * (g**(-0.2))

end function f36

real*8 function f18 (a, z)
      !a: decay factor (eq 42)
      !z: height

      implicit none

      real*8 :: z, a

      f18 = exp(-a * z)

end function f18

real*8 function f42 (dp, ug)
      !dp: particle diameter
      !ug: inlet gas velocity

      implicit none

      real*8 :: dp, ug

      f42 = 315.0 * (dp**0.64) / ug

end function f42

real*8 function f19 (cs, csb)
      !cs: solid concentration in the freeboard (eq 18)
      !csb: solid concentration in the upper limit of the dense bed

      implicit none

      real*8 :: cs, csb

      f19 = 1.0 - 0.75 * (cs / csb)**0.4

end function f19
