!equation list file
!by Murilo Varela

real*8 function f30 (bg, nc, phi)
!bg: stoichiometric coefficient to fully convert CH4 into CO2
!nc: efficiency of combustion
!phi: ratio of oxygen-carrier to fuel

      real*8 :: bg, nc, phi

      f30 = bg * nc * (1.0 /(4.0 * phi))

end function f30

real*8 function f27 (tmr, t)
!tmr: mean residence time of particle in the whole fluidized bed reactor
!t: time

      real*8 :: tmr, t

      f27 = (1.0/tmr) * exp(-t/tmr)

end function f27

real*8 function f28 (xsin, tm, t)
!xsin: mean convertion of particles in
!tm: mean reacting time
!t: time

      real*8 :: xsin, tm, t

      f28 = xsin + t/tm

end function f28

real*8 function f29 (xsout, tm, tmr, t, xsin)
!xsout: mean convertion of particles going out
!tm: mean reacting time
!tmr: mean residence time of particle in the whole fluidized bed reactor
!t: time
!xsin: mean convertion of particles going in

      real*8 :: xsout, tm, tmr, t, xsin

      !f29 = integral method
end function f29
