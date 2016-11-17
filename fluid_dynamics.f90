!Fluid dynamics of the reactor
!by: Murilo Varela

subroutine velocity_min_fluidization (c1, c2, ar, rhog, mdp, mug, umf)

      implicit none

      real*8 :: re, f9_1, c1, c2, ar
      real*8 :: umf, f9_2, rhog, mdp, mug

      re = f9_1 (c1, c2, ar)
      umf = f9_2 (rhog, mdp, mug, re)

end subroutine velocity_min_fluidization

subroutine single_bubble_velocity (u, umf, h, a0, g, ubinf)

      implicit none

      real*8 :: u, umf, h, a0, g, db, ubinf, f36, f15

      db = f36 (u, umf, h, a0, g)
      ubinf = f15 (g, db)

end subroutine single_bubble_velocity

subroutine visible_velocity_bubble_i (psi, u0, umf, sigb, uvis, ubinf)

      implicit none

      real*8 :: psi, u0, umf, sigb, uvis, e, sigb_, ubinf
      real*8 :: f10, f14
      integer :: n

      n = 0
      e = 1.0e-14
      sigb_ = 1.0

      do while (abs(sigb - sigb_) > e)

            sigb = sigb_
            n = n + 1

            uvis = f10 (psi, u0, umf, sigb)
            sigb_ = f14 (uvis, ubinf)

      end do
      sigb = sigb_

end subroutine visible_velocity_bubble_i

subroutine visible_velocity_bubble (u, g, dp, umf, h, a0, u0, sigb, ubinf)

      implicit none

      real*8 :: dp, umf, h, a0, ubinf, g, u
      real*8 :: u0, u0_, sigb, uvis, fb, psi, utf, e
      real*8 :: f8, f11, f12, f13
      integer :: n

      n = 0
      e = 1.0e-14
      u0_ = u

      do while (abs(u0 - u0_) > e)

            u0 = u0_
            n = n + 1

            call single_bubble_velocity (u0, umf, h, a0, g, ubinf)

            fb = f13 (dp, u0, umf)
            psi = f12 (fb, h, a0)

            call visible_velocity_bubble_i (psi, u0, umf, sigb, uvis, ubinf)

            utf = f11 (psi, u0, umf, sigb)
            !print*, fb, psi, utf, u0, u0_
            u0_ = f8 (umf, uvis, utf, sigb)

      end do

      u0 = u0_

end subroutine visible_velocity_bubble
