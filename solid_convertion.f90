!solid convertion subroutines
!by: Murilo Varela

subroutine convertion_assumption (bg, nc, phi, xs1)
      implicit none

      real*8 :: bg, nc, phi
      real*8 :: f30, xs1

      xs1 = f30 (bg, nc, phi)

end subroutine convertion_assumption

subroutine rtd_curve (tmr, t, rtd)
      implicit none

      real*8 :: f27, tmr, t, rtd

      rtd = f27 (tmr, t)

end subroutine rtd_curve

subroutine mean_react_time (xsout, tm, tmr, xsin)
      implicit none

      real*8 :: xsout, tmr, xsin
      real*8 :: f29_i, e, es, tm

      integer :: n, i

      n = 0
      tm = e
      e = 1.0e-6
      es = 1.0e-8
      do while (abs(f29_i (xsout, tm, tmr, xsin) - (1.0 - xsout)) > e)
            tm = tm + es
            n = n + 1
      end do
      print*, tm, n
end subroutine mean_react_time
