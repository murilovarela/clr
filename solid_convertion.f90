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
            real*8 :: sol, sol_

            integer :: n

            n = 0
            e = 1.0e-14
            es = 1.0e-1
            tm = e            !14.733253138492049       1 473 325 263

            do while (abs(f29_i (xsout, tm, tmr, xsin) - (1.0 - xsout)) > e)

                  sol = f29_i (xsout, tm, tmr, xsin) - (1.0 - xsout)
                  sol_ = f29_i (xsout, tm + es, tmr, xsin) - (1.0 - xsout)

                  if (sol/abs(sol) == sol_/abs(sol_)) then
                        sol = f29_i (xsout, tm, tmr, xsin) - (1.0 - xsout)
                        sol_ = f29_i (xsout, tm + es, tmr, xsin) - (1.0 - xsout)

                        tm = tm + es
                        n = n + 1
                  else
                        es = es/10.0
                  end if
            end do
end subroutine mean_react_time
