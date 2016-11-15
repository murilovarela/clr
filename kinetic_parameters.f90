!determination of kinetic parameters
!by: Murilo Varela


subroutine kinetic_parameters (k0, e, r, t, rg, b, vm, c, n, temp, xs, ti)
      implicit none
      real*8 :: k0, e, r, t, rg, b, vm, c, n, i, temp
      real*8 :: f1_1, f1_2, f2, ti, k, xs

      k = f2 (k0, e, r, temp)
      ti = f1_2 (rg, b, vm, k, c, n)
      xs = f1_1 (t, ti)

end subroutine kinetic_parameters
