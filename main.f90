!combustion looping reactor
!by: Murilo Varela

program main

!steps 1
!     solid convertion assumption (eq 30) Xs1
!step 2
!     average solids convertion
!     RTD curve, E(t) (eq 27)
!     solid convertion distribution (eq 28,29)
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
