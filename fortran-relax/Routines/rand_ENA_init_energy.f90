
SUBROUTINE rand_ENA_init_energy( Energy )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Find intial random ENA energy based on paper
! by Zank et al. (2009)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!!!! Outputs
	REAL(KIND=8)		:: Energy

	!!!! Internal
	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: En_low, En_high
	REAL(KIND=8)		:: n, r

	En_low  = 0.1D0
	En_high = 6.0D0

	n = -0.63D0
	r = lfg()
	
	Energy = 1000.0D0*( En_low**n + r*( En_high**n - En_low**n) )**(1.0D0/n)	

END SUBROUTINE

