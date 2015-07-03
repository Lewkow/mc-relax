
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Energy loss routine calculates the energy loss per unit time
! using the universal curve to determine diffusion cross section
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE universal_energy_loss(E, Mp, Mt, EL)
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)			:: E		! Collision energy [eV]
	REAL(KIND=8)			:: Mp		! Projectile mass
	REAL(KIND=8)			:: Mt		! Target Mass

	!! Outputs
	REAL(KIND=8)			:: EL		! Energy Loss [eV]

	!! Internal
	REAL(KIND=8)			:: DFCS, C


	!! get diffusion cross section for given energy
	CALL universal_dfcs_int( E, Mp, Mt, DFCS )	

	!! set mass constant
	C  = 2.0D0*Mp*Mt/(Mp+Mt)**2

	!! Get energy loss
	EL = C*DFCS

END SUBROUTINE universal_energy_loss

