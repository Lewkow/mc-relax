
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Convert Energy [eV] to velocity [m/s]
! using the mass of the particle [amu]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE energy_to_velocity( E, m, v )
	USE physics_constants, ONLY : UTOKG, EVTOJ

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E 	! [eV]
	REAL(KIND=8)		:: m	! [amu]

	!! Outputs	
	REAL(KIND=8)		:: v	! [m/s]

	v = SQRT(2.0D0*E*EVTOJ/(m*UTOKG))	

END SUBROUTINE energy_to_velocity

