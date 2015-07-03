!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine will determine if, after a collision, the 
! target particle is considered hot. If so, it also finds the 
! velocity vector of the secondary hot particle in the planet
! cartesian system. These values are returned so that they can
! be written to an output data file. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE secondary_hot( E1, E2, E_hot, Theta, Phi )
	USE planet, ONLY : E_Therm, MP, MT
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!!! Inputs
	REAL(KIND=8)			:: E1, E2

	!!! Outputs
	REAL(KIND=8)			:: E_hot, Theta, Phi

	!!! Internal
	REAL(KIND=8)			:: C1, C2, M
	REAL(KIND=8)			:: lfg

	!!! Calculate energy of target particle after collision
	E_hot = E1 - E2

	!!! Set secondary particle energy to zero if below thermal temp
	IF ( E_hot .LT. E_Therm ) THEN

		E_hot = 0.0D0	
		Theta = 0.0D0
		Phi   = 0.0D0

	ELSE

		!!! Find scattering angle of target if hot	
		M  = (MP + MT)**2.0D0 / (4.0D0*MP*MT)	
		C1 = E_hot/E1
		C2 = M*C1
		
		Theta = ACOS(SQRT(C2)) 
		Phi   = lfg()*2.0D0*PI

	END IF

END SUBROUTINE secondary_hot



