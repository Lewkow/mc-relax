
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculate the current gravitational acceleration
!! as a function of altitude above Mars
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mars_gravity( z, g )
	USE physics_constants, ONLY : MARS_R, MARS_M, GRAV_G

	IMPLICIT NONE

	!! Inputs	
	REAL(KIND=8)	:: z		! altitude above mars [m]

	!! Outputs
	REAL(KIND=8)	:: g		! gravitational acceleration [m/s^2]

	!! Internal
	REAL(KIND=8)	:: R		! total distance for planet center	

	R    = MARS_R + z
	g    = GRAV_G*MARS_M/(R*R)	

END SUBROUTINE mars_gravity

