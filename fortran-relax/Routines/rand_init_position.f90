
SUBROUTINE rand_init_iso_position( x )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine finds random isotropic initial position vector. 
! x returns a 3 component vector which contains the cartesian
! (x,y,z) representation of the random "isotropic" position vector.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!! Output
	REAL(KIND=8),DIMENSION(3)	:: x

	!! Internal
	REAL(KIND=8)							:: lfg					! Random number generator
	REAL(KIND=8)							:: r1, r2				! Random numbers
	REAL(KIND=8)							:: theta, phi		! Polar angles
	
	r1 = lfg()
	r2 = lfg()
	
	theta = r1*PI
	phi   = r2*2.0D0*PI
	
	x(1)  = SIN(theta)*COS(phi)
	x(2)  = SIN(theta)*SIN(phi)
	x(3)  = COS(theta)	

END SUBROUTINE rand_init_iso_position

