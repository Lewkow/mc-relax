
SUBROUTINE rand_init_velocity( theta, phi, u )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine finds random initial velocity vectors. 
! Theta returns the polar angle which is randomly generated, 
! Phi returns the azimuthal angle which is randomly generated, 
! u returns a 3 component vector which contains the cartesian
! (x,y,z) representation of the random velocity vector.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!! Output
	REAL(KIND=8)							:: theta, phi
	REAL(KIND=8),DIMENSION(3)	:: u

	!! Internal
	REAL(KIND=8)							:: lfg					! Random number generator
	REAL(KIND=8)							:: r1, r2				! Random numbers
	
	r1 = lfg()
	r2 = lfg()
	
	theta = r1*PI
	phi   = r2*2.0D0*PI
	
	u(1)  = SIN(theta)*COS(phi)
	u(2)  = SIN(theta)*SIN(phi)
	u(3)  = COS(theta)	

END SUBROUTINE rand_init_velocity

