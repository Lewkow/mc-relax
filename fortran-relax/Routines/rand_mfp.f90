
SUBROUTINE rand_mfp( Den, TCS, MFP )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Finds a "random" mean free path based on cross sections, 
! and densities of targets for a given homogenous environment. 
! The units of MFP are the same as the density and cross section.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE rand_seed

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: Den		! Density of targets
	REAL(KIND=8)	:: TCS		! Total cross section of targets

	!! Outputs
	REAL(KIND=8)	:: MFP		! Random mean free path

	!! Internal
	REAL(KIND=8)	:: lfg		! Random Number Generator	
	REAL(KIND=8)	:: r			! Random number	
	REAL(KIND=8)	:: L			! Statistical mean free path	

	r = lfg()
	L = 1.0D0/( Den*TCS )
	
	MFP = -L*LOG( 1.0D0 - r )	

END SUBROUTINE rand_mfp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE atmosphere_rand_mfp( TCS, theta, height, R )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Finds a "random" path length based on energy dependent total
! cross sections and position dependent target particle densities. 
! Using the plane parallel path length ds = dz/cos(theta), with 
! dz = 1km, a vector P is constructed to hold the probabilities 
! of traveling a set distance without colliding. This is then
! compared to a "randomly" generated number to see what distance
! the particle actually travels.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE rand_seed
	USE physics_constants, ONLY : PI
	
	IMPLICIT NONE

	!! inputs
	REAL(KIND=8)		:: TCS, theta, height

	!! outputs
	REAL(KIND=8)		:: R

	!! internal
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: P
	INTEGER																:: N_z, i, j
	REAL(KIND=8)													:: z, dz, den, zi, zf, Ab_zi, Ab_zf
	REAL(KIND=8)													:: lfg
	REAL(KIND=8)													:: tot, lamda

	!! set differential for height [m]
	dz = 1.0D3	

	!! if particle traveling down
	IF (theta .LT. PI/2.0D0) THEN
	
		!! Set number of z intervals to height above ground
		N_z = CEILING(height/dz)

		!! Allocate P array
		ALLOCATE(P(N_z))
	
		!! Set inital and final absolute height bounds
		Ab_zi = height
		Ab_zf = 0.0D0	
		
		!! Start loop to find P array	
		DO i=1,N_z
			tot = 0.0D0
			z   = Ab_zi - i*dz
			DO j=0,i
				CALL mars_density(16, z, den)
				lamda = 1.0D0/(den*TCS)
				tot = tot + dz/(COS(theta)*lamda)
			END DO
			P(i) = EXP(-tot)
			WRITE(7,*) P(i)	
		END DO

	END IF

	R = 1.0D0

	!! clean up memory
	DEALLOCATE(P)


END SUBROUTINE atmosphere_rand_mfp
