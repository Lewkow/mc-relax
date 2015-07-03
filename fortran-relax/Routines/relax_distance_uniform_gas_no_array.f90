
SUBROUTINE uniform_relax_distance_no_array(N_MC, E_IN, E_FN)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine takes a number N_MC particles and
! an initial and final energy for each particle 
! and writes the current energy as well as the 
! random scattering angle for that energy to a 
! file in the following format
!
!	Particle_Number	| Energy(eV) | Random_Scatt_Angle(rad)
!
! NO_ARRAY version does not save all histories to arrays
!	to save memory for large MC calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI	
	USE mpi_info

	IMPLICIT NONE

	! Inputs
	INTEGER					:: FRAME				! Which frame to compute random angle in
																	! 0 -> CMF		1 -> LF
	INTEGER					:: N_MC					! Number of Monte Carlo Particles to Trace

	REAL(KIND=8)		:: E_IN					! Initial Energy of all Monte Carlo Particles
	REAL(KIND=8)		:: E_FN					! Final (relaxed) Energy of all MC Particles

	! Internal
	INTEGER					:: MC						! counter
	INTEGER					:: COL					! counter
	INTEGER					:: i						! counter
	INTEGER					:: j						! counter

	INTEGER					:: Max_Coll			! Maximum number of collisions for a single MC particle

	INTEGER					:: My_N_MC			! Number of MC particle for each processor

	INTEGER					:: HARD_SPHERE	! Set to 1 if HS cross section is used instead of quantum

	REAL(KIND=8)		:: HS_R					! Radius to use for HS approximation

	REAL(KIND=8)		:: Enow					! Current energy of particle
	REAL(KIND=8)		:: Enxt					! Next energy of particle after collision
	REAL(KIND=8)		:: SAng 				! random scattering angle
	REAL(KIND=8)		:: L_SAng				! random scattering angle in lab frame

	REAL(KIND=8)		:: Den					! Density of uniform gas in [m^-3]

	REAL(KIND=8)		:: L						! Mean Free Path
	REAL(KIND=8)		:: LL						! Random Free Path 

	REAL(KIND=8)		:: TCS					! Total cross section 

	REAL(KIND=8)		:: lfg					! random number function
	REAL(KIND=8)		:: r1						! random number
	REAL(KIND=8)		:: r2						! random number

	REAL(KIND=8),DIMENSION(3)										:: R_Prev		! Previous Position vector
	REAL(KIND=8),DIMENSION(3)										:: R_Curr		! Current Position vection

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)				:: th				! theta scattering angle array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)				:: ph				! phi scattering angle array

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)			:: u				! Scattering angle array

	REAL(KIND=8)																:: DIS 			! Final Displacment 
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)				:: NOW_DIS 	! Individual Collision Displacment array 
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)				:: NOW_EN	 	! Individual Collision Energy array 

	IF (myid == 0) THEN
		OPEN(UNIT=555, FILE="../Data/MC_ENSEMBLE_UNI_RELAX_FINAL_DISP.dat", ACCESS="APPEND")
		OPEN(UNIT=666, FILE="../Data/MC_ENSEMBLE_UNI_RELAX_ENGY_DISP.dat", ACCESS="APPEND")
	END IF

	! Set HARD_SPHERE to 1 for HS approximation for cross sections
	HARD_SPHERE = 1
	
	! Radius in bohr radii
	HS_R = 0.68D0 
	
	! Set density of bath gas [m^-3]
	! Density equal to 0.1 cm^-3 from Local bubble paper, cloudlet parameters
	Den = 1.0D5

	! Set Maximum number of collisions per MC Particle
	Max_Coll = 1000

	! Set Number of MC particles per processor
	My_N_MC = N_MC/np

	! Allocate arrays
	ALLOCATE( th(Max_Coll), ph(Max_Coll) )
	ALLOCATE( u(Max_Coll, 3) )
	ALLOCATE( NOW_DIS(Max_Coll), NOW_EN(Max_Coll) )

	DO MC=1,My_N_MC
		COL	 = 1 
		Enow = E_IN

		! Set up initial random, isotropic trajectory
		r1    = lfg()	
		r2    = lfg()	
		th(1) = r1*PI 
		ph(1) = r2*2.0D0*PI

		! Apply inital angles to u array
		u(1,1) = SIN(th(1))*COS(ph(1))
		u(1,2) = SIN(th(1))*SIN(ph(1))
		u(1,3) = COS(th(1))

		! Set initial location at origin with energy E_IN
		R_Prev(1) = 0.0D0
		R_Prev(2) = 0.0D0
		R_Prev(3) = 0.0D0

		DO WHILE (Enow .GT. E_FN)

			! Find TCS for current Energy
			IF ( HARD_SPHERE .EQ. 1 ) THEN
				TCS = PI*HS_R*HS_R
			ELSE
				CALL TCS_HeH(Enow,TCS)
			END IF

			! Convert TCS to [m]
			TCS = TCS*(5.29177249D-11)**2.0D0

			! Find MFP for current Energy (in units of [AU])
			L = 1.0D0/(Den*TCS*1.5E11)	

			! Find a random theta scattering angle for current energy
			IF ( HARD_SPHERE .EQ. 1 ) THEN
				r1   = lfg()
				SAng = ACOS(1.0D0 - 2.0D0*r1)
			ELSE
				CALL find_rand_angle(Enow,SAng)
			END IF

			! Find random distance traveled before collision
			r1 = lfg()
			LL = -L*LOG(1.0D0-r1)

			! Find new energy after colliding at SAng angle
			CALL find_new_energy(Enow,SAng,Enxt)

			!	Convert SAng to LF
			CALL angle_to_lab(SAng, L_SAng)
			th(COL+1) = L_SAng

			! Find a random phi scattering angle
			r1 				= lfg()
			ph(COL+1) = r1*2.0D0*PI	

			! Transport particle
			R_Curr(1) = R_Prev(1) + u(COL,1)*LL	
			R_Curr(2) = R_Prev(2) + u(COL,2)*LL	
			R_Curr(3) = R_Prev(3) + u(COL,3)*LL	

			! Rotate axis to update angular u array
    	u(COL+1,1) = u(COL,1)*COS(th(COL+1)) + SIN(th(COL+1))*(u(COL,3)*COS(ph(COL+1))*COS(th(COL)) - SIN(ph(COL+1))*SIN(ph(COL)))
    	u(COL+1,2) = u(COL,2)*COS(th(COL+1)) + SIN(th(COL+1))*(u(COL,3)*COS(ph(COL+1))*SIN(th(COL)) + SIN(ph(COL+1))*COS(ph(COL)))
    	u(COL+1,3) = u(COL,3)*COS(th(COL+1)) - SIN(th(COL+1))*SIN(th(COL))*COS(ph(COL+1))

			Enow = Enxt
		
			NOW_DIS(COL) = SQRT( R_Curr(1)**2.0D0 + R_Curr(2)**2.0D0 + R_Curr(3)**2.0D0 )
			NOW_EN(COL)  = Enow	

			R_Prev(1) = R_Curr(1)
			R_Prev(2) = R_Curr(2)
			R_Prev(3) = R_Curr(3)
	
			COL = COL+1

		END DO ! WHILE

		! Write NOW_DIS and NOW_EN data to file
!		DO i=1,(COL-1)
!			WRITE(666,*) NOW_DIS(i), NOW_EN(i)
!		END DO

		! Calculate total displacment of MC particle
		DIS = SQRT( R_Curr(1)**2.0D0 + R_Curr(2)**2.0D0 + R_Curr(3)**2.0D0 )

		IF (MOD(MC,1000) .EQ. 0) THEN
			WRITE(*,*) MC, " complete with ", COL, " collisions", " Displaced [AU] ", DIS
		END IF

		WRITE(555,*) DIS

	END DO ! MC 

	WRITE(*,*) "Processor ", myid, " Complete"

	CLOSE(555)
	CLOSE(666)

	! free arrays
	DEALLOCATE(th,ph,u,NOW_DIS,NOW_EN)


END SUBROUTINE


