
SUBROUTINE write_energy_rand_angle(FRAME, N_MC, E_IN, E_FN)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	IMPLICIT NONE

	! Inputs
	INTEGER					:: FRAME	! Which frame to compute random angle in
														! 0 -> CMF		1 -> LF
	INTEGER					:: N_MC		! Number of Monte Carlo Particles to Trace

	REAL(KIND=8)		:: E_IN		! Initial Energy of all Monte Carlo Particles
	REAL(KIND=8)		:: E_FN		! Final (relaxed) Energy of all MC Particles

	! Internal
	INTEGER					:: i			! counter
	INTEGER					:: j			! counter

	REAL(KIND=8)		:: Enow		! Current energy of particle
	REAL(KIND=8)		:: Enxt		! Next energy of particle after collision
	REAL(KIND=8)		:: SAng 	! random scattering angle
	REAL(KIND=8)		:: L_SAng	! random scattering angle in lab frame

	OPEN(UNIT=555, FILE="../Data/MC_ENSEMBLE_EN_AN.dat", ACCESS="APPEND")

	DO i=1,N_MC
		j 	 = 0
		Enow = E_IN

		DO WHILE (Enow .GT. E_FN)

			! Find a random scattering angle for current energy
			CALL find_rand_angle(Enow,SAng)

			! Find new energy after colliding at SAng angle
			CALL find_new_energy(Enow,SAng,Enxt)

			!	Convert SAng to LF if FRAME = 1
			IF ( FRAME == 1) THEN
				CALL angle_to_lab(SAng, L_SAng)
				SAng = L_SAng
			END IF

			WRITE(555,*) i, Enow, SAng

			Enow = Enxt

			j = j+1

		END DO ! WHILE

		WRITE(*,*) i, " complete with ", j, " collisions"

	END DO ! i

	CLOSE(555)

END SUBROUTINE

