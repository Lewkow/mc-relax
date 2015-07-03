
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Calculates the average angle for a given energy 
! using several random numbers and writes results to 
! a file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
SUBROUTINE average_scattering_angle
	USE planet
	USE universal

	IMPLICIT NONE

	!! Internal
	CHARACTER(LEN=4)		:: Coll_Type
	REAL(KIND=8)				:: d_En, En, start_En, end_En, ang_tot, ang
	INTEGER							:: N_En, N_MC, i, j, NICK

	!! open files to write to
	OPEN(UNIT=78,FILE='../Data/HeO_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=79,FILE='../Data/HeH_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=80,FILE='../Data/HeHe_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=81,FILE='../Data/HH_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=82,FILE='../Data/HO_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=83,FILE='../Data/HAr_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=84,FILE='../Data/HH2_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=85,FILE='../Data/HN2_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=86,FILE='../Data/HCO_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=87,FILE='../Data/HCO2_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=88,FILE='../Data/HeAr_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=89,FILE='../Data/HeH2_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=90,FILE='../Data/HeN2_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=91,FILE='../Data/HeCO_Average_Scattering_Angles.dat',ACCESS='APPEND')
	OPEN(UNIT=92,FILE='../Data/HeCO2_Average_Scattering_Angles.dat',ACCESS='APPEND')

	!! set energy start, end, number of points
	start_En = 0.01D0
	end_En   = 2000.0D0
	N_En     = 100
	d_En     = (end_En-start_En)/REAL(N_En-1)
		
	!! set number of random angles to collect
	N_MC     = 10000

	NICK = 1
	IF (NICK .EQ. 1) THEN
	N_MC = 10000000
	En   = 5000.0D0
	DO i=1,N_MC
		CALL lin_rand_angle(En,'HeH  ', ang)
		WRITE(640,*) ang
	END DO
	WRITE(*,*) 'NICK, Im done with your test!'
	STOP
	END IF

	!! write info to screen
	WRITE(*,'(A)') 'STARTING AVERAGE SCATTERING CALCULATION'
	WRITE(*,'(A)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	WRITE(*,'(A,I7,A,ES10.2,A,ES10.2,A)') &
	& 'Using ', N_En, ' energies from ', start_En, ' eV to ', end_En, ' eV'
	WRITE(*,'(A,I7,A)') 'Using ', N_MC, ' random scattering events per energy'
	WRITE(*,'(A)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	WRITE(*,*)

	!! do He+O
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HeO  ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(78,*) En, ang
	END DO	

	WRITE(*,*) 'He+O Complete'

	!! do He+H
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HeH  ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(79,*) En, ang
	END DO	

	WRITE(*,*) 'He+H Complete'

	!! do He+He
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HeHe ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(80,*) En, ang
	END DO	

	WRITE(*,*) 'He+He Complete'

	!! do H+H
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HH   ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(81,*) En, ang
	END DO	

	WRITE(*,*) 'H+H Complete'

	PROJ = 'H '
	CALL read_uni_tables
	WRITE(*,*) 'NA H+O: ', NUM_ANG

	!! H+O
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HO   ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(82,*) En, ang
	END DO	

	WRITE(*,*) 'H+O Complete'

	!! H+Ar
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HAr  ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(83,*) En, ang
	END DO	

	WRITE(*,*) 'H+Ar Complete'

	!! H+H2
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HH2  ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(84,*) En, ang
	END DO	

	WRITE(*,*) 'H+H2 Complete'

	!! H+N2
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HN2  ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(85,*) En, ang
	END DO	

	WRITE(*,*) 'H+N2 Complete'

	!! H+CO
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HCO  ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(86,*) En, ang
	END DO	

	WRITE(*,*) 'H+CO Complete'

	!! H+CO2
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HCO2 ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(87,*) En, ang
	END DO	

	WRITE(*,*) 'H+CO2 Complete'

	CALL clean_uni_tables
	PROJ = 'He'
	CALL read_uni_tables

	!! He+Ar
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HeAr ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(88,*) En, ang
	END DO	

	WRITE(*,*) 'He+Ar Complete'

	!! He+H2
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HeH2 ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(89,*) En, ang
	END DO	

	WRITE(*,*) 'He+H2 Complete'

	!! He+N2
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HeN2 ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(90,*) En, ang
	END DO	

	WRITE(*,*) 'He+N2 Complete'

	!! He+CO
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HeCO ', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(91,*) En, ang
	END DO	

	WRITE(*,*) 'He+CO Complete'

	!! He+CO2
	DO j=1,N_En
		En      = start_En + (j-1)*d_En
		ang_tot = 0.0D0
		DO i=1,N_MC		
			CALL lin_rand_angle(En,'HeCO2', ang)
			ang_tot = ang_tot + ang
		END DO
		ang = ang_tot/REAL(N_MC)
		WRITE(92,*) En, ang
	END DO	

	WRITE(*,*) 'He+CO2 Complete'

	!! close files
	CLOSE(78)
	CLOSE(79)
	CLOSE(80)
	CLOSE(81)
	CLOSE(82)
	CLOSE(83)
	CLOSE(84)
	CLOSE(85)
	CLOSE(86)
	CLOSE(87)
	CLOSE(88)
	CLOSE(89)
	CLOSE(90)
	CLOSE(91)
	CLOSE(92)

END SUBROUTINE average_scattering_angle

