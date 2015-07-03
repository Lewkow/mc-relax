
SUBROUTINE ang_prob_data_read
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in Angular Probability data table from file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE tables
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: table_HeO, table_HeH, table_HeHe, table_HH
	REAL,PARAMETER													:: MB_per_real = 8.0D-6
	REAL,PARAMETER													:: GB_per_real = 8.0D-9

	REAL(KIND=8)														:: E0
	REAL(KIND=8)														:: Enow

	INTEGER																	:: NUM_LINES_HeO, NUM_LINES_HeH, NUM_LINES_HeHe
	INTEGER																	:: NUM_LINES_HH
	INTEGER																	:: N
	INTEGER																	:: i 
	INTEGER																	:: j 
	INTEGER																	:: k, s1, s2, s3, s4

	IF ( myid .EQ. 0 ) THEN

		!! open tables files to read from
		OPEN(UNIT=10, FILE="../Tables/AngProb3D_HeO.dat", STATUS="old", ACTION="read")
		OPEN(UNIT=11, FILE="../Tables/AngProb3D_HeH.dat", STATUS="old", ACTION="read")
		OPEN(UNIT=12, FILE="../Tables/AngProb3D_HeHe_NonSym.dat", STATUS="old", ACTION="read")
		OPEN(UNIT=13, FILE="../Tables/AngProb3D_HH_NonSym.dat", STATUS="old", ACTION="read")

		!! read in how many lines are in each file
		READ(10,*) NUM_LINES_HeO
		READ(11,*) NUM_LINES_HeH
		READ(12,*) NUM_LINES_HeHe
		READ(13,*) NUM_LINES_HH

		!! allocate space in tables
		ALLOCATE(table_HeO (NUM_LINES_HeO ,3))
		ALLOCATE(table_HeH (NUM_LINES_HeH ,3))
		ALLOCATE(table_HeHe(NUM_LINES_HeHe,3))
		ALLOCATE(table_HH  (NUM_LINES_HH,3) )

		!! fill in main tables for each collision type
		DO i=1,NUM_LINES_HeO
			READ(10,*) table_HeO(i,1), table_HeO(i,2), table_HeO(i,3)
		END DO

		DO i=1,NUM_LINES_HeH
			READ(11,*) table_HeH(i,1), table_HeH(i,2), table_HeH(i,3)
		END DO

		DO i=1,NUM_LINES_HeHe
			READ(12,*) table_HeHe(i,1), table_HeHe(i,2), table_HeHe(i,3)
		END DO

		DO i=1,NUM_LINES_HH
			READ(13,*) table_HH(i,1), table_HH(i,2), table_HH(i,3)
		END DO

		!! He+O
		E0 		= table_HeO(1,1)
		Enow  = E0	
		N 		= 0 
	
		DO WHILE (Enow == E0)
			N    = N + 1
			Enow = table_HeO(N,1)		
		END DO

		NUM_AN_HeO = N-1
		NUM_EN_HeO = NUM_LINES_HeO/NUM_AN_HeO

		!! He+H
		E0 		= table_HeH(1,1)
		Enow  = E0	
		N 		= 0 
	
		DO WHILE (Enow == E0)
			N    = N + 1
			Enow = table_HeH(N,1)		
		END DO

		NUM_AN_HeH = N-1
		NUM_EN_HeH = NUM_LINES_HeH/NUM_AN_HeH

		!! He+He
		E0 		= table_HeHe(1,1)
		Enow  = E0	
		N 		= 0 
	
		DO WHILE (Enow == E0)
			N    = N + 1
			Enow = table_HeHe(N,1)		
		END DO

		NUM_AN_HeHe = N-1
		NUM_EN_HeHe = NUM_LINES_HeHe/NUM_AN_HeHe

		!! H+H
		E0   = table_HH(1,1)
		Enow = E0
		N    = 0

		DO WHILE (Enow .EQ. E0)
			N    = N+1
			Enow = table_HH(N,1)
		END DO

		NUM_AN_HH = N-1
		NUM_EN_HH = NUM_LINES_HH/NUM_AN_HH

		!! Alocate memoery for energy, angle, PD and spline derivatives for all collision tables
		ALLOCATE(F_ENERGY_HeO(NUM_EN_HeO),F_ANGLE_HeO(NUM_AN_HeO),F_PD_HeO(NUM_EN_HeO,NUM_AN_HeO),SP_Y2_HeO(NUM_AN_HeO))
		ALLOCATE(F_ENERGY_HeH(NUM_EN_HeH),F_ANGLE_HeH(NUM_AN_HeH),F_PD_HeH(NUM_EN_HeH,NUM_AN_HeH),SP_Y2_HeH(NUM_AN_HeH))
		ALLOCATE(F_ENERGY_HeHe(NUM_EN_HeHe),F_ANGLE_HeHe(NUM_AN_HeHe),F_PD_HeHe(NUM_EN_HeHe,NUM_AN_HeHe),SP_Y2_HeHe(NUM_AN_HeHe))
		ALLOCATE(F_ENERGY_HH(NUM_EN_HH),F_ANGLE_HH(NUM_AN_HH),F_PD_HH(NUM_EN_HH,NUM_AN_HH),SP_Y2_HH(NUM_AN_HH))

		s1 = SIZE(F_ENERGY_HeO)  + SIZE(F_ANGLE_HeO)  + SIZE(F_PD_HeO)  + SIZE(SP_Y2_HeO)
		s2 = SIZE(F_ENERGY_HeH)  + SIZE(F_ANGLE_HeH)  + SIZE(F_PD_HeH)  + SIZE(SP_Y2_HeH)
		s3 = SIZE(F_ENERGY_HeHe) + SIZE(F_ANGLE_HeHe) + SIZE(F_PD_HeHe) + SIZE(SP_Y2_HeHe)
		s4 = SIZE(F_ENERGY_HH)   + SIZE(F_ANGLE_HH)   + SIZE(F_PD_HH)   + SIZE(SP_Y2_HH)

		!! File Angle Arrays from Tables
		DO i=1,NUM_AN_HeO
			F_ANGLE_HeO(i)  = table_HeO(i,2)
		END DO

		DO i=1,NUM_AN_HeH
			F_ANGLE_HeH(i)  = table_HeH(i,2)
		END DO

		DO i=1,NUM_AN_HeHe
			F_ANGLE_HeHe(i) = table_HeHe(i,2)
		END DO

		DO i=1,NUM_AN_HH
			F_ANGLE_HH(i)   = table_HH(i,2)
		END DO

		!! File Probability Density Matrix and File Energy Array from Table
		DO i=1,NUM_EN_HeO
			DO j=1,NUM_AN_HeO
				k = (i-1)*NUM_AN_HeO + j
				F_PD_HeO(i,j) = table_HeO(k,3)
			END DO
			F_ENERGY_HeO(i) = table_HeO(k,1)
		END DO

		DO i=1,NUM_EN_HeH
			DO j=1,NUM_AN_HeH
				k = (i-1)*NUM_AN_HeH + j
				F_PD_HeH(i,j) = table_HeH(k,3)
			END DO
			F_ENERGY_HeH(i) = table_HeH(k,1)
		END DO

		DO i=1,NUM_EN_HeHe
			DO j=1,NUM_AN_HeHe
				k = (i-1)*NUM_AN_HeHe + j
				F_PD_HeHe(i,j) = table_HeHe(k,3)
			END DO
			F_ENERGY_HeHe(i) = table_HeHe(k,1)
		END DO

		DO i=1,NUM_EN_HH
			DO j=1,NUM_AN_HH
				k = (i-1)*NUM_AN_HH + j
				F_PD_HH(i,j) = table_HH(k,3)
			END DO
			F_ENERGY_HH(i) = table_HH(k,1)
		END DO

		!! Free memory from tables
		DEALLOCATE(table_HeO, table_HeH, table_HeHe, table_HH)
		CLOSE(10)
		CLOSE(11)
		CLOSE(12)
		CLOSE(13)

	END IF ! root

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_EN_HeO,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_AN_HeO,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_EN_HeH,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_AN_HeH,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_EN_HeHe, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_AN_HeHe, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_EN_HH,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_AN_HH,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

	IF ( myid .NE. 0 ) THEN
    ALLOCATE(F_ENERGY_HeO(NUM_EN_HeO),F_ANGLE_HeO(NUM_AN_HeO),F_PD_HeO(NUM_EN_HeO,NUM_AN_HeO),SP_Y2_HeO(NUM_AN_HeO))
    ALLOCATE(F_ENERGY_HeH(NUM_EN_HeH),F_ANGLE_HeH(NUM_AN_HeH),F_PD_HeH(NUM_EN_HeH,NUM_AN_HeH),SP_Y2_HeH(NUM_AN_HeH))
    ALLOCATE(F_ENERGY_HeHe(NUM_EN_HeHe),F_ANGLE_HeHe(NUM_AN_HeHe),F_PD_HeHe(NUM_EN_HeHe,NUM_AN_HeHe),SP_Y2_HeHe(NUM_AN_HeHe))
    ALLOCATE(F_ENERGY_HH(NUM_EN_HH),F_ANGLE_HH(NUM_AN_HH),F_PD_HH(NUM_EN_HH,NUM_AN_HH),SP_Y2_HH(NUM_AN_HH))
	END IF ! not root

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_ENERGY_HeO, NUM_EN_HeO, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_ANGLE_HeO,  NUM_AN_HeO, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_PD_HeO,NUM_EN_HeO*NUM_AN_HeO,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(F_ENERGY_HeH, NUM_EN_HeH, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_ANGLE_HeH,  NUM_AN_HeH, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_PD_HeH,NUM_EN_HeH*NUM_AN_HeH,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(F_ENERGY_HeHe, NUM_EN_HeHe, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_ANGLE_HeHe,  NUM_AN_HeHe, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_PD_HeHe,NUM_EN_HeHe*NUM_AN_HeHe,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(F_ENERGY_HH, NUM_EN_HH, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_ANGLE_HH,  NUM_AN_HH, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST(F_PD_HH,NUM_EN_HH*NUM_AN_HH,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

	IF (myid == 0) THEN
		34 FORMAT(A, ES10.2)
		35 FORMAT(A, I6, A, I5) 
		WRITE(*,*)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A)')      'He + O TABLE'
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,35)         "N_ENERGIES: ", NUM_EN_HeO, " N_ANGLES: ", NUM_AN_HeO
		WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY_HeO(1)
		WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY_HeO(NUM_EN_HeO)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A)')      'He + H TABLE'
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,35)         "N_ENERGIES: ", NUM_EN_HeH, " N_ANGLES: ", NUM_AN_HeH
		WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY_HeH(1)
		WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY_HeH(NUM_EN_HeH)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A)')      'H + H TABLE'
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,35)         "N_ENERGIES: ", NUM_EN_HH, " N_ANGLES: ", NUM_AN_HH
		WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY_HH(1)
		WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY_HH(NUM_EN_HH)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A)')      'He + He TABLE'
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,35)         "N_ENERGIES: ", NUM_EN_HeHe, " N_ANGLES: ", NUM_AN_HeHe
		WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY_HeHe(1)
		WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY_HeHe(NUM_EN_HeHe)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A,F6.2)') 'Memory Footprint HeO   [MB]: ', np*s1*MB_per_real
		WRITE(*,'(A,F6.2)') 'Memory Footprint HeH   [MB]: ', np*s2*MB_per_real
		WRITE(*,'(A,F6.2)') 'Memory Footprint HH    [MB]: ', np*s4*MB_per_real
		WRITE(*,'(A,F6.2)') 'Memory Footprint HeHe  [MB]: ', np*s3*MB_per_real
		WRITE(*,'(A,F6.2)') 'Total Memory Footprint [MB]: ', np*(s1+s2+s3+s4)*MB_per_real	
		WRITE(*,'(A)')      '------------------------------------------------------'
		WRITE(*,*)
	END IF

END SUBROUTINE ang_prob_data_read

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE pd_data_read
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in DCS data table from file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE tables
	USE mpi_info

	IMPLICIT NONE

	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: table_HeO, table_HeH, table_HeHe
	REAL,PARAMETER													:: MB_per_real = 8.0D-6
	REAL,PARAMETER													:: GB_per_real = 8.0D-9

	REAL(KIND=8)														:: E0
	REAL(KIND=8)														:: Enow

	INTEGER																	:: NUM_LINES_HeO, NUM_LINES_HeH, NUM_LINES_HeHe
	INTEGER																	:: N
	INTEGER																	:: i 
	INTEGER																	:: j 
	INTEGER																	:: k, s1, s2, s3

	OPEN(UNIT=10, FILE="../Tables/Probability_Density_HeO.dat", STATUS="old", ACTION="read")
	OPEN(UNIT=11, FILE="../Tables/Probability_Density_HeH.dat", STATUS="old", ACTION="read")
	OPEN(UNIT=12, FILE="../Tables/Probability_Density_Sym_HeHe.dat", STATUS="old", ACTION="read")

	READ(10,*) NUM_LINES_HeO
	READ(11,*) NUM_LINES_HeH
	READ(12,*) NUM_LINES_HeHe

	ALLOCATE(table_HeO (NUM_LINES_HeO,3))
	ALLOCATE(table_HeH (NUM_LINES_HeH,3))
	ALLOCATE(table_HeHe(NUM_LINES_HeHe,3))

	DO i=1,NUM_LINES_HeO
		READ(10,*) table_HeO(i,1), table_HeO(i,2), table_HeO(i,3)
	END DO

	DO i=1,NUM_LINES_HeH
		READ(11,*) table_HeH(i,1), table_HeH(i,2), table_HeH(i,3)
	END DO

	DO i=1,NUM_LINES_HeHe
		READ(12,*) table_HeHe(i,1), table_HeHe(i,2), table_HeHe(i,3)
	END DO

	E0 		= table_HeO(1,1)
	Enow  = E0	
	N 		= 0 
	
	DO WHILE (Enow == E0)
		N    = N + 1
		Enow = table_HeO(N,1)		
	END DO

	NUM_AN_HeO = N-1
	NUM_EN_HeO = NUM_LINES_HeO/NUM_AN_HeO

	E0 		= table_HeH(1,1)
	Enow  = E0	
	N 		= 0 
	
	DO WHILE (Enow == E0)
		N    = N + 1
		Enow = table_HeH(N,1)		
	END DO

	NUM_AN_HeH = N-1
	NUM_EN_HeH = NUM_LINES_HeH/NUM_AN_HeH

	E0 		= table_HeHe(1,1)
	Enow  = E0	
	N 		= 0 
	
	DO WHILE (Enow == E0)
		N    = N + 1
		Enow = table_HeHe(N,1)		
	END DO

	NUM_AN_HeHe = N-1
	NUM_EN_HeHe = NUM_LINES_HeHe/NUM_AN_HeHe

	!! Alocate memoery for energy, angle, PD and spline derivatives for all collision tables
	ALLOCATE(F_ENERGY_HeO(NUM_EN_HeO),  F_ANGLE_HeO(NUM_AN_HeO),  F_PD_HeO(NUM_EN_HeO,NUM_AN_HeO),  SP_Y2_HeO(NUM_AN_HeO))
	ALLOCATE(F_ENERGY_HeH(NUM_EN_HeH),  F_ANGLE_HeH(NUM_AN_HeH),  F_PD_HeH(NUM_EN_HeH,NUM_AN_HeH),  SP_Y2_HeH(NUM_AN_HeH))
	ALLOCATE(F_ENERGY_HeHe(NUM_EN_HeHe), F_ANGLE_HeHe(NUM_AN_HeHe), F_PD_HeHe(NUM_EN_HeHe,NUM_AN_HeHe), SP_Y2_HeHe(NUM_AN_HeHe))

	s1 = SIZE(F_ENERGY_HeO) + SIZE(F_ANGLE_HeO) + SIZE(F_PD_HeO) + SIZE(SP_Y2_HeO)
	s2 = SIZE(F_ENERGY_HeH) + SIZE(F_ANGLE_HeH) + SIZE(F_PD_HeH) + SIZE(SP_Y2_HeH)
	s3 = SIZE(F_ENERGY_HeHe) + SIZE(F_ANGLE_HeHe) + SIZE(F_PD_HeHe) + SIZE(SP_Y2_HeHe)

	!! File File Angle Arrays from Tables
	DO i=1,NUM_AN_HeO
		F_ANGLE_HeO(i) = table_HeO(i,2)
	END DO

	DO i=1,NUM_AN_HeH
		F_ANGLE_HeH(i) = table_HeH(i,2)
	END DO

	DO i=1,NUM_AN_HeHe
		F_ANGLE_HeHe(i) = table_HeHe(i,2)
	END DO

	!! File File Probability Density Matrix and File Energy Array from Table
	DO i=1,NUM_EN_HeO
		DO j=1,NUM_AN_HeO
			k = (i-1)*NUM_AN_HeO + j
			F_PD_HeO(i,j) = table_HeO(k,3)
		END DO
		F_ENERGY_HeO(i) = table_HeO(k,1)
	END DO

	DO i=1,NUM_EN_HeH
		DO j=1,NUM_AN_HeH
			k = (i-1)*NUM_AN_HeH + j
			F_PD_HeH(i,j) = table_HeH(k,3)
		END DO
		F_ENERGY_HeH(i) = table_HeH(k,1)
	END DO

	DO i=1,NUM_EN_HeHe
		DO j=1,NUM_AN_HeHe
			k = (i-1)*NUM_AN_HeHe + j
			F_PD_HeHe(i,j) = table_HeHe(k,3)
		END DO
		F_ENERGY_HeHe(i) = table_HeHe(k,1)
	END DO

	!! Free memory from tables
	DEALLOCATE(table_HeO, table_HeH, table_HeHe)
	CLOSE(10)
	CLOSE(11)
	CLOSE(12)

	34 FORMAT(A, ES10.2)
	35 FORMAT(A, I6, A, I5) 

	IF (myid == 0) THEN
		WRITE(*,*)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A)')      'He + O TABLE'
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,35)         "N_ENERGIES: ", NUM_EN_HeO, " N_ANGLES: ", NUM_AN_HeO
		WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY_HeO(1)
		WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY_HeO(NUM_EN_HeO)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A)')      'He + H TABLE'
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,35)         "N_ENERGIES: ", NUM_EN_HeH, " N_ANGLES: ", NUM_AN_HeH
		WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY_HeH(1)
		WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY_HeH(NUM_EN_HeH)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A)')      'He + He TABLE'
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,35)         "N_ENERGIES: ", NUM_EN_HeHe, " N_ANGLES: ", NUM_AN_HeHe
		WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY_HeHe(1)
		WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY_HeHe(NUM_EN_HeHe)
		WRITE(*,'(a)')      '******************************************************'
		WRITE(*,'(A,F5.2)') 'Memory Footprint HeO   [MB]: ', np*s1*MB_per_real
		WRITE(*,'(A,F5.2)') 'Memory Footprint HeH   [MB]: ', np*s2*MB_per_real
		WRITE(*,'(A,F5.2)') 'Memory Footprint HeHe  [MB]: ', np*s3*MB_per_real
		WRITE(*,'(A,F5.2)') 'Total Memory Footprint [MB]: ', np*(s1+s2+s3)*MB_per_real	
		WRITE(*,'(A)')      '------------------------------------------------------'
		WRITE(*,*)
	END IF

END SUBROUTINE pd_data_read

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE index_EN(Enow,Table,E_index)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Find the index for the energy for which to use
! the dcs's from the table
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE tables, ONLY : F_ENERGY_HeO, F_ENERGY_HeH, F_ENERGY_HeHe, F_ENERGY_HH
	USE universal, ONLY : F_ENERGY
	
	IMPLICIT NONE
	
	CHARACTER(LEN=4)					:: Table
	REAL(KIND=8)							:: Enow
	REAL(KIND=8)							:: E

	INTEGER										:: E_index
	INTEGER										:: i
	INTEGER										:: z

	z = 0
	i = 0 
	DO WHILE (z == 0)
		i = i + 1

		IF (TRIM(Table) .EQ. 'HeO') THEN
			IF (i .EQ. SIZE(F_ENERGY_HeO))  z = 1
			E = F_ENERGY_HeO(i)
		ELSE IF ( (TRIM(Table) .EQ. 'HeH') .OR. (TRIM(Table) .EQ. 'HHe') ) THEN
			IF (i .EQ. SIZE(F_ENERGY_HeH))  z = 1
			E = F_ENERGY_HeH(i)
		ELSE IF (TRIM(Table) .EQ. 'HeHe') THEN 
			IF (i .EQ. SIZE(F_ENERGY_HeHe)) z = 1
			E = F_ENERGY_HeHe(i)
		ELSE IF (TRIM(Table) .EQ. 'HH') THEN
			IF (i .EQ. SIZE(F_ENERGY_HH))   z = 1
			E = F_ENERGY_HH(i)
		ELSE ! use scaling cross sections
			IF (i .EQ. SIZE(F_ENERGY))      z = 1
			E = F_ENERGY(i)
		END IF

		IF (E .GE. Enow) THEN
			z = 1
		END IF
	END DO 
	
	E_index = i

END SUBROUTINE index_EN

