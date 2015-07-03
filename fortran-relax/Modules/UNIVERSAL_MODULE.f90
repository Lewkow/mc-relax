
MODULE universal

	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE	:: table_X_H, table_X_He, table_X_O, table_X_Ar	
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE	:: table_X_H2, table_X_N2, table_X_CO, table_X_CO2	
	REAL(KIND=8),DIMENSION(:), ALLOCATABLE	:: F_ENERGY,  F_ANGLE
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_O,  F_PD_X_Ar
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_H2, F_PD_X_N2
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_CO, F_PD_X_CO2

	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_ENERGY
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_O,  TCS_X_Ar
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_H2, TCS_X_N2
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_CO, TCS_X_CO2

	INTEGER																	:: NL_TCS
	INTEGER																	:: NL_X_H, NL_X_He, NL_X_O, NL_X_Ar
	INTEGER																	:: NL_X_H2, NL_X_N2, NL_X_CO, NL_X_CO2
	INTEGER																	:: NUM_ENG, NUM_ANG

  REAL,PARAMETER                          :: MB = 8.0D-6
  REAL,PARAMETER                          :: GB = 8.0D-9

	REAL(KIND=8)		:: tau_0, tau_b

CONTAINS


!################################################
!################################################

SUBROUTINE read_uni_tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in Angular Probability data table from file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE planet, ONLY : PROJ
  USE mpi_info

  IMPLICIT NONE

	INCLUDE 'mpif.h'

  REAL(KIND=8)                            :: E0
  REAL(KIND=8)                            :: Enow

  INTEGER                                 :: N
  INTEGER                                 :: i
  INTEGER                                 :: j
  INTEGER                                 :: k, s1, s2, s3, s4

	IF (myid .EQ. 0) THEN
		IF (PROJ .NE. 'H' .AND. PROJ .NE. 'He') WRITE(*,*) 'PROJ ', PROJ, ' NOT KNOWN'
		IF (PROJ .EQ. 'H') THEN
			!! open tables files to read from
 			OPEN(UNIT=10, FILE="../Tables/H_O_ANG.dat", STATUS="old", ACTION="read")
 			OPEN(UNIT=11, FILE="../Tables/H_Ar_ANG.dat", STATUS="old", ACTION="read")
 		 	OPEN(UNIT=12, FILE="../Tables/H_H2_ANG.dat", STATUS="old", ACTION="read")
 		 	OPEN(UNIT=13, FILE="../Tables/H_N2_ANG.dat", STATUS="old", ACTION="read")
 			OPEN(UNIT=14, FILE="../Tables/H_CO_ANG.dat", STATUS="old", ACTION="read")
 			OPEN(UNIT=15, FILE="../Tables/H_CO2_ANG.dat", STATUS="old", ACTION="read")
  		READ(10,*) NL_X_O
  		READ(11,*) NL_X_Ar
  		READ(12,*) NL_X_H2
  		READ(13,*) NL_X_N2
  		READ(14,*) NL_X_CO
  		READ(15,*) NL_X_CO2
  		ALLOCATE(table_X_O(NL_X_O ,3))
  		ALLOCATE(table_X_Ar(NL_X_Ar ,3))
  		ALLOCATE(table_X_H2(NL_X_H2,3))
  		ALLOCATE(table_X_N2(NL_X_N2,3))
  		ALLOCATE(table_X_CO(NL_X_CO ,3))
  		ALLOCATE(table_X_CO2(NL_X_CO2 ,3))
  		DO i=1,NL_X_O
				READ(10,*) table_X_O(i,1), table_X_O(i,2), table_X_O(i,3)
				READ(11,*) table_X_Ar(i,1), table_X_Ar(i,2), table_X_Ar(i,3)
				READ(12,*) table_X_H2(i,1), table_X_H2(i,2), table_X_H2(i,3)
				READ(13,*) table_X_N2(i,1), table_X_N2(i,2), table_X_N2(i,3)
				READ(14,*) table_X_CO(i,1), table_X_CO(i,2), table_X_CO(i,3)
				READ(15,*) table_X_CO2(i,1), table_X_CO2(i,2), table_X_CO2(i,3)
  		END DO
		ELSE IF (PROJ .EQ. 'He') THEN
 			OPEN(UNIT=10, FILE="../Tables/He_Ar_ANG.dat", STATUS="old", ACTION="read")
 		 	OPEN(UNIT=11, FILE="../Tables/He_H2_ANG.dat", STATUS="old", ACTION="read")
 		 	OPEN(UNIT=12, FILE="../Tables/He_N2_ANG.dat", STATUS="old", ACTION="read")
 			OPEN(UNIT=13, FILE="../Tables/He_CO_ANG.dat", STATUS="old", ACTION="read")
 			OPEN(UNIT=14, FILE="../Tables/He_CO2_ANG.dat", STATUS="old", ACTION="read")
  		READ(10,*) NL_X_Ar
  		READ(11,*) NL_X_H2
  		READ(12,*) NL_X_N2
  		READ(13,*) NL_X_CO
  		READ(14,*) NL_X_CO2
  		ALLOCATE(table_X_Ar(NL_X_Ar ,3))
  		ALLOCATE(table_X_H2(NL_X_H2,3))
  		ALLOCATE(table_X_N2(NL_X_N2,3))
  		ALLOCATE(table_X_CO(NL_X_CO ,3))
  		ALLOCATE(table_X_CO2(NL_X_CO2 ,3))
  		DO i=1,NL_X_Ar
				READ(10,*) table_X_Ar(i,1), table_X_Ar(i,2), table_X_Ar(i,3)
				READ(11,*) table_X_H2(i,1), table_X_H2(i,2), table_X_H2(i,3)
				READ(12,*) table_X_N2(i,1), table_X_N2(i,2), table_X_N2(i,3)
				READ(13,*) table_X_CO(i,1), table_X_CO(i,2), table_X_CO(i,3)
				READ(14,*) table_X_CO2(i,1), table_X_CO2(i,2), table_X_CO2(i,3)
  		END DO
		ELSE 
			WRITE(*,*) 'Current Projectile: ', PROJ, ' not currently supported'
		END IF

	 	E0    = table_X_Ar(1,1)
  	Enow  = E0
 	 	N     = 0
 		DO WHILE (Enow == E0)
   		N    = N + 1
    	Enow = table_X_Ar(N,1)
  	END DO

  	NUM_ANG = N-1
  	NUM_ENG = NL_X_Ar/NUM_ANG
		ALLOCATE(F_ENERGY(NUM_ENG), F_ANGLE(NUM_ANG))
  	DO i=1,NUM_ANG
    	F_ANGLE(i)  = table_X_Ar(i,2)
  	END DO

		IF (PROJ .EQ. 'H') THEN
			ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))	
			s1 = SIZE(F_ENERGY) + SIZE(F_ANGLE) + 6*SIZE(F_PD_X_Ar)
  		DO i=1,NUM_ENG
    		DO j=1,NUM_ANG
      		k = (i-1)*NUM_ANG + j
      		F_PD_X_O(i,j) = table_X_O(k,3)
      		F_PD_X_Ar(i,j) = table_X_Ar(k,3)
      		F_PD_X_H2(i,j) = table_X_H2(k,3)
      		F_PD_X_N2(i,j) = table_X_N2(k,3)
      		F_PD_X_CO(i,j) = table_X_CO(k,3)
      		F_PD_X_CO2(i,j) = table_X_CO2(k,3)
    		END DO
    		F_ENERGY(i) = table_X_Ar(k,1)
  		END DO
  		DEALLOCATE(table_X_O, table_X_Ar, table_X_H2, table_X_N2, table_X_CO, table_X_CO2)
  		CLOSE(10)
  		CLOSE(11)
  		CLOSE(12)
  		CLOSE(13)
  		CLOSE(14)
  		CLOSE(15)
		ELSE IF (PROJ .EQ. 'He') THEN
			ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))	
			ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))	
			s1 = SIZE(F_ENERGY) + SIZE(F_ANGLE) + 5*SIZE(F_PD_X_Ar)
  		DO i=1,NUM_ENG
    		DO j=1,NUM_ANG
      		k = (i-1)*NUM_ANG + j
      		F_PD_X_Ar(i,j) = table_X_Ar(k,3)
      		F_PD_X_H2(i,j) = table_X_H2(k,3)
      		F_PD_X_N2(i,j) = table_X_N2(k,3)
      		F_PD_X_CO(i,j) = table_X_CO(k,3)
      		F_PD_X_CO2(i,j) = table_X_CO2(k,3)
    		END DO
    		F_ENERGY(i) = table_X_Ar(k,1)
  		END DO
  		DEALLOCATE(table_X_Ar, table_X_H2, table_X_N2, table_X_CO, table_X_CO2)
  		CLOSE(10)
  		CLOSE(11)
  		CLOSE(12)
  		CLOSE(13)
  		CLOSE(14)
		END IF

  	34 FORMAT(A, ES10.2)
  	35 FORMAT(A, I6, A, I5)

    WRITE(*,'(A)')      'UNIVERSAL TABLES READ'
		WRITE(*,'(a)')      '#################################################'
    WRITE(*,35)         "N_ENERGIES: ", NUM_ENG, " N_ANGLES: ", NUM_ANG
    WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY(1)
    WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY(NUM_ENG)
		WRITE(*,'(a)')      '#################################################'
    WRITE(*,'(A,F6.2)') 'Universal Memory Footprint [MB]: ', np*s1*MB
		WRITE(*,'(a)')      '#################################################'
  END IF ! rank = 0

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(NUM_ENG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(NUM_ANG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	IF (myid .NE. 0) THEN
		ALLOCATE(F_ENERGY(NUM_ENG), F_ANGLE(NUM_ANG))
		IF (PROJ .EQ. 'H') THEN
      ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))
		ELSE IF (PROJ .EQ. 'He') THEN
      ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))
      ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))
		END IF
	END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	CALL MPI_BCAST(F_ENERGY,NUM_ENG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)	
	CALL MPI_BCAST(F_ANGLE,NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)	

	IF (PROJ .EQ. 'H') THEN
		CALL MPI_BCAST(F_PD_X_O,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_Ar,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_H2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_N2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_CO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_CO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	ELSE IF (PROJ .EQ. 'He') THEN
		CALL MPI_BCAST(F_PD_X_Ar,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_H2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_N2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_CO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(F_PD_X_CO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	END IF

END SUBROUTINE read_uni_tables

!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!

SUBROUTINE clean_uni_tables
	USE planet, ONLY : PROJ

	IF (PROJ .EQ. 'H') THEN
  	DEALLOCATE(F_PD_X_O)
  	DEALLOCATE(F_PD_X_Ar)
  	DEALLOCATE(F_PD_X_H2)
  	DEALLOCATE(F_PD_X_N2)
  	DEALLOCATE(F_PD_X_CO)
  	DEALLOCATE(F_PD_X_CO2)
	ELSE IF (PROJ .EQ. 'He') THEN
  	DEALLOCATE(F_PD_X_Ar)
  	DEALLOCATE(F_PD_X_H2)
  	DEALLOCATE(F_PD_X_N2)
  	DEALLOCATE(F_PD_X_CO)
  	DEALLOCATE(F_PD_X_CO2)
	END IF

	DEALLOCATE(F_ENERGY,F_ANGLE)

END SUBROUTINE clean_uni_tables

!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_uni_tcs_tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in Angular Probability data table from file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE planet, ONLY : PROJ
  USE mpi_info

  IMPLICIT NONE

  INCLUDE 'mpif.h'

	INTEGER						:: i
	REAL(KIND=8)			:: dummy

	IF (myid .EQ. 0) THEN
		IF (PROJ .EQ. 'H ') THEN
      OPEN(UNIT=10, FILE="../Tables/H_O_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=11, FILE="../Tables/H_Ar_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=12, FILE="../Tables/H_H2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=13, FILE="../Tables/H_N2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=14, FILE="../Tables/H_CO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=15, FILE="../Tables/H_CO2_TCS.dat", STATUS="old", ACTION="read")
      READ(10,*) NL_TCS
      READ(11,*) 
      READ(12,*) 
      READ(13,*) 
      READ(14,*) 
      READ(15,*) 
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			DO i=1,NL_TCS
				READ(10,*) TCS_ENERGY(i), TCS_X_O(i)
				READ(11,*) dummy,         TCS_X_Ar(i)
				READ(12,*) dummy,         TCS_X_H2(i)
				READ(13,*) dummy,         TCS_X_N2(i)
				READ(14,*) dummy,         TCS_X_CO(i)
				READ(15,*) dummy,         TCS_X_CO2(i)
			END DO ! i
			CLOSE(10)
			CLOSE(11)
			CLOSE(12)
			CLOSE(13)
			CLOSE(14)
			CLOSE(15)
		ELSE IF (PROJ .EQ. 'He') THEN
      OPEN(UNIT=11, FILE="../Tables/He_Ar_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=12, FILE="../Tables/He_H2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=13, FILE="../Tables/He_N2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=14, FILE="../Tables/He_CO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=15, FILE="../Tables/He_CO2_TCS.dat", STATUS="old", ACTION="read")
      READ(11,*) NL_TCS
      READ(12,*) 
      READ(13,*) 
      READ(14,*) 
      READ(15,*) 
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			DO i=1,NL_TCS
				READ(11,*) TCS_ENERGY(i), TCS_X_Ar(i)
				READ(12,*) dummy,         TCS_X_H2(i)
				READ(13,*) dummy,         TCS_X_N2(i)
				READ(14,*) dummy,         TCS_X_CO(i)
				READ(15,*) dummy,         TCS_X_CO2(i)
			END DO ! i
			CLOSE(11)
			CLOSE(12)
			CLOSE(13)
			CLOSE(14)
			CLOSE(15)
		ELSE
			WRITE(*,*) 'Projectile ', PROJ, ' not currently supported in tables'
		END IF ! PROJ
	END IF ! myid = 0

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NL_TCS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF ( myid .NE. 0) THEN
		IF (PROJ .EQ. 'H ') THEN
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
		ELSE IF (PROJ .EQ. 'He') THEN
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
		END IF
	END IF
		
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	CALL MPI_BCAST(TCS_ENERGY, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

	IF (PROJ .EQ. 'H ') THEN
		CALL MPI_BCAST(TCS_X_O,  NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_Ar, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_H2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_N2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	ELSE IF (PROJ .EQ. 'He') THEN	
		CALL MPI_BCAST(TCS_X_Ar, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_H2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_N2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE read_uni_tcs_tables

!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!

SUBROUTINE clean_uni_tcs_tables
	USE planet, ONLY : PROJ

	IMPLICIT NONE

	IF (PROJ .EQ. 'H ') THEN
		DEALLOCATE( TCS_X_O, TCS_X_Ar, TCS_X_H2, TCS_X_N2, TCS_X_CO, TCS_X_CO2 )
	ELSE IF (PROJ .EQ. 'He') THEN
		DEALLOCATE( TCS_X_Ar, TCS_X_H2, TCS_X_N2, TCS_X_CO, TCS_X_CO2 )
	END IF

END SUBROUTINE clean_uni_tcs_tables

END MODULE universal


