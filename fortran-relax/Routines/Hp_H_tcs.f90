
SUBROUTINE test_Hp_H_tcs

	USE tables, ONLY : HpH_E, HpH_tcs, NUM_HpH_tcs

	IMPLICIT NONE

	REAL(KIND=8)		:: E(29), Ein, Efn, dE, tcs
	INTEGER					:: i, N

	N = 29
	
	CALL read_HpH_tcs_table

  OPEN(UNIT=200,FILE="../Data/test_H+_H_elastic_tcs.dat",ACCESS="APPEND")

	E = (/1.0D-2,2.5D-2,5.0D-2,7.5D-2,1.0D-1,2.5D-1,5.0D-1,7.5D-1, &
	&     1.0D0,2.5D0,5.0D0,7.5D0,1.0D1,2.5D1,5.0D1,7.5D1, &
	&     1.0D2,2.5D2,5.0D2,7.5D2,1.0D3,2.5D3,5.0D3,7.5D3, &
	&			1.0D4,2.5D4,5.0D4,7.5D4,1.0D5/)	

	DO i=1,N
		CALL Hp_H_tcs(E(i),tcs)
		WRITE(200,*) E(i), tcs	
	END DO

	CLOSE(200)

	CALL clean_HpH_tcs_table

END SUBROUTINE test_Hp_H_tcs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get H^+ + H elastic total cross section pulling from
! tables provided by Schultz, Krstic, Lee & Raymond 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Hp_H_tcs( E, tcs )
	USE tables, ONLY : HpH_E, HpH_tcs, NUM_HpH_tcs

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E	 ! CM energy [eV]

	!! Outputs
	REAL(KIND=8)		:: tcs ! [m^2]

	!! Internal
	REAL(KIND=8)		:: randy, prob_now, x0, x1, y0, y1, m, alt_now, CC, highest, lowest
	REAL(KIND=8)		:: C1, C2, now
	INTEGER					:: GOING, C

	CC        = 2.80D-17						! [AU^2] -> [cm^2]
	lowest    = HpH_E(1)						! [eV]
	highest   = HpH_E(NUM_HpH_tcs) 	! [eV]
	C 				= 1
	now 			= HpH_E(C)

	DO WHILE ( (E .GT. now) .AND. (C .LT. NUM_HpH_tcs) )
		C 		= C + 1
		now 	= HpH_E(C)
	END DO
	IF (C .GT. 1) THEN
		x0  = HpH_E(C-1)
		y0  = HpH_tcs(C-1)
	ELSE
		x0  = HpH_E(1)
		y0  = HpH_tcs(1)
	END IF
	x1  = HpH_E(C)
	y1  = HpH_tcs(C)
	m   = (y1-y0)/(x1-x0)
	tcs = y0 + m*(E-x0)	
	tcs = tcs*CC
	tcs = tcs*1.0D-4

END SUBROUTINE Hp_H_tcs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE read_HpH_tcs_table
	USE tables, ONLY : HpH_E, HpH_tcs, NUM_HpH_tcs
	USE mpi_info, ONLY : myid, ierr

	IMPLICIT NONE
	
	INCLUDE 'mpif.h'

	REAL(KIND=8)	:: dummy

	INTEGER				:: i

	IF (myid .EQ. 0) THEN
	
		!! open table file for reading
		OPEN(UNIT=72, FILE="../Tables/Hp_H_Elastic_TCS.dat", STATUS="old", ACTION="read")	

		!! read in number of lines in table
		READ(72,*) NUM_HpH_tcs

		ALLOCATE( HpH_E(NUM_HpH_tcs), HpH_tcs(NUM_HpH_tcs) )

		DO i=1,NUM_HpH_tcs
			READ(72,*) HpH_E(i), dummy, HpH_tcs(i)
		END DO

		CLOSE(72)
	
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_HpH_tcs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	IF (myid .NE. 0) THEN
		ALLOCATE( HpH_E(NUM_HpH_tcs), HpH_tcs(NUM_HpH_tcs))	
	END IF
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( HpH_E, NUM_HpH_tcs,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( HpH_tcs, NUM_HpH_tcs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE read_HpH_tcs_table

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE clean_HpH_tcs_table
	USE tables, ONLY : HpH_E, HpH_tcs
	
	IMPLICIT NONE

	DEALLOCATE( HpH_E, HpH_tcs )

END SUBROUTINE clean_HpH_tcs_table


