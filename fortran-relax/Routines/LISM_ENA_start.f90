
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read parameters for ENA start probabilies and 
! save them in the lism module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lism_ena_table_read
	USE lism
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER					:: i, N_H, N_He

	IF (myid .EQ. 0) THEN
		OPEN(UNIT=50,FILE='../Tables/LISM_H_ENA_Start.dat')	
		OPEN(UNIT=60,FILE='../Tables/LISM_He_ENA_Start.dat')	

		READ(50,*) N_H	
		READ(60,*) N_He	

		ALLOCATE( H_ENA_PROB(N_H), H_ENA_DIST(N_H) ) 
		ALLOCATE( He_ENA_PROB(N_He), He_ENA_DIST(N_He) ) 

		DO i=1,N_H
			READ(50,*) H_ENA_PROB(i), H_ENA_DIST(i)	
		END DO

		DO i=1,N_He
			READ(60,*) He_ENA_PROB(i), He_ENA_DIST(i)	
		END DO

		CLOSE(50)
		CLOSE(60)
	END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	CALL MPI_BCAST(N_H,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(N_He,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	IF (myid .NE. 0) THEN
		ALLOCATE( H_ENA_PROB(N_H), H_ENA_DIST(N_H) ) 
		ALLOCATE( He_ENA_PROB(N_He), He_ENA_DIST(N_He) ) 
	END IF

	CALL MPI_BCAST( H_ENA_PROB, N_H, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST( He_ENA_PROB, N_He, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST( H_ENA_DIST, N_H, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST( He_ENA_DIST, N_He, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE lism_ena_table_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Determine starting location and velocity
! for ENA in LISM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lism_ena_start( R, Proj )
	USE physics_constants
	USE lism

	IMPLICIT NONE

	!! Outputs
	REAL(KIND=8)		:: R(3)		! starting location [x,y,z] in [m]
	CHARACTER(Len=2):: Proj   ! projectile name (H or He)

	!! Internal
	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: phi, theta, dist, randy, randy1, randy2
	INTEGER					:: i, N


	IF (ENA_TYPE .EQ. 0) THEN
		Proj   = 'H '
		randy2 = lfg()
		N      = 0
		i      = 1
		DO WHILE (N .EQ. 0)
			IF ( (H_ENA_PROB(i) .GE. randy2) .OR. (i .GE. SIZE(H_ENA_PROB(:))) ) THEN
				dist = H_ENA_DIST(i)
				N    = 1
			END IF
			i = i+1	
		END DO
	ELSE
		Proj   = 'He'
		randy2 = lfg()
		N      = 0
		i      = 1
		DO WHILE (N .EQ. 0)
			IF ( (He_ENA_PROB(i) .GE. randy2) .OR. (i .GE. SIZE(He_ENA_PROB(:))) ) THEN
				dist = He_ENA_DIST(i)
				N    = 1
			END IF
			i = i+1	
		END DO
	END IF


	!! get random location on sphere
	randy1 = lfg()
	randy2 = lfg()
	theta  = randy1*PI
	phi    = randy2*2.0D0*PI
	R(1)   = dist*SIN(theta)*COS(phi)	
	R(2)   = dist*SIN(theta)*SIN(phi)	
	R(3)   = dist*COS(theta)

END SUBROUTINE lism_ena_start

