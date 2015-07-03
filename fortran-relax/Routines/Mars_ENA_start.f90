
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Test random SH init postions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE test_mars_sh_start_alt

	INTEGER				:: N, i
	REAL(KIND=8)	:: z

	N = 100000

	DO i=1,N
		CALL mars_sh_start( z )
		WRITE(666,*) z
	END DO	

END SUBROUTINE test_mars_sh_start_alt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read parameters for ENA start probabilies and 
! save them in the planet module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mars_ena_table_read
	USE planet
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER					:: i, N_H, N_He

	IF (myid .EQ. 0) THEN

		OPEN(UNIT=50,FILE='../Tables/Mars_H_ENA_Start.dat')	
		OPEN(UNIT=60,FILE='../Tables/Mars_He_ENA_Start.dat')	

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

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( N_H, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( N_He, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	IF ( myid .NE. 0 ) THEN
		ALLOCATE( H_ENA_PROB(N_H), H_ENA_DIST(N_H) ) 
		ALLOCATE( He_ENA_PROB(N_He), He_ENA_DIST(N_He) ) 
	END IF
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( H_ENA_PROB, N_H, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( He_ENA_PROB, N_He, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

END SUBROUTINE mars_ena_table_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read parameters for SH atoms/molecules start probabilies and 
! save them in the planet module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mars_sh_table_read
	USE planet
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: Dummy_X, Dummy_Y
	REAL(KIND=8)													:: Tot, dZ, Tot_now
	INTEGER																:: i, N_H, N_He

	IF (myid .EQ. 0) THEN

		IF (ENA_SH .EQ. 1) THEN
			OPEN(UNIT=50,FILE='../Tables/H_Mean_SH_He.dat')
		ELSE IF (ENA_SH .EQ. 4) THEN
			OPEN(UNIT=50,FILE='../Tables/He_Mean_SH_He.dat')
		END IF

		READ(50,*) N_SH_Z0

		ALLOCATE( SH_Z0_Z(N_SH_Z0), SH_Z0_P(N_SH_Z0) ) 
		ALLOCATE( Dummy_X(N_SH_Z0), Dummy_Y(N_SH_Z0) ) 

		DO i=1,N_SH_Z0
			READ(50,*) Dummy_Y(i), Dummy_X(i)	
		END DO
		CLOSE(50)

		dZ 					= Dummy_X(2)-Dummy_X(1)	! [km]
		Dummy_Y(:) 	= Dummy_Y(:)/dZ	

		Tot = 0.0D0
		DO i=1,N_SH_Z0
			Tot = Tot + Dummy_Y(i)*dZ
		END DO

		Dummy_Y(:) 	= Dummy_Y(:)/Tot	
		Tot_now 		= 0.0D0

		DO i=1,N_SH_Z0
			Tot_now = Tot_now + Dummy_Y(i)*dZ
			SH_Z0_Z(i) = Dummy_X(i)
			SH_Z0_P(i) = Tot_now
	!		WRITE(*,*) SH_Z0_Z(i), SH_Z0_P(i)
		END DO
	!	WRITE(*,*) 'TOT: ', Tot_now

		DEALLOCATE( Dummy_X, Dummy_Y )

	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( N_SH_Z0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	IF ( myid .NE. 0 ) THEN
		ALLOCATE( SH_Z0_Z(N_SH_Z0), SH_Z0_P(N_SH_Z0) ) 
	END IF
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( SH_Z0_Z, N_SH_Z0, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( SH_Z0_P, N_SH_Z0, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

END SUBROUTINE mars_sh_table_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Determine starting height for SH atoms/molecules
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mars_sh_start( z )
	USE physics_constants
	USE planet

	IMPLICIT NONE

	!! Outputs
	REAL(KIND=8)		:: z			! starting height above planet in [m]

	!! Internal
	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: phi, theta, dist, randy, randy1, randy2, z1, z2
	INTEGER					:: i, N

	randy2 = lfg()
	N      = 0
	i      = 1
	DO WHILE (N .EQ. 0)
		IF ( (SH_Z0_P(i) .GE. randy2) .OR. (i .GE. SIZE(SH_Z0_P(:))) ) THEN
			z1   = SH_Z0_Z(i)
			DO WHILE ((SH_Z0_Z(i) .EQ. z1) .AND. (i .LT. SIZE(SH_Z0_Z(:))))
				i = i+1
			END DO
			z2   	 	= SH_Z0_Z(i)
			randy1 	= lfg()
			dist 		= z1 + randy1*(z2-z1)
			N    		= 1
		END IF
		i = i+1	
	END DO

	z = dist*1.0D3	! convert form [km] to [m]

END SUBROUTINE mars_sh_start


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Determine starting height for ENAs above Mars
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mars_ena_start( z )
	USE physics_constants
	USE planet

	IMPLICIT NONE

	!! Outputs
	REAL(KIND=8)		:: z			! starting height above planet in [m]

	!! Internal
	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: phi, theta, dist, randy, randy1, randy2
	INTEGER					:: i, N

	IF (Proj .EQ. 'H  ') THEN
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
	ELSE IF (Proj .EQ. 'He ') THEN	
		randy2 = lfg()
		N      = 0
		i      = 1
		DO WHILE (N .EQ. 0)
			IF ( (He_ENA_PROB(i) .GE. randy2) .OR. (i .GE. SIZE(HE_ENA_DIST(:))) ) THEN
				dist = He_ENA_DIST(i)
				N    = 1
			END IF
			i = i+1	
		END DO
	END IF

	z = dist*1.0D3	! convert form [km] to [m]

END SUBROUTINE mars_ena_start

