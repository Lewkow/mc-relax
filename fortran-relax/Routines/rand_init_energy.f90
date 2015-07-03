
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get initial random energy for ENA particle on Mars
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE rand_init_energy( MP, r_energy )
	USE physics_constants, ONLY : UTOKG, EVTOJ

	IMPLICIT NONE
	!! Inputs
	REAL(KIND=8)		:: MP				 ! projectile mass [u]
	
	!! Outputs
	REAL(KIND=8)		:: r_energy	 ! [eV]

	!! Internal
	REAL(KIND=8)		:: r, E, V

	!! get random solar wind velocity
	CALL get_Vsw(V)
	
	!! convert V from km/s to m/s
	V = V*1000.0D0

	!! convert velocity to eV
	E = 0.5D0*MP*UTOKG*V*V/EVTOJ
	r_energy = E

END SUBROUTINE rand_init_energy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read from Vsw arrays and produce a velocity based on 
! a random number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_Vsw( V )
	USE tables, ONLY : Vsw, Vsw_Prob

	IMPLICIT NONE

	!! Outputs
	REAL(KIND=8)		:: V ! [km/s]

	!! Internal
	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: randy, prob_now, x0, x1, y0, y1, m
	INTEGER					:: GOING, C

	randy = lfg()
	
	C 				= 1
	prob_now 	= Vsw_Prob(C)/100.0D0
	DO WHILE (randy .GT. prob_now) 
		C 				= C + 1
		prob_now 	= Vsw_Prob(C)/100.0D0
	END DO
	x0 = Vsw_Prob(C-1)/100.0D0
	x1 = Vsw_Prob(C)/100.0D0
	y0 = Vsw(C-1)
	y1 = Vsw(C)

	m = (y1-y0)/(x1-x0)
	V = y0 + m*(randy-x0)	

END SUBROUTINE get_Vsw

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Parse the Vsw data file and save parameters to Module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_Vsw_table
	USE tables, ONLY : Vsw, Vsw_Prob, NUM_Vsw
	USE mpi_info
	
	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER			:: i

	IF (myid .EQ. 0) THEN
	
		!! open table file for reading
		OPEN(UNIT=72, FILE="../Tables/SW_Velocity.dat", STATUS="old", ACTION="read")	

		!! read in number of lines in table
		READ(72,*) NUM_Vsw

		ALLOCATE( Vsw(NUM_Vsw), Vsw_Prob(NUM_Vsw) )

		DO i=1,NUM_Vsw
			READ(72,*) Vsw_Prob(i), Vsw(i)
		END DO

		CLOSE(72)

	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_Vsw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	IF (myid .NE. 0) ALLOCATE( Vsw(NUM_Vsw), Vsw_Prob(NUM_Vsw) )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Vsw, NUM_Vsw, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Vsw_Prob, NUM_Vsw, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	

END SUBROUTINE read_Vsw_table

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Test the random energy routine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE test_rand_energy
	
	IMPLICIT NONE

	INTEGER 			:: i, N
	REAL(KIND=8)	:: E 

!	CALL read_Vsw_table

	OPEN(UNIT=88, FILE="../Data/Test_Rand_Init_Energy.dat", ACCESS="APPEND")

	N = 1000000

	DO i=1,N
		CALL rand_init_energy( 1.0D0, E )	
		WRITE(88,*) E
	END DO	

	CLOSE(88)

END SUBROUTINE test_rand_energy

