
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read from Vsw arrays and produce a velocity based on 
! a random number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_electron_density( Z, den )
	USE tables, ONLY : electron_den, electron_alt, NUM_e_den

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: Z	 ! [m]

	!! Outputs
	REAL(KIND=8)		:: den ! [1/m^3]

	!! Internal
	REAL(KIND=8)		:: randy, prob_now, x0, x1, y0, y1, m, alt_now, CC, highest, lowest
	REAL(KIND=8)		:: C1, C2
	INTEGER					:: GOING, C

	CC        = 1.0D-3	
	randy 		= Z*CC
	lowest    = electron_alt(1)					! [km]
	highest   = electron_alt(NUM_e_den) ! [km]
	C 				= 1
	alt_now 	= electron_alt(C)

!	write(*,*) 'highest: z', highest, randy

	IF ( randy .LE. highest ) THEN

		DO WHILE ( (randy .GT. alt_now) .AND. (C .LT. SIZE(electron_alt(:))))
			C 				= C + 1
			alt_now 	= electron_alt(C)
		END DO

		IF (C .GT. 1) THEN
			x0  = electron_alt(C-1)
			y0  = electron_den(C-1)
		ELSE
			x0  = electron_alt(1)
			y0  = electron_den(1)
		END IF
		x1  = electron_alt(C)
		y1  = electron_den(C)
		m   = (y1-y0)/(x1-x0)
		den = y0 + m*(randy-x0)	

!		write(*,*) 'low z: n: ', z/1.0d3, den

	END IF

	IF ( randy .GT. highest ) THEN	

		!! density equal to exponential decay fit to last 2 data points of 
		!! electron density table
		!! n = 10^( -(z-C1)/C2 )

		C1  = 45.8D0
		C2  = 376.0D0	
		m   = -(randy - C2)/C1
		den = 10.0D0**m

!		write(*,*) 'hig z: n: ', randy, den

	END IF

	IF ( randy .LE. lowest ) den = electron_den(1)

	den = den*1.0D6	! convert from [1/cm^3] -> [1/m^3]
!	WRITE(*,*) 'x0, x1, y0, y1, m, den ', x0, x1, y0, y1, m, den
!	WRITE(*,*) randy, den

END SUBROUTINE get_electron_density

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Parse the Vsw data file and save parameters to Module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_e_den_table
	USE tables, ONLY : electron_den, electron_alt, NUM_E_den
	USE mpi_info

	IMPLICIT NONE
	
	INCLUDE 'mpif.h'

	INTEGER			:: i

	IF (myid .EQ. 0) THEN
	
		!! open table file for reading
		OPEN(UNIT=72, FILE="../Tables/electron_density.dat", STATUS="old", ACTION="read")	

		!! read in number of lines in table
		READ(72,*) NUM_e_den

		ALLOCATE( electron_alt(NUM_e_den), electron_den(NUM_e_den) )

		DO i=1,NUM_e_den
			READ(72,*) electron_alt(i), electron_den(i)
		END DO

		CLOSE(72)
	
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NUM_E_den, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	IF (myid .NE. 0) THEN
		ALLOCATE( electron_alt(NUM_e_den), electron_den(NUM_e_den))	
	END IF
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( electron_alt(:), NUM_e_den,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( electron_den(:), NUM_e_den, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE read_e_den_table

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Parse the Vsw data file and save parameters to Module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE clean_e_den_table
	USE tables, ONLY : electron_den, electron_alt
	
	IMPLICIT NONE

	DEALLOCATE( electron_alt, electron_den )

END SUBROUTINE clean_e_den_table


