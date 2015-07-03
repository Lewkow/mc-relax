!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module to use for ENA production routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE ENA

	REAL(KIND=8),PARAMETER	:: z_in 		= 2000.0D3	!! Initial height (high) [m]
	REAL(KIND=8),PARAMETER	:: z_fn 		= 100.0D3		!! Final height   (low)  [m]
	REAL(KIND=8),PARAMETER	:: E_sw 		= 1000.0D0  !! Energy of SW ions [eV]
	REAL(KIND=8),PARAMETER	:: V_sw 		= 200.0D3		!! Velocity of SW [m/s]
	REAL(KIND=8),PARAMETER	:: N_sw 		= 5.0D6			!! Density of SW ions at 1 AU [1/m^3]
	REAL(KIND=8),PARAMETER	:: H_sw			= 0.92D0		!! fraction of SW which is protons
	REAL(KIND=8),PARAMETER	:: He_sw		= 0.08D0		!! fraction of SW which is alpha particles
	REAL(KIND=8),PARAMETER	:: B0   		= 0.0D0		!! Strength of LISM magnetic field [T]
!	REAL(KIND=8),PARAMETER	:: B0   		= 2.7D-10		!! Strength of LISM magnetic field [T]
	INTEGER,PARAMETER				:: Nz 			= 1000			!! Number of heights to compute 
	INTEGER									:: L_ion, L_targ			  !! Length of ion and target character array
	REAL(KIND=8),PARAMETER	:: lism_start = 1.5D11	!! start integration at 1 AU [m]
!	REAL(KIND=8),PARAMETER	:: r0       = 1.4959D11	!! Initial ENA production radius [m] 

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: Prod_Height_Alt
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: Prod_Height_Cum
	INTEGER																:: N_Prod_Height

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE read_ENA_prod_height_table
	USE mpi_info
	USE planet

	IMPLICIT NONE

	INTEGER :: i

	INCLUDE 'mpif.h'

	IF ( myid .EQ. 0 ) THEN 
		IF ( TRIM(Proj) .EQ. 'H' ) THEN
			IF ( ATMOSPHERE .EQ. 0 ) THEN
				OPEN(UNIT=55,FILE='../Tables/H_Min_Production_Height.dat',STATUS='old',ACTION='read')
			ELSE IF ( ATMOSPHERE .EQ. 1 ) THEN
				OPEN(UNIT=55,FILE='../Tables/H_Mean_Production_Height.dat',STATUS='old',ACTION='read')
			ELSE IF ( ATMOSPHERE .EQ. 2 ) THEN
				OPEN(UNIT=55,FILE='../Tables/H_Max_Production_Height.dat',STATUS='old',ACTION='read')
			END IF
		ELSE IF ( TRIM(Proj) .EQ. 'He' ) THEN
			IF ( ATMOSPHERE .EQ. 0 ) THEN
				OPEN(UNIT=55,FILE='../Tables/He_Min_Production_Height.dat',STATUS='old',ACTION='read')
			ELSE IF ( ATMOSPHERE .EQ. 1 ) THEN
				OPEN(UNIT=55,FILE='../Tables/He_Mean_Production_Height.dat',STATUS='old',ACTION='read')
			ELSE IF ( ATMOSPHERE .EQ. 2 ) THEN
				OPEN(UNIT=55,FILE='../Tables/He_Max_Production_Height.dat',STATUS='old',ACTION='read')
			END IF
		END IF
		
		READ(55,*) N_Prod_Height

		ALLOCATE( Prod_Height_Alt(N_Prod_Height), Prod_Height_Cum(N_Prod_Height) )
	
		DO i=1,N_Prod_Height
			READ(55,*) Prod_Height_Alt(i), Prod_Height_Cum(i)
		END DO

		CLOSE(55)

	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( N_Prod_Height, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	IF (myid .NE. 0) ALLOCATE( Prod_Height_Alt(N_Prod_Height), Prod_Height_Cum(N_Prod_Height) )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Prod_Height_Alt, N_Prod_Height, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Prod_Height_Cum, N_Prod_Height, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE read_ENA_prod_height_table

!! 

SUBROUTINE clean_ENA_prod_height_table

	DEALLOCATE( Prod_Height_Alt, Prod_Height_Cum )

END SUBROUTINE clean_ENA_prod_height_table

!!

SUBROUTINE rand_ENA_altitude( z )

	IMPLICIT NONE

	!! Outputs
	REAL(KIND=8)	:: z 	!! random altitude [m]

	!! Internal
	REAL(KIND=8)	:: lfg	!! random number generator
	REAL(KIND=8)	:: randy, x1, x2, y1, y2, m, y0
	INTEGER				:: i

	randy = lfg()
	i 		= 1

	DO WHILE ( Prod_Height_Cum(i) .LE. randy )
		i = i + 1
	END DO

	x1 = Prod_Height_Cum(i-1)
	x2 = Prod_Height_Cum(i)
	y1 = Prod_Height_Alt(i-1)
	y2 = Prod_Height_Alt(i)
	m  = (y2-y1)/(x2-x1)
	y0 = y1 - m*x1
	z  = y0 + m*randy
	z  = z*1.0D3

END SUBROUTINE rand_ENA_altitude

END MODULE ENA

