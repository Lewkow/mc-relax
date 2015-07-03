
SUBROUTINE input_reader
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Reads input values and stores them to input module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE inputs
	USE physics_constants
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	IF (myid .EQ. 0) THEN

		OPEN(11, FILE="../Inputs/keys.in", STATUS="old", ACTION="read")

		! Comments
		READ(11,*)
		READ(11,*)
		READ(11,*)
		READ(11,*)
		READ(11,*)
		READ(11,*)

		! LISM calculation ON/OFF
		READ(11,*) DO_LISM

		! planet calculation ON/OFF
		READ(11,*) DO_planet

		! ENA production ON/OFF
		READ(11,*) DO_ENA_Production

		! Calculate average scattering angles
		READ(11,*) DO_Average_Scattering_Angle

		! Do MC SH escape simulation
		READ(11,*) DO_Escape_MC

		! Do transparency integration for SH Escape
		READ(11,*) DO_Escape_Trans

		CLOSE(11)

	END IF ! myid = 0

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(DO_LISM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(DO_planet,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(DO_ENA_Production,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(DO_Average_Scattering_Angle,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(DO_Escape_MC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(DO_Escape_Trans,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

END SUBROUTINE

