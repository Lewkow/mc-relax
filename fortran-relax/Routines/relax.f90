
PROGRAM mc_relax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Main program for the MC relaxation calc
! to calculate the energy of a particle as
! a function of number of collisions using 
! random scattering angles obtained from 
! prob densities from ASTROSCATT. 
! 
! This is for a uniform bath gas for now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants
	USE inputs
	USE rand_seed
	USE tables
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	REAL(KIND=8)		:: PD_Start	! Starting Time 
	REAL(KIND=8)		:: PD_Stop	! Stopping Time 
	REAL(KIND=8)		:: PD_Elap	! Elapsed Time 

	!----------------------
	!	Initialize MPI	
	!----------------------
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
	
	!----------------------
	!	Init and BCast Global 
	! Seed, and ALFG Regs
	!----------------------
	IF (myid == 0) THEN
		CALL gseed
	END IF

	CALL MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL rint(myid)

	!----------------------
	!	Initialize Inputs	
	!----------------------
	CALL cpu_time(PD_Start)

	! read keys.in 
	CALL input_reader

	CALL ang_prob_data_read

	CALL CPU_TIME(PD_Stop)
	PD_Elap = PD_Stop - PD_Start 
	IF (myid == 0) THEN
		WRITE(*,"(A,F5.2,A)") "Inputs read in ", PD_Elap, " SEC"
		WRITE(*,*)
	END IF

	!----------------------
	!	Call different tests	
	!----------------------

	CALL call_tests 

	!----------------------
	!	Finalize MPI	
	!----------------------
	CALL MPI_FINALIZE(ierr)

END PROGRAM



