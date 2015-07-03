!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Module to be used for LISM calculations 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE lism

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: H_ENA_PROB
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: H_ENA_DIST
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: He_ENA_PROB
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: He_ENA_DIST

	REAL(KIND=8),PARAMETER									:: H_COMP  = 0.92D0					! percentage of SW that is H+
	REAL(KIND=8),PARAMETER									:: He_COMP = 1.0D0-H_COMP		! percentage of SW that is He++

  INTEGER		    :: N_Part 
	INTEGER				:: ENA_TYPE	   ! Type of ENA 0 => H | 1 => He
	INTEGER				:: ION_METH    ! Type of ion energy loss
														   ! 0 => none | 1 => Butler | 2 => Bethe
	INTEGER				:: ENG_METH    ! 0 => mono | 1 => Prob Den
	REAL(KIND=8)	:: MONO_ENG    ! Mono Energy [eV]
	INTEGER				:: N_r_zone    ! Number of radial zones
	INTEGER				:: N_Engy_Dist ! Number of energy zones
	REAL(KIND=8)	:: Max_r_zone  ! Maximum radial zone [AU]
	REAL(KIND=8)	:: Engy_Dist_Fn! Maximum energy [eV]
	REAL(KIND=8)	:: dr_zone     ! differential for zone distance [AU]
	REAL(KIND=8)	:: d_Engy_Dist ! differential for energy [eV]
	INTEGER				:: DO_RND_In_EN! Do random initial energy test
	INTEGER				:: DO_RND_AN   ! Do random scattering angle test
	INTEGER				:: QM_ON			 ! Use QM cross sections if = 1
															 ! Use HS cross sections if = 0

CONTAINS

SUBROUTINE lism_read_inputs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read input file ../LISM_keys.in and 
! save parameters to module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE physics_constants

	IMPLICIT NONE

  INCLUDE 'mpif.h'

	INTEGER,PARAMETER	:: N_int  = 7
	INTEGER,PARAMETER	:: N_real = 3
	INTEGER 					:: N_head, i
	INTEGER 					:: D_I(N_int)
	REAL(KIND=8)			:: D_R(N_real)

	IF (myid .EQ. 0) THEN

		OPEN(UNIT=33,FILE='../Inputs/LISM_keys.in',STATUS='OLD')
		N_head = 5
		DO i=1,N_head
			READ(33,*)
		END DO

		READ(33,*) ENA_TYPE 					! Int
		READ(33,*) 
		READ(33,*) N_Part
		READ(33,*) 
		READ(33,*) ION_METH						! Int
		READ(33,*) 
		READ(33,*) ENG_METH						! Int
		READ(33,*) MONO_ENG						! Real
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) N_r_zone						! Int
		READ(33,*) Max_r_zone					! Real
		READ(33,*) 
		READ(33,*) N_Engy_Dist				! Int
		READ(33,*) Engy_Dist_Fn				! Real
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) DO_RND_In_EN				! Int
		READ(33,*) 
		READ(33,*) DO_RND_AN					! Int
		READ(33,*) 

		CLOSE(33)	
	
		CALL lism_inputs_write

		D_I(1) = ENA_TYPE
		D_I(2) = ION_METH	
		D_I(3) = ENG_METH
		D_I(4) = N_R_zone
		D_I(5) = N_Engy_Dist
		D_I(6) = DO_RND_In_EN
		D_I(7) = DO_RND_AN

		D_R(1) = MONO_ENG
		D_R(2) = Max_r_zone
		D_R(3) = Engy_Dist_Fn

	END IF

	!! Broadcast the input arrays to rest of ranks
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( D_I, N_int, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( D_R, N_real, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( N_Part, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	ENA_TYPE 			= D_I(1)
	ION_METH 			= D_I(2)
	ENG_METH 			= D_I(3)
	N_R_zone 			= D_I(4)
	N_ENGY_Dist 	= D_I(5)
	DO_RND_In_EN	= D_I(6)
	DO_RND_AN			= D_I(7)

	MONO_ENG      = D_R(1)
	Max_r_zone    = D_R(2)*6.3D4 ! convert from LY to AU
	Engy_Dist_Fn  = D_R(3)

	dr_zone    		= REAL(Max_r_zone)/REAL(N_r_zone-1)
	d_Engy_Dist   = REAL(Engy_Dist_Fn)/REAL(N_Engy_Dist-1)

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE lism_read_inputs

SUBROUTINE lism_inputs_write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Write parameters from LISM_keys.in
! to screen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info	
	
	IMPLICIT NONE

	IF ( myid .EQ. 0 ) THEN

	WRITE(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'	
	WRITE(*,*) '$$$$$        LISM INPUT PARAMETERS         $$$$$'	
	IF (ENA_TYPE .EQ. 0) THEN
		WRITE(*,*) '$$$$$        ENA TYPE: H                   $$$$$'
	ELSE 
		WRITE(*,*) '$$$$$        ENA TYPE: He                  $$$$$'
	END IF
	IF (ION_METH .EQ. 0) THEN
		WRITE(*,*) '$$$$$        ION LOSS: NONE                $$$$$'
	ELSE IF (ION_METH .EQ. 1) THEN
		WRITE(*,*) '$$$$$        ION LOSS: BUTLER              $$$$$'
	ELSE IF (ION_METH .EQ. 2) THEN
		WRITE(*,*) '$$$$$        ION LOSS: BETHE               $$$$$'
	END IF	
	IF (ENG_METH .EQ. 0) THEN
		WRITE(*,'(A,ES10.2,A)') ' $$$$$        MONO ENERGY: ', MONO_ENG, '       $$$$$'
	ELSE 
		WRITE(*,*) '$$$$$        INIT ENERGY: PROB DIST        $$$$$'
	END IF
	WRITE(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'	

	END IF

END SUBROUTINE lism_inputs_write

SUBROUTINE lism_grid_write
	USE mpi_info
	USE physics_constants

	IMPLICIT NONE

	IF ( myid .EQ. 0 ) THEN
		!! print grid info to screen
    WRITE(*,*)
    WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    WRITE(*,*) '%%%%%              GRID INFO              %%%%%'
    WRITE(*,'(A,ES10.2,A)') ' %%%%%       R_end [LY]: ', Max_r_zone*AUtoLY, '        %%%%%'
    WRITE(*,'(A,ES10.2,A)') ' %%%%%       dR    [LY]: ', dr_zone*AUtoLY, '        %%%%%'
    WRITE(*,'(A,ES10.2,A)') ' %%%%%       E_end [eV]: ', Engy_Dist_Fn, '        %%%%%'
    WRITE(*,'(A,ES10.2,A)') ' %%%%%       dE    [eV]: ', d_Engy_Dist, '        %%%%%'
    WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    WRITE(*,*)
	END IF

END SUBROUTINE lism_grid_write

END MODULE lism

