
MODULE click_data
	USE planet, ONLY : N_Part, MP, h_0, Click_Skip, Max_Click

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: CLICK
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_E_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_H_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_N_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_Ux_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_H_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_H_N
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_dE_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_H_dE
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_H_Ux
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_Ux_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_Ux_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_Ux_O
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_E_O
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_Ux_Ar
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_E_Ar
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_Ux_H2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_E_H2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_Ux_N2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_E_N2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_Ux_CO
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_E_CO
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_Ux_CO2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: C_SHA_H_E_CO2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_E_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_H_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_N_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_Ux_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_H_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_H_N
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_dE_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_H_dE
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)				  :: A_H_Ux
	
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)						:: Esc_Engy_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)						:: Esc_NColl_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)						:: Therm_Time_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)						:: Therm_Height_Dist

	INTEGER																					:: My_Num_Esc

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)						:: ROOT_Esc_Engy_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)						:: ROOT_Esc_NColl_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)						:: ROOT_Therm_Time_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)						:: ROOT_Therm_Height_Dist
	INTEGER																					:: Root_Num_Esc
	INTEGER,ALLOCATABLE,DIMENSION(:)								:: ALL_Num_Esc

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_E_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_H_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_N_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_Ux_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_H_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_H_N
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_dE_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_H_dE
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_C_H_Ux
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_Ux_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_Ux_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_Ux_O
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_E_O
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_Ux_Ar
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_E_Ar
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_Ux_H2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_E_H2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_Ux_N2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_E_N2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_Ux_CO
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_E_CO
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_Ux_CO2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_C_SHA_H_E_CO2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_E_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_H_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_N_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_Ux_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_H_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_H_N
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_dE_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_H_dE
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_A_H_Ux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CLICK holds all data to be written to file
! CLICK( click, histogram_bin, z )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z = 1   time [sec]
! z = 2 	energy [eV]
! z = 3   x position [m]
! z = 4   y position [m]
! z = 5   z position (height) [m]
! z = 6   ux
! z = 7   uy
! z = 8   uz
! z = 9   theta [rad]
! z = 10  phi [rad]
! z = 11  total number of collisions for MC particle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: CLICK_0
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: CLICK_1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CLICK_0 hold data in same way as CLICK but for current click only
! CLICK_1 hold data in same way as CLICK but for next click only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: X_CLICK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X_CLICK hols the reference vectors for histgram data
! X_CLICK( histogram_bin, z )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z = 1  time [sec]
! z = 2  energy [eV]
! z = 3  x position [m]
! z = 4  y position [m]
! z = 5  z position (height) [m]
! z = 6  ux unit velocity
! z = 7  uy unit velocity
! z = 8  uz unit velocity
! z = 9  theta range
! z = 10 phi range
! z = 11 Num Collisions
! z = 12 energy loss [eV]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER,PARAMETER																:: Num_Hist  = 50
	INTEGER																					:: Num_Tot_Clicks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SUBROUTINE allocate_click
	USE planet
	USE mpi_info

	IMPLICIT NONE

  REAL(KIND=8),PARAMETER :: GB_per_real = 8.0D-9

	!! num total clicks plus initial condition of 1 click
	Num_Tot_Clicks = Max_Click/Click_Skip

	ALLOCATE( C_E_T( Num_Tot_Clicks,  Num_Hist, Num_Hist ) )
	ALLOCATE( C_H_T( Num_Tot_Clicks,  Num_Hist, Num_Hist ) )
	ALLOCATE( C_N_T( Num_Tot_Clicks,  Num_Hist, Num_Hist ) )
	ALLOCATE( C_Ux_T( Num_Tot_Clicks, Num_Hist, Num_Hist ) )
	ALLOCATE( C_H_E( Num_Tot_Clicks,  Num_Hist, Num_Hist ) )
	ALLOCATE( C_H_N( Num_Tot_Clicks,  Num_Hist, Num_Hist ) )
	ALLOCATE( C_X_Y( Num_Tot_Clicks,  Num_Hist, Num_Hist ) )
	ALLOCATE( C_dE_T( Num_Tot_Clicks, Num_Hist, Num_Hist ) )
	ALLOCATE( C_H_dE( Num_Tot_Clicks, Num_Hist, Num_Hist ) )
	ALLOCATE( C_H_Ux( Num_Tot_Clicks, Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_Ux_H( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_E_H( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_Ux_He( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_E_He( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_Ux_O( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_E_O( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_Ux_Ar( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_E_Ar( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_Ux_H2( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_E_H2( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_Ux_N2( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_E_N2( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_Ux_CO( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_E_CO( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_Ux_CO2( Num_Hist, Num_Hist ) )
	ALLOCATE( C_SHA_H_E_CO2( Num_Hist, Num_Hist ) )
	ALLOCATE( CLICK_0( N_Part,   11 ) )
	ALLOCATE( CLICK_1( N_Part,   11 ) )
	ALLOCATE( X_CLICK( Num_Hist, 12 ) )
	ALLOCATE( Esc_Engy_Dist( N_Part ) )
	ALLOCATE( Esc_NColl_Dist( N_Part ) )
	ALLOCATE( Therm_Time_Dist( N_Part ) )
	ALLOCATE( Therm_Height_Dist( N_Part ) )

	C_SHA_H_E_H       = 0.0D0	
	C_SHA_H_Ux_H      = 0.0D0	
	C_SHA_H_E_He      = 0.0D0	
	C_SHA_H_Ux_He     = 0.0D0	
	C_SHA_H_E_O       = 0.0D0	
	C_SHA_H_Ux_O      = 0.0D0	
	C_SHA_H_E_Ar      = 0.0D0	
	C_SHA_H_Ux_Ar     = 0.0D0	
	C_SHA_H_E_H2      = 0.0D0	
	C_SHA_H_Ux_H2     = 0.0D0	
	C_SHA_H_E_N2      = 0.0D0	
	C_SHA_H_Ux_N2     = 0.0D0	
	C_SHA_H_E_CO      = 0.0D0	
	C_SHA_H_Ux_CO     = 0.0D0	
	C_SHA_H_E_CO2     = 0.0D0	
	C_SHA_H_Ux_CO2    = 0.0D0	
	Esc_Engy_Dist 		= 0.0D0
	Esc_NColl_Dist 		= 0.0D0
	Therm_Time_Dist 	= 0.0D0
	Therm_Height_Dist = 0.0D0
	My_Num_Esc 				= 0

	IF (WRITE_ALL_E_T  .EQ. 1) THEN
		ALLOCATE(A_E_T(Num_Hist,Num_Hist))
		A_E_T = 0.0D0
	END IF
	IF (WRITE_ALL_H_T  .EQ. 1) THEN
		ALLOCATE(A_H_T(Num_Hist,Num_Hist))
		A_H_T = 0.0D0
	END IF
	IF (WRITE_ALL_N_T  .EQ. 1) THEN
		ALLOCATE(A_N_T(Num_Hist,Num_Hist))
		A_N_T = 0.0D0
	END IF
	IF (WRITE_ALL_Ux_T .EQ. 1) THEN
		ALLOCATE(A_Ux_T(Num_Hist,Num_Hist))
		A_Ux_T = 0.0D0
	END IF
	IF (WRITE_ALL_H_Ux .EQ. 1) THEN
		ALLOCATE(A_H_Ux(Num_Hist,Num_Hist))
		A_H_Ux = 0.0D0
	END IF
	IF (WRITE_ALL_H_N  .EQ. 1) THEN
		ALLOCATE(A_H_N(Num_Hist,Num_Hist))
		A_H_N = 0.0D0
	END IF
	IF (WRITE_ALL_H_E  .EQ. 1) THEN
		ALLOCATE(A_H_E(Num_Hist,Num_Hist))
		A_H_E = 0.0D0
	END IF
	IF (WRITE_ALL_H_dE .EQ. 1) THEN
		ALLOCATE(A_H_dE(Num_Hist,Num_Hist))
		A_H_dE = 0.0D0
	END IF

	IF (myid .EQ. 0) THEN
		WRITE(*,'(A,I7,A)') '####  ', Num_Tot_Clicks, ' Write Clicks Written           ####'
		WRITE(*,'(A,F5.2,A)') '####  ',np*10.0D0*Num_Tot_Clicks*Num_Hist**2*GB_per_real, ' GB used to allocate Click arrays ####'
	END IF

	C_E_T  			= 0.0D0
	C_H_T  			= 0.0D0
	C_N_T  			= 0.0D0
	C_Ux_T 			= 0.0D0
	C_H_E  			= 0.0D0
	C_H_N  			= 0.0D0
	C_X_Y  			= 0.0D0
	C_dE_T 			= 0.0D0
	C_H_dE 			= 0.0D0
	C_H_Ux 			= 0.0D0

END SUBROUTINE allocate_click

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SUBROUTINE clean_click
	USE planet

	DEALLOCATE( C_E_T )
	DEALLOCATE( C_H_T )
	DEALLOCATE( C_N_T )
	DEALLOCATE( C_Ux_T )
	DEALLOCATE( C_H_E )
	DEALLOCATE( C_H_N )
	DEALLOCATE( C_X_Y )
	DEALLOCATE( C_dE_T )
	DEALLOCATE( C_H_dE )
	DEALLOCATE( C_H_Ux )
	DEALLOCATE( C_SHA_H_Ux_H )
	DEALLOCATE( C_SHA_H_E_H )
	DEALLOCATE( C_SHA_H_Ux_He )
	DEALLOCATE( C_SHA_H_E_He )
	DEALLOCATE( C_SHA_H_Ux_O )
	DEALLOCATE( C_SHA_H_E_O )
	DEALLOCATE( C_SHA_H_Ux_Ar )
	DEALLOCATE( C_SHA_H_E_Ar )
	DEALLOCATE( C_SHA_H_Ux_H2 )
	DEALLOCATE( C_SHA_H_E_H2 )
	DEALLOCATE( C_SHA_H_Ux_N2 )
	DEALLOCATE( C_SHA_H_E_N2 )
	DEALLOCATE( C_SHA_H_Ux_CO )
	DEALLOCATE( C_SHA_H_E_CO )
	DEALLOCATE( C_SHA_H_Ux_CO2 )
	DEALLOCATE( C_SHA_H_E_CO2 )
	DEALLOCATE( CLICK_0 )
	DEALLOCATE( CLICK_1 )
	DEALLOCATE( X_CLICK )
	DEALLOCATE( Esc_Engy_Dist )
	DEALLOCATE( Esc_NColl_Dist )
	DEALLOCATE( Therm_Time_Dist )
	DEALLOCATE( Therm_Height_Dist )

	IF (WRITE_ALL_E_T  .EQ. 1) DEALLOCATE(A_E_T)
	IF (WRITE_ALL_H_T  .EQ. 1) DEALLOCATE(A_H_T)
	IF (WRITE_ALL_N_T  .EQ. 1) DEALLOCATE(A_N_T)
	IF (WRITE_ALL_Ux_T .EQ. 1) DEALLOCATE(A_Ux_T)
	IF (WRITE_ALL_H_Ux .EQ. 1) DEALLOCATE(A_H_Ux)
	IF (WRITE_ALL_H_N  .EQ. 1) DEALLOCATE(A_H_N)
	IF (WRITE_ALL_H_E  .EQ. 1) DEALLOCATE(A_H_E)
	IF (WRITE_ALL_H_dE .EQ. 1) DEALLOCATE(A_H_dE)

END SUBROUTINE clean_click

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! set initial value for click data

SUBROUTINE set_init_click_data_Vsw( h, SZA )
	USE physics_constants, ONLY : PI
	USE planet, ONLY : PROJ
	USE ena
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Input
	REAL(KIND=8)		:: h		!! Initial height [m]
	REAL(KIND=8)		:: SZA	!! Solar Zenith Angle [rad]

	!! Internal
	INTEGER					:: i
	REAL(KIND=8)		:: v		!! velocity [m/s]
	REAL(KIND=8)		:: lfg	!! random number generator
	REAL(KIND=8)		:: theta, phi, r_phi, E0, MP, rh

	IF (myid .EQ. 0) THEN
		IF (PROJ .EQ. 'H ') THEN
			MP = 1.0D0
		ELSE IF (PROJ .EQ. 'He') THEN
			MP = 4.0D0
		ELSE
			WRITE(*,*) 'Proj: ', PROJ, ' not known!!'
		END IF
		OPEN(UNIT=65,FILE='../Data/Initial_Energy.dat',ACCESS='APPEND')
		DO i=1,N_Part
			CALL rand_init_energy( MP, E0 )
			CALL rand_ENA_altitude( rh )
			WRITE(65,*) E0
			r_phi 					= lfg()
			phi 						= r_phi*2.0D0*PI
			theta           = PI - SZA
			CLICK_0( i, 1 ) = 0.0D0								! Time [sec]
			CLICK_0( i, 2 ) = E0									! Energy [eV] 		
			CLICK_0( i, 3 ) = 0.0D0								! x0		
			CLICK_0( i, 4 ) = 0.0D0     					! y0
			CLICK_0( i, 5 ) = rh									! z0
			CLICK_0( i, 6 ) = SIN(theta)*COS(phi) ! ux0
			CLICK_0( i, 7 ) = SIN(theta)*SIN(phi) ! uy0
			CLICK_0( i, 8 ) = COS(theta) 					! uz0
			CLICK_0( i, 9 ) = 0.0D0 							! theta scattering angle [rad]
			CLICK_0( i, 10) = 0.0D0								! phi scattering angle [rad]
			CLICK_0( i, 11) = 0.0D0								! number of collisions
		END DO
		CLOSE(65)
	END IF
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(CLICK_0,N_Part*11,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

END SUBROUTINE set_init_click_data_Vsw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! set initial value for click data for
! Nascent SH simulation

SUBROUTINE set_init_click_data_SH
	USE physics_constants, ONLY : PI
	USE planet, ONLY : PROJ
	USE ena
	USE escape_trans
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Internal
	INTEGER					:: i
	REAL(KIND=8)		:: lfg	!! random number generator
	REAL(KIND=8)		:: theta, phi, r_theta, r_phi, E0, rh

	IF (myid .EQ. 0) THEN
		OPEN(UNIT=65,FILE='../Data/Nascent_SH_Initial_Energy.dat',ACCESS='APPEND')
		OPEN(UNIT=66,FILE='../Data/Nascent_SH_Initial_Altitude.dat',ACCESS='APPEND')
		DO i=1,N_Part
			CALL rand_Nascent_SH_Energy( E0 )
			CALL mars_sh_start( rh )
			WRITE(65,*) E0
			WRITE(66,*) rh
			r_phi 					= lfg()
			r_theta 				= lfg()
			phi 						= r_phi*2.0D0*PI
			theta           = r_theta*PI 
			CLICK_0( i, 1 ) = 0.0D0								! Time [sec]
			CLICK_0( i, 2 ) = E0									! Energy [eV] 		
			CLICK_0( i, 3 ) = 0.0D0								! x0		
			CLICK_0( i, 4 ) = 0.0D0     					! y0
			CLICK_0( i, 5 ) = rh									! z0
			CLICK_0( i, 6 ) = SIN(theta)*COS(phi) ! ux0
			CLICK_0( i, 7 ) = SIN(theta)*SIN(phi) ! uy0
			CLICK_0( i, 8 ) = COS(theta) 					! uz0
			CLICK_0( i, 9 ) = 0.0D0 							! theta scattering angle [rad]
			CLICK_0( i, 10) = 0.0D0								! phi scattering angle [rad]
			CLICK_0( i, 11) = 0.0D0								! number of collisions
		END DO
		CLOSE(65)
		CLOSE(66)
	END IF
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(CLICK_0,N_Part*11,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

END SUBROUTINE set_init_click_data_SH

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! set initial value for click data

SUBROUTINE set_init_click_data( E, h, SZA )
	USE physics_constants, ONLY : PI
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Input
	REAL(KIND=8)		:: E		!! Initial energy [eV]
	REAL(KIND=8)		:: h		!! Initial height [m]
	REAL(KIND=8)		:: SZA	!! Solar Zenith Angle [rad]

	!! Internal
	INTEGER					:: i
	REAL(KIND=8)		:: v		!! velocity [m/s]
	REAL(KIND=8)		:: lfg	!! random number generator
	REAL(KIND=8)		:: theta, phi, r_phi

	IF (myid .EQ. 0) THEN
		DO i=1,N_Part
			r_phi 					= lfg()
			phi 						= r_phi*2.0D0*PI
			theta           = PI - SZA
			CLICK_0( i, 1 ) = 0.0D0								! Time [sec]
			CLICK_0( i, 2 ) = E										! Energy [eV] 		
			CLICK_0( i, 3 ) = 0.0D0								! x0		
			CLICK_0( i, 4 ) = 0.0D0     					! y0
			CLICK_0( i, 5 ) = h										! z0
			CLICK_0( i, 6 ) = SIN(theta)*COS(phi) ! ux0
			CLICK_0( i, 7 ) = SIN(theta)*SIN(phi) ! uy0
			CLICK_0( i, 8 ) = COS(theta) 					! uz0
			CLICK_0( i, 9 ) = 0.0D0 							! theta scattering angle [rad]
			CLICK_0( i, 10) = 0.0D0								! phi scattering angle [rad]
			CLICK_0( i, 11) = 0.0D0								! number of collisions
		END DO
	END IF
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(CLICK_0,N_Part*11,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

END SUBROUTINE set_init_click_data

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! set reference vectors for click data

SUBROUTINE click_ref_vectors
	USE physics_constants, ONLY : PI
	USE planet, ONLY : MONO_ENGY, INIT_ENGY, PROJ, high, E_therm
	USE inputs, ONLY : DO_Escape_MC

	IMPLICIT NONE

	INTEGER				:: i
	REAL(KIND=8)	:: TEnd, EEnd, REnd, NEnd, low, dEEnd, highZ

	highZ = 250.0D3
	low  	= 75.0D3
	TEnd 	= 0.75D0	
	REnd 	= 100.0D3
	NEnd 	= Max_Click/4
	dEEnd = 5.0D0

	IF (DO_Escape_MC .EQ. 1) THEN
		EEnd = 5.0D0
	ELSE
		IF (MONO_ENGY .EQ. 1) THEN
			EEnd = INIT_ENGY
		ELSE
			IF (PROJ .EQ. 'H ') THEN
				EEnd = 3.0D3
			ELSE IF (PROJ .EQ. 'He') THEN
				EEnd = 4.0D0*3.0D3
			ELSE
				WRITE(*,*) 'PROJ: ', PROJ, ' not known! '
			END IF
		END IF
	END IF

	DO i=1,Num_Hist
		X_CLICK(i,1)  = 0.0D0 + (i-1)*(TEnd/REAL(Num_Hist-1))				! time [sec]
		X_CLICK(i,2)  = 1.0D0 + (i-1)*(EEnd/REAL(Num_Hist-1))				! energy [eV]
!		X_CLICK(i,2)  = E_therm + (i-1)*(EEnd/REAL(Num_Hist-1))				! energy [eV]
		X_CLICK(i,3)  = -REnd + (i-1)*(2.0D0*REND/REAL(Num_Hist-1))	! x position [m]
		X_CLICK(i,4)  = -REnd + (i-1)*(2.0D0*REND/REAL(Num_Hist-1))	! y position [m]
		X_CLICK(i,5)  = low   + (i-1)*((highZ-low)/REAL(Num_Hist-1))	! z position (height) [m]
		X_CLICK(i,6)  = -1.0D0+ (i-1)*(2.0D0/REAL(Num_Hist-1))			! ux unit velocity
		X_CLICK(i,7)  = -1.0D0+ (i-1)*(2.0D0/REAL(Num_Hist-1))			! uy unit velocity
		X_CLICK(i,8)  = -1.0D0+ (i-1)*(2.0D0/REAL(Num_Hist-1))			! uz unit velocity
		X_CLICK(i,9)  = 0.0D0 + (i-1)*(PI/REAL(Num_Hist-1))					! theta range
		X_CLICK(i,10) = 0.0D0 + (i-1)*(2.0D0*PI/REAL(Num_Hist-1))		! phi range
		X_CLICK(i,11) = 0.0D0 + (i-1)*(NEnd/REAL(Num_Hist-1))				! number of collisions
		X_CLICK(i,12) = 0.0D0 + (i-1)*(dEEnd/REAL(Num_Hist-1))			! energy loss/SHA energy [eV]
	END DO

END SUBROUTINE click_ref_vectors

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! write the reference vectors to files
! in data directory
SUBROUTINE write_click_ref_vectors

	IMPLICIT NONE

	INTEGER		:: i

	OPEN(UNIT=90, FILE="../Data/planet_X_time.dat", ACCESS="APPEND")
	OPEN(UNIT=91, FILE="../Data/planet_X_energy.dat", ACCESS="APPEND")
	OPEN(UNIT=92, FILE="../Data/planet_X_height.dat", ACCESS="APPEND")
	OPEN(UNIT=93, FILE="../Data/planet_X_xy.dat", ACCESS="APPEND")
	OPEN(UNIT=94, FILE="../Data/planet_X_unit_vel.dat", ACCESS="APPEND")
	OPEN(UNIT=95, FILE="../Data/planet_X_theta.dat", ACCESS="APPEND")
	OPEN(UNIT=96, FILE="../Data/planet_X_phi.dat", ACCESS="APPEND")
	OPEN(UNIT=97, FILE="../Data/planet_X_energy_loss.dat", ACCESS="APPEND")
	OPEN(UNIT=98, FILE="../Data/planet_X_Ncoll.dat", ACCESS="APPEND")

	DO i=1,Num_Hist
		WRITE(90,*) X_CLICK(i,1)
		WRITE(91,*) X_CLICK(i,2)
		WRITE(92,*) X_CLICK(i,5)
		WRITE(93,*) X_CLICK(i,3)
		WRITE(94,*) X_CLICK(i,6)
		WRITE(95,*) X_CLICK(i,9)
		WRITE(96,*) X_CLICK(i,10)
		WRITE(97,*) X_CLICK(i,12)
		WRITE(98,*) X_CLICK(i,11)
	END DO

	CLOSE(90)
	CLOSE(91)
	CLOSE(92)
	CLOSE(93)
	CLOSE(94)
	CLOSE(95)
	CLOSE(96)
	CLOSE(97)
	CLOSE(98)

END SUBROUTINE write_click_ref_vectors

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



END MODULE click_data

