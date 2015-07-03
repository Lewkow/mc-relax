
MODULE planet 
	USE physics_constants, ONLY : PI

	CHARACTER(LEN=3)					:: Proj			! Projectile type [H, He3, He4, O]
	CHARACTER(LEN=3)					:: Targ			! Target type [H, He3, He4, O]
	CHARACTER(LEN=2)					:: CM_TYPE	! Collision Method Type [QM, HS]

	REAL(KIND=8)							:: E_0			! Projectile Energy [eV]
	REAL(KIND=8)							:: E_Therm	! Thermal Energy [eV]
	REAL(KIND=8)							:: E_Esc		! Escape Energy [eV]
	REAL(KIND=8)							:: MP 			! Projectile Mass [amu]
	REAL(KIND=8)							:: MT 			! Target Mass [amu]
	REAL(KIND=8)							:: MU 			! Reduced Mass [amu]
	REAL(KIND=8)							:: h_0			! Initial height of projectile [m]
	REAL(KIND=8)							:: SZA			! Solar Zeneith Angle [deg]
	REAL(KIND=8)							:: high  		! height at which particle has escaped [m]

	REAL(KIND=8)							:: mass
	REAL(KIND=8)							:: Escape_Energy
	REAL(KIND=8)							:: Planet_Energy
	REAL(KIND=8)							:: Second_Energy

	REAL(KIND=8)							:: CO2_Ang
	REAL(KIND=8)							:: O_Ang
	REAL(KIND=8)							:: He_Ang
	REAL(KIND=8)							:: H_Ang
	REAL(KIND=8)							:: Ave_Ang

	REAL(KIND=8)							:: ROOT_Ave_Ang
	REAL(KIND=8)							:: ROOT_CO2_Ang
	REAL(KIND=8)							:: ROOT_O_Ang
	REAL(KIND=8)							:: ROOT_He_Ang
	REAL(KIND=8)							:: ROOT_H_Ang
	REAL(KIND=8)							:: INIT_ENGY	! if mono-energetic, starting energy

	INTEGER										:: ATMOSPHERE ! (0 = MIN SA) (1 = Mean SA) (2 = Max SA)
	INTEGER										:: DO_SHA			! DO SHA tracking = 1 	
	INTEGER										:: DO_GRAV		! Gravity transport = 1 	
	INTEGER										:: MONO_ENGY	! Mono-energetic = 1 prob dist = 2	
	INTEGER										:: N_Part 		! Number of particles to run until thermalized
	INTEGER										:: N_Max			! Maximum number of collisions for a particle
	INTEGER										:: atom				! scaling amplitude atom-atom or atom-molecule
	INTEGER										:: Therm_Count, Planet_Count, Escape_Count, Second_Count
	INTEGER										:: H_Count, He_Count, O_Count, CO2_Count, Coll_Count
	INTEGER										:: ROOT_H_Count, ROOT_He_Count, ROOT_O_Count
	INTEGER										:: ROOT_CO2_Count, ROOT_Coll_Count

	INTEGER										:: Click_Skip
	INTEGER										:: Max_Click

	INTEGER										:: N_SH_Z0

	INTEGER										:: WRITE_E_T
	INTEGER										:: WRITE_H_T
	INTEGER										:: WRITE_N_T
	INTEGER										:: WRITE_Ux_T
	INTEGER										:: WRITE_H_Ux
	INTEGER										:: WRITE_H_N
	INTEGER										:: WRITE_H_E
	INTEGER										:: WRITE_H_dE
	INTEGER										:: WRITE_ALL_E_T
	INTEGER										:: WRITE_ALL_H_T
	INTEGER										:: WRITE_ALL_N_T
	INTEGER										:: WRITE_ALL_Ux_T
	INTEGER										:: WRITE_ALL_H_Ux
	INTEGER										:: WRITE_ALL_H_N
	INTEGER										:: WRITE_ALL_H_E
	INTEGER										:: WRITE_ALL_H_dE
	INTEGER										:: WRITE_SHA_H_E
	INTEGER										:: WRITE_SHA_H_Ux

	REAL(KIND=8),PARAMETER		:: HS_R_H 	= 1.20D-10
	REAL(KIND=8),PARAMETER		:: HS_R_He 	= 1.40D-10
	REAL(KIND=8),PARAMETER		:: HS_R_Ar 	= 1.88D-10
	REAL(KIND=8),PARAMETER		:: HS_R_O 	= 1.52D-10
	REAL(KIND=8),PARAMETER		:: HS_R_H2 	= 2.40D-10
	REAL(KIND=8),PARAMETER		:: HS_R_N2 	= 3.10D-10
	REAL(KIND=8),PARAMETER		:: HS_R_CO 	= 3.22D-10
	REAL(KIND=8),PARAMETER		:: HS_R_CO2 = 4.74D-10

	REAL(KIND=8),PARAMETER		:: M_H   = 1.00782503207D0
  REAL(KIND=8),PARAMETER		:: M_H2  = 2.0157D0
  REAL(KIND=8),PARAMETER		:: M_N2  = 28.0134D0
  REAL(KIND=8),PARAMETER		:: M_Ar  = 39.948D0
  REAL(KIND=8),PARAMETER		:: M_He  = 4.00260325415D0
  REAL(KIND=8),PARAMETER		:: M_O   = 15.99491461956D0
	REAL(KIND=8),PARAMETER		:: M_CO2 = 44.00964D0
	REAL(KIND=8),PARAMETER		:: M_CO	 = 28.0147D0

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: H_ENA_PROB		! from tables for starting ENA position
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: H_ENA_DIST
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: He_ENA_PROB
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: He_ENA_DIST

  REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: Nas_ENA_Energy, Nas_ENA_Height

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: r_zone_dE, r_zone_dC
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: r_zone_E, r_zone_C, r_zone
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_r_zone_E, Root_r_zone_C
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_r_zone_dE, Root_r_zone_dC
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: E_zone, SHA_E_zone
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: SH_Z0_Z, SH_Z0_P
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: Root_Term_Ave
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: My_Engy_Dist, Root_Engy_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: My_SHA_Engy_Dist, Root_SHA_Engy_Dist

	INTEGER																	:: N_Atmos
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Z_Atmos
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: H_Atmos, He_Atmos, O_Atmos, Ar_Atmos
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: H2_Atmos, N2_Atmos, CO_Atmos, CO2_Atmos

	INTEGER, PARAMETER											:: ENA_SH = 1 ! 1 = H ENA   4 = He ENA

	CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE mass_finder
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Finds the mass of the target and 
  ! projectile based on the char inputs
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!! Projectile
  IF ( Proj == 'H  ' ) THEN
    MP = M_H 
  ELSE IF ( Proj == 'H2 ' ) THEN
    MP = M_H2 
  ELSE IF ( Proj == 'N2 ' ) THEN
    MP = M_N2
  ELSE IF ( Proj == 'Ar ' ) THEN
    MP = M_Ar
  ELSE IF ( Proj == 'He ' ) THEN
    MP = M_He 
  ELSE IF ( Proj == 'O  ' ) THEN
    MP = M_O 
	ELSE IF ( Proj == 'CO2' ) THEN
    MP = M_CO2 
	ELSE IF ( Proj == 'CO ' ) THEN
    MP = M_CO 
  ELSE
    WRITE(*,*) "Projectile ", Proj, " not known!!!!"
  END IF

  !!! Target
  IF ( TRIM(Targ) == 'H' ) THEN
		MT = M_H
  ELSE IF ( TRIM(Targ) == 'H2' ) THEN
		MT = M_H2
  ELSE IF ( TRIM(Targ) == 'N2' ) THEN
		MT = M_N2
  ELSE IF ( TRIM(Targ) == 'Ar' ) THEN
		MT = M_Ar
  ELSE IF ( TRIM(Targ) == 'He' ) THEN
		MT = M_He
  ELSE IF ( TRIM(Targ) == 'O' ) THEN
		MT = M_O
	ELSE IF ( TRIM(Targ) == 'CO2' ) THEN
		MT = M_CO2
	ELSE IF ( TRIM(Targ) == 'CO' ) THEN
		MT = M_CO
  ELSE
    WRITE(*,*) "Target ", Targ, " not known!!!!"
  END IF

	MU = MP*MT/(MP+MT)

  END SUBROUTINE mass_finder

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	SUBROUTINE read_planet_inputs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Read inputs from planet_keys.in and 
	! save information to planet module
	!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		USE mpi_info
	
		IMPLICIT NONE

		INCLUDE 'mpif.h'

		INTEGER,PARAMETER     :: N_real = 5
		INTEGER,PARAMETER     :: N_int  = 2
		INTEGER,PARAMETER     :: N_char = 1
		REAL(KIND=8) 					:: dum_real(N_real)
		INTEGER								:: dum_int(N_int)
		CHARACTER(LEN=2)			:: dum_char

		IF ( myid .EQ. 0) THEN
			OPEN(21, FILE="../Inputs/planet_keys.in", STATUS="old", ACTION="read")

			!!! Comments
			READ(21,*)	
			READ(21,*)	
			READ(21,*)	
			READ(21,*)	

			READ(21,*) Proj	
			READ(21,*) CM_TYPE
			READ(21,*) ATMOSPHERE
			READ(21,*) E_Therm
			READ(21,*) E_Esc
			READ(21,*) h_0
			READ(21,*) high
			READ(21,*) SZA
			READ(21,*) N_Part
			READ(21,*) DO_SHA	
			READ(21,*) DO_GRAV
			READ(21,*) MONO_ENGY
			READ(21,*) INIT_ENGY
			READ(21,*) Max_click
			READ(21,*) Click_Skip

			READ(21,*)	
			READ(21,*)	
			READ(21,*)	

			READ(21,*) WRITE_E_T	
			READ(21,*) WRITE_H_T	
			READ(21,*) WRITE_N_T	
			READ(21,*) WRITE_Ux_T	
			READ(21,*) WRITE_H_Ux	
			READ(21,*) WRITE_H_N
			READ(21,*) WRITE_H_E
			READ(21,*) WRITE_H_dE

			READ(21,*) WRITE_ALL_E_T	
			READ(21,*) WRITE_ALL_H_T	
			READ(21,*) WRITE_ALL_N_T	
			READ(21,*) WRITE_ALL_Ux_T	
			READ(21,*) WRITE_ALL_H_Ux	
			READ(21,*) WRITE_ALL_H_N
			READ(21,*) WRITE_ALL_H_E
			READ(21,*) WRITE_ALL_H_dE
			
			READ(21,*) WRITE_SHA_H_E
			READ(21,*) WRITE_SHA_H_Ux
	
			65 FORMAT (A,F4.1,A)
			66 FORMAT (A,ES8.2,A)
			67 FORMAT (A,I7,A)
			68 FORMAT (2A)
			69 FORMAT (A)
			70 FORMAT (3A)

!		Targ = Proj
!		CALL mass_finder

			!!! Write inputs to screen
			WRITE(*,*)
			WRITE(*,69) '#################################################'
			WRITE(*,69) '####     PLANETARY INPUTS (planet_keys.in)   ####'
			WRITE(*,69) '#################################################'
			WRITE(*,67) '####  Number of MC Particles   : ', N_Part,'     ####'
			WRITE(*,70) '####  Projectile ENA           : ', PROJ, '         ####'
			IF (CM_TYPE .EQ. 'QM') THEN
				WRITE(*,69) '####  Collision Method         : Quantum     ####'
			END IF
			IF (CM_TYPE .EQ. 'HS') THEN
				WRITE(*,69) '####  Collision Method         : Hard Sphere ####'
			END IF
			IF (ATMOSPHERE .EQ. 0) THEN
				WRITE(*,69) '####  Atmosphere Model         : Min SA      ####'
			ELSE IF (ATMOSPHERE .EQ. 1) THEN
				WRITE(*,69) '####  Atmosphere Model         : Mean SA     ####'
			ELSE IF (ATMOSPHERE .EQ. 2) THEN
				WRITE(*,69) '####  Atmosphere Model         : Max SA      ####'
			END IF
			WRITE(*,66) '####  Escape Height [km]       : ', high/1.0D3, '    ####'
			WRITE(*,66) '####  Escape Energy [eV]       : ', E_Esc, '    ####'
			WRITE(*,66) '####  Thermal Energy [eV]      : ', E_Therm, '    ####'
			WRITE(*,65) '####  Solar Zenith Angle [deg] : ', SZA, '        ####'
			IF (MONO_ENGY .EQ. 1) THEN
				WRITE(*,69) '####  Mono-Energetic Simulation ON           ####'
				WRITE(*,66) '####  Starting Energy [eV]     : ', INIT_ENGY, '    ####'
			ELSE
				WRITE(*,69) '####  Solar Wind Probability Simulation ON   ####'
			END IF
			IF (DO_SHA .EQ. 1) THEN
				WRITE(*,69) '####  Secondary Hot Atom Tracking ON         ####'  
			ELSE
				WRITE(*,69) '####  Secondary Hot Atom Tracking OFF        ####'  
			END IF
			IF (DO_GRAV .EQ. 1) THEN
				WRITE(*,69) '####  Planet Gravity ON                      ####'  
			ELSE
				WRITE(*,69) '####  Planet Gravity OFF                     ####'  
			END IF
			WRITE(*,67) '####  Maximum Number of Clicks : ', Max_click,'     ####'
			WRITE(*,67) '####  Data Write Click Skipped : ', click_skip,'     ####'
			WRITE(*,69) '#################################################'
		END IF ! root

		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

		CALL MPI_BCAST( Proj, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( CM_TYPE, 2, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( E_Therm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( E_Esc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( h_0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( high, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( SZA, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( ATMOSPHERE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( N_Part, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( DO_SHA, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( DO_GRAV, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( MONO_ENGY, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( INIT_ENGY, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( Click_Skip, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( Max_Click, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_E_T, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_H_T, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_N_T, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_Ux_T, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_H_Ux, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_H_N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_H_E, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_H_dE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_ALL_E_T, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_ALL_H_T, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_ALL_N_T, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_ALL_Ux_T, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_ALL_H_Ux, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_ALL_H_N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_ALL_H_E, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_ALL_H_dE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_SHA_H_E, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( WRITE_SHA_H_Ux, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

		!! convert SZA from deg to rad
		SZA = SZA*PI/180.0D0

		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	END SUBROUTINE read_planet_inputs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	SUBROUTINE write_planet_grid( N_r_zone, N_Engy_Dist, Max_r_zone, Engy_Dist_Fn, N, My_N )
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	
	! Write grid information for planet simulation	
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Inputs
	INTEGER				:: N_r_zone, N_Engy_Dist, N, My_N
	REAL(KIND=8)	:: Max_r_zone, Engy_Dist_Fn		! [m], [eV]

    65 FORMAT (A,F4.1,A)
    66 FORMAT (A,ES8.2,A)
    67 FORMAT (A,I7,A)
    68 FORMAT (2A)
   	69 FORMAT (A)
  	70 FORMAT (3A)

	IF (myid .EQ. 0) THEN
		WRITE(*,69) '#################################################'
    WRITE(*,69) '####        PLANETARY GRID INFORMATION       ####'
		WRITE(*,69) '#################################################'
		WRITE(*,67) '####  Number of height zones   : ', N_r_zone, '     ####'
		WRITE(*,67) '####  Number of energy zones   : ', N_Engy_Dist, '     ####'
		WRITE(*,66) '####  Maximum height zone [km] : ', Max_r_zone/1.0D3, '    ####'
		WRITE(*,66) '####  Maximum energy zone [eV] : ', Engy_Dist_Fn, '    ####'
		WRITE(*,69) '#################################################'
		WRITE(*,69) '#################################################'
    WRITE(*,67) '####  Total MC particles       : ', N, '     ####'
    WRITE(*,67) '####  MC particles per rank    : ', My_N, '     ####'
		WRITE(*,69) '#################################################'
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	END SUBROUTINE write_planet_grid

!#############################################
!#############################################

	SUBROUTINE read_mars_density
		USE mpi_info

		IMPLICIT NONE

		INCLUDE 'mpif.h'
	
		INTEGER		:: i

		IF (myid .EQ. 0) THEN
			IF (ATMOSPHERE .EQ. 0) THEN
				OPEN(31, FILE="../Tables/Krasnopolsky_Mars_Min_Density.dat", STATUS="old", ACTION="read")
			ELSE IF (ATMOSPHERE .EQ. 1) THEN
				OPEN(31, FILE="../Tables/Krasnopolsky_Mars_Mean_Density.dat", STATUS="old", ACTION="read")
			ELSE IF (ATMOSPHERE .EQ. 2) THEN
				OPEN(31, FILE="../Tables/Krasnopolsky_Mars_Max_Density.dat", STATUS="old", ACTION="read")
			END IF	
			READ(31,*) N_Atmos
			ALLOCATE(Z_Atmos(N_Atmos))
			ALLOCATE(H_Atmos(N_Atmos))
			ALLOCATE(H2_Atmos(N_Atmos))
			ALLOCATE(He_Atmos(N_Atmos))
			ALLOCATE(O_Atmos(N_Atmos))
			ALLOCATE(Ar_Atmos(N_Atmos))
			ALLOCATE(N2_Atmos(N_Atmos))
			ALLOCATE(CO_Atmos(N_Atmos))
			ALLOCATE(CO2_Atmos(N_Atmos))
			DO i=1,N_Atmos
				READ(31,*) Z_Atmos(i), H_Atmos(i), H2_Atmos(i), He_Atmos(i),  &
				& O_Atmos(i), Ar_Atmos(i), N2_Atmos(i), CO_Atmos(i), CO2_Atmos(i)
			END DO	
			CLOSE(31)
		END IF ! myid = 0	

		CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
		CALL MPI_BCAST(N_Atmos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

		IF (myid .NE. 0) THEN
			ALLOCATE(Z_Atmos(N_Atmos))
			ALLOCATE(H_Atmos(N_Atmos))
			ALLOCATE(H2_Atmos(N_Atmos))
			ALLOCATE(He_Atmos(N_Atmos))
			ALLOCATE(O_Atmos(N_Atmos))
			ALLOCATE(Ar_Atmos(N_Atmos))
			ALLOCATE(N2_Atmos(N_Atmos))
			ALLOCATE(CO_Atmos(N_Atmos))
			ALLOCATE(CO2_Atmos(N_Atmos))
		END IF
			
		CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
		CALL MPI_BCAST(Z_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(H_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(H2_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(He_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(O_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(Ar_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(N2_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(CO_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(CO2_Atmos,N_Atmos,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

	END SUBROUTINE read_mars_density

!#############################################
!#############################################

	SUBROUTINE clean_mars_density

			DEALLOCATE(Z_Atmos)
			DEALLOCATE(H_Atmos)
			DEALLOCATE(H2_Atmos)
			DEALLOCATE(He_Atmos)
			DEALLOCATE(O_Atmos)
			DEALLOCATE(Ar_Atmos)
			DEALLOCATE(N2_Atmos)
			DEALLOCATE(CO_Atmos)
			DEALLOCATE(CO2_Atmos)

	END SUBROUTINE clean_mars_density

!#############################################
!#############################################

	SUBROUTINE mars_table_density(targ, h, den)
		
		IMPLICIT NONE

		!! Inputs 
		CHARACTER(LEN=3)  :: targ			! target atom/molecule
  	REAL(KIND=8)      :: h        ! [m] 

  	!! Outputs
  	REAL(KIND=8)      :: den      ! [1/m^3]

		!! Internal
		REAL(KIND=8)			:: h_km, x1, x2, y1, y2, m, y0, n
		REAL(KIND=8)			:: f_den(N_Atmos)
		INTEGER						:: t1, t2, i, j, NOW

		h_km = h/1.0D3

		IF (TRIM(targ) .EQ. 'H') THEN
			f_den = H_Atmos
		ELSE IF (TRIM(targ) .EQ. 'He') THEN
			f_den = He_Atmos
		ELSE IF (TRIM(targ) .EQ. 'O') THEN
			f_den = O_Atmos
		ELSE IF (TRIM(targ) .EQ. 'Ar') THEN
			f_den = Ar_Atmos
		ELSE IF (TRIM(targ) .EQ. 'H2') THEN
			f_den = H2_Atmos
		ELSE IF (TRIM(targ) .EQ. 'N2') THEN
			f_den = N2_Atmos
		ELSE IF (TRIM(targ) .EQ. 'CO') THEN
			f_den = CO_Atmos
		ELSE IF (TRIM(targ) .EQ. 'CO2') THEN
			f_den = CO2_Atmos
		ELSE
			WRITE(*,*) 'Atmosphere species ', targ, ' not in density file'
		END IF

		IF ( (h_km .GE. Z_Atmos(1)) .AND. (h_km .LE. Z_Atmos(N_Atmos)) ) THEN
			! inside file altitudes
			NOW = 0
			j   = 1
			DO WHILE (NOW .EQ. 0)	
				j = j + 1
				IF ( (h_km .LT. Z_Atmos(j)) .OR. (j .EQ. N_Atmos) ) NOW = 1
			END DO
			t1 = j-1
			t2 = j
			x1 = Z_Atmos(t1)
			x2 = Z_Atmos(t2)
			y1 = LOG(f_den(t1))	
			y2 = LOG(f_den(t2))	
			m  = (y2-y1)/(x2-x1)
			y0 = y2 - m*x2
			n  = EXP(m*h_km + y0)
		ELSE ! outside file altitudes
			IF ( h_km .LT. Z_Atmos(1) ) THEN
				t1 = 1
				t2 = 2
				x1 = Z_Atmos(t1)	
				x2 = Z_Atmos(t2)	
				y1 = LOG(f_den(t1))
				y2 = LOG(f_den(t2))
				m  = (y2-y1)/(x2-x1)
				y0 = y2 - m*x2
				n  = EXP(m*h_km + y0)
			ELSE
				t1 = N_Atmos-1
				t2 = N_Atmos
				x1 = Z_Atmos(t1)	
				x2 = Z_Atmos(t2)	
				y1 = LOG(f_den(t1))
				y2 = LOG(f_den(t2))
				m  = (y2-y1)/(x2-x1)
				y0 = y2 - m*x2
				n  = EXP(m*h_km + y0)
			END IF
		END IF

		den = n*1e6			! convert from 1/cm^3 -> 1/m^3

	END SUBROUTINE mars_table_density

	!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE test_mars_table_density
	
		IMPLICIT NONE

		REAL(KIND=8)		:: zi, zf, dz, z
		REAL(KIND=8)		:: H_den, He_den, O_den, Ar_den, H2_den, N2_den, CO_den, CO2_den
		INTEGER					:: N, i
	
!		CALL read_mars_density

		N  = 1000
		zi = 50.0D3
		zf = 800.0D3
		dz = (zf-zi)/REAL(N-1)

		OPEN(UNIT=61,FILE='../Data/Mars_Density_Test.dat',STATUS='NEW',ACCESS='APPEND')

		DO i=1,N
			z = zi + REAL(i-1)*dz
			CALL mars_table_density( 'H  ', z, H_den   )
			CALL mars_table_density( 'He ', z, He_den  )
			CALL mars_table_density( 'O  ', z, O_den   )
			CALL mars_table_density( 'Ar ', z, Ar_den  )
			CALL mars_table_density( 'H2 ', z, H2_den  )
			CALL mars_table_density( 'N2 ', z, N2_den  )
			CALL mars_table_density( 'CO ', z, CO_den  )
			CALL mars_table_density( 'CO2', z, CO2_den )
			WRITE(61,*) z/1000.0D0, H_den, He_den, O_den, Ar_den, H2_den, N2_den, CO_den, CO2_den
		END DO
	
		CLOSE(61)	

!		CALL clean_mars_density		


	END SUBROUTINE test_mars_table_density
	
END MODULE planet



