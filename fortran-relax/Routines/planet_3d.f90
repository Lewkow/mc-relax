!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D planetary MC relaxation experiment
!!
!! Complete, autonomous routine to calculate
!! planetary relaxation of ENAs including
!! secondary hot atoms. Escape fraction and
!! energy distributions are calculated
!! in a 3D planetary enviornment. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE planet_3d
	USE planet
	USE secondary 
	USE rand_seed
	USE tables
	USE physics_constants
	USE click_data
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER				:: i, j, k, My_N_MC, tot, MC, N_Coll, Transport, GET_TCS
	INTEGER				:: N_r_zone, N_Engy_Dist, N_SHA_Engy_Dist
	INTEGER				:: Root_Escape_Count, Root_Therm_Count, Root_Second_Count

	REAL(KIND=8),PARAMETER :: GB_per_real = 8.0D-9
	REAL(KIND=8),PARAMETER :: MB_per_real = 8.0D-6
	
	REAL(KIND=8)	:: x_0, y_0, z_0, V_0, V, Length
	REAL(KIND=8)	:: x1, y1, z1, x2, y2, z2, dL
	REAL(KIND=8) 	:: E_now, z_now, pos_now, vz_now, vx_now
	REAL(KIND=8) 	:: E_nxt, z_nxt, pos_nxt, vz_nxt, vx_nxt
	REAL(KIND=8)	:: Den_now, theta_prev, phi_prev, ScattAng, TCS	
	REAL(KIND=8)	:: theta_now, phi_now, Lab_ScattAng, E_temp
	REAL(KIND=8)	:: phi, MFP, r, C1, E_targ, Theta_targ, Phi_targ	
	REAL(KIND=8)	:: Second_Ave, Therm_Ave, Planet_Ave, Escape_Ave 
	REAL(KIND=8)	:: HeH_TCS, HeO_TCS, HeHe_TCS, HCO2_TCS
	REAL(KIND=8)	:: Den_Tot, Den_He, Den_O, Den_H, Den_CO2, rand_targ
	REAL(KIND=8)	:: MFP_HeH, MFP_HCO2, MFP_HeO, MFP_HeHe
	REAL(KIND=8)	:: P_H, P_He, P_CO2, P_O, P_tot
	REAL(KIND=8)	:: MR_H, MR_He, MR_CO2, MR_O, h, dh
	REAL(KIND=8)	:: R1_L, R1_U, R2_L, R2_U, R3_L, R3_U
	REAL(KIND=8)	:: Max_r_zone, dr_zone, Engy_Dist_In, Engy_Dist_Fn, d_Engy_Dist
	REAL(KIND=8)	:: Engy_Dist_SHA_In, Engy_Dist_SHA_Fn, d_Engy_Dist_SHA
	REAL(KIND=8)	:: lfg

	REAL(KIND=8),DIMENSION(4) :: Prob_Targ
  REAL(KIND=8),DIMENSION(3) :: x_now, x_nxt
  REAL(KIND=8),DIMENSION(3) :: u_now, u_nxt, u_targ

	!! Read in planetary inputs
	CALL read_planet_inputs

	!! Read Mars electron density tables
!	CALL read_e_den_table

	!! Read Hp+H elastic TCS tables
!	CALL read_HpH_tcs_table

	!! Read in ENA starting location tables
	CALL mars_ena_table_read

	!! Allocate click data array
	CALL allocate_click

	!! Write run info to file if root
	IF ( myid .EQ. 0 ) CALL write_run_info

	!! Set number of height zones for arrays
	N_r_zone   = 100
	Max_r_zone = high ! [m]
	dr_zone		 = Max_r_zone/REAL(N_r_zone-1)

	!! Allocate space for all ranks for height zone arrays
	ALLOCATE( r_zone_E(N_r_zone), r_zone_C(N_r_zone), r_zone(N_r_zone) )
	ALLOCATE( r_zone_dE(N_r_zone), r_zone_dC(N_r_zone) )
	ALLOCATE( r_zone_SHA_E(N_r_zone), r_zone_SHA_C(N_r_zone) )
	ALLOCATE( r_zone_SHA_dE(N_r_zone), r_zone_SHA_dC(N_r_zone) )
	ALLOCATE( r_prod_SHA_E(N_r_zone), r_prod_SHA_C(N_r_zone) )

	!! Allocate space for root zone arrays
	IF (myid .EQ. 0) THEN
		ALLOCATE( Root_r_zone_E(N_r_zone), Root_r_zone_C(N_r_zone), Root_Term_Ave(np,3) )
		ALLOCATE( Root_r_zone_dE(N_r_zone), Root_r_zone_dC(N_r_zone) )
		ALLOCATE( Root_r_zone_SHA_E(N_r_zone), Root_r_zone_SHA_C(N_r_zone) )
		ALLOCATE( Root_r_zone_SHA_dE(N_r_zone), Root_r_zone_SHA_dC(N_r_zone) )
		ALLOCATE( Root_r_prod_SHA_E(N_r_zone), Root_r_prod_SHA_C(N_r_zone) )
	END IF

	!! Set number of energy zones for arrays
	N_Engy_Dist  			= N_r_zone
	Engy_Dist_In 			= 0.0D0
	Engy_Dist_SHA_In 	= 0.0D0
	Engy_Dist_Fn 			= 3.5D3
	Engy_Dist_SHA_Fn  = 1.0D2
	d_Engy_Dist  			= (Engy_Dist_Fn-Engy_Dist_In)/REAL(N_Engy_Dist-1)
	d_Engy_Dist_SHA   = (Engy_Dist_SHA_Fn-Engy_Dist_SHA_In)/REAL(N_Engy_Dist-1)

	!! Allocate space for all ranks for energy zone arrays
	ALLOCATE( My_Engy_Dist(N_r_zone,N_Engy_Dist), E_zone(N_Engy_Dist) )
	ALLOCATE( My_SHA_Engy_Dist(N_r_zone,N_Engy_Dist), SHA_E_Zone(N_Engy_Dist) )
	
	!! Allocate space for root energy dist arrays
	IF (myid .EQ. 0) ALLOCATE( Root_Engy_Dist(N_r_zone,N_Engy_Dist) )	
  IF (myid .EQ. 0) ALLOCATE( Root_SHA_Engy_Dist(N_r_zone,N_Engy_Dist) )

	!! Initialize arrays to zeros
	r_zone_E(:) 					= 0.0D0
	r_zone_C(:) 					= 0.0D0
	r_zone_dE(:)					= 0.0D0	
	r_zone_dC(:)					= 0.0D0	
	r_zone_SHA_E(:) 			= 0.0D0
	r_zone_SHA_C(:) 			= 0.0D0
	r_prod_SHA_E(:) 			= 0.0D0
	r_prod_SHA_C(:) 			= 0.0D0
	r_zone_SHA_dE(:)			= 0.0D0	
	r_zone_SHA_dC(:)			= 0.0D0	
	My_Engy_Dist(:,:) 		= 0.0D0
	My_SHA_Engy_Dist(:,:) = 0.0D0

	!! Fill zone index arrays
	DO i=1,N_r_zone
		r_zone(i) = (i-1)*dr_zone
	END DO

	DO i=1,N_Engy_Dist
		E_zone(i) 		= (i-1)*d_Engy_Dist
		SHA_E_zone(i) = (i-1)*d_Engy_Dist_SHA
	END DO	

  !! Set up number of MC particles per rank
  My_N_MC = N_Part / np

  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  !! Collect My_N_MC to root to make sure full number of particles are being computed
  CALL MPI_REDUCE( My_N_MC, tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  !! Add 1 particle to all ranks N_MC if not enough
  IF (tot .LT. N_Part) THEN
    My_N_MC = My_N_MC + 1
    CALL MPI_REDUCE( My_N_MC, tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BCAST( tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    N_Part = tot
  END IF 

	!! Write grid information to screen
	CALL write_planet_grid( N_r_zone, N_Engy_Dist, Max_r_zone, Engy_Dist_Fn, N_part, My_N_MC )

	!! Write welcome messages to screen
  IF ( myid == 0 ) THEN
    WRITE(*,'(A)') '#################################################'
    WRITE(*,'(A)') '####  STARTING 3D PLANETARY ENA SIMULATION   ####'
    WRITE(*,'(A,I5,A)') '#### ',np, ' Processor(s) Being Used           ####'
		WRITE(*,'(A,F5.2,A)') '####    Total Memory Footprint [GB]: ', &
		& (12.0D0*N_Part + np*(N_r_zone*N_Engy_Dist+N_Engy_Dist) + &
		& N_r_zone*N_Engy_Dist + np*3.0D0*N_r_zone+2.0D0*N_r_zone)*GB_per_real, '   ####'
    WRITE(*,'(A)') '#################################################'
    WRITE(*,*)
    WRITE(*,*)
  END IF

	!! Open file to write to
  IF (myid .EQ. 0) THEN
		OPEN(UNIT=222, FILE="../Data/planet_ENA_Energy_Loss.dat",   ACCESS="APPEND")
		OPEN(UNIT=223, FILE="../Data/planet_SHA_Energy_Loss.dat",   ACCESS="APPEND")
		OPEN(UNIT=333, FILE="../Data/planet_ENA_Energy_Dist.dat",   ACCESS="APPEND")
		OPEN(UNIT=444, FILE="../Data/planet_SHA_Energy_Dist.dat",   ACCESS="APPEND")
		OPEN(UNIT=334, FILE="../Data/planet_ENA_3D_Energy_Dist.dat",ACCESS="APPEND")
		OPEN(UNIT=335, FILE="../Data/planet_SHA_3D_Energy_Dist.dat",ACCESS="APPEND")
		OPEN(UNIT=555, FILE="../Data/planet_SHA_Production.dat",    ACCESS="APPEND")
  END IF
	
	!! Write density profile for atosphere to file
	dh = high/1000.0D0
	DO i=1,1000
		h = i*dh
		CALL kras_mars_density( 44, h, Den_CO2 )
 		CALL kras_mars_density( 1,  h, Den_H   )
 		CALL kras_mars_density( 16, h, Den_O   )
 		CALL kras_mars_density( 4,  h, Den_He  )
	END DO

	!! Init counters
	Planet_Energy = 0.0D0
	Escape_Energy = 0.0D0
	Second_Energy = 0.0D0
	Ave_Ang				= 0.0D0
	CO2_Ang				= 0.0D0
	O_Ang					= 0.0D0
	H_Ang					= 0.0D0
	He_Ang				= 0.0D0
	Therm_Count   = 0
	Planet_Count  = 0
	Escape_Count  = 0
	Second_Count  = 0
	H_Count				= 0
	He_Count			= 0
	O_Count				= 0
	CO2_Count			= 0
	Coll_Count		= 0

	!! Allocate memory for SHA arrays
	CALL allocate_sha

	!! Fill SW velocity arrays
	CALL read_Vsw_table

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Main DO loop for initially hot MC particles
	DO MC=1,My_N_MC

		!! Set SHA counter for current MC particle to 0
		My_SHA_count = 0

		!! Initialize SHA arrays to zero
		CALL clean_sha

		!! Set initial position and height of particle [m]
		CALL mars_ena_start( x_0 )
		IF (x_0 .GE. high) x_0 = high
		y_0 = 0.0D0
		z_0 = 0.0D0

		!! Set energy [eV], position and velocity of hot MC particle
		IF (PROJ .EQ. 'H  ') CALL rand_init_energy( 1.0D0, E_now )
		IF (PROJ .EQ. 'He ') CALL rand_init_energy( 4.0D0, E_now )
		IF (E_now .GT. Engy_Dist_Fn) E_now = Engy_Dist_Fn
!		E_now = 1.0D3

		x_now(1) 		= x_0
		x_now(2) 		= y_0
		x_now(3) 		= z_0
		u_now(1) 		= -COS(SZA)
		u_now(2) 		= SIN(SZA) 
		u_now(3) 		= 0.0D0
		theta_prev 	= ACOS(u_now(3))
		phi_prev   	= ACOS(u_now(1)/SIN(theta_prev))

		CALL planet_transport( E_now, x_now(1), x_now(2), x_now(3), &
		&                      u_now(1), u_now(2), u_now(3), 1 )

!		WRITE(*,*) 'ENA: ', MC, ' complete with ', My_SHA_Count, ' SHAs'

		IF ( DO_SHA .EQ. 1 ) THEN
			!! Loop over SHAs for this particular ENA
			DO i=1,My_SHA_Count
				E_now    		= SHA_E(i)
				x_now(1) 		= SHA_R(i,1)
				x_now(2) 		= SHA_R(i,2)
				x_now(3) 		= SHA_R(i,3)
				u_now(1) 		= SHA_V(i,1)
				u_now(2) 		= SHA_V(i,2)
				u_now(3) 		= SHA_V(i,3)
				theta_prev 	= ACOS(u_now(3))
				phi_prev    = ACOS(u_now(1)/SIN(theta_prev))
				CALL planet_transport( E_now, x_now(1), x_now(2), x_now(3), &
				&                      u_now(1), u_now(2), u_now(3), 2 )
			END DO
		END IF

		IF ( MOD(100.0*REAL(MC)/REAL(My_N_MC),10.0) .EQ. 0 ) THEN 
			WRITE(*,'(A,I4,A,F5.1,A)') "Processor ", myid, " is ", &
			&  100.0*REAL(MC)/REAL(My_N_MC), " % Complete"
		END IF

	END DO ! MC

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr ) 
	!@@@@@@@@@@@@

	!! Free memory for SHAs
	CALL free_sha

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr ) 
	!@@@@@@@@@@@@

	!! Individual ranks averages
	Planet_Ave = Planet_Energy/Planet_Count
	Escape_Ave = Escape_Energy/Escape_Count
	Second_Ave = Second_Energy/Second_Count

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr ) 
	!@@@@@@@@@@@@

	!! Collect termination averages to root
	CALL MPI_REDUCE( Second_Count, Root_Second_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_GATHER( Planet_Ave, 1, MPI_DOUBLE_PRECISION, Root_Term_Ave(:,1), 1, &
	& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_GATHER( Escape_Ave, 1, MPI_DOUBLE_PRECISION, Root_Term_Ave(:,2), 1, &
	& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_GATHER( Second_Ave, 1, MPI_DOUBLE_PRECISION, Root_Term_Ave(:,3), 1, &
	& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	!@@@@@@@@@@@@

	IF (myid .EQ. 0) THEN
		Planet_Ave = 0.0D0
		Escape_Ave = 0.0D0
		Second_Ave = 0.0D0
		DO i=1,np
			Planet_Ave = Planet_Ave + Root_Term_Ave(i,1)
			Escape_Ave = Escape_Ave + Root_Term_Ave(i,2)
			Second_Ave = Second_Ave + Root_Term_Ave(i,3)
		END DO ! np
		Planet_Ave = Planet_Ave/REAL(np)
		Escape_Ave = Escape_Ave/REAL(np)
		Second_Ave = Second_Ave/REAL(np)
	END IF ! root

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	!@@@@@@@@@@@@

	!! Collect all thermalization, planet, and secondary counts
	CALL MPI_Reduce( Therm_Count, Root_Therm_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_Reduce( Escape_Count, Root_Escape_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_Reduce( Second_Count, Root_Second_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	!@@@@@@@@@@@@

	!! Reduce Sum zone parameters to root
	CALL MPI_REDUCE( r_zone_E, root_r_zone_E, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_zone_C, root_r_zone_C, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_zone_SHA_E, root_r_zone_SHA_E, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_zone_SHA_C, root_r_zone_SHA_C, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_prod_SHA_E, root_r_prod_SHA_E, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_prod_SHA_C, root_r_prod_SHA_C, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_zone_dE, root_r_zone_dE, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_zone_dC, root_r_zone_dC, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_zone_SHA_dE, root_r_zone_SHA_dE, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( r_zone_SHA_dC, root_r_zone_SHA_dC, N_r_zone, MPI_DOUBLE_PRECISION, &
	& 							 MPI_SUM, 0, MPI_COMM_WORLD, ierr )


	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	!@@@@@@@@@@@@

	!! Reduce sum energy/height zone array
	CALL MPI_REDUCE( My_Engy_Dist, Root_Engy_Dist, N_Engy_Dist*N_r_zone, &
	& 							 MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( My_SHA_Engy_Dist, Root_SHA_Engy_Dist, N_Engy_Dist*N_r_zone, &
	& 							 MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	!@@@@@@@@@@@@

	!! Reduce sum average scattering angles
	CALL MPI_REDUCE( Ave_Ang, ROOT_Ave_Ang, 1, MPI_DOUBLE_PRECISION, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_REDUCE( CO2_Ang, ROOT_CO2_Ang, 1, MPI_DOUBLE_PRECISION, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_REDUCE( O_Ang, ROOT_O_Ang, 1, MPI_DOUBLE_PRECISION, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_REDUCE( He_Ang, ROOT_He_Ang, 1, MPI_DOUBLE_PRECISION, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_REDUCE( H_Ang, ROOT_H_Ang, 1, MPI_DOUBLE_PRECISION, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			
	!@@@@@@@@@@@@

	!! Reduce sum average scattering angle counters
	CALL MPI_REDUCE( Coll_Count, ROOT_Coll_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_REDUCE( CO2_Count, ROOT_CO2_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_REDUCE( O_Count, ROOT_O_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_REDUCE( He_Count, ROOT_He_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_REDUCE( H_Count, ROOT_H_Count, 1, MPI_INTEGER, &
	& MPI_SUM, 0, MPI_COMM_WORLD, ierr )	

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			
	!@@@@@@@@@@@@

	!! Write zone dependent average energies to file
	IF (myid .EQ. 0) THEN

		DO i=1,N_r_zone
			IF (root_r_zone_C(i) .EQ. 0) THEN
				WRITE(333,*) r_zone(i), 0.0D0, 0.0D0
			ELSE
				WRITE(333,*) r_zone(i), root_r_zone_E(i), root_r_zone_C(i)
			END IF
		END DO

		DO i=1,N_r_zone
			IF (root_r_zone_SHA_C(i) .EQ. 0) THEN
				WRITE(444,*) r_zone(i), 0.0D0, 0.0D0
			ELSE
				WRITE(444,*) r_zone(i), root_r_zone_SHA_E(i), root_r_zone_SHA_C(i)
			END IF
		END DO

		DO i=1,N_r_zone
			IF (root_r_prod_SHA_C(i) .EQ. 0) THEN
				WRITE(555,*) r_zone(i), 0.0D0, 0.0D0
			ELSE
				WRITE(555,*) r_zone(i), root_r_prod_SHA_E(i), root_r_prod_SHA_C(i)
			END IF
		END DO


		DO i=1,N_r_zone
			IF (root_r_zone_dC(i) .EQ. 0) THEN
				WRITE(222,*) r_zone(i), 0.0D0, 0.0D0
			ELSE
				WRITE(222,*) r_zone(i), root_r_zone_dE(i), root_r_zone_dC(i)
			END IF
		END DO

		DO i=1,N_r_zone
			IF (root_r_zone_SHA_dC(i) .EQ. 0) THEN
				WRITE(223,*) r_zone(i), 0.0D0, 0.0D0
			ELSE
				WRITE(223,*) r_zone(i), root_r_zone_SHA_dE(i), root_r_zone_SHA_dC(i)
			END IF
		END DO

		DO i=1,N_r_zone
			DO j=1,N_Engy_Dist
				WRITE(334,*) r_zone(i), E_zone(j), Root_Engy_Dist(i,j)
			END DO
		END DO

		DO i=1,N_r_zone
			DO j=1,N_Engy_Dist
				WRITE(335,*) r_zone(i), SHA_E_zone(j), Root_SHA_Engy_Dist(i,j)
			END DO
		END DO

	END IF

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			
	!@@@@@@@@@@@@

	43 FORMAT (I8,A,F6.1)
	44 FORMAT (F5.1,A,F6.1)
	55 FORMAT (F5.1,A)
	56 FORMAT (F7.3,A,F4.1,A)
	77 FORMAT (A,F4.1,A)

	IF (myid .EQ. 0) THEN
		WRITE(*,*) 
		WRITE(*,'(A)') "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
		WRITE(*,55) 100.0*Root_Therm_Count/REAL(N_Part+Root_Second_Count), &
		&						" % MC particles thermalized "

		IF (Planet_Count .GT. 0) THEN
			WRITE(*,44) 100.0*Planet_Count/REAL(N_Part), &
			& " % MC particles hit the planet with average energy [eV]: ", Planet_Ave
		ELSE
			WRITE(*,44) 100.0*Planet_Count/REAL(N_Part), " % MC particles hit the planet" 
		END IF

		WRITE(*,44) 100.0*Root_Escape_Count/REAL(N_Part), &
		& " % MC particles escaped with average energy [eV]: ", Escape_Ave
		WRITE(*,43) INT(Root_Second_Count/REAL(N_Part)), &
		& " average SHAs per ENA with average energy [eV]: ", Second_Ave 
		WRITE(*,*)
		WRITE(*,77) "Ave lab scatt angle all collisions: ", &
		& Ave_Ang/REAL(Coll_Count), " [deg]"	
		WRITE(*,56) 100.0*REAL(CO2_Count)/REAL(Coll_Count), &
		& " % CO2 collisions, ave lab scatt angle: ", CO2_Ang/REAL(CO2_Count), " [deg]" 
		WRITE(*,56) 100.0*REAL(O_Count)/REAL(Coll_Count), &
		& " % O collisions,   ave lab scatt angle: ", O_Ang/REAL(O_Count), " [deg]" 
		WRITE(*,56) 100.0*REAL(He_Count)/REAL(Coll_Count), &
		& " % He collisions,  ave lab scatt angle: ", He_Ang/REAL(He_Count), " [deg]" 
		WRITE(*,56) 100.0*REAL(H_Count)/REAL(Coll_Count), &
		& " % H collisions,   ave lab scatt angle: ", H_Ang/REAL(H_Count), " [deg]" 
		WRITE(*,'(A)') "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
		WRITE(*,*) 
	END IF ! root

	!! clean up memory for all ranks
	DEALLOCATE( r_zone_E, r_zone_C, r_zone, My_Engy_Dist, E_zone )
	DEALLOCATE( r_zone_dE, r_zone_dC )
	DEALLOCATE( r_zone_SHA_E, r_zone_SHA_C )
	DEALLOCATE( r_prod_SHA_E, r_prod_SHA_C )
	DEALLOCATE( r_zone_SHA_dE, r_zone_SHA_dC )

!	CALL clean_e_den_table
!	CALL clean_HpH_tcs_table
	CALL clean_click


	!! clean up root memory
	IF (myid .EQ. 0) THEN
		DEALLOCATE( Root_r_zone_E, Root_r_zone_C, Root_Engy_Dist ) 
		DEALLOCATE( Root_r_zone_dE, Root_r_zone_dC ) 
		DEALLOCATE( Root_r_zone_SHA_E, Root_r_zone_SHA_C )
		DEALLOCATE( Root_r_zone_SHA_dE, Root_r_zone_SHA_dC )
		DEALLOCATE( Root_r_prod_SHA_E, Root_r_prod_SHA_C )
		DEALLOCATE( Root_Term_Ave )
	END IF

	!! close file handles if root
	IF (myid .EQ. 0) THEN
		CLOSE(222)	
		CLOSE(223)	
		CLOSE(333)	
		CLOSE(444)	
		CLOSE(334)	
		CLOSE(335)	
		CLOSE(555)
	END IF


END SUBROUTINE planet_3d

