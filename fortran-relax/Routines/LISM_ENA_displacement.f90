
SUBROUTINE LISM_ENA_displacement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Experiment to calculate relaxation of ENAs in the LISM which
! orginate throught solar wind charge exchange collisions. 
! Made to be 'main' routine for experiment, in that most 
! computation takes place in other subroutines for ease of 
! addition of higher complexity to the problem in the future. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants
	USE lism
	USE ENA, ONLY : L_ion, L_targ
	USE flux_mapping
	USE planet, ONLY : M_H, M_He
	USE mpi_info	
	
	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Internal
	INTEGER																	:: N_MC, My_N_MC, MC
	INTEGER																	:: q, zz, z, i, j, tot 
	INTEGER																	:: N_hist
	INTEGER																	:: max_collisions 
	INTEGER																	:: N_Coll_zone
	INTEGER																	:: My_10_per_done
	INTEGER,DIMENSION(np) 									:: buff_tot

	CHARACTER(LEN=3)												:: Proj

	REAL,PARAMETER													:: MB_per_real = 8.0D-6
	REAL,PARAMETER													:: GB_per_real = 8.0D-9

	REAL(KIND=8)														:: lfg
	REAL(KIND=8)														:: N_Coll, N_CX_Coll_H, N_Atom_Coll_H, N_Coll_H
	REAL(KIND=8)														:: N_CX_Coll_He, N_Atom_Coll_He, N_Coll_He
	REAL(KIND=8)														:: root_N_CX, root_N_AT, root_N
	REAL(KIND=8)														:: ScattAng, Lab_ScattAng, randy
	REAL(KIND=8)														:: CX_CS, TCS, MFP, Den_now, ion_Den_now
	REAL(KIND=8)														:: E_now, E_nxt, E_temp
	REAL(KIND=8)														:: r1, r2, temp_t, temp_p
	REAL(KIND=8)														:: theta_now, phi_now, theta_prev, phi_prev, phi
	REAL(KIND=8)														:: max_diss, tot_per 
	REAL(KIND=8)														:: Start_t, End_t, Full_Start, Full_End
	REAL(KIND=8)														:: Elapsed, Full_min, Full_sec
	REAL(KIND=8)														:: r_now, r
	REAL(KIND=8)														:: Den_tot, MR_H, MR_Hp, force, ion_path
	REAL(KIND=8)														:: Engy_Dist_In
	REAL(KIND=8)														:: Prob_atom, Prob_ion, Prob_tot, tt, tau

	REAL(KIND=8),DIMENSION(3)								:: u_now, u_nxt
	REAL(KIND=8),DIMENSION(3)								:: x_now, x_nxt
	REAL(KIND=8),DIMENSION(4)								:: XXX

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: My_Data
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: Root_Data
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: My_Engy_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: Root_Engy_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_CX_C_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_CX_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_EL_C_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_EL_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_CX_C_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_CX_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_EL_C_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_EL_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Diss
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Hist_x
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Hist_y 
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: r_zone_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: r_zone_C
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: r_zone
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: E_zone
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_r_zone_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_r_zone_C
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_CX_C_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_CX_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_EL_C_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_EL_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_CX_C_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_CX_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_EL_C_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_EL_E_He

	QM_ON = 1

	!! get inptus from file LISM_keys.in
	CALL lism_read_inputs

	!! set up initial energy tables
	CALL read_Vsw_table

	!! Run tests if input file dictates such
	IF ( DO_RND_IN_EN .EQ. 1) CALL test_rand_energy
	IF ( DO_RND_AN    .EQ. 1) CALL test_lin_rand_angle

	!! set up collision number bins for energy loss and collision numbers
	N_Coll_zone = 50000
	ALLOCATE( My_N_CX_E_H(N_Coll_zone), My_N_CX_C_H(N_Coll_zone) )
	ALLOCATE( My_N_EL_E_H(N_Coll_zone), My_N_EL_C_H(N_Coll_zone) )	
	ALLOCATE( My_N_CX_E_He(N_Coll_zone), My_N_CX_C_He(N_Coll_zone) ) 
	ALLOCATE( My_N_EL_E_He(N_Coll_zone), My_N_EL_C_He(N_Coll_zone) )	
	IF (myid .EQ. 0) THEN
		ALLOCATE( Root_N_CX_E_H(N_Coll_zone), Root_N_CX_C_H(N_Coll_zone) )
		ALLOCATE( Root_N_EL_E_H(N_Coll_zone), Root_N_EL_C_H(N_Coll_zone) )
		ALLOCATE( Root_N_CX_E_He(N_Coll_zone), Root_N_CX_C_He(N_Coll_zone) )
		ALLOCATE( Root_N_EL_E_He(N_Coll_zone), Root_N_EL_C_He(N_Coll_zone) )
	END IF
	My_N_CX_E_H(:)  = 0.0D0
	My_N_CX_C_H(:)  = 0.0D0
	My_N_EL_E_H(:)  = 0.0D0
	My_N_EL_C_H(:)  = 0.0D0
	My_N_CX_E_He(:) = 0.0D0
	My_N_CX_C_He(:) = 0.0D0
	My_N_EL_E_He(:) = 0.0D0
	My_N_EL_C_He(:) = 0.0D0

	ALLOCATE( r_zone_E(N_r_zone), r_zone_C(N_r_zone), r_zone(N_r_zone) )
	IF (myid .EQ. 0) ALLOCATE( Root_r_zone_E(N_r_zone), Root_r_zone_C(N_r_zone) )
	
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	r_zone_E(:) = 0.0D0			
	r_zone_C(:) = 0.0D0			

	!! write grid parameters to screen
	CALL lism_grid_write

	ALLOCATE( My_Engy_Dist(N_r_zone,N_Engy_Dist), E_zone(N_Engy_Dist) )
	IF (myid .EQ. 0) ALLOCATE( Root_Engy_Dist(N_r_zone,N_Engy_Dist) )
	My_Engy_Dist(:,:) = 0.0D0

  DO i=1,N_r_zone
    r_zone(i) = REAL(i-1)*dr_zone
  END DO

	DO i=1,N_Engy_Dist
		E_zone(i) = REAL(i-1)*d_Engy_Dist
	END DO

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Use number of MC particles as specified in keys.in input file
	N_MC = N_Part

	!! Output formats
	11 FORMAT(I7, A)
	12 FORMAT(A, I4, A, I9, A)
	13 FORMAT(A, I4, A, F5.1, A)
	
	!! Open file to write to
	IF (myid .EQ. 0) THEN
		OPEN(UNIT=458, FILE="../Data/LISM_ENA_dE_NColl_EL_H.dat", ACCESS="APPEND")
		OPEN(UNIT=459, FILE="../Data/LISM_ENA_dE_NColl_CX_H.dat", ACCESS="APPEND")
		OPEN(UNIT=558, FILE="../Data/LISM_ENA_dE_NColl_EL_He.dat", ACCESS="APPEND")
		OPEN(UNIT=559, FILE="../Data/LISM_ENA_dE_NColl_CX_He.dat", ACCESS="APPEND")
		OPEN(UNIT=567, FILE="../Data/LISM_ENA_Start_Energy_Distance.dat", ACCESS="APPEND")
		OPEN(UNIT=666, FILE="../Data/LISM_ENA_End_Pos.dat", ACCESS="APPEND")
		OPEN(UNIT=777, FILE="../Data/LISM_ENA_Diss_Hist.dat", ACCESS="APPEND")
		OPEN(UNIT=778, FILE="../Data/LISM_ENA_Diss_Percent.dat", ACCESS="APPEND")
		OPEN(UNIT=888, FILE="../Data/LISM_ENA_Coll_Hist.dat", ACCESS="APPEND")
		OPEN(UNIT=889, FILE="../Data/LISM_ENA_Coll_Percent.dat", ACCESS="APPEND")
		OPEN(UNIT=890, FILE="../Data/LISM_ENA_CX_Coll_Hist.dat", ACCESS="APPEND")
		OPEN(UNIT=891, FILE="../Data/LISM_ENA_Atom_Coll_Hist.dat", ACCESS="APPEND")
		OPEN(UNIT=333, FILE="../Data/LISM_ENA_Energy_Dist.dat", ACCESS="APPEND")
		OPEN(UNIT=555, FILE="../Data/LISM_ENA_3D_Energy_Dist.dat", ACCESS="APPEND")
	END IF

	!! Set up number of MC particles per rank
	My_N_MC 			 = N_MC / np
	My_10_per_done = My_N_MC / 10

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Collect My_N_MC to root to make sure full number of particles are being computed
	CALL MPI_REDUCE( My_N_MC, tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	
	!! Add 1 MC particle to all ranks N_MC if not enough
	IF (tot .LT. N_MC) THEN
		My_N_MC = My_N_MC + 1 
		CALL MPI_REDUCE( My_N_MC, tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		N_MC = tot
	END IF	
	
	!! Allocate Root_Data array for all MC particles
	IF (myid .EQ. 0) THEN
		ALLOCATE( Root_Data(N_MC,8) )
	END IF

	!! Allocate My_Data to hold information for each MC particle for each processor
	ALLOCATE( My_Data(My_N_MC,8) )

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Write welcome messages to screen
	IF ( myid == 0 ) THEN
		WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		WRITE(*,*) 'STARTING LISM ENA SIMULATION'
		WRITE(*,'(I3,A)') np, ' Processors Being Used'
		WRITE(*,'(ES9.2,A)') REAL(N_MC), ' Total Particles Being Simulated'
		WRITE(*,'(A,F5.2)') ' Total Storage Array Memory Footprint [GB]: ', &
		&  (12.0D0*N_MC + np*(N_r_zone*N_Engy_Dist+N_Engy_Dist) + &
		&  N_r_zone*N_Engy_Dist + np*3.0D0*N_r_zone+2.0D0*N_r_zone)*GB_per_real
		WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		WRITE(*,*)
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Write number of MC particles per processor to screen
	IF (myid .EQ. 0) THEN
		WRITE(*,'(A,ES9.2,A)') ' Computing ', REAL(My_N_MC), ' MC Particles per MPI rank' 
		WRITE(*,*)
		Full_Start = MPI_WTIME()
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! read starting location parameter tables
	CALL lism_ena_table_read

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	!! Main DO loop
	DO MC=1,My_N_MC

		!! Find random initial unit velocity vector
		CALL rand_init_velocity( theta_prev, phi_prev, u_now )

		!! Find random initial position [AU]
		CALL lism_ena_start( x_now, Proj )	
		Proj = 'H '

		!! Find random initial energy [eV] at 1 AU
		IF (Proj .EQ. 'H ') CALL rand_init_energy( 1.0D0, E_now )
		IF (Proj .EQ. 'He') CALL rand_init_energy( 4.0D0, E_now )	

		!! Initial distance from star [AU]
		r_now = SQRT( x_now(1)**2 + x_now(2)**2 + x_now(3)**2 )

		!! Get initial temp at creation location
		CALL LISM_temp( r_now, E_temp )

		IF (myid .EQ. 0) WRITE(567,*) E_now, r_now

		!! Initialize the collision counter to zero
		N_Coll					= 0.0D0	
		N_Coll_H      	= 0.0D0
		N_CX_Coll_H   	= 0.0D0
		N_Atom_Coll_H 	= 0.0D0
		N_Coll_He      	= 0.0D0
		N_CX_Coll_He   	= 0.0D0
		N_Atom_Coll_He 	= 0.0D0

		!! Start energy dependent thermalization loop
		DO WHILE ( E_now .GT. E_temp )

			!! get current distance from star
			r_now = SQRT(x_now(1)**2 + x_now(2)**2 + x_now(3)**2)

			!! Get density [1/m^3] and temperature [eV] in current position
			CALL LISM_density( r_now, Den_now )
			CALL LISM_temp( r_now, E_temp )

			!! Find atom-atom TCS for current energy
			IF (Proj .EQ. 'He ') THEN
				IF (QM_ON .EQ. 1) THEN
					IF (E_now .GE. 1.0D4) CALL TCS_HeH( 9.99D3, TCS )
					IF (E_now .LT. 1.0D4) CALL TCS_HeH( E_now, TCS )
				ELSE IF (QM_ON .EQ. 0) THEN
					CALL HS_TCS( M_He, M_H, TCS )
				END IF
				L_targ = 2
			END IF
			IF (Proj .EQ. 'H  ') THEN
				IF (QM_ON .EQ. 1) THEN
					IF (E_now .GE. 1.0D4) CALL TCS_HH( 1.0D4, TCS )
					IF (E_now .LT. 1.0D4) CALL TCS_HH( E_now, TCS )
				ELSE IF (QM_ON .EQ. 0) THEN
					CALL HS_TCS( M_H, M_He, TCS )
				END IF
				L_targ = 2
			END IF

			L_ion = 3 

			!! Convert TCS to [m^2]
			IF (QM_ON .EQ. 1) TCS = TCS*BOHRTOM**2.0D0

			!! Determine to collide with an ion or atom
			CALL LISM_ion_density( r_now, ion_Den_now ) 

			!! Get CX cross sections for proj (H or He) plus H+
			CALL cx_cross_sections( 'H_p  ', Proj, 1, E_now, CX_CS )

			!! Get mixing ratio for neutrals and ions
			Den_tot  = Den_now + ion_Den_now
			MR_H 		 = Den_now/Den_tot
			MR_Hp    = ion_Den_now/Den_tot
			Prob_tot = ion_Den_now*CX_CS + Den_now*TCS

			!! Determine if collision is atom-atom or atom-ion
			Prob_atom = Den_now*TCS/Prob_tot
			Prob_ion  = ion_Den_now*CX_CS/Prob_tot

			!! get random number to deterimine if collision is with atom or ion
			r1 = lfg()

			IF (r1 .LE. prob_atom) THEN
				!#######################################
				!## collision is elastic atom-atom #####
				!#######################################
				!! Get random free path
				CALL rand_mfp( Den_now, TCS, MFP )

				!! Convert MFP from [m] to [AU]
				MFP = MFP * MTOAU 

				!! Find random theta scattering angle
				IF (Proj .EQ. 'H ') THEN
					IF (E_now .GT. 2.5D3) THEN
						IF (QM_ON .EQ. 0) CALL HS_rand_angle( ScattAng )
						IF (QM_ON .EQ. 1) CALL lin_rand_angle( 2.5D3, 'HH   ', ScattAng )	
					ELSE
						IF (QM_ON .EQ. 0) CALL HS_rand_angle( ScattAng )
						IF (QM_ON .EQ. 1) CALL lin_rand_angle( E_now, 'HH   ', ScattAng )	
					END IF
				ELSE
					IF (E_now .GT. 4.0D3) THEN
						IF (QM_ON .EQ. 0) CALL HS_rand_angle( ScattAng )
						IF (QM_ON .EQ. 1) CALL lin_rand_angle( 4.0D3, 'HeH  ', ScattAng )
					ELSE
						IF (QM_ON .EQ. 0) CALL HS_rand_angle( ScattAng )
						IF (QM_ON .EQ. 1) CALL lin_rand_angle( E_now, 'HeH  ', ScattAng )
					END IF
				END IF

				!! Find new energy after collision at angle ScattAng
				CALL find_new_energy( E_now, ScattAng, 1.0D0, E_nxt )

				!! Convert ScattAng to lab frame
				CALL angle_to_lab( ScattAng, 1.0D0, Lab_ScattAng )

				!! Find random phi scattering angle
				r1     = lfg()
				phi    = r1*2.0D0*PI

				!! Update new scattering angles
				theta_now = Lab_ScattAng
				phi_now   = phi 

				!! Transport particle and get new unit velocity and position vectors
				CALL particle_transport( u_now, x_now, theta_prev, phi_prev, &
				&                        theta_now, phi_now, MFP, u_nxt, x_nxt ) 

				IF ( Proj .EQ. 'H ') THEN
					N_Atom_Coll_H = N_Atom_Coll_H + 1.0D0
					IF ( FLOOR(N_Atom_Coll_H) .LT. SIZE(My_N_EL_C_H) ) THEN
						My_N_EL_C_H(FLOOR(N_Atom_Coll_H)) = My_N_EL_C_H(FLOOR(N_Atom_Coll_H)) + 1.0D0				
						My_N_EL_E_H(FLOOR(N_Atom_Coll_H)) = My_N_EL_E_H(FLOOR(N_Atom_Coll_H)) + E_nxt				
					END IF
				ELSE
					N_Atom_Coll_He = N_Atom_Coll_He + 1.0D0
					IF ( FLOOR(N_Atom_Coll_He) .LT. SIZE(My_N_EL_C_He) ) THEN
						My_N_EL_C_He(FLOOR(N_Atom_Coll_He)) = My_N_EL_C_He(FLOOR(N_Atom_Coll_He)) + 1.0D0			
						My_N_EL_E_He(FLOOR(N_Atom_Coll_He)) = My_N_EL_E_He(FLOOR(N_Atom_Coll_He)) + E_nxt			
					END IF
				END IF

			ELSE
				!#######################################
				!## collision is atom-ion CX ###########
				!#######################################
				!! Get random free path of atom before CX collision with ion
				CALL rand_mfp( ion_Den_now, CX_CS, MFP )

				!! Convert MFP from [m] to [AU]
				MFP = MFP * MTOAU
	
				!! Update position location to CX position before transporting ion
				x_now(:) = x_now(:) + u_now(:)*MFP

				!! transport newly created ion 
				CALL ion_transport(1.0D0,1,E_now,CX_CS,x_now/MTOAU,u_now,x_nxt,u_nxt,tt,ion_path)

				!! convert x_nxt from [m] back to [AU]
				x_nxt = x_nxt * MTOAU

				!! find energy loss due to ion transport
				IF ( Proj .EQ. 'H ') THEN

					IF ( ION_METH .EQ. 1) THEN
						CALL butler_drag( 1, ion_path, ion_Den_now, E_temp*11604.505D0, E_now, E_nxt )
					END IF
					IF (ION_METH .EQ. 2) THEN
						CALL bethe_drag( 1, ion_path, ion_Den_now, E_now, E_nxt )
					END IF
					IF (E_nxt .LT. 0.0D0) E_nxt = 0.0D0

				ELSE IF ( Proj .EQ. 'He') THEN

					IF ( ION_METH .EQ. 1) THEN
						CALL butler_drag( 2, ion_path, ion_Den_now, E_temp*11604.505D0, E_now, E_nxt )
					END IF
					IF ( ION_METH .EQ. 2) THEN
						CALL bethe_drag( 2, ion_path, ion_Den_now, E_now, E_nxt )
					END IF
					IF (E_nxt .LT. 0.0D0) E_nxt = 0.0D0

				END IF

				IF ( Proj .EQ. 'H ') THEN
					N_CX_Coll_H = N_CX_Coll_H + 1.0D0
					IF ( FLOOR(N_CX_Coll_H) .LT. SIZE(My_N_CX_C_H) ) THEN
						My_N_CX_C_H(FLOOR(N_CX_Coll_H)) = My_N_CX_C_H(FLOOR(N_CX_Coll_H)) + 1.0D0				
						My_N_CX_E_H(FLOOR(N_CX_Coll_H)) = My_N_CX_E_H(FLOOR(N_CX_Coll_H)) + E_nxt				
					END IF
				ELSE
					N_CX_Coll_He = N_CX_Coll_He + 1.0D0
					IF ( FLOOR(N_CX_Coll_He) .LT. SIZE(My_N_CX_C_He) ) THEN
						My_N_CX_C_He(FLOOR(N_CX_Coll_He)) = My_N_CX_C_He(FLOOR(N_CX_Coll_He)) + 1.0D0				
						My_N_CX_E_He(FLOOR(N_CX_Coll_He)) = My_N_CX_E_He(FLOOR(N_CX_Coll_He)) + E_nxt				
					END IF
				END IF

			END IF

			!! Check for particle crossing mapping surface
			CALL flux_map( x_now(:), x_nxt(:), E_now )


			!! Update values for next itteration
			theta_prev = theta_now
			phi_prev   = phi_now
			x_now(:)   = x_nxt(:)
			u_now(:)   = u_nxt(:)
			E_now      = E_nxt
			N_Coll 		 = N_Coll + 1.0D0

			r = SQRT( x_now(1)**2 + x_now(2)**2 + x_now(3)**2 )
			z = 1
			DO WHILE ( (r .GE. r_zone(z)) .AND. (z .LT. N_r_zone) ) 		
				z = z + 1
			END DO
			IF (z .LE. N_r_zone) THEN
				r_zone_E(z) = r_zone_E(z) + E_now
				r_zone_C(z) = r_zone_C(z) + 1.0D0	
			END IF
	
			zz = 1
			DO WHILE ( (E_now .GE. E_zone(zz)) .AND. (zz .LT. N_Engy_Dist) )
				zz = zz + 1
			END DO
			IF ( (zz .LE. N_Engy_Dist) .AND. (z .LE. N_r_zone) .AND. (E_now .GT. E_temp) ) THEN
				My_Engy_Dist(z,zz) = My_Engy_Dist(z,zz) + 1.0D0
			END IF

			!! update local temp
			CALL LISM_temp( r, E_temp )

		END DO ! WHILE

		!! Save collision data to My_Data array
		My_Data(MC,1) = x_now(1)
		My_Data(MC,2) = x_now(2)
		My_Data(MC,3) = x_now(3) 
		My_Data(MC,4) = N_Coll
		My_Data(MC,5) = N_Atom_Coll_H
		My_Data(MC,6) = N_CX_Coll_H
		My_Data(MC,7) = N_Atom_Coll_He
		My_Data(MC,8) = N_CX_Coll_He

!		WRITE(*,*) 'Finished MC: Proj: EL: CX: TOT: ', MC, PROJ, N_Atom_Coll_H+N_CX_Coll_He, N_CX_Coll_H+N_CX_Coll_He, N_Coll

		!! Write status to screen at set intervals
		IF ( (MOD(MC,My_10_per_done) .EQ. 0) .AND. (N_Part .GE. 10) ) THEN
			WRITE(*,13) 'Rank', myid, ' is ', 100.0*REAL(MC)/REAL(My_N_MC), ' % complete'
		END IF

!		WRITE(*,*) 'RANK ',myid,' finished particle ', MC

	END DO ! MC	

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Write timing data from run to screen
	IF (myid .EQ. 0) THEN
		Full_End = MPI_WTIME()
		Elapsed  = Full_End - Full_Start
		Full_Min = Elapsed/60
		Full_Sec = MOD(Elapsed,60.0)
		WRITE(*,*) 
    WRITE(*,FMT="(A,I5,A,I2,A)") "LISM ENA RELAX complete in ",INT(Full_Min), " min ", INT(Full_Sec)," sec"
    WRITE(*,FMT="(A,I3,A,I5,A,I2,A)") "Time saved by using", np, " processors: ", &
    & INT((Elapsed*(np-1))/60)," min ", INT(MOD((Elapsed*(np-1)),60.0)), " sec"
    WRITE(*,*)
  END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (myid .EQ. 0) WRITE(*,*) 'Gathering data arrays to root'

	!! Gather all Data from processors to Root
	DO i=1,8
!		WRITE(*,*) 'rank: myN: sizes: ', myid, My_N_MC, size(My_Data(:,i)), size(Root_Data(:,i))
		CALL MPI_GATHER( My_Data(:,i), My_N_MC, MPI_DOUBLE_PRECISION,   &
	  &                Root_Data(:,i), My_N_MC, MPI_DOUBLE_PRECISION, &
	  &								 0, MPI_COMM_WORLD, ierr )

		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	END DO

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (myid .EQ. 0) WRITE(*,*) 'Reducing Collision Numbers to Root'

	IF (myid .EQ. 0) THEN	
		!! Communicate collision numbers to root
		root_N_AT = 0.0D0
		root_N_CX = 0.0D0

		DO i=1,size(Root_Data(:,5))
			root_N_AT = root_N_AT + Root_Data(i,5)	
			root_N_CX = root_N_CX + Root_Data(i,6)	
		END DO

		WRITE(*,'(A,ES11.2,A,ES11.2,A,ES11.2)') 'Number Elastic:',root_N_AT, &
		&       '  Number CX:',root_N_CX,'  Number Tot:',(root_N_AT+root_N_CX)
		WRITE(*,'(A)') '--------------------------------------------------------------------------'
		WRITE(*,'(A,F6.2)') '% Elastic Collisions: ', 100.0D0*root_N_AT/(root_N_AT+root_N_CX)
		WRITE(*,'(A,F6.2)') '% CX Collisions:      ', 100.0D0*root_N_CX/(root_N_AT+root_N_CX)
		WRITE(*,*) 
	END IF
	
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (myid .EQ. 0) WRITE(*,*) 'Reducing zone parameters to root'

	!! Reduce Sum zone parameters to root	
	CALL MPI_REDUCE( r_zone_E, root_r_zone_E, N_r_zone, MPI_DOUBLE_PRECISION, &
	&								 MPI_SUM, 0, MPI_COMM_WORLD, ierr )

	CALL MPI_REDUCE( r_zone_C, root_r_zone_C, N_r_zone, MPI_DOUBLE_PRECISION, &
	&								 MPI_SUM, 0, MPI_COMM_WORLD, ierr )

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Reduce Collision number arrays
	CALL MPI_REDUCE( My_N_CX_C_H, Root_N_CX_C_H, N_Coll_zone, MPI_DOUBLE_PRECISION, &
	&									MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( My_N_CX_E_H, Root_N_CX_E_H, N_Coll_zone, MPI_DOUBLE_PRECISION, &
	&									MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( My_N_EL_C_H, Root_N_EL_C_H, N_Coll_zone, MPI_DOUBLE_PRECISION, &
	&									MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( My_N_EL_E_H, Root_N_EL_E_H, N_Coll_zone, MPI_DOUBLE_PRECISION, &
	&									MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( My_N_CX_C_He, Root_N_CX_C_He, N_Coll_zone, MPI_DOUBLE_PRECISION, &
	&									MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( My_N_CX_E_He, Root_N_CX_E_He, N_Coll_zone, MPI_DOUBLE_PRECISION, &
	&									MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( My_N_EL_C_He, Root_N_EL_C_He, N_Coll_zone, MPI_DOUBLE_PRECISION, &
	&									MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( My_N_EL_E_He, Root_N_EL_E_He, N_Coll_zone, MPI_DOUBLE_PRECISION, &
	&									MPI_SUM, 0, MPI_COMM_WORLD, ierr )

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (myid .EQ. 0) WRITE(*,*) 'Reducing energy/distance zone arrays to root'

	!! Reduce Sum energy/distance zone array
	CALL MPI_REDUCE( My_Engy_Dist, Root_Engy_Dist, N_Engy_Dist*N_r_zone, & 
	& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (myid .EQ. 0) WRITE(*,*) 'Root writing data files'

	!! Write data to file if root
	IF (myid .EQ. 0) THEN

		DO i=1,N_Coll_zone
			IF (Root_N_EL_C_H(i) .EQ. 0.0D0) THEN
				WRITE(458,*) i, 0.0D0	
			ELSE	
				WRITE(458,*) i, Root_N_EL_E_H(i)/Root_N_EL_C_H(i) 	
			END IF
			IF (Root_N_CX_C_H(i) .EQ. 0.0D0) THEN
				WRITE(459,*) i, 0.0D0	
			ELSE	
				WRITE(459,*) i, Root_N_CX_E_H(i)/Root_N_CX_C_H(i) 	
			END IF
			IF (Root_N_EL_C_He(i) .EQ. 0.0D0) THEN
				WRITE(558,*) i, 0.0D0	
			ELSE	
				WRITE(558,*) i, Root_N_EL_E_He(i)/Root_N_EL_C_He(i) 	
			END IF
			IF (Root_N_CX_C_He(i) .EQ. 0.0D0) THEN
				WRITE(559,*) i, 0.0D0	
			ELSE	
				WRITE(559,*) i, Root_N_CX_E_He(i)/Root_N_CX_C_He(i) 	
			END IF
		END DO

		DO i=1,N_Part
			WRITE(666,*) Root_Data(i,1), Root_Data(i,2), Root_Data(i,3)
		END DO

		DO i=1,N_r_zone
			IF (root_r_zone_C(i) .EQ. 0) THEN
				WRITE(333,*) r_zone(i), 0
			ELSE
!				WRITE(333,*) r_zone(i), root_r_zone_E(i)*100.0D0/N_MC
				WRITE(333,*) r_zone(i), root_r_zone_E(i)/root_r_zone_C(i)
			END IF
		END DO

		DO i=1,N_r_zone
			DO j=1,N_Engy_Dist
				WRITE(555,*) r_zone(i), E_zone(j), Root_Engy_Dist(i,j)/N_MC
			END DO
		END DO

		!! Do histograms

		!! Determine Distance of Displacement
		ALLOCATE( Diss(N_Part) )
		CALL distance_finder( Root_Data(1:N_Part,1), Root_Data(1:N_Part,2), &
    &                     Root_Data(1:N_Part,3), N_Part, Diss ) 

		!! Find maximum displacment 
		max_diss = MAXVAL(Diss)

		!! Calculate Histograms for Dissplacement 
		!! Set number of Histogram bins
		N_hist = 500
		ALLOCATE( Hist_x(N_hist), Hist_y(N_hist) )
		CALL histogram( Diss, N_Part, 0.0D0, max_diss, N_hist, Hist_x, Hist_y )
		tot_per = 0.0D0	
		DO i=1,N_hist
			WRITE(777,*) Hist_x(i), Hist_y(i)/DBLE(N_Part)
			tot_per = tot_per + Hist_y(i)/DBLE(N_Part)
			WRITE(778,*) Hist_x(i), tot_per
		END DO		

		!! Calculate Histograms for Number of Collisions
		!! Set number of Histogram bins
		N_hist = 100
		max_collisions = MAXVAL(Root_Data(:,4))
		CALL histogram( Root_Data(1:N_Part,4), N_Part, 0.0D0, DBLE(max_collisions), N_hist, Hist_x, Hist_y )	
		tot_per = 0.0D0
		DO i=1,N_hist
			WRITE(888,*) Hist_x(i), Hist_y(i)/DBLE(N_Part)
			tot_per = tot_per + Hist_y(i)/DBLE(N_Part)
			WRITE(889,*) Hist_x(i), tot_per
		END DO

		!! Calculate Atom collision histogram for number of collisions
		N_hist = 100
		max_collisions = MAXVAL(Root_Data(:,5))
		CALL histogram( Root_Data(1:N_Part,5), N_Part, 0.0D0, DBLE(max_collisions), N_hist, Hist_x, Hist_y )	
		DO i=1,N_hist
			WRITE(891,*) Hist_x(i), Hist_y(i)/DBLE(N_Part)
		END DO

		!! Calculate CX collision histogram for number of collisions
		N_hist = 100
		max_collisions = MAXVAL(Root_Data(:,6))
		CALL histogram( Root_Data(1:N_Part,6), N_Part, 0.0D0, DBLE(max_collisions), N_hist, Hist_x, Hist_y )	
		DO i=1,N_hist
			WRITE(890,*) Hist_x(i), Hist_y(i)/DBLE(N_Part)
		END DO

	END IF !root

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (myid .EQ. 0) WRITE(*,*) 'Cleaning memory'
	
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Deallocate My_Data and Root_Data array, close data files
	DEALLOCATE( My_Data )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	DEALLOCATE( r_zone_E, r_zone_C, r_zone )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	DEALLOCATE( My_Engy_Dist, E_zone )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  DEALLOCATE( My_N_CX_E_H, My_N_CX_C_H, My_N_EL_E_H, My_N_EL_C_H )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	DEALLOCATE( My_N_CX_E_He, My_N_CX_C_He, My_N_EL_E_He, My_N_EL_C_He )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (myid .EQ. 0) THEN
		CLOSE( 333 )
		CLOSE( 458 )
		CLOSE( 459 )
		CLOSE( 558 )
		CLOSE( 559 )
		CLOSE( 555 )
		CLOSE( 567 )
		CLOSE( 666 )
		CLOSE( 777 )
		CLOSE( 778 )
		CLOSE( 888 )
		CLOSE( 889 )
		CLOSE( 890 )
		CLOSE( 891 )
		write(*,*) 'root done cleaning'
	END IF


	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE LISM_ENA_displacement

