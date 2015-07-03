!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D planetary MC relaxation experiment
!!
!! Complete, autonomous routine to calculate
!! planetary relaxation of ENAs including
!! secondary hot atoms. Escape fraction and
!! energy distributions are calculated
!! in a 3D planetary enviornment. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE planet_onestep_3d
	USE planet
	USE inputs, ONLY : DO_Escape_MC
	USE escape_trans
	USE secondary 
	USE rand_seed
	USE tables
	USE physics_constants
	USE click_data
	USE universal
	USE ena
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER				:: i, j, k, My_N_MC, tot, MC, N_Coll, Transport, GET_TCS
	INTEGER				:: N_r_zone, N_Engy_Dist, N_SHA_Engy_Dist, WRITE_1D
	INTEGER				:: Root_Escape_Count, Root_Therm_Count, Root_Second_Count
	INTEGER				:: ck, iT, iE, iH, ix, iy, iz
	INTEGEr				:: My_MC_Start, My_MC_End
	INTEGER,ALLOCATABLE,DIMENSION(:)	:: MC_THERM

	REAL(KIND=8),PARAMETER :: GB_per_real = 8.0D-9
	REAL(KIND=8),PARAMETER :: MB_per_real = 8.0D-6
	
	REAL(KIND=8)	:: x_0, y_0, z_0, V_0, V, Length
	REAL(KIND=8)	:: x1, y1, z1, x2, y2, z2, dL, dt
	REAL(KIND=8) 	:: E_now, z_now, pos_now, vz_now, vx_now
	REAL(KIND=8) 	:: E_nxt, z_nxt, pos_nxt, vz_nxt, vx_nxt
	REAL(KIND=8)	:: Den_now, theta_prev, phi_prev, ScattAng, TCS	
	REAL(KIND=8)	:: theta_now, phi_now, Lab_ScattAng, E_temp
	REAL(KIND=8)	:: phi, MFP, r, C1, E_targ, Theta_targ, Phi_targ	
	REAL(KIND=8)	:: Second_Ave, Therm_Ave, Planet_Ave, Escape_Ave 
	REAL(KIND=8)	:: HeH_TCS, HeO_TCS, HeHe_TCS, HCO2_TCS
	REAL(KIND=8)	:: Den_Tot, Den_He, Den_O, Den_H, Den_CO2, rand_targ
	REAL(KIND=8)	:: MFP_HeH, MFP_HCO2, MFP_HeO, MFP_HeHe
	REAL(KIND=8)	:: P_H, P_He, P_CO2, P_O, P_tot, PH_now, PH_nxt, TH_now, TH_nxt
	REAL(KIND=8)	:: MR_H, MR_He, MR_CO2, MR_O, h, dh, Len_U
	REAL(KIND=8)	:: R1_L, R1_U, R2_L, R2_U, R3_L, R3_U
	REAL(KIND=8)	:: Max_r_zone, dr_zone, Engy_Dist_In, Engy_Dist_Fn, d_Engy_Dist
	REAL(KIND=8)	:: Engy_Dist_SHA_In, Engy_Dist_SHA_Fn, d_Engy_Dist_SHA
	REAL(KIND=8)	:: ddE, ddT, ddX, ddY, ddUz, HHH, EEE, TTT, dE1, dE2
	REAL(KIND=8)	:: EE1, EE2, TT1, TT2, HH1, HH2, YY1, YY2, ZZ1, ZZ2 
	REAL(KIND=8)	:: Full_Start, Full_End, Full_Min, Full_Sec, Elapsed
	REAL(KIND=8)	:: Click_Start, Click_End, Click_Elapsed, Click_Hr, Click_Min, Click_Sec
	REAL(KIND=8)	:: lfg

	REAL(KIND=8),DIMENSION(4) :: Prob_Targ
  REAL(KIND=8),DIMENSION(3) :: x_now, x_nxt
  REAL(KIND=8),DIMENSION(3) :: u_now, u_nxt, u_targ

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: my_dummy, root_dummy

	!! Read in planetary inputs
	CALL read_planet_inputs

	!! Read Mars electron density tables
!	CALL read_e_den_table

	!! must have PROJ defined first
	CALL read_uni_tables
	CALL read_uni_tcs_tables
	CALL read_HpH_tcs_table
	CALL read_ENA_prod_height_table
!	CALL write_qm_tcs

	!! read density table for atmosphere
	CALL read_mars_density
 ! CALL test_mars_table_density

	!! Read in ENA starting location tables
	IF (DO_Escape_MC .EQ. 1) THEN
		CALL mars_sh_table_read	
!		CALL test_mars_sh_start_alt
	ELSE
		CALL mars_ena_table_read
	END IF

	!! Write run info to file if root
	IF ( myid .EQ. 0 ) CALL write_run_info

  !! Set up number of MC particles per rank
  My_N_MC = N_Part / np

	!! Set my MC start and MC end counters
	My_MC_Start = myid*My_N_MC + 1
	My_MC_End   = (myid + 1)*My_N_MC

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

	!! Allocate click data array
	CALL allocate_click

	ALLOCATE(MC_THERM(N_Part))
	MC_THERM(:) = 0

	!! Write welcome messages to screen
  IF ( myid == 0 ) THEN
    WRITE(*,'(A)') '#################################################'
    WRITE(*,'(A)') '####  STARTING 3D PLANETARY ENA SIMULATION   ####'
		WRITE(*,'(A)') '####             CLICK-BY-CLICK              ####'
    WRITE(*,'(A,I5,A)') '#### ',np, ' Processor(s) Being Used           ####'
    WRITE(*,'(A)') '#################################################'
    WRITE(*,*)
  END IF

	!! Allocate memory for SHA arrays
	CALL allocate_sha

	!! Fill SW velocity arrays for ENA or SH energy arrays for SH simulations
	IF (DO_Escape_MC .EQ. 1) THEN
		IF (ENA_SH .EQ. 1) THEN
			CALL read_SH_ED_Table( 'He ', 'H ')
		ELSE IF (ENA_SH .EQ. 4) THEN
			CALL read_SH_ED_Table( 'He ', 'He')
		END IF			
	ELSE
		CALL read_Vsw_table
	END IF

	!! Run ion transport for ENA production
!	CALL ion_setup( My_MC_Start, My_MC_End )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! set initial click data for all particles
	IF ( MONO_ENGY .EQ. 1) THEN
		CALL set_init_click_data( INIT_ENGY, h_0, SZA )
	ELSE IF ( MONO_ENGY .EQ. 2) THEN
		!! implament solar wind prob dist for init energy
		!! staring at h_0 with given SZA	
		IF (DO_Escape_MC .EQ. 1) THEN
!			CALL test_SH_Energy_Prob
			CALL set_init_click_data_SH	
		ELSE
			CALL set_init_click_data_Vsw( h_0, SZA )
		END IF
	ELSE 
		WRITE(*,*) 'Energy option: ', MONO_ENGY, ' not known!!'
	END IF

	!! set reference vectors fo click data
	CALL click_ref_vectors

	!! write reference vectors to file 
	IF (myid .EQ. 0) CALL write_click_ref_vectors

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! start simulation timer
	IF (myid .EQ. 0) Full_Start = MPI_WTIME()

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Loop through clicks
	DO ck=1,Max_click
		
		CLICK_1 = 0.0D0

		!! Main DO loop for initially hot MC particles
		DO MC=My_MC_Start,My_MC_End

			!! set current energy
			E_now 		= CLICK_0( MC, 2 )
			x_now(1) 	= CLICK_0( MC, 3 )
			x_now(2) 	= CLICK_0( MC, 4 )
			x_now(3) 	= CLICK_0( MC, 5 )
			u_now(1) 	= CLICK_0( MC, 6 )
			u_now(2) 	= CLICK_0( MC, 7 )
			u_now(3) 	= CLICK_0( MC, 8 )
			TH_now    = CLICK_0( MC, 9 )
			PH_now    = CLICK_0( MC, 10)

			!! Set SHA counter for current MC particle to 0
			My_SHA_count = 0

			!! Initialize SHA arrays to zero
			CALL clean_sha

!			IF ( (E_now .GE. E_Therm) .AND. (x_now(3) .LE. high) ) THEN
			IF ( x_now(3) .LE. high ) THEN
				CALL planet_onestep_transport( E_now, x_now(1), x_now(2), x_now(3), &
				&                      u_now(1), u_now(2), u_now(3), TH_now, PH_now, 2, &
				&                      E_nxt, x_nxt(1), x_nxt(2), x_nxt(3), u_nxt(1),  &
				&											u_nxt(2), u_nxt(3), TH_nxt, PH_nxt, dt )

				IF ( (x_nxt(3) .GE. high) .AND. (E_nxt .GE. E_Esc) ) THEN
					Esc_Engy_Dist( MC ) = E_nxt
					Esc_NColl_Dist( MC ) = CLICK_0( MC, 11 )
				END IF

				Len_U = SQRT( u_nxt(1)**2 + u_nxt(2)**2 + u_nxt(3)**2 )

				CLICK_1( MC, 1 ) = CLICK_0( MC, 1) + dt
				CLICK_1( MC, 2 ) = E_nxt
				CLICK_1( MC, 3 ) = x_nxt(1)
				CLICK_1( MC, 4 ) = x_nxt(2)
				CLICK_1( MC, 5 ) = x_nxt(3)
				CLICK_1( MC, 6 ) = u_nxt(1)/Len_U
				CLICK_1( MC, 7 ) = u_nxt(2)/Len_U
				CLICK_1( MC, 8 ) = u_nxt(3)/Len_U
				CLICK_1( MC, 9 ) = TH_nxt
				CLICK_1( MC, 10) = PH_nxt

				IF ( E_nxt .EQ. E_now ) THEN
					CLICK_1( MC, 11) = CLICK_0( MC, 11) 
				ELSE
					CLICK_1( MC, 11) = CLICK_0( MC, 11) + 1.0D0
				END IF

!				IF ( (CLICK_1(MC,1).LT.0.1D0) ) THEN
!					WRITE(100,*) CLICK_1( MC, 2 )
!				ELSE IF ((CLICK_1(MC,1).GT.0.1D0).AND. &
!				&        (CLICK_1(MC,1).LT.0.5D0)) THEN
!					WRITE(101,*) CLICK_1( MC, 2 )
!				ELSE IF ((CLICK_1(MC,1).GT.0.5D0).AND. &
!				&        (CLICK_1(MC,1).LT.0.7D0)) THEN
!					WRITE(102,*) CLICK_1( MC, 2 )
!				ELSE IF ((CLICK_1(MC,1).GT.0.7D0).AND. &
!				&        (CLICK_1(MC,1).LT.1.0D0)) THEN
!					WRITE(103,*) CLICK_1( MC, 2 )
!				END IF	

				IF ( (E_nxt .LE. E_therm) .AND. (MC_THERM(MC) .EQ. 0) ) THEN
					Therm_Time_Dist( MC )   = CLICK_1( MC, 1 )	
					Therm_Height_Dist( MC ) = CLICK_1( MC, 5 )
					MC_THERM(MC) = 1
				END IF

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
				END IF ! SHA

				!! fill phase space arrays
				!! if click is to be counted
				IF ( MOD(ck,Click_Skip) .EQ. 0 ) THEN
					CALL fill_phase_spaces(ck/Click_Skip,MC)
				END IF

				CALL fill_all_phase_spaces(MC)

			ELSE ! Escaped
				CLICK_1( MC, : ) = CLICK_0( MC, : )
!				IF (x_now(3) .GT. high) CLICK_1( MC, 5 ) = 1.1D0*high
			END IF ! Particle hot or thermal

!			CALL fill_all_phase_spaces(MC)
		
		END DO ! MC

		CLICK_0 = CLICK_1

		IF ( myid .EQ. 0 ) THEN

			IF ( ck .EQ. 10 ) THEN
				Click_End			= MPI_WTIME()
				Click_Elapsed = Max_Click*(Click_End - Full_Start)/10.0D0
				Click_hr      = FLOOR(Click_Elapsed/3600.0D0)
				Click_min     = FLOOR(MOD((Click_Elapsed/60.0D0),60.0))
				Click_sec     = MOD(Click_Elapsed,60.0)
				WRITE(*,*) 'Projected Simulation Time: '
				WRITE(*,'(I4,A,I4,A,I4,A)') INT(Click_hr), ' hours  ', INT(Click_min), ' minutes  ', INT(Click_sec), ' seconds '
				WRITE(*,*) 
			END IF

			IF ( MOD(ck,(Max_click/10)) .EQ. 0) THEN
				IF ( (ck/REAL(Max_click)) .EQ. 0.1 ) THEN
    			WRITE(*,*) 'Percent Complete'
					WRITE(*,*) '----------------'
					WRITE(*,*) '[*         ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.2) THEN
					WRITE(*,*) '[**        ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.3) THEN
					WRITE(*,*) '[***       ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.4) THEN
					WRITE(*,*) '[****      ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.5) THEN
					WRITE(*,*) '[*****     ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.6) THEN
					WRITE(*,*) '[******    ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.7) THEN
					WRITE(*,*) '[*******   ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.8) THEN
					WRITE(*,*) '[********  ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.9) THEN
					WRITE(*,*) '[********* ]'
				ELSE IF ( (ck/REAL(Max_click)) .EQ. 1.0) THEN
					WRITE(*,*) '[**********]'
				END IF
			END IF
		END IF

		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	END DO ! click

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr ) 
	!@@@@@@@@@@@@

	IF (myid .EQ. 0) THEN
		Full_End = MPI_WTIME()
		Elapsed  = Full_End - Full_Start
		Full_min = Elapsed/60
		Full_sec = MOD(Elapsed,60.0)
    WRITE(*,*)
    WRITE(*,FMT="(A,I5,A,I2,A)") "Planet Transport complete in ",INT(Full_Min), " min ", INT(Full_Sec)," sec"
    WRITE(*,FMT="(A,I3,A,I5,A,I2,A)") "Time saved by using", np, " processors: ", &
    & INT((Elapsed*(np-1))/60)," min ", INT(MOD((Elapsed*(np-1)),60.0)), " sec"
    WRITE(*,*)
  END IF

	!! Free memory for SHAs
	CALL free_sha

	!! Free uni memeory
	CALL clean_uni_tables
	CALL clean_uni_tcs_tables

	!! Free density tables
	CALL clean_mars_density

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			
	!@@@@@@@@@@@@

	IF (myid .EQ. 0) THEN
		WRITE(*,*)
		WRITE(*,*) "Start Writing Data Files"
		WRITE(*,*) "Max click: ", Max_click, " skip: ", click_skip, " Tot clicks: ",Num_Tot_Clicks
	END IF

	77 FORMAT (ES15.5)

	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@
	! DO ESC NColl DIST
	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@

	IF (myid .EQ. 0) ALLOCATE( ROOT_ESC_NColl_DIST(N_Part) )

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( Esc_NColl_Dist, ROOT_ESC_NColl_DIST, N_Part, &
	& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	IF (myid .EQ. 0) THEN
		OPEN(UNIT=76,FILE='../Data/Escape_NColl_Distribution.dat',ACCESS='APPEND')
		DO i=1,N_Part
			IF (ROOT_ESC_NColl_DIST(i) .NE. 0.0D0) WRITE(76,*) ROOT_ESC_NColl_DIST(i)	
		END DO
		CLOSE(76)
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@
	! DO ESC ENGY DIST
	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@

	IF (myid .EQ. 0) ALLOCATE( ROOT_ESC_ENGY_DIST(N_Part) )

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( Esc_Engy_Dist, ROOT_ESC_ENGY_DIST, N_Part, &
	& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	IF (myid .EQ. 0) THEN
		OPEN(UNIT=76,FILE='../Data/Escape_Energy_Distribution.dat',ACCESS='APPEND')
		DO i=1,N_Part
			IF (ROOT_ESC_ENGY_DIST(i) .NE. 0.0D0) WRITE(76,*) ROOT_ESC_ENGY_DIST(i)	
		END DO
		CLOSE(76)
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	
	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@
	! DO THERM TIME DIST
	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@

	IF (myid .EQ. 0) ALLOCATE( ROOT_Therm_Time_DIST(N_Part) )

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_REDUCE( Therm_Time_Dist, ROOT_Therm_Time_DIST, N_Part, &
	& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	IF (myid .EQ. 0) THEN
		OPEN(UNIT=76,FILE='../Data/Thermalization_Time_Distribution.dat',ACCESS='APPEND')
		DO i=1,N_Part
			IF (ROOT_Therm_Time_DIST(i) .NE. 0.0D0) WRITE(76,*) ROOT_Therm_Time_DIST(i)	
		END DO
		CLOSE(76)
	END IF
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  !@@@@@@@@@@@@@@@@@@@@@@@@
  !@@@@@@@@@@@@@@@@@@@@@@@@
  ! DO THERM HEIGHT DIST
  !@@@@@@@@@@@@@@@@@@@@@@@@
  !@@@@@@@@@@@@@@@@@@@@@@@@

  IF (myid .EQ. 0) ALLOCATE( ROOT_Therm_Height_DIST(N_Part) )

  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  CALL MPI_REDUCE( Therm_Height_Dist, ROOT_Therm_Height_DIST, N_Part, &
  & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  IF (myid .EQ. 0) THEN
    OPEN(UNIT=76,FILE='../Data/Thermalization_Height_Distribution.dat',ACCESS='APPEND')
    DO i=1,N_Part
      IF (ROOT_Therm_Height_DIST(i) .NE. 0.0D0) WRITE(76,*) ROOT_Therm_Height_DIST(i)
    END DO    
		CLOSE(76)
  END IF  
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@
	! Calculate reflection coefficient
	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@
	CALL reflection

	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@
	! DO CLICK DATA
	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@

	IF ( WRITE_E_T .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_C_E_T( Num_Tot_Clicks, Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( C_E_T(:,:,:), ROOT_C_E_T(:,:,:), Num_Tot_Clicks*Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=41,FILE="../Data/planet_energy_vs_time.dat",ACCESS="APPEND")
			!******************************
			! write E vs T to file	
			!******************************
			DO ck=1,Num_Tot_Clicks
				DO i=1,Num_Hist
					DO j=1,Num_Hist
						WRITE(41,*) ck*click_skip, i, j, ROOT_C_E_T(ck,i,j)	
					END DO ! j
				END DO ! i
			END DO ! ck
			CLOSE(41) 
			DEALLOCATE( ROOT_C_E_T )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_H_T .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_C_H_T( Num_Tot_Clicks, Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( C_H_T(:,:,:), ROOT_C_H_T(:,:,:), Num_Tot_Clicks*Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=42,FILE="../Data/planet_height_vs_time.dat",ACCESS="APPEND")
			!******************************
			! write H vs T to file	
			!******************************
			DO ck=1,Num_Tot_Clicks
				DO i=1,Num_Hist
					DO j=1,Num_Hist
						WRITE(42,*) ck*click_skip, i, j, ROOT_C_H_T(ck,i,j)	
					END DO ! j
				END DO ! i
			END DO ! ck
			CLOSE(42)
			DEALLOCATE( ROOT_C_H_T )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_H_E .EQ. 1) THEN
		IF ( myid .EQ. 0 ) THEN
			ALLOCATE( ROOT_C_H_E( Num_Tot_Clicks, Num_Hist, Num_Hist) )
			ALLOCATE( root_dummy( Num_Hist, Num_Hist ) )
		END IF
		ALLOCATE( my_dummy( Num_Hist, Num_Hist ) )
		DO ck=1,Num_Tot_Clicks
			my_dummy(:,:) = C_H_E(ck,:,:)	
			CALL MPI_REDUCE( my_dummy, root_dummy, Num_Hist*Num_Hist, &
					 & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
			IF (myid .EQ. 0) ROOT_C_H_E(ck,:,:) = root_dummy(:,:)
			CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		END DO
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=43,FILE="../Data/planet_height_vs_energy.dat",ACCESS="APPEND")
			!******************************
			! write H vs E to file	
			!******************************
			DO ck=1,Num_Tot_Clicks
				DO i=1,Num_Hist
					DO j=1,Num_Hist
						WRITE(43,*) ck*click_skip, i, j, ROOT_C_H_E(ck,i,j)	
					END DO ! j
				END DO ! i
			END DO ! ck
			CLOSE(43)
			DEALLOCATE( ROOT_C_H_E )
			DEALLOCATE( root_dummy )
		END IF ! root
		DEALLOCATE( my_dummy )
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_H_dE .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_C_H_dE( Num_Tot_Clicks, Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( C_H_dE(:,:,:), ROOT_C_H_dE(:,:,:), Num_Tot_Clicks*Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=44,FILE="../Data/planet_height_vs_energy_loss.dat",ACCESS="APPEND")
			!******************************
			! write H vs dE to file	
			!******************************
			DO ck=1,Num_Tot_Clicks
				DO i=1,Num_Hist
					DO j=1,Num_Hist
						WRITE(44,*) ck*click_skip, i, j, ROOT_C_H_dE(ck,i,j)	
					END DO ! j
				END DO ! i
			END DO ! ck
			CLOSE(44)
			DEALLOCATE( ROOT_C_H_dE )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_Ux_T .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_C_Ux_T( Num_Tot_Clicks, Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( C_Ux_T(:,:,:), ROOT_C_Ux_T(:,:,:), Num_Tot_Clicks*Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=45,FILE="../Data/planet_vert_vel_vs_time.dat",ACCESS="APPEND")
			!******************************
			! write Ux vs T to file	
			!******************************
			DO ck=1,Num_Tot_Clicks
				DO i=1,Num_Hist
					DO j=1,Num_Hist
						WRITE(45,*) ck*click_skip, i, j, ROOT_C_Ux_T(ck,i,j)	
					END DO ! j
				END DO ! i
			END DO ! ck
			CLOSE(45)
			DEALLOCATE( ROOT_C_Ux_T )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_H_Ux .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_C_H_Ux( Num_Tot_Clicks, Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( C_H_Ux(:,:,:), ROOT_C_H_Ux(:,:,:), Num_Tot_Clicks*Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=46,FILE="../Data/planet_height_vs_vert_vel.dat",ACCESS="APPEND")
			!******************************
			! write H vs Ux to file	
			!******************************
			DO ck=1,Num_Tot_Clicks
				DO i=1,Num_Hist
					DO j=1,Num_Hist
						WRITE(46,*) ck*click_skip, i, j, ROOT_C_H_Ux(ck,i,j)	
					END DO ! j
				END DO ! i
			END DO ! ck
			CLOSE(46)
			DEALLOCATE( ROOT_C_H_Ux )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_N_T .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_C_N_T( Num_Tot_Clicks, Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( C_N_T(:,:,:), ROOT_C_N_T(:,:,:), Num_Tot_Clicks*Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=47,FILE="../Data/planet_Ncoll_vs_time.dat",ACCESS="APPEND")
			!******************************
			! write N vs T to file	
			!******************************
			DO ck=1,Num_Tot_Clicks
				DO i=1,Num_Hist
					DO j=1,Num_Hist
						WRITE(47,*) ck*click_skip, i, j, ROOT_C_N_T(ck,i,j)	
					END DO ! j
				END DO ! i
			END DO ! ck
			CLOSE(47)
			DEALLOCATE( ROOT_C_N_T )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_H_N .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_C_H_N( Num_Tot_Clicks, Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( C_H_N(:,:,:), ROOT_C_H_N(:,:,:), Num_Tot_Clicks*Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=48,FILE="../Data/planet_height_vs_Ncoll.dat",ACCESS="APPEND")
			!******************************
			! write N vs T to file	
			!******************************
			DO ck=1,Num_Tot_Clicks
				DO i=1,Num_Hist
					DO j=1,Num_Hist
						WRITE(48,*) ck*click_skip, i, j, ROOT_C_H_N(ck,i,j)	
					END DO ! j
				END DO ! i
			END DO ! ck
			CLOSE(48)
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@
	! DO ALL CLICK DATA
	!@@@@@@@@@@@@@@@@@@@@@@@@
	!@@@@@@@@@@@@@@@@@@@@@@@@

	IF ( WRITE_ALL_E_T .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_A_E_T( Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( A_E_T(:,:), ROOT_A_E_T(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=41,FILE="../Data/planet_all_energy_vs_time.dat",ACCESS="APPEND")
			!******************************
			! write E vs T to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(41,*) i, j, ROOT_A_E_T(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(41) 
			DEALLOCATE( ROOT_A_E_T )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_ALL_H_T .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_A_H_T( Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( A_H_T(:,:), ROOT_A_H_T(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=42,FILE="../Data/planet_all_height_vs_time.dat",ACCESS="APPEND")
			!******************************
			! write H vs T to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(42,*) i, j, ROOT_A_H_T(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(42)
			DEALLOCATE( ROOT_A_H_T )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_ALL_H_E .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_A_H_E( Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( A_H_E(:,:), ROOT_A_H_E(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=43,FILE="../Data/planet_all_height_vs_energy.dat",ACCESS="APPEND")
			!******************************
			! write H vs E to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(43,*) i, j, ROOT_A_H_E(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(43)
			DEALLOCATE( ROOT_A_H_E )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_ALL_H_dE .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_A_H_dE( Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( A_H_dE(:,:), ROOT_A_H_dE(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=44,FILE="../Data/planet_all_height_vs_energy_loss.dat",ACCESS="APPEND")
			!******************************
			! write H vs dE to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(44,*) i, j, ROOT_A_H_dE(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(44)
			DEALLOCATE( ROOT_A_H_dE )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_ALL_Ux_T .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_A_Ux_T( Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( A_Ux_T(:,:), ROOT_A_Ux_T(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=45,FILE="../Data/planet_all_vert_vel_vs_time.dat",ACCESS="APPEND")
			!******************************
			! write Ux vs T to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(45,*) i, j, ROOT_A_Ux_T(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(45)
			DEALLOCATE( ROOT_A_Ux_T )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_ALL_H_Ux .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_A_H_Ux( Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( A_H_Ux(:,:), ROOT_A_H_Ux(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=46,FILE="../Data/planet_all_height_vs_vert_vel.dat",ACCESS="APPEND")
			!******************************
			! write H vs Ux to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(46,*) i, j, ROOT_A_H_Ux(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(46)
			DEALLOCATE( ROOT_A_H_Ux )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_ALL_N_T .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_A_N_T( Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( A_N_T(:,:), ROOT_A_N_T(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=47,FILE="../Data/planet_all_Ncoll_vs_time.dat",ACCESS="APPEND")
			!******************************
			! write N vs T to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(47,*) i, j, ROOT_A_N_T(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(47)
			DEALLOCATE( ROOT_A_N_T )
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_ALL_H_N .EQ. 1) THEN
		IF ( myid .EQ. 0 ) ALLOCATE( ROOT_A_H_N( Num_Hist, Num_Hist) )
		CALL MPI_REDUCE( A_H_N(:,:), ROOT_A_H_N(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=48,FILE="../Data/planet_all_height_vs_Ncoll.dat",ACCESS="APPEND")
			!******************************
			! write N vs T to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(48,*) i, j, ROOT_A_H_N(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(48)
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_SHA_H_Ux .EQ. 1) THEN
		IF ( myid .EQ. 0 ) THEN
			ALLOCATE( ROOT_C_SHA_H_Ux_H(   Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_Ux_He(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_Ux_O(   Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_Ux_Ar(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_Ux_H2(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_Ux_N2(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_Ux_CO(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_Ux_CO2( Num_Hist, Num_Hist) )
		END IF
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	
		CALL MPI_REDUCE( C_SHA_H_Ux_H(:,:), ROOT_C_SHA_H_Ux_H(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_Ux_He(:,:), ROOT_C_SHA_H_Ux_He(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_Ux_O(:,:), ROOT_C_SHA_H_Ux_O(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_Ux_Ar(:,:), ROOT_C_SHA_H_Ux_Ar(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_Ux_H2(:,:), ROOT_C_SHA_H_Ux_H2(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_Ux_N2(:,:), ROOT_C_SHA_H_Ux_N2(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_Ux_CO(:,:), ROOT_C_SHA_H_Ux_CO(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_Ux_CO2(:,:), ROOT_C_SHA_H_Ux_CO2(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=40,FILE="../Data/planet_SHA_height_vs_Ux_H.dat",ACCESS="APPEND")
			OPEN(UNIT=41,FILE="../Data/planet_SHA_height_vs_Ux_He.dat",ACCESS="APPEND")
			OPEN(UNIT=42,FILE="../Data/planet_SHA_height_vs_Ux_O.dat",ACCESS="APPEND")
			OPEN(UNIT=43,FILE="../Data/planet_SHA_height_vs_Ux_Ar.dat",ACCESS="APPEND")
			OPEN(UNIT=44,FILE="../Data/planet_SHA_height_vs_Ux_H2.dat",ACCESS="APPEND")
			OPEN(UNIT=45,FILE="../Data/planet_SHA_height_vs_Ux_N2.dat",ACCESS="APPEND")
			OPEN(UNIT=46,FILE="../Data/planet_SHA_height_vs_Ux_CO.dat",ACCESS="APPEND")
			OPEN(UNIT=47,FILE="../Data/planet_SHA_height_vs_Ux_CO2.dat",ACCESS="APPEND")
			!******************************
			! write SHA H vs Ux to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(40,*) X_CLICK(i,5), X_CLICK(j,8), ROOT_C_SHA_H_Ux_H(i,j)	
					WRITE(41,*) X_CLICK(i,5), X_CLICK(j,8), ROOT_C_SHA_H_Ux_He(i,j)	
					WRITE(42,*) X_CLICK(i,5), X_CLICK(j,8), ROOT_C_SHA_H_Ux_O(i,j)	
					WRITE(43,*) X_CLICK(i,5), X_CLICK(j,8), ROOT_C_SHA_H_Ux_Ar(i,j)	
					WRITE(44,*) X_CLICK(i,5), X_CLICK(j,8), ROOT_C_SHA_H_Ux_H2(i,j)	
					WRITE(45,*) X_CLICK(i,5), X_CLICK(j,8), ROOT_C_SHA_H_Ux_N2(i,j)	
					WRITE(46,*) X_CLICK(i,5), X_CLICK(j,8), ROOT_C_SHA_H_Ux_CO(i,j)	
					WRITE(47,*) X_CLICK(i,5), X_CLICK(j,8), ROOT_C_SHA_H_Ux_CO2(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(40)
			CLOSE(41)
			CLOSE(42)
			CLOSE(43)
			CLOSE(44)
			CLOSE(45)
			CLOSE(46)
			CLOSE(47)
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF ( WRITE_SHA_H_E .EQ. 1) THEN
		IF ( myid .EQ. 0 ) THEN
			ALLOCATE( ROOT_C_SHA_H_E_H(   Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_E_He(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_E_O(   Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_E_Ar(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_E_H2(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_E_N2(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_E_CO(  Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_C_SHA_H_E_CO2( Num_Hist, Num_Hist) )
		END IF
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_E_H(:,:), ROOT_C_SHA_H_E_H(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_E_He(:,:), ROOT_C_SHA_H_E_He(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_E_O(:,:), ROOT_C_SHA_H_E_O(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_E_Ar(:,:), ROOT_C_SHA_H_E_Ar(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_E_H2(:,:), ROOT_C_SHA_H_E_H2(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_E_N2(:,:), ROOT_C_SHA_H_E_N2(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_E_CO(:,:), ROOT_C_SHA_H_E_CO(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_REDUCE( C_SHA_H_E_CO2(:,:), ROOT_C_SHA_H_E_CO2(:,:), Num_Hist*Num_Hist, &
		& MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
		IF ( myid .EQ. 0 ) THEN
			OPEN(UNIT=50,FILE="../Data/planet_SHA_height_vs_energy_H.dat",ACCESS="APPEND")
			OPEN(UNIT=51,FILE="../Data/planet_SHA_height_vs_energy_He.dat",ACCESS="APPEND")
			OPEN(UNIT=52,FILE="../Data/planet_SHA_height_vs_energy_O.dat",ACCESS="APPEND")
			OPEN(UNIT=53,FILE="../Data/planet_SHA_height_vs_energy_Ar.dat",ACCESS="APPEND")
			OPEN(UNIT=54,FILE="../Data/planet_SHA_height_vs_energy_H2.dat",ACCESS="APPEND")
			OPEN(UNIT=55,FILE="../Data/planet_SHA_height_vs_energy_N2.dat",ACCESS="APPEND")
			OPEN(UNIT=56,FILE="../Data/planet_SHA_height_vs_energy_CO.dat",ACCESS="APPEND")
			OPEN(UNIT=57,FILE="../Data/planet_SHA_height_vs_energy_CO2.dat",ACCESS="APPEND")
			!******************************
			! write SHA H vs Ux to file	
			!******************************
			DO i=1,Num_Hist
				DO j=1,Num_Hist
					WRITE(50,*) X_CLICK(i,5), X_CLICK(j,12), ROOT_C_SHA_H_E_H(i,j)	
					WRITE(51,*) X_CLICK(i,5), X_CLICK(j,12), ROOT_C_SHA_H_E_He(i,j)	
					WRITE(52,*) X_CLICK(i,5), X_CLICK(j,12), ROOT_C_SHA_H_E_O(i,j)	
					WRITE(53,*) X_CLICK(i,5), X_CLICK(j,12), ROOT_C_SHA_H_E_Ar(i,j)	
					WRITE(54,*) X_CLICK(i,5), X_CLICK(j,12), ROOT_C_SHA_H_E_H2(i,j)	
					WRITE(55,*) X_CLICK(i,5), X_CLICK(j,12), ROOT_C_SHA_H_E_N2(i,j)	
					WRITE(56,*) X_CLICK(i,5), X_CLICK(j,12), ROOT_C_SHA_H_E_CO(i,j)	
					WRITE(57,*) X_CLICK(i,5), X_CLICK(j,12), ROOT_C_SHA_H_E_CO2(i,j)	
				END DO ! j
			END DO ! i
			CLOSE(50)
			CLOSE(51)
			CLOSE(52)
			CLOSE(53)
			CLOSE(54)
			CLOSE(55)
			CLOSE(56)
			CLOSE(57)
		END IF ! root
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			

	IF (myid .EQ. 0) THEN
		WRITE(*,*)
		WRITE(*,*) 'data files written'
	END IF

!	CALL clean_e_den_table
	CALL clean_HpH_tcs_table
	CALL clean_click
	CALL clean_ENA_prod_height_table
	
	DEALLOCATE(MC_THERM)

	!@@@@@@@@@@@@
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )			
	!@@@@@@@@@@@@

END SUBROUTINE planet_onestep_3d

