!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! Transport particle in planetary atmosphere
!! 
!! Step-by-step methods used in this transport 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE planet_onestep_transport( E0, x0, y0, z0, ux0, uy0, uz0, t0, p0, SHA_ON, E1, x1, y1, z1, ux1, uy1, uz1, t1, p1, dt ) 
	USE planet
	USE physics_constants
	USE secondary
	USE click_data
	USE rand_seed
	USE mpi_info

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)				:: E0								!! initial energy [eV]
	REAL(KIND=8)				:: x0								!! initial position [m]
	REAL(KIND=8)				:: y0								!! initial position [m]
	REAL(KIND=8)				:: z0								!! initial position [m]
	REAL(KIND=8)				:: ux0							!! initial unit velocity [m/s]
	REAL(KIND=8)				:: uy0							!! initial unit velocity [m/s]
	REAL(KIND=8)				:: uz0							!! initial unit velocity [m/s]
	REAL(KIND=8)				:: t0								!! initial theta angle [rad]
	REAL(KIND=8)				:: p0								!! initial phi angle [rad]
	INTEGER							:: SHA_ON						!! SHA tracking on if = 1 off if = 0
																					!! if SHA = 2 write data to SHA arrays
																					!! instead of ENA arrays
	!! Outputs
	REAL(KIND=8)				:: E1								!! final energy [eV]
	REAL(KIND=8)				:: x1								!! final position [m]
	REAL(KIND=8)				:: y1								!! final position [m]
	REAL(KIND=8)				:: z1								!! final position [m]
	REAL(KIND=8)				:: ux1							!! final unit velocity [m/s]
	REAL(KIND=8)				:: uy1							!! final unit velocity [m/s]
	REAL(KIND=8)				:: uz1							!! final unit velocity [m/s]
	REAL(KIND=8)				:: t1								!! final theta angle [rad]
	REAL(KIND=8)				:: p1								!! final phi angle [rad]
	REAL(KIND=8)				:: dt								!! time elapsed in lab frame during one step [sec]

	!! Internal
	REAL(KIND=8),DIMENSION(3)	:: x_now, u_now, u_targ, u_nxt, x_nxt, dum_x
	REAL(KIND=8),DIMENSION(8)	:: Prob_Targ
	REAL(KIND=8)							:: lfg
	REAL(KIND=8)							:: theta_now, phi_now, elec_den, HpHTCS
	REAL(KIND=8)							:: theta_prev, phi_prev, E_now, E_nxt, E_targ
	REAL(KIND=8)							:: HCO2_TCS, HeH_TCS, HeO_TCS, HeHe_TCS
	REAL(KIND=8)							:: TCS_X_CO2, TCS_X_CO, TCS_X_H2, TCS_X_N2 
	REAL(KIND=8)							:: TCS_X_O, TCS_X_He, TCS_X_H, TCS_X_Ar 
	REAL(KIND=8)							:: Den_CO2, Den_H2, Den_N2, Den_CO
	REAL(KIND=8)							:: Den_H, Den_O, Den_He, Den_Ar, Den_Tot
	REAL(KIND=8)							:: MFP_X_CO2, MFP_X_CO, MFP_X_H2, MFP_X_N2 
	REAL(KIND=8)							:: MFP_X_O, MFP_X_He, MFP_X_H, MFP_X_Ar, MFP 
	REAL(KIND=8)							:: MR_CO2, MR_CO, MR_H2, MR_N2, phi
	REAL(KIND=8)							:: MR_O, MR_He, MR_H, MR_Ar
	REAL(KIND=8)							:: P_Tot, rand_targ, Theta_targ, Phi_targ
	REAL(KIND=8)							:: C1, r, dL, Length, ScattAng, Lab_ScattAng
	REAL(KIND=8)							:: dE, U_tot, velocity, gravity, Target_Theta
	REAL(KIND=8)							:: Target_Phi, r_theta, r_phi, Target_ux0
	REAL(KIND=8)							:: Target_uy0, Target_uz0, dum_dx
	REAL(KIND=8)							:: pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8
	INTEGER										:: N_Coll, Transport, GET_TCS	
	INTEGER										:: i, j, xyz, NOW, ii, jj
	
	!! Set steplength for transport [m]
	dL          = 1.0D3	

	E_now       = E0
	x_now(1) 		= x0
	x_now(2) 		= y0
	x_now(3) 		= z0

	u_now(1) 		= ux0
	u_now(2) 		= uy0
	u_now(3) 		= uz0

	U_tot       = SQRT(ux0**2 + uy0**2 + uz0**2)
	ux0         = ux0/U_tot
	uy0         = uy0/U_tot
	uz0         = uz0/U_tot

	theta_prev 	= t0
	phi_prev 		= p0

	!! Init collision counter to zero
	N_Coll = 0

	!! Set transport/TCS to ON
	GET_TCS   = 1

	!! Find total cross section [m^2] if energy has changed 
	IF ( GET_TCS .EQ. 1) THEN
		IF (CM_TYPE .EQ. 'QM') THEN	
			CALL get_qm_tcs(E_now,TCS_X_CO2,TCS_X_CO,TCS_X_N2,TCS_X_H2,TCS_X_O,TCS_X_Ar,TCS_X_He,TCS_X_H)
		ELSE IF (CM_TYPE .EQ. 'HS') THEN
			CALL get_hs_tcs(TCS_X_CO2,TCS_X_CO,TCS_X_N2,TCS_X_H2,TCS_X_O,TCS_X_Ar,TCS_X_He,TCS_X_H)
		END IF ! HS or QM TCS
	END IF ! if GET_TCS
	
	CALL mars_table_density( 'H  ', x_now(3), Den_H   )
	CALL mars_table_density( 'He ', x_now(3), Den_He  )
	CALL mars_table_density( 'O  ', x_now(3), Den_O   )
	CALL mars_table_density( 'Ar ', x_now(3), Den_Ar  )
	CALL mars_table_density( 'H2 ', x_now(3), Den_H2  )
	CALL mars_table_density( 'N2 ', x_now(3), Den_N2  )
	CALL mars_table_density( 'CO ', x_now(3), Den_CO  )
	CALL mars_table_density( 'CO2', x_now(3), Den_CO2 )

	MFP_X_CO2 = Den_CO2 *TCS_X_CO2
	MFP_X_CO  = Den_CO  *TCS_X_CO
	MFP_X_H2  = Den_H2  *TCS_X_H2
	MFP_X_N2  = Den_N2  *TCS_X_N2
	MFP_X_O   = Den_O   *TCS_X_O
	MFP_X_He  = Den_He  *TCS_X_He
	MFP_X_Ar  = Den_Ar  *TCS_X_Ar
	MFP_X_H   = Den_H   *TCS_X_H

	Den_Tot = Den_CO2 + Den_H + Den_O + Den_He + Den_CO + Den_H2 + Den_N2 + Den_Ar

	!! Total mean free path for all target species
	MFP = 1.0D0/(MFP_X_CO2+MFP_X_CO+MFP_X_N2+MFP_X_H2+MFP_X_O+MFP_X_Ar+MFP_X_He+MFP_X_H)

	!! if 20% of mfp is less than 1km, use that as steplength
	IF (0.2D0*MFP .LT. 1.0D3) THEN
		dL = 0.2D0*MFP
	ELSE
		dL = 1.0D3
	END IF

	C1  = EXP( -dL/MFP )
	r   = lfg()		

	!! Determine if collision occurs in current steplength  
	IF ( r .GT. C1 ) THEN
		!!-------------
		!! Collision
		!!-------------
		!! Find exact location of collision
	!	WRITE(*,*) myid, ' COLLISION ', x_now(3)
		Length = -MFP*LOG(r)

		CALL energy_to_velocity( E_now, MP, velocity )
		dt = Length/velocity

		!! Find mixing ratios for target species at collision altitude
		MR_CO2 = Den_CO2/Den_Tot
		MR_CO  = Den_CO /Den_Tot
		MR_He  = Den_He /Den_Tot
		MR_Ar  = Den_Ar /Den_Tot
		MR_H   = Den_H  /Den_Tot
		MR_H2  = Den_H2 /Den_Tot
		MR_N2  = Den_N2 /Den_Tot
		MR_O   = Den_O  /Den_Tot

		!! Find probability array for target species collisions 
		pp1 		= MR_CO2*TCS_X_CO2 
		pp2 		= MR_CO*TCS_X_CO 
		pp3 		= MR_H2*TCS_X_H2 
		pp4 		= MR_N2*TCS_X_N2 
		pp5 		= MR_O*TCS_X_O 
		pp6 		= MR_He*TCS_X_He 
		pp7 		= MR_Ar*TCS_X_Ar 
		pp8 		= MR_H*TCS_X_H
		P_Tot 	= pp1 + pp2 + pp3 + pp4 + pp5 + pp6 + pp7 + pp8
		Prob_Targ(1) = pp1/P_Tot
		Prob_Targ(2) = Prob_Targ(1) + pp2/P_Tot
		Prob_Targ(3) = Prob_Targ(2) + pp3/P_Tot
		Prob_Targ(4) = Prob_Targ(3) + pp4/P_Tot
		Prob_Targ(5) = Prob_Targ(4) + pp5/P_Tot
		Prob_Targ(6) = Prob_Targ(5) + pp6/P_Tot
		Prob_Targ(7) = Prob_Targ(6) + pp7/P_Tot
		Prob_Targ(8) = Prob_Targ(7) + pp8/P_Tot

		!! Find random target atom
		rand_targ = lfg()
		IF (rand_targ .LT. Prob_Targ(1)) THEN
			!! Target is CO2
			Targ = 'CO2'
			atom = 2
		ELSE IF ( (rand_targ .GT. Prob_Targ(1)) .AND. (rand_targ .LT. Prob_Targ(2)) ) THEN
			!! Target is CO
			Targ = 'CO '
			atom = 2
		ELSE IF ( (rand_targ .GT. Prob_Targ(2)) .AND. (rand_targ .LT. Prob_Targ(3)) ) THEN
			!! Target is H2
			Targ = 'H2 '
			atom = 2
		ELSE IF ( (rand_targ .GT. Prob_Targ(3)) .AND. (rand_targ .LT. Prob_Targ(4)) ) THEN
			!! Target is N2 
			Targ = 'N2 '
			atom = 2
		ELSE IF ( (rand_targ .GT. Prob_Targ(4)) .AND. (rand_targ .LT. Prob_Targ(5)) ) THEN
			!! Target is O
			Targ = 'O  '
			atom = 1 
		ELSE IF ( (rand_targ .GT. Prob_Targ(5)) .AND. (rand_targ .LT. Prob_Targ(6)) ) THEN
			!! Target is He 
			Targ = 'He '
			atom = 1 
		ELSE IF ( (rand_targ .GT. Prob_Targ(6)) .AND. (rand_targ .LT. Prob_Targ(7)) ) THEN
			!! Target is Ar 
			Targ = 'Ar '
			atom = 1 
		ELSE IF ( (rand_targ .GT. Prob_Targ(7)) .AND. (rand_targ .LT. Prob_Targ(8)) ) THEN
			!! Target is H 
			Targ = 'H  '
			atom = 1 
		END IF

		!! Find mass of target atom
		CALL mass_finder

		!! Find reduced mass [kg]
		mass = MU*TOAMU*AMUTOKG

		!! Update collision counter
		Coll_Count = Coll_Count + 1

		!!!!!!!!!!!
		!! Get random scattering angle in CM frame for given collision
		!!!!!!!!!!!
		!! use 1 keV for all energies > 1keV until tables complete
		!!!!!!!!!!!
		IF (CM_TYPE .EQ. 'QM') THEN
			IF (TRIM(PROJ) .EQ. 'H') THEN
				IF (TRIM(Targ) .EQ. 'H') THEN
					CALL lin_rand_angle( E_now, 'HH   ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'He') THEN
					CALL lin_rand_angle( E_now, 'HHe  ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'O') THEN
					CALL lin_rand_angle( E_now, 'HO   ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'Ar') THEN
					CALL lin_rand_angle( E_now, 'HAr  ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'H2') THEN
					CALL lin_rand_angle( E_now, 'HH2  ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'N2') THEN
					CALL lin_rand_angle( E_now, 'HN2  ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'CO') THEN
					CALL lin_rand_angle( E_now, 'HCO  ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'CO2') THEN
					CALL lin_rand_angle( E_now, 'HCO2 ', ScattAng )
				END IF
			ELSE IF (TRIM(PROJ) .EQ. 'He') THEN
				IF (TRIM(Targ) .EQ. 'H') THEN
					CALL lin_rand_angle( E_now, 'HeH  ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'He') THEN
					CALL lin_rand_angle( E_now, 'HeHe ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'O') THEN
					CALL lin_rand_angle( E_now, 'HeO  ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'Ar') THEN
					CALL lin_rand_angle( E_now, 'HeAr ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'H2') THEN
					CALL lin_rand_angle( E_now, 'HeH2 ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'N2') THEN
					CALL lin_rand_angle( E_now, 'HeN2 ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'CO') THEN
					CALL lin_rand_angle( E_now, 'HeCO ', ScattAng )
				ELSE IF (TRIM(Targ) .EQ. 'CO2') THEN
					CALL lin_rand_angle( E_now, 'HeCO2', ScattAng )
				END IF
			END IF
		ELSE IF (CM_TYPE .EQ. 'HS') THEN
			CALL HS_rand_angle( ScattAng )
		END IF

		!! Find new projectile energy
		CALL find_new_energy( E_now, ScattAng, MP/MT, E_nxt )
		dE = E_now - E_nxt

		IF ( (WRITE_SHA_H_Ux .EQ. 1) .AND. (dE .GE. E_Therm) ) THEN
			Target_Theta = (PI-ScattAng)/2.0D0 ! target scattering angle in lab frame [rad]
			Target_Phi   = lfg()*2.0D0*PI
			CALL particle_dir_cos_transport(u_now,x_now,Target_Theta,Target_Phi,dum_dx,u_targ,dum_x)
			Target_ux0   = u_targ(1)
			Target_uy0   = u_targ(2)
			Target_uz0   = u_targ(3)
			NOW = 0
			ii  = 0
			DO WHILE ( NOW .EQ. 0 ) 
				ii = ii + 1
				IF ( (x_now(3) .LE. X_CLICK(ii,5)) .OR. (ii .GE. Num_Hist) ) THEN
					NOW = 1	
				END IF 	
			END DO	
			NOW = 0
			jj  = 0
			DO WHILE ( NOW .EQ. 0 ) 
				jj = jj + 1
				IF ( (Target_uz0 .LE. X_CLICK(jj,8)) .OR. (jj .GE. Num_Hist) ) THEN
					NOW = 1	
				END IF 	
			END DO
			IF (TRIM(Targ) .EQ. 'H') THEN
				C_SHA_H_Ux_H(ii,jj) = C_SHA_H_Ux_H(ii,jj) + 1.0D0/N_Part	
			ELSE IF (TRIM(Targ) .EQ. 'He') THEN
				C_SHA_H_Ux_He(ii,jj) = C_SHA_H_Ux_He(ii,jj) + 1.0D0/N_Part	
			ELSE IF (TRIM(Targ) .EQ. 'O') THEN
				C_SHA_H_Ux_O(ii,jj) = C_SHA_H_Ux_O(ii,jj) + 1.0D0/N_Part	
			ELSE IF (TRIM(Targ) .EQ. 'Ar') THEN
				C_SHA_H_Ux_Ar(ii,jj) = C_SHA_H_Ux_Ar(ii,jj) + 1.0D0/N_Part	
			ELSE IF (TRIM(Targ) .EQ. 'H2') THEN
				C_SHA_H_Ux_H2(ii,jj) = C_SHA_H_Ux_H2(ii,jj) + 1.0D0/N_Part	
			ELSE IF (TRIM(Targ) .EQ. 'N2') THEN
				C_SHA_H_Ux_N2(ii,jj) = C_SHA_H_Ux_N2(ii,jj) + 1.0D0/N_Part	
			ELSE IF (TRIM(Targ) .EQ. 'CO') THEN
				C_SHA_H_Ux_CO(ii,jj) = C_SHA_H_Ux_CO(ii,jj) + 1.0D0/N_Part	
			ELSE IF (TRIM(Targ) .EQ. 'CO2') THEN
				C_SHA_H_Ux_CO2(ii,jj) = C_SHA_H_Ux_CO2(ii,jj) + 1.0D0/N_Part	
			ELSE
				WRITE(*,*) 'UNKNOWN TARGET IN H vs Uz SHA, TARGET: ', Targ
			END IF
		END IF ! write_sha_H_UX on

		IF ( (WRITE_SHA_H_E .EQ. 1) .AND. (dE .GE. E_Therm) ) THEN
			NOW = 0
			ii  = 0
			DO WHILE ( NOW .EQ. 0 )
				ii = ii + 1
				IF ( (x_now(3) .LE. X_CLICK(ii,5)) .OR. (ii .GE. Num_Hist)) THEN
					NOW = 1
				END IF
			END DO
			NOW = 0
			jj  = 0
			DO WHILE ( NOW .EQ. 0 )
				jj = jj + 1
				IF ( (dE .LE. X_CLICK(jj,12) ) .OR. (jj .GE. Num_Hist) )	THEN
					NOW = 1
				END IF
			END DO
			IF (TRIM(Targ) .EQ. 'H') THEN
				C_SHA_H_E_H(ii,jj) = C_SHA_H_E_H(ii,jj) + 1.0D0/N_Part
			ELSE IF (TRIM(Targ) .EQ. 'He') THEN
				C_SHA_H_E_He(ii,jj) = C_SHA_H_E_He(ii,jj) + 1.0D0/N_Part
			ELSE IF (TRIM(Targ) .EQ. 'O') THEN
				C_SHA_H_E_O(ii,jj) = C_SHA_H_E_O(ii,jj) + 1.0D0/N_Part
			ELSE IF (TRIM(Targ) .EQ. 'Ar') THEN
				C_SHA_H_E_Ar(ii,jj) = C_SHA_H_E_Ar(ii,jj) + 1.0D0/N_Part
			ELSE IF (TRIM(Targ) .EQ. 'H2') THEN
				C_SHA_H_E_H2(ii,jj) = C_SHA_H_E_H2(ii,jj) + 1.0D0/N_Part
			ELSE IF (TRIM(Targ) .EQ. 'N2') THEN
				C_SHA_H_E_N2(ii,jj) = C_SHA_H_E_N2(ii,jj) + 1.0D0/N_Part
			ELSE IF (TRIM(Targ) .EQ. 'CO') THEN
				C_SHA_H_E_CO(ii,jj) = C_SHA_H_E_CO(ii,jj) + 1.0D0/N_Part
			ELSE IF (TRIM(Targ) .EQ. 'CO2') THEN
				C_SHA_H_E_CO2(ii,jj) = C_SHA_H_E_CO2(ii,jj) + 1.0D0/N_Part
			ELSE
				WRITE(*,*) 'UNKNOWN TARGET IN H vs E SHA, TARGET: ', Targ
			END IF
		END IF ! write_sha_H_E on	

		!! Run Secondary routines if SHA_ON = 1
		IF (SHA_ON .EQ. 1) THEN
			!! Find info on secondary hot particle
			CALL secondary_hot( E_now, E_nxt, E_targ, Theta_targ, Phi_targ )

			!! Write secondary hot info to file if hot enough
			IF (E_targ .GT. E_Therm) THEN
!				WRITE(*,*) 'SHA Energy ', E_targ
				CALL velocity_unit_vector( u_now, theta_prev, phi_prev, & 
				&													Theta_targ, Phi_targ, u_targ )
				Second_Count          = Second_Count  + 1
				Second_Energy         = Second_Energy + E_targ
				My_SHA_Count          = My_SHA_Count  + 1
				SHA_E(My_SHA_Count)   = E_targ
				SHA_R(My_SHA_Count,1) = x_now(1)
				SHA_R(My_SHA_Count,2) = x_now(2)
				SHA_R(My_SHA_Count,3) = x_now(3)
				SHA_V(My_SHA_Count,1) = u_targ(1)
				SHA_V(My_SHA_Count,2) = u_targ(2)
				SHA_V(My_SHA_Count,3) = u_targ(3)
				i = 1
				IF (My_SHA_Count .GE. MAX_SHA) THEN
					WRITE(*,*) 'SHA COUNT FOR ENA EXCEEDS MAX!!!'
				END IF
			END IF ! E_targ > E_therm
		END IF ! SHA_ON

		!! Convert ScattAng to lab frame
		CALL angle_to_lab( ScattAng, MP/MT, Lab_ScattAng )

!    WRITE(*,*) 'Collision with : ', Targ, ' theta_L : ', Lab_ScattAng, ' dE: ', dE

		!! Update target atom angle counters
		IF (TRIM(Targ) .EQ. 'CO2') THEN
			CO2_Ang = CO2_Ang + Lab_ScattAng*180.0D0/PI
		ELSE IF (TRIM(Targ) .EQ. 'O'  ) THEN
			O_Ang   = O_Ang   + Lab_ScattAng*180.0D0/PI
		ELSE IF (TRIM(Targ) .EQ. 'H'  ) THEN
			H_Ang   = H_Ang   + Lab_ScattAng*180.0D0/PI
		ELSE IF (TRIM(Targ) .EQ. 'He' ) THEN
			He_Ang  = He_Ang  + Lab_ScattAng*180.0D0/PI
		END IF

		!! Update average scattering angle counter
		Ave_Ang = Ave_Ang + Lab_ScattAng*180.0D0/PI

		!! Find random phi scattering angle
		r   = lfg()
		phi = r*2.0D0*PI

		!! Update new scattering angles
		theta_now = Lab_ScattAng
		phi_now   = phi

		!! Get new unit velocity and position vectors
		CALL particle_dir_cos_transport( u_now, x_now, theta_now, phi_now, Length, u_nxt, x_nxt )

		IF ( DO_GRAV .EQ. 1 ) THEN
			!! gravity current altitude
			CALL mars_gravity(x_now(1), gravity)

			WRITE(*,*) 'CO Zvel Gvel: ', velocity*u_now(1), gravity*dt
			!! update z component of velocity to include gravity  
			u_nxt(1) = velocity*u_nxt(1) - gravity*dt

			!! update unit velocity vector
			U_tot    = SQRT( (u_nxt(1))**2 + (velocity*u_nxt(2))**2 + (velocity*u_nxt(3))**2 )
			u_nxt(1) = u_nxt(1)/U_tot
			u_nxt(2) = u_nxt(2)/U_tot
			u_nxt(3) = u_nxt(3)/U_tot
		END IF

		!! Update values for next itteration
		theta_prev = theta_now
		phi_prev   = phi_now
		x_now(:)   = x_nxt(:)
		u_now(:)   = u_nxt(:)
		E_now      = E_nxt
		N_Coll     = N_Coll + 1
		GET_TCS    = 1

		E1  = E_nxt
		x1  = x_nxt(1)
		y1  = x_nxt(2)
		z1  = x_nxt(3)
		ux1 = u_nxt(1)
		uy1 = u_nxt(2)
		uz1 = u_nxt(3)
		t1  = theta_now
		p1  = phi_now

!		WRITE(*,*) 'COLL ux uy uz: ', ux1, uy1, uz1

	ELSE
		!!!!!!!!!!!!!!!!
		!! No Collision
		!!!!!!!!!!!!!!!!

!		WRITE(*,*) myid, ' NO COLLISION ', x_now(3)

		dE = 0.0D0

   	CALL energy_to_velocity( E_now, MP, velocity )
    dt = dL/velocity

		IF ( DO_GRAV .EQ. 1 ) THEN
			!! gravity current altitude
			CALL mars_gravity(x_now(1), gravity)

			WRITE(*,*) 'NC Zvel Gvel: ', velocity*u_now(1), gravity*dt
			!! update z component of velocity to include gravity	
			u_nxt(1) = velocity*u_now(1) - gravity*dt
			u_nxt(2) = velocity*u_now(2)
			u_nxt(2) = velocity*u_now(2)

			!! update unit velocity vector
			U_tot    = SQRT( u_nxt(1)**2 + u_nxt(2)**2 + u_nxt(3)**2 )
			u_nxt(1) = u_nxt(1)/U_tot
			u_nxt(2) = u_nxt(2)/U_tot
			u_nxt(3) = u_nxt(3)/U_tot
		ELSE
			U_tot    = SQRT( u_now(1)**2 + u_now(2)**2 + u_now(3)**2 )
			u_nxt(1) = u_now(1)/U_tot
			u_nxt(2) = u_now(2)/U_tot
			u_nxt(3) = u_now(3)/U_tot
		END IF

		!! Transport particle in straight line trajectory 
		x_nxt(1) = x_now(1) + dL*u_now(1)
		x_nxt(2) = x_now(2) + dL*u_now(2)
		x_nxt(3) = x_now(3) + dL*u_now(3)

		E1  = E0
		x1  = x_nxt(1)
		y1  = x_nxt(2)
		z1  = x_nxt(3)
		ux1 = u_nxt(1)
		uy1 = u_nxt(2)
		uz1 = u_nxt(3)
		t1  = t0
		p1  = p0

!		WRITE(*,*) 'NO COLL ux uy uz: ', ux1, uy1, uz1

		GET_TCS = 0
			
	END IF ! transport step

END SUBROUTINE planet_onestep_transport

