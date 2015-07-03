
SUBROUTINE ion_setup( MC_Start, MC_End )
	USE planet
	USE physics_constants, ONLY : PI
	
	IMPLICIT NONE

	!! Inputs
	INTEGER				:: MC_Start, MC_End

	!! Implicit
	INTEGER				:: i, NOW, CX, NC, Coll, TOO_HIGH
	REAL(KIND=8)	:: E0, E00, E1, x0, y0, z0, x1, y1, z1
	REAL(KIND=8)	:: ux0, uy0, uz0, ux1, uy1, uz1, t0, t1
	REAL(KIND=8)	:: p0, p1, dt, r_phi, phi, theta, MassP
	REAL(KIND=8)	:: lfg

	ALLOCATE( Nas_ENA_Energy(N_Part) )	
	ALLOCATE( Nas_ENA_Height(N_Part) )	

	Nas_ENA_Energy = 0.0D0
	Nas_ENA_Height = 0.0D0

	IF (PROJ .EQ. 'H ') MassP = 1.0D0
	IF (PROJ .EQ. 'He') MassP = 4.0D0

	DO i=MC_Start,MC_End	
		NC = 0
		IF (MONO_ENGY .EQ. 1) THEN
			E0 = E_0
		ELSE
			CALL rand_init_energy( MassP, E00 )
!			E00 = 2000.0D0
			WRITE(500,*) E00
		END IF	
		r_phi = lfg()
		phi		= r_phi*2.0D0*PI
		theta = PI - SZA
		ux0   = SIN(theta)*COS(phi)
		uy0   = SIN(theta)*SIN(phi)
		uz0   = COS(theta)
		x0    = 0.0D0	
		y0    = 0.0D0	
		z0    = 800.0D3 ! [m]
		E0    = E00
		NOW   = 1	
		TOO_HIGH = 0
		DO WHILE (NOW .EQ. 1) 
			CALL planet_onestep_ion_transport( E0, x0, y0, z0, ux0, uy0, uz0, t0, p0, &
			&                E1, x1, y1, z1, ux1, uy1, uz1, t1, p1, dt, Coll, CX )
			E0  = E1
			x0  = x1
			y0  = y1
			z0  = z1
			ux0 = ux1
			uy0 = uy1
			uz0 = uz1
			IF (CX .EQ. 1) NOW = 0	
			IF (Coll .EQ. 1) NC  = NC + 1
			IF (z1 .GT. 1000.0D3) THEN
				TOO_HIGH 	= 1
				NOW				= 0
			END IF	
		END DO
		IF (TOO_HIGH .EQ. 0) THEN	
			Nas_ENA_Energy(i) = E0
			WRITE(600,*) E0
			Nas_ENA_Height(i) = z0
			WRITE(700,*) z0
			WRITE(800,*) NC
			WRITE(*,'(A,I7,3ES12.2,I7)') 'MC: E0: Ef: Z0: NC: ', i, E00, E0, z0/1000.0D0, NC
		END IF
	END DO ! i

END SUBROUTINE ion_setup

!###################
!###################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! Transport particle in planetary atmosphere
!! 
!! Step-by-step methods used in this transport 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE planet_onestep_ion_transport(E0,x0,y0,z0,ux0,uy0,uz0,t0,p0,E1,x1,y1,z1,ux1,uy1,uz1,t1,p1,dt,Coll,CX) 
	USE planet
	USE physics_constants
	USE secondary
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
	INTEGER							:: CX								!! if a CX collision occurs CX = 1, else CX = 0
	INTEGER							:: Coll							!! if an elastic collision occurs, Coll = 1

	!! Internal
	REAL(KIND=8),DIMENSION(3)	:: x_now, u_now, u_targ, u_nxt, x_nxt, dum_x
	REAL(KIND=8),DIMENSION(16):: Prob_Targ
	REAL(KIND=8)							:: lfg
	REAL(KIND=8)							:: theta_now, phi_now, elec_den, HpHTCS
	REAL(KIND=8)							:: theta_prev, phi_prev, E_now, E_nxt, E_targ
	REAL(KIND=8)							:: HCO2_TCS, HeH_TCS, HeO_TCS, HeHe_TCS
	REAL(KIND=8)							:: CX_X_CO2, CX_X_CO, CX_X_H2, CX_X_N2 
	REAL(KIND=8)							:: CX_X_O, CX_X_He, CX_X_H, CX_X_Ar 
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
	REAL(KIND=8)							:: pp9, pp10, pp11, pp12, pp13, pp14, pp15, pp16
	INTEGER										:: N_Coll, Transport, GET_TCS	
	INTEGER										:: i, j, xyz, click, NOW, ii, jj
	
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

	!! Find total cross section [m^2] 
	CALL get_cx_tcs(E_now,CX_X_CO2,CX_X_CO,CX_X_N2,CX_X_H2,CX_X_O,CX_X_Ar,CX_X_He,CX_X_H)
	CALL get_ion_atom_tcs(E_now,TCS_X_CO2,TCS_X_CO,TCS_X_N2,TCS_X_H2,TCS_X_O,TCS_X_Ar,TCS_X_He,TCS_X_H)

	CALL mars_table_density( 'H  ', x_now(3), Den_H   )
	CALL mars_table_density( 'He ', x_now(3), Den_He  )
	CALL mars_table_density( 'O  ', x_now(3), Den_O   )
	CALL mars_table_density( 'Ar ', x_now(3), Den_Ar  )
	CALL mars_table_density( 'H2 ', x_now(3), Den_H2  )
	CALL mars_table_density( 'N2 ', x_now(3), Den_N2  )
	CALL mars_table_density( 'CO ', x_now(3), Den_CO  )
	CALL mars_table_density( 'CO2', x_now(3), Den_CO2 )

	MFP_X_CO2 = Den_CO2 *(TCS_X_CO2 + CX_X_CO2)
	MFP_X_CO  = Den_CO  *(TCS_X_CO  + CX_X_CO )
	MFP_X_H2  = Den_H2  *(TCS_X_H2  + CX_X_H2 )
	MFP_X_N2  = Den_N2  *(TCS_X_N2  + CX_X_N2 )
	MFP_X_O   = Den_O   *(TCS_X_O   + CX_X_O  )
	MFP_X_He  = Den_He  *(TCS_X_He  + CX_X_He )
	MFP_X_Ar  = Den_Ar  *(TCS_X_Ar  + CX_X_Ar )
	MFP_X_H   = Den_H   *(TCS_X_H   + CX_X_H  )
	
!	WRITE(40,*) x_now(3)/1000.0D0, Den_H, Den_He, Den_O, Den_Ar, Den_H2, Den_N2, Den_CO, Den_CO2
!	WRITE(41,*) x_now(3)/1000.0D0, TCS_X_H, TCS_X_He, TCS_X_O, TCS_X_Ar, TCS_X_H2, TCS_X_N2, TCS_X_CO, TCS_X_CO2
!	WRITE(42,*) x_now(3)/1000.0D0, CX_X_H, CX_X_He, CX_X_O, CX_X_Ar, CX_X_H2, CX_X_N2, CX_X_CO, CX_X_CO2

!	WRITE(*,*) 'Z:   ', x_now(3)/1000.0D0
!	WRITE(*,*) 'MOL: ', MFP_X_CO2, MFP_X_CO, MFP_X_H2, MFP_X_N2
!	WRITE(*,*) 'ATOM:', MFP_X_O, MFP_X_He, MFP_X_Ar, MFP_X_H

	Den_Tot = Den_CO2 + Den_H + Den_O + Den_He + Den_CO + Den_H2 + Den_N2 + Den_Ar

	!! Total mean free path for all target species
	MFP = 1.0D0/(MFP_X_CO2+MFP_X_CO+MFP_X_N2+MFP_X_H2+MFP_X_O+MFP_X_Ar+MFP_X_He+MFP_X_H)
!	WRITE(*,*) '----------------'
!	WRITE(*,*) 'E: z: ', E0, z0/1000.0D0
!	WRITE(*,'(9ES10.2)') MFP, MFP_X_H, MFP_X_He, MFP_X_O, MFP_X_Ar, MFP_X_H2, MFP_X_N2, MFP_X_CO, MFP_X_CO2

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
		pp9 		= MR_CO2*CX_X_CO2 
		pp10 		= MR_CO*CX_X_CO 
		pp11 		= MR_H2*CX_X_H2 
		pp12 		= MR_N2*CX_X_N2 
		pp13 		= MR_O*CX_X_O 
		pp14 		= MR_He*CX_X_He 
		pp15 		= MR_Ar*CX_X_Ar 
		pp16 		= MR_H*CX_X_H
		P_Tot 	= pp1+pp2+pp3+pp4+pp5+pp6+pp7+pp8 &
		&        +pp9+pp10+pp11+pp12+pp13+pp14+pp15+pp16
		Prob_Targ(1)  = pp1/P_Tot
		Prob_Targ(2)  = Prob_Targ(1)  + pp2/P_Tot
		Prob_Targ(3)  = Prob_Targ(2)  + pp3/P_Tot
		Prob_Targ(4)  = Prob_Targ(3)  + pp4/P_Tot
		Prob_Targ(5)  = Prob_Targ(4)  + pp5/P_Tot
		Prob_Targ(6)  = Prob_Targ(5)  + pp6/P_Tot
		Prob_Targ(7)  = Prob_Targ(6)  + pp7/P_Tot
		Prob_Targ(8)  = Prob_Targ(7)  + pp8/P_Tot
		Prob_Targ(9)  = Prob_Targ(8)  + pp9/P_Tot
		Prob_Targ(10) = Prob_Targ(9)  + pp10/P_Tot
		Prob_Targ(11) = Prob_Targ(10) + pp11/P_Tot
		Prob_Targ(12) = Prob_Targ(11) + pp12/P_Tot
		Prob_Targ(13) = Prob_Targ(12) + pp13/P_Tot
		Prob_Targ(14) = Prob_Targ(13) + pp14/P_Tot
		Prob_Targ(15) = Prob_Targ(14) + pp15/P_Tot
		Prob_Targ(16) = Prob_Targ(15) + pp16/P_Tot

		WRITE(*,*)
		WRITE(*,*) 'Z: ', z0/1000.0D0
		WRITE(*,*) 'EL: ', Prob_Targ(8)
		WRITE(*,*) 'CX: ', 1.0D0-Prob_Targ(8)
		WRITE(*,*)

		!! Find random target atom
		rand_targ = lfg()
		IF (rand_targ .LT. Prob_Targ(1)) THEN
			!! Target is CO2
			Targ = 'CO2'
			atom = 2
			CX   = 0
			Coll = 1
		ELSE IF ( (rand_targ .GT. Prob_Targ(1)) .AND. (rand_targ .LT. Prob_Targ(2)) ) THEN
			!! Target is CO
			Targ = 'CO '
			atom = 2
			CX   = 0
			Coll = 1
		ELSE IF ( (rand_targ .GT. Prob_Targ(2)) .AND. (rand_targ .LT. Prob_Targ(3)) ) THEN
			!! Target is H2
			Targ = 'H2 '
			atom = 2
			CX   = 0
			Coll = 1
		ELSE IF ( (rand_targ .GT. Prob_Targ(3)) .AND. (rand_targ .LT. Prob_Targ(4)) ) THEN
			!! Target is N2 
			Targ = 'N2 '
			atom = 2
			CX   = 0
			Coll = 1
		ELSE IF ( (rand_targ .GT. Prob_Targ(4)) .AND. (rand_targ .LT. Prob_Targ(5)) ) THEN
			!! Target is O
			Targ = 'O  '
			atom = 1 
			CX   = 0
			Coll = 1
		ELSE IF ( (rand_targ .GT. Prob_Targ(5)) .AND. (rand_targ .LT. Prob_Targ(6)) ) THEN
			!! Target is He 
			Targ = 'He '
			atom = 1 
			CX   = 0
			Coll = 1
		ELSE IF ( (rand_targ .GT. Prob_Targ(6)) .AND. (rand_targ .LT. Prob_Targ(7)) ) THEN
			!! Target is Ar 
			Targ = 'Ar '
			atom = 1 
			CX   = 0
			Coll = 1
		ELSE IF ( (rand_targ .GT. Prob_Targ(7)) .AND. (rand_targ .LT. Prob_Targ(8)) ) THEN
			!! Target is H 
			Targ = 'H  '
			atom = 1 
			CX   = 0
			Coll = 1
		ELSE IF ( (rand_targ .GT. Prob_Targ(8)) .AND. (rand_targ .LT. Prob_Targ(9)) ) THEN
			!! Target is CO2 CX 
			Targ = 'CO2'
			atom = 1 
			CX   = 1
			Coll = 0
		ELSE IF ( (rand_targ .GT. Prob_Targ(9)) .AND. (rand_targ .LT. Prob_Targ(10)) ) THEN
			!! Target is CO CX
			Targ = 'CO '
			atom = 1 
			CX   = 1
			Coll = 0
		ELSE IF ( (rand_targ .GT. Prob_Targ(10)) .AND. (rand_targ .LT. Prob_Targ(11)) ) THEN
			!! Target is H2 CX
			Targ = 'H2 '
			atom = 1 
			CX   = 1
			Coll = 0
		ELSE IF ( (rand_targ .GT. Prob_Targ(11)) .AND. (rand_targ .LT. Prob_Targ(12)) ) THEN
			!! Target is N2
			Targ = 'N2 '
			atom = 1 
			CX   = 1
			Coll = 0
		ELSE IF ( (rand_targ .GT. Prob_Targ(12)) .AND. (rand_targ .LT. Prob_Targ(13)) ) THEN
			!! Target is O CX
			Targ = 'O  '
			atom = 1 
			CX   = 1 
			Coll = 0
		ELSE IF ( (rand_targ .GT. Prob_Targ(13)) .AND. (rand_targ .LT. Prob_Targ(14)) ) THEN
			!! Target is He CX 
			Targ = 'He '
			atom = 1 
			CX   = 1 
			Coll = 0
		ELSE IF ( (rand_targ .GT. Prob_Targ(14)) .AND. (rand_targ .LT. Prob_Targ(15)) ) THEN
			!! Target is Ar CX
			Targ = 'Ar '
			atom = 1 
			CX   = 1 
			Coll = 0
		ELSE IF ( (rand_targ .GT. Prob_Targ(15)) .AND. (rand_targ .LT. Prob_Targ(16)) ) THEN
			!! Target is H CX
			Targ = 'H  '
			atom = 1 
			CX   = 1
			Coll = 0
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
					CALL lin_rand_angle( E_now, 'HeH  ', ScattAng )
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
		WRITE(*,*) 'Theta: Enow: Enxt: dE: ', ScattAng*PI/180.0D0, E_now, E_nxt, dE

		!! Convert ScattAng to lab frame
		CALL angle_to_lab( ScattAng, MP/MT, Lab_ScattAng )

		!! Find random phi scattering angle
		r   = lfg()
		phi = r*2.0D0*PI

		!! Update new scattering angles
		theta_now = Lab_ScattAng
		phi_now   = phi

		!! Get new unit velocity and position vectors
		CALL particle_dir_cos_transport( u_now, x_now, theta_now, phi_now, Length, u_nxt, x_nxt )

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

	ELSE
		!!!!!!!!!!!!!!!!
		!! No Collision
		!!!!!!!!!!!!!!!!
		CX   = 0
		Coll = 0
		dE   = 0.0D0

   	CALL energy_to_velocity( E_now, MP, velocity )
    dt = dL/velocity

		U_tot    = SQRT( u_now(1)**2 + u_now(2)**2 + u_now(3)**2 )
		u_nxt(1) = u_now(1)/U_tot
		u_nxt(2) = u_now(2)/U_tot
		u_nxt(3) = u_now(3)/U_tot

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

	END IF ! transport step

END SUBROUTINE planet_onestep_ion_transport

