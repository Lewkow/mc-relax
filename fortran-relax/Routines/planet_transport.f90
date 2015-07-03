!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! Transport particle in planetary atmosphere
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE planet_transport( E0, x0, y0, z0, ux0, uy0, uz0, SHA_ON ) 
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
	INTEGER							:: SHA_ON						!! SHA tracking on if = 1 off if = 0
																					!! if SHA = 2 write data to SHA arrays
																					!! instead of ENA arrays

	!! Outputs


	!! Internal
	REAL(KIND=8),DIMENSION(3)	:: x_now, u_now, u_targ, u_nxt, x_nxt
	REAL(KIND=8),DIMENSION(4)	:: Prob_Targ
	REAL(KIND=8)							:: lfg
	REAL(KIND=8)							:: theta_now, phi_now, elec_den, HpHTCS
	REAL(KIND=8)							:: theta_prev, phi_prev, E_now, E_nxt, E_targ
	REAL(KIND=8)							:: HCO2_TCS, HeH_TCS, HeO_TCS, HeHe_TCS
	REAL(KIND=8)							:: Den_CO2, Den_H, Den_O, Den_He, Den_Tot
	REAL(KIND=8)							:: MFP_HCO2, MFP_HeH, MFP_HeO, MFP_HeHe, MFP
	REAL(KIND=8)							:: MR_CO2, MR_He, MR_H, MR_O, phi
	REAL(KIND=8)							:: P_Tot, rand_targ, Theta_targ, Phi_targ
	REAL(KIND=8)							:: C1, r, dL, Length, ScattAng, Lab_ScattAng
	REAL(KIND=8)							:: dE, U_tot, XXX(5)
	INTEGER										:: N_Coll, Transport, GET_TCS	
	INTEGER										:: i, j, xyz, click
	
	!! Set steplength for transport [m]
	dL          = 1.0D3	

	E_now       = E0
	click				= 0 
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

	theta_prev 	= ACOS(uz0)
	phi_prev = ACOS(ux0/SIN(theta_prev))

	!! Init collision counter to zero
	N_Coll = 0

	!! Set transport/TCS to ON
	Transport = 1
	GET_TCS   = 1

	DO WHILE (Transport .EQ. 1)
		!! Find total cross section [m^2] if energy has changed 
		IF ( GET_TCS .EQ. 1) THEN
			IF (CM_TYPE .EQ. 'QM') THEN	
				CALL get_qm_tcs( E_now, HCO2_TCS, HeH_TCS, HeO_TCS, HeHe_TCS )
			ELSE IF (CM_TYPE .EQ. 'HS') THEN
				CALL get_hs_tcs( HCO2_TCS, HeH_TCS, HeO_TCS, HeHe_TCS )	
			END IF ! HS or QM TCS
		END IF ! if GET_TCS

		CALL kras_mars_density( 44, x_now(1), Den_CO2 )
		CALL kras_mars_density( 1,  x_now(1), Den_H   )
		CALL kras_mars_density( 16, x_now(1), Den_O   )
		CALL kras_mars_density( 4,  x_now(1), Den_He  )

!		CALL get_electron_density( x_now(1), elec_den )	
!		CALL Hp_H_tcs( E_now, HpHTCS )

		MFP_HCO2 = Den_CO2*HCO2_TCS
		MFP_HeH  = Den_H  *HeH_TCS
		MFP_HeO  = Den_O  *HeO_TCS
		MFP_HeHe = Den_He *HeHe_TCS

		Den_Tot = Den_CO2 + Den_H + Den_O + Den_He

		!! Total mean free path for all target species
		MFP = 1.0D0/( MFP_HCO2 + MFP_HeH + MFP_HeO + MFP_HeHe )
		C1  = EXP( -dL/MFP )
		r   = lfg()		

		!! Determine if collision occurs in current steplength  
		IF ( r .GT. C1 ) THEN
			!!-------------
			!! Collision
			!!-------------
			!! Find exact location of collision
			!! CHECK THIS FOR CORRECTNESS
			!! EITHER MFP OR DL
!			Length = -MFP*LOG(r)
			Length = -dL*LOG(r)

			!! Find mixing ratios for target species at collision altitude
			MR_CO2 = Den_CO2/Den_Tot
			MR_He  = Den_He /Den_Tot
			MR_H   = Den_H  /Den_Tot
			MR_O   = Den_O  /Den_Tot

			!! Find probability array for target species collisions 
			P_Tot        = MR_CO2*HCO2_TCS + MR_He*HeHe_TCS + MR_H*HeH_TCS + MR_O*HeO_TCS
			Prob_Targ(1) = MR_CO2*HCO2_TCS/P_Tot
			Prob_Targ(2) = Prob_Targ(1) + MR_He*HeHe_TCS/P_Tot
			Prob_Targ(3) = Prob_Targ(2) + MR_H *HeH_TCS/P_Tot
			Prob_Targ(4) = Prob_Targ(3) + MR_O *HeO_TCS/P_Tot

			!! Find random target atom
			rand_targ = lfg()
			IF (rand_targ .LT. Prob_Targ(1)) THEN
				!! Target is CO2
				Targ      = 'CO2'
				CO2_Count = CO2_Count + 1
			ELSE IF ( (rand_targ .GT. Prob_Targ(1)) .AND. (rand_targ .LT. Prob_Targ(2)) ) THEN
				!! Target is He
				Targ     = 'He'
				He_Count = He_Count + 1
			ELSE IF ( (rand_targ .GT. Prob_Targ(2)) .AND. (rand_targ .LT. Prob_Targ(3)) ) THEN
				!! Target is H
				Targ    = 'H'
				H_Count = H_Count + 1
			ELSE IF ( rand_targ .GT. Prob_Targ(3) ) THEN
				!! Target is O
				Targ    = 'O'
				O_Count = O_Count + 1
			END IF

			!! Find mass of target atom
			CALL mass_finder

			!! Find reduced mass [kg]
			mass = MU*TOAMU*AMUTOKG

			!! Update collision counter
			Coll_Count = Coll_Count + 1

			!! Find random scattering angle
			IF (Targ .EQ. 'CO2') THEN
				atom = 2
			ELSE
				atom = 1
			END IF

			!!!!!!!!!!!
			!! Get random scattering angle in CM frame for given collision
			!!!!!!!!!!!
			!! use 1 keV for all energies > 1keV until tables complete
			!!!!!!!!!!!
			IF (CM_TYPE .EQ. 'QM') THEN
				IF (Targ .EQ. 'CO2') CALL lin_rand_angle( E_now, 'HeO ', ScattAng )
				IF (Targ .EQ. 'H'  ) CALL lin_rand_angle( E_now, 'HeH ', ScattAng )
				IF (Targ .EQ. 'O'  ) CALL lin_rand_angle( E_now, 'HeO ', ScattAng )
				IF (Targ .EQ. 'He' ) CALL lin_rand_angle( E_now, 'HeHe', ScattAng )
			ELSE IF (CM_TYPE .EQ. 'HS') THEN
				CALL HS_rand_angle( ScattAng )
			END IF

			!! Find new projectile energy
			CALL find_new_energy( E_now, ScattAng, MP/MT, E_nxt )
			dE = E_now - E_nxt

			!! Run Secondary routines if SHA_ON = 1
			IF (SHA_ON .EQ. 1) THEN
				!! Find info on secondary hot particle
				CALL secondary_hot( E_now, E_nxt, E_targ, Theta_targ, Phi_targ )

				!! Write secondary hot info to file if hot enough
				IF (E_targ .GT. E_Therm) THEN
!					WRITE(*,*) 'SHA Energy ', E_targ
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
      		DO WHILE ( (x_now(1) .GE. r_zone(i)) .AND. (i .LT. SIZE(r_zone(:))) )
        		i = i + 1
      		END DO
					r_prod_SHA_E(i) = r_prod_SHA_E(i) + E_now	
					r_prod_SHA_C(i) = r_prod_SHA_C(i) + 1.0D0 	
					IF (My_SHA_Count .GE. MAX_SHA) THEN
						WRITE(*,*) 'SHA COUNT FOR ENA EXCEEDS MAX!!!'
						EXIT
					END IF
				END IF ! E_targ > E_therm
			END IF ! SHA_ON

			!! Convert ScattAng to lab frame
			CALL angle_to_lab( ScattAng, MP/MT, Lab_ScattAng )

			!! Update target atom angle counters
			IF (Targ .EQ. 'CO2') THEN
				CO2_Ang = CO2_Ang + Lab_ScattAng*180.0D0/PI
			ELSE IF (Targ .EQ. 'O'  ) THEN
				O_Ang   = O_Ang   + Lab_ScattAng*180.0D0/PI
			ELSE IF (Targ .EQ. 'H'  ) THEN
				H_Ang   = H_Ang   + Lab_ScattAng*180.0D0/PI
			ELSE IF (Targ .EQ. 'He' ) THEN
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
			CALL particle_transport( u_now, x_now, theta_prev, phi_prev, &
			&                        theta_now, phi_now, Length, u_nxt, x_nxt )

			!! Update values for next itteration
			theta_prev = theta_now
			phi_prev   = phi_now
			x_now(:)   = x_nxt(:)
			u_now(:)   = u_nxt(:)
			E_now      = E_nxt
			N_Coll     = N_Coll + 1
			GET_TCS    = 1

!			IF (SHA_ON .EQ. 2) WRITE(*,*) '  COLL z: ', x_now(1)/1.0D3

		ELSE
			!!!!!!!!!!!!!!!!
			!! No Collision
			!!!!!!!!!!!!!!!!

			dE = 0.0D0

			!! Transport particle in straight line trajectory 
			x_now(1) = x_now(1) + dL*u_now(1)
			x_now(2) = x_now(2) + dL*u_now(2)
			x_now(3) = x_now(3) + dL*u_now(3)

			GET_TCS = 0
		
!			IF (SHA_ON .EQ. 2) WRITE(*,*) '  NO COLL z: ', x_now(1)/1.0D3
			
		END IF ! transport step

!		IF (SHA_ON .EQ. 2) WRITE(*,*) '     E: z: ', E_now, x_now(1)/1.0D3

		!! Exit if any of the following critera are met
		IF ( x_now(1) .LE. 0.0D0) THEN
			Planet_Energy = Planet_Energy + E_now
			Planet_Count  = Planet_Count  + 1
			Transport 		= 0
		END IF
	
		IF ( (x_now(1) .GE. high) .AND. (E_now .GE. E_esc) ) THEN
			Escape_Energy = Escape_Energy + E_now
			Escape_Count  = Escape_Count  + 1
			Transport 		= 0
		END IF

		IF ( E_now .LE. E_Therm)  THEN
			Therm_Count = Therm_Count + 1
			Transport 	= 0
!			WRITE(420,*) x_now(1)
		END IF

		IF ( x_now(1) .GE. 2.0D0*high ) THEN
			Transport = 0
		END IF

		!! Update individual ranks arrays
		IF ( x_now(1) .LE. high ) THEN
			i = 1
			DO WHILE ( x_now(1) .GE. r_zone(i) )
				i = i + 1
			END DO
			IF (SHA_ON .EQ. 2) THEN ! update SHA arrays
				r_zone_SHA_E(i)  = r_zone_SHA_E(i)  + E_now
				r_zone_SHA_C(i)  = r_zone_SHA_C(i)  + 1.0D0
				r_zone_SHA_dE(i) = r_zone_SHA_dE(i) + dE
				r_zone_SHA_dC(i) = r_zone_SHA_dC(i) + 1.0D0
				j = 1
				DO WHILE ((E_now .GE. SHA_E_zone(j)) .AND. (j .LT. SIZE(SHA_E_zone(:))))
					j = j+1
				END DO
				My_SHA_Engy_Dist(i,j) = My_SHA_Engy_Dist(i,j) + 1.0D0/REAL(N_Part)
			ELSE !! update ENA arrays
				r_zone_E(i)  = r_zone_E(i)  + E_now
				r_zone_C(i)  = r_zone_C(i)  + 1.0D0
				r_zone_dE(i) = r_zone_dE(i) + dE
				r_zone_dC(i) = r_zone_dC(i) + 1.0D0
				j = 1
				DO WHILE ((E_now .GE. E_zone(j)) .AND. (j .LT. SIZE(E_zone(:))))
					j = j+1
				END DO
				My_Engy_Dist(i,j)     = My_Engy_Dist(i,j) + 1.0D0/REAL(N_Part)
			END IF ! SHA or ENA
		END IF ! height less than high
	click = click + 1
	END DO ! Transport = 1

	WRITE(456,*) click

END SUBROUTINE planet_transport

