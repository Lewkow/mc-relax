
SUBROUTINE ion_planet_simulation
	USE mpi_info
	USE tables
	USE rand_seed
	USE physics_constants
	USE ENA
	
	IMPLICIT NONE

	INTEGER,PARAMETER	:: N_Targs = 5
	REAL(KIND=8)			:: z0, zn, E0, En, Ex, mfp, HH_el_tcs, HH_cx_tcs, H_den, tot_den, tot_tcs
	REAL(KIND=8)			:: HH_el_mfp,HH_cx_mfp,C,He_den,O_den,CO2_den,He_tcs,O_tcs,CO2_tcs
	REAL(KIND=8)			:: Ar_den, H2_den
	REAL(KIND=8)			:: Prob_Targ(N_Targs), P_Tot, HH_MR, e_den, dEe, dL, H_mfp, HpH_mfp
	REAL(KIND=8)			:: lt, r, r_targ, Length, zx, He_el_tcs, O_el_tcs, CO2_el_tcs
	REAL(KIND=8)			:: He_cx_tcs, O_cx_tcs, CO2_cx_tcs, CO2_MR, HE_MR, O_MR
	REAL(KIND=8)			:: dummy(4)
	REAL(KIND=8)			:: lfg
	INTEGER						:: N, ion, i
	CHARACTER					:: ionn*10, targ*10

	dummy(:) = 0.0D0

	CALL read_e_den_table
	CALL read_HpH_tcs_table
	CALL read_Vsw_table

	OPEN(UNIT=101,FILE='../Data/ENA_Prod_Height.dat', ACCESS='APPEND')

	!! Set number of projectiles
	N = 100000	

	WRITE(*,*) 'Working...'
	WRITE(*,*) '[          ]'

	!! start loop over number number of particles
	DO i=1,N

		!! Set initial height
		z0 = 800.0D3		! [m]
		zn = z0

		!! Set initial energy
		CALL rand_init_energy( 1.0D0, E0 )
!		E0 = 3000.0D0		! [eV]
		En = E0

		!! start as ion
		ion = 1

		!! total displacement
		lt = 0.0D0	

		!! Start while loop while projectile is ion
		DO WHILE ( ion .EQ. 1 )

			!! find density of neutrals at current height [1/m^3]
			CALL kras_mars_density( 1, zn, H_den )
			CALL kras_mars_density( 4, zn, He_den )
			CALL kras_mars_density( 16, zn, O_den )
			CALL kras_mars_density( 16, zn, Ar_den )
			CALL kras_mars_density( 44, zn, CO2_den )

			!! get electron density at current height [1/m^3]
			CALL get_electron_density( zn, e_den )

			!! find total density of neutrals at current height
			tot_den = H_den + He_den + O_den + CO2_den

			!! find total cross section
			CALL Hp_H_tcs( En, HH_el_tcs )
			CALL TCS_HeH( dummy(:), En, He_el_tcs )
			
			!! convert tcs from [a0^2] -> [m^2]
			He_el_tcs = He_el_tcs*BOHRTOM**2

			ionn   = 'H_p'
			L_ion  = LEN_TRIM(ionn)
			targ   = 'H'
			L_targ = LEN_TRIM(targ)
			CALL cx_cross_sections( TRIM(ionn), TRIM(targ), 1, En, HH_cx_tcs ) ! [m^2]
			targ   = 'He'
			L_targ = LEN_TRIM(targ)
			CALL cx_cross_sections( TRIM(ionn), TRIM(targ), 1, En, He_cx_tcs ) ! [m^2]
			targ   = 'O'
			L_targ = LEN_TRIM(targ)
			CALL cx_cross_sections( TRIM(ionn), TRIM(targ), 1, En, O_cx_tcs ) ! [m^2]
			targ   = 'CO2'
			L_targ = LEN_TRIM(targ)
			CALL cx_cross_sections( TRIM(ionn), TRIM(targ), 1, En, CO2_cx_tcs ) ! [m^2]

			!! get mean free paths	[1/m]	
			O_el_tcs   = HH_el_tcs
			CO2_el_tcs = HH_el_tcs
			HH_el_mfp  = H_den*HH_el_tcs + He_den*He_el_tcs + O_den*O_el_tcs + CO2_den*CO2_el_tcs
			HH_cx_mfp  = H_den*HH_cx_tcs + He_den*He_cx_tcs + O_den*O_cx_tcs + CO2_den*CO2_cx_tcs

!			WRITE(*,*) H_den, HH_el_tcs, HH_cx_tcs
!			WRITE(*,*) 'MFP: EL | CX ', HH_el_mfp, HH_cx_mfp

			!! total mean free path [m]
			mfp = 1.0D0/( HH_el_mfp + HH_cx_mfp )

!			WRITE(*,*) 'MFP EL: CX: TOT: ', HH_el_mfp, HH_cx_mfp, mfp

			!! set step length as 0.25*mfp [m]
!			dL = mfp/4.0D0
			dL = 1.0D3

			!! find probability for colliding within steplength
			C  = EXP( -dL/mfp )

			!! get random number to decide if collision occurred
			r  = lfg()

!			WRITE(*,*) 'mfp: C: r: ', mfp, C, r

			IF ( r .GT. C ) THEN
				!!!!!!!!!!!!	
				!! Collision
				!!!!!!!!!!!!	
		
				!! exact location of collision [m]	
				Length = -mfp*LOG(r)

				!! get height now
				zx     = zn - Length
!				WRITE(*,*) 'Collision: zx: zn: ', zx/1.0d3, zn/1.0d3		

				!! get mixing ratios for target species at collision height
				HH_MR  = H_den / tot_den
				He_MR  = He_den / tot_den
				O_MR   = O_den / tot_den
				CO2_MR = CO2_den / tot_den

				!! get total collision probability
				P_Tot  				= HH_MR*HH_el_tcs + HH_MR*HH_cx_tcs
				Prob_Targ(1) 	= HH_MR*HH_cx_tcs/P_Tot
				Prob_Targ(2)  = Prob_Targ(1) + HH_MR*HH_el_tcs/P_TOT
				Prob_Targ(3)  = Prob_Targ(2) + He_MR*He_el_tcs/P_TOT
				Prob_Targ(4)  = Prob_Targ(3) + O_MR*O_el_tcs/P_TOT
				Prob_Targ(5)  = Prob_Targ(4) + CO2_MR*CO2_el_tcs/P_TOT

				!! get random number for random target
				r_targ = lfg()
				
				IF ( r_targ .LT. Prob_Targ(1) ) THEN
					!! Target is H and collision is CX

					!! get energy loss due to electron drag
					CALL butler_drag(1, Length, e_den, 300.0D0, En, Ex)
					ion = 0
					lt  = lt + Length
					zx  = zn - Length
					Ex  = Ex - 0.5d0
					WRITE(101,*) zn, En		
!					WRITE(*,*) 'H CX ', zx, En 
				ELSE IF ( (r_targ .GT. Prob_Targ(1)) .AND. (r_targ .LT. Prob_Targ(2)) ) THEN
					!! get energy loss due to electron drag
					CALL butler_drag(1, Length, e_den, 300.0D0, En, Ex)
					Ex = Ex - 0.5d0
					zx = zn - Length
					lt = lt + Length			
!					WRITE(*,*) 'H EL ', zx, Ex 
				ELSE IF ( (r_targ .GT. Prob_Targ(2)) .AND. (r_targ .LT. Prob_Targ(3)) ) THEN
					!! Target is He and collision is CX
					CALL butler_drag(1, Length, e_den, 300.0D0, En, Ex)
					ion = 0
					lt  = lt + Length
					zx  = zn - Length
					Ex  = Ex - 0.5d0
					WRITE(101,*) zn, En		
!					WRITE(*,*) 'He CX ', zx, En 
				ELSE IF ( (r_targ .GT. Prob_Targ(3)) .AND. (r_targ .LT. Prob_Targ(4)) ) THEN
					!! Target is O and collision is CX
					CALL butler_drag(1, Length, e_den, 300.0D0, En, Ex)
					ion = 0
					lt  = lt + Length
					zx  = zn - Length
					Ex  = Ex - 0.5d0
					WRITE(101,*) zn, En		
!					WRITE(*,*) 'O CX ', zx, En 
				ELSE IF ( (r_targ .GT. Prob_Targ(4)) .AND. (r_targ .LT. Prob_Targ(5)) ) THEN
					!! Target is CO2 and collision is CX
					CALL butler_drag(1, Length, e_den, 300.0D0, En, Ex)
					ion = 0
					lt  = lt + Length
					zx  = zn - Length
					Ex  = Ex - 0.5d0
					WRITE(101,*) zn, En		
!					WRITE(*,*) 'CO2 CX ', zx, En 
				END IF 

			ELSE
				!!!!!!!!!!!
				!! No Collision	
				!!!!!!!!!!!

				!! get energy loss due to electron drag
				CALL butler_drag(1, dL, e_den, 300.0D0, En, Ex)

				!! update total length
				lt = lt + dL 
				zx = zn - dL

!				WRITE(*,*) 'No Collision zx: zn: dE: ', zx/1.0d3, zn/1.0d3, En-Ex

			END IF

			zn = zx
			En = Ex

			IF ( (zn .LT. 0.0D0) .OR. (En .LT. 1.0D-3) ) THEN
				zn  = 0.0D0
				En  = 0.0D0
				ion = 0
			END IF

		END DO ! while ion

		IF (100.0d0*i/REAL(N) .EQ. 10 ) WRITE(*,*) '[*         ]'
		IF (100.0d0*i/REAL(N) .EQ. 20 ) WRITE(*,*) '[**        ]'
		IF (100.0d0*i/REAL(N) .EQ. 30 ) WRITE(*,*) '[***       ]'
		IF (100.0d0*i/REAL(N) .EQ. 40 ) WRITE(*,*) '[****      ]'
		IF (100.0d0*i/REAL(N) .EQ. 50 ) WRITE(*,*) '[*****     ]'
		IF (100.0d0*i/REAL(N) .EQ. 60 ) WRITE(*,*) '[******    ]'
		IF (100.0d0*i/REAL(N) .EQ. 70 ) WRITE(*,*) '[*******   ]'
		IF (100.0d0*i/REAL(N) .EQ. 80 ) WRITE(*,*) '[********  ]'
		IF (100.0d0*i/REAL(N) .EQ. 90 ) WRITE(*,*) '[********* ]'
		IF (100.0d0*i/REAL(N) .EQ. 100 )WRITE(*,*) '[**********]'



	END DO ! i


CLOSE(101)


END SUBROUTINE ion_planet_simulation

