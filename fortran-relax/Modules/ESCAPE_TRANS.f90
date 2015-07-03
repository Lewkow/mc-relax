
MODULE escape_trans
	USE planet

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: SH_Prod_Alt
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: SH_Prod_Rate
	INTEGER																:: N_SH_Prod

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: SH_ED_Energy
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: SH_ED_Prob
	INTEGER																:: N_SH_ED

CONTAINS

!!!!!!!!

SUBROUTINE get_SH_Energy_Prob( E, prob )
	IMPLICIT NONE

	REAL(KIND=8)	:: E, prob

	INTEGER 			:: i, keep_going
	REAL(KIND=8)	:: x0, x1, y0, y1, Enow

	keep_going = 1
	i = 2 
	prob = 0.0D0
	DO WHILE (keep_going .EQ. 1)
		Enow = SH_ED_Energy(i)
!		WRITE(*,*) 'i: znow: ztable: ptable: ', i, znow, SH_Prod_Alt(i), SH_Prod_Rate(i)
		IF ((E .LT. Enow) .OR. ((i+1) .GE. SIZE(SH_ED_Energy(:)))) THEN
			x0 					= SH_ED_Energy(i-1)
			x1 					= Enow	
			y0 					= SH_ED_Prob(i-1)
			y1 					= SH_ED_Prob(i)
			prob 				= y0 + ( (E-x0)*y1 - (E-x0)*y0 )/(x1-x0)
			keep_going 	= 0
		ELSE IF (E .GT. SH_ED_Energy(N_SH_ED)) THEN
			prob = 0.0D0
			keep_going = 0
		ELSE
			i = i+1
		END IF
	END DO	

END SUBROUTINE get_SH_Energy_Prob

!!!!!!!!

SUBROUTINE test_SH_Energy_Prob

	IMPLICIT NONE

	REAL(KIND=8)	:: E, dE, P, Ei, Ef
	INTEGER				:: i, N

	N  = 1000	
	Ei = 0.0D0
	Ef = 4.0D0
	dE = (Ef-Ei)/REAL(N-1)

	DO i=1,N
		E = Ei + REAL(i-1)*dE
		CALL get_SH_Energy_Prob(E,P)
	END DO

END SUBROUTINE test_SH_Energy_Prob

!!!!!!!!

SUBROUTINE rand_Nascent_SH_Energy( E0 )

	IMPLICIT NONE

	! Outputs
	REAL(KIND=8)	:: E0	! [eV]

	! Internal
	REAL(KIND=8)	:: lfg
	REAL(KIND=8)	:: randy, dE, E_now, Tot, P_now, TTot, Ei, Ef
	INTEGER				:: i, NE

	randy = lfg()
	Ei    = 0.0D0
	Ef    = 4.0D0
	NE		= 1000
	dE    = (Ef-Ei)/REAL(NE-1)
	Tot   = 0.0D0

	TTot = 0.0D0
	DO i=1,NE
		E_now = Ei + REAL(i-1)*dE
		CALL get_SH_Energy_Prob( E_now, P_now )
		TTot = TTot + P_now
	END DO

	i = 0	
	DO WHILE (Tot .LT. randy)
		E_now = Ei + REAL(i-1)*dE
		CALL get_SH_Energy_Prob( E_now, P_now )
		Tot 	= Tot + P_now/TTot
		i 		= i+1
	END DO

	E0 = E_now		

END SUBROUTINE rand_Nascent_SH_Energy

!!!!!!!!

SUBROUTINE get_Alt_Yield( z, yield )

	IMPLICIT NONE

	! inputs
	REAL(KIND=8)	:: z 		 ! [km]

	! outputs
	REAL(KIND=8)	:: yield ! [1/cm]

	INTEGER				:: i, keep_going
	REAL(KIND=8)	:: znow, dz, km_to_cm

	km_to_cm    = 1.0D5
	keep_going 	= 1
	i 					= 1
	yield 			= 0.0D0
	dz    			= SH_Prod_Alt(2) - SH_Prod_Alt(1)
	dz					= dz*km_to_cm
	
	DO WHILE (keep_going .EQ. 1)
		znow = SH_Prod_Alt(i)
		IF ( (i+1) .GE. SIZE(SH_PROD_Alt(:))) THEN
			yield = 0.0D0
			keep_going = 0
		ELSE IF (z .LT. znow) THEN
			yield = SH_Prod_Rate(i)/dz	
			keep_going = 0
		ELSE
			i = i+1
		END IF
	END DO	

END SUBROUTINE get_Alt_Yield

!!!!!!!!

SUBROUTINE trans_integral( E, z, MSH, T )

	IMPLICIT NONE

	! inputs
	REAL(KIND=8)	:: E   ! projectile energy [eV]
	REAL(KIND=8)	:: z   ! initial altitude  [cm]
	REAL(KIND=8)	:: MSH ! projectile mass   [amu]

	! outputs
	REAL(KIND=8)	:: T ! []

	! internal
	REAL(KIND=8)	:: dz, z_esc, km_to_cm, cm_to_m, a0_to_cm, znow, C1, C2, tot
	REAL(KIND=8)	:: m3_to_cm3
	REAL(KIND=8)	:: n_H, n_He, n_O, n_Ar, n_H2, n_N2, n_CO, n_CO2
	REAL(KIND=8)	:: c_H, c_He, c_O, c_Ar, c_H2, c_N2, c_CO, c_CO2
	INTEGER				:: N_int, i

	! conversion factors
	km_to_cm  = 1.0D5
	cm_to_m   = 1.0D-2
	a0_to_cm  = 5.29D-9
	m3_to_cm3 = 1.0D-6

	! number of integrations steps
	N_int = 100

	! escape altitude [cm] 
	z_esc = 700.0D0*km_to_cm

	! height integration chunk [cm]
	dz = (z_esc - z)/REAL(N_int-1)

	tot = 0.0D0

	DO i=1,N_int
		znow = z + REAL(i-1)*dz
	
		! get densities at current altitude
		CALL mars_table_density('H  ', znow*cm_to_m, n_H  )	
		CALL mars_table_density('He ', znow*cm_to_m, n_He )	
		CALL mars_table_density('O  ', znow*cm_to_m, n_O  )	
		CALL mars_table_density('Ar ', znow*cm_to_m, n_Ar )	
		CALL mars_table_density('H2 ', znow*cm_to_m, n_H2 )	
		CALL mars_table_density('N2 ', znow*cm_to_m, n_N2 )	
		CALL mars_table_density('CO ', znow*cm_to_m, n_CO )	
		CALL mars_table_density('CO2', znow*cm_to_m, n_CO2)	

		n_H = n_H * m3_to_cm3
		n_He = n_He * m3_to_cm3
		n_O = n_O * m3_to_cm3
		n_Ar = n_Ar * m3_to_cm3
		n_H2 = n_H2 * m3_to_cm3
		n_N2 = n_N2 * m3_to_cm3
		n_CO = n_CO * m3_to_cm3
		n_CO2 = n_CO2 * m3_to_cm3

		! get diffusion cross sections for current energy
		CALL universal_dfcs_int( E, MSH, M_H,   1, c_H  )
		CALL universal_dfcs_int( E, MSH, M_He,  1, c_He )
		CALL universal_dfcs_int( E, MSH, M_O,   1, c_O  )
		CALL universal_dfcs_int( E, MSH, M_Ar,  1, c_Ar )
		CALL universal_dfcs_int( E, MSH, M_H2,  1, c_H2 )
		CALL universal_dfcs_int( E, MSH, M_N2,  1, c_N2 )
		CALL universal_dfcs_int( E, MSH, M_CO,  1, c_CO )
		CALL universal_dfcs_int( E, MSH, M_CO2, 1, c_CO2)

		c_H   = c_H  * a0_to_cm**2
		c_He  = c_He * a0_to_cm**2
		c_O   = c_O  * a0_to_cm**2
		c_Ar  = c_Ar * a0_to_cm**2
		c_H2  = c_H2 * a0_to_cm**2
		c_N2  = c_N2 * a0_to_cm**2
		c_CO  = c_CO * a0_to_cm**2
		c_CO2 = c_CO2* a0_to_cm**2

!		WRITE(66,*) n_H, n_He, n_O, n_Ar, n_H2, n_N2, n_CO, n_CO2
!		WRITE(77,*) c_H, c_He, c_O, c_Ar, c_H2, c_N2, c_CO, C_CO2

		C1 	= n_H*c_H+n_He*c_He+n_O*c_O+n_Ar*c_Ar+n_H2*c_H2+n_N2*c_N2+n_CO*c_CO+n_CO2*c_CO2
		C2 	= C1*dz
		tot = tot + C2
!		WRITE(*,*) 'z: dz: C1: C2: t: ', znow, dz, C1, C2, tot

	END DO

	T = tot

END SUBROUTINE trans_integral

!!!!!!!!

SUBROUTINE read_SH_Prod_Alt_Table( SH, ENA )

	IMPLICIT NONE

	! inputs
	CHARACTER(LEN=3)	:: SH		! SH atom/molecule to find escape trans
	CHARACTER(LEN=2)	:: ENA	! ENA which creates SH atom/molecule

	INTEGER						:: i

	IF (SH .EQ. 'H  ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_H.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_H.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'He ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_He.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_He.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'O  ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_O.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_O.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'Ar ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_Ar.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_Ar.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'H2 ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_H2.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_H2.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'N2 ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_N2.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_N2.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'CO ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_CO.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_CO.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'CO2') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_CO2.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_CO2.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE
		WRITE(*,*) 'SH is UNKNOWN!'
	END IF

	READ(88,*) N_SH_Prod
	ALLOCATE( SH_Prod_Alt( N_SH_Prod ), SH_Prod_Rate( N_SH_Prod ) )
	
	WRITE(*,*) 'Production Altitudes for SH: ENA: ', SH, ENA	
	DO i=1,N_SH_Prod
		READ(88,*) SH_Prod_Rate(i), SH_Prod_Alt(i)
		WRITE(*,*) SH_Prod_Alt(i), SH_Prod_Rate(i)
	END DO

	CLOSE(88)

END SUBROUTINE read_SH_Prod_Alt_Table 

!!!!!!!!

SUBROUTINE read_SH_ED_Table( SH, ENA )

	IMPLICIT NONE

	! inputs
	CHARACTER(LEN=3)	:: SH		! SH atom/molecule to find escape trans
	CHARACTER(LEN=2)	:: ENA	! ENA which creates SH atom/molecule

	INTEGER						:: i

	IF (SH .EQ. 'H  ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_ED_H.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_ED_H.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'He ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_ED_He.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_ED_He.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'O  ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_ED_O.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_ED_O.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'Ar ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_ED_Ar.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_ED_Ar.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'H2 ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_ED_H2.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_ED_H2.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'N2 ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_ED_N2.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_ED_N2.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'CO ') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_ED_CO.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_ED_CO.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE IF (SH .EQ. 'CO2') THEN
		IF ( ENA .EQ. 'H ') THEN
			OPEN(UNIT=88, FILE="../Tables/H_Mean_SH_ED_CO2.dat", STATUS="old", ACTION="read")
		ELSE
			OPEN(UNIT=88, FILE="../Tables/He_Mean_SH_ED_CO2.dat", STATUS="old", ACTION="read")
		END IF	
	ELSE
		WRITE(*,*) 'SH is UNKNOWN!'
	END IF

	READ(88,*) N_SH_ED	
	ALLOCATE( SH_ED_Energy( N_SH_ED ), SH_ED_Prob( N_SH_ED ) )

	WRITE(*,*) 'Energy Dist for SH: ENA: ', SH, ENA	
	DO i=1,N_SH_ED
		READ(88,*) SH_ED_Energy(i), SH_ED_Prob(i)
		WRITE(*,*) SH_ED_Energy(i), SH_ED_Prob(i)
	END DO

	CLOSE(88)

END SUBROUTINE read_SH_ED_Table

!!!!!!!!

SUBROUTINE clean_SH_ED_Table

	DEALLOCATE(SH_ED_Energy)
	DEALLOCATE(SH_ED_Prob)

END SUBROUTINE clean_SH_ED_Table

!!!!!!!!

SUBROUTINE clean_SH_Prod_Alt_Table

	DEALLOCATE(SH_Prod_Alt)
	DEALLOCATE(SH_Prod_Rate)

END SUBROUTINE clean_SH_Prod_Alt_Table

!!!!!!!!

SUBROUTINE all_trans

	CHARACTER(LEN=2)		:: ENA
	CHARACTER(LEN=3)		:: SH

	INTEGER							:: ALL_ON

	ALL_ON = 0

	! determine which atmosphere model to use
	! 0 = min SA; 1 = mean SA; 2 = max SA
	ATMOSPHERE = 1

	! set up atmosphere density table
	CALL read_mars_density

	IF (ALL_ON .EQ. 0) THEN
	ENA = 'H '
	SH  = 'H  '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table
	ENA = 'He'
	SH  = 'H  '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table
	ELSE
	ENA = 'H '
	SH  = 'He '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'H '
	SH  = 'O  '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'H '
	SH  = 'H2 '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'H '
	SH  = 'N2 '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'H '
	SH  = 'CO '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'He'
	SH  = 'He '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'He'
	SH  = 'O  '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'He'
	SH  = 'H2 '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'He'
	SH  = 'N2 '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	ENA = 'He'
	SH  = 'CO '
	CALL read_SH_ED_Table(SH,ENA)
	CALL read_SH_Prod_Alt_Table(SH,ENA)
	CALL trans(SH,ENA)
	CALL clean_SH_ED_Table
	CALL clean_SH_Prod_Alt_Table

	END IF

	! clean up density tables
	CALL clean_mars_density

END SUBROUTINE all_trans

!!!!!!!!

SUBROUTINE trans( SH, ENA )

	USE planet

	IMPLICIT NONE

	! inputs
	CHARACTER(LEN=3)	:: SH		! SH atom/molecule to find escape trans
	CHARACTER(LEN=2)	:: ENA	! ENA which creates SH atom/molecule

	! internal
	REAL(KIND=8)			:: SH_mass, Esc_E, Esc_P, E, prob, z, yield, tot, T, alt_tot
	REAL(KIND=8)			:: km_to_cm, z_min, z_max, E_max, dE, dz, E_now, zm_now, zi_now, z_now
	REAL(KIND=8)			:: u_now, u_tot, TT, E_tot, du
	INTEGER						:: i, N_u, N_z, N_E, jz, iz, iE, iu, YIELD_TEST, OLD_WAY, j

	WRITE(*,*) 'DOING ESCAPE TRANS --> ENA: SHA: ', ENA, SH

	km_to_cm = 1.0D5

	IF (SH .EQ. 'H  ') THEN
		SH_mass =	M_H
		Esc_E   = 0.1D0
		IF ( ENA .EQ. 'H ') THEN
			Esc_P = 1.0D0
			OPEN(UNIT=555, FILE="../Data/H_ENA_SH_H_Escape.dat", ACCESS="APPEND")
		ELSE
			Esc_P = 1.0D0
			OPEN(UNIT=555, FILE="../Data/He_ENA_SH_H_Escape.dat", ACCESS="APPEND")
		END IF	
	ELSE IF (SH .EQ. 'He ') THEN
		SH_mass =	M_He
		Esc_E   = 0.4D0
		IF ( ENA .EQ. 'H ') THEN
			Esc_P = 0.801802D0
			OPEN(UNIT=555, FILE="../Data/H_ENA_SH_He_Escape.dat", ACCESS="APPEND")
		ELSE
			Esc_P = 0.9375D0
			OPEN(UNIT=555, FILE="../Data/He_ENA_SH_He_Escape.dat", ACCESS="APPEND")
		END IF	
	ELSE IF (SH .EQ. 'O  ') THEN
		SH_mass =	M_O
		Esc_E   = 1.7D0
		IF ( ENA .EQ. 'H ') THEN
			Esc_P = 0.01944D0
			OPEN(UNIT=555, FILE="../Data/H_ENA_SH_O_Escape.dat", ACCESS="APPEND")
		ELSE
			Esc_P = 0.033415D0
			OPEN(UNIT=555, FILE="../Data/He_ENA_SH_O_Escape.dat", ACCESS="APPEND")
		END IF	
	ELSE IF (SH .EQ. 'Ar ') THEN
		SH_mass =	M_Ar
		Esc_E   = 4.3D0
		IF ( ENA .EQ. 'H ') THEN
			Esc_P = 0.0D0
			OPEN(UNIT=555, FILE="../Data/H_ENA_SH_Ar_Escape.dat", ACCESS="APPEND")
		ELSE
			Esc_P = 0.0D0
			OPEN(UNIT=555, FILE="../Data/He_ENA_SH_Ar_Escape.dat", ACCESS="APPEND")
		END IF	
	ELSE IF (SH .EQ. 'H2 ') THEN
		SH_mass =	M_H2
		Esc_E   = 0.2D0
		IF ( ENA .EQ. 'H ') THEN
			Esc_P = 0.935484D0
			OPEN(UNIT=555, FILE="../Data/H_ENA_SH_H2_Escape.dat", ACCESS="APPEND")
		ELSE
			Esc_P = 1.0D0
			OPEN(UNIT=555, FILE="../Data/He_ENA_SH_H2_Escape.dat", ACCESS="APPEND")
		END IF	
	ELSE IF (SH .EQ. 'N2 ') THEN
		SH_mass =	M_N2
		Esc_E   = 3.0D0
		IF ( ENA .EQ. 'H ') THEN
			Esc_P = 0.000131D0
			OPEN(UNIT=555, FILE="../Data/H_ENA_SH_N2_Escape.dat", ACCESS="APPEND")
		ELSE
			Esc_P = 0.000215D0
			OPEN(UNIT=555, FILE="../Data/He_ENA_SH_N2_Escape.dat", ACCESS="APPEND")
		END IF	
	ELSE IF (SH .EQ. 'CO ') THEN
		SH_mass =	M_CO
		Esc_E   = 3.0D0
		IF ( ENA .EQ. 'H ') THEN
			Esc_P = 0.000594D0
			OPEN(UNIT=555, FILE="../Data/H_ENA_SH_CO_Escape.dat", ACCESS="APPEND")
		ELSE
			Esc_P = 0.000576D0
			OPEN(UNIT=555, FILE="../Data/He_ENA_SH_CO_Escape.dat", ACCESS="APPEND")
		END IF	
	ELSE IF (SH .EQ. 'CO2') THEN
		SH_mass =	M_CO2
		Esc_E   = 4.8D0
		IF ( ENA .EQ. 'H ') THEN
			Esc_P = 0.0D0
			OPEN(UNIT=555, FILE="../Data/H_ENA_SH_CO2_Escape.dat", ACCESS="APPEND")
		ELSE
			Esc_P = 0.0D0
			OPEN(UNIT=555, FILE="../Data/He_ENA_SH_CO2_Escape.dat", ACCESS="APPEND")
		END IF	
	ELSE
		WRITE(*,*) 'SH is UNKNOWN!'
	END IF

	! minimum altitude [cm]
	z_min = 70.0D0*km_to_cm

	! maximum altitude [cm]
	z_max = 700.0D0*km_to_cm

	! maximum energy [eV]
	E_max = 5.0D0

	! number of energy/altitude integrations
	N_z   = 100
	N_E		= 100
	N_u   = 100

	! differential energy [eV]
	dE = (E_max-Esc_E)/REAL(N_E-1)

	! differential altitude [cm]
	dz = (z_max-z_min)/REAL(N_z-1)

	tot = 0.0D0	

	! OLD_WAY turns on original, 1D transparency
	! OLD_WAY off calculates, upward escape velocity cone transparency
	OLD_WAY = 0

	IF (OLD_WAY .EQ. 1) THEN
		DO iz=1,N_z
			z_now = z_min + REAL(iz-1)*dz
			CALL get_Alt_Yield( z_now/km_to_cm, yield )
			alt_tot = 0.0D0
			IF (yield .GT. 0.0D0) THEN
				DO iE=1,N_E
					E_now = Esc_E + REAL(iE-1)*dE
					CALL get_SH_Energy_Prob( E_now, prob )
					CALL trans_integral( E_now, z_now, SH_mass, T )
					TT = EXP(-T)
					WRITE(*,'(5ES10.2)') z_now/km_to_cm, E_now, yield, prob, TT
					alt_tot = alt_tot + yield*prob*TT*dE*dz	
					tot 		= tot + yield*prob*TT*dE*dz	
				END DO
			END IF
!			WRITE(69,*) z_now, alt_tot
			WRITE(555,*) z_now/km_to_cm, alt_tot
			WRITE(*,*) 'z: Esc/Inc: ', z_now, alt_tot
		END DO
	ELSE
		tot = 0.0D0
		DO jz=1,N_z
			zm_now = z_min + REAL(jz-1)*dz
			CALL get_Alt_Yield( zm_now/km_to_cm, yield )
			alt_tot = 0.0D0
			IF (yield .GT. 0.0D0) THEN
				E_tot = 0.0D0
				DO iE=1,N_E
					E_now = Esc_E + REAL(iE-1)*dE
					CALL get_SH_Energy_Prob( E_now, prob )
					du 		= (1.0D0-SQRT(Esc_E/E_now))/REAL(N_u-1)				
					u_tot = 0.0D0
					CALL trans_integral( E_now, zm_now, SH_mass, T )
					DO iu=1,N_u
						u_now = SQRT(Esc_E/E_now)+REAL(iu-1)*du
						TT 		= EXP((-1.0D0/u_now)*T)
						u_tot = u_tot + TT*du	
						WRITE(*,*) 'T: TT: u_tot: ', T, TT, u_tot
					END DO ! iu		
					E_tot = E_tot + dE*prob*0.5D0*u_tot
					WRITE(*,*) 'E_tot: u_tot: ', E_tot, u_tot
				END DO ! iE	
				alt_tot = alt_tot + dz*yield*E_tot
				tot 		= tot     + dz*yield*E_tot
			END IF ! yield > 0
			WRITE(555,*) zm_now/km_to_cm, alt_tot
			WRITE(*,*) 'z: Esc/Inc: ', zm_now/km_to_cm, alt_tot
		END DO ! jz
	END IF

	WRITE(*,*) 'Total Esc/Inc = ', tot

!	DO i=1,100
!		z = 70.0D0 + REAL(i-1)*(300.0D0-70.0D0)/100.0D0
!		CALL trans_integral( 3.0D0, z*km_to_cm, SH_mass, T )
!		WRITE(*,*) 'z: T: ', z, T	
!	END DO

	YIELD_TEST = 0
	IF (YIELD_TEST .EQ. 1) THEN
		tot = 0.0D0
		DO i=1,100
			z = 70.0D0 + REAL(i-1)*(250.0D0-70.0D0)/100.0D0
			CALL get_Alt_Yield( z, yield )
			tot = tot + yield	
			WRITE(*,*) 'z: yield: tot: ', z, yield, tot
		END DO
	END IF

!	DO i=1,100
!		E = REAL(i-1)*5.0D0/100.0D0
!		CALL get_SH_Energy_Prob( E, prob )
!		WRITE(120,*) E, prob
!	END DO

	CLOSE(555)

END SUBROUTINE trans




END MODULE escape_trans


