
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Return charge exchange cross sections for given
! input ion-atom/molecule collision which results in 
! energetic neutral atoms. 
!
! Parameters for H^+ + (H,O,N2,O2) given by 
! Lindsay and Stebbings, 2005, JGR, 110. Table 1. 
! Units of Fittings: [keV], [cm^2]
!
! Parameters for H^+   + He -> H    + He^+
!								 He^+  + H  -> He   + H^+
!								 He^+  + He -> He   + He^+
!                He^2+ + H  -> He^+ + H^+
!								 He^2+ + He -> He^+ + He^+
!                He^2+ + He -> He   + He^2+
! given by Barnett, CF et al., 1990, 
!          Atomic Data Fusion Vol. 1, ORNL. 
! Units of Fittings: [eV], [cm^2]
!
! H^+ + CO  -> H + CO^+ 
! H^+ + CO2 -> H + CO2^+ given by 
! Kusakabe, Buenker and Kimura, Charge Transfer Processes
! in Collisions of H+ ions with H2, D2, CO, CO2, CH4
! C2H2, H2H6 adn C3H8 below 10 keV
!
! He^++ + CO2 -> He + CO2^++ given by 
! from Greenwood, Chutjian, and Smith, 
! (2000), ApJ, 529:605-609
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE cx_cross_sections( ion, targ, ne, Et, CS ) 
	USE ENA, ONLY : L_ion, L_targ
	USE physics_constants, ONLY : Mass_p, Mass_n, Mass_e

	IMPLICIT NONE

	!! Inputs
	CHARACTER(Len=5)	:: ion					!! H_p, He_p, He_pp
	CHARACTER(Len=3)  :: targ					!! H, He, O, Ar, H2, N2, CO, CO2
	INTEGER						:: ne						!! Number of electrons transferred during collision
	REAL(KIND=8)			:: Et						!! Lab Frame [eV]

	!! Outputs
	REAL(KIND=8)			:: CS						!! [m^2]

	!! Internal	
	REAL(KIND=8)			:: amu, Emin, Emax, X, a0, a1, a2, a3, a4, a5, a6, a7, a8
	REAL(KIND=8)			:: E, C1, C2, C3, C4, C5, C6, C7, C8, C
	REAL(KIND=8)			:: T1, T2, T3, T4, T5, T6, T7, T8
	REAL(KIND=8)			:: ma, mb, mu, A, B, x1, y1, x2, y2, m, y0

	!! Convert cm^2 to m^2
	C = 1.0D0/10000.0D0

	!! H^+ projectile
	IF (TRIM(ion) .EQ. 'H_p') THEN

		!! H^+ + H2 -> H + H2^+
		IF (TRIM(targ) .EQ. 'H2') THEN
			E  = Et/1000.0D0
			y0 = 0.19D0; 
			m  = 3.7D0; 
			CS = 1.0D-16*(y0+m*E)
			CS = CS*C	
		END IF

		!! H^+ + Ar -> H + Ar^+
		IF (TRIM(targ) .EQ. 'Ar') THEN
			CS = 1.0D-15
			CS = CS*C
		END IF

		!! H^+ + CO -> H + CO^+
		IF (TRIM(targ) .EQ. 'CO') THEN
			E  = Et/1000.0D0
			x1 = 0.1D0
			y1 = 2.0D-15 
			x2 = 10.0D0
			y2 = 1.0D-15
			m  = (y2-y1)/(x2-x1)
			y0 = y1 - m*x1
			CS = y0 + m*E
			CS = CS*C
		END IF

		!! H^+ + N2 -> H + N2^+
		IF (TRIM(targ) .EQ. 'N2') THEN
			E  = Et/1000.0D0
			a1 = 12.5D0
			a2 = 1.52D0
			a3 = 3.97D0
			a4 = 0.36D0
			a5 = -1.20D0
			a6 = 0.208D0
			a7 = 0.741D0
			CS = a1*(EXP(-((LOG(E)-a2)**2)/a3))*(1.0D0-EXP(-E/a4))**2 &
			&    + ((a5-a6*LOG(E))**2)*(1.0D0-EXP(-a7/E))**2
			CS = CS*1e-16*C
		END IF
		
		!! H^+ + CO2 -> H + CO2^+
		IF (TRIM(targ) .EQ. 'CO2') THEN
			A  = -0.22544D0
			B  = -14.065D0	
			CS = 10.0D0**(A*LOG10(Et)+B)
			CS = CS*C
		END IF

		!! H^+ + H -> H + H^+
		IF (TRIM(targ) .EQ. 'H') THEN
			E  = Et/1000.0D0
			a1 = 4.15D0
			a2 = 0.531D0
			a3 = 67.3D0
			C1 = a1 - a2*LOG(E)
			C2 = 1.0D0 - EXP(-a3/E)
			C3 = C1*C1
			C4 = C2**4.5
			CS = C3*C4*1.0D-16
!WRITE(*,'(6ES10.2)') E, C1, C2, C3, C4, CS
			CS = CS*C
		END IF

		!! H^+ + O -> H + O^+
		IF (TRIM(targ) .EQ. 'O') THEN
			E  = Et/1000.0D0
			a1 = 2.91D0
			a2 = 0.0886D0
			a3 = 50.9D0
			a4 = 4.73D0
			a5 = -0.862D0
			a6 = 0.0306D0
			C1 = a1 - a2*LOG(E)
			C2 = 1.0D0 - EXP(-a3/E)
			C3 = C1*C1
			C4 = C2*C2
			C5 = a4 - a5*LOG(E)
			C6 = 1.0D0 - EXP(-a6/E)
			C7 = C6*C6
			CS = C3*C4 + C5*C7
			CS = CS*C*1.0D-16
		END IF

		!! H^+ + O2 -> H + O2^+
		IF (TRIM(targ) .EQ. 'O2') THEN
			E  = Et/1000.0D0
			a1 = 1.83D0
			a2 = -0.545D0
			a3 = 15.8D0
			a4 = 6.35D0
			a5 = -0.801D0
			a6 = 0.240D0
			C1 = a1 - a2*LOG(E)
			C2 = 1.0D0 - EXP(-a3/E)
			C3 = a4 - a5*LOG(E)
			C4 = 1.0D0 - EXP(-a6/E)
			C5 = C1*C1
			C6 = C2**1.5
			C7 = C3*C3
			CS = C5*C6 + C7*C4
			CS = CS*C*1.0D-16
		END IF

		!! H^+ + He -> H + He^+
		IF (TRIM(targ) .EQ. 'He') THEN
			amu  = Mass_p
			E    = Et/amu
			Emin = 9.9D1
			Emax = 1.1D7
			IF (E .LT. Emin) E = Emin
			IF (E .GT. Emax) E = Emax
			X    = ((LOG(E)-LOG(Emin))-(LOG(Emax)-LOG(E)))/(LOG(Emax)-LOG(Emin)) 
			a0   = -92.7819D0
			a1   = -4.69073D0
			a2   = -8.99666D0
			a3   = -0.490484D0
			a4   =  0.935276D0
			a5   = -0.0504853D0
			a6   = -0.0410750D0
			a7   = -0.0730136D0
			a8   =  0.0189320D0
			CALL chebyshev(1,X,T1)
			CALL chebyshev(2,X,T2)
			CALL chebyshev(3,X,T3)
			CALL chebyshev(4,X,T4)
			CALL chebyshev(5,X,T5)
			CALL chebyshev(6,X,T6)
			CALL chebyshev(7,X,T7)
			CALL chebyshev(8,X,T8)
			C1 = a0/2.0D0 + a1*T1+a2*T2+a3*T3+a4*T4+a5*T5+a6*T6+a7*T7+a8*T8
			CS = EXP(C1)
			CS = CS*C
!WRITE(*,'(8ES10.2)') T1, T2, T3, T4, T5, T6, T7, T8
!WRITE(*,'(5ES10.2)') C1, CS
		END IF
	
	!! He^+ projectile
	ELSE IF (TRIM(ion) .EQ. 'He_p') THEN

		!! He^+ + H -> He + H^+
		IF (TRIM(targ) .EQ. 'H') THEN
			amu  = 2.0D0*Mass_p + 2.0D0*Mass_n + Mass_e
			E    = Et/amu
			Emin = 2.5D2
			Emax = 7.0D5
!			IF (E .LT. Emin) E = Emin
!			IF (E .GT. Emax) E = Emax
			X    = ((LOG(E)-LOG(Emin))-(LOG(Emax)-LOG(E)))/(LOG(Emax)-LOG(Emin)) 
			a0   = -79.2449D0	
			a1   = -3.43445D0
			a2   = -3.38701D0
			a3   = -0.953151D0
			a4   =  0.0800233D0
			a5   =  0.153935D0	
			a6   = -0.198503D0
			a7   = -0.126958D0	
			a8   =  0.0492936D0
			CALL chebyshev(1,X,T1)
			CALL chebyshev(2,X,T2)
			CALL chebyshev(3,X,T3)
			CALL chebyshev(4,X,T4)
			CALL chebyshev(5,X,T5)
			CALL chebyshev(6,X,T6)
			CALL chebyshev(7,X,T7)
			CALL chebyshev(8,X,T8)
			C1 = a0/2.0D0 + a1*T1+a2*T2+a3*T3+a4*T4+a5*T5+a6*T6+a7*T7+a8*T8
			CS = EXP(C1)
			CS = CS*C
		END IF

		!! He^+ + He -> He + He^+
		IF (TRIM(targ) .EQ. 'He') THEN
			amu  = 2.0D0*Mass_p + 2.0D0*Mass_n + Mass_e
			E    = Et/amu
			Emin = 4.0D-2
			Emax = 5.3D5
			IF (E .LT. Emin) E = Emin
			IF (E .GT. Emax) E = Emax
			X    = ((LOG(E)-LOG(Emin))-(LOG(Emax)-LOG(E)))/(LOG(Emax)-LOG(Emin)) 
			a0   = -71.6915D0	
			a1   = -3.72702D0
			a2   = -1.97026D0
			a3   = -1.16040D0
			a4   = -0.605858D0 
			a5   = -0.259050D0	
			a6   = -0.0892322D0
			a7   = -0.0165873D0	
			a8   =  9.75783D-3 
			CALL chebyshev(1,X,T1)
			CALL chebyshev(2,X,T2)
			CALL chebyshev(3,X,T3)
			CALL chebyshev(4,X,T4)
			CALL chebyshev(5,X,T5)
			CALL chebyshev(6,X,T6)
			CALL chebyshev(7,X,T7)
			CALL chebyshev(8,X,T8)
			C1 = a0/2.0D0 + a1*T1+a2*T2+a3*T3+a4*T4+a5*T5+a6*T6+a7*T7+a8*T8
			CS = EXP(C1)
			CS = CS*C
		END IF	

	!! He^2+ projectile
	ELSE IF (TRIM(ion) .EQ. 'He_pp') THEN 

		!! He^2+ + CO2 -> He + CO2^2+
		IF (TRIM(targ) .EQ. 'CO2') THEN
			C1 = -6.747D-24
			C2 = 1.7632D-19	
			C3 = -2.496D-17
			CS = C1*Et*Et + C2*Et + C3
			CS = CS*C	
		END IF

		!! He^2+ + H -> He^+ + H^+
		IF (TRIM(targ) .EQ. 'H') THEN
			amu  = 2.0D0*Mass_p + 2.0D0*Mass_n
			E    = Et/amu
			Emin = 7.0D1
			Emax = 5.0D5
			IF (E .LT. Emin) E = Emin
			IF (E .GT. Emax) E = Emax
			X    = ((LOG(E)-LOG(Emin))-(LOG(Emax)-LOG(E)))/(LOG(Emax)-LOG(Emin)) 
			a0   = -78.3976D0
			a1   =  0.773499D0
			a2   = -5.12346D0
			a3   = -0.125575D0
			a4   = -0.232331D0
			a5   =  0.261317D0
			a6   = -0.0479652D0
			a7   =  0.189343D0
			a8   = -0.138708D0
			CALL chebyshev(1,X,T1)
			CALL chebyshev(2,X,T2)
			CALL chebyshev(3,X,T3)
			CALL chebyshev(4,X,T4)
			CALL chebyshev(5,X,T5)
			CALL chebyshev(6,X,T6)
			CALL chebyshev(7,X,T7)
			CALL chebyshev(8,X,T8)
			C1 = a0/2.0D0 + a1*T1+a2*T2+a3*T3+a4*T4+a5*T5+a6*T6+a7*T7+a8*T8
			CS = EXP(C1)
			CS = CS*C
		END IF

		!! He^2+ + He -> He^+ + He^+
		IF ( (TRIM(targ) .EQ. 'He') .AND. (ne .EQ. 1) ) THEN
			amu  = 2.0D0*Mass_p + 2.0D0*Mass_n
			E    = Et/amu
			Emin = 2.4D-4
			Emax = 2.0D6
			IF (E .LT. Emin) E = Emin
			IF (E .GT. Emax) E = Emax
			X    = ((LOG(E)-LOG(Emin))-(LOG(Emax)-LOG(E)))/(LOG(Emax)-LOG(Emin)) 
			a0   = -83.5447D0
			a1   =  0.758421D0
			a2   = -0.673803D0
			a3   = -3.68088D0
			a4   = -2.02719D0
			a5   = -0.148166D0
			a6   =  0.0237951D0
			a7   = -0.0738396D0
			a8   =  0.243993D0
			CALL chebyshev(1,X,T1)
			CALL chebyshev(2,X,T2)
			CALL chebyshev(3,X,T3)
			CALL chebyshev(4,X,T4)
			CALL chebyshev(5,X,T5)
			CALL chebyshev(6,X,T6)
			CALL chebyshev(7,X,T7)
			CALL chebyshev(8,X,T8)
			C1 = a0/2.0D0 + a1*T1+a2*T2+a3*T3+a4*T4+a5*T5+a6*T6+a7*T7+a8*T8
			CS = EXP(C1)
			CS = CS*C
		END IF

		!! He^2+ + He -> He + He^2+
		IF ( (TRIM(targ) .EQ. 'He') .AND. (ne .EQ. 2) ) THEN
			amu  = 2.0D0*Mass_p + 2.0D0*Mass_n
			E    = Et/amu
			Emin = 7.0D0 
			Emax = 4.0D5
			IF (E .LT. Emin) E = Emin
			IF (E .GT. Emax) E = Emax
			X    = ((LOG(E)-LOG(Emin))-(LOG(Emax)-LOG(E)))/(LOG(Emax)-LOG(Emin)) 
			a0   = -75.3410D0
			a1   = -3.70255D0
			a2   = -2.16467D0
			a3   = -1.23660D0
			a4   = -0.590825D0
			a5   = -0.179613D0
			a6   =  0.0396927D0
			a7   =  0.102768D0
			a8   =  0.0774058D0
			CALL chebyshev(1,X,T1)
			CALL chebyshev(2,X,T2)
			CALL chebyshev(3,X,T3)
			CALL chebyshev(4,X,T4)
			CALL chebyshev(5,X,T5)
			CALL chebyshev(6,X,T6)
			CALL chebyshev(7,X,T7)
			CALL chebyshev(8,X,T8)
			C1 = a0/2.0D0 + a1*T1+a2*T2+a3*T3+ &
      &    a4*T4+a5*T5+a6*T6+a7*T7+a8*T8
			CS = EXP(C1)
			CS = CS*C
		END IF

	END IF! projectiles

END SUBROUTINE cx_cross_sections

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE write_cx_cs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Write charge exchange cross sections
! to file when called for all species
!
! ALL CX CS are in [m^2] as is the 
! output from cx_cross_sections routine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE ENA, ONLY : L_ion, L_targ

  IMPLICIT NONE

  REAL(KIND=8)  :: E_in, E_fn, E, dE, E_H, E_He
  REAL(KIND=8)  :: Hp_CO2, Hp_O2, Hp_CO, Hp_O, Hp_H, Hp_He, Hep_H, Hp_N2
  REAL(KIND=8)  :: Hep_He, Hepp_CO2, Hepp_1H, Hepp_1He, Hepp_2He, Hp_H2, Hp_Ar
  INTEGER       :: NE, i
  CHARACTER     :: ion*10, targ*10

  OPEN(UNIT=55,FILE="../Data/H_CX_CS.dat",ACCESS="APPEND")
  OPEN(UNIT=56,FILE="../Data/He_CX_CS.dat",ACCESS="APPEND")

  E_in = 1.0D2    !! [eV]
  E_fn = 2.0D3    !! [eV]
  NE   = 50     !! Number of energies
  dE   = (E_fn-E_in)/REAL(NE-1)

  DO i=1,NE
    E_H    = E_in + (i-1)*dE
    E_He   = E_H*4.0D0
    ion    = 'H_p'
    L_ion  = LEN_TRIM(ion)
    targ   = 'H2'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_H2)
    targ   = 'Ar'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_Ar)
    targ   = 'CO'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_CO)
    targ   = 'CO2'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_CO2)
    targ   = 'O2'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_O2)
    targ   = 'N2'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_N2)
    targ   = 'O'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_O)
    targ   = 'He'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_He)
    targ   = 'H'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_H,Hp_H)
    ion    = 'He_p'
    L_ion  = LEN_TRIM(ion)
    targ   = 'H'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_He,Hep_H)
    targ   = 'He'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_He,Hep_He)
    ion    = 'He_pp'
    L_ion  = LEN_TRIM(ion)
    targ   = 'CO2'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,2,E_He,Hepp_CO2)
    targ   = 'H'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_He,Hepp_1H)
    targ   = 'He'
    L_targ = LEN_TRIM(targ)
    CALL cx_cross_sections(ion,targ,1,E_He,Hepp_1He)
    CALL cx_cross_sections(ion,targ,2,E_He,Hepp_2He)

    WRITE(55,'(10ES11.3)') E_H, Hp_H, Hp_He, Hp_O, Hp_Ar, Hp_H2, Hp_N2, Hp_CO, Hp_O2, Hp_CO2
    WRITE(56,'(7ES11.3)') E_He, Hep_H, Hep_He, Hepp_1H, Hepp_1He, Hepp_2He, Hepp_CO2

  END DO

  CLOSE(55)
  CLOSE(56)

END SUBROUTINE write_cx_cs


