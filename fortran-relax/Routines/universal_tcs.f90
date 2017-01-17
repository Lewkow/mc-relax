SUBROUTINE test_pd

	IMPLICIT NONE

	REAL(KIND=8)	:: Mp, Mt, E
	INTEGER		:: atom

	Mp = 4.0d0
	Mt = 44.0d0
	E  = 500.0d0
	
	atom = 2

	CALL universal_pd( E, Mp, Mt, atom )

END SUBROUTINE test_pd

!#########################################################################
!#########################################################################

SUBROUTINE test_universal_tcs

	IMPLICIT NONE

	REAL(KIND=8)	:: TCS, E_in, E_fn, E, dE, U_TCS, U_DFCS, Q_TCS, Mp, Mt, U_EL
	INTEGER		:: NE, i, atom

	Mp   = 4.0D0
	Mt   = 1.0D0

	atom = 1

	NE   = 100
	E_in = 10.0D0
	E_fn = 1000.0D0
	dE   = (E_fn-E_in)/REAL(NE-1)

	DO i=1,NE
		E = E_in + (i-1)*dE
		CALL universal_tcs_int( E, Mp, Mt, atom, TCS )
		WRITE(22,*) E, TCS
	END DO

END SUBROUTINE test_universal_tcs

!#########################################################################
!#########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Obtain universal tcs from tables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE universal_table_tcs( E, Targ, TCS )
	USE planet, ONLY : PROJ
	USE universal, ONLY : TCS_ENERGY, TCS_X_O, TCS_X_Ar, TCS_X_H2, TCS_X_N2, TCS_X_CO, TCS_X_CO2, NL_TCS

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E		! Center of Mass frame collision energy [eV]
	CHARACTER(LEN=3)	:: Targ ! Target atom

	!! Outputs
	REAL(KIND=8)		:: TCS	! Total cross section [a0^2]

	!! Internal
	REAL(KIND=8)		:: TCS_A(NL_TCS)
	REAL(KIND=8)		:: x1, x2, y1, y2, y0, m
	INTEGER			:: i, z

	IF (PROJ .EQ. 'H ') THEN
		IF (TRIM(Targ) .EQ. 'O')   TCS_A = TCS_X_O		
		IF (TRIM(Targ) .EQ. 'Ar')  TCS_A = TCS_X_Ar		
		IF (TRIM(Targ) .EQ. 'H2')  TCS_A = TCS_X_H2		
		IF (TRIM(Targ) .EQ. 'N2')  TCS_A = TCS_X_N2		
		IF (TRIM(Targ) .EQ. 'CO')  TCS_A = TCS_X_CO		
		IF (TRIM(Targ) .EQ. 'CO2') TCS_A = TCS_X_CO2		
	ELSE IF (PROJ .EQ. 'He') THEN
		IF (TRIM(Targ) .EQ. 'Ar')  TCS_A = TCS_X_Ar		
		IF (TRIM(Targ) .EQ. 'H2')  TCS_A = TCS_X_H2		
		IF (TRIM(Targ) .EQ. 'N2')  TCS_A = TCS_X_N2		
		IF (TRIM(Targ) .EQ. 'CO')  TCS_A = TCS_X_CO		
		IF (TRIM(Targ) .EQ. 'CO2') TCS_A = TCS_X_CO2		
	ELSE 
		WRITE(*,*) 'Projectile ', PROJ, ' not know from universal table tcs'	
	END IF

	z = 0
	i = 0
	DO WHILE (z .EQ. 0) 
		i = i + 1
		IF ( (TCS_ENERGY(i) .GE. E) .OR. (i .EQ. SIZE(TCS_ENERGY)) ) THEN
			z = 1
			IF (i .EQ. 1) i = 2
		END IF
	END DO	
	x1 = TCS_ENERGY(i-1)
	x2 = TCS_ENERGY(i)
	y1 = TCS_A(i-1)
	y2 = TCS_A(i)
	m  = (y2-y1)/(x2-x1)
	y0 = y1 - m*x1
	
	TCS = y0 + m*E	

END SUBROUTINE universal_table_tcs

!#########################################################################
!#########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gives the total elastic cross section calculated from 
! the universal amplitude curve by numeric integration. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE universal_tcs_int( E, Mp, Mt, atom, TCS )
	USE physics_constants, ONLY : PI 

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E			! CM Energy [eV]
	REAL(KIND=8)		:: Mp			! Projectile mass [amu]
	REAL(KIND=8)		:: Mt			! Target mass [amu]
	INTEGER			:: atom			! 1 if atom-atom, 2 if atom-molecule

	!! Outputs
	REAL(KIND=8)		:: TCS			! Total Cross Section [a0^2]

	!! Internal
	REAL(KIND=8)		:: CC, A, B, C, lam, Mu, T_min, T_max, dT, T, ret, DCS
	INTEGER			:: N, i

	CC = PI/180.0D0

	!! calculate reduced mass
	Mu = Mp*Mt/(Mp+Mt)

	!! set minimum and maximum integration angles cm [deg]
	T_min = 0.01D0			
	T_max = 170.0D0	

	!! set number of integration points
	N  = 10000
	dT = (T_max-T_min)/REAL(N-1) 

	ret = 0.0D0

	DO i=1,N
		T   = T_min + (i-1)*dT
		CALL universal_dcs( E, Mu, T, atom, DCS )
		ret = ret + 2.0D0*PI*SIN(T*CC)*DCS*dT*CC
	END DO	
	
	TCS = ret	

END SUBROUTINE universal_tcs_int

!#############################################################
!#############################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Probability Density Routine for Universal Curve
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE universal_pd( E, Mp, Mt, atom )
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E		! [eV]
	REAL(KIND=8)	:: Mp		! [u]
	REAL(KIND=8)	:: Mt		! [u]
	INTEGER		:: atom		! 1 if atom-atom, 2 if atom-molecule

	!! Internal
	REAL(KIND=8)	:: PD, C, T, dT, tot, tcs, dcs, Mu
	INTEGER		:: i, N

	C = 180.0D0/PI

	Mu = Mp*Mt/(Mp+Mt)

	!! Number of integration points
	N  = 50000

	dT = PI/REAL(N-1)

	!! Get total cross section
	CALL universal_tcs_int( E, Mp, Mt, atom, tcs )

	tot = 0.0D0	

	DO i=1,N
		t   = REAL(i-1)*dT
		CALL universal_dcs_core( E, Mu, t*C, atom, dcs )
		PD  = 2.0D0*PI*SIN(T)*dcs*dT/tcs
		tot = tot + PD
	END DO	

END SUBROUTINE universal_pd

!#############################################################
!#############################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gives the total elastic cross section calculated from 
! the universal amplitude curve. 
!
! Analyitic solution using error function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE universal_tcs( E, Mp, Mt, TCS )
	USE physics_constants, ONLY : PI 

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E			! Energy [eV]
	REAL(KIND=8)		:: Mp			! Projectile mass [amu]
	REAL(KIND=8)		:: Mt			! Target mass [amu]

	!! Outputs
	REAL(KIND=8)		:: TCS			! Total Cross Section [a0^2]

	!! Internal
	REAL(KIND=8)		:: Mu, Tau_Min, Tau_Max, Theta_min, Theta_max, Z_min, Z_max
	REAL(KIND=8)		:: AA, A, B, C, gam, lam, CC, Err_min, Err_max

	!! Quadtratic fitting parameters
	A   = -0.13D0
	AA  = ABS(A)
	B   = 1.0D0
	C   = 2.7D0
	lam = 1.0D0		!! lam = 1.5 for atom-diatomic molecule collision, lam = 1.0 atom-atom collision
	gam = C + LOG(lam)
	CC  = PI*SQRT(PI/AA)*EXP(gam - (B*B)/(4.0D0*A))

	!! calculate reduced mass
	Mu = Mp*Mt/(Mp+Mt)

	!! calculate minimum and maximum scattering angles
  CALL uni_min_angle( Mp, Mt, 0.999D0, Theta_min )	
	Theta_min = Theta_min*180.0D0/PI
	IF (Mp .GE. Mt) Theta_max = ASIN(Mt/Mp)*180.0D0/PI
	IF (Mp .LT. Mt)	Theta_max = 180.0D0
	Tau_min   = E*Theta_min/Mu
	Tau_max   = E*Theta_max/Mu
	WRITE(*,'(A,4ES10.2)') "E: Theta_Min: Tau_Min: Tau_Max: ", E, Theta_min, Tau_min, Tau_max
!	WRITE(77,*) E*Theta_min/Mu, E

	!! calculate arguments for error functions
	Z_min = ( B + 2.0D0*A*LOG(Tau_min) ) / ( 2.0D0*SQRT(AA) )
	Z_max = ( B + 2.0D0*A*LOG(Tau_max) ) / ( 2.0D0*SQRT(AA) )

	!! calculate min and max error functions
	CALL errf(Z_min,Err_min)
	CALL errf(Z_max,Err_max)

!WRITE(*,*) Z_min, Z_max, Err_min, Err_max

	TCS = CC*(Err_min - Err_max)

END SUBROUTINE universal_tcs

!#############################################################
!#############################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gives the diffusion (transport) cross section calculated from 
! the universal amplitude curve by numeric integration. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE universal_dfcs_int( E, Mp, Mt, atom, TCS )
	USE physics_constants, ONLY : PI 

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E			! Energy [eV]
	REAL(KIND=8)		:: Mp			! Projectile mass [amu]
	REAL(KIND=8)		:: Mt			! Target mass [amu]
	INTEGER			:: atom			! 1 if atom-atom, 2 if atom-molecule

	!! Outputs
	REAL(KIND=8)		:: TCS			! Total Cross Section [a0^2]

	!! Internal
	REAL(KIND=8)		:: CC, A, B, C, lam, Mu, T_min, T_max, dT, T, ret, DCS
	INTEGER			:: N, i

	CC = PI/180.0D0

	!! calculate reduced mass
	Mu = Mp*Mt/(Mp+Mt)

	!! integration angles [deg] cm frame
	T_min = 0.01D0
	T_max = 170.0D0

	!! set number of integration points
	N  = 10000
	dT = (T_max-T_min)/REAL(N-1) 

	ret = 0.0D0

	DO i=1,N
		T   = T_min + (i-1)*dT
		CALL universal_dcs( E, Mu, T, atom, DCS )
		ret = ret + 2.0D0*PI*(1.0D0-COS(T*CC))*SIN(T*CC)*DCS*dT*CC
	END DO	
	
	TCS = ret	

END SUBROUTINE universal_dfcs_int

