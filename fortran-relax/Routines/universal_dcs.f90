SUBROUTINE test_universal_dcs

	IMPLICIT NONE

	REAL(KIND=8)	:: E, T, dT, DCS, Mp, Mt, Mu
	INTEGER				:: N, i, atom

	N  = 50000
	Mp = 4.0D0
	Mt = 16.0D0	
	Mu = Mp*Mt/(Mp+Mt)
	dT = 180.0D0/REAL(N)

	atom = 1 

	E = 500.0D0

	DO i=1,N
		T = i*dT
		CALL universal_dcs_core(E,Mu,T,atom,DCS)
		WRITE(99,*) T, DCS
	END DO	

END SUBROUTINE test_universal_dcs 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate elastic differential cross section
! using universal curve fit with quadratic eqn
! core fit with analyitic form
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE universal_dcs( E, Mu, theta, atom, DCS )
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)			:: E		!! Energy in [eV]
	REAL(KIND=8)			:: Mu 	!! Reduced mass in [amu]
	REAL(KIND=8)			:: theta!! Angle [deg]
	INTEGER						:: atom !! 1 if atom-atom, 2 if atom-molecule

	!! Outputs
	REAL(KIND=8)			:: DCS	!! Differential Cross Section [a0^2]

	!! Internal
	REAL(KIND=8)			:: CC, x0, y0, m, A, B, C, Z, lam, G, T, ZZ, C0

	CC = PI/180.0D0

	IF (atom .EQ. 1) THEN
		lam = 1.0D0
	ELSE
		lam = 1.4D0
	END IF

	A = -0.13
	B = 1.0
	C = 2.7
	G = C + log(lam)
	T = E*theta/mu
	Z = 1/(theta*sin(theta*pi/180))

	x0 = 50.12
	m  = -0.02574
	y0 = 2.03811
	C0 = 32.3416
	ZZ = 10

	IF (T > x0) THEN
  	DCS = Z*exp(A*(log(T))**2 + B*log(T) + G)
	ELSE
  	CC  = ZZ*exp(y0 + m*log(T)) + C0
  	DCS = lam*Z*CC
	END IF

END SUBROUTINE universal_dcs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate elastic differential cross section
! using universal curve fit with quadratic eqn
! with core for small angles using guassian
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE universal_dcs_core( E, Mu, T, atom, DCS )
	USE universal
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)			:: E		!! Energy in [eV]
	REAL(KIND=8)			:: Mu 	!! Reduced mass in [amu]
	REAL(KIND=8)			:: T 		!! Angle [deg]
	INTEGER						:: atom !! 1 if atom-atom, 2 if atom-molecule

	!! Outputs
	REAL(KIND=8)			:: DCS	!! Differential Cross Section [a0^2]

	!! Internal
	REAL(KIND=8)			:: A, B, C, lam, gam, TT, Tau, t0, tb

	IF (atom .EQ. 1) 									 lam = 1.0D0
	IF (atom .EQ. 2) 									 lam = 1.4D0
	IF (atom .NE. 1 .AND. atom .NE. 2) lam = 1.0D0

	A   = -0.13D0
	B   = 1.0D0
	C   = 2.7D0
	gam = C + LOG(lam)

	!! Get core parameters
	CALL uc_tcs( E )	
!	CALL uc_ael( E )	
!	CALL uc_dEdt( E )	

	TT  = T*E/Mu

	t0  = tau_0*Mu/E	!! [deg]
	tb  = tau_b*Mu/E	!! [deg]

	IF ( TT .LE. tau_0 ) THEN
		Tau = tau_0
		DCS = (EXP(A*(LOG(Tau)*LOG(Tau))+B*LOG(Tau)+gam)/((Mu*Tau/E)*SIN(Tau*PI*Mu/(180.0D0*E)))) &
		&     *EXP( (Tau**2-TT**2)/(2.0D0*tau_b**2) )
	ELSE
		Tau = TT 
		DCS = EXP(A*(LOG(Tau)*LOG(Tau))+B*LOG(Tau)+gam)/((Mu*TT/E)*SIN(Mu*TT*PI/(180.0D0*E)))
	END IF

END SUBROUTINE universal_dcs_core

