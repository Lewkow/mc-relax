SUBROUTINE uni_ang_test

	IMPLICIT NONE

	INTEGER				:: N, NE, i, j, atom
	REAL(KIND=8)	:: C, tot1, tot2, E_min, E_max, dE, E, Mp, Mt, T1, T2, LT2
	REAL(KIND=8)	:: TCS1, TCS2, p, dE1, dE2, rat, per_rat

	C = 180.0D0/3.14D0

	N      = 5000
	NE 		 = 50
	E_min  = 10.0D0
	E_max  = 1000.0D0
	dE     = (E_max-E_min)/REAL(NE-1)	
	Mp 		 = 4.0D0
	Mt     = 1.0D0
	p      = Mp/Mt
	atom   = 1 

	OPEN(UNIT=55,FILE='../Data/Uni_test.dat',ACCESS='APPEND')

	DO j=1,NE
		E 	 = E_min + (j-1)*dE
		tot1 = 0.0D0
		tot2 = 0.0D0
		DO i=1,N
			CALL universal_rand_angle_int( E, Mp, Mt, atom, T1 )
			tot1 = tot1+ T1
			CALL lin_rand_angle( E, 'HeH ', T2 )
			tot2 = tot2+ T2
		END DO
		T1  = tot1/REAL(N)
		T2  = tot2/REAL(N)
		dE1 = E*((2.0D0*p)/(1.0D0+p)**2)*(1.0D0-COS(T1))
		dE2 = E*((2.0D0*p)/(1.0D0+p)**2)*(1.0D0-COS(T2))
		CALL universal_tcs_int(E,Mp,Mt,atom,TCS1)
		CALL TCS_HeO(E,TCS2)
		rat = (TCS2*dE2)/(TCS1*dE1)
		per_rat = ABS(1.0D0-rat)*100.0D0
!		WRITE(55,*) E, T1*C
!		WRITE(*,*) E, T1*C
		WRITE(55,'(9ES10.2)') E, T1*C, T2*C, dE1, dE2, TCS1, TCS2, rat, per_rat
		WRITE(*,'(9ES10.2)') E, T1*C, T2*C, dE1, dE2, TCS1, TCS2, rat, per_rat
	END DO

	CLOSE(55)

END SUBROUTINE uni_ang_test

!##########################################################################
!##########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gives a random scattering angle using the universal
! amplitude curve fit with a quadradic. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE universal_rand_angle( E, Mp, Mt, Theta )
	USE physics_constants, ONLY : PI 

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E			! Energy [eV] 
	REAL(KIND=8)		:: Mp			! Projectile mass [amu]
	REAL(KIND=8)		:: Mt			! Target mass [amu]

	!! Outputs
	REAL(KIND=8)		:: Theta	! Lab Scattering Angle [deg]

	!! Internal
	REAL(KIND=8)		:: Z_min, Z_max, Mu, Tau_min, Theta_min, r, lfg, A, B
	REAL(KIND=8)		:: Z, T, Err, Err_min, Err_max, Theta_max, dT, tot
	INTEGER					:: i

	!! Quadtratic fitting parameters
	A = -0.13D0
	B = 1.0D0

	!! calculate reduced mass
	Mu = Mp*Mt/(Mp+Mt)

	!! calculate minimum and maximum scattering angles
	Tau_min   = 1.0D-2					! Minimum tau from fitting function
	Theta_min = Tau_min*Mu/E

	!! set angular increments for integration
	dT = 0.001									! [deg]

	IF (Mp .GE. Mt) THEN
		Theta_max = ASIN(Mt/Mp)*180.0D0/PI
	ELSE
		Theta_max = 180.0D0
	END IF

	!! calculate arguments for error functions
	Z_min = (B + 2.0D0*A*LOG(E*Theta_min/Mu))/(2.0D0*SQRT(ABS(A)))
	Z_max = (B + 2.0D0*A*LOG(E*Theta_max/Mu))/(2.0D0*SQRT(ABS(A)))

	!! calculate min and max error functions
	CALL errf(Z_min,Err_min)
	CALL errf(Z_max,Err_max)

	!! get random number for scattering angle
	r = lfg()
	
	!! set initial tot
	tot = 0.0D0
	i   = 0

	DO WHILE (tot .LT. r)
		T = Theta_min + i*dT
		Z = (B + 2.0D0*A*LOG(E*T/Mu))/(2.0D0*SQRT(ABS(A)))
		CALL errf(Z,Err)	
		tot = (Err - Err_min)/(Err_max - Err_min)
		i   = i+1
	END DO

	!! Convert angle from deg to rad
	Theta = T*PI/180.0D0

END SUBROUTINE universal_rand_angle

!##########################################################################
!##########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Universal random angle from integration
!
! This routine finds a random angle from the 
! universal curve using numeric integration
! for both the total cross section and 
! for the probability density. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE universal_rand_angle_int( E, Mp, Mt, atom, SAng )
	USE physics_constants, ONLY : PI
	
	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E		! Collision Energy [eV] Lab Frame
	REAL(KIND=8)		:: Mp		! Projectile Mass
	REAL(KIND=8)		:: Mt 	! Target Mass
	INTEGER					:: atom	! 1 if atom-atom, 2 if atom-molecule

	!! Outputs
	REAL(KIND=8)		:: SAng	! Random Scattering Angle [rad] CM Frame

	!! Internal
	REAL(KIND=8)		:: lfg	! Random number generator function
	REAL(KIND=8)		:: Cang, tcs, dcs, tot, t, dt, t_min, t_max, r, Mu, cc  
	REAL(KIND=8)		:: E_CM, Y1, Y2, X1, X2, Y0, m, pd
	INTEGER					:: N, i, C

	CC = 180.0D0/PI

	!! Convert collision energy to CM
	E_CM = mt*E/(mp+mt)

	!! Set reduced mass
	Mu = mt*mp/(mt+mp)

	!! Get tcs for current energy
	CALL universal_tcs_int( E_CM, Mp, Mt, atom, tcs )

	!! Set number of integration points for integration over prob density 
	N = 1000

	!! Set minimum and maximum angles used for integration
	!! [rad] in CM Frame
	t_min = 0.02D0/CC
	t_max = 170.0D0/CC 

	!! Set differential for angular integration
	!! [rad] in CM Frame
	dt = (t_max-t_min)/REAL(N-1)

	!! Get random number associated with scattering angle
	r = lfg()

	!! Initialize and zero integration sum 
	tot = 0.0D0

	C = 1
	i = 1

	DO WHILE (C .EQ. 1)
		t   = t_min + (i-1)*dt
		CALL universal_dcs( E_CM, Mu, t*CC, atom, dcs )		
		pd  = 2.0D0*PI*SIN(t)*dcs*dt/tcs
		tot = tot + pd
		IF (tot .GT. r) THEN
			X1 		= tot - pd
			X2 		= tot
			Y1 		= REAL(i)*dt
			Y2 		= REAL(i+1)*dt
			m  		= (Y2-Y1)/(X2-X1)
			Y0 		= Y2 - m*X2
			Cang 	= Y0 + m*r
			C    	= 0
		ELSE
			i = i+1
		END IF
	END DO	

	!! Convert scattering angle from CM to Lab Frame
!	CALL angle_to_lab(Cang, mp/mt, Sang)	

	IF ( ( Cang .LT. 0.0D0 ) .OR. ( Cang .GT. PI ) ) THEN
		SAng = PI/10.0D0
	ELSE
		SAng = Cang
	END IF

END SUBROUTINE universal_rand_angle_int

