SUBROUTINE ion_transport_test
	
	IMPLICIT NONE

	REAL(KIND=8)			:: tt, amu, E, CX, R0(3), V0(3), rout(3), vout(3), path
	INTEGER						:: Z, i, N

	OPEN(UNIT=45,FILE='../Data/ion_trans_test_end_r.dat',ACCESS='APPEND')
	OPEN(UNIT=46,FILE='../Data/ion_trans_test_end_v.dat',ACCESS='APPEND')

	N   	= 10000
	Z   	= 2
	amu 	= 4
	E   	= 1.0D3
	CX  	= 1.0D-19
	R0(1) = 0.0D0
	R0(2) = 0.0D0
	R0(3) = 0.0D0
	V0(1) = 0.0D0
	V0(2) = 0.0D0
	V0(3) = 1.0D0

	DO i=1,N
		CALL ion_transport( amu, Z, E, CX, R0, V0, rout, vout, tt, path )
		write(45,*) rout(1), rout(2), rout(3)
		write(46,*) vout(1), vout(2), vout(3)
	END DO

	CLOSE(45)
	CLOSE(46)

END SUBROUTINE ion_transport_test

!#################################################################
!#################################################################

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! ion transport routine
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ion_transport( amu, Z, E, CX, R0, V0, rout, vout, tt, ion_path )
	USE ENA
	USE physics_constants, ONLY : PI
	USE mpi_info 

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Inputs
	REAL(KIND=8)			:: amu			! mass of ion [amu]
	INTEGER						:: Z				! positive charge of ion
	REAL(KIND=8)			:: E				! energy of ion [eV]
	REAL(KIND=8)			:: CX				! charge exchange cross section for ion with neutral H [m^2]
	REAL(KIND=8)			:: R0(3)		! initial position vector [x y z] [m]
	REAL(KIND=8)			:: V0(3)		! initial unit velocity vector [vx vy vz]

	!! Outputs
	REAL(KIND=8)			:: rout(3) 	! final position vector [x y z] [m]
	REAL(KIND=8)			:: vout(3) 	! final unit velocity vector [vx vy vz]
	REAL(KIND=8)			:: tt				! time during transport
	REAL(KIND=8)			:: ion_path	! path_length of ion transport [m] 

	!! Internal
	INTEGER						:: GOOD, NumL, i, j
	REAL(KIND=8)			:: M, M_p, q, JeV, Bt, Bp, B(3), eps, dt, L, lam, Mvout
	REAL(KIND=8)			:: v_perp, v_parr, VV0(3), v, pathL, d, RR, den, r, theta
	REAL(KIND=8)			:: AU, Rad, w0, v0_perp(3), v0_parr(3), C(3), n1(3), n2(3), MC
	REAL(KIND=8)			:: RC(3), rnow(3), rnxt(3), dL, t, P, nowL, rt(3), rp(3)
	REAL(KIND=8)			:: lfg
	REAL(KIND=8)			:: t_start, t_end, t_elap

	AU  = 1.5D11					! astronomical unit [m]
	M_p = 1.672D-27				! mass of proton [kg]
	q   = -1.602D-19			! charge of electron [C]
	JeV = 6.24150974D18		! Joules to eV  1 [J] = JeV [eV]

	!!! Get absolute distance from sun
	RR = SQRT(R0(1)**2 + R0(2)**2 + R0(3)**2)

	!!! Get magnetic field vector
	Bt   = 38.0D0*PI/180.0D0	
	Bp   = 23.0D0*PI/180.0D0
	B(1) = B0*COS(Bp)*COS(Bt)
	B(2) = B0*COS(Bp)*SIN(Bt)
	B(3) = B0*SIN(Bp)
!	Bt   = lfg()*PI	
!	Bp   = 2.0D0*PI*lfg()
!	B(1) = B0*SIN(Bt)*COS(Bp)
!	B(2) = B0*SIN(Bt)*SIN(Bp)
!	B(3) = B0*COS(Bt)

	!!! Get mass of ion
	M    = amu*M_p 
!WRITE(*,*) 'mass: ', M, amu, M_p

	!!! Percent of rad freq per time step
	eps  = 0.01D0

	!!! Time step
	dt   = eps*M/(Z*ABS(q)*B0)

	!!! Step length for transport [m]
	L    = AU

	!!! Mean free path [m]	
	CALL LISM_density(RR,den)
!	CALL LB_density(RR,den)

	lam  = 1.0D0/(den*CX) 
	
	!!! Get random path length
	P    = EXP(-L/lam)
	GOOD = 0
	NumL = 0
	
	DO WHILE (GOOD .EQ. 0)
		r = lfg()
		IF ( r .GE. P ) THEN 	!! CX collision back to neutral
			d    = -lam*LOG(r)
			GOOD = 1	
		ELSE 									!! No CX collision, keep transporting
			NumL = NumL + 1
		END IF
	END DO 	

	pathL = L*NumL + d	
	ion_path = pathL
! WRITE(*,*) "MFP: Rand FP: [AU] ", lam/AU, pathL/AU 

	!!! get initial velocity
	v   = SQRT(2.0D0*E/(M*JeV))		! magnitude of velocity [m/s]
	VV0 = V0*v										! initial velocity vector [m/s]
!WRITE(*,'(A,7ES10.2)') 'v: V0: VV0: ', v, V0(1), V0(2), V0(3), VV0(1), VV0(2), VV0(3)

	!!! get angle between B and VV0
	theta = ACOS( (B(1)*VV0(1)+B(2)*VV0(2)+B(3)*VV0(3))/(v*B0) )
!WRITE(*,*) 'B: ', B(1), B(2), B(3)
!WRITE(*,*) 'V: ', VV0(1), VV0(2), VV0(3)
!WRITE(*,*) 'theta params: ', B(1)*VV0(1), B(2)*VV0(2), B(3)*VV0(3)
!WRITE(*,*) 'theta: ', (B(1)*VV0(1)+B(2)*VV0(2)+B(3)*VV0(3))/(V*B0), theta
	
	!!! get perp and parr velocity vectors
	v_perp  = v*SIN(theta)
!WRITE(*,*) 'v: t: S(t): vp: ', v, theta*180.0D0/PI, SIN(theta), v_perp
	v_parr  = v*COS(theta)	
	v0_perp = VV0 - ( (VV0(1)*B(1)+VV0(2)*B(2)+VV0(3)*B(3))/(B0*B0) )*B
	v0_parr =  ( (VV0(1)*B(1)+VV0(2)*B(2)+VV0(3)*B(3))/(B0*B0) )*B

	!!! get natural frequency and radius of oscillation
	w0      = Z*q*B0/M
	Rad     = M*v_perp/ABS(Z*q*B0)	
	
	!!! get unit vector for perp plane
	C(1)    = v0_perp(2)*B(3) - v0_perp(3)*B(2)
	C(2)    = v0_perp(3)*B(1) - v0_perp(1)*B(3)
	C(3)    = v0_perp(1)*B(2) - v0_perp(2)*B(1)
	MC      = SQRT( C(1)**2 + C(2)**2 + C(3)**2 )
	n1      = C/MC
	C(1)    = n1(2)*B(3) - n1(3)*B(2)
	C(2)    = n1(3)*B(1) - n1(1)*B(3)
	C(3)    = n1(1)*B(2) - n1(2)*B(1)
	MC      = SQRT( C(1)**2 + C(2)**2 + C(3)**2 )
	n2      = C/MC

	!!! get constant for position vector
	RC      = R0 - Rad*n1
	
	!!! get position and velocity at time t	
	rp   = v0_parr*t
	rt   = Rad*COS(w0*t)*n1 + Rad*SIN(w0*t)*n2
	rout = RC + rp + rt
	vout = v0_parr + Rad*w0*( -SIN(w0*t)*n1 + COS(w0*t)*n2 )

	Mvout = SQRT( vout(1)**2 + vout(2)**2 + vout(3)**2 )	

	!!! get time at which particle has traveled pathL
	t       = pathL/Mvout

!	write(*,'(A,3ES12.2)') 'mfp: v: t: ', pathL, Mvout, t

	vout(1) = vout(1)/Mvout	
	vout(2) = vout(2)/Mvout	
	vout(3) = vout(3)/Mvout	

	tt      = t

!WRITE(*,*)
!WRITE(*,*) '+++++++++++++++++++++++'
!WRITE(*,*) 'V, vxyz: ', Mvout, vout(1), vout(2), vout(3)
!WRITE(*,*) M, v_perp, ABS(Z*q*B0) 
!WRITE(*,*) Rad
!WRITE(*,*) R0(1), R0(2), R0(3)
!WRITE(*,*) n1(1), n1(2), n1(3)
!WRITE(*,*) RC(1), RC(2), RC(3)
!WRITE(*,*) rp(1), rp(2), rp(3)
!WRITE(*,*) rt(1), rt(2), rt(3)
!WRITE(*,*) '+++++++++++++++++++++++'

END SUBROUTINE ion_transport

!#################################################################
!#################################################################
SUBROUTINE butler_drag( z, d, n, T, E, E_out )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates energy loss due to ion-electron interaction
! within the ISM plasma using Butler & Buckingham 1962
!
! E(t) = E0*exp(-2*t/Tau0)
! E0   = initial energy 
! t    = time exposed to plasma [sec]
! T0   = 3*mi*me*ve^3/(16*sqrt(pi)*k^2*z^2*e^4*ne*log(lam)
! lam  = 2/t_m
! t_m  = pi-2*acot(k*q1*q2/(2*E*DL))
! DL   = sqrt(eps*kb*T/(n*e^2))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI, EVTOJ, eps0, mass_e, mass_p, mass_n, qe, UTOKG, KB

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E		! initial energy of ion [eV]
	REAL(KIND=8)	:: d		! distance traveled as ion [m]
	REAL(KIND=8)	:: n		! plasma density [m^-3]
	REAL(KIND=8)	:: T		! plasma temperature [K]
	INTEGER				:: z		! unit of charge for ion

	!! Outputs
	REAL(KIND=8)	:: E_out	! energy after transport [eV]

	!! Internal
	REAL(KIND=8)	:: vi, ve, DL, tt, tm, lam, T0

	DL = SQRT(eps0*kb*T/(n*qe*qe))		! Debye length [m]

	IF (z .EQ. 1) vi = SQRT(2.0D0*E*EVTOJ/(mass_p*UTOKG))
	IF (z .EQ. 2) vi = SQRT(2.0D0*E*EVTOJ/(4.0D0*mass_p*UTOKG))
	
	ve  = SQRT(KB*T/(mass_e*UTOKG))
	tt  = d/vi ! [s]
	tm  = PI - 2.0D0*ATAN(2.0D0*E*EVTOJ*DL/((1.0D0/(4.0E0*PI*eps0))*qe*qe))
	lam = 2.0D0/tm
	T0  = 3.0D0*mass_p*UTOKG*mass_e*UTOKG*ve*ve*ve/ &
	&     (16.0D0*SQRT(PI)*(1.0D0/(4.0D0*PI*eps0))**2*REAL(z)**2*qe**4*n*LOG(lam))

	E_out = E*EXP(-2.0D0*tt/T0)

END SUBROUTINE butler_drag

!#################################################################
!#################################################################

SUBROUTINE bethe_drag( z, d, n, E, E_out )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates energy loss due to ion-electron interaction
! within the ISM plasma
!
! dE/dx = - (4*pi*n*z^2/(me*v^2))*(e^2/(4*pi*eps0))^2*
!           log( 2*me*v^2/I )
!
! I = mean excitation potential
! I(z=1) = 19 eV
! I(z=2) = 42 eV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI, EVTOJ, eps0, mass_e, mass_p, mass_n, qe, UTOKG

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E		! initial energy of ion [eV]
	REAL(KIND=8)	:: d		! distance traveled as ion [m]
	REAL(KIND=8)	:: n		! ion density [m^-3]
	INTEGER				:: z		! unit of charge for ion

	!! Outputs
	REAL(KIND=8)	:: E_out	! energy after transport [eV]

	!! Internal
	REAL(KIND=8)	:: C1, C2, m, me, v, k, EJ, I, drag

	me = mass_e*UTOKG

	IF (z .EQ. 1) THEN
		I = 19.0D0*EVTOJ
		m = mass_p*UTOKG
	ELSE IF (z .EQ. 2) THEN
		I = 42.0D0*EVTOJ
		m = 2.0D0*(mass_p+mass_n)*UTOKG
	END IF

	k  = 1.0D0/(4.0D0*PI*eps0)
	EJ = E*EVTOJ
	v  = SQRT( 2.0D0*EJ/m )
!	WRITE(*,*) E, EJ, m, v/1000.0

	C1   = -((4.0D0*PI*n*REAL(z*z))/(me*v*v))*((qe*qe*k)**2)/EVTOJ
	C2   = LOG(2.0D0*me*v*v/I)
	drag = C1*C2

!	WRITE(520,*) v/1000.0D0, C1, C2, drag

	IF (drag .GT. E) THEN
		E_out = 0.0D0
	ELSE
		E_out = E - d*drag
	END IF	

END SUBROUTINE bethe_drag

!#################################################################
!#################################################################

SUBROUTINE ion_drag( r, E, Proj, drag )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Ion Drag coefficient such that
! E = E0*exp(-t/drag) after t seconds of ion transport
!
! drag = 1/3 sqrt(2*pi) n * vt * rho^2 * CL
! vt	 = sqrt( KT/m )
! rho  = (1/4/pi/eps) * q^2/K/T
! CL   = int_{0}^{inf} dx * 2 * exp( -x ) * 
!					log( (2*DB*x + rho)/(2*a*x+rho) ) 
! DB   = sqrt( eps*K*T/n/q^2 )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI, KB, EVTOK, EVTOJ

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: r, E
	CHARACTER(len=2)	:: Proj

	!! Outputs
	REAL(KIND=8)	:: drag

	!! Internal	
	REAL(KIND=8)	:: eps, qe, DB, CL, rho, vt, dx, tot, sml
	REAL(KIND=8)	:: now, nxt, dff, T, x, a, den, m, mp, vel, FC
	INTEGER				:: GO, i

	FC  = 1.602D-19
	eps = 8.8541878D-12								! [C^2/N/m^2]
	qe  = 1.6021765D-19								! [C]
	mp  = 1.67262178D-27							! [kg]
	!! get density and temperture
	CALL LISM_ion_density( r, den )
	CALL LISM_temp( r, T )
!	CALL LB_density( r, den )
!	CALL LB_temp( r, T )

	!! get ion mass
	IF (PROJ .EQ. 'H ') THEN
		m = mp
	ELSE
		m = 4.0D0*mp
	END IF

	T   = T*EVTOK
	DB  = SQRT(eps*kb*T/(den*qe**2)) 		! [m]
	rho = qe**2/(kb*T*4.0D0*PI*eps)		! [m]
	dx  = 1.0D-2	 										! [m]
	a   = 1.0D-10											! [m] "dust" size
	vt  = SQRT(kb*T/m)								
	tot = 0.0D0 
	nxt = 0.0D0
	sml = 1.0D-6 
	GO  = 1 
	i   = 0
	vel = SQRT(2.0D0*E*EVTOJ/m)
 
	DO WHILE (GO .EQ. 1)
		x   = REAL(i*dx)
		now = 2.0D0*dx*EXP(-x)*LOG( (2.0D0*DB*x+rho)/(2.0D0*a*x+rho) )
		tot = tot + now
		IF ( (dff .LE. sml) .AND. (i .GT. 0) ) THEN
			GO = 0
		ELSE	
			i = i+1
		END IF
	END DO

	CL   = tot

!	drag = (2.0D0/3.0D0)*SQRT(2.0D0*PI)*den*m*vt*rho*rho*CL*vel
	drag = (2.0D0/3.0D0)*SQRT(2.0D0*PI)*den*m*vt*rho*rho*CL*vel/FC

!	drag = (1.0D0/3.0D0)*SQRT(2.0D0*PI)*den*vt*rho*rho*CL	

!	WRITE(*,'(7ES12.2)') T, DB, rho, vt, CL, vel, drag

END SUBROUTINE ion_drag

