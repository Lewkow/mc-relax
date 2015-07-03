
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Parameter routines for core fitting of the
! universal curve. 
!
! for t<=t0
!
! |f(E,t)|^2 = |f(E,t0)|^2 * 
!							 exp( (t0^2 - t^2)/(2*tb^2) )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!########################################################
!########################################################

SUBROUTINE uc_tcs( E )
	USE universal

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E			! [eV]

	!! Internal
	REAL(KIND=8)	:: p1_0, p2_0, p1_b, p2_b

	p1_0 =  0.126d0
	p2_0 = -0.288d0

	p1_b =  0.0525d0
	p2_b = -2.3d0

	tau_0 = p1_0*E + p2_0
	tau_b = p1_b*E + p2_b

END SUBROUTINE uc_tcs

!########################################################
!########################################################

SUBROUTINE uc_ael( E )
	USE universal	

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E			! [eV]

	!! Internal
	REAL(KIND=8)	:: p1_0, p2_0, p1_b, p2_b

	p1_0 =  0.138d0
	p2_0 =  2.27d0

	p1_b =  0.0644d0
	p2_b = -3.76d0

	tau_0 = p1_0*E + p2_0
	tau_b = p1_b*E + p2_b

END SUBROUTINE uc_ael

!########################################################
!########################################################

SUBROUTINE uc_dEdt( E )
	USE universal	

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E			! [eV]

	!! Internal
	REAL(KIND=8)	:: p1_0, p2_0, p1_b, p2_b

	p1_0 =  0.0119d0
	p2_0 =  12.0d0

	p1_b =  0.0027d0
	p2_b =  2.5d0

	tau_0 = p1_0*E + p2_0
	tau_b = p1_b*E + p2_b

END SUBROUTINE uc_dEdt

