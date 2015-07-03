
SUBROUTINE angle_to_lab(theta_cm, rat, theta_lab)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Converts scattering angle from center of mass (CM)
! to lab frame. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!------------------------
	! Inputs
	!------------------------

	REAL(KIND=8)				:: theta_cm	!! [rad]
	REAL(KIND=8)				:: rat			!! Ratio of masses rat = MP/MT

	!------------------------
	! Outputs
	!------------------------

	REAL(KIND=8)				:: theta_lab

	!------------------------
	! Internal 
	!------------------------

	REAL(KIND=8)				:: ret

	ret       = ( COS(theta_cm) + rat )/( SQRT(1.0D0 + 2.0D0*rat*COS(theta_cm) + rat*rat) )
	theta_lab = ACOS(ret)
	
!	ret 			= SIN(theta_cm)/(rat + COS(theta_cm))
!	theta_lab = DATAN(ret) 

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

