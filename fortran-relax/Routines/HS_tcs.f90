!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! HS total cross sections for all 
!! collision types
!!
!! TCS = pi*R^2
!! R   = rp + rt
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HS_tcs( RP, RT, tcs )
	USE physics_constants, ONLY : PI	

	IMPLICIT NONE 

	!! Outputs
	REAL(KIND=8)		:: tcs


	!! Inputs
	REAL(KIND=8)		:: RP, RT

	!! Internal
	REAL(KIND=8)		:: r	

	r   = RP + RT 
	tcs = PI*r*r

END SUBROUTINE HS_tcs

