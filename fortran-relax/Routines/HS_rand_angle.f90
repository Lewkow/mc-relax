!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! Generate random scattering angles using
!! hard sphere cross sections
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HS_rand_angle( ScattAng )
	USE rand_seed

	IMPLICIT NONE

	!! Outputs
	REAL(KIND=8)			:: ScattAng

	!! Internal
	REAL(KIND=8)			:: r
	REAL(KIND=8)			:: lfg

	r        = lfg()
	ScattAng = ACOS(1.0D0-2.0D0*r)

END SUBROUTINE HS_rand_angle


