
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate TCS for elastic H+CO2 collisions
! using same methods as Kallos et al. 2000. 
!
! TCS(E) = a0*E^a
!	E: Collision energy in keV
! a0: 2.1D-16	[cm^2]
! a: -.666 
!
! Input E in [eV] and output TCS in [m^2]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TCS_HCO2( E, TCS )

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E

	!! Outputs
	REAL(KIND=8)		:: TCS

	!! Internal
	REAL(KIND=8)		:: a0, a, C

	a0 = 2.1D-16
	a  = -0.666D0
	C  = 1.0D0/100.0D0
	C  = C*C
	a0 = a0*C

	TCS = a0*(E/1000)**a

END SUBROUTINE
