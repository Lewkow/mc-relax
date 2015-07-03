!############################################################
!
! Calculate the total cross sections for H+H collisions
! non-symmetric collisions
!
! TCS for energies between 0.1 eV - 100 eV were confirmed
!	when compared to
!
! Consistent definitions for, and relationships among, 
! cross sections for elastic scattering of hydrogen ion, 
! atoms, and molecules, Krstic and Schultz, 1999, PRA, 60, 3. 
!
!############################################################
SUBROUTINE TCS_HH( E, TCS )
	
	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: E 	! eV
	
	!! Outputs
	REAL(KIND=8)		:: TCS ! a0^2

	!! Internal	
	REAL(KIND=8)		:: L(12), R(12), alp(12), sig(12)
	INTEGER					:: i, j

	R(1)  = 0.07D0
	R(2)  = 0.25D0
	R(3)  = 5.0D0
	R(4)  = 10.0D0
	R(5)  = 50.0D0
	R(6)  = 300.0D0
	R(7)  = 420.0D0
	R(8)  = 870.0D0
	R(9)  = 1000.0D0
	R(10) = 3000.0D0
	R(11) = 5680.0D0
	R(12) = 11000.0D0

	L(1)  = 0.0D0
	L(2)  = R(1)
	L(3)  = R(2)
	L(4)  = R(3)
	L(5)  = R(4)
	L(6)  = R(5)
	L(7)  = R(6)
	L(8)  = R(7)
	L(9)  = R(8)
	L(10) = R(9)
	L(11) = R(10)
	L(12) = R(11)

	alp(1)  = 0.13445D0
	alp(2)  = 0.10158D0
	alp(3)  = 0.11143D0
	alp(4)  = 0.13537D0
  alp(5)  = 0.16816D0
  alp(6)  = 0.19399D0
	alp(7)  = 0.14963D0
	alp(8)  = 0.28865D0
	alp(9)  = 0.40730D0
	alp(10) = 0.51443D0
	alp(11) = 0.75148D0
	alp(12) = 0.73176D0

	sig(1)  = 61.275D0
	sig(2)  = 85.378D0
	sig(3)  = 78.972D0
	sig(4)  = 68.640D0
  sig(5)  = 59.743D0
	sig(6)  = 55.041D0
	sig(7)  = 58.450D0
	sig(8)  = 51.792D0
  sig(9)  = 50.963D0
	sig(10) = 51.425D0
	sig(11) = 67.559D0
	sig(12) = 65.714D0

	DO i=1,12
		IF ( ( E .GE. L(i) ) .AND. ( E .LE. R(i) ) ) THEN
			TCS = sig(i)*(1000.0D0/E)**alp(i)
		END IF
	END DO

END SUBROUTINE TCS_HH



