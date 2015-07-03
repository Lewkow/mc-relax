!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate Chebyshev polynomials up to order 8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE chebyshev(n,X,T)

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: X
	INTEGER					:: n
	
	!! Outputs	
	REAL(KIND=8)		:: T

	!! Internal
	REAL(KIND=8)		:: X2, X3, X4, X5, X6, X7, X8

	X2 = X*X
	X3 = X2*X
	X4 = X2*X2
	X5 = X4*X
	X6 = X5*X
	X7 = X6*X
	X8 = X4*X4

	IF (n .EQ. 1) THEN
		T = X	
	ELSE IF (n .EQ. 2) THEN
		T = 2.0D0*X2 - 1.0D0
	ELSE IF (n .EQ. 3) THEN
		T = 4.0D0*X3 - 3.0D0*X
	ELSE IF (n .EQ. 4) THEN
		T = 8.0D0*X4 - 8.0D0*X2 + 1.0D0
	ELSE IF (n .EQ. 5) THEN
		T = 16.0D0*X5 - 20.0D0*X3 + 5.0D0*X
	ELSE IF (n .EQ. 6) THEN
		T = 32.0D0*X6 - 48.0D0*X4 + 18.0D0*X2 - 1.0D0
	ELSE IF (n .EQ. 7) THEN
		T = 64.0D0*X7 - 112.0D0*X5 + 56.0D0*X3 - 7.0D0*X
	ELSE IF (n .EQ. 8) THEN
		T = 128.0D0*X8 - 256.0D0*X6 + 160.0D0*X4 - 32.0D0*X2 + 1.0D0
	END IF

END SUBROUTINE chebyshev
