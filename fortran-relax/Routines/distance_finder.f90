
SUBROUTINE distance_finder( x, y, z, L, d )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Distance finder takes 3 arrays, x, y, z, of
! length L and determines the dissplacment
! array d such that
! d(i) = SQRT( x(i)^2 + y(i)^2 + z(i)^2 )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!! Inputs
	INTEGER										:: L
	REAL(KIND=8),DIMENSION(L)	:: x
	REAL(KIND=8),DIMENSION(L)	:: y
	REAL(KIND=8),DIMENSION(L)	:: z

	!! Outputs 
	REAL(KIND=8),DIMENSION(L)	:: d
	INTEGER										:: i

	DO i=1,L
		d(i) = SQRT( x(i)**2.0D0 + y(i)**2.0D0 + z(i)**2.0D0 )
	END DO

END SUBROUTINE distance_finder
