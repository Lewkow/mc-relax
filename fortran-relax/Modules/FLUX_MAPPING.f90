
MODULE flux_mapping

	REAL(KIND=8),PARAMETER :: R1_map = 1.0D0 ! [AU]



CONTAINS

SUBROUTINE flux_map( x_now, x_nxt, E )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Determine if particle crosses mapping surface
! and has sunward velocity component
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 	IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)  :: x_now(3), x_nxt(3)
	REAL(KIND=8)	:: E

	!! Internal
	REAL(KIND=8)	:: R1, R2, R0, t0
	REAL(KIND=8)	:: Rm(3)


	R1 = SQRT(x_now(1)**2+x_now(2)**2+x_now(3)**2)	
	R2 = SQRT(x_nxt(1)**2+x_nxt(2)**2+x_nxt(3)**2)	
 
	IF (R1 .GT. R1_map .AND. R2 .LT. R1_map) THEN
		R0 		= SQRT( (x_nxt(1)-x_now(1))**2 + (x_nxt(2)-x_now(2))**2 + (x_nxt(3)-x_now(3))**2 )	
		t0 		= R1_map/R0
		Rm(1) = t0*(x_nxt(1)-x_now(1))
		Rm(2) = t0*(x_nxt(2)-x_now(2))
		Rm(3) = t0*(x_nxt(3)-x_now(3))
		WRITE(420,*) Rm(1), Rm(2), Rm(3), E
	END IF

END SUBROUTINE flux_map


END MODULE flux_mapping

