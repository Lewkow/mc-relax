
SUBROUTINE disp_relax_info(N, P, T)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Write information about particle
! relaxation to the screen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!-------------------
	! INPUTS	
	!-------------------

	INTEGER					:: N	! Number of collisions to thermalize
	INTEGER					:: P 	! MC Particle Number

	REAL(KIND=8)		:: T	! Time that P particle took to calculate relaxation

	IF (MOD(P,1000) == 0) THEN
		WRITE(*,*) "| Particle ", P, "| Num Collisions ", N, "| Calc Time ", T
	END IF

END SUBROUTINE

