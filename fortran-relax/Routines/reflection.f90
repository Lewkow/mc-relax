
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculate reflection coefficient for current
!! ensemble and write to file
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE reflection
	USE click_data, ONLY : ROOT_ESC_ENGY_DIST
	USE planet, ONLY : N_Part, high, INIT_ENGY, h_0
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER				:: i, St, Fn
	REAL(KIND=8)	:: N_Esc, tot, ref, root_tot

	tot	= 0.0D0

	IF (myid .EQ. 0) THEN
	
		DO i=1,N_Part
			IF (ROOT_ESC_ENGY_DIST(i) .NE. 0.0D0) tot = tot + 1.0D0
		END DO

		Ref = tot/REAL(N_Part)

		39 FORMAT (A,ES8.2)
		40 FORMAT (A,ES8.2,A)
		41 FORMAT (A)

		WRITE(*,*)
		WRITE(*,40) 'Initial height: ', h_0/1.0D3, ' [km]'
		WRITE(*,40) 'Initial energy: ', INIT_ENGY, ' [eV]'
		WRITE(*,41) '----------------'
		WRITE(*,39) 'Reflection Coeff: ', Ref
		WRITE(*,*)

	END IF



END SUBROUTINE reflection

