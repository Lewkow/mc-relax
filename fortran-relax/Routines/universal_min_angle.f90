
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates the minimal angle to use for Universal
! curve as determined by the minimum percent energy
! loss input into function and masses of projectile
! and target particles. 
!
! 
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE uni_min_angle( Mp, Mt, X, T )
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)		:: Mp			! Mass of projectile
	REAL(KIND=8)		:: Mt			! Mass of target
	REAL(KIND=8)		:: X			! Minimum percentage energy loss 
														! such that E = E0*X
														! E  -> Final Energy  
														! E0 -> Initial Energy  

	!! Outputs
	REAL(KIND=8)		:: T			! Minimum lab frame scattering
														! angle for universal curve

	!! Internal
	REAL(KIND=8)		:: C1, C2, C, Mu, T_CM

	Mu = Mp*Mt/(Mp+Mt)
	C  = Mu/Mp

	C1 = (0.5D0/C)*(X-1.0D0) - C + 1.0D0
	C2 = 1.0D0 - C

	T_CM = ACOS(C1/C2)

	!! convert scattering angle from cm to lab frame
	CALL angle_to_lab(T_CM, Mp/Mt, T)	
!WRITE(*,*) "C1: C2: T_CM: T_LB:", C1, C2, T_CM*180.0D0/PI, T*180.0D0/PI

END SUBROUTINE uni_min_angle

