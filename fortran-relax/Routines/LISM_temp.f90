
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Tempeture fucntion to be called for LISM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LISM_temp( r, temp )
	USE physics_constants, ONLY : PCTOAU

  IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)      :: r    !! radial distance from sol [m]

  !! Outputs
  REAL(KIND=8)      :: temp  !! density of LISM neutral H [1/m^3]

  !! if in local cloud
  IF (r .LE. 4.5D0*PCTOAU) THEN
    temp = 0.6D0   	!! [eV] = 7000 K 
  !! else in local bubble
  ELSE
    temp = 86.0D0 	!! [eV] = 10^6 K
  END IF

END SUBROUTINE LISM_temp

!############################################
!############################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Tempeture fucntion to be called for Local Bubble
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LB_temp( r, temp )
	USE physics_constants, ONLY : PCTOAU

  IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)      :: r    !! radial distance from sol [m]

  !! Outputs
  REAL(KIND=8)      :: temp  !! density of LISM neutral H [1/m^3]

  temp = 86.0D0 	!! [eV] 

END SUBROUTINE LB_temp


