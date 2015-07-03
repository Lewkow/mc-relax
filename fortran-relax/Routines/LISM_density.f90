
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Density fucntion to be called for LISM
!
! Neutral density in LIC taken from Table 1 Frisch et al. 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LISM_density( r, den )
	USE physics_constants, ONLY : PCTOAU
	
	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)			:: r		!! radial distance from sol [m]
	
	!! Outputs
	REAL(KIND=8)			:: den	!! density of LISM neutral H [1/m^3]

	!! if in local cloud
	IF (r .LE. 4.5D0*PCTOAU) THEN
		den = 0.2D6		!! [1/m^3]
	!! else in local bubble
	ELSE 
		den = 0.0D0	!! [1/m^3] 
	END IF

END SUBROUTINE LISM_density

!############################################
!############################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Proton Density fucntion to be called for LISM
!
! Ion density in LIC taken from Table 1 Frisch et al. 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LISM_ion_density( r, den )
	USE physics_constants, ONLY : PCTOAU
	
	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)			:: r		!! radial distance from sol [m]
	
	!! Outputs
	REAL(KIND=8)			:: den	!! density of LISM neutral H [1/m^3]

	!! if in local cloud
	IF (r .LE. 4.5D0*PCTOAU) THEN
		den = 0.018D6		!! [1/m^3]
	!! else in local bubble
	ELSE 
		den = 0.005D6	!! [1/m^3] 
	END IF

END SUBROUTINE LISM_ion_density

!############################################
!############################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hydrogen density in local bubble, constant and taken as 
! 5% of total density as designated in Welsh et al. 2009. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LB_density( r, den )

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)			:: r		!! radial distance from star

	!! Outputs
	REAL(KIND=8)			:: den	!! density of LB hydrogen

	den = 0.005D6 ! [1/m^3]

END SUBROUTINE LB_density

!############################################
!############################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Proton density in local bubble, constant and taken as 
! 95% of total density as designated in Welsh et al. 2009. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LB_ion_density( r, den )

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)			:: r		!! radial distance from star

	!! Outputs
	REAL(KIND=8)			:: den	!! density of LB protons

	den = 0.095D6 ! [1/m^3]

END SUBROUTINE LB_ion_density

