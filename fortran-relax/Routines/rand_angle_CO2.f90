!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gives and random angle for H+CO2 as
! given by Kollis et al. (2000). 
!
! T_rand(RND) = ( RND*(T_mx^a - T_mn^a) + T_mn^a )^(1/a)
! RND: random number [0-1]
! T_mx: maximum scattering angle: PI for H+CO2
! T_mn: minimum scattering angle: 1 deg for H+CO2 (Kollis et al. )
! a: -0.666
!
! T_rand is given in the cm frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE rand_angle_CO2( angle )
USE physics_constants, ONLY : PI
USE planet, ONLY : MT, MP

	IMPLICIT NONE

	!! Output
	REAL(KIND=8)		:: angle

	!! Internal
	REAL(KIND=8)		:: a, T_mn, T_mx, r, p, ang, C, C1, C2, C3
	REAL(KIND=8)		:: lfg

	C    = PI/180.0D0
	a    = -0.666D0
	T_mn = PI/180.0D0
	T_mx = PI
!	T_mn = 1.0D0 
!	T_mx = 180.0D0
	r    = lfg()

	C1   = T_mn**a
	C2   = T_mx**a - C1
	C3   = r*C2 + C1
!	ang  = C3**(1.0D0/a)

	ang  = ( r*(T_mx**a - T_mn**a) + T_mn**a )**(1.0D0/a)
!	ang  = ang*C

	!! convert angle from lab to CM frame
	p = MP/MT
	
	angle = ACOS( -p + p*COS(ang)**2.0 + COS(ang)*SQRT(	1.0D0 - p*p*SIN(ang)**2.0) )


END SUBROUTINE rand_angle_CO2
