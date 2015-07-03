
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! kras_mars_density is neutral density profile
! given by krasnopolsky JGR, 107, E12, 5128 (2002)
! figure 2a for solar mean conditions. data points
! from the figure were fit with a quadratic and 
! high altitudes were fit with line. intersection
! of two points takes place when values are equal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE kras_mars_density(atom, h, den)

	IMPLICIT NONE

	!! Inputs
	INTEGER						:: atom
	REAL(KIND=8)			:: h				! [m]	

	!! Outputs
	REAL(KIND=8)			:: den			! [1/m^3]

	!! Internal
	REAL(KIND=8)			:: q, h_km, PL1, PL2, PL3, PH1, PH2, z0, A, B, C

	h_km = h/1000.0D0

	! CO2
	IF (atom .EQ. 44) THEN
		PL1 =  1.0648D-4
		PL2 = -7.4214D-2	
		PL3 =  18.493D0
		PH1 = -3.4949D-2
		PH2 = 15.786D0
	! H
	ELSE IF (atom .EQ. 1) THEN
		PL1 = 4.9289D-5
		PL2 = -2.6595D-2
		PL3 = 8.5637D0
		PH1 = -3.0103D-3
		PH2 = 5.9031D0
	! He
	ELSE IF (atom .EQ. 4) THEN	
		PL1 = 9.0656D-5
		PL2 = -4.5914D-2
		PL3 = 10.698D0
		PH1 = -4.7712D-3
		PH2 = 6.4314D0
	ELSE IF (atom .EQ. 16) THEN
		PL1 = 9.0065D-5
		PL2 = -5.292D-2
		PL3 = 14.208D0
		PH1 = -1.165D-2
		PH2 = 9.7959D0	
	ELSE
		WRITE(*,*) "atom of mass: ", atom, &
		&				" is unknown to me and thus has no density profile on Mars!!!"
		STOP	
	END IF

	A  = PL1
	B  = PL2 - PH1
	C  = PL3 - PH2
	z0 = (-B + SQRT(B*B - 4.0D0*A*C))/(2.0D0*A)

	IF (h_km < z0) THEN
		q   = PL1*h_km*h_km + PL2*h_km + PL3
		den = 10.0D0**q
	ELSE
		q   = PH1*h_km + PH2 
		den = 10.0D0**q
	END IF	

	!! convert density from [1/cm^3] to [1/m^3]
	den = den*1.0D6

END SUBROUTINE kras_mars_density

!########################################################
!########################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! mars_density holds density data for Mars as given
! by Figure 3 Fox et al. JGR (2009). 
! atom gives a string of the particular atom which
! which density is requested. 
! h is the height above the Martian surface which 
! the density is being requested at. 
! den returns the density in SI [m^-3]
! 
! Currently all densities are for Low Solar Activity!!! 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE mars_density(atom, h, den)

	IMPLICIT NONE

	!! Inputs
	INTEGER						:: atom		! atomic mass [integer]
	REAL(KIND=8)			:: h			! [m]

	!! Outputs
	REAL(KIND=8)			:: den		! [1/m^3]

	!! Internal
	INTEGER						:: i
	REAL(KIND=8)			:: n0, h0, HH, h_km  

	h_km = h/1000.0D0

	IF (atom .EQ. 4) THEN
		!! Helium Denisty
		n0 = 1.0D6 !! [1/m^3]
		IF (h_km .GE. 125.0D0) THEN
			h0 = 1876.9D0 !! [km] 
			HH = -124.27D0
		ELSE
			h0 = 279.99D0 !! [km]
			HH = -10.857D0
		END IF
		den = n0*EXP( (h_km - h0)/HH )

	ELSE IF (atom .EQ. 1) THEN
		!! Hydrogen Density
		n0 = 1.0D6 !! [1/m^3]	
		IF (h_km .GE. 120.0D0) THEN
			h0 = 6420.8D0 !! [km] 
			HH = -493.26D0
		ELSE
			h0 = 300.0D0	!! [km]
			HH = -13.029D0
		END IF 
		den = n0*EXP( (h_km - h0)/HH )

	ELSE IF (atom .EQ. 16) THEN
		!! Oxygen Density
		n0 = 1.0D6 !! [1/m^3]
		IF (h_km .GE. 140.0D0) THEN
			h0 = 790.01D0	!! [km]
			HH = -34.744D0
		ELSE
			h0 = 310.0D0	!! [km]
			HH = -8.6859D0
		END IF
		den = n0*EXP( (h_km - h0)/HH )

	ELSE IF (atom .EQ. 44) THEN
		n0 = 1.0D12		!! [cm^-3]
		h0 = 102.0D0	!! [km]
		HH = 7.8D0		!! [km]
		n0 = n0*1.0D6 !! convert to [m^-3]
		h0 = h0*1.0D3 !! convert to [m]
		HH = HH*1.0D3 !! convert to [m]
		den = n0*EXP( -(h-h0)/HH )

	ELSE
		WRITE(*,*) "atom of mass: ", atom, &
		&				" is unknown and thus has no density profile on Mars!!!"
		STOP	
	END IF
	

END SUBROUTINE mars_density

