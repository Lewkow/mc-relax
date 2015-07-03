!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calculate ena production through transparency integration of 
! atmosphere of the planet Mars
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mars_ena_trans
	use ena
	use planet

	implicit none

	real(kind=8)		:: alt_h, alt_l, f, z, dz, den
	integer					:: n, i

	ATMOSPHERE = 2
	WRITE(*,*) 'CALCULATING MARS ENA PRODUCTION with atmosphere model: ', ATMOSPHERE
	CALL read_mars_density

	n  	  = 1000			!! number of altitude bins
	alt_l = 60.0d3		!! low altitude limit
	alt_h	= 800.0d3	!! high altitude limit
	dz   	= (alt_h - alt_l)/real(n-1)

	open(unit=100,file="../Data/ena_production_hp.dat",access="append")
	open(unit=200,file="../Data/ena_production_hepp.dat",access="append")
	
	!! write charge exchange cross sections to file
	call write_cx_cs
	
	do i=1,n
		z = alt_l + (i-1)*dz
		call sw_hp_ena_prod(z,f)
		write(100,*) z/1000.0d0, f
		call sw_hepp_ena_prod(z,f)
		write(200,*) z/1000.0d0, f
	end do

	close(100)
	close(200)

end subroutine mars_ena_trans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calculate ena production through transparency integration of 
! the local interstellar medium 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lism_ena_trans
	use ENA

	implicit none

	real(kind=8)		:: r_start, r_end, f, r, dr, den, AU
	integer					:: n, i

	AU    	= 1.5D11				!! AU in [m]
	n  	  	= 50000					!! number of altitude bins
	r_start = AU						!! integration starting point
	r_end		= 1000.0D0*AU		!! high altitude limit
	dr   		= (r_end-r_start)/real(n-1)

	open(unit=110,file="../Data/lism_production_hp.dat",access="append")
	open(unit=210,file="../Data/lism_production_hep.dat",access="append")
	
	!! write charge exchange cross sections to file
!	call write_cx_cs
	
	do i=1,n
		r = r_start + (i-1)*dr
		call lism_sw_hp_ena_prod(r,f)
		write(110,*) r/AU, f
		call lism_sw_hep_ena_prod(r,f)
		write(210,*) r/AU, f
	end do

	close(110)
	close(210)

end subroutine lism_ena_trans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Proton Component of solar wind
!
! Calculates SW flux using the following equation
! 
! j_SW = N_SW * V_SW * exp[ - sum_i int_{inf}^{r} (n_{i} * CS_{i} ]
!
! N_SW   = Density of SW at r [1/m^3]
! V_SW   = Velocity of SW [m/s] 
! n_{i}  = Density of ith neutral gas at r [1/m^3]
! CS_{i} = Charge exchange cross section for energy of SW [m^2]
!
! From Kallio et al. JGR, 1997, 102, A10. 
!
! Inputs: r -> distance from star [m]
!
! Outputs: f -> SW flux with CX losses at distance r [1/m^2/s]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lism_SW_Hp_ENA_prod(r,f)
  USE ENA
  USE physics_constants

  IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)  :: r 	!! distance from star [m]

  !! Outputs
  REAL(KIND=8)  :: f  !! Density of newly created neutrals [1/m^3/s]

  !! Internal
  REAL(KIND=8)	:: den_H, den_O, den_O2, den_He, den_CO2
  REAL(KIND=8)  :: C1, c2, C3, E, cx_H, cx_O, cx_CO2, cx_O2, cx_He
  REAL(KIND=8)  :: AU, now, tot, dr, x, r_start, bot
  CHARACTER     :: ion*10, targ*10
  INTEGER       :: N, i, j

	IF (r .LT. lism_start) WRITE(*,*) "Need to increase lism_start to complete calculation"


	AU      = 1.5D11					!! AU [m]
  N   		= 10000     			!! Number of integration steps for total interval
  r_start = lism_start      !! starting distance from star [m]
  dr  		= (r-r_start)/REAL(N-1)

	!! get density of hydrogen [m]
	CALL LISM_density( r, den_H )	

  ion   = 'H_p'
  L_ion = LEN_TRIM(ion)
  E     = 0.5D0*Mass_p*TOAMU*AMUTOKG*V_sw*V_sw
	E     = E/EVTOJ !! convert from J to eV

  !! Get CX CS for each collision at current energy
  targ   = 'H'
  L_targ = LEN_TRIM(targ)
  CALL cx_cross_sections( ion, targ, 1, E, cx_H )

  tot = 0.0D0

  DO i=1,N
    now = ( den_H*cx_H )*dr
    tot = tot + now
  END DO

	C1 = H_sw*N_sw*V_sw*(AU/r)**2
	C2 = EXP(-tot)
	C3 =  den_H*cx_H 
	f  = C1*C2*C3

END SUBROUTINE lism_SW_Hp_ENA_prod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Proton Component of solar wind
!
! Calculates SW flux using the following equation
! 
! j_SW = N_SW * V_SW * exp[ - sum_i int_{inf}^{r} (n_{i} * CS_{i} ]
!
! N_SW   = Density of SW at r [1/m^3]
! V_SW   = Velocity of SW [m/s] 
! n_{i}  = Density of ith neutral gas at r [1/m^3]
! CS_{i} = Charge exchange cross section for energy of SW [m^2]
!
! From Kallio et al. JGR, 1997, 102, A10. 
!
! Inputs: r -> distance from star [m]
!
! Outputs: f -> SW flux with CX losses at distance r [1/m^2/s]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lism_SW_Hep_ENA_prod(r,f)
  USE ENA
  USE physics_constants

  IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)  :: r 	!! distance from star [m]

  !! Outputs
  REAL(KIND=8)  :: f  !! Density of newly created neutrals [1/m^3/s]

  !! Internal
  REAL(KIND=8)	:: den_H, den_O, den_O2, den_He, den_CO2
  REAL(KIND=8)  :: C1, c2, C3, E, cx_H, cx_O, cx_CO2, cx_O2, cx_He
  REAL(KIND=8)  :: AU, now, tot, dr, x, r_start, bot
  CHARACTER     :: ion*10, targ*10
  INTEGER       :: N, i, j

	IF (r .LT. lism_start) WRITE(*,*) "Need to increase lism_start to complete calculation"

	AU      = 1.5D11					!! AU [m]
  N   		= 10000     			!! Number of integration steps for total interval
  r_start = lism_start      !! starting distance from star [m]
  dr  		= (r-r_start)/REAL(N-1)

	!! get density of hydrogen [m]
	CALL LISM_density( r, den_H )	

  ion   = 'He_p'
  L_ion = LEN_TRIM(ion)
  E     = 2.0D0*Mass_p*TOAMU*AMUTOKG*V_sw*V_sw
	E     = E/EVTOJ !! convert from J to eV

  !! Get CX CS for each collision at current energy
  targ   = 'H'
  L_targ = LEN_TRIM(targ)
  CALL cx_cross_sections( ion, targ, 1, E, cx_H )

  tot = 0.0D0

  DO i=1,N
    now = ( den_H*cx_H )*dr
    tot = tot + now
  END DO

	C1 = He_sw*N_sw*V_sw*(AU/r)**2
	C2 = EXP(-tot)
	C3 = den_H*cx_H 
	f  = C1*C2*C3

END SUBROUTINE lism_SW_Hep_ENA_prod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Proton Component of solar wind
!
! Calculates SW flux using the following equation
! 
! j_SW = N_SW * V_SW * exp[ - sum_i int_{inf}^{r} (n_{i} * CS_{i} ]
!
! N_SW   = Density of SW at Mars [1/m^3]
! V_SW   = Velocity of SW [m/s] 
! n_{i}  = Density of ith neutral gas at r [1/m^3]
! CS_{i} = Charge exchange cross section for energy of SW [m^2]
!
! From Kallio et al. JGR, 1997, 102, A10. 
!
! Inputs: z -> altitude above planet [m]
!
! Outputs: f -> SW flux with CX losses at altitude r [1/m^2/s]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SW_Hp_ENA_prod(z,f)
  USE ENA
  USE physics_constants
	USE planet

  IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)  :: z 	!! Altitude above surface [m]

  !! Outputs
  REAL(KIND=8)  :: f  !! Density of newly created neutrals [1/m^3/s]

  !! Internal
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: den_H, den_O, den_O2, den_He, den_CO2
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: den_H2, den_N2, den_CO, den_Ar
  REAL(KIND=8)  :: C1, c2, C3, E, cx_H, cx_O, cx_CO2, cx_O2, cx_He, cx_CO, cx_Ar
  REAL(KIND=8)  :: now, tot, dx, x, r, top, bot, VSW, cx_H2, cx_N2
  CHARACTER     :: ion*10
  INTEGER       :: N, i, j

	IF (z .GT. z_in) WRITE(*,*) "Need to increase z_in to complete calculation"

  N   = 10000     !! Number of integration steps for total interval
  top = z_in      !! Top of atmosphere [m]
  r   = z   			!! Current height in atmosphere [m]
  dx  = (top-r)/REAL(N-1)

  ALLOCATE( den_H(N), den_O(N), den_He(N), den_CO2(N) )
	ALLOCATE( den_H2(N), den_N2(N), den_CO(N), den_Ar(N) )

  ion   = 'H_p'
  L_ion = LEN_TRIM(ion)
	E     = 1.0D3
  VSW   = SQRT(E*EVTOJ/(2.0D0*Mass_p*TOAMU*AMUTOKG))

  !! Get CX CS for each collision at current energy
  CALL cx_cross_sections( ion, 'H  ', 1, E, cx_H )
  CALL cx_cross_sections( ion, 'O  ', 1, E, cx_O )
  CALL cx_cross_sections( ion, 'He ', 1, E, cx_He )
  CALL cx_cross_sections( ion, 'CO2', 1, E, cx_CO2 )
  CALL cx_cross_sections( ion, 'H2 ', 1, E, cx_H2 )
  CALL cx_cross_sections( ion, 'CO2', 1, E, cx_N2 )
  CALL cx_cross_sections( ion, 'CO ', 1, E, cx_CO )
  CALL cx_cross_sections( ion, 'Ar ', 1, E, cx_Ar )

!	WRITE(*,'(A,8ES10.2)') 'H He O Ar H2 N2 CO CO2: ', cx_H, cx_He, cx_O, cx_Ar, cx_H2, cx_N2, cx_CO, cx_CO2

  DO i=1,N
    x = top - (i-1)*dx
    CALL mars_table_density( 'H  ', x, den_H(i) )
    CALL mars_table_density( 'He ', x, den_He(i) )
    CALL mars_table_density( 'O  ', x, den_O(i) )
    CALL mars_table_density( 'Ar ', x, den_Ar(i) )
    CALL mars_table_density( 'H2 ', x, den_H2(i) )
    CALL mars_table_density( 'N2 ', x, den_N2(i) )
    CALL mars_table_density( 'CO ', x, den_CO(i) )
    CALL mars_table_density( 'CO2', x, den_CO2(i) )
  END DO

  tot = 0.0D0

  DO i=1,N
    now = ( den_H(i)*cx_H + den_O(i)*cx_O + den_He(i)*cx_He + den_CO2(i)*cx_CO2 + &
		&       den_Ar(i)*cx_Ar + den_H2(i)*cx_H2 + den_N2(i)*cx_N2 + den_CO(i)*cx_CO)*dx
!		WRITE(*,*) 'H   - den: cx: ', den_H(i), cx_H
!		WRITE(*,*) 'He  - den: cx: ', den_He(i), cx_He
!		WRITE(*,*) 'O   - den: cx: ', den_O(i), cx_O
!		WRITE(*,*) 'Ar  - den: cx: ', den_Ar(i), cx_Ar
!		WRITE(*,*) 'H2  - den: cx: ', den_H2(i), cx_H2
!		WRITE(*,*) 'N2  - den: cx: ', den_N2(i), cx_N2
!		WRITE(*,*) 'CO  - den: cx: ', den_CO(i), cx_CO
!		WRITE(*,*) 'CO2 - den: cx: ', den_CO2(i), cx_CO2
!		WRITE(*,*) 'now: tot: ', now, tot
    tot = tot + now
  END DO


  CALL mars_table_density( 'H  ', r, den_H(1) )
  CALL mars_table_density( 'He ', r, den_He(1) )
  CALL mars_table_density( 'O  ', r, den_O(1) )
  CALL mars_table_density( 'Ar ', r, den_Ar(1) )
  CALL mars_table_density( 'H2 ', r, den_H2(1) )
  CALL mars_table_density( 'N2 ', r, den_N2(1) )
  CALL mars_table_density( 'CO ', r, den_CO(1) )
  CALL mars_table_density( 'CO2', r, den_CO2(1) )

	C1 = H_sw*N_sw*VSW
	C2 = EXP(-tot)
	C3 = (den_H(1)*cx_H + den_O(1)*cx_O + den_He(1)*cx_He + den_CO2(1)*cx_CO2 + &
  &     den_H2(1)*cx_H2 + den_N2(1)*cx_N2 + den_CO(1)*cx_CO + den_Ar(1)*cx_Ar)
!	WRITE(*,*) 'z: C1: C2: C3: ', r, C1, C2, C3
	f  = C1*C2*C3

  DEALLOCATE( den_H, den_O, den_He, den_CO2 )
	DEALLOCATE( den_H2, den_N2, den_CO, den_Ar ) 

END SUBROUTINE SW_Hp_ENA_prod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE rand_SW_Hp_ENA_prod(E,z)
  USE ENA
  USE physics_constants
	USE planet, ONLY : mars_table_density

  IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)  :: E 	!! Energy of proton [eV]

  !! Outputs
  REAL(KIND=8)  :: z 	!! Altitude above surface [m]

  !! Internal
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: den_H, den_O, den_O2, den_He, den_CO2
  REAL(KIND=8)  		:: C1, c2, C3, cx_H, cx_O, cx_CO2, cx_O2, cx_He
  REAL(KIND=8)  		:: now, tot, dx, x, r, top, bot, f
  CHARACTER(Len=10) :: ion, targ
  INTEGER       		:: N, i, j

	IF (z .GT. z_in) WRITE(*,*) "Need to increase z_in to complete calculation"

  N   = 10000     !! Number of integration steps for total interval
  top = z_in      !! Top of atmosphere [m]
  r   = z   			!! Current height in atmosphere [m]
  dx  = (top-r)/REAL(N-1)

  ALLOCATE( den_H(N), den_O(N), den_He(N), den_CO2(N) )

  ion   = 'H_p'
  L_ion = LEN_TRIM(ion)
	E     = 1.0D3
!  E     = 0.5D0*Mass_p*TOAMU*AMUTOKG*V_sw*V_sw
!	E     = E/EVTOJ !! convert from J to eV

  !! Get CX CS for each collision at current energy
  targ   = 'H'
  L_targ = LEN_TRIM(targ)
  CALL cx_cross_sections( ion, targ, 1, E, cx_H )
  targ   = 'O'
  L_targ = LEN_TRIM(targ)
  CALL cx_cross_sections( ion, targ, 1, E, cx_O )
  targ   = 'He'
  L_targ = LEN_TRIM(targ)
  CALL cx_cross_sections( ion, targ, 1, E, cx_He )
	!! Use O2 cross section with CO2 density for now
  targ   = 'CO2'
	L_targ = LEN_TRIM(targ)
  CALL cx_cross_sections( ion, targ, 1, E, cx_CO2 )

  DO i=1,N
    x = top - (i-1)*dx
		CALL mars_table_density( 'H  ', x, den_H(i)   )
		CALL mars_table_density( 'He ', x, den_He(i)  )
		CALL mars_table_density( 'O  ', x, den_O(i)   )
		CALL mars_table_density( 'CO2', x, den_CO2(i) )
  END DO

  tot = 0.0D0

  DO i=1,N
    now = ( den_H(i)*cx_H + den_O(i)*cx_O + den_He(i)*cx_He + den_CO2(i)*cx_CO2 )*dx
    tot = tot + now
  END DO

	CALL mars_table_density( 'H  ', r, den_H(1)   )
	CALL mars_table_density( 'He ', r, den_He(1)  )
	CALL mars_table_density( 'O  ', r, den_O(1)   )
	CALL mars_table_density( 'CO2', r, den_CO2(1) )

	C1 = H_sw*N_sw*V_sw
	C2 = EXP(-tot)
	C3 = (den_H(1)*cx_H + den_O(1)*cx_O + den_He(1)*cx_He + den_CO2(1)*cx_CO2)
	f  = C1*C2*C3

  DEALLOCATE( den_H, den_O, den_He, den_CO2 )

END SUBROUTINE rand_SW_Hp_ENA_prod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Alpha Particle Component of solar wind
!
! Calculates SW flux using the following equation
! 
! j_SW = N_SW * V_SW * exp[ - sum_i int_{inf}^{r} (n_{i} * CS_{i} ]
!
! N_SW   = Density of SW at Mars [1/m^3]
! V_SW   = Velocity of SW [m/s] 
! n_{i}  = Density of ith neutral gas at r [1/m^3]
! CS_{i} = Charge exchange cross section for energy of SW [m^2]
!
! From Kallio et al. JGR, 1997, 102, A10. 
!
! Inputs: z -> altitude above planet [m]
!
! Outputs: f -> SW flux with CX losses at altitude r [1/m^2/s]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SW_Hepp_ENA_prod(z,f)
  USE ENA
  USE physics_constants
	USE planet

  IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)  :: z 	!! Altitude above surface [m]

  !! Outputs
  REAL(KIND=8)  :: f  !! Density of newly created neutrals [1/m^3/s]

  !! Internal
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: den_H, den_He, den_O, den_Ar, den_CO2
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: den_O2, den_N2, den_CO, den_H2
  REAL(KIND=8)  :: C1, c2, C3, E, cx_H, cx_Ar, cx_O, cx_CO2, cx_O2, cx_He
	REAL(KIND=8)	:: cx_H2, cx_N2, cx_CO
  REAL(KIND=8)  :: now, tot, dx, x, r, top, bot, VSW
  CHARACTER     :: ion*10
  INTEGER       :: N, i, j

	IF (z .GT. z_in) WRITE(*,*) "Need to increase z_in to complete calculation"

  N   = 10000     !! Number of integration steps for total interval
  top = z_in      !! Top of atmosphere [m]
  r   = z   			!! Current height in atmosphere [m]
  dx  = (top-r)/REAL(N-1)

  ALLOCATE( den_H(N), den_O(N), den_He(N), den_CO2(N) )
  ALLOCATE( den_H2(N), den_CO(N), den_Ar(N), den_N2(N) )

  ion   = 'He_pp'
  L_ion = LEN_TRIM(ion)
	E     = 4.0D3
	VSW   = SQRT(E*EVTOJ/(2.0D0*Mass_p*TOAMU*AMUTOKG))
!  E     = (Mass_n+Mass_p)*TOAMU*AMUTOKG*V_sw*V_sw
!	E     = E/EVTOJ !! convert from J to eV

  !! Get CX CS for each collision at current energy
  CALL cx_cross_sections( ion, 'He ', 2, E, cx_He )
	CALL cx_cross_sections( ion, 'CO2', 2, E, cx_CO2 ) 

	!! for now set these CX CS equal as we do not have a He^+2 + O CX CS
	cx_O  = cx_CO2
	cx_Ar = cx_CO2
	cx_N2 = cx_CO2
	cx_H2 = cx_CO2
	cx_CO = cx_CO2

  DO i=1,N
    x = top - (i-1)*dx
		CALL mars_table_density( 'He ', x, den_He(i) )
		CALL mars_table_density( 'O  ', x, den_O(i) )
		CALL mars_table_density( 'Ar ', x, den_Ar(i) )
		CALL mars_table_density( 'H2 ', x, den_H2(i) )
		CALL mars_table_density( 'N2 ', x, den_N2(i) )
		CALL mars_table_density( 'CO ', x, den_CO(i) )
		CALL mars_table_density( 'CO2', x, den_CO2(i) )
  END DO

  tot = 0.0D0

	!! For now use He^2+ + He -> He + He^2+ for all CX CS

  DO i=1,N
    now = ( den_O(i)*cx_O + den_He(i)*cx_He + den_CO2(i)*cx_CO2 + &
		& den_H2(i)*cx_H2 + den_N2(i)*cx_N2 + den_CO(i)*cx_CO + den_Ar(i)*cx_Ar)*dx
    tot = tot + now
  END DO

	CALL mars_table_density( 'He ', r, den_He(1) )
	CALL mars_table_density( 'O  ', r, den_O(1) )
	CALL mars_table_density( 'Ar ', r, den_Ar(1) )
	CALL mars_table_density( 'H2 ', r, den_H2(1) )
	CALL mars_table_density( 'N2 ', r, den_N2(1) )
	CALL mars_table_density( 'CO ', r, den_CO(1) )
	CALL mars_table_density( 'CO2', r, den_CO2(1) )

	C1 = He_sw*N_sw*VSW
	C2 = EXP(-tot)
	C3 = (den_O(1)*cx_O + den_He(1)*cx_He + den_CO2(1)*cx_CO2 + &
	&     den_H2(1)*cx_H2 + den_N2(1)*cx_N2 + den_CO(1)*cx_CO + den_Ar(1)*cx_Ar)
	f  = C1*C2*C3

  DEALLOCATE( den_H, den_O, den_He, den_CO2 )
	DEALLOCATE( den_Ar, den_CO, den_N2, den_H2 ) 

END SUBROUTINE SW_Hepp_ENA_prod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

