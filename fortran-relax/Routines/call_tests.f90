!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Specify which MC tests to run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE call_tests
	USE inputs, ONLY : DO_LISM, DO_planet, DO_ENA_production, DO_Average_Scattering_Angle, DO_Escape_Trans, DO_Escape_MC
	USE planet, ONLY : test_mars_table_density
	USE escape_trans

	IMPLICIT NONE

	IF (DO_LISM   .EQ. 1)         					CALL LISM_ENA_displacement
	IF (DO_planet .EQ. 1 .OR. DO_Escape_MC .EQ. 1) CALL planet_onestep_3d
	IF (DO_Average_Scattering_Angle .EQ. 1) CALL average_scattering_angle
	IF (DO_Escape_Trans .EQ. 1) 						CALL all_trans
	IF (DO_ENA_production .EQ. 1) THEN
		CALL lism_ENA_trans
		CALL mars_ENA_trans
	END IF

CALL test_Hp_H_tcs

!CALL write_cx_cs
	
!CALL test_mars_table_density

!CALL ion_planet_simulation

!CALL ion_transport_test
!CALL test_pd
!CALL test_universal_dcs
!CALL test_universal_tcs
!CALL uni_ang_test

END SUBROUTINE call_tests



