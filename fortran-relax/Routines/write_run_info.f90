
SUBROUTINE write_run_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Writes to a file in ../Data/ which holds all the info
! on a given run to make data analysis easier. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE planet

	IMPLICIT NONE

	OPEN(UNIT=200,FILE="../Data/planet_3d_run_info.dat",ACCESS="APPEND")

	WRITE(200,'(A,A)')      "Projectile Type:     ", Proj
	WRITE(200,'(A,I10)')    "Number MC Particles: ", N_Part 
	WRITE(200,'(A,ES10.2)') "Thermal Energy [eV]: ", E_Therm
	WRITE(200,'(A,ES10.2)') "Initial height [km]: ", h_0/1000.0D0 
	WRITE(200,'(A,ES10.2)') "Escape  height [km]: ", high 
	WRITE(200,'(A,ES10.2)') "Solar Zeneith Angle: ", SZA 
	

	CLOSE(200)

END SUBROUTINE write_run_info

