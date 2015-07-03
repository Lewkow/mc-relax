
MODULE inputs

	REAL(KIND=8)				:: EI														! Initial Energy
	REAL(KIND=8)				:: EF														! Final Energy
	REAL(KIND=8)				:: DE														! Difference in Energy for each Calc 
	REAL(KIND=8)				:: MP 													! Projectile Mass
	REAL(KIND=8)				:: MT 													! Target Mass
	REAL(KIND=8)				:: MU 													! Reduced Mass

	INTEGER							:: N_Max												! Maximum number of collisions for a particle
	INTEGER							:: NE														! Number of Energies to Compute
	INTEGER							:: N_Part

	INTEGER							:: DO_LISM											! Perform LISM calculation if = 1
	INTEGER							:: DO_planet										! Perform planet calculation if = 1
	INTEGER							:: DO_ENA_production						! Perform ENA production calculation if = 1
	INTEGER							:: DO_Average_Scattering_Angle	! Perform Average Scattering Angle Calculation if = 1
	INTEGER							:: DO_Escape_Trans
	INTEGER							:: DO_Escape_MC
	
END MODULE inputs

