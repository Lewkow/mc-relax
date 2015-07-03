
MODULE physics_constants

	REAL(KIND=8), PARAMETER		:: TOK      = 315776.177845143					! AU to K 
	REAL(KIND=8), PARAMETER		:: TOAMU    = 1822.888479779						! AU to amu
	REAL(KIND=8), PARAMETER		:: AMUTOKG  = 9.10938291E-31						! amu to kg
	REAL(KIND=8), PARAMETER		:: UTOKG		= 1.66053886E-27						! AU to kg
	REAL(KIND=8), PARAMETER		:: TOEV     = 27.21138386								! AU to eV
	REAL(KIND=8), PARAMETER		:: EVTOJ    = 1.602176487E-19						! eV to J
	REAL(KIND=8), PARAMETER		:: EVTOK		= 11604.505									! eV to K
	REAL(KIND=8), PARAMETER		:: PI       = 3.141592653589793238462
	REAL(KIND=8), PARAMETER		:: TODEG    = 180.0D0/PI								! Rad to Deg
	REAL(KIND=8), PARAMETER		:: TOA      = 0.529177249								! Bohr to Angstroms 
	REAL(KIND=8), PARAMETER		:: BOHRTOM  = 5.29177249D-11						! Bohr to Meters
	REAL(KIND=8), PARAMETER		:: CMTOBOHR = 188972598.85789						! CM to Bohr
	REAL(KIND=8), PARAMETER		:: CMTOM		= 1.0D0/100.0D0							! CM to M 
!	REAL(KIND=8), PARAMETER		:: KB       = 8.6173324E-5/TOEV					! Boltzmann constant in AU/K 
	REAL(KIND=8), PARAMETER		:: KB       = 1.3806488E-23							! Boltzmann constant in J/K 
	REAL(KIND=8), PARAMETER		:: MTOAU		= 6.68458134E-12						! Meters to Astronomical Units	
	REAL(KIND=8), PARAMETER		:: PCTOM		= 3.08567758E16							! Parsecs to Meters
	REAL(KIND=8), PARAMETER		:: PCTOAU		= PCTOM*MTOAU								! Parsecs to Meters
	REAL(KIND=8), PARAMETER		:: MARS_R   = 3396.0D3	  							! Radius of Mars [m]
	REAL(KIND=8), PARAMETER		:: MARS_M   = 6.41693E23								! Mass of Mars [kg]
	REAL(KIND=8), PARAMETER		:: GRAV_G   = 6.673E-11									! Gravitational Constant [m^3/kg/s^2]
	REAL(KIND=8), PARAMETER		:: Mass_p   = 1.007276466812D0					! Mass of proton [amu]
	REAL(KIND=8), PARAMETER		:: Mass_n		= 1.00866491600D0						! Mass of neutron [amu]
	REAL(KIND=8), PARAMETER   :: Mass_e		= 5.4857990946D-4						! Mass of electron [amu]
	REAL(KIND=8), PARAMETER		:: AUtoLY   = 15.813D-6									! AU to LY converstion
	REAL(KIND=8), PARAMETER		:: eps0     = 8.854187817D-12						! permiativity of free space
	REAL(KIND=8), PARAMETER   :: qe       = 1.60217657D-19						! elementrary charge [C]

END MODULE physics_constants
