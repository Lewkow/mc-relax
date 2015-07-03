
MODULE rand_seed

	INTEGER(KIND=4),PARAMETER		:: modulus   = 2147483647
	INTEGER(KIND=4),PARAMETER		:: small_lag = 97
	INTEGER(KIND=4),PARAMETER		:: large_lag = 127
	INTEGER(KIND=4),PARAMETER		:: cbit      = 22

	INTEGER(KIND=4)							:: register(large_lag)
	INTEGER(KIND=4)							:: lag1 
	INTEGER(KIND=4)							:: lag2
	INTEGER(KIND=4)							:: seed

	REAL(KIND=8)								:: rmult

END MODULE

