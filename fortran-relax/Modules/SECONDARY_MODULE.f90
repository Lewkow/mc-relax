
MODULE secondary

	INTEGER, PARAMETER											:: MAX_SHA  = 2000		! Maximum number of SHAs per ENA
	
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: SHA_E							! Energy array for SHAs
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: SHA_V							! Unit velocity array for SHAs
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: SHA_R							! Location array for SHAs
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: r_zone_SHA_E     	! height SHA E array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: r_zone_SHA_C     	! height SHA count array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: r_prod_SHA_E     	! height SHA E array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: r_prod_SHA_C     	! height SHA prod count array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: r_zone_SHA_dE     	! height SHA E array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: r_zone_SHA_dC     	! height SHA count array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: Root_r_zone_SHA_E	! root height SHA E array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: Root_r_zone_SHA_C  ! root height SHA count array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: Root_r_prod_SHA_E	! root height SHA E array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: Root_r_prod_SHA_C  ! root height SHA count array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: Root_r_zone_SHA_dE	! root height SHA E array
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)   :: Root_r_zone_SHA_dC ! root height SHA count array

	INTEGER																	:: My_SHA_count

	CONTAINS
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	SUBROUTINE allocate_sha

		ALLOCATE( SHA_E(MAX_SHA), SHA_V(MAX_SHA,3), SHA_R(MAX_SHA,3) )

	END SUBROUTINE allocate_sha

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE clean_sha

		SHA_E(:) 	 	= 0.0D0
		SHA_V(:,1) 	= 0.0D0
		SHA_V(:,2) 	= 0.0D0
		SHA_V(:,3) 	= 0.0D0
		SHA_R(:,1) 	= 0.0D0
		SHA_R(:,2) 	= 0.0D0
		SHA_R(:,3) 	= 0.0D0

	END SUBROUTINE clean_sha

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE free_sha

		DEALLOCATE( SHA_E, SHA_V, SHA_R )

	END SUBROUTINE free_sha

END MODULE secondary

