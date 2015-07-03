SUBROUTINE test_lin_rand_angle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Test lin_rand_angle routine by finding the
! average scattering angle and standard dev
! for all collision species and over all
! collision energies
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	INTEGER				:: N_En, N_An, i, j, k
	REAL(KIND=8)	:: Ang_tot, dE, E_in, E_fn
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: Ang, Ang_mean, Ang_std, E


	N_En = 1000
	N_An = 1000
	E_in = 0.01D0
	E_fn = 3.5D3
	dE   = (E_fn-E_in)/REAL(N_En-1)

	OPEN(UNIT=456, FILE="../Data/Average_Collision_Angles.dat", ACCESS="APPEND")
	ALLOCATE(Ang(N_An), Ang_mean(N_En), Ang_std(N_En), E(N_En))

	DO i=1,N_En
		E(i)		= E_in + (i-1)*dE
		Ang_tot = 0.0D0
		DO j=1,N_An
			CALL lin_rand_angle(E(i),'HeO  ',Ang(j))
			Ang_tot = Ang_tot + Ang(j)*180.0D0/PI
		END DO
		Ang_mean(i) = Ang_tot/REAL(N_An)
		Ang_tot = 0.0D0
		DO j=1,N_An
			Ang_tot = Ang_tot + (Ang(j)-Ang_mean(i))**2	
		END DO		
		Ang_std(i) = SQRT( (1.0D0/N_An)*Ang_tot )	

		WRITE(456,*) E(i), Ang_mean(i), Ang_std(i)	

	END DO

	DEALLOCATE(Ang,Ang_mean,Ang_std,E)
	CLOSE(456)

END SUBROUTINE test_lin_rand_angle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE lin_rand_angle(Enow,Coll_Type,SAng)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Find random angle from AngProb3D tables 
! | Energy (eV) | Angle (deg) | Probability |
! Use linear interpolation inbetween angle
! points. 
!
! Energy in Lab Frame!!!
! SAng   in CM  Frame!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants, ONLY : PI
	USE rand_seed
	USE tables
	USE universal
	USE inputs, ONLY : DO_Average_Scattering_Angle
	USE planet, ONLY : MP, MT
	
	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)			:: Enow 				!! [eV] in lab frame
	CHARACTER(LEN=5)	:: Coll_Type

	!! Outputs
	REAL(KIND=8)			:: SAng				 !! [rad] in CM frame

	!! Internal
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: Ang_now
	REAL(KIND=8)													:: lfg
	REAL(KIND=8)													:: M, r, Pnow, x0, y0, x1, y1 
	REAL(KIND=8)													:: rat, CM_theta, LB_theta, Ecm
	INTEGER																:: NA, DONE, j, E_index, k

	!-------------------------	
	! Get index for energy	
	!-------------------------
	!! convert energy from lab to CM
	IF (DO_Average_Scattering_Angle .EQ. 1) THEN 
		Ecm = Enow
	ELSE
		Ecm = Enow*MT/(MP+MT)	
	END IF

	CALL index_EN(Ecm, Coll_Type, E_index)

	IF ( TRIM(Coll_Type) .EQ. 'HeO') THEN
		NA = NUM_AN_HeO
	ELSE IF ( TRIM(Coll_Type) .EQ. 'HeH') THEN
		NA = NUM_AN_HeH
	ELSE IF ( TRIM(Coll_Type) .EQ. 'HeHe') THEN
		NA = NUM_AN_HeHe
	ELSE IF ( TRIM(Coll_Type) .EQ. 'HH') THEN
		NA = NUM_AN_HH
	ELSE ! Use Scaling CS
		NA = NUM_ANG
	END IF

	!! Allocate PD_now based on Number of Angles for a given collision type
	ALLOCATE( PD_now(NA), Ang_now(NA) )

	IF ( TRIM(Coll_Type) .EQ. 'HeO') THEN
		PD_now(:)  = F_PD_HeO(E_index,:)
		Ang_now(:) = F_ANGLE_HeO(:)
	END IF
	IF ( (TRIM(Coll_Type) .EQ. 'HeH') .OR. (TRIM(Coll_Type) .EQ. 'HHe')) THEN
		PD_now(:)  = F_PD_HeH(E_index,:)
		Ang_now(:) = F_ANGLE_HeH(:)
	END IF
	IF ( TRIM(Coll_Type) .EQ. 'HeHe') THEN
		PD_now(:)  = F_PD_HeHe(E_index,:)
		Ang_now(:) = F_ANGLE_HeHe(:)
	END IF
	IF ( TRIM(Coll_Type) .EQ. 'HH') THEN
		PD_now(:)  = F_PD_HH(E_index,:)
		Ang_now(:) = F_ANGLE_HH(:)
	END IF
	IF ( TRIM(Coll_Type) .EQ. 'HO') THEN
		PD_now(:)  = F_PD_X_O(E_index,:)
		Ang_now(:) = F_ANGLE(:)
	END IF
	IF ( (TRIM(Coll_Type) .EQ. 'HAr') .OR. (TRIM(Coll_Type) .EQ. 'HeAr') ) THEN
		PD_now(:)  = F_PD_X_Ar(E_index,:)
		Ang_now(:) = F_ANGLE(:)
	END IF
	IF ( (TRIM(Coll_Type) .EQ. 'HH2') .OR. (TRIM(Coll_Type) .EQ. 'HeH2') ) THEN
		PD_now(:)  = F_PD_X_H2(E_index,:)
		Ang_now(:) = F_ANGLE(:)
	END IF
	IF ( (TRIM(Coll_Type) .EQ. 'HN2') .OR. (TRIM(Coll_Type) .EQ. 'HeN2') ) THEN
		PD_now(:)  = F_PD_X_N2(E_index,:)
		Ang_now(:) = F_ANGLE(:)
	END IF
	IF ( (TRIM(Coll_Type) .EQ. 'HCO') .OR. (TRIM(Coll_Type) .EQ. 'HeCO') ) THEN
		PD_now(:)  = F_PD_X_CO(E_index,:)
		Ang_now(:) = F_ANGLE(:)
	END IF
	IF ( (TRIM(Coll_Type) .EQ. 'HCO2') .OR. (TRIM(Coll_Type) .EQ. 'HeCO2') ) THEN
		PD_now(:)  = F_PD_X_CO(E_index,:)
		Ang_now(:) = F_ANGLE(:)
	END IF

	!! Get max probability
	M = PD_now(NA)

	!! Get random number for angle
	r    = M*lfg()
	DONE = 0
	j    = 1
	
	DO WHILE (DONE .EQ. 0)
		Pnow = PD_now(j)
		IF ( (r .LT. Pnow) .OR. ((j+1) .GE. SIZE(PD_now(:))) ) THEN
			x0   			= PD_now(j-1)
			x1   			= Pnow
			y0   			= Ang_now(j-1)
			y1   			= Ang_now(j)
			CM_theta 	= y0 + ( (r-x0)*y1 - (r-x0)*y0 )/(x1-x0)
			DONE 			= 1
		ELSE
			j = j+1	
		END IF
	END DO

	DEALLOCATE(PD_now, Ang_now)

	!! Convert angle from deg to rads
	SAng = CM_theta*PI/180.0D0

END SUBROUTINE lin_rand_angle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE find_rand_angle(Enow,Coll_Type,SAng)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Find random angle from prob density from
! ASTROSCATT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants , ONLY : PI
	USE rand_seed
	USE tables
	USE mpi_info

	IMPLICIT NONE

	!-----------------
	! INPUTS	
	!-----------------

	REAL(KIND=8)		:: Enow
	CHARACTER(LEN=4):: Coll_Type

	!-----------------
	! OUTPUTS	
	!-----------------

	REAL(KIND=8)		:: SAng

	!-----------------
	! INTERNAL	
	!-----------------

	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: r1
	REAL(KIND=8)		:: r2
	REAL(KIND=8)		:: x
	REAL(KIND=8)		:: y
	REAL(KIND=8)		:: yr
	REAL(KIND=8)		:: M

	INTEGER					:: E_index
	INTEGER					:: k
	INTEGER					:: good_rand
	INTEGER					:: NA

	!-------------------------	
	! Get index for energy	
	!-------------------------	
	CALL index_EN(Enow, Coll_Type, E_index)

	IF ( Coll_Type .EQ. 'HeO ') NA = NUM_AN_HeO
	IF ( Coll_Type .EQ. 'HeH ') NA = NUM_AN_HeH
	IF ( Coll_Type .EQ. 'HeHe') NA = NUM_AN_HeHe

	!! Allocate PD_now based on Number of Angles for a given collision type
	ALLOCATE( PD_now(NA) )

	DO k=1,NA
		IF ( Coll_Type .EQ. 'HeO ') PD_now(k) = F_PD_HeO(E_index,k)
		IF ( Coll_Type .EQ. 'HeH ') PD_now(k) = F_PD_HeH(E_index,k)
		IF ( Coll_Type .EQ. 'HeHe') PD_now(k) = F_PD_HeHe(E_index,k)
	END DO

	!-------------------------	
	! Calculate M for Energy	
	!-------------------------	
	CALL max_pd(M,NA)

	!-------------------------	
	! Calculate Splines for Enow	
	!-------------------------	
	IF ( Coll_Type .EQ. 'HeO ') CALL spline(F_ANGLE_HeO,  PD_now, NA, 1.0E31, 1.0E31, SP_Y2_HeO)
	IF ( Coll_Type .EQ. 'HeH ') CALL spline(F_ANGLE_HeH,  PD_now, NA, 1.0E31, 1.0E31, SP_Y2_HeH)
	IF ( Coll_Type .EQ. 'HeHe') CALL spline(F_ANGLE_HeHe, PD_now, NA, 1.0E31, 1.0E31, SP_Y2_HeHe)

	good_rand = 0

	!-------------------------	
	! Pull from PDF for random angle	
	!-------------------------	
	DO WHILE (good_rand == 0)
		r1 = lfg()
		r2 = lfg()
		x  = r1*PI	

		!-------------------------	
		! Calculate Splines for y(x)	
		!-------------------------	
		IF ( Coll_Type .EQ. 'HeO ' ) CALL splint(F_ANGLE_HeO,  PD_now, SP_Y2_HeO, NA, x*180.0D0/PI, y)
		IF ( Coll_Type .EQ. 'HeH ' ) CALL splint(F_ANGLE_HeH,  PD_now, SP_Y2_HeH, NA, x*180.0D0/PI, y)
		IF ( Coll_Type .EQ. 'HeHe' ) CALL splint(F_ANGLE_HeHe, PD_now, SP_Y2_HeHe, NA, x*180.0D0/PI, y)
		
		yr = r2*M
	
		IF (yr .LE. y) THEN
			good_rand = 1
		END IF

	END DO

	DEALLOCATE(PD_now)

	SAng = x 

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE max_pd(M,NUM_AN)
	USE tables, ONLY : PD_now

	IMPLICIT NONE

	REAL(KIND=8)			:: M
	INTEGER						:: NUM_AN

	INTEGER						:: i

	M = 0.0

	DO i=1,NUM_AN
		IF (PD_now(i) > M) THEN
			M = PD_now(i)
		END IF
	END DO

	M = M + 1

END SUBROUTINE

