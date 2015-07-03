
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file houses the different MC experiments to be done
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE average_angle(Energy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This expt takes an ensemble of P particles
! and finds the looks at the distribution
! of scattering angles at a given energy. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE inputs

	IMPLICIT NONE

	REAL(KIND=8)				:: Energy
	REAL(KIND=8)				:: SAng 

	INTEGER							:: P

	CHARACTER(LEN=21)		:: Fname
	
	WRITE(UNIT=Fname, FMT="(A, ES7.1, A)") "AngDist_", Energy, "eV.dat"

	OPEN(UNIT=33, FILE=Fname, ACCESS="APPEND")
	
	WRITE(*,*) "Starting average angle for Energy ", Energy, " eV"	
	
	DO P=1,N_Part
		CALL find_rand_angle(Energy,SAng)			
		WRITE(33,*) SAng
	END DO

	CLOSE(33)

END SUBROUTINE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE ensemble_every_N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This expt calculates the averag energy
! loss for an ensemble of P particles, all
! starting at an energy E0. The average
! energy after one collision is then the 
! next E0 for the ensemble. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE inputs
	USE rand_seed

	IMPLICIT NONE

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Enxt_Ave
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Enxt_Array

	REAL(KIND=8)														:: Enow
	REAL(KIND=8)														:: Enxt
	REAL(KIND=8)														:: SAng 
	REAL(KIND=8)														:: Enow_Sum 
	REAL(KIND=8)														:: Var_Sum 
	REAL(KIND=8)														:: t_stop 
	REAL(KIND=8)														:: t_start 
	REAL(KIND=8)														:: Elap 
	REAL(KIND=8)														:: C 

	INTEGER																	:: N
	INTEGER																	:: P
	INTEGER																	:: N_N

	OPEN(UNIT=666, FILE="ENSEMB_PER_N.dat", ACCESS="APPEND")

	N_N  = 100
	Enow = EI

	ALLOCATE(Enxt_Ave(N_N))
	ALLOCATE(Enxt_Array(N_Part))

	DO N=1,N_N

		Enow_Sum = 0.0		
		CALL CPU_TIME(t_start)

		DO P=1,N_Part
			
			CALL find_rand_angle(Enow,SAng)
			CALL find_new_energy(Enow,SAng,Enxt)
			
			Enxt_Array(P) = Enxt	
			Enow_Sum 			= Enow_Sum + Enxt	

		END DO ! P

		CALL CPU_TIME(t_stop)
		Elap        = t_stop - t_start
		Enxt_Ave(N) = Enow_Sum/DBLE(N_Part)
		Enow        = Enxt_Ave(N)
		Var_Sum     = 0.0

		DO P=1,N_Part
			C 			= (Enxt_Array(P) - Enxt_Ave(N))**2.0
			Var_Sum = Var_Sum + C	
		END DO

		Var_Sum = Var_Sum/DBLE(N_Part-1)

		WRITE(*,*) N, "Collision Complete in ", Elap, " sec. Average Energy ", Enow, " Var ", Var_Sum	
		WRITE(666,*) N, Enow, Var_Sum

	END DO ! N

	CLOSE(666)

	DEALLOCATE(Enxt_Ave)	
	DEALLOCATE(Enxt_Array)	

END SUBROUTINE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE ensemble_prop_Nset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This expt starts an ensemble of particles at an 
! energy EI and then watches as the propogate for a 
! set number of collisions. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE inputs
	USE rand_seed

	IMPLICIT NONE

	REAL(KIND=8)						:: N_TOT	
	REAL(KIND=8)						:: Enow	
	REAL(KIND=8)						:: T_Start
	REAL(KIND=8)						:: F_START
	REAL(KIND=8)						:: SAng
	REAL(KIND=8)						:: Enxt
	REAL(KIND=8)						:: T_Stop
	REAL(KIND=8)						:: F_Stop
	REAL(KIND=8)						:: T_Elap 
	REAL(KIND=8)						:: F_Elap 
	REAL(KIND=8)						:: Ave_N 

	INTEGER									:: z
	INTEGER									:: P
	INTEGER									:: N 
	
 	OPEN(UNIT=444, FILE="Ensemble_N_Coll_Relax.dat", ACCESS="APPEND")

  DO z=1,NE
    EF    = EI - z*DE
    N_TOT = 0

    CALL CPU_TIME(F_START)

    DO P=1,N_Part

      Enow = EI
      CALL cpu_time(T_Start)

      DO N=1,100
        !----------------------
        ! Find Random Angle 
        !----------------------
        CALL find_rand_angle(Enow,SAng)

        !----------------------
        ! Find New Energy 
        !----------------------
        CALL find_new_energy(Enow,SAng,Enxt)

        Enow   = Enxt

        WRITE(444,*) N, Enow

      END DO !N !! Enow

      CALL cpu_time(T_Stop)
      T_Elap     = T_Stop-T_Start

      !----------------------
      ! Write Part Relax Info 
      !----------------------
      CALL disp_relax_info(N-1,P,T_Elap)

    END DO ! P

    CALL CPU_TIME(F_STOP)
    F_ELAP = F_STOP-F_START

    !----------------------
    ! Calc Ave Num Coll Ave_N 
    !----------------------
    Ave_N = DBLE(N_TOT)/DBLE(N_Part)
    WRITE(*,*) "T ", F_ELAP, "E ", EF, "<N> ", Ave_N
  END DO ! z

  !----------------------
  ! Clean Up Memory 
  !----------------------

  CLOSE(444)

END SUBROUTINE
