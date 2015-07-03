
SUBROUTINE TCS_HeH(E, TCS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Uses fitted TCS fucntions to obtain
! the total elastic cross section
! for the He+H collision at energy E
!
! Input  E   [eV]
! Output TCS [a0^2]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	
	IMPLICIT NONE

	! Inputs
	REAL(KIND=8)		:: E

	! Outputs
	REAL(KIND=8)		:: TCS

	! Internal
	REAL(KIND=8)		:: L1
	REAL(KIND=8)		:: L2
	REAL(KIND=8)		:: L3
	REAL(KIND=8)		:: L4
	REAL(KIND=8)		:: L5
	REAL(KIND=8)		:: L6
	REAL(KIND=8)		:: L7
	REAL(KIND=8)		:: L8
	REAL(KIND=8)		:: L9
	REAL(KIND=8)		:: L10

	REAL(KIND=8)		:: R1
	REAL(KIND=8)		:: R2
	REAL(KIND=8)		:: R3
	REAL(KIND=8)		:: R4
	REAL(KIND=8)		:: R5
	REAL(KIND=8)		:: R6
	REAL(KIND=8)		:: R7
	REAL(KIND=8)		:: R8
	REAL(KIND=8)		:: R9
	REAL(KIND=8)		:: R10

	REAL(KIND=8)		:: Sig 
	REAL(KIND=8)		:: Alh 

	L1 = 0.00D0
  L2 = 0.06D0
  L3 = 0.36D0
  L4 = 2.50D0
  L5 = 19.0D0
  L6 = 70.0D0
  L7 = 470.0D0
  L8 = 1000.0D0
  L9 = 3100.0D0
  L10 = 8400.0D0

  R1 = 0.06D0
  R2 = 0.36D0
  R3 = 2.5D0
  R4 = 19.0D0
  R5 = 70.0D0
  R6 = 470.0D0
  R7 = 1000.0D0
  R8 = 3100.0D0
  R9 = 8400.0D0
  R10 = 10000.0D0

  IF      (E .LT. R1) THEN
    sig = 61.896216
    alh = 0.087133
  ELSE IF ( (L2 .LE. E) .AND. (E .LT. R2) ) THEN
    sig = 70.309258
    alh = 0.074753
  ELSE IF ( (L3 .LE. E) .AND. (E .LT. R3) ) THEN
    sig = 55.363928
    alh = 0.105018
  ELSE IF ( (L4 .LE. E) .AND. (E .LT. R4) ) THEN
    sig = 44.212133
    alh = 0.143602
  ELSE IF ( (L5 .LE. E) .AND. (E .LT. R5) ) THEN
    sig = 38.393605
    alh = 0.179375
  ELSE IF ( (L6 .LE. E) .AND. (E .LT. R6) ) THEN
    sig = 33.227513
    alh = 0.233673
  ELSE IF ( (L7 .LE. E) .AND. (E .LT. R7) ) THEN
    sig = 30.921883
    alh = 0.329190
  ELSE IF ( (L8 .LE. E) .AND. (E .LT. R8) ) THEN
    sig = 30.978307
    alh = 0.447127
  ELSE IF ( (L9 .LE. E) .AND. (E .LT. R9) ) THEN
    sig = 34.323398
    alh = 0.538563
  ELSE IF ( E .GT. L10 ) THEN
    sig = 38.446181
    alh = 0.592398
 	ELSE 
		WRITE(*,*) 'Rank: ', myid
    WRITE(*,*) "Your Energy ", E, " is out of range for this function He+H!!"
    sig = 0.0D0
!    sig = 38.446181
    alh = 0.592398
		TCS = sig*(1.0D3/E)**alh
		E   = 1.0D0		
		WRITE(*,*) 'Here is the TCS: ', TCS
 	END IF 

  TCS = sig*(1000.0D0/E)**alh
	
END SUBROUTINE

!##############################################################################
!##############################################################################
!##############################################################################

SUBROUTINE TCS_HeO(E, TCS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Uses fitted TCS fucntions to obtain
! the total elastic cross section
! for the He+O collision at energy E
!
! Input  E   [eV]
! Output TCS [a0^2]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	
	IMPLICIT NONE

	! Inputs
	REAL(KIND=8)		:: E

	! Outputs
	REAL(KIND=8)		:: TCS

	! Internal
	REAL(KIND=8)		:: L1
	REAL(KIND=8)		:: L2
	REAL(KIND=8)		:: L3
	REAL(KIND=8)		:: L4
	REAL(KIND=8)		:: L5
	REAL(KIND=8)		:: L6
	REAL(KIND=8)		:: L7
	REAL(KIND=8)		:: L8
	REAL(KIND=8)		:: L9
	REAL(KIND=8)		:: L10

	REAL(KIND=8)		:: R1
	REAL(KIND=8)		:: R2
	REAL(KIND=8)		:: R3
	REAL(KIND=8)		:: R4
	REAL(KIND=8)		:: R5
	REAL(KIND=8)		:: R6
	REAL(KIND=8)		:: R7
	REAL(KIND=8)		:: R8
	REAL(KIND=8)		:: R9
	REAL(KIND=8)		:: R10

	REAL(KIND=8)		:: Sig 
	REAL(KIND=8)		:: Alh 


	L1 = 0.01D0
  L2 = 0.04D0
  L3 = 0.15D0
  L4 = 1.50D0
  L5 = 15.0D0
  L6 = 40.0D0
  L7 = 300.0D0
  L8 = 1050.0D0
  L9 = 3000.0D0

  R1 = 0.04D0
  R2 = 0.15D0
  R3 = 1.5D0
  R4 = 15.0D0
  R5 = 40.0D0
  R6 = 300.0D0
  R7 = 1050.0D0
  R8 = 3000.0D0
  R9 = 10000.0D0

  IF      (E .LT. R1) THEN
    sig = 7.441D0 
    alh = 0.325D0
  ELSE IF ( (L2 .LE. E) .AND. (E .LT. R2) ) THEN
    sig = 33.009D0
    alh = 0.178D0
  ELSE IF ( (L3 .LE. E) .AND. (E .LT. R3) ) THEN
    sig = 77.212D0
    alh = 0.082D0
  ELSE IF ( (L4 .LE. E) .AND. (E .LT. R4) ) THEN
    sig = 78.369D0
    alh = 0.080D0
  ELSE IF ( (L5 .LE. E) .AND. (E .LT. R5) ) THEN
    sig = 73.408D0
    alh = 0.096D0
  ELSE IF ( (L6 .LE. E) .AND. (E .LT. R6) ) THEN
    sig = 68.502D0
    alh = 0.117D0
  ELSE IF ( (L7 .LE. E) .AND. (E .LT. R7) ) THEN
    sig = 67.371D0
    alh = 0.131D0
  ELSE IF ( (L8 .LE. E) .AND. (E .LT. R8) ) THEN
    sig = 67.588D0
    alh = 0.153D0
	ELSE IF ( E .GE. L9 ) THEN
    sig = 70.206D0
    alh = 0.188D0
 	ELSE 
		WRITE(*,*) 'Rank: ', myid
    WRITE(*,*) "Your Energy ", E, " is out of range for this function He+O!!"
    sig = 0.0D0
!    sig = 70.206D0
    alh = 0.188D0
		TCS = sig*(1.0D3/E)**alh
		E   = 1.0D0		
		WRITE(*,*) 'Here is the TCS: ', TCS
 	END IF 

  TCS = sig*(1000.0D0/E)**alh

END SUBROUTINE

!##############################################################################
!##############################################################################
!##############################################################################

SUBROUTINE TCS_HeHe(E, TCS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Uses fitted TCS fucntions to obtain
! the total elastic cross section
! for the He+He collision at energy E
!
! Input  E   [eV]
! Output TCS [a0^2]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	
	IMPLICIT NONE

	! Inputs
	REAL(KIND=8)		:: E

	! Outputs
	REAL(KIND=8)		:: TCS

	! Internal
	REAL(KIND=8)		:: L1
	REAL(KIND=8)		:: L2
	REAL(KIND=8)		:: L3
	REAL(KIND=8)		:: L4
	REAL(KIND=8)		:: L5
	REAL(KIND=8)		:: L6
	REAL(KIND=8)		:: L7
	REAL(KIND=8)		:: L8
	REAL(KIND=8)		:: L9
	REAL(KIND=8)		:: L10

	REAL(KIND=8)		:: R1
	REAL(KIND=8)		:: R2
	REAL(KIND=8)		:: R3
	REAL(KIND=8)		:: R4
	REAL(KIND=8)		:: R5
	REAL(KIND=8)		:: R6
	REAL(KIND=8)		:: R7
	REAL(KIND=8)		:: R8
	REAL(KIND=8)		:: R9
	REAL(KIND=8)		:: R10

	REAL(KIND=8)		:: Sig 
	REAL(KIND=8)		:: Alh 


	L1 = 0.01D0
  L2 = 0.06D0
  L3 = 0.35D0
  L4 = 2.50D0
  L5 = 15.0D0
  L6 = 70.0D0
  L7 = 500.0D0
  L8 = 1000.0D0

  R1 = 0.06D0
  R2 = 0.35D0
  R3 = 2.5D0
  R4 = 15.0D0
  R5 = 70.0D0
  R6 = 500.0D0
  R7 = 1000.0D0
  R8 = 10000.0D0

  IF      (E .LT. R1) THEN
    sig = 60.9D0
    alh = 0.067D0
  ELSE IF ( (L2 .LE. E) .AND. (E .LT. R2) ) THEN
    sig = 67.8D0
    alh = 0.056D0
  ELSE IF ( (L3 .LE. E) .AND. (E .LT. R3) ) THEN
    sig = 57.971D0
    alh = 0.075D0
  ELSE IF ( (L4 .LE. E) .AND. (E .LT. R4) ) THEN
    sig = 49.208D0
    alh = 0.103D0
  ELSE IF ( (L5 .LE. E) .AND. (E .LT. R5) ) THEN
    sig = 45.143D0
    alh = 0.123D0
  ELSE IF ( (L6 .LE. E) .AND. (E .LT. R6) ) THEN
    sig = 41.931D0
    alh = 0.150D0
  ELSE IF ( (L7 .LE. E) .AND. (E .LT. R7) ) THEN
    sig = 41.178D0
    alh = 0.177D0
  ELSE IF ( (L8 .LE. E) .AND. (E .LT. R8) ) THEN
    sig = 41.175D0
    alh = 0.210D0
 	ELSE 
		WRITE(*,*) 'Rank: ', myid
    WRITE(*,*) "Your Energy ", E, " is out of range for this function He+He!!"
		sig = 0.0D0
!    sig = 41.175D0
    alh = 0.210D0
		TCS = sig*(1.0D3/E)**alh
		E   = 1.0D0		
		WRITE(*,*) 'Here is the TCS: ', TCS
 	END IF 

  TCS = sig*(1000.0D0/E)**alh

END SUBROUTINE


