
!! Fill all phase spaces which will
!! be written to files for data
!! analysis
!! ALL CLICKS AVERAGE

SUBROUTINE fill_all_phase_spaces( MC )
	USE click_data
	USE planet

	IMPLICIT NONE

	!! Input
	INTEGER					:: MC		! current MC particle

	!! Internal
	INTEGER					:: i, j, iX, iY
	REAL(KIND=8)		:: ddT, ddE, ddN, dEL, ddXY, ddU, ddZ, X, Y, dX, dY, X1, X2, Y1, Y2

	ddT  = (X_CLICK(2,1) - X_CLICK(1,1))/2.0D0
	ddE  = (X_CLICK(2,2) - X_CLICK(1,2))/2.0D0
	ddXY = (X_CLICK(2,3) - X_CLICK(1,3))/2.0D0
	ddZ  = (X_CLICK(2,5) - X_CLICK(1,5))/2.0D0
	ddU  = (X_CLICK(2,6) - X_CLICK(1,6))/2.0D0
	dEL  = (X_CLICK(2,12) - X_CLICK(1,12))/2.0D0
	ddN  = (X_CLICK(2,11) - X_CLICK(1,11))/2.0D0

	IF (WRITE_ALL_E_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <E> vs T
		! Energy vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 2 ! energy
		iX = 1 ! time
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddE
		dX = ddT
		DO i=1,Num_Hist
			X1 = X_CLICK(i,iX)-dX
			X2 = X_CLICK(i,iX)+dX
			DO j=1,Num_Hist
				Y1 = X_CLICK(j,iY)-dY
				Y2 = X_CLICK(j,iY)+dY
				IF ((X.GE.X1).AND.(X.LE.X2)) THEN
					IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
						A_E_T(i,j) = A_E_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
					END IF ! Y
				END IF ! X 
			END DO ! j 
		END DO ! i
	END IF

	IF (WRITE_ALL_H_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <H> vs T
		! Height vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 5 ! height
		iX = 1 ! time
 		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddZ
		dX = ddT
		DO i=1,Num_Hist
 	  	X1 = X_CLICK(i,iX)-dX
 	   	X2 = X_CLICK(i,iX)+dX
    	DO j=1,Num_Hist
     		Y1 = X_CLICK(j,iY)-dY
      	Y2 = X_CLICK(j,iY)+dY
      	IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        	IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          	A_H_T(i,j) = A_H_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
        	END IF ! Y
      	END IF ! X 
    	END DO ! j 
  	END DO ! i
	END IF

	IF (WRITE_ALL_N_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <N> vs T
		! Number Collisions vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 11! number collisions
		iX = 1 ! time
 		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddN
		dX = ddT
 	 	DO i=1,Num_Hist
    	X1 = X_CLICK(i,iX)-dX
    	X2 = X_CLICK(i,iX)+dX
    	DO j=1,Num_Hist
     		Y1 = X_CLICK(j,iY)-dY
      	Y2 = X_CLICK(j,iY)+dY
      	IF ((X.GT.X1).AND.(X.LE.X2)) THEN
        	IF ((Y.GT.Y1).AND.(Y.LE.Y2)) THEN
          	A_N_T(i,j) = A_N_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
        	END IF ! Y
      	END IF ! X 
    	END DO ! j 
  	END DO ! i
	END IF

	IF (WRITE_ALL_Ux_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y   vs x
		! <Uz> vs T
		! Vert Vel vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 8 ! vert vel
		iX = 1 ! time
 		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddU
		dX = ddT
 		DO i=1,Num_Hist
    	X1 = X_CLICK(i,iX)-dX
    	X2 = X_CLICK(i,iX)+dX
    	DO j=1,Num_Hist
     		Y1 = X_CLICK(j,iY)-dY
      	Y2 = X_CLICK(j,iY)+dY
      	IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        	IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          	A_Ux_T(i,j) = A_Ux_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
        	END IF ! Y
      	END IF ! X 
    	END DO ! j 
  	END DO ! i
	END IF

	IF (WRITE_ALL_H_Ux .EQ. 1) THEN
	  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y   vs x
		! <Uz> vs H
 		! Vert Velocity vs Height
  	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iX = 5 ! height
		iY = 8 ! vert vel
  	X  = CLICK_1(MC,iX)
  	Y  = CLICK_1(MC,iY)
  	dX = ddZ
  	dY = ddU
  	DO i=1,Num_Hist
    	X1 = X_CLICK(i,iX)-dX
    	X2 = X_CLICK(i,iX)+dX
    	DO j=1,Num_Hist
     		Y1 = X_CLICK(j,iY)-dY
      	Y2 = X_CLICK(j,iY)+dY
      	IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        	IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          	A_H_Ux(i,j) = A_H_Ux(i,j) + 1.0D0/REAL(N_Part*Max_Click)
        	END IF ! Y
      	END IF ! X 
    	END DO ! j 
  	END DO ! i
	END IF

	IF (WRITE_ALL_H_E .EQ. 1) THEN
	  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <E> vs H
 		! Energy vs Height Prob Den
  	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 2 ! energy
		iX = 5 ! height
  	X  = CLICK_1(MC,iX)
  	Y  = CLICK_1(MC,iY)
  	dY = ddE
  	dX = ddZ
!		WRITE(*,'(A,2ES12.2)') 'Z: E: ', X, Y
  	DO i=1,(Num_Hist-1)
    	X1 = X_CLICK(i,iX)
    	X2 = X_CLICK(i+1,iX)
    	DO j=1,(Num_Hist-1)
      	Y1 = X_CLICK(j,iY)
      	Y2 = X_CLICK(j+1,iY)
      	IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        	IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          	A_H_E(i,j) = A_H_E(i,j) + 1.0D0/REAL(N_Part*Max_Click)
        	END IF ! Y
      	END IF ! X 
    	END DO ! j 
  	END DO ! i
	END IF

	IF (WRITE_ALL_H_N .EQ. 1) THEN
	  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <N> vs H
 		! Collisions Number vs Height Prob Den
	  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iX = 5 ! height
		iY = 11! collision number
  	X  = CLICK_1(MC,iX)
  	Y  = CLICK_1(MC,iY)
  	dX = ddZ
  	dY = ddN
  	DO i=1,Num_Hist
    	X1 = X_CLICK(i,iX)-dX
    	X2 = X_CLICK(i,iX)+dX
    	DO j=1,Num_Hist
      	Y1 = X_CLICK(j,iY)-dY
      	Y2 = X_CLICK(j,iY)+dY
      	IF ((X.GT.X1).AND.(X.LE.X2)) THEN
        	IF ((Y.GT.Y1).AND.(Y.LE.Y2)) THEN
          	A_H_N(i,j) = A_H_N(i,j) + 1.0D0/REAL(N_Part*Max_Click)
        	END IF ! Y
      	END IF ! X 
    	END DO ! j 
  	END DO ! i
	END IF

	IF (WRITE_ALL_H_dE .EQ. 1) THEN
	  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y   vs x
		! <dE> vs H
 		! Energy Loss vs Height Prob Den
  	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iX = 5  ! height
		iY = 12 ! dE for X ref
  	X  = CLICK_1(MC,iX)
  	Y  = CLICK_0(MC,2)-CLICK_1(MC,2)
  	dX = ddZ
  	dY = dEL
  	DO i=1,Num_Hist
    	X1 = X_CLICK(i,iX)-dX
    	X2 = X_CLICK(i,iX)+dX
    	DO j=1,Num_Hist
      	Y1 = X_CLICK(j,iY)-dY
      	Y2 = X_CLICK(j,iY)+dY
      	IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        	IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          	A_H_dE(i,j) = A_H_dE(i,j) + 1.0D0/REAL(N_Part*Max_Click)
        	END IF ! Y
      	END IF ! X 
    	END DO ! j 
  	END DO ! i
	END IF

END SUBROUTINE fill_all_phase_spaces

!#######################################
!#######################################
!#######################################

!! Fill all phase spaces which will
!! be written to files for data
!! analysis

SUBROUTINE fill_phase_spaces( ck, MC )

	USE click_data
	USE planet, ONLY : N_Part

	IMPLICIT NONE

	!! Input
	INTEGER					:: ck		! click
	INTEGER					:: MC		! current MC particle

	!! Internal
	INTEGER					:: i, j, iX, iY
	REAL(KIND=8)		:: ddT, ddE, ddN, dEL, ddXY, ddU, ddZ, X, Y, dX, dY, X1, X2, Y1, Y2

	ddT  = (X_CLICK(2,1) - X_CLICK(1,1))/2.0D0
	ddE  = (X_CLICK(2,2) - X_CLICK(1,2))/2.0D0
	ddXY = (X_CLICK(2,3) - X_CLICK(1,3))/2.0D0
	ddZ  = (X_CLICK(2,5) - X_CLICK(1,5))/2.0D0
	ddU  = (X_CLICK(2,6) - X_CLICK(1,6))/2.0D0
	dEL  = (X_CLICK(2,12) - X_CLICK(1,12))/2.0D0
	ddN  = (X_CLICK(2,11) - X_CLICK(1,11))/2.0D0

	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	! Energy vs Time Prob Den
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	iX = 2 ! energy
	iY = 1 ! time
	X  = CLICK_1(MC,iX)
	Y  = CLICK_1(MC,iY)
	dX = ddE
	dY = ddT
	DO i=1,Num_Hist
		X1 = X_CLICK(i,iX)-dX
		X2 = X_CLICK(i,iX)+dX
		DO j=1,Num_Hist
			Y1 = X_CLICK(j,iY)-dY
			Y2 = X_CLICK(j,iY)+dY
			IF ((X.GE.X1).AND.(X.LE.X2)) THEN
				IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
					C_E_T(ck,i,j) = C_E_T(ck,i,j) + 1.0D0/REAL(N_Part)
				END IF ! Y
			END IF ! X 
		END DO ! j 
	END DO ! i

	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	! Height vs Time Prob Den
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	iX = 5 ! height
	iY = 1 ! time
 	X  = CLICK_1(MC,iX)
	Y  = CLICK_1(MC,iY)
	dX = ddZ
	dY = ddT
  DO i=1,Num_Hist
    X1 = X_CLICK(i,iX)-dX
    X2 = X_CLICK(i,iX)+dX
    DO j=1,Num_Hist
      Y1 = X_CLICK(j,iY)-dY
      Y2 = X_CLICK(j,iY)+dY
      IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          C_H_T(ck,i,j) = C_H_T(ck,i,j) + 1.0D0/REAL(N_Part)
        END IF ! Y
      END IF ! X 
    END DO ! j 
  END DO ! i

	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	! Number Collisions vs Time Prob Den
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	iX = 11! number collisions
	iY = 1 ! time
 	X  = CLICK_1(MC,iX)
	Y  = CLICK_1(MC,iY)
	dX = ddN
	dY = ddT
  DO i=1,Num_Hist
    X1 = X_CLICK(i,iX)-dX
    X2 = X_CLICK(i,iX)+dX
    DO j=1,Num_Hist
      Y1 = X_CLICK(j,iY)-dY
      Y2 = X_CLICK(j,iY)+dY
      IF ((X.GT.X1).AND.(X.LE.X2)) THEN
        IF ((Y.GT.Y1).AND.(Y.LE.Y2)) THEN
          C_N_T(ck,i,j) = C_N_T(ck,i,j) + 1.0D0/REAL(N_Part)
        END IF ! Y
      END IF ! X 
    END DO ! j 
  END DO ! i

	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	! Vert Vel vs Time Prob Den
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	iX = 8 ! vert vel
	iY = 1 ! time
 	X  = CLICK_1(MC,iX)
	Y  = CLICK_1(MC,iY)
	dX = ddU
	dY = ddT
  DO i=1,Num_Hist
    X1 = X_CLICK(i,iX)-dX
    X2 = X_CLICK(i,iX)+dX
    DO j=1,Num_Hist
      Y1 = X_CLICK(j,iY)-dY
      Y2 = X_CLICK(j,iY)+dY
      IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          C_Ux_T(ck,i,j) = C_Ux_T(ck,i,j) + 1.0D0/REAL(N_Part)
        END IF ! Y
      END IF ! X 
    END DO ! j 
  END DO ! i

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! Height vs Vert Velocity
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	iX = 5 ! height
	iY = 8 ! vert vel
  X  = CLICK_1(MC,iX)
  Y  = CLICK_1(MC,iY)
  dX = ddZ
  dY = ddU
  DO i=1,Num_Hist
    X1 = X_CLICK(i,iX)-dX
    X2 = X_CLICK(i,iX)+dX
    DO j=1,Num_Hist
      Y1 = X_CLICK(j,iY)-dY
      Y2 = X_CLICK(j,iY)+dY
      IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          C_H_Ux(ck,i,j) = C_H_Ux(ck,i,j) + 1.0D0/REAL(N_Part)
        END IF ! Y
      END IF ! X 
    END DO ! j 
  END DO ! i

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! Height vs Energy Prob Den
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	iX = 5 ! height
	iY = 2 ! energy
  X  = CLICK_1(MC,iX)
  Y  = CLICK_1(MC,iY)
  dX = ddZ
  dY = ddE
  DO i=1,Num_Hist
    X1 = X_CLICK(i,iX)-dX
    X2 = X_CLICK(i,iX)+dX
    DO j=1,Num_Hist
      Y1 = X_CLICK(j,iY)-dY
      Y2 = X_CLICK(j,iY)+dY
      IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          C_H_E(ck,i,j) = C_H_E(ck,i,j) + 1.0D0/REAL(N_Part)
        END IF ! Y
      END IF ! X 
    END DO ! j 
  END DO ! i

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! Height vs Collisions Number Prob Den
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	iX = 5 ! height
	iY = 11! collision number
  X  = CLICK_1(MC,iX)
  Y  = CLICK_1(MC,iY)
  dX = ddZ
  dY = ddN
  DO i=1,Num_Hist
    X1 = X_CLICK(i,iX)-dX
    X2 = X_CLICK(i,iX)+dX
    DO j=1,Num_Hist
      Y1 = X_CLICK(j,iY)-dY
      Y2 = X_CLICK(j,iY)+dY
      IF ((X.GT.X1).AND.(X.LE.X2)) THEN
        IF ((Y.GT.Y1).AND.(Y.LE.Y2)) THEN
          C_H_N(ck,i,j) = C_H_N(ck,i,j) + 1.0D0/REAL(N_Part)
        END IF ! Y
      END IF ! X 
    END DO ! j 
  END DO ! i

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! Height vs Energy Loss Prob Den
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	iX = 5  ! height
	iY = 12 ! dE for X ref
  X  = CLICK_1(MC,iX)
  Y  = CLICK_0(MC,2)-CLICK_1(MC,2)
  dX = ddZ
  dY = dEL
  DO i=1,Num_Hist
    X1 = X_CLICK(i,iX)-dX
    X2 = X_CLICK(i,iX)+dX
    DO j=1,Num_Hist
      Y1 = X_CLICK(j,iY)-dY
      Y2 = X_CLICK(j,iY)+dY
      IF ((X.GE.X1).AND.(X.LE.X2)) THEN
        IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
          C_H_dE(ck,i,j) = C_H_dE(ck,i,j) + 1.0D0/REAL(N_Part)
        END IF ! Y
      END IF ! X 
    END DO ! j 
  END DO ! i



END SUBROUTINE fill_phase_spaces

