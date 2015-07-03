SUBROUTINE vec_histogram( A, L, XA, XL, hist_y )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Histogram takes a data array A, of length L (Y-data)
! and x-vector XA and length XL (X-data)
! and outputs hist_y which contain the histogram
! data from A and has length XL. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!! Inputs
	INTEGER											:: L
	INTEGER											:: XL
	REAL(KIND=8),DIMENSION(L)		:: A
	REAL(KIND=8),DIMENSION(XL)	:: XA

	!! Outputs
	REAL(KIND=8),DIMENSION(XL)	:: hist_y

!! Internal
	INTEGER										:: i
	INTEGER										:: j
	INTEGER										:: tot, tot_y 
	REAL(KIND=8)							:: dx
	REAL(KIND=8)							:: left 
	REAL(KIND=8)							:: right 

	tot_y = 0.0D0

	!! Bin Width	
	dx = XA(2) - XA(1)

	!! Loop over all bins	
	DO i=1,XL
		!! bin left value
		IF (i .EQ. 1) THEN
			left = -1.0D20
		ELSE
			left = XA(i) - dx/2.0D0
		END IF
		!! bin right value
		IF (i .EQ. XL) THEN
			right = 1.0D20
		ELSE
			right = XA(i) + dx/2.0D0
		END IF
		!! initialize counter for bin to 0
		tot = 0
		!! Loop over entire array A
		DO j=1,L
			!! Increment counter if value is within bin left and right
			IF ( A(j) .GT. left .AND. A(j) .LT. right ) tot = tot+1
		END DO
		!! Set hist_y to bin counter
		hist_y(i) = DBLE(tot)	
		tot_y			= tot_y + DBLE(tot)
	END DO

	DO i=1,XL
		hist_y(i) = hist_y(i)/tot_y
	END DO

END SUBROUTINE vec_histogram



