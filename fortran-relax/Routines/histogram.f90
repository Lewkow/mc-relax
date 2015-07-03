
SUBROUTINE histogram( A, L, x_start, x_end, N, hist_x, hist_y )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Histogram takes an array A, of length L, 
! starting and ending x values, x_start, x_end, 
! and the number of bins, N, and outputs
! 2 arrays, hist_x and hist_y which contain the histogram
! data from A. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!! Inputs
	INTEGER										:: L
	INTEGER										:: N 
	REAL(KIND=8)							:: x_start	
	REAL(KIND=8)							:: x_end
	REAL(KIND=8),DIMENSION(L)	:: A

	!! Outputs
	REAL(KIND=8),DIMENSION(N)	:: hist_x
	REAL(KIND=8),DIMENSION(N)	:: hist_y

	!! Internal
	INTEGER										:: i
	INTEGER										:: j
	INTEGER										:: tot 
	REAL(KIND=8)							:: dx
	REAL(KIND=8)							:: left 
	REAL(KIND=8)							:: right 

	!! Bin Width	
	dx = (x_end - x_start)/DBLE(N)

	!! Loop over all bins	
	DO i=1,N
		!! Set up hist_x, evenly spaced bin array with values of the bin centers
		hist_x(i) = x_start + (i-1)*dx + dx/2.0D0	

		!! bin left value
		left      = hist_x(i) - dx/2.0D0
		!! bin right value
		right     = hist_x(i) + dx/2.0D0
		!! initialize counter for bin to 0
		tot       = 0
		!! Loop over entire array A
		DO j=1,L
			!! Increment counter if value is within bin left and right
			IF ( A(j) .GT. left .AND. A(j) .LT. right ) tot = tot+1
		END DO
		!! Set hist_y to bin counter
		hist_y(i) = DBLE(tot)	
	END DO

END SUBROUTINE histogram

