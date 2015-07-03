
SUBROUTINE find_new_energy(Enow,SAng,rat,Enxt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Find new energy based on scattering angle
!	SAng is CM scattering angle
!	Enow and Enxt are both LF energies
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!--------------------
	!	INPUTS	
	!--------------------

	REAL(KIND=8)		:: rat			! rat = MP/MT ratio of masses	
	REAL(KIND=8)		:: Enow
	REAL(KIND=8)		:: SAng 

	!--------------------
	!	OUTPUTS	
	!--------------------

	REAL(KIND=8)		:: Enxt

	!--------------------
	!	INTERNAL	
	!--------------------

	REAL(KIND=8)		:: C 

	C    = (1.0D0 + 2.0D0*rat*COS(SAng) + rat*rat)/((1.0D0 + rat)**2.0)
	IF (C .GT. 1.0D0) THEN
!		WRITE(*,*) 'Ang: rat: C: ', SAng*180.0D0/3.14D0, rat, C
		C = 1.0D0
	END IF
	Enxt = Enow*C
	IF (Enxt .GE. Enow) Enxt = Enow

	IF (Enxt .GT. Enow) WRITE(*,*) C, SAng*180.0D0/3.14D0, Enow, Enxt, abs(Enow-Enxt)

!	WRITE(*,*) 'Ang: Enow: Enxt: ', SAng*180.0D0/3.14D0, Enow, Enxt

END SUBROUTINE

