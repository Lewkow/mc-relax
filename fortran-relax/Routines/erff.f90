!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates error function using numerical recipe
! special functions from chapter 6: 
! Incomplete Gamma Function, Error Function, 
! Chi-Square Probability Function, Cumulative
! Poisson Function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE errf(x,ret)
	!!
	!! Uses gammp
	!!

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: x

	!! Outputs
	REAL(KIND=8)	:: ret

	!! Internal 
	REAL(KIND=8)	:: gammp

	IF (x .LT. 0.0D0) THEN
		ret = -gammp(0.5d0,x**2)
	ELSE
		ret =  gammp(0.5d0,x**2)
	END IF

END SUBROUTINE errf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

FUNCTION gammp(a,x)
	!!
	!! Uses gser, gcf
	!!

	IMPLICIT NONE

	REAL(KIND=8)		:: a, gammp, x
	REAL(KIND=8)		:: gammcf, gamser, gln
	IF (x.LT.0. .OR. a .LE. 0.) WRITE(*,*) 'bad arguments in gammp'
	IF (x.LT.a+1.)THEN
		CALL gser(gamser,a,x,gln)
		gammp = gamser
	ELSE
		CALL gcf(gammcf,a,x,gln)
		gammp=1.-gammcf
	END IF
	RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

SUBROUTINE gser(gamser,a,x,gln)
	!!
	!! Uses gammln
	!!

	IMPLICIT NONE

	INTEGER 			:: ITMAX
	REAL(KIND=8)	:: a, gamser, gln, x, EPS
	PARAMETER (ITMAX=100, EPS=3.e-7)
	INTEGER				:: n
	REAL(KIND=8)	:: ap, del, sum, gammln
	
	gln = gammln(a)
	IF (x.LE.0.)THEN
		IF (x.LT.0.) WRITE(*,*) 'x < 0 in gser'
		gamser=0.
		RETURN
	END IF
	ap = a
	sum = 1./a
	del = sum
	DO n=1,ITMAX
		ap = ap+1.
		del = del*x/ap
		sum = sum+del
		if(abs(del).lt.abs(sum)*EPS)goto 1
	END DO
	WRITE(*,*) 'a too large, ITMAX too small in gser'
1 gamser=sum*exp(-x+a*log(x)-gln)
	RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

SUBROUTINE gcf(gammcf,a,x,gln)
	!!
	!! Uses gammln
	!!

	IMPLICIT NONE

	INTEGER 			:: ITMAX
	REAL(KIND=8)	:: a, gammcf, gln, x, EPS, FPMIN
	PARAMETER (ITMAX=100, EPS=3.e-7, FPMIN=1.e-30)
	INTEGER 			:: i
	REAL(KIND=8)	:: an, b, c, d, del, h, gammln
	
	gln = gammln(a)
	b = x+1.-a
	c = 1./FPMIN
	d = 1./b
	h = d
	do i=1,ITMAX
		an = -i*(i-a)
		b = b+2.
		d = an*d+b
		if(abs(d).lt.FPMIN) d=FPMIN
		c = b+an/c
		if(abs(c).lt.FPMIN) c=FPMIN
		d = 1./d
		del = d*c
		h = h*del
		if(abs(del-1.).lt.EPS)goto 1
	end do
	WRITE(*,*) 'a too large, ITMAX too small in gcf'
1 gammcf = exp(-x+a*log(x)-gln)*h
	RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

FUNCTION gammln(xx)
  IMPLICIT NONE

	REAL(KIND=8)		:: gammln, xx
	INTEGER					:: j, i
	REAL(KIND=8)		:: ser, stp, tmp, x, y

  DATA stp /2.5066282746310005d0/
  REAL(KIND=8), DIMENSION(6) :: coef = (/76.18009172947146d0,&
    -86.50532032941677d0,24.01409824083091d0,&
    -1.231739572450155d0,0.1208650973866179d-2,&
    -0.5395239384953d-5/)
	SAVE coef, stp

  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do i=1,size(coef)
    y=y+1.0d0
    ser=ser+coef(i)/y
  end do
  gammln=tmp+log(stp*ser/x)
	RETURN
END FUNCTION gammln

