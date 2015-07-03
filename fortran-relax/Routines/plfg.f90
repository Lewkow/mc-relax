
!*********************SUBROUTINE RINT*******************
!													   !
! Subroutine that initializes canonical rectangle of   !
! LFG(127,97,31) by employing the Galois binary shift  !
! register		                                       !
!													   !
! Input variable:                                      !
!													   !
! myid - integer process id used to seed Galois 	   !
!        register									   !
!													   !
!                                                      !
!-------------------------------------------------------
subroutine rint(myid)
!-------------------------------------------------------
	USE rand_seed

	IMPLICIT NONE
	
integer*4 bit, myid, galois
INTEGER		:: i
INTEGER		:: j
!
!
!
!-----------------------------------------------------
!PHASE I - INITIALIZE VARIABLES 
!-----------------------------------------------------
!Mask seed with process ID
seed=seed+myid
!-----------------------------------------------------
lag1=1
lag2=large_lag-small_lag+1
register=0
rmult = 1.0d0/modulus
!-----------------------------------------------------
!
!
!-----------------------------------------------------
!PHASE II - INITIALIZE CANONICAL RECTANGLE
!-----------------------------------------------------
!Run-up the Galois register
!-----------------------------------------------------
do i=1,1000
 bit=galois()
end do
!stop
!-----------------------------------------------------
!Set canonical bit in LFG register
!-----------------------------------------------------
register(cbit)=ibset(register(cbit),0)
!-----------------------------------------------------
!Initialize canonical rectangle
!-----------------------------------------------------
do i=1,30
 do j=1,(large_lag-1)                    
  bit = galois()
  if (bit==1) then
	register(j)=ibset(register(j),i)
  end if
 end do
end do
!-----------------------------------------------------
!Purge transient flat-spot in the ALFG register 
!-----------------------------------------------------
do i=1,1500
 register(lag1)=register(lag1)+register(lag2)
 register(lag1)=iand(register(lag1),modulus)
 lag1=lag1+1
 lag2=lag2+1
 if (lag2 > large_lag) then
	lag2=1
 elseif (lag1>large_lag) then
	lag1=1
 end if
end do
!-----------------------------------------------------

return
end

!**********************FUNCTION LFG****************************
!															  !
! Function that returns a random real number from             !
! LFG(127,97,31) output is real*8		                      !
!                                                             !
!--------------------------------------------------------------
function  lfg()
!--------------------------------------------------------------

	USE rand_seed

	IMPLICIT NONE

integer*4    bit, numtemp, galois
real*8       lfg
!
!
!
!
!-----------------------------------------------------
!PHASE I - GENERATE RANDOM INTEGER
!-----------------------------------------------------
register(lag1)=register(lag1)+register(lag2)
register(lag1)=iand(register(lag1),modulus)
numtemp=register(lag1)
!-----------------------------------------------------
!
!
!-----------------------------------------------------
!PHASE II - REPAIR DEFECTIVE BIT
!-----------------------------------------------------
!This phase replaces the least significant bit of the 
!random integer with output from the galois register
!fixes the defect in the lsb.
!-----------------------------------------------------
numtemp=ibclr(numtemp,0)
bit = galois()
if (bit==1) then
  numtemp=ibset(numtemp,0)
end if		
!-----------------------------------------------------
!
!
!-----------------------------------------------------
!PHASE III - OBTAIN RANDOM REAL
!-----------------------------------------------------
!Obtain random real number and advance the register
!forward
!-----------------------------------------------------
lfg=numtemp*rmult
lag1=lag1+1
lag2=lag2+1
if (lag2 > large_lag) then
  lag2=1
elseif (lag1>large_lag) then
  lag1=1
end if
!-----------------------------------------------------

return 
end
	

!*********************SUBROUTINE GALOIS****************************
!																  !
!  Subroutine to accomplish the Galois shift register.            !
!  This is an efficient and effective way to initialize 		  !
!  the ALFG which corresponds to a primitive trinomial.			  !
!  Output is integer*4.                                            ! 
!																  !
!------------------------------------------------------------------
function galois()
!-----------------------------------------------------
	USE rand_seed

	IMPLICIT NONE

integer*4  galois, mask

! Galois mask corresponding to primitive polynomial
! p(x) = 1 + x^25 + x^27 + x^29 + x^30 + x^31 + x^32
!mask = 2#10000000000000000000000001010111 !CAREFUL -- Not portable across non-MS-DEVELOPER studio!!
data mask /B'10000000000000000000000001010111'/    !USE this instead on non MS-DEV platforms


!OBTAIN THE LEAST SIGNIFICANT BIT
!------------------------------------
galois=0
if (btest(seed,0)) then		             
	galois=1
end if
!------------------------------------

!RIGHT SHIFT THE SEED AND XOR IT WITH (LSB*MASK)
!------------------------------------
seed=xor(ishft(seed,-1),galois*mask)	
!------------------------------------

return

end
!
!*******************SUBROUTINE GSEED*******************************
!																  !
!      Subroutine to initialize global seed						  !
!																  !
!******************************************************************
!
subroutine gseed

	USE rand_seed

	IMPLICIT NONE

integer itimes(8),iyr,imon,iday,ihr,imin,isec,imil
integer ndays,iyrs
integer jdays(12)
INTEGER :: i

!
!
data (jdays(i),i=1,12)/0,31,59,90,120,151,181,212,243,273,304,334/ 
!
!---------------------------------------------------------------------
!
!
!
! get initial seed
!
call date_and_time (values = itimes) 
iyr  = itimes(1)
imon = itimes(2)
iday = itimes(3)
ihr  = itimes(5)
imin = itimes(6)
isec = itimes(7)
imil = itimes(8)
!
! corrections for leap year
!
ndays = iday+jdays(imon)
if (mod(iyr,4) .eq. 0 .and. mod(iyr/100,4) .ne. 0 .and. imon .gt. 2) then
  ndays = ndays + 1
end if
iyrs = mod(iyr,32)
seed=iyrs+2**5*(imon+2**4*(iday+2**5*(ihr+2**5*(imin+2**6*(isec))))) 

return

end
