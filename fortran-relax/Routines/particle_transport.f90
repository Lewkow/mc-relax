
SUBROUTINE particle_dir_cos_transport( u_now, x_now, theta, phi, dx, u_nxt, x_nxt ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates transport parameters for particle collisions. 
! Given a nascent particle with direction cosines u(i) at location x(i)
! as well as a random free path dx, the particle may be transported
! to the next collision point x(i+1) at which point, given a
! scattering angle theta [0-pi] and a random angle phi [0-2pi] 
! the new direction cosine u(i+1) may be found.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8),DIMENSION(3)	:: u_now			! Unit velocity vector currently
	REAL(KIND=8),DIMENSION(3)	:: x_now			! Position vector currently
	REAL(KIND=8)							:: theta			! Current theta scattering angle
	REAL(KIND=8)							:: phi				! Current phi scattering angle
	REAL(KIND=8)							:: dx 				! Random free path

	!! Outputs
	REAL(KIND=8),DIMENSION(3)	:: u_nxt			! Next unit velocity vector
	REAL(KIND=8),DIMENSION(3)	:: x_nxt			! Next position vector

	!! Internal
	INTEGER										:: i					! Index counter
	REAL(KIND=8)							:: C1


	!!!! Transport particle from x_now to x_nxt
	!!!! Before update directional cosines
	DO i=1,3
		x_nxt(i) = x_now(i) + u_now(i)*dx
	END DO

!	WRITE(*,*) 'dL: ', dx
!	WRITE(*,*) 'x0: x1: ux: ', x_now(1), x_nxt(1), u_now(1)
!	WRITE(*,*) 'y0: y1: uy: ', x_now(2), x_nxt(2), u_now(2)
!	WRITE(*,*) 'z0: z1: uz: ', x_now(3), x_nxt(3), u_now(3)

	IF (u_now(3) .EQ. 1.0D0) THEN
		u_nxt(1) = SIN(theta)*COS(phi)
		u_nxt(2) = SIN(theta)*SIN(phi)
		u_nxt(3) = COS(theta)
	ELSE IF (u_now(3) .EQ. -1.0D0) THEN
		u_nxt(1) = SIN(theta)*COS(phi)
		u_nxt(2) = -SIN(theta)*SIN(phi)
		u_nxt(3) = -COS(theta)
	ELSE
		C1				= SQRT(1.0D0-u_now(3)**2)
		u_nxt(1) 	= (SIN(theta)*(u_now(1)*u_now(3)*COS(phi)-u_now(2)*SIN(phi)))/C1 + u_now(1)*COS(theta)
		u_nxt(2) 	= (SIN(theta)*(u_now(2)*u_now(3)*COS(phi)+u_now(1)*SIN(phi)))/C1 + u_now(2)*COS(theta)
		u_nxt(3) 	= -C1*SIN(theta)*COS(phi)+u_now(3)*COS(theta)
	END IF	

END SUBROUTINE particle_dir_cos_transport

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE particle_transport( u_now, x_now, theta_prev, phi_prev, theta_now, phi_now, dx, u_nxt, x_nxt ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates transport parameters for particle collisions. 
! Given a nascent particle with unit velocity u(i) at location x(i)
! as well as a random free path dx, the particle may be transported
! to the next collision point x(i+1) at which point, given a
! scattering angle theta [0-pi] and a random angle phi [0-2pi] 
! the new unit velocity u(i+1) may be found.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8),DIMENSION(3)	:: u_now			! Unit velocity vector currently
	REAL(KIND=8),DIMENSION(3)	:: x_now			! Position vector currently
	REAL(KIND=8)							:: theta_prev	! Previous theta scattering angle
	REAL(KIND=8)							:: phi_prev		! Previous phi scattering angle
	REAL(KIND=8)							:: theta_now	! Current theta scattering angle
	REAL(KIND=8)							:: phi_now		! Current phi scattering angle
	REAL(KIND=8)							:: dx 				! Random free path

	!! Outputs
	REAL(KIND=8),DIMENSION(3)	:: u_nxt			! Next unit velocity vector
	REAL(KIND=8),DIMENSION(3)	:: x_nxt			! Next position vector

	!! Internal
	INTEGER										:: i					! Index counter

	!!!! Transport particle from x_now to x_nxt
	DO i=1,3
		x_nxt(i) = x_now(i) + u_now(i)*dx
	END DO

!	WRITE(*,'(A,6ES10.2)') 'zn: zx: tn: tx: pn: px: ', x_now(1), x_nxt(1), theta_prev, theta_now, phi_prev, phi_now

	!!!! Rotate velocity vector for new unit velocity
	u_nxt(1) = u_now(1)*COS(theta_now) + SIN(theta_now) * (u_now(3)*COS(phi_now)*COS(theta_prev) - SIN(phi_now)*SIN(phi_prev))
	u_nxt(2) = u_now(2)*COS(theta_now) + SIN(theta_now) * (u_now(3)*COS(phi_now)*SIN(theta_prev) + SIN(phi_now)*COS(phi_prev))
	u_nxt(3) = u_now(3)*COS(theta_now) - SIN(theta_now)*SIN(theta_prev)*COS(phi_now)


END SUBROUTINE 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE velocity_unit_vector( u_now, theta_prev, phi_prev, theta_now, phi_now, u_nxt ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates new unit velocity vector in lab frame given 
! the previous unit vector, scattering angles and the 
! current scattering angles. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8),DIMENSION(3)	:: u_now			! Unit velocity vector currently
	REAL(KIND=8)							:: theta_prev	! Previous theta scattering angle
	REAL(KIND=8)							:: phi_prev		! Previous phi scattering angle
	REAL(KIND=8)							:: theta_now	! Current theta scattering angle
	REAL(KIND=8)							:: phi_now		! Current phi scattering angle

	!! Outputs
	REAL(KIND=8),DIMENSION(3)	:: u_nxt			! Next unit velocity vector

	!! Internal
	INTEGER										:: i					! Index counter

	!!!! Rotate velocity vector for new unit velocity
	u_nxt(1) = u_now(1)*COS(theta_now) + SIN(theta_now) * (u_now(3)*COS(phi_now)*COS(theta_prev) - SIN(phi_now)*SIN(phi_prev))
	u_nxt(2) = u_now(2)*COS(theta_now) + SIN(theta_now) * (u_now(3)*COS(phi_now)*SIN(theta_prev) + SIN(phi_now)*COS(phi_prev))
	u_nxt(3) = u_now(3)*COS(theta_now) - SIN(theta_now)*SIN(theta_prev)*COS(phi_now)


END SUBROUTINE 

