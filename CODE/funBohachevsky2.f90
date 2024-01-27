!====================================
! Bohachevsky Function 2
! n = 2
! -100 <= x <= 100
! fstar = 0
!====================================
subroutine funct(n,xx,f)
	implicit none
	integer		:: n
	real*8		:: xx(n),x(n), f
	real*8, PARAMETER :: PI = 3.1415923
    x = Xx + 5.d0*dsqrt(2.d0)
	f = x(1)**2.d0 + 2.d0*x(2)**2.d0 
	f = f - 0.3d0*dcos(3.d0*PI*x(1))*dcos(4.d0*PI*x(2)) + 0.3d0

	return
end subroutine funct

subroutine grad(n,y,g)
	implicit none
	integer		:: n
	real*8		:: x(n), y(n), g(n)
	real*8, PARAMETER :: PI = 3.1415923

	x = y + 5.d0*sqrt(2.d0)

	g(1) = 2.d0*x(1) + 0.9d0*PI*dsin(3.d0*PI*x(1))*dcos(4.d0*PI*x(2))
	g(2) = 4.d0*x(2) + 0.3d0*4.d0*PI*dcos(3.d0*PI*x(1))*dsin(4.d0*PI*x(2))

	return
end subroutine grad
