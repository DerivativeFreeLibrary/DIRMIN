!====================================
! Bohachevsky Function 1
! n = 2
! -100 <= x <= 100
! fstar = 0
!====================================
subroutine funct(n,x,f)
	implicit none
	integer		:: n
	real*8		:: x(n), f,xx(n)
	real*8, PARAMETER :: PI = 3.1415923
    xx = X + 5.d0*dsqrt(2.d0)
	f = xx(1)**2.d0 + 2.d0*xx(2)**2.d0 - 0.3d0*dcos(3.d0*PI*xx(1))
	f = f - 0.4d0*dcos(4.d0*PI*xx(2)) + 0.7d0

	return
end subroutine funct

subroutine grad(n,x,g)
	implicit none
	integer		:: n
	real*8		:: x(n), g(n), xx(n)
	real*8, PARAMETER :: PI = 3.1415923

    xx = X + 5.d0*dsqrt(2.d0)

	g(1) = 2.d0*xx(1) + 0.9d0*PI*dsin(3.d0*PI*xx(1))
	g(2) = 4.d0*( xx(2) + 0.4d0*PI*dsin(4.d0*PI*xx(2)) )

	return
end subroutine grad
