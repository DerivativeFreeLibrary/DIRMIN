!====================================
! Booth function
! n = 2
! -10 <= x <= 10
! fstar = 0
!====================================
subroutine funct(n,x,f)
	implicit none
	integer		:: n
	real*8		:: x(n), f, A, B

	A = x(1) + 2.d0*x(2) - 7.d0
	B = 2.d0*x(1) + x(2) - 5.d0
	
	f = A**2.d0 + B**2.d0

	return
end subroutine funct

subroutine grad(n,x,g)
	implicit none
	integer		:: n
	real*8		:: x(n), g(n), A, B

	A = x(1) + 2.d0*x(2) - 7.d0
	B = 2.d0*x(1) + x(2) - 5.d0
	g(1) = 2.d0*A + 4.d0*B
	g(2) = 4.d0*A + 2.d0*B

	return
end subroutine grad
