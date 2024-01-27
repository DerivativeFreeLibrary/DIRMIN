!=========================================
!
! FUNZIONE DI GRIEWANK (intervallo ??)
!
!=========================================

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F, xaux(n)

      INTEGER          :: I

	  xaux = X + 5.d0*dsqrt(2.d0)

      F = 1.D0

      DO I=1,N
         F = F*DCOS(xaux(I)/DSQRT(DBLE(I)))
      END DO

      F = -F+1.D0

      DO I=1,N
         F = F + xaux(I)*xaux(I)/4.0D+3
      END DO

      RETURN

END


subroutine grad(n,x,g)

implicit none

integer				::n, i, j
double precision	::x(n), g(n), xaux(n), prod(n)

xaux = X + 5.d0*dsqrt(2.d0)

prod = 1.0d0

do i=1,n
	do j= 1,n
		if (i.ne.j) then
			prod(i) = prod(i)*dcos(xaux(j)/dsqrt(dble(j))) 
		endif
	enddo
end do


do i =1,n

	g(i) = (2.0d0*xaux(i))/4.0d3+(1/dsqrt(dble(i)))*dsin(xaux(i)/dsqrt(dble(i)))*prod(i)

enddo


return

end subroutine
