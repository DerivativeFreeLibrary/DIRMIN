!=========================================
!
! FUNZIONE DI TRECCANI
! 
!  X(1) .GE. -2.5  ;  X(2) .LE.  2.5
!
!=========================================

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F, xaux(n)

      DOUBLE PRECISION :: X1Q

      Xaux(1) = X(1) + 0.25d0*dsqrt(2.d0)
      Xaux(2) = X(2) + 0.25d0*dsqrt(2.d0)
	
      F = 0.D0
      X1Q = Xaux(1)*Xaux(1)
      F = X1Q*X1Q+4.D0*X1Q*Xaux(1)+4.D0*X1Q+Xaux(2)*Xaux(2)

      RETURN

END


subroutine grad(n,x,g)

integer				::n
double precision	::x(n),g(n), Xaux(n)

Xaux(1) = X(1) + 0.25d0*dsqrt(2.d0)
Xaux(2) = X(2) + 0.25d0*dsqrt(2.d0)


g(1) = 4.0d0*Xaux(1)**3+12.0d0*Xaux(1)**2+8.0d0*Xaux(1)

g(2) = 2.0d0*Xaux(2)

return

end subroutine grad