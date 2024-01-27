!=========================================
!
! Funzione Esponenziale
!
!=========================================

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F, xaux(n)

      DOUBLE PRECISION :: F1
      INTEGER          :: I

	  XAUX = X + 0.5d0*dsqrt(2.d0)

      F1 = 0.0D0

      DO I=1,N
         F1 = F1+XAUX(I)**2
      END DO

      F = -DEXP(-0.5D0*F1)

      RETURN

END


subroutine grad(n,x,g)

implicit none

integer			::n,i
double precision::x(n),g(n),xaux(n),f1

xaux = x + 0.5d0*dsqrt(2.d0)

f1 = 0.0d0

do i=1,n
 f1 = f1+xaux(i)**2
end do

do i=1,n

	g(i) = dexp(-0.5d0*f1)*xaux(i)

enddo

return

end subroutine grad
