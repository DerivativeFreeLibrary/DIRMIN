!=========================================
!
! Funzione di Goldstein e Price (n=2)
!
!=========================================

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      F = (1.0d0+(X(1)+X(2)+1.0d0)**2 *(19.0D0-14.0D0*X(1)+3.0D0*(X(1)**2) &
          -14.0D0*X(2)+6.0D0*X(1)*X(2)+3.0D0*(X(2)**2) ))  * (30.0D0+(2.0D0*X(1) &
          -3.0D0*X(2))**2*(18.0D0-32.0D0*X(1)+12.0D0*(X(1)**2)+48.0D0*X(2) &
          -36.0D0*X(1)*X(2)+27.0D0*(X(2)**2)))

      RETURN

END



subroutine grad(n,x,g)


implicit none

integer				::n
double precision	::x(n),g(n)
double precision	::A, B, C, D, E, F


C = X(1)+X(2)+1.0d0

D =19.0D0-14.0D0*X(1)+3.0D0*(X(1)**2)-14.0D0*X(2)+6.0D0*X(1)*X(2)+3.0D0*(X(2)**2)

E = 2.0D0*X(1)-3.0D0*X(2)

F = 18.0D0-32.0D0*X(1)+12.0D0*(X(1)**2)+48.0D0*X(2) &
          -36.0D0*X(1)*X(2)+27.0D0*(X(2)**2)
A = 1.0d0+(C**2)*D

B = 30.0d0+(E**2)*F

g(1) = (2.0d0*(C*D)+(C**2)*(-14.0d0 + 6.0d0*x(1) + 6.0d0*x(2) ) ) * B + &
	A*(4.0d0*(E*F)+(E**2)*(-32.0d0+24.0d0*x(1)-36.0d0*x(2)))


g(2) = (2.0d0*C*D+C**2*(-14.0d0+6.0d0*x(1)+6.0d0*x(2)))*B+&
	A*(-6.0d0*E*F+(E**2)*(48.0d0-36.0d0*x(1)+54.0d0*x(2)))


return

end subroutine grad