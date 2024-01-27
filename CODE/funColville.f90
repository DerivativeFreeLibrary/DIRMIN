!================================================================
!
! Funzione Colville (n=4) 
!
!================================================================
!				VARIABLE BOUNDS
!================================================================
!       -5.D0                5.D0          Limits of variable 1
!       -5.D0                5.D0          Limits of variable 2
!       -5.D0                5.D0          Limits of variable 3
!       -5.D0                5.D0          Limits of variable 4
! fstar = 0
!================================================================
SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      INTEGER          :: I

      F = 0.0D0
	  F = 100.0D0*(X(1)**2-X(2))**2 + (X(1)-1.0D0)**2 + (X(3)-1.0D0)**2
	  F = F + 90.0D0*(X(3)**2-X(4))**2 + 10.1D0*((X(2)-1.0D0)**2+(X(4)-1.0D0)**2)
	  F = F + 19.8D0*(X(2)-1.0D0)*(X(4)-1.0D0)
      RETURN

END

SUBROUTINE GRAD(N,X,G)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), G(N)

	G(1) =  400.d0*x(1)*(x(1)**2 - x(2)) + 2.d0*(x(1)-1.d0)
	G(2) = -200.d0*(x(1)**2 - x(2)) + 20.2d0*(x(2)-1.d0) + 19.8d0*(x(4)-1.d0)
	G(3) =    2.d0*(x(3)-1.d0) + 180.d0*2.d0*x(3)*(x(3)**2-x(4))
	G(4) = -180.d0*(x(3)**2-x(4)) + 20.2d0*(x(4)-1.d0) + 19.8d0*(x(2)-1.d0)

      RETURN
END
