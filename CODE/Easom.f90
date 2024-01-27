!=========================================
! Easom function
!
! n=2 (default)
!   
! molti minimi locali       
!
! valore ottimo della funzione obiettivo = -1.0
!
! punto ottimo = (pi , pi), 
!
!         -100.0 =< x(i) =< 100.0   
!====================================================    

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      DOUBLE PRECISION, PARAMETER :: PI = 3.1415923D0

		!write(*,*) 'funct:',x
      F = -dcos(x(1))*dcos(x(2))*dexp(-(x(1)-pi)**2-(x(2)-pi)**2)
		!write(*,*) 'fine funct'
      RETURN
	   
END

SUBROUTINE DATA(N,L,U,FOTT)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: L(N),U(N),FOTT
      INTEGER          :: I

      DO I=1,N
         L(I) = -100.d0
         U(I) =  100.d0
      ENDDO

      FOTT = -1.D0

      RETURN

END

subroutine grad(n,x,g)


implicit none

integer				::n
double precision	::x(n),g(n)
double precision, parameter :: pi = 3.1415923d0

		!write(*,*) 'grad:',x

g(1) = dsin(x(1))*dcos(x(2))*dexp(-(x(1)-pi)**2-(x(2)-pi)**2)-&
	dcos(x(1))*dcos(x(2))*dexp(-(x(1)-pi)**2-(x(2)-pi)**2)*(-2.0d0*(x(1)-pi))

g(2) = dcos(x(1))*dsin(x(2))*dexp(-(x(1)-pi)**2-(x(2)-pi)**2)-&
	dcos(x(1))*dcos(x(2))*dexp(-(x(1)-pi)**2-(x(2)-pi)**2)*(-2.0d0*(x(2)-pi))


		!write(*,*) 'fine grad'

return

end subroutine grad