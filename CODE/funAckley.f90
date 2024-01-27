!=========================================
! Ackley function
!
! n=2 (default)
!   
! molti minimi locali       
!
! valore ottimo della funzione obiettivo = 0.0
!
! punto ottimo = 0
!
!         -15.0 =< x(i) =< 30.0   
!====================================================    

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      DOUBLE PRECISION, PARAMETER :: PI = 3.1415923D0

      DOUBLE PRECISION :: a,b,c,s1,s2
      INTEGER          :: I

      a = 20.d0
	  b = 0.2d0
	  c = 2.d0*pi

      s1 = 0.d0 
	  s2 = 0.d0

      DO I=1,N
         s1 = s1+x(i)*x(i)
         s2 = s2+dcos(c*x(i))
      ENDDO

      F = -a*dexp(-b*sqrt((1.d0/dble(n))*s1))-dexp((1.d0/dble(n))*s2)+a+dexp(1.d0)

      RETURN

END



subroutine grad(n,x,g)

implicit none

integer			::n,i
double precision::x(n), g(n)
double precision :: a,b,c,s1,s2
double precision, parameter :: pi = 3.1415923d0

a = 20.d0
b = 0.2d0
c = 2.d0*pi

s1 = 0.d0 
s2 = 0.d0

do i=1,n
 s1 = s1+x(i)*x(i)
 s2 = s2+dcos(c*x(i))
enddo

do i = 1,n
     
	g(i) = a*b*dexp(-b*sqrt((1.d0/dble(n))*s1)) * (1.d0/(dsqrt((1.d0/dble(n))*s1)))&
	 *x(i)*1.d0/dble(n) + dexp((1.d0/dble(n))*s2)*dsin(c*x(i))*c*(1.d0/dble(n))
	
enddo
return

end subroutine grad


SUBROUTINE DATA(N,L,U,FOTT)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: L(N),U(N),FOTT
      INTEGER          :: I

      DO I=1,N
         L(I) = -15.d0
         U(I) =  30.d0
      ENDDO

      FOTT = 0.D0

      RETURN

END