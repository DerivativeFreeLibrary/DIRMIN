!=========================================
! Michalewics function
!
! n= variabile
!   
! molti minimi locali       
!
! valore ottimo della funzione obiettivo
!   
!       n=2, f(x*) = -1.8013.
!       n=5, f(x*) = -4.687658.
!       n=10, f(x*) = -9.66015.
!
! punto ottimo = ?
!
!         -pi =< x(i) =< pi   
!====================================================    

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      DOUBLE PRECISION, PARAMETER :: PI = 3.1415923D0

      DOUBLE PRECISION :: s
      INTEGER          :: I,M

      m = 10

      s = 0.d0
      DO i = 1,n
         s = s+dsin(x(i))*(dsin(dble(i)*x(i)**2/pi))**(2*m)
      ENDDO

      F = -s

      RETURN

END

SUBROUTINE DATA(N,L,U,FOTT)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: L(N),U(N),FOTT
      INTEGER          :: I
      DOUBLE PRECISION, PARAMETER :: PI = 3.1415923D0

      DO I=1,N
         L(I) = -PI
         U(I) =  PI
      ENDDO

      IF(N.EQ.2) FOTT = -1.8013D0
	  IF(N.EQ.5) FOTT = -4.687658D0
	  IF(N.EQ.10) FOTT = -9.66015D0

	  IF((N.NE.2).AND.(N.NE.5).AND.(N.NE.10)) THEN
	    WRITE(*,*) ' funct. Michalewics non definito funct. ott.'
		STOP
      ENDIF

      RETURN

END


subroutine grad(n,x,g)


implicit none

integer				::n,i,m
double precision	::x(n), g(n)
double precision, parameter :: pi = 3.1415923d0

m = 10
do i=1,n

	g(i)=-dcos(x(i))*(dsin(dble(i)*x(i)**2/pi))**(2*m) &
	-dsin(x(i))*4.0d0*dble(m)*dble(i)*x(i)*(dsin(dble(i)*(x(i)**2)/pi))**(2*m-1)*dcos(dble(i)*(x(i)**2)/pi)/pi

enddo


return


end subroutine