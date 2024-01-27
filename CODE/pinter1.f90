!===============================================================
!
! FUNZIONE DI PINTER (n=2,5,10,20,....)
! HANDBOOK OF GLOBAL OPTIMIZATION VOLUME 2 CAPITOLO 15 DI PINTER
!
!===============================================================
!			VARIABLE BOUNDS
!===============================================================
!       -5.D+0               5.D+0        Limits of variable i

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: XSTAR(N), X(N), F, P1, P2

      INTEGER          :: I

	  DO I=1,N
		IF (MOD(I,2) == 0) THEN
			XSTAR(I) = (DBLE(I)/DBLE(N))*4.0D0-0.5
		ELSE
			XSTAR(I) = -(DBLE(I)/DBLE(N))*4.0D0+0.5
		ENDIF

	  ENDDO
	  

      F = 0.D0

      DO I=1,N
         F = F+0.025D0*DBLE(N)*(X(I)-XSTAR(I))**2
      END DO

      P1 = 0.0D0
      P2 = 0.0D0
      DO I=1,N
         P1 = P1 + (X(I)-XSTAR(I)) + (X(I)-XSTAR(I))**2
         P2 = P2 + (X(I)-XSTAR(I)) 
	  END DO

	  F = F+DSIN(P1)**2+DSIN(P2)**2
	  
	  	
      RETURN

END

subroutine grad(n,x,g)


implicit none


integer				::n, i
double precision	::x(n), g(n), s, xstar(n), p, pq

do i=1,n
if (mod(i,2) == 0) then
	xstar(i) = (dble(i)/dble(n))*4.0d0-0.5
else
	xstar(i) = -(dble(i)/dble(n))*4.0d0+0.5
endif

enddo

p = 0.0d0
pq = 0.0d0

do i = 1,n
	
	p = p+x(i)-xstar(i)
	pq = pq + (x(i)-xstar(i))**2

enddo

s = 0.025d0*dble(n)
do i=1,n

	g(i) = 2.0d0*s*(x(i)-xstar(i))+2.0d0*dsin(p+pq)*dcos(p+pq)*(1.0d0+2.0d0*(x(i)-xstar(i))) &
			+2.0d0*dsin(p)*dcos(p)

enddo


return

end subroutine grad


