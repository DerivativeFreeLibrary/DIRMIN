!=========================================
! Dixon & Price function
!
! n= variabile (25 (default))
!        
! valore ottimo della funzione obiettivo = 0.0
!
!         -10.0 =< x(i) =< 10.0   
!====================================================    

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N,I
      DOUBLE PRECISION :: X(N),F

	  F=0.D0
      DO I = 2,N
         F = F+dble(I)*(2.D0*x(I)**2-x(I-1))**2    
      ENDDO
      F = F+(x(1)-1.D0)**2

      RETURN

END

SUBROUTINE DATA(N,L,U,FOTT)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: L(N),U(N),FOTT
      INTEGER          :: I

      DO I=1,N
         L(I) = -10.d0
         U(I) =  10.d0
      ENDDO

      FOTT = 0.D0

      RETURN

END


subroutine grad(n,x,g)


implicit none

integer				::n, i
double precision	::x(n), g(n)

do i=1,n

if (i==1) then

	g(1) = 2.0d0*(x(i)-1.0d0)-2.0d0*dble(i+1)*(2.0d0*(x(2)**2)-x(1))
else
	if (i< n) then
		
		g(i) = 8.0d0*dble(i)*x(i)*(2.0d0*(x(i)**2)-x(i-1))-2.0d0*dble(i+1)*(2.0d0*(x(i+1)**2)-x(i)) 

	else
		
		g(i) = 8.0d0*dble(i)*x(i)*(2.0d0*(x(i)**2)-x(i-1))

	endif
endif

enddo	

return

end subroutine grad


