!=========================================
!
! Funzione di HARTMAN (n=3)
!
!=========================================

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      DOUBLE PRECISION :: A(4,3), C(4), P(4,3), FA
      INTEGER          :: I, J

        A(1,1)=3.D0
        A(2,1)=0.1D0
        A(3,1)=3.D0
        A(4,1)=0.1D0

        A(1,2)=10.D0
        A(2,2)=10.D0
        A(3,2)=10.D0
        A(4,2)=10.D0

        A(1,3)=30.D0
        A(2,3)=35.D0
        A(3,3)=30.D0
        A(4,3)=35.D0

        C(1)=1.D0
        C(2)=1.2D0
        C(3)=3.D0
        C(4)=3.2D0

        P(1,1)=0.3689D0
        P(2,1)=0.4699D0
        P(3,1)=0.1091D0
        P(4,1)=0.03815D0

        P(1,2)=0.1170D0
        P(2,2)=0.4387D0
        P(3,2)=0.8732D0
        P(4,2)=0.5743D0

        P(1,3)=0.2673D0
        P(2,3)=0.7470D0
        P(3,3)=0.5547D0
        P(4,3)=0.8828D0

        F  = 0.D0
        FA = 0.D0

        DO I=1,4
           DO J=1,3
              FA = FA + A(I,J)*((X(J)-P(I,J))**2)
           END DO
           F  = F - C(I)*DEXP(-FA)
           FA = 0.D0
        END DO

      RETURN

END



subroutine grad(n,x,g)

implicit none

integer			::n, i, j
double precision::x(n), g(n)
double precision :: a(4,3), c(4), p(4,3), fa(4)

a(1,1)=3.d0
a(2,1)=0.1d0
a(3,1)=3.d0
a(4,1)=0.1d0

a(1,2)=10.d0
a(2,2)=10.d0
a(3,2)=10.d0
a(4,2)=10.d0

a(1,3)=30.d0
a(2,3)=35.d0
a(3,3)=30.d0
a(4,3)=35.d0

c(1)=1.d0
c(2)=1.2d0
c(3)=3.d0
c(4)=3.2d0

p(1,1)=0.3689d0
p(2,1)=0.4699d0
p(3,1)=0.1091d0
p(4,1)=0.03815d0

p(1,2)=0.1170d0
p(2,2)=0.4387d0
p(3,2)=0.8732d0
p(4,2)=0.5743d0

p(1,3)=0.2673d0
p(2,3)=0.7470d0
p(3,3)=0.5547d0
p(4,3)=0.8828d0

!write(*,*) 'x=',x
if(isnan(x(1)).or.isnan(x(2)).or.isnan(x(3))) pause
fa = 0.0d0
g = 0.d0
do i=1,4
   do j=1,3
      fa(i) = fa(i) + a(i,j)*((x(j)-p(i,j))**2)
   end do
end do
!write(*,*) 'fa=',fa
do i = 1,4
	
	do j = 1,3
		g(j) = g(j)+c(i)*dexp(-fa(i))*2.0d0*a(i,j)*(x(j)-p(i,j)) 
	enddo

enddo

return

end subroutine grad