!=========================================
!
! Funzione di Shekel 10 minimi (n=4) 
!
!=========================================

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      DOUBLE PRECISION :: A(10,4), C(10), FA
      INTEGER          :: I, J

      DO I=1,4
         A(1,I)=4.D0
         A(2,I)=1.D0
         A(3,I)=8.D0
         A(4,I)=6.D0
      END DO
      DO I=1,2
         A(5,2*(I-1)+1)=3.D0
         A(5,2*I)=7.D0
         A(6,2*(I-1)+1)=2.D0
         A(6,2*I)=9.D0
         A(7,I)=5.D0
         A(7,I+2)=3.D0
         A(8,2*(I-1)+1)=8.D0
         A(8,2*I)=1.D0
         A(9,2*(I-1)+1)=6.D0
         A(9,2*I)=2.D0
         A(10,2*(I-1)+1)=7.D0
         A(10,2*I)=3.6D0
      END DO

      C(1)=0.1D0
      C(2)=0.2D0
      C(3)=0.2D0
      C(4)=0.4D0
      C(5)=0.4D0
      C(6)=0.6D0
      C(7)=0.3D0
      C(8)=0.7D0
      C(9)=0.5D0
      C(10)=0.5D0

      F  = 0.0D0
      FA = 0.0D0

      DO I=1,10
         DO J=1,4
            FA = FA +(X(J)-A(I,J))**2
         END DO
         IF ((FA+C(I)).EQ.0.D0) THEN
            F=1.D25
            RETURN
         ENDIF
         F  = F -1.0D0/(FA+C(I))
         FA = 0.0D0
      END DO

      RETURN

END


subroutine grad(n,x,g)


implicit none


integer			::n, i,j
double precision::x(n), g(n)

double precision :: a(10,4), c(10), fa


do i=1,4
 a(1,i)=4.d0
 a(2,i)=1.d0
 a(3,i)=8.d0
 a(4,i)=6.d0
end do
do i=1,2
 a(5,2*(i-1)+1)=3.d0
 a(5,2*i)=7.d0
 a(6,2*(i-1)+1)=2.d0
 a(6,2*i)=9.d0
 a(7,i)=5.d0
 a(7,i+2)=3.d0
 a(8,2*(i-1)+1)=8.d0
 a(8,2*i)=1.d0
 a(9,2*(i-1)+1)=6.d0
 a(9,2*i)=2.d0
 a(10,2*(i-1)+1)=7.d0
 a(10,2*i)=3.6d0
end do

c(1)=0.1d0
c(2)=0.2d0
c(3)=0.2d0
c(4)=0.4d0
c(5)=0.4d0
c(6)=0.6d0
c(7)=0.3d0
c(8)=0.7d0
c(9)=0.5d0
c(10)=0.5d0

fa = 0.0d0
g  = 0.0d0

do i=1,10

 do j=1,4
    fa = fa +(x(j)-a(i,j))**2
 end do

do j = 1,4
 g(j) = g(j)+ 2.0d0*(x(j)-a(i,j))/(fa+c(i))**2 
enddo

fa = 0.0d0


end do

return


end subroutine grad