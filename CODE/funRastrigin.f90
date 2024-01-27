!===============================================================
!
! Funzione Rastrigin (n variabile) consigliata n=5
! ALI-TORN POPULATION SET BASED GO ALGORITHMS: SOME MODIFICATIONS AND NUMERICAL STUDIES
! COMPUTERS AND OPERATION RESEARCH 
!===============================================================
!			VARIABLE BOUNDS
!===============================================================
!       -5.12D0              5.12D0          Limits of variable i
! FGLOB = 0.0D0 XSTAR = 0.0D0

subroutine funct(n,x,f)

      implicit none

      integer          :: n
      double precision :: x(n), f, xaux(n)

      double precision, parameter :: pi = 3.1415923d0

      integer          :: i

	  xaux = x + 2.56d0*dsqrt(2.d0)

      f = 10.d0*dble(n)

	  do i = 1,n
		f = f + xaux(i)**2.d0 -10.d0*dcos(2.d0*pi*xaux(i))
	  enddo

      return

end


subroutine grad(n,x,g)

implicit none

integer			::n, i
double precision::x(n), g(n), xaux(n)
double precision, parameter :: pi = 3.1415923d0

xaux = x + 2.56d0*dsqrt(2.d0)

do i=1,n

	g(i) = 2.0d0*xaux(i)+20.0d0*pi*dsin(2.d0*pi*xaux(i))

enddo

return

end subroutine grad
