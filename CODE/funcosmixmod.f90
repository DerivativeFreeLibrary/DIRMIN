!=========================================
!
! Funzione Cosine Mixture
!
!=========================================

subroutine funct(n,x,f)

      implicit none

      integer          :: n
      double precision :: x(n), f, xaux(n)

      double precision, parameter :: pi = 3.1415923

      integer          :: i

	  xaux = x + 0.5d0*dsqrt(2.d0)

      f = 0.0d0

      do i=1,n
         f = f + 0.1d0*dcos(5.0d0*pi*xaux(i))-xaux(i)**2
      end do
      
      f = - f

      return

end subroutine funct



subroutine grad(n,x,g)

implicit none

integer				::n,i
double precision	::x(n), g(n), xaux(n)

double precision, parameter :: pi = 3.1415923

 xaux = x + 0.5d0*dsqrt(2.d0)

do i = 1,n

	g(i) =0.1d0*sin(5.0d0*pi*xaux(i))*5.0d0*pi+2.0d0*xaux(i)

enddo

return

end subroutine grad
