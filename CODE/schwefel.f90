! schwefel function
! -500 <= x <= 500
! fstar = 0
!--------------------
!     SCHWEFEL FUNCTION
!--------------------
      subroutine funct(n,x,f)
      implicit none
      integer k,n
      real*8 f,x(n)
      f=418.9828873d0*n
      do k=1,n
       f=f-x(k)*dsin(dsqrt(dabs(x(k))))
      end do
      return
      end
!-------------
!     GRADIENT
!-------------
      subroutine grad(n,x,g)
      implicit none
      integer k,n
      real*8 x(n),g(n),rd
      do k=1,n
       rd=dsqrt(dabs(x(k)))
       g(k)=-dsin(rd)-0.5d0*rd*dcos(rd)
      end do
      return
      end
