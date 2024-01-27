! powell function
! -4 <= x <= 5
! fstar = 0
!--------------------
!     POWELL FUNCTION
!--------------------
      subroutine funct(n,x,f)
      implicit none
      integer j,k,n,m
      real*8 f,x(n),y(n/4,4)
      m=n/4
      f=0.d0
      do k=1,m
       do j=1,4
        y(k,j)=x(j+4*(k-1))
       end do
       f=f+(y(k,1)+10.d0*y(k,2))**2+5.d0*(y(k,3)-y(k,4))**2  &
          +(y(k,2)-2.d0*y(k,3))**4+10.d0*(y(k,1)-y(k,4))**4
      end do
      return
      end
!-------------
!     GRADIENT
!-------------
      subroutine grad(n,x,g)
      implicit none
      integer i,j,k,n,m
      real*8 x(n),g(n),y(n/4,4),gy(4)
      m=n/4
      do k=1,m
       do j=1,4
        y(k,j)=x(j+4*(k-1))
       end do
       gy(1)=2.d0*(y(k,1)+10.d0*y(k,2))+40.d0*(y(k,1)-y(k,4))**3
       gy(2)=20.d0*(y(k,1)+10.d0*y(k,2))+4.d0*(y(k,2)-2.d0*y(k,3))**3
       gy(3)=10.d0*(y(k,3)-y(k,4))-8.d0*(y(k,2)-2.d0*y(k,3))**3
       gy(4)=-10.d0*(y(k,3)-y(k,4))-40.d0*(y(k,1)-y(k,4))**3
       do j=1,4
        g(j+4*(k-1))=gy(j)
       end do
      end do
      return
      end
