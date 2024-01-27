! power sum function
! 0 <= x <= n
! fstar = 0
!--------------------
!     POWER SUM FUNCTION
!--------------------
      subroutine funct(n,x,f)
      implicit none
      integer j,k,n
      real*8 f,x(n),s,b(4)
	b(1)=8.d0
	b(2)=18.d0
	b(3)=44.d0
	b(4)=114.d0
      f=0.d0
      do k=1,n
       s=0.d0
       do j=1,n
        s=s+x(j)**k
       end do
       f=f+(s-b(k))**2
      end do
      return
      end
!-------------
!     GRADIENT
!-------------
      subroutine grad(n,x,g)
      implicit none
      integer i,j,k,n
      real*8 x(n),g(n),s1,s,b(4)
	b(1)=8.d0
	b(2)=18.d0
	b(3)=44.d0
	b(4)=114.d0
      s1=0.d0
      do j=1,n
       s1=s1+x(j)
      end do
      do i=1,n
       g(i)=2.d0*(s1-b(1))
       do k=2,n
        s=0.d0
        do j=1,n
         s=s+x(j)**k
        end do
        g(i)=g(i)+2.d0*(s-b(k))*k*x(i)**(k-1)
       end do
      end do
      return
      end
