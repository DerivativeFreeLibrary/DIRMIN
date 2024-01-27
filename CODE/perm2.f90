! Perm function 2
! -n <= x <= n
! fstar = 0
!--------------------
!     PERM2 FUNCTION
!--------------------
      subroutine funct(n,x,f)
      implicit none
      integer j,k,n
      real*8 f,x(n),s
      f=0.d0
      do k=1,n
       s=0.d0
       do j=1,n
        s=s+(j+10.d0)*(x(j)**k-(1.d0/j)**k)
       end do
       f=f+s**2
      end do
      return
      end
!-------------
!     GRADIENT
!-------------
      subroutine grad(n,x,g)
      implicit none
      integer i,j,k,n
      real*8 x(n),g(n),s1,s
      s1=0.d0
      do j=1,n
       s1=s1+(j+10.d0)*(x(j)-1.d0/j)
      end do
      do i=1,n
       g(i)=2.d0*s1*(i+10.d0)
       do k=2,n
        s=0.d0
        do j=1,n
         s=s+(j+10.d0)*(x(j)**k-(1.d0/j)**k)
        end do
        g(i)=g(i)+2.d0*s*k*(i+10.d0)*x(i)**(k-1)
       end do
      end do
      return
      end
