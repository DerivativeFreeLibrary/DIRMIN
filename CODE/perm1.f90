! Perm function 1
! -n <= x <= n
! fstar = 0
!--------------------
!     PERM1 FUNCTION
!--------------------
      subroutine funct(n,x,f)
      implicit none
      integer j,k,n
      real*8 f,x(n),s
      f=0.d0
      do k=1,n
       s=0.d0
       do j=1,n
	   !write(*,*)j,k,(x(j)/j)**k
        s=s+(dble(j)**dble(k)+0.5d0)*((x(j)/j)**k-1.d0)
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
       s1=s1+(j+0.5d0)*((x(j)/j)-1.d0)
      end do
      do i=1,n
       g(i)=2.d0*s1*(i+0.5d0)/i
       do k=2,n
        s=0.d0
        do j=1,n
         s=s+(j**k+0.5d0)*((x(j)/j)**k-1.d0)
        end do
        g(i)=g(i)+2.d0*s*k*(i**k+0.5d0)*((x(i)/i)**(k-1)/i)
       end do
      end do
      return
      end
