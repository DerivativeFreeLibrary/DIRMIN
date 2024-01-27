
subroutine funct(n,x,f)

use m_testCE ! modulo che contiene la matrice di rotazione

implicit none

integer			::n, i
double precision::x(n),f, z(n), u(n), summ
double precision, parameter :: pi = 3.1415923d0


do i=1,n
	if (mod(i,2)==0) then
		u(i) = x(i)-1.0d0  !atan(1.0d0)
	else
		u(i) = x(i)-1.0d0  !dexp(1.0d0)
	endif
enddo

z = matmul(M_griew,u)


summ = 0.0d0
do i = 1,n

	summ = summ +z(i)**2-10.0d0*dcos(2.0d0*pi*z(i))+10.0d0

enddo

f = summ-330.0d0


return

end subroutine funct

subroutine grad(n,x,g)
use m_testCE

implicit none

integer			::n,i
double precision::x(n), g(n),  z(n), u(n), gaux(n)
double precision, parameter :: pi = 3.1415923d0

do i=1,n
	if (mod(i,2)==0) then
		u(i) = x(i)-1.0d0 !atan(1.0d0)
	else
		u(i) = x(i)-1.0d0 !dexp(1.0d0)
	endif
enddo

z = matmul(M_griew,u)


do i=1,n

	gaux(i) = 2.0d0*z(i)+20.0d0*pi*dsin(2.d0*pi*z(i))

enddo

g = matmul(transpose(M_griew), gaux)



return

end subroutine grad



