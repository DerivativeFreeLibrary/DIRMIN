
subroutine funct(n,x,f)

use m_testCE ! modulo che contiene la matrice di rotazione

implicit none

integer			::n, i
double precision::x(n),f, z(n), u(n), summ, prod

do i=1,n
	if (mod(i,2)==0) then
		u(i) = x(i)-atan(1.0d0)
	else
		u(i) = x(i)-dexp(1.0d0)
	endif
enddo

!write(*,*) 'qui',M_griew
z = matmul(M_griew,u)

summ = 0.0d0
prod = 1.0d0
do i = 1,n

	summ = summ +z(i)**2/4.0d3
	prod = prod*dcos(z(i)/dsqrt(dble(i)))

enddo

f = summ-prod+1.0d0-180.0d0


return

end subroutine funct

subroutine grad(n,x,g)

use m_testCE
implicit none

integer				::n,i,j
double precision	::g(n),x(n), xaux(n), u(n), prod(n), gaux(n)


do i=1,n
	if (mod(i,2)==0) then
		u(i) = x(i)-atan(1.0d0)
	else
		u(i) = x(i)-dexp(1.0d0)
	endif
enddo

xaux = matmul(M_griew,u)

prod = 1.0d0

do i=1,n
	do j= 1,n
		if (i.ne.j) then
			prod(i) = prod(i)*dcos(xaux(j)/dsqrt(dble(j))) 
		endif
	enddo
end do


do i =1,n

	gaux(i) = (2.0d0*xaux(i))/4.0d3+(1/dsqrt(dble(i)))*dsin(xaux(i)/dsqrt(dble(i)))*prod(i)

enddo

g = matmul(transpose(M_griew), gaux)

return

end subroutine






