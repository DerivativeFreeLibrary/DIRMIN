
!module Test_Function

!integer							:: n, l2, l3, n_func_eval, n_func_grad_eval, dim, dim_m
!double precision				:: h, k, a2, c1, c2, funzione
!double precision, allocatable	:: alpha0(:,:),alpha2(:,:),alpha3(:,:), beta0(:,:), beta2(:,:), beta3(:,:),&
!							   df(:,:),du(:,:), f(:), u(:), kk(:), aa1(:), dw(:,:),w(:), gradiente(:), &
!							   o1(:), do1(:), add(:), dadd(:), add_g(:), dadd_g(:), dg(:), dtmp(:), go_point(:)
 
!integer, allocatable			:: m(:), p(:,:)

!integer, parameter				:: MAX_N=1000
!double precision, parameter		:: A1=12.0/25.0
!double precision, parameter		:: AA2=16.0/125.0

!end module Test_Function

!!!!!!!!!!!!!!!!!!!!!!!!attenzione x il globale direct type: il modulo potrebbe già essere incluso nel progetto dentro modules.f90!!!!!!!!!!!!!!!!!!!


double precision function GammaComputation(x,level3, grad)

use Test_Function

implicit none

double precision	:: x(dim),g, aux, aux1, aux2, aux3, grad(dim)
integer				:: level3,i,j,t

!write(*,*) 'level3 = ',level3

do i=1,n
	w(i) = 0.0d0
	do j=1,n
		
		w(i) = w(i)+dw(i,j)*x(j)
	enddo
	!write(*,*) 'w ',i,'=',w(i)
enddo

do i=1,n
	
	o1(i) = -h*cos(aa1(i)*(w(i)-c1))+h
	do1(i) = h*aa1(i)*sin(aa1(i)*(w(i)-c1))

enddo

do i=n+1,n+dim_m-2

	add(i-n) = (x(i)-2.5d0)*(x(i)-2.5d0)-h*cos(a2*(x(i)-2.5d0))+h
	dadd(i-n) = 2.0d0*(x(i)-2.5d0)+h*a2*sin(a2*(x(i)-2.5d0))
	!write(*,*) ' add',i-n,'=', add(i-n)
	!write(*,*) 'dadd',i-n,'=',dadd(i-n)
	!pause
enddo
!write(*,*) 'm=',m
do i=1,dim_m

	f(i)=0.0d0
	do j=1,i-2
		f(i) = f(i)+add(j)
		df(i,n+j) = dadd(j)
		!write(*,*) 'i=',i,' j=',j
		!pause
	enddo
	do j=1,m(i)
		
		if (w(j)<0.0d0) then
			
			f(i) = f(i)+alpha0(level3,j)+alpha2(level3,j)*(w(j)-c1)*(w(j)-c1)&
			+alpha3(level3,j)*(w(j)-c1)*(w(j)-c1)*(w(j)-c1)+o1(j)
			dtmp(j) = 2.0d0*alpha2(level3,j)*(w(j)-c1)+3.0d0*alpha3(level3,j)*(w(j)-c1)*(w(j)-c1)+do1(j)

		else

			f(i) = f(i)+beta0(level3,j)+beta2(level3,j)*(w(j)-c2)*(w(j)-c2)&
			+beta3(level3,j)*(w(j)-c2)*(w(j)-c2)*(w(j)-c2)+o1(j)
			dtmp(j) = 2.0d0*beta2(level3,j)*(w(j)-c2)+3.0d0*beta3(level3,j)*(w(j)-c2)*(w(j)-c2)+do1(j)
			
		endif
		!write(*,*) 'dtmp ',j,'=',dtmp(j)

	enddo
	do j=m(i)+1,n

		if (p(level3,j)==1) then
			
			f(i) = f(i)+0.5d0*(w(j)-c1)*(w(j)-c1)+o1(j)+2.0d0
			dtmp(j) = (w(j)-c1)+do1(j)
		else
			
			f(i) = f(i)+0.5d0*(w(j)-c2)*(w(j)-c2)+o1(j)+2.0d0
			dtmp(j) = (w(j)-c2)+do1(j)
		
		endif

		!write(*,*) 'dtmp ',j,'=',dtmp(j)
	enddo
	do j=1,n

		df(i,j) = dtmp(1)*dw(1,j)
		do t=2,n

			df(i,j) = df(i,j)+dtmp(t)*dw(t,j)
		
		enddo
			!write(*,*)'df(',i,',',j,')=',df(i,j)
!			pause

	enddo
	!write(*,*) 'f(',i,')=',f(i)
	
enddo
!pause

do i=1,n

	dg(i) = df(1,i)
	!write(*,*) 'dg[',i,']=',dg(i)

enddo

g=f(1)

do i=2,dim_m
!	write(*,*) 'x[',n+i-1,']=',x(n+i-1)
!	write(*,*) 'g=',g,' a2=',a2,' a1=',a1
!	write(*,*) 'f(',i,')=',f(i)
	if (x(n+i-1).le.0.0d0) then

		aux = (g+2.0d0*f(i))*(x(n+i-1)+2.5d0)
		aux1 = -(g+f(i))*cos(a2*(x(n+i-1)+2.5d0))+f(i)+g
		aux2 = (g+f(i))*a2*sin(a2*(x(n+i-1)+2.5d0))
		g = g+a1*aux*(x(n+i-1)+2.5d0)-AA2*aux*(x(n+i-1)+2.5d0)*(x(n+i-1)+2.5d0)+aux1
		dg(n+i-1) = 2.0d0*a1*aux-3.0d0*AA2*aux*(x(n+i-1)+2.50d0)+aux2
		do j=1,n+i-1-1
			
			aux3 = (dg(j)+2.0d0*df(i,j))*(x(n+i-1)+2.50d0)
			dg(j) = -AA2*aux3*(x(n+i-1)+2.50d0)*(x(n+i-1)+2.50d0)+a1*aux3*(x(n+i-1)+2.50d0)&
				+dg(j)+(df(i,j)+dg(j))*(1.0d0-cos(a2*(x(n+i-1)+2.50d0)))

		enddo

	else

		aux = (2.0d0*g+f(i))*(x(n+i-1)-2.5d0)
		aux1 = -(g+f(i))*cos(a2*(x(n+i-1)-2.5d0))+f(i)+g
		aux2 = (g+f(i))*a2*sin(a2*(x(n+i-1)-2.5d0))
		g = f(i)+a1*aux*(x(n+i-1)-2.5d0)+AA2*aux*(x(n+i-1)-2.5d0)*(x(n+i-1)-2.5d0)+aux1
		dg(n+i-1) = 2.0d0*a1*aux+3.0d0*AA2*aux*(x(n+i-1)-2.50d0)+aux2
		do j=1,n+i-1-1
			
			aux3 = (2.0d0*dg(j)+df(i,j))*(x(n+i-1)-2.50d0)
			dg(j) = AA2*aux3*(x(n+i-1)-2.50d0)*(x(n+i-1)-2.50d0)+a1*aux3*(x(n+i-1)-2.50d0)&
				+df(i,j)+(df(i,j)+dg(j))*(1.0d0-cos(a2*(x(n+i-1)-2.50d0)))

		enddo

	endif

enddo


do i=1,n+dim_m-1
	
	grad(i) = dg(i)
!	write(*,*) 'grad ',i,'=',grad(i)
enddo

GammaComputation = g

return

end function GammaComputation


subroutine funct_grad(x)

use Test_Function

implicit none


double precision		:: x(dim), aux, aux1, aux2, aux3, GammaComputation
integer					::i,j

n_func_eval = n_func_eval+1

n_func_grad_eval = n_func_grad_eval+1

do i = 1,l3

	u(i) = GammaComputation(x,i,du(i,1:dim))
!	write(*,*)'u(',i,')=', u(i)
!	do j = 1,dim
!		write(*,*)'du(',i,j,')=',du(i,j)	
!	enddo
enddo

!pause

do i=n+dim_m+1,n+dim_m+l3-2

	add_g(i-n-dim_m)=(x(i-1)-2.5d0)*(x(i-1)-2.5d0)-h*cos(a2*(x(i-1)-2.5d0))+h
	dadd_g(i-n-dim_m)=2.0d0*(x(i-1)-2.5d0)+h*a2*sin(a2*(x(i-1)-2.5d0))
	!write(*,*)'add_g(i-n-dim_m)=', add_g(i-n-dim_m) 
	!write(*,*)'dadd_g(i-n-dim_m)=', dadd_g(i-n-dim_m) 
enddo

do i=3,l3

	do j=1,i-2
		
		du(i,n+dim_m-1+j) = dadd_g(j)
		u(i) = u(i)+add_g(j)
		
	enddo

enddo

do i=1,n+dim_m-1

	gradiente(i) = du(1,i)
	!write(*,*)'du(1,i)= ', du(1,i)

enddo

funzione = u(1)

do i=2,l3
	
	if (x(n+dim_m+i-2).le.0.0d0) then
		
		aux = (funzione +2.0d0*u(i)+2.0d0-(1.0d0/(dble(l3))))*(x(n+dim_m+i-2)+2.50d0)
		aux1 = -(u(i)+funzione)*cos(a2*(x(n+dim_m+i-2)+2.50d0))+(u(i)+funzione)
		aux2 = (u(i)+funzione)*a2*sin(a2*(x(n+dim_m+i-2)+2.50d0))
		funzione = funzione+(1.0d0/dble(l3))+a1*aux*(x(n+dim_m+i-2)+2.50d0)&
		-AA2*aux*(x(n+dim_m+i-2)+2.50d0)*(x(n+dim_m+i-2)+2.50d0)+aux1
		gradiente(n+dim_m+i-2) = 2.0d0*a1*aux-3.0d0*AA2*aux*(x(n+dim_m+i-2)+2.50d0)+aux2
		do j=1,n+dim_m+i-3
			
			aux3 = (gradiente(j)+2.0d0*du(i,j))*(x(n+dim_m+i-2)+2.50d0)
			gradiente(j)=gradiente(j)-AA2*aux3*(x(n+dim_m+i-2)+2.50d0)*(x(n+dim_m+i-2)+2.50d0)&
			+a1*aux3*(x(n+dim_m+i-2)+2.50d0)+(gradiente(j)+du(i,j))*(1.0d0-cos(a2*(x(n+dim_m+i-2)+2.50d0)))
			
		enddo

	else
		
		aux = (2.0d0*funzione+u(i)+2.0d0)*(x(n+dim_m+i-2)-2.50d0)
		aux1 = -(u(i)+funzione)*cos(a2*(x(n+dim_m+i-2)-2.50d0))+(u(i)+funzione)
		aux2 = (u(i)+funzione)*a2*sin(a2*(x(n+dim_m+i-2)-2.50d0))
		funzione = u(i)+a1*aux*(x(n+dim_m+i-2)-2.50d0)&
		+AA2*aux*(x(n+dim_m+i-2)-2.50d0)*(x(n+dim_m+i-2)-2.50d0)+aux1
		gradiente(n+dim_m+i-2) = 2.0d0*a1*aux+3.0d0*AA2*aux*(x(n+dim_m+i-2)-2.50d0)+aux2
		do j=1,n+dim_m+i-3
			
			aux3 = (2.0d0*gradiente(j)+du(i,j))*(x(n+dim_m+i-2)-2.50d0)
			gradiente(j)=du(i,j)+AA2*aux3*(x(n+dim_m+i-2)-2.50d0)*(x(n+dim_m+i-2)-2.50d0)&
			+a1*aux3*(x(n+dim_m+i-2)-2.50d0)+(gradiente(j)+du(i,j))*(1.0d0-cos(a2*(x(n+dim_m+i-2)-2.50d0)))
			
		enddo

	endif

enddo


return

end subroutine funct_grad




subroutine funct(nn,x,f_ret)
use Test_function
implicit none

integer:: nn

double precision::x(n), f_ret

call funct_grad(x)

f_ret = funzione

return




end subroutine funct


subroutine grad(nn,x,g_ret)
use Test_function
implicit none

integer:: nn

double precision::x(n), g_ret(n)

call funct_grad(x)

g_ret = gradiente

return




end subroutine grad

