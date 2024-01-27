module mod_type
	type intervallo
		real*8, allocatable			:: cent(:), dimen(:), lbs(:), ubs(:), xbars(:)
		real*8					:: fint, diam, maxdim, der
		logical					:: flagloc, flagdiv, flagcon, flagopt
		integer					:: id
		type(intervallo), pointer	:: next, pred
	end type intervallo
	type vertice
		type(intervallo), pointer	:: int
		type(vertice),    pointer	:: next
	end type vertice
	type colonna
		real*8						:: diam
		type(colonna),    pointer	:: next, pred
		type(intervallo), pointer	:: int
	end type colonna
	type fpunt
		real*8						:: f
		type(intervallo), pointer	:: punt
	end type fpunt
end module mod_type

module mod_globale
	real*8, allocatable				:: globxbest(:)
	real*8						:: globfbest
	integer						:: globnftot, globnf, globmaxnf
end module mod_globale

module mod_box
	real*8, allocatable				:: lb(:), ub(:), xbar(:)
	real*8, allocatable				:: lbs(:), ubs(:)
	real*8, allocatable				:: xtemp(:), ytemp(:)
	real*8							:: ampiezza, fattore  !metà dell'ampiezza
end module mod_box

module mod_suddividi
	real*8,  allocatable			:: vetf1(:), vetf2(:)
	real*8,  allocatable			:: xsud(:), ysud(:)
	logical, allocatable			:: mask(:)
end module mod_suddividi

module mod_mem
	logical							:: memerror
	integer							:: num_el_L
end module mod_mem

module Test_Function

integer							:: n, l2, l3, n_func_eval, n_func_grad_eval, dim, dim_m
double precision				:: h, k, a2, c1, c2, funzione
double precision, allocatable	:: alpha0(:,:),alpha2(:,:),alpha3(:,:), beta0(:,:), beta2(:,:), beta3(:,:),&
								   df(:,:),du(:,:), f(:), u(:), kk(:), aa1(:), dw(:,:),w(:), gradiente(:), &
								   o1(:), do1(:), add(:), dadd(:), add_g(:), dadd_g(:), dg(:), dtmp(:), go_point(:)
 
integer, allocatable			:: m(:), p(:,:)

integer, parameter				:: MAX_N=1000
double precision, parameter		:: A1=12.0/25.0
double precision, parameter		:: AA2=16.0/125.0

end module Test_Function

module m_testCE

double precision, allocatable:: M_griew(:,:)


end module m_testCE
