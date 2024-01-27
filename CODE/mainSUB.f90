!============================================================================================
!    DIRMIN - DIRECT Algorithm with derivative-free local minimizations.  
!    A Derivative-Free algorithm for bound  constrained global optimization problems 
!    proposed by G.Liuzzi, S.Lucidi, V.Piccialli (see Ref. below)
!
!    Copyright (C) 2014  G.Liuzzi, S.Lucidi, V.Piccialli
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    G.Liuzzi, S.Lucidi, V.Piccialli: Exploiting derivative-free local searches in 
!    DIRECT-type algorithms for global optimization. 
!    Computational Optimization and Applications, to appear (2014)
!
!==============================================================================================
subroutine nuovoglobale(n,lbb,ubb,xbest,fbest,nf,nminloc,fglob,maxint,maxnf,	      &
			iprint,trovato,nint,mindiam,maxdiam,maxL,xmin,fmin,imin,jmin,	&
			alfa_stop,tolglob,trigLS)
	use mod_type
	use mod_box
	use mod_globale
	use mod_suddividi 
	use m_testCE 
	use mod_mem
	implicit none

	interface
		subroutine aggiorna_struttura(root,Ltilde,fdir,nint)
			use mod_type
			implicit none

			type(colonna),pointer		:: root
			real*8				:: Ltilde,fdir
			integer			:: nint
		end subroutine aggiorna_struttura
		subroutine genera_partizione(root,n,tol,iprint)
			use mod_type
			implicit none

			type(colonna),pointer		:: root
			integer			:: n,iprint
			real*8				:: tol
		end subroutine genera_partizione
		subroutine elimina_colonna(currcol,num)
			use mod_type
			implicit none
			type(colonna),pointer		:: currcol
			integer			:: num
		end subroutine elimina_colonna
		subroutine ricintervallo_dx(root,convexhull,nconv,iprint,Ltilde,eps,eps_cvx,fmin)
			use mod_type
			use mod_mem
			implicit none
 
			type(colonna),pointer		:: root
			type(vertice),   pointer	:: convexhull
			real*8				:: Ltilde, eps, eps_cvx
			real*8				:: fmin
			integer			:: iprint, nconv
		end subroutine ricintervallo_dx
		subroutine ricintervallo(root,convexhull,nconv,iprint,maxL,eps_cvx)
			use mod_type
			implicit none

			type(colonna),pointer		:: root
			type(vertice),pointer		:: convexhull
			integer			:: nconv,iprint
			real*8				:: maxL, eps_cvx
		end subroutine ricintervallo
		subroutine riduciconvexhull(convexhull,nelim,eps,toldiam,fmin)
			use mod_type
			implicit none

			type(vertice),pointer		:: convexhull
			integer			:: nelim
			real*8				:: eps, fmin, toldiam
		end subroutine riduciconvexhull
		subroutine dealloca_struct(root)
			use mod_type
			implicit none

			type(colonna), pointer	:: root
		end subroutine dealloca_struct
		subroutine suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
			use mod_type
			use mod_suddividi
			implicit none

			type(vertice), pointer	:: currch
			type(colonna),    pointer	:: root
			integer			:: n, nf, nint, cont
			real*8				:: xdir(n), fdir, xdir_unscaled(n)
			real*8				:: maxL,Ltilde,tol
		end subroutine suddividi
		subroutine stampa_intervalli(n,root)
			use mod_type
			implicit none

			type(colonna),   pointer	:: root
			integer			:: n
		end subroutine stampa_intervalli
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8				:: tol
		end subroutine find_colonna
	end interface

	type(colonna),    pointer		:: root, currcol, tempcol
	type(intervallo), target		:: primo
	type(intervallo), pointer		:: curr
	type(vertice),    pointer		:: convexhull, currch, currch1, currchmax
	type(vertice),	  pointer		:: filter
	type(fpunt)				:: ottimo

	external				:: funct

	integer				:: n, iprint
	real*8					:: alfa_stop, tolglob, trigLS
 
	integer				:: nf, ng, nint, num, i, nelim, nconv, k ,ktot, maxiter, numnf, numng
	integer				:: maxnf, maxint, iexit, nminloc, cont, maxcont
	logical				:: halt, direct_puro, trovato, lista_vuota, minric !minric=true fatte min loc per ric
	logical                         	:: vicino
	real*8					:: norma, maxdiam, mindiam, toldiam, basdiam, fglob
	real*8					:: eps, eps_cvx
	real*8					:: lbb(n),ubb(n)
	real*8					:: xbest(n), fbest, ftemp, tmpder, minder, maxL, Ltilde, gg(n), bl(n), bu(n), blsd(n), busd(n)
	real*8					:: xdir(n),  xx(n), ff, fdir, tol, xdir_unscaled(n), fminloc,  xtemp_sc(n)
	real*8					:: xmin(n,100), fmin(100)
	real*8					:: Lr, epsglob, tau
	integer				:: imin, jmin

	integer				:: seed(2), istop
	real					:: rr		

	! parametri ed inizializzazione

	cont        = 0
	maxcont     = 0
	tol         = 1.d-12
	toldiam     = 1.d+1*sqrt(dble(n))/2.d0 
	basdiam     = 0.0d0
	eps         = 1.d-4
	eps_cvx     = 0.d0
	minder      = 0.0d0
	maxL        = 0.0d0
	Ltilde      = 1.d+60
	halt        = .false.
	trovato     = .false.
	lista_vuota = .false.
	memerror    = .false.
	minric      = .false.
	nf = 0
	!questo è il max numero di iterazioni del locale, attenzione è nel modulo!

	lb  = lbb
	ub  = ubb

	allocate(root)
	nullify(root%next)
	nullify(root%pred)

	allocate(root%int)
	call alloca_intervallo(root%int,n)
	root%int%cent    = 0.5d0
	root%int%dimen   = 1.d0
	root%int%maxdim  = 1.d0
	root%int%der     = 0.d0
	root%int%xbars   = 0.5d0
	root%int%lbs     = 0.d0
	root%int%ubs     = 1.d0

	root%int%diam    = norma(n,root%int%dimen)/2.d0
	root%diam        = norma(n,root%int%dimen)/2.d0
	root%int%flagloc = .false.
	root%int%flagdiv = .true.
	root%int%flagcon = .false.
	root%int%flagopt = .false.
	root%int%id      = 1
	nullify(root%int%next)
	nullify(root%int%pred)

	call unscalevars(n,root%int%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
    	ytemp = xtemp

	call funct(n,xtemp,root%int%fint)

	nf          = 1
	nint        = 1

	fdir        = root%int%fint
	xdir        = root%int%cent

	Lr          = 1.d0
	tau         = 0.05d0
	epsglob     = 0.1d0*tau*dble(n)

	write(*,*) 'record del miglior punto iniziale ...'

	call system('del best.txt')

	ng          = 0
	nconv       = 1
	nelim       = 0
	nminloc     = 0

	xdir_unscaled = xbest
	call stampa_ottimo(n,xbest,fbest,nint)

	direct_puro = .true.

	do while (.not.halt) 

		cont    = cont + 1
		maxcont = maxcont + 1

		currcol => root

		if(iprint > 0) then
			write(*,*) 'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo'
			write(*,*) '  id          diam           fint'
		endif
		
		do while (associated(currcol))
			if(currcol%diam < basdiam) then
				!------------------------------------------------
				! Elimina tutta la colonna corrispondente
				!------------------------------------------------
				call elimina_colonna(currcol,num)
				nint = nint - num
				if(.not.associated(currcol%next)) then
					lista_vuota = .true.
					deallocate(currcol)
					nullify(currcol)
					nullify(root)
					exit
				else
					currcol => currcol%next
					deallocate(root)
					root => currcol
					nullify(root%pred)
				endif
			else
				if(.not.associated(currcol%next)) then
					maxdiam = currcol%diam
					! non ci vuole un exit?					
					! exit
				endif
				currcol => currcol%next
			endif
		enddo

!------------------------------------------------------------------------------

		if(tau*maxdiam**2.d0 < epsglob) epsglob = tau*maxdiam**2.d0

		if(associated(root)) then
			nullify(root%pred)
			if((.not.associated(root%int)).and.(.not.associated(root%next))) lista_vuota = .true.
		else
			lista_vuota = .true.
		endif
		if(lista_vuota) then
			halt = .true.
			write(*,*) 'lista vuota'
			exit
		endif

		mindiam = root%diam

		!--------------------------
		! stopping criterion
		!--------------------------

		trovato = (abs(fbest-fglob)/dmax1(1.0d0,abs(fglob)) < tolglob)
		if((nf >= maxnf).or.(nint >= maxint).or.trovato) then
			halt = .true.
			cycle
		endif

		!---------------------------------------
		! ric. intervallo potenzialmente ottimo
		!---------------------------------------

		write(*,*) 'inizio ric. potenzialmente ottimi'
		!call ricintervallo_dx(root,convexhull,nconv,-1,Ltilde,eps,eps_cvx,fdir)
		call ricintervallo   (root,convexhull,nconv,-1,Ltilde,eps_cvx)
		if (memerror) then
			write(*,*) '  fine memoria disponibile'
			halt = .true.
		endif
		write(*,*) '  fine ric. potenzialmente ottimi'

		!write(*,*) 'nconv=',nconv

		!----------------------------------------------
		! riduci il convex hull con il criterio su
		! fmin
		!----------------------------------------------

		if (direct_puro) toldiam = 0.0d0
		
		write(*,*) 'inizio riduzione convex hull'
!	        call riduciconvexhull(convexhull,nelim,eps,toldiam,fdir)
		call riduciconvexhull(convexhull,nelim,eps,toldiam,fbest)
		write(*,*) '  fine riduzione convex hull'

		currch => convexhull
		
		if(iprint > 0) write(*,*) 'nconv = ',nconv,' nelim = ',nelim

		do i = 1,nelim
			currch => currch%next
		enddo

		!-----------------------------------------------------
		! Fai le minimizzazioni locali dai centroidi di
		! intervalli pot. ottimi da cui non sei gia partito
		!-----------------------------------------------------

		if(associated(currch)) then
			!write(*,*) 'currch e'' associato dopo eliminazione'
			minder = 0.d0

			write(*,*) 'inizio min. locali'

	    		currch1   => currch
			do while (associated(currch1))
				if(nf >= maxnf) exit
				if(nf < trigLS*maxnf) exit
				if( (.not.currch1%int%flagloc) ) then
				!if( (.not.currch1%int%flagloc).and.	&
				!    (currch1%int%fint - fbest + epsglob - Lr*currch1%int%diam*2.d0) < 0.d0 ) then

					maxiter = 10000

					call unscalevars(n, currch1%int%cent, ytemp, xbar, lbs, ubs)
					call unscalevars_direct(n,ytemp,xtemp)

					numnf = 0
					globnf = nf

					call sd_box(n,xtemp,fminloc,lb,ub,currch1%int%dimen*(ub-lb)/2.d0,alfa_stop,min(100*n,maxnf-nf),numnf,iprint-1,istop)

					globnf = globnftot

					nf = nf + numnf
					nminloc = nminloc + 1


					if(fminloc < fbest) then
						if (fminloc < fbest-1.d-0) then
							imin  = 0
							jmin  = 0
							write(*,*)'ho riazzerato la lista'
						else
							write(*,*)'ff = ', fminloc, ' fbest =', fbest
						endif
						fbest = fminloc
						xbest = xtemp

					elseif((fminloc-fbest <= 1.d1).and.(imin < 100)) then
						if(maxval(abs(xbest-xtemp)) > 1.d-1) then
							vicino = .false.
							do i = 1,imin
								if(maxval(abs(xmin(1:n,i)-xtemp)) <= 1.d-1) then
									vicino = .true.
									exit
								endif
							enddo
							if(.not.vicino) then
								imin           = imin+1
								xmin(1:n,imin) = xtemp
								fmin(imin)     = fminloc
							endif
						endif 
					endif

					currch1%int%flagloc = .true.

				endif
				currch1 => currch1%next
			enddo

			write(*,*) '  fine min. locali'

		else
			toldiam = 1.d-1*toldiam
			write(*,*) 'toldiam riduz. ',toldiam
		endif


		!-----------------------------------------------------
		! Se non hai esaurito il budget di calcoli di funzione
		! prosegui DIRECT con le suddivisioni
		!-----------------------------------------------------

		if(nf < maxnf) then

			!-----------------------------------------------
			! rimuovo i rimanenti intervalli sul convexhull
			! dalla struttura dati per inserirli in seguito
			! nella posizione corretta
			!-----------------------------------------------
			if(iprint > 0) then
				write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
				write(*,*) '  id          diam           fint'
				write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
				write(24,*) '     diam           fint',fdir
			endif

			currch1 => currch
			do while (associated(currch1))
				if(iprint > 0) then
					write(*,*) currch1%int%id, currch1%int%diam, currch1%int%fint
					write(24,*) currch1%int%diam, currch1%int%fint
				endif

				call find_colonna(root,currch1%int,tol,tempcol)
				if(tempcol%int%id == currch1%int%id) then
					if(associated(currch1%int%next)) then
						tempcol%int => currch1%int%next
						nullify(tempcol%int%pred)
					else
						!write(*,*) 'la colonna si e'' svuotata quindi la elimino'
						!--------------------------------
						! la colonna si e' svuotata
						! quindi la elimino
						!--------------------------------
						if((.not.associated(tempcol%pred)) .and. (.not.associated(tempcol%next))) then
							!--------------------------------
							! E' l'unica colonna
							!--------------------------------
							nullify(tempcol)
							deallocate(root)
							nullify(root)
							!write(*,*) 'la col e'' unica e la elimino',associated(root)
							exit
						elseif(.not.associated(tempcol%pred)) then
							!write(*,*) 'la col e'' la prima ma non unica'
							!--------------------------------
							! E' la prima ma non l'unica
							!--------------------------------
							tempcol => root
							root => root%next
							deallocate(tempcol)
							nullify(root%pred)
						elseif(.not.associated(tempcol%next)) then
							!--------------------------------
							! E' l'ultima ma non l'unica
							!--------------------------------
							!write(*,*) 'la col e'' l''ultima ma non unica',tempcol%int%id,root%int%id
							nullify(tempcol%pred%next)
							deallocate(tempcol)
						else
							!--------------------------------
							! E' in mezzo
							!--------------------------------
							!write(*,*) 'la col e'' in mezzo',tempcol%int%id,root%int%id
							tempcol%pred%next => tempcol%next
							tempcol%next%pred => tempcol%pred
							deallocate(tempcol)
						endif 
					endif
				else
					write(*,*) 'shit shit !! il primo non e'' il primo',tempcol%int%id,currch1%int%id
					stop
				endif
				currch1 => currch1%next
			enddo

			if(iprint > 0) then
				write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
				write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			endif


			if(iprint > 0) write(*,*) 'ciclo principale: dopo esegui'

			if(fdir < fbest) then
				fbest = fdir
				xbest = xdir_unscaled
				call stampa_ottimo(n,xbest,fbest,nint)
			endif

			if(associated(root)) then
				nullify(root%pred)
				if((.not.associated(root%int)).and.(.not.associated(root%next))) lista_vuota = .true.
			else
				lista_vuota = .true.
			endif
			if(lista_vuota) then
				write(*,*) 'lista vuota'
				lista_vuota = .false.
			endif

			write(*,*) 'inizio suddivisioni'

			do while (associated(currch).and. (.not. memerror))
				!----------------------------------------------------
				! suddivisione dell'intervallo potenzialmente ottimo
				!----------------------------------------------------
	!			if(currch%int%flagdiv) call suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
				if(currch%int%flagdiv) call suddividi(currch,root,n,nf,nint,xdir,fbest,xdir_unscaled,maxL,Ltilde,cont,tol)
				if(cont == 0) maxcont = 0
				if (memerror) then
					write(*,*)'fine memoria, esco al momento della suddivisione '
					halt = .true.
					exit
				endif

				currch => currch%next
			enddo

			write(*,*) '  fine suddivisioni'

			if(iprint > 0) write(*,*) 'ciclo principale: dopo suddividi'


		endif  ! if(nf < maxnf) then





		!-------------------------------------------
		! dealloca la lista che memorizza gli
		! intervalli sul convex hull
		!-------------------------------------------

		currch => convexhull
		do while (associated(currch))
			!if(associated(currch%int)) currch%int%flagcon = .false.
			currch => currch%next
			nullify(convexhull%int)
			nullify(convexhull%next)
			deallocate(convexhull)
			convexhull => currch
		enddo

		write(*,*) '----- ',num_el_L

		if(associated(root)) then
			nullify(root%pred)
			if((.not.associated(root%int)).and.(.not.associated(root%next))) lista_vuota = .true.
		else
			lista_vuota = .true.
		endif
		if(lista_vuota) then
			halt = .true.
			write(*,*) 'lista vuota'
			exit
		endif

		if(fdir < fbest) then
			fbest = fdir
			xbest = xdir_unscaled
			call stampa_ottimo(n,xbest,fbest,nint)
!			call aggiorna_struttura(root,Ltilde,fdir,nint)
		endif

		write(*,800) fbest,fdir,nf,nelim,nconv-nelim,mindiam,maxdiam
		write(*,810) maxL*maxval(ubb-lbb),Ltilde,cont,nf-nint,nint

		!call stampa_intervalli(n,root)
	enddo

	if(fdir < fbest) then
		fbest = fdir
		xbest = xdir_unscaled
		call stampa_ottimo(n,xbest,fbest,nint)
	endif

	call dealloca_struct(root)

	return

800		FORMAT('  fmin=',es11.4,'    fDIR=',es11.4,'      nf=',i11,   '  nelim=',i11,' ncnv-nelim=',i8,' diam=',es11.4,' DIAM=',es11.4)
810		FORMAT('  maxL=',es11.4,'  Ltilde=',es11.4,'    cont=',I11,   ' num_el=',I11,    '   nint=',i11)
830		FORMAT(' fbest=',es11.4)
end subroutine nuovoglobale

subroutine stampa_ottimo(n,x,f,nint)
	implicit none
	integer :: n, nint
	real*8  :: x(n), f
	integer :: i

	open(12,file='best.txt',status='unknown',access='append')
	do i = 1,n
		write(12,100) i,x(i)
	enddo
	write(12,110) f,nint
	close(12)

	return

100 format(1x,'x(',i2,') = ',es17.10)
110 format(1x,'fbest = ',es17.10,' nint = ',i12)

end subroutine stampa_ottimo
