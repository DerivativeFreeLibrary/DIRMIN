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
subroutine riduciconvexhull(convexhull,nelim,eps,toldiam,fmin)
	use mod_type
	implicit none

	type(vertice),pointer		:: convexhull
	type(vertice),pointer		:: currch
	integer						:: nelim, nugua
	real*8						:: eps, toldiam, fmin, L
	logical						:: halt
	
	nelim = 0
	nugua = 0
	halt = .false.
	currch => convexhull
	do while (.not.halt)
		!write(*,*) '1'
		if(associated(currch)) then
			if(currch%int%diam < toldiam) then
				nelim = nelim + 1
				currch => currch%next
				cycle
			endif
		else
			halt = .true.
			exit
		endif

		!write(*,*) '2'
		if(associated(currch%next)) then
		!write(*,*) '3'
			if((currch%next%int%diam - currch%int%diam) > 0.d0) then
			  
		!write(*,*) '4'
				L = (currch%next%int%fint - currch%int%fint) / (currch%next%int%diam - currch%int%diam)
				!if( currch%int%fint - L*currch%int%diam >  fmin - eps*(1.d0 + abs(fmin))) then
				!write(*,*) currch%int%diam, toldiam
				!write(*,*) currch%int%fint,L,currch%int%diam,fmin,eps,abs(fmin)
				!write(*,*) currch%int%fint - L*currch%int%diam , fmin-eps*abs(fmin)

				if( currch%int%fint - L*currch%int%diam >  fmin - eps*max(abs(fmin),1.d-6)  ) then
					
					!if(associated(currch%next)) then
						nelim = nelim + 1 + nugua
						nugua = 0
					!endif
					currch => currch%next
					!currch => convexhull
					!convexhull => convexhull%next
					!nullify(currch%int)
					!nullify(currch%next)
					!deallocate(currch)

				else
					halt = .true.
				endif
			else
				nugua = nugua + 1
				currch => currch%next
			endif
		else
			halt = .true.
		endif
	enddo

	!read(*,*)

	nullify(currch)

	return
end subroutine riduciconvexhull

subroutine ricintervallo_dx(root,convexhull,nconv,iprint,Ltilde,eps,eps_cvx,fmin)
	use mod_type
	use mod_mem
	implicit none
 
	type(colonna),pointer		:: root, currcol, tempcol, ultimacol
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo
	real*8				:: Ltilde, stimaL, fh, dh, eps, eps_cvx, maxL, maxLtemp
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam, fmin
	logical				:: halt
	integer				:: nconv, iprint, sv

	!search the column with max diameter
	ultimacol => root
	do while (associated(ultimacol%next))
		ultimacol => ultimacol%next
	enddo

	!ultimacol points to the column with maximum diameter
	!record the first interval (top, right in the CVX HULL)
	primo => ultimacol%int

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	convexhull%int => primo
	primo%flagcon = .true.
   	currch => convexhull
	currcol => root

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '  id          diam           fint'
		write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fint
		write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*) '     diam           fint',fmin
		write(24,*) convexhull%int%diam, convexhull%int%fint
	endif

	nconv = 1

	maxL = -1.d+10
	currcol => ultimacol

	do while (associated(currcol%pred))
		if(maxL < (ultimacol%int%fint - currcol%pred%int%fint)/(ultimacol%diam - currcol%pred%diam)) then
		  maxL = (ultimacol%int%fint - currcol%pred%int%fint)/(ultimacol%diam - currcol%pred%diam)
		endif
		currcol => currcol%pred
	enddo


	halt = .false.
	if (associated(ultimacol%pred)) then
		currcol => ultimacol%pred
	else
		nullify(currch%next)
		if(iprint > 0) then
			write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			write(24,*)
		endif
		return
	endif
	
	do while (associated(currcol%pred))
		tempcol => currcol%pred
		maxLtemp = -1.d+10
		!do while ((maxLtemp <= maxL).and.(associated(tempcol)))
		do while (associated(tempcol))
			!write(*,*) '=====',maxLtemp,(currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam),maxL,tempcol%int%fint
			if(maxLtemp < (currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam)) then
				maxLtemp = (currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam)
			endif
			tempcol => tempcol%pred
		enddo
		!write(*,*) '[[[[[[[[[[[[',currcol%diam,currcol%int%fint,maxLtemp

		if(maxLtemp <= maxL) then

			!write(*,*) currcol%int%fint,fmin,eps,currcol%diam,maxL
			!write(*,*) (currcol%int%fint - fmin + eps*abs(fmin))/currcol%diam, maxL
			!write(*,*)
			if((currcol%int%fint - fmin + eps*abs(fmin))/currcol%diam <= maxL) then
				allocate(currch%next, stat = sv)
				if(sv.ne.0) then
					memerror = .true.
					return
				endif
				currch%next%int => currcol%int
				currch => currch%next
				if(iprint > 0) then
					write(*,*) currch%int%id, currch%int%diam, currch%int%fint
					write(24,*) currch%int%diam, currch%int%fint
				endif
				nconv  = nconv + 1
				maxL = maxLtemp
			else
				exit
			endif

		endif
		currcol => currcol%pred
		
	enddo

	!write(*,*) root%int%fint,fmin,eps,root%diam,maxL
	!write(*,*) (root%int%fint - fmin + eps*abs(fmin))/root%diam, maxL
	!write(*,*)
	if((root%int%fint - fmin + eps*abs(fmin))/root%diam <= maxL) then
		allocate(currch%next, stat = sv)
		if(sv.ne.0) then
			memerror = .true.
			return
		endif
		currch%next%int => root%int
		currch => currch%next
		if(iprint > 0) then
			write(*,*) currch%int%id, currch%int%diam, currch%int%fint
			write(24,*) currch%int%diam, currch%int%fint
		endif
		nconv  = nconv + 1
	endif

	nullify(currch%next)

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*)
	endif

	return
end subroutine ricintervallo_dx

subroutine ricintervallo(root,convexhull,nconv,iprint,Ltilde,eps_cvx)
	use mod_type
	use mod_mem
	implicit none
 
	type(colonna),pointer		:: root, currcol, aux
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo
	real*8				:: Ltilde, stimaL, fh, dh, eps_cvx
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam
	logical				:: halt
	integer				:: nconv, iprint, sv

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	currcol  => root
	aux      => root
	minfunc  = root%int%fint
	maxdiam  = root%diam
	primo    => root%int

	!----------------------------------------
	! search for the interval with max. diam.
	! and min. objective function
	!----------------------------------------
	do while (associated(currcol%next))
		if( ((currcol%next%int%fint < minfunc).or.		&
		    ((currcol%next%diam > maxdiam).and.(currcol%next%int%fint - minfunc <= 1.d-9))).and. &
			(currcol%next%diam > 1.d-8) ) then
			primo   =>  currcol%next%int
			aux     =>  currcol%next	
			maxdiam =  primo%diam
			minfunc =  primo%fint
		endif
		currcol => currcol%next
	enddo

	currcol => aux%next

	!--------------------------------------
	! record the first interval belonging
	! to the convex hull so far identified
	!--------------------------------------
	convexhull%int => primo
	primo%flagcon = .true.
    	currch => convexhull

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '  id          diam           fint'
		write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fint
		write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*) '      diam           fint'
		write(24,*) convexhull%int%diam, convexhull%int%fint
	endif

	nconv = 1

	halt = .false.
	
	do while (.not.halt)
		!-------------------------------------
		! among those intervals in the upper
		! right region, find the one with
		! maximum cosine with vector (1,0)
		!-------------------------------------
		maxcos = -1.d0
		stimaL = 0.d0
		do while (associated(currcol))
			norma = dsqrt((currcol%diam - currch%int%diam)**2.d0+(currcol%int%fint - currch%int%fint)**2.d0)
			coseno = (currcol%diam - currch%int%diam) / norma			
			if(coseno > maxcos) then
				!write(*,*) 'coseno ',coseno,maxcos
				stimaL = (currcol%int%fint - currch%int%fint)/(currcol%diam - currch%int%diam)
				maxcos = coseno
				primo => currcol%int
				aux   => currcol
			endif
			currcol => currcol%next
		enddo
		currcol => aux%next
		if(stimaL > Ltilde) exit
		if(maxcos > 0.d0) then
			allocate(currch%next, stat = sv)
			if(sv.ne.0) then
				memerror = .true.
				return
			endif
			currch%next%int => primo
			currch => currch%next
			nconv  = nconv + 1
			primo%flagcon = .true.
			if (iprint > 0)	then
				write(*,*) currch%int%id, currch%int%diam, currch%int%fint
				write(24,*) currch%int%diam, currch%int%fint
			endif
			if(.false.) then
				fh = primo%fint
				dh = primo%diam
				primo => primo%next
				do while (associated(primo))
					if((primo%fint - fh)/dh <= EPS_CVX) then
						allocate(currch%next, stat = sv)
						if(sv.ne.0) then
							memerror = .true.
							return
						endif
						currch%next%int => primo
						currch => currch%next
						nconv  = nconv + 1
						primo%flagcon = .true.
						if (iprint > 0)	then
							write(*,*) currch%int%id, currch%int%diam, currch%int%fint
							write(24,*) currch%int%diam, currch%int%fint
						endif
						primo => primo%next
					else
						exit
					endif				
				enddo
            endif
		else
			halt = .true.
		endif
	enddo

	nullify(currch%next)
	
	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*)
	endif

	return
end subroutine ricintervallo

subroutine suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
	use mod_type
	use mod_globale
	use mod_suddividi
	use mod_mem
	use mod_box
	implicit none

	type(vertice), pointer		:: currch
	type(intervallo), pointer	:: curr
	type(colonna),    pointer	:: root, temp
	integer						:: n, nf, nint, cont
	real*8						:: xdir(n), fdir, xdir_unscaled(n)

	integer						:: i
	!real*8						:: lbs(n), ubs(n)
	!real*8						:: vetf1(n), vetf2(n)
	!logical					:: mask(n)
	integer						:: numtrue, ind1(1), ind2(1)
	real*8						:: maxL,Ltilde,tol !,ytemp(n)
	logical						:: flder

	interface
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8						:: tol
		end subroutine find_colonna
		subroutine insert_intervallo(curcol,int)
			use mod_type
			implicit none

			type(colonna),    pointer	:: curcol
			type(intervallo), pointer	:: int
		end subroutine insert_intervallo
		subroutine triplica(primo,root,n,ind,f1,f2,nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flag,cont,tol)
			use mod_type
			use mod_box
			use mod_mem
			implicit none

			type(intervallo), pointer	:: primo
			type(colonna),    pointer	:: root
			integer						:: n, ind, nint, cont
			real*8						:: xdir(n), fdir,xdir_unscaled(n)
			real*8						:: f1, f2, maxL, Ltilde, tol
			logical						:: flag
		end subroutine triplica
	end interface
	
	curr => currch%int
	numtrue = 0
	do i = 1,n
		if(curr%maxdim == curr%dimen(i)) then
			ysud = curr%cent 
			ysud(i) = curr%cent(i) + 1.d0*curr%dimen(i)/3.d0
			call unscalevars(n,ysud,ytemp,xbar,lbs,ubs)
			call unscalevars_direct(n,ytemp,xsud)
			call funct(n,xsud,vetf1(i))
		      if((vetf1(i) < globfbest).and.(nf <= globmaxnf)) then
		 	globxbest = xsud
			globfbest = vetf1(i)
			globnftot = nf
			globnf    = nf
		      endif

			ysud(i) = curr%cent(i) - 1.d0*curr%dimen(i)/3.d0
			call unscalevars(n,ysud,ytemp,xbar,lbs,ubs)
			call unscalevars_direct(n,ytemp,xsud)
			call funct(n,xsud,vetf2(i))
		      if((vetf2(i) < globfbest).and.(nf <= globmaxnf)) then
		 	globxbest = xsud
			globfbest = vetf2(i)
			globnftot = nf
			globnf    = nf
		      endif
			mask(i) = .true.
			numtrue = numtrue + 1

			nf = nf+2
		else
			vetf1(i) = 1.d+30
			vetf2(i) = 1.d+30
			mask(i)  = .false.
		endif
	enddo

	curr%der = 0.d0
	flder    = .false.

	do i = 1,numtrue
		ind1 = minloc(vetf1,mask)
		ind2 = minloc(vetf2,mask)
		if(vetf1(ind1(1)) < vetf2(ind2(1))) then
			mask(ind1(1)) = .false.
			!write(*,*) '======= suddividi: chiamo triplica', numtrue
			call triplica(curr,root,n,ind1(1),vetf1(ind1(1)),vetf2(ind1(1)),nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flder,cont,tol)
			if (memerror) then
				write(*,*)'fine memoria disponibile mentre triplico'
				return
			endif
		else
			mask(ind2(1)) = .false.
			!write(*,*) '======= suddividi: chiamo triplica', numtrue
			call triplica(curr,root,n,ind2(1),vetf1(ind2(1)),vetf2(ind2(1)),nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flder,cont,tol)
			if (memerror) then
				write(*,*)'fine memoria disponibile mentre triplico'
				return
			endif
		endif
		nint = nint + 2
	enddo

	!write(*,*) 'controllo il centrale'
	if(curr%fint - Ltilde*curr%diam/2.d0 <= fdir) then
	!if(.true.) then
		curr%flagcon = .false.
		!write(*,*) '-------insert   primo--------',associated(root),curr%id,curr%diam
		!if(associated(root)) write(*,*) '-------insert   primo--------',associated(root),associated(root%pred),associated(root%next),root%diam
		call find_colonna(root,curr,tol,temp)
		call insert_intervallo(temp,curr)
		!nullify(curr)
	else
		!write(*,*) '-------elimino  primo--------',curr%id

		deallocate(curr%cent,curr%dimen)
		nullify(curr%next)
		nullify(curr%pred)
		deallocate(curr)
		!deallocate(currch%int)
		!nullify(currch%int)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif

	return
end subroutine suddividi

subroutine triplica(primo,root,n,ind,f1,f2,nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flag,cont,tol)
	use mod_type
	use mod_box
	use mod_mem
	implicit none

	interface
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8						:: tol
		end subroutine find_colonna
		subroutine insert_intervallo(curcol,int)
			use mod_type
			implicit none

			type(colonna),    pointer	:: curcol
			type(intervallo), pointer	:: int
		end subroutine insert_intervallo
	end interface

	type(intervallo), pointer	:: primo
	type(colonna),    pointer	:: root, temp
	integer						:: n, ind, nint, cont, sv
	real*8						:: xdir(n), fdir, xdir_unscaled(n)
	real*8						:: f1, f2, norma, deltax, maxL, Ltilde, tol
	real*8						:: g(n)
	real*8						:: newder
	logical						:: flag

	type(intervallo), pointer	:: secondo
	type(intervallo), pointer	:: terzo

	allocate(secondo, stat = sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif
	allocate(terzo, stat = sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	call alloca_intervallo(secondo,n)
	if (memerror) then
		write(*,*)'fine memoria disponibile mentre triplico'
		return
	endif
	
	call alloca_intervallo(terzo,n)
	if (memerror) then
		write(*,*)'fine memoria disponibile mentre triplico'
		return
	endif
 
	!allocate(secondo%cent(n),secondo%dimen(n))
	!allocate(terzo%cent(n),  terzo%dimen(n))

	secondo%cent = primo%cent
	terzo%cent   = primo%cent

	secondo%cent(ind) = secondo%cent(ind) + 1.d0*primo%dimen(ind)/3.d0
	terzo%cent(ind)   = terzo%cent(ind)   - 1.d0*primo%dimen(ind)/3.d0

	secondo%dimen = primo%dimen
	terzo%dimen   = primo%dimen

	primo%dimen(ind)   = primo%dimen(ind)/3.d0
	secondo%dimen(ind) = secondo%dimen(ind)/3.d0
	terzo%dimen(ind)   = terzo%dimen(ind)/3.d0

	primo%maxdim = maxval(primo%dimen)
	primo%diam   = norma(n,primo%dimen)/2.d0

	secondo%maxdim = maxval(secondo%dimen)
	secondo%diam   = norma(n,secondo%dimen)/2.d0
	secondo%xbars  = primo%xbars
	secondo%lbs    = primo%lbs
	secondo%ubs    = primo%ubs

	terzo%maxdim = maxval(terzo%dimen)
	terzo%diam   = norma(n,terzo%dimen)/2.d0
    terzo%xbars  = primo%xbars
	terzo%lbs    = primo%lbs
	terzo%ubs    = primo%ubs
	secondo%fint = f1

	if(f1 < fdir) then
		fdir = f1
		xdir = secondo%cent
		call unscalevars(n,xdir,ytemp,xbar,lbs,ubs)
		call unscalevars_direct(n,ytemp,xdir_unscaled)
		cont = 0
	endif

	terzo%fint = f2

	if(f2 < fdir) then
		fdir = f2
		xdir = terzo%cent
		call unscalevars(n,xdir,ytemp,xbar,lbs,ubs)
		call unscalevars_direct(n,ytemp,xdir_unscaled)
		cont = 0
	endif

	secondo%flagopt = .false.
	terzo%flagopt   = .false.
	primo%flagopt   = .false.

	secondo%flagloc = .false.
	terzo%flagloc   = .false.

	secondo%flagdiv = .true.
	terzo%flagdiv   = .true.

	secondo%flagcon = .false.
	terzo%flagcon   = .false.

	secondo%id      = nint+1
	terzo%id        = nint+2

	call unscalevars(n,secondo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = xtemp(ind)

	call unscalevars(n,primo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = deltax - xtemp(ind)
	secondo%der		= abs(f1 - primo%fint)/abs(deltax)
!	call grad(xtemp,n,g)
!	secondo%der = norma(n,g)

	if (maxL < secondo%der	) then
		maxL = secondo%der
		cont = 0
	endif

	call unscalevars(n,terzo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = xtemp(ind)

	call unscalevars(n,primo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = deltax - xtemp(ind)
	terzo%der		= abs(primo%fint - f2)/abs(deltax)

	if (maxL < terzo%der ) then
		maxL = terzo%der
		cont = 0
	endif

	primo%der = max(primo%der,abs(secondo%der),abs(terzo%der))

	if(secondo%fint - Ltilde*secondo%diam/2.d0 <= fdir) then
		!write(*,*) '-------insert  secondo--------',associated(root)
		call find_colonna(root,secondo,tol,temp)
		call insert_intervallo(temp,secondo)
		!write(*,*) '-------insert  secondo--------',associated(root),associated(root%pred),associated(root%next)
	else
		!write(*,*) '-------elimino secondo--------'
		deallocate(secondo%cent,secondo%dimen)
		deallocate(secondo)
		nullify(secondo)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif

	if(terzo%fint - Ltilde*terzo%diam/2.d0 <= fdir) then
		!write(*,*) '-------insert   terzo--------',associated(root)
		call find_colonna(root,terzo,tol,temp)
		call insert_intervallo(temp,terzo)
		!write(*,*) '-------insert   terzo--------',associated(root),associated(root%pred),associated(root%next)
	else
		!write(*,*) '-------elimino  terzo--------'
		deallocate(terzo%cent,terzo%dimen)
		deallocate(terzo)
		nullify(terzo)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif
		


	return
end subroutine triplica

subroutine scalevars(n,x,y,xbar,lbs,ubs)
	implicit none

	!-----------------------------------------------------------------
	! dato [lbs , ubs] subset [0 , 1] e xbar in [lbs , ubs]
	! dato x in [lbs , ubs]
	! restituisce y trasformato ricentrando xbar in (lbs+ubs)/2
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n), xbar(n), lbs(n), ubs(n)
	real*8						:: cent(n)
	integer						:: i
	
	cent = (lbs+ubs) / 2.d0

	do i = 1,n
		if(x(i) <= xbar(i)) then
			y(i) = ( (x(i)-lbs(i)) / (xbar(i) - lbs(i)) ) * ( cent(i) -lbs(i) ) + lbs(i)
		else
			y(i) = ( (x(i) - xbar(i)) / (ubs(i) - xbar(i)) ) * ( ubs(i) - cent(i) ) + cent(i)
		endif
	enddo

	return

end subroutine scalevars

subroutine unscalevars(n,y,x,xbar,lbs,ubs)
	implicit none

	!-----------------------------------------------------------------
	! dato [lbs , ubs] subset [0 , 1] e xbar in [lbs , ubs]
	! dato y in [lbs , ubs] trasformato
	! restituisce x in modo che (lbs+ubs)/2 va in xbar
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n), xbar(n), lbs(n), ubs(n)
	real*8						:: cent(n)
	integer						:: i
	
	cent = (lbs+ubs) / 2.d0

	do i = 1,n
		if(y(i) <= cent(i)) then
			x(i) = ( (y(i)-lbs(i)) / (cent(i) - lbs(i)) ) * ( xbar(i) -lbs(i) ) + lbs(i)
		else
			x(i) = ( (y(i) - cent(i)) / (ubs(i) - cent(i)) ) * ( ubs(i) - xbar(i) ) + xbar(i)
		endif
	enddo

	return

end subroutine unscalevars

subroutine scalevars_direct(n,x,y)
	use mod_box
	implicit none

	!-----------------------------------------------------------------
	! dato [lb , ub] e x in [lb , ub]
	! restituisce y in [0 , 1]
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n)
	integer						:: i

	y = (x-lb)/(ub-lb)

	return

end subroutine scalevars_direct

subroutine unscalevars_direct(n,y,x)
	use mod_box
	implicit none

	!-----------------------------------------------------------------
	! dato [lb , ub] e y in [0 , 1]
	! restituisce x in [lb , ub]
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n)
	integer						:: i

	x = lb + y*(ub-lb)
	
	return
end subroutine unscalevars_direct

real*8 function norma(n,x)
	implicit none
	
	integer						:: n
	real*8						:: x(n)
	
	norma = dsqrt(dot_product(x,x))

	return 
end function norma

subroutine gradnum(n,x,g)
	implicit none
	integer			:: n
	real*8			:: x(n),g(n),z(n)
	real*8			:: eps, fp, fm
	integer			:: i
	
	eps = 1.d-6
	z = x
	do i = 1,n
		z(i) = x(i) + eps
		call funct(n,z,fp)
		z(i) = x(i) - eps
		call funct(n,z,fm)
		z(i) = x(i)
		g(i) = (fp-fm)/(2.d0*eps)
	enddo

	return
end subroutine gradnum

subroutine SetParam(file_input)

use Test_Function
implicit none

integer,parameter				:: maxn=1000

integer							:: r, input_iniz, in, i, j, l, t, rr, icheck, k1, kappa
double precision				:: x(dim), db 
double precision, parameter		::pi = 3.14159265358979323846

integer random

character ( len = 64)			::file_input
character ( len = 100)			::line
character (len = 1)				::char
character (len = 2)				::char2
 
input_iniz = 1

open(unit = input_iniz, file = file_input, status = 'old', iostat = icheck)

read(input_iniz,*)n
write(*,*)n

if ((n > maxn).or.(n<1)) then

	write(*,*)'error: the value for n must be in the interval 1-',maxn
	stop

endif

read(input_iniz,*)l2
write(*,*)l2

db=2.0d0**(n+1)-1.0d0

if ((l2 > (db)).or.(l2<1)) then

	write(*,*)'error: the value for l2 must be in the interval 1-',db
	stop

endif

read(input_iniz,*)l3
write(*,*)l3
if ((l3 > floor(dsqrt(dble(n)))).or.(l3<1)) then

	write(*,*)'error: the value for l3 must be in the interval 1-',floor(dsqrt(dble(n)))
	stop

endif

read(input_iniz,*)h
write(*,*)h

if ((h > 30.0d0).or.(h<10.0d0)) then

	write(*,*)'error: the value for h must be in the interval 1-30'
	stop

endif

read(input_iniz,*) random

if ((random /= 0).and.(random /= 1)) then

	write(*,*)'error: the value for random must be 0 or 1'
	stop

endif

allocate(kk(n),aa1(n))

if (random == 0) then

	read(input_iniz,*) k

	if ((k<10).or.(k>20)) then

		write(*,*)'error: the value for k must be in 10-20'
		stop

	endif

	kk=k

else
	
	k=0
	read(input_iniz,*) (kk(i),i=1,n)
	do i=1,n

		k=k+kk(i)

	enddo
	k=k/dble(n)

	!write(*,*)kk

endif


allocate(m(floor(log(dble(l2))/log(2.0d0))+1))

r = l2
kappa = 0
t=0

do while (r>0)

	if (mod(r,2)==1) then

		kappa = kappa+1
		m(kappa)=t
				
	endif
	t=t+1
	r=r/2

enddo

write(*,*)'n= ', n, 'k= ', kappa, 'L3= ',l3 

dim = n+kappa+l3-2
dim_m = kappa
r = m(dim_m)
write(*,*)'dim= ', dim, 'dim_m= ', dim_m, 'r= ',r 
!pause

!write(*,*)'floor(log(dble(l2))/log(2.0d0))+1', floor(log(dble(l2))/log(2.0d0))+1, 'dim_m ', dim_m


allocate(gradiente(dim), go_point(dim), o1(n), do1(n))
allocate(p(l3,n))
allocate(alpha0(l3,r))
allocate(alpha2(l3,r))
allocate(alpha3(l3,r))
allocate(beta0(l3,r))
allocate(beta2(l3,r))
allocate(beta3(l3,r))

if (dim_m>2) then
	
	allocate(add(dim_m-2),dadd(dim_m-2))

endif

allocate(w(n), dtmp(n+dim_m-1), dw(n,n), f(dim_m), df(dim_m,n+dim_m-1), dg(n+dim_m-1))

if (l3 >2) then

	allocate(add_g(l3-2), dadd_g(l3-2))

endif

allocate(u(l3),du(l3,dim))

read(input_iniz,*)c1

write(*,*)'c1', c1

if ((c1<-3.5d0).or.(c1>-2.0d0)) then
	
		write(*,*)'error: the value for c1 must be in [-3.5,-2]'
		stop

endif

read(input_iniz,*)c2

if ((c2<2.d0).or.(c2>3.5d0)) then
	
		write(*,*)'error: the value for c2 must be in [2,3.5]'
		stop

endif

write(*,*)'c2', c2

do i=1,l3

	read(input_iniz,*) (p(i,j),j=1,n)

enddo


do i=1,n

	read(input_iniz,*) (dw(i,j),j=1,n)

enddo

!do i=1,n

!	write(*,*)dw(i,:)

!enddo


do i=1,n

	aa1(i) = (2.0d0*pi/(c2-c1))*ceiling(kk(i)*(c2-c1)/10.0d0)

enddo

a2 = (2.0d0*pi/5.0d0)*ceiling(k*0.5d0)


do i=1,l3

	do j=1,r

		alpha0(i,j) = dble(p(i,j))
		alpha2(i,j) = -3.0d0*(p(i,j)-5.0d0)/((-c1)*(-c1))
		alpha3(i,j) = 2.0d0*(p(i,j)-5.0d0)/((-c1)*(-c1)*(-c1))
		beta0(i,j) = 1.0d0-dble(p(i,j))
		beta2(i,j) = 3.0d0*(dble(p(i,j))+4.0d0)/((-c2)*(-c2))
		beta3(i,j) = -2.0d0*(dble(p(i,j))+4.0d0)/((-c2)*(-c2)*(-c2))
		!write(*,*)'alpha0', i, j,'=',alpha0(i,j)
		!write(*,*)'alpha2', i, j,'=',alpha2(i,j)
		!write(*,*)'alpha3', i, j,'=',alpha3(i,j)
		!write(*,*)'beta0', i, j,'=',beta0(i,j)
		!write(*,*)'beta2', i, j,'=',beta2(i,j)
		!write(*,*)'beta3', i, j,'=',beta3(i,j)
		!pause

	enddo

enddo

do j=1,m(dim_m)

	if (p(l3,j)==1) then
		
		go_point(j) = c2

	else

		go_point(j) = c1
		
	endif

enddo

do j= m(dim_m)+1,n

	if (p(l3,j)==1) then
		
		go_point(j) = c1

	else

		go_point(j) = c2
		
	endif
	

enddo


do j=n+1,dim

	go_point(j)=2.5d0

enddo

write(*,*)'go_point = ', go_point
return

end subroutine SetParam

subroutine rot_matrix(n,c)

use m_testCE
implicit none

integer ::i,j,n, seed(12), info
double precision:: A(n,n),rr, tau(n), P(n,n), Q1(n,n), &
 Q(n,n), Q2(n,n), v(n,1), Id(n,n), u(n), c, const
double precision:: work(2*n)
real r

seed(1) = 1967408543         
seed(2) = 499
call random_seed(put=seed)

Id=0.0d0
do i=1,n

	Id(i,i) = 1.0d0

enddo

do i=1,n
	do j=1,n
		call random_number(r)
		rr=dble(r)
		A(i,j)=-1.0d0+rr*2.0d0
	enddo
enddo


!CALL DGEQRF(n,n,A,n, tau,work,2*n,info )


P=Id
do i=1,n
	v = 0.0d0
	v(i,1) = 1.0d0
	v(i+1:n,1) = A(i+1:n,i)
	Q1 = P
	Q2 = Id-tau(i)*matmul(v,transpose(v)) 
	P = MATMUL(Q1,Q2)
	
enddo




do i=1,n
	do j=1,n
		call random_number(r)
		rr=dble(r)
		A(i,j)=-1.0d0+rr*2.0d0
	enddo
enddo


!CALL DGEQRF(n,n,A,n, tau,work,2*n,info )

Q = Id
do i=1,n
	v = 0.0d0
	v(i,1) = 1.0d0
	v(i+1:n,1) = A(i+1:n,i)
	Q1 = Q
	Q2 = Id-tau(i)*matmul(v,transpose(v)) 
	Q = MATMUL(Q1,Q2)

enddo


do i = 1,n
	call random_number(r)
	rr=dble(r)
	v(i,1) = rr

enddo
const = (maxval(v(:,1))-minval(v(:,1)))

do i=1,n
	u(i)=c*( (v(i,1)-minval(v(:,1)) )/const)
enddo



do i = 1,n

	Id(i,i) = u(i)

enddo

M_griew=matmul(matmul(P,Id),Q);

!do i=1,n

!	write(*,*)M_griew(i,:)
!	pause


!enddo


return


end subroutine rot_matrix

subroutine read_M(n,file_input)

use m_testCE
implicit none

integer::i, j, icheck,n
character ( len = 64)			::file_input


OPEN(unit=10, FILE=file_input, STATUS='OLD', IOSTAT=ICHECK)
!write(*,*)'n = ',n
do i=1,n

	READ(10,*) (M_griew(i,j),j=1,n)
!	write(*,*)M_griew(i,:)
!	pause

enddo

close(10)

!WRITE(*,*) M_griew
!pause

end subroutine read_M

