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
program main
	call mainNuovoGLOB()
end program main
subroutine mainNuovoGLOB()

	use mod_box
	use mod_suddividi
	use mod_mem
	use mod_globale
	use m_testCE
	
	implicit none
	include 'parametri.fi'

	integer		:: n
	character*40		:: nomefun
	character*10		:: fun_name
	character*11		:: tempstr
	character ( len = 64)	:: file_input
	real*8			:: lbb(nn), ubb(nn), lb0(nn), ub0(nn), fglob, xaux(nn), yaux(nn)
	real*8			:: xbest(nn), fbest, ftemp, mindist, maxdist
	real*8			:: mindiam,maxdiam,maxL, c
	real			:: rr,time_begin, time_end, timetot, f_max, f_min
	integer		:: imin, jmin
	integer		:: seed(2), iprint, i, nf, numnf, numng, k, ktot, j, ng
	integer		:: maxint, nint, maxnf
	integer*8		:: ninttot, nftot, nminloc
	logical		:: trovato
	real*8			:: alfa_stop, tolglob, trigLS
	real*8, allocatable	:: xmin(:,:), fmin(:)
	real, allocatable	:: lbins(:), ubins(:), xr(:)
	external		:: funct, grad
 
	memerror = .false.

	n = nn

	write(*,*) 'n=',n

	!----------------------------------------
	! Allocate storage arrays
	!----------------------------------------
	allocate(lb(n),ub(n),xtemp(n),ytemp(n),xbar(n),lbs(n),ubs(n))
	allocate(vetf1(n),vetf2(n),xsud(n),ysud(n),mask(n))
	allocate(xmin(n,100),fmin(100))
	allocate(lbins(n),ubins(n),xr(n))
	allocate(globxbest(n))

	!----------------------------------------
	! Set bounds, nomefun and other problem
	! dependent parameters
	!----------------------------------------
	open(2,file='nomefun',status='OLD')
	read(2,100) nomefun
	close(2)
 	write(*,*) nomefun
	if (nomefun==' Griewrot2') then
		allocate(M_griew(n,n))
		if (n ==2) then
			file_input = 'matrixM_2.txt'
			call read_M(n,file_input)
		endif
		if (n ==10) then
			file_input = 'matrixM_10.txt'
			call read_M(n,file_input)
		endif
		if (n ==30) then
			file_input = 'matrixM_30.txt'
			call read_M(n,file_input)
		endif
		if (n ==50) then
			file_input = 'matrixM_50.txt'
			call read_M(n,file_input)

		endif
		if ((n.ne.10).and.(n.ne. 2).and.(n.ne.30).and.(n.ne.50)) then
			c = 3.0d0
			write(*,*) c
			call rot_matrix(n,c)
			do i=1,n
				do j=1,n
					call random_number(rr)
					M_griew(i,j)=M_griew(i,j)*(1.0d0+0.3d0*dble(rr))
				enddo
			enddo
		endif
	endif
	if (nomefun ==' ShiftRotRastr') then
		allocate(M_griew(n,n))
		if (n ==2) then
			file_input = 'matrixM_Rastr_2.txt'
			call read_M(n,file_input)
		endif
		if (n ==10) then
			file_input = 'matrixM_Rastr_10.txt'
			call read_M(n,file_input)
		endif
		if (n ==30) then
			file_input = 'matrixM_Rastr_30.txt'
			call read_M(n,file_input)
		endif
		if (n ==50) then
			file_input = 'matrixM_Rastr_50.txt'
			call read_M(n,file_input)
		endif

	endif
	open(2,file='PROBDAT.d',status='OLD')
	read(2,*) fglob

	lbs = 0.d0
	ubs = 1.d0
	do i = 1,n
		read(2,*) lb0(i), ub0(i)
	lbb(i)=lb0(i)
		ubb(i)=ub0(i)
	enddo
	close(2)

	xbest    = (ubb+lbb)/2.d0
	xbar     = (ubs+lbs)/2.d0
	fattore  = 1.0d0
	ampiezza = 1.0d0

	call funct(n,xbest,fbest)
	write(*,*) 'fbest = ',fbest
	write(*,*) 'xbest = ',xbest
	write(*,*) '  ubb = ',ubb
	write(*,*) '  lbb = ',lbb

	!-----------------------------------------------------------------------------
	! Read customizing parameters from file: DIRMIN_params.txt
	!-----------------------------------------------------------------------------
	!
	! iprint    (MUST BE an integer) = the printing level
	!
	! rmaxnf    (MUST BE an integer) = relative maximum num. of fun. evaluations
	!           the maximum number of function evaluations is:  
	!              maxnf = rmaxnf*N  where N is the problem dimension
	!
	! maxint    (MUST BE an integer) = maximum number of intervals.
	!           if maxint < 0, then it is set to 1000*maxnf
	!
	! tolglob   (MUST BE a  real)    = tolerance in the stopping condition
	!           code stops when:
	!             abs(fbest-fglob)/max(1.0d0,abs(fglob)) < tolglob
	!
	! alfa_stop (MUST BE a  real)    = required accuracy for the LineSearches
	!           if (alfa_stop < 1.e-6) then alfa_stop is set to 1.e-6 
	!
	! trigLS    (MUST BE a  real)    = percentage of maxnf function evaluations
	!           performed BEFORE first Linesearch is started
	!           trigLS must be between 0.0 and 1.0
	! 
	!-----------------------------------------------------------------------------
	open(2,file='DIRMIN_params.txt',status='OLD')
	read(2,*) tempstr, iprint
	read(2,*) tempstr, maxnf
	maxnf = maxnf*n
	read(2,*) tempstr, maxint
	if(maxint < 0) then
		maxint = min(1000000,abs(1000*maxnf))
	endif
	read(2,*) tempstr, tolglob
	read(2,*) tempstr, alfa_stop
	if(alfa_stop < 1.d-6) alfa_stop = 1.d-6
	read(2,*) tempstr, trigLS
	if(trigLS < 0.d0) trigLS = 0.d0
	if(trigLS > 1.d0) trigLS = 1.d0
	close(2)
	
	!iprint    = 0
	!maxnf     = 10*n
	!maxint    = 1000*maxnf !200000000
	!tolglob   = 0.d0 !1.d-4 
	!alfa_stop = 1.d-3
	!trigLS    = 0.1d0

	nftot     = 1
	ninttot   = 0
	nminloc   = 0
	imin      = 0
	jmin      = 0

	globxbest = xbest
	globfbest = fbest
	globnf    = 1
	globnftot = 1
	globmaxnf = maxnf

	call cpu_time(time_begin)

	call nuovoglobale(n,lbb,ubb,xbest,fbest,nf,ng,fglob,maxint,maxnf,iprint,trovato, nint,mindiam,maxdiam,maxL, &
					  xmin,fmin,imin,jmin,alfa_stop,tolglob,trigLS)

	call cpu_time(time_end)
	timetot = time_end - time_begin

	nftot   = nftot + nf
	ninttot = ninttot + nint
	nminloc = nminloc + ng

	write(*,*)
	write(*,*) '---------- sommario risultati -----------'
	write(*,*) '  cpu_time = ',timetot, ' secondi'
	write(*,*) '        nf = ',nftot
	write(*,*) '        ng = ',ng
	write(*,*) '     fbest = ',fbest
	write(*,*) '     xbest = ',xbest
	write(*,*) '-----------------------------------------'

	write(1,*)
	write(1,*) '---------- sommario risultati -----------'
	write(1,*) '  cpu_time = ',timetot, ' secondi'
	write(1,*) '        nf = ',nftot
	write(1,*) '        ng = ',ng
	write(1,*) '     fbest = ',fbest
	write(1,*) '     xbest = ',xbest
	write(1,*) '-----------------------------------------'


	if(((globfbest-fglob)/dmax1(1.0d0,abs(fglob)) .ge. 1.d-4)) then

		!write(2,800) nomefun, n, globfbest, globnftot, nftot, nminloc, ninttot, mindiam,maxdiam ,maxL
		write(2,900) nomefun, n, globfbest, globnftot, maxnf, nftot

	else

		!write(2,801) nomefun, n, globfbest, globnftot, nftot, nminloc, ninttot, mindiam,maxdiam ,maxL
		write(2,901) nomefun, n, globfbest, globnftot, maxnf, nftot

	endif

	deallocate(lb,ub,xtemp,ytemp,xbar,lbs,ubs)
	deallocate(vetf1,vetf2,xsud,ysud,mask)
	deallocate(xmin,fmin)
	deallocate(lbins,ubins,xr)
	deallocate(globxbest)

100 format(a40)

800 FORMAT(a20,' & ',i4,' &  \bf ', es16.8,4(' & ',i11),' & ', es10.2,' & ', es10.2,' & ', es10.2,' \\')
801 FORMAT(a20,' & ',i4,' &      ', es16.8,4(' & ',i11),' & ', es10.2,' & ', es10.2,' & ', es10.2,' \\')
802 FORMAT(a20,' & ',i4,' &   *  ', es16.8,4(' & ',i11),' & ', es10.2,' & ', es10.2,' & ', es10.2,' \\')

900 FORMAT(a20,' & ',i4,' &  \bf ', es16.8,3(' & ',i11),' \\')
901 FORMAT(a20,' & ',i4,' &      ', es16.8,3(' & ',i11),' \\')

end subroutine mainNuovoGLOB

