!***********************************************************************************
!
! REALIZZAZIONE SEQUENZIALE DEL METODO SENZA DERIVATE
! S.Lucidi M. Sciandrone, "A DERIVATIVE-FREE ALGORITHM FOR BOUND CONSTRAINED
! OPTIMIZATION", Computational Optim. and Appl. Vol. 21, pp 119-142 (2002) 
!
!***********************************************************************************
PROGRAM seqsdbox
	USE mpi
	IMPLICIT NONE
	INTEGER:: n, i, k, h
	INTEGER:: maxiter, nf
	LOGICAL:: flag_dir
	INCLUDE 'parametri.f'
	PARAMETER(maxiter=100000)
	DOUBLE PRECISION:: x(n), LB(n), UB(n), d(n), z(n), gamma, theta
	DOUBLE PRECISION:: alpha_tilde(n), f, fz, f_ext, alpha_max, alpha, tol
	!-----------------------------------------------------------------------
	!Parametri per MPI
	INTEGER:: rank, nproc, ierr, tag
	LOGICAL:: stop
	INTEGER:: status(MPI_STATUS_SIZE)
	!-----------------------------------------------------------------------
	! TIMING
	DOUBLE PRECISION:: tstart, tend
	
	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
	call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)

	!INIZIALIZZAZIONI
	nf = 0
	alpha_tilde = 1.d0
	d = 1.d0
	gamma = 1.d-6	
	theta = 0.5d0
	tol = 1.d-6
	LB = -1.d+24
	UB =  1.d+24

	if(rank == 0) then	
	call startp(n,x)
	call funct(n,x,f)
	nf = nf+1
	f_ext = f
	!write(*,*) 'f.val:',f

	! TIMING
	tstart = MPI_Wtime()
	
	! CICLO PRINCIPALE
	do k=1, maxiter
		WRITE(*,10) f, maxval(alpha_tilde),k,nf
		z = x
		if(maxval(alpha_tilde) <= tol) exit
		do i=1,n
			!pause
			flag_dir = .true.
			!--------------
			!STEP 1--> di
			!--------------
			f = f_ext
			if(d(i) > 0.d0) then
				alpha_max = UB(i)-x(i)
			else
				alpha_max = x(i)-LB(i)
			end if
			alpha = min(alpha_tilde(i),alpha_max)
			z(i) = x(i)+alpha*d(i)
			call funct(n,z,fz)
			nf = nf+1
			!write(*,*) 'Step1:',alpha, fz, f-gamma*alpha**2.d0
			if((alpha > 0.d0).and.(fz <= f-gamma*alpha**2.d0)) then
				!--------------
				!STEP 3(1)
				!--------------
				call linesearch(i,n,x,z,d,alpha_max,alpha,gamma,f,fz,f_ext,nf,nproc)
				alpha_tilde(i) = alpha
				flag_dir = .false.
			end if
			!--------------
			!STEP 2--> -di
			!--------------
			if(flag_dir) then
				if(d(i) > 0.d0) then
					alpha_max = x(i)-LB(i)
				else
					alpha_max = UB(i)-x(i)
				end if
				alpha = min(alpha_tilde(i),alpha_max)
				z(i) = x(i)-alpha*d(i)
				call funct(n,z,fz)
				nf = nf+1
				!write(*,*) 'Step2:',alpha, fz, f-gamma*alpha**2.d0
				if((alpha > 0.d0).and.(fz <= f-gamma*alpha**2.d0)) then
					!--------------
					!STEP 3(2)
					!--------------
					d(i) = -d(i)
					call linesearch(i,n,x,z,d,alpha_max,alpha,gamma,f,fz,f_ext,nf,nproc)
					alpha_tilde(i) = alpha
				else
					alpha = 0.d0
					alpha_tilde(i) = theta*min(alpha_tilde(i),alpha_max)
				end if
			end if
			!--------------
			!STEP 4
			!--------------
			x(i) = x(i) + alpha*d(i)
			z(i) = x(i)
			f = f_ext
			!write(*,*) 'x:',x
		end do
		! ORACOLO: TROVA x_new t.c. f(x_new)f_ext < f(x)
		call ORACOLO(n,x,f,z,fz)
		x = z 
		f = fz
	end do   
	
	! TIMING
	tend = MPI_Wtime()

	WRITE(2,20) f,k,nf,tend-tstart
	WRITE(*,*) 'tempo totale di esecuzione : ',tend-tstart
	WRITE(*,10) f, maxval(alpha_tilde),k,nf

	else	!rank /= 0
		do while (.true.)
			CALL MPI_RECV(stop,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,status,ierr)
			if(stop) exit
			CALL MPI_RECV(tag,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierr)
			CALL MPI_RECV(x,n,MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
			CALL funct(n,x,f)
			CALL MPI_SEND(f,1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,ierr)
		enddo
	endif

	20&
	FORMAT(1x,' & ',es11.4,' & ',i6,' & ',i6,' & ',f7.2,' \\\hline')

	10&
	FORMAT('obj.f:'es11.4,' max alpha:',d10.4,' k:',i6,' nf:',i6)
		
END PROGRAM

SUBROUTINE ORACOLO(n,x,f,z,fz)
	IMPLICIT NONE
	INTEGER:: n
	DOUBLE PRECISION:: x(n), z(n), f,fz
	
	z = x
	fz = f
	RETURN  

END SUBROUTINE ORACOLO

SUBROUTINE linesearch(i,n,x,z,d,alpha_max,alpha,gamma,f,fz,f_ext,nf,nproc)
	USE mpi	
	IMPLICIT NONE
	!INCLUDE 'mpif.h'
	INTEGER:: n, i, nf, nproc, ierr, j, maxj, tag
	INTEGER:: status(MPI_STATUS_SIZE)
	DOUBLE PRECISION:: x(n), z(n), f, fz, gamma, d(n),f_ext
	DOUBLE PRECISION:: alpha_max, alpha, delta, alpha_tilde
	DOUBLE PRECISION:: alphain
	DOUBLE PRECISION:: job(nproc-1), res(nproc-1), temp(nproc-1)
	LOGICAL:: lancia_altre, stop, iflag, occupato(nproc-1)

	!INIZIALIZZAZIONI
	delta = 0.5d0
	f_ext = fz
	alphain = alpha
	alpha_tilde = min(alpha_max,alphain/delta)
	tag = 1
	stop = .false.
	occupato = .false.

	do while (.true.)
!		write(*,*) occupato
!		CALL FLUSH()
		do j = 1,nproc-1 
			maxj = j
			job(j) = x(i)+alpha_tilde*d(i)
			!----------------------------------
			!e lo manda al processo j-esimo
			!----------------------------------
			CALL MPI_ISEND(stop,1,MPI_LOGICAL,j,0,MPI_COMM_WORLD,status,ierr)
			CALL MPI_ISEND(tag,1,MPI_INTEGER,j,0,MPI_COMM_WORLD,status,ierr)
			CALL MPI_ISEND(job(j),n,MPI_DOUBLE_PRECISION,j,tag,MPI_COMM_WORLD,status,ierr)
			occupato(j) = .true.
			if(alpha_tilde < alpha_max) then
				alpha = alpha_tilde
				alpha_tilde = min(alpha_max,alpha/delta)
			else
				exit
			endif
		enddo
		do while (.true.)
!			write(*,*) occupato
!			READ(*,*)
!			CALL FLUSH()
			lancia_altre = .false.
			do j = 1,maxj-1
				CALL MPI_IPROBE(j, tag, MPI_COMM_WORLD, iflag, status, ierr)
				if(iflag) then
					CALL MPI_RECV(res(j),1,MPI_DOUBLE_PRECISION,j,tag,MPI_COMM_WORLD,status,ierr)
					occupato(j) = .false.
					if((j==1).AND.(res(j) > f_ext - gamma*alphain**2.d0)) then
						alpha = alphain
						return
					endif
					CALL MPI_IPROBE(j+1, tag, MPI_COMM_WORLD, iflag, status, ierr)
					if(iflag) then
						CALL MPI_RECV(res(j+1),1,MPI_DOUBLE_PRECISION,j+1,tag,MPI_COMM_WORLD,status,ierr)
						occupato(j+1) = .false.
						if((j+1==maxj).AND.(res(j+1) <= f_ext - gamma*(alphain/(delta**(maxj)))**2.d0)) then
							tag = tag + 1
							lancia_altre = .true.
							exit
						elseif((res(j)   <= f_ext - gamma*(alphain/(delta**(j  )))**2.d0).AND. &
									 (res(j+1) >  f_ext - gamma*(alphain/(delta**(j+1)))**2.d0)) then
							alpha = alphain/(delta**j)
							return
						endif
					endif
				endif
			enddo
			if(lancia_altre) then
				do j = 1,maxj
					if(occupato(j)) then
!						write(*,*) j
						CALL MPI_IRECV(temp(j),1,MPI_DOUBLE_PRECISION,j,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
						occupato(j) = .false.
					endif	
				enddo
				alpha_tilde = min(alpha_max,alphain/(delta**(maxj+1)))
				exit
			endif
		enddo
	enddo
	
	!z(i) = x(i)+alpha_tilde*d(i) 
	!call funct(n,z,fz)
	!nf = nf+1
	!do while((alpha_tilde < alpha_max ).and.(fz<=f-gamma*alpha_tilde**2.d0))
	!	f_ext = fz
	!	alpha = alpha_tilde
	!	alpha_tilde = min(alpha_max,alpha/delta)
	!	z(i) = x(i)+alpha_tilde*d(i)
	!	call funct(n,z,fz)
	!	nf = nf+1
	!end do
END SUBROUTINE linesearch

SUBROUTINE delay()
	IMPLICIT NONE
	DOUBLE PRECISION:: tstart, tend	

	return
	
	CALL CPU_TIME(tstart)
	CALL CPU_TIME(tend)
	DO WHILE (tend - tstart < 0.01)
		CALL CPU_TIME(tend)
	END DO
	
	RETURN
END SUBROUTINE delay
