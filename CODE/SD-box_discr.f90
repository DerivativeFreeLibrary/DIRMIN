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
      subroutine sd_box(n,x,f,bl,bu,alfa_iniz,alfa_stop,nf_max,nf,iprint,istop)
      use mod_globale
      implicit none
	  logical :: cambio_eps
      integer :: n,i,j,i_corr,nf,ni,nf_max,index_int(n)
      integer :: num_fal,istop
      integer :: iprint,i_corr_fall
	  integer :: flag_fail(n)

      real*8 :: x(n),z(n),d(n),alfa_iniz(n)
      real*8 :: alfa_d(n),alfa,alfa_max, alfa_d_old
      real*8 :: f,fz , eta
	  real*8 :: bl(n),bu(n),alfa_stop,maxeps,scale_int(n) 
	  logical:: discr_change

!     vettore dei valori di f sui punti di un simplesso n+1 dim.

      real*8 :: fstop(n+1)
	  index_int = 0
	  scale_int = 0.0d0

!     num_fal rappresenta il numero di fallimenti consecutivi

!     i_corr rappresenta l'indice della direzione corrente


!     inizializzazione

      nf = 0
      ni = 0
	  discr_change = .false. 

	  eta = 1.d-6

      flag_fail=0

	  num_fal=0

      istop = 0

      fstop=0.d0


      if(iprint.ge.1)  then
         do i=1,n
              write(*,*) ' bl(',i,')=',bl(i)
              write(*,*) ' bu(',i,')=',bu(i)
			  write(*,*) ' x(',i,')=',x(i)
           
		 end do
      endif
  

!     ---- scelta iniziale dei passi lungo le direzioni --------

      do i=1,n

        if(index_int(i).eq.0) then
        
           !alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i))))
		   !alfa_d(i)=1.d0
		   alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,alfa_iniz(i)))

           if(iprint.ge.1) then
              write(*,*) ' alfainiz(',i,')=',alfa_d(i)
              write(1,*) ' alfainiz(',i,')=',alfa_d(i)
           endif
		else

		   alfa_d(i)=dmax1(scale_int(i),dmin1(2.d0*scale_int(i),dabs(x(i))))
      
           if(iprint.ge.1) then
              write(*,*) ' alfainiz(',i,')=',alfa_d(i)
              write(1,*) ' alfainiz(',i,')=',alfa_d(i)
           endif
		end if
      end do
!     -----------------------------------------------------------

!     ---- scelta iniziale delle direzioni ----------------------

      do i=1,n      
        d(i)=1.d0 
      end do
!     -----------------------------------------------------------  

      x = max(bl,min(bu,x))	  
     
      call funct(n,x,f)
      if((f < globfbest).and.(globnf+nf <= globmaxnf)) then
 	globxbest = x
	globfbest = f
	globnftot = globnf + nf
      endif

!	  nf=nf+1

	  i_corr=1

      fstop(i_corr)=f

      do i=1,n
	    z(i)=x(i)
      end do

      if(iprint.ge.2) then
        write(*,*) ' ----------------------------------'
        write(1,*) ' ----------------------------------'
        write(*,*) ' finiz =',f
        write(1,*) ' finiz =',f
        do i=1,n
          write(*,*) ' xiniz(',i,')=',x(i)
          write(1,*) ' xiniz(',i,')=',x(i)
        enddo
      endif

!---------------------------   
!     ciclo principale
!---------------------------

      do 

         if(iprint.ge.1) then
           write(*,*) '----------------------------------------------'
           write(1,*) '----------------------------------------------'
           write(*,100) ni,nf,f,alfa_max
           write(1,100) ni,nf,f,alfa_max
100        format(' ni=',i4,'  nf=',i5,'   f=',d12.5,'   alfamax=',d12.5)
         endif
         if(iprint.ge.2) then
	       do i=1,n
                write(*,*) ' x(',i,')=',x(i)
                write(1,*) ' x(',i,')=',x(i)
            enddo
         endif
!-------------------------------------
!    campionamento lungo asse i_corr
!-------------------------------------
         if(index_int(i_corr).eq.0) then 
 
                call linesearchbox_cont(n,scale_int,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                           alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)

         else
                alfa_d_old=alfa_d(i_corr)
                call linesearchbox_discr(n,eta,index_int,scale_int,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                           alfa_max,i_corr_fall,iprint,bl,bu,ni,nf,discr_change,flag_fail)

         endif

         if(dabs(alfa).ge.1.d-12) then
		    
			flag_fail(i_corr)=0
		               
            x(i_corr) = x(i_corr)+alfa*d(i_corr)
            f=fz
 	        fstop(i_corr)=f
			     
            num_fal=0
            ni=ni+1
      
         else

			flag_fail(i_corr)=1

			if ((index_int(i_corr).eq.1).and.(alfa_d_old.gt.scale_int(i_corr)))  flag_fail(i_corr)=0

	        if(i_corr_fall.lt.2) then 

		      fstop(i_corr)=fz         

              num_fal=num_fal+1
              ni=ni+1

	        endif

	     end if

		 z(i_corr) = x(i_corr)

        !write(*,*) 'icorr=',i_corr
         if(i_corr.lt.n) then
            i_corr=i_corr+1
         else
		    if(.not.discr_change) then
				do i = 1,n
					if((index_int(i).eq.1).and.(alfa_d(i) > 1)) then
						discr_change = .true.
						exit
					endif
				enddo
				if(.not.discr_change) then
					eta = eta/2.d0
				endif
			endif
            i_corr=1
	   	    discr_change = .false. 
         end if 

         call stop(n,index_int,scale_int,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max,flag_fail)

         if (istop.ge.1) exit

         !------------------------------------------------
         ! Aggiornamento parametro di smoothing eps
         !------------------------------------------------

!	  call funct(n,x,f)
	     
!         if(max(0.d0,maxval(constr)).gt.0.d0) then


      enddo
      return
    


      end
        


!     #######################################################

      subroutine stop(n,index_int,scale_int, alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max, flag_fail)

      implicit none
      
      integer :: n,istop,i,nf,ni,nf_max
	  integer :: index_int(n), flag_fail(n)

      real*8 :: alfa_d(n),alfa_max,fstop(n+1),ffstop,ffm,f,alfa_stop
	  real*8 :: scale_int(n) 

	  logical :: test

      istop=0

      alfa_max=0.0d0


      do i=1,n				
	    if (index_int(i) == 0) then
          if(alfa_d(i) > alfa_max) then
            alfa_max=alfa_d(i)
          endif
		endif
      enddo
      !write(*,*) 'alfa_d',alfa_d
      !write(*,*) 'index_int',index_int
    ! write(*,*) 'alfamax = ',alfa_max
!     ffm=1.d+30
!      do i=1,n+1
!        if(fstop(i).lt.ffm) ffm=fstop(i)
!      enddo
      if(ni.ge.(n+1)) then
        ffm=f
        do i=1,n
          ffm=ffm+fstop(i)
        enddo
        ffm=ffm/dfloat((n+1))

 !       ffstop=(f-ffm)*(f-ffm)
 !       do i=1,n
 !          ffstop=ffstop+(fstop(i)-ffm)*(fstop(i)-ffm)
 !       enddo
 
!        ffstop=dsqrt(ffstop/dfloat(n+1))

!        if(ffstop.le.alfa_stop) then
!         istop = 1
!       end if

	  endif

!      write(*,*) ffstop, ' n=',n
!	  write(*,*) fstop
        
      if(alfa_max <= alfa_stop) then
	    test=.true.
		do i=1,n
		
         if (index_int(i).eq.1) then 
        
		  if((alfa_d(i).ne.scale_int(i)).or.(flag_fail(i).ne.1)) then
!		    write(1,*) ' scale_int(',i,')=',scale_int(i),'  flag_fail(',i,')=',flag_fail(i)
!			write(*,*) ' scale_int(',i,')=',scale_int(i),'  flag_fail(',i,')=',flag_fail(i)
		    test=.false.
	      end if
        
		 end if
	  

		end do
        if (test.eqv..true.) then
		   istop = 1
		end if
        
	  end if
      


      if(nf.gt.nf_max) then
        istop = 2
      end if

      if(ni.gt.nf_max) then
        istop = 3
      end if

      return

      end




!     *********************************************************
!     *         
!     *         linesearch lungo le variabili continue
!     *
!     ********************************************************
           
 
      subroutine linesearchbox_cont(n,scale_int,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                                 alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)
      use mod_globale
      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: ni,num_fal
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n),scale_int(n)
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma, gamma_int
      real*8 :: delta,delta1,fpar,fzdelta

!      if(ni.eq.0) gamma=1.d-6*dmax1(1.d-9,dmin1(dabs(f),1.d0)) 
	  
	  gamma=1.d-6      !-6

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0

!     indice della direzione corrente

      j=i_corr

	  if(iprint.ge.1) then
			write(*,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
			write(1,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
	  endif


	  if(dabs(alfa_d(j)).le.1.d-3*dmin1(1.d0,alfa_max)) then
			alfa=0.d0
			if(iprint.ge.1) then
				 write(*,*) '  alfa piccolo'
				 write(1,*) '  alfa piccolo'
				 write(*,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
				 write(1,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
			endif
			return
	  endif
      
!     ciclo su verso della direzione

	  do ielle=1,2

		 if(d(j).gt.0.d0) then

		     if((alfa_d(j)-(bu(j)-x(j))).lt.(-1.d-6)) then                 
   			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
			    alfa=bu(j)-x(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' punto espan. sulla front. *',alfa
					   write(1,*) ' punto espan. sulla front. *',alfa
				endif
			 endif

		  else

			 if((alfa_d(j)-(x(j)-bl(j))).lt.(-1.d-6)) then
			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
				alfa=x(j)-bl(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' punto espan. sulla front. *',alfa
					   write(1,*) ' punto espan. sulla front. *',alfa
				endif
			 endif

		  endif

		  if(dabs(alfa).le.1.d-3*dmin1(1.d0,alfa_max)) then
  
			 d(j)=-d(j)
			 i_corr_fall=i_corr_fall+1
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta per alfa piccolo'
				   write(1,*) ' direzione opposta per alfa piccolo'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   write(1,*) ' j =',j,'    d(j) =',d(j)
				   write(*,*) ' alfa=',alfa,'    alfamax=',alfa_max
				   write(1,*) ' alfa=',alfa,'    alfamax=',alfa_max
			  endif
			  alfa=0.d0
			  cycle

		  endif

		  alfaex=alfa

		  z(j) = x(j)+alfa*d(j)
    
		  call funct(n,z,fz)
		  
		  nf=nf+1

	      if((fz < globfbest).and.(globnf+nf <= globmaxnf)) then
	 	globxbest = z
		globfbest = fz
		globnftot = globnf + nf
	      endif

		  if(iprint.ge.1) then
				write(*,*) ' fz =',fz,'   alfa =',alfa
				write(1,*) ' fz =',fz,'   alfa =',alfa
		  endif
		  if(iprint.ge.2) then
			  do i=1,n
				  write(*,*) ' z(',i,')=',z(i)
				  write(1,*) ' z(',i,')=',z(i)
			  enddo
		  endif

		  fpar= f-gamma*alfa*alfa

!         test sulla direzione

		  if(fz.lt.fpar) then

!            espansione

			 do

!				 if((ifront.eq.1).or.(num_fal.gt.n-1)) then
				  if((ifront.eq.1)) then

			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
				         write(1,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
			         endif
				     alfa_d(j)=delta*alfa

				     return

				 end if

				 if(d(j).gt.0.d0) then
							
					 if((alfa/delta1-(bu(j)-x(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=bu(j)-x(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
							write(1,*) ' punto espan. sulla front.'
						 endif
					 end if

				 else

					 if((alfa/delta1-(x(j)-bl(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=x(j)-bl(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
							write(1,*) ' punto espan. sulla front.'
						 endif
					 end if

				 endif
						 
				 z(j) = x(j)+alfaex*d(j) 
				   
				call funct(n,z,fzdelta)
				
				 nf=nf+1

			      if((fzdelta < globfbest).and.(globnf+nf <= globmaxnf)) then
			 	globxbest = z
				globfbest = fzdelta
				globnftot = globnf + nf
			      endif

				 if(iprint.ge.1) then
					  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
					  write(1,*) ' fzex=',fzdelta,'  alfaex=',alfaex
				 endif
				 if(iprint.ge.2) then
					  do i=1,n
						 write(*,*) ' z(',i,')=',z(i)
						 write(1,*) ' z(',i,')=',z(i)
					  enddo
				 endif

				 fpar= f-gamma*alfaex*alfaex

				 if(fzdelta.lt.fpar) then

					 fz=fzdelta
					 alfa=alfaex

				 else               
					 alfa_d(j)=delta*alfa
			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto fz =',fz,'   alfa =',alfa
				         write(1,*) ' accetta punto fz =',fz,'   alfa =',alfa
			         endif
					 return
				 end if

		     enddo ! ciclo estrapolazione

		  else   ! direzione opposta     

			 d(j)=-d(j)
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta'
				   write(1,*) ' direzione opposta'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   write(1,*) ' j =',j,'    d(j) =',d(j)
			 endif

		  endif       ! test sulla direzione
			  
	  enddo       ! ciclo sul verso della direzione

	  if(i_corr_fall.eq.2) then
			 alfa_d(j)=alfa_d(j)
	  else
			 alfa_d(j)=delta*alfa_d(j)
	  end if

	  alfa=0.d0

	  if(iprint.ge.1) then
			write(*,*) ' fallimento direzione'
			write(1,*) ' fallimento direzione'
	  endif

	  return      
	  
      end

!     *********************************************************
!     *         
!     *         linesearch lungo le variabili discrete
!     *
!     ********************************************************
           
 
      subroutine linesearchbox_discr(n,eta,index_int,scale_int,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                                 alfa_max,i_corr_fall,iprint,bl,bu,ni,nf,discr_change,flag_fail)
      use mod_globale
      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: ni,num_fal
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
	  integer :: index_int(n),flag_fail(n)
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n),scale_int(n)
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma, gamma_int,eta
      real*8 :: delta,delta1,fpar,fzdelta
	  logical:: discr_change, test

!      if(ni.eq.0) gamma=1.d-6*dmax1(1.d-9,dmin1(dabs(f),1.d0)) 
	  
      gamma_int=1.d-0 

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0

!     indice della direzione corrente

      j=i_corr


	  if(iprint.ge.1) then
			   write(*,*) 'variabile discreta  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
			   write(1,*) 'variabile discreta  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
	  endif

	  test=.true.

      if(alfa_d(i_corr).eq.scale_int(i_corr)) then  
        
		  do i=1,n
		  
		    if((flag_fail(i).eq.0)) then

		      test=.false.
			  exit

	        end if

		  enddo

          if(test) then

		     alfa=0.d0

		     if(iprint.ge.1) then
			    write(*,*) ' direzione gia analizzata'
			    write(1,*) ' direzione gia analizzata'
		     endif

             return
          endif
                  
	   end if
      
!      ciclo per il verso della direzione	 
	 
	   do ielle=1,2

		   if(d(j).gt.0.d0) then

				if( ( ( bu(j)-x(j)-alfa_d(j) ) ).lt.0.d0 ) then  
!				   alfa=floor(bu(j)-x(j))
				   alfa=      bu(j)-x(j)
				   ifront=1
				   if (alfa.eq.0.d0) then 
				      d(j)=-d(j)
					  ifront=0
				      cycle
				   endif
                else
				   alfa=alfa_d(j)
				end if		

		   else

				if( ((x(j)-alfa_d(j)-bl(j))).lt.0.0d0 ) then
!				   alfa=floor(x(j)-bl(j))
				   alfa=      x(j)-bl(j)
				   if(iprint .gt. 1) then
					   write(*,*) 'alfa =',alfa
				   endif
				   ifront=1
				   if (alfa.eq.0.d0) then
				      d(j)=-d(j)
					  ifront=0
				      cycle
				   endif
				 else
				   alfa=alfa_d(j)
				endif

		   endif

           alfaex=alfa

		   z(j) = x(j)+alfa*d(j)
    
		   call funct(n,z,fz)
		   
		   nf=nf+1

	      if((fz < globfbest).and.(globnf+nf <= globmaxnf)) then
	 	globxbest = z
		globfbest = fz
		globnftot = globnf + nf
	      endif

		   if(iprint.ge.1) then
				write(*,*) ' fz =',fz,'   alfa =',alfa
				write(1,*) ' fz =',fz,'   alfa =',alfa
		   endif
		   if(iprint.ge.2) then
			   do i=1,n
				  write(*,*) ' z(',i,')=',z(i)
				  write(1,*) ' z(',i,')=',z(i)
			   enddo
		   endif

		   fpar= f-gamma_int*eta

!          test sulla direzione

		   if(fz.lt.fpar) then
			  
			  discr_change = .true.

!             espansione

			  do 
                  if(ifront.eq.1) then 

                     if(iprint.ge.1) then
				              write(*,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
				              write(1,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
			         endif

				     return

				  endif

				  if(d(j).gt.0.d0) then
							
					 if((bu(j)-x(j)-2.0d0*alfa ).lt.(0.0d0)) then

!					    alfaex=floor(bu(j)-x(j))
					    alfaex=      bu(j)-x(j)
				        ifront=1

				        if (alfaex.le.alfa) then
						 
						   alfa_d(j)=max(scale_int(j),max(dble(floor((alfa/2.0d0)/scale_int(j)+0.5d0)),1.d0)*scale_int(j))

			               if(iprint.ge.1) then
				              write(*,*) ' accetta punto quasi frontiera fz =',fz,'   alfa =',alfa
				              write(1,*) ' accetta punto quasi frontiera fz =',fz,'   alfa =',alfa
			               endif

						   return

					    endif
						 
					  else

					     alfaex=alfa*2.0d0						
					
					  end if

				   else

					  if(( x(j)-2.0d0*alfa-bl(j) ).lt.(0.0d0)) then

!					      alfaex=floor(x(j)-bl(j))
					      alfaex=      x(j)-bl(j)
						  ifront=1

						  if (alfaex.le.alfa) then
						 
						     alfa_d(j)=max(scale_int(j),max(dble(floor((alfa/2.0d0)/scale_int(j)+0.5d0)),1.d0)*scale_int(j))

			                 if(iprint.ge.1) then
				               write(*,*) ' accetta punto quasi frontiera fz =',fz,'   alfa =',alfa
				               write(1,*) ' accetta punto quasi frontiera fz =',fz,'   alfa =',alfa
			                 endif

						     return

						  endif
						 
					  else

					      alfaex=alfa*2.0d0						
					
					  end if

				   endif
						 
				   z(j) = x(j)+alfaex*d(j) 
				   
     
				  call funct(n,z,fzdelta)
				
				   nf=nf+1

			      if((fzdelta < globfbest).and.(globnf+nf <= globmaxnf)) then
			 	globxbest = z
				globfbest = fzdelta
				globnftot = globnf + nf
			      endif

				   if(iprint.ge.1) then
					  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
					  write(1,*) ' fzex=',fzdelta,'  alfaex=',alfaex
				   endif
				   if(iprint.ge.2) then
					  do i=1,n
						 write(*,*) ' z(',i,')=',z(i)
						 write(1,*) ' z(',i,')=',z(i)
					  enddo
				   endif

				   fpar= f-gamma_int*eta

				   if(fzdelta.lt.fpar) then

					  fz=fzdelta
					  alfa=alfaex

				   else               
					   alfa_d(j)=max(scale_int(j),max(dble(floor((alfa/2.0d0)/scale_int(j)+0.5d0)),1.d0)*scale_int(j))
			           if(iprint.ge.1) then
				         write(*,*) ' accetta punto  fz =',fz,'   alfa =',alfa
				         write(1,*) ' accetta punto  fz =',fz,'   alfa =',alfa
			           endif

					  return
				   end if

				enddo

			 else 

				d(j)=-d(j)
				ifront=0

				if(iprint.ge.1) then
				   write(*,*) ' direzione opposta'
				   write(1,*) ' direzione opposta'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   write(1,*) ' j =',j,'    d(j) =',d(j)
				endif

			 endif
			  
		  enddo

		  alfa_d(j)=max(scale_int(j),max(dble(floor((alfa/2.0d0)/scale_int(j)+0.5d0)),1.d0)*scale_int(j))

		  alfa=0.d0
		  if(iprint.ge.1) then
			  write(*,*) ' fallimento direzione'
			  write(1,*) ' fallimento direzione'
		  endif

		  return
	  	     
      end


