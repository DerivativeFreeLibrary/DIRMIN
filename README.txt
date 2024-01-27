-----------------------------------------------------------
 How to run the derivative-free optimizer DIRMIN for 
 on a collection of bound constrained global optimization 
 problems
-----------------------------------------------------------
 The package provides a Fortran90 version of the DIRMIN algorithm.

0- Gunzip and untar the archive in a folder on your computer by
   issuing in a directory of your choice (ex. curdir) the command

   $> tar -zxvf DIRMIN.tar.gz

1- Edit the parameters file DIRMIN_params.txt by specifying:

	iprint    (MUST BE an integer) = the printing level
	
	rmaxnf    (MUST BE an integer) = relative maximum num. of fun. evaluations
	          the maximum number of function evaluations is:  
	          maxnf = rmaxnf*N  where N is the problem dimension
	
	maxint    (MUST BE an integer) = maximum number of intervals.
	          if maxint < 0, then it is set to 1000*maxnf
	
	tolglob   (MUST BE a  real)    = tolerance in the stopping condition
	          code stops when:
	             abs(fbest-fglob)/max(1.0d0,abs(fglob)) < tolglob
	          where fglob is the known global minimum value and
	          fbest is the best value found so far
	
	alfa_stop (MUST BE a  real)    = required accuracy for the LineSearches
	          if (alfa_stop < 1.e-6) then alfa_stop is set to 1.e-6 
	
	trigLS    (MUST BE a  real)    = percentage of maxnf function evaluations
	          performed BEFORE first Linesearch is started
	          trigLS must be between 0.0 and 1.0

IF YOU WANT TO EXECUTE DIRMIN ON THE PROVIDED COLLECTION OF TEST 
PROBLEMS FOR BOUND CONSTRAINED GLOBAL OPTIMIZATION, THEN PROCEED TO STEP 2.
OTHERWISE, IF YOU WANT TO EXECUTE DIRMIN ON A SINGLE USER-SPECIFIED PROBLEM,
PROCEED TO STEP 3
   
2- At command prompt in curdir execute

   $> ./run_dirmin.sh

   This will run DIRMIN on a collection of 74 bound constrained global 
   optimization problems and output the results in file ris.tex

   The result file ris.tex is a formatted text file having 6 columns.
   In particular:
   column 1 gives the PROBLEM NAME
   column 2 gives the NUMBER OF VARIABLES
   column 3 gives the best value found (f*) within MAXNF f.evals.
   column 4 gives the number of f.evals (nf) required to obtain f*
   column 5 gives the maximum number of f.evals MAXNF
   column 6 gives the actual number of f.evals (anf) performed

   Note that nf <= maxnf <= anf and f* is updated until nf <= maxnf.

3- Edit file curdir/parametri.fi to define the number of variables (nn) 

4- Edit file curdir/PROBDAT.d to define 
   fglob, that is the global minimum value (if known) or an underestimate to it
   lower and upper bounds on the variables

5- Edit file curdir/nomefun to define a name for the problem to be solved 

6- Provide a file myfun.f90 with objective function definition. The file MUST
   contain a subroutine definition as the following example

	SUBROUTINE FUNCT(N,X,F)
		IMPLICIT NONE
		INTEGER          :: N,I
		DOUBLE PRECISION :: X(N),F

		F=0.D0
		DO I = 2,N
			F = F+dble(I)*(2.D0*x(I)**2-x(I-1))**2    
		ENDDO
		F = F+(x(1)-1.D0)**2

		RETURN
	END

2- At command prompt in curdir execute 

     $> make FUNCT=myfun
 
   which will create the executable 'dirmin'

4- execute

     $> ./dirmin

   and you are done!
