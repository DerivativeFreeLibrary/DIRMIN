!*************************************************
!	FUNZIONE ZIRILLI2
!*************************************************
      	SUBROUTINE FUNCT(N,X,F)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N)
	F=(1./4.)*X(1)**4-(1./2.)*X(1)**2+.1*X(1)+(1./2.)*X(2)**2
       	RETURN
	END
!*************************************************
!	GRADIENTE
!*************************************************
	SUBROUTINE GRAD(N,X,G)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N),G(N)
	G(1)=X(1)**3-X(1)+.1
	G(2)=X(2)
	RETURN
	END
!*************************************************	
!	HESSIANO
!*************************************************
	SUBROUTINE HESS(N,X,H)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N),H(N,N)
	H(1,1)=3.*X(1)**2-1.
	H(1,2)=0.
	H(2,1)=H(1,2)
	H(2,2)=1.
	RETURN
	END
!*************************************************
!	HESSIANO x DIREZIONE
!*************************************************
	SUBROUTINE HESDIR(N,X,D,HD)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N),HD(N),D(N),H(2,2)
      	CALL HESS(N,X,H)
        DO 2 I=1,N
      	HD(I)=0.	 
      	DO 2 J=1,N
    2  	HD(I)=HD(I)+H(I,J)*D(J)					
     	RETURN
	END
