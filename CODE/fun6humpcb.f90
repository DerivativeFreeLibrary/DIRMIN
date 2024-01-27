!*************************************************
!	FUNZIONE CAMEL
!*************************************************
      	SUBROUTINE FUNCT(N,X,F)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N)
	F=(4.-2.1*X(1)**2+(X(1)**4)/3.)*X(1)**2+X(1)*X(2) &
     	+(-4.+4.*X(2)**2)*X(2)**2
       	RETURN
	END
!*************************************************
!	GRADIENTE
!*************************************************
	SUBROUTINE GRAD(N,X,G)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N),G(N)
	G(2)=X(1)-8.*X(2)+16.*X(2)**3
	G(1)=X(2)+8.*X(1)-8.4*X(1)**3+2.*X(1)**5
	RETURN
	END
!*************************************************	
!	HESSIANO
!*************************************************
	SUBROUTINE HESS(N,X,H)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(N),H(N,N)
	H(1,1)=8.-25.2*X(1)**2+10*X(1)**4
	H(1,2)=1.
	H(2,1)=H(1,2)
	H(2,2)=-8.+48.*X(2)**2
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
