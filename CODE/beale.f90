! n = 2
! -4.5 <= x <= 4.5
! fstar = 0
!**************************************************
!	FUNZIONE BEALE
!**************************************************

      SUBROUTINE STARTP(N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N)

      X(1)=1.0D0
      X(2)=1.0D0

      RETURN
      END
!**************************************************
!	FUNZIONE BEALE
!**************************************************
	SUBROUTINE FUNCT(N,X,F)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (M=3)
	DIMENSION X(N),Y(M)
        Y(1)=1.5D0
        Y(2)=2.25D0
        Y(3)=2.625D0
        IF(DABS(1.D0-X(2)).GT.1.D+16) GO TO 1
        IF(DABS(X(1)).GT.1.D+5) GO TO 1
        IF(DABS(X(2)).GT.1.D+5) GO TO 1
	F=(Y(1)-X(1)*(1.D0-X(2)))**2
        IF(DABS(1.D0-X(2)**2).GT.1.D+16) GO TO 1
	F=F+(Y(2)-X(1)*(1.D0-X(2)**2))**2
        IF(DABS(1.D0-X(2)**3).GT.1.D+16) GO TO 1
	F=F+(Y(3)-X(1)*(1.D0-X(2)**3))**2
	RETURN
    1   F=1.D+34
        RETURN
	END
!**************************************************
!	GRADIENTE
!**************************************************
	SUBROUTINE GRAD(N,X,G)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (M=3)
	DIMENSION X(N),Y(M),G(N)
        Y(1)=1.5D0
        Y(2)=2.25D0
        Y(3)=2.625D0
	G(1)=-2.D0*(Y(1)-X(1)*(1.D0-X(2)))*(1.D0-X(2))
	G(1)=G(1)-2.D0*(Y(2)-X(1)*(1.D0-X(2)**2))*(1.D0-X(2)**2)
	G(1)=G(1)-2.D0*(Y(3)-X(1)*(1.D0-X(2)**3))*(1.D0-X(2)**3)

	G(2)=2.D0*(Y(1)-X(1)*(1.D0-X(2)))*X(1)
	G(2)=G(2)+2.D0*(Y(2)-X(1)*(1.D0-X(2)**2))*(2.D0*X(1)*X(2))
	G(2)=G(2)+2.D0*(Y(3)-X(1)*(1.D0-X(2)**3))*(3.D0*X(1)*X(2)**2)
	RETURN
	END
!**************************************************	
!	HESSIANO
!**************************************************
	SUBROUTINE HESS(N,X,H)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (M=3)
	DIMENSION X(N),Y(M),H(N,N)
        Y(1)=1.5D0
        Y(2)=2.25D0
        Y(3)=2.625D0
	H(1,1)=2.D0*(1.D0-X(2))*(1.D0-X(2))
	H(1,1)=H(1,1)+2.D0*(1.D0-X(2)**2)*(1.D0-X(2)**2)
	H(1,1)=H(1,1)+2.D0*(1.D0-X(2)**3)*(1.D0-X(2)**3)
      
	H(2,2)=2.D0*X(1)*X(1)
	H(2,2)=H(2,2)+2.D0*(2.D0*X(1)*X(2))**2 +2.D0*(Y(2)-X(1)*(1.D0-X(2)**2))*(2.D0*X(1))
	H(2,2)=H(2,2)+2.D0*(3.D0*X(1)*X(2)**2)**2 +2.D0*(Y(3)-X(1)*(1.D0-X(2)**3))*(6.D0*X(1)*X(2))

	H(1,2)=-2.D0*(1.D0-X(2))*X(1) +2.D0*(Y(1)-X(1)*(1.D0-X(2)))
	H(1,2)=H(1,2)-2.D0*(1.D0-X(2)**2)*(2.D0*X(1)*X(2)) +2.D0*(Y(2)-X(1)*(1.D0-X(2)**2))*(2.D0*X(2))
	H(1,2)=H(1,2)-2.D0*(1.D0-X(2)**3)*(3.D0*X(1)*X(2)**2) +2.D0*(Y(3)-X(1)*(1.D0-X(2)**3))*(3.D0*X(2)**2)

	H(2,1)=H(1,2)

	RETURN
	END
!**************************************************
!	HESSIANO x DIREZIONE
!**************************************************
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




