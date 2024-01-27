!*************************************************
!	FUNZIONE LEVY 5
!*************************************************
      SUBROUTINE FUNCT(N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N)

      PG=DACOS(-1.D0)
      CT=25.D-2
      RK=10.D0
      RA=1.D0
      NM=N-1

      F=RK*DSIN(PG*(1.D0+CT*(X(1)-1.D0)))**2

      DO 1 I=1,NM
      IP=I+1
      F=F+(((1.D0+CT*(X(I)-1.D0))-RA)**2)*(1.D0+RK*DSIN(PG* &
     (1.D0+CT*(X(I+1)-1.D0)))**2)
    1 CONTINUE

      F=F+((1.D0+CT*(X(N)-1.D0))-RA)**2

      F=(PG/N)*F

      RETURN
      END

!*************************************************
!	GRADIENTE
!*************************************************
      SUBROUTINE GRAD(N,X,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),G(N)

      PG=DACOS(-1.D0)
      CT=25.D-2
      RK=10.D0
      RA=1.D0
      NM=N-1

      G(1)=2.*PG*RK*DSIN(PG*(1.D0+CT*(X(1)-1.D0)))&
     *DCOS(PG*(1.D0+CT*(X(1)-1.D0)))
      
      G(1)=G(1)+2.*((1.D0+CT*(X(1)-1.D0))-RA)*(1.D0+&
     RK*DSIN(PG*(1.D0+CT*(X(2)-1.D0)))**2)
      G(1)=(PG/N)*CT*G(1)

      DO 1 I=2,NM
      IM=I-1
      IP=I+1

      G(I)=2.*RK*PG*(((1.D0+CT*(X(IM)-1.D0))-RA)**2)*&
     DSIN(PG*(1.D0+CT*(X(I)-1.D0)))*DCOS(PG*(1.D0+CT*(X(I)-1.D0)))
      G(I)=G(I)+2.*((1.D0+CT*(X(I)-1.D0))-RA)*(1.D0+RK*&
     DSIN(PG*(1.D0+CT*(X(IP)-1.D0)))**2)
      G(I)=(PG/N)*CT*G(I)

    1 CONTINUE

      G(N)=2.*RK*PG*(((1.D0+CT*(X(N-1)-1.D0))-RA)**2)*&
     DSIN(PG*(1.D0+CT*(X(N)-1.D0)))*DCOS(PG*(1.D0+CT*(X(N)-1.D0)))
      G(N)=G(N)+2.*((1.D0+CT*(X(N)-1.D0))-RA)
      G(N)=(PG/N)*CT*G(N)

      RETURN
      END
!*************************************************
!	HESSIANO
!*************************************************
      
      SUBROUTINE HESS(N,X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),H(N,N)

      PG=DACOS(-1.D0)
      CT=25.D-2
      RK=10.D0
      RA=1.D0
      NM=N-1

      DO 1 I=1,N
      DO 1 J=1,N
    1 H(I,J)=0.

      H(1,1)=2.*RK*(PG**2)*(DCOS(PG*(1.D0+CT*(X(1)-1.D0)))**2&
     -DSIN(PG*(1.D0+CT*(X(1)-1.D0)))**2)
      H(1,1)=H(1,1)+2.*(1.D0+RK*DSIN(PG*(1.D0+CT*(X(2)-1.D0)))**2)
      H(1,1)=(PG/N)*CT*CT*H(1,1)

      H(1,2)=4.*RK*PG*((1.D0+CT*(X(1)-1.D0))-RA)*&
     DSIN(PG*(1.D0+CT*(X(2)-1.D0)))*DCOS(PG*(1.D0+CT*(X(2)-1.D0)))

      H(1,2)=(PG/N)*CT*CT*H(1,2)

      H(2,1)=H(1,2)

      DO 2 I=2,NM
      IP=I+1
      IM=I-1
      
      
      H(I,I)=2.*RK*(PG**2)*(((1.D0+CT*(X(IM)-1.D0))-RA)**2)&
	  *(DCOS(PG*(1.D0+CT*(X(I)-1.D0)))**2-DSIN(PG*(1.D0+CT*(X(I)-1.D0)))**2)
      H(I,I)=H(I,I)+2.*(1.D0+RK*DSIN(PG*(1.D0+CT*(X(IP)-1.D0)))**2)
      H(I,I)=(PG/N)*CT*CT*H(I,I)

      H(I,IP)=4.*RK*PG*((1.D0+CT*(X(I)-1.D0))-RA)&
     *DSIN(PG*(1.D0+CT*(X(IP)-1.D0)))*DCOS(PG*(1.D0+CT*(X(IP)-1.D0)))
      H(I,IP)=(PG/N)*CT*CT*H(I,IP)

      H(IP,I)=H(I,IP)

    2 CONTINUE

      H(N,N)=2.*RK*(PG**2)*(((1.D0+CT*(X(N-1)-1.D0))-RA)**2)&
     *(DCOS(PG*(1.D0+CT*(X(N)-1.D0)))**2-DSIN(PG*(1.D0+CT*(X(N)-1.D0)))**2)+2.

      H(N,N)=(PG/N)*CT*CT*H(N,N)

      RETURN
      END
!*************************************************
!	HESSIANO x DIREZIONE
!*************************************************
      SUBROUTINE HESDIR(N,X,D,HD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),HD(N),D(N)
      NM=N-1
      DO 1 I=1,N
      HD(I)=0.
    1 CONTINUE
      C1=2.+1200.*X(1)**2-400.*X(2)
      C2=-400.*X(1)
      HD(1)=HD(1)+C1*D(1)+C2*D(2)
      HD(2)=HD(2)+C2*D(1)
      DO 2 I=2,NM
      IP=I+1
      C1=202.+1200.*X(I)**2-400.*X(IP)
      C2=-400.*X(I)
      HD(I)=HD(I)+C1*D(I)+C2*D(I+1)
      HD(I+1)=HD(I+1)+C2*D(I)
    2 CONTINUE
      HD(N)=HD(N)+200.*D(N)
      RETURN
      END
