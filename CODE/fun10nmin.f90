!*************************************************
!	FUNZIONE LEVY 8
!*************************************************
      SUBROUTINE FUNCT(N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N)

      PG=DACOS(-1.D0)
      RK=10.D0
      RA=1.D0
      NM=N-1

      F=RK*DSIN(PG*X(1))**2

      DO 1 I=1,NM
      IP=I+1
      F=F+((X(I)-RA)**2)*(1.D0+RK*DSIN(PG*X(I+1))**2)
    1 CONTINUE

      F=F+(X(N)-RA)**2

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
      RK=10.D0
      RA=1.D0
      NM=N-1

      G(1)=2.*PG*RK*DSIN(PG*X(1))*DCOS(PG*X(1))
      G(1)=G(1)+2.*(X(1)-RA)*(1.D0+RK*DSIN(PG*X(2))**2)
      G(1)=(PG/N)*G(1)

      DO 1 I=2,NM
      IM=I-1
      IP=I+1

      G(I)=2.*RK*PG*((X(IM)-RA)**2)*DSIN(PG*X(I))*DCOS(PG*X(I))
      G(I)=G(I)+2.*(X(I)-RA)*(1.D0+RK*DSIN(PG*X(IP))**2)
      G(I)=(PG/N)*G(I)

    1 CONTINUE

      G(N)=2.*RK*PG*((X(N-1)-RA)**2)*DSIN(PG*X(N))*DCOS(PG*X(N))
      G(N)=G(N)+2.*(X(N)-RA)
      G(N)=(PG/N)*G(N)

      RETURN
      END
!*************************************************
!	HESSIANO
!*************************************************
      
      SUBROUTINE HESS(N,X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),H(N,N)

      PG=DACOS(-1.D0)
      RK=10.D0
      RA=1.D0
      NM=N-1

      DO 1 I=1,N
      DO 1 J=1,N
    1 H(I,J)=0.

      H(1,1)=2.*RK*(PG**2)*(DCOS(PG*X(1))**2-DSIN(PG*X(1))**2)
      H(1,1)=H(1,1)+2.*(1.D0+RK*DSIN(PG*X(2))**2)
      H(1,1)=(PG/N)*H(1,1)

      H(1,2)=4.*RK*PG*(X(1)-RA)*DSIN(PG*X(2))*DCOS(PG*X(2))

      H(1,2)=(PG/N)*H(1,2)

      H(2,1)=H(1,2)

      DO 2 I=2,NM
      IP=I+1
      IM=I-1
      
      
      H(I,I)=2.*RK*(PG**2)*((X(IM)-RA)**2)*(DCOS(PG*X(I))**2-&
     DSIN(PG*X(I))**2)
      H(I,I)=H(I,I)+2.*(1.D0+RK*DSIN(PG*X(IP))**2)
      H(I,I)=(PG/N)*H(I,I)

      H(I,IP)=4.*RK*PG*(X(I)-RA)*DSIN(PG*X(IP))*DCOS(PG*X(IP))
      H(I,IP)=(PG/N)*H(I,IP)

      H(IP,I)=H(I,IP)

    2 CONTINUE

      H(N,N)=2.*RK*(PG**2)*((X(N-1)-RA)**2)*(DCOS(PG*X(N))**2-DSIN(PG*X(N))**2)+2.

      H(N,N)=(PG/N)*H(N,N)

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
