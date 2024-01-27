!*************************************************
!	FUNZIONE ZIRILLI 14
!*************************************************
      SUBROUTINE FUNCT(N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N)

      PG=DACOS(-1.D0)
      COST=.1
      NM=N-1

      F=DSIN(3.D0*PG*X(1))**2

      DO 1 I=1,NM
      IP=I+1
      F=F+((X(I)-1.D0)**2)*(1.D0+DSIN(3.D0*PG*X(IP))**2)
    1 CONTINUE

      F=F+((X(N)-1.D0)**2)*(1.D0+DSIN(2.D0*PG*X(N))**2)
      
      F=COST*F

      RETURN
      END

!*************************************************
!	GRADIENTE
!*************************************************
      SUBROUTINE GRAD(N,X,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),G(N)

      PG=DACOS(-1.D0)
      COST=.1
      NM=N-1

      G(1)=6.*DSIN(3.D0*PG*X(1))*DCOS(3.D0*PG*X(1))*PG
      G(1)=G(1)+2.*(X(1)-1.D0)*(1.D0+DSIN(3.D0*PG*X(2))**2)
      G(1)=COST*G(1)

      DO 1 I=2,NM
      IM=I-1
      IP=I+1

      G(I)=6.*((X(IM)-1.)**2)*DSIN(3.D0*PG*X(I))*DCOS(3.D0*PG*X(I))*PG
      G(I)=G(I)+2.*(X(I)-1.D0)*(1.D0+DSIN(3.D0*PG*X(IP))**2)
      G(I)=COST*G(I)

    1 CONTINUE

      G(N)=6.*((X(N-1)-1.)**2)*DSIN(3.D0*PG*X(N))*DCOS(3.D0*PG*X(N))*PG
      G(N)=G(N)+2.*(X(N)-1.D0)*(1.D0+DSIN(2.D0*PG*X(N))**2)
      G(N)=G(N)+4.*((X(N)-1.)**2)*DSIN(2.D0*PG*X(N))*DCOS(2.D0*PG*X(N))*PG
      G(N)=COST*G(N)

      RETURN
      END
!*************************************************
!	HESSIANO
!*************************************************
      SUBROUTINE HESS(N,X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),H(N,N)

      PG=DACOS(-1.D0)
      COST=.1
      NM=N-1

      DO 1 I=1,N
      DO 1 J=1,N
    1 H(I,J)=0.

      H(1,1)=18.*(PG**2)*(DCOS(3.D0*PG*X(1))**2-DSIN(3.D0*PG*X(1))**2)
      H(1,1)=H(1,1)+2.*(1.D0+DSIN(3.D0*PG*X(2))**2)
      H(1,1)=COST*H(1,1)
      
      H(1,2)=12.*PG*(X(1)-1.D0)*DSIN(3.D0*PG*X(2))*DCOS(3.D0*PG*X(2))
      H(1,2)=COST*H(1,2)

      H(2,1)=H(1,2)

      DO 2 I=2,NM
      IP=I+1
      IM=I-1

      H(I,I)=18.*(PG**2)*((X(IM)-1.)**2)*(DCOS(3.D0*PG*X(I))**2&
      -DSIN(3.D0*PG*X(I))**2)
      H(I,I)=H(I,I)+2.*(1.D0+DSIN(3.D0*PG*X(IP))**2)
      H(I,I)=COST*H(I,I)

      H(I,IP)=12.*PG*(X(I)-1.)*DSIN(3.D0*PG*X(IP))*DCOS(3.D0*PG*X(IP))
      H(I,IP)=COST*H(I,IP)
      H(IP,I)=H(I,IP)

    2 CONTINUE

      H(N,N)=18.*(PG**2)*((X(N-1)-1.)**2)*(DCOS(3.D0*PG*X(N))**2 &
     -DSIN(3.D0*PG*X(N))**2)+2.*(1.D0+DSIN(2.D0*PG*X(N))**2) &
     +8.*(X(N)-1.)*DSIN(2.D0*PG*X(N))*DCOS(2.D0*PG*X(N))*PG &
     +8.*(X(N)-1.)*DSIN(2.D0*PG*X(N))*DCOS(2.D0*PG*X(N))*PG &
     +8.*(PG**2)*((X(N)-1.)**2)*(DCOS(2.D0*PG*X(N))**2-DSIN(2.D0*PG*X(N))**2)

      H(N,N)=COST*H(N,N)

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
