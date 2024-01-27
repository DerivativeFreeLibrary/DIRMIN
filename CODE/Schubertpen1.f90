!**************************************************
!	FUNZIONE SCHUBERT
!**************************************************
      SUBROUTINE FUNCT(N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N)

      F1=0.D0
      F2=0.D0

      DO 1 I=1,5
      F1=F1+I*DCOS((I+1.D0)*X(1)+I)
      F2=F2+I*DCOS((I+1.D0)*X(2)+I)
    1 CONTINUE

      F=F1*F2
      
      F=F+((X(1)+1.42513)**2+(X(2)+0.80032)**2)/2.D0

      RETURN
      END

!**************************************************
!	GRADIENTE
!**************************************************
      SUBROUTINE GRAD(N,X,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),G(N)

      F1=0.D0
      F2=0.D0

      DO 1 I=1,5
      F1=F1+I*DCOS((I+1.D0)*X(1)+I)
      F2=F2+I*DCOS((I+1.D0)*X(2)+I)
    1 CONTINUE


      G(1)=0.D0
      G(2)=0.D0

      DO 2 I=1,5

      G(1)=G(1)-F2*(I+1.)*I*DSIN((I+1.D0)*X(1)+I)
      G(2)=G(2)-F1*(I+1.)*I*DSIN((I+1.D0)*X(2)+I)

    2 CONTINUE

      G(1)=G(1)+X(1)+1.42513
      G(2)=G(2)+X(2)+0.80032

      RETURN
      END
!**************************************************
!	HESSIANO
!**************************************************
      
      SUBROUTINE HESS(N,X,H)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),H(N,N)

      F1=0.D0
      F2=0.D0

      DO 1 I=1,5
      F1=F1+I*DCOS((I+1.D0)*X(1)+I)
      F2=F2+I*DCOS((I+1.D0)*X(2)+I)
    1 CONTINUE

      DO 3 I=1,N
      DO 3 J=1,N
      H(I,J)=0.D0
    3 CONTINUE

      DO 2 I=1,5

      H(1,1)=H(1,1)-F2*((I+1.)**2)*I*DCOS((I+1.D0)*X(1)+I)
      H(2,2)=H(2,2)-F1*((I+1.)**2)*I*DCOS((I+1.D0)*X(2)+I)

    2 CONTINUE

      H(1,1)=H(1,1)+1.
      H(2,2)=H(2,2)+1.


      F11=0.D0
      F22=0.D0

      DO 4 I=1,5
      F11=F11-I*(I+1.)*DSIN((I+1.D0)*X(1)+I)
      F22=F22-I*(I+1.)*DSIN((I+1.D0)*X(2)+I)
    4 CONTINUE

      H(1,2)=F11*F22
      H(2,1)=H(1,2)

      RETURN
      END
!**************************************************
!	HESSIANO x DIREZIONE
!**************************************************
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
