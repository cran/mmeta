C
C       ============================================================
C       Purpose: This program computes the hypergeometric function 
C                F(a,b,c,x) using subroutine HYGFX
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument ( x � 1 )
C       Output:  HF --- F(a,b,c,x)
C       Example:
C              b = 3.30,  c = 6.70
C              a     F(a,b,c,.25)     F(a,b,c,.55)    F(a,b,c,.85)
C            ------------------------------------------------------
C            -2.5   .72356129D+00    .46961432D+00   .29106096D+00
C            -0.5   .93610145D+00    .85187390D+00   .75543187D+00
C             0.5   .10689695D+01    .11795358D+01   .13510497D+01
C             2.5   .14051563D+01    .23999063D+01   .57381566D+01
C
C              a = 3.30,  b = 6.70
C              c     F(a,b,c,.25)     F(a,b,c,.55)    F(a,b,c,.85)
C            ------------------------------------------------------
C            -5.5   .15090670D+05    .10170778D+11   .58682088D+19
C            -0.5  -.21631479D+04   -.30854772D+07  -.10217370D+13
C             0.5   .26451677D+03    .11967860D+06   .92370648D+10
C             4.5   .41946916D+01    .58092729D+02   .20396914D+05
C       ============================================================
C
        SUBROUTINE HYGFX(A,B,C,X,HF,IERR)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C                IERR --- Error flag
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        DIMENSION G(26)
        
        DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
        
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IERR=0
        IF (C.EQ.0.0D0.OR.C.EQ.(-INT(DABS(C)))) THEN
           IERR=1
           RETURN
        ENDIF
        L0=C.EQ.INT(C).AND.C.LT.0.0D0
        L1=DABS(1.0D0-X).LT.1.0D-15
        L2=DABS(X+1.0D0).LT.1.0D-15
        L3=DABS(X).LT.1.0D-15
        L4=DABS(X-1.0D0).LT.1.0D-15
        L5=DABS(X+1.0D0).LT.1.0D-15
        IF (L0.OR.L1) THEN
           IERR=2
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95D0) EPS=1.0D-8
        IF (X.EQ.0.0D0.OR.A.EQ.0.0D0.OR.B.EQ.0.0D0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (L3.OR.A.EQ.C.OR.B.EQ.C) THEN
           HF=1.0D0
           RETURN
        ELSE IF (L4.AND.C-A-B.GT.0.0D0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (L5.AND.C-A+B.LE.0.0D0) THEN
           A1=A
           B1=C-B
           C1=A+1.0D0-C
           CALL GAMMA(A1,GA1)
           CALL GAMMA(B1,GB1)
           CALL GAMMA(C1,GC1)
           CALL GAMMA(A1+B1,GAB1)
           HF=GA1*GB1/(GC1*GAB1)*((-X)**(-A))
           RETURN
        ELSE IF (L5.AND.C-A+B.GT.0.0D0) THEN
           A1=A
           B1=C-B
           C1=A+1.0D0-C
           CALL GAMMA(A1,GA1)
           CALL GAMMA(B1,GB1)
           CALL GAMMA(C1,GC1)
           CALL GAMMA(A1+B1,GAB1)
           HF=GA1*GB1/(GC1*GAB1)*((-X)**(-A))
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0D0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B).LT.1.0D-15) THEN
              M=INT(C-A-B+EPS)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO J=1,ABS(M)-1
                 GM=GM*J
              END DO
              RM=1.0D0
              DO J=1,ABS(M)
                 RM=RM*J
              END DO
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=GM*GC/(GA*GB)*(1.0D0-X)**(M)
                 DO K=1,M
                    R0=R0*(A+K-1.0D0)*(B+K-1.0D0)/
     &                 (K*(A+B+K-1.0D0))*(1.0D0-X)
                 END DO
                 F0=C0*R0
                 DO K=1,M
                    SP0=SP0+1.0D0/(A+K-1.0D0)+1.0D0/(B+K-1.0D0)-
     &                  1.0D0/K-1.0D0/(A+B+K-1.0D0)
                 END DO
                 F1=C1*(PA+PB+2.0D0*EL+SP0+DLOG(1.0D0-X))
                 DO K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0D0))+
     &                 (1.0D0-B)/(K*(B+K-1.0D0))
                    SM=0.0D0
                    DO J=1,M
                       SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0D0))+
     &                    1.0D0/(B+J+K-1.0D0)
                    END DO
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0D0)/
     &                 (K*(A+B+M+K-1.0D0))*(1.0D0-X)
                    F1=F1+C1*R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 15
                    HW=F1
                 END DO
 15              HF=F0+F1
              ELSE
                 M=-M
                 C0=GM*GC/(GA*GB)*(1.0D0-X)**(M)
                 C1=GM*GC/(GAM*GBM)
                 DO K=1,M
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0D0)/
     &                 (K*(A+B-M+K-1.0D0))*(1.0D0-X)
                 END DO
                 F0=C0*R0
                 DO K=1,M
                    SP0=SP0+1.0D0/(A-M+K-1.0D0)+1.0D0/(B-M+K-1.0D0)
     &                  -1.0D0/K-1.0D0/(A+B-M+K-1.0D0)
                 END DO
                 F1=C1*(PA+PB+2.0D0*EL+SP0+DLOG(1.0D0-X))
                 DO K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0D0))+(1.0D0-B)/
     &                 (K*(B+K-1.0D0))
                    SM=0.0D0
                    DO J=1,M
                       SM=SM+1.0D0/(A+J+K-1.0D0)+1.0D0/
     &                    (B+J+K-1.0D0)
                    END DO
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0D0)/(K*(A+B+K-1.0D0))
     &                 *(1.0D0-X)
                    F1=F1+C1*R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 25
                    HW=F1
                 END DO
 25              HF=F0+F1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0D0)/
     &              (K*(A+B-C+K))*(1.0D0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0D0)/
     &              (K*(C-A-B+K))*(1.0D0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 35
                 HW=HF
              END DO
 35           HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.GT.B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 45
              HW=HF
           END DO
 45        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) IERR=3
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function (x)
C       Input :  x  --- Argument of gamma function
C                       ( x is not equal to 0,-1,-2,...)
C       Output:  GA --- (x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        
        DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
        
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO K=2,M1
                 GA=GA*K
              END DO
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO K=1,M
                 R=R*(Z-K)
              END DO
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           GR=G(26)
           DO K=25,1,-1
              GR=GR*Z+G(K)
           END DO
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE PSI(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=IDINT(XA)
           DO K=1,N-1
              S=S+1.0D0/K
           END DO
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=IDINT(XA-.5)
           DO K=1,N
              S=S+1.0/(2.0D0*K-1.0D0)
           END DO
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO K=0,N-1
                 S=S+1.0D0/(XA+K)
              END DO
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+A6)*X2+A5)*X2
     &        +A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
