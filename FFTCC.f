C   SUBROUTINE FFTCC (A,N,IWK,WK)
C./ ADD NAME=FFTCC
C./ NUMBER NEW1=00000010,INCR=10
C   IMSL ROUTINE NAME   - FFTCC
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - IBM77/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A
C                           COMPLEX VALUED SEQUENCE
C
C   USAGE               - CALL FFTCC (A,N,IWK,WK)
C
C   ARGUMENTS    A      - COMPLEX VECTOR OF LENGTH N. ON INPUT A
C                           CONTAINS THE COMPLEX VALUED SEQUENCE TO BE
C                           TRANSFORMED. ON OUTPUT A IS REPLACED BY THE
C                           FOURIER TRANSFORM.
C                N      - INPUT NUMBER OF DATA POINTS TO BE
C                           TRANSFORMED. N MAY BE ANY POSITIVE
C                           INTEGER.
C                IWK    - INTEGER WORK VECTOR OF LENGTH 6*N+150.
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
C                WK     - REAL WORK VECTOR OF LENGTH 6*N+150.
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  FFTCC COMPUTES THE FOURIER TRANSFORM, X, ACCORDING
C                TO THE FOLLOWING FORMULA;
C
C                  X(K+1) = SUM FROM J = 0 TO N-1 OF
C                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))
C                  FOR K=0,1,...,N-1 AND PI=3.1415...
C
C                NOTE THAT X OVERWRITES A ON OUTPUT.
C            2.  FFTCC CAN BE USED TO COMPUTE
C
C                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF
C                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))
C                  FOR K=0,1,...,N-1 AND PI=3.1415...
C
C                BY PERFORMING THE FOLLOWING STEPS;
C
C                     DO 10 I=1,N
C                        A(I) = CONJG(A(I))
C                  10 CONTINUE
C                     CALL FFTCC (A,N,IWK,WK)
C                     DO 20 I=1,N
C                        A(I) = CONJG(A(I))/N
C                  20 CONTINUE
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FFTCC (A,N,IWK,WK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IWK(1)
      DOUBLE PRECISION   WK(1)
      COMPLEX*16         A(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IAM,IAP,IBM,IBP,IC,ICC,ICF,ICK,ID,IDM1,II,
     1                   IJA,IKB,IKT,ILL,IM,IRD,ISF,ISK,ISP,ISS,ITA,ITB,
     2                   J,JA,JF,JJ,JK,K,K0,K1,K2,K3,KA,KB,KD2,KF,KH,KN,
     3                   KT,KTP,L,L1,M,MM,MM1,MP
      DOUBLE PRECISION   CM,SM,C1,C2,C3,S1,S2,S3,C30,RAD,A0,A1,A4,B4,
     1                   A2,A3,B0,B1,B2,B3,ZERO,HALF,ONE,TWO,Z0(2),
     2                   Z1(2),Z2(2),Z3(2),Z4(2)
      COMPLEX*16         ZA0,ZA1,ZA2,ZA3,ZA4,AK2
      EQUIVALENCE        (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)),
     1                   (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),(A1,Z1(1)),
     2                   (B1,Z1(2)),(A2,Z2(1)),(B2,Z2(2)),(A3,Z3(1)),
     3                   (B3,Z3(2)),(ZA4,Z4(1)),(Z4(1),A4),(Z4(2),B4)
      DATA               RAD/6.283185307179586D0/,
     1                   C30/.8660254037844386D0/
      DATA               ZERO,HALF,ONE,TWO/0.0D0,0.5D0,1.0D0,2.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .EQ. 1) GO TO 9005
      K = N
      M = 0
      J = 2
      JJ = 4
      JF = 0
C                                  DETERMINE THE SQUARE FACTORS OF N
      IWK(1) = 1
    5 I = K/JJ
      IF (I*JJ .NE. K) GO TO 10
      M = M+1
      IWK(M+1) = J
      K = I
      GO TO 5
   10 J = J + 2
      IF (J .EQ. 4) J = 3
      JJ = J * J
      IF (JJ .LE. K) GO TO 5
      KT = M
C                                  DETERMINE THE REMAINING FACTORS OF N
      J = 2
   15 I = K / J
      IF (I*J .NE. K) GO TO 20
      M = M + 1
      IWK(M+1) = J
      K = I
      GO TO 15
   20 J = J + 1
      IF (J .EQ. 3) GO TO 15
      J = J + 1
      IF(J.LE.K) GO TO 15
      K = IWK(M+1)
      IF (IWK(KT+1) .GT. IWK(M+1)) K = IWK(KT+1)
      IF(KT.LE.0) GO TO 30
      KTP = KT + 2
      DO 25  I = 1,KT
         J = KTP - I
         M = M+1
         IWK(M+1) = IWK(J)
   25 CONTINUE
   30 MP = M+1
      IC = MP+1
      ID = IC+MP
      ILL = ID+MP
      IRD = ILL+MP+1
      ICC = IRD+MP
      ISS = ICC+MP
      ICK = ISS+MP
      ISK = ICK+K
      ICF = ISK+K
      ISF = ICF+K
      IAP = ISF+K
      KD2 = (K-1) / 2 + 1
      IBP = IAP + KD2
      IAM = IBP + KD2
      IBM = IAM + KD2
      MM1 = M-1
      I = 1
   35 L = MP - I
      J = IC - I
      IWK(ILL+L) = 0
      IF ((IWK(J-1) + IWK(J)) .EQ. 4) IWK(ILL+L) = 1
      IF (IWK(ILL+L) .EQ. 0) GO TO 40
      I = I + 1
      L = L - 1
      IWK(ILL+L) = 0
   40 I = I + 1
      IF(I.LE.MM1) GO TO 35
      IWK(ILL+1) = 0
      IWK(ILL+MP) = 0
      IWK(IC) = 1
      IWK(ID) = N
      DO 45  J = 1,M
         K = IWK(J+1)
         IWK(IC+J) = IWK(IC+J-1) * K
         IWK(ID+J) = IWK(ID+J-1) / K
         WK(IRD+J) = RAD/IWK(IC+J)
         C1 = RAD/K
         IF (K .LE. 2) GO TO 45
         WK(ICC+J) = DCOS(C1)
         WK(ISS+J) = DSIN(C1)
   45 CONTINUE
      MM = M
      IF (IWK(ILL+M) .EQ. 1) MM = M - 1
      IF (MM .LE. 1) GO TO 50
      SM = IWK(IC+MM-2) * WK(IRD+M)
      CM = DCOS(SM)
      SM = DSIN(SM)
   50 KB = 0
      KN = N
      JJ = 0
      I = 1
      C1 = ONE
      S1 = ZERO
      L1 = 1
   55 IF (IWK(ILL+I+1) .EQ. 1) GO TO 60
      KF = IWK(I+1)
      GO TO 65
   60 KF = 4
      I = I+1
   65 ISP = IWK(ID+I)
      IF (L1 .EQ. 1) GO TO 70
      S1 = JJ * WK(IRD+I)
      C1 = DCOS(S1)
      S1 = DSIN(S1)
C                                  FACTORS OF 2, 3, AND 4 ARE
C                                  HANDLED SEPARATELY.
   70 IF (KF .GT. 4) GO TO 140
      GO TO (75,75,90,115), KF
   75 K0 = KB + ISP
      K2 = K0 + ISP
      IF (L1 .EQ. 1) GO TO 85
   80 K0 = K0 - 1
      IF (K0 .LT. KB) GO TO 190
      K2 = K2 - 1
      ZA4 = A(K2+1)
      A0 = A4*C1-B4*S1
      B0 = A4*S1+B4*C1
      A(K2+1) = A(K0+1)-ZA0
      A(K0+1) = A(K0+1)+ZA0
      GO TO 80
   85 K0 = K0 - 1
      IF (K0 .LT. KB) GO TO 190
      K2 = K2 - 1
      AK2 = A(K2+1)
      A(K2+1) = A(K0+1)-AK2
      A(K0+1) = A(K0+1)+AK2
      GO TO 85
   90 IF (L1 .EQ. 1) GO TO 95
      C2 = C1 * C1 - S1 * S1
      S2 = TWO * C1 * S1
   95 JA = KB + ISP - 1
      KA = JA + KB
      IKB = KB+1
      IJA = JA+1
      DO 110 II = IKB,IJA
         K0 = KA - II + 1
         K1 = K0 + ISP
         K2 = K1 + ISP
         ZA0 = A(K0+1)
         IF (L1 .EQ. 1) GO TO 100
         ZA4 = A(K1+1)
         A1 = A4*C1-B4*S1
         B1 = A4*S1+B4*C1
         ZA4 = A(K2+1)
         A2 = A4*C2-B4*S2
         B2 = A4*S2+B4*C2
         GO TO 105
  100    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
  105    A(K0+1) = DCMPLX(A0+A1+A2,B0+B1+B2)
         A0 = -HALF * (A1+A2) + A0
         A1 = (A1-A2) * C30
         B0 = -HALF * (B1+B2) + B0
         B1 = (B1-B2) * C30
         A(K1+1) = DCMPLX(A0-B1,B0+A1)
         A(K2+1) = DCMPLX(A0+B1,B0-A1)
  110 CONTINUE
      GO TO 190
  115 IF (L1 .EQ. 1) GO TO 120
      C2 = C1 * C1 - S1 * S1
      S2 = TWO * C1 * S1
      C3 = C1 * C2 - S1 * S2
      S3 = S1 * C2 + C1 * S2
  120 JA = KB + ISP - 1
      KA = JA + KB
      IKB = KB+1
      IJA = JA+1
      DO 135 II = IKB,IJA
         K0 = KA - II + 1
         K1 = K0 + ISP
         K2 = K1 + ISP
         K3 = K2 + ISP
         ZA0 = A(K0+1)
         IF (L1 .EQ. 1) GO TO 125
         ZA4 = A(K1+1)
         A1 = A4*C1-B4*S1
         B1 = A4*S1+B4*C1
         ZA4 = A(K2+1)
         A2 = A4*C2-B4*S2
         B2 = A4*S2+B4*C2
         ZA4 = A(K3+1)
         A3 = A4*C3-B4*S3
         B3 = A4*S3+B4*C3
         GO TO 130
  125    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
         ZA3 = A(K3+1)
  130    A(K0+1) = DCMPLX(A0+A2+A1+A3,B0+B2+B1+B3)
         A(K1+1) = DCMPLX(A0+A2-A1-A3,B0+B2-B1-B3)
         A(K2+1) = DCMPLX(A0-A2-B1+B3,B0-B2+A1-A3)
         A(K3+1) = DCMPLX(A0-A2+B1-B3,B0-B2-A1+A3)
  135 CONTINUE
      GO TO 190
  140 JK = KF - 1
      KH = JK/2
      K3 = IWK(ID+I-1)
      K0 = KB + ISP
      IF (L1 .EQ. 1) GO TO 150
      K = JK - 1
      WK(ICF+1) = C1
      WK(ISF+1) = S1
      DO 145 J = 1,K
         WK(ICF+J+1) = WK(ICF+J) * C1 - WK(ISF+J) * S1
         WK(ISF+J+1) = WK(ICF+J) * S1 + WK(ISF+J) * C1
  145 CONTINUE
  150 IF (KF .EQ. JF) GO TO 160
      C2 = WK(ICC+I)
      WK(ICK+1) = C2
      WK(ICK+JK) = C2
      S2 = WK(ISS+I)
      WK(ISK+1) = S2
      WK(ISK+JK) = -S2
      DO 155 J = 1,KH
         K = JK - J
         WK(ICK+K) = WK(ICK+J) * C2 - WK(ISK+J) * S2
         WK(ICK+J+1) = WK(ICK+K)
         WK(ISK+J+1) = WK(ICK+J) * S2 + WK(ISK+J) * C2
         WK(ISK+K) = -WK(ISK+J+1)
  155 CONTINUE
  160 K0 = K0 - 1
      K1 = K0
      K2 = K0 + K3
      ZA0 = A(K0+1)
      A3 = A0
      B3 = B0
      DO 175 J = 1,KH
         K1 = K1 + ISP
         K2 = K2 - ISP
         IF (L1 .EQ. 1) GO TO 165
         K = KF - J
         ZA4 = A(K1+1)
         A1 = A4*WK(ICF+J)-B4*WK(ISF+J)
         B1 = A4*WK(ISF+J)+B4*WK(ICF+J)
         ZA4 = A(K2+1)
         A2 = A4*WK(ICF+K)-B4*WK(ISF+K)
         B2 = A4*WK(ISF+K)+B4*WK(ICF+K)
         GO TO 170
  165    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
  170    WK(IAP+J) = A1 + A2
         WK(IAM+J) = A1 - A2
         WK(IBP+J) = B1 + B2
         WK(IBM+J) = B1 - B2
         A3 = A1 + A2 + A3
         B3 = B1 + B2 + B3
  175 CONTINUE
      A(K0+1) = DCMPLX(A3,B3)
      K1 = K0
      K2 = K0 + K3
      DO 185 J = 1,KH
         K1 = K1 + ISP
         K2 = K2 - ISP
         JK = J
         A1 = A0
         B1 = B0
         A2 = ZERO
         B2 = ZERO
         DO 180  K = 1,KH
            A1 = A1 + WK(IAP+K) * WK(ICK+JK)
            A2 = A2 + WK(IAM+K) * WK(ISK+JK)
            B1 = B1 + WK(IBP+K) * WK(ICK+JK)
            B2 = B2 + WK(IBM+K) * WK(ISK+JK)
            JK = JK + J
            IF (JK .GE. KF) JK = JK - KF
  180    CONTINUE
         A(K1+1) = DCMPLX(A1-B2,B1+A2)
         A(K2+1) = DCMPLX(A1+B2,B1-A2)
  185 CONTINUE
      IF (K0 .GT. KB) GO TO 160
      JF = KF
  190 IF ( I .GE. MM ) GO TO 195
      I = I + 1
      GO TO 55
  195 I = MM
      L1 = 0
      KB = IWK(ID+I-1) + KB
      IF (KB .GE. KN) GO TO 215
  200 JJ = IWK(IC+I-2) + JJ
      IF (JJ .LT. IWK(IC+I-1)) GO TO 205
      I = I - 1
      JJ = JJ - IWK(IC+I)
      GO TO 200
  205 IF (I .NE. MM) GO TO 210
      C2 = C1
      C1 = CM * C1 - SM * S1
      S1 = SM * C2 + CM * S1
      GO TO 70
  210 IF (IWK(ILL+I) .EQ. 1) I = I + 1
      GO TO 55
  215 I = 1
      JA = KT - 1
      KA = JA + 1
      IF(JA.LT.1) GO TO 225
      DO 220  II = 1,JA
         J = KA - II
         IWK(J+1) = IWK(J+1) - 1
         I = IWK(J+1) + I
  220 CONTINUE
C                                  THE RESULT IS NOW PERMUTED TO
C                                  NORMAL ORDER.
  225 IF (KT .LE. 0) GO TO 270
      J = 1
      I = 0
      KB = 0
  230 K2 = IWK(ID+J) + KB
      K3 = K2
      JJ = IWK(IC+J-1)
      JK = JJ
      K0 = KB + JJ
      ISP = IWK(IC+J) - JJ
  235 K = K0 + JJ
  240 ZA4 = A(K0+1)
      A(K0+1) = A(K2+1)
      A(K2+1) = ZA4
      K0 = K0 + 1
      K2 = K2 + 1
      IF (K0 .LT. K) GO TO 240
      K0 = K0 + ISP
      K2 = K2 + ISP
      IF (K0 .LT. K3) GO TO 235
      IF (K0 .GE. K3 + ISP) GO TO 245
      K0 = K0 - IWK(ID+J) + JJ
      GO TO 235
  245 K3 = IWK(ID+J) + K3
      IF (K3 - KB .GE. IWK(ID+J-1)) GO TO 250
      K2 = K3 + JK
      JK = JK + JJ
      K0 = K3 - IWK(ID+J) + JK
      GO TO 235
  250 IF (J .GE. KT) GO TO 260
      K = IWK(J+1) + I
      J = J + 1
  255 I = I + 1
      IWK(ILL+I) = J
      IF (I .LT. K) GO TO 255
      GO TO 230
  260 KB = K3
      IF (I .LE. 0) GO TO 265
      J = IWK(ILL+I)
      I = I - 1
      GO TO 230
  265 IF (KB .GE. N) GO TO 270
      J = 1
      GO TO 230
  270 JK = IWK(IC+KT)
      ISP = IWK(ID+KT)
      M = M - KT
      KB = ISP/JK-2
      IF (KT .GE. M-1 ) GO TO 9005
      ITA = ILL+KB+1
      ITB = ITA+JK
      IDM1 = ID-1
      IKT = KT+1
      IM = M+1
      DO 275 J = IKT,IM
         IWK(IDM1+J) = IWK(IDM1+J)/JK
  275 CONTINUE
      JJ = 0
      DO 290 J = 1,KB
         K = KT
  280    JJ = IWK(ID+K+1) + JJ
         IF (JJ .LT. IWK(ID+K)) GO TO 285
         JJ = JJ - IWK(ID+K)
         K = K + 1
         GO TO 280
  285    IWK(ILL+J) = JJ
         IF (JJ .EQ. J) IWK(ILL+J) = -J
  290 CONTINUE
C                                  DETERMINE THE PERMUTATION CYCLES
C                                  OF LENGTH GREATER THAN OR EQUAL
C                                  TO TWO.
      DO 300  J = 1,KB
         IF (IWK(ILL+J) .LE. 0) GO TO 300
         K2 = J
  295    K2 = IABS(IWK(ILL+K2))
         IF (K2 .EQ. J) GO TO 300
         IWK(ILL+K2) = -IWK(ILL+K2)
         GO TO 295
  300 CONTINUE
C                                  REORDER A FOLLOWING THE
C                                  PERMUTATION CYCLES
      I = 0
      J = 0
      KB = 0
      KN = N
  305 J = J + 1
      IF (IWK(ILL+J) .LT. 0) GO TO 305
      K = IWK(ILL+J)
      K0 = JK * K + KB
  310 ZA4 = A(K0+I+1)
      WK(ITA+I) = A4
      WK(ITB+I) = B4
      I = I + 1
      IF (I .LT. JK) GO TO 310
      I = 0
  315 K = -IWK(ILL+K)
      JJ = K0
      K0 = JK * K + KB
  320 A(JJ+I+1) = A(K0+I+1)
      I = I + 1
      IF (I .LT. JK) GO TO 320
      I = 0
      IF (K .NE. J) GO TO 315
  325 A(K0+I+1) = DCMPLX(WK(ITA+I),WK(ITB+I))
      I = I + 1
      IF (I .LT. JK) GO TO 325
      I = 0
      IF (J .LT. K2) GO TO 305
      J = 0
      KB = KB + ISP
      IF (KB .LT. KN) GO TO 305
 9005 RETURN
      END SUBROUTINE FFTCC

