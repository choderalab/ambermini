c-----------------------------------------------------------------------
      SUBROUTINE LIGENB(A,VEC,EIG,IA,N,NDIM)
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
C     IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- A GIVENS HOUSHOLDER MATRIX DIAGONALIZATION
C           ROUTINE SAME AS EIGEN BUT WORKS WITH A
C           LINEAR ARRAY -----
C
      DIMENSION A(*),VEC(NDIM,*),EIG(*),IA(*)
      DIMENSION W(50),GAMA(50),BETA(50),BETASQ(50)
      DIMENSION P(50),Q(50),IPOSV(50),IVPOS(50),IORD(50)
      EQUIVALENCE (P(1),Q(1)),(P(1),BETA(1)),(P(1),IVPOS(1))
      EQUIVALENCE (BETASQ(1),IORD(1)),(IPOSV(1),GAMA(1))
C
      DATA ZERO,PT5,ONE,TWO /0.0E+00,0.5E+00,1.0E+00,2.0E+00/
      DATA RHOSQ /1.0E-22/
C
      IF(N .EQ. 0) GO TO 900
      N1 = N-1
      N2 = N-2
      GAMA(1) = A(1)
      IF(N2) 360,340,100
  100 DO 320 NR = 1,N2
      IK = IA(NR+1)+NR
      B = A(IK)
      S = ZERO
      DO 120 I = NR,N2
      IJ = IA(I+2)+NR
  120 S = S+A(IJ)**2
C
C     ----- PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION -----
C
      A(IK) = ZERO
      IF(S .LE. ZERO) GO TO 300
      S = S+B*B
      SGN = +ONE
      IF(B .GE. ZERO) GO TO 140
      SGN = -ONE
  140 SQRTS = SQRT(S)
      D = SGN/(SQRTS+SQRTS)
      TEMP = SQRT(PT5+B*D)
      W(NR) = TEMP
      A(IK) = TEMP
      D = D/TEMP
      B = -SGN*SQRTS
C
C     ----- -D- IS FACTOR OF PROPORTIONALITY. NOW
C           COMPUTE AND SAVE -W- VECTOR. EXTRA SINGLY
C           SUBSCRIPTED -W- VECTOR FOR SPEED -----
C
      DO 160 I = NR,N2
      IJ = IA(I+2)+NR
      TEMP = D*A(IJ)
      W(I+1) = TEMP
  160 A(IJ) = TEMP
C
C     ----- PREMULTIPLY VECTOR -W- BY MATRIX -A- TO
C           OBTAIN -P- VECTOR. SIMULTANEOUSLY ACCUMULATE
C           DOT PRODUCT -WP- -- SCALR -K- -----
C
      WTAW = ZERO
      DO 240 I = NR,N1
      SUM = ZERO
      II = IA(I+1)
      DO 180 J = NR,I
      IJ = II+J+1
  180 SUM = SUM+A(IJ)*W(J)
      I1 = I+1
      IF(N1 .LT. I1) GO TO 220
      DO 200 J = I1,N1
      IJ = IA(J+1)+I+1
  200 SUM = SUM+A(IJ)*W(J)
  220 P(I) = SUM
  240 WTAW = WTAW+SUM*W(I)
      DO 260 I = NR,N1
  260 Q(I) = P(I)-WTAW*W(I)
C
C     ----- NOW FORM -PAP- MATRIX, REQUIRED PART -----
C
      DO 280 J = NR,N1
      QJ = Q(J)
      WJ = W(J)
      JJ = J+1
      DO 280 I = J,N1
      IJ = IA(I+1)+JJ
  280 A(IJ) = A(IJ)-TWO*(W(I)*QJ+WJ*Q(I))
  300 BETA(NR) = B
      BETASQ(NR) = B*B
      IL = IK+1
  320 GAMA(NR+1) = A(IL)
  340 IJ = IA(N)+N-1
      B = A(IJ)
      BETA(N-1) = B
      BETASQ(N-1) = B*B
      IJ = IJ+1
      GAMA(N) = A(IJ)
  360 BETASQ(N) = ZERO
C
C     ----- ADJOIN AN IDENTYTY MATRIX TO BE POST-
C           MULTIPLIED BY ROTATIONS -----
C
      DO 400 I = 1,N
      DO 380 J = 1,N
  380 VEC(I,J) = ZERO
  400 VEC(I,I) = ONE
      M = N
      SUM = ZERO
      NPAS = 1
      GO TO 600
  420 SUM = SUM+SHIFT
      COSA = ONE
      G = GAMA(1)-SHIFT
      PP = G
      PPBS = PP*PP+BETASQ(1)
      PPBR = SQRT(PPBS)
      DO 540 J = 1,M
      COSAP = COSA
      IF(PPBS .NE. ZERO) GO TO 440
      SINA = ZERO
      SINA2 = ZERO
      COSA = ONE
      GO TO 500
  440 SINA = BETA(J)/PPBR
      SINA2 = BETASQ(J)/PPBS
      COSA = PP/PPBR
C
C     ----- POSTMULTIPLY IDENTITY BY -P- TRANSPOSE -----
C
      NT = J+NPAS
      IF(NT .LT. N) GO TO 460
      NT = N
  460 CONTINUE
      DO 480 I = 1,NT
      TEMP = COSA*VEC(I,J)+SINA*VEC(I,J+1)
      VEC(I,J+1) = -SINA*VEC(I,J)+COSA*VEC(I,J+1)
  480 VEC(I,J) = TEMP
  500 DIA = GAMA(J+1)-SHIFT
      U = SINA2*(G+DIA)
      GAMA(J) = G+U
      G = DIA-U
      PP = DIA*COSA-SINA*COSAP*BETA(J)
      IF(J .NE. M) GO TO 520
      BETA(J) = SINA*PP
      BETASQ(J) = SINA2*PP*PP
      GO TO 560
  520 PPBS = PP*PP+BETASQ(J+1)
      PPBR = SQRT(PPBS)
      BETA(J) = SINA*PPBR
  540 BETASQ(J) = SINA2*PPBS
  560 GAMA(M+1) = G
C
C     ----- TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT -----
C
      NPAS = NPAS+1
      IF(BETASQ(M) .GT. RHOSQ) GO TO 620
  580 EIG(M+1) = GAMA(M+1)+SUM
  600 BETA(M) = ZERO
      BETASQ(M) = ZERO
      M = M-1
      IF(M .EQ. 0) GO TO 660
      IF(BETASQ(M) .LE. RHOSQ) GO TO 580
C
C     ----- TAKE ROOT OF CORMER 2 BY 2 NEAREST TO
C           LOWER DIAGONAL IN VALUE AS ESTIMATE OF
C           EIGENVALUE TO USE FOR SHIFT -----
C
  620 A2 = GAMA(M+1)
      R2 = PT5*A2
      R1 = PT5*GAMA(M)
      R12 = R1+R2
      DIF = R1-R2
      TEMP = SQRT(DIF*DIF+BETASQ(M))
      R1 = R12+TEMP
      R2 = R12-TEMP
      DIF = ABS(A2-R1)- ABS(A2-R2)
      IF(DIF .LT. ZERO) GO TO 640
      SHIFT = R2
      GO TO 420
  640 SHIFT = R1
      GO TO 420
  660 EIG(1) = GAMA(1)+SUM
      DO 680 J = 1,N
      IPOSV(J) = J
      IVPOS(J) = J
  680 IORD(J) = J
      M = N
      GO TO 740
  700 DO 720 J = 1,M
      IF(EIG(J) .LE. EIG(J+1)) GO TO 720
      TEMP = EIG(J)
      EIG(J) = EIG(J+1)
      EIG(J+1) = TEMP
      ITEMP = IORD(J)
      IORD(J) = IORD(J+1)
      IORD(J+1) = ITEMP
  720 CONTINUE
  740 M = M-1
      IF(M .NE. 0) GO TO 700
      IF(N1 .EQ. 0) GO TO 800
      DO 780 L = 1,N1
      NV = IORD(L)
      NP = IPOSV(NV)
      IF(NP .EQ. L) GO TO 780
      LV = IVPOS(L)
      IVPOS(NP) = LV
      IPOSV(LV) = NP
      DO 760 I = 1,N
      TEMP = VEC(I,L)
      VEC(I,L) = VEC(I,NP)
  760 VEC(I,NP) = TEMP
  780 CONTINUE
C
C     ----- BACK TRANSFORM THE VECTORS OF THE TRIPLE
C           DIAGONAL MATRIX -----
C
  800 DO 880 NRR = 1,N
      K = N1
  820 K = K-1
      IF(K .LE. 0) GO TO 880
      SUM = ZERO
      DO 840 I = K,N1
      IJ = IA(I+1)+K
  840 SUM = SUM+VEC(I+1,NRR)*A(IJ)
      SUM = SUM+SUM
      DO 860 I = K,N1
      IJ = IA(I+1)+K
  860 VEC(I+1,NRR) = VEC(I+1,NRR)-SUM*A(IJ)
      GO TO 820
  880 CONTINUE
  900 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
