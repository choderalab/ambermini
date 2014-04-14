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
c-----------------------------------------------------------------------
      SUBROUTINE matinv(A,N,D,L,M)
      implicit double precision (a-h,o-z)
C
C     ----- STANDARD IBM MATRIX INVERSION ROUTINE -----
c
c     arguments:
c       a:  square matrix of dimension nxn
c       d:  resultant determinant
c       l:  work vector of length n
c       m:  work vector of length n
C
      DIMENSION A(*),L(*),M(*)
C
C     ----- SEARCH FOR LARGEST ELEMENT -----
C
      D = 1.0d0
      NK = -N
      DO 80 K = 1,N
      NK = NK+N
      L(K) = K
      M(K) = K
      KK = NK+K
      BIGA = A(KK)
      DO 20 J = K,N
      IZ = N*(J-1)
      DO 20 I = K,N
      IJ = IZ+I
      IF( ABS(BIGA)- ABS(A(IJ))) 15,20,20
   15 BIGA = A(IJ)
      L(K) = I
      M(K) = J
   20 CONTINUE
C
C     ----- INTERCHANGE ROWS -----
C
      J = L(K)
      IF(J-K) 35,35,25
   25 KI = K-N
      DO 30 I = 1,N
      KI = KI+N
      HOLD = -A(KI)
      JI = KI-K+J
      A(KI) = A(JI)
   30 A(JI) = HOLD
C
C     ----- INTERCHANGE COLUMNS -----
C
   35 I = M(K)
      IF(I-K) 45,45,38
   38 JP = N*(I-1)
      DO 40 J = 1,N
      JK = NK+J
      JI = JP+J
      HOLD = -A(JK)
      A(JK) = A(JI)
   40 A(JI) = HOLD
C
C     ----- DIVIDE COLUMN BY MINUS PIVOT -----
C
   45 IF(BIGA) 48,46,48
   46 D = 0.0d0
      GO TO 150
   48 DO 55 I = 1,N
      IF(I-K) 50,55,50
   50 IK = NK+I
      A(IK) = A(IK)/(-BIGA)
   55 CONTINUE
C
C     ----- REDUCE MATRIX -----
C
      DO 65 I = 1,N
      IK = NK+I
      HOLD = A(IK)
      IJ = I-N
      DO 65 J = 1,N
      IJ = IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ = IJ-I+K
      A(IJ) = HOLD*A(KJ)+A(IJ)
   65 CONTINUE
C
C     ----- DIVIDE ROW BY PIVOT -----
C
      KJ = K-N
      DO 75 J = 1,N
      KJ = KJ+N
      IF(J-K) 70,75,70
   70 A(KJ) = A(KJ)/BIGA
   75 CONTINUE
C
C     ----- PRODUCT OF PIVOTS -----
C
      D = D*BIGA
C
C     ----- REPLACE PIVOT BY RECIPROCAL -----
C
      A(KK) = 1.0d0/BIGA
   80 CONTINUE
C
C     ----- FINAL ROW AND COLUMN INTERCHANGE -----
C
      K = N
  100 K = (K-1)
      IF(K) 150,150,105
  105 I = L(K)
      IF(I-K) 120,120,108
  108 JQ = N*(K-1)
      JR = N*(I-1)
      DO 110 J = 1,N
      JK = JQ+J
      HOLD = A(JK)
      JI = JR+J
      A(JK) = -A(JI)
  110 A(JI) = HOLD
  120 J = M(K)
      IF(J-K) 100,100,125
  125 KI = K-N
      DO 130 I = 1,N
      KI = KI+N
      HOLD = A(KI)
      JI = KI-K+J
      A(KI) = -A(JI)
  130 A(JI) = HOLD
      GOTO 100
  150 RETURN
      END
