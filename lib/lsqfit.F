      SUBROUTINE LSQFIT(ISTEP,FIT,FITI,FIT2,LOUT)
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
      implicit double precision (a-h,o-z)
c     Rev 17-May-89: hardwired the use of brief output.
      LOGICAL LOUT
C
C     ----- ROUTINE TO MAKE A LEAST SQUARE FIT TO A STRAIGHT LINE ----
C
c  ---dac change 2/90: protect against divide by zero:
c
      if (istep.lt.3) return
c
      FIT = FIT/ISTEP
      FITI = FITI/ISTEP
      FIT2 = FIT2/ISTEP
      IM1 = ISTEP-1
      IP1 = ISTEP+1
      IPM1 = IP1*IM1
      ALINE = (2.d0*FITI-IP1*FIT)*6.d0/IPM1
      BLINE = FIT-ALINE*IP1/2.d0
      RLINE = 0.0d0
      DUM = FIT2-FIT*FIT
      IF(ABS(DUM).GT.1.0d-06) RLINE = IPM1/(DUM*12.0d0)
      IF(RLINE.LT.0.d0) RLINE = 0.d0
      RLINE = ALINE* SQRT(RLINE)
      CHILIN = FIT2-2.d0*ALINE*FITI-2.d0*BLINE*FIT+BLINE**2+
     +             ALINE**2*IP1*(IP1+IM1+1)/6.d0+ALINE*BLINE*IP1
      IF(CHILIN.LT.0.d0) CHILIN = 0.d0
      CHILIN =  SQRT(CHILIN)
      if (lout) WRITE(6,9118) ISTEP,ALINE,BLINE
 9118 FORMAT(/5X,' RESULT OF LEAST SQUARE FIT OVER',I5,' STEPS',
     +       /5X,' ENERGY DRIFT PER STEP =',F14.6,5X,'ETOT(AT X=0) =',
     +       1PE11.3)
      RETURN
      END
