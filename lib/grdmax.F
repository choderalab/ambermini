      subroutine grdmax(n,g,iatmax,fdmax)
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
      implicit double precision (a-h,o-z)
c     Rev A mods:  converted this routine from function to subrt.
c                  added iatmax = atom number of max gradient to args.
      DIMENSION G(*)
      DUM = 0.0D0
      iatmax = 1
      DO 100 I = 1,N
           GI = ABS(G(I))
           IF (GI.GT.DUM) then
                DUM = GI
                iatmax = i
           endif
  100 CONTINUE
      fdmax = DUM
      RETURN
      END
