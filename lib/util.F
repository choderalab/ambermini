c
c    util.f - non-precision-dependent, simple routines
c
c-----------------------------------------------------------------------
      SUBROUTINE IWIND(NF)
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
      REWIND NF
      RETURN
      END
      SUBROUTINE CLOSC(NF)
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
#ifdef IBM_VM_CMS
c     vm/cms do not allow close of unit 6
      if (nf.eq.6) return
#endif
      CLOSE(UNIT=NF)
      RETURN
      END
c-----------------------------------------------------------------------
