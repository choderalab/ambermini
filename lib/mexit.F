      subroutine mexit(output_unit, status)
C
C  mexit() - machine-dependent exit() procedure, designed to return an 
C            appropriate (success/failure) value to the operating system.
c            This is non-mpi version; parallel codes should create their
c            own equivalents!
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
C                Regents of the University of California               **
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
      implicit none
      integer output_unit  ! close this unit if greater than zero
      integer status       ! exit status; error if non-zero


      if (output_unit.gt.0) then
          close(unit=output_unit)
      endif

#if XLF90 || IBM3090 || F2C
      if (status.ne.0) then
          stop 1
      else
          stop 0
      endif
#else
      call exit(status)
#endif
      end



