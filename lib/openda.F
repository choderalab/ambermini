c-----------------------------------------------------------------------
      SUBROUTINE OPENDA(IREST,NAMF)
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
c
c     Rev A mod:  (G. Seibel) removed common/iofile/ (not used)
c                 cleaned up and made as portable as possible.
c                 This routine is machine dependent in that the method
c                 of specification of record length of a direct access
c                 file is not part of the ANSI Fortran standard.  Most
c                 byte-oriented machines want the length in bytes, except
c                 dec vms, which wants it in 'longwords', which are 4 byte
c                 words.  The 64bit Cray wants it in bytes, but remember
c                 that a Cray word is 8 bytes.
c
      CHARACTER*80 NAMF
#ifdef IBM_VM_CMS
      integer*4 irc
#endif
      COMMON/DAFIL/IDAF,NREC,NAV,IODA(510,3),HH(10240)
      COMMON/DAIOLN/IRECLN,IRECST,IFILEN(510)
c
#ifdef IBM_VM_CMS
      call vname(namf, 80)
#endif
      IDAF = 15
      IRECLN = 512
c
c     --- set up the "wordsize" ---
c
      iwrdsz = 4
#ifdef VMS
      iwrdsz = 1
#endif
#ifdef CRAYFISH
      iwrdsz = 8
#endif
#ifdef T3D
      iwrdsz = 8
#endif
      lenrec = irecln * iwrdsz
c
      IRECST = 1
      NREC = 0
      DO 10 I = 1, 510
           IODA(I,3) = -1
   10 CONTINUE
#ifdef IBM_VM_CMS
c
c     VM/CMS - need to set extent of file fairly large for
c     making db4.dat
c
      call fileinf(irc, 'maxrec', 1000)
      if (irc .ne. 0) then
         write(6,'(/2x,a,i6)') 'Error on fileinf: ',irc
         call mexit(6, 1)
      endif
#endif
c
      if (irest .lt. 1) then
c          --- creating new database (PREP only) ---
           IF (IREST .eq. 0) then
             OPEN(UNIT=IDAF, FILE=NAMF, STATUS='new', ACCESS='DIRECT',
     +          FORM='UNFORMATTED', recl=lenrec, iostat=ios)
           else if (irest .eq. -1) then
             OPEN(UNIT=IDAF, FILE=NAMF, STATUS='unknown',
     .          ACCESS='DIRECT', FORM='UNFORMATTED', recl=lenrec, 
     .          iostat=ios)
           endif
           if (ios .ne. 0) then
                write(6,'(/2x,a,a)') 'Error on open: ',namf
                call mexit(6, 1)
           endif
           IRECST = IRECST + 4
           WRITE(IDAF,REC=1) NREC,IRECST,IFILEN
           WRITE(IDAF,REC=2) (IODA(J,1),J=1,510)
           WRITE(IDAF,REC=3) (IODA(J,2),J=1,510)
           WRITE(IDAF,REC=4) (IODA(J,3),J=1,510)
      else 
c         --- reading & maybe appending ---
          if (irest.eq.1) then
c            --- old dbase - appending by prep ---
             OPEN(UNIT=IDAF, file=namf, STATUS='OLD', ACCESS='DIRECT',
     +          FORM='UNFORMATTED', recl=lenrec, iostat=ios)
             if (ios .ne. 0) then
                write(6,'(/2x,a,a)') 'Error on dbase open: ', namf
                call mexit(6, 1)
             endif
          else
c             --- link - read_only ---
#ifdef VMS
c             --- link uses unit # ---
              OPEN(UNIT=IDAF, STATUS='OLD', ACCESS='DIRECT',
     +          FORM='UNFORMATTED', recl=lenrec, iostat=ios,
     +          readonly)
              if (ios .ne. 0) then
                write(6,'(/2x,a,i6)') 'Error on open: unit ',idaf
                write(6,'(20x,a)') '(Unit must be assigned for VMS)'
                call mexit(6, 1)
              endif
#else
              OPEN(UNIT=IDAF, FILE=NAMF, STATUS='OLD', ACCESS='DIRECT',
     +          FORM='UNFORMATTED', recl=lenrec, iostat=ios)
              if (ios .ne. 0) then
                write(6,'(/2x,a,a)') 'Error on open: ',namf
                call mexit(6, 1)
              endif
#endif
           endif
           READ(IDAF,REC=1) NREC,IRECST,IFILEN
           READ(IDAF,REC=2) (IODA(J,1),J=1,510)
           READ(IDAF,REC=3) (IODA(J,2),J=1,510)
           READ(IDAF,REC=4) (IODA(J,3),J=1,510)
      endif
#ifdef IBM_VM_CMS
      call fileinf(irc)
      if (irc .ne. 0) then
         write(6,'(/2x,a,i6)') 'Error on fileinf reset: ',irc
         call mexit(6, 1)
      endif
#endif
      RETURN
      END
c-----------------------------------------------------------------------

