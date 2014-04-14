      subroutine amopen(lun,fname,fstat,fform,facc)
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
#ifdef VMS
      implicit none
#endif

c     INPUT:
c
      integer lun
c        ... logical unit number
      character*(*) fname
c        ... file name (not used in VAX/VMS implementation)
      character*1 fstat
c        ... status code: "N", "O", or "U" = new, old, unk.
      character*1 fform
c        ... format code: "U", "F" = unform., form.
      character*1 facc
c        ... access code: "R", "W", "A" = read, read/write, append
C The amopen subroutines are ifdefd for:
C
C       VMS
C       Unix/generic
C
#if VMS || IBM_VM_CMS
c-----------------------------------------------------------------------
c     THIS IS VAX/VMS and IBM VM/CMS VERSION
c     (grouped because both require external assignment of units)
c     Author: George Seibel
c
c
c     INTERNAL:
c
      character*7 stat
c        ... status keyword
      character*11 kform
c        ... form keyword
      integer ios
c        ... i/o status variable
c
#ifdef IBM_VM_CMS
c     VM/CMS: unit 6 is automatically opened so leave it alone
      if (lun.eq.6) return
c
#endif
      if (fstat .eq. 'N') then
           stat = 'NEW'
      elseif (fstat .eq. 'O') then
           stat = 'OLD'
      elseif (fstat .eq. 'U') then
           stat = 'UNKNOWN'
      else
           write(6,'(/,2x,a,i4)')
     +                'amopen: bogus fstat, unit ', lun
           call mexit(6, 1)
      endif
c
      if (fform .eq. 'U') then
           kform = 'UNFORMATTED'
      elseif (fform .eq. 'F') then
           kform = 'FORMATTED'
      else
           write(6,'(/,2x,a,i4)')
     +                'amopen: bogus fform, unit ', lun
           call mexit(6, 1)
      endif
c
#ifdef VMS
      if (facc .eq. 'R') then
           if (fform .eq. 'U') then
                open(unit=lun,status=stat,form=kform,readonly,
     +               iostat=ios)
           elseif (fform .eq. 'F') then
                open(unit=lun,status=stat,carriagecontrol='list',
     +               form=kform,readonly,iostat=ios)
           endif
      elseif (facc .eq. 'W' .or. facc .eq. 'A') then
           if (fform .eq. 'U') then
                open(unit=lun,status=stat,form=kform,iostat=ios)
           elseif (fform .eq. 'F') then
                open(unit=lun,status=stat,
     +               carriagecontrol='list',form=kform,iostat=ios)
           endif
      else
           write(6,'(/,2x,a,i4)')
     +                'amopen: bogus facc, unit ', lun
           call mexit(6, 1)
      endif
#else
      open(unit=lun,status=stat,form=kform,iostat=ios)
#endif
c
      if (ios .ne. 0) then
#ifdef IBM_VM_CMS
           if (lun.ne.6) close(unit=lun)
#else
           close(unit=lun)
#endif
           if (lun .eq. 6)  then
c               this is the only place outside of mexit() where a stop
c               stmt should occur
                stop 'Error on Unit 6 OPEN'
           else
                write(6,'(/,2x,a,i4)')
     +                'Error on OPEN, Fortran unit ', lun
                call mexit(6, 1)
           endif
      endif
      return
      end

#else
c-----------------------------------------------------------------------
c     THIS IS UNIX VERSION
c     Author: George Seibel
c     Rev 13-Jun-90:  add rewind after open.
 
c     INTERNAL:
 
      character*7 stat
c        ... status keyword
      character*11 kform
c        ... form keyword
      integer ios
c        ... i/o status variable
 
      if (fstat .eq. 'N') then
           stat = 'NEW'
      elseif (fstat .eq. 'O') then
           stat = 'OLD'
      elseif (fstat .eq. 'U') then
           stat = 'UNKNOWN'
      else
           write(6,'(/,2x,a,i4)')
     $           'amopen: bogus fstat, unit ', lun
           call mexit(6, 1)
      endif
c
      if (fform .eq. 'U') then
           kform = 'UNFORMATTED'
      elseif (fform .eq. 'F') then
           kform = 'FORMATTED'
      else
           write(6,'(/,2x,a,i4)')
     $           'amopen: bogus fform, unit', lun
           call mexit(6, 1)
      endif
c
      open(unit=lun,file=fname,status=stat,form=kform,iostat=ios)
c
      if (ios .ne. 0) then
           if (lun .eq. 6) then
#ifndef DUMB
                write(0,'(/,2x,a,i4,a,a)') 'Unit ', lun, 
     $                                    ' Error on OPEN: ',fname
#endif
           else
#ifndef DUMB
                write(0,'(/,2x,a,i4,a,a)') 'Unit ', lun, 
     $                 ' Error on OPEN: ',fname(1:len_trim(fname))
#endif
                write(6,'(/,2x,a,i4,a,a)') 'Unit ', lun, 
     $                 ' Error on OPEN: ',fname(1:len_trim(fname))
                close(unit=6)
           endif
           call mexit(6, 1)
      endif
      rewind(lun)
      return
      end
#endif
