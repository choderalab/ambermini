c-----------------------------------------------------------------------
      SUBROUTINE RGROUP(NATOM,NATC,NRES,NGRP,IPRES,LBRES,IGRAPH,ISYMBL,
     +                  ITREE,IGROUP,WEIT,XC,KONST,BELLY,nfu)
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
      !implicit double precision (a-h,o-z)
      implicit none
c
c     Mods for Rev A by gls:
c     - cpp selectable precision
c     - wildcard functionality fixed. (see findp)
c     - format change to print less garbage
c     - unreferenced statement labels removed
c     - changed title write from 12 to 19 A4. (not 20 to avoid wraps)
c     - changed to read from unit nf -- dac 12/90
c
      logical :: konst, misc, belly
C
C     ----- READ IN GROUPS EACH ATOM IS IN -----
C
C           WILL KEEP READING IN GROUPS OF CARDS UNTIL A BLANK CARD IS
C           READ
C           ALL ATOMS IN THIS CARD GROUP WILL BE PUT IN THE SAME GROUP
C           THE ONLY EXCEPTION IS A RES CARD WHICH STARTS WITH A -
C           RES NO. THE RESIDUES IN THIS GROUP WILL ALL BE PUT IN INDIV.
C           GROUPS ANY OTHER GROUPS IN THIS SECTION MUST ALSO START
C           WITH A  I.E. THEY MUST ALSO BE INDIV. GROUPS
C           ANY TIME A "FIND" CARD IS READ, ONLY THOSE ATOMS THAT MEET
C           THE SPECIFICATIONS ON AT LEAST 1 OF THE FIND CARDS
C           WILL BE INCLUDED IN THE GROUP
C           THIS OPTION MAY BE ENDED BY READING IN "NFIND",TERMINATING
C           THE INPUT FOR THE GROUPS, OR BY READING IN ANOTHER "FIND"
C
      character(len=4) jgraph,jresnm,jsymbl,jtree
      integer :: isrch
      COMMON/PROPF/JGRAPH(11),JRESNM(11),JSYMBL(11),JTREE(11),ISRCH
      character(len=4) ititl1(20)
      character(len=4) katn, ifind, nfind, iiend, ksear, ihol(20), ktypg
      character(len=4) igraph(*),isymbl(*),itree(*),lbres(*)
      
      integer :: i, i1, i2, iend, ifld(20), igroup(*), igrp(8), iii,
     +           istart, ipres(*), itime, ivar(20), izero,
     +           j, jfld(20), jgrp(8),
     +           k,
     +           lfind, lsign,
     +           natc, natmg, natom, nf, nfu, ngrp, nres
      
      double precision :: fvar(20),
     +                    xc(*),
     +                    weit(*), wt

      DATA IFLD/ 2*0, 13*2 , 5*0 /
      DATA JFLD/ 4*1, 16*0 /
      DATA KATN /'ATOM'/
      DATA IFIND /'FIND'/
      DATA NFIND /'NFIN'/
      DATA IIEND /'END '/
      DATA KSEAR /'SEAR'/
C
C     ----- RES CARD LISTS RESIDUE GROUPS IGRP(I) TO JGRP(I) ----
C
C           IF 1ST VALUE IS NEGATIVE, THEN EVERY RESIDUE FROM IGRP(I)
C           TO JGRP(I) IS CONSIDERED A SEPARATE GROUP
C           IF 2ND VALUE   = 0, THEN SUBGROUP CONSISTS OF 1 RESIDUE
C           ATOM CARD READS IN GROUPS OF ATOMS
C           IF 2ND ATOM IN EACH PAIR  = 0 , THEN SUBGROUP CONSISTS OF
C           1 ATOM RES AND ATOM CARDS MAY BE READ IN ANY ORDER
C           END INPUT WITH AN "END " CARD
C
C           ZERO NGRP BEFORE CALLING THIS ROUTINE
C           ROUTINE WILL RETURN WITH THE NUMBER OF THE LAST GROUP READ
C
      NGRP = 0
      nf = nfu
      if (nf.le.0) nf=5
C
C     ----- INITIALISE THE GROUP ARRAY -----
C
      DO 100 I = 1,NATOM
        IGROUP(I) = 0
        IF (KONST) WEIT(I) = 0.0d0
  100 CONTINUE
   22 CONTINUE
        ITIME = 0
        LSIGN = 0
        IZERO = 0
        ISRCH = 0
        NATMG = 0
C
C       ----- READ DIFFERENT GROUPS -----
C
        READ(nf,9208) (ITITL1(K),K=1,19)
        IF(ITITL1(1).EQ.IIEND) then
          write(6, '(4x,a)') '----- END OF GROUP READ -----'
          GO TO 900
        endif
        NGRP = NGRP+1
        write(6,'(4x,a,i5,a)') '----- READING GROUP ',
     +                       ngrp, '; TITLE:'
        WRITE(6,9218) (ITITL1(K),K=1,19)
C
C       ----- IF CONSTRAINED GROUPS READ THE WEIGHT FOR EACH GROUP -----
C
        IF (KONST) then
          IFLD(1) = 3
          IFLD(2) = 0
          CALL RFREE(IFLD,IHOL,IVAR,FVAR,nf,6)
          WT = FVAR(1)
          WRITE(6,9018) NGRP,WT
        endif
   10   CONTINUE
C
C       ----- READ THE GROUP CARDS -----
C
        IFLD(1) = 1
        IFLD(2) = 2
        CALL RFREE(IFLD,IHOL,IVAR,FVAR,nf,6)
        KTYPG = IHOL(1)
        K = 1
        DO 120 I = 1,7
          IGRP(I) = IVAR(K)
          JGRP(I) = IVAR(K+1)
          K = K+2
  120   CONTINUE
        IF (KTYPG.EQ.IIEND) GO TO 16
        IF (KTYPG.EQ.NFIND) then
          WRITE(6,199)
          ISRCH = 0
          GO TO 10
        endif
        IF (KTYPG.eq.IFIND) then
C
C         ----- FIND OPTION ... READ THE ATOM SPECIFICS -----
C
          WRITE(6,200)
          DO 64 III = 1,12
            CALL RFREE(JFLD,IHOL,IVAR,FVAR,nf,6)
            JGRAPH(III) = IHOL(1)
            JSYMBL(III) = IHOL(2)
            JTREE(III) = IHOL(3)
            JRESNM(III) = IHOL(4)
            IF(JGRAPH(III).EQ.KSEAR) GO TO 65
            WRITE(6,202) JGRAPH(III),JSYMBL(III),JTREE(III),JRESNM(III)
   64     CONTINUE
   65     CONTINUE
C
          ISRCH = III-1
          IF(ISRCH.GT.10) WRITE(6,66) ISRCH
C
C         ----- NOW READ IN RES AND ATOMS TO BE SEARCHED -----
C
          GO TO 10
        endif
        ITIME = ITIME+1
        IF (KTYPG.ne.KATN) then
C
C         ----- 1ST GROUP OF 1ST CARD MUST BE - IF ANY - NUMBERS ARE
C               FOUND -----
C
          IF(ITIME.EQ.1.AND.IGRP(1).LT.0) LSIGN = 1
          DO 12 I = 1,7
            I1 = IGRP(I)
            IF (I1.EQ.0) GO TO 10
            I2 = JGRP(I)
            IF(I2.GT.NRES) I2 = NRES
            IF(I1.GT.0.AND.LSIGN.EQ.1) GO TO 36
            IF(I1.LT.0.AND.LSIGN.EQ.0) GO TO 36
            IF(I1.LT.0) I1 = -I1
            IF(I2.LE.0) I2 = I1
            IF(LSIGN.EQ.0) WRITE(6,14) NGRP,I1,I2
            DO 13 J = I1,I2
              ISTART = IPRES(J)
              IEND = IPRES(J+1)-1
              DO 45 K = ISTART,IEND
                IF (ISRCH.GT.0) 
     +            CALL FINDP(K,J,LFIND,NRES,IPRES,LBRES,ISYMBL,
     +                       ITREE,IGRAPH)
                IF (ISRCH.GT.0.AND.LFIND.EQ.0) GO TO 45
                IGROUP(K) = NGRP
                IF(KONST) WEIT(K) = WT
                NATMG = NATMG+1
   45         CONTINUE
              IF(LSIGN.EQ.1) WRITE(6,46) NGRP,J
              IF(LSIGN.EQ.1) NGRP = NGRP+1
   13       CONTINUE
   12     CONTINUE
          GO TO 10
        endif
C
C       ----- ATOM TYPE CONSTRAINTS -----
C
        IF(LSIGN.EQ.1) GO TO 36
        WRITE(6,51) NGRP
        DO 17 I = 1,7
          I1 = IGRP(I)
          IF(I1.LT.0) GO TO 36
          IF(I1.EQ.0) GO TO 10
          I2 = JGRP(I)
          IF(I2.GT.NATOM) I2 = NATOM
          IF(I2.LE.0) I2 = I1
          WRITE(6,52) I1,I2
          DO 18 J = I1,I2
            IF(ISRCH.GT.0) 
     +        CALL FINDP(J,IZERO,LFIND,NRES,IPRES,LBRES,ISYMBL,
     +                   ITREE,IGRAPH)
            IF(ISRCH.GT.0.AND.LFIND.EQ.0) GO TO 18
            NATMG = NATMG+1
            IGROUP(J) = NGRP
            IF(KONST) WEIT(J) = WT
   18     CONTINUE
   17   CONTINUE
        GO TO 10
C
   16   IF(LSIGN.EQ.1) NGRP = NGRP-1
        IF(ITIME.EQ.0) NGRP = NGRP-1
c       IF(ISRCH.GT.0) WRITE(6,199)
        WRITE(6,222) NATMG
      GO TO 22
C
   36 CONTINUE
      WRITE(6,127) KTYPG,(IGRP(I),JGRP(I),I = 1,7)
      GO TO 10
C
C     ----- ALL GROUPS ARE READ RETURN -----
C
  900 CONTINUE
      IF (KONST) then
C
C       ----- GATHER ALL THE CONSTRAINED ATOMS TOGETHER -----
C
        NATC = 0
        DO 920 I = 1,NATOM
          IF(IGROUP(I).LE.0) GO TO 920
          NATC = NATC+1
          IGROUP(NATC) = I
          WEIT(NATC) = WEIT(I)
  920   CONTINUE
C
C       ----- do not PRINT THE HISTORY OF CONSTRAINED ATOMS -----
C
#ifdef debug
        WRITE(6,9108)
        DO 940 I = 1,NATC
          J = IGROUP(I)
          J3 = 3*J-3
          WRITE(6,9118) I,J,IGRAPH(J),ISYMBL(J),ITREE(J),
     +                  (XC(J3+K),K=1,3), WEIT(I)
  940   CONTINUE
#endif
      else if (.not.BELLY) then
C
C       ----- PUT THE ATOMS WHICH ARE NOT IN THE DEFINED GROUPS
C             AS THE LAST GROUP -----
C
        MISC = .FALSE.
        DO 820 I = 1,NATOM
          IF(IGROUP(I).NE.0) GO TO 820
          MISC = .TRUE.
          IGROUP(I) = NGRP+1
  820   CONTINUE
        IF(MISC) NGRP = NGRP+1
c       IF(MISC) WRITE(6,9308) NGRP
      endif
*  11 FORMAT(A4,1X,7(2I5))
  199 FORMAT(1H ,5X,'END OF ATOM SPECIFICATION',/)
  200 FORMAT(1H ,5X,'ALL ATOMS THAT MEET 1 OF THE FOLLOWING',
     +    ' SPECIFICATIONS WILL BE INCLUDED IN GROUP BELOW',/)
*  61 FORMAT(A4,A2,A1,A4)
  202 FORMAT(1H ,5X,'GRAPH NAME  = ',A4,2X,'SYMBOL  = ',A2,4X,
     +        'TREE SYMBOL  = ',A1,5X,'RESIDUE TYPE  = ',A4,/)
   66 FORMAT(1H ,5X,'**** NUMBER OF FIND CARDS  = ',I5,2X,
     +        'IS TOO BIG ******',/)
   14 FORMAT(' GRP',I5,' RES',I5,' TO ',I5)
   46 FORMAT(1H ,5X,'GROUP',I5,3X,'CONSISTS OF RESIDUE',I5,/)
   51 FORMAT(1H ,5X,'GROUP',I5,3X,'CONSISTS OF ATOMS -',/)
   52 FORMAT(1H ,34X,I5,2X,'TO',I5)
  222 FORMAT(1H ,5X,'Number of atoms in this group  = ',I5)
  127 FORMAT(1H ,5X,'***PROBLEMS WITH GROUP',A4,14I5,'*******',/)
 9018 FORMAT(/5X,'GROUP',I5,' HAS HARMONIC CONSTRAINTS',F12.5)
 9118 FORMAT(I5,I6,1X,A4,1X,2A4,3F10.4,F12.5)
 9108 FORMAT(/ /10X,'HISTORY OF CONSTRAINED ATOMS',/ /)
 9208 FORMAT(20A4)
 9218 FORMAT(1X,20A4)
c9308 FORMAT(/5X,'THE GROUP ',I4, ' CONTAINS ALL ATOMS NOT DEFINED',
c    +       ' AS GROUPS BY THE INPUT',/)
      RETURN
      END SUBROUTINE RGROUP
c-----------------------------------------------------------------------
      SUBROUTINE FINDP(IATOM,IRES,ILOC,NRES,IPRES,LBRES,ISYMBL,ITREE,
     +                 IGRAPH)
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
      !implicit double precision (a-h,o-z)
      implicit none
      
c     Rev A mod (G. Seibel, Apr 89)
c     Changed iblank, jblank to iwild, jwild, to give the wildcard
c     functionality promised by the doc.
c     isymbl() dimensioned. (this was long-standing bug)
C
      character(len=4) jgraph,jresnm,jsymbl,jtree
      integer :: isrch
      COMMON/PROPF/JGRAPH(11),JRESNM(11),JSYMBL(11),JTREE(11),ISRCH
      
      character(len=4), intent(in) :: igraph(*), isymbl(*), itree(*),
     +                                lbres(*)
      character(len=4) :: iwild, jwild

      integer, intent(in) :: iatom, ipres(*),
     +                       nres
      integer, intent(inout) :: iloc, ires
      
      integer :: i
C
C     ----- CHECKS IF A GIVEN ATOM HAS CERTAIN CHARACTERISTICS -----
C
      DATA iwild /'*   '/
      DATA jwild /'*   '/
C
      ILOC = 0
      IF(IRES.EQ.0) CALL FINDRS(IATOM,IRES,NRES,IPRES)
      DO 10 I = 1,ISRCH
        IF((JRESNM(I).NE.iwild).AND.(JRESNM(I).NE.LBRES(IRES)))
     +                GO TO 10
        IF((JGRAPH(I).NE.iwild).AND.(JGRAPH(I).NE.IGRAPH(IATOM)))
     +                GO TO 10
        IF((JTREE(I).NE.jwild).AND.(JTREE(I).NE.ITREE(IATOM)))
     +                GO TO 10
        IF((JSYMBL(I).NE.jwild).AND.(JSYMBL(I).NE.ISYMBL(IATOM)))
     +                GO TO 10
        ILOC = 1
        GO TO 20
   10 CONTINUE
   20 CONTINUE
      RETURN
      END SUBROUTINE FINDP
c-----------------------------------------------------------------------
      subroutine findrs(numa,ires,nres,ipres)
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
      implicit none
      
      integer, intent(in) :: ipres(*),
     +                       nres, numa
      integer, intent(inout) :: ires
      
      integer :: i, im
      
      if (numa .gt. 0) then
         if (numa .ge. ipres(nres)) then
            ires = nres
            return
         else
            im = nres - 1
            do i=1,im
               if (numa .ge. ipres(i) .and. numa .lt. ipres(i+1)) then
                  ires = i
                  return
               end if
            end do
         end if
      end if
      
      write(6,100) numa
      call mexit(6, 1)
      
  100 format(/2x,'PROBLEMS FINDING RESIDUE OF ATOM ',i5)
      
      end subroutine findrs
