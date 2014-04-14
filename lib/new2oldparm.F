C Program NEW2OLDPARM
C
C This program converts a new format PARM file, with headers before
C each section, to an old format PARM file without headers (and
C possibly with only a subset of the data). 

C  Works as a filter:   new2oldparm < prmtop-new > prmtop-old
C
C Author: David Pearlman
C Date:   9/00
C
      PARAMETER (MXPRM = 1000000)

      CHARACTER*80 INFIL,OUTFIL,FMT
      CHARACTER*80 TYPE

      DIMENSION IV(MXPRM),IW(MXPRM),IX(MXPRM),IY(MXPRM),IZ(MXPRM)
      double precision XX(MXPRM)
  
C

      CALL NXTSEC(5,0,1,' ','TITLE',FMT,IOK)

C If IOK = -1, then this is not a 'new' format parm file. Stop with error.

      IF (IOK.EQ.-1) THEN
         WRITE(6,10)
   10    FORMAT('ERROR: Input file is not a new format PARM file; ',/,
     *          '       No %VERSION line found')
         STOP
      END IF

      READ(5,FMT) (IX(I),I=1,20)
      WRITE(6,'(20A4)') (IX(I),I=1,20)

      TYPE = 'POINTERS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ(5,FMT)        NATOM,NTYPES, NBONH, MBONA,NTHETH,MTHETA, 
     *                   NPHIH, MPHIA,NHPARM, NPARM,   NNB,  NRES,
     *                   NBONA,NTHETA, NPHIA,NUMBND,NUMANG, NPTRA,
     *                   NATYP,  NPHB,IFPERT, NBPER, NGPER, NDPER,
     *                   MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP,
     *                   NUMEXTRA
      WRITE(6,'(12I6)')  NATOM,NTYPES, NBONH, MBONA,NTHETH,MTHETA, 
     *                   NPHIH, MPHIA,NHPARM, NPARM,   NNB,  NRES,
     *                   NBONA,NTHETA, NPHIA,NUMBND,NUMANG, NPTRA,
     *                   NATYP,  NPHB,IFPERT, NBPER, NGPER, NDPER,
     *                   MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP,
     *                   NUMEXTRA

C Make sure scratch arrays are big enough to hold required data

      IIMAX = MAX(NATOM,NTYPES*NTYPES,NATYP,NBONH,NBONA,NTHETH)
      iimax = max(iimax,NTHETA,NPHIH,NPHIA,NNB,2*NBPER,2*NGPER,2*NDPER)
      IF (MXPRM.LT.IIMAX) GO TO 9001

      TYPE = 'ATOM_NAME'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NATOM)
      WRITE(6,'(20A4)') (IX(I),I=1,NATOM)

      TYPE = 'CHARGE'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NATOM)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NATOM)

      TYPE = 'MASS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NATOM)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NATOM)

      TYPE = 'ATOM_TYPE_INDEX'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NATOM)
      WRITE(6,'(12I6)') (IX(I),I=1,NATOM)

      TYPE = 'NUMBER_EXCLUDED_ATOMS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NATOM)
      WRITE(6,'(12I6)') (IX(I),I=1,NATOM)

      TYPE = 'NONBONDED_PARM_INDEX'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NTYPES*NTYPES)
      WRITE(6,'(12I6)') (IX(I),I=1,NTYPES*NTYPES)

      TYPE = 'RESIDUE_LABEL'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NRES)
      WRITE(6,'(20A4)') (IX(I),I=1,NRES)

      TYPE = 'RESIDUE_POINTER'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NRES)
      WRITE(6,'(12I6)') (IX(I),I=1,NRES)

      TYPE = 'BOND_FORCE_CONSTANT'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NUMBND)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NUMBND)

      TYPE = 'BOND_EQUIL_VALUE'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NUMBND)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NUMBND)

      TYPE = 'ANGLE_FORCE_CONSTANT'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NUMANG)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NUMANG)

      TYPE = 'ANGLE_EQUIL_VALUE'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NUMANG)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NUMANG)

      TYPE = 'DIHEDRAL_FORCE_CONSTANT'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NPTRA) 
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NPTRA)

      TYPE = 'DIHEDRAL_PERIODICITY'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NPTRA) 
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NPTRA) 

      TYPE = 'DIHEDRAL_PHASE'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NPTRA) 
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NPTRA) 

      TYPE = 'SOLTY'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NATYP)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NATYP) 

      TYPE = 'LENNARD_JONES_ACOEF'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NTYPES*(NTYPES+1)/2)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NTYPES*(NTYPES+1)/2) 

      TYPE = 'LENNARD_JONES_BCOEF'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NTYPES*(NTYPES+1)/2)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NTYPES*(NTYPES+1)/2) 

      TYPE = 'BONDS_INC_HYDROGEN'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),IY(I),IZ(I),I=1,NBONH)
      WRITE(6,'(12I6)') (IX(I),IY(I),IZ(I),I=1,NBONH)

      TYPE = 'BONDS_WITHOUT_HYDROGEN'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),IY(I),IZ(I),I=1,NBONA)
      WRITE(6,'(12I6)') (IX(I),IY(I),IZ(I),I=1,NBONA)

      TYPE = 'ANGLES_INC_HYDROGEN'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),IY(I),IZ(I),IV(I),I=1,NTHETH)
      WRITE(6,'(12I6)') (IX(I),IY(I),IZ(I),IV(I),I=1,NTHETH)

      TYPE = 'ANGLES_WITHOUT_HYDROGEN'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),IY(I),IZ(I),IV(I),I=1,NTHETA)
      WRITE(6,'(12I6)') (IX(I),IY(I),IZ(I),IV(I),I=1,NTHETA)

      TYPE = 'DIHEDRALS_INC_HYDROGEN'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),IY(I),IZ(I),IV(I),IW(I),I=1,NPHIH)
      WRITE(6,'(12I6)') (IX(I),IY(I),IZ(I),IV(I),IW(I),I=1,NPHIH)

      TYPE = 'DIHEDRALS_WITHOUT_HYDROGEN'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),IY(I),IZ(I),IV(I),IW(I),I=1,NPHIA)
      WRITE(6,'(12I6)') (IX(I),IY(I),IZ(I),IV(I),IW(I),I=1,NPHIA)

      TYPE = 'EXCLUDED_ATOMS_LIST'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NNB)
      WRITE(6,'(12I6)') (IX(I),I=1,NNB)

      TYPE = 'HBOND_ACOEF'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NPHB)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NPHB)

      TYPE = 'HBOND_BCOEF'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NPHB)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NPHB)

      TYPE = 'HBCUT'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (XX(I),I=1,NPHB)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NPHB)

      TYPE = 'AMBER_ATOM_TYPE'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NATOM)
      WRITE(6,'(20A4)') (IX(I),I=1,NATOM)

      TYPE = 'TREE_CHAIN_CLASSIFICATION'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NATOM)
      WRITE(6,'(20A4)') (IX(I),I=1,NATOM)

      TYPE = 'JOIN_ARRAY'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NATOM)
      WRITE(6,'(12I6)') (IX(I),I=1,NATOM)

      TYPE = 'IROTAT'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NATOM)
      WRITE(6,'(12I6)') (IX(I),I=1,NATOM)

      IF (IFBOX.LE.0) GO TO 50

      TYPE = 'SOLVENT_POINTERS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      IPTRES,  NSPM,NSPSOL
      WRITE(6,'(12I6)') IPTRES,  NSPM,NSPSOL

      TYPE = 'ATOMS_PER_MOLECULE'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      (IX(I),I=1,NSPM)
      WRITE(6,'(12I6)') (IX(I),I=1,NSPM)

      TYPE = 'BOX_DIMENSIONS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      BETA,BOX1,BOX2,BOX3
      WRITE(6,'(1P5E16.8)') BETA,BOX1,BOX2,BOX3

   50 IF (IFCAP.LE.0) GO TO 100

      TYPE = 'CAP_INFO'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT)      NATCAP,CUTCAP,XCAP,YCAP,ZCAP
      WRITE(6,'(I6,/,1P5E16.8)') NATCAP,CUTCAP,XCAP,YCAP,ZCAP

  100 IF (IFPERT.LE.0) GO TO 150

      TYPE = 'PERT_BOND_ATOMS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),IY(I),I=1,NBPER)
      WRITE(6,'(12I6)') (IX(I),IY(I),I=1,NBPER)

      TYPE = 'PERT_BOND_PARAMS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),I=1,2*NBPER)
      WRITE(6,'(12I6)') (IX(I),I=1,2*NBPER)

      TYPE = 'PERT_ANGLE_ATOMS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),IY(I),IZ(I),I=1,NGPER)
      WRITE(6,'(12I6)') (IX(I),IY(I),IZ(I),I=1,NGPER)

      TYPE = 'PERT_ANGLE_PARAMS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),I=1,2*NGPER)
      WRITE(6,'(12I6)') (IX(I),I=1,2*NGPER)

      TYPE = 'PERT_DIHEDRAL_ATOMS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),IY(I),IZ(I),IV(I),I=1,NDPER)
      WRITE(6,'(12I6)') (IX(I),IY(I),IZ(I),IV(I),I=1,NDPER)

      TYPE = 'PERT_DIHEDRAL_PARAMS'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),I=1,2*NDPER)
      WRITE(6,'(12I6)') (IX(I),I=1,2*NDPER)

      TYPE = 'PERT_RESIDUE_NAME'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),I=1,NRES)
      WRITE(6,'(20A4)') (IX(I),I=1,NRES)

      TYPE = 'PERT_ATOM_NAME'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),I=1,NATOM)
      WRITE(6,'(20A4)') (IX(I),I=1,NATOM)

      TYPE = 'PERT_ATOM_SYMBOL'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),I=1,NATOM)
      WRITE(6,'(20A4)') (IX(I),I=1,NATOM)

      TYPE = 'ALMPER'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (XX(I),I=1,NATOM)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NATOM)

      TYPE = 'IAPER'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),I=1,NATOM)
      WRITE(6,'(12I6)') (IX(I),I=1,NATOM)

      TYPE = 'PERT_ATOM_TYPE_INDEX'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (IX(I),I=1,NATOM)
      WRITE(6,'(12I6)') (IX(I),I=1,NATOM)

      TYPE = 'PERT_CHARGE'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (XX(I),I=1,NATOM)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NATOM)
c
C Following includes polarizability info. This will be included if
C it exists in the PARM file we are translating, or not (with info
C message to user) if it's not there. Since polarizability info
C appeared at end of standard old parm file, it won't hurt if we
C don't include it...
c
  150 CONTINUE

      TYPE = 'POLARIZABILITY'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) THEN
         WRITE(0,151)
  151    FORMAT('INFO: No polarizability info found in "New" file')
         GO TO 200
      END IF
      READ (5,FMT) (XX(I),I=1,NATOM)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NATOM)

      IF (IFPERT.LE.0) GO TO 200

      TYPE = 'PERT_POLARIZABILITY'
      CALL NXTSEC(5,0,1,' ',TYPE,FMT,IOK)
      IF (IOK.EQ.-2) GO TO 9000
      READ (5,FMT) (XX(I),I=1,NATOM)
      WRITE(6,'(1P5E16.8)') (XX(I),I=1,NATOM)
c
C All done:
c
  200 CONTINUE
      GO TO 1000
 
      
C Errors:

 9000 WRITE(0,9500)  TYPE
 9500 FORMAT('ERROR: Required data section not found. Data title:',/,
     *       A)
      STOP

 9001 WRITE(0,9501)  
 9501 FORMAT('ERROR: Scratch array not large enough. Increase MXPRM')
      STOP
 1000 END
