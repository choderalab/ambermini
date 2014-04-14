/*
 *  File:   tripos.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *  Description:
 *      Code to read TRIPOS files from SYBYL.
 */

/*        Modifications 
 *        Mason Louchart  (2011)
 *        Universite de Picardie - Jules Verne, Amiens
 *        http://q4md-forcefieldtools.org
 *
 *        mol3 file format
 */

#include    "basics.h"
#include    "classes.h"
#include    "tripos.h"
#include    "parmLib.h"

/*
 *  uTriposReadUnit
 *
 *  Author: Christian Schafmeister (1991)
 *
 *  Read a single UNIT from the TRIPOS file.
 *  Return the UNIT, or NULL if nothing was read.
 */
UNIT uTriposReadUnit(fIn)
FILE *fIn;
{
    STRING sLine;
    int iAtoms;
    VARARRAY vaAtoms;
    VARARRAY vaResidues;
    UNIT uUnit;
    STRING sUnitName;
    int iRet;
    double dX, dY, dZ;//, dDummy;
    STRING sType;
    int iResidue;
    STRING sTemp;
    VECTOR vPos;
    int iIndex, iA, iB;
    ATOM aA, aB;
    int iOrder;
    ATOM aAtom;
    BOOL bGotIt;
    int iFileAtoms = 0;
    int iFileBonds = 0;
    int iFileSubstructures = 0;
    int iTemp, i;
    double dCharge;
    RESIDUE rRes;
    STRING sName;
    STRING sOrder;
    int iRootAtom;
    //STRING sSize;
    //STRING sDescriptor, sDesc;
    STRING sCmd;
    int iElement;//, iTag, iDummy;
    //PARMSET psSet;
    STRING sMol3="none";

#define TRIPOS_MOLECULE     "@<TRIPOS>MOLECULE"
#define TRIPOS_ATOM     "@<TRIPOS>ATOM"
#define TRIPOS_BOND     "@<TRIPOS>BOND"
#define TRIPOS_SUBSTRUCTURE "@<TRIPOS>SUBSTRUCTURE"
#define TRIPOS_HEADTAIL  "@<TRIPOS>HEADTAIL\n"
#define TRIPOS_RESIDUECONNECT  "@<TRIPOS>RESIDUECONNECT\n"

/*
fgets a line, skipping whitespace only lines and converting control
characters to spaces; then scan that line; goto fail on eof or error.
see below for usage examples.
*/
#define T_FSCANF( iret,f,s,ss,fail) { \
              char *r ;\
              int  isallwhite = 1 ;\
              while(1) { ;\
                  r = fgets( s, sizeof(s), fIn ) ;\
                  if (r == NULL) goto fail ;\
                  for( ; *r != '\0' ; ++r ) { ;\
                      if ( isspace(*r) ) continue ;\
                      if ( isprint(*r) ) { isallwhite=0; continue; } ;\
                      if ( iscntrl(*r) ) { ;\
                          *r = ' ' ;\
                          VP0 ( ( "Warning: control character read; " ) ) ;\
                          VP0 ( ( "it was converted to a space!\n" ) ) ;\
                      } ;\
                  } ;\
                  if ( ! isallwhite ) break ;\
              } ;\
              iret=sscanf ss ;\
              if (feof(f)) goto fail; }


    vaAtoms = NULL;
    vaResidues = NULL;

    /* Search for the TRIPOS_MOLECULE string */

    bGotIt = FALSE;
    while (TRUE) {
        T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
        if (strcmp(sCmd, TRIPOS_MOLECULE) == 0) {
            bGotIt = TRUE;
            break;
        }
    }
    if (!bGotIt)
        goto FAIL;

    /* Read the UNITs name */

    uUnit = (UNIT) oCreate(UNITid);
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sUnitName), FAIL);
    VP0 ( ( "Reading MOLECULE named %s", sLine ) ) ;
    ContainerSetName(uUnit, sUnitName);
    MESSAGE(("Reading unit: %s\n", sUnitName));
    T_FSCANF(iRet, fIn, sLine, (sLine, "%d %d %d %d %d", &iFileAtoms,
                                &iFileBonds, &iFileSubstructures,
                                &iTemp, &iTemp), FAIL);
    MESSAGE(("Number of atoms: %d\n", iFileAtoms));
    MESSAGE(("Number of bonds: %d\n", iFileBonds));
    MESSAGE(("Number of substructures: %d\n", iFileSubstructures));

    /* Search for the TRIPOS_ATOM string */

    bGotIt = FALSE;
    while (TRUE) {
        T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
        if (strcmp(sCmd, TRIPOS_ATOM) == 0) {
            bGotIt = TRUE;
            break;
        }
    }
    if (!bGotIt)
        goto FAIL;

#if 0
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sSize), FAIL);
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sDescriptor), FAIL);
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sTemp), FAIL);
    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sTemp), FAIL);

    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
    if (strcmp(sCmd, TRIPOS_ATOM) != 0)
        goto FAIL;
#endif

    /* Read in the ATOM records */

    vaAtoms = vaVarArrayCreate(sizeof(ATOM));
    for (i = 0; i < iFileAtoms; i++) {
        dCharge = 0.0;
        T_FSCANF(iRet, fIn, sLine,
                 (sLine, "%d %s %lf %lf %lf %s %d %s %lf",
                  &iIndex, sName,
                  &dX, &dY, &dZ, sType, &iResidue, sTemp, &dCharge), FAIL);
        if (iRet != 9)
            break;

        /* don't allow the residue number to be greater than no. of
           substructures, nor less than 1:   */
        iResidue = iResidue > iFileSubstructures ? iFileSubstructures : iResidue;
        iResidue = iResidue < 1 ? 1 : iResidue;

        MESSAGE((" Atom: %s\n", sName));
        aAtom = (ATOM) oCreate(ATOMid);
        ContainerSetName(aAtom, sName);
        AtomSetTempInt(aAtom, iResidue);
        VectorDef(&vPos, dX, dY, dZ);
        AtomSetPosition(aAtom, vPos);

/*  iElement = iElementNumberFromAmber(sType);   */
/*    dac change: just use the first letter: how else to easily tell
      CA from calcium???                                              */

        sTemp[0] = cUpper(sType[0]);
        sTemp[1] = '\0';
        iElement = iElementNumberFromAmber(sTemp);

#if 0
/*  initial attempt to use the atom map, but how do we get that read in?  */
        if (bParmLibDefaultExists()) {
            fprintf(stderr, "have the library\n");
            PARMLIB_DEFAULT_LOOP(psSet,
                                 (iTag = iParmSetFindAtom(psSet, sType)));
            if (iTag != PARM_NOT_FOUND) {
                ParmSetAtom(psSet, iTag, sTemp,
                            &dDummy, &dDummy, &dDummy, &dDummy, &dDummy,
                            &dDummy, &iElement, &iDummy, sDesc);
            }
        }
#endif

        AtomSetElement(aAtom, iElement);
        AtomSetType(aAtom, sType);
        AtomSetCharge(aAtom, dCharge);
        VarArrayAdd(vaAtoms, (GENP) & aAtom);
    }

    T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
    if (strcmp(sCmd, TRIPOS_BOND) != 0)
        goto FAIL;

    iAtoms = iVarArrayElementCount(vaAtoms);
    for (i = 0; i < iFileBonds; i++) {
        T_FSCANF(iRet, fIn, sLine,
                 (sLine, "%d %d %d %s", &iIndex, &iA, &iB, sOrder), FAIL);
        if (iRet != 4)
            break;
        MESSAGE((" Bond %d - %d\n", iA, iB));
        if (iA > iAtoms || iB > iAtoms) {
            printf("Cannot form bond between atoms %d and %d\n", iA, iB);
        } else {
            aA = *PVAI(vaAtoms, ATOM, iA - 1);
            aB = *PVAI(vaAtoms, ATOM, iB - 1);
            StringLower(sOrder);
            iOrder = BONDSINGLE;
            if (strcmp(sOrder, "1") == 0) {
                iOrder = BONDSINGLE;
            } else if (strcmp(sOrder, "2") == 0) {
                iOrder = BONDDOUBLE;
            } else if (strcmp(sOrder, "3") == 0) {
                iOrder = BONDTRIPLE;
            } else if (strcmp(sOrder, "ar") == 0) {
                iOrder = BONDAROMATIC;
            } else if (strcmp(sOrder, "am") == 0) {
                iOrder = BONDAROMATIC;
            }
            AtomBondToOrder(aA, aB, iOrder);
        }
    }

    /* Read the RESIDUE stuff */

    if( iFileSubstructures ){
        T_FSCANF(iRet, fIn, sLine, (sLine, "%s", sCmd), FAIL);
        if (strcmp(sCmd, TRIPOS_SUBSTRUCTURE) != 0) goto FAIL;

        vaResidues = vaVarArrayCreate(sizeof(RESIDUE));
        for (i = 0; i < iFileSubstructures; i++) {
            T_FSCANF(iRet, fIn, sLine,
                 (sLine, "%d %s %d %s",
                  &iIndex, sName, &iRootAtom, sType), FAIL);
            if (iRet != 4) break;

            MESSAGE((" Substructure: %s\n", sName));
            rRes = (RESIDUE) oCreate(RESIDUEid);
            ContainerSetName(rRes, sName);
            ContainerSetSequence(rRes, iIndex);
            VarArrayAdd(vaResidues, (GENP) & rRes);
        }
    } else {

        /* create a single residue with the same name as the unit name:  */

        vaResidues = vaVarArrayCreate(sizeof(RESIDUE));
        rRes = (RESIDUE) oCreate(RESIDUEid);
        ContainerSetName(rRes, sUnitName);
        ContainerSetSequence(rRes, iIndex);
        VarArrayAdd(vaResidues, (GENP) & rRes);
    }

    /* Connect everything together */
    MESSAGE((" Building the UNIT\n" ));

    /* Put the ATOMs within the RESIDUES */

    for (i = 0; i < iVarArrayElementCount(vaAtoms); i++) {
        aAtom = *PVAI(vaAtoms, ATOM, i);
        iResidue = iAtomTempInt(aAtom);
        rRes = *PVAI(vaResidues, RESIDUE, iResidue - 1);
        ContainerAdd((CONTAINER) rRes, (OBJEKT) aAtom);
    }

/*_____________________________________________________________________________
|                                                                              |
|       Author: Mason Louchart (2011)                                          |
|       http://q4md-forcefieldtools.org                                        |
|       Universite de Picardie - Jules Verne, Amiens                           |
|                                                                              |
|       Tutorial available at                                                  |
|       http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php                 |
|                                                                              |
|       If the TRIPOS_HEADTAIL tag is read, the rest of the information is     |
|       read (mol3File: head and tail)                                         |
|_____________________________________________________________________________*/

    fgets(sMol3,25,fIn);
    if ( strcmp(sMol3,TRIPOS_HEADTAIL) == 0 ){
        int LineNumber,AtomNumber;
        int HeadTailParentNumber,ResSequenceNumber;
        STRING HeadTailName;   		

        /* Set the Head and Tail */
        for (LineNumber=1; LineNumber<=2; LineNumber++){
            fscanf(fIn, "%s %d", HeadTailName, &HeadTailParentNumber);
            ATOM aShortLived = (ATOM) oCreate(ATOMid);
            for (AtomNumber = 0; AtomNumber < iVarArrayElementCount(vaAtoms); AtomNumber++){
                aShortLived = *PVAI(vaAtoms, ATOM, AtomNumber);
                ResSequenceNumber = iContainerSequence(cContainerWithin(aShortLived));
                if ( ResSequenceNumber == HeadTailParentNumber ){
                    if ( strcmp(aShortLived->cHeader.sName, HeadTailName) == 0 ){
                        if (LineNumber == 1){
                            UnitSetHead(uUnit, aShortLived);
                        }
                        else {
                            UnitSetTail(uUnit, aShortLived);
                        }
                    }
                }
            }
        }
        fseek(fIn, 1, SEEK_CUR);
    }

/*_____________________________________________________________________________
|                                                                              |
|       Author: Mason Louchart (2011)                                          |
|       http://q4md-forcefieldtools.org                                        |
|       Universite de Picardie - Jules Verne, Amiens                           |
|                                                                              |
|       Tutorial available at                                                  |
|       http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php                 |
|                                                                              |
|       If the TRIPOS_RESIDUECONNECT tag is read, the rest of the information  |
|       is read (mol3File: residue connects)                                   |
|_____________________________________________________________________________*/

    if ( strcmp(sMol3,TRIPOS_HEADTAIL) == 0 ){
        fgets(sMol3,25,fIn);
        if ( strcmp(sMol3,TRIPOS_RESIDUECONNECT) == 0 ){
            STRING ConnectName[6];
            int j,AtomNumber,ResNumber,ResSequenceNumber;

            /* Set residues connects */
            for (i=0; i < iVarArrayElementCount(vaResidues); i++){
                fscanf(fIn, "%d %s %s %s %s %s %s",&ResNumber, ConnectName[0], ConnectName[1],
                    ConnectName[2], ConnectName[3], ConnectName[4], ConnectName[5]);
                RESIDUE rRes = *PVAI(vaResidues, RESIDUE, i);
                ATOM aShortLived = (ATOM) oCreate(ATOMid);
                for (j=0; j<=5; j++){
                    if(ConnectName[j]!=0){
                        for (AtomNumber=0; AtomNumber < iVarArrayElementCount(vaAtoms); AtomNumber++){
                            aShortLived = *PVAI(vaAtoms, ATOM, AtomNumber);
                            ResSequenceNumber = iContainerSequence(cContainerWithin(aShortLived));
                            if ( ResSequenceNumber==ResNumber ){
                                if ( strcmp(aShortLived->cHeader.sName, ConnectName[j]) == 0 ){
                                    ResidueSetConnectAtom(rRes,j,aShortLived);	
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /* Put the RESIDUES within the UNIT */

    for (i = 0; i < iVarArrayElementCount(vaResidues); i++) {
        rRes = *PVAI(vaResidues, RESIDUE, i);
        ContainerAdd((CONTAINER) uUnit, (OBJEKT) rRes);
    }
    VarArrayDestroy(&vaAtoms);
    VarArrayDestroy(&vaResidues);

    return (uUnit);

  FAIL:
    VP0 ( ( "Fatal Error:  last line read: %s", sLine ) ) ;
    /* Destroy everything */
    if (vaAtoms) {
        for (i = 0; i < iVarArrayElementCount(vaAtoms); i++) {
            aAtom = *PVAI(vaAtoms, ATOM, i);
            Destroy((OBJEKT *) & aAtom);
        }
        VarArrayDestroy(&vaAtoms);
    }
    if (vaResidues) {
        for (i = 0; i < iVarArrayElementCount(vaResidues); i++) {
            rRes = *PVAI(vaResidues, RESIDUE, i);
            Destroy((OBJEKT *) & rRes);
        }
        VarArrayDestroy(&vaResidues);
    }
    return (NULL);
}
