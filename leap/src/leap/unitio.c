/*
 *        File:        unitio.newparm.c
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
 *        Description:
 *                Input/output routines for UNITs
 *                this has been seperated out to make unit.c smaller
 *                and because there is ALOT of complex code here
 *                that doesn't have alot to do with the day to day operation
 *                of UNITs.
 *
 *
 *      Class:
 *              UNIT
 *      Superclass:
 *              CONTAINER, LOOP
 *
 *      Description:
 *              A UNIT is a subclass of CONTAINER.
 *              UNITS can contain molecules, residues, and atoms.
 *              UNITS are used to contain the entire molecular
 *              system plus other information about the system.
 *
 *        NOTE:        OFF files are used to write UNITs to files.
 *                In OFF files there is no implicit ordering
 *                of ANY DATA WHATSOEVER.  The PROGRAMMER is
 *                to assume that there is no order regardless
 *                of the ordering that the code in this file
 *                generates.
 *
 *        NOTE2:        When setting up tables for writing to OFF files,
 *                Indices are FORTRAN indices, where the first
 *                element has index=1.
 *
 */

/*      Modifications induced by the implementation of the savemol2 command
*       Christine Cezard (2007) 
*       Universite de Picardie - Jules Verne
*       http://q4md-forcefieldtools.org
*       zbUnitIOIndexBondParameters and zUnitDoAtoms are now "extern functions" 
*/ 

#include <time.h>

#include        "basics.h"
#include        "vector.h"
#include        "classes.h"
#include        "restraint.h"
#include        "bag.h"
#include        "dictionary.h"
#include        "database.h"
#include        "parmLib.h"
#include        "avl.h"
#include        "defaults.h"
#include        "tools.h"
#include        "fortran.h"
#include        "mathop.h"
#include        "sort.h"
#include        "cmap.h"
#ifdef BINTRAJ
#  include      "netcdf.h"
#endif

int iFatal;

/*
 *        Private data types
 *
 */


#define PERTURBED       0x00000001
#define BOUNDARY        0x00000002
#define JSBFAC          0.89089872 /* = 1/(2^(1/6)), needed for Jayaram et al. (M)GB radii */

typedef struct {
    CONTAINERNAMEt sName;
    CONTAINERNAMEt sPertName;
    ATOMTYPEt sType;
    ATOMTYPEt sPertType;
    int iTypeIndex;
    int iPertTypeIndex;
    int iElement;
    int iPertElement;
    double dCharge;
    double dPertCharge;
    int iResidueIndex;
    VECTOR vPos;
    VECTOR vVelocity;
    int iSequence;
    FLAGS fFlags;
    ATOM aAtom;
} SAVEATOMt;

typedef struct {
    int iType;
    FLAGS fFlags;
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iAtom4;
    double dKx;
    double dX0;
    double dN;
    int iParmIndex;                /* This is filled in when */
    /* the parameters are written */
    /* to the file */
} SAVERESTRAINTt;

typedef struct {
    int iAtom1;
    int iAtom2;
    int iParmIndex;
    int iPertParmIndex;
    FLAGS fFlags;
} SAVEBONDt;

typedef struct {
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iParmIndex;
    int iPertParmIndex;
    FLAGS fFlags;
} SAVEANGLEt;

typedef struct {
    BOOL bProper;
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iAtom4;
    int iParmIndex;
    BOOL bCalc14;
    int iPertParmIndex;
    BOOL bPertCalc14;
    FLAGS fFlags;
} SAVETORSIONt;


typedef struct {
    int iAtom1;
    int iAtom2;
    FLAGS fFlags;
} SAVECONNECTIVITYt;


typedef struct {
    STRING sName;
    int iSequenceNumber;
    int iNextChildSequence;
    MOLECULE mMolecule;
} SAVEMOLECULEt;

typedef struct {
    char sAboveType[2];
    int iAboveIndex;
    char sBelowType[2];
    int iBelowIndex;
} SAVEHIERARCHYt;

typedef struct {
    int iConnect;
} SAVECONNECTt;


typedef struct {
    BOOL bCapableOfHBonding;
    double dE;
    double dR;
    double dE14;
    double dR14;
    typeStr sType;
} NONBONDt;

typedef struct {
  typeStr         sType1;
  typeStr         sType2;
  double          dEI;
  double          dEJ;
  double          dRI;
  double          dRJ;
  double          dA;
  double          dC;
  DESCRIPTION     sDesc;
} NBEDITt;

typedef struct {
    double dA;
    double dC;
    double dA14;
    double dC14;
} NONBONDACt;

                        /* Used to save the bounding box info */
                        /* All of this is stored in one OFF entry */
                        /* dUseBox is greater than zero if the bounding */
                        /* box is used, and not greater than zero if it is */
                        /* not */

typedef struct {
    double dUseBox;
    double dBeta;
    double dXWidth;
    double dYWidth;
    double dZWidth;
} SAVEBOXt;

typedef struct {
    double dUseCap;
    double dX;
    double dY;
    double dZ;
    double dRadius;
} SAVECAPt;


typedef struct {
    int iGroupIndex;
    int iIndexAtom;
} SAVEGROUPSt;


/*
static void debugtypes(UNIT uUnit, char *str)
{
    int i;
    ATOM aAtom;

    fprintf(stderr, "--------------------- %s\n", str);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
        fprintf(stderr, "\t%d\t%s\t%s\n",
                iAtomId(aAtom), sAtomName(aAtom), sAtomType(aAtom));
    }
    fprintf(stderr, "--------------------------------\n");
}
*/

static void SetTreeType(ATOM aAtom, int iCoordLeft)
{
    switch (iCoordLeft) {
    case 0:
        AtomSetTempDouble(aAtom, (double) 'E');
        break;
    case 1:
        AtomSetTempDouble(aAtom, (double) 'S');
        break;
    case 2:
        AtomSetTempDouble(aAtom, (double) 'B');
        break;
    case 3:
        AtomSetTempDouble(aAtom, (double) '3');
        break;
    case 4:
        VP0(("  (%s: unexpected coordination 4 - %s)\n",
             sAtomName(aAtom), "using nonstandard tree type '4'"));
        AtomSetTempDouble(aAtom, (double) '4');
        break;
    case 5:
        AtomSetTempDouble(aAtom, (double) '5');
        VP0(("  (%s: unexpected coordination 5 - %s)\n",
             sAtomName(aAtom), "using nonstandard tree type '5'"));
        break;
    case 6:
        AtomSetTempDouble(aAtom, (double) '6');
        VP0(("  (%s: unexpected coordination 6 - %s)\n",
             sAtomName(aAtom), "using nonstandard tree type '6'"));
        break;
    default:
        AtomSetTempDouble(aAtom, (double) 'X');
        VP0(("  (%s: unexpected coordination %d - %s)\n",
             sAtomName(aAtom), iCoordLeft, "using tree type 'X'"));
        break;
    }
}

/*
 *      zUnitLoadTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Load the tables which can be used to construct the UNIT
 *      from a DATABASE.
 */
BOOL zbUnitIOLoadTables(UNIT uUnit, DATABASE db)
{
    int iSize, iCount, iAtomCount, iType;
    STRING sName;
    SAVEATOMt *saPAtom;
    SAVECONNECTIVITYt *scPConnectivity;
    SAVERESTRAINTt *srPRestraint;
    SAVERESIDUEt *srPResidue, *srPResTemp;
    SAVEMOLECULEt *smPMolecule;
    SAVEHIERARCHYt *shPHierarchy;
    SAVEGROUPSt *sgPGroupAtom;
    int iBondCount, iRestraintCount, iSequence;
    BOOL bGotOne;
    SAVEBOXt sbBox;
    SAVECAPt scCap;
    int iTemp, iLen, i;


    DBPushPrefix(db, "unit.");

    if (!bDBGetValue(db, "name", &iLen, (GENP) sName, sizeof(sName))) {
        bGotOne = FALSE;
        goto DONE;
    }
    ContainerSetName(uUnit, sName);

    bDBGetValue(db, "childsequence", &iLen, (GENP) & iSequence, 0);
    ContainerSetNextChildsSequence(uUnit, iSequence);

    /* Construct an array of atoms with names and types */

    uUnit->vaAtoms = vaVarArrayCreate(sizeof(SAVEATOMt));
    bDBGetType(db, "atoms", &iType, &iAtomCount);
    VarArraySetSize((uUnit->vaAtoms), iAtomCount);
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    iSize = sizeof(SAVEATOMt);
    bDBGetTable(db, "atoms", &iAtomCount,
                3, (char *) &(saPAtom->iTypeIndex), iSize,
                4, (char *) &(saPAtom->iResidueIndex), iSize,
                5, (char *) &(saPAtom->fFlags), iSize,
                6, (char *) &(saPAtom->iSequence), iSize,
                7, (char *) &(saPAtom->iElement), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                8, (char *) &(saPAtom->dCharge), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) saPAtom->sName, iSize,
                2, (char *) saPAtom->sType, iSize,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    bDBGetTable(db, "atomspertinfo", &iAtomCount,
                3, (char *) &(saPAtom->iPertTypeIndex), iSize,
                4, (char *) &(saPAtom->iPertElement), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                5, (char *) &(saPAtom->dPertCharge), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) saPAtom->sPertName, iSize,
                2, (char *) saPAtom->sPertType, iSize,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Get the atom positions */

    bDBGetTable(db, "positions", &iAtomCount,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) &dVX(&(saPAtom->vPos)), sizeof(SAVEATOMt),
                2, (char *) &dVY(&(saPAtom->vPos)), sizeof(SAVEATOMt),
                3, (char *) &dVZ(&(saPAtom->vPos)), sizeof(SAVEATOMt),
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Get the atoms velocities */

    bDBGetTable(db, "velocities", &iAtomCount,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) &dVX(&(saPAtom->vVelocity)), sizeof(SAVEATOMt),
                2, (char *) &dVY(&(saPAtom->vVelocity)), sizeof(SAVEATOMt),
                3, (char *) &dVZ(&(saPAtom->vVelocity)), sizeof(SAVEATOMt),
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Load the bounding box information */

    bDBGetValue(db, "boundbox", &iLen, (GENP) & sbBox,
                sizeof(sbBox.dUseBox));
    UnitSetUseBox(uUnit, (sbBox.dUseBox > 0.0));
    UnitSetBeta(uUnit, sbBox.dBeta);
    UnitSetBox(uUnit, sbBox.dXWidth, sbBox.dYWidth, sbBox.dZWidth);
    ToolSanityCheckBox(uUnit);

    /* Load the cap information */

    bDBGetValue(db, "solventcap", &iLen, (GENP) & scCap,
                sizeof(scCap.dUseCap));
    UnitSetUseSolventCap(uUnit, (scCap.dUseCap > 0.0));
    UnitSetSolventCap(uUnit, scCap.dX, scCap.dY, scCap.dZ, scCap.dRadius);

    /* Load the connectivity information */

    uUnit->vaConnectivity = vaVarArrayCreate(sizeof(SAVECONNECTIVITYt));
    bDBGetType(db, "connectivity", &iType, &iBondCount);
    VarArraySetSize((uUnit->vaConnectivity), iBondCount);
    if (iBondCount) {
        scPConnectivity =
            PVAI(uUnit->vaConnectivity, SAVECONNECTIVITYt, 0);
        iSize = sizeof(SAVECONNECTIVITYt);
        bDBGetTable(db, "connectivity", &iBondCount,
                    1, (char *) &(scPConnectivity->iAtom1), iSize,
                    2, (char *) &(scPConnectivity->iAtom2), iSize,
                    3, (char *) &(scPConnectivity->fFlags), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

        /* Load the restraints information */
    }
    uUnit->vaRestraints = vaVarArrayCreate(sizeof(SAVERESTRAINTt));
    bDBGetType(db, "restraints", &iType, &iRestraintCount);
    VarArraySetSize((uUnit->vaRestraints), iRestraintCount);
    if (iRestraintCount) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        iSize = sizeof(SAVERESTRAINTt);
        bDBGetTable(db, "restraints", &iBondCount,
                    1, (char *) &(srPRestraint->iType), iSize,
                    2, (char *) &(srPRestraint->fFlags), iSize,
                    3, (char *) &(srPRestraint->iAtom1), iSize,
                    4, (char *) &(srPRestraint->iAtom2), iSize,
                    5, (char *) &(srPRestraint->iAtom3), iSize,
                    6, (char *) &(srPRestraint->iAtom4), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    7, (char *) &(srPRestraint->dKx), iSize,
                    8, (char *) &(srPRestraint->dX0), iSize,
                    9, (char *) &(srPRestraint->dN), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);
    }
    /* Load the UNIT connect atoms */

    uUnit->vaConnect = vaVarArrayCreate(sizeof(int));
    VarArraySetSize((uUnit->vaConnect), 2);
    bDBGetValue(db, "connect", &iLen,
                (GENP) PVAI(uUnit->vaConnect, int, 0), sizeof(int));

    /* Load an array of residues */

    uUnit->vaResidues = vaVarArrayCreate(sizeof(SAVERESIDUEt));
    bDBGetType(db, "residues", &iType, &iCount);
    VarArraySetSize((uUnit->vaResidues), iCount);
    srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
    iSize = sizeof(SAVERESIDUEt);
    if (MAXCONNECT != 6)
        DFATAL(("MAXCONNECT has been changed, update UnitLoadTables"));
    bDBGetTable(db, "residues", &iCount,
                2, (char *) &(srPResidue->iSequenceNumber), iSize,
                3, (char *) &(srPResidue->iNextChildSequence), iSize,
                4, (char *) &(srPResidue->iAtomStartIndex), iSize,
                6, (char *) &(srPResidue->iImagingAtomIndex), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) srPResidue->sName, iSize,
                5, (char *) srPResidue->sResidueType, iSize,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Load the RESIDUE PDB sequence numbers from an */
    /* extra table to remain compatible with old OFF files */

    if (bDBGetType(db, "residuesPdbSequenceNumber", &iType, &iTemp)) {
        bDBGetValue(db, "residuesPdbSequenceNumber", &iTemp,
                    (GENP) & (srPResidue->iPdbSequenceNumber), iSize);
    } else {

        /* If no sequence numbers are defined then create some */

        srPResTemp = srPResidue;
        for (i = 0; i < iCount; i++) {
            srPResTemp->iPdbSequenceNumber = i + 1;
            srPResTemp++;
        }
    }

    bDBGetTable(db, "residueconnect", &iCount,
                1, (char *) &(srPResidue->iaConnectIndex[0]), iSize,
                2, (char *) &(srPResidue->iaConnectIndex[1]), iSize,
                3, (char *) &(srPResidue->iaConnectIndex[2]), iSize,
                4, (char *) &(srPResidue->iaConnectIndex[3]), iSize,
                5, (char *) &(srPResidue->iaConnectIndex[4]), iSize,
                6, (char *) &(srPResidue->iaConnectIndex[5]), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Construct an array of molecules */

    uUnit->vaMolecules = vaVarArrayCreate(sizeof(SAVEMOLECULEt));
    bDBGetType(db, "molecules", &iType, &iCount);
    VarArraySetSize((uUnit->vaMolecules), iCount);
    if (iCount) {
        smPMolecule = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        iSize = sizeof(SAVEMOLECULEt);
        bDBGetTable(db, "molecules", &iCount,
                    2, (char *) &(smPMolecule->iSequenceNumber), iSize,
                    3, (char *) &(smPMolecule->iNextChildSequence), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    1, (char *) smPMolecule->sName, iSize,
                    0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);
    }
    /* Construct an array which describes the hierarchy */

    uUnit->vaHierarchy = vaVarArrayCreate(sizeof(SAVEHIERARCHYt));
    bDBGetType(db, "hierarchy", &iType, &iCount);
    VarArraySetSize((uUnit->vaHierarchy), iCount);
    if (iCount) {
        shPHierarchy = PVAI(uUnit->vaHierarchy, SAVEHIERARCHYt, 0);
        iSize = sizeof(SAVEHIERARCHYt);
        bDBGetTable(db, "hierarchy", &iCount,
                    2, (char *) &(shPHierarchy->iAboveIndex), iSize,
                    4, (char *) &(shPHierarchy->iBelowIndex), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    1, (char *) shPHierarchy->sAboveType, iSize,
                    3, (char *) shPHierarchy->sBelowType, iSize,
                    0, NULL, 0, 0, NULL, 0, 0, NULL, 0);
    }
    /* Load the atom groups */

    if (bDBGetType(db, "groupNames", &iType, &iCount)) {
        uUnit->vaGroupNames = vaVarArrayCreate(sizeof(STRING));
        uUnit->vaGroupAtoms = vaVarArrayCreate(sizeof(SAVEGROUPSt));
        VarArraySetSize((uUnit->vaGroupNames), iCount);
        if (iCount) {
            bDBGetValue(db, "groupNames", &iCount,
                        (GENP) PVAI(uUnit->vaGroupNames, STRING, 0),
                        sizeof(STRING));
            if (bDBGetType(db, "groupAtoms", &iType, &iCount)) {
                VarArraySetSize(uUnit->vaGroupAtoms, iCount);
                if (iCount) {
                    sgPGroupAtom = PVAI(uUnit->vaGroupAtoms, SAVEGROUPSt,
                                        0);
                    iSize = sizeof(SAVEGROUPSt);
                    bDBGetTable(db, "groupAtoms", &iCount,
                                1, (char *) &(sgPGroupAtom->iGroupIndex),
                                iSize, 2,
                                (char *) &(sgPGroupAtom->iIndexAtom),
                                iSize, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0,
                                0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0,
                                NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL,
                                0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0,
                                NULL, 0, 0, NULL, 0);
                }
            }
        }
    }

    /* 
     *  DELETED: uUnit->psParameters = psParmSetLoad( db ); 
     *  since saved params confuse the issue of precedence.
     */

    bGotOne = TRUE;

  DONE:
    DBPopPrefix(db);
    return (bGotOne);
}

/*
 *      zUnitIOSaveTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Save the tables which can be used to construct the UNIT
 *      from a DATABASE.
 */
void zUnitIOSaveTables(UNIT uUnit, DATABASE db)
{
    int iSize, iCount, iSequence;
    SAVEATOMt *saPAtom;
    SAVEBONDt *sbPBond;
    SAVECONNECTIVITYt *scPConnectivity;
    SAVEANGLEt *saPAngle;
    SAVETORSIONt *stPTorsion;
    SAVERESTRAINTt *srPRestraint;
    SAVERESIDUEt *srPResidue;
    SAVEMOLECULEt *smPMolecule;
    SAVEHIERARCHYt *shPHierarchy;
    SAVEGROUPSt *sgPGroupAtom;
    SAVEBOXt sbBox;
    SAVECAPt scCap;


    if (!iVarArrayElementCount(uUnit->vaAtoms)) {
        VP0(("\tUnit has no atoms!\n"));
        return;
    }
    DBPushPrefix(db, "unit.");

    DBPutValue(db, "name", ENTRYSINGLE | ENTRYSTRING, 1,
               (GENP) sContainerName(uUnit), 0);

    iSequence = iContainerNextChildsSequence(uUnit);
    DBPutValue(db, "childsequence", ENTRYSINGLE | ENTRYINTEGER, 1,
               (GENP) & iSequence, 0);

    /* Save the array of atoms with names and types */

    iCount = iVarArrayElementCount(uUnit->vaAtoms);
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    iSize = sizeof(SAVEATOMt);
    DBPutTable(db, "atoms", iCount,
               3, "typex", (char *) &(saPAtom->iTypeIndex), iSize,
               4, "resx", (char *) &(saPAtom->iResidueIndex), iSize,
               5, "flags", (char *) &(saPAtom->fFlags), iSize,
               6, "seq", (char *) &(saPAtom->iSequence), iSize,
               7, "elmnt", (char *) &(saPAtom->iElement), iSize,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               8, "chg", (char *) &(saPAtom->dCharge), iSize,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               1, "name", (char *) saPAtom->sName, iSize,
               2, "type", (char *) saPAtom->sType, iSize,
               0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);

    DBPutTable(db, "atomspertinfo", iCount,
               3, "ptypex", (char *) &(saPAtom->iPertTypeIndex), iSize,
               4, "pelmnt", (char *) &(saPAtom->iPertElement), iSize,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               5, "pchg", (char *) &(saPAtom->dPertCharge), iSize,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               1, "pname", (char *) saPAtom->sPertName, iSize,
               2, "ptype", (char *) saPAtom->sPertType, iSize,
               0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);

    /* Get the atom positions */

    DBPutTable(db, "positions", iCount,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               1, "x", (char *) &dVX(&(saPAtom->vPos)), sizeof(SAVEATOMt),
               2, "y", (char *) &dVY(&(saPAtom->vPos)), sizeof(SAVEATOMt),
               3, "z", (char *) &dVZ(&(saPAtom->vPos)), sizeof(SAVEATOMt),
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);

    /* Get the atom velocities */

    DBPutTable(db, "velocities", iCount,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               1, "x", (char *) &dVX(&(saPAtom->vVelocity)),
               sizeof(SAVEATOMt), 2, "y",
               (char *) &dVY(&(saPAtom->vVelocity)), sizeof(SAVEATOMt), 3,
               "z", (char *) &dVZ(&(saPAtom->vVelocity)),
               sizeof(SAVEATOMt), 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
               NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
               NULL, 0);

    /* Save the bounding box information */

    ToolSanityCheckBox(uUnit);
    sbBox.dBeta = dUnitBeta(uUnit);
    UnitGetBox(uUnit, &(sbBox.dXWidth),
               &(sbBox.dYWidth), &(sbBox.dZWidth));
    if (bUnitUseBox(uUnit))
        sbBox.dUseBox = 1.0;
    else
        sbBox.dUseBox = -1.0;
    DBPutValue(db, "boundbox", ENTRYDOUBLE | ENTRYARRAY,
               sizeof(sbBox) / sizeof(sbBox.dUseBox),        /* # of doubles */
               (GENP) & sbBox, sizeof(sbBox.dUseBox));

    /* Save the cap information */

    if (bUnitUseSolventCap(uUnit))
        scCap.dUseCap = 1.0;
    else
        scCap.dUseCap = -1.0;
    UnitGetSolventCap(uUnit, &(scCap.dX), &(scCap.dY), &(scCap.dZ),
                      &(scCap.dRadius));
    DBPutValue(db, "solventcap", ENTRYDOUBLE | ENTRYARRAY,
               sizeof(scCap) / sizeof(scCap.dUseCap),
               (GENP) & scCap, sizeof(scCap.dUseCap));

    /* Save the connectivity information */

    if ((iCount = iVarArrayElementCount(uUnit->vaConnectivity))) {
        scPConnectivity =
            PVAI(uUnit->vaConnectivity, SAVECONNECTIVITYt, 0);
        iSize = sizeof(SAVECONNECTIVITYt);
        DBPutTable(db, "connectivity", iCount,
                   1, "atom1x", (char *) &(scPConnectivity->iAtom1), iSize,
                   2, "atom2x", (char *) &(scPConnectivity->iAtom2), iSize,
                   3, "flags", (char *) &(scPConnectivity->fFlags), iSize,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save the restraints information */

    if ((iCount = iVarArrayElementCount(uUnit->vaRestraints))) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        iSize = sizeof(SAVERESTRAINTt);
        DBPutTable(db, "restraints", iCount,
                   1, "type", (char *) &(srPRestraint->iType), iSize,
                   2, "flags", (char *) &(srPRestraint->fFlags), iSize,
                   3, "atom1x", (char *) &(srPRestraint->iAtom1), iSize,
                   4, "atom2x", (char *) &(srPRestraint->iAtom2), iSize,
                   5, "atom3x", (char *) &(srPRestraint->iAtom3), iSize,
                   6, "atom4x", (char *) &(srPRestraint->iAtom4), iSize,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   7, "kx", (char *) &(srPRestraint->dKx), iSize,
                   8, "x0", (char *) &(srPRestraint->dX0), iSize,
                   9, "n", (char *) &(srPRestraint->dN), iSize,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save the UNIT connect atoms */

    DBPutValue(db, "connect", ENTRYARRAY | ENTRYINTEGER, 2,
               (GENP) PVAI(uUnit->vaConnect, int, 0), sizeof(int));

    /* Save the array of residues */

    if ((iCount = iVarArrayElementCount(uUnit->vaResidues))) {
        srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        iSize = sizeof(SAVERESIDUEt);
        if (MAXCONNECT != 6)
            DFATAL(("MAXCONNECT has been changed, update UnitSaveTables"));

        DBPutTable(db, "residues", iCount,
                   2, "seq", (char *) &(srPResidue->iSequenceNumber),
                   iSize, 3, "childseq",
                   (char *) &(srPResidue->iNextChildSequence), iSize, 4,
                   "startatomx", (char *) &(srPResidue->iAtomStartIndex),
                   iSize, 6, "imagingx",
                   (char *) &(srPResidue->iImagingAtomIndex), iSize, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 1, "name",
                   (char *) srPResidue->sName, iSize, 5, "restype",
                   (char *) srPResidue->sResidueType, iSize, 0, NULL, NULL,
                   0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);

        if (iVarArrayElementCount(uUnit->vaConnect)) {
            DBPutValue(db, "connect", ENTRYARRAY | ENTRYINTEGER, 2,
                       (GENP) PVAI(uUnit->vaConnect, int, 0), sizeof(int));
        }
        DBPutValue(db, "residuesPdbSequenceNumber",
                   ENTRYARRAY | ENTRYINTEGER, iCount,
                   (GENP) & srPResidue->iPdbSequenceNumber, iSize);

        DBPutTable(db, "residueconnect", iCount,
                   1, "c1x", (char *) &(srPResidue->iaConnectIndex[0]),
                   iSize, 2, "c2x",
                   (char *) &(srPResidue->iaConnectIndex[1]), iSize, 3,
                   "c3x", (char *) &(srPResidue->iaConnectIndex[2]), iSize,
                   4, "c4x", (char *) &(srPResidue->iaConnectIndex[3]),
                   iSize, 5, "c5x",
                   (char *) &(srPResidue->iaConnectIndex[4]), iSize, 6,
                   "c6x", (char *) &(srPResidue->iaConnectIndex[5]), iSize,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save an array of molecules */

    if ((iCount = iVarArrayElementCount(uUnit->vaMolecules))) {
        smPMolecule = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        iSize = sizeof(SAVEMOLECULEt);
        DBPutTable(db, "molecules", iCount,
                   2, "seqnum", (char *) &(smPMolecule->iSequenceNumber),
                   iSize, 3, "childseq",
                   (char *) &(smPMolecule->iNextChildSequence), iSize, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 1, "name", (char *) smPMolecule->sName,
                   iSize, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save an array which describes the hierarchy */

    if ((iCount = iVarArrayElementCount(uUnit->vaHierarchy))) {
        shPHierarchy = PVAI(uUnit->vaHierarchy, SAVEHIERARCHYt, 0);
        iSize = sizeof(SAVEHIERARCHYt);
        DBPutTable(db, "hierarchy", iCount,
                   2, "abovex", (char *) &(shPHierarchy->iAboveIndex),
                   iSize, 4, "belowx",
                   (char *) &(shPHierarchy->iBelowIndex), iSize, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 1, "abovetype",
                   (char *) shPHierarchy->sAboveType, iSize, 3,
                   "belowtype", (char *) shPHierarchy->sBelowType, iSize,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save the group information */

    if (uUnit->vaGroupNames) {
        DBPutValue(db, "groupNames", ENTRYSTRING | ENTRYARRAY,
                   iVarArrayElementCount(uUnit->vaGroupNames),
                   (GENP) PVAI(uUnit->vaGroupNames, STRING, 0),
                   iVarArrayElementSize(uUnit->vaGroupNames));
        sgPGroupAtom = PVAI(uUnit->vaGroupAtoms, SAVEGROUPSt, 0);
        iSize = iVarArrayElementSize(uUnit->vaGroupAtoms);
        DBPutTable(db, "groupAtoms",
                   iVarArrayElementCount(uUnit->vaGroupAtoms),
                   1, "groupIndex", (char *) &(sgPGroupAtom->iGroupIndex),
                   iSize, 2, "atomIndex",
                   (char *) &(sgPGroupAtom->iIndexAtom), iSize, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }

    /*
     *  DELETED: if ( uUnit->psParameters != NULL )
     *               ParmSetSave( uUnit->psParameters, db );
     *  since saving these confuses the precedence of params
     */

    /* Save the bond information */

    if ((iCount = iVarArrayElementCount(uUnit->vaBonds))) {
        sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        iSize = sizeof(SAVEBONDt);
        DBPutTable(db, "bonds", iCount,
                   1, "atom1x", (char *) &(sbPBond->iAtom1), iSize,
                   2, "atom2x", (char *) &(sbPBond->iAtom2), iSize,
                   3, "parmx", (char *) &(sbPBond->iParmIndex), iSize,
                   4, "pertparmx", (char *) &(sbPBond->iPertParmIndex),
                   iSize, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0);
    }

    /* Save the angles information */

    if ((iCount = iVarArrayElementCount(uUnit->vaAngles))) {

        saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
        iSize = sizeof(SAVEANGLEt);
        DBPutTable(db, "angles", iCount,
                   1, "atom1x", (char *) &(saPAngle->iAtom1), iSize,
                   2, "atom2x", (char *) &(saPAngle->iAtom2), iSize,
                   3, "atom3x", (char *) &(saPAngle->iAtom3), iSize,
                   4, "parmx", (char *) &(saPAngle->iParmIndex), iSize,
                   5, "pertparmx", (char *) &(saPAngle->iPertParmIndex),
                   iSize, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0);
    }


    /* Save the torsions information */

    if ((iCount = iVarArrayElementCount(uUnit->vaTorsions))) {

        stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        iSize = sizeof(SAVETORSIONt);
        DBPutTable(db, "torsions", iCount,
                   1, "atom1x", (char *) &(stPTorsion->iAtom1), iSize,
                   2, "atom2x", (char *) &(stPTorsion->iAtom2), iSize,
                   3, "atom3x", (char *) &(stPTorsion->iAtom3), iSize,
                   4, "atom4x", (char *) &(stPTorsion->iAtom4), iSize,
                   5, "parmx", (char *) &(stPTorsion->iParmIndex), iSize,
                   6, "pertparmx", (char *) &(stPTorsion->iPertParmIndex),
                   iSize, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }

    DBPopPrefix(db);
}





/*
 *      zUnitIOTableAddAtom
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Add an atom to the UNITs tables.
 *        Return FALSE in (*bPFailed) if the atom could not be added,
 *        this can happen when the type is unknown.
 */
static void
zUnitIOTableAddAtom(UNIT uUnit, ATOM aAtom, int i, PARMLIB plParameters,
                    BOOL * bPFailed, BOOL bPert)
{
    SAVEATOMt *saPAtom;
    int iIndex, iElement, iHybridization, iTemp;
    double dMass, dPolar, dDepth, dRStar, dDepth14, dRStar14, dScreenF;
    STRING sType, sDesc, sTemp;
    PARMSET psTemp;

    *bPFailed = FALSE;

    /* Define the ATOM index number in the SAVEATOMt array */


    ContainerSetTempInt(aAtom, i + 1);
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);


    saPAtom->aAtom = aAtom;
    REF(aAtom);

    strcpy(saPAtom->sName, sContainerName(aAtom));
    if (strlen(sAtomPertName(aAtom)) != 0) {
        strcpy(saPAtom->sPertName, sAtomPertName(aAtom));
    } else {
        MESSAGE((" no pert atom name set for %s - using nonpert\n",
                 sContainerFullDescriptor((CONTAINER) aAtom, sTemp)));
        strcpy(saPAtom->sPertName, sContainerName(aAtom));
    }
    strcpy(saPAtom->sType, sAtomType(aAtom));
    if (strlen(sAtomPertType(aAtom)) != 0) {
        strcpy(saPAtom->sPertType, sAtomPertType(aAtom));
    } else {
        MESSAGE((" no pert atom type set for %s - using nonpert\n",
                 sContainerFullDescriptor((CONTAINER) aAtom, sTemp)));
        strcpy(saPAtom->sPertType, sAtomType(aAtom));
    }
    saPAtom->dCharge = dAtomCharge(aAtom);
    saPAtom->dPertCharge = dAtomPertCharge(aAtom);
    saPAtom->iSequence = iContainerSequence(aAtom);
    saPAtom->iElement = iAtomElement(aAtom);
    saPAtom->iPertElement = iAtomPertElement(aAtom);

    /* If the atom is contained by a residue then set */
    /* which residue it is, otherwise set it to zero */

    if (iObjectType(cContainerWithin(aAtom)) == RESIDUEid)
        saPAtom->iResidueIndex =
            iContainerTempInt(cContainerWithin(aAtom));
    else
        saPAtom->iResidueIndex = 0;
    saPAtom->vPos = vAtomPosition(aAtom);
    saPAtom->vVelocity = vAtomVelocity(aAtom);
    saPAtom->fFlags = fAtomFlags(aAtom);
    saPAtom->iTypeIndex = 0;
    saPAtom->iPertTypeIndex = 0;

    /* If parameters should be generated then lookup the atom */
    /* in the PARMLIB, and add it to the UNITs PARMSET */
    /* then point the atoms TypeIndex and PertTypeIndex's to */
    /* the entry in the PARMSET */

    if (strlen(sAtomType(aAtom)) == 0) {
        /*
         *  may not matter, since UnitCheck catches
         *      this for saving parm, & not a problem 
         *      for saveoff
         */
        /* try to mess things up.. bPFailed seems to be ignored */
        iIndex = -9999;
    } else if (plParameters != NULL) {
        iIndex = iParmSetFindAtom(uUnit->psParameters, sAtomType(aAtom));
        if (iIndex == PARM_NOT_FOUND) {

            PARMLIB_LOOP(plParameters, psTemp,
                         (iTemp = iParmSetFindAtom(psTemp,
                                                   sAtomType(aAtom))));
            if (iTemp != PARM_NOT_FOUND) {
                ParmSetAtom(psTemp, iTemp, sType,
                            &dMass, &dPolar, &dDepth, &dRStar,
			    &dDepth14, &dRStar14, &dScreenF, &iElement,
			    &iHybridization, sDesc);
                iIndex =
                    iParmSetAddAtom(uUnit->psParameters, sAtomType(aAtom),
                                    dMass, dPolar, dDepth, dRStar,
                                    dDepth14, dRStar14, dScreenF, iElement,
                                    iHybridization, sDesc);
            } else {
                iIndex = 0;
                VP0(("For atom: %s Could not find vdW (or other) parameters for type: %s\n",
                     sContainerFullDescriptor((CONTAINER) aAtom, sTemp),
                     sAtomType(aAtom)));
                *bPFailed = TRUE;
            }
        } else {
            ParmSetAtom(uUnit->psParameters, iIndex, sType,
                        &dMass, &dPolar, &dDepth, &dRStar, &dDepth14,
                        &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc);
        }
        saPAtom->iTypeIndex = iIndex + 1;
        if(iElement != NOELEMENT)
            saPAtom->iElement = iElement;

        if (bPert && bAtomFlagsSet(aAtom, ATOMPERTURB)) {
            iIndex = iParmSetFindAtom(uUnit->psParameters,
                                      saPAtom->sPertType);
            if (iIndex == PARM_NOT_FOUND) {
                PARMLIB_LOOP(plParameters, psTemp,
                             (iTemp = iParmSetFindAtom(psTemp,
                                                       saPAtom->
                                                       sPertType)));
                if (iTemp != PARM_NOT_FOUND) {
                    ParmSetAtom(psTemp, iTemp, sType,
                                &dMass, &dPolar, &dDepth, &dRStar,
                                &dDepth14, &dRStar14, &dScreenF, &iElement,
                                &iHybridization, sDesc);
                    iIndex =
                        iParmSetAddAtom(uUnit->psParameters,
                                        saPAtom->sPertType, dMass, dPolar,
                                        dDepth, dRStar, dDepth14, dRStar14,
					dScreenF,
                                        iElement, iHybridization, sDesc);
                } else {
                    iIndex = 0;
                    VP0(("For atom: %s type %s\n\t%s %s\n",
                         sContainerFullDescriptor((CONTAINER) aAtom,
                                                  sTemp), sType,
                         "Could not find perturbed type: ",
                         saPAtom->sPertType));
                    *bPFailed = TRUE;
                }
            } else {
                ParmSetAtom(uUnit->psParameters, iIndex, sType,
                            &dMass, &dPolar, &dDepth, &dRStar, &dDepth14,
                            &dRStar14, &dScreenF, &iElement, &iHybridization,
			    sDesc);

            }
        }
        if (bPert)
            saPAtom->iPertTypeIndex = iIndex + 1;
    }

}


/*
 *  avl tree related declarations - the 1-4 ones are here
 *        for sharing between routines, the impropers to keep
 *        the 1-4's company.
 */
static IX_REC e14, *ePImp = NULL;
static IX_DESC scr14_index, improper_index;
static int *Pint1, *Pint2;


/*
 *      zUnitIOSetCalc14Flags
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      This routine sets the flags that determine whether
 *      or not 1-4 interactions are calculated for the two
 *      atoms across the torsion.  It does this
 *      by checking the bond tables, and angle tables to see
 *      if the two outer atoms in the torsion are represented
 *      as having a bond or bond angle between them, if this is
 *      TRUE then there is no 1-4 interaction.  Also check
 *      to make sure that no 1-4 interaction has already been
 *      set to be calculated.
 *
 *        The *bPCalc14, *bPCalcPert14 flags are used to determine
 *        if the Calc14 flag has already been set for this torsion.
 *        For each torsion, the first time this routine is called,
 *        both must be set to TRUE.  The first term that has an 
 *        interaction will have it's Calc14 flag set.
 *
 */
static void
zUnitIOSetCalc14Flags(SAVETORSIONt * stPTorsion, BOOL * bPCalc14,
                      BOOL * bPCalcPert14)
{
    BOOL bCheck;

    /*
     *  order the 1st, last pointers by atom #
     */
    if (stPTorsion->iAtom1 < stPTorsion->iAtom4) {
        *Pint1 = stPTorsion->iAtom1;
        *Pint2 = stPTorsion->iAtom4;
    } else {
        *Pint1 = stPTorsion->iAtom4;
        *Pint2 = stPTorsion->iAtom1;
    }

    stPTorsion->bCalc14 = FALSE;
    stPTorsion->bPertCalc14 = FALSE;

    /* Check if we have to check torsions to see if the 1-4 has been set */

    bCheck = FALSE;
    if (stPTorsion->iParmIndex != 0 && *bPCalc14)
        bCheck = TRUE;
    if (stPTorsion->iPertParmIndex != 0 && *bPCalcPert14)
        bCheck = TRUE;

    if (!bCheck)
        return;

    switch (find_key(&e14, &scr14_index)) {
    case IX_OK:
        if (e14.recptr != NULL) {
            /*
             *  bond or angle
             */
            return;
        }
        /*
         *  duplicates a previous torsion's 1-4
         */
        if (stPTorsion->iParmIndex != 0)
            *bPCalc14 = FALSE;

        if (stPTorsion->iPertParmIndex != 0)
            *bPCalcPert14 = FALSE;
        break;
    case IX_FAIL:

        /*
         *  a new torsion - add to index
         */
        e14.recptr = (void *) NULL;
        if (!add_key(&e14, &scr14_index))
            DFATAL(("1-4 torsions: cannot add key %d %d\n",
                    *Pint1, *Pint2));
        break;
    default:
        DFATAL(("unexpected index error\n"));
    }
    if (stPTorsion->iParmIndex != 0) {
        stPTorsion->bCalc14 = *bPCalc14;
        *bPCalc14 = FALSE;
    }
    if (stPTorsion->iPertParmIndex != 0) {
        stPTorsion->bPertCalc14 = *bPCalcPert14;
        *bPCalcPert14 = FALSE;
    }

}




/*
 *        zUnitIndexBondParameters
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *        For all of the bonds, search for their parameters
 *        within the UNITs PARMSET.  If they are
 *        found then set the index to the entry, otherwise
 *        add the bond parameter to the PARMSET and set the index.
 *        If either of the ATOMs in the bond are to be perturbed then
 *        do the same with the perturbation parameter.
 *
 *        Return TRUE if there was a problem generating parameters.
 */
BOOL
zbUnitIOIndexBondParameters(PARMLIB plLib, UNIT uUnit, BOOL bPert)
{
    int iCount, iIndex, iTemp;
    LOOP lTemp;
    SAVEBONDt *sbPBond;
    ATOM aAtom1, aAtom2;
    BOOL bFailedGeneratingParameters;
    double dKb, dR0;
    STRING sAtom1, sAtom2, sDesc;
    PARMSET psTemp;
#ifdef  DEBUG
    STRING sTemp1, sTemp2;
#endif


    bFailedGeneratingParameters = FALSE;

    if (uUnit->vaBonds != NULL) {
        VP0(("Rebuilding bond parameters.\n"));
        VarArrayDestroy(&(uUnit->vaBonds));
    } else
        VP0(("Building bond parameters.\n"));

    uUnit->vaBonds = vaVarArrayCreate(sizeof(SAVEBONDt));
    iCount = 0;
    lTemp = lLoop((OBJEKT) uUnit, BONDS);
    while (oNext(&lTemp) != NULL)
        iCount++;
    VarArraySetSize((uUnit->vaBonds), iCount);
    if (iCount) {
        lTemp = lLoop((OBJEKT) uUnit, BONDS);
        sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        if (sbPBond == NULL)
            DFATAL((" ?? null\n"));
        for (; oNext(&lTemp) != NULL; sbPBond++) {
            LoopGetBond(&lTemp, &aAtom1, &aAtom2);
            sbPBond->iAtom1 = iContainerTempInt(aAtom1);
            sbPBond->iAtom2 = iContainerTempInt(aAtom2);
            sbPBond->iParmIndex = 0;
            sbPBond->iPertParmIndex = 0;
            sbPBond->fFlags = 0;
            strcpy(sAtom1, sAtomType(aAtom1));
            strcpy(sAtom2, sAtomType(aAtom2));
            iIndex = iParmSetFindBond(uUnit->psParameters, sAtom1, sAtom2);
            if (iIndex == PARM_NOT_FOUND) {
                PARMLIB_LOOP(plLib, psTemp,
                             (iTemp = iParmSetFindBond(psTemp,
                                                       sAtom1, sAtom2)));
                if (iTemp != PARM_NOT_FOUND) {
                    ParmSetBond(psTemp, iTemp, sAtom1, sAtom2, &dKb, &dR0,
                                sDesc);
                    iIndex =
                        iParmSetAddBond(uUnit->psParameters, sAtom1,
                                        sAtom2, dKb, dR0, sDesc);
                } else {
                    bFailedGeneratingParameters = TRUE;
                    iIndex = 0;
                    VP0(("Could not find bond parameter for: %s - %s\n",
                         sAtom1, sAtom2));
                }
            }
            sbPBond->iParmIndex = iIndex + 1;

            if (bPert &&
                (bAtomFlagsSet(aAtom1, ATOMPERTURB) ||
                 bAtomFlagsSet(aAtom2, ATOMPERTURB))) {

                /* Note that the bond is perturbed and whether or */
                /* not it is on the boundary between perturbed and */
                /* non-perturbed */
                MESSAGE(("Pert interaction between: %s-%s\n",
                         sContainerFullDescriptor((CONTAINER) aAtom1,
                                                  sTemp1),
                         sContainerFullDescriptor((CONTAINER) aAtom2,
                                                  sTemp2)));
                MESSAGE(
                        ("Indexes = %d-%d\n", sbPBond->iAtom1,
                         sbPBond->iAtom2));
                sbPBond->fFlags |= PERTURBED;
                if (!(bAtomFlagsSet(aAtom1, ATOMPERTURB) &&
                      bAtomFlagsSet(aAtom2, ATOMPERTURB))) {
                    sbPBond->fFlags |= BOUNDARY;
                    MESSAGE(("-Boundary\n"));
                }

                if (bAtomFlagsSet(aAtom1, ATOMPERTURB)) {
                    if (strlen(sAtomPertType(aAtom1)) != 0)
                        strcpy(sAtom1, sAtomPertType(aAtom1));
                    else
                        strcpy(sAtom1, sAtomType(aAtom1));
                } else
                    strcpy(sAtom1, sAtomType(aAtom1));
                if (bAtomFlagsSet(aAtom2, ATOMPERTURB)) {
                    if (strlen(sAtomPertType(aAtom2)) != 0)
                        strcpy(sAtom2, sAtomPertType(aAtom2));
                    else
                        strcpy(sAtom2, sAtomType(aAtom2));
                } else
                    strcpy(sAtom2, sAtomType(aAtom2));
                iIndex = iParmSetFindBond(uUnit->psParameters,
                                          sAtom1, sAtom2);
                if (iIndex == PARM_NOT_FOUND) {
                    PARMLIB_LOOP(plLib, psTemp,
                                 (iTemp = iParmSetFindBond(psTemp,
                                                           sAtom1,
                                                           sAtom2)));
                    if (iTemp != PARM_NOT_FOUND) {
                        ParmSetBond(psTemp, iTemp, sAtom1, sAtom2,
                                    &dKb, &dR0, sDesc);
                        iIndex = iParmSetAddBond(uUnit->psParameters,
                                                 sAtom1, sAtom2, dKb, dR0,
                                                 sDesc);
                    } else {
                        bFailedGeneratingParameters = TRUE;
                        iIndex = 0;
                        VP0(("No bond parameter for: %s - %s\n",
                             sAtom1, sAtom2));
                    }
                }
                sbPBond->iPertParmIndex = iIndex + 1;
            }
        }
    }

    return (bFailedGeneratingParameters);
}


/*
 *        zUnitIndexAngleParameters
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *        For all of the angles, search for their parameters
 *        within the UNITs PARMSET.  If they are
 *        found then set the index to the entry, otherwise
 *        add the angle parameter to the PARMSET and set the index.
 *        If any of the ATOMs in the angle  are to be perturbed then
 *        do the same with the perturbation parameter.
 *
 *        Return TRUE if there was a problem generating parameters.
 */
static BOOL
zbUnitIOIndexAngleParameters(PARMLIB plLib, UNIT uUnit, BOOL bPert)
{
    LOOP lTemp;
    SAVEANGLEt saAngle;
    ATOM aAtom1, aAtom2, aAtom3;
    int iIndex, iTemp;
    BOOL bFailedGeneratingParameters;
    STRING sAtom1, sAtom2, sAtom3, sDesc;
    PARMSET psTemp;
    double dKt, dT0, dTkub, dRkub;

    bFailedGeneratingParameters = FALSE;

    /* Now generate the ANGLE table */

    if (uUnit->vaAngles != NULL) {
        VP0(("Rebuilding angle parameters.\n"));
        VarArrayDestroy(&(uUnit->vaAngles));
    } else
        VP0(("Building angle parameters.\n"));

    uUnit->vaAngles = vaVarArrayCreate(sizeof(SAVEANGLEt));

    lTemp = lLoop((OBJEKT) uUnit, ANGLES);
    while (oNext(&lTemp) != NULL) {
        LoopGetAngle(&lTemp, &aAtom1, &aAtom2, &aAtom3);

        saAngle.iAtom1 = iContainerTempInt(aAtom1);
        saAngle.iAtom2 = iContainerTempInt(aAtom2);
        saAngle.iAtom3 = iContainerTempInt(aAtom3);
        saAngle.iParmIndex = 0;
        saAngle.iPertParmIndex = 0;
        saAngle.fFlags = 0;

        strcpy(sAtom1, sAtomType(aAtom1));
        strcpy(sAtom2, sAtomType(aAtom2));
        strcpy(sAtom3, sAtomType(aAtom3));

/* TODO:Fix this UNGODLY HACK, the PARMSET type has to be modified to */
/* TODO:Allow the user to specify interactions that should be ignored */
/* TODO:Like HW-HW-OW angles for TIP3 waters */

        if (zbUnitIgnoreAngle(sAtom1, sAtom2, sAtom3))
            goto IGNORE1;

        iIndex = iParmSetFindAngle(uUnit->psParameters,
                                   sAtom1, sAtom2, sAtom3);
        if (iIndex == PARM_NOT_FOUND) {
            PARMLIB_LOOP(plLib, psTemp,
                         (iTemp = iParmSetFindAngle(psTemp, sAtom1,
                                                    sAtom2, sAtom3)));
            if (iTemp != PARM_NOT_FOUND) {
                ParmSetAngle(psTemp, iTemp,
                             sAtom1, sAtom2, sAtom3,
                             &dKt, &dT0, &dTkub, &dRkub, sDesc);
                iIndex = iParmSetAddAngle(uUnit->psParameters,
                                          sAtom1, sAtom2, sAtom3,
                                          dKt, dT0, dTkub, dRkub, sDesc);
            } else {
                bFailedGeneratingParameters = TRUE;
                iIndex = 0;
                VP0(("Could not find angle parameter: %s - %s - %s\n",
                     sAtom1, sAtom2, sAtom3));
            }
        }
        saAngle.iParmIndex = iIndex + 1;

      IGNORE1:
        ;

        if (bPert &&
            (bAtomFlagsSet(aAtom1, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom2, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom3, ATOMPERTURB))) {

            /* Note that the angle is perturbed and whether or */
            /* not it is on the boundary between perturbed and */
            /* non-perturbed */
            saAngle.fFlags |= PERTURBED;
            if (!(bAtomFlagsSet(aAtom1, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom2, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom3, ATOMPERTURB)))
                saAngle.fFlags |= BOUNDARY;
            if (bAtomFlagsSet(aAtom1, ATOMPERTURB)) {
                if (strlen(sAtomPertType(aAtom1)) != 0)
                    strcpy(sAtom1, sAtomPertType(aAtom1));
                else
                    strcpy(sAtom1, sAtomType(aAtom1));
            } else
                strcpy(sAtom1, sAtomType(aAtom1));
            if (bAtomFlagsSet(aAtom2, ATOMPERTURB)) {
                if (strlen(sAtomPertType(aAtom2)) != 0)
                    strcpy(sAtom2, sAtomPertType(aAtom2));
                else
                    strcpy(sAtom2, sAtomType(aAtom2));
            } else
                strcpy(sAtom2, sAtomType(aAtom2));
            if (bAtomFlagsSet(aAtom3, ATOMPERTURB)) {
                if (strlen(sAtomPertType(aAtom3)) != 0)
                    strcpy(sAtom3, sAtomPertType(aAtom3));
                else
                    strcpy(sAtom2, sAtomType(aAtom2));
            } else
                strcpy(sAtom3, sAtomType(aAtom3));

/* TODO:Fix this UNGODLY HACK, the PARMSET type has to be modified to */
/* TODO:Allow the user to specify interactions that should be ignored */
/* TODO:Like HW-HW-OW angles for TIP3 waters */

            if (zbUnitIgnoreAngle(sAtom1, sAtom2, sAtom3))
                goto IGNORE2;

            iIndex = iParmSetFindAngle(uUnit->psParameters,
                                       sAtom1, sAtom2, sAtom3);
            if (iIndex == PARM_NOT_FOUND) {
                PARMLIB_LOOP(plLib, psTemp,
                             (iTemp = iParmSetFindAngle(psTemp, sAtom1,
                                                        sAtom2, sAtom3)));
                if (iTemp != PARM_NOT_FOUND) {
                    ParmSetAngle(psTemp, iTemp,
                                 sAtom1, sAtom2, sAtom3,
                                 &dKt, &dT0, &dTkub, &dRkub, sDesc);
                    iIndex = iParmSetAddAngle(uUnit->psParameters,
                                              sAtom1, sAtom2, sAtom3,
                                              dKt, dT0, dTkub, dRkub,
                                              sDesc);
                } else {
                    bFailedGeneratingParameters = TRUE;
                    iIndex = 0;
                    VP0(("Can't find angle parameter: %s - %s - %s\n",
                         sAtom1, sAtom2, sAtom3));
                }
            }
            saAngle.iPertParmIndex = iIndex + 1;
        }
      IGNORE2:

        if (saAngle.iPertParmIndex == 0 && saAngle.iParmIndex == 0)
            continue;

        /* Only add the angle interactions if there is a normal interaction */
        /* or a perturbed one */

        VarArrayAdd(uUnit->vaAngles, (GENP) & saAngle);
    }
    return (bFailedGeneratingParameters);
}


/*
 *  BoilTorsions() - reduce list of torsion params to
 *        unique numerical ones, updating param pointers
 *        in the topological list.
 *
 *        Bill Ross, May 1996
 */
static void
BoilTorsions(VARARRAY * vaPParms, int iParmOffset,
             VARARRAY vaTorsions, int iTorsionOffset)
{
    VARARRAY vaB;
    TORSIONPARMt *tpA, *tpB, tC;
    int iIndex, iA, iB, iParmCount, iTorsionCount;

    strcpy(tC.sDesc, "reduced params");
    strcpy(tC.sType1, "__");
    strcpy(tC.sType2, "__");
    strcpy(tC.sType3, "__");
    strcpy(tC.sType4, "__");
    strcpy(tC.sOrder, "___");

    iParmCount = iVarArrayElementCount(*vaPParms);
    iTorsionCount = iVarArrayElementCount(vaTorsions);

    iIndex = iParmOffset;
    vaB = vaVarArrayCreate(sizeof(TORSIONPARMt));
    tpA = PVAI(*vaPParms, TORSIONPARMt, 0);
    for (iA = 0; iA < iParmCount; iA++, tpA++) {
        int i, iOldIndex;
        SAVETORSIONt *stP;

        if (!strcmp(tpA->sType1, "__"))
            continue;
        /*
         *  torsion hasn't been marked as 'superfluous'
         *      so add to new array
         */
        tC.dKp = tpA->dKp;
        tC.iN = tpA->iN;
        tC.dP0 = tpA->dP0;
	tC.dScEE = tpA->dScEE;
	tC.dScNB = tpA->dScNB;
        VarArrayAdd(vaB, (GENP) & tC);
        iIndex++;
        iOldIndex = iParmOffset + iA + 1;

        /*
         *  update any affected torsions
         */
        if (iIndex != iOldIndex) {
            stP = PVAI(vaTorsions, SAVETORSIONt, iTorsionOffset);
            for (i = iTorsionOffset; i < iTorsionCount; i++, stP++) {
                if (stP->iParmIndex == iOldIndex)
                    stP->iParmIndex = iIndex;
                if (stP->iPertParmIndex == iOldIndex)
                    stP->iPertParmIndex = iIndex;
            }
        }

        /*
         *  mark any subsequent duplicates 'superfluous'
         *      and update indexes into the array
         */
        for (tpB = tpA + 1, iB = iA + 1; iB < iParmCount; iB++, tpB++) {
            if (!strcmp(tpB->sType1, "__"))
                continue;
            if (tpB->iN != tpA->iN)
                continue;
            if (tpB->dKp != tpA->dKp)
                continue;
            if (tpB->dP0 != tpA->dP0)
                continue;
	    if (tpB->dScEE != tpA->dScEE)
		continue;
	    if (tpB->dScNB != tpA->dScNB)
		continue;

            /*
             *  B is a duplicate of A
             */
            strcpy(tpB->sType1, "__");
            iOldIndex = iParmOffset + iB + 1;
            stP = PVAI(vaTorsions, SAVETORSIONt, iTorsionOffset);
            for (i = iTorsionOffset; i < iTorsionCount; i++, stP++) {
                if (stP->iParmIndex == iOldIndex)
                    stP->iParmIndex = iIndex;
                if (stP->iPertParmIndex == iOldIndex)
                    stP->iPertParmIndex = iIndex;
            }

        }
    }

    /*
     *  throw away the old parms array & put the new one in place
     */
    VarArrayDestroy(vaPParms);
    *vaPParms = vaB;
}




/*
 *        zbUnitIOIndexTorsionParameters
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *        For all of the angles, search for their parameters
 *        within the UNITs PARMSET.  If they are
 *        found then set the index to the entry, otherwise
 *        add the angle parameter to the PARMSET and set the index.
 *        If any of the ATOMs in the angle  are to be perturbed then
 *        do the same with the perturbation parameter.
 *
 *        (bProper) is TRUE if generating parameters for PROPER torsions,
 *        otherwise IMPROPER torsions.
 *
 *        Return TRUE if there was a problem generating parameters.
 *
 *      Omitting the cacheing function of the UNIT's PARMSET, here
 *      is the intended algorithm:
 *
 *        foreach dihedral in molecule
 *          foreach parmset, working back from the most recently loaded
 *            foreach dihedral in parmset
 *              if an exact match is found
 *                throw away any previous wild card matches
 *                save term (and continue collecting terms in this parmset)
 *              endif
 *              if a wild card match is found and nothing else found yet
 *                  save it
 *              endif
 *            end/parmset_dihedrals
 *            *if an exact match found*
 *              don't search any more parmsets
 *            endif
 *          end/parmsets
 *          ..plug params in for molecule
 *        end/molecule_dihedrals
 */
extern int itest;
static BOOL
zbUnitIOIndexTorsionParameters(PARMLIB plLib, UNIT uUnit,
                               BOOL bProper, BOOL bPert )
{
    LOOP lTemp;
    SAVETORSIONt stTorsion;
    ATOM aAtom1, aAtom2, aAtom3, aAtom4;
    STRING sAtom1, sAtom2, sAtom3, sAtom4;
    STRING sPert1, sPert2, sPert3, sPert4;
    STRING sOrigAtom1, sOrigAtom2, sOrigAtom3, sOrigAtom4;
    STRING sOrigPert1, sOrigPert2, sOrigPert3, sOrigPert4;
    TORSION tTorsion, tPertTorsion;
    BOOL bPerturbTorsion;
    PARMSET psTemp;
    int iTerm, iPertTerm;
    BOOL bDone, bUse, bUsePert, bCopy, bCopyPert, bEnd, bPertEnd;
    int iN, iPertIndex, iPertN, iLastN, iLastPertN;
    double dKp, dP0, dPertKp, dPertP0;
    double dScEE, dScNB, dPScEE, dPScNB;
    BOOL bCalc14, bCalcPert14;
#ifdef  DEBUG2
    STRING s1, s2, s3, s4;
    int iTParm, iTmp;
    double dTK, dTP;
    STRING sT1, sT2, sT3, sT4, sTemp;
#endif
        STRING sDesc;
    BOOL bFailedGeneratingParameters;
    int iTN, iTNPert = 0;
    int iaIndexes[4];
    char *cPaTypes[4];
    int iImproper = 0, iParmOffset = 0, iTorsionOffset = 0;
    int iCount = 0, i, iIndex;

#define                MAX_N                9999

    bFailedGeneratingParameters = FALSE;

    VP0(("Building %s torsion parameters.\n",
         (bProper ? "proper" : "improper")));
    if (!bProper) {
        iParmOffset = iParmSetProperCount(uUnit->psParameters);
        iTorsionOffset = iVarArrayElementCount(uUnit->vaTorsions);
    }
    /* 
     *  NOTE: In order for the 1-4 interactions to be calculated 
     *  properly, add constraint bonds and angles AFTER 
     *  the proper torsions are added.  This allows the 1-4 
     *  interaction checker to use the bond lists and angle 
     *  lists to check connectivity
     */
    if (bProper) {
        int iMax;

        /*
         *  set up 1-4 checking stuff by initializing index with bond
         *  and angle atom pairs
         *
         *  create index, non-duplicate keys, key length == 2 ints
         */
        create_index(&scr14_index, 0, 2 * sizeof(int));

        /*
         *  set up convenience pointers for stuffing key
         *      with 2 ints
         */
        Pint1 = (int *) &e14.key;
        Pint2 = Pint1 + 1;

        /* 
         *  Using the index pointer as a flag; 1==bond,angle, 0==torsion.
         *      The bonds should all be unique pairs. An angle could
         *      duplicate a bond if there is a 'triangle', and could
         *      duplicate an angle if there is a 'square'.
         */
        e14.recptr = (RECPOS) 1;

        /* plop in bonded pairs if any */
        if ((iMax = iVarArrayElementCount(uUnit->vaBonds))) {
            SAVEBONDt *sbPBondT = PVAI(uUnit->vaBonds, SAVEBONDt, 0);

            for (i = 0; i < iMax; i++, sbPBondT++) {
                if (sbPBondT->iAtom1 < sbPBondT->iAtom2) {
                    *Pint1 = sbPBondT->iAtom1;
                    *Pint2 = sbPBondT->iAtom2;
                } else {
                    *Pint1 = sbPBondT->iAtom2;
                    *Pint2 = sbPBondT->iAtom1;
                }
                if (add_key(&e14, &scr14_index) != IX_OK)
                    DFATAL(
                           ("1-4: cannot add bond %d %d\n%s", *Pint1, *Pint2,
                            "This may be caused by duplicate bond "
                            "specifications;\n"
                            "for example, explicit bond commands in addition "
                            "to PDB conect records.\n"));
            }
        }
        /* plop in angle pairs */
        if ((iMax = iVarArrayElementCount(uUnit->vaAngles))) {
            SAVEANGLEt *saPAngleT = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
            for (i = 0; i < iMax; i++, saPAngleT++) {
                if (saPAngleT->iAtom1 < saPAngleT->iAtom3) {
                    *Pint1 = saPAngleT->iAtom1;
                    *Pint2 = saPAngleT->iAtom3;
                } else {
                    *Pint1 = saPAngleT->iAtom3;
                    *Pint2 = saPAngleT->iAtom1;
                }
                if (add_key(&e14, &scr14_index) != IX_OK)
                    VP0(("1-4: angle %d %d %s %s\n",
                         *Pint1, *Pint2,
                         "duplicates bond ('triangular' bond)",
                         "or angle ('square' bond)\n"));
            }
        }
    } else {
        /*
         *  Impropers - use an index to track instantiations of
         *      the pure types templated in the ff to actual
         *      residue & atom names. For this, we use string-based
         *      indexing to count duplicate instantiations, as opposed 
         *      to the fixed-length (integer) records used for
         *      1-4 tracking in the proper torsions.
         *      
         *  create index, 'count'-style duplicate keys, key length == string
         *  create record w/ key portion longer than the default
         */
        create_index(&improper_index, 2, 0);
        MALLOC(ePImp, IX_REC *, sizeof(IX_REC) + 80);
        /*
         *  set recptr to avoid mem check complaint 
         *      (since this index doesn't use its 
         *      recptr)        
         */
        ePImp->recptr = NULL;
    }

    lTemp = lLoop((OBJEKT) uUnit, (bProper ? PROPERS : IMPROPERS));

    while (oNext(&lTemp) != NULL) {
        if (bProper) {
            LoopGetTorsion(&lTemp, &aAtom1, &aAtom2, &aAtom3, &aAtom4);
        } else {
            LoopGetImproper(&lTemp, &aAtom1, &aAtom2, &aAtom3, &aAtom4);
        }

/*
        fprintf( stderr, "%s torsion %s %s %s %s\n", 
                        (bProper ? "Proper" : "Improper"),
                        sAtomName(aAtom1), sAtomName(aAtom2),
                        sAtomName(aAtom3), sAtomName(aAtom4) );
*/

        MESSAGE(("%s torsion %s %s %s %s\n",
                 (bProper ? "Proper" : "Improper"),
                 sAtomName(aAtom1), sAtomName(aAtom2),
                 sAtomName(aAtom3), sAtomName(aAtom4)));

        stTorsion.bProper = bProper;
        stTorsion.bCalc14 = FALSE;
        stTorsion.bPertCalc14 = FALSE;


        /* iContainerTempInt contains atom indices into SAVEATOMt arrays */

        stTorsion.iAtom1 = iContainerTempInt(aAtom1);
        stTorsion.iAtom2 = iContainerTempInt(aAtom2);
        stTorsion.iAtom3 = iContainerTempInt(aAtom3);
        stTorsion.iAtom4 = iContainerTempInt(aAtom4);
        stTorsion.fFlags = 0;

        /* Define the names of the unperturbed ATOMs */

        strcpy(sAtom1, sAtomType(aAtom1));
        strcpy(sAtom2, sAtomType(aAtom2));
        strcpy(sAtom3, sAtomType(aAtom3));
        strcpy(sAtom4, sAtomType(aAtom4));
        strcpy(sOrigAtom1, sAtom1);
        strcpy(sOrigAtom2, sAtom2);
        strcpy(sOrigAtom3, sAtom3);
        strcpy(sOrigAtom4, sAtom4);

        /* skip if this torsion refers to an extra point:  */
        if( GDefaults.iDeleteExtraPointAngles ){
            if( strcmp( sAtom1, "EP" ) == 0 || strcmp( sAtom4, "EP" ) == 0 )
            continue;
        }

        /* Check if the torsion is to be perturbed, if it */
        /* is then set flags saying so, and create a TORSION for */
        /* the perturbation */
        bPerturbTorsion = FALSE;
        if (bPert &&
            (bAtomFlagsSet(aAtom1, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom2, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom3, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom4, ATOMPERTURB))) {
            bPerturbTorsion = TRUE;


            /* Note that the torsion is perturbed and whether or */
            /* not it is on the boundary between perturbed and */
            /* non-perturbed */

            stTorsion.fFlags |= PERTURBED;
            if (!(bAtomFlagsSet(aAtom1, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom2, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom3, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom4, ATOMPERTURB))) {
                stTorsion.fFlags |= BOUNDARY;
                MESSAGE(("Boundary torsion: %s-%s-%s-%s\n",
                         sAtom1, sAtom2, sAtom3, sAtom4));
            }

            /* Define the names of the perturbed atoms */
            if (bAtomFlagsSet(aAtom1, ATOMPERTURB) &&
                strlen(sAtomPertType(aAtom1)) != 0)
                strcpy(sPert1, sAtomPertType(aAtom1));
            else
                strcpy(sPert1, sAtom1);
            if (bAtomFlagsSet(aAtom2, ATOMPERTURB) &&
                strlen(sAtomPertType(aAtom2)) != 0)
                strcpy(sPert2, sAtomPertType(aAtom2));
            else
                strcpy(sPert2, sAtom2);
            if (bAtomFlagsSet(aAtom3, ATOMPERTURB) &&
                strlen(sAtomPertType(aAtom3)) != 0)
                strcpy(sPert3, sAtomPertType(aAtom3));
            else
                strcpy(sPert3, sAtom3);
            if (bAtomFlagsSet(aAtom4, ATOMPERTURB) &&
                strlen(sAtomPertType(aAtom4)) != 0)
                strcpy(sPert4, sAtomPertType(aAtom4));
            else
                strcpy(sPert4, sAtom4);
            strcpy(sOrigPert1, sPert1);
            strcpy(sOrigPert2, sPert2);
            strcpy(sOrigPert3, sPert3);
            strcpy(sOrigPert4, sPert4);
        }

        /* First search through the UNITs PARMSET for the */
        /* torsion parameters that we need. */

        tTorsion = tParmSetTORSIONCreate();

        if (bProper) {
            iTN = iParmSetFindProperTerms(uUnit->psParameters,
                                          tTorsion, TRUE,
                                          sAtom1, sAtom2, sAtom3, sAtom4);
        } else {
            iTN = iParmSetFindImproperTerms(uUnit->psParameters,
                                            tTorsion, TRUE,
                                            sAtom1, sAtom2,
                                            sAtom3, sAtom4);
        }

        if (iTN != PARM_NOT_FOUND) {
            MESSAGE(("Found existing %s terms in UNIT PARMSET.\n",
                     (bProper ? "PROPER" : "IMPROPER")));
        }
        if (bPerturbTorsion) {
            tPertTorsion = tParmSetTORSIONCreate();
            if (bProper) {
                iTNPert = iParmSetFindProperTerms(uUnit->psParameters,
                                                  tPertTorsion, TRUE,
                                                  sPert1, sPert2,
                                                  sPert3, sPert4);
            } else {
                iTNPert = iParmSetFindImproperTerms(uUnit->psParameters,
                                                    tPertTorsion, TRUE,
                                                    sPert1, sPert2,
                                                    sPert3, sPert4);
            }
            if (iTNPert != PARM_NOT_FOUND) {
                MESSAGE(("Found existing %s terms in UNIT PARMSET.\n",
                         (bProper ? "PROPER" : "IMPROPER")));
            }
        }

        if (iTN != PARM_FOUND_EXACT ||
            (bPerturbTorsion && iTNPert != PARM_FOUND_EXACT)) {

            /* search through the PARMLIBs */
            /* put 1st torsion into tTorsion since
               hopefully we're searching back from most
               recent parm loaded */

            PARMLIB_LOOP_ALL(plLib, psTemp) {
                if (iTN != PARM_FOUND_EXACT) {
                    if (bProper) {
                        iTN = iParmSetFindProperTerms(psTemp,
                                                      tTorsion, FALSE,
                                                      sAtom1, sAtom2,
                                                      sAtom3, sAtom4);
                    } else {
                        iTN = iParmSetFindImproperTerms(psTemp,
                                                        tTorsion, FALSE,
                                                        sAtom1, sAtom2,
                                                        sAtom3, sAtom4);
                    }
                }
                if (bPerturbTorsion && iTNPert != PARM_FOUND_EXACT) {
                    if (bProper) {
                        iTNPert = iParmSetFindProperTerms(psTemp,
                                                          tPertTorsion,
                                                          FALSE, sPert1,
                                                          sPert2, sPert3,
                                                          sPert4);
                    } else {
                        iTNPert = iParmSetFindImproperTerms(psTemp,
                                                            tPertTorsion,
                                                            FALSE, sPert1,
                                                            sPert2, sPert3,
                                                            sPert4);
                    }
                }
                if (iTN == PARM_FOUND_EXACT) {
                    if (bPerturbTorsion) {
                        if (iTNPert == PARM_FOUND_EXACT)
                            break;
                    } else
                        break;
                }
            }
        }
#ifdef        DEBUG2
        MESSAGE(("%s %s-%s-%s-%s found %d terms\n",
                 (bProper ? "PROPER" : "IMPROPER"),
                 sAtom1, sAtom2, sAtom3, sAtom4,
                 iParmSetTORSIONTermCount(tTorsion)));
        for (i = 0; i < iParmSetTORSIONTermCount(tTorsion); i++) {
            ParmSetTORSIONTerm(tTorsion, i,
                               &iTParm,
                               sT1, sT2, sT3, sT4,
                               &iTmp, &dTK, &dTP, sTemp);
            MESSAGE(("Term %3d  %d %s-%s-%s-%s  %d  %lf  %lf\n",
                     i, iTParm, sT1, sT2, sT3, sT4, iTmp, dTK, dTP));
        }
        if (bPerturbTorsion) {
            MESSAGE(("Pert%s %s-%s-%s-%s found %d terms\n",
                     (bProper ? "PROPER" : "IMPROPER"),
                     sPert1, sPert2, sPert3, sPert4,
                     iParmSetTORSIONTermCount(tPertTorsion)));
            for (i = 0; i < iParmSetTORSIONTermCount(tPertTorsion); i++) {
                ParmSetTORSIONTerm(tPertTorsion, i,
                                   &iTParm,
                                   sT1, sT2, sT3, sT4,
                                   &iTmp, &dTK, &dTP, sTemp);
                MESSAGE(("Term %3d  %d %s-%s-%s-%s  %d  %lf  %lf\n",
                         i, iTParm, sT1, sT2, sT3, sT4, iTmp, dTK, dTP));
            }
        }
#endif


        /* Now we have a complete TORSION in tTorsion, and */
        /* if the torsion is being perturbed we have another */
        /* complete TORSION in tPertTorsion */
        /* The (stTorsion) terms must now be created */

        iTerm = 0;
        iPertTerm = 0;
        bDone = FALSE;
        iIndex = PARM_NOT_FOUND;
        iPertIndex = PARM_NOT_FOUND;
        bEnd = FALSE;
        bPertEnd = FALSE;
        iN = MAX_N;
        iPertN = MAX_N;
        bCalc14 = TRUE;
        bCalcPert14 = TRUE;

        /*
         *  get 1st term
         */
        if (iParmSetTORSIONTermCount(tTorsion) != 0) {
            ParmSetTORSIONTerm(tTorsion, iTerm,
                               &iIndex,
                               sAtom1, sAtom2, sAtom3, sAtom4,
                               &iN, &dKp, &dP0, &dScEE, &dScNB,
			       sDesc);
            MESSAGE(("First non-perturbed multiplicity: %d\n", iN));
        } else {
            if (bProper) {
                VP0((" ** No torsion terms for  %-s-%-s-%-s-%-s\n",
                     sOrigAtom1, sOrigAtom2, sOrigAtom3, sOrigAtom4));
                bFailedGeneratingParameters = TRUE;
            } else if ( iAtomHybridization(aAtom3) == 2 ){
                VP1((" ** Warning: No sp2 improper torsion term for  %-s-%-s-%-s-%-s\n",
                     sOrigAtom1, sOrigAtom2, sOrigAtom3, sOrigAtom4));
                        VP1(( "        atoms are: %s %s %s %s\n", 
                        sAtomName(aAtom1), sAtomName(aAtom2),
                        sAtomName(aAtom3), sAtomName(aAtom4) ));
            }
            bEnd = TRUE;
        }
        if (bPerturbTorsion) {
            if (iParmSetTORSIONTermCount(tPertTorsion) != 0) {
                ParmSetTORSIONTerm(tPertTorsion, iPertTerm,
                                   &iPertIndex,
                                   sPert1, sPert2, sPert3, sPert4,
                                   &iPertN, &dPertKp, &dPertP0, &dPScEE,
				   &dPScNB, sDesc);
                MESSAGE(("First perturbed multiplicity: %d\n", iPertN));
            } else
                bPertEnd = TRUE;
        }
        if (!bPerturbTorsion)
            bDone = bEnd;
        else
            bDone = bEnd && bPertEnd;

        /* If the interaction is an improper then reorder */
        /* the atoms to get them in the same order as */
        /* was in the original parameter set */

        if (!bProper && iParmSetTORSIONTermCount(tTorsion) > 0) {
            cPaTypes[0] = sAtomType(aAtom1);
            cPaTypes[1] = sAtomType(aAtom2);
            cPaTypes[2] = sAtomType(aAtom3);
            cPaTypes[3] = sAtomType(aAtom4);
            iaIndexes[0] = stTorsion.iAtom1;
            iaIndexes[1] = stTorsion.iAtom2;
            iaIndexes[2] = stTorsion.iAtom3;
            iaIndexes[3] = stTorsion.iAtom4;
            MESSAGE(("Old order: %d %d %d %d\n",
                     stTorsion.iAtom1,
                     stTorsion.iAtom2,
                     stTorsion.iAtom3, stTorsion.iAtom4));

            ParmSetImproperOrderAtoms(tTorsion, 0, cPaTypes, iaIndexes);

            stTorsion.iAtom1 = iaIndexes[0];
            stTorsion.iAtom2 = iaIndexes[1];
            stTorsion.iAtom3 = iaIndexes[2];
            stTorsion.iAtom4 = iaIndexes[3];
            MESSAGE(("New order: %d %d %d %d\n",
                     stTorsion.iAtom1,
                     stTorsion.iAtom2,
                     stTorsion.iAtom3, stTorsion.iAtom4));
        }

        /* Loop over all of the terms */

        while (!bDone) {
            bUse = FALSE;
            bUsePert = FALSE;
            bCopy = FALSE;
            bCopyPert = FALSE;
            stTorsion.iParmIndex = 0;
            stTorsion.iPertParmIndex = 0;

            /* Get the next term, the one with the */
            /* lowest multiplicity, and advance to the */
            /* next multiplicity */

            if (!bPerturbTorsion) {
                /* Advance to the next term within the TORSION */
                bUse = TRUE;
                if (iIndex == PARM_NOT_FOUND)
                    bCopy = TRUE;
            } else {

                /* If the multiplicity is the same for */
                /* both the nonperturbed and perturbed */
                /* term then advance them both */
                if (iPertN == iN) {
                    bUse = TRUE;
                    bUsePert = TRUE;
                    if (iIndex == PARM_NOT_FOUND)
                        bCopy = TRUE;
                    if (iPertIndex == PARM_NOT_FOUND)
                        bCopyPert = TRUE;
                } else if (iN < iPertN) {
                    bUse = TRUE;
                    if (iIndex == PARM_NOT_FOUND)
                        bCopy = TRUE;
                } else {
                    bUsePert = TRUE;
                    if (iPertIndex == PARM_NOT_FOUND)
                        bCopyPert = TRUE;
                }
            }

            MESSAGE(
                    ("Flags:  bUse:%d bCopy:%d  bUsePert:%d bCopyPert:%d\n",
                     bUse, bCopy, bUsePert, bCopyPert));

            /* Now save the terms into the UNITs PARMSET if */
            /* they are not already there */

            if (bCopy) {
                if (bProper)
                    iIndex = iParmSetAddProperTerm(uUnit->psParameters,
                                                   sAtom1, sAtom2, sAtom3,
                                                   sAtom4, iN, dKp, dP0,
						   dScEE, dScNB,
                                                   sDesc);
/*                else if ( !GDefaults.iCharmm )    ???---should I do this????     */
                else
                    iIndex = iParmSetAddImproperTerm(uUnit->psParameters,
                                                     sAtom1, sAtom2,
                                                     sAtom3, sAtom4, iN,
                                                     dKp, dP0, dScEE, dScNB,sDesc);
            }
            if (bCopyPert) {
                if (bProper) {
                    iPertIndex = iParmSetAddProperTerm(uUnit->psParameters,
                                                       sPert1, sPert2,
                                                       sPert3, sPert4,
                                                       iPertN, dPertKp,
                                                       dPertP0, dScEE,
						       dScNB, sDesc);
                } else {
                    iPertIndex =
                        iParmSetAddImproperTerm(uUnit->psParameters,
                                                sPert1, sPert2, sPert3,
                                                sPert4, iPertN, dPertKp,
                                                dPertP0,dScEE,  dScNB,  sDesc);
                }
                MESSAGE(("iPertIndex = %d\n", iPertIndex));
            }

            /* Save the multiplicity */

            iLastN = iN;
            iLastPertN = iPertN;
            if (bUse) {
                stTorsion.iParmIndex = iParmOffset + iIndex + 1;
                iTerm++;
                if (iTerm >= iParmSetTORSIONTermCount(tTorsion)) {
                    bEnd = TRUE;
                    iN = MAX_N;
                } else {
                    ParmSetTORSIONTerm(tTorsion, iTerm,
                                       &iIndex,
                                       sAtom1, sAtom2, sAtom3, sAtom4,
                                       &iN, &dKp, &dP0, &dScEE, &dScNB,
				       sDesc);
                }
                MESSAGE(
                        ("Advancing non-perturbed multiplicity to %d\n",
                         iN));
            }
            if (bUsePert) {
                stTorsion.iPertParmIndex = iParmOffset + iPertIndex + 1;
                iPertTerm++;
                if (iPertTerm >= iParmSetTORSIONTermCount(tPertTorsion)) {
                    bPertEnd = TRUE;
                    iPertN = MAX_N;
                } else {
                    ParmSetTORSIONTerm(tPertTorsion, iPertTerm,
                                       &iPertIndex,
                                       sPert1, sPert2, sPert3, sPert4,
                                       &iPertN, &dPertKp, &dPertP0, &dScEE,
				       &dScNB, sDesc);
                }
                MESSAGE(
                        ("Advancing perturbed multiplicity to %d\n",
                         iPertN));
            }
            if (bProper)
                zUnitIOSetCalc14Flags(&stTorsion, &bCalc14, &bCalcPert14);
            /* 
             *  If we are writing a perturbation topology file,
             *  for every unperturbed term for a given torsion,
             *  there must be a perturbed term (and vice-versa).  
             *  This is because multiple torsional potentials
             *  may apply to a single torsion, and each is perturbed
             *  individually in gibbs.
             */

            /* At this point iLastN, iLastPertN should be the */
            /* multiplicity of the torsion term */

            if (bPerturbTorsion) {
                if (stTorsion.iParmIndex == 0 ||
                    stTorsion.iPertParmIndex == 0) {
                    bFailedGeneratingParameters = TRUE;
                    VP0(("*** %s torsion parameters missing ***\n",
                         (bProper ? "Proper" : "Improper")));
                    VP0((" atom names: %-s-%-s-%-s-%-s\n",
                         sAtomName(aAtom1), sAtomName(aAtom2),
                         sAtomName(aAtom3), sAtomName(aAtom4)));

                    VP0(
                        (" atom types: %-s-%-s-%-s-%-s  =pert=>  %-s-%-s-%-s-%-s\n",
                         sOrigAtom1, sOrigAtom2, sOrigAtom3, sOrigAtom4,
                         sOrigPert1, sOrigPert2, sOrigPert3, sOrigPert4));
                    /* note: missing multiplicity is for the _other_ state */
                    VP0(
                        ("Please add a dummy parameter of multiplicity %d\n",
                         (stTorsion.iParmIndex !=
                          0 ? iLastN : iLastPertN)));
                    VP0(
                        ("for the %spert types to your parameter set.\n",
                         (stTorsion.iParmIndex == 0 ? "non-" : "")));
                    VP0(
                        (" - e.g. %-s-%-s-%-s-%-s  %s    0.0     0.       %d.\n",
                         (stTorsion.iParmIndex ==
                          0 ? sOrigAtom1 : sOrigPert1),
                         (stTorsion.iParmIndex ==
                          0 ? sOrigAtom2 : sOrigPert2),
                         (stTorsion.iParmIndex ==
                          0 ? sOrigAtom3 : sOrigPert3),
                         (stTorsion.iParmIndex ==
                          0 ? sOrigAtom4 : sOrigPert4),
                         (bProper ? "1" : ""),
                         (stTorsion.iParmIndex !=
                          0 ? iLastN : iLastPertN)));

                    VP0(("%s %s\n%s %s\n", "(This is because multiple",
                         "torsional potentials may apply to a",
                         "single torsion, and each is perturbed",
                         "individually in gibbs.)"));
                }
            }

            /* Add the term to the table */
            VarArrayAdd(uUnit->vaTorsions, (GENP) & stTorsion);
            iCount++;
            if (!bProper) {
                STRING sDesc1, sDesc2, sDesc3, sDesc4;

                AtomDescStr(aAtom1, FALSE, sDesc1);
                AtomDescStr(aAtom2, FALSE, sDesc2);
                AtomDescStr(aAtom3, FALSE, sDesc3);
                AtomDescStr(aAtom4, FALSE, sDesc4);
                sprintf(ePImp->key, "%s - %s - %s - %s",
                        sDesc1, sDesc2, sDesc3, sDesc4);
                if (add_key(ePImp, &improper_index) != IX_OK)
                    DFATAL(("add_key() impropers\n"));

                iImproper++;
            }
#ifdef        DEBUG2
            MESSAGE(("Adding %s : %s - %s - %s - %s\n",
                     (bProper ? "PROPER" : "IMPROPER"),
                     sContainerFullDescriptor((CONTAINER) aAtom1, s1),
                     sContainerFullDescriptor((CONTAINER) aAtom2, s2),
                     sContainerFullDescriptor((CONTAINER) aAtom3, s3),
                     sContainerFullDescriptor((CONTAINER) aAtom4, s4)));
            MESSAGE(("Perturbed: %s\n", sBOOL(bPerturbTorsion)));
#endif
            if (!bPerturbTorsion) {
                bDone = bEnd;
            } else {
                bDone = bEnd && bPertEnd;
            }
        }

        /* Release the memory used by the TORSIONs */

        ParmSetTORSIONDestroy(&tTorsion);

        if (bPerturbTorsion) {
            ParmSetTORSIONDestroy(&tPertTorsion);
        }
    }

    if (bProper) {

        /*
         *  throw away 1-4 scratch index
         */
        destroy_index(&scr14_index);

        /*
         *  'boil down' numerically redundant params
         */
        BoilTorsions(&uUnit->psParameters->vaTorsions, 0,        /* parm table */
                     uUnit->vaTorsions, 0);        /* atom lists */

    } else {
        RESIDUE rRes;
        IMPROPERt *IP;
        int iImproper2 = 0, iPrep = 0;

        BoilTorsions(&uUnit->psParameters->vaImpropers,        /* parm table */
                     iParmOffset, uUnit->vaTorsions,        /* atom lists */
                     iTorsionOffset);

        lTemp = lLoop((OBJEKT) uUnit, RESIDUES);
        while ((rRes = (RESIDUE) oNext(&lTemp))) {
            if (rRes->vaImpropers) {
                if (iPrep++ == 0)
                    VP0(("old PREP-specified impropers:\n"));
            }
            if ((iCount = iVarArrayElementCount(rRes->vaImpropers))) {
                sContainerDescriptor((CONTAINER) rRes, sDesc);
                IP = PVAI(rRes->vaImpropers, IMPROPERt, 0);
                for (i = 0; i < iCount; i++, IP++) {
                    iImproper2++;
                    VP0((" %s:  %.4s %.4s %.4s %.4s\n", sDesc + 1,
                         IP->sName1, IP->sName2, IP->sName3, IP->sName4));
                }
            }
        }
        if (iImproper && GiVerbosityLevel > 0) {

            VP0(("--Impropers:\n"));
            first_key(&improper_index);
            while (next_key(ePImp, &improper_index) == IX_OK)
                VP1(("  %d\t%s\n", ePImp->count, ePImp->key));
        }
        VP0((" total %d improper torsion%s applied\n",
             iImproper, (iImproper != 1 ? "s" : "")));
        if (iPrep)
            VP0((" %d improper torsions in old prep form\n", iImproper2));
        destroy_index(&improper_index);
        FREE(ePImp);
    }


    return (bFailedGeneratingParameters);
}


/*
 *        zUnitIOFindAndCountMolecules
 *
 *        The caller must supply a VARARRAY whose elements are (int)s.
 *
 *        Loop through the residues list and count the number
 *        of molecules.
 *        For each molecule, count the number of ATOMs and place
 *        that count in a VARARRAY.  Return the VARARRAY, the number
 *        of molecules counted, and the index of the first solvent molecule.
 */
static void
zUnitIOFindAndCountMolecules(UNIT uUnit, VARARRAY * vaPMolecules,
                             int *iPFirstSolvent)
{
    SAVERESIDUEt *srPRes;
    int i, iResidues, iAtom, iCount;
    LOOP lSpanning;
    BOOL bSeenFirstSolvent;
    ATOM aAtom;

    bSeenFirstSolvent = FALSE;

    /* Clear the ATOMTOUCHED flag on all the ATOMs */

    ContainerResetAllAtomsFlags((CONTAINER) uUnit, ATOMTOUCHED);

    /* Get the first RESIDUE */

    srPRes = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
    iResidues = iVarArrayElementCount(uUnit->vaResidues);

    /* Loop over all RESIDUES */

    for (i = 0; i < iResidues; i++, srPRes++) {

        /* Search for the next RESIDUE whose first ATOM has not */
        /* been touched */

        iAtom = srPRes->iAtomStartIndex - 1;
        if (iAtom < 0) {
            /* skip empty residue */
            continue;
        }
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iAtom)->aAtom;
        if (bAtomFlagsSet(aAtom, ATOMTOUCHED)) {
            continue;
        }

        /* Touch all of the ATOMs within the molecule that */
        /* contains the current RESIDUE */

        iCount = 0;
        lSpanning = lLoop((OBJEKT) aAtom, SPANNINGTREE);
        FOREACH(aAtom, ATOM, lSpanning) {
            AtomSetFlags(aAtom, ATOMTOUCHED);
            iCount++;
        }

        /* Add the molecule to the molecule ATOM count array */

        VarArrayAdd(*vaPMolecules, (GENP) & iCount);
        if (!bSeenFirstSolvent) {
            if (cResidueType(srPRes->rResidue) == RESTYPESOLVENT) {
                bSeenFirstSolvent = TRUE;
                *iPFirstSolvent = iVarArrayElementCount(*vaPMolecules) - 1;
            }
        }

    }

    if (!bSeenFirstSolvent) {
        *iPFirstSolvent = iVarArrayElementCount(*vaPMolecules);
    }
}



/*
 *        zUnitDoResidue
 *
 *        Build a RESIDUE entry into the uUnit->vaResidues table.
 */
static void 
zUnitDoResidue(UNIT uUnit, RESIDUE rRes, int *iPPos)
{
    SAVERESIDUEt *srPResidue;

    ContainerSetTempInt(rRes, (*iPPos) + 1);

    srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, *iPPos);
    srPResidue->rResidue = rRes;
    REF(rRes);
    strcpy(srPResidue->sName, sContainerName(rRes));
    srPResidue->sResidueType[0] = cResidueType(rRes);
    srPResidue->sResidueType[1] = '\0';
    srPResidue->iSequenceNumber = iContainerSequence(rRes);
    srPResidue->iNextChildSequence = iContainerNextChildsSequence(rRes);
    srPResidue->iPdbSequenceNumber = iResiduePdbSequence(rRes);
    (*iPPos)++;

}


/*
 *        zUnitDoAtoms
 *
 *        Loop over all of the ATOMs within the RESIDUE and add them
 *        to the UNITs tables.
 *
 *        If the RESIDUE is of type RESTYPESOLVENT then the first
 *        ATOM will be the imaging atom, this is to maintain compatibility
 *        with AMBER parm files.
 */
void
zUnitDoAtoms(UNIT uUnit, PARMLIB plParameters, RESIDUE rRes, int *iPPos,
             BOOL * bPFailed, BOOL bPert)
{
    SAVERESIDUEt *srPResidue;
    LOOP lTemp;
    ATOM aAtom, aIgnore;
    BOOL bFailed;

    srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt,
                      iContainerTempInt(rRes) - 1);
    srPResidue->iAtomStartIndex = 0;

    /* Check if it is a solvent RESIDUE */
    /* If it is, put the imaging ATOM into the table first */

    aIgnore = NULL;
    if (cResidueType(rRes) == RESTYPESOLVENT) {
        if (aResidueImagingAtom(rRes) != NULL) {
            aIgnore = aResidueImagingAtom(rRes);
            zUnitIOTableAddAtom(uUnit, aIgnore, *iPPos,
                                plParameters, &bFailed, bPert);
            (*iPPos)++;
            *bPFailed |= bFailed;
            if (srPResidue->iAtomStartIndex == 0) {
                srPResidue->iAtomStartIndex = iContainerTempInt(aIgnore);
            }
        }
    }

    /* Now put the rest of the ATOMs */

    lTemp = lLoop((OBJEKT) rRes, DIRECTCONTENTSBYSEQNUM);
    while ((aAtom = (ATOM) oNext(&lTemp)) != NULL) {

        /* Ignore the imaging ATOM if there is one */

        if (aAtom == aIgnore)
            continue;

        zUnitIOTableAddAtom(uUnit, aAtom, *iPPos, plParameters, &bFailed,
                            bPert);
        (*iPPos)++;
        *bPFailed |= bFailed;
        if (srPResidue->iAtomStartIndex == 0) {
            srPResidue->iAtomStartIndex = iContainerTempInt(aAtom);
        }
    }
}


/* Put RESIDUEs in the following order: */
/* non-solvent residues */
/* solvent residues not in solvent cap */
/* solvent residues in solvent cap */
/* THIS IS ONLY DONE BECAUSE AMBER PARM FILES */
/* REQUIRE IT!!!!!!  PROGRAMS SHOULD ASSUME */
/* THAT RESIDUES HAVE ARBITRARY ORDERING WITHIN */
/* OBJECT FILE FORMAT FILES AND SORT */
/* THEM THEMSELVES TO PREVENT INCOMPATIBILITIES */
/* WITH FUTURE RELEASES!!!!!!!!!!!!!!!!!!!!!!!! */

int
zUnitIOAmberOrderResidues( UNIT uUnit )
{
        LOOP    lResidues;
        RESIDUE rRes;
        int     i = 0;
        int     iResidueCount = 0;

        /*
        **  count residues
        */
        lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
        while (oNext(&lResidues) != NULL)
                iResidueCount++;

        if (iResidueCount == 0)
                return(0);

        /*
        **  allocate array for residues
        */
        uUnit->vaResidues = vaVarArrayCreate(sizeof(SAVERESIDUEt));
        VarArraySetSize((uUnit->vaResidues), iResidueCount);

        /* Loop through solvent RESIDUEs and set the temp */
        /* flag saying if they are in the cap or not */

        lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
        while ((rRes = (RESIDUE) oNext(&lResidues))) {
            if (cResidueType(rRes) == RESTYPESOLVENT) {
                if (bUnitCapContainsContainer(uUnit, (CONTAINER) rRes)) {
                    ResidueSetFlags(rRes, RESIDUEINCAP);
                } else {
                    ResidueResetFlags(rRes, RESIDUEINCAP);
                }
            }
        }
        
        if ( !GDefaults.reorder_residues ) {
          printf("\"order_residues\" off: keep input residue order. This is risky: unkown behavior may ensue!\n");
          lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
          while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
              if (iObjectType(rRes) != RESIDUEid)
                  continue;
              zUnitDoResidue(uUnit, rRes, &i);
          }
        }

        else {
          //~ printf("\"reorder_residues\" on: move solvent residues to end.\n");
          /*
          **  solute residues
          */
          lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
          while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
              if (iObjectType(rRes) != RESIDUEid)
                  continue;
              if (cResidueType(rRes) == RESTYPESOLVENT)
                  continue;
              zUnitDoResidue(uUnit, rRes, &i);
          }
          /*
          **  solvent non-cap residues
          */
          lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
          while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
              if (iObjectType(rRes) != RESIDUEid)
                  continue;
              if (cResidueType(rRes) != RESTYPESOLVENT)
                  continue;
              if (bResidueFlagsSet(rRes, RESIDUEINCAP))
                  continue;
              zUnitDoResidue(uUnit, rRes, &i);
          }
          /*
          **  solvent cap residues
          */
          lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
          while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
              if (iObjectType(rRes) != RESIDUEid)
                  continue;
              if (cResidueType(rRes) != RESTYPESOLVENT)
                  continue;
              if (bResidueFlagsSet(rRes, RESIDUEINCAP))
                  zUnitDoResidue(uUnit, rRes, &i);
          }
        }
        return(iResidueCount);
}

/*
 *      zUnitIOBuildTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Build a table representation of the UNIT.
 *        If (bPert) then build tables for perturbation run.
 *
 *        NOTE: Programmers should assume that the order of entries
 *        NOTE: within tables and OFF files is COMPLETELY ARBITRARY.
 */
void
zUnitIOBuildTables(UNIT uUnit, PARMLIB plParameters,
                   BOOL * bPGeneratedParameters, BOOL bPert, BOOL bCheck)
{
    SAVECONNECTIVITYt *scPCon;
    SAVEMOLECULEt *smPMolecule;
    SAVERESTRAINTt *srPRestraint;
    SAVERESIDUEt *srPResidue;
    SAVEHIERARCHYt *shPHierarchy;
    SAVEGROUPSt sgAtomGroup;
    BAGLOOP blRestraint;
    RESTRAINT rRest;
    MOLECULE mMol;
    RESIDUE rRes;
    ATOM aAtom1, aAtom2, aAtom3, aAtom4;
    LOOP lMolecules;
    OBJEKT oAbove, oBelow;
    STRING sAtom1, sAtom2, sDesc;
    BOOL bGenerateParameters, bFailedGeneratingParameters;
    PARMSET psTemp;
    LOOP lTemp, lResidues;
    DICTLOOP dlGroups;
    LISTLOOP llAtoms;
    LIST lGroup;
    int i, j, iAtomCount, iMoleculeCount, iResidueCount;
    int iIndex, iCount, iGroup, iErrors = 0, iWarnings = 0;
    double dKx, dX0, dN, dA, dB, dEI, dEJ, dRI, dRJ;

    if (bCheck) {
        VP0(("Checking Unit.\n"));
        iFatal = 0;
        UnitCheck(uUnit, &iErrors, &iWarnings);
        if (iFatal) {
            VP0(("Failed to generate parameters\n"));
            *bPGeneratedParameters = FALSE;
            return;
        }
        if (iErrors || iWarnings) {
            /* just for fun, grammar */
            VP0(("\n -- ignoring the %s%s%s%s%s.\n\n",
                 (iErrors ? "error" : ""),
                 (iErrors > 1 ? "s" : ""),
                 (iErrors && iWarnings ? " and " : ""),
                 (iWarnings ? "warning" : ""),
                 (iWarnings > 1 ? "s" : "")));
        }
    }

    VP0(("Building topology.\n"));

    bFailedGeneratingParameters = FALSE;
    if (iUnitMode(uUnit) != UNITNORMAL) {
        DFATAL(("The UNIT must be in NORMAL mode!"));
    }

    UnitSetMode(uUnit, UNITTABLES);

    bGenerateParameters = FALSE;
    if (plParameters != NULL) {
        uUnit->psParameters = (PARMSET) oCreate(PARMSETid);
        bGenerateParameters = TRUE;
    }


    /* Build the molecule information */

    iMoleculeCount = 0;
    lMolecules = lLoop((OBJEKT) uUnit, MOLECULES);
    while (oNext(&lMolecules) != NULL)
        iMoleculeCount++;

    if (iMoleculeCount) {
        uUnit->vaMolecules = vaVarArrayCreate(sizeof(SAVEMOLECULEt));
        VarArraySetSize((uUnit->vaMolecules), iMoleculeCount);

        lMolecules = lLoop((OBJEKT) uUnit, MOLECULES);
        smPMolecule = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        for (i = 0; (mMol = (MOLECULE) oNext(&lMolecules));
             i++, smPMolecule++) {
            ContainerSetTempInt(mMol, i + 1);
            smPMolecule->mMolecule = mMol;
            REF(mMol);
            strcpy(smPMolecule->sName, sContainerName(mMol));
            smPMolecule->iSequenceNumber = iContainerSequence(mMol);
            smPMolecule->iNextChildSequence =
                iContainerNextChildsSequence(mMol);
        }
    }


    /* Build the residue information */

    /* zUnitIOAmberOrderResidues puts residues in arbitrary amber order */
    iResidueCount = zUnitIOAmberOrderResidues( uUnit );

    if (iResidueCount) {

        /* Build the array for the atoms */

        VP0(("Building atom parameters.\n"));

        iAtomCount = 0;
        lTemp = lLoop((OBJEKT) uUnit, ATOMS);
        while (oNext(&lTemp) != NULL)
            iAtomCount++;

        if (iAtomCount) {
            uUnit->vaAtoms = vaVarArrayCreate(sizeof(SAVEATOMt));
            VarArraySetSize((uUnit->vaAtoms), iAtomCount);
            i = 0;

            /* Put solvent ATOMs after other residues */
            /* THIS IS ONLY DONE BECAUSE AMBER PARM FILES */
            /* REQUIRE IT!!!!!!  PROGRAMS SHOULD ASSUME */
            /* THAT ATOMS HAVE ARBITRARY ORDERING AND SORT */
            /* THEM THEMSELVES TO PREVENT INCOMPATIBILITIES */
            /* WITH FUTURE RELEASES!!!!!!!!!!!!!!!!!!!!!!!! */

          if ( !GDefaults.reorder_residues ) {
            lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
            uUnit->iCapTempInt=0;
            while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
              if ( !uUnit->iCapTempInt ) {
                if (cResidueType(rRes) == RESTYPESOLVENT &&
                    bResidueFlagsSet(rRes, RESIDUEINCAP))
                  uUnit->iCapTempInt = i;
              }
              zUnitDoAtoms(uUnit, plParameters, rRes, &i, 
                           &bFailedGeneratingParameters, bPert);                                 
            }
          }
          else {
            lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
            while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
                if (cResidueType(rRes) != RESTYPESOLVENT)
                    zUnitDoAtoms(uUnit, plParameters, rRes,
                                 &i, &bFailedGeneratingParameters, bPert);
            }
            lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
            while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
                if (cResidueType(rRes) == RESTYPESOLVENT &&
                    !bResidueFlagsSet(rRes, RESIDUEINCAP))
                    zUnitDoAtoms(uUnit, plParameters, rRes,
                                 &i, &bFailedGeneratingParameters, bPert);
            }

            /* Set the uUnit->iCapTempInt integer to point to */
            /* the last ATOM that is not CAP solvent */

            uUnit->iCapTempInt = i;
            lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
            while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
                if (cResidueType(rRes) == RESTYPESOLVENT &&
                    bResidueFlagsSet(rRes, RESIDUEINCAP))
                    zUnitDoAtoms(uUnit, plParameters, rRes,
                                 &i, &bFailedGeneratingParameters, bPert);
            }
          }
        }
    }

/*
 *. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
 *
 *        Now the ATOMs have their indices stored in ContainerTempInt
 *
 */

    /* Generate the table storing the RESTRAINTs */

    if (iBagSize(uUnit->bRestraints)) {
        uUnit->vaRestraints = vaVarArrayCreate(sizeof(SAVERESTRAINTt));
        VarArraySetSize((uUnit->vaRestraints),
                        iBagSize(uUnit->bRestraints));
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        blRestraint = blBagLoop(uUnit->bRestraints);
        while ((rRest = (RESTRAINT) PBagNext(&blRestraint))) {
            srPRestraint->iType = iRestraintType(rRest);
            srPRestraint->fFlags = fRestraintFlags(rRest);
            switch (iRestraintType(rRest)) {
            case RESTRAINTBOND:
                RestraintBondGet(rRest, &aAtom1, &aAtom2, &dKx, &dX0);
                srPRestraint->iAtom1 = iContainerTempInt(aAtom1);
                srPRestraint->iAtom2 = iContainerTempInt(aAtom2);
                srPRestraint->iAtom3 = 0;
                srPRestraint->iAtom4 = 0;
                srPRestraint->dKx = dKx;
                srPRestraint->dX0 = dX0;
                srPRestraint->dN = 0.0;
                break;
            case RESTRAINTANGLE:
                RestraintAngleGet(rRest, &aAtom1, &aAtom2, &aAtom3,
                                  &dKx, &dX0);
                srPRestraint->iAtom1 = iContainerTempInt(aAtom1);
                srPRestraint->iAtom2 = iContainerTempInt(aAtom2);
                srPRestraint->iAtom3 = iContainerTempInt(aAtom3);
                srPRestraint->iAtom4 = 0;
                srPRestraint->dKx = dKx;
                srPRestraint->dX0 = dX0;
                srPRestraint->dN = 0.0;
                break;
            case RESTRAINTTORSION:
                RestraintTorsionGet(rRest, &aAtom1, &aAtom2,
                                    &aAtom3, &aAtom4, &dKx, &dX0, &dN);
                srPRestraint->iAtom1 = iContainerTempInt(aAtom1);
                srPRestraint->iAtom2 = iContainerTempInt(aAtom2);
                srPRestraint->iAtom3 = iContainerTempInt(aAtom3);
                srPRestraint->iAtom4 = iContainerTempInt(aAtom4);
                srPRestraint->dKx = dKx;
                srPRestraint->dX0 = dX0;
                srPRestraint->dN = dN;
                break;
            default:
                DFATAL(("Illegal restraint type!"));
            }
            srPRestraint++;
        }
    }


    /* Generate the indices for the RESIDUE connect atoms */
    /* and Imaging ATOMs */

    if (iVarArrayElementCount(uUnit->vaResidues)) {
        srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (i = 0; i < iResidueCount; srPResidue++, i++) {
            srPResidue->iImagingAtomIndex = 0;
            if (aResidueImagingAtom(srPResidue->rResidue) != NULL) {
                srPResidue->iImagingAtomIndex =
                    iContainerTempInt(aResidueImagingAtom
                                      (srPResidue->rResidue));
            }
            for (j = 0; j < MAXCONNECT; j++) {
                if (bResidueConnectUsed(srPResidue->rResidue, j)) {
                    srPResidue->iaConnectIndex[j] =
                        iContainerTempInt(srPResidue->rResidue->
                                          aaConnect[j]);
                } else
                    srPResidue->iaConnectIndex[j] = 0;
            }
        }
    }

    /* Now generate the connectivity table */

    iCount = 0;
    lTemp = lLoop((OBJEKT) uUnit, BONDS);
    while (oNext(&lTemp) != NULL)
        iCount++;
    uUnit->vaConnectivity = vaVarArrayCreate(sizeof(SAVECONNECTIVITYt));
    VarArraySetSize((uUnit->vaConnectivity), iCount);
    if (iCount) {
        i = 0;
        lTemp = lLoop((OBJEKT) uUnit, BONDS);
        scPCon = PVAI(uUnit->vaConnectivity, SAVECONNECTIVITYt, 0);
        while (oNext(&lTemp) != NULL) {
            LoopGetBond(&lTemp, &aAtom1, &aAtom2);
            scPCon->fFlags = fAtomFindBondFlags(aAtom1, aAtom2);
            scPCon->iAtom1 = iContainerTempInt(aAtom1);
            scPCon->iAtom2 = iContainerTempInt(aAtom2);
            i++;
            scPCon++;
        }
    }


    /* Build a table for the Hierarchy information */

    iCount = 0;
    lTemp = lLoop((OBJEKT) uUnit, CONTAINERS);
    while (oNext(&lTemp) != NULL)
        iCount++;
    if (iCount) {
        uUnit->vaHierarchy = vaVarArrayCreate(sizeof(SAVEHIERARCHYt));
        VarArraySetSize((uUnit->vaHierarchy), iCount);
        lTemp = lLoop((OBJEKT) uUnit, CONTAINERS);
        i = 0;
        shPHierarchy = PVAI(uUnit->vaHierarchy, SAVEHIERARCHYt, 0);
        while ((oBelow = oNext(&lTemp)) != NULL) {
            oAbove = (OBJEKT) cContainerWithin(oBelow);
            shPHierarchy->sAboveType[0] = iObjectType(oAbove);
            shPHierarchy->sAboveType[1] = '\0';
            shPHierarchy->sBelowType[0] = iObjectType(oBelow);
            shPHierarchy->sBelowType[1] = '\0';
            if (iObjectType(oAbove) == UNITid) {
                shPHierarchy->iAboveIndex = 0;
            } else {
                shPHierarchy->iAboveIndex = iContainerTempInt(oAbove);
            }
            shPHierarchy->iBelowIndex = iContainerTempInt(oBelow);
            i++;
            shPHierarchy++;
        }
    }

    /* Build a table for the UNIT Connect information */

    uUnit->vaConnect = vaVarArrayCreate(sizeof(int));
    VarArraySetSize((uUnit->vaConnect), 2);
    if (bUnitHeadUsed(uUnit))
        iIndex = iContainerTempInt(aUnitHead(uUnit));
    else
        iIndex = 0;
    *PVAI(uUnit->vaConnect, int, 0) = iIndex;
    if (bUnitTailUsed(uUnit))
        iIndex = iContainerTempInt(aUnitTail(uUnit));
    else
        iIndex = 0;
    *PVAI(uUnit->vaConnect, int, 1) = iIndex;



    /* Build tables for the UNIT groups */

    if (iDictionaryElementCount(uUnit->dAtomGroups)) {
        uUnit->vaGroupNames = vaVarArrayCreate(sizeof(STRING));
        VarArraySetSize(uUnit->vaGroupNames,
                        iDictionaryElementCount(uUnit->dAtomGroups));

        /* assuming there are atoms; the array is grown on the fly */
        uUnit->vaGroupAtoms = vaVarArrayCreate(sizeof(SAVEGROUPSt));
        iGroup = 0;
        dlGroups = ydlDictionaryLoop(uUnit->dAtomGroups);
        while ((lGroup =
               (LIST) yPDictionaryNext(uUnit->dAtomGroups, &dlGroups))) {
            strcpy((char *) (PVAI(uUnit->vaGroupNames, STRING, iGroup)),
                   sDictLoopKey(dlGroups));
            llAtoms = llListLoop(lGroup);
            for (i = 0; i < iListSize(lGroup); i++) {
                sgAtomGroup.iGroupIndex = iGroup + 1;
                sgAtomGroup.iIndexAtom =
                    iContainerTempInt((ATOM) oListNext(&llAtoms));
                VarArrayAdd(uUnit->vaGroupAtoms, (GENP) & sgAtomGroup);
            }
            iGroup++;
        }
    }

/* The rest of the code should be executed ONLY if PARAMETERS are being */
/* generated */

    if (bFailedGeneratingParameters == FALSE && bGenerateParameters) {


        /* Now generate the BOND table */

        bFailedGeneratingParameters |=
            zbUnitIOIndexBondParameters(plParameters, uUnit, bPert);


        /* Now generate the ANGLE table */

        bFailedGeneratingParameters |=
            zbUnitIOIndexAngleParameters(plParameters, uUnit, bPert);


        /* Now generate the TORSION and IMPROPER tables */
        /* Place them in the SAME VARARRAY !!!!! */

        if (uUnit->vaTorsions != NULL) {
            VP0(("Regenerating proper and improper torsions.\n"));
            VarArrayDestroy(&(uUnit->vaTorsions));
        }
        uUnit->vaTorsions = vaVarArrayCreate(sizeof(SAVETORSIONt));
        bFailedGeneratingParameters |=
            zbUnitIOIndexTorsionParameters(plParameters, uUnit, TRUE, bPert);
        bFailedGeneratingParameters |=
            zbUnitIOIndexTorsionParameters(plParameters, uUnit, FALSE, bPert);


        /* Generate the potential H-Bond parameters             */
        /* Do this by looping through ALL H-Bond parameters     */
        /* looking for those where both atoms are defined       */
        /* within this UNITs atom list                          */
        /* If there is one then add it to this UNITs PARMSET    */

        VP0(("Building H-Bond parameters.\n"));
        PARMLIB_LOOP_ALL(plParameters, psTemp) {
            for (i = 0; i < iParmSetTotalHBondParms(psTemp); i++) {
                ParmSetHBond(psTemp, i, sAtom1, sAtom2, &dA, &dB, sDesc);
                if ((iParmSetFindAtom(uUnit->psParameters, sAtom1)
                     != PARM_NOT_FOUND) &&
                    (iParmSetFindAtom(uUnit->psParameters, sAtom2)
                     != PARM_NOT_FOUND)) {
                    if (iParmSetFindHBond
                        (uUnit->psParameters, sAtom1,
                         sAtom2) == PARM_NOT_FOUND) {
                        iParmSetAddHBond(uUnit->psParameters, sAtom1,
                                         sAtom2, dA, dB, sDesc);
                    }
                }
            }
        }

	/* Carry LJ edits into the unit topology */

	VP0(("Incorporating Non-Bonded adjustments.\n"));
	PARMLIB_LOOP_ALL(plParameters, psTemp) {
	  for (i = 0; i < iParmSetTotalNBEdits(psTemp); i++) {
	    ParmSetNBEdit(psTemp, i, sAtom1, sAtom2, &dEI, &dEJ, &dRI, &dRJ,
			  sDesc);
	    if ((iParmSetFindAtom(uUnit->psParameters, sAtom1)
		 != PARM_NOT_FOUND) &&
		(iParmSetFindAtom(uUnit->psParameters, sAtom2)
		 != PARM_NOT_FOUND)) {
	      if (iParmSetFindNBEdit(uUnit->psParameters, sAtom1, sAtom2) ==
		  PARM_NOT_FOUND) {
		iParmSetAddNBEdit(uUnit->psParameters, sAtom1, sAtom2, dEI,
				  dEJ, dRI, dRJ, sDesc);
	      }
	    }
	  }
	}
    }
    if (bFailedGeneratingParameters == FALSE)
        *bPGeneratedParameters = TRUE;
    else
        *bPGeneratedParameters = FALSE;
}






/*
 *      zUnitIOBuildFromTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Build a UNIT from its tables.
 */
void zUnitIOBuildFromTables(UNIT uUnit)
{
    SAVEATOMt *saPAtom;
    SAVECONNECTIVITYt *scPCon;
    SAVERESTRAINTt *srPRestraint;
    SAVERESIDUEt *srPResidue;
    SAVEMOLECULEt *smPMolecule;
    SAVEHIERARCHYt *shPHierarchy;
    SAVEGROUPSt *sgPGroupAtom;
    RESTRAINT rRest;
    ATOM aAtom, aAtom1, aAtom2, aAtom3, aAtom4;
    RESIDUE rRes;
    MOLECULE mMol;
    OBJEKT oObj1, oObj2;
    int i, j, iMax, iIndex, iChildSeq, iSeq, iGroup;
    STRING sGroup;
    LIST lGroup;

    /* Build the atoms */

    iMax = iVarArrayElementCount(uUnit->vaAtoms);
    if (!iMax)
        return;

    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    for (i = 0; i < iMax; i++, saPAtom++) {
        aAtom = (ATOM) oCreate(ATOMid);
        ContainerSetName(aAtom, saPAtom->sName);
        ContainerSetSequence(aAtom, saPAtom->iSequence);
        AtomSetElement(aAtom, saPAtom->iElement);
        AtomSetPertElement(aAtom, saPAtom->iPertElement);
        AtomSetPertName(aAtom, saPAtom->sPertName);
        AtomSetType(aAtom, saPAtom->sType);
        AtomSetPertType(aAtom, saPAtom->sPertType);
        AtomSetCharge(aAtom, saPAtom->dCharge);
        AtomSetPertCharge(aAtom, saPAtom->dPertCharge);
        AtomSetPosition(aAtom, saPAtom->vPos);
        AtomSetVelocity(aAtom, saPAtom->vVelocity);
        AtomDefineFlags(aAtom, saPAtom->fFlags);
        /*
         *  the ATOM objekt has 1 reference at this point..
         */
        saPAtom->aAtom = aAtom;
    }

    /* Build the Residues */

    srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
    iMax = iVarArrayElementCount(uUnit->vaResidues);
    for (i = 0; i < iMax; i++, srPResidue++) {
        rRes = (RESIDUE) oCreate(RESIDUEid);
        ContainerSetName(rRes, srPResidue->sName);
        ContainerSetSequence(rRes, srPResidue->iSequenceNumber);
        ContainerSetNextChildsSequence(rRes,
                                       srPResidue->iNextChildSequence);
        ResidueSetPdbSequence(rRes, srPResidue->iPdbSequenceNumber);
        ResidueSetType(rRes, srPResidue->sResidueType[0]);

        /* Define the imaging ATOM */

        iIndex = srPResidue->iImagingAtomIndex;
        if (iIndex != 0)
            aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iIndex - 1)->aAtom;
        else
            aAtom = NULL;
        ResidueSetImagingAtom(rRes, aAtom);

        for (j = 0; j < MAXCONNECT; j++) {
            iIndex = srPResidue->iaConnectIndex[j];
            if (iIndex != 0) {
                aAtom = (PVAI(uUnit->vaAtoms, SAVEATOMt,
                              iIndex - 1))->aAtom;
                ResidueSetConnectAtom(rRes, j, aAtom);
            }
        }
        /*
         *  the RESIDUE objekt has 1 reference at this point..
         */
        srPResidue->rResidue = rRes;
    }

    /* Build the Molecules */

    if ((iMax = iVarArrayElementCount(uUnit->vaMolecules))) {
        smPMolecule = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        for (i = 0; i < iMax; i++, smPMolecule++) {
            mMol = (MOLECULE) oCreate(MOLECULEid);
            ContainerSetName(mMol, smPMolecule->sName);
            ContainerSetSequence(mMol, smPMolecule->iSequenceNumber);
            ContainerSetNextChildsSequence(rRes,
                                           smPMolecule->
                                           iNextChildSequence);
            smPMolecule->mMolecule = mMol;
        }
    }

    /* Create bonds between atoms */

    if ((iMax = iVarArrayElementCount(uUnit->vaConnectivity))) {
        for (i = 0; i < iMax; i++) {
            scPCon = PVAI(uUnit->vaConnectivity, SAVECONNECTIVITYt, i);
            aAtom1 =
                PVAI(uUnit->vaAtoms, SAVEATOMt, scPCon->iAtom1 - 1)->aAtom;
            aAtom2 =
                PVAI(uUnit->vaAtoms, SAVEATOMt, scPCon->iAtom2 - 1)->aAtom;
            AtomBondToFlags(aAtom1, aAtom2, scPCon->fFlags);
        }
    }

    /* Build the hierarchy of Unit/Molecules/Residues/Atoms */

    if ((iMax = iVarArrayElementCount(uUnit->vaHierarchy))) {
        shPHierarchy = PVAI(uUnit->vaHierarchy, SAVEHIERARCHYt, 0);
        for (i = 0; i < iMax; i++, shPHierarchy++) {

            switch (shPHierarchy->sAboveType[0]) {
            case UNITid:
                oObj1 = (OBJEKT) uUnit;
                break;
            case MOLECULEid:
                oObj1 = (OBJEKT)
                    (PVAI(uUnit->vaMolecules, SAVEMOLECULEt,
                          shPHierarchy->iAboveIndex - 1))->mMolecule;
                break;
            case RESIDUEid:
                oObj1 = (OBJEKT)
                    (PVAI(uUnit->vaResidues, SAVERESIDUEt,
                          shPHierarchy->iAboveIndex - 1))->rResidue;
                break;
            case ATOMid:
                oObj1 = (OBJEKT)
                    (PVAI(uUnit->vaAtoms, SAVEATOMt,
                          shPHierarchy->iAboveIndex - 1))->aAtom;
                break;
            }
            switch (shPHierarchy->sBelowType[0]) {
            case UNITid:
                oObj2 = (OBJEKT) uUnit;
                break;
            case MOLECULEid:
                oObj2 = (OBJEKT)
                    (PVAI(uUnit->vaMolecules, SAVEMOLECULEt,
                          shPHierarchy->iBelowIndex - 1))->mMolecule;
                break;
            case RESIDUEid:
                oObj2 = (OBJEKT)
                    (PVAI(uUnit->vaResidues, SAVERESIDUEt,
                          shPHierarchy->iBelowIndex - 1))->rResidue;
                break;
            case ATOMid:
                oObj2 = (OBJEKT)
                    (PVAI(uUnit->vaAtoms, SAVEATOMt,
                          shPHierarchy->iBelowIndex - 1))->aAtom;
                break;
            }
            iChildSeq = iContainerNextChildsSequence(oObj1);
            iSeq = iContainerSequence(oObj2);
            /*
             *  ContainerAdd() increments Obj2's reference counter
             */
            ContainerAdd((CONTAINER) oObj1, oObj2);
            ContainerSetNextChildsSequence(oObj1, iChildSeq);
            ContainerSetSequence(oObj2, iSeq);
        }
    }

    /* Set up the RESTRAINTs */

    if ((iMax = iVarArrayElementCount(uUnit->vaRestraints))) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        for (i = 0; i < iMax; i++, srPRestraint++) {
            rRest = rRestraintCreate();
            RestraintDefineFlags(rRest, srPRestraint->fFlags);
            switch (srPRestraint->iType) {
            case RESTRAINTBOND:
                aAtom1 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom1 - 1)->aAtom;
                aAtom2 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom2 - 1)->aAtom;
                RestraintBondSet(rRest, aAtom1, aAtom2,
                                 srPRestraint->dKx, srPRestraint->dX0);
                break;
            case RESTRAINTANGLE:
                aAtom1 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom1 - 1)->aAtom;
                aAtom2 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom2 - 1)->aAtom;
                aAtom3 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom3 - 1)->aAtom;
                RestraintAngleSet(rRest, aAtom1, aAtom2, aAtom3,
                                  srPRestraint->dKx, srPRestraint->dX0);
                break;
            case RESTRAINTTORSION:
                aAtom1 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom1 - 1)->aAtom;
                aAtom2 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom2 - 1)->aAtom;
                aAtom3 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom3 - 1)->aAtom;
                aAtom4 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom4 - 1)->aAtom;
                RestraintTorsionSet(rRest, aAtom1, aAtom2, aAtom3, aAtom4,
                                    srPRestraint->dKx, srPRestraint->dX0,
                                    srPRestraint->dN);
                break;
            default:
                DFATAL(("Illegal RESTRAINT type loaded"));
            }
            UnitAddRestraint(uUnit, rRest);
        }
    }

    /* Set up the UNIT Connect information */

    iIndex = *PVAI(uUnit->vaConnect, int, 0);
    aAtom = NULL;
    if (iIndex != 0)
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iIndex - 1)->aAtom;
    UnitSetHead(uUnit, aAtom);

    iIndex = *PVAI(uUnit->vaConnect, int, 1);
    aAtom = NULL;
    if (iIndex != 0)
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iIndex - 1)->aAtom;
    UnitSetTail(uUnit, aAtom);

    /* Set up the ATOM groups information */

    if ((iMax = iVarArrayElementCount(uUnit->vaGroupNames))) {
        for (i = 0; i < iMax; i++) {
            strcpy(sGroup,
                   (char *) (PVAI(uUnit->vaGroupNames, STRING, i)));
            bUnitGroupCreate(uUnit, sGroup);
        }
    }
    if ((iMax = iVarArrayElementCount(uUnit->vaGroupAtoms))) {
        sgPGroupAtom = PVAI(uUnit->vaGroupAtoms, SAVEGROUPSt, 0);
        for (i = 0; i < iMax; i++, sgPGroupAtom++) {
            iGroup = sgPGroupAtom->iGroupIndex - 1;
            iIndex = sgPGroupAtom->iIndexAtom;
            lGroup = lUnitGroup(uUnit,
                                (char *) PVAI(uUnit->vaGroupNames, STRING,
                                              iGroup));
            aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iIndex - 1)->aAtom;
            ListAddUnique(lGroup, (GENP) aAtom);
        }
    }

}






/*
 *      zUnitIODestroyTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Destroy all of the tables associated with the UNIT.
 */
void zUnitIODestroyTables(UNIT uUnit)
{
    int i, iCount;

    /*
     *  since the SAVEATOMt's have pointers to the ATOMs,
     *      need to deref the ATOMs so that they will
     *      be freed when the last pointer to them is
     *      deleted. Same goes for RESIDUEs..
     */
    if (uUnit->vaAtoms != NULL) {
        SAVEATOMt *aP = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
        iCount = iVarArrayElementCount(uUnit->vaAtoms);
        for (i = 0; i < iCount; i++, aP++)
            DEREF(aP->aAtom);
        VarArrayDestroy(&(uUnit->vaAtoms));
    }
    if (uUnit->vaResidues != NULL) {
        SAVERESIDUEt *rP = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        iCount = iVarArrayElementCount(uUnit->vaResidues);
        for (i = 0; i < iCount; i++, rP++)
            DEREF(rP->rResidue);
        VarArrayDestroy(&(uUnit->vaResidues));
    }

    if (uUnit->vaMolecules != NULL) {
        SAVEMOLECULEt *mP = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        iCount = iVarArrayElementCount(uUnit->vaMolecules);
        for (i = 0; i < iCount; i++, mP++)
            DEREF(mP->mMolecule);
        VarArrayDestroy(&(uUnit->vaMolecules));
    }
/* %%% */
    if (uUnit->vaHierarchy != NULL)
        VarArrayDestroy(&(uUnit->vaHierarchy));
    if (uUnit->vaBonds != NULL)
        VarArrayDestroy(&(uUnit->vaBonds));
    if (uUnit->vaAngles != NULL)
        VarArrayDestroy(&(uUnit->vaAngles));
    if (uUnit->vaTorsions != NULL)
        VarArrayDestroy(&(uUnit->vaTorsions));
    if (uUnit->vaRestraints != NULL)
        VarArrayDestroy(&(uUnit->vaRestraints));
    if (uUnit->vaConnect != NULL)
        VarArrayDestroy(&(uUnit->vaConnect));
    if (uUnit->vaConnectivity != NULL)
        VarArrayDestroy(&(uUnit->vaConnectivity));
    if (uUnit->vaGroupNames != NULL)
        VarArrayDestroy(&(uUnit->vaGroupNames));
    if (uUnit->vaGroupAtoms != NULL)
        VarArrayDestroy(&(uUnit->vaGroupAtoms));

    UnitSetMode(uUnit, UNITNORMAL);
}

/*
 *      CheckTypeNBEdit
 *
 *        Author:       David S. Cerutti (2013)
 *
 *      Check to see whether this atom type (which is imminently to be
 *      compacted in the nonbonded parameters array) is mentioned in any
 *      nonbonded pair potential adjustments.  If so, then it has to remain
 *      its own type.
 */
static int
CheckTypeNBEdit(typeStr sType, VARARRAY vaPNBEdits)
{
  int i;

  for (i = 0; i < vaPNBEdits->count; i++) {
    if (strcmp(sType, PVAI(vaPNBEdits, NBEDITt, i)->sType1) == 0 ||
	strcmp(sType, PVAI(vaPNBEdits, NBEDITt, i)->sType2) == 0) {
      return 1;
    }
  }

  return 0;
}

/*
 *      CheckAgainstNBEdits
 *
 *        Author:       David S. Cerutti (2013)
 *
 *      Check to see whether an atom type, when paired with any other atom
 *      type, has any Lennard-Jones sigma and epsilon pairs that do not
 *      conform to the standard combining rules.
 */
static void
CheckAgainstNBEdits(VARARRAY vaPNBEdits, typeStr tI, typeStr tJ,
		    double *dA, double *dC)
{
  int i;

  for (i = 0; i < vaPNBEdits->count; i++) {
    if ((strcmp(PVAI(vaPNBEdits, NBEDITt, i)->sType1, tI) == 0 &&
	 strcmp(PVAI(vaPNBEdits, NBEDITt, i)->sType2, tJ) == 0) ||
	(strcmp(PVAI(vaPNBEdits, NBEDITt, i)->sType1, tJ) == 0 &&
	 strcmp(PVAI(vaPNBEdits, NBEDITt, i)->sType2, tI) == 0)) {

      /* This is an edited pair interaction */
      MathOpConvertNonBondToAC(PVAI(vaPNBEdits, NBEDITt, i)->dEI,
			       PVAI(vaPNBEdits, NBEDITt, i)->dRI,
			       PVAI(vaPNBEdits, NBEDITt, i)->dEJ,
			       PVAI(vaPNBEdits, NBEDITt, i)->dRJ, dA, dC);
      break;
    }
  }
}



/*
 *      zUnitIOBuildNonBondArrays
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Build the two NON-BOND arrays.  One is NxN where N is the total
 *      number of types, and the other is Mx(M+1)/2 where M is the
 *      number of UNIQUE non-bond types.
 *      The NxN (vaNBIndexMatrix) array is square matrix containing 
 *      integer indices+1 into the Mx(M+1)/2 (vaNBParameters) array
 *      which contains the NON-BOND A and C parameters.
 *      (vaPNBIndex) will return a pointer to a VarArray which contains
 *      indices into the (vaPNonBonds) array which contains the indices
 *      for the unique atom types (all those with unique parameters).
 */
static void
zUnitIOBuildNonBondArrays(UNIT uUnit, VARARRAY * vaPNBIndexMatrix,
                          VARARRAY * vaPNBParameters,
                          VARARRAY * vaPNBIndex, VARARRAY * vaPNonBonds)
{
    VARARRAY vaNBIndex, vaNonBonds, vaPNBEdits;
    int i, j, iNonBonds, iNBIndices, iTemp, iI, iJ, iX, iY;
    int iIndex, iHBondIndex, iNBIndex, iElement, iHybridization;
    double dMass, dPolar, dE, dR, dE14, dR14, dA, dC, dEI, dRI, dEJ, dRJ;
    double dScreenF;
    STRING sType, sTypeI, sTypeJ, sDesc;
    NONBONDt *nbPFirst, *nbPCur, *nbPTemp, *nbPLast;

    /* Build the NON-BONDED arrays */
    /* First reduce the non-bond parameters to the absolute */
    /* minimum, by rolling up parameters with duplicate values */
    /* maintaining a many to one mapping into the */
    /* rolled up non-bond parameter array */


    vaNBIndex = vaVarArrayCreate(sizeof(int));
    vaNonBonds = vaVarArrayCreate(sizeof(NONBONDt));
    iNonBonds = iParmSetTotalAtomParms(uUnit->psParameters);
    iNBIndices = iNonBonds;
    VarArraySetSize(vaNonBonds, iNonBonds);
    VarArraySetSize(vaNBIndex, iNonBonds);
    vaPNBEdits = uUnit->psParameters->vaNBEdits;
    for (i = 0; i < iNonBonds; i++) {
        ParmSetAtom(uUnit->psParameters, i,
                    sType, &dMass, &dPolar, &dE, &dR, &dE14, &dR14,
                    &dScreenF, &iElement, &iHybridization, sDesc);
        PVAI(vaNonBonds, NONBONDt, i)->bCapableOfHBonding =
            bParmSetCapableOfHBonding(uUnit->psParameters, sType);
        PVAI(vaNonBonds, NONBONDt, i)->dE = dE;
        PVAI(vaNonBonds, NONBONDt, i)->dR = dR;
        PVAI(vaNonBonds, NONBONDt, i)->dE14 = dE14;
        PVAI(vaNonBonds, NONBONDt, i)->dR14 = dR14;
        strcpy(PVAI(vaNonBonds, NONBONDt, i)->sType, sType);
        *PVAI(vaNBIndex, int, i) = i;
    }

    /* Now roll up the equivalent NON-BOND parameters */

    if (!GDefaults.iCharmm) {
        nbPLast = PVAI(vaNonBonds, NONBONDt, iNonBonds - 1);
        for (iNBIndex = 0; iNBIndex < iVarArrayElementCount(vaNBIndex);
             iNBIndex++) {
            nbPFirst = PVAI(vaNonBonds, NONBONDt, 0);
            nbPCur =
                PVAI(vaNonBonds, NONBONDt, *PVAI(vaNBIndex, int, iNBIndex));
            if (!nbPCur->bCapableOfHBonding) {
                while (nbPFirst < nbPCur) {
                    if (nbPFirst->dE == nbPCur->dE &&
                        nbPFirst->dR == nbPCur->dR &&
                        nbPFirst->dE14 == nbPCur->dE14 &&
                        nbPFirst->dR14 == nbPCur->dR14 &&
			CheckTypeNBEdit(nbPFirst->sType, vaPNBEdits) == 0 &&
			CheckTypeNBEdit(nbPCur->sType, vaPNBEdits) == 0) {
                        for (nbPTemp = nbPCur; nbPTemp < nbPLast;
                             nbPTemp++) *nbPTemp = *(nbPTemp + 1);
                       iNonBonds--;
                        nbPLast--;

                        /* Update the indices into the NON-BOND array */
                        /* Making sure that indices into parameters that */
                        /* follow the one at j will be moved back on */

                        j = nbPFirst - PVAI(vaNonBonds, NONBONDt, 0);
                        *PVAI(vaNBIndex, int, iNBIndex) = j;
                        for (iTemp = iNBIndex + 1; iTemp < iNBIndices;
                             iTemp++) {
                            if (*PVAI(vaNBIndex, int, iTemp) > j)
                                 (*PVAI(vaNBIndex, int, iTemp))--;
                        }
                        break;
                    }
                    nbPFirst++;
                }
            }
        }
    }

    /* Now the first iNonBonds entries of the array vaNonBonds */
    /* contain unique NON-BOND parameters, and the array vaNBIndex */
    /* contains indices into the vaNonBonds array for every */
    /* NON-BOND type */
    /* Change the size of the vaNonBonds array to the number of */
    /* valid non-bond parameters */
    VarArraySetSize(vaNonBonds, iNonBonds);

    /* Build the NxN integer index array */

    *vaPNBIndexMatrix = vaVarArrayCreate(sizeof(int));
    VarArraySetSize((*vaPNBIndexMatrix), iNonBonds * iNonBonds);
    for (i = 0; i < iNBIndices; i++) {
        for (j = 0; j < iNBIndices; j++) {

            /* Calculate the position of the parameters for the */
            /* non-bond interaction i-j within the vaPNBParameters */
            /* array */

            iI = *PVAI(vaNBIndex, int, i);
            iJ = *PVAI(vaNBIndex, int, j);
            iX = MIN(iI, iJ);
            iY = MAX(iI, iJ);
            iIndex = iY * (iY + 1) / 2 + iX + 1;        /* +1 because they are FORTRAN */
            /* style arrays !!!!! */

            /* Check if there is an H-Bond parameter for this */
            /* interaction, if there is, make iIndex = -iHBondIndex */
            /* -ve to signify that the index is into the HBOND tables */

            ParmSetAtom(uUnit->psParameters, i, sTypeI, &dMass, &dPolar,
                        &dE, &dR, &dE14, &dR14, &dScreenF, &iElement,
		       	&iHybridization, sDesc);
            ParmSetAtom(uUnit->psParameters, j, sTypeJ, &dMass, &dPolar,
                        &dE, &dR, &dE14, &dR14, &dScreenF, &iElement,
			&iHybridization, sDesc);
            iHBondIndex =
                iParmSetFindHBond(uUnit->psParameters, sTypeI, sTypeJ);
            if (iHBondIndex != PARM_NOT_FOUND)
                iIndex = -(iHBondIndex + 1);
            *PVAI(*vaPNBIndexMatrix, int, iI * iNonBonds + iJ) = iIndex;
        }
    }

    /* Now calculate the A,C parameters for all unique */
    /* NON-BOND interactions */

    *vaPNBParameters = vaVarArrayCreate(sizeof(NONBONDACt));
    VarArraySetSize((*vaPNBParameters), iNonBonds * (iNonBonds + 1) / 2);
    for (j = 0; j < iNonBonds; j++) {
        for (i = 0; i <= j; i++) {
            iX = i;
            iY = j;
            iIndex = iY * (iY + 1) / 2 + iX;
            dEI = PVAI(vaNonBonds, NONBONDt, i)->dE;
            dRI = PVAI(vaNonBonds, NONBONDt, i)->dR;
            dEJ = PVAI(vaNonBonds, NONBONDt, j)->dE;
            dRJ = PVAI(vaNonBonds, NONBONDt, j)->dR;
            MathOpConvertNonBondToAC(dEI, dRI, dEJ, dRJ, &dA, &dC);
	    CheckAgainstNBEdits(vaPNBEdits,
				PVAI(vaNonBonds, NONBONDt, i)->sType,
				PVAI(vaNonBonds, NONBONDt, j)->sType,
				&dA, &dC);
            PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dA = dA;
            PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dC = dC;
            if (GDefaults.iCharmm) {
                dEI = PVAI(vaNonBonds, NONBONDt, i)->dE14;
                dRI = PVAI(vaNonBonds, NONBONDt, i)->dR14;
                dEJ = PVAI(vaNonBonds, NONBONDt, j)->dE14;
                dRJ = PVAI(vaNonBonds, NONBONDt, j)->dR14;
                MathOpConvertNonBondToAC(dEI, dRI, dEJ, dRJ, &dA, &dC);
		CheckAgainstNBEdits(vaPNBEdits,
				    PVAI(vaNonBonds, NONBONDt, i)->sType,
				    PVAI(vaNonBonds, NONBONDt, j)->sType,
				    &dA, &dC);
                PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dA14 = dA;
                PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dC14 = dC;
            }
        }
    }

    /* Return the arrays that refer to the atoms */
    *vaPNonBonds = vaNonBonds;
    *vaPNBIndex = vaNBIndex;
    /* The other two arrays, vaPNBIndexMatrix and vaPNBParameters, */
    /* have already been assigned their return values. */

}

/*
 *  MarkSideChain() - depth-first descent of side chain
 */
static void MarkSideChain(ATOM aParentAtom, ATOM aAtom, int *iP)
{
    int i, j, k, iChildTag = *iP;
    ATOM aChildAtom, aNbrs[MAXBONDS];

    /*
     *  figure tree type - count 'downstream', 
     *      undesignated neighbors, i.e. omits
     *      parent and anything else already 
     *      encountered. Mark all claimed atoms
     *      as belonging to this one for the
     *      benefit of recursion looping back
     *      around to this point before getting
     *      to the atom later in this routine
     */
    j = 0;
    for (i = 0; i < iAtomCoordination(aAtom); i++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, i);
        if (aChildAtom == aParentAtom)
            continue;
        if (dAtomTemp(aChildAtom) == (double) 'x' &&
            iAtomTempInt(aChildAtom) == -1) {
            AtomSetTempInt(aChildAtom, iChildTag);
            j++;
        }
    }
    SetTreeType(aAtom, j);

    AtomSetSeenId(aAtom, *iP);
    (*iP)++;

    /*
     *  put eligible children in order for readability -
     *      get obvious 'E' types 1st
     */
    k = 0;
    for (j = 0; j < iAtomCoordination(aAtom); j++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, j);
        if (aChildAtom == aParentAtom)
            continue;
        if (dAtomTemp(aChildAtom) != (double) 'x')
            continue;
        if (iAtomCoordination(aChildAtom) == 1)
            aNbrs[k++] = aChildAtom;
    }
    /*
     *  bubble sort them into alphabetical order
     */
    while (1) {
        int iMore = 0;

        for (j = 0; j < k - 1; j++) {
            if (strcmp(sAtomName(aNbrs[j]), sAtomName(aNbrs[j + 1])) > 0) {
                ATOM aTmp;
                aTmp = aNbrs[j + 1];
                aNbrs[j + 1] = aNbrs[j];
                aNbrs[j] = aTmp;
                iMore++;
            }
        }
        if (!iMore)
            break;
    }
    /*
     *  print
     */
    for (j = 0; j < k; j++)
        MarkSideChain(aAtom, aNbrs[j], iP);

    /*
     *  get remaining eligible children
     */
    k = 0;
    for (j = 0; j < iAtomCoordination(aAtom); j++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, j);
        if (aChildAtom == aParentAtom)
            continue;
        if (iAtomCoordination(aChildAtom) == 1)        /* done above */
            continue;

        if (dAtomTemp(aChildAtom) == (double) 'x' &&
            iAtomTempInt(aChildAtom) == iChildTag) {
            aNbrs[k++] = aChildAtom;
        }
    }
    /*
     *  bubble sort
     */
    while (1) {
        int iMore = 0;

        for (j = 0; j < k - 1; j++) {
            if (strcmp(sAtomName(aNbrs[j]), sAtomName(aNbrs[j + 1])) > 0) {
                ATOM aTmp;
                aTmp = aNbrs[j + 1];
                aNbrs[j + 1] = aNbrs[j];
                aNbrs[j] = aTmp;
                iMore++;
            }
        }
        if (!iMore)
            break;
    }
    /*
     *  print
     */
    for (j = 0; j < k; j++)
        MarkSideChain(aAtom, aNbrs[j], iP);
}

static int MarkMainChainAtoms(RESIDUE rRes, int complain)
{
    int iAtomCount, i, j, iLevel, iMin, iMax, iNext, ierr = 0;
    ATOM aAtom, aAtom0, aAtom1, aChildAtom;
    VARARRAY vaAtoms;
    LOOP lTemp;
    char *cPResName;
    cPResName = sContainerName(rRes);

    /*
     *  set up for breadth-1st search: count atoms,
     *      setting up array of atom pointers & initializing depth
     */
    iAtomCount = 0;
    lTemp = lLoop((OBJEKT) rRes, ATOMS);
    while ((aAtom = (ATOM) oNext(&lTemp)) != NULL) {
        AtomSetTempInt(aAtom, -1);
        AtomSetTempDouble(aAtom, (double) 'x');
        iAtomCount++;
    }
    if (!iAtomCount) {
        VP0(("  %s: no atoms\n", cPResName));
        return (0);
    }
    if (iAtomCount == 1) {
        /*
         *  call it a main chain whether connected or not
         */
        lTemp = lLoop((OBJEKT) rRes, ATOMS);
        aAtom = (ATOM) oNext(&lTemp);
        AtomSetTempDouble(aAtom, (double) 'M');
        return (1);
    }

    /*
     *  check connect atoms
     */
    aAtom0 = (ATOM) rRes->aaConnect[0];
    aAtom1 = (ATOM) rRes->aaConnect[1];
    if (aAtom0 == NULL) {
        if (complain)
            VP0(("  %s:  connect0 not defined\n", cPResName));
        ierr++;
    }
    if (aAtom1 == NULL) {
        if (complain)
            VP0(("  %s:  connect1 not defined\n", cPResName));
        ierr++;
    }
    if (ierr)
        return (-iAtomCount);
/*
fprintf(stderr," connect %s .. %s total %d\n",
sAtomName( aAtom0 ), sAtomName( aAtom1), iAtomCount);
*/
    /*
     *  prepare 'stack' for tracking depth-1st search
     */
    vaAtoms = vaVarArrayCreate(sizeof(ATOM));
    VarArraySetSize(vaAtoms, iAtomCount + 1);
    *PVAI(vaAtoms, ATOM, 0) = aAtom0;
    iLevel = 0;
    iMin = 0;
    iMax = 1;

    /*
     *  loop over each level in tree starting at aAtom0
     *      until target atom found, marking atoms w/ tree level
     */
    while (1) {
        /*
         *  prepare to accumulate next level down
         */
        iNext = iMax;
        /*
         *  loop over atoms in  current level
         */
/*
fprintf(stderr,"--- lev %d:  %d .. %d\n", iLevel, iMin, iMax);
*/
        for (i = iMin; i < iMax; i++) {
            /*
             *  1st time, i=0 and 0th elt of vaAtoms is aAtom0
             */
            aAtom = *PVAI(vaAtoms, ATOM, i);
/*
fprintf(stderr,"     down %s lev %d coord %d\n",
sAtomName(aAtom), iLevel, iAtomCoordination(aAtom));
*/
            /* 
             *  mark level and see if this is target
             */
            AtomSetTempInt(aAtom, iLevel);
            if (aAtom == aAtom1)
                goto FOUND;

            /*
             *  put all unseen 'children' on
             *      next level stack
             */
            for (j = 0; j < iAtomCoordination(aAtom); j++) {
                aChildAtom = aAtomBondedNeighbor(aAtom, j);
/*
fprintf(stderr,"         child %s add %c\n",
sAtomName(aChildAtom),
(iAtomTempInt( aChildAtom ) == -1 ? 'Y' : 'N'));
*/
                /*
                 *  stay within residue
                 */
                if (cContainerWithin(aChildAtom) !=
                    cContainerWithin(aAtom)) continue;
                if (iAtomTempInt(aChildAtom) == -1) {
                    *PVAI(vaAtoms, ATOM, iNext++) = aChildAtom;
                    /*
                     *  re-flag atom to ensure
                     *      that it doesn't get included
                     *      again as a child of another
                     *      atom at this level
                     */
                    AtomSetTempInt(aChildAtom, -2);
                }
            }

        }
        iMin = iMax;
        iMax = iNext;
        iLevel++;
    }
  FOUND:
    /*
     *  Traverse marked tree from connect1 back up to
     *      connect0, marking main chain
     */
    aAtom = aAtom1;
    while (1) {
        ATOM aSelected = (ATOM) NULL;

        AtomSetTempDouble(aAtom, (double) 'M');
        if (aAtom == aAtom0)
            break;
        for (j = 0; j < iAtomCoordination(aAtom); j++) {
            aChildAtom = aAtomBondedNeighbor(aAtom, j);
            /*
             *  stay within residue
             */
            if (cContainerWithin(aChildAtom) != cContainerWithin(aAtom))
                continue;
            if (iAtomTempInt(aChildAtom) == iLevel - 1) {
                /*
                 *  choose lesser atom if >1, for readability
                 */
                if (aSelected == NULL) {
                    aSelected = aChildAtom;
                } else if (strcmp(sAtomName(aChildAtom),
                                  sAtomName(aSelected)) < 0) {
                    aSelected = aChildAtom;
                }
            }
        }
        aAtom = aSelected;
        iLevel--;
    }
    /*
     *  clean up
     */
    VarArrayDestroy(&vaAtoms);

    return (iAtomCount);
}

/*
 *  MarkSideChains() - used for prmtop
 *
 *        NOTE: MarkMainChainAtoms() must be run 1st
 */
static void MarkSideChains(RESIDUE rRes)
{
    ATOM aAtom, aParentAtom, aAtom0, aAtom1;
    LOOP lTemp;
    int iCount;

    aAtom0 = (ATOM) rRes->aaConnect[0];
    aAtom1 = (ATOM) rRes->aaConnect[1];
    if (aAtom0 == NULL)                /* 1-atom residue */
        return;

    /*
     *  traverse tree from top, following main chain
     *
     *      reset atom temp ints to use in marking 'seen' atoms
     */
    lTemp = lLoop((OBJEKT) rRes, ATOMS);
    while ((aAtom = (ATOM) oNext(&lTemp)) != NULL) {
        AtomSetTempInt(aAtom, -1);
        AtomSetSeenId(aAtom, -1);
    }

    aAtom = aAtom0;
    aParentAtom = NULL;
    iCount = 0;
    while (1) {
        ATOM aNextMain;
        int j;

        AtomSetSeenId(aAtom, iCount);
        /*
         *  find next main chain down in this residue, 
         *      marking side chains
         */
        aNextMain = NULL;
        for (j = 0; j < iAtomCoordination(aAtom); j++) {
            ATOM aChildAtom = aAtomBondedNeighbor(aAtom, j);
            if (aChildAtom == aParentAtom)
                continue;
            if (cContainerWithin(aChildAtom) != cContainerWithin(aAtom))
                continue;
            if (dAtomTemp(aChildAtom) == (double) 'M') {
                aNextMain = aChildAtom;
            } else if (dAtomTemp(aChildAtom) == (double) 'x') {
                /*
                 *  side chain not seen before
                 */
                MarkSideChain(aAtom, aChildAtom, &iCount);
            }
        }
        if (aAtom == aAtom1)
            break;
        if (aNextMain == NULL)        /* 1 'M' atom in residue */
            break;
        aParentAtom = aAtom;
        aAtom = aNextMain;
    }
}


/*
 *      zUnitIOSaveAmberParmFormat
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Save the UNIT in the AMBER PARM file format.
 *      This requires that the UNIT tables be built and that
 *      the UNIT contain a parameter set.
 *      The iContainerTempInt(atom) should still return the index
 *      of the atom within the vaAtoms array.
 *      Atom coordinates are written to the file (fCrd).
 *
 *NOTE:        This routine depends on the order of the RESIDUEs in
 *        vaResidues being such that solvent residues follow
 *        all other RESIDUEs.  I know that this is going
 *        against the philosophy that the data written to
 *        OFF files has NO implicit order, and outside of
 *        this program that is how they should be handled.
 *        But it was SO convenient to sort the RESIDUEs
 *        as they are put into the table that I could
 *        not resist.
 *
 *TODO: Add RESTRAINT code
 *TODO: Add CAP information
 */
#define AMBERINDEX(i)   3*(i-1)
#define INTFORMAT       "%8d"
#define DBLFORMAT       "%16.8lE"
#define LBLFORMAT       "%-4s"
#define IDFORMAT       "%-8s"
#define ELECTRONTOKCAL  18.2223

        /* RESTRAINTLOOP is used to loop over the RESTRAINTs */
        /* for adding constants to tables of constants */
#define        RESTRAINTLOOP( type, field, indexStart ) { \
int        ii, iiMax, jj = 0; \
    if ( (iiMax = iVarArrayElementCount( uUnit->vaRestraints )) ) { \
            srPRestraint = PVAI( uUnit->vaRestraints, SAVERESTRAINTt, 0 ); \
            for ( ii=0; ii<iiMax; ii++, srPRestraint++ ) { \
                if ( srPRestraint->iType == type ) { \
                            FortranWriteDouble( srPRestraint->field ); \
                            srPRestraint->iParmIndex = indexStart+jj; \
                            jj++; \
                } \
            } \
    } \
}
#define        bPERT_BOND(bp,a1,a2)        (bp && (bAtomFlagsSet(a1,ATOMPERTURB)\
                                || bAtomFlagsSet(a2,ATOMPERTURB)))
#define        bPERT_ANGLE(bp,a1,a2,a3) (bp && (bAtomFlagsSet(a1,ATOMPERTURB) \
                                || bAtomFlagsSet(a2,ATOMPERTURB)\
                                || bAtomFlagsSet(a3,ATOMPERTURB)))
#define        bPERT_TORSION(bp,a1,a2,a3,a4)        (bp && (bAtomFlagsSet(a1,ATOMPERTURB) \
                                || bAtomFlagsSet(a2,ATOMPERTURB)\
                                || bAtomFlagsSet(a3,ATOMPERTURB)\
                                || bAtomFlagsSet(a4,ATOMPERTURB)))


typedef struct {
SAVETORSIONt *tp;
} SAVETORSIONtp;

static void
get4atoms(UNIT u, SAVETORSIONt *pt, SAVEATOMt *sa4[])
{
            sa4[0] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom1 - 1);
            sa4[1] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom2 - 1);
            sa4[2] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom3 - 1);
            sa4[3] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom4 - 1);
}
                                  
static int
copyatoms(int atoms[], SAVETORSIONt *sa4, SAVETORSIONt *sb4, int k1, int k2)
{
int atma[4], atmb[4], i;

  if ( k1 > 0 ) {
    atma[0] = AMBERINDEX( sa4->iAtom1 );
    atma[1] = AMBERINDEX( sa4->iAtom2 );
    atma[2] = AMBERINDEX( sa4->iAtom3 );
    atma[3] = AMBERINDEX( sa4->iAtom4 );
    } else {
    atma[3] = AMBERINDEX( sa4->iAtom1 );
    atma[2] = AMBERINDEX( sa4->iAtom2 );
    atma[1] = AMBERINDEX( sa4->iAtom3 );
    atma[0] = AMBERINDEX( sa4->iAtom4 );
    }

  if ( k2 > 0 ) {
    atmb[0] = AMBERINDEX( sb4->iAtom1 );
    atmb[1] = AMBERINDEX( sb4->iAtom2 );
    atmb[2] = AMBERINDEX( sb4->iAtom3 );
    atmb[3] = AMBERINDEX( sb4->iAtom4 );
    } else {
    atmb[3] = AMBERINDEX( sb4->iAtom1 );
    atmb[2] = AMBERINDEX( sb4->iAtom2 );
    atmb[1] = AMBERINDEX( sb4->iAtom3 );
    atmb[0] = AMBERINDEX( sb4->iAtom4 );
    }
// check overlap
   for (i=0; i < 3; i++) if (atmb[i] != atma[i+1])
    {
      VP0(("Atom indices mismatch?? in copyatoms %i %i\n",atma[i+1],atmb[i]));
      return -1;
    } else {
      atoms[i] = atma[i];
      atoms[i+2] = atmb[i+1];
    }
return 4;

}

static int
cmpresname1(UNIT u, SAVEATOMt *sa4, WRD reslist[], int nres)
{
int i,j;
int l;
char *sname;

   sname = PVAI(u->vaResidues,SAVERESIDUEt,sa4->iResidueIndex - 1)->sName;
   for (i = 0; i < nres; i++) {
      if ( strcmp(sname, reslist[i]) == 0 ) return (i+1);
   }

return -1;

}

static int
cmpresname4(UNIT u, SAVEATOMt *sa4[], WRD reslist[], int nres)
{
int i,j;
int l;
char *sname;

  for (j=0; j < 4; j++) {
      if ((i=cmpresname1(u, sa4[j], reslist, nres)) > 0 ) return (i+1); }

return -1;

}

static int
cmp4vs4(SAVEATOMt *sa4[], WRD atm4[])
{
int i, l1, l2;
l1=0;
l2=0;
for (i = 0; i< 4; i++) if (strcmp(sa4[i]->sName,   atm4[i]) == 0)   l1++ ; 
if (l1 == 4 ) return 4;
for (i = 0; i< 4; i++) if (strcmp(sa4[3-i]->sName, atm4[i]) == 0)   l2++ ; 
if (l2 == 4 ) return -4;

return 0;

}

static int
cmp_residx(SAVEATOMt *sa4[], SAVEATOMt *sb4[], int *residx)
{
int i, l1, l2;
int l, idx0, idx1;

idx0 = sa4[0]->iResidueIndex;
idx1 = residx[0];
l1=0;
l2=0;
for (i = 0; i< 4; i++) {
     if ((sa4[i]->iResidueIndex - idx0) == (residx[i] - idx1 )) 
         l1++ ;  }
for (i = 0; i< 4; i++) {
     if ((sb4[i]->iResidueIndex - idx0) == (residx[i+1] - idx1 )) 
         l2++ ;  }

if (l1 == 4 && l2 == 4 ) return 1;

return 0;

}

static void
SaveAmberParmCMAP(UNIT uUnit, FILE * fOut)
{
    // 
    // CMAP parameters, Mengjuei Hsieh and Yong Duan
    //
        int i,j,k,l;
        int iterCMAP, iterRes, mapid, mapcount, maptypes;
        int *mapflag, *mapidx;
        int iNumDIH, iNumRes;
        int *prospect, nprospect, ires;
        ATOM aE, aF, aG, aH;
        SAVEATOMt *saA, *saB, *saC, *saD, *saE, *saF, *saG, *saH;
        SAVEATOMt *sa4[4], *sb4[4];
        SAVETORSIONt *stPTorsion2, *stPTorsion, *stPTorsiont, *stPTorsion2t;
        SAVETORSIONtp *stdptt, *stdpt0;
        CMNT *cmntt;
        STRING sTmp;
        PHIPSI *phipsi;
        CMAPLST *cmaplstt;
        int k1, k2;
        CMAP *cmap;
        int maxmap;
        
    if (!GDefaults.iCMAP) return;
    if (mapnum <= 0) return;

        mapcount = 0;
        mapflag = (int *) malloc(sizeof(int) * (mapnum+1) );
        mapidx  = (int *) malloc(sizeof(int) * (mapnum+1) );

        i = 0;
        for (i=0; i <= mapnum; i++) {
          mapflag[i] = 0;
          mapidx [i] = 0;
        }
        
    for (cmaplstt = cmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next)
      {
           cmap=cmaplstt->cmap;
      }
        iNumDIH=iVarArrayElementCount(uUnit->vaTorsions);
        if (iNumDIH > 0){
            stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        } else {
            return;
        }
        
        iNumRes=iVarArrayElementCount(uUnit->vaResidues);

        maxmap = iNumDIH;
        if (maxmap > 0){
            //srPRes = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
            phipsi = (PHIPSI *) malloc(sizeof(PHIPSI) * (maxmap));
        } else {
            exit(1);
        }

        stdpt0 = (SAVETORSIONtp *) malloc(sizeof(SAVETORSIONtp) * (iNumDIH));
        nprospect = 0;
        mapcount  = 0;
        // Loop over dihedral list
// pre-filter removes the irrelevant torsions first ...
        for (i = 0; i < iNumDIH; i++, stPTorsion++) {

//            if (stPTorsion->bCalc14 == 0) continue; //cycle

            get4atoms(uUnit, stPTorsion, sa4);

            for (cmaplstt = cmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next)
               {
                  cmap=cmaplstt->cmap;
                  if (cmpresname4(uUnit, sa4, cmap->reslist, cmap->nres) > 0)  // potential match
                  {
                  if ( abs(cmp4vs4(sa4, (WRD *) (&cmap->atmname[0]))) == 4 || 
                       abs(cmp4vs4(sa4, (WRD *) (&cmap->atmname[1]))) == 4 )
                  {   // keep in the list
                  stdpt0[nprospect].tp = stPTorsion;
                  nprospect ++;
                  break;
                  }
                 } else if ( cmap->termmap > 0 )  {
                  if (cmpresname4(uUnit, sa4, cmap->creslist, cmap->nres) > 0)  {
                  if ( abs(cmp4vs4(sa4, (WRD *) (&cmap->catmname[0]))) == 4 || 
                       abs(cmp4vs4(sa4, (WRD *) (&cmap->catmname[1]))) == 4 ) 
                  {
                  stdpt0[nprospect].tp = stPTorsion;
                  nprospect ++;
                  break;
                  } }
                  else if (cmpresname4(uUnit, sa4, cmap->nreslist, cmap->nres) > 0) 
                  if ( abs(cmp4vs4(sa4, (WRD *) (&cmap->natmname[0]))) == 4 || 
                       abs(cmp4vs4(sa4, (WRD *) (&cmap->natmname[1]))) == 4 ) 
                  {
                  stdpt0[nprospect].tp = stPTorsion;
                  nprospect ++;
                  break;
                  }
                 }
               }
            }
//
        k = 0;
        int k1c,k2c,k1n,k2n;
        for (i = 0; i < nprospect; i++) {
            stPTorsion = stdpt0[i].tp;
            get4atoms(uUnit, stPTorsion, sa4);
            mapid = 0;
            for (cmaplstt = cmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next)
              {  cmap=cmaplstt->cmap;
                 mapid ++;
                 k1=cmp4vs4(sa4, (WRD *) (&cmap->atmname[0])); // check atmnames 0..3
                 k1c=cmp4vs4(sa4, (WRD *) (&cmap->catmname[0])); // check atmnames 0..3
                 k1n=cmp4vs4(sa4, (WRD *) (&cmap->natmname[0])); // check atmnames 0..3

                 if (cmap->termmap == 0) { k1c=0; k1n=0; }

                 int iresc, iresn;

                 if (abs(k1) == 4 || abs(k1c) == 4 || abs(k1n) == 4 ) { // match the atom names 
                   // "0" in residx[] marks "present residue"
                 for ( l = 0; l < 5 && cmap->residx[l] != 0; l++);
                 if ( l < 4 ) { ires = cmpresname1(uUnit, sa4[l], cmap->reslist, cmap->nres);
                                iresc= cmpresname1(uUnit, sa4[l], cmap->creslist, cmap->nres);
                                iresn= cmpresname1(uUnit, sa4[l], cmap->nreslist, cmap->nres); 

                 if (abs(k1) != 4 ) ires =0;
                 if (abs(k1c) != 4 ) iresc =0;
                 if (abs(k1n) != 4 ) iresn =0; }

                 if ( ires > 0 || iresc > 0 || iresn > 0 ||  // found on the reslist of the cmap
                        l == 4 ) { // or the "present residue" pointer is the last one.
                 for (j = 0; j < nprospect; j++) {
                     stPTorsion2 = stdpt0[j].tp;
                     get4atoms(uUnit, stPTorsion2, sb4);
                     if ( l == 4 ) 
                        { ires = cmpresname1(uUnit, sb4[3], cmap->reslist, cmap->nres);
                          iresc= cmpresname1(uUnit, sb4[3], cmap->creslist, cmap->nres);
                          iresn= cmpresname1(uUnit, sb4[3], cmap->nreslist, cmap->nres); }

                     if (cmap->termmap == 0) {iresc = 0; iresn = 0; }

                     k2 = cmp4vs4(sb4, (WRD *) (&cmap->atmname[1])); // check atmnames 1..4
                     k2c= cmp4vs4(sb4, (WRD *) (&cmap->catmname[1])); // check atmnames 1..4
                     k2n= cmp4vs4(sb4, (WRD *) (&cmap->natmname[1])); // check atmnames 1..4

                     if (abs(k2) == 4 && abs(k1) == 4 && ires > 0 ) 
                      // check residue indicies, 0: present, -1: residue before, +1: residue after
                           if ( cmp_residx(sa4, sb4, cmap->residx) )  // we found two torsions and the cmap
                              {   copyatoms(phipsi[k].atoms, stPTorsion, stPTorsion2, k1, k2);
                                  phipsi[k].mapid = mapid;
                                  // mark as used
                                  mapflag[mapid-1] = 1;
                                  k ++; }

                     if (abs(k2c) == 4 && abs(k1c) == 4 && iresc > 0 ) 
                           if ( cmp_residx(sa4, sb4, cmap->cresidx) )
                              {   copyatoms(phipsi[k].atoms, stPTorsion, stPTorsion2, k1c, k2c);
                                  phipsi[k].mapid = mapid;
                                  // mark as used
                                  mapflag[mapid-1] = 1;
                                  k ++; }

                     if (abs(k2n) == 4 && abs(k1n) == 4 && iresn > 0 ) 
                           if ( cmp_residx(sa4, sb4, cmap->nresidx) )
                              {   copyatoms(phipsi[k].atoms, stPTorsion, stPTorsion2, k1n, k2n);
                                  phipsi[k].mapid = mapid;
                                  // mark as used
                                  mapflag[mapid-1] = 1;
                                  k ++; }
                           }
                       }
                   }
                }
            }

        int mk=0;
        for (i=0; i<k; i++) {
            int flag = 1;

// Remove the redundant torsions.

            int l;
            if (i>0) for (l=0; l<i; l++) {
            if (
                phipsi[i].atoms[0] == phipsi[l].atoms[0] &&
                phipsi[i].atoms[1] == phipsi[l].atoms[1] &&
                phipsi[i].atoms[2] == phipsi[l].atoms[2] &&
                phipsi[i].atoms[3] == phipsi[l].atoms[3] &&
                phipsi[i].atoms[4] == phipsi[l].atoms[4] &&
                phipsi[i].mapid == phipsi[l].mapid )
                flag=-1;
        }
//
            for (j=0; j<5; j++) if (phipsi[i].atoms[j] < 0) flag=-1;
            if (phipsi[i].mapid <= 0) flag=-1;
            if (flag > 0) mk++;
    }

        //wmap
        maptypes = 0;
  
        for (i=0; i < mapnum; i++) {
           if (mapflag[i] == 1) {
              mapidx[i] = maptypes;
              maptypes ++;
            }
        }
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG CMAP_COUNT");
        //FortranWriteString("%COMMENT");
//        for (cmntt=cmnt0; cmntt->next != NULL; cmntt=cmntt->next){
//            sTmp[0]='\0';
//            strcat(sTmp,"%COMMENT");
//            strcat(sTmp,cmntt->record);
//            FortranWriteString(sTmp);
            //fprintf(fpout,"%%COMMENT%s",cmntt->record);
//        }
        FortranWriteString("%FORMAT(2I8)");
//        sprintf(sTmp,"%8d%8d",mapcount,mapnum);
//        sprintf(sTmp,"%8d%8d",k,maptypes);
        sprintf(sTmp,"%8d%8d",mk,maptypes);
        FortranWriteString(sTmp);
        //FortranEndLine();
        
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG CMAP_RESOLUTION");
        FortranWriteString("%FORMAT(20I4)");
        FortranFormat(20, "%4d");
        mapid = 0;
        for (cmaplstt = cmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next)
           {  cmap=cmaplstt->cmap;
              if (mapflag[mapid]) FortranWriteInt(cmap->resolution) ;
              mapid ++;
           }
        FortranEndLine();
        
        mapid = 0;
        for (cmaplstt = cmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next)
           {  cmap=cmaplstt->cmap;
            STRING str;
            int msize;
            if (mapflag[mapid] == 1) {
            sprintf(str,"0%1i",mapidx[mapid]+1);
            if (mapidx[mapid]+1>9) sprintf(str,"%1i",mapidx[mapid]+1);
            sprintf(sTmp,"%%FLAG CMAP_PARAMETER_%2s",str);
            FortranFormat(1, "%-80s");
            sTmp[79]='\0';
            FortranWriteString(sTmp);
            sprintf(sTmp,"%%COMMENT  %s",cmap->title);
            FortranWriteString(sTmp);
            FortranWriteString("%FORMAT(8F9.5)");
            FortranFormat(8, "%9.5lf");
            msize = cmap->resolution*cmap->resolution;
            for (j=0; j<msize; j+=8){
                int l;
                for (l=j; l<msize && l<j+8; l++){
                    FortranWriteDouble(cmap->map[l]);
                }
            }
            FortranEndLine();
            }
              mapid ++;
        }
        
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG CMAP_INDEX");
        FortranWriteString("%FORMAT(6I8)");
        FortranFormat(6, "%8d");
        for (i=0; i<k; i++) {
            int flag = 1;

// Remove the redundant torsions.

            int l;
            if (i>0) for (l=0; l<i; l++) {
            if (
                phipsi[i].atoms[0] == phipsi[l].atoms[0] && 
                phipsi[i].atoms[1] == phipsi[l].atoms[1] && 
                phipsi[i].atoms[2] == phipsi[l].atoms[2] && 
                phipsi[i].atoms[3] == phipsi[l].atoms[3] && 
                phipsi[i].atoms[4] == phipsi[l].atoms[4] && 
                phipsi[i].mapid == phipsi[l].mapid ) 
                flag=-1;
               }
//
            for (j=0; j<5; j++) if (phipsi[i].atoms[j] < 0) flag=-1;
            if (phipsi[i].mapid <= 0) flag=-1;
            if (flag > 0) {
                for (j=0; j<5; j++) FortranWriteInt(phipsi[i].atoms[j]/3+1);
                FortranWriteInt(mapidx[phipsi[i].mapid-1]+1);
            }
        }
        FortranEndLine();
    }

/*
 *      zUnitIOSaveAmberNetcdf
 * 
 *      Author: Robin Betz (2013)
 * 
 *      Writes the coordinates in UNIT to a coordinate file in netCdf format.
 *      This is written based from NetcdfFile.cpp in cpptraj and AmberNetcdf.F90
 *      in pmemd/sander
 * 
 *      Arguments:
 *          uUnit     - UNIT to save
 *          filename  - name of netcdf file to write
 */
void
zUnitIOSaveAmberNetcdf( UNIT uUnit, char *filename )
{
#   ifndef BINTRAJ
    VP0(( "Error: Compiled without NETCDF support. Recompile with -DBINTRAJ\n"));
#   else

  int ncid;                             // netcdf file handle
  int did_spatial, did_atom, did_frame; // dimension IDs
  int vid_spatial, vid_coord; // variable IDs
  int did_cell_spatial, did_cell_angular, did_label;
  int vid_cell_spatial, vid_cell_angular, vid_cell_length, vid_cell_angle, vid_time;
  int dimensionID[NC_MAX_VAR_DIMS];
  
  size_t start[2], count[2];
  
  // Get number of atoms T_T
  int iAtomCount = iVarArrayElementCount(uUnit->vaAtoms);
  printf("There are %i atoms\n", iAtomCount);
  
  // Create the file
  if (nc_create( filename, NC_64BIT_OFFSET, &ncid ) != NC_NOERR ) {
    VP0(( "%s: Error creating file\n", filename ));
  }
  // Spatial dimension and variable
  if ( nc_def_dim( ncid, "spatial", 3, &did_spatial) != NC_NOERR ) {
    VP0(( "%s: Error defining spatial dimension\n", filename ));
  }
  dimensionID[0] = did_spatial;
  if ( nc_def_var( ncid, "spatial", NC_CHAR, 1, dimensionID, &vid_spatial ) != NC_NOERR ) {
    VP0(( "%s: Error defining spatial variable\n", filename ));
  }
  // Atom dimension
  if ( nc_def_dim( ncid, "atom", iAtomCount, &did_atom)
    != NC_NOERR ) {
    VP0(( "%s: Error defining atom dimension\n", filename ));
    }
    // Time dimension and variable
    if ( nc_def_var(ncid, "time", NC_DOUBLE, 0, dimensionID, &vid_time) != NC_NOERR ) {
      VP0(( "%s: Error defining time variable\n", filename ));
    }
    if ( nc_put_att_text(ncid, vid_time, "units", 10, "picosecond") != NC_NOERR ) {
      VP0(( "%s: Error setting time units to picosecond\n", filename ));
    }
    dimensionID[0] = did_atom;
    dimensionID[1] = did_spatial;
    // Coord variable and attribute text
    if ( nc_def_var( ncid, "coordinates", NC_DOUBLE, 2, dimensionID, &vid_coord ) != NC_NOERR ) {
      VP0(( "%s: Error defining coordinate variable\n", filename ));
    }
    if ( nc_put_att_text( ncid, vid_coord, "units", 8, "angstrom") != NC_NOERR ) {
      VP0(( "%s: Error setting coordinate units to angstrom\n", filename ));
    }
    
    // Define box if it exists
    if ( bUnitUseBox(uUnit) == TRUE ) {
      printf("Using the unit box\n");
      // Cell spatial
      if ( nc_def_dim( ncid, "cell_spatial", 3, &did_cell_spatial) != NC_NOERR ){
        VP0(( "%s: Error defining cell spatial dimension\n", filename ));
      }
      dimensionID[0] = did_cell_spatial;
      if (nc_def_var(ncid, "cell_spatial", NC_CHAR, 1, dimensionID, &vid_cell_spatial) 
        != NC_NOERR) {
        VP0(( "%s: Error defining cell spatial variable\n", filename ));
        }
        // Cell angular
        if ( nc_def_dim(ncid, "label", 5, &did_label) != NC_NOERR) {
          VP0(( "%s: Error defining label dimension\n", filename ));
        }
        if ( nc_def_dim( ncid, "cell_angular", 3, &did_cell_angular) != NC_NOERR) {
          VP0(( "%s: Error defining cell angular dimension\n", filename ));
        }
        dimensionID[0] = did_cell_angular;
        dimensionID[1] = did_label;
        if ( nc_def_var( ncid, "cell_angular", NC_CHAR, 2, dimensionID, &vid_cell_angular) 
          != NC_NOERR ) {
            VP0(( "%s: Error defining cell angular variable\n", filename ));
          }
          // Box dimensions
          dimensionID[0] = did_cell_spatial;
        if ( nc_def_var( ncid, "cell_lengths", NC_DOUBLE, 1, dimensionID, &vid_cell_length )
          != NC_NOERR ) {
            VP0(( "%s: Error defining cell lengths\n", filename ));
          }
          if ( nc_put_att_text(ncid, vid_cell_length, "units", 8, "angstrom") != NC_NOERR) {
            VP0(( "%s: Error setting cell length units to angstrom\n", filename ));
          }
          dimensionID[0] = did_cell_angular;
        if (nc_def_var(ncid, "cell_angles",NC_DOUBLE,1,dimensionID, &vid_cell_angle) 
          != NC_NOERR) {
            VP0(( "%s: Error defining cell angles variable\n", filename ));
          }
          if (nc_put_att_text(ncid,vid_cell_angle,"units",6,"degree") != NC_NOERR) {
            VP0(( "%s: Error setting cell angle units to degree\n", filename ));
          }
    }
    
    // Conventions and file attributes
    if (nc_put_att_text( ncid,NC_GLOBAL,"title", strlen(sContainerName(uUnit)),
      sContainerName(uUnit) ) != NC_NOERR ) {
      VP0(( "%s: Error writing title\n", filename ));
      }
      if (nc_put_att_text(ncid,NC_GLOBAL,"application",5,"AMBER") != NC_NOERR ) {
        VP0(( "%s: Error writing application string\n", filename ));
      }
      if (nc_put_att_text(ncid,NC_GLOBAL,"program",4,"leap") != NC_NOERR ){
        VP0(( "%s: Error writing program string\n", filename ));
      }
      if (nc_put_att_text(ncid,NC_GLOBAL,"programVersion",3,"1.0") != NC_NOERR){
        VP0(( "%s: Error writing program version string\n", filename ));
      }
      if (nc_put_att_text(ncid,NC_GLOBAL,"Conventions",12,"AMBERRESTART") != NC_NOERR) {
        VP0(( "%s: Error writing conventions\n", filename ));
      }
      if (nc_put_att_text(ncid,NC_GLOBAL,"ConventionVersion",3,"1.0") != NC_NOERR) {
        VP0(( "%s: Error writing conventions version\n", filename ));
      }
      
      // Set fill mode and end definitions
      if (nc_set_fill(ncid, NC_NOFILL, dimensionID) != NC_NOERR) {
        VP0(( "%s: Error setting fill mode\n", filename ));
      }
      if (nc_enddef(ncid) != NC_NOERR) {
        VP0(( "%s: NetCDF error on ending definitions\n", filename ));
      }
      
      // Spatial dimension labels
      start[0] = 0; count[0] = 3;
      char xyz[3];
      xyz[0]='x'; xyz[1]='y'; xyz[2]='z';
      if (nc_put_vara_text(ncid, vid_spatial, start, count, xyz) != NC_NOERR) {
        VP0(( "%s: Error writing spatial labels\n", filename ));
      }
      if (bUnitUseBox(uUnit) == TRUE) {
        xyz[0]='a';xyz[1]='b',xyz[2]='c';
        if (nc_put_vara_text(ncid,vid_cell_spatial,start,count,xyz) != NC_NOERR) {
          VP0(( "%s: Error writing cell spatial labels\n", filename ));
        }
        char abc[15] = { 'a', 'l', 'p', 'h', 'a',
        'b', 'e', 't', 'a', ' ',
        'g', 'a', 'm', 'm', 'a' };
        start[0]=0; start[1]=0; count[0]=3; count[1]=5;
        if (nc_put_vara_text(ncid,vid_cell_angular, start,count,abc) != NC_NOERR) {
          VP0(( "%s: Error writing cell angular labels\n", filename ));
        }
      }

// Create data array
double *data = (double*)malloc(3*iAtomCount*sizeof(double));
int counter = 0;
VECTOR vPos; ATOM aAtom;
LOOP lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );

// Write time = 0
double time = 0.0;
if (nc_put_var_double(ncid, vid_time, &time) != NC_NOERR) {
  VP0(( "%s: Error writing start time\n", filename ));
}

// Calculate box shift if there's a box and we're centering
double dX, dY, dZ;
double dX2, dY2, dZ2;
if ( bUnitUseBox(uUnit) == TRUE
  && GDefaults.nocenter == 0 ) {
  UnitGetBox(uUnit, &dX, &dY, &dZ);
dX2 = dX * 0.5;
dY2 = dY * 0.5;
dZ2 = dZ * 0.5;
  } else {
    dX2 = dY2 = dZ2 = 0.0;
  }
  
  // Fill the coordinate array with calculated shift (none if nobox)
  int i;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    vPos = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->vPos;
    data[counter] = dVX(&vPos) + dX2; ++counter;
    data[counter] = dVY(&vPos) + dY2; ++counter;
    data[counter] = dVZ(&vPos) + dZ2; ++counter;
  }
  
  // Write the coordinate array to the file
  start[0] = start[1] =  0;
  count[0] = iAtomCount; count[1] = 3;
  if (nc_put_vara_double(ncid, vid_coord, start, count, data) != NC_NOERR) {
    free(data);
    VP0(( "%s: Error writing coordinate data\n", filename ));
  }
  
  // Write box lengths and angles if necessary
  if ( bUnitUseBox(uUnit) == TRUE ) {
    count[0] = 3; count[1] = 0;
    double lengths[3];
    lengths[0] = dX;
    lengths[1] = dY;
    lengths[2] = dZ;
    if ( nc_put_vara_double(ncid, vid_cell_length, start, count, lengths) != NC_NOERR) {
      free(data);
      VP0(( "%s: Error writing cell lengths\n", filename ));
    }
    lengths[0] = lengths[1] = lengths[2] = dUnitBeta(uUnit) / DEGTORAD;
    if ( nc_put_vara_double(ncid, vid_cell_angle, start, count, lengths) != NC_NOERR) {
      free(data);
      VP0(( "%s: Error writing cell angles\n", filename ));
    }
  }
  
  // Close the file
  if (nc_close(ncid) != NC_NOERR) {
    free(data);
    VP0(( "%s: Error closing file\n", filename ));
  }
  free(data);
  printf("Successfully saved NetCDF inpcrd file \"%s\"\n", filename);
#   endif
}

void
zUnitIOSaveAmberParmFormat(UNIT uUnit, FILE * fOut, char* crdName,
                           BOOL bPolar, BOOL bPert, BOOL bNetcdf)
{
    int i, iMax, iIndex;
    LOOP lTemp, lSpan;
    ATOM aAtom, aAtomA, aA, aB, aC, aD;
    int iCount, iBondWith, iBondWithout;
    int iAngleWith, iAngleWithout;
    int iTorsionWith, iTorsionWithout;
    int iNumExtra;
    int iResidueIndex;
    VARARRAY vaExcludedAtoms, vaExcludedCount;
    VARARRAY vaNBIndexMatrix, vaNBParameters;
    VARARRAY vaNBIndex, vaNonBonds;
    int iCountPerturbed, iCountBondPerturbed, iCountBondBoundary;
    int iCountAnglePerturbed, iCountAngleBoundary;
    int iCountTorsionPerturbed, iCountTorsionBoundary;
    SAVEBONDt *sbPBond;
    SAVEANGLEt *saPAngle;
    SAVEATOMt *saPAtom;
    SAVETORSIONt *stPTorsion;
    SAVERESTRAINTt *srPRestraint;
    double dMass, dPolar, dR, dKb, dR0, dKt, dT0, dTkub, dRkub, dKp, dP0,
        dC, dD, dTemp;
    double dScEE, dScNB;
    double dScreenF, dSceeScaleFactor;
    STRING sAtom1, sAtom2, sAtom3, sAtom4, sType1, sType2;
    int iN, iAtoms, iMaxAtoms, iTemp, iAtom, iCalc14, iProper;
    int iElement, iHybridization, iStart, iFirstSolvent;
    RESIDUE rRes;
    BOOL bFoundSome;
    VECTOR vPos;
    char *cPTemp;
    double dX, dY, dZ, dEpsilon, dRStar, dEpsilon14, dRStar14;
    STRING sDesc, sType;
    VARARRAY vaMolecules;
    IX_REC *ePResEnt;
    IX_DESC iResIx;
    char sVersionHeader[81];
    time_t tp;
    double dGBrad, dGBscreen;

    // Open the coordinate file
    FILE *fCrd = FOPENCOMPLAIN( crdName, "w" );
    if ( fCrd == NULL ) {
      VP0(( "%s: Could not open file: %s\n", crdName ));
    }
    
// sanity check

if (bPolar && GDefaults.iIPOL <= 0) // NOT allowed to save IPOL=0 for polarizable parm
  {
    VP0(("  Conflict: polarizable prmtop can not have IPOL <= 0.\n"));
    VP0(("  Please change IPOL in frcmod/parmxx.dat or set default IPOL.\n"));
    return;
  } else if ( ! bPolar && GDefaults.iIPOL > 0 ) {
    VP0(("  Conflict: non-polarizable prmtop can not have IPOL > 0.\n"));
    VP0(("  Please change IPOL in frcmod/parmxx.dat or set default IPOL.\n"));
    return;
  }

    /* Build the excluded atom list */


    MESSAGE(("Building the excluded atom list\n"));
    vaExcludedCount = vaVarArrayCreate(sizeof(int));
    vaExcludedAtoms = vaVarArrayCreate(sizeof(int));

    iCountPerturbed = 0;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;

        if (bAtomFlagsSet(aAtom, ATOMPERTURB))
            iCountPerturbed++;

        lSpan = lLoop((OBJEKT) aAtom, SPANNINGTREE);
        iCount = 0;
        bFoundSome = FALSE;
        iStart = iVarArrayElementCount(vaExcludedAtoms);
        while ((aA = (ATOM) oNext(&lSpan))) {

            if (aA == aAtom)
                continue;

            /* If the atom is more than three away from the first atom */
            /* then it is not in the excluded atom list */

            if (iAtomBackCount(aA) >= 4)
                break;

            if (iContainerTempInt(aA) > iContainerTempInt(aAtom)) {
                VarArrayAdd(vaExcludedAtoms,
                            (GENP) & iContainerTempInt(aA));
                bFoundSome = TRUE;
                iCount++;
            }
        }
        if (!bFoundSome) {
            iAtoms = 0;
            VarArrayAdd(vaExcludedAtoms, (GENP) & iAtoms);
            iCount++;
        } else {

            /* Sort the part of the VARARRAY just added so that */
            /* the excluded ATOMs are in ascending order by index */

            SortByInteger((GENP) PVAI(vaExcludedAtoms, int, iStart),
                          iCount,
                          sizeof(int),
                          (GENP) PVAI(vaExcludedAtoms, int, iStart), TRUE);
        }

        VarArrayAdd(vaExcludedCount, (GENP) & iCount);
    }

    /*
     *  mark main chain atoms where possible, noting the 
     *  number of atoms in the largest residue. keep
     *  track of residues which can't be marked.
     */
    VP0(("Not Marking per-residue atom chain types.\n"));
    iMaxAtoms = 0;

    create_index(&iResIx, 2, 0);
    MALLOC(ePResEnt, IX_REC *, sizeof(IX_REC) + 8);
    ePResEnt->recptr = NULL;        /* for Purify */

    VP0(("Marking per-residue atom chain types.\n"));
    iMaxAtoms = 0;
    lTemp = lLoop((OBJEKT) uUnit, RESIDUES);
    while ((rRes = (RESIDUE) oNext(&lTemp))) {
        int iAtoms = MarkMainChainAtoms(rRes, 0);
        if (iAtoms > 0)
            (void) MarkSideChains(rRes);
        if (iAtoms < 0) {
            iAtoms = -iAtoms;
            /*
             *  couldn't mark main chains
             */

            strcpy(ePResEnt->key, rRes->cHeader.sName);
            if (add_key(ePResEnt, &iResIx) != IX_OK)
                DFATAL(("add_key() residue chain\n"));
        }
        if (iAtoms > iMaxAtoms)
            iMaxAtoms = iAtoms;
    }
    /*
     *  print warnings
     */
    first_key(&iResIx);
    i = 1;
    while (next_key(ePResEnt, &iResIx) == IX_OK) {
        if (i) {
            VP0(("  (Residues lacking connect0/connect1 - \n"));
            VP0(("   these don't have chain types marked:\n\n"));
            VP0(("\tres\ttotal affected\n\n"));
            i = 0;
        }
        VP0(("\t%s\t%d\n", ePResEnt->key, ePResEnt->count));
    }
    if (!i)
        VP0(("  )\n"));
    destroy_index(&iResIx);
    FREE(ePResEnt);

    /* Build the NON-BOND arrays that AMBER needs */

    zUnitIOBuildNonBondArrays(uUnit, &vaNBIndexMatrix, &vaNBParameters,
                              &vaNBIndex, &vaNonBonds);
    FortranFile(fOut);

#if 0
    /*
     *---------------------------------------------------------
     *
     *      Turn on debugging of fortran format output file
     *      by sticking comments into the file.
     */

    FortranDebugOn();
#endif


    /* -1- Save the title of the UNIT */
    FortranDebug("-1-");
    MESSAGE(("Saving the name of the UNIT\n"));
    FortranFormat(1, "%-80s");
    time(&tp);
    strftime(sVersionHeader, 81,
             "%%VERSION  VERSION_STAMP = V0001.000  DATE = %m/%d/%y  %H:%M:%S\0",
             localtime(&tp));
    FortranWriteString(sVersionHeader);
    FortranWriteString("%FLAG TITLE");
    FortranWriteString("%FORMAT(20a4)");
    FortranWriteString(sContainerName(uUnit));

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG POINTERS");
    FortranWriteString("%FORMAT(10I8)");

    /* -2- Save control information */
    FortranDebug("-2-");
    MESSAGE(("Saving all the main control variables\n"));
    FortranFormat(10, INTFORMAT);

     /*NTOTAT*/ FortranWriteInt(iVarArrayElementCount(uUnit->vaAtoms));
     /*NTYPES*/ FortranWriteInt(iVarArrayElementCount(vaNonBonds));

    /* Count the number of bonds with hydrogens, and without */

    iBondWith = 0;
    iBondWithout = 0;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++) {
        sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
        aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1 - 1)->aAtom;
        aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2 - 1)->aAtom;
        if (bPERT_BOND(bPert, aA, aD))
            continue;
        if (iAtomElement(aA) == HYDROGEN || iAtomElement(aD) == HYDROGEN)
            iBondWith++;
/*
 *  dac: just a syntax check for potential change in functionality:
 *      if (sAtomName(aA)[0] == 'H' || sAtomName(aD)[0] == 'H')
 *          iBondWith++;
 */
        else
            iBondWithout++;
    }
     /*NBONH*/ FortranWriteInt(iBondWith);
     /*NBONA*/ FortranWriteInt(iBondWithout);

    /* Count the number of angles with hydrogens, and without */

    iAngleWith = 0;
    iAngleWithout = 0;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++) {
        saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
        aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1 - 1)->aAtom;
        aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2 - 1)->aAtom;
        aD = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3 - 1)->aAtom;
        if (bPERT_ANGLE(bPert, aA, aB, aD))
            continue;
        if ( iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aD) == HYDROGEN )
            iAngleWith++;
        else
            iAngleWithout++;
    }
     /*NTHETH*/ FortranWriteInt(iAngleWith);
     /*NTHETA*/ FortranWriteInt(iAngleWithout);

    /* Count the number of torsions with hydrogens, and without */

    iTorsionWith = 0;
    iTorsionWithout = 0;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++) {
        stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
        aA =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1 - 1)->aAtom;
        aB =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2 - 1)->aAtom;
        aC =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3 - 1)->aAtom;
        aD =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4 - 1)->aAtom;
        if (bPERT_TORSION(bPert, aA, aB, aC, aD))
            continue;
        if (iAtomElement(aA) == HYDROGEN ||
            iAtomElement(aB) == HYDROGEN ||
            iAtomElement(aC) == HYDROGEN || iAtomElement(aD) == HYDROGEN)
            iTorsionWith++;
        else
            iTorsionWithout++;
    }
     /*NPHIH*/ FortranWriteInt(iTorsionWith);
     /*NPHIA*/ FortranWriteInt(iTorsionWithout);

     /*JHPARM*/ FortranWriteInt(0);
     /*JPARM*/ FortranWriteInt(0);

    /* Write the number of excluded atoms */

     /*NEXT*/ FortranWriteInt(iVarArrayElementCount(vaExcludedAtoms));
     /*NTOTRS*/ FortranWriteInt(iVarArrayElementCount(uUnit->vaResidues));

    /* Write the number of bonds/angles/torsions without hydrogens */
    /* PLUS the number of RESTRAINT bonds/angles/torsions */

     /*MBONA*/
        FortranWriteInt(iBondWithout +
                        iUnitRestraintTypeCount(uUnit, RESTRAINTBOND));
     /*MTHETA*/
        FortranWriteInt(iAngleWithout +
                        iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE));
     /*MPHIA*/
        FortranWriteInt(iTorsionWithout +
                        iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION));

    /* Write the number of unique bond types, angle types, torsion types */
    /* Add in the number of RESTRAINT bonds/angles/torsion because */
    /* they will have new parameters */

     /*MUMBND*/
        FortranWriteInt(iParmSetTotalBondParms(uUnit->psParameters) +
                        iUnitRestraintTypeCount(uUnit, RESTRAINTBOND));
     /*MUMANG*/
        FortranWriteInt(iParmSetTotalAngleParms(uUnit->psParameters) +
                        iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE));
     /*NPTRA*/
        FortranWriteInt(iParmSetTotalTorsionParms(uUnit->psParameters) +
                        iParmSetTotalImproperParms(uUnit->psParameters) +
                        iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION));

    /* TODO - have different arrays for different restraint types */
    if (iVarArrayElementCount(uUnit->vaRestraints))
        VP0((" Restraints:  Bond %d  Angle %d  Torsion %d\n",
             iUnitRestraintTypeCount(uUnit, RESTRAINTBOND),
             iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE),
             iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION)));
    else
        VP0((" (no restraints)\n"));

     /*NATYP*/
    FortranWriteInt(iParmSetTotalAtomParms(uUnit->psParameters));
     /*NHB*/
    FortranWriteInt(iParmSetTotalHBondParms(uUnit->psParameters));
     /*IFPERT*/
    if (bPert)
        FortranWriteInt(1);
    else
        FortranWriteInt(0);



    /* Count the number of bonds to be perturbed, and those across the */
    /* perturbation/non-perturbed boundary */

    iCountBondPerturbed = 0;
    iCountBondBoundary = 0;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++) {
        sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
        if ((sbPBond->fFlags & PERTURBED) != 0) {
            iCountBondPerturbed++;
            if ((sbPBond->fFlags & BOUNDARY) != 0) {
                MESSAGE(("Boundary pert bond %d-%d\n",
                         sbPBond->iAtom1, sbPBond->iAtom2));
                iCountBondBoundary++;
            }
        }
    }

    MESSAGE(("Perturbed bonds: %d\n", iCountBondPerturbed));
    MESSAGE(("Perturbed boundary bonds: %d\n", iCountBondBoundary));

    /* Count the number of angles to be perturbed, and those on the */
    /* boundary */

    iCountAnglePerturbed = 0;
    iCountAngleBoundary = 0;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++) {
        saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
        if ((saPAngle->fFlags & PERTURBED) != 0)
            iCountAnglePerturbed++;
        if ((saPAngle->fFlags & BOUNDARY) != 0)
            iCountAngleBoundary++;
    }

    /* Count the number of torsions and impropers to be perturbed */
    /* and those on the boundary */

    iCountTorsionPerturbed = 0;
    iCountTorsionBoundary = 0;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++) {
        stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
        if ((stPTorsion->fFlags & PERTURBED) != 0)
            iCountTorsionPerturbed++;
        if ((stPTorsion->fFlags & BOUNDARY) != 0)
            iCountTorsionBoundary++;
    }

     /*NBPER*/ FortranWriteInt(iCountBondPerturbed);
     /*NGPER*/ FortranWriteInt(iCountAnglePerturbed);
     /*NDPER*/ FortranWriteInt(iCountTorsionPerturbed);
     /*MBPER*/ FortranWriteInt(iCountBondPerturbed - iCountBondBoundary);
     /*MGPER*/ FortranWriteInt(iCountAnglePerturbed - iCountAngleBoundary);
     /*MDPER*/
        FortranWriteInt(iCountTorsionPerturbed - iCountTorsionBoundary);

    /* Save flag for periodic boundary conditions */

     /*IFBOX*/
    if (bUnitUseBox(uUnit)) {
        if (bUnitBoxOct(uUnit))
            FortranWriteInt(2);
        else
            FortranWriteInt(1);
    } else
        FortranWriteInt(0);

    /* Save the number of atoms in the largest residue */

    /*NMXRS*/
    FortranWriteInt(iMaxAtoms);

    /* Save flag for cap information */

    /*IFCAP*/
    if (bUnitUseSolventCap(uUnit))
        FortranWriteInt(1);
    else
        FortranWriteInt(0);

    /*NUMEXTRA*/
    iNumExtra = 0;
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
                cPTemp = sAtomType( PVAI(uUnit->vaAtoms, SAVEATOMt, i )->aAtom );
                if( !strncmp( cPTemp, "EP", 2 )) iNumExtra++;
    }
    FortranWriteInt( iNumExtra );

    FortranEndLine();


    /* -3-  write out the names of the atoms */
    FortranDebug("-3-");

    MESSAGE(("Writing the names of the atoms\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ATOM_NAME");
    FortranWriteString("%FORMAT(20a4)");
    FortranFormat(20, LBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        cPTemp = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sName;
        if (strlen(cPTemp) > 4)
            cPTemp += (strlen(cPTemp) - 4);
        FortranWriteString(cPTemp);
    }
    FortranEndLine();

    /* -4- write out the atomic charges */
    FortranDebug("-4-");

    MESSAGE(("Writing the atomic charges\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CHARGE");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        FortranWriteDouble(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->dCharge *
                           ELECTRONTOKCAL);
    }
    FortranEndLine();

    MESSAGE(("Writing the atomic numbers\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ATOMIC_NUMBER");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        FortranWriteInt(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iElement);
    }
    FortranEndLine();

    /* -5- write out the atomic masses */
    FortranDebug("-5-");

    MESSAGE(("Writing the atomic masses\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG MASS");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
        iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
        ParmSetAtom(uUnit->psParameters, iIndex, sType,
                    &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                    &dRStar14, &dScreenF, &iElement, &iHybridization,
		    sDesc);
        FortranWriteDouble(dMass);
    }
    FortranEndLine();

    /* -6- write out the atomic types */
    FortranDebug("-6-");

    MESSAGE(("Writing the atomic types\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ATOM_TYPE_INDEX");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        iAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iTypeIndex - 1;
        iTemp = *PVAI(vaNBIndex, int, iAtom);
        FortranWriteInt(iTemp + 1);
    }
    FortranEndLine();

    /* -7- write out the starting index into the excluded atom list */
    FortranDebug("-7-");

    MESSAGE(("Writing the starting index into the excluded atom list\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG NUMBER_EXCLUDED_ATOMS");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        FortranWriteInt(*PVAI(vaExcludedCount, int, i));
    }
    FortranEndLine();

    /* -8- Write the index for the position of the non bond type */
    /* of each type */
    FortranDebug("-8-");

    MESSAGE(("writing position of the non bond type of each type\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG NONBONDED_PARM_INDEX");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(vaNBIndexMatrix); i++) {
        FortranWriteInt(*PVAI(vaNBIndexMatrix, int, i));
    }
    FortranEndLine();

    /* -9- Residue labels */
    /* Trim the string down to at most 3 characters by */
    /* taking the last three characters if it is too long */
    FortranDebug("-9-");

    MESSAGE(("Writing the residue labels\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG RESIDUE_LABEL");
    FortranWriteString("%FORMAT(20a4)");
    FortranFormat(20, LBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
        cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sName;
        if (strlen(cPTemp) > 3)
            cPTemp += (strlen(cPTemp) - 3);
        FortranWriteString(cPTemp);
    }
    FortranEndLine();

        if( GDefaults.iHaveResIds && GDefaults.iUseResIds ){
                MESSAGE(("Writing the residue ids\n"));
                FortranFormat(1, "%-80s");
                FortranWriteString("%FLAG RESIDUE_ID");
                FortranWriteString("%FORMAT(10a8)");
                FortranFormat(10, IDFORMAT);
                for (i = 1; i <= iVarArrayElementCount(uUnit->vaResidues); i++) {
                        FortranWriteString(GDefaults.sResidueId[i]);
                }
                FortranEndLine();
        }

    /* -10- Pointer list for all the residues */
    FortranDebug("-10-");

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG RESIDUE_POINTER");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
        FortranWriteInt(PVAI(uUnit->vaResidues,
                             SAVERESIDUEt, i)->iAtomStartIndex);
    }
    FortranEndLine();

    /* -11- Force constants for bonds */
    FortranDebug("-11-");

    MESSAGE(("Writing bond force constants\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG BOND_FORCE_CONSTANT");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalBondParms(uUnit->psParameters); i++) {
        ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0,
                    sDesc);
        FortranWriteDouble(dKb);
    }
    /* Write the RESTRAINT constants AND set the index */
    /* for where the interaction can find its constants */
    RESTRAINTLOOP(RESTRAINTBOND, dKx, i + 1);
    FortranEndLine();

    /* -12- Equilibrium bond lengths */
    FortranDebug("-12-");

    MESSAGE(("Writing equilibrium bond lengths\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG BOND_EQUIL_VALUE");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalBondParms(uUnit->psParameters); i++) {
        ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0,
                    sDesc);
        FortranWriteDouble(dR0);
    }
    /* Write the bond RESTRAINT constants AND set the index */
    /* for where the interaction can find its constants */
    RESTRAINTLOOP(RESTRAINTBOND, dX0, i + 1);
    FortranEndLine();

    /* -13- Force constants for angles */
    FortranDebug("-13-");

    MESSAGE(("Writing angle force constants\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ANGLE_FORCE_CONSTANT");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
        ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                     &dKt, &dT0, &dTkub, &dRkub, sDesc);
        FortranWriteDouble(dKt);
    }
    /* Write the angle RESTRAINT constants AND set the index */
    /* for where the interaction can find its constants */
    RESTRAINTLOOP(RESTRAINTANGLE, dKx, i + 1);
    FortranEndLine();

    /* -14- Equilibrium angle values */
    FortranDebug("-14-");

    MESSAGE(("Writing equilibrium angle values\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ANGLE_EQUIL_VALUE");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
        ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                     &dKt, &dT0, &dTkub, &dRkub, sDesc);
        FortranWriteDouble(dT0);
    }
    /* Write the angle RESTRAINT constants AND set the index */
    /* for where the interaction can find its constants */
    RESTRAINTLOOP(RESTRAINTANGLE, dX0, i + 1);
    FortranEndLine();

    /* -15- Force constants for torsions */
    FortranDebug("-15-");

    MESSAGE(("Writing torsional force constants\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG DIHEDRAL_FORCE_CONSTANT");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    MESSAGE(("There are %d torsions and %d impropers\n",
             iParmSetTotalTorsionParms(uUnit->psParameters),
             iParmSetTotalImproperParms(uUnit->psParameters)));
    for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
        ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2,
                       sAtom3, sAtom4, &iN, &dKp, &dP0, &dScEE,
		       &dScNB, sDesc);
        MESSAGE(("Torsion %d  %s-%s-%s-%s %d %lf %lf\n",
                 i, sAtom1, sAtom2, sAtom3, sAtom4, iN, dKp, dP0));
        FortranWriteDouble(dKp);
    }
    for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
        ParmSetImproper(uUnit->psParameters, i, sAtom1, sAtom2,
                        sAtom3, sAtom4, &iN, &dKp, &dP0, sDesc);
        MESSAGE(("Improper %d  %s-%s-%s-%s %d %lf %lf\n",
                 i, sAtom1, sAtom2, sAtom3, sAtom4, iN, dKp, dP0));
        FortranWriteDouble(dKp);
    }
    /* Write the torsion RESTRAINT constants AND set the index */
    /* for where the interaction can find its constants */
    RESTRAINTLOOP(RESTRAINTTORSION, dKx, i + 1);
    FortranEndLine();

    /* -16- Periodicity for the dihedral angles */
    FortranDebug("-16-");

    MESSAGE(("Writing periodicity of torsion interaction\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG DIHEDRAL_PERIODICITY");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
        ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2,
                       sAtom3, sAtom4, &iN, &dKp, &dP0, &dScEE,
		       &dScNB, sDesc);
        dTemp = iN;
        FortranWriteDouble(dTemp);
    }
    for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
        ParmSetImproper(uUnit->psParameters, i, sAtom1, sAtom2,
                        sAtom3, sAtom4, &iN, &dKp, &dP0, sDesc);
        dTemp = iN;
        FortranWriteDouble(dTemp);
    }
    /* Write the torsion RESTRAINT constants AND set the index */
    /* for where the interaction can find its constants */
    RESTRAINTLOOP(RESTRAINTTORSION, dX0, i + 1);
    FortranEndLine();

    /* -17- Phase for torsions */
    FortranDebug("-17-");

    MESSAGE(("Writing phase for torsion interactions\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG DIHEDRAL_PHASE");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
        ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2,
                       sAtom3, sAtom4, &iN, &dKp, &dP0, &dScEE,
		       &dScNB, sDesc);
        FortranWriteDouble(dP0);
    }
    for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
        ParmSetImproper(uUnit->psParameters, i, sAtom1, sAtom2,
                        sAtom3, sAtom4, &iN, &dKp, &dP0, sDesc);
        FortranWriteDouble(dP0);
    }
    /* Write the torsion RESTRAINT constants AND set the index */
    /* for where the interaction can find its constants */
    RESTRAINTLOOP(RESTRAINTTORSION, dN, i + 1);
    FortranEndLine();

/*    if (GDefaults.dSceeScaleFactor > 0.0) {  */
//    1-4 NB and 1-4 EEL scaling factors are written regardless.
        /* -17B- */
        FortranDebug("-17B-");
        
        MESSAGE(("Writing SCEE_SCALE_FACTOR torsion\n"));
        //dSceeScaleFactor = 0.0;
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG SCEE_SCALE_FACTOR");
        FortranWriteString("%FORMAT(5E16.8)");
        FortranFormat(5, DBLFORMAT);
        for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
            ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2,
                           sAtom3, sAtom4, &iN, &dKp, &dP0, &dScEE,
			   &dScNB, sDesc);
            if ( dScEE < 0.0 ) dScEE = GDefaults.dSceeScaleFactor;
            FortranWriteDouble(dScEE);
        }
        for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
            FortranWriteDouble(0.0);
        }
        FortranEndLine();
        
        /* -17C- */
        FortranDebug("-17C-");
        
        MESSAGE(("Writing SCNB_SCALE_FACTOR torsion\n"));
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG SCNB_SCALE_FACTOR");
        FortranWriteString("%FORMAT(5E16.8)");
        FortranFormat(5, DBLFORMAT);
        for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
            ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2,
                           sAtom3, sAtom4, &iN, &dKp, &dP0, &dScEE,
			   &dScNB, sDesc);
	    if ( dScNB < 0.0 ) dScNB = GDefaults.dScnbScaleFactor;
            FortranWriteDouble(dScNB);
        }

        for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
            FortranWriteDouble(0.0);
        }
        FortranEndLine();
/*
    }
 */

    /* -18- Not used, reserved for future use, uses NATYP */
    /* Corresponds to the AMBER SOLTY array */
    FortranDebug("-18-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG SOLTY");
    FortranWriteString("%FORMAT(5E16.8)");

    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalAtomParms(uUnit->psParameters); i++) {
        FortranWriteDouble(0.0);
    }
    FortranEndLine();

    /* -19- Lennard jones r**12 term for all possible interactions */
    /* CN1 array */
    FortranDebug("-19-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG LENNARD_JONES_ACOEF");
    FortranWriteString("%FORMAT(5E16.8)");

    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
        FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dA);
    }
    FortranEndLine();

    /* -20- Lennard jones r**6 term for all possible interactions */
    /* CN2 array */
    FortranDebug("-20-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG LENNARD_JONES_BCOEF");
    FortranWriteString("%FORMAT(5E16.8)");

    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
        FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dC);
    }
    FortranEndLine();

    /* -21- Write the bond interactions that include hydrogen */
    /* Write the two indices into the atom table, then the index */
    /* into the interaction table */
    FortranDebug("-21-");

    MESSAGE(("Writing the bond interactions with hydrogens\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG BONDS_INC_HYDROGEN");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++) {
        sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
        aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1 - 1)->aAtom;
        aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2 - 1)->aAtom;
        if (bPERT_BOND(bPert, aA, aD))
            continue;
        if (iAtomElement(aA) == HYDROGEN || iAtomElement(aD) == HYDROGEN) {
            FortranWriteInt(AMBERINDEX(sbPBond->iAtom1));
            FortranWriteInt(AMBERINDEX(sbPBond->iAtom2));
            FortranWriteInt(sbPBond->iParmIndex);
        }
    }
    FortranEndLine();

    /* -22- Write the bond interactions that dont include hydrogen */
    /* Write the two indices into the atom table, then the index */
    /* into the interaction table */
    FortranDebug("-22-");

    MESSAGE(("Writing the bond interactions without hydrogens\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG BONDS_WITHOUT_HYDROGEN");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++) {
        sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
        aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1 - 1)->aAtom;
        aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2 - 1)->aAtom;
        if (bPERT_BOND(bPert, aA, aD))
            continue;
        if (!(iAtomElement(aA) == HYDROGEN ||
              iAtomElement(aD) == HYDROGEN)) {
            FortranWriteInt(AMBERINDEX(sbPBond->iAtom1));
            FortranWriteInt(AMBERINDEX(sbPBond->iAtom2));
            FortranWriteInt(sbPBond->iParmIndex);
        }
    }
    /* Write out the (bond without H) RESTRAINT interactions */
    /* The iParmIndex field is set in RESTRAINTLOOP */
    if ((iMax = iVarArrayElementCount(uUnit->vaRestraints))) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        for (i = 0; i < iMax; i++, srPRestraint++) {
            if (srPRestraint->iType == RESTRAINTBOND) {
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom1));
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom2));
                FortranWriteInt(srPRestraint->iParmIndex);
            }
        }
    }
    FortranEndLine();

    /* -23- Write the angle interactions that include hydrogen */
    /* Write the three indices into the atom table, then the index */
    /* into the interaction table */
    FortranDebug("-23-");

    MESSAGE(("Writing the angle interactions with hydrogens\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ANGLES_INC_HYDROGEN");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++) {
        saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
        aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1 - 1)->aAtom;
        aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2 - 1)->aAtom;
        aD = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3 - 1)->aAtom;
        if (bPERT_ANGLE(bPert, aA, aB, aD))
            continue;
        if ( iAtomElement(aA) == HYDROGEN
                || iAtomElement(aB) == HYDROGEN
                || iAtomElement(aD) == HYDROGEN ) {
            FortranWriteInt(AMBERINDEX(saPAngle->iAtom1));
            FortranWriteInt(AMBERINDEX(saPAngle->iAtom2));
            FortranWriteInt(AMBERINDEX(saPAngle->iAtom3));
            FortranWriteInt(saPAngle->iParmIndex);
        }
    }
    FortranEndLine();

    /* -24- Write the angle interactions that dont include hydrogen */
    /* Write the three indices into the atom table, then the index */
    /* into the interaction table */
    FortranDebug("-24-");

    MESSAGE(("Writing the angle interactions without hydrogens\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ANGLES_WITHOUT_HYDROGEN");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++) {
        saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
        aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1 - 1)->aAtom;
        aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2 - 1)->aAtom;
        aD = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3 - 1)->aAtom;
        if (bPERT_ANGLE(bPert, aA, aB, aD))
            continue;
        if ( !(iAtomElement(aA) == HYDROGEN
                || iAtomElement(aB) == HYDROGEN
                || iAtomElement(aD) == HYDROGEN) ) {
            FortranWriteInt(AMBERINDEX(saPAngle->iAtom1));
            FortranWriteInt(AMBERINDEX(saPAngle->iAtom2));
            FortranWriteInt(AMBERINDEX(saPAngle->iAtom3));
            FortranWriteInt(saPAngle->iParmIndex);
        }
    }
    /* Write out the RESTRAINT interactions */
    /* The iParmIndex field is set in RESTRAINTLOOP */
    if ((iMax = iVarArrayElementCount(uUnit->vaRestraints))) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        for (i = 0; i < iMax; i++, srPRestraint++) {
            if (srPRestraint->iType == RESTRAINTANGLE) {
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom1));
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom2));
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom3));
                FortranWriteInt(srPRestraint->iParmIndex);
            }
        }
    }
    FortranEndLine();

    /* -25- Write the torsion interactions that include hydrogen */
    /* Write the three indices into the atom table, then the index */
    /* into the interaction table */
    FortranDebug("-25-");

    MESSAGE(("Writing the torsion interactions with hydrogens\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG DIHEDRALS_INC_HYDROGEN");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++) {
        stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
        aA =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1 - 1)->aAtom;
        aB =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2 - 1)->aAtom;
        aC =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3 - 1)->aAtom;
        aD =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4 - 1)->aAtom;
        if (bPERT_TORSION(bPert, aA, aB, aC, aD))
            continue;
        if (iAtomElement(aA) == HYDROGEN
            || iAtomElement(aB) == HYDROGEN
            || iAtomElement(aC) == HYDROGEN
            || iAtomElement(aD) == HYDROGEN) {
            if ((AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                (AMBERINDEX(stPTorsion->iAtom4) == 0)) {
                MESSAGE(
                        ("Had to turn torsion around to avoid K,L == 0\n"));
                MESSAGE(
                        ("Outer atoms: %s --- %s\n", sContainerName(aA),
                         sContainerName(aD)));
                MESSAGE(
                        ("Old order %d %d %d %d\n", stPTorsion->iAtom1,
                         stPTorsion->iAtom2, stPTorsion->iAtom3,
                         stPTorsion->iAtom4));
                SWAP(stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp);
                SWAP(stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp);
                MESSAGE(("New order %d %d %d %d\n",
                         stPTorsion->iAtom1,
                         stPTorsion->iAtom2,
                         stPTorsion->iAtom3, stPTorsion->iAtom4));
            }
            if (stPTorsion->bProper)
                iProper = 1;
            else
                iProper = -1;
            if (stPTorsion->bCalc14)
                iCalc14 = 1;
            else
                iCalc14 = -1;
            if (GDefaults.iCharmm && iProper == -1) {
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1) * iCalc14);
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
            } else {
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3) * iCalc14);
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
            }
            FortranWriteInt(stPTorsion->iParmIndex);
        }
    }
    FortranEndLine();

    /* -26- Write the torsion interactions that dont include hydrogen */
    /* Write the three indices into the atom table, then the index */
    /* into the interaction table */
    FortranDebug("-26-");

    MESSAGE(("Writing the torsion interactions without hydrogens\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG DIHEDRALS_WITHOUT_HYDROGEN");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++) {
        stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
        aA =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1 - 1)->aAtom;
        aB =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2 - 1)->aAtom;
        aC =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3 - 1)->aAtom;
        aD =
            PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4 - 1)->aAtom;
        if (bPERT_TORSION(bPert, aA, aB, aC, aD))
            continue;
        if (!(iAtomElement(aA) == HYDROGEN
              || iAtomElement(aB) == HYDROGEN
              || iAtomElement(aC) == HYDROGEN
              || iAtomElement(aD) == HYDROGEN)) {
            if ((AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                (AMBERINDEX(stPTorsion->iAtom4) == 0)) {
                MESSAGE(("Had to turn torsion to avoid K,L == 0\n"));
                MESSAGE(("Outer atoms: %s --- %s\n",
                         sContainerName(aA), sContainerName(aD)));
                SWAP(stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp);
                SWAP(stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp);
            }
            if (stPTorsion->bCalc14)
                iCalc14 = 1;
            else
                iCalc14 = -1;
            if (stPTorsion->bProper)
                iProper = 1;
            else
                iProper = -1;
            if (GDefaults.iCharmm && iProper == -1) {
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1) * iCalc14);
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
            } else {
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3) * iCalc14);
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
            }
            FortranWriteInt(stPTorsion->iParmIndex);
        }
    }
    /* Write out the RESTRAINT interactions */
    /* The iParmIndex field is set in RESTRAINTLOOP */
    if ((iMax = iVarArrayElementCount(uUnit->vaRestraints))) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        for (i = 0; i < iMax; i++, srPRestraint++) {
            if (srPRestraint->iType == RESTRAINTTORSION) {
                if ((AMBERINDEX(srPRestraint->iAtom3) == 0) ||
                    (AMBERINDEX(srPRestraint->iAtom4) == 0)) {
                    MESSAGE(
                            ("Had to turn RESTRAINT torsion around to avoid\n"));
                    MESSAGE(("K,L == 0\n"));
                    SWAP(srPRestraint->iAtom1, srPRestraint->iAtom4,
                         iTemp);
                    SWAP(srPRestraint->iAtom2, srPRestraint->iAtom3,
                         iTemp);
                }
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom1));
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom2));
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom3));
                FortranWriteInt(AMBERINDEX(srPRestraint->iAtom4));
                FortranWriteInt(srPRestraint->iParmIndex);
            }
        }
    }
    FortranEndLine();

    /* -27- Write the excluded atom list */
    FortranDebug("-27-");

    MESSAGE(("Writing the excluded atom list\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG EXCLUDED_ATOMS_LIST");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(vaExcludedAtoms); i++) {
        FortranWriteInt(*PVAI(vaExcludedAtoms, int, i));
    }
    FortranEndLine();

    /* -28- Write the R^12 term for the Hydrogen bond equation */
    FortranDebug("-28-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG HBOND_ACOEF");
    FortranWriteString("%FORMAT(5E16.8)");

    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
        ParmSetHBond(uUnit->psParameters, i, sType1, sType2, &dC, &dD,
                     sDesc);
        FortranWriteDouble(dC);
    }
    FortranEndLine();

    /* -29- Write the R^10 term for the Hydrogen bond equation */
    FortranDebug("-29-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG HBOND_BCOEF");
    FortranWriteString("%FORMAT(5E16.8)");

    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
        ParmSetHBond(uUnit->psParameters, i, sType1, sType2, &dC, &dD,
                     sDesc);
        FortranWriteDouble(dD);
    }
    FortranEndLine();

    /* -30- No longer used, but stored */
    FortranDebug("-30-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG HBCUT");
    FortranWriteString("%FORMAT(5E16.8)");

    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
        FortranWriteDouble(0.0);
    }
    FortranEndLine();

    /* -31- List of atomic symbols */
    FortranDebug("-31-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG AMBER_ATOM_TYPE");
    FortranWriteString("%FORMAT(20a4)");

    FortranFormat(20, LBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        FortranWriteString(sAtomType
                           (PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom));
    }
    FortranEndLine();

    /* -32- List of tree symbols */
    FortranDebug("-32-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG TREE_CHAIN_CLASSIFICATION");
    FortranWriteString("%FORMAT(20a4)");

    FortranFormat(20, LBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
        if (dAtomTemp(aAtom) == (double) 'M')
            FortranWriteString("M  ");
        else if (dAtomTemp(aAtom) == (double) 'E')
            FortranWriteString("E  ");
        else if (dAtomTemp(aAtom) == (double) 'S')
            FortranWriteString("S  ");
        else if (dAtomTemp(aAtom) == (double) 'B')
            FortranWriteString("B  ");
        else if (dAtomTemp(aAtom) == (double) '3')
            FortranWriteString("3  ");
        else if (dAtomTemp(aAtom) == (double) '4')
            FortranWriteString("4  ");
        else if (dAtomTemp(aAtom) == (double) '5')
            FortranWriteString("5  ");
        else if (dAtomTemp(aAtom) == (double) '6')
            FortranWriteString("6  ");
        else if (dAtomTemp(aAtom) == (double) 'X')
            FortranWriteString("X  ");
        else
            FortranWriteString("BLA");
    }
    FortranEndLine();

    /* -33- Tree Joining information !!!!!!! Add support for this !!!!! */
    FortranDebug("-33-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG JOIN_ARRAY");
    FortranWriteString("%FORMAT(10I8)");

    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        FortranWriteInt(0);
    }
    FortranEndLine();

    /* -34- Who knows, something to do with rotating atoms */
    FortranDebug("-34-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG IROTAT");
    FortranWriteString("%FORMAT(10I8)");

    FortranFormat(10, INTFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        FortranWriteInt(0);
    }
    FortranEndLine();

    /* -35A- The last residue before "solvent" */
    /* Number of molecules */
    /* Index of first molecule that is solvent */

    if (bUnitUseBox(uUnit)) {
        FortranDebug("-35A-");

        /* Find the index of the first solvent RESIDUE */

        for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
            if (PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->
                sResidueType[0] == RESTYPESOLVENT) break;
        }
        iTemp = i;

        /* 
         *  Find the molecules and return the number of ATOMs in each 
         *  molecule, along with the index of the first solvent molecule
         */

        vaMolecules = vaVarArrayCreate(sizeof(int));
        zUnitIOFindAndCountMolecules(uUnit, &vaMolecules, &iFirstSolvent);

        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG SOLVENT_POINTERS");
        FortranWriteString("%FORMAT(3I8)");
        FortranFormat(3, INTFORMAT);
        FortranWriteInt(iTemp);
        FortranWriteInt(iVarArrayElementCount(vaMolecules));
        FortranWriteInt(iFirstSolvent + 1);        /* FORTRAN index */

        FortranEndLine();

        /* -35B- The number of ATOMs in the Ith RESIDUE */

        FortranDebug("-35B-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG ATOMS_PER_MOLECULE");
        FortranWriteString("%FORMAT(10I8)");
        FortranFormat(10, INTFORMAT);
        for (i = 0; i < iVarArrayElementCount(vaMolecules); i++) {
            FortranWriteInt(*PVAI(vaMolecules, int, i));
        }
        FortranEndLine();

        /* -35C- BETA, (BOX(I), I=1,3 ) */

        FortranDebug("-35C-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG BOX_DIMENSIONS");
        FortranWriteString("%FORMAT(5E16.8)");
        FortranFormat(4, DBLFORMAT);
        FortranWriteDouble(dUnitBeta(uUnit) / DEGTORAD);
        UnitGetBox(uUnit, &dX, &dY, &dZ);
        FortranWriteDouble(dX);
        FortranWriteDouble(dY);
        FortranWriteDouble(dZ);
        FortranEndLine();
    }

    /* -35D- NATCAP */

    if (bUnitUseSolventCap(uUnit)) {
        FortranDebug("-35D-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG CAP_INFO");
        FortranWriteString("%FORMAT(10I8)");
        FortranFormat(1, INTFORMAT);
        FortranWriteInt(uUnit->iCapTempInt);
        FortranEndLine();

        /* -35E- CUTCAP, XCAP, YCAP, ZCAP */
        FortranDebug("-35E-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG CAP_INFO2");
        FortranWriteString("%FORMAT(5E16.8)");
        FortranFormat(4, DBLFORMAT);
        UnitGetSolventCap(uUnit, &dX, &dY, &dZ, &dR);
        FortranWriteDouble(dR);
        FortranWriteDouble(dX);
        FortranWriteDouble(dY);
        FortranWriteDouble(dZ);
        FortranEndLine();
    }

    /* write out the GB radii */

    MESSAGE(("Writing the GB radii\n"));

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG RADIUS_SET");
    FortranWriteString("%FORMAT(1a80)");
    switch (GDefaults.iGBparm) {
        case 0: FortranWriteString("Bondi radii (bondi)"); break;
        case 1: FortranWriteString("amber6 modified Bondi radii (amber6)"); break;
        case 2: FortranWriteString("modified Bondi radii (mbondi)"); break;
#if 0
        case 3: FortranWriteString("Huo and Kollman optimized radii (pbamber)"); break;
#endif
        case 6: FortranWriteString("H(N)-modified Bondi radii (mbondi2)"); break;
        case 7: FortranWriteString("Parse radii (parse)"); break;
        case 8: FortranWriteString("ArgH and AspGluO modified Bondi2 radii (mbondi3)"); break;
        default: FortranWriteString("Unknown radius set (leap needs to be modified!)"); break;
    }  
    FortranEndLine();
    FortranWriteString("%FLAG RADII");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    iResidueIndex = -1;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
        iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
        ParmSetAtom(uUnit->psParameters, iIndex, sType,
                    &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                    &dRStar14, &dScreenF, &iElement, &iHybridization,
		    sDesc);

        if( GDefaults.iGBparm < 3  || GDefaults.iGBparm == 6 || GDefaults.iGBparm == 8 ) {
            /* Bondi or modified Bondi radii */
            switch( iElement ) {
                case  1: dGBrad = 1.2; 
                    /* make the modifications that hydrogen radii
                       depend upon the atoms they are bonded to.  
                       iGBparm=1 corresponds to Amber 6, JACS 122:2489 (2000);
                       iGBparm=2 adds the update of Biopolymers 56: 275 (2001)   
                     */
                    if( iAtomCoordination(saPAtom->aAtom) > 0 ) {
                        /* For multiply bonded Hydrogen atoms use the first
                         * bond for determining modified GB radii.
                         * WAT contains multiply bonded Hydrogen atoms 
                         * so do not emit a warning.
                         */

                        aAtomA = aAtomBondedNeighbor(saPAtom->aAtom, 0);
                        if( GDefaults.iGBparm == 1 || GDefaults.iGBparm == 2 ) {
                            if( sAtomType(aAtomA)[0] == 'C' ||
                                    sAtomType(aAtomA)[0] == 'c' ) dGBrad = 1.3;
                            if( sAtomType(aAtomA)[0] == 'O' ||
                                    sAtomType(aAtomA)[0] == 'o' ) dGBrad = 0.8;
                            if( sAtomType(aAtomA)[0] == 'S' ||
                                    sAtomType(aAtomA)[0] == 's' ) dGBrad = 0.8;
                            if( (sAtomType(aAtomA)[0] == 'N' ||
                                    sAtomType(aAtomA)[0] == 'n')  &&
                                    GDefaults.iGBparm == 2) dGBrad = 1.3;
                            if( (sAtomType(aAtomA)[0] == 'H' ||
                                    sAtomType(aAtomA)[0] == 'h')  &&
                                (sAtomType(aAtomA)[1] == 'W' ||
                                    sAtomType(aAtomA)[1] == 'w')) dGBrad = 0.8; 
                        }
                        else if( GDefaults.iGBparm == 6 || GDefaults.iGBparm == 8 ) { 
                            /* try Alexey's scheme */
                            if( sAtomType(aAtomA)[0] == 'N' ||
                                     sAtomType(aAtomA)[0] == 'n' ) {
                                dGBrad = 1.3;
                                if ( GDefaults.iGBparm == 8 ) {
                                    // update residue as appropriate
                                    if( saPAtom->iResidueIndex != iResidueIndex ) {
                                        iResidueIndex = saPAtom->iResidueIndex;
                                        cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt,
                                                      iResidueIndex-1)->sName;
                                        if (strlen(cPTemp) > 3)
                                            cPTemp += (strlen(cPTemp) - 3);
                                    }
                                    // adjust Arg HH and HE
                                    //VP0(( "%3d. Is %s ARG?\n", i, cPTemp ));
                                    if( !strcmp(cPTemp, "ARG") &&
                                        !(strncmp(sAtomName(saPAtom->aAtom),"HH",2) &&
                                          strcmp(sAtomName(saPAtom->aAtom),"HE"))) {
                                        dGBrad = 1.17;
                                    }
                                }
                            }
                        }
                    }
                    else {
                        VP0(( "WARNING: Unbonded Hydrogen atom %s in %s.\n"
                                " Cannot determine the requested GB radius"
                                " for this atom.\n"
                                " Writing the unmodified Bondi GB radius.\n",
                                saPAtom->aAtom->cHeader.sName,
                                saPAtom->aAtom->cHeader.cContainedBy->sName ));
                    }
                    break;
                case  6:
                    /* Use the mass of the carbon atom. We are testing for
                     * carbons here. C1 == CH, C2 == CH2, C3 == CH3. UA carbons
                     * have a larger radius (2.2), so we want to make sure that
                     * the C1, C2, and C3 atom types _really_ correspond to UA
                     * UA carbons. C1 atoms should have a mass of 12.01 + 1.01,
                     * C2 should be 12.01 + 2.02, and C3 should be 12.01 + 3.03.
                     * This mneumonic will not work for 13C named "C1". This is
                     * a (hopefully) temporary kludge.
                     *              --JMS 11/4/2013
                     */
                    if ( strncmp(sType,"C1",2) && strncmp(sType,"C2",2) && strncmp (sType,"C3",2) )
                       dGBrad = 1.7; 
                    else if ( !strncmp(sType, "C1", 2) && dMass < 13.0 )
                       dGBrad = 1.7;
                    else if ( !strncmp(sType, "C2", 2) && dMass < 14.0 )
                       dGBrad = 1.7;
                    else if ( !strncmp(sType, "C3", 2) && dMass < 15.0 )
                       dGBrad = 1.7;
                    else
                       dGBrad = 2.2;
                    break;
                case  7: dGBrad = 1.55; break;
                case  8:
                    dGBrad = 1.5;
                    if ( GDefaults.iGBparm == 8 ) {
                        // update residue as appropriate
                        if( saPAtom->iResidueIndex != iResidueIndex ) {
                            iResidueIndex = saPAtom->iResidueIndex;
                            cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt,
                                          iResidueIndex-1)->sName;
                            if (strlen(cPTemp) > 3)
                                cPTemp += (strlen(cPTemp) - 3);
                        }
                        // adjust Asp OD and Glu OE, and terminal OXT
                        if( !(strcmp(cPTemp, "ASP") || strncmp(sAtomName(saPAtom->aAtom),"OD",2)) ||
                            !(strcmp(cPTemp, "AS4") || strncmp(sAtomName(saPAtom->aAtom),"OD",2)) ||
                            !(strcmp(cPTemp, "GLU") || strncmp(sAtomName(saPAtom->aAtom),"OE",2)) ||
                            !(strcmp(cPTemp, "GL4") || strncmp(sAtomName(saPAtom->aAtom),"OE",2)) ||
                            (!strcmp(sAtomName(saPAtom->aAtom),"OXT") ||
                             (i+1<iVarArrayElementCount(uUnit->vaAtoms)&&
                              !strcmp(sAtomName(PVAI(uUnit->vaAtoms, SAVEATOMt, i+1)->aAtom),
                                      "OXT")))){
                            dGBrad = 1.4;
                        }
                    }
                    break;
                case  9: dGBrad = 1.5; break;
                case 14: dGBrad = 2.1; break;
                case 15: dGBrad = 1.85; break;
                case 16: dGBrad = 1.8; break;
                case 17: dGBrad = 1.7; break;
                default: dGBrad = 1.5; break;
            }
        } else if ( GDefaults.iGBparm == 3 ) {  /* radii from Huo & Kollman */
            switch( iElement ) {
                case  1: dGBrad = 1.15; break;
                case  6: dGBrad = 1.85; break;
                case  7: dGBrad = 1.64; break;
                case  8: dGBrad = 1.53; break;
                case  9: dGBrad = 1.53; break;
                case 15: dGBrad = 2.02; break;
                case 16: dGBrad = 2.00; break;
                case 17: dGBrad = 1.97; break;
                case 35: dGBrad = 2.03; break;
                case 53: dGBrad = 2.10; break;
                default: dGBrad = 1.5; break;  /* DAC made this up!  */
            }
        } else if ( GDefaults.iGBparm == 7 ) { /* Parse radii */
            switch( iElement ) {
                case  1: dGBrad = 1.00; break;
                case  6: dGBrad = 1.70; break;
                case  7: dGBrad = 1.50; break;
                case  8: dGBrad = 1.40; break;
                case 16: dGBrad = 1.85; break;
                default: dGBrad = 1.50; break;  /* Radii from J. Phys. Chem. 1994, 98, 1978-1988 */
            }
        }
        FortranWriteDouble(dGBrad);
    }
    FortranEndLine();

    /* write out the GB screening parameters */

    MESSAGE(("Writing the GB screening parameters\n"));
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG SCREEN");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
        iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
        ParmSetAtom(uUnit->psParameters, iIndex, sType,
                    &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                    &dRStar14, &dScreenF, &iElement, &iHybridization,
		    sDesc);
        if( GDefaults.iGBparm < 4 || GDefaults.iGBparm == 6 || GDefaults.iGBparm == 8 ){
            /* for now, hardwire the Bondi radii  */
            switch( iElement ){
                case 1:  dGBscreen = 0.85; break;
                case 6:  dGBscreen = 0.72; break;
                case 7:  dGBscreen = 0.79; break;
                case 8:  dGBscreen = 0.85; break;
                case 9:  dGBscreen = 0.88; break;
                case 15: dGBscreen = 0.86; break;
                case 16: dGBscreen = 0.96; break;
                default: dGBscreen = 0.8; break;  /* or should fail?? */
            }
        } else if( GDefaults.iGBparm == 4 ){ /* param for Jayaram et al. 'GB' */
            switch( iElement ){
                case 1:  dGBscreen = 0.8461; break;
                case 6:  dGBscreen = 0.9615; break;
                case 7:  dGBscreen = 0.9343; break;
                case 8:  dGBscreen = 1.0088; break;
                case 11: dGBscreen = 1.0000; break; 
                case 12: dGBscreen = 1.0000; break; /* set by HG */
                case 15: dGBscreen = 1.0700; break; 
                case 16: dGBscreen = 1.1733; break; 
                default: dGBscreen = 0.8000; break; /* set by HG */
            }
        } else if ( GDefaults.iGBparm == 5 ){  /* param for Jayaram et al. 'MGB' */
            switch( iElement ){
                case 1:  dGBscreen = 0.8846; break;
                case 6:  dGBscreen = 0.9186; break;
                case 7:  dGBscreen = 0.8733; break;
                case 8:  dGBscreen = 0.8836; break;
                case 11: dGBscreen = 1.0000; break;
                case 12: dGBscreen = 1.0000; break; /* set by HG */
                case 15: dGBscreen = 0.9604; break;
                case 16: dGBscreen = 0.9323; break;
                default: dGBscreen = 0.8000; break; /* set by HG */
            }
        }
        FortranWriteDouble(dGBscreen);
    }
    FortranEndLine();

//
// write IPOL near the end of prmtop
//
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG IPOL");
    FortranWriteString("%FORMAT(1I8)");
    FortranFormat(1, INTFORMAT);
    FortranWriteInt(GDefaults.iIPOL);
    FortranEndLine();


    /* Write the perturbation information */

    if (bPert) {

        /* -36A- Bonds that are to be perturbed */
        /* Totally perturbed bonds first, */
        /* boundary second */
        FortranDebug("-36A-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_BOND_ATOMS");
        FortranWriteString("%FORMAT(10I8)");

        FortranFormat(10, INTFORMAT);
        if (iVarArrayElementCount(uUnit->vaBonds))
            sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds);
             i++, sbPBond++) {
            if (((sbPBond->fFlags & PERTURBED) != 0)
                && ((sbPBond->fFlags & BOUNDARY) == 0)) {
                FortranWriteInt(AMBERINDEX(sbPBond->iAtom1));
                FortranWriteInt(AMBERINDEX(sbPBond->iAtom2));
            }
        }
        if (iVarArrayElementCount(uUnit->vaBonds))
            sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds);
             i++, sbPBond++) {
            if (((sbPBond->fFlags & PERTURBED) != 0)
                && ((sbPBond->fFlags & BOUNDARY) != 0)) {
                FortranWriteInt(AMBERINDEX(sbPBond->iAtom1));
                FortranWriteInt(AMBERINDEX(sbPBond->iAtom2));
            }
        }
        FortranEndLine();

        /* -36B- Index into bond interaction arrays */
        FortranDebug("-36B-");

        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_BOND_PARAMS");
        FortranWriteString("%FORMAT(10I8)");
        FortranFormat(10, INTFORMAT);


        /* First LAMBDA = 0 */
        
        if (iVarArrayElementCount(uUnit->vaBonds))
            sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds);
             i++, sbPBond++) {
            if (((sbPBond->fFlags & PERTURBED) != 0)
                && ((sbPBond->fFlags & BOUNDARY) == 0)) {
                FortranWriteInt(sbPBond->iParmIndex);
            }
        }
        if (iVarArrayElementCount(uUnit->vaBonds))
            sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds);
             i++, sbPBond++) {
            if (((sbPBond->fFlags & PERTURBED) != 0)
                && ((sbPBond->fFlags & BOUNDARY) != 0)) {
                FortranWriteInt(sbPBond->iParmIndex);
            }
        }

        /* Then LAMBDA = 1 */

        if (iVarArrayElementCount(uUnit->vaBonds))
            sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds);
             i++, sbPBond++) {
            if (((sbPBond->fFlags & PERTURBED) != 0)
                && ((sbPBond->fFlags & BOUNDARY) == 0)) {
                FortranWriteInt(sbPBond->iPertParmIndex);
            }
        }
        if (iVarArrayElementCount(uUnit->vaBonds))
            sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds);
             i++, sbPBond++) {
            if (((sbPBond->fFlags & PERTURBED) != 0)
                && ((sbPBond->fFlags & BOUNDARY) != 0)) {
                FortranWriteInt(sbPBond->iPertParmIndex);
            }
        }
        FortranEndLine();


        /* -36C- Angles that are to be perturbed */
        FortranDebug("-36C-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_ANGLE_ATOMS");
        FortranWriteString("%FORMAT(10I8)");

        FortranFormat(10, INTFORMAT);
        if (iVarArrayElementCount(uUnit->vaAngles))
            saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles);
             i++, saPAngle++) {
            if (((saPAngle->fFlags & PERTURBED) != 0)
                && ((saPAngle->fFlags & BOUNDARY) == 0)) {
                FortranWriteInt(AMBERINDEX(saPAngle->iAtom1));
                FortranWriteInt(AMBERINDEX(saPAngle->iAtom2));
                FortranWriteInt(AMBERINDEX(saPAngle->iAtom3));
            }
        }
        if (iVarArrayElementCount(uUnit->vaAngles))
            saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles);
             i++, saPAngle++) {
            if (((saPAngle->fFlags & PERTURBED) != 0)
                && ((saPAngle->fFlags & BOUNDARY) != 0)) {
                FortranWriteInt(AMBERINDEX(saPAngle->iAtom1));
                FortranWriteInt(AMBERINDEX(saPAngle->iAtom2));
                FortranWriteInt(AMBERINDEX(saPAngle->iAtom3));
            }
        }
        FortranEndLine();

        /* -36D- Index into angle interaction arrays */
        FortranDebug("-36D-");

        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_ANGLE_PARAMS");
        FortranWriteString("%FORMAT(10I8)");

        FortranFormat(10, INTFORMAT);

        /* First LAMBDA = 0 */

        if (iVarArrayElementCount(uUnit->vaAngles))
            saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles);
             i++, saPAngle++) {
            if (((saPAngle->fFlags & PERTURBED) != 0)
                && ((saPAngle->fFlags & BOUNDARY) == 0)) {
                FortranWriteInt(saPAngle->iParmIndex);
            }
        }
        if (iVarArrayElementCount(uUnit->vaAngles))
            saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles);
             i++, saPAngle++) {
            if (((saPAngle->fFlags & PERTURBED) != 0)
                && ((saPAngle->fFlags & BOUNDARY) != 0)) {
                FortranWriteInt(saPAngle->iParmIndex);
            }
        }

        /* Then LAMBDA = 1 */

        if (iVarArrayElementCount(uUnit->vaAngles))
            saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles);
             i++, saPAngle++) {
            if (((saPAngle->fFlags & PERTURBED) != 0)
                && ((saPAngle->fFlags & BOUNDARY) == 0)) {
                FortranWriteInt(saPAngle->iPertParmIndex);
            }
        }
        if (iVarArrayElementCount(uUnit->vaAngles))
            saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles);
             i++, saPAngle++) {
            if (((saPAngle->fFlags & PERTURBED) != 0)
                && ((saPAngle->fFlags & BOUNDARY) != 0)) {
                FortranWriteInt(saPAngle->iPertParmIndex);
            }
        }
        FortranEndLine();

        /* -36E- Torsions that are to be perturbed */
        FortranDebug("-36E-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_DIHEDRAL_ATOMS");
        FortranWriteString("%FORMAT(10I8)");

        FortranFormat(10, INTFORMAT);
        if (iVarArrayElementCount(uUnit->vaTorsions))
            stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions);
             i++, stPTorsion++) {
            if (((stPTorsion->fFlags & PERTURBED) != 0)
                && ((stPTorsion->fFlags & BOUNDARY) == 0)) {

                if ((AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                    (AMBERINDEX(stPTorsion->iAtom4) == 0)) {
                    MESSAGE(
                            ("Had to turn torsion around to avoid K,L == 0\n"));
                    MESSAGE(
                            ("Outer atoms: %s --- %s\n",
                             sContainerName(aA), sContainerName(aD)));
                    SWAP(stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp);
                    SWAP(stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp);
                }
                if (stPTorsion->bProper)
                    iProper = 1;
                else
                    iProper = -1;
                if (stPTorsion->bCalc14)
                    iCalc14 = 1;
                else
                    iCalc14 = -1;
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3) * iCalc14);
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
            }
        }
        if (iVarArrayElementCount(uUnit->vaTorsions))
            stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions);
             i++, stPTorsion++) {
            if (((stPTorsion->fFlags & PERTURBED) != 0)
                && ((stPTorsion->fFlags & BOUNDARY) != 0)) {

                if ((AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                    (AMBERINDEX(stPTorsion->iAtom4) == 0)) {
                    MESSAGE(
                            ("Had to turn torsion around to avoid K,L == 0\n"));
                    MESSAGE(
                            ("Outer atoms: %s --- %s\n",
                             sContainerName(aA), sContainerName(aD)));
                    SWAP(stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp);
                    SWAP(stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp);
                }
                if (stPTorsion->bProper)
                    iProper = 1;
                else
                    iProper = -1;
                if (stPTorsion->bCalc14)
                    iCalc14 = 1;
                else
                    iCalc14 = -1;
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3) * iCalc14);
                FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
            }
        }
        FortranEndLine();

        /* -36F- Index into torsion interaction arrays */
        FortranDebug("-36F-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_DIHEDRAL_PARAMS");
        FortranWriteString("%FORMAT(10I8)");

        FortranFormat(10, INTFORMAT);

        /* First LAMBDA = 0 */

        if (iVarArrayElementCount(uUnit->vaTorsions))
            stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions);
             i++, stPTorsion++) {
            if (((stPTorsion->fFlags & PERTURBED) != 0)
                && ((stPTorsion->fFlags & BOUNDARY) == 0)) {
                FortranWriteInt(stPTorsion->iParmIndex);
            }
        }
        if (iVarArrayElementCount(uUnit->vaTorsions))
            stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions);
             i++, stPTorsion++) {
            if (((stPTorsion->fFlags & PERTURBED) != 0)
                && ((stPTorsion->fFlags & BOUNDARY) != 0)) {
                FortranWriteInt(stPTorsion->iParmIndex);
            }
        }

        /* Then LAMBDA = 1 */

        if (iVarArrayElementCount(uUnit->vaTorsions))
            stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions);
             i++, stPTorsion++) {
            if (((stPTorsion->fFlags & PERTURBED) != 0)
                && ((stPTorsion->fFlags & BOUNDARY) == 0)) {
                FortranWriteInt(stPTorsion->iPertParmIndex);
            }
        }
        if (iVarArrayElementCount(uUnit->vaTorsions))
            stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions);
             i++, stPTorsion++) {
            if (((stPTorsion->fFlags & PERTURBED) != 0)
                && ((stPTorsion->fFlags & BOUNDARY) != 0)) {
                FortranWriteInt(stPTorsion->iPertParmIndex);
            }
        }
        FortranEndLine();

        /* -36G- Residue labels at LAMBDA = 1 */
        /* Just write the labels at LAMBDA = 0 */

        /* Trim the string down to at most 3 characters by */
        /* taking the last three characters if it is too long */
        FortranDebug("-36G-");

        MESSAGE(("Writing the residue labels\n"));
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_RESIDUE_NAME");
        FortranWriteString("%FORMAT(20a4)");
        FortranFormat(20, LBLFORMAT);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
            cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sName;
            if (strlen(cPTemp) > 3)
                cPTemp += (strlen(cPTemp) - 3);
            FortranWriteString(cPTemp);
        }
        FortranEndLine();

        /* -36H- Atom names at LAMBDA = 0 */
        FortranDebug("-36H-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_ATOM_NAME");
        FortranWriteString("%FORMAT(20a4)");

        FortranFormat(20, LBLFORMAT);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
            cPTemp = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sPertName;
            if (strlen(cPTemp) == 0)
                cPTemp = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sName;
            if (strlen(cPTemp) > 4)
                cPTemp += (strlen(cPTemp) - 4);
            FortranWriteString(cPTemp);
        }
        FortranEndLine();

        /* -36I- List of atomic symbols (atom types??????) */
        FortranDebug("-36I-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_ATOM_SYMBOL");
        FortranWriteString("%FORMAT(20a4)");

        FortranFormat(20, LBLFORMAT);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
            cPTemp =
                sAtomPertType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
            if (strlen(cPTemp) == 0)
                cPTemp =
                    sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
            if (strlen(cPTemp) > 3)
                cPTemp += (strlen(cPTemp) - 3);
            FortranWriteString(cPTemp);
        }
        FortranEndLine();

        /* -36J- Value of LAMBDA for each ATOM ????????? */
        /* TODO: Figure out what the hell this is */
        FortranDebug("-36J-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG ALMPER");
        FortranWriteString("%FORMAT(5E16.8)");

        FortranFormat(5, DBLFORMAT);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
            FortranWriteDouble(0.0);
        }
        FortranEndLine();

        /* -36K- Flag to tell whether the atom is perturbed */
        FortranDebug("-36K-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG IAPER");
        FortranWriteString("%FORMAT(10I8)");

        FortranFormat(10, INTFORMAT);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
            if (bAtomPerturbed(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom))
                FortranWriteInt(1);
            else
                FortranWriteInt(0);
        }
        FortranEndLine();

        /* -36L- List of atom types - IACPER */
        FortranDebug("-36L-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_ATOM_TYPE_INDEX");
        FortranWriteString("%FORMAT(10I8)");

        FortranFormat(10, INTFORMAT);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
            iAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iPertTypeIndex - 1;
            iTemp = *PVAI(vaNBIndex, int, iAtom);
            FortranWriteInt(iTemp + 1);
        }
        FortranEndLine();

        /* -36M- Perturbed charges */
        FortranDebug("-36M-");
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG PERT_CHARGE");
        FortranWriteString("%FORMAT(5E16.8)");

        FortranFormat(5, DBLFORMAT);
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
            aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
                        if( GDefaults.iGibbs ){
                                if( bAtomPerturbed(aAtom))
                                        FortranWriteDouble(ELECTRONTOKCAL *
                                   dAtomPertCharge(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom));
                                else
                                        FortranWriteDouble(ELECTRONTOKCAL *
                                   dAtomCharge(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom));
                        } else {
                                FortranWriteDouble(ELECTRONTOKCAL *
                                   (dAtomCharge(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom) +
                                   dAtomPertCharge(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom)));
                        }
        }
        FortranEndLine();

    }

    /*
     *  polarizabilities
     */
    if (bPolar) {
        iCount = 0;
        iCountPerturbed = 0;
        iMax = iVarArrayElementCount(uUnit->vaAtoms);
        MESSAGE(("Writing the atomic polarizabilities\n"));
        FortranFormat(1, "%-80s");
        FortranWriteString("%FLAG POLARIZABILITY");
        FortranWriteString("%FORMAT(5E16.8)");
        FortranFormat(5, DBLFORMAT);
        saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
        for (i = 0; i < iMax; i++, saPAtom++) {
            iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
            ParmSetAtom(uUnit->psParameters, iIndex, sType,
                        &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                        &dRStar14, &dScreenF, &iElement, &iHybridization,
                        sDesc);
            if (dPolar == -1.0) {
                dPolar = 0.0;
                iCount++;
            }
            FortranWriteDouble(dPolar);
        }
        if (iCount > 0)
            VP0(("Total atoms with default polarization=0.0: %d of %d\n",
                 iCount, iMax));
        FortranEndLine();
        
        
            FortranFormat(1, "%-80s");
            FortranWriteString("%FLAG DIPOLE_DAMP_FACTOR");
            FortranWriteString("%FORMAT(5E16.8)");
            FortranFormat(5, DBLFORMAT);
            saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
            for (i = 0; i < iMax; i++, saPAtom++) {
                iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
                ParmSetAtom(uUnit->psParameters, iIndex, sType,
                            &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                            &dRStar14, &dScreenF, &iElement, &iHybridization,
                            sDesc);
                if (dScreenF == 0.0) {
                    dScreenF = GDefaults.dDipoleDampFactor;
                }
                FortranWriteDouble(dScreenF);
            }                       
        FortranEndLine();
        
        if (bPert) {
            int iPertTot = 0;
            FortranFormat(1, "%-80s");
            FortranWriteString("%FLAG PERT_POLARIZABILITY");
            FortranWriteString("%FORMAT(5E16.8)");
            FortranFormat(5, DBLFORMAT);
            saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
            for (i = 0; i < iMax; i++, saPAtom++) {
                BOOL bTmp = bAtomPerturbed(saPAtom->aAtom);
                
                if (bTmp) {
                    iIndex = iParmSetFindAtom(uUnit->psParameters,
                                              saPAtom->sPertType);
                    iPertTot++;
                } else
                    iIndex = iParmSetFindAtom(uUnit->psParameters,
                                              saPAtom->sType);
                ParmSetAtom(uUnit->psParameters, iIndex, sType,
                            &dMass, &dPolar, &dEpsilon, &dRStar,
                            &dEpsilon14, &dRStar14, &dScreenF, &iElement,
                            &iHybridization, sDesc);
                if (dPolar == -1.0) {
                    dPolar = 0.0;
                    if (bTmp)
                        iCountPerturbed++;
                }
                FortranWriteDouble(dPolar);
            }
            FortranEndLine();
            if (iCountPerturbed > 0)
                VP0(
                    ("Total pert atoms with default polarization=0.0: %d of %d\n",
                     iCountPerturbed, iPertTot));
        }
    }

    /*  Charmm-style parameters  */

    if (GDefaults.iCharmm) {
        /* -19- Lennard jones r**12 term for all 14 interactions */
        /* CN114 array */
        FortranDebug("-19-");

        FortranFormat(5, DBLFORMAT);
        for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
            FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dA14);
        }
        FortranEndLine();

        /* -20- Lennard jones r**6 term for all 14 interactions */
        /* CN214 array */
        FortranDebug("-20-");

        FortranFormat(5, DBLFORMAT);
        for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
            FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dC14);
        }
        FortranEndLine();

        /* -13- Force constants for Urey-Bradley */
        FortranDebug("-13-");

        FortranFormat(5, DBLFORMAT);
        for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
            ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                         &dKt, &dT0, &dTkub, &dRkub, sDesc);
            FortranWriteDouble(dTkub);
        }
        FortranEndLine();

        /* -14- Equilibrium distances for Urey-Bradley */
        FortranDebug("-14-");

        FortranFormat(5, DBLFORMAT);
        for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
            ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                         &dKt, &dT0, &dTkub, &dRkub, sDesc);
            FortranWriteDouble(dRkub);
        }
        FortranEndLine();

    }

    // 
    // CMAP parameters, Mengjuei Hsieh and Yong Duan
    //

    SaveAmberParmCMAP(uUnit, fOut);

        /********************************************************/
        /* Write the coordinate file                            */
        /********************************************************/
    if ( bNetcdf == TRUE ) {
      zUnitIOSaveAmberNetcdf(uUnit, crdName);
      // Clean up arrays and return, netcdf file will be written later
      VarArrayDestroy(&vaNBIndexMatrix);
      VarArrayDestroy(&vaNBParameters);
      VarArrayDestroy(&vaExcludedAtoms);
      VarArrayDestroy(&vaExcludedCount);
      VarArrayDestroy(&vaNBIndex);
      VarArrayDestroy(&vaNonBonds);
      return;
    }
    FortranFile(fCrd);

    FortranFormat(1, "%s");
    FortranWriteString(sContainerName(uUnit));
    FortranEndLine();

    FortranFormat(1, "%6d");
    FortranWriteInt(iVarArrayElementCount(uUnit->vaAtoms));
    FortranEndLine();

    FortranFormat(6, "%12.7lf");
    if (bUnitUseBox(uUnit)) {
        double dX2, dY2, dZ2;

        if( GDefaults.nocenter == 0 ){
          UnitGetBox(uUnit, &dX, &dY, &dZ);
          dX2 = dX * 0.5;
          dY2 = dY * 0.5;
          dZ2 = dZ * 0.5;
        } else {
          dX2 = 0.0;
          dY2 = 0.0;
          dZ2 = 0.0;
        }

        /*
        *  shift box to Amber spot; later, add a cmd opt or environment
        *      var to switch between 0,0,0 center (spasms) or corner
        */
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
            vPos = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->vPos;
            FortranWriteDouble(dVX(&vPos) + dX2);
            FortranWriteDouble(dVY(&vPos) + dY2);
            FortranWriteDouble(dVZ(&vPos) + dZ2);
        }
        FortranEndLine();
        FortranWriteDouble(dX);
        FortranWriteDouble(dY);
        FortranWriteDouble(dZ);
        FortranWriteDouble(dUnitBeta(uUnit) / DEGTORAD);
        FortranWriteDouble(dUnitBeta(uUnit) / DEGTORAD);
        FortranWriteDouble(dUnitBeta(uUnit) / DEGTORAD);
        FortranEndLine();
    } else {
        for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
            vPos = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->vPos;
            FortranWriteDouble(dVX(&vPos));
            FortranWriteDouble(dVY(&vPos));
            FortranWriteDouble(dVZ(&vPos));
        }
        FortranEndLine();
    }

    VarArrayDestroy(&vaNBIndexMatrix);
    VarArrayDestroy(&vaNBParameters);
    VarArrayDestroy(&vaExcludedAtoms);
    VarArrayDestroy(&vaExcludedCount);
    VarArrayDestroy(&vaNBIndex);
    VarArrayDestroy(&vaNonBonds);
    fclose( fCrd );

}


/*
 *      zUnitIOSaveAmberParmFormat_old
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save the UNIT in the AMBER PARM file format.
 *      This requires that the UNIT tables be built and that
 *      the UNIT contain a parameter set.
 *      The iContainerTempInt(atom) should still return the index
 *      of the atom within the vaAtoms array.
 *      Atom coordinates are written to the file (fCrd).
 *
 *NOTE: This routine depends on the order of the RESIDUEs in
 *      vaResidues being such that solvent residues follow
 *      all other RESIDUEs.  I know that this is going
 *      against the philosophy that the data written to
 *      OFF files has NO implicit order, and outside of
 *      this program that is how they should be handled.
 *      But it was SO convenient to sort the RESIDUEs
 *      as they are put into the table that I could
 *      not resist.
 *
 *TODO: Add RESTRAINT code
 *TODO: Add CAP information
 */
#undef INTFORMAT
//#define AMBERINDEX(i)   3*(i-1)
#define INTFORMAT       "%6d"
//#define DBLFORMAT       "%16.8lE"
//#define LBLFORMAT       "%-4s"
//#define ELECTRONTOKCAL  18.2223

void
zUnitIOSaveAmberParmFormat_old( UNIT uUnit, FILE *fOut, char *crdName, 
        BOOL bPolar, BOOL bPert )
{
int             i, iMax, iIndex;
LOOP            lTemp, lSpan;
ATOM            aAtom, aA, aB, aC, aD;
int             iCount, iBondWith, iBondWithout;
int             iAngleWith, iAngleWithout;
int             iTorsionWith, iTorsionWithout;
VARARRAY        vaExcludedAtoms, vaExcludedCount;
VARARRAY        vaNBIndexMatrix, vaNBParameters;
VARARRAY        vaNBIndex, vaNonBonds;
int             iCountPerturbed, iCountBondPerturbed, iCountBondBoundary;
int             iCountAnglePerturbed, iCountAngleBoundary;
int             iCountTorsionPerturbed, iCountTorsionBoundary;
int             iNumExtra;
SAVEBONDt       *sbPBond;
SAVEANGLEt      *saPAngle;
SAVEATOMt       *saPAtom;
SAVETORSIONt    *stPTorsion;
SAVERESTRAINTt  *srPRestraint;
double          dMass, dPolar, dR, dKb, dR0, dKt, dT0, dTkub, dRkub, dKp, dP0, dC, dD, dTemp;
double		dScEE, dScNB, dScreenF;
STRING          sAtom1, sAtom2, sAtom3, sAtom4, sType1, sType2;
int             iN, iAtoms, iMaxAtoms, iTemp, iAtom, iCalc14, iProper;
int             iElement, iHybridization, iStart, iFirstSolvent;
RESIDUE         rRes;
BOOL            bFoundSome;
VECTOR          vPos;
char            *cPTemp;
double          dX, dY, dZ, dEpsilon, dRStar, dEpsilon14, dRStar14;
STRING          sDesc, sType;
VARARRAY        vaMolecules;
IX_REC          *ePResEnt;
IX_DESC         iResIx;

  // Open the coordinate file
  FILE *fCrd = FOPENCOMPLAIN( crdName, "w" );
  if ( fOut == NULL ) {
    VP0(( "Could not open file: %s\n", crdName ));
  }

                /* Build the excluded atom list */

    MESSAGE(( "Building the excluded atom list\n" ));
    vaExcludedCount = vaVarArrayCreate( sizeof(int) );
    vaExcludedAtoms = vaVarArrayCreate( sizeof(int) );

    iCountPerturbed = 0;
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        aAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->aAtom;

        if ( bAtomFlagsSet( aAtom, ATOMPERTURB ) )
                iCountPerturbed++;

        lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
        iCount = 0;
        bFoundSome = FALSE;
        iStart = iVarArrayElementCount( vaExcludedAtoms );
        while ( (aA = (ATOM)oNext(&lSpan)) ) {

            if ( aA == aAtom ) continue;
            
                /* If the atom is more than three away from the first atom */
                /* then it is not in the excluded atom list */

            if ( iAtomBackCount(aA) >= 4 ) break;

            if ( iContainerTempInt(aA) > iContainerTempInt(aAtom) ) { 
                VarArrayAdd( vaExcludedAtoms, (GENP)&iContainerTempInt(aA) );
                bFoundSome = TRUE;
                iCount++;
            }
        }
        if ( !bFoundSome ) {
            iAtoms = 0;
            VarArrayAdd( vaExcludedAtoms, (GENP)&iAtoms );
            iCount++;
        } else {

                /* Sort the part of the VARARRAY just added so that */
                /* the excluded ATOMs are in ascending order by index */

            SortByInteger( (GENP) PVAI( vaExcludedAtoms, int, iStart ),
                                iCount,
                                sizeof(int),
                               (GENP) PVAI( vaExcludedAtoms, int, iStart ),
                                TRUE );
        } 

        VarArrayAdd( vaExcludedCount, (GENP)&iCount );
    }

    /*
     *  mark main chain atoms where possible, noting the 
     *  number of atoms in the largest residue. keep
     *  track of residues which can't be marked.
     */
    VP0(( "Not Marking per-residue atom chain types.\n" ));
    iMaxAtoms = 0;
    create_index( &iResIx, 2, 0 );
    MALLOC( ePResEnt, IX_REC *, sizeof(IX_REC)+8 );
    ePResEnt->recptr = NULL;    /* for Purify */

    VP0(( "Marking per-residue atom chain types.\n" ));

    iMaxAtoms = 0;
    lTemp = lLoop( (OBJEKT)uUnit, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lTemp)) ) {
        int     iAtoms = MarkMainChainAtoms( rRes, 0 );
        if ( iAtoms > 0 )
                (void)MarkSideChains( rRes );
        if ( iAtoms < 0 ) {
                iAtoms = -iAtoms;
                /*
                 *  couldn't mark main chains
                 */

                strcpy( ePResEnt->key, rRes->cHeader.sName );
                if ( add_key( ePResEnt, &iResIx ) != IX_OK )
                        DFATAL(( "add_key() residue chain\n" ));
        }
        if ( iAtoms > iMaxAtoms ) 
                iMaxAtoms = iAtoms;
    }
    /*
     *  print warnings
     */
    first_key( &iResIx );
    i = 1;
    while ( next_key( ePResEnt, &iResIx ) == IX_OK ) {
        if ( i ) {
                VP0(( "  (Residues lacking connect0/connect1 - \n" ));
                VP0(( "   these don't have chain types marked:\n\n" ));
                VP0(( "\tres\ttotal affected\n\n" ));
                i = 0;
        }
        VP0(( "\t%s\t%d\n", ePResEnt->key, ePResEnt->count));
    }
    if (!i)
        VP0(( "  )\n" ));
    destroy_index( &iResIx );
    FREE( ePResEnt );

                /* Build the NON-BOND arrays that AMBER needs */

    zUnitIOBuildNonBondArrays( uUnit, &vaNBIndexMatrix, &vaNBParameters,
                                        &vaNBIndex, &vaNonBonds );
 

    FortranFile( fOut );

#if 0
        /*
         *---------------------------------------------------------
         *
         *      Turn on debugging of fortran format output file
         *      by sticking comments into the file.
         */

    FortranDebugOn();
#endif

    
        /* -1- Save the title of the UNIT */
    FortranDebug( "-1-" );
    MESSAGE(( "Saving the name of the UNIT\n" ));
    FortranFormat( 1, "%-80s" );
    FortranWriteString( sContainerName(uUnit) );

        /* -2- Save control information */
    FortranDebug( "-2-" );
    MESSAGE(( "Saving all the main control variables\n" ));
    FortranFormat( 12, INTFORMAT );

/*NTOTAT*/      
    FortranWriteInt( iVarArrayElementCount( uUnit->vaAtoms ) );
/*NTYPES*/      
    FortranWriteInt( iVarArrayElementCount( vaNonBonds ) );
    
        /* Count the number of bonds with hydrogens, and without */

    iBondWith = 0;
    iBondWithout = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaBonds ); i++ ) {
        sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1 )->aAtom;
        if ( bPERT_BOND(bPert,aA,aD) ) 
                continue;
        if ( iAtomElement(aA) == HYDROGEN ||
                iAtomElement(aD) == HYDROGEN ) iBondWith++;
        else iBondWithout++;
    }
/*NBONH*/       
    FortranWriteInt( iBondWith );
/*NBONA*/       
    FortranWriteInt( iBondWithout );

        /* Count the number of angles with hydrogens, and without */

    iAngleWith = 0;
    iAngleWithout = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaAngles ); i++ ) {
        saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1 )->aAtom;
        if ( bPERT_ANGLE(bPert,aA,aB,aD) ) 
                continue;
        if ( iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aD) == HYDROGEN )
            iAngleWith++;
        else
            iAngleWithout++;
    }
/*NTHETH*/      
    FortranWriteInt( iAngleWith );
/*NTHETA*/      
    FortranWriteInt( iAngleWithout );

        /* Count the number of torsions with hydrogens, and without */

    iTorsionWith = 0;
    iTorsionWithout = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaTorsions ); i++ ) {
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2-1 )->aAtom;
        aC = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4-1 )->aAtom;
        if ( bPERT_TORSION(bPert,aA,aB,aC,aD) ) 
            continue;
        if ( iAtomElement(aA) == HYDROGEN || 
                iAtomElement(aB) == HYDROGEN ||
                iAtomElement(aC) == HYDROGEN ||
                iAtomElement(aD) == HYDROGEN ) 
            iTorsionWith++;
        else 
            iTorsionWithout++;
    }
/*NPHIH*/       
    FortranWriteInt( iTorsionWith );
/*NPHIA*/       
    FortranWriteInt( iTorsionWithout );

/*JHPARM*/      
    FortranWriteInt( 0 );
/*JPARM*/       
    FortranWriteInt( 0 );

        /* Write the number of excluded atoms */

/*NEXT*/        
    FortranWriteInt( iVarArrayElementCount(vaExcludedAtoms) );
/*NTOTRS*/      
    FortranWriteInt( iVarArrayElementCount(uUnit->vaResidues) );

        /* Write the number of bonds/angles/torsions without hydrogens */
        /* PLUS the number of RESTRAINT bonds/angles/torsions */

/*MBONA*/       
    FortranWriteInt( iBondWithout+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTBOND ) );
/*MTHETA*/      
    FortranWriteInt( iAngleWithout+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTANGLE ) );
/*MPHIA*/       
    FortranWriteInt( iTorsionWithout+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTTORSION ) );

        /* Write the number of unique bond types, angle types, torsion types */
        /* Add in the number of RESTRAINT bonds/angles/torsion because */
        /* they will have new parameters */

/*MUMBND*/      
    FortranWriteInt( iParmSetTotalBondParms(uUnit->psParameters)+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTBOND ) );
/*MUMANG*/      
    FortranWriteInt( iParmSetTotalAngleParms(uUnit->psParameters)+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTANGLE ) );
/*NPTRA*/       
    FortranWriteInt( iParmSetTotalTorsionParms(uUnit->psParameters) +
                     iParmSetTotalImproperParms(uUnit->psParameters) +
                     iUnitRestraintTypeCount( uUnit, RESTRAINTTORSION ) );

                /* TODO - have different arrays for different restraint types*/
    if ( iVarArrayElementCount( uUnit->vaRestraints ) )
        VP0(( " Restraints:  Bond %d  Angle %d  Torsion %d\n",
                iUnitRestraintTypeCount( uUnit, RESTRAINTBOND ),
                iUnitRestraintTypeCount( uUnit, RESTRAINTANGLE ),
                iUnitRestraintTypeCount( uUnit, RESTRAINTTORSION ) ));
    else
        VP0(( " (no restraints)\n" ));

        /* The next parameter corresponds to NATYP in AMBER */
        /* I don't know what it does, and Dave Spellmeyer says that */
        /* he only uses it to skip over the SOLTY array */

/*NATYP*/       
    FortranWriteInt( iParmSetTotalAtomParms(uUnit->psParameters) );
/*NHB*/         
    FortranWriteInt( iParmSetTotalHBondParms(uUnit->psParameters) );
/*IFPERT*/
    if ( bPert )
        FortranWriteInt( 1 );
    else
        FortranWriteInt( 0 );



        /* Count the number of bonds to be perturbed, and those across the */
        /* perturbation/non-perturbed boundary */

    iCountBondPerturbed = 0;
    iCountBondBoundary = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaBonds ); i++ ) {
        sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, i );
        if ( (sbPBond->fFlags&PERTURBED) != 0 ) {
            iCountBondPerturbed++;
            if ( (sbPBond->fFlags&BOUNDARY) != 0 ) {
                MESSAGE(( "Boundary pert bond %d-%d\n", 
                                sbPBond->iAtom1, sbPBond->iAtom2 ));
                iCountBondBoundary++;
            }
        }
    }

    MESSAGE(( "Perturbed bonds: %d\n", iCountBondPerturbed ));
    MESSAGE(( "Perturbed boundary bonds: %d\n", iCountBondBoundary ));

        /* Count the number of angles to be perturbed, and those on the */
        /* boundary */

    iCountAnglePerturbed = 0;
    iCountAngleBoundary = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaAngles ); i++ ) {
        saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, i );
        if ( (saPAngle->fFlags&PERTURBED) != 0 ) iCountAnglePerturbed++;
        if ( (saPAngle->fFlags&BOUNDARY) != 0 ) iCountAngleBoundary++;
    }

        /* Count the number of torsions and impropers to be perturbed */
        /* and those on the boundary */

    iCountTorsionPerturbed = 0;
    iCountTorsionBoundary = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaTorsions ); i++ ) {
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, i );
        if ( (stPTorsion->fFlags&PERTURBED) != 0 ) iCountTorsionPerturbed++;
        if ( (stPTorsion->fFlags&BOUNDARY) != 0 ) iCountTorsionBoundary++;
    }

    /*NBPER*/   
    FortranWriteInt( iCountBondPerturbed );
    /*NGPER*/   
    FortranWriteInt( iCountAnglePerturbed );
    /*NDPER*/   
    FortranWriteInt( iCountTorsionPerturbed );
    /*MBPER*/   
    FortranWriteInt( iCountBondPerturbed-iCountBondBoundary );
    /*MGPER*/   
    FortranWriteInt( iCountAnglePerturbed-iCountAngleBoundary );
    /*MDPER*/   
    FortranWriteInt( iCountTorsionPerturbed-iCountTorsionBoundary );

        /* Save flag for periodic boundary conditions */

    /*IFBOX*/
    if ( bUnitUseBox(uUnit) ) {
        if ( bUnitBoxOct(uUnit) )
                FortranWriteInt( 2 );
        else
                FortranWriteInt( 1 );
    } else
        FortranWriteInt( 0 );

        /* Save the number of atoms in the largest residue */

    /*NMXRS*/
    printf("iMaxAoms (2) %i\n",iMaxAtoms);
    FortranWriteInt( iMaxAtoms );

        /* Save flag for cap information */

    /*IFCAP*/
    if ( bUnitUseSolventCap(uUnit) )
        FortranWriteInt( 1 );
    else
        FortranWriteInt( 0 );

    /*NUMEXTRA*/
    iNumExtra = 0;
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
                cPTemp = sAtomType( PVAI(uUnit->vaAtoms, SAVEATOMt, i )->aAtom );
                if( !strncmp( cPTemp, "EP", 2 )) iNumExtra++;
    }
    FortranWriteInt( iNumExtra );

    FortranEndLine();


        /* -3-  write out the names of the atoms */
    FortranDebug( "-3-" );

    MESSAGE(( "Writing the names of the atoms\n" ));
    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        cPTemp = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->sName;
        if ( strlen( cPTemp ) > 4 ) cPTemp += ( strlen(cPTemp)-4 );
        FortranWriteString( cPTemp );
    }
    FortranEndLine();
    
        /* -4- write out the atomic charges */
    FortranDebug( "-4-" );
        
    MESSAGE(( "Writing the atomic charges\n" ));
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteDouble( PVAI( uUnit->vaAtoms, SAVEATOMt, i )->dCharge *
                                ELECTRONTOKCAL );
    }
    FortranEndLine();

        /* -5- write out the atomic masses */
    FortranDebug( "-5-" );

    MESSAGE(( "Writing the atomic masses\n" ));         
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        saPAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i );
        iIndex = iParmSetFindAtom( uUnit->psParameters, saPAtom->sType );
        ParmSetAtom( uUnit->psParameters, iIndex, sType,
                        &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, 
			&dScreenF,
            &iElement, &iHybridization, sDesc );
        FortranWriteDouble( dMass );
    }
    FortranEndLine();

        /* -6- write out the atomic types */
    FortranDebug( "-6-" );

    MESSAGE(( "Writing the atomic types\n" ));          
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        iAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->iTypeIndex-1;
        iTemp = *PVAI( vaNBIndex, int, iAtom );
        FortranWriteInt( iTemp+1 );
    }
    FortranEndLine();

        /* -7- write out the starting index into the excluded atom list */
    FortranDebug( "-7-" );

    MESSAGE(( "Writing the starting index into the excluded atom list\n" ));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteInt( *PVAI( vaExcludedCount, int, i ) );
    }
    FortranEndLine();

        /* -8- Write the index for the position of the non bond type */
                /* of each type */
    FortranDebug( "-8-" );

    MESSAGE(( "writing position of the non bond type of each type\n"));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(vaNBIndexMatrix); i++ ) {
        FortranWriteInt( *PVAI( vaNBIndexMatrix, int, i ) );
    }
    FortranEndLine();

        /* -9- Residue labels */
                /* Trim the string down to at most 3 characters by */
                /* taking the last three characters if it is too long */
    FortranDebug( "-9-" );

    MESSAGE(( "Writing the residue labels\n" ));
    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaResidues); i++ ) {
        cPTemp = PVAI( uUnit->vaResidues, SAVERESIDUEt, i )->sName;
        if ( strlen( cPTemp ) > 3 ) cPTemp += ( strlen(cPTemp)-3 );
        FortranWriteString( cPTemp );
    }
    FortranEndLine();

        /* -10- Pointer list for all the residues */
    FortranDebug( "-10-" );

    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaResidues); i++ ) {
        FortranWriteInt( PVAI( uUnit->vaResidues, 
                                SAVERESIDUEt, i )->iAtomStartIndex );
    }
    FortranEndLine();

        /* -11- Force constants for bonds */
    FortranDebug( "-11-" );

    MESSAGE(( "Writing bond force constants\n" ));
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalBondParms(uUnit->psParameters); i++ ) {
       ParmSetBond( uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, sDesc );
       FortranWriteDouble( dKb );
    }
                /* Write the RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTBOND, dKx, i+1 );
    FortranEndLine();

        /* -12- Equilibrium bond lengths */
    FortranDebug( "-12-" );

    MESSAGE(( "Writing equilibrium bond lengths\n" ));
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalBondParms(uUnit->psParameters); i++ ) {
        ParmSetBond( uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, sDesc );
        FortranWriteDouble( dR0 );
    }
                /* Write the bond RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTBOND, dX0, i+1 );
    FortranEndLine();

        /* -13- Force constants for angles */
    FortranDebug( "-13-" );

    MESSAGE(( "Writing angle force constants\n" ));
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalAngleParms(uUnit->psParameters); i++ ) {
        ParmSetAngle( uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                        &dKt, &dT0, &dTkub, &dRkub, sDesc );
        FortranWriteDouble( dKt );
    }
                /* Write the angle RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTANGLE, dKx, i+1 );
    FortranEndLine();

        /* -14- Equilibrium angle values */
    FortranDebug( "-14-" );

    MESSAGE(( "Writing equilibrium angle values\n" ));
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalAngleParms(uUnit->psParameters); i++ ) {
        ParmSetAngle( uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                        &dKt, &dT0, &dTkub, &dRkub, sDesc );
        FortranWriteDouble( dT0 );
    }
                /* Write the angle RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTANGLE, dX0, i+1 );
    FortranEndLine();

        /* -15- Force constants for torsions */
    FortranDebug( "-15-" );

    MESSAGE(( "Writing torsional force constants\n" ));
    FortranFormat( 5, DBLFORMAT );
    MESSAGE(( "There are %d torsions and %d impropers\n", 
                iParmSetTotalTorsionParms(uUnit->psParameters),
                iParmSetTotalImproperParms(uUnit->psParameters) ));
    for ( i=0; i<iParmSetTotalTorsionParms(uUnit->psParameters); i++ ) {
        ParmSetTorsion( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc );
        MESSAGE(( "Torsion %d  %s-%s-%s-%s %d %lf %lf\n",
                        i, sAtom1, sAtom2, sAtom3, sAtom4,
                        iN, dKp, dP0 ));
        FortranWriteDouble( dKp );
    }
    for ( i=0; i<iParmSetTotalImproperParms(uUnit->psParameters); i++ ) {
        ParmSetImproper( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, sDesc );
        MESSAGE(( "Improper %d  %s-%s-%s-%s %d %lf %lf\n",
                        i, sAtom1, sAtom2, sAtom3, sAtom4,
                        iN, dKp, dP0 ));
        FortranWriteDouble( dKp );
    }
                /* Write the torsion RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTTORSION, dKx, i+1 );
    FortranEndLine();

        /* -16- Division factor for the dihedral angles */
    FortranDebug( "-16-" );

    MESSAGE(( "Writing multiplicity of torsion interaction\n" ));
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalTorsionParms(uUnit->psParameters); i++ ) {
        ParmSetTorsion( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc );
        dTemp = iN;
        FortranWriteDouble( dTemp );
    }
    for ( i=0; i<iParmSetTotalImproperParms(uUnit->psParameters); i++ ) {
        ParmSetImproper( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, sDesc );
        dTemp = iN;
        FortranWriteDouble( dTemp );
    }
                /* Write the torsion RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTTORSION, dX0, i+1 );
    FortranEndLine();

        /* -17- Phase for torsions */
    FortranDebug( "-17-" );

    MESSAGE(( "Writing phase for torsion interactions\n" ));
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalTorsionParms(uUnit->psParameters); i++ ) {
        ParmSetTorsion( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc );
        FortranWriteDouble( dP0 );
    }
    for ( i=0; i<iParmSetTotalImproperParms(uUnit->psParameters); i++ ) {
        ParmSetImproper( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, sDesc );
        FortranWriteDouble( dP0 );
    }
                /* Write the torsion RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTTORSION, dN, i+1 );
    FortranEndLine();

        /* -18- Not used, reserved for future use, uses NATYP */
                /* Corresponds to the AMBER SOLTY array */
    FortranDebug( "-18-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalAtomParms(uUnit->psParameters); i++ ) {
        FortranWriteDouble( 0.0 );
    }
    FortranEndLine();

        /* -19- Lennard jones r**12 term for all possible interactions */
                /* CN1 array */
    FortranDebug( "-19-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(vaNBParameters); i++ ) {
        FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dA );
    }
    FortranEndLine();
    
        /* -20- Lennard jones r**6 term for all possible interactions */
                /* CN2 array */
    FortranDebug( "-20-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(vaNBParameters); i++ ) {
        FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dC );
    }
    FortranEndLine();
 
        /* -21- Write the bond interactions that include hydrogen */
                /* Write the two indices into the atom table, then the index */
                /* into the interaction table */
    FortranDebug( "-21-" );

    MESSAGE(( "Writing the bond interactions with hydrogens\n" ));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaBonds ); i++ ) {
        sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1 )->aAtom;
        if ( bPERT_BOND(bPert,aA,aD) ) 
            continue;
        if ( iAtomElement(aA) == HYDROGEN ||
                iAtomElement(aD) == HYDROGEN ) {
            FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
            FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
            FortranWriteInt( sbPBond->iParmIndex );
        }
    }
    FortranEndLine();

        /* -22- Write the bond interactions that dont include hydrogen */
                /* Write the two indices into the atom table, then the index */
                /* into the interaction table */
    FortranDebug( "-22-" );

    MESSAGE(( "Writing the bond interactions without hydrogens\n" ));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaBonds ); i++ ) {
        sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1 )->aAtom;
        if ( bPERT_BOND(bPert,aA,aD) ) 
            continue;
        if ( !(iAtomElement(aA) == HYDROGEN ||
                iAtomElement(aD) == HYDROGEN) ) {
            FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
            FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
            FortranWriteInt( sbPBond->iParmIndex );
        }
    }
        /* Write out the (bond without H) RESTRAINT interactions */
        /* The iParmIndex field is set in RESTRAINTLOOP */
    if ( (iMax = iVarArrayElementCount( uUnit->vaRestraints )) ) {
        srPRestraint = PVAI( uUnit->vaRestraints, SAVERESTRAINTt, 0 );
        for ( i=0; i<iMax; i++, srPRestraint++ ) {
                if ( srPRestraint->iType == RESTRAINTBOND ) {
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom1) );
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom2) );
                        FortranWriteInt( srPRestraint->iParmIndex );
                }
        }
    }
    FortranEndLine();

        /* -23- Write the angle interactions that include hydrogen */
                /* Write the three indices into the atom table, then the index*/
                /* into the interaction table */
    FortranDebug( "-23-" );

    MESSAGE(( "Writing the angle interactions with hydrogens\n" ));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaAngles ); i++ ) {
        saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1 )->aAtom;
        if ( bPERT_ANGLE(bPert,aA,aB,aD) ) 
            continue;
        if ( iAtomElement(aA) == HYDROGEN
                || iAtomElement(aB) == HYDROGEN
                || iAtomElement(aD) == HYDROGEN ) {
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
            FortranWriteInt( saPAngle->iParmIndex );
        }
    }
    FortranEndLine();

        /* -24- Write the angle interactions that dont include hydrogen */
                /* Write the three indices into the atom table, then the index*/
                /* into the interaction table */
    FortranDebug( "-24-" );

    MESSAGE(( "Writing the angle interactions without hydrogens\n" ));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaAngles ); i++ ) {
        saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1 )->aAtom;
        if ( bPERT_ANGLE(bPert,aA,aB,aD) ) 
            continue;
        if ( !(iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aD) == HYDROGEN) ) {
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
            FortranWriteInt( saPAngle->iParmIndex );
        }
    }
        /* Write out the RESTRAINT interactions */
        /* The iParmIndex field is set in RESTRAINTLOOP */
    if ( (iMax = iVarArrayElementCount( uUnit->vaRestraints )) ) {
        srPRestraint = PVAI( uUnit->vaRestraints, SAVERESTRAINTt, 0 );
        for ( i=0; i<iMax; i++, srPRestraint++ ) {
                if ( srPRestraint->iType == RESTRAINTANGLE ) {
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom1) );
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom2) );
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom3) );
                        FortranWriteInt( srPRestraint->iParmIndex );
                }
        }
    }
    FortranEndLine();

        /* -25- Write the torsion interactions that include hydrogen */
               /* Write the three indices into the atom table, then the index*/
               /* into the interaction table */
    FortranDebug( "-25-" );

    MESSAGE(( "Writing the torsion interactions with hydrogens\n" ));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaTorsions ); i++ ) {
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2-1 )->aAtom;
        aC = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4-1 )->aAtom;
        if ( bPERT_TORSION(bPert,aA,aB,aC,aD) ) 
            continue;
        if ( iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aC) == HYDROGEN
                || iAtomElement(aD) == HYDROGEN ) {
            if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                 (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                MESSAGE(( "Had to turn torsion around to avoid K,L == 0\n" ));
                MESSAGE(( "Outer atoms: %s --- %s\n", 
                                sContainerName(aA), sContainerName(aD) ));
                MESSAGE(( "Old order %d %d %d %d\n",
                                stPTorsion->iAtom1,
                                stPTorsion->iAtom2,
                                stPTorsion->iAtom3,
                                stPTorsion->iAtom4 ));
                SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
                MESSAGE(( "New order %d %d %d %d\n",
                                stPTorsion->iAtom1,
                                stPTorsion->iAtom2,
                                stPTorsion->iAtom3,
                                stPTorsion->iAtom4 ));
            }
            if ( stPTorsion->bProper )  iProper = 1;
            else                        iProper = -1;
            if ( stPTorsion->bCalc14 )  iCalc14 = 1;
            else                        iCalc14 = -1;
            if( GDefaults.iCharmm && iProper == -1 ){
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1)*iCalc14 );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            } else {
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3)*iCalc14 );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            }
            FortranWriteInt( stPTorsion->iParmIndex );
        }
    }
    FortranEndLine();

        /* -26- Write the torsion interactions that dont include hydrogen */
                /* Write the three indices into the atom table, then the index*/
                /* into the interaction table */
    FortranDebug( "-26-" );

    MESSAGE(( "Writing the torsion interactions without hydrogens\n" ));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaTorsions ); i++ ) {
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2-1 )->aAtom;
        aC = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4-1 )->aAtom;
        if ( bPERT_TORSION(bPert,aA,aB,aC,aD) ) 
            continue;
        if ( !(iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aC) == HYDROGEN
                || iAtomElement(aD) == HYDROGEN) ) {
            if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                 (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                MESSAGE(( "Had to turn torsion to avoid K,L == 0\n" ));
                MESSAGE(( "Outer atoms: %s --- %s\n", 
                                sContainerName(aA), sContainerName(aD) ));
                SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
            }
            if ( stPTorsion->bCalc14 )  iCalc14 = 1;
            else                        iCalc14 = -1;
            if ( stPTorsion->bProper )  iProper = 1;
            else                        iProper = -1;
            if( GDefaults.iCharmm && iProper == -1 ){
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1)*iCalc14 );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            } else {
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3)*iCalc14 );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            }
            FortranWriteInt( stPTorsion->iParmIndex );
        }
    }
        /* Write out the RESTRAINT interactions */
        /* The iParmIndex field is set in RESTRAINTLOOP */
    if ( (iMax = iVarArrayElementCount( uUnit->vaRestraints )) ) {
        srPRestraint = PVAI( uUnit->vaRestraints, SAVERESTRAINTt, 0 );
        for ( i=0; i<iMax; i++, srPRestraint++ ) {
            if ( srPRestraint->iType == RESTRAINTTORSION ) {
                if ( (AMBERINDEX(srPRestraint->iAtom3) == 0 ) ||
                     (AMBERINDEX(srPRestraint->iAtom4) == 0 ) ) {
                    MESSAGE(( "Had to turn RESTRAINT torsion around to avoid\n" ));
                    MESSAGE(( "K,L == 0\n" ));
                    SWAP( srPRestraint->iAtom1, srPRestraint->iAtom4, iTemp );
                    SWAP( srPRestraint->iAtom2, srPRestraint->iAtom3, iTemp );
                }
                FortranWriteInt( AMBERINDEX(srPRestraint->iAtom1) );
                FortranWriteInt( AMBERINDEX(srPRestraint->iAtom2) );
                FortranWriteInt( AMBERINDEX(srPRestraint->iAtom3) );
                FortranWriteInt( AMBERINDEX(srPRestraint->iAtom4) );
                FortranWriteInt( srPRestraint->iParmIndex );
            }
        }
    }
    FortranEndLine();

        /* -27- Write the excluded atom list */
    FortranDebug( "-27-" );

    MESSAGE(( "Writing the excluded atom list\n" ));
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( vaExcludedAtoms ); i++ ) {
        FortranWriteInt( *PVAI( vaExcludedAtoms, int, i ) );
    }
    FortranEndLine();

        /* -28- Write the R^12 term for the Hydrogen bond equation */
    FortranDebug( "-28-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalHBondParms(uUnit->psParameters); i++ ) {
        ParmSetHBond( uUnit->psParameters, i, sType1, sType2, &dC, &dD, sDesc );
        FortranWriteDouble( dC );
    }
    FortranEndLine();

        /* -29- Write the R^10 term for the Hydrogen bond equation */
    FortranDebug( "-29-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalHBondParms(uUnit->psParameters); i++ ) {
        ParmSetHBond( uUnit->psParameters, i, sType1, sType2, &dC, &dD, sDesc );
        FortranWriteDouble( dD );
    }
    FortranEndLine();

        /* -30- No longer used, but stored */
    FortranDebug( "-30-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalHBondParms(uUnit->psParameters); i++ ) {
        FortranWriteDouble( 0.0 );
    }
    FortranEndLine();

        /* -31- List of atomic symbols */
    FortranDebug( "-31-" );

    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteString( 
                sAtomType( PVAI(uUnit->vaAtoms, SAVEATOMt, i )->aAtom ) );
    }
    FortranEndLine();

        /* -32- List of tree symbols */
    FortranDebug( "-32-" );

    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        aAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->aAtom;
        if ( dAtomTemp( aAtom ) == (double)'M' )
                FortranWriteString( "M  " );
        else if ( dAtomTemp( aAtom ) == (double)'E' )
                FortranWriteString( "E  " );
        else if ( dAtomTemp( aAtom ) == (double)'S' )
                FortranWriteString( "S  " );
        else if ( dAtomTemp( aAtom ) == (double)'B' )
                FortranWriteString( "B  " );
        else if ( dAtomTemp( aAtom ) == (double)'3' )
                FortranWriteString( "3  " );
        else if ( dAtomTemp( aAtom ) == (double)'4' )
                FortranWriteString( "4  " );
        else if ( dAtomTemp( aAtom ) == (double)'5' )
                FortranWriteString( "5  " );
        else if ( dAtomTemp( aAtom ) == (double)'6' )
                FortranWriteString( "6  " );
        else if ( dAtomTemp( aAtom ) == (double)'X' )
                FortranWriteString( "X  " );
        else
                FortranWriteString( "BLA" );
    }
    FortranEndLine();

        /* -33- Tree Joining information !!!!!!! Add support for this !!!!! */
    FortranDebug( "-33-" );

    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteInt( 0 );
    }
    FortranEndLine();

        /* -34- Who knows, something to do with rotating atoms */
    FortranDebug( "-34-" );

    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteInt( 0 );
    }
    FortranEndLine();

        /* -35A- The last residue before "solvent" */
                /* Number of molecules */
                /* Index of first molecule that is solvent */

    if ( bUnitUseBox(uUnit) ) {
        FortranDebug( "-35A-" );

                /* Find the index of the first solvent RESIDUE */

        for ( i=0; i<iVarArrayElementCount(uUnit->vaResidues); i++ ) {
            if ( PVAI(uUnit->vaResidues,SAVERESIDUEt,i)->sResidueType[0] ==
                RESTYPESOLVENT ) break;
        }
        iTemp = i;

        /* 
         *  Find the molecules and return the number of ATOMs in each 
         *  molecule, along with the index of the first solvent molecule
         */

        vaMolecules = vaVarArrayCreate(sizeof(int));
        zUnitIOFindAndCountMolecules( uUnit, &vaMolecules, &iFirstSolvent );

        FortranFormat( 3, INTFORMAT );
        FortranWriteInt( iTemp );
        FortranWriteInt( iVarArrayElementCount(vaMolecules) );
        FortranWriteInt( iFirstSolvent+1 );     /* FORTRAN index */

        FortranEndLine();

                /* -35B- The number of ATOMs in the Ith RESIDUE */

        FortranDebug( "-35B-" );
        FortranFormat( 12, INTFORMAT );
        for ( i=0; i<iVarArrayElementCount(vaMolecules); i++ ) {
            FortranWriteInt( *PVAI(vaMolecules,int,i) );
        }
        FortranEndLine();

                /* -35C- BETA, (BOX(I), I=1,3 ) */

        FortranDebug( "-35C-" );
        FortranFormat( 4, DBLFORMAT );
        FortranWriteDouble( dUnitBeta(uUnit)/DEGTORAD );
        UnitGetBox( uUnit, &dX, &dY, &dZ );
        FortranWriteDouble( dX );
        FortranWriteDouble( dY );
        FortranWriteDouble( dZ );
        FortranEndLine();
    }

        /* -35D- NATCAP */

    if ( bUnitUseSolventCap(uUnit) ) {
        FortranDebug( "-35D-" );
        FortranFormat( 1, INTFORMAT );
        FortranWriteInt( uUnit->iCapTempInt );
        FortranEndLine();
        
        /* -35E- CUTCAP, XCAP, YCAP, ZCAP */
        FortranDebug( "-35E-" );

        FortranFormat( 4, DBLFORMAT );
        UnitGetSolventCap( uUnit, &dX, &dY, &dZ, &dR );
        FortranWriteDouble( dR );
        FortranWriteDouble( dX );
        FortranWriteDouble( dY );
        FortranWriteDouble( dZ );
        FortranEndLine();
    }

                /* Write the perturbation information */

    if ( bPert ) {

                /* -36A- Bonds that are to be perturbed */
                        /* Totally perturbed bonds first, */
                        /* boundary second */
        FortranDebug( "-36A-" );

        FortranFormat( 12, INTFORMAT );
        if ( iVarArrayElementCount(uUnit->vaBonds) )
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                 && ((sbPBond->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
                FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaBonds) )
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                 && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
                FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
            }
        }
        FortranEndLine();

                /* -36B- Index into bond interaction arrays */
        FortranDebug( "-36B-" );

        FortranFormat( 12, INTFORMAT );


                        /* First LAMBDA = 0 */

        if ( iVarArrayElementCount(uUnit->vaBonds) )
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( sbPBond->iParmIndex );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaBonds) )
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( sbPBond->iParmIndex );
            }
        }

                        /* Then LAMBDA = 1 */

        if ( iVarArrayElementCount(uUnit->vaBonds) )
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( sbPBond->iPertParmIndex );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaBonds) )
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( sbPBond->iPertParmIndex );
            }
        }
        FortranEndLine();
        

                /* -36C- Angles that are to be perturbed */
        FortranDebug( "-36C-" );

        FortranFormat( 12, INTFORMAT );
        if ( iVarArrayElementCount(uUnit->vaAngles) )
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaAngles) )
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
            }
        }
        FortranEndLine();

                /* -36D- Index into angle interaction arrays */
        FortranDebug( "-36D-" );

        FortranFormat( 12, INTFORMAT );


                        /* First LAMBDA = 0 */

        if ( iVarArrayElementCount(uUnit->vaAngles) )
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( saPAngle->iParmIndex );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaAngles) )
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( saPAngle->iParmIndex );
            }
        }

                        /* Then LAMBDA = 1 */

        if ( iVarArrayElementCount(uUnit->vaAngles) )
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( saPAngle->iPertParmIndex );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaAngles) )
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( saPAngle->iPertParmIndex );
            }
        }
        FortranEndLine();

                /* -36E- Torsions that are to be perturbed */
        FortranDebug( "-36E-" );

        FortranFormat( 12, INTFORMAT );
        if ( iVarArrayElementCount(uUnit->vaTorsions) )
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions); 
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)==0) ) {

                if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                     (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                    MESSAGE(( "Had to turn torsion around to avoid K,L == 0\n" ));
                    MESSAGE(( "Outer atoms: %s --- %s\n", 
                                    sContainerName(aA), sContainerName(aD) ));
                    SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                    SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
                }
                if ( stPTorsion->bProper )  iProper = 1;
                else                        iProper = -1;
                if ( stPTorsion->bCalc14 )  iCalc14 = 1;
                else                        iCalc14 = -1;
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1) );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3)*iCalc14 );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaTorsions) )
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions); 
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {

                if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                     (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                    MESSAGE(( "Had to turn torsion around to avoid K,L == 0\n" ));
                    MESSAGE(( "Outer atoms: %s --- %s\n", 
                                    sContainerName(aA), sContainerName(aD) ));
                    SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                    SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
                }
                if ( stPTorsion->bProper )  iProper = 1;
                else                        iProper = -1;
                if ( stPTorsion->bCalc14 )  iCalc14 = 1;
                else                        iCalc14 = -1;
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1) );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3)*iCalc14 );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            }
        }
        FortranEndLine();

                /* -36F- Index into torsion interaction arrays */
        FortranDebug( "-36F-" );

        FortranFormat( 12, INTFORMAT );

                        /* First LAMBDA = 0 */

        if ( iVarArrayElementCount(uUnit->vaTorsions) )
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions);
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( stPTorsion->iParmIndex );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaTorsions) )
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions); 
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( stPTorsion->iParmIndex );
            }
        }

                        /* Then LAMBDA = 1 */

        if ( iVarArrayElementCount(uUnit->vaTorsions) )
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions); 
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( stPTorsion->iPertParmIndex );
            }
        }
        if ( iVarArrayElementCount(uUnit->vaTorsions) )
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions);
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( stPTorsion->iPertParmIndex );
            }
        }
        FortranEndLine();

                /* -36G- Residue labels at LAMBDA = 1 */
                        /* Just write the labels at LAMBDA = 0 */

                /* Trim the string down to at most 3 characters by */
                /* taking the last three characters if it is too long */
        FortranDebug( "-36G-" );

        MESSAGE(( "Writing the residue labels\n" ));
        FortranFormat( 20, LBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaResidues); i++ ) {
            cPTemp = PVAI( uUnit->vaResidues, SAVERESIDUEt, i )->sName;
            if ( strlen( cPTemp ) > 3 ) cPTemp += ( strlen(cPTemp)-3 );
            FortranWriteString( cPTemp );
        }
        FortranEndLine();

                /* -36H- Atom names at LAMBDA = 0 */
        FortranDebug( "-36H-" );

        FortranFormat( 20, LBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            cPTemp = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->sPertName;
            if ( strlen( cPTemp ) == 0 )
                cPTemp = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->sName;
            if ( strlen( cPTemp ) > 4 ) cPTemp += ( strlen(cPTemp)-4 );
            FortranWriteString( cPTemp );
        }
        FortranEndLine();

                /* -36I- List of atomic symbols (atom types??????) */
        FortranDebug( "-36I-" );

        FortranFormat( 20, LBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            cPTemp = sAtomPertType(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom);
            if ( strlen(cPTemp) == 0 )
                cPTemp = sAtomType(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom);
            if ( strlen( cPTemp ) > 3 ) cPTemp += ( strlen(cPTemp)-3 );
            FortranWriteString( cPTemp );
        }
        FortranEndLine();

                /* -36J- Value of LAMBDA for each ATOM ????????? */
                /* TODO: Figure out what the hell this is */
        FortranDebug( "-36J-" );

        FortranFormat( 5, DBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            FortranWriteDouble( 0.0 );
        }
        FortranEndLine();

                /* -36K- Flag to tell whether the atom is perturbed */
        FortranDebug( "-36K-" );

        FortranFormat( 12, INTFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            if ( bAtomPerturbed(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom) ) 
                FortranWriteInt( 1 );
            else FortranWriteInt( 0 );
        }
        FortranEndLine();

                /* -36L- List of atom types - IACPER */
        FortranDebug( "-36L-" );

        FortranFormat( 12, INTFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            iAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->iPertTypeIndex-1;
            iTemp = *PVAI( vaNBIndex, int, iAtom );
            FortranWriteInt( iTemp+1 );
        }
        FortranEndLine();

                /* -36M- Perturbed charges */
        FortranDebug( "-36M-" );

        FortranFormat( 5, DBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            aAtom = PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom;
            if ( bAtomPerturbed( aAtom ) )
                FortranWriteDouble( ELECTRONTOKCAL *
                    dAtomPertCharge(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom));
            else
                FortranWriteDouble( ELECTRONTOKCAL *
                    dAtomCharge(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom));
        }
        FortranEndLine();

    }

    /*
     *  polarizabilities
     */
    if ( bPolar ) {
        iCount = 0;
        iCountPerturbed = 0;
        iMax = iVarArrayElementCount(uUnit->vaAtoms);
        MESSAGE(( "Writing the atomic polarizabilities\n" ));         
        FortranFormat( 5, DBLFORMAT );
        saPAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, 0 );
        for ( i=0; i<iMax; i++, saPAtom++ ) {
            iIndex = iParmSetFindAtom( uUnit->psParameters, saPAtom->sType );
            ParmSetAtom( uUnit->psParameters, iIndex, sType,
                        &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14, &dRStar14,
			&dScreenF,
            &iElement, &iHybridization, sDesc );
            if ( dPolar == -1.0 ) {
                dPolar = 0.0;
                iCount++;
            }
            FortranWriteDouble( dPolar );
        }
        if ( iCount > 0 )
                VP0(( "Total atoms with default polarization=0.0: %d of %d\n",
                                                        iCount, iMax));
                                                
	/*
        FortranEndLine();
	if (GDefaults.dDipoleDampFactor > 1.0) {
           FortranFormat(1, "%-80s");
           FortranWriteString("%FLAG DIPOLE_DAMP_FACTOR");
           FortranWriteString("%FORMAT(5E16.8)");
           FortranFormat(5, DBLFORMAT);
           for (i = 0; i < iMax; i++, saPAtom++) {
               iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
               ParmSetAtom(uUnit->psParameters, iIndex, sType,
                           &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                           &dRStar14, &dScreenF, &iElement, &iHybridization,
			   sDesc); 
	       if (dScreenF == 0.0) {
		  dScreenF = GDefaults.dDipoleDampFactor;
	       }
               FortranWriteDouble(dScreenF);
	   }                       
	}
	*/

        FortranEndLine();
        if ( bPert ) {
                int     iPertTot = 0;
                saPAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, 0 );
                for ( i=0; i<iMax; i++, saPAtom++ ) {
                        BOOL    bTmp = bAtomPerturbed( saPAtom->aAtom );

                        if ( bTmp ) {
                                iIndex = iParmSetFindAtom( uUnit->psParameters, 
                                                        saPAtom->sPertType );
                                iPertTot++;
                        } else
                                iIndex = iParmSetFindAtom( uUnit->psParameters, 
                                                        saPAtom->sType );
                        ParmSetAtom( uUnit->psParameters, iIndex, sType,
                                        &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                    &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc );
                        if ( dPolar == -1.0 ) {
                                dPolar = 0.0;
                                if ( bTmp ) iCountPerturbed++;
                        }
                        FortranWriteDouble( dPolar );
                }
                if ( iCountPerturbed > 0 )
                    VP0(( 
                     "Total pert atoms with default polarization=0.0: %d of %d\n",
                                        iCountPerturbed, iPertTot ));
        }
    }

        /*  Charmm-style parameters  */

        if( GDefaults.iCharmm ){
        /* -19- Lennard jones r**12 term for all 14 interactions */
                /* CN114 array */
        FortranDebug( "-19-" );

        FortranFormat( 5, DBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(vaNBParameters); i++ ) {
                FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dA14 );
        }
        FortranEndLine();
    
        /* -20- Lennard jones r**6 term for all 14 interactions */
                /* CN214 array */
        FortranDebug( "-20-" );

        FortranFormat( 5, DBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(vaNBParameters); i++ ) {
                FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dC14 );
        }
        FortranEndLine();

        /* -13- Force constants for Urey-Bradley */
                FortranDebug( "-13-" );

                FortranFormat( 5, DBLFORMAT );
                for ( i=0; i<iParmSetTotalAngleParms(uUnit->psParameters); i++ ) {
                        ParmSetAngle( uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                                &dKt, &dT0, &dTkub, &dRkub, sDesc );
                        FortranWriteDouble( dTkub );
                }
                FortranEndLine();

        /* -14- Equilibrium distances for Urey-Bradley*/
                FortranDebug( "-14-" );

                FortranFormat( 5, DBLFORMAT );
                for ( i=0; i<iParmSetTotalAngleParms(uUnit->psParameters); i++ ) {
                        ParmSetAngle( uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                                                &dKt, &dT0, &dTkub, &dRkub, sDesc );
                        FortranWriteDouble( dRkub );
                }
                FortranEndLine();
 
        }


        /********************************************************/
        /* Write the coordinate file                            */
        /********************************************************/

    FortranFile( fCrd );

    FortranFormat( 1, "%s" );
    FortranWriteString( sContainerName( uUnit ) );
    FortranEndLine();

    FortranFormat( 1, "%5d" );
    FortranWriteInt( iVarArrayElementCount( uUnit->vaAtoms ) );
    FortranEndLine();

    FortranFormat( 6, "%12.7lf" );
    if ( bUnitUseBox(uUnit) ) {
        double  dX2, dY2, dZ2;

        UnitGetBox( uUnit, &dX, &dY, &dZ );
        dX2 = dX * 0.5;
        dY2 = dY * 0.5;
        dZ2 = dZ * 0.5;

        /*
         *  shift box to Amber spot; later, add a cmd opt or environment
         *      var to switch between 0,0,0 center (spasms) or corner
         */
        for ( i = 0; i<iVarArrayElementCount( uUnit->vaAtoms ); i++ ) {
            vPos = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->vPos;
            FortranWriteDouble( dVX(&vPos) + dX2 );
            FortranWriteDouble( dVY(&vPos) + dY2 );
            FortranWriteDouble( dVZ(&vPos) + dZ2 );
        }
        FortranEndLine();
        FortranWriteDouble( dX );
        FortranWriteDouble( dY );
        FortranWriteDouble( dZ );
        FortranWriteDouble( dUnitBeta( uUnit ) / DEGTORAD );
        FortranWriteDouble( dUnitBeta( uUnit ) / DEGTORAD );
        FortranWriteDouble( dUnitBeta( uUnit ) / DEGTORAD );
        FortranEndLine();
    } else {
        for ( i = 0; i<iVarArrayElementCount( uUnit->vaAtoms ); i++ ) {
            vPos = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->vPos;
            FortranWriteDouble( dVX(&vPos) );
            FortranWriteDouble( dVY(&vPos) );
            FortranWriteDouble( dVZ(&vPos) );
        }
        FortranEndLine();
    }

    VarArrayDestroy( &vaNBIndexMatrix );
    VarArrayDestroy( &vaNBParameters );
    VarArrayDestroy( &vaExcludedAtoms );
    VarArrayDestroy( &vaExcludedCount );
    VarArrayDestroy( &vaNBIndex );
    VarArrayDestroy( &vaNonBonds );
    fclose( fCrd );

}
static char *prepfmt = "   %-3d %-4s  %-4s  %c      %f  %f  %f    %f\n";

/*
 *  PrintSideChain() - depth-first descent of side chain,
 *        printing atoms, noting loops
 */
static void
PrintSideChain(FILE * fOut, ATOM aParentAtom, ATOM aAtom,
               int *iP, VARARRAY vaLoopAtoms)
{
    VECTOR vPos;
    int i, j, k, iChildTag = *iP;
    ATOM aChildAtom, aNbrs[MAXBONDS];

    /*
     *  figure tree type - count 'downstream', 
     *      undesignated neighbors, i.e. omits
     *      parent and anything else already 
     *      encountered. Mark all claimed atoms
     *      as belonging to this one for the
     *      benefit of recursion looping back
     *      around to this point before getting
     *      to the atom later in this routine
     */
/*
fprintf(stderr, "-- parent %s -> %s coord %d\n", 
sAtomName(aParentAtom),sAtomName(aAtom), iAtomCoordination(aAtom));
*/
    j = 0;
    for (i = 0; i < iAtomCoordination(aAtom); i++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, i);
        if (aChildAtom == aParentAtom)
            continue;
/*
fprintf(stderr, "   child %s type %c\n", 
sAtomName(aChildAtom), (char)dAtomTemp( aChildAtom ));
*/
        if (dAtomTemp(aChildAtom) == (double) 'x' &&
            iAtomTempInt(aChildAtom) == -1) {
            AtomSetTempInt(aChildAtom, iChildTag);
            j++;
        }
    }
    SetTreeType(aAtom, j);

    /*
     *  print
     */
    vPos = vAtomPosition(aAtom);
    fprintf(fOut, prepfmt, *iP,
            sAtomName(aAtom),
            sAtomType(aAtom),
            (char) dAtomTemp(aAtom),
            vPos.dX, vPos.dY, vPos.dZ, dAtomCharge(aAtom));
    AtomSetSeenId(aAtom, *iP);
    (*iP)++;

    /*
     *  put eligible children in order for readability -
     *      get obvious 'E' types 1st
     */
    k = 0;
    for (j = 0; j < iAtomCoordination(aAtom); j++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, j);
        if (aChildAtom == aParentAtom)
            continue;
        if (dAtomTemp(aChildAtom) != (double) 'x')
            continue;
        if (iAtomCoordination(aChildAtom) == 1)
            aNbrs[k++] = aChildAtom;
    }
    /*
     *  bubble sort them into alphabetical order
     */
    while (1) {
        int iMore = 0;

        for (j = 0; j < k - 1; j++) {
            if (strcmp(sAtomName(aNbrs[j]), sAtomName(aNbrs[j + 1])) > 0) {
                ATOM aTmp;
                aTmp = aNbrs[j + 1];
                aNbrs[j + 1] = aNbrs[j];
                aNbrs[j] = aTmp;
                iMore++;
            }
        }
        if (!iMore)
            break;
    }
    /*
     *  print
     */
    for (j = 0; j < k; j++)
        PrintSideChain(fOut, aAtom, aNbrs[j], iP, vaLoopAtoms);

    /*
     *  get remaining eligible children
     */
    k = 0;
    for (j = 0; j < iAtomCoordination(aAtom); j++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, j);
        if (aChildAtom == aParentAtom)
            continue;
        if (iAtomCoordination(aChildAtom) == 1)        /* done above */
            continue;

        if (dAtomTemp(aChildAtom) == (double) 'x' &&
            iAtomTempInt(aChildAtom) == iChildTag) {
            aNbrs[k++] = aChildAtom;
        } else {
            STRING sTemp;

            /*
             *  atom seen already: may have found a new 
             *      loop closing bond
             */
            if (strcmp(sAtomName(aAtom), sAtomName(aChildAtom)) < 0)
                sprintf(sTemp, "%-4s %-4s",
                        sAtomName(aAtom), sAtomName(aChildAtom));
            else
                sprintf(sTemp, "%-4s %-4s",
                        sAtomName(aChildAtom), sAtomName(aAtom));
            for (i = 0; i < iVarArrayElementCount(vaLoopAtoms); i++) {
                if (!strcmp(sTemp, PVAI(vaLoopAtoms, char, i)))
                     break;
            }
            if (i == iVarArrayElementCount(vaLoopAtoms))
                VarArrayAdd(vaLoopAtoms, (GENP) sTemp);
        }
    }
    /*
     *  bubble sort
     */
    while (1) {
        int iMore = 0;

        for (j = 0; j < k - 1; j++) {
            if (strcmp(sAtomName(aNbrs[j]), sAtomName(aNbrs[j + 1])) > 0) {
                ATOM aTmp;
                aTmp = aNbrs[j + 1];
                aNbrs[j + 1] = aNbrs[j];
                aNbrs[j] = aTmp;
                iMore++;
            }
        }
        if (!iMore)
            break;
    }
    /*
     *  print
     */
    for (j = 0; j < k; j++)
        PrintSideChain(fOut, aAtom, aNbrs[j], iP, vaLoopAtoms);
}
static void WritePrepRes(RESIDUE rRes, FILE * fOut)
{
    int i, j, iMax;
    LOOP lTemp;
    ATOM aAtom, aAtom0, aAtom1, aChildAtom, aParentAtom;
    VARARRAY vaLoopAtoms;
    char *cPResName;

    cPResName = sContainerName(rRes);
    VP0(("  saving prep, residue %s\n", cPResName));

    /*
     *  mark main chain atoms; worry about side chains later
     */

    if (MarkMainChainAtoms(rRes, 1) < 1)
        return;

    aAtom0 = (ATOM) rRes->aaConnect[0];
    aAtom1 = (ATOM) rRes->aaConnect[1];


    /*
     *  write residue header & dummy atoms
     */
    if (strlen(cPResName) > 3)
        cPResName += strlen(cPResName) - 3;
    fprintf(fOut, " leap-generated prep residue\n");
    fprintf(fOut, "%s.res\n", cPResName);
    fprintf(fOut, "%-3s  INT     0\n", cPResName);
    fprintf(fOut, "CHANGE   NOMIT DU   BEG \n");
    fprintf(fOut, "   0.00000\n");
    fprintf(fOut, "   %s\n   %s\n   %s\n",
            "1   DUMM  DU    M      0.000000  0.000000  0.000000  0.0",
            "2   DUMM  DU    M      1.000000  0.000000  0.000000  0.0",
            "3   DUMM  DU    M      1.000000  1.000000  0.000000  0.0");

    /*
     *  traverse tree from top, following main chain
     *      and ordering so that branches precede further
     *      main chain, printing atoms
     *
     *      reset atom temp ints to use in marking 'seen' atoms
     *      and SeenId for tracking order in file for ordering
     *      atoms in impropers; dAtomTemps were set in 
     *      MarkMainChainAtoms()
     */
    lTemp = lLoop((OBJEKT) rRes, ATOMS);
    while ((aAtom = (ATOM) oNext(&lTemp)) != NULL) {
        AtomSetTempInt(aAtom, -1);
        AtomSetSeenId(aAtom, -1);
    }

    vaLoopAtoms = vaVarArrayCreate(2 * ATOMTYPELEN + 1);
    i = 4;
    aAtom = aAtom0;
    aParentAtom = NULL;
    while (1) {
        ATOM aNextMain;
        VECTOR vPos;

        vPos = vAtomPosition(aAtom);
        fprintf(fOut, prepfmt, i,
                sAtomName(aAtom),
                sAtomType(aAtom),
                (char) dAtomTemp(aAtom),
                vPos.dX, vPos.dY, vPos.dZ, dAtomCharge(aAtom));
        AtomSetSeenId(aAtom, i++);
        /*
         *  find next main chain down in this residue, 
         *      marking side chains
         */
        for (j = 0; j < iAtomCoordination(aAtom); j++) {
            aChildAtom = aAtomBondedNeighbor(aAtom, j);
            if (aChildAtom == aParentAtom)
                continue;
            /*
             *  stay within residue
             */
            if (cContainerWithin(aChildAtom) != cContainerWithin(aAtom))
                continue;
            if (dAtomTemp(aChildAtom) == (double) 'M') {
                aNextMain = aChildAtom;
            } else if (dAtomTemp(aChildAtom) == (double) 'x') {
                /*
                 *  side chain not seen before
                 */
                PrintSideChain(fOut, aAtom, aChildAtom, &i, vaLoopAtoms);
            }
        }
        if (aAtom == aAtom1)
            break;
        aParentAtom = aAtom;
        aAtom = aNextMain;
    }

    /*
     *  atoms written; write extra line to indicate no
     *      internal coords follow
     */
    fprintf(fOut, "\n");

    /*
     *  write any LOOP closing bonds
     */
    iMax = iVarArrayElementCount(vaLoopAtoms);
    if (iMax) {
        fprintf(fOut, "\nLOOP\n");
        for (i = 0; i < iMax; i++) {
            fprintf(fOut, " %s\n", PVAI(vaLoopAtoms, char, i));
        }
    }

    /*
     *  write list of potential IMPROPERs
     */
    lTemp = lLoop((OBJEKT) rRes, IMPROPERS);
    for (i = 0; oNext(&lTemp);) {
        ATOM aAtomA, aAtomB, aAtomC, aAtomD;
        int iSame;

        LoopGetImproper(&lTemp, &aAtomA, &aAtomB, &aAtomC, &aAtomD);

        if (aAtomC == aAtom0 || aAtomC == aAtom1)
            continue;
        if (iAtomCoordination(aAtomC) != 3)
            continue;

        if (!i++)
            fprintf(fOut, "\nIMPROPER\n");

        /*
         *  sort peripheral atoms by type
         */
        iSame = strcmp(sAtomType(aAtomA), sAtomType(aAtomB));
        if (iSame > 0 ||
            (!iSame && iAtomSeenId(aAtomA) > iAtomSeenId(aAtomB))) {
            aChildAtom = aAtomB;
            aAtomB = aAtomA;
            aAtomA = aChildAtom;
        }
        iSame = strcmp(sAtomType(aAtomB), sAtomType(aAtomD));
        if (iSame > 0 ||
            (!iSame && iAtomSeenId(aAtomB) > iAtomSeenId(aAtomD))) {
            aChildAtom = aAtomD;
            aAtomD = aAtomB;
            aAtomB = aChildAtom;
        }
        iSame = strcmp(sAtomType(aAtomA), sAtomType(aAtomB));
        if (iSame > 0 ||
            (!iSame && iAtomSeenId(aAtomA) > iAtomSeenId(aAtomB))) {
            aChildAtom = aAtomB;
            aAtomB = aAtomA;
            aAtomA = aChildAtom;
        }
        fprintf(fOut, " %-4s %-4s %-4s %-4s\n",
                sAtomName(aAtomA), sAtomName(aAtomB),
                sAtomName(aAtomC), sAtomName(aAtomD));
    }
    if (iAtomCoordination(aAtom0) == 2) {
        ATOM aAtomA, aAtomB;
        if (!i++)
            fprintf(fOut, "\nIMPROPER\n");
        aAtomA = aAtomBondedNeighbor(aAtom0, 0);
        aAtomB = aAtomBondedNeighbor(aAtom0, 1);
        if (strcmp(sAtomType(aAtomA), sAtomType(aAtomB)) > 0) {
            aChildAtom = aAtomB;
            aAtomB = aAtomA;
            aAtomA = aChildAtom;
        }
        fprintf(fOut, " %-4s %-4s %-4s %-4s\n",
                "-M", sAtomName(aAtomA), sAtomName(aAtom0),
                sAtomName(aAtomB));
    }
    if (iAtomCoordination(aAtom1) == 2) {
        ATOM aAtomA, aAtomB;
        if (!i++)
            fprintf(fOut, "\nIMPROPER\n");
        aAtomA = aAtomBondedNeighbor(aAtom1, 0);
        aAtomB = aAtomBondedNeighbor(aAtom1, 1);
        if (strcmp(sAtomType(aAtomA), sAtomType(aAtomB)) > 0) {
            aChildAtom = aAtomB;
            aAtomB = aAtomA;
            aAtomA = aChildAtom;
        }
        fprintf(fOut, " %-4s %-4s %-4s %-4s\n",
                sAtomName(aAtomA), "+M", sAtomName(aAtom1),
                sAtomName(aAtomB));
    }

    /*
     *  done
     */
    fprintf(fOut, "\nDONE\n");


    /*
     *  clean up
     */
    VarArrayDestroy(&vaLoopAtoms);

}

void UnitIOSaveAmberPrep(UNIT uUnit, FILE * fOut)
{
    LOOP lResidues;
    RESIDUE rRes;

    /*
     *  write a default beginning-of-prep
     */
    fprintf(fOut, " 0 0 0\n\n");

    /*
     *  put each residue in
     */
    lResidues = lLoop((OBJEKT) uUnit, RESIDUES);
    while ((rRes = (RESIDUE) oNext(&lResidues)))
        WritePrepRes(rRes, fOut);
    /*
     *  terminate file
     */
    fprintf(fOut, "STOP\n");

}
