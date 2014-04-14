/*
 *      File:  mol2File.c based on pdbFile.c 
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
 *   mol2File.c      based on pdbFile.c
 * 
 *      Christine Cezard (2007)
 *
 *     Universite de Picardie - Jules Verne, Amiens
 *     http://q4md-forcefieldtools.org
 *
 *      Description:
 *              Writes a Mol2 file
 */

#include        <limits.h>

#include        "basics.h"

#include        "classes.h"
#include        "dictionary.h"
#include        "database.h"
#include        "pdb.h"
#include        "varArray.h"
#include        "avl.h"
#include        "sort.h"
#include        "matrix.h"
#include        "model.h"
#include        "minimizer.h"
#include        "leap.h"
#include        "defaults.h"
#include        "atom.h"
#include        "tripos.h"
#include        "stdio.h"
#include        "parmLib.h"
#include        "unitio.h"



                /* MOL2WRITEt is used to store data for writing Mol2 files */

typedef struct  {
        FILE    *fPdbFile;
        int     iRecordNumber;
        int     iResidueSeq;
        STRING  sResidueName;
} MOL2WRITEt;

typedef struct  {
        char    sName[10];
        int     iTerminator;
        int     iPdbSequence;
} RESIDUENAMEt;

typedef struct  {
        BOOL    bUsed;
        MATRIX  mTransform;
} PDBMATRIXt;

typedef struct {
    int iAtom1;
    int iAtom2;
    FLAGS fFlags;
} SAVECONNECTIVITYt;

typedef struct {
    int iAtom1;
    int iAtom2;
    int iParmIndex;
    int iPertParmIndex;
    FLAGS fFlags;
} SAVEBONDt;

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

typedef struct  {
        FILE            *fPdbFile;
        VARARRAY        vaUnits;                
        VARARRAY        vaResidues;
        VARARRAY        vaResidueSeq;
		VARARRAY        vaAtoms;
	    VARARRAY        vaMatrices;
		UNIT            uUnit;
		RESIDUE         rRes;
		int             iPdbSequence;
        int             iNextUnit;
        int             iMaxSerialNum;
} PDBREADt;


/*
----------------------------------------------------------------------

        Static variables
        
*/

//static  DICTIONARY      SdResidueNameMap = NULL;
//static  DICTIONARY      SdAtomNameMap = NULL;



#define SEQUENCE_NUMBER_STEPS   10000

extern int zUnitIOAmberOrderResidues( UNIT );

static void
zMol2FileBegin( MOL2WRITEt *pwPFile, FILE *fOut )
{
    pwPFile->fPdbFile = fOut;
    pwPFile->iRecordNumber= 1;
    pwPFile->iResidueSeq = 1;
}

void
zMol2FileWriteAtomRecord( MOL2WRITEt *pwPFile, ATOM aAtom, int choice )
{
pdb_record      p;
#ifdef DEBUG
STRING          sElement;
#endif 
STRING          sName, sTemp, sTypeAt;

    p.pdb.atom.serial_num = pwPFile->iRecordNumber++;
    strcpy( p.pdb.atom.residue.name, pwPFile->sResidueName );

	strcpy( sName, sAtomType(aAtom) );
	strcpy( sTemp, sContainerName((CONTAINER)aAtom) );
	    
    
	if (choice == 0) {
	  strcpy( sTypeAt, sTemp);
	  sTypeAt[1] = '\0';
	  strcpy( p.pdb.atom.type_at, sTypeAt );
	  }
	if (choice == 1) {
	  strcpy( sTypeAt, sName);
	  strcpy( p.pdb.atom.type_at, sTypeAt );
      }	  
	  
    strcpy( p.pdb.atom.name, sTemp ); 

	MESSAGE(( "Element: |%s|   pdb_name=|%s|\n", sElement, sName ));

    p.pdb.atom.residue.chain_id = ' ';
    p.pdb.atom.residue.seq_num = pwPFile->iResidueSeq;
    p.pdb.atom.residue.insert_code=' ';
/*    p.pdb.atom.alt_loc = ' ' ; */
    p.pdb.atom.x = dVX(&vAtomPosition(aAtom));
    p.pdb.atom.y = dVY(&vAtomPosition(aAtom));
    p.pdb.atom.z = dVZ(&vAtomPosition(aAtom));
    p.pdb.atom.occupancy = 1.0;
    p.pdb.atom.temp_factor = dAtomCharge(aAtom);
    p.pdb.atom.ftnote_num = 0;
    p.record_type = MOL2_ATOM;
    pdb_write_record( pwPFile->fPdbFile, &p, NULL, 0 );

}


static void
zMol2FileWriteContainer( MOL2WRITEt *pwPFile, CONTAINER cCont, int choice )
{
//const char      C_TERMINAL_PREFIX = 'C';
//const char      N_TERMINAL_PREFIX = 'N';
const int       RESIDUE_NAME_LENGTH = 4;
LOOP            lContents;
ATOM            aAtom;
char            *cPTemp;

   cPTemp = sContainerName(cCont);
        strncpy( pwPFile->sResidueName, cPTemp, RESIDUE_NAME_LENGTH );
      
        pwPFile->sResidueName[ RESIDUE_NAME_LENGTH ] = '\0';
       if ( strlen(cPTemp) > RESIDUE_NAME_LENGTH ) {
            VP0(( " Truncating residue name for PDB format: %s -> %s\n", 
                        sContainerName(cCont), pwPFile->sResidueName ));
        } 
 /*   } 
*/
    if ( iObjectType(cCont) == ATOMid ) {
        zMol2FileWriteAtomRecord( pwPFile, (ATOM)cCont, choice);
        pwPFile->iResidueSeq++;
    } else if ( iObjectType(cCont) == RESIDUEid ) {
        RESIDUE rRes = (RESIDUE) cCont;
        pwPFile->iResidueSeq = rRes->iTemp;
        if ( pwPFile->iResidueSeq == 0 )
                pwPFile->iResidueSeq = 1;
        lContents = lLoop( (OBJEKT)cCont, DIRECTCONTENTSBYSEQNUM );
        while ( (aAtom = (ATOM)oNext(&lContents)) ) {
            zMol2FileWriteAtomRecord( pwPFile, aAtom, choice);
			}
    } else if ( iObjectType(cCont) == MOLECULEid ) {
        lContents = lLoop( (OBJEKT)cCont, ATOMS );
        while ( (aAtom = (ATOM)oNext(&lContents)) ) {
            zMol2FileWriteAtomRecord( pwPFile, aAtom, choice );
        }
        pwPFile->iResidueSeq++;
    }
}

static char *
zMol2FileWriteResidueContainer( MOL2WRITEt *pwPFile, CONTAINER cCont )
{
//const char      C_TERMINAL_PREFIX = 'C';
//const char      N_TERMINAL_PREFIX = 'N';
const int       RESIDUE_NAME_LENGTH = 3;
//LOOP            lContents;
char            *cPTemp;
RESIDUE         rRes;

     cPTemp = sContainerName(cCont);
     strncpy( pwPFile->sResidueName, cPTemp, RESIDUE_NAME_LENGTH );
        /* The intentional side effect is to truncate long names. */
     pwPFile->sResidueName[ RESIDUE_NAME_LENGTH ] = '\0';
     if ( strlen(cPTemp) > RESIDUE_NAME_LENGTH ) {
        VP0(( " Truncating residue name for PDB format: %s -> %s\n", 
                   sContainerName(cCont), pwPFile->sResidueName ));
        } 
     rRes = (RESIDUE) cCont;

     return(sContainerName(cCont));	
}


void
Mol2Write( FILE *fOut, UNIT uUnit, int choice )
{

//      int             i, iResidueCount, iAtomCount = 0, iAtom, iBondCount, iCount;
        int             i, iResidueCount, iAtomCount = 0, iBondCount, iCount;
//`		int             j,k,l,m,n;
		int             j,k;
        LOOP            lContents,lTemp,lResidues;
        SAVERESIDUEt    *srPResidue;
		SAVECONNECTIVITYt *scPCon;
		MOL2WRITEt      pwFile;
		STRING     sTemp;
		char     *sName; 
		ATOM       aAtom1,aAtom2;
		RESIDUE rRes1;
        BOOL   bPert,bFailedGeneratingParameters;
//		STRING sAtom1, sAtom2, sDesc;

		
zMol2FileBegin( &pwFile, fOut );

/* @<TRIPOS>MOLECULE Bloc */
fprintf(fOut, "@<TRIPOS>MOLECULE\n") ;
iResidueCount = zUnitIOAmberOrderResidues( uUnit );	


strcpy( sTemp, sContainerName((CONTAINER) uUnit));

iCount = 0;
bPert = FALSE;
lTemp = lLoop((OBJEKT) uUnit, BONDS);
 while (oNext(&lTemp) != NULL)
        iCount++;

iAtomCount = 0;
        lTemp = lLoop((OBJEKT) uUnit, ATOMS);
        while (oNext(&lTemp) != NULL)
            iAtomCount++;		
		
iBondCount = iCount;

fprintf(fOut, "%s\n", sTemp) ;	
fprintf(fOut, "%5d %5d %5d     0     1 \n", iAtomCount,iBondCount,iResidueCount) ;
fprintf(fOut, "SMALL\n");
fprintf(fOut, "USER_CHARGES\n");



/* @<TRIPOS>ATOM Bloc */
fprintf(fOut, "@<TRIPOS>ATOM\n") ;	

srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (i = 0; i < iResidueCount; srPResidue++, i++) {
                RESIDUE rRes = srPResidue->rResidue;
				rRes->iTemp = i + 1;
                zMol2FileWriteContainer( &pwFile, (CONTAINER) rRes, choice); 
		}
fprintf(fOut, "@<TRIPOS>BOND\n") ;	

 /* Now generate the connectivity table */
 /* Inspired by zUnitIOBuildTables in unitio.c */

 /*zbUnitIOIndexBondParameters(plParameters, uUnit, bPert); */
 
        iAtomCount = 0;
        lTemp = lLoop((OBJEKT) uUnit, ATOMS);
        while (oNext(&lTemp) != NULL)
            iAtomCount++;

        if (iAtomCount) {
            uUnit->vaAtoms = vaVarArrayCreate(sizeof(SAVEATOMt));
            VarArraySetSize((uUnit->vaAtoms), iAtomCount);
            i = 0;
     

           lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTSBYSEQNUM);
            while ((rRes1 = (RESIDUE) oNext(&lResidues)) != NULL) {
            zUnitDoAtoms(uUnit, NULL, rRes1, &i, &bFailedGeneratingParameters, bPert);
            }
         }  
    
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
            /*LoopGetBond(&lTemp, &aAtom1, &aAtom2); */
			 (*(&aAtom1))=(ATOM)(&lTemp)->oaObj[0] ;
			 (*(&aAtom2))=(ATOM)(&lTemp)->oaObj[1] ;
            i++;
            scPCon++;
			fprintf(fOut, "%5d %5d %5d 1\n", i,(((CONTAINER)(aAtom1))->iTempInt),(((CONTAINER)(aAtom2))->iTempInt));
        }
    }
 

fprintf(fOut, "@<TRIPOS>SUBSTRUCTURE\n") ;	
srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
k=0;
for (i = 0; i < iResidueCount; srPResidue++, i++) {
				ATOM aAtom;
                RESIDUE rRes = srPResidue->rResidue;
                j=0;
				sName = zMol2FileWriteResidueContainer( &pwFile, (CONTAINER) rRes); 
                lContents = lLoop((OBJEKT)rRes, DIRECTCONTENTSBYSEQNUM );
				while((aAtom = (ATOM)oNext(&lContents)) )
				j=j+1 ;
				k=k+j;
				fprintf(fOut, "%7d %4s %14d ****               0 ****  **** \n", i+1, sName, k-j+1); 
						
		} 
}
 

