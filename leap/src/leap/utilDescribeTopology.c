/*
 *      File:   utilDescTop.c
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
 *      Description:
 *              This program dumps a description of all the interactions
 *              in an AMBER topology file.
 *              It can be used for debugging topology file generators.
 */




/* #include        "classes.h" */

#ifdef sgi 
#include <getopt.h> 
#endif

#include	"basics.h"

#include	"vector.h"

#include        "varArray.h"
#include        "dictionary.h"
#include        "fortran.h"

typedef char    LABELTYPE[90];

usage(arg0)
char	*arg0;
{
        fprintf( stderr, 
	  "Usage: %s [options] {topologyFile} {coordinateFile}\n", arg0 );
        fprintf( stderr, "Options:\n" );
	fprintf( stderr, "  D  - Describe the load\n" );
        fprintf( stderr, "  n  - Include atom index numbers in names.\n" );
        fprintf( stderr, "  a  - Include atom labels in names.\n" );
        fprintf( stderr, "  r  - Include residue names in atom names.\n" );
        fprintf( stderr, "  t  - Include atom types in names.\n" );
	fprintf( stderr, "  I  - Print indices for some interactions.\n" );
        fprintf( stderr, "  e  - Print those energies that are calculated.\n" );
	fprintf( stderr, "  i  - Ignore non-bonds (This speeds things up alot).\n" );
	fprintf( stderr, "  s  - Sort torsion atom names\n" );
        exit(2);
}

typedef	struct	{
		int		iIntI;
		int		iIntJ;
		int		iIntK;
		int		iIntL;
		GENP		PData;
		} BLOCKt;



/*
 *      Global variables
 */

static	BOOL    SbIncludeAtomNumberInName = FALSE;
static	BOOL    SbIncludeAtomTypeInName = FALSE;
static	BOOL    SbIncludeAtomName = FALSE;
static	BOOL    SbIncludeResidueName = FALSE;
static	BOOL    SbIncludeIndices = FALSE;
static	BOOL    SbPrintEnergy = FALSE;
static	BOOL	SbIgnoreNonbonds = FALSE;
static	BOOL	SbDescLoad = FALSE;
static	BOOL	SbSortTorsions = FALSE;

int             ntotat, ntypes, nbonh, nbona, ntheth, ntheta, nphih, nphia;
int             jhparm, jparm, next, ntotrs, mbona, mtheta, mphia, mumbnd;
int             mumang, mptra, natyp, nhb, ifpert, nbper, ngper, ndper, mbper;
int             mgper, mdper, ifbox, nmxrs, ifcap;
VARARRAY        vaIgraph;
VARARRAY        vaChrg;
VARARRAY        vaAmass;
VARARRAY        vaIac;
VARARRAY        vaNumex;
VARARRAY        vaNno;
VARARRAY        vaLabres;
VARARRAY        vaIpres;
VARARRAY        vaRk;
VARARRAY        vaReq, vaTk, vaTeq, vaPk, vaPn, vaPhase, vaSolty;
VARARRAY        vaCn1, vaCn2;
VARARRAY        vaIbh, vaJbh, vaIcbh;
VARARRAY        vaIb, vaJb, vaIcb;
VARARRAY        vaIth, vaJth, vaKth, vaIcth;
VARARRAY        vaIt, vaJt, vaKt, vaIct;
VARARRAY        vaIph, vaJph, vaKph, vaLph, vaIcph;
VARARRAY        vaIp, vaJp, vaKp, vaLp, vaIcp;
VARARRAY        vaNatex;
VARARRAY        vaAg, vaBg, vaHbcut, vaIsymbl, vaItree, vaJoin, vaIrotat;
VARARRAY	vaCoords;

int		iptres, nspm, nspsol;
VARARRAY	vaNsp;
double		beta, box1, box2, box3;

int		iNatcap;
double		dCutcap, dXcap, dYcap, dZcap;

VARARRAY	vaIbper, vaJbper;
VARARRAY	vaIcbper;

double		cutcap;
double		xcap;
double		ycap;
double		zcap;
VARARRAY	vaItper;
VARARRAY	vaJtper;
VARARRAY	vaKtper;
VARARRAY	vaIctper;
VARARRAY	vaIpper;
VARARRAY	vaJpper;
VARARRAY	vaKpper;
VARARRAY	vaLpper;
VARARRAY	vaIcpper;
VARARRAY	vaLabper;
VARARRAY	vaIgrper;
VARARRAY	vaIsmper;
VARARRAY	vaAlmper;
VARARRAY	vaIaper;
VARARRAY	vaIacper;
VARARRAY	vaCgper;
int		natcap;


#define GETNAME(s,vaS,vaI,i) strcpy(s,PVAI(vaS,LABELTYPE,\
                                abs(*PVAI(vaI,int,i))/3));




/*
 *      sGetName
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the residue number and name of the atom.
 */
char*   sGetName( iIndexIndex, vaIndex, sStr )
int             iIndexIndex;
VARARRAY        vaIndex;
char*           sStr;
{
STRING          sTemp, sNumber, sType;
int             iRes, iAtom;
int		iFirst;
STRING		sLabel;
STRING		sName;

    strcpy( sName, "" );

    if ( vaIndex == NULL ) iAtom = iIndexIndex;
    else iAtom = abs(*PVAI(vaIndex,int,iIndexIndex)/3);
    for ( iRes = 1; iRes<ntotrs; iRes++ ) {
	iFirst = *PVAI( vaIpres, int, iRes );
        if ( iFirst >(iAtom+1) ) break;
    }

    if ( SbIncludeResidueName ) {
	sprintf( sTemp, "%2d(%3s)", iRes, PVAI(vaLabres,LABELTYPE,iRes-1 ));
	strcat( sName, sTemp );
    }
    if ( SbIncludeAtomName ) {
	strcat( sName, PVAI(vaIgraph,LABELTYPE,iAtom) );
    }
    if ( SbIncludeAtomNumberInName ) {
	sprintf( sTemp, "[%3d", iAtom+1 );
	strcat( sName, sTemp );
     }
    if ( SbIncludeAtomTypeInName ) {
        sprintf( sType, "(%s", PVAI( vaIsymbl, LABELTYPE, iAtom ) );
	strcat( sName, sType );
    }

    sRemoveSpaces( sName, sStr );
    return(sStr);
}





/*
 *      FillInt
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the VARARRAY with integer values.
 */
void    FillInt( vaPArray, iElements )
VARARRAY*       vaPArray;
int             iElements;
{
int*    	iPCur;
int     	i;
STRING		sTemp;


    *vaPArray = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( (*vaPArray), iElements );
    if ( !iElements ) {
		FortranSkipLine();
		return;
    }
    iPCur = PVAI( *vaPArray, int, 0 );
    for ( i=0; i<iElements; i++ ) {
        *iPCur = iFortranReadInt();
	if ( SbDescLoad ) {
	    fprintf( stderr, "INT:  %d\n", *iPCur );
	}
        iPCur++;
    }
}




/*
 *      FillDbl
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the VARARRAY with double values.
 */
void    FillDbl( vaPArray, iElements )
VARARRAY*       vaPArray;
int             iElements;
{
double*         dPCur;
int             i;
STRING		sTemp;


    *vaPArray = vaVarArrayCreate( sizeof(double) );
    VarArraySetSize( (*vaPArray), iElements );
    if ( !iElements ) {
		FortranSkipLine();
		return;
    }
    dPCur = PVAI( *vaPArray, double, 0 );
    for ( i=0; i<iElements; i++ ) {
        *dPCur = dFortranReadDouble();
	if ( SbDescLoad ) {
	    fprintf( stderr, "DBL: %lf\n", *dPCur );
	}
        dPCur++;
    }
}


/*
 *      FillLabel
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the VARARRAY with labels.
 */
void    FillLabel( vaPArray, iElements )
VARARRAY*       vaPArray;
int             iElements;
{
char*           cPCur;
int             i;
LABELTYPE       lTemp;


    MESSAGE(( "FillLabel #elements=%d\n", iElements ));
    *vaPArray = vaVarArrayCreate( sizeof(LABELTYPE) );
    VarArraySetSize( (*vaPArray), iElements );
    if ( !iElements ) {
		FortranSkipLine();
		return;
    }
    cPCur = PVAI( *vaPArray, char, 0 );
    for ( i=0; i<iElements; i++ ) {
        strcpy( cPCur, sFortranReadLabel( lTemp ) );
	if ( SbDescLoad ) {
	    fprintf( stderr, "Label: %s\n", cPCur );
	}
        cPCur += sizeof(LABELTYPE);
    }
}





/*
 *      PrintInteraction
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Print the constants for the interaction between the two atoms
 *      and the type of interaction.
 */
void    PrintInteraction( iA, iB, iTotal, vaIndex, vaCn1, vaCn2, vaAg, vaBg )
int             iA;
int             iB;
int             iTotal;
VARARRAY        vaIndex;
VARARRAY        vaCn1;
VARARRAY        vaCn2;
VARARRAY        vaAg;
VARARRAY        vaBg;
{
int             iIndexIndex;
int             iIndex;


    iIndexIndex = (iA-1)*iTotal+(iB-1);
    iIndex = *PVAI( vaIndex, int, iIndexIndex );
    if ( iIndex > 0 ) {
        printf( "NonBond   CN1= %9.3lE   CN2= %9.3lE",
                        *PVAI( vaCn1, double, iIndex-1 ),
                        *PVAI( vaCn2, double, iIndex-1 ) );
    } else {
        printf( "HBond  A= %9.3lE  B= %9.3lE",
                        *PVAI( vaAg, double, (-iIndex-1) ),
                        *PVAI( vaBg, double, (-iIndex-1) ) );
    }
}

/*
 *      AddBondInteraction
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a bond interaction to the DICTIONARY
 */
void    AddBondInteraction( dBonds, i, vaIgraph, vaIb, vaJb, vaIcb, bRest )
DICTIONARY      dBonds;
int		i;
VARARRAY        vaIgraph;
VARARRAY        vaIb;
VARARRAY        vaJb;
VARARRAY        vaIcb;
BOOL		bRest;
{
STRING          sI, sJ, sTemp;
STRING          sName;
BLOCKt*		bPBlock;


    sGetName( i, vaIb, sI );
    sGetName( i, vaJb, sJ );

                /* Order the atoms so that they are in ascending order */
    if ( strcmp( sI, sJ ) > 0 ) {
        strcpy( sTemp, sI ); strcpy( sI, sJ ); strcpy( sJ, sTemp );
    }
    if ( bRest ) strcpy( sName, "RESTRAIN " );
    else         strcpy( sName, "Normal   " );
    strcat( sName, sI ); strcat( sName, "-" );
    strcat( sName, sJ ); 

    MALLOC( bPBlock, BLOCKt*, sizeof(BLOCKt) );
    bPBlock->iIntI = *PVAI(vaIb,int,i);
    bPBlock->iIntJ = *PVAI(vaJb,int,i);
    bPBlock->PData = (GENP)PVAI(vaIcb,int,i);
    DictionaryAdd( dBonds, sName, (GENP)bPBlock ); 
}



/*
 *      AddAngleInteraction
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add an angle interaction to the DICTIONARY
 */
void    AddAngleInteraction( dAngles, i, vaIgraph, vaIt, vaJt, vaKt, vaIct,
				bRest )
DICTIONARY      dAngles;
int		i;
VARARRAY        vaIgraph;
VARARRAY        vaIt;
VARARRAY        vaJt;
VARARRAY        vaKt;
VARARRAY        vaIct;
BOOL		bRest;
{
STRING          sI, sJ, sK, sTemp;
STRING          sName;


    sGetName( i, vaIt, sI );
    sGetName( i, vaJt, sJ );
    sGetName( i, vaKt, sK );

                /* Order the atoms so that they are in ascending order */
    if ( strcmp( sI, sK ) > 0 ) {
        strcpy( sTemp, sI ); strcpy( sI, sK ); strcpy( sK, sTemp );
    }
    if ( bRest ) strcpy( sName, "RESTRAIN " );
    else         strcpy( sName, "Normal   " );
    strcat( sName, sI ); strcat( sName, "-" );
    strcat( sName, sJ ); strcat( sName, "-" );
    strcat( sName, sK ); 

    DictionaryAdd( dAngles, sName, PVAI( vaIct, int, i ) ); 
}



/*
 *      AddTorsion14Interaction
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a torsion interaction to the dictionary and
 *      if there is a 1-4 add it to the 14 dictionary.
 */
void    AddTorsion14Interaction( dTorsions, d14s, i, vaIgraph,
                                vaIp, vaJp, vaKp, vaLp, vaIcp, bRest )
DICTIONARY      dTorsions;
DICTIONARY      d14s;
int             i;
VARARRAY        vaIgraph;
VARARRAY        vaIp;
VARARRAY        vaJp;
VARARRAY        vaKp;
VARARRAY        vaLp;
VARARRAY        vaIcp;
BOOL		bRest;
{
BOOL            bImproper, bCalc14;
STRING          sI, sJ, sK, sL, sTemp, sName;
STRING          s1, s2, s3;
BLOCKt*		bPBlock;


        sGetName( i, vaIp, sI );
        sGetName( i, vaJp, sJ );
        bCalc14 = (*PVAI(vaKp,int,i)>=0);
        sGetName( i, vaKp, sK );
        bImproper = (*PVAI(vaLp,int,i)<0);
        sGetName( i, vaLp, sL );

        if ( bImproper ) bCalc14 = FALSE;

                /* Order the four atoms so that the outer ones are */
                /* in ascending order, if it is a PROPER torsion */
                /* Otherwise sort the first, second and fourth atoms */
                /* do a quick and dirty sort */

	if ( SbSortTorsions ) {
	    if ( !bImproper ) {
		if ( strcmp( sI, sL ) > 0 ) {
		    strcpy( sTemp, sI ); strcpy( sI, sL ); strcpy( sL, sTemp );
		    strcpy( sTemp, sJ ); strcpy( sJ, sK ); strcpy( sK, sTemp );
		}
	    } else {
		    /* If sI is the smallest */
		if ( strcmp( sI, sJ )<0 && strcmp( sI, sL )<0 ) {
		    strcpy( s1, sI );
		    if ( strcmp( sJ, sL )< 0 ) {
			strcpy( s2, sJ );
			strcpy( s3, sL );
		    } else {
			strcpy( s2, sL );
			strcpy( s3, sJ );
		    }
		    /* If sJ is the smallest */
		} else if ( strcmp(sJ,sI)<0 && strcmp(sJ,sL)<0 ) {
		    strcpy( s1, sJ );
		    if ( strcmp( sI, sL )< 0 ) {
			strcpy( s2, sI );
			strcpy( s3, sL );
		    } else {
			strcpy( s2, sL );
			strcpy( s3, sI );
		    }
		    /* Otherwise sL is the smallest */
		} else {
		    strcpy( s1, sL );
		    if ( strcmp( sI, sJ )< 0 ) {
			strcpy( s2, sI );
			strcpy( s3, sJ );
		    } else {
			strcpy( s2, sJ );
			strcpy( s3, sI );
		    }
		}
		strcpy( sI, s1 );
		strcpy( sJ, s2 );
		strcpy( sL, s3 );
	    }       
	}

        if ( bRest ) strcpy( sName, "RSTR " );
        else         strcpy( sName, "NRML " );
        if ( bImproper) strcat( sName, "I " );
        else            strcat( sName, "P " );
        strcat( sName, sI ); strcat( sName, "-" );
        strcat( sName, sJ ); strcat( sName, "-" );
        strcat( sName, sK ); strcat( sName, "-" );
        strcat( sName, sL );

	MALLOC( bPBlock, BLOCKt*, sizeof(BLOCKt) );    
	bPBlock->iIntI = *PVAI(vaIp,int,i);
	bPBlock->iIntJ = *PVAI(vaJp,int,i);
	bPBlock->iIntK = *PVAI(vaKp,int,i);
	bPBlock->iIntL = *PVAI(vaLp,int,i);
	bPBlock->PData = (GENP)PVAI(vaIcp,int,i);
        DictionaryAdd( dTorsions, sName, bPBlock ); 

        if ( bCalc14 ) {
            strcpy( sName, sI );
            strcat( sName, "-" );
            strcat( sName, sL );

                /* The pointer to the data does not mean anything */
                /* in this case, it just has to be non NULL */

            DictionaryAdd( d14s, sName, PVAI( vaIcp, int, i ) ); 
        }
}









/*
 *      AddAllPertBondInteractions
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a perturbed bond interaction to the DICTIONARY
 */
void    AddAllPertBondInteractions( dPBonds )
DICTIONARY*      dPBonds;
{
STRING          sI, sJ, sTemp;
STRING          sName;
BLOCKt*		bPBlock;
int		i;
double		dReq, dRk;
char*		cPNew;
DICTIONARY	dBonds;
int		iParm;


    dBonds = dDictionaryCreate();

    for ( i=0; i<nbper; i++ ) {
	
        sGetName( i, vaIbper, sI );
        sGetName( i, vaJbper, sJ );

                /* Order the atoms so that they are in ascending order */
        if ( strcmp( sI, sJ ) > 0 ) {
            strcpy( sTemp, sI ); strcpy( sI, sJ ); strcpy( sJ, sTemp );
        }
        if ( i < (mbper) ) {
	    strcpy( sName, "PERT " );
        } else {
	    strcpy( sName, "EDGE " );
        }
        strcat( sName, sI );
        strcat( sName, "-" );
        strcat( sName, sJ );
        strcat( sName, " " );
        iParm = *PVAI( vaIcbper, int, i );
        dRk = *PVAI( vaRk, double, iParm-1 );
        dReq = *PVAI( vaReq, double, iParm-1 );
        sprintf( sTemp, "%9.3lE %9.3lE >< ", dRk, dReq );
        strcat( sName, sTemp );

        iParm = *PVAI( vaIcbper, int, i+nbper );
        dRk = *PVAI( vaRk, double, iParm-1 );
        dReq = *PVAI( vaReq, double, iParm-1 );
        sprintf( sTemp, "%9.3lE %9.3lE", dRk, dReq );
        strcat( sName, sTemp );

        DictionaryAdd( dBonds, sName, sName ); 
    }
    *dPBonds = dBonds;
}










/*
 *      AddAllPertAngleInteractions
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add all perturbed angle interactions to the DICTIONARY
 */
void    AddAllPertAngleInteractions( dPAngles )
DICTIONARY*      dPAngles;
{
STRING          sI, sJ, sK, sTemp;
STRING          sName;
BLOCKt*		bPBlock;
int		i;
double		dReq, dRk;
char*		cPNew;
DICTIONARY	dAngles;
int		iParm;
double		dTk, dTeq;


    dAngles = dDictionaryCreate();

    for ( i=0; i<ngper; i++ ) {
	
        sGetName( i, vaItper, sI );
        sGetName( i, vaJtper, sJ );
        sGetName( i, vaKtper, sK );

                /* Order the atoms so that they are in ascending order */
        if ( strcmp( sI, sK ) > 0 ) {
            strcpy( sTemp, sI ); strcpy( sI, sK ); strcpy( sK, sTemp );
        }
        if ( i < (mgper) ) {
	    strcpy( sName, "PERT " );
        } else {
	    strcpy( sName, "EDGE " );
        }
        strcat( sName, sI );
        strcat( sName, "-" );
        strcat( sName, sJ );
        strcat( sName, "-" );
        strcat( sName, sK );
        strcat( sName, " " );
        iParm = *PVAI( vaIctper, int, i );
        dTk = *PVAI( vaTk, double, iParm-1 ),
        dTeq= *PVAI( vaTeq, double, iParm-1 ); 
        sprintf( sTemp, "%9.3lE %9.3lE >< ", dTk, dTeq );
        strcat( sName, sTemp );

        iParm = *PVAI( vaIctper, int, i+ngper );
        dTk = *PVAI( vaTk, double, iParm-1 );
        dTeq= *PVAI( vaTeq, double, iParm-1 ); 
        sprintf( sTemp, "%9.3lE %9.3lE >< ", dTk, dTeq );
        strcat( sName, sTemp );

        DictionaryAdd( dAngles, sName, sName ); 
    }
    *dPAngles = dAngles;
}





/*
 *      AddAllPertTorsionInteractions
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add all perturbed Torsion interactions to the DICTIONARY
 */
void    AddAllPertTorsionInteractions( dPTorsions )
DICTIONARY*      dPTorsions;
{
STRING          sI, sJ, sK, sL, sTemp;
STRING          sName;
BLOCKt*		bPBlock;
int		i;
double		dReq, dRk;
char*		cPNew;
DICTIONARY	dTorsions;
int		iParm;
double		dKp, dN, dPhase;
BOOL		bImproper;
int		iN;


    dTorsions = dDictionaryCreate();

    for ( i=0; i<ndper; i++ ) {

        bImproper = (*PVAI(vaLpper,int,i) < 0 );
        sGetName( i, vaIpper, sI );
        sGetName( i, vaJpper, sJ );
        sGetName( i, vaKpper, sK );
        sGetName( i, vaLpper, sL );

                /* Order the atoms so that they are in ascending order */
        if ( strcmp( sI, sL ) > 0 ) {
            strcpy( sTemp, sI ); strcpy( sI, sL ); strcpy( sL, sTemp );
            strcpy( sTemp, sJ ); strcpy( sJ, sK ); strcpy( sK, sTemp );
        }
        if ( i < (mdper) ) {
	    strcpy( sName, "Pe" );
        } else {
	    strcpy( sName, "Ed" );
        }
	if ( bImproper ) {
	    strcat( sName, "Im" );
	} else {
	    strcat( sName, "Pr" );
        }
        strcat( sName, sI );
        strcat( sName, "-" );
        strcat( sName, sJ );
        strcat( sName, "-" );
        strcat( sName, sK );
        strcat( sName, "-" );
        strcat( sName, sL );
        strcat( sName, " " );
        iParm = *PVAI( vaIcpper, int, i );
	if ( iParm == 0 ) {
	    strcat( sName, "NONE|" );
	} else {
	    dKp = *PVAI(vaPk,double,iParm-1);
	    dN  = fabs(*PVAI(vaPn,double,iParm-1));
	    dPhase = *PVAI(vaPhase,double,iParm-1);
	    iN = (int)dN;
	    sprintf( sTemp, "%9.3lE %d", dKp, iN );
	    strcat( sName, sTemp );
	    if ( SbIncludeIndices ) {
		sprintf( sTemp, "#%d", iParm );
		strcat( sName, sTemp );
	    }
	    strcat( sName, "|" );
	}
        iParm = *PVAI( vaIcpper, int, i+ndper );
	if ( iParm == 0 ) {
	    strcat( sName, "NONE" );
        } else {
	    dKp = *PVAI(vaPk,double,iParm-1);
	    dN  = fabs(*PVAI(vaPn,double,iParm-1));
	    dPhase = *PVAI(vaPhase,double,iParm-1);
	    iN = (int)dN;
	    sprintf( sTemp, "%9.3lE %d", dKp, iN );
	    strcat( sName, sTemp );
	    if ( SbIncludeIndices ) {
		sprintf( sTemp, "#%d", iParm );
		strcat( sName, sTemp );
	    }
	}

        DictionaryAdd( dTorsions, sName, sName ); 
    }
    *dPTorsions = dTorsions;
}






void main( argc, argv )
int             argc;
char*           argv[];
{
/*types*/
int             iIndex;
FILE*           fIn;
STRING          sTitle;
VARARRAY        vaTorsions;
int             i, j;
DICTLOOP        dlLoop, dlLoop2;
VARARRAY        vaAtoms;
DICTIONARY      dTypes;
DICTIONARY      dBonds;
DICTIONARY      dAngles;
DICTIONARY      dTorsions;
DICTIONARY      d14s;
DICTIONARY	dNonBonds;
DICTIONARY	dPertBonds;
STRING          sLine, sEntry, sTemp, sFileName;
STRING		sTemp1, sTemp2;
int		iExclude, iExcludeNum;
int		iTemp;
int		iLeft, iRight;
int		iI, iJ, iK, iL;
VECTOR		vI, vJ, vK, vL;
double		dAngle, dN, dKp, dPhase, dEnergy;
STRING		sCoords;
FILE*		fCrd;
char*		cPStart;
BLOCKt*		bPBlock;
double		dTotalEnergy, dTotalBondEnergy;
DICTIONARY	dPertAngles, dPertTorsions;
int		cc;
extern  int     optind;
extern  char    *optarg;



    BasicsInitialize();

    fprintf( stderr, "AMBER topology file reader\n" );
    if ( argc < 3 ) {
	fprintf(stderr, "need file names\n");
	usage(argv[0]);
    }

    SbSortTorsions = FALSE;
    while ((cc = getopt(argc, argv, "sDnartIei")) != EOF) {
	switch( cc ) {
	    case 's':
		SbSortTorsions = TRUE;
		break;
	    case 'r':
		SbIncludeResidueName = TRUE;
		break;
	    case 'I':
		SbIncludeIndices = TRUE;
		break;
	    case 'a':
		SbIncludeAtomName = TRUE;
		break;
            case 'n': 
                SbIncludeAtomNumberInName = TRUE;
                break;
            case 't':
                SbIncludeAtomTypeInName = TRUE;
                break;
	    case 'i':
		SbIgnoreNonbonds = TRUE;
		break;
	    case 'e':
		SbPrintEnergy = TRUE;
		break;
	    case 'D':
		SbDescLoad = TRUE;
		break;
            default:
		fprintf(stderr, "unprogrammed option: [%c] %d\n", cc, cc);
                usage(argv[0]);
        }
    }
    if (argc - optind != 2) {
	fprintf(stderr, "incorrect number of args %d %d\n", argc, optind);
	usage(argv[0]);
    }
    strcpy( sFileName, argv[argc-2] );
    strcpy( sCoords, argv[argc-1] );

    fIn = FOPENCOMPLAIN( sFileName, "r" );
    if ( fIn == NULL ) {
        printf( "Could not open file: %s\n", sFileName );
        exit(1);
    }

 
    FortranFile( fIn );

    FortranSkipLine();

    fprintf( stderr, "Reading topology file\n" );

    ntotat      = iFortranReadInt();
    ntypes      = iFortranReadInt();
    nbonh       = iFortranReadInt();
    nbona       = iFortranReadInt();
    ntheth      = iFortranReadInt();
    ntheta      = iFortranReadInt();
    nphih       = iFortranReadInt();
    nphia       = iFortranReadInt();
    jhparm      = iFortranReadInt();
    jparm       = iFortranReadInt();
    next        = iFortranReadInt();
    ntotrs      = iFortranReadInt();
    mbona       = iFortranReadInt();
    mtheta      = iFortranReadInt();
    mphia       = iFortranReadInt();
    mumbnd      = iFortranReadInt();
    mumang      = iFortranReadInt();
    mptra       = iFortranReadInt();
    natyp       = iFortranReadInt();
    nhb         = iFortranReadInt();
    ifpert      = iFortranReadInt();
    nbper       = iFortranReadInt();
    ngper       = iFortranReadInt();
    ndper       = iFortranReadInt();
    mbper       = iFortranReadInt();
    mgper       = iFortranReadInt();
    mdper       = iFortranReadInt();
    ifbox       = iFortranReadInt();
    nmxrs       = iFortranReadInt();
    ifcap       = iFortranReadInt();

        /* -3- */
    if ( SbDescLoad ) fprintf( stderr, "-3-\n" );
    FillLabel( &vaIgraph, ntotat );

        /* -4- */
    if ( SbDescLoad ) fprintf( stderr, "-4-\n" );
    FillDbl( &vaChrg, ntotat );

        /* -5- */
    if ( SbDescLoad ) fprintf( stderr, "-5-\n" );
    FillDbl( &vaAmass, ntotat );

        /* -6- */
    if ( SbDescLoad ) fprintf( stderr, "-6-\n" );
    FillInt( &vaIac, ntotat );

        /* -7- */
    if ( SbDescLoad ) fprintf( stderr, "-7-\n" );
    FillInt( &vaNumex, ntotat );

        /* -8- */
    if ( SbDescLoad ) fprintf( stderr, "-8-\n" );
    FillInt( &vaNno, ntypes*ntypes );

        /* -9- */
    if ( SbDescLoad ) fprintf( stderr, "-9-\n" );
    FillLabel( &vaLabres, ntotrs );

        /* -10- */
    if ( SbDescLoad ) fprintf( stderr, "-10-\n" );
    FillInt( &vaIpres, ntotrs );

        /* -11- */
    if ( SbDescLoad ) fprintf( stderr, "-11-\n" );
    FillDbl( &vaRk, mumbnd );

        /* -12- */
    if ( SbDescLoad ) fprintf( stderr, "-12-\n" );
    FillDbl( &vaReq, mumbnd );

        /* -13- */
    if ( SbDescLoad ) fprintf( stderr, "-13-\n" );
    FillDbl( &vaTk, mumang );

        /* -14- */
    if ( SbDescLoad ) fprintf( stderr, "-14-\n" );
    FillDbl( &vaTeq, mumang );

        /* -15- */
    if ( SbDescLoad ) fprintf( stderr, "-15-\n" );
    FillDbl( &vaPk, mptra );

        /* -16- */
    if ( SbDescLoad ) fprintf( stderr, "-16-\n" );
    FillDbl( &vaPn, mptra );

        /* -17- */
    if ( SbDescLoad ) fprintf( stderr, "-17-\n" );
    FillDbl( &vaPhase, mptra );

        /* -18- */
    if ( SbDescLoad ) fprintf( stderr, "-18-\n" );
    FillDbl( &vaSolty, natyp );

        /* -19- */
    if ( SbDescLoad ) fprintf( stderr, "-19-\n" );
    FillDbl( &vaCn1, ntypes*(ntypes+1)/2 );

        /* -20- */
    if ( SbDescLoad ) fprintf( stderr, "-20-\n" );
    FillDbl( &vaCn2, ntypes*(ntypes+1)/2 );

        /* -21- */

    if ( SbDescLoad ) fprintf( stderr, "-21-\n" );
    vaIbh = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIbh, nbonh );
    vaJbh = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaJbh, nbonh );
    vaIcbh = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIcbh, nbonh );
    for ( i=0; i<nbonh; i++ ) {
        *PVAI( vaIbh, int, i ) = iFortranReadInt();
        *PVAI( vaJbh, int, i ) = iFortranReadInt();
        *PVAI( vaIcbh, int, i ) = iFortranReadInt();
	if ( SbDescLoad ) {
	    fprintf( stderr, "%d %d %d\n",
        	*PVAI( vaIbh, int, i ),
        	*PVAI( vaJbh, int, i ),
        	*PVAI( vaIcbh, int, i ) );
	}
    }

        /* -22- */

    if ( SbDescLoad ) fprintf( stderr, "-22-\n" );
    vaIb = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIb, mbona );
    vaJb = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaJb, mbona );
    vaIcb = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIcb, mbona );
    for ( i=0; i<mbona; i++ ) {
        *PVAI( vaIb, int, i ) = iFortranReadInt();
        *PVAI( vaJb, int, i ) = iFortranReadInt();
        *PVAI( vaIcb, int, i ) = iFortranReadInt();
	if ( SbDescLoad ) {
	    fprintf( stderr, "%d %d %d\n",
        	*PVAI( vaIb, int, i ),
        	*PVAI( vaJb, int, i ),
        	*PVAI( vaIcb, int, i ) );
	}
    }


        /* -23- */

    if ( SbDescLoad ) fprintf( stderr, "-23-\n" );
    vaIth = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIth, ntheth );
    vaJth = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaJth, ntheth );
    vaKth = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaKth, ntheth );
    vaIcth = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIcth, ntheth );
    for ( i=0; i<ntheth; i++ ) {
        *PVAI( vaIth, int, i ) = iFortranReadInt();
        *PVAI( vaJth, int, i ) = iFortranReadInt();
        *PVAI( vaKth, int, i ) = iFortranReadInt();
        *PVAI( vaIcth, int, i ) = iFortranReadInt();
	if ( SbDescLoad ) {
	    fprintf( stderr, "%d %d %d %d\n",
        	*PVAI( vaIth, int, i ),
        	*PVAI( vaJth, int, i ),
        	*PVAI( vaKth, int, i ),
        	*PVAI( vaIcth, int, i ) );
	}
    }

        /* -24- */

    if ( SbDescLoad ) fprintf( stderr, "-24-\n" );
    vaIt = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIt, mtheta );
    vaJt = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaJt, mtheta );
    vaKt = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaKt, mtheta );
    vaIct = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIct, mtheta );
    for ( i=0; i<mtheta; i++ ) {
        *PVAI( vaIt, int, i ) = iFortranReadInt();
        *PVAI( vaJt, int, i ) = iFortranReadInt();
        *PVAI( vaKt, int, i ) = iFortranReadInt();
        *PVAI( vaIct, int, i ) = iFortranReadInt();
	if ( SbDescLoad ) {
	    fprintf( stderr, "%d %d %d %d\n",
        	*PVAI( vaIt, int, i ),
        	*PVAI( vaJt, int, i ),
        	*PVAI( vaKt, int, i ),
        	*PVAI( vaIct, int, i ) );
	}
    }

        /* -25- Torsions with Hydrogens */

    if ( SbDescLoad ) fprintf( stderr, "-25-\n" );
    if ( SbDescLoad ) {
	fprintf( stderr, "Torsions with hydrogens\n" );
    }
    vaIph = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIph, nphih );
    vaJph = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaJph, nphih );
    vaKph = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaKph, nphih );
    vaLph = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaLph, nphih );
    vaIcph = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIcph, nphih );
    for ( i=0; i<nphih; i++ ) {
        *PVAI( vaIph, int, i ) = iFortranReadInt();
        *PVAI( vaJph, int, i ) = iFortranReadInt();
        *PVAI( vaKph, int, i ) = iFortranReadInt();
        *PVAI( vaLph, int, i ) = iFortranReadInt();
        *PVAI( vaIcph, int, i ) = iFortranReadInt();
	if ( SbDescLoad ) {
	    fprintf( stderr, "%d %d %d %d %d\n",
        	*PVAI( vaIph, int, i ),
        	*PVAI( vaJph, int, i ),
        	*PVAI( vaKph, int, i ),
        	*PVAI( vaLph, int, i ),
        	*PVAI( vaIcph, int, i ) );
	}
    }

        /* -26- Torsions without hydrogens */

    if ( SbDescLoad ) fprintf( stderr, "-26-\n" );
    if ( SbDescLoad ) {
	fprintf( stderr, "Torsions without hydrogens\n" );
    }
    vaIp = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIp, mphia );
    vaJp = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaJp, mphia );
    vaKp = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaKp, mphia );
    vaLp = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaLp, mphia );
    vaIcp = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaIcp, mphia );
    for ( i=0; i<mphia; i++ ) {
        *PVAI( vaIp, int, i ) = iFortranReadInt();
        *PVAI( vaJp, int, i ) = iFortranReadInt();
        *PVAI( vaKp, int, i ) = iFortranReadInt();
        *PVAI( vaLp, int, i ) = iFortranReadInt();
        *PVAI( vaIcp, int, i ) = iFortranReadInt();
	if ( SbDescLoad ) {
	    fprintf( stderr, "%d %d %d %d %d\n",
        	*PVAI( vaIp, int, i ),
        	*PVAI( vaJp, int, i ),
        	*PVAI( vaKp, int, i ),
        	*PVAI( vaLp, int, i ),
        	*PVAI( vaIcp, int, i ) );
	}
    }


        /* -27- */

    if ( SbDescLoad ) fprintf( stderr, "-27-\n" );
    FillInt( &vaNatex, next );

        /* -28- */

    if ( SbDescLoad ) fprintf( stderr, "-28-\n" );
    FillDbl( &vaAg, nhb );

        /* -29- */

    if ( SbDescLoad ) fprintf( stderr, "-29- vaBg\n" );
    FillDbl( &vaBg, nhb );

        /* -30- */
    if ( SbDescLoad ) fprintf( stderr, "-30- vaHbcut\n" );
    FillDbl( &vaHbcut, nhb );

        /* -31- */
    if ( SbDescLoad ) fprintf( stderr, "-31- vaIsymbl\n" );
    FillLabel( &vaIsymbl, ntotat );

        /* -32- */
    if ( SbDescLoad ) fprintf( stderr, "-32- vaItree\n" );
    FillLabel( &vaItree, ntotat );

        /* -33- */
    if ( SbDescLoad ) fprintf( stderr, "-33-\n" );
    FillInt( &vaJoin, ntotat );

        /* -34- */
    if ( SbDescLoad ) fprintf( stderr, "-34-\n" );
    FillInt( &vaIrotat, ntotat );

	/* -35A- */
    if ( ifbox == 0 ) goto OVERBOX;

    if ( SbDescLoad ) fprintf( stderr, "-35A-\n" );
    iptres = iFortranReadInt();
    nspm   = iFortranReadInt();
    nspsol = iFortranReadInt();
    if ( SbDescLoad ) {
	fprintf( stderr, "INT: %d\n", iptres );
	fprintf( stderr, "INT: %d\n", nspm );
	fprintf( stderr, "INT: %d\n", nspsol );
    }

	/* -35B- */
    if ( SbDescLoad ) fprintf( stderr, "-35B-\n" );
    FillInt( &vaNsp, nspm );

	/* -35C- */

    if ( SbDescLoad ) fprintf( stderr, "-35C-\n" );
    beta = dFortranReadDouble();
    box1 = dFortranReadDouble();
    box2 = dFortranReadDouble();
    box3 = dFortranReadDouble();
    if ( SbDescLoad ) {
	fprintf( stderr, "DBL: %lf\n", beta );
	fprintf( stderr, "DBL: %lf\n", box1 );
	fprintf( stderr, "DBL: %lf\n", box2 );
	fprintf( stderr, "DBL: %lf\n", box3 );
    }
   

OVERBOX:


    if ( ifcap == 0 ) goto OVERCAP;

	/* -35D- */

    if ( SbDescLoad ) fprintf( stderr, "-35D-\n" );
    natcap = iFortranReadInt();

	/* -35E- */

    if ( SbDescLoad ) fprintf( stderr, "-35E-\n" );
    cutcap = dFortranReadDouble();
    xcap = dFortranReadDouble();
    ycap = dFortranReadDouble();
    zcap = dFortranReadDouble();

OVERCAP:

    if ( ifpert == 0 ) goto OVERPERT;

	/* -36A- */

    if ( SbDescLoad ) fprintf( stderr, "-36A-\n" );
    vaIbper = vaVarArrayCreate(sizeof(int));
    vaJbper = vaVarArrayCreate(sizeof(int));
    for ( i=0; i<nbper; i++ ) {
	iTemp = iFortranReadInt();
	VarArrayAdd( vaIbper, (GENP)&iTemp );
	iTemp = iFortranReadInt();
	VarArrayAdd( vaJbper, (GENP)&iTemp );
    }

	/* -36B- */

    if ( SbDescLoad ) fprintf( stderr, "-36B-\n" );
    FillInt( &vaIcbper, 2*nbper );



	/* -36C- */

    if ( SbDescLoad ) fprintf( stderr, "-36C-\n" );
    vaItper = vaVarArrayCreate(sizeof(int));
    vaJtper = vaVarArrayCreate(sizeof(int));
    vaKtper = vaVarArrayCreate(sizeof(int));
    for ( i=0; i<ngper; i++ ) {
	iTemp = iFortranReadInt();
	VarArrayAdd( vaItper, (GENP)&iTemp );
	iTemp = iFortranReadInt();
	VarArrayAdd( vaJtper, (GENP)&iTemp );
	iTemp = iFortranReadInt();
	VarArrayAdd( vaKtper, (GENP)&iTemp );
    }

	/* -36D- */

    if ( SbDescLoad ) fprintf( stderr, "-36D-\n" );
    FillInt( &vaIctper, 2*ngper );



	/* -36E- */

    if ( SbDescLoad ) fprintf( stderr, "-36E-\n" );
    vaIpper = vaVarArrayCreate(sizeof(int));
    vaJpper = vaVarArrayCreate(sizeof(int));
    vaKpper = vaVarArrayCreate(sizeof(int));
    vaLpper = vaVarArrayCreate(sizeof(int));
    for ( i=0; i<ndper; i++ ) {
	iTemp = iFortranReadInt();
	VarArrayAdd( vaIpper, (GENP)&iTemp );
	iTemp = iFortranReadInt();
	VarArrayAdd( vaJpper, (GENP)&iTemp );
	iTemp = iFortranReadInt();
	VarArrayAdd( vaKpper, (GENP)&iTemp );
	iTemp = iFortranReadInt();
	VarArrayAdd( vaLpper, (GENP)&iTemp );
    }

	/* -36F- */

    if ( SbDescLoad ) fprintf( stderr, "-36F-\n" );
    FillInt( &vaIcpper, 2*ndper );


	/* -36G- */

    if ( SbDescLoad ) fprintf( stderr, "-36G-\n" );
    FillLabel( &vaLabper, ntotrs );

	/* -36H- */

    if ( SbDescLoad ) fprintf( stderr, "-36H-\n" );
    FillLabel( &vaIgrper, ntotat );

	/* -36I- */

    if ( SbDescLoad ) fprintf( stderr, "-36I-\n" );
    FillLabel( &vaIsmper, ntotat );

	/* -36J- */

    if ( SbDescLoad ) fprintf( stderr, "-36J-\n" );
    FillDbl( &vaAlmper, ntotat );

	/* -36K- */

    if ( SbDescLoad ) fprintf( stderr, "-36K-\n" );
    FillInt( &vaIaper, ntotat );

	/* -36L- */

    if ( SbDescLoad ) fprintf( stderr, "-36L-\n" );
    FillInt( &vaIacper, ntotat );

	/* -36M- */

    if ( SbDescLoad ) fprintf( stderr, "-36M-\n" );
    FillDbl( &vaCgper, ntotat );

OVERPERT:



		/* Read the coordinate file */

    fCrd = fopen( sCoords, "r" );
    fgets( sLine, sizeof(STRING), fCrd );
    fgets( sLine, sizeof(STRING), fCrd );

    FortranFile( fCrd );

    if ( SbDescLoad ) fprintf( stderr, "-COORDINATES-\n" );
    FillDbl( &vaCoords, ntotat*3 );


        /*
         *---------------------------------------------------------------
         */


        /* BEGIN processing data */

    fprintf( stderr, "Building atom type information\n" );

        /* Construct a dictionary of atom names, types, and charges */
    vaAtoms = vaVarArrayCreate( sizeof(STRING) );
    for ( i=0; i<ntotat; i++ ) {
        sprintf( sEntry, "%-20s %8.3lf %4s",
                        sGetName( i, NULL, sTemp ),
                        *PVAI( vaChrg, double, i ),
                        PVAI( vaIsymbl, LABELTYPE, i ) );
        VarArrayAdd( vaAtoms, (GENP)sEntry );
    }
    SortByString( PVAI(vaAtoms,STRING,0), ntotat, sizeof(STRING),
                        PVAI(vaAtoms,STRING,0), TRUE );

        /* Construct a dictionary of atom type names to atom type IDs */

    dTypes = dDictionaryCreate();
    for ( i=0; i<ntotat; i++ ) {
        if ( yPDictionaryFind( dTypes, PVAI( vaIsymbl, LABELTYPE, i ))==NULL ) {
            DictionaryAdd( dTypes, PVAI( vaIsymbl, LABELTYPE, i ),
                                PVAI( vaIac, int, i ) );
        }
    }

        /* Construct a dictionary of bond interactions */

    fprintf( stderr, "Building bond information\n" );
    dBonds = dDictionaryCreate();
    for ( i=0; i<nbonh; i++ ) {
        AddBondInteraction( dBonds, i, vaIgraph, vaIbh, vaJbh, vaIcbh, FALSE );
    }
    for ( i=0; i<mbona; i++ ) {
        AddBondInteraction( dBonds, i, vaIgraph, vaIb, vaJb, vaIcb,
		(i>=nbona) );
    }


        /* Construct a dictionary of angle interactions */

    fprintf( stderr, "Building angle information\n" );
    dAngles = dDictionaryCreate();
    for ( i=0; i<ntheth; i++ ) {
        AddAngleInteraction( dAngles, i, vaIgraph, 
                                vaIth, vaJth, vaKth, vaIcth, FALSE );
    }
    for ( i=0; i<mtheta; i++ ) {
        AddAngleInteraction( dAngles, i, vaIgraph, 
                                vaIt, vaJt, vaKt, vaIct, (i>=ntheta) );
    }

        /* Construct a dictionary of torsion interactions */

    fprintf( stderr, "Building torsion and 1-4 information\n" );
    dTorsions = dDictionaryCreate();
    d14s = dDictionaryCreate();
    for ( i=0; i<nphih; i++ ) {
        AddTorsion14Interaction( dTorsions, d14s, i, vaIgraph,
                                vaIph, vaJph, vaKph, vaLph, vaIcph, FALSE );
    }
    for ( i=0; i<mphia; i++ ) {
        AddTorsion14Interaction( dTorsions, d14s, i, vaIgraph,
                                vaIp, vaJp, vaKp, vaLp, vaIcp, (i>=nphia) );
    }


	/* Construct a dictionary of non-bond interactions */

    if ( !SbIgnoreNonbonds ) {
        fprintf( stderr, "Building non-bond information\n" );
	iExclude = 0;
	dNonBonds = dDictionaryCreate();
	for ( i=0; i<ntotat-1; i++ ) {
	    iExcludeNum = *PVAI( vaNumex, int, i );
	    for ( j=i+1; j<ntotat; j++ ) {
	        if ( iExcludeNum ) {
		    iTemp = *PVAI( vaNatex, int, iExclude );
		    if ( j+1 == iTemp ) {
		        iExcludeNum--;
		        iExclude++;
		        continue;
		    }
	        }

	        if ( strcmp( (char*)(PVAI(vaIgraph,LABELTYPE,i)), 
			(char*)(PVAI(vaIgraph,LABELTYPE,j)) ) < 0 ) {
		    iLeft = i;
		    iRight = j;
	        } else {
		    iLeft = j;
		    iRight = i;
	        }
	        sprintf( sEntry, "%s[%s] -- %s[%s]",
			sGetName( iLeft, NULL, sTemp1 ),
			PVAI( vaIsymbl, LABELTYPE, iLeft ),
			sGetName( iRight, NULL, sTemp2 ),
			PVAI( vaIsymbl, LABELTYPE, iRight ) );
	        DictionaryAdd( dNonBonds, sEntry, sEntry );
	    }
	}
    }


    if ( ifpert != 0 ) {
	AddAllPertBondInteractions( &dPertBonds );
	AddAllPertAngleInteractions( &dPertAngles );
	AddAllPertTorsionInteractions( &dPertTorsions );
    }







        /*
         *--------------------------------------------------------------
         */

        /* BEGIN output of data */

    fprintf( stderr, "Writing output\n" );

    	/* Print box information if there is any */

    if ( ifbox != 0 ) {
	printf( "Box\n" );
	printf( "Beta: %9.3lE\n", beta );
	printf( "Box dimensions: %9.3lE, %9.3lE, %9.3lE\n",
		box1, box2, box3 );
	printf( "Last residue not solvent: %d\n", iptres );
	printf( "Number of molecules: %d\n", nspm );
	printf( "First solvent molecule: %d\n", nspsol );
	printf( "Number of atoms per molecule: %d\n" );
	printf( "%-20s %-20s\n", "Molecule", "#Atoms" );
	for ( i=0; i<nspm; i++ ) {
	    printf( "%20d %20d\n", i+1, *PVAI(vaNsp,int,i) );
        }
    } else {
	printf( "NO box.\n" );
    }

    if ( ifcap != 0 ) {
	printf( "Cap\n" );
	printf( "Number of atoms in cap: %d\n", natcap );
	printf( "Cutcap: %9.3lE\n", cutcap );
	printf( "xcap, ycap, zcap = %9.3lE, %9.3lE, %9.3lE\n", xcap, ycap, zcap );
    }

    printf( "=Residue names\n" );
    for ( i=0; i<ntotrs; i++ ) {
        printf( "%5d  %s\n", i+1, PVAI( vaLabres, LABELTYPE, i ) );
    }

    printf( "Res.  Atom     Charges Type\n" );
    for ( i=0; i<ntotat; i++ ) {
        printf( "%s\n", PVAI( vaAtoms, STRING, i ) );
    }

 
    printf( "=Types of Non-bond interactions.\n" );
    dlLoop = ydlDictionaryLoop(dTypes);
    while ( yPDictionaryNext( dTypes, &dlLoop ) ) {
        dlLoop2 = dlLoop;
        do {
            printf( "%5s -%5s: ", sDictLoopKey( dlLoop ),
                                  sDictLoopKey( dlLoop2 ) );
            PrintInteraction( *(int*)PDictLoopData(dlLoop),
                              *(int*)PDictLoopData(dlLoop2),
                              ntypes, vaNno, vaCn1, vaCn2, vaAg, vaBg );
            printf( "\n" );
            yPDictionaryNext( dTypes, &dlLoop2 );
        } while ( dlLoop2 != NULL );
    }


    if ( !SbIgnoreNonbonds ) {
    	printf( "=Non-bond interactions by atom [type].\n" );
    	dlLoop = ydlDictionaryLoop(dNonBonds);
    	while ( yPDictionaryNext( dNonBonds, &dlLoop ) ) {
	    printf( "%s\n", sDictLoopKey(dlLoop) );
	}
    }

  
    dTotalBondEnergy = 0.0; 
    printf( "=Bond interactions.\n" );
    dlLoop = ydlDictionaryLoop(dBonds);
    while ( yPDictionaryNext( dBonds, &dlLoop ) ) {
	bPBlock = (BLOCKt*)PDictLoopData(dlLoop);
	iIndex = (*(int*)bPBlock->PData)-1;
	    iI = bPBlock->iIntI;
	    iJ = bPBlock->iIntJ;
	    VectorDef( &vI, *PVAI(vaCoords,double,iI),
			*PVAI(vaCoords,double,iI+1),
			*PVAI(vaCoords,double,iI+2) );
	    VectorDef( &vJ, *PVAI(vaCoords,double,iJ),
			*PVAI(vaCoords,double,iJ+1),
			*PVAI(vaCoords,double,iJ+2) );

	    dAngle = dVectorAtomLength( &vI, &vJ );
	    dKp = *PVAI(vaRk,double,iIndex);
	    dPhase = *PVAI(vaReq,double,iIndex);
	    dEnergy = dKp*(dAngle - dPhase)*(dAngle-dPhase);
	    dTotalBondEnergy += dEnergy;

        printf( "%s %9.3lE %9.3lE",
                sDictLoopKey(dlLoop),
                *PVAI( vaRk, double, iIndex ),
                *PVAI( vaReq, double, iIndex )); 
	if ( SbPrintEnergy ) printf( " E=%9.3lE\n", dEnergy );
	else printf( "\n" );
    }
    if ( SbPrintEnergy ) {
        printf( "Total bond energy = %9.3lE\n", dTotalBondEnergy );
    }

    printf( "=Angle interactions.\n" );
    dlLoop = ydlDictionaryLoop(dAngles);
    while ( yPDictionaryNext( dAngles, &dlLoop ) ) {
        iIndex = *(int*)PDictLoopData(dlLoop)-1;
        printf( "%s %9.3lE %9.3lE\n",
                sDictLoopKey(dlLoop),
                *PVAI( vaTk, double, iIndex ),
                *PVAI( vaTeq, double, iIndex ) ); 
    }

    printf( "=Torsion interactions.\n" );

    dTotalEnergy = 0.0;
    vaTorsions = vaVarArrayCreate( sizeof(STRING) );
    dlLoop = ydlDictionaryLoop(dTorsions);
    while ( yPDictionaryNext( dTorsions, &dlLoop )) {

	bPBlock = (BLOCKt*)PDictLoopData(dlLoop);
	iIndex = (*(int*)bPBlock->PData)-1;

                /* Keep printing interactions until Pn(iIndex)>0.0 */
                /* If it is less then zero it means that there is */
                /* another interaction for this dihedral at iIndex+1 */

                /* If K<0 then don't calculate 1-4 interactions */
                /* In topology files where the torsions are 'rolled up' */
                /* Don't calculate 1-4 interactions on all subsequent */
                /* torsion interactions */

        do {
	    iI = bPBlock->iIntI;
	    iJ = bPBlock->iIntJ;
	    iK = abs(bPBlock->iIntK);
	    iL = abs(bPBlock->iIntL);
	    VectorDef( &vI, *PVAI(vaCoords,double,iI),
			*PVAI(vaCoords,double,iI+1),
			*PVAI(vaCoords,double,iI+2) );
	    VectorDef( &vJ, *PVAI(vaCoords,double,iJ),
			*PVAI(vaCoords,double,iJ+1),
			*PVAI(vaCoords,double,iJ+2) );
	    VectorDef( &vK, *PVAI(vaCoords,double,iK),
			*PVAI(vaCoords,double,iK+1),
			*PVAI(vaCoords,double,iK+2) );
	    VectorDef( &vL, *PVAI(vaCoords,double,iL),
			*PVAI(vaCoords,double,iL+1),
			*PVAI(vaCoords,double,iL+2) );

	    dAngle = dVectorAtomTorsion( &vI, &vJ, &vK, &vL );
	    dKp = *PVAI(vaPk,double,iIndex);
	    dN  = fabs(*PVAI(vaPn,double,iIndex));
	    dPhase = *PVAI(vaPhase,double,iIndex);
	    dEnergy = dKp*(1+cos(dN*dAngle - dPhase));
	    dTotalEnergy += dEnergy;

	    if ( SbPrintEnergy ) {
		sprintf( sLine, "%s %1.1lf %9.3lE %5.1lf phi= %9.1lf E=%9.2lE",
                	sDictLoopKey(dlLoop),
                	fabs(*PVAI( vaPn, double, iIndex )),
                	*PVAI( vaPk, double, iIndex ),
               		*PVAI( vaPhase, double, iIndex )/DEGTORAD,
			dAngle/DEGTORAD,
			dEnergy );
	    } else {
		sprintf( sLine, "%s %1.1lf %9.3lE %5.1lf",
                	sDictLoopKey(dlLoop),
                	fabs(*PVAI( vaPn, double, iIndex )),
                	*PVAI( vaPk, double, iIndex ),
               		*PVAI( vaPhase, double, iIndex )/DEGTORAD );
	    }

            VarArrayAdd( vaTorsions, (GENP)sLine );
            iIndex++;
        } while ( *PVAI( vaPn, double, iIndex-1 ) < 0.0 );
    }
    printf( "Please wait, sorting . . .\n" );
    SortByString( PVAI( vaTorsions, STRING, 0 ), 
                iVarArrayElementCount(vaTorsions),
                sizeof(STRING), PVAI( vaTorsions, STRING, 0 ), TRUE );
    for ( i=0; i<iVarArrayElementCount(vaTorsions); i++ ) {
        printf( "%s\n", PVAI( vaTorsions, STRING, i ) );
    }
    if ( SbPrintEnergy ) printf( "Total torsion energy = %9.3lE\n", dTotalEnergy );


    printf( "=1-4 non-bond and electrostatic interactions.\n" );
    dlLoop = ydlDictionaryLoop(d14s);
    while ( yPDictionaryNext( d14s, &dlLoop ) ) {
        iIndex = *(int*)PDictLoopData(dlLoop)-1;
        printf( "%s\n", sDictLoopKey(dlLoop) );
    }

    if ( ifpert != 0 ) {
	printf( "Perturbation stuff.\n" );
	printf( "=Perturbed bonds: PERT: %d, EDGE: %d.\n", (mbper), (nbper-mbper) );
	dlLoop = ydlDictionaryLoop( dPertBonds );
	while ( yPDictionaryNext( dPertBonds, &dlLoop ) ) {
	    printf( "%s\n", sDictLoopKey(dlLoop) );
	}
	printf( "=Perturbed angles: PERT: %d, EDGE: %d.\n", (mgper), (ngper-mgper) );
	dlLoop = ydlDictionaryLoop( dPertAngles );
	while ( yPDictionaryNext( dPertAngles, &dlLoop ) ) {
	    printf( "%s\n", sDictLoopKey(dlLoop) );
	}
	printf( "=Perturbed torsions: PERT: %d, EDGE: %d.\n", 
		    (mdper), (ndper-mdper) );
	dlLoop = ydlDictionaryLoop( dPertTorsions );
	while ( yPDictionaryNext( dPertTorsions, &dlLoop ) ) {
	    printf( "%s\n", sDictLoopKey(dlLoop) );
	}
    }

    fclose( fIn );

    

    exit(0);
}
