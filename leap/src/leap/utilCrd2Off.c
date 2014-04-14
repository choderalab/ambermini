/*
 *	File:	utilCrd2Off.c
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
 *	Description:
 *		Convert an AMBER coordinate file to an OFF
 *		file.
 */



#include	"basics.h"

#include	"vector.h"

#include        "varArray.h"
#include        "dictionary.h"
#include        "fortran.h"

#include	"database.h"

typedef char    LABELTYPE[90];


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
static	BOOL    SbPrintEnergy = FALSE;
static	BOOL	SbIgnoreNonbonds = FALSE;
static	BOOL	SbPrintNumbers = FALSE;

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

VARARRAY	vaIbper, vaJpber;
VARARRAY	vaIcbper;



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


    if ( vaIndex == NULL ) iAtom = iIndexIndex;
    else iAtom = abs(*PVAI(vaIndex,int,iIndexIndex)/3);
    for ( iRes = 1; iRes<ntotrs; iRes++ ) {
	iFirst = *PVAI( vaIpres, int, iRes );
        if ( iFirst >(iAtom+1) ) break;
    }

    strcpy( sTemp, (char*)(PVAI( vaIgraph, LABELTYPE, iAtom ) ));
    if ( SbIncludeAtomNumberInName ) sprintf( sNumber, "[%3d]", iAtom+1 );
    else                             strcpy( sNumber, "" );
    if ( SbIncludeAtomTypeInName ) 
        sprintf( sType, "(%s)", PVAI( vaIsymbl, LABELTYPE, iAtom ) );
    else strcpy( sType, "" );
    sprintf( sStr, "%2d(%3s):%s%s%s", iRes, PVAI(vaLabres,LABELTYPE,iRes-1),
                 sTemp, sType, sNumber );
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
    iPCur = PVAI( *vaPArray, int, 0 );
    for ( i=0; i<iElements; i++ ) {
        *iPCur = iFortranReadInt();
        iPCur++;
    }
    if ( iElements == 0 ) {
	FortranSkipLine();
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
    dPCur = PVAI( *vaPArray, double, 0 );
    for ( i=0; i<iElements; i++ ) {
        *dPCur = dFortranReadDouble();
        dPCur++;
    }
    if ( iElements == 0 ) {
	FortranSkipLine();
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
    cPCur = PVAI( *vaPArray, char, 0 );
    for ( i=0; i<iElements; i++ ) {
        strcpy( cPCur, sFortranReadLabel( lTemp ) );
        cPCur += sizeof(LABELTYPE);
    }
    if ( iElements == 0 ) {
	FortranSkipLine();
    }
}
void main( argc, argv )
int             argc;
char*           argv[];
{
int             iIndex;
FILE*           fIn;
FILE*           fCrd;
int		iCount;
STRING		sName;
STRING          sTitle;
int             i, j;
VARARRAY        vaAtoms;
STRING          sLine, sEntry, sTemp, sFileName, sCoords, sOff;
DATABASE	dbOff;
VARARRAY	vaResidueSequence, vaAtomSequence;
int		iAtom, iResidue;
STRING		sTemp1, sTemp2;


    BasicsInitialize();

    fprintf( stderr, "Convert AMBER crd file to OFF\n" );
    if ( argc == 1 ) {
        fprintf( stderr, "Usage: %s {topologyFile} {crdFile} {OffFile}\n", 
			argv[0] );
        exit(0);
    }

    if ( argc == 4 ) {
        strcpy( sFileName, argv[1] );
        strcpy( sCoords, argv[2] );
        strcpy( sOff, argv[3] );
    } else {
        fprintf( stderr, "Illegal number of arguments.\n" );
        exit(1);
    }

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
    FillLabel( &vaIgraph, ntotat );

        /* -4- */
    FillDbl( &vaChrg, ntotat );

        /* -5- */
    FillDbl( &vaAmass, ntotat );

        /* -6- */
    FillInt( &vaIac, ntotat );

        /* -7- */
    FillInt( &vaNumex, ntotat );

        /* -8- */
    FillInt( &vaNno, ntypes*ntypes );

        /* -9- */
    FillLabel( &vaLabres, ntotrs );

        /* -10- */
    FillInt( &vaIpres, ntotrs );

        /* -11- */
    FillDbl( &vaRk, mumbnd );

        /* -12- */
    FillDbl( &vaReq, mumbnd );

        /* -13- */
    FillDbl( &vaTk, mumang );

        /* -14- */
    FillDbl( &vaTeq, mumang );

        /* -15- */
    FillDbl( &vaPk, mptra );

        /* -16- */
    FillDbl( &vaPn, mptra );

        /* -17- */
    FillDbl( &vaPhase, mptra );

        /* -18- */
    FillDbl( &vaSolty, natyp );

        /* -19- */

    FillDbl( &vaCn1, ntypes*(ntypes+1)/2 );

        /* -20- */

    FillDbl( &vaCn2, ntypes*(ntypes+1)/2 );

        /* -21- */

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
    }

        /* -22- */

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
    }


        /* -23- */

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
    }

        /* -24- */

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
    }

        /* -25- Torsions with Hydrogens */

    if ( SbPrintNumbers ) {
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
	if ( SbPrintNumbers ) {
	    fprintf( stderr, "%d %d %d %d %d\n",
        	*PVAI( vaIph, int, i ),
        	*PVAI( vaJph, int, i ),
        	*PVAI( vaKph, int, i ),
        	*PVAI( vaLph, int, i ),
        	*PVAI( vaIcph, int, i ) );
	}
    }

        /* -26- Torsions without hydrogens */

    if ( SbPrintNumbers ) {
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
	if ( SbPrintNumbers ) {
	    fprintf( stderr, "%d %d %d %d %d\n",
        	*PVAI( vaIp, int, i ),
        	*PVAI( vaJp, int, i ),
        	*PVAI( vaKp, int, i ),
        	*PVAI( vaLp, int, i ),
        	*PVAI( vaIcp, int, i ) );
	}
    }


        /* -27- */

    FillInt( &vaNatex, next );

        /* -28- */

    FillDbl( &vaAg, nhb );

        /* -29- */

    FillDbl( &vaBg, nhb );

        /* -30- */
    FillDbl( &vaHbcut, nhb );

        /* -31- */
    FillLabel( &vaIsymbl, ntotat );

        /* -32- */
    FillLabel( &vaItree, ntotat );

        /* -33- */
    FillInt( &vaJoin, ntotat );

        /* -34- */
    FillInt( &vaIrotat, ntotat );

	/* -35A- */
    if ( ifbox == 0 ) goto OVERBOX;

    iptres = iFortranReadInt();
    nspm   = iFortranReadInt();
    nspsol = iFortranReadInt();

	/* -35B- */
    FillInt( &vaNsp, nspm );

	/* -35C- */

    beta = dFortranReadDouble();
    box1 = dFortranReadDouble();
    box2 = dFortranReadDouble();
    box3 = dFortranReadDouble();

OVERBOX:


		/* Write the table of residue/atom sequence number */

    dbOff = dbDBSeqOpen( sOff, OPENREADWRITE );

    vaResidueSequence = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaResidueSequence, ntotat );
    vaAtomSequence = vaVarArrayCreate( sizeof(int) );
    VarArraySetSize( vaAtomSequence, ntotat );

    iResidue = 1;
    for ( i=1; i<=ntotat; i++ ) {
	if ( iResidue<ntotrs ) {
	    if ( i>=*PVAI(vaIpres,int,iResidue) ) iResidue++;
	}
	iAtom = i-*PVAI(vaIpres,int,iResidue-1)+1;
	*PVAI( vaResidueSequence, int, i-1 ) = iResidue;
	*PVAI( vaAtomSequence, int, i-1 ) = iAtom;
    }

    DBPutTable( dbOff, "crd.atom.names", ntotat,
		1, "resseq", PVAI(vaResidueSequence,int,0), sizeof(int),
		2, "atomseq", PVAI(vaAtomSequence,int,0), sizeof(int),
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
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0 );


		/* Read the coordinate file */

    fCrd = fopen( sCoords, "r" );
    fgets( sLine, sizeof(STRING), fCrd );
    fgets( sLine, sizeof(STRING), fCrd );

    FortranFile( fCrd );

    iCount = 1;
    while ( !feof(fCrd) ) {
	FillDbl( &vaCoords, ntotat*3 );
	if ( feof(fCrd) ) break;
	sprintf( sName, "crd.atoms.crd.%d", iCount ); 
	iCount++;
        DBPutTable( dbOff, sName, ntotat,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,

		1, "x", PVAI(vaCoords,double,0), sizeof(double)*3,
		2, "y", PVAI(vaCoords,double,1), sizeof(double)*3,
		3, "z", PVAI(vaCoords,double,2), sizeof(double)*3,
		0, NULL, NULL, 0,

		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0 );
    }

    DBClose( &dbOff );
    fclose( fIn );

    

    exit(0);
}
