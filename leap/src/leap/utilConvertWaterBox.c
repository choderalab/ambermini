/*
 *	File:	convertWat.c
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
 *		Convert AMBER waterbox to a LEaP UNIT.
 */







#include	"basics.h"

#include	"classes.h"

#include	"dictionary.h"
#include	"database.h"
#include	"library.h"




main( argc, argv )
int		argc;
char*		argv[];
{
FILE*		fWat;
LIBRARY		lWater;
double		daCoords[1000][3][4];
UNIT		uUnit;
RESIDUE		rRes;
ATOM		aH1, aH2, aO;
VECTOR		v;
char		sLine[1000];
int		iMol, i, j, k, iTotal;
double		dCoord;
MOLECULE	mMol;
STRING		sName;
double		dWidth;


    BasicsInitialize();

    printf( "Convert AMBER formatted waterbox to LEaP library.\n" );
    if ( argc != 3 ) {
    	printf( "Usage: %s {waterbox file} {library}\n", argv[0] );
    	exit(1);
    }

    printf( "Reading water box coordinates from: %s into LIBRARY: %s\n",
    			argv[1], argv[2] );

    fWat = FOPENCOMPLAIN( argv[1], "r" );

	/* Skip the first line */

    fgets( sLine, sizeof(sLine), fWat );
    sscanf( sLine, "%d %lf", &iMol, &dWidth );

    iTotal = 0;
    for ( k=0; k<4; k++ ) {
	for ( j=0; j<3; j++ ) {
	    for ( i=0; i<iMol; i++ ) {
	        if ( fscanf( fWat, "%lf", &dCoord ) <= 0 ) goto DONE;
	        daCoords[i][j][k] = dCoord;
	        iTotal++;
	    }
	}
    }

DONE:
    fclose( fWat );

    printf( "Read %d numbers\n", iTotal );

    uUnit = (UNIT)oCreate(UNITid);
    ContainerSetName( uUnit, "WATBOX216" );
    UnitSetUseBox( uUnit, TRUE );
    UnitSetBeta( uUnit, 90.0 );			/* Beta is stored in degrees */
    UnitSetBox( uUnit, dWidth, dWidth, dWidth );

    for ( i=0; i<iMol; i++ ) {
	rRes = (RESIDUE)oCreate(RESIDUEid);
	ContainerSetName( rRes, "WAT" );
	ResidueSetType( rRes, RESTYPESOLVENT );
	aH1 = (ATOM)oCreate(ATOMid);
	ContainerSetName( aH1, "H1" );
	AtomSetType( aH1, "HW" );
	AtomSetElement( aH1, HYDROGEN );
	VectorDef( &v, daCoords[i][0][1], 
			daCoords[i][1][1], 
			daCoords[i][2][1] );
	AtomSetPosition( aH1, v );
	AtomSetCharge( aH1, 0.417 );

	aH2 = (ATOM)oCreate(ATOMid);
	ContainerSetName( aH2, "H2" );
	AtomSetType( aH2, "HW" );
	AtomSetElement( aH2, HYDROGEN );
	VectorDef( &v, daCoords[i][0][2], 
			daCoords[i][1][2], 
			daCoords[i][2][2] );

	AtomSetPosition( aH2, v );
	AtomSetCharge( aH2, 0.417 );

	aO = (ATOM)oCreate(ATOMid);
	ContainerSetName( aO, "O" );
	AtomSetType( aO, "OW" );
	AtomSetElement( aO, OXYGEN );
	VectorDef( &v, daCoords[i][0][0], 
			daCoords[i][1][0], 
			daCoords[i][2][0] );
	AtomSetPosition( aO, v );
	AtomSetCharge( aO, -0.834 );
	
	AtomBondTo( aH1, aO );
	AtomBondTo( aH2, aO );
	AtomBondTo( aH1, aH2 );
	
	ContainerAdd( (CONTAINER)rRes, (OBJEKT)aO );
	ContainerAdd( (CONTAINER)rRes, (OBJEKT)aH1 );
	ContainerAdd( (CONTAINER)rRes, (OBJEKT)aH2 );

		/* Define the main atom of the solvent */
	ResidueSetImagingAtom( rRes, aO );
       
	ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rRes );
    }
    
    lWater = lLibraryOpen( argv[2], OPENREADWRITE );
    LibrarySave( lWater, "WATBOX218", uUnit, NULL );
    LibraryClose( &lWater );
  
    exit(0);
}
