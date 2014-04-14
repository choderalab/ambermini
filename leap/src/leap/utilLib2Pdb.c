/*
 *	File:	utilLib2Pdb.c
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
 *		A program that quickly generates a pdb
 *		file from a LIBRARY file.
 *		A pdb file is generated for each UNIT in the
 *		LIBRARY.
 */


#include	"basics.h"

#include	"classes.h"

#include	"library.h"

#include	"pdbFile.h"


main( argc, argv )
int		argc;
char*		argv[];
{
STRING		sFile;
LIBRARY		lUnits;
UNIT		uUnit;
STRING		sOutput;
FILE*		fPdb;
char*		cPUnit;


    BasicsInitialize();
    
    if ( argc != 2 ) goto BADARG;
    
    strcpy( sFile, argv[1] );
    lUnits = lLibraryOpen( sFile, OPENREADONLY );
    if ( lUnits == NULL ) {
    	fprintf( stderr, "Could not open file: %s\n", sFile );
	exit(1);
    }
    
    LibraryLoop( lUnits );
    while ( cPUnit = sLibraryNext(lUnits) ) {
    	uUnit = (UNIT)oLibraryLoad( lUnits, cPUnit );
	if ( uUnit == NULL ) {
	    fprintf( stderr, "Could not find UNIT: %s\n", cPUnit );
	    exit(1);
	}
        strcpy( sOutput, sContainerName(uUnit) );
	strcat( sOutput, ".pdb" );
	fPdb = fopen( sOutput, "w" );
	if ( fPdb == NULL ) {
	    fprintf( stderr, "Could not open file: %s\n", sOutput );
	    exit(1);
	}
	fprintf( stderr, "Writing PDB file: %s\n", sOutput );
	PdbWrite( fPdb, uUnit );
	fclose(fPdb);
    }
    LibraryClose( &lUnits );
    exit(0);

BADARG:
    fprintf( stderr, "Usage: %s LIBRARY\n", argv[0] );
    exit(1);
    
}
