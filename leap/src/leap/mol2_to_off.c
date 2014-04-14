/*
 *	File:	utilTripos.c
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
 *		Read a TRIPOS file and write OFF files.
 */



#include	"basics.h"

#include	"classes.h"

#include	"library.h"

#include	"tripos.h"




/*
 *	main
 *
 *	Author:	Christian Schafmeister (1991)
 */
void main( argc, argv )
int             argc;
char*           argv[];
{
FILE*           fIn;
UNIT		uUnit;
STRING		sName, sLibrary;
LIBRARY		lLibrary;


    BasicsInitialize();

    if ( argc != 2 ) {
        fprintf( stderr, "Usage: %s TRIPOS_FILE\n", 
			argv[0] );
        exit(0);
    }

    fIn = FOPENCOMPLAIN( argv[1], "r" );
    if ( fIn == NULL ) {
        printf( "Could not open file: %s\n", argv[1] );
        exit(1);
    }


	/* Read a TRIPOS UNIT */

    while ( uUnit = uTriposReadUnit(fIn) ) {

		/* Rename the UNIT so that leap doesnt choke on it */

	strcpy( sName, "u" );
	strcat( sName, sContainerName(uUnit) );
	ContainerSetName( uUnit, sName );
	printf( "Read unit: %s\n", sContainerName(uUnit) );

	strcpy( sLibrary, sContainerName(uUnit) );
	strcat( sLibrary, ".lib" );

	lLibrary = lLibraryOpen( sLibrary, OPENREADWRITE );
	LibrarySave( lLibrary, sContainerName(uUnit), (OBJEKT) uUnit, NULL );
	LibraryClose( &lLibrary );
    }

    fclose( fIn );

}
