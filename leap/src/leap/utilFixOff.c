/*
 *	File:	utilFixOff.c
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
 *		Fix OFF files.
 */


#include	<stdio.h>



typedef	char	STRING[500];

#define	TEMPFILE	"temp1234"
#define	COPY		0
#define	GET_INDEX	1

STRING	saTitles[200];
int	iTitles;


void	toLower( cPStr )
char*		cPStr;
{
    while ( *cPStr ) {
	if ( *cPStr>='A' && *cPStr <='Z' ) {	
	    *cPStr = *cPStr +32;
	}
	cPStr++;
    }
}


void	splitHeader( cPHeader, cPFront, cPName, cPTail )
char*		cPHeader;
char*		cPFront;
char*		cPName;
char*		cPTail;
{
int	iFirst, iSecond;

    iFirst = strlen( "!entry." ) - 1;
    for ( iSecond = iFirst+1; iSecond < strlen(cPHeader); iSecond++ ) {
	if ( cPHeader[iSecond] == '.' ) break;
    }
    if ( iSecond == strlen(cPHeader) ) {
	printf( "Illegal header: %s\n", cPHeader );
	exit(1);
    }

    strncpy( cPFront, cPHeader, iFirst+1 );
    strncpy( cPName, &(cPHeader[iFirst+1]), iSecond-iFirst-1 );
    strcpy( cPTail, &(cPHeader[iSecond]) );

}

	


void correctHeader( cPLine )
char*		cPLine;
{
STRING		sTemp;
STRING		sHead, sName, sTail;
int		i;

    if ( strncmp( cPLine, "!entry", strlen("!entry") ) == 0 ) {
	splitHeader( cPLine, sHead, sName, sTail );
	for ( i=0; i<iTitles; i++ ) {
	    strcpy( sTemp, saTitles[i] );
	    toLower( sTemp );
	    if ( strcmp( sName, sTemp )== 0 ) {
		strcpy( sName, saTitles[i] );
		break;
	   }
	}
	if ( i == iTitles ) {
	    printf( "Could not find name: %s in header: %s\n", sName, cPLine );
	    printf( "sHead=%s\n", sHead );
	    printf( "sName=%s\n", sName );
	    printf( "sTail=%s\n", sTail );
	}
	strcpy( cPLine, sHead );
	strcat( cPLine, sName );
	strcat( cPLine, sTail );
    }

}


 

main( argc, argv )
int		argc;
char*		argv[];
{
int		iState, i;
STRING		sLine;
STRING		sEntry;
FILE*		fIn;
FILE*		fOut;

    iTitles = 0;

    if ( argc != 2 ) {
	printf( "Usage: %s OFF_FILE\n", argv[0] );
	exit(1);
    }

    fIn = fopen( argv[1], "r" );
    if ( fIn == NULL ) {
	printf( "Could not open file: %s\n", argv[1] );
	exit(1);
    }
    fOut = fopen( TEMPFILE, "w" );
    if ( fOut == NULL ) {
	printf( "Could not open temporary file: %s\n", TEMPFILE );
	exit(1);
    }

    iState = COPY;

    while ( !feof(fIn) ) {
	fgets( sLine, sizeof(sLine), fIn );
	if ( feof(fIn) ) break;

		/* Change state */

	switch ( iState ) {
	    case COPY:
		if ( strncmp( sLine, "!!index", strlen("!!index") ) == 0 ) {
		    iState = GET_INDEX;
		}
		break;
	    case GET_INDEX:
		if ( strncmp( sLine, "!entry", strlen("!entry") ) == 0 ) {
		    iState = COPY;
		}
		break;
	}


		/* act on the current state */

	switch ( iState ) {
	    case GET_INDEX:
		if ( sLine[0] == '!' ) break;
		for ( i=0; i<strlen(sLine)-4; i++ ) {
		    sEntry[i] = sLine[i+2];
		}
		sEntry[i] = '\0';
		strcpy( saTitles[iTitles], sEntry );
		printf( "Got title: %s\n", sEntry );
		iTitles++;
		break;
	    case COPY:
		correctHeader( sLine );
		break;
	}
	fprintf( fOut, "%s", sLine );


    }

    fclose( fIn );
    fclose( fOut );

	/* Rename the files */
    sprintf( sLine, "rm %s", argv[1] );
    system( sLine );

    sprintf( sLine, "mv %s %s", TEMPFILE, argv[1] );
    system( sLine );

}


