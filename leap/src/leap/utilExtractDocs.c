/*
 *      File:   utilExtractDocs.c
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
 */


#include	"basics.h"


#define	CMDBEGIN	"\\beginCcmd"
#define	CMDEND		"\\endCcmd"
#define	CMDLABEL	"\\Ccmd"


#define MAX_HELPS	100

typedef	struct	LINEs {
		STRING		sText;
		struct LINEs 	*lNext;
		} LINEt;

typedef	LINEt*	LINE;

typedef	struct	{
		STRING		sKeyword;
		LINE		lFirst;
		LINE		lLast;
		int		iLine;
		} HELPt;





/*
 *	zUEDAppendHelpLine
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Append a line of text to the help.
 */
void	zUEDAppendHelpLine( hPHelp, sLine )
HELPt*		hPHelp;
char*		sLine;
{
LINE		lNew;


    MALLOC( lNew, LINE, sizeof(LINEt) );
    lNew->lNext = NULL;
    strcpy( lNew->sText, sLine );
    if ( hPHelp->lFirst == NULL ) {
	hPHelp->lFirst = lNew;
	hPHelp->lLast = lNew;
    } else {
	hPHelp->lLast->lNext = lNew;
	hPHelp->lLast = lNew;
    }
    MESSAGE(( "Appending: %s\n", sLine ));
}



/*
 *	zUEDDefineKeyword
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Define the keyword for the help.
 */
void	zUEDDefineKeyword( hPHelp, sKeyword )
HELPt*		hPHelp;
STRING		sKeyword;
{


    strcpy( hPHelp->sKeyword, sKeyword );

}


main( argc, argv )
int		argc;
char*		argv[];
{
int		iKeyword;
STRING          sLine, sTemp;
FILE*		fIn;
FILE*		fOut;
char*		cPKeyword;
HELPt		haHelp[MAX_HELPS];
int		i;
LINE		lCur;
int		iLine;

#define	FGETS(s,l,f) { fgets(s,l,f); \
if ( s[strlen(s)-1] == '\n' ) {s[strlen(s)-1] = '\0';} }


    BasicsInitialize();
    if ( argc != 3 ) {
	printf( "Usage: %s {source} {TeXfile}\n", argv[0] );
	exit(1);
    }
    fIn = fopen( argv[1], "r" );
    fOut = fopen( argv[2], "w" );
    iKeyword = 0; 
    iLine = 0;
    MESSAGE(( "Source file: %s\n", argv[1] ));
    MESSAGE(( "TeX file   : %s\n", argv[2] ));
    while ( !feof(fIn) ) { 
        strcpy( sLine, "" ); 
        FGETS( sLine, sizeof(STRING), fIn ); 
	MESSAGE(( "Line: %s\n", sLine ));
	MESSAGE(( "CMP : %s\n", CMDBEGIN ));
	if ( strncmp( sLine, CMDBEGIN, strlen(CMDBEGIN) ) == 0 ) { 
	    MESSAGE(( "Starting help\n" ));
	    haHelp[iKeyword].lLast = NULL; 
	    haHelp[iKeyword].lFirst = NULL;
	    haHelp[iKeyword].iLine = ++iLine;
	    strcpy( haHelp[iKeyword].sKeyword, "" );
	    zUEDAppendHelpLine( &(haHelp[iKeyword]), sLine );
	    while ( !feof(fIn) ) {
		strcpy( sLine, "" );
		FGETS( sLine, sizeof(STRING), fIn );
		iLine++;
		zUEDAppendHelpLine( &(haHelp[iKeyword]), sLine );
	        if ( strncmp( sLine, CMDEND, strlen(CMDEND) ) == 0 ) break;
		if ( strncmp( sLine, CMDLABEL, strlen(CMDLABEL) ) == 0 ) {
		    for ( i=0; i<strlen(sLine); i++ ) {
			if ( sLine[i] == '{' ) cPKeyword = &(sLine[i]);
			if ( sLine[i] == '}' ) sLine[i] = '\0';
		    }
		    zUEDDefineKeyword( &(haHelp[iKeyword]), cPKeyword );
		}
	    }
	    if ( feof(fIn) ) {
		DFATAL(( "Missing closing docs starting on line: %d\n",
				haHelp[iKeyword].iLine ));
	    }
	    iKeyword++;
	}
    }

		/* Now sort the entries */

    SortByString( haHelp, iKeyword, 
    			sizeof(haHelp[0]), haHelp[0].sKeyword, TRUE );
			
		/* Then write them out */

    for ( i=0; i<iKeyword; i++ ) {
	lCur = haHelp[i].lFirst;
	while ( lCur ) {
	    fprintf( fOut, "%s\n", lCur->sText );
	    lCur = lCur->lNext;
	}
	fprintf( fOut, "\n\n\n\n" );
    }
    fprintf( fOut, "\n\n" );


    exit(0);
}
