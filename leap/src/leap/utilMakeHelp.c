/*
 *      File:   utilMakeHelp.c
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
 *		This program generates a c program that is
 *              used to initialize a help database.
 *              The help text is obtained from a special text
 *              file.
 *              This program takes text from the standard input
 *              and writes a 'c' program to the standard output.
 *              The standard input has the form:
 *			@@beginHelp: (keyWord)
 *				(lines of text)
 *			@@endHelp
 */


#include	"basics.h"

BOOL	GbGraphicalEnvironment;	/* HACK for lazy Imakefile */


#define	BEGINHELP	"@@beginHelp:"
#define	ENDHELP		"@@endHelp"


#define MAXKEYWORDS     100

#define MAXTEXTLINES    200
#define MAXTEXTLEN      80


typedef STRING  TEXTBLOCK[MAXTEXTLINES];



/*
 *      sOutputString
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a string that has quotes (") converted to (\")
 */
char *
sOutputString( char *sStr, char *sOut )
{
char	*sTemp;


    sTemp = sOut;
    while ( *sStr != '\0' ) {
        if ( *sStr == '"' ) {
            *sTemp = '\\';
            sTemp++;
        }
        *sTemp = *sStr;
        sTemp++;
        sStr++;
    }
    *sTemp = '\0';
    return(sOut);
}

    


/*
 *      Define global variables
 *
 *      This will take up a huge amount of storage according
 *      to the above defines.  About 1Meg.
 *
 */

STRING          saKeywords[MAXKEYWORDS];
int             iaLines[MAXKEYWORDS];
TEXTBLOCK       taText[MAXKEYWORDS];


int main( int argc, char *argv[] )
{
int             iKeyword, iLine, i, j;
STRING          sLine, sTemp;
FILE		*fIn, *fOut;
BOOL		bBlank;

#define	FGETS(s,l,f) { fgets(s,l,f); \
if ( s[strlen(s)-1] == '\n' ) {s[strlen(s)-1] = '\0';} }


    BasicsInitialize();
    if ( argc != 3 ) {
	printf( "Usage: %s {helptext} {helpcode}\n", argv[0] );
	exit(1);
    }
    if ((fIn = fopen( argv[1], "r" )) == NULL) {
	perror( argv[1] );
	exit(1);
    }
    if ((fOut = fopen( argv[2], "w" )) == NULL) {
	perror( argv[2] );
	exit(1);
    }
    iKeyword = 0;
    while ( !feof(fIn) ) {
	strcpy( sLine, "" );
        FGETS( sLine, sizeof(STRING), fIn );
	if ( !*sLine ) continue;
	if ( *sLine == '\f' ) continue;
	if ( strncmp( sLine, BEGINHELP, strlen(BEGINHELP) ) != 0 ) continue;
	sscanf( sLine, "%s %s", sTemp, saKeywords[iKeyword] );
	fflush(stdout);
        if ( feof(fIn) ) break;
        iLine = 0;
	bBlank = FALSE;
        while ( TRUE ) {
            FGETS( sLine, sizeof(STRING), fIn );
	    if ( !*sLine && !bBlank ) {
		bBlank = TRUE;
		continue;
	    }
	    bBlank = FALSE;
	    if ( *sLine == '\f' ) continue;
	    if ( strncmp( sLine, ENDHELP, strlen(ENDHELP) ) == 0 ) break;
	    strcpy( taText[iKeyword][iLine], sLine );
            iLine++;
	    if ( iLine > MAXTEXTLINES )
	    	DFATAL(( "Too many text lines, maximum: %d\n", MAXTEXTLINES ));
        }
        iaLines[iKeyword] = iLine;
        iKeyword++;
	if ( iKeyword > MAXKEYWORDS ) 
	    DFATAL(( "Too many keywords, maximum: %d\n", MAXKEYWORDS ));
    }

    fprintf( fOut, "/* This file was generated by utilMakeHelp.c */\n" );
    fprintf( fOut, "\n\n\n" );
    fprintf( fOut, "extern  void    HelpAdd();\n" );
    fprintf( fOut, "\n\n\n" );
    fprintf( fOut, "void HTInit()\n" );
    fprintf( fOut, "{\n" );
    for ( i=0; i<iKeyword; i++ ) {
        fprintf( fOut, "    HelpAdd( \"%s\", \"\\\n", 
                        saKeywords[i]);
        for ( j=0; j<iaLines[i]; j++ ) {
            fprintf( fOut, "%s\\n\\\n", sOutputString( taText[i][j], sTemp ) );
        }
        fprintf( fOut, "\\n\" );\n" );
    }
    fprintf( fOut, "}\n" );

    exit(0);
}
