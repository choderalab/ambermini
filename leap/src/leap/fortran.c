/*
 *      File:   fortran.c
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
 *
 *              This file contains routines for performing
 *              FORTRAN style formatted output.
 *              It contains a routine that allows the caller
 *              to specify the number of entries that should
 *              appear on each line and the format of each
 *              entry on the line.
 *              Then there is a routine that the user uses
 *              to write each subsequent piece of data.
 *
 *		Debugging assistance is given with 
 *			FortranDebugOn()
 *			FortranDebug(str)
 *			GsFortranDebug
 *
 *		When FortranDebugOn is called, subsequent
 *		calls to FortranDebugStr write strings into
 *		the output file that are always ignored by
 *		these fortran input routines, but can
 *		be used by the programmer to get their
 *		bearings within these huge ugly files.
 *		Whenever a FortranDebug string is found
 *		by the reading routines, it is copied into
 *		the global variable GsFortranDebug.
 *
 */



#include	"basics.h"

#include        "stringExtra.h"


/*
 *--------------------------------------------------------------------
 *
 *      Static variables 
 *
 */


#define	FORTRAN_DEBUG_COMMENT_CHAR	')'


static  STRING  SsFormat;
static  FILE*   SfFile;
static  int     SiPerLine;
static  int     SiOnLine = 0;   /* Stores the # of entries already on line */
static  BOOL    SbWroteNothing;

static  STRING  SsInput;
static  BOOL    SbNeedInput = TRUE;

static	BOOL	SbFortranDebug = FALSE;
static	STRING	GsFortranDebug;



/*
 *----------------------------------------------------------
 *
 *	Private routines
 *
 */

/* 
 *	zFortranGetInputLine
 *
 *	Read the next line into the STATIC variable SsInput.
 */
static void
zFortranGetInputLine()
{
   while ( !feof(SfFile) ) {
	SsInput[0] = '\0';
	fgets( SsInput, sizeof(SsInput), SfFile );
	if ( SsInput[0] == FORTRAN_DEBUG_COMMENT_CHAR ) {
	    strcpy( GsFortranDebug, &(SsInput[1]) );
	    MESSAGE(( "FORTRAN DEBUG:%s", GsFortranDebug ));
	    continue;
	}
	break;
   }
}



/*
 *==========================================================
 *
 *	Public routines
 *
 */


/*
 *      FortranFile
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Save the file where output is to be placed.
 */
void
FortranFile( FILE *fOut )
{
    SbNeedInput = TRUE;
    SfFile = fOut;
}




/*
 *      FortranFormat
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Save the format and number of entries per line.
 */
void
FortranFormat( int iPerLine, char *sFormat )
{
    SiPerLine = iPerLine;
    strcpy( SsFormat, sFormat );
    SiOnLine = 0;
    SbWroteNothing = TRUE;
}



/*
 *      FortranWriteInt
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Write an integer out to the file, increment the number of entries on 
 *      the line, if it is equal to SiPerLine, then move to the next line.
 */
void
FortranWriteInt( int iVal )
{
    fprintf( SfFile, SsFormat, iVal );
    SiOnLine++;
    SbWroteNothing = FALSE;
    if ( SiOnLine >= SiPerLine ) {
        fprintf( SfFile, "\n" );
        SiOnLine = 0;
    }
}

  



/*
 *      FortranWriteDouble
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Write a double out to the file, increment the number of entries on 
 *      the line, if it is equal to SiPerLine, then move to the next line.
 */
void
FortranWriteDouble( double dVal )
{
    fprintf( SfFile, SsFormat, dVal );
    SiOnLine++;
    SbWroteNothing = FALSE;
    if ( SiOnLine >= SiPerLine ) {
        fprintf( SfFile, "\n" );
        SiOnLine = 0;
    }
}




/*
 *      FortranWriteString
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Write a String out to the file, increment the number of entries on 
 *      the line, if it is equal to SiPerLine, then move to the next line.
 */
void
FortranWriteString( char *sVal )
{
    fprintf( SfFile, SsFormat, sVal );
    SiOnLine++;
    SbWroteNothing = FALSE;
    if ( SiOnLine >= SiPerLine ) {
        fprintf( SfFile, "\n" );
        SiOnLine = 0;
    }
}




/*
 *      FortranEndLine
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      If the number of objects on the line is not zero then print
 *      and end of line.
 */
void
FortranEndLine()
{
    if ( SbWroteNothing || SiOnLine != 0 ) 
	fprintf( SfFile, "\n" );
    SbWroteNothing = TRUE;
    SiOnLine = 0;
}


/*
 *      sFortranReadString
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Read the next string from the file
 */
char *
sFortranReadString( char *sString )
{
BOOL            bGotOne;
int             iLen;

    bGotOne = FALSE;
    strcpy( sString, "" );
    while ( !bGotOne ) {
        if ( SbNeedInput ) {
            if ( feof(SfFile) ) return(NULL);
            strcpy( SsInput, "" );
	    zFortranGetInputLine();
            if ( (iLen=strlen(SsInput)) > 0 ) SsInput[iLen-1] = '\0';
        }
        SbNeedInput = FALSE;
        sRemoveLeadingSpaces( SsInput );
        sRemoveFirstString( SsInput, sString );
        if ( strlen(sString)!=0 ) bGotOne = TRUE;
        SbNeedInput = ( strlen(SsInput)==0 );
   }
   return(sString);
}




/*
 *      sFortranReadLabel
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Read the next label from the file
 */
char *
sFortranReadLabel( char *sString )
{
BOOL            bGotOne;
int             iLen;

    bGotOne = FALSE;
    strcpy( sString, "" );
    while ( !bGotOne ) {
        if ( SbNeedInput ) {
            if ( feof(SfFile) ) return(NULL);
	    zFortranGetInputLine();
            if ( (iLen=strlen(SsInput)) > 0 ) SsInput[iLen-1] = '\0';
        }
        SbNeedInput = FALSE;
        strncpy( sString, SsInput, 4 );
        sString[4] = '\0';
        strcpy( SsInput, SsInput+4 );
        if ( strlen(sString)!=0 ) bGotOne = TRUE;
        SbNeedInput = ( strlen(SsInput)==0 );
   }
   return(sString);
}



/*
 *      iFortranReadInt
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Read the next string and return it as an integer.
 */
int
iFortranReadInt( )
{
STRING          sInt;
int             iVal;

    sFortranReadString( sInt );
    sscanf( sInt, "%d", &iVal );
    return(iVal);
}




/*
 *      dFortranReadDouble
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Read the next string and return it as a double.
 */
double  
dFortranReadDouble( )
{
STRING          sDouble;
double          dVal;

    sFortranReadString( sDouble );
    sscanf( sDouble, "%lf", &dVal );
    return(dVal); 
}



/*
 *      FortranSkipLine
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Skip a line of input
 */
void
FortranSkipLine()
{
    if ( feof(SfFile) ) 
	return;
    zFortranGetInputLine();
}




/*
 *	FortranDebugOn
 *
 *	Turn fortran debugging ON.  Default is OFF.
 *
 *	When on, FortranDebug calls write single strings
 *	Into the fortran format file which will be ignored
 *	by fortran reading routines, but can be
 *	used by the programmer to get their bearings within
 *	the nasty looking files.
 */
void
FortranDebugOn()
{
    SbFortranDebug = TRUE;
}

/*
 *	FortranDebug
 *
 *	If SbFortranDebug is TRUE then write the
 *	string to the output file with a comment character
 *	prefixed to it.
 */
void
FortranDebug( char *sStr )
{
    if ( SbFortranDebug ) {
        fprintf( SfFile, "%c%s\n", FORTRAN_DEBUG_COMMENT_CHAR, sStr );
    }
}
