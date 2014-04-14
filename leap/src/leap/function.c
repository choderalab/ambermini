/*
 *	File:	function.c
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
 *		Store filename/start-stop line numbers/function name
 *		and index them by integer.
 *
 */


#include	"basics.h"

#include	"varArray.h"

#include	"function.h"



VARARRAY	GvaFunctions = NULL;
VARARRAY	GvaFiles = NULL;


/*
 *	FunctionAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a function definition to the VARARRAY.
 */
void
FunctionAdd( char *sFilename, int iStart, int iStop, char *sFunction )
{
FUNCTIONt	fFunction;

    strcpy( fFunction.sFilename, sFilename );
    fFunction.iStart = iStart;
    fFunction.iStop = iStop;
    strcpy( fFunction.sFunction, sFunction );
    fFunction.bPrint = FALSE;

    if ( GvaFunctions == NULL ) {
	GvaFunctions = vaVarArrayCreate( sizeof(FUNCTIONt) );
    }
    VarArrayAdd( GvaFunctions, (GENP)&fFunction );
}



/*
 *	iFunctionFindWithFilenameLine
 *
 *	Find and return the index of the function using
 *	the functions filename and a line number within the
 *	function.
 */
int	
iFunctionFindWithFilenameLine( char *sFilename, int iLine )
{
int		iCur;
FUNCTIONt	*fPFunc;

    fPFunc = PVAI(GvaFunctions,FUNCTIONt,0);
    for ( iCur = 0; iCur<iVarArrayElementCount(GvaFunctions); iCur++ ) {
	if ( fPFunc->iStart <= iLine && iLine <= fPFunc->iStop ) {
	    if ( strcmp( sFilename, fPFunc->sFilename ) == 0 ) {
		return(iCur);
	    }
	}
	fPFunc++;
    }
    return(NO_FUNCTION);
}



/*
 *	iFunctionFindWithFunction
 *
 *	Find and return the index of the function using
 *	the functions name.
 */
int	
iFunctionFindWithFunction( char *sFunction )
{
int		iCur;
FUNCTIONt	*fPFunc;

    fPFunc = PVAI(GvaFunctions,FUNCTIONt,0);
    for ( iCur = 0; iCur<iVarArrayElementCount(GvaFunctions); iCur++ ) {
	if ( strcmp( sFunction, fPFunc->sFunction ) == 0 ) {
	    return(iCur);
	}
	fPFunc++;
    }
    DFATAL(( "Cannot find function: %s\n", sFunction ));
}




/*
 *	FunctionFileAdd
 *
 *	Add a filename to the file list.
 */
void	
FunctionFileAdd( char *sName )
{
STRING		sFilename;

    if ( GvaFiles == NULL ) {
	GvaFiles = vaVarArrayCreate(sizeof(STRING));
    }
    strcpy( sFilename, sName );
    VarArrayAdd( GvaFiles, (GENP)sFilename );
}






/*
 *	iFunctionFilenameFind
 *
 *	Find the filename in the filename array.
 */
int
iFunctionFilenameFind( char *sFilename )
{
int		i;


    for ( i=0; i<iVarArrayElementCount(GvaFiles); i++ ) {
	if ( strcmp( sFilename, PVAI(GvaFiles,STRING,i) ) == 0 ) {
	    return(i);
	}
    }
    return(NO_FUNCTION);
}




