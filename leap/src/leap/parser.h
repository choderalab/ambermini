/*
 *	File:	parser.h
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
 *		Define typedefs for RESULTt.
 *		This also defines GLOBAL variables that are allocated
 *		within the parser and are used by several other files
 *		like 'tools.c'.
 */


#ifndef	PARSER_H
#define	PARSER_H



#ifndef	PARMLIB_H
# include	"parmLib.h"
#endif


/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 *
 *	Result structure
 *
 *	The RESULTt type is used to pass commands and information
 *	from the backend/parser part of LEaP to the frontend.
 *	The result codes are currently slanted towards commands
 *	to the graphics system.
 *
 *	The result codes are:
 *		CNONE	=	No command.
 *		CQUIT	=	Quit LEaP.
 *		CVIEW	=	Display the UNIT whose variable name
 *				is passed in the sVariable part of
 *				RESULTt on the graphics display.
 *		CRESETVIEW =	Reset the viewing matrix for the
 *				current UNIT.
 *
 */		

#define	CNONE		0
#define	CQUIT		1
#define	CVIEW		2
#define	CEDIT		3
#define CVERBOSITY      4

typedef	struct	{
	int	iCommand;
	STRING	sVariable;
	OBJEKT	oObject;    /* Davids Changes:  Was 'UNIT uUnit',*/
		                    /* but needed to be changed to be    */
		                    /* used by parmsets                  */
} RESULTt;

extern RESULTt		GrMainResult;
extern PARMLIB		GplAllParameters;
extern BOOL		GbGraphicalEnvironment;
extern STRING		GsProgramName;


/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *
 *	The following is used by the parser to maintain a stack of
 *	files where input is received from.  If the file is NULL
 *	then input is received from the main program in the form
 *	of BLOCKS.  The main program is then responsible for reading
 *	the stdin ( in the command line interface ) or for gathering
 *	keypress events from X-Windows ( in the graphical interface ).
 *
 */

#define	MAXINPUTFILES	10		/* Maximum 10 input files can be */
					/* open at once */

extern int	GiInputFileStackPos;
extern FILE	*GfaInputFileStack[MAXINPUTFILES];
extern FILE	*fINPUTFILE();


#define INPUTPUSHFILE( f )      ( GfaInputFileStack[++GiInputFileStackPos] = f )
#define INPUTPOPFILE()   	( GiInputFileStackPos-- )
#define iINPUTSTACKDEPTH()      ( GiInputFileStackPos )
#define bINPUTMAXDEPTH()        ( GiInputFileStackPos == MAXINPUTFILES-1 )


#endif
