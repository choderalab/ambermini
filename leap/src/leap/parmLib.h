/*
 *      File:   parmLib.h
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
 *              A PARMLIB maintains an ordered list of PARMLIBs
 *              which are searched sequentially for parameters.
 */

#ifndef PARMLIB_H
#define PARMLIB_H

# include	"classes.h"


typedef struct  {
	LIST            lParmSets;
	LISTLOOP        llParmSetLoop;
} PARMLIBt;

typedef PARMLIBt	*PARMLIB;


	/* Maintain a default PARMLIB */
extern	PARMLIB	GplDefaultParmLib;


extern PARMLIB	plParmLibCreate();
extern void	ParmLibDestroy(PARMLIB *plPLib);
extern void	ParmLibAddParmSet(PARMLIB plLib, PARMSET psSet);

extern void	ParmLibParmSetLoop(PARMLIB plLib);
extern BOOL	bParmLibNextParmSet(PARMLIB plLib, PARMSET *psPSet);

#define iParmLibSize( p )	iListSize(((PARMLIB)(p))->lParmSets)


		/* Manage a default PARMLIB */

#define	ParmLibDefineDefault( pl )	( GplDefaultParmLib = pl )
#define	bParmLibDefaultExists()		( GplDefaultParmLib != NULL )
#define	plParmLibDefault()		GplDefaultParmLib


/*
 *	PARMLIB_LOOP, PARMLIB_DEFAULT_LOOP
 *
 *	Creates a loop over all of the PARMSETs in the
 *	PARMLIB, testing each for the condition (cond)
 *	if (cond) != PARM_NOT_FOUND then the loop will stop.
 *
 *	(cond) Must be of the form: ( id = iParmSetFindxxxx() )
 *	the caller can then test id to see if the parameter was 
 *	found.
 */


#define	PARMLIB_LOOP( pl, ps, cond ) \
	ParmLibParmSetLoop(pl); \
	while ( bParmLibNextParmSet( pl, &ps ) ) { \
    		if ( (cond) != PARM_NOT_FOUND ) break; \
	}


/*
 *	PARMLIB_LOOP_ALL can be used to loop over all PARMSETs
 */
#define	PARMLIB_LOOP_ALL( pl, ps ) \
	ParmLibParmSetLoop(pl); \
	while ( bParmLibNextParmSet( pl, &ps ) ) 


#define	PARMLIB_DEFAULT_LOOP( ps, cond ) \
	ParmLibParmSetLoop(GplDefaultParmLib); \
	while ( bParmLibNextParmSet( GplDefaultParmLib, &ps ) ) { \
    		if ( cond != PARM_NOT_FOUND ) break; \
	}


/*
 *	PARMLIB_DEFAULT_LOOP_ALL can be used to loop over all PARMSETs
 */
#define	PARMLIB_DEFAULT_LOOP_ALL( ps ) \
	ParmLibParmSetLoop(GplDefaultParmLib); \
	while ( bParmLibNextParmSet( GplDefaultParmLib, &ps ) ) 


#endif /* PARMLIB_H */
