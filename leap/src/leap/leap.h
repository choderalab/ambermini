/*
 *      File:   leap.h
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
 *      Contains all definitions that are used in LEaP
 */
 

#ifndef	LEAP_H
#define	LEAP_H

#ifndef	PARMLIB_H
# include        "parmLib.h"
#endif

#ifndef	AMBER_H
# include        "amber.h"
#endif




/* 
 *----------------------------------------------------
 *
 *      Global variables for LEaP
 */

#define LEAPRC  "leaprc"                /* Startup filename */

extern	BOOL	GbGraphicalEnvironment;	/* TRUE if in a graphical environment*/







/*
 *      variables.c
 */
extern void	VariablesInit();
extern void	VariablesList();
extern void	VariableSet(char *sName, OBJEKT oObj);
extern void	VariableRemove( char *sName );
extern OBJEKT	oVariable(char *sName);
extern void	VariablesDestroy();
extern DICTIONARY dVariablesDictionary();



/*
 *      set.c
 */
extern  void    SetAttribute();         /* ( OBJEKT, STRING, OBJEKT ) */


#endif
