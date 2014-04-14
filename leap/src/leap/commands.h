/*
 *      File:   commands.h
 * B
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
 *      External declarations for ALL command handlers.
 *
 *      Each command handler has it's name prefixed with 'Cmd' and
 *      has two arguments, an integer representing the number
 *      of arguments and an array of OBJEKTS which are the arguments
 *
 */



#ifndef	COMMANDS_H
#define	COMMANDS_H

#ifndef	VARARRAY_H
#include	"varArray.h"
#endif

typedef OBJEKT  (*FUNCTION)();

typedef struct {
        STRING          sName;
        FUNCTION        fCallback;
} COMMANDt;

extern  COMMANDt        cCommands[];


typedef struct {
        STRING  sName;
        STRING  sCommand;
} ALIASt;

typedef ALIASt	*ALIAS;

extern	VARARRAY	GvaAlias;

OBJEKT
oCmd_charge( int, ASSOC* );

OBJEKT
oCmd_check( int, ASSOC* );
#endif
