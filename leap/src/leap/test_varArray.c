/*
 *	File:	test_varArray.c
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
 */

#include	"basics.h"

#include	"varArray.h"


main()
{
VARARRAY	vaTest;

int		i;

    BasicsInitialize();

    vaTest = vaVarArrayCreate(sizeof(int));
    i=1;
    VarArrayInsertBefore( vaTest, 0, &i );
    i=3;
    VarArrayInsertBefore( vaTest, 1, &i );
    i=2;
    VarArrayInsertBefore( vaTest, 1, &i );
    i = 0;
    VarArrayInsertBefore( vaTest, 0, &i );

    for ( i=0; i<iVarArrayElementCount(vaTest); i++ ) {
	printf( "varArray[%d] = %d\n", i, *PVAI(vaTest,int,i) );
    }
}

