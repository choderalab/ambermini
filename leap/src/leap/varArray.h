#ifndef  VARARRAY_H
#define  VARARRAY_H
 
/*
 *      File: varArray.h
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
 *              A VARARRAY is a regular array whose size can increase.
 *
 *
 *       Body of varArray.[ch] modified by Vladimir Romanovski (1994).
 *
 */

#include "basics.h"

typedef struct {
        int     count;  /* real count of element in array  */ 
        int     size;   /* all elements have the same size */
        int     slot;   /* max available count of element before new realloc*/
        char    *data;
} HeaderStruct, *VARARRAY;

extern int      iVarArrayElementSize(VARARRAY header);
extern int      iVarArrayElementCount(VARARRAY header);
extern int      iVarArrayPointerToIndex(VARARRAY header, char *data);
extern char     *PVarArrayIndex( VARARRAY header, int pos);

extern VARARRAY vaVarArrayCreate(int size);
extern void     VarArrayDestroy(VARARRAY *header);
extern void     VarArrayAdd(VARARRAY header, GENP data);
extern VARARRAY vaVarArrayCopy(VARARRAY header);
extern VARARRAY vaVarArrayCopy2(VARARRAY header1, VARARRAY header2);
extern void     VarArraySetSize(VARARRAY header, int ncount);
extern GENP     PVarArrayDebugIndex(VARARRAY header, int pos, 
                        char *file, int line);

extern void     VarArrayInsertBefore(VARARRAY header, int pos, GENP data);
extern void     VarArrayInsertBeforeMore(VARARRAY header, int pos, int num);
     
#define VarArrayDelete(t,i) VarArrayDeleteMore(t,i,1)

#define PVAI(va,tc,i) ((tc*)PVarArrayIndex(va,(i)))


#endif  /* VARARRAY_H */




































