/*
 *	File:	hash.h
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
 *		Handle hash tables.
 *
 */

#ifndef HASH_H
#define	HASH_H

typedef	struct	{
	char	*cPKey;
	GENP	PData;
} HASH_ENTRYt;


typedef	struct HASH_LISTs {
	HASH_ENTRYt		heEntry;
	struct HASH_LISTs	*hlPNext;
} HASH_LISTt;


typedef	struct	{
	HASH_LISTt	**hlPaTable;
	int		iTableSize;
	int		iElementCount;
} HASH_TABLEt;


typedef	HASH_TABLEt	*HASH_TABLE;


extern HASH_TABLE	htHTCreate(int iTableSize);
extern void		HTDestroy(HASH_TABLE *htPTable);
extern void		HTDescribe(HASH_TABLE htTable);

extern void		HTAdd(HASH_TABLE htTable, char *sKey, GENP PData);
extern GENP		PHTFind(HASH_TABLE htTable, char *sKey);
extern BOOL		bHTDelete(HASH_TABLE htTable, char *sKey, 
				char **cPPKey, GENP *PPData);

extern void		HTWalk(HASH_TABLE htTable, GENP PCallback, 
				VFUNCTION VFunc);


#define	iHTTableSize(h)		((h)->iTableSize)
#define	iHTElementCount(h)	((h)->iElementCount)


#endif /* HASH_H */
