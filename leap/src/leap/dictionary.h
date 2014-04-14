/*
 *      File: dictionary.h
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
 *              A DICTIONARY is an object that maintains
 *              a collection of entries that are made up of two
 *              parts: a key and a pointer to an object.
 *              The name can be a string up to 255 characters and
 *              is case sensitive.
 *
 *              The DICTIONARY allows the caller to:
 *                      create new dictionaries
 *                      destroy dictionaries
 *                      add new objects
 *                      return objects using the key
 *                      delete objects using the key
 *                      loop over the contents of the dictionary
 *				in ascending alphabetical order.
 */


#ifndef	DICTIONARY_H
#define	DICTIONARY_H


#include        "basics.h"
#include	"hash.h"
#include	"varArray.h"


typedef	struct	{
	char	*cPKey;
	GENP	PData;
} DICTLOOPt;

typedef DICTLOOPt	*DICTLOOP;


typedef	struct	{
	GENP		PKlass;
	HASH_TABLE	htEntries;
	VARARRAY	vaSort;		/* DICTIONARY_SORTt */
	DICTLOOP	dlLastSort;
} DICTIONARYt;

typedef	DICTIONARYt	*DICTIONARY;


		/* Dictionary CLASS structure */

extern	KLASSt	GkTDictionary;


/*
------------------------------------------------------------------

        PUBLIC methods.
*/

extern DICTIONARY	dDictionaryCreate();    
extern void		DictionaryDestroy(DICTIONARY *dPDict);
extern void		DictionaryAdd(DICTIONARY dDict, char *sKey, GENP PData);
extern GENP		yPDictionaryFind(DICTIONARY dDict, char *sKey);
extern GENP		yPDictionaryDelete(DICTIONARY dDict, char *sKey);
extern DICTLOOP		ydlDictionaryLoop(DICTIONARY dDict);
extern GENP		yPDictionaryNext(DICTIONARY dDict, DICTLOOP *dlPLoop);
extern void		yDictionaryDescribe(DICTIONARY dDict);

/* TODO: Implement a dlDictonaryFastLoop/PDictionaryFastNext that does */
/* TODO: NOT sort the list */


#define iDictionaryElementCount(d) \
		(iHTElementCount((d)->htEntries))
#define sDictLoopKey(dl)	((dl)->cPKey)
#define PDictLoopData(dl)	((dl)->PData)



#endif  /* ifndef DICTIONARY_H */
