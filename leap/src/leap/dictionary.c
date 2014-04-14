/*
 *      File: dictionary.c
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


#include	"basics.h"

#include	"hash.h"

#include	"dictionary.h"


/*
 *--------------------------------------------------
 *
 *	Define the DICTIONARY class structure
 */

KLASSt	GkTDictionary = {
	"DICTIONARY",
	NULL
};



#define	DICTIONARY_HASH_TABLE_SIZE	1009 /* Make this bigger if needed */


/*
 *	zDictionarySortCleanup
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	If there is a vaSort VARARRAY defined then
 *	destroy it.
 */
static void
zDictionarySortCleanup(DICTIONARY dDict)
{
    if ( dDict->vaSort ) {
	VarArrayDestroy( &(dDict->vaSort) );
    }
}


/*
 *	zDictonaryDestroyKey
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	This is a callback for HTWalk which
 *	is passed NULL, a key, and a data
 *	element.  This routine FREEs the memory
 *	for the key.
 */
static void
zDictionaryDestroyKey( DICTIONARY dDict, char *cPKey, GENP PData )
{
    FREE(cPKey);
}


/*
 *	zDictonaryFillSortVarArray
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	This is a callback for HTWalk which
 *	is passed the DICTIONARY, and
 *	a (char*) to a DICTIONARY entry key
 *	and a GENP to a DICTIONARY entry data.
 *	Add both pieces of data to the DICTIONARY sort
 *	VARARRAY in dDict->dlLastSort and increment
 *	dDict->dlLastSort.
 */
static void
zDictionaryFillSortVarArray( DICTIONARY dDict, char *cPKey, GENP PData )
{
    (dDict->dlLastSort)->cPKey = cPKey;
    (dDict->dlLastSort)->PData = PData;
    dDict->dlLastSort++;
}



/*
 *	ziDictionaryCompareDictLoopEntries
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	This is a callback for the Sort routine.
 *	Compare to DICTLOOP entries and
 *	return -1 if the first is less than the
 *	second, +1 if it is greator, and zero
 *	if they are the same.
 */
static int
ziDictionaryCompareDictLoopEntries( DICTLOOP dlA, DICTLOOP dlB )
{
    return(strcmp( dlA->cPKey, dlB->cPKey ));
}
 

/*
 *-------------------------------------------------------
 *
 *	Public routines.
 *
 */


/*
 *	dDictionaryCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create a DICTIONARY.
 */
DICTIONARY
dDictionaryCreate()
{
DICTIONARY	dNew;

    MALLOC( dNew, DICTIONARY, sizeof(DICTIONARYt) );
		/* Create a HASH_TABLE */

		/* Define the class */

    dNew->PKlass = (GENP)&GkTDictionary;
    dNew->htEntries = htHTCreate(DICTIONARY_HASH_TABLE_SIZE);

		/* Initialize the sort array */

    dNew->vaSort = NULL;
    dNew->dlLastSort = NULL;

    return(dNew);
}




/*
 *	DictionaryDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the DICTIONARY.
 */
void	
DictionaryDestroy( DICTIONARY *dPDict )
{

    zDictionarySortCleanup(*dPDict);

		/* First FREE all of the keys */

    HTWalk( (*dPDict)->htEntries, NULL, zDictionaryDestroyKey );

		/* Destroy the HASH_TABLE */
    HTDestroy( &((*dPDict)->htEntries) );

    (*dPDict)->PKlass = NULL;
    FREE((*dPDict));
    *dPDict = NULL;

}



/*
 *	DictionaryAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add an entry to the DICTIONARY.
 */
void
DictionaryAdd( DICTIONARY dDict, char *sKey, GENP PData )
{
int		iLen;
char		*cPKey;

    zDictionarySortCleanup(dDict);

		/* Allocate memory for the key */

    iLen = strlen(sKey);
    MALLOC( cPKey, char*, iLen+1 );
    strcpy( cPKey, sKey );

    HTAdd( dDict->htEntries, cPKey, PData );
}



/*
 *	yPDictionaryFind
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the data associated with the key.
 */
GENP
yPDictionaryFind( DICTIONARY dDict, char *sKey )
{
		/* Dont try to clean up the sort arrays */
		/* because we need SPEED!!! */

    return( PHTFind( dDict->htEntries, sKey ) );
}


/*
 *	yPDictionaryDelete
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Delete a DICTIONARY entry and return
 *	the data that was associated with the entry.
 */
GENP
yPDictionaryDelete( DICTIONARY dDict, char *sKey )
{
char		*cPKey;
GENP		PData;
BOOL		bFound;


    bFound = bHTDelete( dDict->htEntries, sKey, &cPKey, &PData );

    if ( bFound ) {
		/* Free the memory for the Key */
	FREE(cPKey);
	return(PData);
    }
    return(NULL);
}





/*
 *	ydlDictionaryLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Sort the DICTIONARY, store the sorted entries
 *	in the VARARRAY (vaSort) and return a DICTLOOP
 *	that can be passed to PDictionaryNext for
 *	each successive entry in the sorted list.
 */
DICTLOOP
ydlDictionaryLoop( DICTIONARY dDict )
{
		/* Cleanup an old sorted array */

    zDictionarySortCleanup(dDict);
    if ( iHTElementCount(dDict->htEntries) == 0 ) return(NULL);

		/* Setup a static variable for the */
		/* sort VARARRAY and walk through the */
		/* HASH_TABLE copying elements into the */
		/* VARARRAY and then sort them */

    dDict->vaSort = vaVarArrayCreate(sizeof(DICTLOOPt));
    VarArraySetSize( dDict->vaSort, iHTElementCount(dDict->htEntries) );
    dDict->dlLastSort = PVAI( dDict->vaSort, DICTLOOPt, 0 );

    HTWalk( dDict->htEntries, (GENP)dDict, zDictionaryFillSortVarArray );

		/* Now sort the VARARRAY */

    qsort( PVAI(dDict->vaSort,DICTLOOPt,0),
		iVarArrayElementCount(dDict->vaSort),
		iVarArrayElementSize(dDict->vaSort),
		(int (*) (const void *, const void *) )ziDictionaryCompareDictLoopEntries );

		/* Now return a DICTLOOP that can be */
		/* used by PDictionaryNext to get successive */
		/* sorted elements */

    return(NULL);
}


/*
 *	yPDictionaryNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return successive elements of the sorted
 *	DICTIONARY.
 *
 * TODO:	Handle nested DICTIONARY loops.
 *
 */
GENP
yPDictionaryNext( DICTIONARY dDict, DICTLOOP *dlPLoop )
{
    if ( *dlPLoop == NULL ) {
	if ( dDict->vaSort == NULL ) return(NULL);
	*dlPLoop = PVAI( dDict->vaSort, DICTLOOPt, 0 );
    } else {
	(*dlPLoop)++;
	if ( *dlPLoop >= dDict->dlLastSort ) {
	    *dlPLoop = NULL;
	    return(NULL);
	} 
    }
    return((*dlPLoop)->PData);
}




/*
 *	DictionaryDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Describe the DICTIONARY.
 */
void
DictionaryDescribe( DICTIONARY dDict )
{
DICTLOOP	dlLoop;

    VP0(( "Dictionary with %d elements.\n", 
		iDictionaryElementCount(dDict) ));
    dlLoop = ydlDictionaryLoop(dDict);
    while ( yPDictionaryNext( dDict, &dlLoop ) ) {
	VP0(( "%40s = %lX\n", 
		sDictLoopKey(dlLoop), 
		PDictLoopData(dlLoop) ));
    }
}


