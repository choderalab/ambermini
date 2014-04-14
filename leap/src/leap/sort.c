/*
 *      File: sort.c
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
 *              This is an interface to the general sorting routine
 *              'qsort'.  It allows for sorting by string, integer, and
 *              double keys.
 */



#include	"basics.h"


static  int     SiAscending;
static  int     SiOffset;

/*
 * ---------------------------------------------------------------
 *
 *      Comparison routines
 *
 */


/*
 *	iCompareInteger
 *
 *	Author:	Christian Schafmeister (1991)
 */
static int
iCompareInteger( GENP PA, GENP PB )
{
int	*iPA, *iPB;

    iPA = (int*)( ((char*)PA) + SiOffset );
    iPB = (int*)( ((char*)PB) + SiOffset );

    if ( *iPA < *iPB ) return(-1*SiAscending);
    if ( *iPA == *iPB ) return(0);
    return(SiAscending);
}
 

/*
 *	iCompareDouble
 *
 *	Author:	Christian Schafmeister (1991)
 */
static int
iCompareDouble( GENP PA, GENP PB )
{
double	*dPA, *dPB;

    dPA = (double*)( ((char*)PA) + SiOffset );
    dPB = (double*)( ((char*)PB) + SiOffset );

    if ( *dPA < *dPB ) return(-1*SiAscending);
    if ( *dPA == *dPB ) return(0);
    return(SiAscending);
}
 


/*
 *	iCompareString
 *
 *	Author:	Christian Schafmeister (1991)
 */
static int
iCompareString( GENP PA, GENP PB )
{
char	*cPA, *cPB;

    cPA = (char*)( ((char*)PA) + SiOffset );
    cPB = (char*)( ((char*)PB) + SiOffset );

    return(strcmp( cPA, cPB )*SiAscending);
}



/* 
 *      SortByInteger
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Sort the array by into order by integer keys.
 *      The first integer in the first record is pointed to by
 *      PFirst.
 */
void
SortByInteger( GENP PStart, int iElements, int iSize, GENP PFirst, 
		BOOL bAscending )
{
    SiAscending = bAscending ? 1 : -1;
    SiOffset = (char*)PFirst - (char*)PStart;
    qsort( PStart, iElements, iSize, (int (*) (const void *, const void *) )iCompareInteger );
}

/* 
 *      SortByDouble
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Sort the array by into order by double keys.
 *      The first double in the first record is pointed to by
 *      PFirst.
 */
void
SortByDouble( GENP PStart, int iElements, int iSize, GENP PFirst, 
		BOOL bAscending )
{

    SiAscending = bAscending ? 1 : -1;
    SiOffset = (char*)PFirst - (char*)PStart;
    qsort( PStart, iElements, iSize, (int (*) (const void *, const void *) )iCompareDouble );
}




/* 
 *      SortByString
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Sort the array by into order by string keys.
 *      The first string in the first record is pointed to by
 *      PFirst.
 */
void
SortByString( GENP PStart, int iElements, int iSize, 
		GENP PFirst, BOOL bAscending )
{
    SiAscending = bAscending ? 1 : -1;
    SiOffset = (char*)PFirst - (char*)PStart;
    qsort( PStart, iElements, iSize, (int (*) (const void *, const void *) )iCompareString );
}


/*
 *	Sift
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Sift the array according to an arbitrary criteria.
 *	A sift seperates elements of an array into two groups
 *	one where the criteria is true and the other where
 *	it is false.  The group where the criteria is true
 *	is placed before the group where the criteria is
 *	false.
 */

typedef	BOOL	(*SIFTFUNCTION)();

void
Sift( GENP PData, int iElementSize, int iElements, SIFTFUNCTION bFCriteria, 
	int *iPFirstFalse )
{
GENP		PSwapBuffer, PCur, PTop;
int		iTopTrue, iBottomFalse;

	/* Temporarily allocate a swap buffer */


    MALLOC( PSwapBuffer ,GENP ,iElementSize );

    iTopTrue = -1;
    iBottomFalse = iElements;

    PCur = PData;
    PTop = (char*)PData + iElementSize*(iElements-1);
	
    while ( iTopTrue < iBottomFalse-1 ) {

		/* Test current element */

	if ( bFCriteria( PCur ) ) {

		/* If true then leave it alone and go on to the next */

	    PCur = (char*)PCur + iElementSize;
	    iTopTrue++;
	} else {

		/* If false then swap it with the bottom of the false */
		/* group, and leave the current pointer where it is  */

	    memmove( PSwapBuffer, PCur, iElementSize );
	    memmove( PCur, PTop, iElementSize );
	    memmove( PTop, PSwapBuffer, iElementSize );

	    PTop = (char*)PTop - iElementSize;
	    iBottomFalse--;
	}
    }
    FREE(PSwapBuffer);
    *iPFirstFalse = iBottomFalse;
}


