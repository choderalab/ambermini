/*
 *	File:	ring.h
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
 *		Maintain a RING structure.  A doubly linked list of
 *		nodes where insertions can be made relative to other
 *		elements within the RING.  Each node of the RING
 *		has a pointer to an arbitrary data type.
 *
 */



#include	"basics.h"

#include	"ring.h"



/*
 *-------------------------------------------------------------------
 *
 *	Private routines
 */


/*
 *	zrnPRingCreateNode
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create and return a RING node.
 */
static RINGNODEt *
zrnPRingCreateNode()
{
RINGNODEt	*rnPNew;

    MALLOC( rnPNew, RINGNODEt*, sizeof(RINGNODEt) );

    rnPNew->PData = NULL;
    rnPNew->rnPPrev = NULL;
    rnPNew->rnPNext = NULL;

    return(rnPNew);
}


/*
 *	zRingNodeDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the RING node and set the pointer to NULL.
 */
static void
zRingNodeDestroy( RINGNODEt **rnPPNode )
{
    FREE( (*rnPPNode) );
    *rnPPNode = NULL;
}




/*
 *	zrnPRingFind
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Find the node that has the data the caller wants.
 *	If PData is NULL then just return the first node.
 */
static RINGNODEt *
zrnPRingFind( RING rRing, GENP PData )
{
RINGLOOP	rlSearch;
GENP		PCur;

    if ( PData == NULL ) 
	return(rRing->rnPFirst);

    rlRingLoop(rRing, &rlSearch);
    while ( (PCur = PRingNext(&rlSearch)) != NULL ) {
	if ( PCur == PData ) 
		return(rnPRingNode(rlSearch));
    }
    return(NULL);
}




/*
 *=================================================================
 *
 *	Public routines
 *
 */



/*
 *	rRingCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create and return a RING.
 */
RING
rRingCreate()
{
RING	rNew;

    MALLOC( rNew, RING, sizeof(RINGt) );

    rNew->iElements = 0;
    rNew->rnPFirst = NULL;
    return(rNew);
}




/*
 *	RingDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the RING.
 *	The caller is expected to take care of the contents.
 */
void
RingDestroy( RING *rPRing )
{
RINGLOOP	rlElements;

	/* Destroy all of the nodes of the RING */

    rlRingLoop(*rPRing, &rlElements);
    while ( PRingNext(&rlElements) != NULL ) {
	FREE( rnPRingNode(rlElements) );
    }

    FREE( *rPRing );
}





/*
 *	RingAfterAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a node to the RING directly after the node with
 *	the caller specified data in (PAfter).
 *	If PAfter is NULL then the data will be added after the
 *	first node.
 */
void
RingAfterAdd( RING rRing, GENP PAfter, GENP PData )
{
RINGNODEt	*rnPCur, *rnPNew;

    rnPNew = zrnPRingCreateNode();
    rnPNew->PData = PData;

		/* If the RING is empty then just add the node */

    if ( bRingEmpty(rRing) ) {
        rRing->rnPFirst = rnPNew;
    } else {
	rnPCur = zrnPRingFind( rRing, PAfter );
	if ( rnPCur == NULL )
		DFATAL(( "RingAfterAdd() no PAfter\n" ));

	rnPNew->rnPPrev = rnPCur;
	rnPNew->rnPNext = rnPCur->rnPNext;

	rnPCur->rnPNext = rnPNew;
	if (rnPNew->rnPNext != NULL)
		rnPNew->rnPNext->rnPPrev = rnPNew;
    }
    rRing->iElements++;
}


	


/*
 *	RingBeforeAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a node to the RING directly before the node with
 *	the caller specified data in (PBefore).
 */
void
RingBeforeAdd( RING rRing, GENP PBefore, GENP PData )
{
RINGNODEt	*rnPCur, *rnPNew;

    rnPNew = zrnPRingCreateNode();
    rnPNew->PData = PData;

		/* If the RING is empty then just add the node */

    if ( bRingEmpty(rRing) ) {
        rRing->rnPFirst = rnPNew;
    } else {
	rnPCur = zrnPRingFind( rRing, PBefore );
	if ( rnPCur == NULL )
		DFATAL(( "RingAfterAdd() no PBefore\n" ));

	rnPNew->rnPNext = rnPCur;
	rnPNew->rnPPrev = rnPCur->rnPPrev;

	rnPCur->rnPPrev = rnPNew;
	if (rnPNew->rnPPrev != NULL)
		rnPNew->rnPPrev->rnPNext = rnPNew;
	if (rRing->rnPFirst == rnPCur)
        	rRing->rnPFirst = rnPNew;
    }
    rRing->iElements++;
}

	


/*
 *	rlRingLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Set up a loop over the elements of the RING.
 */
void
rlRingLoop( RING rRing, RINGLOOP *rlP )
{
    rlP->iCount = iRingSize(rRing);
    rlP->rnPNext = rRing->rnPFirst;
    rlP->rnPCur  = NULL;
}




/*
 *	PRingNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a pointer to the next piece of data in the RING.
 */
GENP
PRingNext( RINGLOOP *rlPLoop )
{
    rlPLoop->iCount--;
    if ( rlPLoop->iCount < 0 ) 
	return(NULL);
    rlPLoop->rnPCur = rlPLoop->rnPNext;
    rlPLoop->rnPNext = rlPLoop->rnPNext->rnPNext;
    return(rlPLoop->rnPCur->PData);
}



/*
 *	bRingRemove
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Remove the node with the data PData from the RING.
 */
BOOL
bRingRemove( RING rRing, GENP PData )
{
RINGNODEt	*rnPNode;

	/* Look for the data in the RING */

    rnPNode = zrnPRingFind( rRing, PData );
    if ( rnPNode == NULL ) 
	return(FALSE);

	/* Remove the node from the RING */

    rRing->iElements--;
    if ( rRing->iElements == 0 ) {
	zRingNodeDestroy( &rnPNode );
	rRing->rnPFirst = NULL;
    } else {
	if (rnPNode->rnPPrev != NULL)
		rnPNode->rnPPrev->rnPNext = rnPNode->rnPNext;
	if (rnPNode->rnPNext != NULL)
		rnPNode->rnPNext->rnPPrev = rnPNode->rnPPrev;
	if (rRing->rnPFirst == rnPNode)
		rRing->rnPFirst = rnPNode->rnPNext;
	zRingNodeDestroy( &rnPNode );
    }
    return(TRUE);
}



