/*
 *	File:	bag.c
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
 *		Bag is a doubly linked bag of GENP.
 */




#include	"basics.h"

#include	"bag.h"


/*
 *-------------------------------------------------------------------------
 *
 *       Private routines
*/



/*
 *	searchBag
 *
 *      Find the node that points to the object pointed to by the
 *      argument.  Return a pointer to the pointer to the node
 *      which points to the object.
 */
static BNODEP *
searchBag( BNODEP *bnPPBag, GENP PObj )
{
BNODEP   bnPCur, *bnPPPrev;

        bnPPPrev = bnPPBag;
        bnPCur = (*bnPPPrev);
                /* search for the node that contains the object */
        while ( bnPCur != NULL ) {
                if ( bnPCur->PObject == PObj ) break;
                bnPPPrev = &(bnPCur->bnPNextNode);
                bnPCur = (*bnPPPrev);
        }

        if ( bnPCur == NULL ) return(NULL);
        return(bnPPPrev);
}

        


/*
 *-----------------------------------------------------------------
 *
 *       Public routines
 *       
 */


/*
 *        bBagCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *        Return a new bag.
 */
BAG
bBagCreate()
{
BAG    bBag;

        MALLOC( bBag, BAG, sizeof(BAGt) );
	bBag->iElementCount = 0;
        bBag->bnPFirstNode = NULL;
        bBag->bnPLastNode = NULL;
        return(bBag);
}




/*
 *      BagDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the bag.  The caller is responsible for destroying the
 *      contents of the bag.
 */
void
BagDestroy( BAG *bPBag )
{
BNODEP   bnPTempPtr, bnPFree;

        /* Destroy the BAG nodes */

    bnPTempPtr = (*bPBag)->bnPFirstNode;
    while ( bnPTempPtr != NULL ) {
        bnPFree = bnPTempPtr;
        bnPTempPtr = bnPTempPtr->bnPNextNode;
        
                /* Free the node memory. */
        FREE(bnPFree);
    }

        /* Now destroy the BAG itself */

    FREE( *bPBag );
    *bPBag = NULL;
}



/*
 *        BagAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *        Add an element to the head of a bag.
 *        Simply return if the user tries to add NULL since
 *        NULLs take up no space.
 */
void
BagAdd( BAG bBag, GENP PObj )
{
BNODEP   bnPTempBag;

    if ( PObj == NULL ) 
	return;

                /* Insert the new node at the head of the bag */

    bnPTempBag = bBag->bnPFirstNode;
    MALLOC( bBag->bnPFirstNode, BNODEP, sizeof(BNODE) );

                /* Keep up to date the last node of the bag */

    if ( bnPTempBag == NULL ) {
        bBag->bnPLastNode = bBag->bnPFirstNode;
    }

    (bBag->bnPFirstNode)->PObject = PObj;
    (bBag->bnPFirstNode)->bnPNextNode = bnPTempBag;
    bBag->iElementCount++;
}





/*
 *      BagAddToEnd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add an element to the tail of a bag.
 *      Simply return if the user tries to add NULL since
 *      NULLs take up no space.
 */
void
BagAddToEnd( BAG bBag, GENP PObj )
{
BNODEP   bnPTempBag;

    if ( PObj == NULL ) 
	return;

                /* Insert the new node at the tail of the bag */

    bnPTempBag = bBag->bnPLastNode;
    MALLOC( bBag->bnPLastNode, BNODEP, sizeof(BNODE) );

                /* Keep up to date the last node of the bag */

    if ( bnPTempBag == NULL ) {
        bBag->bnPFirstNode = bBag->bnPLastNode;
    } else {
        bnPTempBag->bnPNextNode = bBag->bnPLastNode;
    }

    (bBag->bnPLastNode)->PObject = PObj;
    (bBag->bnPLastNode)->bnPNextNode = NULL;
    bBag->iElementCount++;
}




/*
 *      BagDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the string with a description of the bag.
 *      It is the callers responsibility to ensure that
 *      the string is long enough to store the description.
 */
void
BagDescribe( BAG bBag )
{
    VP0(( "Bag elt count=%d\n", bBag->iElementCount ));
}




        
/*
 *      bBagRemove
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Remove an object from the bag.
 *      
 *      Return:
 *              FALSE if the element is not in the bag.
 */
BOOL
bBagRemove( BAG bBag, GENP PPtr )
{
BNODEP   *bnPPPrev;
BNODEP   bnPNode;

                /* If the bag is empty then return */

    if ( bBag->bnPFirstNode==NULL ) return(FALSE);

                /* Search the bag for the object */
    bnPPPrev = searchBag( &(bBag->bnPFirstNode), PPtr );
    if ( bnPPPrev == NULL ) return(FALSE);

                /* Now actually remove the object from the bag */
    bnPNode = *bnPPPrev;
    (*bnPPPrev) = (*bnPPPrev)->bnPNextNode;

                /* If the node is the last node then update the last node */

    if ( bnPNode == bBag->bnPLastNode ) {
        bBag->bnPLastNode = *bnPPPrev;
    }

                /* Destroy the node itself */
    FREE(bnPNode);

    bBag->iElementCount--;

    return(TRUE);
}





        
/*
 *      bBagContains
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Find out if a bag contains a particular object.
 *      
 *      Return:
 *              FALSE if the element is not in the bag.
 */
BOOL
bBagContains( BAG bBag, GENP PPtr )
{
BNODEP   *bnPPPrev;

                /* If the bag is empty then return */

    if ( bBag->bnPFirstNode==NULL ) return(FALSE);

                /* Search the bag for the object */
    bnPPPrev = searchBag( &(bBag->bnPFirstNode), PPtr );
    if ( bnPPPrev == NULL ) return(FALSE);
    return(TRUE);
}


/*
 *	BagUniqueAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add an element to the head of a bag.
 *	First check if it is already in the BAG and
 *	simply return if it is.
 *	Simply return if the user tries to add NULL since
 *	NULLs take up no space.
 */
void
BagUniqueAdd( BAG bBag, GENP PObj )
{

    if ( PObj == NULL ) 
	return;

    if ( bBagContains( bBag, PObj ) ) 
	return;

    BagAdd( bBag, PObj );
}







/*
 *      blBagLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Prepare the bag to be looped over.
 *      Calls to BagNext will return the next object in the bag.
 */
BAGLOOP
blBagLoop( BAG bBag )
{
BAGLOOP		blNew;

    blNew.bnPCur = bBag->bnPFirstNode;
    if ( blNew.bnPCur == NULL ) blNew.bnPNext = NULL;
    else			 blNew.bnPNext = (blNew.bnPCur)->bnPNextNode;
    return(blNew);
}




/*
 *      oBagNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the next element in the bag.
 *
 */
GENP
PBagNext( BAGLOOP *blPBagLoop )
{
BNODEP   bnPNode;

    if ( blPBagLoop->bnPCur == NULL ) return(NULL);

		/* Get the current element */

    bnPNode = blPBagLoop->bnPCur;

		/* Advance to the next element */

    blPBagLoop->bnPCur = blPBagLoop->bnPNext;
    if ( blPBagLoop->bnPNext == NULL ) blPBagLoop->bnPNext = NULL;
    else	blPBagLoop->bnPNext =(blPBagLoop->bnPNext)->bnPNextNode;

		/* Return the current element */

    return(bnPNode->PObject);
}







 
