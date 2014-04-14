/*
 *	File:	list.c
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
 *		List is a doubly linked list of OBJEKTS.
 *
 *	WARNING:
 *		NEVER NEVER NEVER remove OBJEKTs from
 *		a LIST that you are currently looping over.
 *		This can lead to inconsistancies.
 *		The solution is to change the structure 
 *		LISTLOOP so that it contains a pointer
 *		to the next element in the list instead
 *		of just to the current one which is to
 *		be Removed from the list.
 *TODO:	Change LISTLOOP to contain the pointer to the next
 *TODO:	element in the list.
 */




#include	"basics.h"

#include	"classes.h"


/*
 *-------------------------------------------------------------------------
 *
 *	Private routines
 */



/*
 *	searchList
 *
 *	Find the node that points to the object pointed to by the
 *	argument.  Return a pointer to the pointer to the node
 *	which points to the object.
 */
static NODEP *
searchList( NODEP *nPPList, OBJEKT oObj, NODEP *nPPPrevious )
{
NODEP   nPCur, *nPPPrev;

        nPPPrev = nPPList;
        nPCur = (*nPPPrev);
	*nPPPrevious = NULL;
                /* search for the node that contains the object */
        while ( nPCur != NULL ) {
                if ( nPCur->PObject == oObj ) break;
		*nPPPrevious = nPCur;
                nPPPrev = &(nPCur->nPNextNode);
                nPCur = (*nPPPrev);
        }

        if ( nPCur == NULL ) return(NULL);
        return(nPPPrev);
}

        


/*
 *-----------------------------------------------------------------
 *
 *        Public routines
 *        
 */


/*
 *	lListCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a new list.
 */
LIST
lListCreate()
{
LIST    lList;

        MALLOC( lList, LIST, sizeof(LISTt) );
        CollectionSetSize( lList, 0 );
        lList->nPFirstNode = NULL;
        lList->nPLastNode = NULL;
        return(lList);
}


/*
 *      ListDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the list.  The caller is responsible for destroying the
 *      contents of the list.
 */
void
ListDestroy( LIST *lPList )
{
NODEP   nPTempPtr, nPFree;

        /* Destroy the LIST nodes */

    nPTempPtr = (*lPList)->nPFirstNode;
    while ( nPTempPtr != NULL ) {
        nPFree = nPTempPtr;
        nPTempPtr = nPTempPtr->nPNextNode;
        
                /* Free the node memory. */
        FREE(nPFree);
    }

        /* Now destroy the LIST itself */

    FREE( *lPList );
    *lPList = NULL;

}



/*
 *        ListAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *        Add an element to the head of a list.
 *        Simply return if the user tries to add NULL since
 *        NULLs take up no space.
 */
void
ListAdd( LIST lList, OBJEKT oObj )
{
NODEP   nPTempList;

    if ( oObj == NULL ) 
	return;

                /* Insert the new node at the head of the list */

    nPTempList = lList->nPFirstNode;
    MALLOC( lList->nPFirstNode, NODEP, sizeof(NODE) );

                /* Keep up to date the last node of the list */

    if ( nPTempList == NULL ) {
        lList->nPLastNode = lList->nPFirstNode;
    }

    (lList->nPFirstNode)->PObject = oObj;
    (lList->nPFirstNode)->nPNextNode = nPTempList;
    CollectionSetSize( lList, iCollectionSize(lList)+1 );

                /* Reference the object because it is being added */
    REF( oObj );
}





/*
 *      ListAddToEnd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add an element to the tail of a list.
 *      Simply return if the user tries to add NULL since
 *      NULLs take up no space.
 */
void
ListAddToEnd( LIST lList, OBJEKT oObj )
{
NODEP   nPTempList;

    if ( oObj == NULL ) 
	return;

                /* Insert the new node at the tail of the list */

    nPTempList = lList->nPLastNode;
    MALLOC( lList->nPLastNode, NODEP, sizeof(NODE) );

                /* Keep up to date the last node of the list */

    if ( nPTempList == NULL ) {
        lList->nPFirstNode = lList->nPLastNode;
    } else {
        nPTempList->nPNextNode = lList->nPLastNode;
    }

    (lList->nPLastNode)->PObject = oObj;
    (lList->nPLastNode)->nPNextNode = NULL;
    CollectionSetSize( lList, iCollectionSize(lList)+1 );

                /* Reference the object because it is being added */
    REF( oObj );
}




/*
 *	ListAddUnique
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add the object to the LIST only if it isn't already in
 *	there.
 */
void
ListAddUnique( LIST lList, GENP PData )
{
    if ( bListContains( lList, PData ) ) 
	return;
    ListAddToEnd( lList, (OBJEKT)PData );
}





/*
 *	ListConcat
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add the contents of lList2 to lList1.
 *	On returning, lList2 will be empty and should be destroyed.
 */
void
ListConcat( LIST lList1, LIST lList2 )
{
    VERIFYOBJEKT( lList1, LISTid );
    VERIFYOBJEKT( lList2, LISTid );

    MESSAGE(( "Before concat List1 size = %d,    List2 size = %d\n",
		iCollectionSize(lList1),
		iCollectionSize(lList2) ));

    if ( iListSize(lList2) == 0 ) return;

	/* Make the LIST larger */

    CollectionSetSize( lList1, iCollectionSize(lList1)+
				iCollectionSize(lList2) );

	/* Point the last node of lList1 to the first node of lList2 */

    if ( lList1->nPLastNode != NULL ) {
	lList1->nPLastNode->nPNextNode = lList2->nPFirstNode;
    } else {
	lList1->nPFirstNode = lList2->nPFirstNode;
    }

    lList1->nPLastNode = lList2->nPLastNode;

	/* Make List2 empty */

    lList2->nPFirstNode = NULL;
    lList2->nPLastNode = NULL;
    CollectionSetSize( lList2, 0 );

    MESSAGE(( "After concat List1 size = %d,    List2 size = %d\n",
		iCollectionSize(lList1),
		iCollectionSize(lList2) ));

    
}




/*
 *      ListDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the string with a description of the list.
 *      It is the callers responsibility to ensure that
 *      the string is long enough to store the description.
 */
void
ListDescribe( LIST lList )
{
OBJEKT          oObj;
int             iSize;
LISTLOOP        llL;

    iSize = iCollectionSize(lList);
    VP0(( "List size=%d\n", iSize ));
    llL = (LISTLOOP)PCollectionLoop( (COLLECTION)lList );
    if ( llL == NULL ) 
	return;
    while ( (oObj = oCollectionNext( (COLLECTION)lList, (GENP *)&llL )) 
								!= NULL ) {
        Describe( oObj );
    }
    VP0(( "--End of list\n" ));
}




        
/*
 *      bListRemove
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Remove an object from the list.
 *      DEREF the object, possible Destroying it if there are
 *      no more references to it.
 *      
 *      Return:
 *              FALSE if the element is not in the list.
 */
BOOL
bListRemove( LIST lList, GENP PPtr )
{
NODEP   *nPPPrev;
NODEP   nPNode;
NODEP	nPPrevious;

                /* If the list is empty then return */

    if ( lList->nPFirstNode==NULL ) return(FALSE);

                /* Search the list for the object */
    nPPPrev = searchList( &(lList->nPFirstNode), (OBJEKT)PPtr, &nPPrevious );
    if ( nPPPrev == NULL ) return(FALSE);

                /* Now actually remove the object from the list */
    nPNode = *nPPPrev;
    (*nPPPrev) = (*nPPPrev)->nPNextNode;

                /* If the node is the last node then update the last node */

    if ( nPNode == lList->nPLastNode ) {
	if ( *nPPPrev == NULL ) lList->nPLastNode = nPPrevious;
	else 			lList->nPLastNode = *nPPPrev;
    }

                /* Destroy the node itself */
    FREE(nPNode);

    CollectionSetSize( lList, iCollectionSize(lList) - 1 );

                /* DEREF the object */

    DEREF( PPtr );

    return(TRUE);
}





        
/*
 *      bListContains
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Find out if a list contains a particular object.
 *      
 *      Return:
 *              FALSE if the element is not in the list.
 */
BOOL
bListContains( LIST lList, GENP PPtr )
{
NODEP   *nPPPrev;
NODEP	nPPrevious;

                /* If the list is empty then return */

    if ( lList->nPFirstNode==NULL ) return(FALSE);

                /* Search the list for the object */
    nPPPrev = searchList( &(lList->nPFirstNode), (OBJEKT)PPtr, &nPPrevious );
    if ( nPPPrev == NULL ) return(FALSE);
    return(TRUE);
}




/*
 *      llListLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Prepare the list to be looped over.
 *      Calls to ListNext will return the next object in the list.
 */
LISTLOOP
llListLoop( LIST lList )
{
    if ( lList == NULL )
	DFATAL(( "llListLoop called with NULL list\n" ));
    return(lList->nPFirstNode);
}




/*
 *      oListNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the next element in the list.
 *
 */
OBJEKT
oListNext( LISTLOOP *llPListLoop )
{
NODEP   nPNode;

    if ( *llPListLoop == NULL ) 
	return(NULL);
    nPNode = *llPListLoop;
    *llPListLoop = (*llPListLoop)->nPNextNode;
    return(nPNode->PObject);
}







/*
 *      lListDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a duplicate of the list.
 *	Should only be called by the collectionduplicate
 *	routine (TODO - delete collection), which is
 *	called by oObjectDuplicate(),
 *	which sets objekt attributes.
 */
LIST
lListDuplicate( LIST lOld )
{
OBJEKT          lNew, oObj, oNew;
LISTLOOP        llLoop;

    lNew = oCreate(LISTid);
    llLoop = llListLoop(lOld);
    while ( ( oObj = oListNext(&llLoop) ) != NULL ) {
        oNew = oObjectDuplicate(oObj);
        ListAdd( (LIST)lNew, oNew );
	oNew->iReferences = 1;	/* since ListAdd() increments */
    }
    return((LIST)lNew);
}
    
