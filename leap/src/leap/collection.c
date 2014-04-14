/*
 *      File: collection.c
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
 *      Class: 
 *              COLLECTION
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              A COLLECTION is an abstract class whose sub-classes
 *              have the following properties.
 *              They have an element count.
 *              Their members can be accessed one by one.
 *
 *              Note: this is an ABSTRACT class.
 *                      Objects of type COLLECTION should
 *                      never be created, only use sub-classes
 *                      of COLLECTION.
 *
 *
 */





#include	"basics.h"

#include	"classes.h"

/*
*******************************************************************

        Define public message routines here.
*/



/*
 *      cCollectionCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a new collection.
 */
COLLECTION
cCollectionCreate( int iType )
{
COLLECTION      cCollect;

    cCollect = NULL;
    switch ( iType ) {
        case LISTid:
            cCollect = (COLLECTION)lListCreate();
            break;
        default:
            DFATAL( ("Unknown COLLECTION sub-class\n") );
            break;
    }
    
    return(cCollect);
}



   
/*
 *      CollectionDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy a collection.
 */
void
CollectionDestroy( COLLECTION *cPCollect )
{

    switch ( iObjectType(*cPCollect) ) {
        case LISTid:
            ListDestroy((LIST *)cPCollect);
            break;
        default:
            DFATAL( ("Unknown COLLECTION sub-class\n") );
            break;
    }
}


   
/*
 *      CollectionDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe a collection.
 */
void
CollectionDescribe( COLLECTION cCollect )
{
    switch ( iObjectType(cCollect) ) {
        case LISTid:
            ListDescribe( (LIST)cCollect );
            break;
        default:
            DFATAL( ("Unknown COLLECTION sub-class\n") );
            break;
    }
}



/*
 *      PCollectionLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a pointer to an object that can be used to
 *      loop over the contents of the COLLECTION.
 */
GENP
PCollectionLoop( COLLECTION cCollect )
{
GENP    PLoop;     

    PLoop = NULL;
    switch ( iObjectType( cCollect ) ) {
        case LISTid:
            PLoop = (GENP)llListLoop( (LIST)cCollect );
            break;
        default:
            DFATAL( ("Attempting to loop over invalid collection id: %d",
                        iObjectType(cCollect) ) );
    }
    return(PLoop);
}


/*
 *      oCollectionNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the next object in the collection loop.
 */
OBJEKT
oCollectionNext( COLLECTION cCollect, GENP *PPNode )
{
OBJEKT  oObj;

    oObj = NULL;
    switch ( iObjectType(cCollect) ) {
        case LISTid:
            oObj = oListNext((LISTLOOP *)PPNode);
            break;
        default:
            DFATAL( ("Unknown COLLECTION sub-class\n") );
            break;
    }
    return(oObj);
}






/*
 *      cCollectionDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Send messages to sub-classes to duplicate themselves.
 *	Should only be called by oObjectDuplicate(),
 *	which sets objekt attributes.
 */
COLLECTION
cCollectionDuplicate( COLLECTION cCol )
{
COLLECTION      cNew;

    cNew = NULL;
    switch ( iObjectType(cCol) ) {
        case LISTid:
            cNew = (COLLECTION) lListDuplicate((LIST)cCol);
            break;
        default:
            DFATAL( ("Cannot duplicate unknown COLLECTION subclass id: %d",
                        iObjectType(cCol) ) );
            break;
    }
    
    return(cNew);
}
