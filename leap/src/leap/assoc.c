/*
 *      File: assoc.c
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
 *              ASSOC
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              An ASSOC (short for association) is an object that
 *              stores a short name and an OBJEKT.
 *              It associates the name with the OBJEKT.
 *              This is used by the command line interface to
 *              store objects with the names that are bound to them.
 *
 */



#include	"basics.h"

#include        "classes.h"



/*
-------------------------------------------------------------------

        Define static variables here.
*/




/*
===================================================================

        Define private routines here.
*/


/*
*******************************************************************

        Define public message routines here.
*/




/*
 *      aAssocCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an ASSOC.
 */
ASSOC   
aAssocCreate()
{
ASSOC   aNew;

    MALLOC( aNew, ASSOC, sizeof(ASSOCt) );
    AssocSetName( aNew, "" );
    AssocSetObject( aNew, NULL );
    return(aNew);
}



/*
 *      AssocDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy an ASSOC.
 *      DEREF the OBJEKT within it.
 */
void
AssocDestroy( ASSOC *aPAssoc )
{
    DEREF( oAssocObject(*aPAssoc) );
    FREE( *aPAssoc );
    *aPAssoc = NULL;
}



/*
 *      AssocDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe the ASSOC.
 */
void
AssocDescribe( ASSOC aAssoc )
{
	if ( strlen(sAssocName(aAssoc)) == 0)
		VP0(( " <assoc>: " ));
	else
		VP0(( "%s: ", sAssocName(aAssoc) ));
	Describe( oAssocObject(aAssoc) );
}
