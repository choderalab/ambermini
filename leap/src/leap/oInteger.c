/*
 *      File: oInteger.c
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
 *              OINTEGER        
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              OINTEGER is an object that has the same properties
 *              as normal integers.  It is only declared to
 *              make passing arguments to functions as homogenous as
 *              possible.  It should not be used for any other purpose. 
 *
 *
 */
 

#include	"basics.h"

#include        "classes.h"

/*
 *      oiOIntegerCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an OINTEGER
 */
OINTEGER
oiOIntegerCreate()
{
OINTEGER        oiNew;

    MALLOC( oiNew, OINTEGER, sizeof(OINTEGERt) );

    return(oiNew);
}




/*
 *      OIntegerDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the OINTEGER object.
 */
void
OIntegerDestroy( OINTEGER *oiPInt )
{
    FREE( *oiPInt );
    *oiPInt = NULL;
}




/* 
 *      OIntegerDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe the OINTEGER object.
 */
void
OIntegerDescribe( OINTEGER oiInt )
{
    VP0(( "%d\n", iOInteger(oiInt) ));
}






/*
 *      oiOIntegerDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a duplicate of the object.
 *	Should only be called by oObjectDuplicate(),
 *	which sets objekt attributes.
 */
OINTEGER
oiOIntegerDuplicate( OINTEGER oi )
{
OINTEGER        oiNew;

    oiNew = (OINTEGER)oCreate(OINTEGERid);
    OIntegerSet( oiNew, iOInteger(oi) );
    return(oiNew);
}



