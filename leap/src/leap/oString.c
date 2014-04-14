/*
 *      File: oString.c
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
 *              OSTRING
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              An OSTRING is a object which contains strings.
 *              It has all of the same properties as normal C strings
 *              but must be accessed as an OBJEKT.  OSTRINGs are 
 *              not supposed to replace C strings, they are only provided
 *              to maintain consistancy in the command line interface.
 *              Since the command line interpreter will be dealing with
 *              MOLECULEs, UNITs, RESIDUEs, ATOMs, LISTs, and also
 *              integers, doubles, and strings there should be
 *              a consistant way of handling these things.
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
******************************************************************

        Define Public routines here.

*/

/*
 *      osOStringCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create and initialize a new string.
 */
OSTRING
osOStringCreate()
{
OSTRING osNew;

    MALLOC( osNew, OSTRING, sizeof(OSTRINGt) );
    osNew->iMaxLen = 0;
    osNew->sString = NULL;
    return(osNew);
}




/*
 *      OStringDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy a OSTRING.
 */
void
OStringDestroy( OSTRING *osPStr )
{

        /* Deallocate the string if there is one */

    if ( (*osPStr)->sString != NULL ) 
	FREE( (*osPStr)->sString );

    FREE( *osPStr );
    *osPStr = NULL;
}



/*
 *      OStringDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe the OSTRING.
 */
void
OStringDescribe( OSTRING osStr )
{
    VP0(( "STRING (with no reference): '%s'\n", sOString(osStr) ));
}





/*
 *      OStringCopy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Copy the contents of one OSTRING into another.
 */
void
OStringCopy( OSTRING osDest, OSTRING osSrc )
{
        /* First deallocate the string for osDest */


    if ( osDest->sString != NULL ) FREE( osDest->sString );
    osDest->sString = NULL;
    osDest->iMaxLen = 0;

        /* Then MALLOC the same amount of space as osSrc */

 
    if ( osSrc->sString != NULL ) {
        if ( strlen(osSrc->sString) >0 ) {
            MALLOC( osDest->sString, char*, strlen(osSrc->sString)+1 );
            osDest->iMaxLen = strlen(osSrc->sString)+1;
            strcpy( osDest->sString, osSrc->sString );
        }
    }
}





/*
 *      OStringDefine
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Define the contents of an OSTRING using a regular string.
 */
void
OStringDefine( OSTRING osStr, char *sStr )
{
        /* First throw out what was in the OSTRING from before */


    if ( osStr->sString != NULL ) {
        FREE( osStr->sString );
        osStr->sString = NULL;
        osStr->iMaxLen = 0;
    }

        /* MALLOC memory and copy in sStr */

    if ( sStr != NULL ) {
        if ( strlen(sStr) != 0 ) {
            MALLOC( osStr->sString, char*, strlen(sStr)+1 );
            strcpy( osStr->sString, sStr );
            osStr->iMaxLen = strlen(sStr)+1;
        }
    }
}



/*
 *      OStringConcat
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Concatenate two OSTRINGs together.
 */
void
OStringConcat( OSTRING osA, OSTRING osB )
{
int             iSize;

        /* If the first string is empty then just Define it */

    if ( osA->sString == NULL ) {
        OStringDefine( osA, osB->sString );
        return;
    }

        /* If the second string is empty then just return */

    if ( osB->sString == NULL ) return;

        /* Figure out how big the new string will have to be. */

    iSize = strlen( osA->sString ) + strlen( osB->sString ) + 1;
    REALLOC( osA->sString, char*, osA->sString, iSize ); 
    osA->iMaxLen = iSize;
    strcat( osA->sString, osB->sString );
}





/*
 *      osOStringDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a duplicate of the object.
 *	Should only be called by oObjectDuplicate(),
 *	which sets objekt attributes.
 */
OSTRING
osOStringDuplicate( OSTRING os )
{
OSTRING osNew;


    osNew = (OSTRING)oCreate(OSTRINGid);
    OStringDefine( osNew, sOString(os) );
    return(osNew);
}



