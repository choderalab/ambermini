/*
 *      File:   parmLib.c
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
 *              A PARMLIB maintains an ordered list of PARMSETs
 *              which are searched sequentially for parameters.
 */





#include	"basics.h"

#include        "classes.h"

#include        "parmLib.h"



PARMLIB	GplDefaultParmLib = NULL;



/*
 *================================================================
 *
 *	Private routines
 *
 */

/*
 *	zbParmLibOK
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Check if the variable passed is a valid PARMLIB.
 */
static BOOL
zbParmLibOk( PARMLIB pl )
{
BOOL		bOk;

    bOk = ( pl != NULL );
#if 0
    if ( !bOk ) {
	VP0(( "There are no PARMSETs loaded!\n" ));
    }
#endif
    return(bOk);
}




/*
 *---------------------------------------------------------------
 *
 *	Public routines.
 */



/*
 *      plParmLibCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a PARMLIB
 */
PARMLIB
plParmLibCreate()
{
PARMLIB     plTemp;

    MALLOC( plTemp, PARMLIB, sizeof(PARMLIBt) );

    plTemp->lParmSets = (LIST)oCreate(LISTid);
    return(plTemp);
}


/*
 *      ParmLibDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy a PARMLIB
 */
void
ParmLibDestroy( PARMLIB *plPLib )
{
PARMSET		psSet;

    ParmLibParmSetLoop( *plPLib );
    while( bParmLibNextParmSet( *plPLib, &psSet ) ) {
	ParmSetDestroy( &psSet );
    }
    Destroy( (OBJEKT *)&((*plPLib)->lParmSets) );
    FREE( *plPLib );
    *plPLib = NULL;
}


/*
 *      ParmLibAddParmSet
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a PARMSET to the PARMLIB
 */
void    
ParmLibAddParmSet( PARMLIB plLib, PARMSET psSet )
{
    if ( !zbParmLibOk(plLib) ) 
	return;
    ListAdd( plLib->lParmSets, (OBJEKT)psSet );
}




/*
 *      ParmLibParmSetLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a loop over all PARMSETs in the PARMLIB.
 */
void    
ParmLibParmSetLoop( PARMLIB plLib )
{
    if ( !zbParmLibOk(plLib) ) 
	return;

    plLib->llParmSetLoop = llListLoop( plLib->lParmSets );
}





/*
 *      bParmLibNextParmSet
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return each PARMSET one at a time.
 *      Return FALSE if there are no more PARMSETs.
 */
BOOL
bParmLibNextParmSet( PARMLIB plLib, PARMSET *psPSet )
{
    if ( !zbParmLibOk(plLib) ) 
	return(FALSE);
    *psPSet = (PARMSET)oListNext(&(plLib->llParmSetLoop));
    if ( *psPSet == NULL ) 
	return(FALSE);
    return(TRUE);
}


