/*
 *	File:	restraint.c
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
 *		Maintain RESTRAINTs.
 */



#include	"basics.h"

#include	"classes.h"
#include	"restraint.h"



/*
 *	rRestraintCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a pointer to a RESTRAINT of undefined type.
 */
RESTRAINT
rRestraintCreate()
{
RESTRAINT	rNew;

    MALLOC( rNew, RESTRAINT, sizeof(RESTRAINTt) );
    rNew->iType = RESTRAINTNONE;
    RestraintDefineFlags( rNew, RESTRAINTALWAYS );
    return(rNew);
}




/*
 *	RestraintDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	FREE the memory allocated to the RESTRAINT.
 */
void
RestraintDestroy( RESTRAINT *rPRes )
{
    FREE( *rPRes );
    *rPRes = NULL;
}





/*
 *	bRestraintBondMatchAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return TRUE if the RESTRAINT bond has the
 *	same atoms as the user passed in.
 */
BOOL
bRestraintBondMatchAtoms( RESTRAINT rRes, ATOM aAtom1, ATOM aAtom2 )
{
    if ( (rRes->rType.rbBond.aAtom1 == aAtom1 &&
	  rRes->rType.rbBond.aAtom2 == aAtom2 ) ||
	 (rRes->rType.rbBond.aAtom1 == aAtom2 &&
	  rRes->rType.rbBond.aAtom2 == aAtom1 ) ) return(TRUE);
    return(FALSE);
}

	
/*
 *	bRestraintAngleMatchAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return TRUE if the RESTRAINT angle has the
 *	same atoms as the user passed in.
 */
BOOL
bRestraintAngleMatchAtoms( RESTRAINT rRes, 
	ATOM aAtom1, ATOM aAtom2, ATOM aAtom3 )
{

    if ( (rRes->rType.raAngle.aAtom1 == aAtom1 &&
          rRes->rType.raAngle.aAtom2 == aAtom2 &&
	  rRes->rType.raAngle.aAtom3 == aAtom3 ) ||
	 (rRes->rType.raAngle.aAtom1 == aAtom3 &&
	  rRes->rType.raAngle.aAtom2 == aAtom2 &&
	  rRes->rType.raAngle.aAtom3 == aAtom1 ) ) return(TRUE);
    return(FALSE);
}



/*	
 *	bRestraintTorsionMatchAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return TRUE if the RESTRAINT torsion has the
 *	same atoms as the user passed in.
 */
BOOL
bRestraintTorsionMatchAtoms( RESTRAINT rRes, 
	ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, ATOM aAtom4 )
{
    if ( (rRes->rType.rtTorsion.aAtom1 == aAtom1 &&
          rRes->rType.rtTorsion.aAtom2 == aAtom2 &&
	  rRes->rType.rtTorsion.aAtom3 == aAtom3 &&
	  rRes->rType.rtTorsion.aAtom4 == aAtom4 ) ||
	 (rRes->rType.rtTorsion.aAtom1 == aAtom4 &&
	  rRes->rType.rtTorsion.aAtom2 == aAtom3 &&
	  rRes->rType.rtTorsion.aAtom3 == aAtom2 &&
	  rRes->rType.rtTorsion.aAtom4 == aAtom1 ) ) return(TRUE);
    return(FALSE);
}



/*
 *	bRestraintContainsAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Test to see if the atom is one of the atoms of
 *	the RESTRAINT.
 *	If it is return TRUE.
 */
BOOL
bRestraintContainsAtom( RESTRAINT rRes, ATOM aAtom )
{
    switch ( rRes->iType ) {
	case RESTRAINTTORSION:
	    if ( rRes->rType.rtTorsion.aAtom1 == aAtom ||
		 rRes->rType.rtTorsion.aAtom2 == aAtom ||
		 rRes->rType.rtTorsion.aAtom3 == aAtom ||
		 rRes->rType.rtTorsion.aAtom4 == aAtom ) return(TRUE);
	    return(FALSE);
	case RESTRAINTANGLE:
	    if ( rRes->rType.raAngle.aAtom1 == aAtom ||
		 rRes->rType.raAngle.aAtom2 == aAtom ||
		 rRes->rType.raAngle.aAtom3 == aAtom ) return(TRUE);
	    return(FALSE);
	case RESTRAINTBOND:
	    if ( rRes->rType.rbBond.aAtom1 == aAtom ||
		 rRes->rType.rbBond.aAtom2 == aAtom ) return(TRUE);
	    return(FALSE);
	default:
	    DFATAL(( "Illegal RESTRAINT type!" ));
    }
    return(FALSE);	/* for lint */
}



