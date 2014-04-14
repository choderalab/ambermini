/*
 *	File:	select.c
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
 *		Select different parts of a UNIT.
 */




#include	"basics.h"

#include	"classes.h"

#include	"minimizer.h"

#include	"select.h"

#include        "graphUtil.h"

#include        "build.h"

/*
 *	SelectAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Just select the ATOM.
 */
void
SelectAtom( ATOM aAtom, BOOL bOn )
{
    if ( bOn )
	AtomSetFlags( aAtom, ATOMSELECTED );
    else
	AtomResetFlags( aAtom, ATOMSELECTED );
}



/*
 *	bSelectRingWithAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Select the RINGs that contain the ATOM.
 *	If there are no RINGs then return FALSE, otherwise TRUE.
 */
BOOL
bSelectRingWithAtom( UNIT uUnit, ATOM aAtom, BOOL bOn )
{
BOOL		bFoundOne;
ATOM		aCur;
INTERNAL	iRing;
LOOP		lInternals, lAtoms;

    bFoundOne = FALSE;
    GraphUtilFindAllSmallestRings( uUnit );

	/* Now check if the ATOM contains ring INTERNALs */

    lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
    while ( (iRing = (INTERNAL)oNext(&lInternals)) ) {
	if ( iInternalType(iRing) == INTERNALRING ) {
	    bFoundOne = TRUE;
	    InternalRingLoopAtoms(iRing);
	    while ( (aCur = aInternalRingNextAtom(iRing)) ) {
		if ( bOn ) AtomSetFlags( aCur, ATOMSELECTED );
		else	AtomResetFlags( aCur, ATOMSELECTED );
	    }
	}
    }

	/* Destroy the INTERNALS */

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    BuildDestroyInternals(&lAtoms);

    return(bFoundOne);
}





/*
 *	SelectResidueWithAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Select the RESIDUE that contains the ATOM.
 */
void	
SelectResidueWithAtom( UNIT uUnit, ATOM aAtom, BOOL bOn )
{
LOOP		lAtoms;
ATOM		aCur;

    lAtoms = lLoop( (OBJEKT)cContainerWithin((CONTAINER)aAtom), ATOMS );
    while ( (aCur = (ATOM)oNext(&lAtoms)) ) {
	if ( bOn ) AtomSetFlags( aCur, ATOMSELECTED );
	else	AtomResetFlags( aCur, ATOMSELECTED );
    }
}


/*
 *	SelectMoleculeWithAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Select the molecule that contains the ATOM.
 *	Do this by generating a spanning tree.
 */
void
SelectMoleculeWithAtom( UNIT uUnit, ATOM aAtom, BOOL bOn )
{
LOOP		lSpan;

    lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
    while ( (aAtom = (ATOM)oNext(&lSpan)) ) {
	if ( bOn ) AtomSetFlags( aAtom, ATOMSELECTED );
	else	AtomResetFlags( aAtom, ATOMSELECTED );
    }
}



/*
 *	SelectEverything
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Select all of the ATOMs in the UNIT.
 */
void
SelectEverything( UNIT uUnit, BOOL bOn )
{
LOOP		lAtoms;
ATOM		aAtom;

    MESSAGE(( "Selecting everything.  Select=%s\n", sBOOL(bOn) ));
    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
	if ( bOn ) AtomSetFlags( aAtom, ATOMSELECTED );
	else	AtomResetFlags( aAtom, ATOMSELECTED );
    }
}




/*
 *	bSelectChainBetween
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Select the shortest chain of ATOMs between aA and aB.
 *	Return FALSE if no chain was found.
 */
BOOL
bSelectChainBetween( UNIT uUnit, ATOM aA, ATOM aB, BOOL bOn )
{
LOOP		lSpan;
ATOM		aCur;

		/* Span out from aA looking for aB */

    lSpan = lLoop( (OBJEKT)aA, SPANNINGTREE );
    while ( (aCur = (ATOM)oNext(&lSpan)) ) {
	if ( aCur == aB ) break;
    }
    if ( aCur == NULL ) return(FALSE);

		/* Select the chain back to aA */

    for ( aCur = aB; aCur; aCur = aAtomBackSpan(aCur) ) {
	if ( bOn ) AtomSetFlags( aCur, ATOMSELECTED );
	else	AtomResetFlags( aCur, ATOMSELECTED );
    }
    return(TRUE);
}




 
 
/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 *
 *	Act on the selected ATOMs
 *
 */







/*
 *	SelectRelaxInFramework
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Relax all of the selected ATOMs within the framework of
 *	the un-selected ATOMs that they are bound to.
 *	Use the MINIMIZER provided by the caller.
 */
void
SelectRelaxInFramework( UNIT uUnit, MINIMIZER mMinimizer )
{
	/* First set the ATOMNEEDSMINIMIZER flag for all selected atoms */


    ContainerWithFlagsSetAtomFlags( (CONTAINER)uUnit, ATOMSELECTED,
					ATOMNEEDSMINIMIZER );

	/* Relax strain */

    BuildRelaxInFramework( uUnit, mMinimizer );

	/* Clean up */

    ContainerWithFlagsResetAtomFlags( (CONTAINER)uUnit, ATOMSELECTED,
					ATOMNEEDSMINIMIZER );

}


