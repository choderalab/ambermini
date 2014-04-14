/*
 *	File:	chirality.c
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
 *		Routines to calculate chirality around an ATOM, and
 *		to calculate the orientation of an ATOM with respect to
 *		a plane formed by three other ATOMs.
 *
 *		Chirality is measured using at least three ATOMs
 *
 *			B     A
 *			 \   /
 *			  \ /
 *			   C
 *			   |
 *			   |
 *			   D	Out or into the plane?
 *
 *		The ATOM (D) can be placed either above or below the plane
 *		defined by (A)-(C)-(B).  The chirality of (C) is calculated
 *		by ordering the ATOMs around it and then crossing (A)-(C) into
 *		(B)-(C) and doing a dot product with (D)-(C).  If the result
 *		is positive then the chirality is +1.0 otherwise -1.0.
 *		The ordering of the ATOMs is determined by sorting them
 *		using the ATOM ID as the key.
 *
 *		Often it is important to calculate the orientation of
 *		an ATOM with respect to three others.  This can be done
 *		by transforming the chirality calculated with respect
 *		to the above ordering into an -orientation- with respect
 *		to a different ordering.
 *
 *		Chirality is defined as the value obtained when the ordering
 *			of the ATOMs is with respect to increasing ATOM ID.
 *		Orientation is defined as the value obtained when the
 *			ordering is arbitrary.
 */




#include	"basics.h"

#include	"classes.h"

#include	"chirality.h"



/*
 *	zChiralityTransformOrientation
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Transform the orientation that has been measured with 
 *	respect to the ordering in aaOrig[4] to the ordering
 *	in aaNew[4].  Return the result.
 *
 *	The transformation is done by swapping ATOMs in aaOrig until
 *	the order matches that of aaNew, each time two ATOMs are swapped,
 *	flip the sign of the orientation.
 *
 *	SIDE EFFECT:   The order in aaOrig is changed.
 */
static void
zChiralityTransformOrientation( double dOrig, ATOM aaOrig[4], 
				double *dPNew, ATOM aaNew[4] )
{
int	i, j;
ATOM	aTemp;

	*dPNew = dOrig;
/*
for ( i=0; i<4; i++ ) {
VP0((" -- orig/new %s %s\n", sContainerName(aaOrig[i]), sContainerName(aaNew[i])
));
}
*/
	for ( i=0; i<4; i++ ) {
		for ( j=i; j<4; j++ ) {
			if ( aaOrig[j] == aaNew[i] ) 
				break;
		}
		if ( j >= 4 ) {
			STRING sOrigDesc[4];
			STRING sNewDesc[4];
			VP0(( "ERROR: Comparing atoms\n"
                  "        %s, \n"
                  "        %s, \n"
                  "        %s, and \n"
                  "        %s \n"
                  "       to atoms\n"
                  "        %s, \n"
                  "        %s, \n"
                  "        %s, and \n"
                  "        %s \n"
			      "       This error may be due to faulty Connection atoms.\n",
				sContainerFullDescriptor( (CONTAINER) aaOrig[0], sOrigDesc[0] ),
				sContainerFullDescriptor( (CONTAINER) aaOrig[1], sOrigDesc[1] ),
				sContainerFullDescriptor( (CONTAINER) aaOrig[2], sOrigDesc[2] ),
				sContainerFullDescriptor( (CONTAINER) aaOrig[3], sOrigDesc[3] ),
				sContainerFullDescriptor( (CONTAINER) aaNew[0],  sNewDesc[0] ),
				sContainerFullDescriptor( (CONTAINER) aaNew[1],  sNewDesc[1] ),
				sContainerFullDescriptor( (CONTAINER) aaNew[2],  sNewDesc[2] ),
				sContainerFullDescriptor( (CONTAINER) aaNew[3],  sNewDesc[3] )  ));
/*
			Describe( cContainerWithin( aaOrig[0] ) ) ;
			Describe( cContainerWithin( aaNew[i] ) ) ;
*/
			DFATAL(( "Atom named %s from %s did not match !\n",
				sContainerName( aaNew[i] ),
				sContainerName( (CONTAINER) cContainerWithin(aaNew[i]) ) ));
		}
		/* Swap elements and flip sign */
		if ( j != i ) {
			SWAP( aaOrig[j], aaOrig[i], aTemp );
			*dPNew = -(*dPNew);
		}
	}
}
                




/*
 *      zaChiralityFindLeastLargerThan
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the atom whose name is canonically less than
 *      all others but larger than aLast.
 */
static ATOM
zaChiralityFindLeastLargerThan( ATOM aAtom, ATOM aLast )
{
ATOM    aCur, aSmall;
int     i;

	aSmall = NULL;
	for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
		aCur = aAtomBondedNeighbor( aAtom, i );
		if ( aLast != NULL ) {
			if ( iAtomId(aLast) >= iAtomId(aCur) ) 
				continue;
		}
		if ( aSmall == NULL ) 
			aSmall = aCur;
		else if ( iAtomId(aCur) < iAtomId(aSmall) ) 
			aSmall = aCur;
	}
	return(aSmall);
}


     


/*
 *      ChiralityOrderNeighbors
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Sort the neighbors of the atom by name and return the first
 *      four.  This is used by the Builder to calculate chirality and 
 *      where to place neighbors.
 */
void
ChiralityOrderNeighbors( ATOM aAtom, 
		ATOM *aPAtomA, ATOM *aPAtomB, ATOM *aPAtomC, ATOM *aPAtomD )
{
int	iCoor;

	*aPAtomA = NULL;
	*aPAtomB = NULL;
	*aPAtomC = NULL;
	*aPAtomD = NULL;

	iCoor = iAtomCoordination(aAtom);
	if ( iCoor == 0 ) 
		return;

	*aPAtomA = zaChiralityFindLeastLargerThan( aAtom, NULL );
	MESSAGE(( "Order atom A=%s\n", sContainerName(*aPAtomA) ));
	if ( iCoor <= 1 ) 
		return;

	*aPAtomB = zaChiralityFindLeastLargerThan( aAtom, *aPAtomA );
	MESSAGE(( "Order atom B=%s\n", sContainerName(*aPAtomB) ));
	if ( iCoor <= 2 ) 
		return;

	*aPAtomC = zaChiralityFindLeastLargerThan( aAtom, *aPAtomB );
	MESSAGE(( "Order atom C=%s\n", sContainerName(*aPAtomC) ));
	if ( iCoor <= 3 ) 
		return;

	*aPAtomD = zaChiralityFindLeastLargerThan( aAtom, *aPAtomC );
	MESSAGE(( "Order atom D=%s\n", sContainerName(*aPAtomD) ));
}






/*
 *      dChiralityForAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the chirality of the atom.
 *      Either +1.0, -1.0,   or 0.0 if not known or undefined.
 *      Currently the criteria for chirality is absolute orientation
 *      of the vectors joining this atom to its neighbors.
 *      The neighbors are assigned an order (currently by sorting their
 *      atom names) eg: A,B,C  and then the vectors to the neighbors
 *      are calculated and operated on to determine the chirality.
 */
double
dChiralityForAtom( ATOM aAtom )
{
ATOM	aAtomA, aAtomB, aAtomC, aAtomD;
double	dChi;
BOOL	bKnowA, bKnowB, bKnowC, bKnowD;
VECTOR	*vPA = NULL, *vPB = NULL, *vPC = NULL, *vPD = NULL;

	if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) ) 
		return(0.0);
    
                /* Only atoms with 3,4 neighbors have chirality */
                
	dChi = 0.0;
	if ( iAtomCoordination(aAtom) == 3 ||
	     iAtomCoordination(aAtom) == 4    ) {
        
		ChiralityOrderNeighbors( aAtom, &aAtomA, &aAtomB, &aAtomC, 
								&aAtomD );

		bKnowA = ( aAtomA!=NULL && 
			   bAtomFlagsSet(aAtomA,ATOMPOSITIONKNOWN) );
		bKnowB = ( aAtomB!=NULL && 
			   bAtomFlagsSet(aAtomB,ATOMPOSITIONKNOWN) );
		bKnowC = ( aAtomC!=NULL && 
			   bAtomFlagsSet(aAtomC,ATOMPOSITIONKNOWN) );
		bKnowD = ( aAtomD!=NULL && 
			   bAtomFlagsSet(aAtomD,ATOMPOSITIONKNOWN) );

		if ( bKnowA ) 
			vPA = &vAtomPosition(aAtomA);
		if ( bKnowB ) 
			vPB = &vAtomPosition(aAtomB);
		if ( bKnowC ) 
			vPC = &vAtomPosition(aAtomC);
		if ( bKnowD ) 
			vPD = &vAtomPosition(aAtomD);

		dChi = dVectorAtomNormalizedChirality( &vAtomPosition(aAtom),
						iAtomCoordination(aAtom),
						vPA, bKnowA,
						vPB, bKnowB,
						vPC, bKnowC,
						vPD, bKnowD );
	}

	return(dChi);
}





/*
 *	dChiralityToOrientation
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Transform the chirality which has been measured with
 *	respect to ATOM ID ordering to an arbitrary ordering.
 *
 */
double
dChiralityToOrientation( double dChirality, ATOM aCenter, 
			ATOM aAtomA, ATOM aAtomB, ATOM aAtomC, ATOM aAtomD )
{
ATOM		aaOrig[4], aaNew[4];
double		dOrient;
BOOL		bNewNull, bOrigNull, bFound;
int		i, j;

    if ( iDoubleCompare(dChirality,0.0)==0 ) return(0.0);

    ChiralityOrderNeighbors( aCenter, &aaOrig[0], &aaOrig[1], 
				&aaOrig[2], &aaOrig[3] );

    aaNew[0] = aAtomA;
    aaNew[1] = aAtomB;
    aaNew[2] = aAtomC;
    aaNew[3] = aAtomD;

		/* aaNew must have the same number of ATOMs as aaOrig */
		/* If there is a missing ATOM in aaNew */
		/* then find which one is missing from the */
		/* aaOrig list */

    bNewNull = ( aaNew[3] == NULL );
    bOrigNull = ( aaOrig[3] == NULL );
    if ( bNewNull && !bOrigNull ) {
	for ( i=0; i<4; i++ ) {
	    bFound = FALSE;
	    for ( j=0; j<3; j++ ) bFound |= (aaOrig[i] == aaNew[j]); 
	    if ( !bFound ) {
		aaNew[3] = aaOrig[i];
		goto CONT;
	    }
	}
    } else if ( !bNewNull && bOrigNull ) {
	DFATAL(( "Only three neighbors around: %s, but orientation has 4\n",
			sContainerName(aCenter) ));
    }

CONT:
    zChiralityTransformOrientation( dChirality, aaOrig, &dOrient, aaNew );

    return(dOrient);
}



/*
 *	dChiralityFromOrientation
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Transform the chirality which has been measured with
 *	respect to some arbitrary ordering to an ordering
 *	based on ATOM ID's.
 *
 */
double
dChiralityFromOrientation( double dOrient, ATOM aCenter, 
			ATOM aAtomA, ATOM aAtomB, ATOM aAtomC, ATOM aAtomD )
{
ATOM		aaOrig[4], aaNew[4];
double		dChirality;
BOOL		bNewNull, bOrigNull, bFound;
int		i, j;

    if ( iDoubleCompare(dOrient,0.0) == 0 ) return(0.0);

    aaOrig[0] = aAtomA;
    aaOrig[1] = aAtomB;
    aaOrig[2] = aAtomC;
    aaOrig[3] = aAtomD;

    ChiralityOrderNeighbors( aCenter, &aaNew[0], &aaNew[1], 
				&aaNew[2], &aaNew[3] );

		/* aaOrig must have the same number of ATOMs as aaNew */
		/* If there is a missing ATOM in aaOrig */
		/* then find which one is missing from the */
		/* aaNew list */

    bNewNull = ( aaNew[3] == NULL );
    bOrigNull = ( aaOrig[3] == NULL );
    if ( !bNewNull && bOrigNull ) {
	for ( i=0; i<4; i++ ) {
	    bFound = FALSE;
	    for ( j=0; j<3; j++ ) bFound |= (aaNew[i] == aaOrig[j]); 
	    if ( !bFound ) {
		aaOrig[3] = aaNew[i];
		goto CONT;
	    }
	}
    } else if ( bNewNull && !bOrigNull ) {
	DFATAL(( "Only three neighbors around: %s, but orientation has 4\n",
			sContainerName(aCenter) ));
    }
CONT:
    zChiralityTransformOrientation( dOrient, aaOrig, &dChirality, aaNew );

    return(dChirality);
}


