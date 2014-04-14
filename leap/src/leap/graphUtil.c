/*
 *	File:	graphUtil.c
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
 *		This file contains functions that perform
 *		various operations on graphs.  Some of the
 *		functions are finding the shortest paths between
 *		two vertices, finding all of the smallest cycles,
 *		so on and so on.
 */





#include	"basics.h"

#include	"classes.h"

#include	"graphUtil.h"

#include        "sort.h"


/*
 *	bGraphUtilFindShortestPath
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Find the shortest path between two atoms within a UNIT
 *	by generating a spanning tree from the first atom and move
 *	out until you hit the second atom.  Return TRUE if a path
 *	was found, otherwise FALSE.  Call the function fPCallback
 *	for each atom found.
 */
BOOL
bGraphUtilFindShortestPath( UNIT uUnit, ATOM aStart, ATOM aStop, 
		void (*fPCallback)())
{
LOOP		lSpanning;
ATOM		aCur;
BOOL		bFound;

    bFound = FALSE;
    lSpanning = lLoop( (OBJEKT)aStart, SPANNINGTREE );
    while ( (aCur = (ATOM)oNext(&lSpanning)) ) {
	if ( aCur == aStop ) {
	    bFound = TRUE;
	    break;
	}
    }

	/* Now if bFound is TRUE then we found a path and use the */
	/* back pointers that were setup in the spanning tree to */
	/* find the shortest path */

    if ( bFound ) {
	aCur = aStop;
	while ( aCur != aStart ) {
	    fPCallback( aCur );
	    aCur = aAtomBackSpan(aCur);
	}
	fPCallback( aStart );
    }

    return bFound;
}




/*
 *-------------------------------------------------------------------
 *
 *	Code to find rings.
 *
 */


		/* BROKENBONDt is used to keep track of bonds temporarily */
		/* while they are broken during the search for rings */

typedef	struct	{
		ATOM		aAtom1;
		ATOM		aAtom2;
		FLAGS		fFlags;
		} BROKENBONDt;

typedef	struct	{
		ATOM		aAtom;
		BOOL		bInSmallRing;
		} RINGOVERLAPt;

typedef	struct	{
		INTERNAL	inRing;
		int		iSize;
		} RINGSORTt;



/*
 *	zGraphUtilDescribeRing
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Print a description of the ring.
 */
void
zGraphUtilDescribeRing( INTERNAL inRing )
{
ATOM		aCur;
STRING		sTemp;

    PRINTF(( "Ring with %d atoms: ", iInternalRingSize(inRing) ));
    InternalRingLoopAtoms(inRing);
    while ( (aCur = aInternalRingNextAtom(inRing)) ) {
	PRINTF(( "       %s\n", sContainerFullDescriptor((CONTAINER)aCur,sTemp) ));
    }
}


/*
 *	zbGraphUtilSeparate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Seperate two rings.  First check to make sure that one
 *	ring contains a set of ATOMs that is a subset of the
 *	other ring.
 *	This routine will only work if the rings passed are
 *	in order by size and the difference between iBig and iSmall
 *	is a minimal ring.
 *	On exit the ring in iBig will be minimal.
 */
BOOL	
zbGraphUtilSeparate( INTERNAL iBig, int iBigIndex, 
		INTERNAL iSmall, int iSmallIndex )
{
RINGOVERLAPt	*roPFirst, *roPLast, *roPBefore, *roPAfter, *roPCur;
RINGOVERLAPt	*roPA, *roPB, *roPAPrevNotInSmall, *roPAtom;
LOOP		lSpanning;
ATOM		aCur, aAtom;
VARARRAY	vaAtoms;
int		i, iOverlap;
FLAGS		fBond;
#ifdef	DEBUG
STRING		sTemp, sA, sB;
#endif



    MESSAGEEXECUTE( {
			MESSAGE(( "Seperating the following rings:\n" ));
			zGraphUtilDescribeRing(iBig);
			zGraphUtilDescribeRing(iSmall);
			MESSAGE(( "------------\n" ));
		    } );

	/* Allocate an array to keep track of ATOMs in the larger ring */
	/* and whether or not they are also in the smaller ring */

    vaAtoms = vaVarArrayCreate( sizeof(RINGOVERLAPt) );
    VarArraySetSize( vaAtoms, iInternalRingSize(iBig) );

	/* Copy the ATOMs in the larger ring into the VARARRAY */

    InternalRingLoopAtoms(iBig);
    roPAtom = PVAI( vaAtoms, RINGOVERLAPt, 0 );
    for ( i=0; i<iInternalRingSize(iBig); i++ ) {
	aAtom = aInternalRingNextAtom(iBig);
	roPAtom->aAtom = aAtom;
	roPAtom->bInSmallRing = FALSE;
	roPAtom++;
    }

	/* Mark which ATOMs in the big ring are in the small ring */
	/* if there are less than two ATOMs then no seperation */
	/* is required */

    iOverlap = 0;
    InternalRingLoopAtoms(iSmall);
    while ( (aAtom = aInternalRingNextAtom(iSmall)) ) {
	roPAtom = PVAI( vaAtoms, RINGOVERLAPt, 0 );
	for ( i=0; i<iVarArrayElementCount(vaAtoms); i++ ) {
	    if ( aAtom == roPAtom->aAtom ) {
		roPAtom->bInSmallRing = TRUE;
		iOverlap++;
	    }
	    roPAtom++;
	}
    }
    if ( iOverlap < 2 ) return(FALSE);

	/* Now find the two atoms that define the boundary where the */
	/* larger ring seperates from the smaller ring */

    roPFirst  = PVAI(vaAtoms,RINGOVERLAPt,0);
    roPLast   = roPFirst + iVarArrayElementCount(vaAtoms);
    roPBefore = roPLast-1;
    roPCur    = roPFirst;
    roPAfter  = roPFirst+1;
    while ( roPCur < roPLast ) {

		/* Check if ATOM at (i) is one of the boundary ATOMs */

	if ( roPCur->bInSmallRing ) {
	    if ( !roPBefore->bInSmallRing ) {
		roPA = roPCur;
		roPAPrevNotInSmall = roPBefore;
	    }
	    if ( !roPAfter->bInSmallRing )  roPB = roPCur;
	}
		
		/* Increment (i), (iBefore), (iAfter) */

	roPBefore++;
	if ( roPBefore == roPLast ) roPBefore = roPFirst;
	roPCur++;
	roPAfter++;
	if ( roPAfter == roPLast ) roPAfter = roPFirst;
    }

	/* Find the ATOM that immediatly follows roPA->aAtom in the */
	/* ring.  If that ATOM is one of the ATOMs in the shortest */
	/* path between roPA->aAtom and roPB->aAtom then the rings */
	/* do not require seperation. */

    roPAfter = roPA+1;
    if ( roPAfter == roPLast ) roPAfter = roPFirst;

    if ( !roPAfter->bInSmallRing ) {
	DFATAL(( "The ATOM after the first fused ATOM is not fused" ));
    }

    MESSAGE(( "The boundary atoms are: %s and %s\n",
			sAtomName(roPA->aAtom), sAtomName(roPB->aAtom) ));

	/* Now find the shortest path between the two boundary ATOMs */
	/* Generate a path from roPA->aAtom to roPB->aAtom */
	/* Make sure that you only look at ATOMs that are in the small ring */
	/* do this by breaking */

    fBond = fAtomFindBondFlags( roPA->aAtom, roPAPrevNotInSmall->aAtom );
    AtomRemoveBond( roPA->aAtom, roPAPrevNotInSmall->aAtom );
    lSpanning = lLoop( (OBJEKT)roPB->aAtom, SPANNINGTREE );
    MESSAGE(( "Searching for shortest path from %s to %s\n",
		sContainerFullDescriptor( (CONTAINER)roPB->aAtom, sB ),
		sContainerFullDescriptor( (CONTAINER)roPA->aAtom, sA ) ));
    while ( (aCur = (ATOM)oNext(&lSpanning)) ) {
	MESSAGE(( "    found %s\n", 
			sContainerFullDescriptor( (CONTAINER)aCur, sA ) ));
	if ( aCur == roPA->aAtom ) break;
    }
    AtomBondToFlags( roPA->aAtom, roPAPrevNotInSmall->aAtom, fBond );

	/* If the second ATOM in the shortest path is the same as */
	/* roPAfter->aAtom then the rings do not need to be separated */

    if ( roPAfter->aAtom == aAtomBackSpan(roPA->aAtom) ) return(FALSE);

	/* Now remove the ATOMs that are not on the boundary but are */
	/* in the smaller ring from the bigger ring */

    roPCur = roPFirst;
    while ( roPCur < roPLast ) {

		/* Remove all ATOMs that are in the smaller ring */
		/* but not boundary ATOMs from the bigger ring */

	if ( roPCur->bInSmallRing &&
	     roPCur != roPA &&
	     roPCur != roPB ) {
	    if ( bInternalRingRemoveAtom( iBig, roPCur->aAtom ) ) {
		MESSAGE(( "Removed %s from the big ring\n", 
				sContainerFullDescriptor((CONTAINER)roPCur->aAtom, sTemp ) ));
	    } else {
		DFATAL(( "Atom: %s was not in big ring", 
			sAtomName(roPCur->aAtom) ));
	    }
	}
	roPCur++;
    }

	/* Follow the path back and add the new ATOMs to the big ring */
	/* Add the new ATOMs after the */

    aCur = aAtomBackSpan(roPA->aAtom);
    while ( aCur != roPB->aAtom ) {
	MESSAGE(( "Adding %s back to the big ring\n", 
		sContainerFullDescriptor( (CONTAINER)aCur, sTemp ) ));
	InternalRingAddAtomAfter( iBig, aCur, roPB->aAtom );
	aCur = aAtomBackSpan(aCur);
    }

    MESSAGEEXECUTE( {
			MESSAGE(( "Result of seperation:\n" ));
			zGraphUtilDescribeRing(iBig);
			zGraphUtilDescribeRing(iSmall);
			MESSAGE(( "=============\n" ));
		    } );

    return(TRUE);
}




		
	



/*
 *	ziGraphUtilJoinRingGroups
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Join the two ring groups indexed by iGroup1 and iGroup2.
 *	Put the rings of the smaller group into the larger group and
 *	return the new group.  Set the VARARRAY entry for the smaller
 *	group to NULL.
 */
int	
ziGraphUtilJoinRingGroups( VARARRAY vaRingGroup, int iGroup1, int iGroup2 )
{
int		iLarger, iSmaller;
LISTLOOP	llRings;
INTERNAL	inRing;
ATOM		aAtom;
LIST		lGroup1, lGroup2, lLarger, lSmaller;

		/* Find the larger group, the one with more rings */
		/* is more likely to have more ATOMs which all have */
		/* to be changed to the new ring group */

    lGroup1 = *PVAI(vaRingGroup,LIST,iGroup1);
    lGroup2 = *PVAI(vaRingGroup,LIST,iGroup2);

    if ( iListSize(lGroup1) < iListSize(lGroup2) ) {
	iLarger = iGroup2;
	iSmaller = iGroup1;
	lLarger = lGroup2;
	lSmaller = lGroup1;
    } else {
	iLarger = iGroup1;
	iSmaller = iGroup2;
	lLarger = lGroup1;
	lSmaller = lGroup2;
    }

    MESSAGE(( "Concatenating ring group #%d to ring group #%d\n",
		iSmaller, iLarger ));
    MESSAGE(( "Ring group #%d has %d rings.\n",
		iLarger, iListSize(lLarger) ));
    MESSAGE(( "Ring group #%d has %d rings.\n",
		iSmaller, iListSize(lSmaller) ));

		/* Change all of the ATOMs in the smaller ring */
		/* group by putting them in the larger ring group */

    llRings = llListLoop(lSmaller);
    while ( (inRing = (INTERNAL)oListNext(&llRings)) ) {
	InternalRingLoopAtoms(inRing);
	while ( (aAtom = aInternalRingNextAtom(inRing)) ) {
	    AtomSetTempInt( aAtom, iLarger );
	}
    }
    ListConcat( lLarger, lSmaller );
    Destroy( (OBJEKT *) &lSmaller );
    *PVAI(vaRingGroup,LIST,iSmaller) = NULL;

    return(iLarger);
}






/*
 *	GraphUtilFindAllSmallestRings
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Find all of the rings within the UNIT and create INTERNALs
 *	that describe the rings.  The INTERNALs are added to the contents
 *	of the ATOM CONTAINERs that form them.
 *
 *	The rings are found by building a spanning tree from an arbitrary
 *	ATOMs within the seperate molecules of the UNIT.
 *	The molecules are found by building spanning trees.
 *	The rings are found by building a spanning tree and then looking
 *	for bonds that are not part of the spanning tree.  These bonds
 *	are ring closing bonds.  The ring closing bonds are broken
 *	and the shortest paths between the ring closing atoms are found
 *	These paths combined with the ring closing bonds are the rings.
 *
 *	The rings found by the above method are not always the smallest
 *	rings possible, as fused ring assemblies tend to come up not as
 *	the set of smallest rings, but as collections of large rings
 *	which contain smaller rings.  These rings must then be teased
 *	apart.  The seperation of the rings is done by finding all
 *	of the ring groups.  Groups of rings that are fused together.
 *	The fused rings are then seperated by finding the shortest 
 *	path between the two atoms that define the junction between the
 *	two fused rings and taking the ATOMs that the two rings have in
 *	common out of the larger ring and putting the shortest path between
 *	the rings into the larger ring.
 *
 *	Return a VARARRAY that contains LISTs of INTERNALs (rings) where
 *	each LIST contains an INTERNAL ring for each system of
 *	fused rings.
 *
 *
 */
void	
GraphUtilFindAllSmallestRingsAndRingGroups( UNIT uUnit, VARARRAY *vaPRingGroups)
{
LOOP		lAtoms, lTemp;
VARARRAY	vaBrokenBonds;
ATOM		aFirst, aAtom, aA, aB;
BROKENBONDt	bbTemp;
FLAGS		fFlags;
int		h, i, j, iPut;
#ifdef DEBUG
STRING		sA, sB;
#endif
INTERNAL	inRing;
VARARRAY	vaRingGroups;
LIST		lRingGroup;
int		iRingGroupIndex;
BOOL		bNewRingGroup;
LISTLOOP	llRings;
ATOM		aLast;
VARARRAY	vaRingSort;
RINGSORTt*	rsPCur;
INTERNAL	inRingBig, inRingSmall;

		/* Turn off the ATOMTOUCHED flag for all ATOMs in UNIT */
		/* ATOMTOUCHED will be used to determine whether the ATOM */
		/* has been used yet in a spanning tree */


    ContainerResetAllAtomsFlags( (CONTAINER)uUnit, ATOMTOUCHED|ATOMTOUCHED2 );

		/* Maintain VARARRAYs of broken bonds and ring groups */

    vaBrokenBonds = vaVarArrayCreate(sizeof(BROKENBONDt));
    vaRingGroups = vaVarArrayCreate(sizeof(LIST));

		/* Build a spanning tree over the atoms of the UNIT */

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aFirst = (ATOM)oNext(&lAtoms)) ) {
	if ( bAtomFlagsSet( aFirst, ATOMTOUCHED ) ) continue;
	MESSAGE(( "Starting a new molecule at atom: %s\n",
			sContainerFullDescriptor((CONTAINER)aFirst,sA) ));

			/* First build a spanning tree to build all of the */
			/* back pointers, also set the ATOMTOUCHED flag */
			/* to note that the atom has been used in a spanning */
			/* tree */

	lTemp = lLoop( (OBJEKT)aFirst, SPANNINGTREE );
	while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
	    MESSAGE(( "Molecule part: %s   ", 
		sContainerFullDescriptor((CONTAINER)aAtom,sA) ));
	    if ( aAtomBackSpan(aAtom) == NULL ) {
		MESSAGE(( "Points back to nowhere\n" ));
	    } else {
		MESSAGE(( "Points back to: %s\n", 
			sAtomName(aAtomBackSpan(aAtom)) ));
	    }
	    AtomSetFlags( aAtom, ATOMTOUCHED );
	}
    }

			/* Now look through all the bonds, searching */
			/* for those where there are no back pointers */
			/* connecting the two atoms.  These bonds will */
			/* be loop closing bonds.  Save the bonds in the */
			/* vaBrokenBonds array and break the bonds */

    lTemp = lLoop( (OBJEKT)uUnit, BONDS );
    while ( oNext(&lTemp ) ) {
	LoopGetBond( &lTemp, &aA, &aB );
	if ( ! ( aAtomBackSpan(aA) == aB ||
		 aAtomBackSpan(aB) == aA ) ) {
	    MESSAGE(( "Loop closing bond %s - %s\n", 
			sAtomName(aA), sAtomName(aB) ));
	    bbTemp.aAtom1 = aA;
	    bbTemp.aAtom2 = aB;
	    bbTemp.fFlags = fAtomFindBondFlags( aA, aB );
	    VarArrayAdd( vaBrokenBonds, (GENP)&bbTemp );
	}
    }

		/* Break the bonds that were loop closing bonds */

    for ( i=0; i<iVarArrayElementCount(vaBrokenBonds); i++ ) {
	aA = PVAI(vaBrokenBonds,BROKENBONDt,i)->aAtom1;
	aB = PVAI(vaBrokenBonds,BROKENBONDt,i)->aAtom2;

	MESSAGE(( "Breaking bond between %s and %s\n",
			sContainerFullDescriptor((CONTAINER)aA,sA),
			sContainerFullDescriptor((CONTAINER)aB,sB) ));
	AtomRemoveBond( aA, aB );
    }

	/* Now for all of the broken bonds, build a spanning tree */
	/* to search from the first atom of the broken bond to its */
	/* old neighbor */
	/* The path between them will contain atoms that are parts */
	/* of rings, but are not minimal rings, the have to be seperated */
	/* from each other */

    MESSAGE(( "There are %d broken bonds\n", 
		iVarArrayElementCount(vaBrokenBonds) ));
    for ( i=0; i<iVarArrayElementCount(vaBrokenBonds); i++ ) {
	aA = PVAI( vaBrokenBonds, BROKENBONDt, i )->aAtom1;
	aB = PVAI( vaBrokenBonds, BROKENBONDt, i )->aAtom2;

		/* Loop over the atoms in the ring, the shortest bath */
		/* between aA and aB */
		/* This will set the backpointers from (aB) back to (aA) */

	MESSAGE(( "Tracing shortest path between: %s and %s\n",
			sContainerFullDescriptor((CONTAINER)aA,sA), 
			sContainerFullDescriptor((CONTAINER)aB,sB) ));

	lTemp = lLoop( (OBJEKT)aA, SPANNINGTREE );
	while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
	    if ( aAtom == aB ) break;
	}

		/* First pull all of the rings that are fused to this */
		/* ring together into one ring group */
		/* Loop over the ATOMs in the shortest path by following */
		/* the backspan pointers until no more are found */

	iRingGroupIndex = -1;
	for ( aAtom = aB; aAtom != NULL; aAtom = aAtomBackSpan(aAtom) ) {
	    if ( bAtomFlagsSet( aAtom, ATOMTOUCHED2 ) ) {
		if ( iRingGroupIndex != -1 ) {

				/* If we have found a ring that is */
				/* not part of the current ring group */
				/* then join the two ring groups */

		    if ( iAtomTempInt(aAtom) != iRingGroupIndex ) {
			iRingGroupIndex = ziGraphUtilJoinRingGroups( 
						vaRingGroups,
						iAtomTempInt(aAtom),
						iRingGroupIndex );
		    }
		} else iRingGroupIndex = iAtomTempInt(aAtom);
	    }
	}

		/* If no RingGroup was found then create a new one */

	if ( iRingGroupIndex == -1 ) {
	    iRingGroupIndex = iVarArrayElementCount(vaRingGroups);
	    lRingGroup = (LIST)oCreate(LISTid);
	    bNewRingGroup = TRUE;
	} else {
	    lRingGroup = *PVAI(vaRingGroups,LIST,iRingGroupIndex);
	    bNewRingGroup = FALSE;
	}
	    
		/* Now put the atoms into the ring in order of how */
		/* they appear around the ring and put the ring into */
		/* the ring group. */

	inRing = iInternalRing();
	aLast = NULL;
	for ( aAtom = aB; aAtom != NULL; aAtom = aAtomBackSpan(aAtom) ) {
	    MESSAGE(( "Adding %s to ring\n",
				sContainerFullDescriptor((CONTAINER)aAtom,sA) ));
	    InternalRingAddAtomAfter( inRing, aAtom, aLast );
	    aLast = aAtom;
	    AtomSetFlags( aAtom, ATOMTOUCHED2 );
	    AtomSetTempInt( aAtom, iRingGroupIndex );
	}

		/* Add the ring to the ring group */
		/* and if the ring group is a new one then add it to */
		/* the VARARRAY of ring groups */

	ListAdd( lRingGroup, (OBJEKT)inRing );
	if ( bNewRingGroup ) 
		VarArrayAdd( vaRingGroups, (GENP)&lRingGroup );
    }

	/* Rebuild the broken bonds */

    MESSAGE(( "About to rejoin bonds\n" ));
    for ( i=0; i<iVarArrayElementCount(vaBrokenBonds); i++ ) {
	aA = PVAI( vaBrokenBonds, BROKENBONDt, i )->aAtom1;
	aB = PVAI( vaBrokenBonds, BROKENBONDt, i )->aAtom2;
	MESSAGE(( "Rejoining bond between: %s - %s\n",
			sAtomName(aA), sAtomName(aB) ));
	fFlags = PVAI( vaBrokenBonds, BROKENBONDt, i )->fFlags;
	AtomBondToFlags( aA, aB, fFlags );
    }

	/* For debugging purposes print up a summary of the ring groups */

#ifdef DEBUG
    MESSAGE(( "There are total %d ring groups\n", 
		iVarArrayElementCount(vaRingGroups) ));
    for ( i=0; i<iVarArrayElementCount(vaRingGroups); i++ ) {
	if ( *PVAI(vaRingGroups,LIST,i) == NULL ) {
	    MESSAGE(( "Ring group #%d is EMPTY\n", i ));
	} else {
	    lRingGroup = *PVAI(vaRingGroups,LIST,i);
	    MESSAGE(( "Ring group #%d has %d rings\n",
			i, iListSize(lRingGroup) ));
	}
    }
#endif

		/* Tease apart the rings within the ring groups */

    for ( h=0; h<iVarArrayElementCount(vaRingGroups); h++ ) {
	lRingGroup = *PVAI(vaRingGroups,LIST,h);
	if ( lRingGroup != NULL ) {

	    MESSAGE(( "About to sort rings for group: %d\n", h ));

			/* First sort the rings within the ring group */
			/* so that we separate them in order */

	    vaRingSort = vaVarArrayCreate(sizeof(RINGSORTt));
	    VarArraySetSize( vaRingSort, iListSize(lRingGroup) );
	    llRings = llListLoop(lRingGroup);
	    rsPCur = PVAI( vaRingSort, RINGSORTt, 0 );
	    while ( (inRing = (INTERNAL)oListNext(&llRings)) ) {
		rsPCur->inRing = inRing;
		rsPCur->iSize = iInternalRingSize(inRing);
		rsPCur++;
	    }
	    SortByInteger( PVAI(vaRingSort,RINGSORTt,0),
				iVarArrayElementCount(vaRingSort),
				sizeof(RINGSORTt),
				&(PVAI(vaRingSort,RINGSORTt,0)->iSize),
				FALSE );

	    MESSAGE(( "About to separate the rings.\n" ));

	    for ( i=0; i<iVarArrayElementCount(vaRingSort); i++ ) {
		inRingBig = PVAI( vaRingSort, RINGSORTt, i )->inRing;
		for ( j=i+1; j<iVarArrayElementCount(vaRingSort); j++ ) {
		    inRingSmall = PVAI( vaRingSort, RINGSORTt, j )->inRing;
		    if ( zbGraphUtilSeparate( inRingBig, i, inRingSmall, j ) ) break;
		}
	    }
	    VarArrayDestroy( &vaRingSort );
	}
    }

		/* Remove the empty ring groups */
		/* Do this by compacting the (vaRingGroups) VARARRAY */

    iPut = 0;
    for ( i=0; i<iVarArrayElementCount(vaRingGroups); i++ ) {
	if ( !*PVAI(vaRingGroups,LIST,i) ) continue;
	*PVAI(vaRingGroups,LIST,iPut) = *PVAI(vaRingGroups,LIST,i);
	iPut++;
    }
    VarArraySetSize( vaRingGroups, iPut );

    *vaPRingGroups = vaRingGroups;
    VarArrayDestroy(&vaBrokenBonds);

}





/*
 *	GraphUtilDestroyRingGroupVarArray
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the ring group VARARRAY, generated by:
 *	GraphUtilFindAllSmallestRingsAndRingGroups
 */
void	
GraphUtilDestroyRingGroupVarArray( VARARRAY *vaPRingGroup )
{
int		i;
LIST		lRingGroup;

		/* Now destroy the lists, the ring groups */

    for ( i=0; i<iVarArrayElementCount(*vaPRingGroup); i++ ) {
	lRingGroup = *PVAI(*vaPRingGroup,LIST,i);
	if ( lRingGroup != NULL ) 
		Destroy( (OBJEKT *)&lRingGroup );
    }

    VarArrayDestroy(vaPRingGroup);

}




/*
 *	GraphUtilFindAllSmallestRings
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Find all of the smallest rings, and add the
 *	ring INTERNALs to the ATOMs that make up the RINGs.
 */
void
GraphUtilFindAllSmallestRings( UNIT uUnit )
{
VARARRAY	vaRingGroups;

		/* Find the rings and the ring groups */

    GraphUtilFindAllSmallestRingsAndRingGroups( uUnit, &vaRingGroups );

		/* Destroy the ring groups */

    GraphUtilDestroyRingGroupVarArray( &vaRingGroups );
}


