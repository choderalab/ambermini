/*
 *      File:   octree.c
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
 *
 *  octree.c by Bill Ross
 *
 *	Notes
 *
 *	What is a an octree? Take a cube, divide it into 8 cubes, divide
 *	these into 8, down to the desired limit of resolution. Octrees 
 *	are a simple way to deal with space using minimal resources to 
 *	achieve the desired resolution - the subdivision of a given cube 
 *	can be avoided if the cube can be shown to be homogenous.
 *
 *	Recursive programming can be used to elegantly traverse tree-type 
 *	data structures such as octrees. The basic concept is that a
 *	routine for processing a node calls itself for processing any
 *	sub-nodes. So that recursion does not extend infinitely, a
 *	recursive routine will test some condition or conditions to
 *	determine if the end of recursion is reached. This technique
 *	hasn't been used in Fortran very much because Fortran did not
 *	originally allow multiple calls to a single routine in the
 *	stack, because local variables were not stack-based.
 */

#include "classes.h"
#include "octree.h"
#include <math.h>
#include <float.h>
#include <sys/types.h>
#include <signal.h>

time_t	time_start;

/*	status of octree node	*/

#define OCT_UNKNOWN	-1
#define OCT_INCLUDED	1
#define OCT_EXCLUDED	2
#define OCT_PARTIAL	3

/* 
 *  globals used for octree evaluation
 */

double	dGridSize;
double	dShellRadius;
int	iColorMethod;

int	iMaxDepth, iTreeGridPoints, iNodeCount;
float	fVolume;
int	*PiDensities;
double	*PdHalfEdges, *PdHalfDiagonals;

int	iDistanceCharge;
int	iChargeAtoms;
double	dCutRadius;
float	*PfCharges;
ATOM	*PaChargeAtoms;
VECTOR	*PvAtomCrds;
double	*PdAtomChgs;
float	fMinCharge, fMaxCharge;
VECTOR	vMinCharge, vMaxCharge;

VECTOR	vNewPoint;
double	dDeleteRadius;
float	fNewCharge;

ATOM	aClosestAtom;
double	dClosestDistance;

FILE	*fChargeFile;

double
dDistanceSq( VECTOR *Pv1, VECTOR *Pv2 )
{
double  x, y, z;
        x = Pv1->dX - Pv2->dX;
        x = x * x;
 
        y = Pv1->dY - Pv2->dY;
        y = y * y;
 
        z = Pv1->dZ - Pv2->dZ;
        z = z * z;
 
        return(x + y + z);
}

double
dDistance( VECTOR *Pv1, VECTOR *Pv2 )
{
double  x, y, z;

        x = Pv1->dX - Pv2->dX;
        x = x * x;
 
        y = Pv1->dY - Pv2->dY;
        y = y * y;
 
        z = Pv1->dZ - Pv2->dZ;
        z = z * z;
 
        return(sqrt(x + y + z));
}

static void
dumpNode( OCTNODEt *PonNode )
{
int	i;
	fprintf(stderr, "@");
	for (i=0; i<PonNode->iDepth; i++)
		fprintf(stderr, " ");
	fprintf(stderr, "%d %d %f %f %f\n", 
		PonNode->iDepth, PonNode->iStatus,
		PonNode->vCorner.dX, PonNode->vCorner.dY, PonNode->vCorner.dZ);
	if ( PonNode->iDepth == iMaxDepth  &&  PonNode->iStatus == OCT_PARTIAL) {
		fprintf(stderr, " bad node\n");
		kill(0,5);
	}
	if ( PonNode->iStatus == OCT_PARTIAL  &&  PonNode->PonChildren == NULL ) {
		fprintf(stderr, "partial w/out children\n");
		kill(0,5);
	}
	if ( PonNode->PonChildren != NULL ) {
		if ( PonNode->iDepth >= iMaxDepth ) {
			fprintf(stderr, "children too deep\n");
			kill(0,5);
		}
		for (i=0; i<8; i++) 
			dumpNode( &PonNode->PonChildren[i] );
	}
}

void
dumpOctree( OCTREE octTree )
{
	iMaxDepth = octTree->iMaxDepth;
	dumpNode( &octTree->onHead );
}

/*************************************************************************
 *************************************************************************/
static OCTNODEt *
PonMakeChildren( VECTOR *PvCorner, int iDepth, int iStatus )
{
OCTNODEt	*PonChildren;
int		i;
double		dHalfEdge = PdHalfEdges[iDepth-1];
	/*
	 *  Make nodes & set simple stuff by initializing 1st node
	 *	and copying it.
	 */
	MALLOC( PonChildren, OCTNODEt *, 8 * sizeof(OCTNODEt) );
	memset(PonChildren, 0, sizeof(OCTNODEt));
	PonChildren->iStatus = iStatus;
	PonChildren->iDepth = iDepth;
	PonChildren->PonChildren = NULL;
	PonChildren->PaAtomList = NULL;
	PonChildren->PfCharges = NULL;
	PonChildren->iAtoms = 0;
	PonChildren->vCorner = *PvCorner;
	PonChildren->iNodeNum = iNodeCount++;
	for (i=1; i<8; i++) {
		PonChildren[i] = PonChildren[0];
		PonChildren[i].iNodeNum = iNodeCount++;
	}

	/* 
	 *  Differentiate the non-0-corner children by giving them
	 *	seperate corners of the octant. The order is, do
	 *	Z-stacks along the Y-axis, forming slices on the X-axis.
	 *	(This corresponds to the traversal with inner loop of Z,
	 *	middle loop of Y, and outer loop of X). This order is
	 *	important later for charge redistribution when subdividing
	 *	an existing node to delete a new atom's volume from the tree.
	 */
	PonChildren[1].vCorner.dZ += dHalfEdge;		/*      0 + Z 	   */
	PonChildren[2].vCorner.dY += dHalfEdge;		/*      0 + Y 	   */

	PonChildren[3].vCorner = PonChildren[2].vCorner;
	PonChildren[3].vCorner.dZ += dHalfEdge;		/*    0 + Y + Z    */

	PonChildren[4].vCorner.dX += dHalfEdge;		/*      0 + X 	   */

	PonChildren[5].vCorner = PonChildren[4].vCorner;
	PonChildren[5].vCorner.dZ += dHalfEdge;		/*    0 + X + Z	   */
	
	PonChildren[6].vCorner = PonChildren[4].vCorner;
	PonChildren[6].vCorner.dY += dHalfEdge;		/*    0 + X + Y	   */

	PonChildren[7].vCorner = PonChildren[6].vCorner;
	PonChildren[7].vCorner.dZ += dHalfEdge;		/*  0 + X + Y + Z  */

	return( PonChildren );
}

static int
FinalCheck( OCTNODEt *PonNode, int iAtoms, ATOM *PaAtomList )
{
int	i;
ATOM	*PaAtom;
#ifdef OCTDEBUG
	depth[PonNode->iDepth]++;
#endif
	/*
	 *  Only need to check corner point. Must be in at least
	 *	one outer radius and in no inner radius.
	 */
	PonNode->iStatus = OCT_EXCLUDED;
	PaAtom = PaAtomList;
	for (i=0; i<iAtoms; i++, PaAtom++) {
		double	d;
		d = dDistance( &PonNode->vCorner, &vAtomPosition( *PaAtom ) );
		if ( d < dShellRadius )
			PonNode->iStatus = OCT_INCLUDED;
		if ( d < dAtomTemp( *PaAtom ) ) {
#ifdef OCTDEBUG
			insex[oct->depth]++;
#endif
			PonNode->iStatus = OCT_EXCLUDED;
			break;
		}
	}
#ifdef OCTDEBUG
	if (PonNode->iStatus == OCT_INCLUDED)
		included[oct->depth]++;
	else
		outside[oct->depth]++;
#endif
	return(PonNode->iStatus);
}

/*
 *  DestroyOctant() - destroys any type of octree
 */
static void
DestroyOctant( OCTNODEt *PonNode, int iStatus )
{

#ifdef DBG2
	fprintf(stderr, "@@@ destroy node 0x%x level %d type %d -> %d", 
			PonNode, PonNode->iLevel,  PonNode->iStatus, iStatus);
	if ( PonNode->PaAtomList != NULL )
		fprintf(stderr, " atomlist 0x%x", PonNode->PaAtomList);
	else
		fprintf(stderr, " atomlist null");
	if ( PonNode->PonChildren != NULL )
		fprintf(stderr, " children 0x%x\n",  PonNode->PonChildren);
	else
		fprintf(stderr, " children null\n");
#endif

	if ( PonNode->PaAtomList != NULL )
		FREE( PonNode->PaAtomList );

	if ( PonNode->PonChildren != NULL ) {
		int		i;

		/*
		 *  Destroy children. Recursive destroys completely
		 *	delete nodes, so mark OCT_UNKNOWN to avoid
		 *	confusion..
		 */
		for (i=0; i<8; i++)
			DestroyOctant( &PonNode->PonChildren[i], OCT_UNKNOWN );
		FREE( PonNode->PonChildren );
	}

	PonNode->iStatus = iStatus;

}
static int
iBuildShellOctant( OCTNODEt *PonNode, int iAtoms, ATOM *PaAtomList )
{
int		i, iPartialIn, iPartialOut, iIncluded, iExcluded, iNewAtoms;
ATOM		*PaAtom, *PaNewAtoms;
double		d, dHalfEdge, dHalfDiagonal;
VECTOR		vCenter;
OCTNODEt	*PonChildren;


	PonNode->PaAtomList = NULL;
	PonNode->iAtoms = 0;

	/*
	 *  End recursion if the limit of resolution (octree max depth)
	 *	has been reached.
	 */

	if ( PonNode->iDepth == iMaxDepth ) 
		return( FinalCheck( PonNode, iAtoms, PaAtomList ) );

	/*
	 *  At less than maximum depth. Evaluate whole box with respect 
	 *	to atoms by considering overlaps with single atoms' inner 
	 *	and outer radii. Use center of box +- half diagonal for test.
	 *	Optimally, whole box can be included or excluded. Otherwise,
	 *	it will be divided into 8 and the parts evaluated.
	 */

	dHalfEdge = PdHalfEdges[PonNode->iDepth];
	vCenter = PonNode->vCorner;
	vCenter.dX += dHalfEdge;
	vCenter.dY += dHalfEdge;
	vCenter.dZ += dHalfEdge;
	dHalfDiagonal = PdHalfDiagonals[PonNode->iDepth];

	/*
	 *  Get memory for list of atoms overlapping this cube, which 
	 *	is necessarily <= size of parent's.
	 */

	MALLOC( PaNewAtoms, ATOM *, iAtoms * sizeof(ATOM) );
	memset( PaNewAtoms, 0, iAtoms * sizeof(ATOM) );
	PonNode->PaAtomList = PaNewAtoms;
	iNewAtoms = 0;
#ifdef DBG2
	fprintf(stderr, "@@@ node 0x%x atoms 0x%x\n", PonNode, PaNewAtoms);
#endif

	/*
	 *  Consider all parent's atoms in relation to this box.
	 */
	iIncluded = 0;
	iPartialIn = 0;
	iPartialOut = 0;
	PaAtom = PaAtomList;
	for (i=0; i<iAtoms; i++, PaAtom++) {

		d = dDistance( &vCenter, &vAtomPosition( *PaAtom ) );

		/*
		 *  Consider inner radius first, since complete
		 *	inclusion in it means the box is excluded
		 *	and the outer radius doesn't need to be looked at.
		 */
		if ( d + dHalfDiagonal < dAtomTemp( *PaAtom ) ) {
			/*
			 *  Complete inclusion in inner radius of
			 *	this atom: no need to ever look at 
			 *	this box again.
			 */
#ifdef OCTDEBUG
			insex[oct->depth]++;
			depth[oct->depth]++;
#endif
#ifdef DBG2
			fprintf(stderr, "@@@ node 0x%x free atoms 0x%x\n", 
					PonNode, PonNode->PaAtomList);
#endif
			FREE( PaNewAtoms );
			PonNode->PaAtomList = NULL;
			PonNode->iStatus = OCT_EXCLUDED;
			return(OCT_EXCLUDED);
		} 

		/*
		 *  Consider outer radius
		 */
		if ( d - dHalfDiagonal < dShellRadius ) { 
			/* 
			 *  Box is at least partially inside atom's outer 
			 *	radius: add to new list. Even if box
			 *	is completely included in an outer radius, 
			 *	later discovery of partial overlap with
			 *	some atom's inner radius will force deeper
			 *	evaluation, so the list must be complete.
			 */
			PaNewAtoms[iNewAtoms++] = *PaAtom;

			if ( d + dHalfDiagonal < dShellRadius )
				iIncluded++; /* completely inside outer */
			else
				iPartialOut++; /* partially inside outer */
			if ( d - dHalfDiagonal < dAtomTemp( *PaAtom ) ) {
				/* 
				 *  Partial inclusion in inner radius of
				 *	this atom.
				 */
				iPartialIn++;
			}
		}
	}

	/*
	 *  Box is not completely inside some inner radius. If it
	 *	is partially inside one or more inner radii, it will
	 *	be evaluated deeper. But first, if box is fully included 
	 *	in target zone or completely outside outer radius, finish 
	 *	it off.
	 */
	if (!iPartialIn) {
		if (iIncluded) {
			/*
		 	 *  Box is completely in 1 or more outer radii and
		 	 *	does not overlap any inner radii.
		 	 */
#ifdef OCTDEBUG
			included[oct->depth]++;
			depth[oct->depth]++;
#endif
			PonNode->iStatus = OCT_INCLUDED;
			PonNode->iAtoms = iNewAtoms;
			return(OCT_INCLUDED);
		} 
		if ( !iPartialOut ) {
			/*
			 *  Box is completely outside outer radii
			 */
#ifdef OCTDEBUG
			outside[oct->depth]++;
			depth[oct->depth]++;
#endif
#ifdef DBG2
			fprintf(stderr, "@@@ node 0x%x free atoms 0x%x\n", 
					PonNode, PonNode->PaAtomList);
#endif
			PonNode->iStatus = OCT_EXCLUDED;
			FREE( PaNewAtoms );
			PonNode->PaAtomList = NULL;
			return(OCT_EXCLUDED);
		}
	}
	PonNode->iAtoms = iNewAtoms;

	/*
	 *  Partial inclusion in target zone (part(s) of box overlap
	 *	inner radius or outer).
	 *  Create sub-boxes for recursion.
	 */
	PonChildren = PonMakeChildren( &PonNode->vCorner, 
					PonNode->iDepth + 1, OCT_UNKNOWN );
	PonNode->PonChildren = PonChildren;

	/*
	 *  Recurse on sub-boxes, keeping track of their status.
	 */
	iIncluded = 0;
	iExcluded = 0;
	for (i=0; i<8; i++) {
		switch (iBuildShellOctant( &PonChildren[i], 
					  iNewAtoms, PaNewAtoms ) ) {
		    case OCT_INCLUDED:
			iIncluded++;
			break;
		    case OCT_EXCLUDED:
			iExcluded++;
			break;
		    case OCT_PARTIAL:
			break;
		    default: 
			fprintf(stderr, "bad type\n");
			exit(1);
		}
	}

	/*
	 *  See if whole box is inside or outside zone of interest
	 *	as a sum of subbox-singleatom considerations.
	 */

	if ( iExcluded == 8 ) {

		/*
		 *  Whole box is within vdws of multiple atoms.
		 */

#ifdef OCTDEBUG
		multex[oct->depth]++;
#endif
		DestroyOctant( PonNode, OCT_EXCLUDED );
		return( OCT_EXCLUDED );
	} 
	if ( iIncluded == 8 ) {

		for (i=0; i<8; i++) {
			if (PonChildren[i].PaAtomList != NULL) {
#ifdef DBG2
				fprintf(stderr, 
				  "@@@ not null incl depth %d free 0x%x\n", 
				  PonChildren[i].iDepth, PonChildren[i].PaAtomList);
#endif
				FREE ( PonChildren[i].PaAtomList );
			}
		}
		/*
		 *  Whole box included via outer radii of multiple atoms.
		 */
#ifdef OCTDEBUG
		multwhole[oct->depth]++;
#endif
		FREE( PonChildren );
		PonNode->PonChildren = NULL;
		PonNode->iStatus = OCT_INCLUDED;
		return(OCT_INCLUDED);
	}

	/*
	 *  Remaining case: some children included, others not.
	 */
	fVolume += iIncluded * dHalfEdge * dHalfEdge * dHalfEdge;
	iTreeGridPoints += iIncluded * PiDensities[PonNode->iDepth+1]; 
	PonNode->iStatus = OCT_PARTIAL;
	return(OCT_PARTIAL);
}

/*
 *  For an interior octant, the point is to keep atom lists at
 *	the bottom nodes of the tree.
 */
static int
iBuildInteriorOctant( OCTNODEt *PonNode, int iAtoms, ATOM *PaAtomList )
{
int		i, iIncluded, iExcluded; 
ATOM		*PaAtom;
ATOM		*PaNewAtoms;
int		iNewAtoms;
double		d, dCenterRadiusSq, dHalfEdge;
VECTOR		vCenter;
OCTNODEt	*PonChildren;


#ifdef OCTDEBUG
		depth[PonNode->iDepth]++;
#endif

	/*
	 *  Evaluate whole box with respect to atoms by considering
	 *	overlaps with single atoms' inner and outer radii.
	 *	Use center of box + half length of diagonal for test.
	 */

	dHalfEdge = PdHalfEdges[PonNode->iDepth];
	vCenter = PonNode->vCorner;
	vCenter.dX += dHalfEdge;
	vCenter.dY += dHalfEdge;
	vCenter.dZ += dHalfEdge;
	dCenterRadiusSq = PdHalfDiagonals[PonNode->iDepth];
	dCenterRadiusSq += 3; /* TODO - arbitrary.. */
	/* dCenterRadiusSq *= dCenterRadiusSq; */

	/*
	 *  Get memory for list of atoms inside sphere enclosing this cube, which 
	 *	must be <= size of parent's.
	 */
	MALLOC( PaNewAtoms, ATOM *, iAtoms * sizeof(ATOM) );
	memset( PaNewAtoms, 0, iAtoms * sizeof(ATOM) );
	PonNode->PaAtomList = PaNewAtoms;
#ifdef DBG2
	fprintf(stderr, "@@@ inode 0x%x atoms 0x%x\n", PonNode, PonNode->PaAtomList);
#endif
	iNewAtoms = 0;

	/*
	 *  Consider all parent's atoms in relation to this box.
	 */
	iIncluded = 0;
	PaAtom = PaAtomList;
	for (i=0; i<iAtoms; i++, PaAtom++) {
		d = dDistance( &vCenter, &vAtomPosition( *PaAtom ) );
		if ( d < dCenterRadiusSq ) {
			/* 
			 *  It's w/in larger sphere
			 */
			PaNewAtoms[iNewAtoms++] = *PaAtom;
		}
	}
	if (!iNewAtoms) {
		/*
		 *  No atoms w/in box-enclosing sphere
		 */
#ifdef OCTDEBUG
		outside[oct->depth]++;
		depth[oct->depth]++;
#endif
		FREE( PonNode->PaAtomList );
		PonNode->PaAtomList = NULL;
		PonNode->iStatus = OCT_EXCLUDED;
		return(OCT_EXCLUDED);
	}

	PonNode->iAtoms = iNewAtoms;

	/*
	 *  End recursion if the limit of resolution (octree max depth)
	 *	has been reached or if the density is adequate for 
	 *	linear search.
	 */
	if ( PonNode->iDepth == iMaxDepth || ( iNewAtoms && iNewAtoms < 10 )) {
		iTreeGridPoints += PiDensities[PonNode->iDepth]; 
		return(OCT_INCLUDED);
	}

	/*
	 *  Create sub-boxes for recursion.
	 */
	PonChildren = PonMakeChildren( &PonNode->vCorner, 
					PonNode->iDepth + 1, OCT_UNKNOWN );
	PonNode->PonChildren = PonChildren;

	/*
	 *  recurse on sub-boxes
	 */
	iIncluded = 0;
	iExcluded = 0;
	for (i=0; i<8; i++) {
		switch (iBuildInteriorOctant( &PonChildren[i], 
					  iNewAtoms, PaNewAtoms ) ) {
		  case OCT_INCLUDED:
			iIncluded++;	/* for debug tracking */
			break;
		  case OCT_EXCLUDED:
			iExcluded++;
			break;
		  default: 
			fprintf(stderr, "bad type\n");
			exit(1);
		}
	}

	/*
	 *  
	 */
	if ( iExcluded == 8 ) {
		FREE( PonNode->PonChildren );
		PonNode->PonChildren = NULL;
	} else {
		FREE( PonNode->PaAtomList );
#ifdef DBG2
		fprintf(stderr, "@@@ node 0x%x free atoms 0x%x\n", PonNode, PonNode->PaAtomList);
#endif
		PonNode->PaAtomList = NULL;
	}

	PonNode->iStatus = OCT_INCLUDED;
	return(OCT_INCLUDED);
}
OCTREE
octOctTreeCreate( UNIT uUnit, int iType, double dGridSpace, double dAddExtent, double dShellExtent,
			int iIncludeSolvent )
{
OCTREE		octTree;
VECTOR		vMinCorner, vMaxCorner, vAtom;
VARARRAY	vaAtoms;
ATOM		aAtom, *PaAtoms;
LOOP		lAtoms, lRes;
RESIDUE		rRes;
int		i, j, iAtoms, iDefaultedRadius;
double		dMaxRadius, dCharge;
double		dTx, dTy, dTz, dTmax, dTmp;

	if ( !uUnit )
		return(NULL);

	/*
	 *  Set globals.
	 */
	dGridSize = dGridSpace;
	dShellRadius = dShellExtent;	/* shell: clearance beyond inner shell*/

	/*
	 *  Create the octree "object" and initialize
	 */
	MALLOC( octTree, OCTREE, sizeof(OCTREEt) );
	octTree->iType = iType;
	octTree->dGridSize = dGridSpace;
	octTree->PfCharges = NULL;
	octTree->PdHalfEdges = NULL;
	octTree->PdHalfDiagonals = NULL;

	/*
	 *  Make array for atom pointers
	 */
	vaAtoms = vaVarArrayCreate( sizeof(ATOM) );

	/*
	 *  Fill array of atom pointers according to octree type, initializing
	 *	the atoms' temporary floating-point values w/ their sizes.
	 */

	iDefaultedRadius = 0;
	lRes = lLoop( (OBJEKT)uUnit, RESIDUES );
	switch ( iType ) {
	    case OCT_SHELL:
	    case OCT_INTERIOR_SOLUTE:
		while ((rRes = (RESIDUE) oNext(&lRes))) {
			if ( iIncludeSolvent  ||  
			     cResidueType( rRes ) != RESTYPESOLVENT ) {
    	    			lAtoms = lLoop( (OBJEKT)rRes, ATOMS );
    	    			while ((aAtom = (ATOM) oNext(&lAtoms))) {
					VarArrayAdd( vaAtoms, (GENP)&aAtom );
					iDefaultedRadius +=
						iAtomSetTmpRadius( aAtom );
				}
	    		}
		}
		break;

	    case OCT_INTERIOR_SOLVENT: {
		int nowarning = 1;
		while ((rRes = (RESIDUE) oNext(&lRes))) {
			if ( cResidueType( rRes ) == RESTYPESOLVENT ) {
				i = 0;
    	    			lAtoms = lLoop( (OBJEKT)rRes, ATOMS );
    	    			while ((aAtom = (ATOM) oNext(&lAtoms))) {
					if (!i)
						VarArrayAdd( vaAtoms, (GENP)&aAtom );
					iDefaultedRadius +=
						iAtomSetTmpRadius( aAtom );
					i++;
				}
				if ( i > 3  &&  nowarning ) {
					VP0(( "Warning: non-water solvent %s\n",
						"may lead to steric problems"));
					nowarning = 0;
				}
	    		}
		}
		}
		break;
	    default:
		VP0(("Octree type not implemented\n"));
		exit(1);
	}

	/*
	 *  If no atoms, nothing to build.
	 */
	iAtoms = iVarArrayElementCount( vaAtoms );
/*
VP0(("atoms: %d\n", iAtoms));
dTmax = 20.0;
for (dTmp=dTmax, iMaxDepth=0; dTmp>dGridSize; iMaxDepth++, dTmp/=2.0);
octTree->iMaxDepth = iMaxDepth;
for (dTmax=dGridSize, i=0; i<iMaxDepth; i++, dTmax*=2.0);
MALLOC( PdHalfEdges, double *, (iMaxDepth+1) * sizeof(double) );
octTree->PdHalfEdges = PdHalfEdges;
MALLOC( PdHalfDiagonals, double *, (iMaxDepth+1) * sizeof(double) );
octTree->PdHalfDiagonals = PdHalfDiagonals;
return(octTree);
*/

	if ( !iAtoms ) {
		VarArrayDestroy( &vaAtoms );
		FREE( octTree );
		return(NULL);
	}
	octTree->vaAtoms = vaAtoms;

	if ( iDefaultedRadius )
		VP0(( "Used default radius %5.2f for %d atoms\n", 
				ATOM_DEFAULT_RADIUS, iDefaultedRadius ));

	/*
	 *  Find bounding box of atom centers and max radius, while
	 *	adding the additional extents to the temporary atom radii.
	 *	For purposes of finding max radius, initialize 1st atom
	 *	then loop over others, comparing.
	 */
	PaAtoms = PVAI( vaAtoms, ATOM, 0 );
	dCharge = dAtomCharge(  *PaAtoms );
	dMaxRadius = dAtomTemp( *PaAtoms );
	vMinCorner = vMaxCorner = vAtomPosition( *PaAtoms );
	AtomTempDoubleIncrement( *PaAtoms, dAddExtent );
	for (i=1,PaAtoms++; i<iAtoms; i++, PaAtoms++) {
		dCharge += dAtomCharge( *PaAtoms );
		dMaxRadius = MAX( dMaxRadius, dAtomTemp( *PaAtoms ) );
		AtomTempDoubleIncrement( *PaAtoms, dAddExtent );
		vAtom = vAtomPosition( *PaAtoms );
		if (vAtom.dX < vMinCorner.dX)
			vMinCorner.dX = vAtom.dX;
		else if (vAtom.dX > vMaxCorner.dX)
			vMaxCorner.dX = vAtom.dX;
		if (vAtom.dY < vMinCorner.dY)
			vMinCorner.dY = vAtom.dY;
		else if (vAtom.dY > vMaxCorner.dY)
			vMaxCorner.dY = vAtom.dY;
		if (vAtom.dZ < vMinCorner.dZ)
			vMinCorner.dZ = vAtom.dZ;
		else if (vAtom.dZ > vMaxCorner.dZ)
			vMaxCorner.dZ = vAtom.dZ;
	}

	/*
	 *  Expand bounding box if neccessary.
	 */
	switch ( iType ) {
	    case OCT_SHELL:
		/*
		 *  Offset shell by max solute atom radius
		 */
		dShellRadius += dMaxRadius + dAddExtent;

		/*
		 *  Enclose atom centers with shell
		 */
		vMinCorner.dX -= dShellRadius;
		vMinCorner.dY -= dShellRadius;
		vMinCorner.dZ -= dShellRadius;
		vMaxCorner.dX += dShellRadius;
		vMaxCorner.dY += dShellRadius;
		vMaxCorner.dZ += dShellRadius;
		break;

	    case OCT_INTERIOR_SOLUTE:
	    case OCT_INTERIOR_SOLVENT:
		/*
		 *  Enclose atom centers with radius of largest atom
		 */
		vMinCorner.dX -= dMaxRadius;
		vMinCorner.dY -= dMaxRadius;
		vMinCorner.dZ -= dMaxRadius;
		vMaxCorner.dX += dMaxRadius;
		vMaxCorner.dY += dMaxRadius;
		vMaxCorner.dZ += dMaxRadius;
		break;
	    default:
		VP0(("Octtree type not implemented\n"));
		exit(1);
	}
	VP1(( "Total solute charge:  %5.2f  Max atom radius:  %5.2f\n", 
					dCharge, dMaxRadius ));
	if ( iType == OCT_SHELL )
		VP0(( "Grid extends from solute vdw + %.2f  to  %.2f\n",
				dAddExtent, dShellRadius ));
	VP1(( "Box:\n" ));
	VP1(( "   enclosing:  %5.2f %5.2f %5.2f   %5.2f %5.2f %5.2f\n",
				vMinCorner.dX, vMinCorner.dY, vMinCorner.dZ,
				vMaxCorner.dX, vMaxCorner.dY, vMaxCorner.dZ));

	/*  
	 *  Find longest edge of box
	 */
	dTx = vMaxCorner.dX - vMinCorner.dX;
	dTy = vMaxCorner.dY - vMinCorner.dY;
	dTz = vMaxCorner.dZ - vMinCorner.dZ;

	dTmax = MAX( dTx, dTy );
	dTmax = MAX( dTmax, dTz );

	/*  
	 *  Determine octree descent depth required to obtain 
	 *	requested grid resolution by dividing longest edge in
	 *	half until it is less than dGridSize.
	 */
	for (dTmp=dTmax, iMaxDepth=0; dTmp>dGridSize; iMaxDepth++, dTmp/=2.0);
	octTree->iMaxDepth = iMaxDepth;

	/*  
	 *  Scale box up to modulo grid size 
	 */
	for (dTmax=dGridSize, i=0; i<iMaxDepth; i++, dTmax*=2.0);

	/*
	 *  Just for fun, calc the new max corner
	 */
	vMaxCorner.dX += dTmax - dTx;
	vMaxCorner.dY += dTmax - dTy;
	vMaxCorner.dZ += dTmax - dTz;

	/*  
	 *  Set up arrays (indexed by depth) of edge size & half-diagonals 
	 *	to save recalculation
	 */
	MALLOC( PdHalfEdges, double *, (iMaxDepth+1) * sizeof(double) );
	octTree->PdHalfEdges = PdHalfEdges;
	MALLOC( PdHalfDiagonals, double *, (iMaxDepth+1) * sizeof(double) );
	octTree->PdHalfDiagonals = PdHalfDiagonals;

	for (dTmp=dTmax/2.0, i=0; i<=iMaxDepth; dTmp/=2.0, i++) {
		PdHalfEdges[i] = dTmp;
		/* diagonal = sqrt( 3 * side^2 ) */
		PdHalfDiagonals[i] = sqrt( 3.0 * (dTmp * dTmp) );
	}

	/*  
	 *  Set up array (indexed by depth) of points per box, used
	 *	for rejuggling charge array when splitting boxes.
	 */
	MALLOC( PiDensities, int *, (iMaxDepth+1) * sizeof(int) );
	octTree->PiDensities = PiDensities;
	for (i=iMaxDepth, j=1; i>-1; i--, j*=8)
		PiDensities[i] = j;

	VP1(( "   sized:\t\t\t      %5.2f %5.2f %5.2f\n",
			vMaxCorner.dX, vMaxCorner.dY, vMaxCorner.dZ));
	VP1(( "   edge:        %5.2f\n", dTmax ));
	VP0(( "Resolution:     %5.2f Angstrom.\n", dGridSize ));
	VP1(( "Tree depth: %d\n", iMaxDepth));

	/*
	 *  Build head node w/ all atoms in list.
	 */
	octTree->onHead.iStatus = OCT_UNKNOWN;
	octTree->onHead.iDepth = 0;
	octTree->onHead.vCorner = vMinCorner;
	/* 
	 *  The master atom list is part of a vararray under the OCTREE
	 *	so, unlike the seperately allocated atom lists in the
	 *	lower nodes on the tree, the master list is not freed:
	 *	rather, its vararray is destroyed. Reset the ATOM* ptr.
	 */
	octTree->onHead.PaAtomList = NULL;	
	PaAtoms = PVAI( vaAtoms, ATOM, 0 );

	/*
	 *  Build octree recursively.
	 */
	iTreeGridPoints = 0;
	fVolume = 0.0;
	iNodeCount = 0;
        time_start = time((time_t *) 0);

	switch ( iType ) {
	  case OCT_SHELL:
		iBuildShellOctant( &octTree->onHead, iAtoms, PaAtoms );
		break;
	  case OCT_INTERIOR_SOLUTE:
	  case OCT_INTERIOR_SOLVENT:
		iBuildInteriorOctant( &octTree->onHead, iAtoms, PaAtoms );
		break;
	  default:
		VP0(("bad switch\n"));
		exit(1);
	}
	octTree->iTreePoints = iTreeGridPoints;	/* global */

	VP0(( "grid build: %ld sec\n", 
			time((time_t *)0) - time_start ));
	VP1(( "Volume = %5.2f%% of box, grid points %d\n", 
			100 * fVolume / ( dTmax * dTmax * dTmax ),
			iTreeGridPoints ));
#ifdef OCTDEBUG
	fprintf(stderr, "depth  r_inc   r_inex  r_out  multin  multout \n");
	for (i=0;i<=iMaxDepth; i++) {
		fprintf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\n", 
			depth[i], 
			included[i],
			insex[i], 
			outside[i], 
			multwhole[i],
			multex[i]
			);
	}
#endif
	/*
	 *  Reset atom radii.
	 */
	PaAtoms = PVAI( vaAtoms, ATOM, 0 );
	for (i=0; i<iAtoms; i++, PaAtoms++)
		AtomTempDoubleIncrement( *PaAtoms, -dAddExtent );
	return(octTree);
}

/*************************************************************************

BuildShellOctant: Build octree for shell around a list of atoms. 

	"Included" means in the shell.
	"Excluded" means either inside an atom or beyond the shell.

	Lists of atoms enclosing or intersecting a box with outer radius
	are built and passed to the children to reduce their evaluation 
	time. For the shell octree, these lists are not saved.

 *************************************************************************/





/*************************************************************************
 *************************************************************************/

void
OctTreeDestroy( OCTREE *PoctTree )
{
	/*
	 *  Free top-level things.
	 */
	if ( (*PoctTree)->vaAtoms )
		VarArrayDestroy( &(*PoctTree)->vaAtoms );

	if ( (*PoctTree)->PfCharges != NULL )
		FREE( (*PoctTree)->PfCharges );

	FREE( (*PoctTree)->PdHalfEdges );
	FREE( (*PoctTree)->PdHalfDiagonals );
	FREE( (*PoctTree)->PiDensities );

	/*
	 *  Set global and recursively delete the tree.
	 */
	iMaxDepth = (*PoctTree)->iMaxDepth;
	DestroyOctant( &(*PoctTree)->onHead, OCT_UNKNOWN );
	FREE( *PoctTree );
	*PoctTree = NULL;
	return;
}


/* @@@ */
/*************************************************************************
 *************************************************************************/

static void
OctNodeInitCharges( OCTNODEt *PonNode )
{
int	i, j, k, l;
int 	ct = iMaxDepth - PonNode->iDepth + 1;
int	iCompCharge;
VECTOR	vPoint;
ATOM	*PaAtom;
double	d;

	if ( PonNode->iStatus == OCT_PARTIAL ) {
		for (i=0; i<8; i++)
			OctNodeInitCharges( &PonNode->PonChildren[i] );
		return;
	}
	if ( PonNode->iStatus == OCT_EXCLUDED )
		return;

	/*
	 *  Node is included: calculate charges.
	 * 	Set node pointer into master charge array 
	 */
	PonNode->PfCharges = PfCharges;

	/* 
	 *  For each point x,y,z in this quadrant, 
	 *	loop over atoms, accumulating charge 
	 */

	vPoint.dX = PonNode->vCorner.dX;
	for (i=0; i<ct; i++, vPoint.dX+=dGridSize) {
	    vPoint.dY = PonNode->vCorner.dY;
	    for (j=0; j<ct; j++, vPoint.dY+=dGridSize) {
		vPoint.dZ = PonNode->vCorner.dZ;
		for (k=0; k<ct; k++, vPoint.dZ+=dGridSize) {
			/*
			 *  Got point: loop over atoms, accumulating charge.
			 */
			*PfCharges = 0.0;
			PaAtom = PaChargeAtoms;
			iCompCharge = 1;
			for (l=0; l<iChargeAtoms; l++, PaAtom++) {
				double 	dX, dY, dZ;

				dX = vPoint.dX - vAtomPosition(*PaAtom).dX;
				dX = dX * dX;
				dY = vPoint.dY - vAtomPosition(*PaAtom).dY;
				dY = dY * dY;
				dZ = vPoint.dZ - vAtomPosition(*PaAtom).dZ;
				dZ = dZ * dZ;
				d = dX + dY + dZ;
				if ( iDistanceCharge )
					d = sqrt(d);
				*PfCharges += dAtomCharge( *PaAtom ) / d;
				if ( d < dAtomTemp( *PaAtom ) )
					iCompCharge = 0;
			}
			if ( iCompCharge ) {
				/*
				 *  Keep track of max, min charges and 
				 *	their locations.
				 */
				if ( *PfCharges > fMaxCharge ) {
					fMaxCharge = *PfCharges;
					vMaxCharge = vPoint;
				} else if ( *PfCharges < fMinCharge ) {
					fMinCharge = *PfCharges;
					vMinCharge = vPoint;
				}
			}
			PfCharges++;
		}
	    }
	}
	return;
}
void
OctTreeInitCharges( OCTREE octTree, int iAtomOption, int iDielectric, double dCutDist,
							VECTOR *PvMin, VECTOR *PvMax )
{
int	i;
ATOM	*PaAtom;

	if ( octTree->iType != OCT_SHELL ) {
		VP0(( "InitCharges: wrong tree type\n" ));
		exit(1);
	}
	if ( !octTree->iTreePoints ) {
		VP0(( "InitCharges: no grid (?)\n" ));
		exit(1);
	}
	VP0(( "Calculating grid charges\n" ));
        time_start = time((time_t *) 0);

	/*
	 *  Get the memory for the charges and note the dielectric.
	 */

	MALLOC( PfCharges, float *, octTree->iTreePoints * sizeof(float) );
	octTree->PfCharges = PfCharges;
	octTree->iDielectric = iDielectric;

	/*
	 *  Set up globals for node descent.
	 */
	if ( iAtomOption == AT_OCTREE ) {
		iChargeAtoms = iVarArrayElementCount( octTree->vaAtoms );
		PaChargeAtoms = PVAI( octTree->vaAtoms, ATOM, 0 );
	} else
		exit(1);
	fMaxCharge = -FLT_MAX;
	fMinCharge = FLT_MAX;
	iMaxDepth = octTree->iMaxDepth;
	dGridSize = octTree->dGridSize;

	if ( iDielectric == DIEL_R2 ) {
		iDistanceCharge = 0;
	} else
		iDistanceCharge = 1;

	/*
	 *  Add new radius to atoms.
	 */
	if ( iDielectric == DIEL_R2 ) { 
		PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
		for (i=0; i<iChargeAtoms; i++, PaAtom++) {
/* TODO - chargeatoms not necc whole set */
			AtomTempDoubleIncrement( *PaAtom, dCutDist );
			AtomTempDoubleSquare( *PaAtom );
		}
	} else {
		PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
		for (i=0; i<iChargeAtoms; i++, PaAtom++)
			AtomTempDoubleIncrement( *PaAtom, dCutDist );
	}
	/*
	 *  Descend the octree.
	 */
	OctNodeInitCharges( &octTree->onHead );
	*PvMin = vMinCharge;
	*PvMax = vMaxCharge;

	/*
	 *  Restore atoms' radii.
	 */
	if ( iDielectric == DIEL_R2 ) { 
		PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
		for (i=0; i<iChargeAtoms; i++, PaAtom++) {
			AtomTempDoubleSquareRoot( *PaAtom );
			AtomTempDoubleIncrement( *PaAtom, -dCutDist );
		}
	} else {
		PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
		for (i=0; i<iChargeAtoms; i++, PaAtom++)
			AtomTempDoubleIncrement( *PaAtom, -dCutDist );
	}
	VP0(( "charges: %ld sec\n", time((time_t *) 0) - time_start));
	return;
}


/*************************************************************************
 *************************************************************************/

static void
OctNodePrintGrid( OCTNODEt *PonNode )
{
VECTOR	vPoint;
int	i, j, k, ct;
float	*PfCharge;
	iNodeCount++;
	if ( PonNode->iStatus == OCT_PARTIAL ) {
	  	for (i=0; i<8; i++) 
			OctNodePrintGrid( &PonNode->PonChildren[i] );
		return;
	}
	if ( PonNode->iStatus == OCT_EXCLUDED )
		return;

	/*
	 *  Included - print the box.
	 */

	ct = iMaxDepth - PonNode->iDepth + 1;
	PfCharge = PonNode->PfCharges;
	vPoint.dX = PonNode->vCorner.dX;
	for (i=0; i<ct; i++, vPoint.dX+=dGridSize) {
	    vPoint.dY = PonNode->vCorner.dY;
	    for (j=0; j<ct; j++, vPoint.dY+=dGridSize) {
	    	vPoint.dZ = PonNode->vCorner.dZ;
		for (k=0; k<ct; k++, vPoint.dZ+=dGridSize) {
			switch ( iColorMethod ) {
			  case COLOR_RANGE:
				fprintf(fChargeFile, ".color %d\n",
			  		(int)floor( 5 + 60 *
					(*PfCharge-fMinCharge) / 
					(fMaxCharge-fMinCharge) ) );
				fprintf(fChargeFile, ".dot %f %f %f\n",
					vPoint.dX, vPoint.dY, vPoint.dZ);
				break;
			  case COLOR_CUT:
				if ( *PfCharge < -0.1 )
					fprintf(fChargeFile, ".color yellow\n");
				else if ( *PfCharge > 0.1 )
					fprintf(fChargeFile, ".color cyan\n");
				else
					fprintf(fChargeFile, ".color black\n");
				fprintf(fChargeFile, ".dot %f %f %f\n",
					vPoint.dX, vPoint.dY, vPoint.dZ);
				break;
			  case COLOR_DEPTH:
				if ( PonNode->iDepth == iMaxDepth )
					fprintf(fChargeFile, ".color white\n");
				else if ( PonNode->iDepth % 2 )
					fprintf(fChargeFile, ".color red\n");
				else
					fprintf(fChargeFile, ".color cyan\n");
				fprintf(fChargeFile, ".dot %f %f %f\n",
					vPoint.dX, vPoint.dY, vPoint.dZ);
				break;
			  case COLOR_NONE:
				fprintf(fChargeFile, "%f   %f %f %f\n",
					*PfCharge, 
					vPoint.dX, vPoint.dY, vPoint.dZ);
				break;
			  default:
				fprintf(fChargeFile, ".color white\n");
			}
			PfCharge++;
		}
	    }
	}
}

void
OctTreePrintGrid( OCTREE octTree, char *sFileName, int iColor )
{
	if ( iColor != COLOR_DEPTH  &&  !octTree->PfCharges ) {
		VP0(( "charge coloring but no charges\n" ));
		return;
	}
	fChargeFile = fopen(sFileName, "w");
	if ( fChargeFile == NULL ) {
		perror(sFileName);
		exit(1);
	}
	
	/*
	 *  Set up globals for this octree.
	 */

	iColorMethod = iColor;
	iMaxDepth = octTree->iMaxDepth;
	dGridSize = octTree->dGridSize;
	iNodeCount = 0;

	/*
	 *  Print the nodes recursively.
	 */

	OctNodePrintGrid( &octTree->onHead );

	VP1(( "total nodes in octree %d\n", iNodeCount ));
	fclose( fChargeFile );
	return;
}
	

/*************************************************************************
 *************************************************************************/

/* global for speed */

static int	boxct[8];

static void
SplitIncludedNode( OCTNODEt *PonNode )
{
int		i, j, k, nchild, ct, ct2;
float		*PfTmpCharges, *PfCharge;
OCTNODEt	*PonChildren;

	if ( PonNode->PonChildren != NULL )
		DFATAL(( "Programming error\n" ));

	/*
	 *  Subdivide this node: set up children array
	 */		
	PonChildren = PonMakeChildren( &PonNode->vCorner, 
					PonNode->iDepth + 1, OCT_INCLUDED );
	PonNode->PonChildren = PonChildren;

	/*
	 *  copy the parent's charges to a temp array
	 */
	MALLOC( PfTmpCharges, float *, 
			sizeof(float) * PiDensities[PonNode->iDepth]);
	memcpy( PfTmpCharges, PonNode->PfCharges, 
			sizeof(float) * PiDensities[PonNode->iDepth]);

	/*
	 *  divide the parent's actual charge array space 
	 *	between the children (just assigning pointers)
	 */
	for (i=0; i<8; i++) {
		PonChildren[i].PfCharges =
				PonNode->PfCharges + 
				i * PiDensities[PonNode->iDepth+1];
	}

	/*
	 *  Distribute charges to the children:
	 *	Reshuffle from the order in the current node
	 *	to the orders in the subboxes. I.e. while scanning
	 *	the copy of the parent's charge array in its order,
	 *	distribute the charges to the proper places in the 
	 *	childrens' arrays. 
	 *
	 *  For this to be relatively simple, the order of the children
	 *	in space (in the octree definition) mimics the order of 
	 *	evaluation of the parent's charges - see PonMakeChildren():
	 *
	 *		inner loop on Z axis (flagged by low order bit = '1')
	 *		next loop on Y axis (second bit = '2') 
	 *		outer loop on X axis (3rd bit = '4').
	 *
	 *	Use of the bits as described orders the subboxes properly.
	 */


	/*
	 *  initialize per-child point counts for indexing
	 *	childrens' point arrays
	 */
	for (i=0; i<8; i++)
		boxct[i] = 0;

	/*
	 *  nchild indicates child that current point of parent
	 *	will go to; 1st child has all bits off (000)
	 */
	nchild = 0;

	/*
	 *  total points in parent node
	 */
	ct = iMaxDepth - PonNode->iDepth + 1;
	ct2 = ct / 2;

	/*
	 *  loop over points in parent node
	 */
	PfCharge = PfTmpCharges;
	for (i=0; i<ct; i++) {
		/* loop over X axis */
		if (i==ct2)		/* 2nd half of X axis */
			nchild |= 4;	/* turn on X bit for 2nd slice*/
		nchild &= ~2;	/* turn off Y bit */
		for (j=0; j<ct; j++) {
			/* loop over Y axis */
			if (j==ct2)		/* 2nd half of Y axis */
				nchild |= 2;	/* Y bit on for 2nd slice */
			nchild &= ~1;	/* turn off Z bit */
			for (k=0; k<ct; k++) {
				/* loop over Z axis */
				if (k==ct2)		/* 2nd half of Z axis */
					nchild |= 1;	/* Z bit on for 2nd slice */
		      
				PonChildren[nchild].PfCharges[boxct[nchild]++] = 
							*PfCharge;
				PfCharge++;
			}
		}
	}
	FREE( PfTmpCharges );
}
static int
OctNodeDeleteSphere( OCTNODEt *PonNode )
{
int		i, iIncluded, iExcluded;
double		d, dHalfEdge, dHalfDiagonal;
OCTNODEt	*PonChildren;
VECTOR		vCenter;

	/*
	 *  If already excluded, quit. 
	 */
	if ( PonNode->iStatus == OCT_EXCLUDED )
		return(OCT_EXCLUDED);

	/*
	 *  If resolution reached, check corner only.
	 */
	if ( PonNode->iDepth == iMaxDepth ) {
		d = dDistance( &PonNode->vCorner, &vNewPoint );
		if ( d < dDeleteRadius ) {
#ifdef OCTDEBUG
			fprintf(subtractions, ".dot %f %f %f\n", 
				PonNode->vCorner.dX, 
				PonNode->vCorner.dY, 
				PonNode->vCorner.dZ);
#endif
			PonNode->iStatus = OCT_EXCLUDED;
		}
		return(PonNode->iStatus);
	}

	/*
	 *  consider whole box from center
	 */
	dHalfEdge = PdHalfEdges[PonNode->iDepth];
	vCenter = PonNode->vCorner;
	vCenter.dX += dHalfEdge;
	vCenter.dY += dHalfEdge;
	vCenter.dZ += dHalfEdge;
	d = dDistance( &vCenter, &vNewPoint );
	dHalfDiagonal = PdHalfDiagonals[PonNode->iDepth];

	/*
	 *  if box completely outside sphere, leave it alone.
	 */
	if ( d - dHalfDiagonal > dDeleteRadius )
		return(OCT_INCLUDED);

	/*
	 *  if box completely inside sphere, clean up its contents.
	 */
	if ( d + dHalfDiagonal < dDeleteRadius ) {
#ifdef SUB_DBG
		XYZ	point;
		int     j, k, ct = iMaxDepth - oct->depth + 1;

		point.X = corner->X;
		for (i=0; i<ct; i++) {
		  point.Y = corner->Y;
		  for (j=0; j<ct; j++) {
		    point.Z = corner->Z;
		    for (k=0; k<ct; k++) {
			fprintf(subtractions, ".dot %f %f %f\n",
				point.X, point.Y, point.Z);
			point.Z += dGridSize;
		    }
		    point.Y += dGridSize;
		  }
		  point.X += dGridSize;
		}
#endif
		DestroyOctant( PonNode, OCT_EXCLUDED );
		return(OCT_EXCLUDED);
	}

	/*
	 *  recurse on sub-boxes
	 */

	if ( PonNode->iStatus == OCT_INCLUDED )
		SplitIncludedNode( PonNode );

	PonChildren = PonNode->PonChildren;

	/*
	 *  Evaluate children.
	 */

	iIncluded = 0;
	iExcluded = 0;
	for (i=0; i<8; i++) {
		switch ( OctNodeDeleteSphere( &PonChildren[i] ) ) {
		    case OCT_INCLUDED:
			iIncluded++;	/* for debug tracking */
			break;
		    case OCT_EXCLUDED:
			iExcluded++;
			break;
		    case OCT_PARTIAL:
			break;
		    default: 
			fprintf(stderr, "bad type\n");
			exit(1);
		}
	}
	if ( iExcluded == 8 ) {

		/*
		 *  Whole box is within vdw of multiple atoms.
		 */
		DestroyOctant( PonNode, OCT_EXCLUDED );
		return(OCT_EXCLUDED);
	} 
	/*
	 * (don't fold into 1 node if iIncluded == 8 to avoid complications 
	 *	in regrouping charges from subnodes - so even if technically 
	 *	OCT_INCLUDED, it is evaluated as an OCT_PARTIAL)
	 */
	PonNode->iStatus = OCT_PARTIAL;
	return(OCT_PARTIAL);
}
void
OctTreeDeleteSphere( OCTREE octTree, VECTOR *PvPoint, double dRadius )
{

	/*
	 *  Set up globals for this octree.
	 */

	iMaxDepth = octTree->iMaxDepth;
	PiDensities = octTree->PiDensities;
	vNewPoint = *PvPoint;
	dDeleteRadius = dRadius;
	dGridSize = octTree->dGridSize;

	/*
	 *  Delete recursively.
	 */

	OctNodeDeleteSphere( &octTree->onHead );
	return;
}



/*************************************************************************
 *************************************************************************/

static void
OctNodeUpdateCharge( OCTNODEt *PonNode, int iParentAtoms, ATOM *PaParentAtoms )
{
VECTOR	vPoint;
double	d;
int	i, j, k, l, ct;
int	iCompOk;
float	*PfCharge;
ATOM	*PaAtom, *PaAtoms;
int	iAtoms;

	/*
	 *  If partial, recurse.
	 */
	if ( PonNode->iStatus == OCT_PARTIAL ) {
		for (i=0; i<8; i++)
			OctNodeUpdateCharge( &PonNode->PonChildren[i],
				PonNode->iAtoms, PonNode->PaAtomList );
		return;
	}

	/*
	 *  If already excluded, quit.
	 */
	if ( PonNode->iStatus == OCT_EXCLUDED )
		return;

	/*
	 *  Included, do whole thing
	 */

	if ( PonNode->PaAtomList != NULL ) {
		PaAtoms = PonNode->PaAtomList;
		iAtoms = PonNode->iAtoms;
	} else {
		PaAtoms = PaParentAtoms;
		iAtoms = iParentAtoms;
	}
	PfCharge = PonNode->PfCharges;
	ct = iMaxDepth - PonNode->iDepth + 1;

	vPoint.dX = PonNode->vCorner.dX;
	for (i=0; i<ct; i++, vPoint.dX+=dGridSize) {
	    vPoint.dY = PonNode->vCorner.dY;
	    for (j=0; j<ct; j++, vPoint.dY+=dGridSize) {
		vPoint.dZ = PonNode->vCorner.dZ;
		for (k=0; k<ct; k++, vPoint.dZ+=dGridSize) {
			double 	dX, dY, dZ;

			dX = vNewPoint.dX - vPoint.dX;
			dX = dX * dX;
			dY = vNewPoint.dY - vPoint.dY;
			dY = dY * dY;
			dZ = vNewPoint.dZ - vPoint.dZ;
			dZ = dZ * dZ;
			d = dX + dY + dZ;
			if ( iDistanceCharge )
				d = sqrt(d);
			*PfCharge += fNewCharge / d;
			/*
			 *  Look at neighboring atoms to make sure
			 *	this point is valid for min/max check.
			 */
			iCompOk = 1;
			PaAtom = PaAtoms;
			for (l=0; l<iAtoms; l++, PaAtom++) {
				dX = vAtomPosition( *PaAtom ).dX - vPoint.dX;
				dX = dX * dX;
				dY = vAtomPosition( *PaAtom ).dY - vPoint.dY;
				dY = dY * dY;
				dZ = vAtomPosition( *PaAtom ).dZ - vPoint.dZ;
				dZ = dZ * dZ;
				d = dX + dY + dZ;
				if ( d < dAtomTemp( *PaAtom ) ) {
					iCompOk = 0;
					break;
				}
			}
			if ( iCompOk ) {
				if (*PfCharge > fMaxCharge ) {
					fMaxCharge = *PfCharge;
					vMaxCharge = vPoint;
				} else if (*PfCharge < fMinCharge ) {
					fMinCharge = *PfCharge;
					vMinCharge = vPoint;
				}
			} else {
				/* HACK to ensure printgrid coloring ok */
				*PfCharge = 0.0;
			}
			PfCharge++;
		}
	    }
	}
}
void
OctTreeUpdateCharge( OCTREE octTree, VECTOR *PvNewPoint, float fCharge, double dCutDist,
							VECTOR *PvMax, VECTOR *PvMin )
{
ATOM	*PaAtom;
int	i;
	if ( octTree->iType != OCT_SHELL ) {
		VP0(( "UpdateCharge: wrong tree type\n" ));
		exit(1);
	}
	if ( !octTree->PfCharges ) {
		VP0(( "UpdateCharge: charges not initted\n" ));
		exit(1);
	}

	/*
	 *  Set up globals for octree.
	 */
	vNewPoint = *PvNewPoint;
	fNewCharge = fCharge;
	fMaxCharge = -FLT_MAX;
	fMinCharge = FLT_MAX;
	iMaxDepth = octTree->iMaxDepth;
	dGridSize = octTree->dGridSize;
	if ( octTree->iDielectric == DIEL_R2 )
		iDistanceCharge = 0;
	else
		iDistanceCharge = 1;

	/*
	 *  Add new radius to atoms.
	 */
	PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
	for (i=0; i<iChargeAtoms; i++, PaAtom++) {
		AtomTempDoubleIncrement( *PaAtom, dCutDist );
		AtomTempDoubleSquare( *PaAtom );
	}

	/*
	 *  Descend octree.
	 */
	OctNodeUpdateCharge( &octTree->onHead, iChargeAtoms, 
				PVAI( octTree->vaAtoms, ATOM, 0 ) );
	*PvMin = vMinCharge;
	*PvMax = vMaxCharge;
	/*
	 *  Restore atoms' radii.
	 */
	PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 );
	for (i=0; i<iChargeAtoms; i++, PaAtom++) {
		AtomTempDoubleSquareRoot( *PaAtom );
		AtomTempDoubleIncrement( *PaAtom, -dCutDist );
	}
	return;
}


/*************************************************************************
 *************************************************************************/
/*
 *  Solvent check: see if a solvent overlaps a new ion. 
 *	Descend tree nodes which contain the new point inside their
 *	box-including sphere till a final node is reached, then
 *	check the atom list for that node.
 *
 *  NOTE - this code assumes that atom lists are updated as solvent
 *	is deleted. In fact, this is not done, so this code is left
 *	for adaptation in finding vdw pairs.
 */

static int
OctNodeCheckSolvent( OCTNODEt *PonNode )
{
int	i;
ATOM	*PaAtom;
VECTOR	vCenter;
double	d, dHalfEdge, dHalfDiagonal;

	if ( PonNode->iStatus != OCT_INCLUDED )
		return(0);

	/*
	 *  See if newpoint is in minimally box-enclosing sphere
	 */

	dHalfEdge = PdHalfEdges[PonNode->iDepth];
	vCenter = PonNode->vCorner;
	vCenter.dX += dHalfEdge;
	vCenter.dY += dHalfEdge;
	vCenter.dZ += dHalfEdge;
	d = dDistance( &vCenter, &vNewPoint );
	dHalfDiagonal = PdHalfDiagonals[PonNode->iDepth];

	if ( d > dHalfDiagonal ) {
		return(0);
	}
/*
VP0(("lev %d within %f %f %f  %f %f %f\n", 
PonNode->iDepth,
PonNode->vCorner.dX,
PonNode->vCorner.dY,
PonNode->vCorner.dZ,
PonNode->vCorner.dX + 2 * PdHalfEdges[PonNode->iDepth],
PonNode->vCorner.dY + 2 * PdHalfEdges[PonNode->iDepth],
PonNode->vCorner.dZ + 2 * PdHalfEdges[PonNode->iDepth]));
*/
	/*
	 *  If at last depth w/ atomlist, check them.
	 */
	if ( PonNode->PonChildren == NULL ) {
/*
VP0(("final? \n"));
*/
		PaAtom = PonNode->PaAtomList;
		for (i=0; i<PonNode->iAtoms; i++, PaAtom++) {
			d = dDistanceSq( &vNewPoint,
					&vAtomPosition( *PaAtom ) );
/*
VP0(("d= %f \n", d));
*/
			if ( d < dClosestDistance ) {
				aClosestAtom = *PaAtom;
				dClosestDistance = d;
			}
		}
		return(1);
	}

	/*
	 *  partial: subdivide
	 */
	for (i=0; i<8; i++) 
		if( OctNodeCheckSolvent( &PonNode->PonChildren[i] ) )
			return(1);

	return(0);
}



RESIDUE
rOctTreeCheckSolvent( OCTREE octTree, VECTOR *PvPoint )
{
ATOM	*PaAtom;
	if ( octTree->iType != OCT_INTERIOR_SOLVENT ) {
		VP0(( "CheckSolvent: wrong octree type\n" ));
		exit(1);
	}
	/*
	 *  Set up globals for octree.
	 */
	vNewPoint = *PvPoint;
	iMaxDepth = octTree->iMaxDepth;
	dGridSize = octTree->dGridSize;

	/*
	 *  Initialize closest atom for comparison.
	 */

	PaAtom = PVAI( octTree->vaAtoms, ATOM, 0 ); 
	aClosestAtom = *PaAtom;
	dClosestDistance = dDistanceSq( &vNewPoint,
					&vAtomPosition( aClosestAtom ) );

	/*
	 *  Descend octree.
	 */
	if (!OctNodeCheckSolvent( &octTree->onHead ) ) {
		VP0(("Completely out of solvent bounding area\n"));
		return(NULL);
	}

	/*
	 *  Check if closest atom overlaps at all.
	 */
	if ( dClosestDistance < 9.0 ) {	/* HACK - approx ion+wat */
		*PvPoint = vNewPoint;
		return( (RESIDUE) cContainerWithin( aClosestAtom ) );
	}
	VP0(("No overlap w/ solvent\n"));
	return(NULL);
}

