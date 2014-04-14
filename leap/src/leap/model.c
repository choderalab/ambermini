int		this = 0;
/*
 *      File:   model.c
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
 *              The MODEL object is an object with a single instance
 *              that interacts with the BUILD object and generates
 *              model bond lengths, angles, torsions etc which
 *              are determined from atom types.
 */




#include	"basics.h"

#include        "classes.h"

#include	"chirality.h"

#include        "model.h"

#include        "sort.h"

#include        "atom.h"

#include        "zMatrix.h"



/*
--------------------------------------------------------------------

        Private routines

*/

		/* MODELTORSIONt is passed to the modelling routines */
		/* to build INTERNAL torsions */
 
typedef	struct	{
		ATOM		aAtom;
		BOOL		bPosKnown;
		VECTOR		vPos;
		BOOL		bBuildInternals;
		} MODELATOMt;


typedef	struct	{
		int		iXFirstUnknown;
		int		iXBonds;
		MODELATOMt	maaXBonds[MAXBONDS];
		double		dXOrientation;
		int		iXHybridization;
		MODELATOMt	maX;
		double		dYOrientation;
		int		iYHybridization;
		MODELATOMt	maY;
		int		iYFirstUnknown;
		int		iYBonds;
		MODELATOMt	maaYBonds[MAXBONDS];
		double		dAbsolute;
		} MODELTORSIONt;





/*
 *	ziModelAtomWeight
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the 'weight' of the atom.  This is used for
 *	forcing the 'heaviest' atoms around a torsion trans to
 *	each other.
 *	The 'weight' of an atom is currently just its element
 *	number, unless the atom is CARBON, then it is 1000,
 *	making it the 'heaviest' atom.
 */
static int
ziModelAtomWeight( ATOM aAtom )
{
    if ( iAtomElement(aAtom) == CARBON ) 
	return(1000);
    return(iAtomElement(aAtom));
}




/*
 *	ModelTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 */
static void
ModelTorsion( MODELTORSIONt *mtPTorsions, int iBondX, int iBondY, double dVal )
{
MODELATOMt	*maPA, *maPX, *maPY, *maPD;
#ifdef DEBUG
STRING		s1, s2, s3, s4;
#endif

    if ( iBondX >= mtPTorsions->iXBonds ||
	 iBondY >= mtPTorsions->iYBonds ) {
	 return;
    }
/*
if (this)
fprintf(stderr, "ibond  %d %d\n", iBondX, iBondY);
*/

    maPA = &(mtPTorsions->maaXBonds[iBondX]);
    maPX = &(mtPTorsions->maX);
    maPY = &(mtPTorsions->maY);
    maPD = &(mtPTorsions->maaYBonds[iBondY]);
    if ( !(maPA->bBuildInternals || maPD->bBuildInternals) ) {
/*
if (this)
fprintf(stderr, " %s needs internals\n",
sAtomName(maPA->aAtom));
*/
	return;
    }
    
	/* If the coordinates for the atoms are defined then */
	/* measure the torsion angle between them and use that for */
	/* the internal */

    if ( maPA->bPosKnown &&
		maPX->bPosKnown &&
		maPY->bPosKnown &&
		maPD->bPosKnown ) {
/*
if (this)
fprintf(stderr, " %s replacing dval\n",
sAtomName(maPA->aAtom));
*/
	dVal = dVectorAtomTorsion( &(maPA->vPos),
				   &(maPX->vPos),
				   &(maPY->vPos),
				   &(maPD->vPos) );
    }

	/* Create an INTERNAL only if there isn't one already */

    if ( !iInternalFindTorsion( maPA->aAtom, 
				maPX->aAtom, 
				maPY->aAtom, 
				maPD->aAtom ) ) {
/*
if (this)
fprintf(stderr, " %s adding internal %f\n", sAtomName(maPA->aAtom),
dVal/DEGTORAD);
*/
	iInternalTorsion( maPA->aAtom,
			  maPX->aAtom,
			  maPY->aAtom,
			  maPD->aAtom,
			  dVal );
        MESSAGE(( "++++Torsion INTERNAL: %lf to %s - %s - %s - %s\n",
		dVal/DEGTORAD,
		sContainerFullDescriptor((CONTAINER)maPA->aAtom,s1),
		sContainerFullDescriptor((CONTAINER)maPX->aAtom,s2),
		sContainerFullDescriptor((CONTAINER)maPY->aAtom,s3),
		sContainerFullDescriptor((CONTAINER)maPD->aAtom,s4) ));
    } else {
	MESSAGE(( "Torsional INTERNAL already exists\n" ));
    }
}



/*
 *	ModelCreateSp3Sp3Torsions
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create the torsion INTERNALs for all torsions
 *	contained within the mtPTorsions record.
 *
 *	The orientation of the ATOMs around the central
 *	ATOMs defines which torsions should be used to
 *	assign INTERNALs to those torsions that are not
 *	yet defined.
 */
static void
ModelCreateSp3Sp3Torsions( MODELTORSIONt *mtPTorsions )
{
double		dADOffset, d180, dm60, d60;

	/* First twist the torsion so that the AD torsion has */
	/* the same absolute angle that is measured */
	/* and twist all the others with it */

    dADOffset = mtPTorsions->dAbsolute - (180.0*DEGTORAD);
    d180 = 180.0*DEGTORAD + dADOffset;
    dm60 = -60.0*DEGTORAD + dADOffset;
    d60  =  60.0*DEGTORAD + dADOffset;
/*
if (this)
fprintf(stderr, "dXOrientation %f  dYOrientation %f  offset %f\n",
mtPTorsions->dXOrientation, mtPTorsions->dYOrientation, dADOffset);
*/

    if ( mtPTorsions->dXOrientation > 0.0 ) {
	if ( mtPTorsions->dYOrientation > 0.0 ) {
		ModelTorsion( mtPTorsions, 0, 0, d180 );
		ModelTorsion( mtPTorsions, 0, 1, dm60 );
		ModelTorsion( mtPTorsions, 0, 2, d60 );
		ModelTorsion( mtPTorsions, 1, 0, dm60 );
		ModelTorsion( mtPTorsions, 1, 1, d60 );
		ModelTorsion( mtPTorsions, 1, 2, d180 );
		ModelTorsion( mtPTorsions, 2, 0, d60 );
		ModelTorsion( mtPTorsions, 2, 1, d180 );
		ModelTorsion( mtPTorsions, 2, 2, dm60 );
	} else {
		ModelTorsion( mtPTorsions, 0, 0,  d180 );
		ModelTorsion( mtPTorsions, 0, 1,  d60 );
		ModelTorsion( mtPTorsions, 0, 2,  dm60 );
		ModelTorsion( mtPTorsions, 1, 0,  dm60 );
		ModelTorsion( mtPTorsions, 1, 1,  d180 );
		ModelTorsion( mtPTorsions, 1, 2,  d60 );
		ModelTorsion( mtPTorsions, 2, 0,  d60 );
		ModelTorsion( mtPTorsions, 2, 1,  dm60 );
		ModelTorsion( mtPTorsions, 2, 2,  d180 );
	}
    } else {
	if ( mtPTorsions->dYOrientation > 0.0 ) {
		ModelTorsion( mtPTorsions, 0, 0, d180 );
		ModelTorsion( mtPTorsions, 0, 1, dm60 );
		ModelTorsion( mtPTorsions, 0, 2, d60 );
		ModelTorsion( mtPTorsions, 1, 0, d60 );
		ModelTorsion( mtPTorsions, 1, 1, d180 );
		ModelTorsion( mtPTorsions, 1, 2, dm60 );
		ModelTorsion( mtPTorsions, 2, 0, dm60 );
		ModelTorsion( mtPTorsions, 2, 1, d60 );
		ModelTorsion( mtPTorsions, 2, 2, d180 );
	} else {
		ModelTorsion( mtPTorsions, 0, 0, d180 );
		ModelTorsion( mtPTorsions, 0, 1, d60 );
		ModelTorsion( mtPTorsions, 0, 2, dm60 );
		ModelTorsion( mtPTorsions, 1, 0, d60 );
		ModelTorsion( mtPTorsions, 1, 1, dm60 );
		ModelTorsion( mtPTorsions, 1, 2, d180 );
		ModelTorsion( mtPTorsions, 2, 0, dm60 );
		ModelTorsion( mtPTorsions, 2, 1, d180 );
		ModelTorsion( mtPTorsions, 2, 2, d60 );
	}
    }
}



/*
 *	ModelCreateSp3Sp2Torsions
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create the torsion INTERNALs for an SP3-SP2 linkage.
 *
 *	The orientation of the ATOMs around the central
 *	ATOMs defines which torsions should be used to
 *	assign INTERNALs to those torsions that are not
 *	yet defined.
 */
static void
ModelCreateSp3Sp2Torsions( MODELTORSIONt *mtPTorsions )
{
double	d180, dm60, d60, d120, dm120, d0;
double	dADOffset;

	/* First twist the torsion so that the AD torsion has */
	/* the same absolute angle that is measured */
	/* and twist all the others with it */


    MESSAGE(( "dADAbsolute = %lf\n", 
			mtPTorsions->dAbsolute/DEGTORAD ));

    dADOffset = mtPTorsions->dAbsolute - (180.0*DEGTORAD);
    d180 =  180.0*DEGTORAD + dADOffset;
    dm60 =  -60.0*DEGTORAD + dADOffset;
    d60  =   60.0*DEGTORAD + dADOffset;
    dm120= -120.0*DEGTORAD + dADOffset;
    d120 =  120.0*DEGTORAD + dADOffset;
    d0   =    0.0*DEGTORAD + dADOffset;

    if ( mtPTorsions->dXOrientation > 0.0 ) {
	ModelTorsion( mtPTorsions, 0, 0, d180 );
	ModelTorsion( mtPTorsions, 0, 1, d0 );
	ModelTorsion( mtPTorsions, 1, 0, dm60 );
	ModelTorsion( mtPTorsions, 1, 1, d120 );
	ModelTorsion( mtPTorsions, 2, 0, d60 );
	ModelTorsion( mtPTorsions, 2, 1, dm120 );
    } else {
	ModelTorsion( mtPTorsions, 0, 0, d180 );
	ModelTorsion( mtPTorsions, 0, 1, d0 );
	ModelTorsion( mtPTorsions, 1, 0, d60 );
	ModelTorsion( mtPTorsions, 1, 1, dm120 );
	ModelTorsion( mtPTorsions, 2, 0, dm60 );
	ModelTorsion( mtPTorsions, 2, 1, d120 );
    }
}






/*
 *	ModelCreateSp2Sp2Torsions
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create the torsion INTERNALs for an SP3-SP2 linkage.
 */
static void
ModelCreateSp2Sp2Torsions( MODELTORSIONt *mtPTorsions )
{
double	d180, d0, dADOffset;

	/* First twist the torsion so that the AD torsion has */
	/* the same absolute angle that is measured */
	/* and twist all the others with it */

    dADOffset = mtPTorsions->dAbsolute - (180.0*DEGTORAD);
    d180 =  180.0*DEGTORAD + dADOffset;
    d0   =    0.0*DEGTORAD + dADOffset;

    ModelTorsion( mtPTorsions, 0, 0, d180 );
    ModelTorsion( mtPTorsions, 0, 1, d0 );
    ModelTorsion( mtPTorsions, 1, 0, d0 );
    ModelTorsion( mtPTorsions, 1, 1, d180 );
}


/*
 *	dModelCalculateOrientation
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	The ATOMs in maPA, maPY, maPB have been ordered
 *	by zModelOrderAtoms.
 *	Calculate the orientation of the maPB ATOM with
 *	respect to the triangle (maPA - maPX - maPY).
 *	This orientation tells the CreateSp3Sp3 and CreateSp3Sp2 routines
 *	which torsion values to use.
 */
static double
dModelCalculateOrientation( MODELATOMt *maPX, MODELATOMt *maPA, 
		MODELATOMt *maPY, MODELATOMt *maPB )
{
double		dOrientation, dChirality;
INTERNAL	iChirality;


    if ( maPX->aAtom == NULL ||
	 maPA->aAtom == NULL ||
	 maPY->aAtom == NULL ||
	 maPB->aAtom == NULL ) return(1.0);

    if ( maPX->bPosKnown &&
	 maPA->bPosKnown &&
	 maPY->bPosKnown &&
	 maPB->bPosKnown ) {
	dOrientation = dVectorAtomChirality( &(maPX->vPos),
						&(maPA->vPos),
						&(maPY->vPos),
						&(maPB->vPos) );
    } else {

		/* If the chirality cannot be obtained from */
		/* the coordinates then get it from */
		/* the INTERNAL if it is defined. */
		/* If the INTERNAL isn't defined then use */
		/* an arbitrary chirality of 1.0 */

	iChirality = iInternalFindChirality(maPX->aAtom);
	if ( iChirality != NULL ) {
	    dChirality = dInternalValue(iChirality);
	} else {
/*
if (this)
fprintf(stderr, "default chirality on %s\n",
sAtomName(maPX->aAtom));
*/
	    dChirality = 1.0;
	}

		/* The orientation just calculated is wrt the ATOM */
		/* ordering scheme in 'atom.c'.  Transform it */
		/* to the ordering scheme maPA, maPY, maPB */

		/* Do the transformation by swapping ATOMs in */
		/* one ordering scheme until we get the second */
		/* ordering scheme, each time flipping the orientation */
		/* value */

	dOrientation = dChiralityToOrientation( dChirality,
				maPX->aAtom,
				maPA->aAtom,
				maPY->aAtom,
				maPB->aAtom,
				NULL );

    }

    return(dOrientation);
}





/*
 *	zbModelTrueIfPosKnown
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return TRUE if the atom's position is known.
 *	An ATOMs position is known if the SfAtomPositionKnownFlag
 *	is set.
 */
static BOOL
zbModelTrueIfPosKnown( MODELATOMt *maPAtom )
{
BOOL	bRes;

    bRes = maPAtom->bPosKnown;
    return(bRes);
}






/*
 *	zModelOrderAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Order the atoms bonded to aX (other than aY) into
 *	the array aaX[].  Sift the atoms in aaX by whether or not
 *	they have their coordinates defined.  Those with coordinates
 *	defined are first in the array and those without are last.
 *
 *	Put the number of atoms in (iX)
 *	and the index of the first atom with undefined coordinates in
 *	*iPFirstXUnknown.  Within each group of atoms place the
 *	atom with the largest 'weight' at the start of the group.
 *	This forces the 'heaviest' atoms to be trans to each other.
 */
static void
zModelOrderAtoms( MODELATOMt *maPX, MODELATOMt *maPY, MODELATOMt maaXBonds[],
		int iX, int *iPFirstXUnknown )
{
int		i, iHighest, iPos, iWeight;
MODELATOMt	maTemp;

		/* Sift the atoms by whether or not they have defined */
		/* coordinates */


    Sift( maaXBonds, sizeof(maaXBonds[0]), iX, 
		zbModelTrueIfPosKnown, iPFirstXUnknown );

		/* Find the 'heaviest' atom in the second half of the */
		/* list and put it at the front */

    iHighest = 0;
    for ( i=0; i<*iPFirstXUnknown; i++ ) {
	iWeight = ziModelAtomWeight(maaXBonds[i].aAtom);
	if ( iHighest < iWeight ) {
	    iHighest = iWeight;
	    iPos = i;
	}
    }
    if (iHighest!=0) SWAP( maaXBonds[0], maaXBonds[iPos], maTemp );

		/* Find the 'heaviest' atom in the second half of the */
		/* list and put it at the front */

    iHighest = 0;
    for ( i=*iPFirstXUnknown; i<iX; i++ ) {
	iWeight = ziModelAtomWeight(maaXBonds[i].aAtom);
	if ( iHighest < iWeight ) {
	    iHighest = iWeight;
	    iPos = i;
	}
    }
    if ( iHighest != 0 ) 
	SWAP( maaXBonds[*iPFirstXUnknown], maaXBonds[iPos], maTemp );
}







/*
 *	zModelBuildMockExternals
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build mock external coordinates using the torsional INTERNALs
 *	in iaTorsions.  Build the coordinates for the ATOMs around
 *	(a1) and (a2) and set the flag (fKnownSet) for the ATOMs
 *	that got mock coordinates, reset it for the rest.
 */
static void
zModelBuildMockExternals( MODELTORSIONt *mtPTorsions, 
		INTERNAL iaTorsions[], int iTorsions )
{
VECTOR		vTemp;
INTERNAL	iInt;
MODELATOMt	*maPNew, *maPCur;
MODELATOMt	*maPC1, *maPC2, *maPC3;
MODELATOMt	*maPTemp1, *maPTemp2, *maPTemp3, *maPTemp4;
BOOL		bGotOne;
int		i, j, iLeft;
#ifdef DEBUG
STRING		s1, s2, s3, s4;
#endif

		/* Define coordinates for the central ATOMs */

    VectorDef( &vTemp, 0.0, 0.0, 0.0 );
    mtPTorsions->maX.vPos = vTemp;
    mtPTorsions->maX.bPosKnown = TRUE;

    VectorDef( &vTemp, 1.0, 0.0, 0.0 );
    mtPTorsions->maY.vPos = vTemp;
    mtPTorsions->maY.bPosKnown = TRUE;

		/* Tell the outer ATOMs that they dont have */
		/* positions defined */

    for ( i=0; i<mtPTorsions->iXBonds; i++ ) {
	mtPTorsions->maaXBonds[i].bPosKnown = FALSE;
    }
    for ( i=0; i<mtPTorsions->iYBonds; i++ ) {
	mtPTorsions->maaYBonds[i].bPosKnown = FALSE;
    }

		/* Place the first ATOM in the XY plane */

    iLeft = iTorsions;
    iInt = iaTorsions[0];
    VectorDef( &vTemp, 1.0, 1.0, 0.0 );
    maPCur = ((MODELATOMt*)PAtomTempPtr(aInternalAtom1(iInt)));
    maPCur->vPos = vTemp;
    maPCur->bPosKnown = TRUE;

    MESSAGE(( "=======  Started mock coords from: %s\n",
		sContainerFullDescriptor((CONTAINER)aInternalAtom1(iInt),s1) ));

		/* Now start looping through the torsions, looking for */
		/* those that have one ATOM defined and then build the */
		/* coordinate for the other ATOM and clear the INTERNAL */
		/* from the ARRAY */

MESSAGEEXECUTE( {
    MESSAGE(( "========  %d Torsions to build mock coords from:\n",
		iTorsions ));
    for ( i=0; i<iTorsions; i++ ) {
	MESSAGE(( "------- Known torsion: %s - %s - %s - %s\n",
		sContainerFullDescriptor((CONTAINER)aInternalAtom1(iaTorsions[i]),s1),
		sContainerFullDescriptor((CONTAINER)aInternalAtom2(iaTorsions[i]),s2),
		sContainerFullDescriptor((CONTAINER)aInternalAtom3(iaTorsions[i]),s3),
		sContainerFullDescriptor((CONTAINER)aInternalAtom4(iaTorsions[i]),s4) ));
    }
 		} );

    for ( i=0; i<iTorsions; i++ ) {
	bGotOne = FALSE;
	for ( j=0; j<iTorsions; j++ ) {
	    iInt = iaTorsions[j];
	    if ( iInt == NULL ) continue;
	    maPTemp1 = (MODELATOMt*)PAtomTempPtr(aInternalAtom1(iInt));
	    maPTemp2 = (MODELATOMt*)PAtomTempPtr(aInternalAtom2(iInt));
	    maPTemp3 = (MODELATOMt*)PAtomTempPtr(aInternalAtom3(iInt));
	    maPTemp4 = (MODELATOMt*)PAtomTempPtr(aInternalAtom4(iInt));
	    if ( maPTemp4->bPosKnown ) {
		maPNew= maPTemp1;
		maPC1 = maPTemp2;
		maPC2 = maPTemp3;
		maPC3 = maPTemp4;
	    } else if ( maPTemp1->bPosKnown ) {
		bGotOne = TRUE;
		maPNew= maPTemp4;
		maPC1 = maPTemp3;
		maPC2 = maPTemp2;
		maPC3 = maPTemp1;
	    } else continue;
	    bGotOne = TRUE;

	    MESSAGE(( "======= Building mock coord for: %s\n",
			sContainerFullDescriptor((CONTAINER)maPNew->aAtom,s1) ));
	    MESSAGE(( "======= Using torsion: %s - %s - %s - %s\n",
		sContainerFullDescriptor((CONTAINER)aInternalAtom1(iaTorsions[j]),s1),
		sContainerFullDescriptor((CONTAINER)aInternalAtom2(iaTorsions[j]),s2),
		sContainerFullDescriptor((CONTAINER)aInternalAtom3(iaTorsions[j]),s3),
		sContainerFullDescriptor((CONTAINER)aInternalAtom4(iaTorsions[j]),s4) ));

			/* Now build the coordinate for aNew */

	    ZMatrixBondAngleTorsion( &(maPNew->vPos),
				     &(maPC1->vPos),
				     &(maPC2->vPos),
				     &(maPC3->vPos),
				     1.0,
				     90.0*DEGTORAD,
				     dInternalValue(iInt) );
	    maPNew->bPosKnown = TRUE;
	    iaTorsions[j] = NULL;
	    iLeft--;
	    break;
	}
	if ( !bGotOne ) {
	    DFATAL(( "There are %d torsions left over for mock coords\n",
			iLeft ));
	}
    }
}




/*
 *===================================================================
 *
 *        Public routines
 *
 */





/*
 *	ModelAssignTorsionsAround
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Assign torsion angles around the two atoms.
 *
 *NOTE:	This routine is VERY complicated.
 *
 *	The routine works on three cases:
 *
 *1)	If no coordinates are defined for (aX), (aY) and the
 *	ATOMs around them, then it simply generates reasonable
 *	torsion INTERNALs for all of the torsions around the
 *	central ATOMs.  The torsions are generated based on
 *	the hybridization of (aX) and (aY).
 *
 *2)	If there are coordinates defined for (aX) and (aY) and
 *	there are a few ATOMs around (aX) and (aY) with positions
 *	defined then torsions are calculated for the remaining
 *	ATOMs based on the hybridization, and chirality of
 *	(aX) and (aY).
 *
 *3)	The most difficult case.  If there are some INTERNALs
 *	defined for torsions across (aX) and (aY) then
 *	temporary coordinates are calculated for the ATOMs 
 *	connected to (aX) and (aY) that have the torsion
 *	INTERNALs (calculated from the INTERNAL torsions and
 *	chiralities of (aX), (aY)) and then case (2) is applied.
 *
 *	Reasonable torsion angles are found using the following
 *	recipe:
 *		For SP3-SP3, the central atoms are called X - Y
 *		the atoms around X are called A, B, C and
 *		those around Y are called D, E, F.
 *		There are four ways to orient the atoms around
 *		X and Y depending on the chirality of X and Y.
 *		All four ways are represented in the function
 *		(zModelCreateSp3Sp3Torsions).  This routine
 *		determines which of the four the current
 *		torsion is by calculating the chiralities of X
 *		and Y and then imposing the torsions for that
 *		pair of chiralities.  The chiralities are calculated
 *		using the external coordinates of the ATOMs
 *		if they are defined, or by building mock
 *		external coordinates from the internal coordinates
 *		if they are defined, or by randomly assigning
 *		+1, +1 as the chiralities if nothing is
 *		defined.
 *
 *NOTE:	THIS ROUTINE AND THOSE THAT ARE CALLED BY THIS ROUTINE
 *	USE THE ATOM TEMPORARY POINTER.
 *
 */
void
ModelAssignTorsionsAround( ATOM aX, ATOM aY, ATOM aOnly )
{
int		iHX, iHY;
BOOL		bKnownX, bKnownY;
int		iTorsions, i, iTemp;
INTERNAL	iaTorsions[MAXTORSIONSAROUNDBOND];
ATOM		aTemp;
MODELTORSIONt	mtTorsions;
MODELATOMt	*maPAtom;
STRING		s1;
#ifdef DEBUG
STRING s2;
#endif


    iHX = iAtomHybridization(aX);
    iHY = iAtomHybridization(aY);

		/* Make sure that the hybridization of (aX) is */
		/* greater than that of (aY) because routines */
		/* called by this one require this */

    if ( iHX < iHY ) {
	SWAP( aX, aY, aTemp );
	SWAP( iHX, iHY, iTemp );
    }

if (!strcmp(sAtomName(aX), "C1") && !strcmp(sAtomName(aY), "O3"))
this++;

		/* First check what is going to be done to fit the */
		/* torsions generated with the existing structure around */
		/* the central atoms */

		/* Check if there is at least one ATOM on either side */
		/* of the central pair that has ATOMPOSITIONKNOWN */
		/* then use the fixed coordinates to deduce how to fit */
		/* the new torsions in.  If there are INTERNAL torsions */
		/* defined around the central pair then use those to generate */
		/* mock coordinates and set the global variable that tells */
		/* the rest of the 'model' code to use ATOMMOCKPOSITIONKNOWN */
		/* to figure out how to fit the new torsions in */

		/* Also set up the MODELTORSIONt record */

    bKnownX = FALSE;
    mtTorsions.iXBonds = iAtomCoordination(aX) - 1;
    maPAtom = &(mtTorsions.maaXBonds[0]);
    for ( i=0; i<iAtomCoordination(aX); i++ ) {
	aTemp = aAtomBondedNeighbor(aX,i);
	if ( aTemp != aY ) {
	    AtomSetTempPtr( aTemp, maPAtom );
	    maPAtom->aAtom = aTemp;
	    maPAtom->vPos = vAtomPosition(aTemp);
	    if ( (maPAtom->bPosKnown = 
		  bAtomFlagsSet(aTemp,ATOMPOSITIONKNOWN)) ) bKnownX = TRUE;
	    if ( aOnly == NULL ) maPAtom->bBuildInternals = TRUE;
	    else 		 maPAtom->bBuildInternals = ( aOnly==aTemp);
	    maPAtom++;
	}
    }
    if ( iAtomCoordination(aX) != 
          1 + maPAtom - &(mtTorsions.maaXBonds[0]) ) {
	VP0(( "Error: Atom %s has force field coordination %i\n"
	      "       but only %i bonded neighbors.\n"
	      "       The cause may be an incorrect atom type, and\n"
	      "       the effect may be a crash very soon.\n",
	      sContainerFullDescriptor((CONTAINER)aX,s1),
	      iAtomCoordination(aX),
	      1 + maPAtom - &(mtTorsions.maaXBonds[0]) ));
    }
    for ( i=iAtomCoordination(aX); i<MAXBONDS; i++ ) {
	maPAtom->aAtom = NULL;
	maPAtom++;
    }

    mtTorsions.maX.aAtom = aX;
    mtTorsions.maX.vPos = vAtomPosition(aX);
    mtTorsions.maX.bPosKnown = bAtomFlagsSet(aX,ATOMPOSITIONKNOWN);
    mtTorsions.iXHybridization = iAtomHybridization(aX);
    AtomSetTempPtr( aX, &(mtTorsions.maX) );

    mtTorsions.maY.aAtom = aY;
    mtTorsions.maY.vPos = vAtomPosition(aY);
    mtTorsions.maY.bPosKnown = bAtomFlagsSet(aY,ATOMPOSITIONKNOWN);
    mtTorsions.iYHybridization = iAtomHybridization(aY);
    AtomSetTempPtr( aY, &(mtTorsions.maY) );

    bKnownY = FALSE;
    mtTorsions.iYBonds = iAtomCoordination(aY) - 1;
    maPAtom = &(mtTorsions.maaYBonds[0]);
    for ( i=0; i<iAtomCoordination(aY); i++ ) {
	aTemp = aAtomBondedNeighbor(aY,i);
	if ( aTemp != aX ) {
	    AtomSetTempPtr( aTemp, maPAtom );
	    maPAtom->aAtom = aTemp;
	    maPAtom->vPos = vAtomPosition(aTemp);
	    if ( (maPAtom->bPosKnown = 
		  bAtomFlagsSet(aTemp,ATOMPOSITIONKNOWN)) ) bKnownY = TRUE;
	    if ( aOnly == NULL ) maPAtom->bBuildInternals = TRUE;
	    else 		 maPAtom->bBuildInternals = ( aOnly==aTemp);
	    maPAtom++;
	}
    }
    if ( iAtomCoordination(aY) != 
          1 + maPAtom - &(mtTorsions.maaYBonds[0]) ) {
	VP0(( "Error: Atom %s has force field coordination %i\n"
	      "       but only %i bonded neighbors.\n"
	      "       The cause may be an incorrect atom type, and\n"
	      "       the effect may be a crash very soon.\n",
	      sContainerFullDescriptor((CONTAINER)aY,s1),
	      iAtomCoordination(aY),
	      1 + maPAtom - &(mtTorsions.maaYBonds[0]) ));
    }
    for ( i=iAtomCoordination(aY); i<MAXBONDS; i++ ) {
	maPAtom->aAtom = NULL;
	maPAtom++;
    }

		/* If they are not both known then try looking for INTERNALs */
		/* that will define the outer ATOM coordinates */

    if ( !( bKnownX && bKnownY ) ) {
	iTorsions = iInternalFindAllTorsionInternalsAround( aX, aY, 
							iaTorsions );
/*
if(this)
fprintf(stderr, "not both %d \n", iTorsions);
*/
	if ( iTorsions != 0 ) {

	    MESSAGE(( "Using INTERNALs to fit new torsions around: %s - %s\n",
			sContainerFullDescriptor((CONTAINER)aX,s1),
			sContainerFullDescriptor((CONTAINER)aY,s2) ));

			/* There are torsion INTERNALs around the */
			/* central ATOMs, build mock external coordinates */
			/* using them.  This will allow the routine to */
			/* order the ATOMs around (aX), (aY) and then */
			/* assign the proper torsion INTERNALs for the */
			/* ATOMs that lack them */

	    zModelBuildMockExternals( &mtTorsions, iaTorsions, iTorsions );

	} else {
	    MESSAGE(( "Completely free in assigning new torsions for: %s - %s\n",
			sContainerFullDescriptor((CONTAINER)aX,s1),
			sContainerFullDescriptor((CONTAINER)aY,s2) ));
	}
    } else {

/*
if(this)
fprintf(stderr, "using ext\n");
*/
	MESSAGE(( "Using externals to fit new torsions around: %s - %s\n",
			sContainerFullDescriptor((CONTAINER)aX,s1),
			sContainerFullDescriptor((CONTAINER)aY,s2) ));

    }
		/* Order the ATOMs around (aX) */
		/* Sift them so that the ATOMs with positions known */
		/* appear at the front of the list */

    zModelOrderAtoms( &(mtTorsions.maX),
		      &(mtTorsions.maY),
		      mtTorsions.maaXBonds,
		      mtTorsions.iXBonds,
		      &(mtTorsions.iXFirstUnknown) );

    zModelOrderAtoms( &(mtTorsions.maY),
		      &(mtTorsions.maX),
		      mtTorsions.maaYBonds,
		      mtTorsions.iYBonds,
		      &(mtTorsions.iYFirstUnknown) );

	/* Now calculate the chirality around each atom X,Y */

    mtTorsions.dXOrientation = 0.0;
    mtTorsions.dYOrientation = 0.0;
    if ( mtTorsions.iXHybridization == HSP3 ) {
	mtTorsions.dXOrientation = dModelCalculateOrientation(
					&(mtTorsions.maX),
					&(mtTorsions.maaXBonds[0]),
					&(mtTorsions.maY),
					&(mtTorsions.maaXBonds[1]) );
    }
    if ( mtTorsions.iYHybridization == HSP3 ) {
	mtTorsions.dYOrientation = dModelCalculateOrientation(
					&(mtTorsions.maY),
					&(mtTorsions.maaYBonds[0]),
					&(mtTorsions.maX),
					&(mtTorsions.maaYBonds[1]) );
    }

MESSAGEEXECUTE( {
    MESSAGE(( "Orientation around: %s = %lf\n", 
		sContainerName(mtTorsions.maX.aAtom),
		mtTorsions.dXOrientation ));
    for ( i=0; i<mtTorsions.iXBonds; i++ ) {
	MESSAGE(( "Atom %d: = %s\n", i, 
		sContainerName(mtTorsions.maaXBonds[i].aAtom) ));
    }
    MESSAGE(( "Orientation around: %s = %lf\n", 
		sContainerName(mtTorsions.maY.aAtom),
		mtTorsions.dYOrientation ));
    for ( i=0; i<mtTorsions.iYBonds; i++ ) {
	MESSAGE(( "Atom %d: = %s\n", i, 
		sContainerName(mtTorsions.maaYBonds[i].aAtom) ));
    }
		} );

	/* Calculate the actual torsion angle between A-X-Y-D */

    if ( mtTorsions.maaXBonds[0].bPosKnown &&
	 mtTorsions.maX.bPosKnown &&
	 mtTorsions.maY.bPosKnown &&
	 mtTorsions.maaYBonds[0].bPosKnown ) {
/*
if (this)
fprintf(stderr, "all known - %s %s %s %s\n",
sAtomName(mtTorsions.maaXBonds[0].aAtom),
sAtomName(mtTorsions.maX.aAtom),
sAtomName(mtTorsions.maY.aAtom),
sAtomName(mtTorsions.maaYBonds[0].aAtom));
*/
	mtTorsions.dAbsolute = dVectorAtomTorsion( 
					&(mtTorsions.maaXBonds[0].vPos),
					&(mtTorsions.maX.vPos),
					&(mtTorsions.maY.vPos),
					&(mtTorsions.maaYBonds[0].vPos) );
    } else {
	mtTorsions.dAbsolute = 180.0*DEGTORAD;
/*
if (this)
fprintf(stderr, "all NOT known - %s %c  %s %c  %s %c  %s %c\n",
sAtomName(mtTorsions.maaXBonds[0].aAtom),
(mtTorsions.maaXBonds[0].bPosKnown ? 'y' : 'n'),
sAtomName(mtTorsions.maX.aAtom),
(mtTorsions.maX.bPosKnown ? 'y' : 'n'),
sAtomName(mtTorsions.maY.aAtom),
(mtTorsions.maY.bPosKnown ? 'y' : 'n'),
sAtomName(mtTorsions.maaYBonds[0].aAtom),
(mtTorsions.maaYBonds[0].bPosKnown ? 'y' : 'n')
);
*/

    }

		/* Build new INTERNALs */

    iHX = mtTorsions.iXHybridization;
    iHY = mtTorsions.iYHybridization;
    if ( iAtomCoordination(aX) > 1 &&
	 iAtomCoordination(aY) > 1 ) {
	if ( iHX == HSP3 && iHY == HSP3 ) {
	    ModelCreateSp3Sp3Torsions( &mtTorsions );
	} else if ( iHX==HSP3 && iHY==HSP2 ) {
	    ModelCreateSp3Sp2Torsions( &mtTorsions );
	} else if ( iHX==HSP2 && iHY==HSP2 ) {
	    ModelCreateSp2Sp2Torsions( &mtTorsions ); 
	} else {
	    PRINTF(( "Currently only Sp3-Sp3/Sp3-Sp2/Sp2-Sp2 are supported\n" ));
	    PRINTF(( "---Tried to superimpose torsions for: *-%s-%s-*\n",
			sContainerName(aX), sContainerName(aY) ));
	    PRINTF(( "--- With Sp%d - Sp%d\n", iHX, iHY ));
	    PRINTF(( "--- Sp0 probably means a new atom type is involved\n" ));
	    PRINTF(( "--- which needs to be added via addAtomTypes\n" ));
	}
    }
this = 0;
}









/*
 *      dModelBondLength
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return an ideal bond length for the bond between the two
 *      atoms.
 */
double
dModelBondLength( ATOM aAtom1, ATOM aAtom2 )
{
double		dValue, dK;
STRING		s1, s2, sDesc;
PARMSET		psTemp;
int		iTag;

    if ( bParmLibDefaultExists() ) {
    	PARMLIB_DEFAULT_LOOP( psTemp,
		( iTag = iParmSetFindBond( psTemp, 
						sAtomType(aAtom1),
						sAtomType(aAtom2) ) ) );
	if ( iTag != PARM_NOT_FOUND ) {
	    ParmSetBond( psTemp, iTag, s1, s2, &dK, &dValue, sDesc );
	    return(dValue);
	}
    }

		/* If there is no parameter then return a model parameter */
	
    ModelBondParm( aAtom1, aAtom2, &dK, &dValue );

    MESSAGE(( "Model bond length for: %s - %s = %lf\n",
		sContainerFullDescriptor((CONTAINER)aAtom1,s1), 
		sContainerFullDescriptor((CONTAINER)aAtom2,s2),
		dValue ));
    return(dValue);
}



/*
 *      dModelBondAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return an ideal bond angle for the angle between the three
 *      atoms.
 */
double
dModelBondAngle( ATOM aAtom1, ATOM aAtom2, ATOM aAtom3 )
{
double		dValue, dK, dTkub, dRkub;
STRING		s1, s2, s3, sDesc;
PARMSET		psTemp;
int		iTag;

	/* First look up the atoms in the default PARMLIBRARY */


    if ( bParmLibDefaultExists() ) {
    	PARMLIB_DEFAULT_LOOP( psTemp,
		( iTag = iParmSetFindAngle( psTemp, 
						sAtomType(aAtom1),
						sAtomType(aAtom2),
						sAtomType(aAtom3) ) ) );
	if ( iTag != PARM_NOT_FOUND ) {
	    ParmSetAngle( psTemp, iTag, s1, s2, s3, &dK, &dValue, 
			&dTkub, &dRkub, sDesc );
	    return(dValue);
	}
    }

    ModelAngleParm( aAtom1, aAtom2, aAtom3, &dK, &dValue );

    MESSAGE(( "Model bond angle for: %s - %s - %s  = %lf\n",
		sContainerFullDescriptor((CONTAINER) aAtom1, s1 ),
		sContainerFullDescriptor((CONTAINER) aAtom2, s2 ),
		sContainerFullDescriptor((CONTAINER) aAtom3, s3 ),
		dValue/DEGTORAD ));

    return(dValue);
}





/*
 *	ModelAddHydrogens
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add Hydrogens to all atoms that are lacking them.
 *	The number of hydrogens required is determined by calculating
 *	the number of hydrogens that would be required to fill the
 *	bonding requirements of the different atoms depending
 *	on the element type and hybridization.
 */
void
ModelAddHydrogens( UNIT uUnit )
{
LOOP		lAtoms;
ATOM		aAtom, aNew;
int		i, iNeeded;
STRING		sName;

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {

	MESSAGE(( "Looking at atom: %s\n", sContainerName(aAtom) ));

		/* How many protons does each element take when SP3 */

	switch ( iAtomElement(aAtom) ) {
	    case CARBON:
		iNeeded = 4;
		break;
	    case NITROGEN:
		iNeeded = 3;
		break;
	    case OXYGEN:
		iNeeded = 2;
		break;
	    case FLOURINE:
	    case CHLORINE:
	    case BROMINE:
	    case HYDROGEN:
		iNeeded = 1;
		break;
	    default:
		iNeeded = 0;
		break;
	}

		/* Delete the number of bonds there already are */

	iNeeded -= iAtomCoordination(aAtom);

		/* Take account the hybridization state */

	switch ( iAtomHybridization(aAtom) ) {
	    case 0:
		VP0(( "(%s: hybridization for type %s unknown)\n",
					sContainerName(aAtom), 
					sAtomType( aAtom ) ));
		break;
	    case HSP3:
		iNeeded -= 0;
		break;
	    case HSP2:
		iNeeded -= 1;
		break;
	    case HSP1:
		iNeeded -= 1;
		break;
	    default:
		DFATAL(( "%s: Illegal hybridization (%d; bonds needed: %d)",
					sContainerName(aAtom), 
					iAtomHybridization(aAtom), 
					iNeeded ));
		break;
	}

	MESSAGE(( "Number of protons to add: %d\n", iNeeded ));

	if ( iNeeded > 0 ) {

		/* Create the protons */

	    for ( i=0; i<iNeeded; i++ ) {
		aNew = (ATOM)oCreate(ATOMid);
		AtomDefineFlags( aNew,
				ATOMNEEDSBUILD|ATOMPOSITIONDRAWN );
		AtomSetElement( aNew, HYDROGEN );
		AtomBondTo( aAtom, aNew );
		ContainerAdd( cContainerWithin(aAtom), (OBJEKT)aNew );
		sprintf( sName, "H%d", iContainerSequence(aNew) );
		ContainerSetName( aNew, sName );
	    }
	}
    }
}






/*
 *------------------------------------------------------------------
 *
 *	Return parameters for interactions between the ATOMs
 *
 *	First look in the default PARMLIBRARY, if nothing is found
 *	then return default values based on the hybridization.
 */



typedef	struct	{
		int	iHybrid2, iHybrid3;
		int	iN;
		double	dK;
		double	dE;
		double  dScEE;
		double  dScNB;
		} H_PROPERPARMt;

typedef	struct	{
		int	iHybrid2;
		double	dK;
		double	dE;
		} H_ANGLEPARMt;

typedef	struct	{
		int	iHybrid1, iHybrid2;
		double	dK;
		double	dE;
		} H_BONDPARMt;


	/* Keep iHybrid2 <= iHybrid3 */

#define	TFORCE	20.0
static	H_PROPERPARMt	SppaPropers[] = {
{	HSP3,	HSP3,	3,	1.0,		0.0,	1.2,	2.0 },	/* Non bond */
{	HSP2,	HSP3,	6,	-2.0,		0.0,	1.2,	2.0 },	/* Non bond */
{	HSP2,	HSP2,	2,	-4.0,		0.0,	1.2,	2.0 },	/* Pi bond overlap */
{	HSP1,	HSP1,	1,	0.0,		0.0,	1.2,	2.0 }	/* Not interesting */
};

#define	AFORCE	100.0
#define	A109	109.5*DEGTORAD
#define	A120	120.0*DEGTORAD
#define	A180	180.0*DEGTORAD
static	H_ANGLEPARMt	SapaAngles[] = {
{	HSP3,	AFORCE,	A109 },
{	HSP2,	AFORCE,	A120 },
{	HSP1,	AFORCE,	A180 },
};


	/* Keep iHybrid1 <= iHybrid2 */
#define	BFORCE	400.0
static	H_BONDPARMt	SbpaBonds[] = {
{	HSP3,	HSP3,	BFORCE,	1.5 },
{	HSP2,	HSP3,	BFORCE,	1.4 },
{	HSP1,	HSP3,	BFORCE,	1.3 },
{	HSP2,	HSP2,	BFORCE,	1.35 },
{	HSP1,	HSP2,	BFORCE,	1.3 },
{	HSP1,	HSP1,	BFORCE,	1.1 }
};





/*
 *	ModelProperTerms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Search for model terms for the proper torsion and
 *	add them to the TORSION.
 */
void
ModelProperTerms( ATOM a1, ATOM a2, ATOM a3, ATOM a4, TORSION tTorsion )
{
int		i, iHybrid2, iHybrid3;
STRING		sDesc;

		/* It may no longer be needed to look for special */
		/* proper terms */
		
    iHybrid2 = iAtomHybridization(a2);
    iHybrid3 = iAtomHybridization(a3);

		/* Sort the hybridizations so that iHybrid2<iHybrid3 */
		/* like in the table */

    if ( iHybrid2>iHybrid3) SWAP( iHybrid2, iHybrid3, i );

    for ( i=0; i<sizeof(SppaPropers)/sizeof(SppaPropers[0]); i++ ) {
	if ( iHybrid2 == SppaPropers[i].iHybrid2 &&
	     iHybrid3 == SppaPropers[i].iHybrid3 ) {
	     bParmSetTORSIONAddProperTerm( tTorsion,
	     				sAtomName(a1), sAtomName(a2),
					sAtomName(a3), sAtomName(a3),
					SppaPropers[i].iN,
					SppaPropers[i].dK,
					SppaPropers[i].dE,
					SppaPropers[i].dScEE,
					SppaPropers[i].dScNB,
					sDesc );
	}
    }
}







/*
 *	ModelAngleParm
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a parameter for the angle involving the three ATOMs.
 */
void
ModelAngleParm( ATOM a1, ATOM a2, ATOM a3, double *dPK, double *dPE )
{
int		i, iHybrid2;

    *dPK = 0.0;
    *dPE = 0.0;
    iHybrid2 = iAtomHybridization(a2);

    for ( i=0; i<sizeof(SapaAngles)/sizeof(SapaAngles[0]); i++ ) {
	if ( iHybrid2 == SapaAngles[i].iHybrid2 ) {
	    *dPK = SapaAngles[i].dK;
	    *dPE = SapaAngles[i].dE;
	    return;
	}
    }
}






/*
 *	ModelBondParm
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a parameter for the bond involving the two ATOMs.
 */
void
ModelBondParm( ATOM a1, ATOM a2, double *dPK, double *dPE )
{
int		i, iHybrid1, iHybrid2;

    /*
     *  treat H's as a special case
     */
    if ( iAtomElement(a1) == HYDROGEN  ||  iAtomElement(a2) == HYDROGEN  ) {
	*dPK = BFORCE;
	*dPE = 1.0;
	return;
    }
    *dPK = 0.0;
    *dPE = 0.0;

    iHybrid1 = iAtomHybridization(a1);
    iHybrid2 = iAtomHybridization(a2);

		/* Sort the hybridizations so that iHybrid2<iHybrid3 */
		/* like in the table */

    if ( iHybrid1>iHybrid2) SWAP( iHybrid1, iHybrid2, i );

    for ( i=0; i<sizeof(SbpaBonds)/sizeof(SbpaBonds[0]); i++ ) {
	if ( iHybrid1 == SbpaBonds[i].iHybrid1 &&
	     iHybrid2 == SbpaBonds[i].iHybrid2 ) {
	    *dPK = SbpaBonds[i].dK;
	    *dPE = SbpaBonds[i].dE;
	    return;
	}
    }
    VP0(( "failed to find default bond length %s-%s, types %s-%s\n", 
		sAtomName(a1), sAtomName(a2), 
		sAtomType(a1), sAtomType(a2) ));
}


