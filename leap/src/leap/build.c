/*
 *      File:   build.c
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
 *
 *		This file contains routines which are used to
 *		make changes to the coordinates of ATOMs.
 *
 *		It contains code to generate internal coordinates
 *		from externals, make changes to internals, and
 *		then generate externals from internals.
 *		
 */




#include	"basics.h"

#include        "vector.h"
#include	"matrix.h"
#include        "zMatrix.h"
#include        "classes.h"

#include	"chirality.h"

#include	"minimizer.h"
#include	"parmLib.h"
#include        "model.h"
#include        "build.h"
#include        "sort.h"
#include        "graphUtil.h"

#define         MAXNEWTONSTEPS  20

/*
------------------------------------------------------------------

        Private routines
        
*/



/*
 *	zaBuildFindNextAtomWithNameButNot
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Look over all the neighbors of aA for an ATOM with
 *	the name sName.  Ignore an ATOM if it is the
 *	same ATOM as aB.  If you find one, return it, otherwise
 *	return NULL.
 */
static ATOM
zaBuildFindNextAtomWithNameButNot( ATOM aA, char *sName, ATOM aB )
{
int		i;
ATOM		aNeighbor;

    for ( i=0; i<iAtomCoordination(aA); i++ ) {
	aNeighbor = aAtomBondedNeighbor( aA, i );
	if ( aNeighbor == aB ) continue;
	if ( strcmp(sName,sContainerName(aNeighbor)) == 0 ) {
	    return(aNeighbor);
	}
    }
    return(NULL);
}




/*
 *      BuildExternalWithAll
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinate for the atom when
 *      a bond length, bond angle, and torsion angle is specified.
 */
static void
BuildExternalWithAll( ATOM aAtom, INTERNAL iTorsion, INTERNAL iAngle, 
			INTERNAL iBond )
{
ATOM            aAtom2, aAtom3, aAtom4;
VECTOR          vAtom2, vAtom3, vAtom4;
double          dRadius, dAngle, dTorsion;
VECTOR          vNew;
#ifdef	DEBUG
STRING		s1, s2, s3, s4;
#endif

    if ( aAtom == aInternalAtom1(iTorsion)) {
        aAtom2 = aInternalAtom2(iTorsion);
        aAtom3 = aInternalAtom3(iTorsion);
        aAtom4 = aInternalAtom4(iTorsion);
    } else {
        aAtom2 = aInternalAtom3(iTorsion);
        aAtom3 = aInternalAtom2(iTorsion);
        aAtom4 = aInternalAtom1(iTorsion);
    }

    vAtom2 = vAtomPosition(aAtom2);
    vAtom3 = vAtomPosition(aAtom3);
    vAtom4 = vAtomPosition(aAtom4);
    dRadius = dInternalValue(iBond);
    dAngle = dInternalValue(iAngle);
    dTorsion = dInternalValue(iTorsion);
 
    MESSAGE(( "Building atom %s using torsion/angle/bond\n",
		sContainerFullDescriptor((CONTAINER)aAtom,s1) ));
    MESSAGE(( "Using - %s - %s - %s\n",
		sContainerFullDescriptor((CONTAINER)aAtom2,s2),
                sContainerFullDescriptor((CONTAINER)aAtom3,s3), 
		sContainerFullDescriptor((CONTAINER)aAtom4,s4) ));
    MESSAGE(( "Torsion = %lf\n", dTorsion/DEGTORAD ));
    MESSAGE(( "Angle   = %lf\n", dAngle/DEGTORAD ));
    MESSAGE(( "Bond    = %lf\n", dRadius ));

    ZMatrixBondAngleTorsion( &vNew, &vAtom2, &vAtom3, &vAtom4, 
                                dRadius, dAngle, dTorsion );
    AtomSetPosition( aAtom, vNew );
}







/*
 *      zdBuildCalculateOrientation
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Calculate the orientation of aNew ( above or below the plane )
 *      required to satisfy the chirality at aAtom.
 */
static double
zdBuildCalculateOrientation( ATOM aAtom, ATOM aNew, double dChirality )
{
ATOM    aAtomA, aAtomB, aAtomC, aAtomD;
double  dOrient;

    ChiralityOrderNeighbors( aAtom, &aAtomA, &aAtomB, &aAtomC, &aAtomD );

    dOrient = 0.0;
    if ( aNew == aAtomA ) dOrient = dChirality;
    else if ( aNew == aAtomB ) dOrient = - dChirality;
    else if ( aNew == aAtomC ) dOrient = dChirality;
    else if ( aNew == aAtomD ) dOrient = - dChirality;
    return(dOrient);
}
   


/*
 *      BuildExternalWithTwoAngles
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinate for the atom when
 *      the chirality, a bond length, and two angles are specified.
 *      If the chirality is not known then calculate it.  If
 *      it is still unknown then use positive chirality.
 */
static void
BuildExternalWithTwoAngles( ATOM aAtom, INTERNAL iAngleA, INTERNAL iAngleB, 
				INTERNAL iBond )
{
ATOM            aAtomC, aAtomA, aAtomB;
VECTOR          vAtomC, vAtomA, vAtomB;
double          dBond, dAngleA, dAngleB;
INTERNAL        iChirality;
LOOP            lInternals;
double          dChirality, dOrient;
VECTOR          vNew;


#ifdef	DEBUG
STRING		s1, s2, s3;
#endif

#define SWAPATOM( AA, BB ) { ATOM CC; CC = AA; AA = BB; BB = CC; }


    if ( aAtom == aInternalAtom1(iAngleA)) {
        aAtomC = aInternalAtom2(iAngleA);
        aAtomA = aInternalAtom3(iAngleA);
    } else {
        aAtomC = aInternalAtom2(iAngleA);
        aAtomA = aInternalAtom1(iAngleA);
    }

    if ( aAtom == aInternalAtom1(iAngleB)) aAtomB = aInternalAtom3(iAngleB);
    else                                   aAtomB = aInternalAtom1(iAngleB);

    MESSAGE(( "Building atom %s using two angles\n",
		sContainerFullDescriptor((CONTAINER)aAtom,s1) ));
    MESSAGE(( "Using first-center-second %s - %s - %s\n",
                sContainerFullDescriptor((CONTAINER)aAtomA,s1), 
		sContainerFullDescriptor((CONTAINER)aAtomC,s2),
                sContainerFullDescriptor((CONTAINER)aAtomB,s3) ));

    vAtomC = vAtomPosition(aAtomC);
    vAtomA = vAtomPosition(aAtomA);
    vAtomB = vAtomPosition(aAtomB);
    dBond = dInternalValue(iBond);
    dAngleA = dInternalValue(iAngleA);
    dAngleB = dInternalValue(iAngleB);
    MESSAGE(( "AngleA  = %lf\n", dAngleA/DEGTORAD ));
    MESSAGE(( "AngleB  = %lf\n", dAngleB/DEGTORAD ));
    MESSAGE(( "Bond    = %lf\n", dBond ));
    
                /* Calculate the chirality of aAtomC */
                /* If it already has a chirality then make sure the new */
                /* coordinate is in keeping with that chirality */
                /* Otherwise check if there is an INTERNAL chirality */
                /* defined */ 
    
    dChirality = dChiralityForAtom( aAtomC );
    if ( dChirality != 0.0 ) {
	MESSAGE(( "Got EXTERNAL chirality: %lf\n", dChirality ));
        dOrient = zdBuildCalculateOrientation( aAtomC, aAtom, dChirality );
    } else {
	dChirality = 1.0;
                /* Check to see if aAtomC has an internal chirality defined */

        lInternals = lLoop( (OBJEKT)aAtomC, INTERNALS );
        while ( (iChirality = (INTERNAL)oNext(&lInternals)) ) {
            if ( iInternalType(iChirality) == INTERNALCHIRALITY ) {
                dChirality = dInternalValue(iChirality);
		MESSAGE(( "Got INTERNAL chirality: %lf\n", 
				dChirality ));
	    }
        }
    }

		/* Calculate the orientation of aAtomD with respect */
		/* to aAtomA - aAtomC - aAtomB */

    dOrient = dChiralityToOrientation( dChirality, aAtomC, aAtomA, aAtomB,
						aAtom, NULL );

    MESSAGE(( "The chirality of the ATOM to build is: %lf\n", dChirality ));
    MESSAGE(( "The orientation of the atom to build is: %lf\n", dOrient ));

                /* Now that the orientation is calculated calculate the */
                /* coordinate */
                
    ZMatrixBondTwoAnglesOrientation( &vNew, &vAtomC, &vAtomA, &vAtomB,
                                        dBond, dAngleA, dAngleB, dOrient );
    AtomSetPosition( aAtom, vNew );
}






/*
 *      BuildExternalWithAngleAndBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinate for the atom when
 *      only a bond length and bond angle is specified.
 *      This is done by constraining the coordinates of the atoms
 *      to lie in the plane formed by the Y axis and the 
 *      vector joining the two specified atoms.
 */
static void
BuildExternalWithAngleAndBond( ATOM aAtom, INTERNAL iAngle, INTERNAL iBond )
{
ATOM            aAtom2, aAtom3;
VECTOR          vAtom2, vAtom3;
double          dRadius, dAngle;
VECTOR          vNew;

#ifdef	DEBUG
STRING		s1, s2, s3;
#endif

    if ( aAtom == aInternalAtom1(iBond) ) 
	aAtom2 = aInternalAtom2(iBond);
    else
	aAtom2 = aInternalAtom1(iBond);
    if ( aAtom == aInternalAtom1(iAngle)) 
	aAtom3 = aInternalAtom3(iAngle);
    else
	aAtom3 = aInternalAtom1(iAngle);

    MESSAGE(( "Building atom %s using angle/bond\n",
		sContainerFullDescriptor((CONTAINER)aAtom,s1) ));
    MESSAGE(( "Using - %s - %s\n",
		sContainerFullDescriptor((CONTAINER)aAtom2,s2),
                sContainerFullDescriptor((CONTAINER)aAtom3,s3) ));

    vAtom2 = vAtomPosition(aAtom2);
    vAtom3 = vAtomPosition(aAtom3);
    dRadius = dInternalValue(iBond);
    dAngle = dInternalValue(iAngle);

    MESSAGE(( "Angle    = %lf\n", dAngle/DEGTORAD ));
    MESSAGE(( "Bond     = %lf\n", dRadius ));

    ZMatrixBondAngle( &vNew, &vAtom2, &vAtom3, dRadius, dAngle );
    AtomSetPosition( aAtom, vNew );
}





/*
 *      BuildExternalWithBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinate for the atom when
 *      only a bond length is specified.
 *      This is done by placing the atom a distance along the X axis
 *      as specified by the internal coordinate.
 */
static void
BuildExternalWithBond( ATOM aAtom, INTERNAL iBond )
{
ATOM            aAtom2;
VECTOR          vAtom2, vNew;
double          dBond;

#ifdef	DEBUG
STRING		s1, s2;
#endif

    if ( aAtom == aInternalAtom1(iBond)  ) aAtom2 = aInternalAtom2(iBond);
    else                                   aAtom2 = aInternalAtom1(iBond);

    MESSAGE(( "Building atom %s using bond\n",
		sContainerFullDescriptor((CONTAINER)aAtom,s1) ));
    MESSAGE(( "Using - %s\n",
		sContainerFullDescriptor((CONTAINER)aAtom2,s2) ));

    vAtom2 = vAtomPosition(aAtom2);
    dBond = dInternalValue(iBond);

    MESSAGE(( "Bond     = %lf\n", dBond ));

    ZMatrixBond( &vNew, &vAtom2, dBond );
    AtomSetPosition( aAtom, vNew );
}





/*
 *	BuildInternalsForOneAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build internal coordinates for all bonds, angles, and torsions
 *	that terminate on this ATOM.
 */
static void
BuildInternalsForOneAtom( ATOM aAtom )
{
#ifdef DEBUG
STRING		sTemp;
#endif
INTERNAL	iaTorsions[MAXTORSIONSAROUNDBOND];
ATOM		a1, a2, a3, aB, aC, aTemp;
#ifdef DEBUG
STRING		s1, s2, s3;
#endif
int		i, j, iTorsions, iShouldBe;
double		dValue;
LOOP		lTemp;

    MESSAGE(( "Building internals for: %s\n", 
		sContainerFullDescriptor((CONTAINER)aAtom,sTemp) ));

                /* Build torsion angles */
		/* These have to be built a special way, because */
		/* torsion INTERNALs cannot be simply looked up in a table */
		/* The model builder has code to add torsion INTERNALs */
		/* around central pairs of ATOMs taking into account any */
		/* external coordinates that may be defined already */

    for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
	aB = aAtomBondedNeighbor( aAtom, i );
	for ( j=0; j<iAtomCoordination(aB); j++ ) {
	    aC = aAtomBondedNeighbor( aB, j );
	    if ( aC == aAtom ) continue;
	    MESSAGE(( "Building torsion INTERNALs for: %s  around: %s - %s\n",
			sContainerFullDescriptor((CONTAINER)aAtom,s1),
			sContainerFullDescriptor((CONTAINER)aB,s2),
			sContainerFullDescriptor((CONTAINER)aC,s3) ));
	    if ( aC==aAtom ) continue;
	    iTorsions = iInternalFindAllTorsionInternalsAround( 
					aB, aC, iaTorsions );
	    iShouldBe = (iAtomCoordination(aB)-1)*(iAtomCoordination(aC)-1);
	    if ( iShouldBe != iTorsions ) {
	        ModelAssignTorsionsAround( aB, aC, aAtom );
	    }
	}
    }
		/* Build bond angles */

    LOOPOVERALL( aAtom, ANGLES|ALLOWDUPLICATES, aTemp, ATOM, lTemp ) {
	LoopGetAngle( &lTemp, &a1, &a2, &a3 );
	MESSAGE(( "Building angle INTERNAL for: %s - %s - %s\n",
			sContainerFullDescriptor( (CONTAINER)a1, s1 ),
			sContainerFullDescriptor( (CONTAINER)a2, s2 ),
			sContainerFullDescriptor( (CONTAINER)a3, s3 ) ));
	if ( !iInternalFindAngle( a1, a2, a3 ) ) {
	    if ( bAtomFlagsSet( a1, ATOMPOSITIONKNOWN ) &&
		 bAtomFlagsSet( a2, ATOMPOSITIONKNOWN ) &&
		 bAtomFlagsSet( a3, ATOMPOSITIONKNOWN ) ) {
		dValue = dVectorAtomAngle( &vAtomPosition(a1),
					   &vAtomPosition(a2),
					   &vAtomPosition(a3) );
		MESSAGE(( "Got bond angle from externals\n" ));
	    } else {
		dValue = dModelBondAngle( a1, a2, a3 );
		MESSAGE(( "Got bond angle from model builder\n" ));
	    }
	    MESSAGE(( "++++Angle INTERNAL: %lf  for %s - %s - %s\n", 
		dValue/DEGTORAD,
		sContainerFullDescriptor((CONTAINER)a1,s1),
		sContainerFullDescriptor((CONTAINER)a2,s2),
		sContainerFullDescriptor((CONTAINER)a3,s3) ));
	    iInternalAngle( a1, a2, a3, dValue );
	} else {
	    MESSAGE(( "Angle INTERNAL was already defined\n" ));
	}
    }

		/* Build bond internals */

    LOOPOVERALL( aAtom, BONDS|ALLOWDUPLICATES, aTemp, ATOM, lTemp ) {
	LoopGetBond( &lTemp, &a1, &a2 );
	MESSAGE(( "Building bond INTERNAL for: %s - %s\n",
			sContainerFullDescriptor( (CONTAINER)a1, s1 ),
			sContainerFullDescriptor( (CONTAINER)a2, s2 ) ));
	if ( !iInternalFindBond( a1, a2 ) ) {
            if ( bAtomFlagsSet( a1, ATOMPOSITIONKNOWN ) &&
                 bAtomFlagsSet( a2, ATOMPOSITIONKNOWN ) ) {
		dValue = dVectorAtomLength( &vAtomPosition(a1),
					    &vAtomPosition(a2) );
	        MESSAGE(( "Got bond length from externals\n" ));
	    } else {
                dValue = dModelBondLength( a1, a2 );
		MESSAGE(( "Got bond length from the model builder\n" ));
	    }
            iInternalBond( a1, a2, dValue );
	    MESSAGE(( "++++Bond INTERNAL: %lf  for %s - %s\n", 
		dValue,
		sContainerFullDescriptor((CONTAINER)a1,s1),
		sContainerFullDescriptor((CONTAINER)a2,s2) ));
	} else {
	    MESSAGE(( "Bond length INTERNAL already defined\n" ));
	}
    }

                /* Build chirality */
                /* Chirality is only defined for atoms with 3 or 4 bonds */
                /* Chirality is the only INTERNAL that does not */
                /* need to be specified */

    if ( iInternalFindChirality(aAtom) == NULL ) {
	dValue = dChiralityForAtom(aAtom);
	if ( dValue != 0.0 ) {
	    iInternalChirality( aAtom, dValue );
	    MESSAGE(( "Got chirality from external coordinates\n" ));
	    MESSAGE(( "++++Chirality INTERNAL: %lf  for %s\n", 
		dValue,
		sContainerFullDescriptor((CONTAINER)a1,s1) ));
	} else {
	    MESSAGE(( "Left chirality undefined\n" ));
	}
    } else {
	MESSAGE(( "Chirality is already defined\n" ));
    }
}








/*
 *      BuildExternalForOneAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      If the coordinates of the atom are not fixed then look
 *      for internal coordinates where this atom is a terminal
 *      atom and all other atoms in the internal coordinate are fixed.
 *      Then build the external coordinate for the atom from this
 *      and set the atom flag that the atoms coordinate is BUILD
 */
static void
BuildExternalForOneAtom( ATOM aAtom )
{
LOOP            lInternals;
INTERNAL        iTemp, iBond, iAngle, iTorsion;
INTERNAL        iAngle1, iAngle2;
VECTOR          vPos;
#ifdef	DEBUG
STRING		s1;
#endif

    if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONFIXED ) ) {
	MESSAGE(( "Building one atom\n" ));
        iBond = NULL;
        iAngle = NULL;
        iAngle1 = NULL;
        iAngle2 = NULL;
        iTorsion = NULL;
        
                /* Loop over all internals involving this atom and look */
                /* for a torsion that has this atom as a terminal atom and */
                /* the other three with their coordinates specified. */
                /* Then find the angle and bond that share this atom and */
                /* that lie on the torsion */
                
        lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
        while ( (iTemp = (INTERNAL)oNext(&lInternals)) ) {
            if ( iInternalType(iTemp)==INTERNALTORSION ) {
                if ( bInternalGoodTorsion( iTemp, aAtom,
					&iTorsion, &iAngle, &iBond )) break;
            }
        }
        if ( iTorsion != NULL ) {
            BuildExternalWithAll( aAtom, iTorsion, iAngle, iBond );
            return;
        }
        
                /* Loop over all internals involving this atom and look */
                /* for an angle that has this atom as a terminal atom and */
                /* the other two with their coordinates specified. */
                /* Then look for another angle which has this atom as */
                /* the terminal and shares the center atom with the first */
                /* angle and also has all of its coordinates specified */
        
        lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
        while ( (iTemp = (INTERNAL)oNext(&lInternals)) ) {
            if ( iInternalType(iTemp)==INTERNALANGLE ) {
                if ( bInternalGoodPairOfAngles( iTemp, aAtom, &iAngle1,
                                                &iAngle2, &iBond )) break;
            }
        }
        if ( iAngle2 != NULL ) {
            BuildExternalWithTwoAngles( aAtom, iAngle1, iAngle2, iBond );
            return;
        }
        
                /* Loop over all internals involving this atom and look */
                /* for an angle that has this atom as a terminal atom and */
                /* the other two with their coordinates specified. */
                /* Then find the bond that shares */
                /* this atom and that lies on the angle */
        
        lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
        while ( (iTemp = (INTERNAL)oNext(&lInternals)) ) {
            if ( iInternalType(iTemp)==INTERNALANGLE ) {
                if ( bInternalGoodAngle( iTemp, aAtom, &iAngle, &iBond ))
                    break;
            }
        }
        if ( iAngle != NULL ) {
            BuildExternalWithAngleAndBond( aAtom, iAngle, iBond );
            return;
        }
        
        
                /* Loop over all internals involving this atom and look */
                /* for a bond that has this atom as a terminal atom*/
                
        lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
        while ( (iTemp = (INTERNAL)oNext(&lInternals)) ) {
            if ( iInternalType(iTemp)==INTERNALBOND ) {
                if ( bInternalGoodBond( iTemp, aAtom, &iBond )) break;
            }
        }
        if ( iBond != NULL ) {
            BuildExternalWithBond( aAtom, iBond );
            return;
        }
        
                /* No internal coordinates can be used to determine the */
                /* external coordinate for this atom, so by default */
                /* place it at the origin */
        MESSAGE(( "--Building with nothing, placing %s at origin\n",
			sContainerFullDescriptor( (CONTAINER)aAtom, s1 ) ));

        ZMatrixNothing( &vPos );
        AtomSetPosition( aAtom, vPos );

    }
}




/*
 *	zBuildVerifyInternals
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Loop through all bonds, angles, torsions, and
 *	make sure that there is an INTERNAL for every one.
 */
#if 0
static void
zBuildVerifyInternals( LOOP lAtoms )
{
LOOP		lTemp;
ATOM		a1, a2, a3, a4;
//STRING		s1, s2, s3, s4;
LOOP		lInternals;
INTERNAL	iInt;
int		iCount;
ATOM		aAtom;

		/* First check bonds */

  while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {

    lTemp = lLoop( (OBJEKT)aAtom, BONDS );
    while ( oNext(&lTemp) ) {
	LoopGetBond( &lTemp, &a1, &a2 );
	iCount = 0;
	lInternals = lLoop( (OBJEKT)a1, INTERNALS );
	while ( (iInt = (INTERNAL)oNext(&lInternals)) ) {
	    if ( iInternalType(iInt) == INTERNALBOND && 
		 ( aInternalAtom1(iInt) == a2 ||
		   aInternalAtom2(iInt) == a2 ) ) iCount++;
	}
	if ( iCount != 1 ) {
	    MESSAGE(( "!!!! %d INTERNALBONDs between: %s - %s\n",
			iCount,
			sContainerFullDescriptor((CONTAINER)a1,s1),
			sContainerFullDescriptor((CONTAINER)a2,s2) ));
	}
    }

		/* Now check angles */

    lTemp = lLoop( (OBJEKT)aAtom, ANGLES );
    while ( oNext(&lTemp) ) {
	LoopGetAngle( &lTemp, &a1, &a2, &a3 );
	iCount = 0;
	lInternals = lLoop( (OBJEKT)a1, INTERNALS );
	while ( (iInt = (INTERNAL)oNext(&lInternals)) ) {
	    if ( iInternalType(iInt) == INTERNALANGLE && 
		 aInternalAtom2(iInt) == a2 && 
		( aInternalAtom1(iInt) == a3 ||
		  aInternalAtom3(iInt) == a3 ) ) iCount++;
	}
	if ( iCount != 1 ) {
	    MESSAGE(( "!!!! %d INTERNALANGLEs between: %s-%s-%s\n",
			iCount,
			sContainerFullDescriptor((CONTAINER)a1,s1),
			sContainerFullDescriptor((CONTAINER)a2,s2),
			sContainerFullDescriptor((CONTAINER)a3,s3) ));
	}
    }

		/* Now check torsions */

    lTemp = lLoop( (OBJEKT)aAtom, PROPERS );
    while ( oNext(&lTemp) ) {
	LoopGetTorsion( &lTemp, &a1, &a2, &a3, &a4 );
	iCount = 0;
	lInternals = lLoop( (OBJEKT)a1, INTERNALS );
	while ( (iInt = (INTERNAL)oNext(&lInternals)) ) {
	    if ( iInternalType(iInt) == INTERNALTORSION &&
		( ( aInternalAtom2(iInt) == a2 &&
		    aInternalAtom3(iInt) == a3 &&
		    aInternalAtom4(iInt) == a4 ) ||
		  ( aInternalAtom3(iInt) == a2 &&
		    aInternalAtom2(iInt) == a3 &&
		    aInternalAtom1(iInt) == a4 ) ) ) iCount++; 
	}
	if ( iCount != 1 ) {
	    MESSAGE(( "!!!! %d INTERNALTORSIONs between: %s-%s-%s-%s\n",
			iCount,
			sContainerFullDescriptor((CONTAINER)a1,s1),
			sContainerFullDescriptor((CONTAINER)a2,s2),
			sContainerFullDescriptor((CONTAINER)a3,s3),
			sContainerFullDescriptor((CONTAINER)a4,s4) ));
	}
    }
  }
}
#endif




/*
 *	zbBuildTorsionInternalsForRingGroupMaybe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	If the ring groups satisfies certian criteria.
 *	If the ring group contains a single ring and in
 *	each ring at most one ATOM has ATOMPOSITIONKNOWN
 *	and all ATOMs are SP2 or SP3 then build
 *	INTERNAL torsions for BENZENE or CYCLOHEXANE.
 */
static BOOL
zbBuildTorsionInternalsForRingGroupMaybe( LIST lRingGroup )
{
int		iHybridization, i;
BOOL		bAllSame;
ATOM		aAtom;
ATOM		aaA[6];
int		iKnown;
LISTLOOP	llList;
INTERNAL	iRing;

    if ( iListSize(lRingGroup) != 1 ) return(FALSE);

		/* Get the single ring in the ring group */

    llList = llListLoop(lRingGroup);
    iRing = (INTERNAL)oListNext(&llList);

		/* Loop over all ATOMs, find the hybridization */
		/* and count the number of ATOMs with */
		/* ATOMPOSITIONKNOWN */

    iKnown = 0;
    iHybridization = HUNDEFINED;
    bAllSame = TRUE;
    i = 0;

    InternalRingLoopAtoms(iRing);
    while ( (aAtom = aInternalRingNextAtom(iRing)) ) {
	MESSAGE(( "Looking at atom: %s\n", sContainerName(aAtom) ));
	aaA[i] = aAtom;
	i++;
	if ( iHybridization == HUNDEFINED ) {
	    iHybridization = iAtomHybridization(aAtom);
	}
	if ( iHybridization != iAtomHybridization(aAtom) ) bAllSame = FALSE;
	if ( bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) ) iKnown++;
    }

    if ( iKnown > 1 ) return(FALSE);
    if ( !bAllSame ) return(FALSE);
    if ( i != 6 ) return(FALSE);

		/* If SP2 ring then create INTERNALs for a BENZENE skeleton */

    if ( iHybridization == HSP2 ) {
	MESSAGE(( "Building BENZENE skeleton\n" ));
	iInternalTorsion( aaA[0], aaA[1], aaA[2], aaA[3], 0.0 );
	iInternalTorsion( aaA[1], aaA[2], aaA[3], aaA[4], 0.0 );
	iInternalTorsion( aaA[2], aaA[3], aaA[4], aaA[5], 0.0 );
	iInternalTorsion( aaA[3], aaA[4], aaA[5], aaA[0], 0.0 );
	iInternalTorsion( aaA[4], aaA[5], aaA[0], aaA[1], 0.0 );
	iInternalTorsion( aaA[5], aaA[0], aaA[1], aaA[2], 0.0 );

        return(TRUE);
    }

		/* If SP3 ring then create INTERNALs for a CYCLOHEXANE */

    if ( iHybridization == HSP3 ) {
	MESSAGE(( "Building CYCLOHEXANE skeleton\n" ));

	iInternalTorsion( aaA[0], aaA[1], aaA[2], aaA[3],  60.0*DEGTORAD );
	iInternalTorsion( aaA[1], aaA[2], aaA[3], aaA[4], -60.0*DEGTORAD );
	iInternalTorsion( aaA[2], aaA[3], aaA[4], aaA[5],  60.0*DEGTORAD );
	iInternalTorsion( aaA[3], aaA[4], aaA[5], aaA[0], -60.0*DEGTORAD );
	iInternalTorsion( aaA[4], aaA[5], aaA[0], aaA[1],  60.0*DEGTORAD );
	iInternalTorsion( aaA[5], aaA[0], aaA[1], aaA[2], -60.0*DEGTORAD );

	return(TRUE);
    }
    return(FALSE);
}
	


/*
 *====================================================================
 *====================================================================
 *====================================================================    
 *
 *        Public routines.
 *
 */






/*
 *	BuildInternalsUsingFlags
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build INTERNALs for all ATOMs defined by the LOOP.
 *	Build INTERNALs only for ATOMs that have the flags
 *	(fForSet) set and (fForReset) reset.  Each ATOM that
 *	has INTERNALs build will have (fSetFlags) set and (fResetFlags)
 *	reset.
 */
void	
BuildInternalsUsingFlags( LOOP *lPAtoms, FLAGS fForSet, FLAGS fForReset, 
				FLAGS fSetFlags, FLAGS fResetFlags )
{
ATOM	aAtom;

    LoopUseMemory(lPAtoms);

    while ( (aAtom = (ATOM)oNext(lPAtoms)) ) {
	if ( bAtomFlagsSet( aAtom, fForSet ) &&
		bAtomFlagsReset( aAtom, fForReset ) ) {
	    BuildInternalsForOneAtom( aAtom );
	}
    }

		/* Change the flags */

    LoopRewindMemory(lPAtoms);
    while ( (aAtom=(ATOM)oNext(lPAtoms)) ) {
	if ( bAtomFlagsSet( aAtom, fForSet ) &&
		bAtomFlagsReset( aAtom, fForReset ) ) {
	    AtomSetFlags( aAtom, fSetFlags );
	    AtomResetFlags( aAtom, fResetFlags );
	}
    }

    LoopDestroyMemory(lPAtoms);
}



/*
 *      BuildExternalsUsingFlags
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinates from the internal coordinates.
 *      ALL internal coordinates must be specified, and reasonable,
 *	and preferably internally consistant.
 *
 *      The caller passes a LOOP structure that embodies
 *      a spanning tree over the atoms that are to be built.
 *	The spanning tree is passed instead of the atom to start
 *	on so that the caller can set LOOP flags that make certian
 *	types of ATOMs invisible to the builder, preventing coordinates
 *	being built for those atoms.
 *
 *	The builder works by generating a Z-Matrix for the atoms 
 *	defined in the LOOP.  It traverses the spanning tree
 *	and for each atom searches for a torsion/angle/bond that
 *	contains other atoms whose positions are known.  It
 *	then uses these internal coordinates to build the external
 *	coordinate for that atom.
 *
 *	Only ATOMs that have flags (fForSet) set and (fForReset) reset
 *	will have externals build.  As each external is build,
 *	(fSetFlags) will be set and (fResetFlags) will be reset.
 */
void
BuildExternalsUsingFlags( LOOP *lPLoop, FLAGS fForSet, FLAGS fForReset,
					FLAGS fSetFlags, FLAGS fResetFlags,
					int *iPAddH, int *iPAddHeavy, 
					int *iPAddUnk, BOOL bMsg )
{
ATOM	aAtom;
STRING	sTemp;

                /* Go through each atom and build its coordinates */

    LoopUseMemory(lPLoop); 

    while ( (aAtom=(ATOM)oNext(lPLoop))!= NULL ) {
	if ( bAtomFlagsSet( aAtom, fForSet ) &&
		bAtomFlagsReset( aAtom, fForReset ) ) {
            BuildExternalForOneAtom( aAtom );
	    if ( iAtomElement( aAtom ) == NOELEMENT )
		(*iPAddUnk)++;
	    else if ( iAtomElement( aAtom ) <= HYDROGEN )
		(*iPAddH)++;
	    else {
		(*iPAddHeavy)++;
		if ( bMsg )
		    VP0(("  Added missing heavy atom: %s\n",
			sContainerFullDescriptor( (CONTAINER)aAtom, sTemp ) ));
	    }
	}
    }

		/* Change the flags */

    LoopRewindMemory(lPLoop);
    while ( (aAtom=(ATOM)oNext(lPLoop)) ) {
	if ( bAtomFlagsSet( aAtom, fForSet ) &&
		bAtomFlagsReset( aAtom, fForReset ) ) {
	    MESSAGE(( "Updating flags for: %s\n",
			sContainerFullDescriptor((CONTAINER)aAtom,sTemp) ));
	    AtomSetFlags( aAtom, fSetFlags );
	    AtomResetFlags( aAtom, fResetFlags );
	}
    }

    LoopDestroyMemory(lPLoop);
}




/*
 *	BuildInternalsForContainer
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build INTERNALs for all internals within the container.
 *	Set all of the (fSet) flags on all the ATOMs and reset all of
 *	the (fReset) flags.
 *
 *	This routine is much faster than BuildInternalsWithFlags
 *	because it does not worry about duplicating INTERNALs
 *	because it loops over all internals in the CONTAINER,
 *	not all internals on each ATOM.
 */
void	
BuildInternalsForContainer( CONTAINER cCont, FLAGS fSet, FLAGS fReset )
{
ATOM		aTemp;
ATOM		a1, a2, a3;
LOOP		lTemp;
#ifdef DEBUG
STRING		s1, s2, s3;
#endif
double		dValue;

    MESSAGE(( "Building internals for: %s\n", 
		sContainerName(cCont) ));

                /* Build torsion angles */
		/* These have to be built a special way, because */
		/* torsion INTERNALs cannot be simply looked up in a table */
		/* The model builder has code to add torsion INTERNALs */
		/* around central pairs of ATOMs taking into account any */
		/* external coordinates that may be defined already */

    LOOPOVERALL( cCont, BONDS, aTemp, ATOM, lTemp ) {
	LoopGetBond( &lTemp, &a2, &a3 );
	if ( iAtomCoordination(a2) <= 1 ) continue;
	if ( iAtomCoordination(a3) <= 1 ) continue;

		/* Assign INTERNALs for ALL torsions around (a2), (a3) */

	MESSAGE(( "Building torsion INTERNALs around: %s - %s\n",
			sContainerFullDescriptor( (CONTAINER)a2, s2 ),
			sContainerFullDescriptor( (CONTAINER)a3, s3 ) ));

	ModelAssignTorsionsAround( a2, a3, NULL );

    }
		/* Build bond angles */

    LOOPOVERALL( cCont, ANGLES, aTemp, ATOM, lTemp ) {
	LoopGetAngle( &lTemp, &a1, &a2, &a3 );
	MESSAGE(( "Building angle INTERNAL for: %s - %s - %s\n",
			sContainerFullDescriptor( (CONTAINER)a1, s1 ),
			sContainerFullDescriptor( (CONTAINER)a2, s2 ),
			sContainerFullDescriptor( (CONTAINER)a3, s3 ) ));
	if ( bAtomFlagsSet( a1, ATOMPOSITIONKNOWN ) &&
	     bAtomFlagsSet( a2, ATOMPOSITIONKNOWN ) &&
	     bAtomFlagsSet( a3, ATOMPOSITIONKNOWN ) ) {
	    dValue = dVectorAtomAngle( &vAtomPosition(a1),
					   &vAtomPosition(a2),
					   &vAtomPosition(a3) );
	    MESSAGE(( "Got bond angle from externals\n" ));
	} else {
	    dValue = dModelBondAngle( a1, a2, a3 );
	    MESSAGE(( "Got bond angle from model builder\n" ));
	}
	MESSAGE(( "++++Angle INTERNAL: %lf  for %s - %s - %s\n", 
		dValue/DEGTORAD,
		sContainerFullDescriptor((CONTAINER)a1,s1),
		sContainerFullDescriptor((CONTAINER)a2,s2),
		sContainerFullDescriptor((CONTAINER)a3,s3) ));
	iInternalAngle( a1, a2, a3, dValue );
    }

		/* Build bond internals */

    LOOPOVERALL( cCont, BONDS, aTemp, ATOM, lTemp ) {
	LoopGetBond( &lTemp, &a1, &a2 );
	MESSAGE(( "Building bond INTERNAL for: %s - %s\n",
			sContainerFullDescriptor( (CONTAINER)a1, s1 ),
			sContainerFullDescriptor( (CONTAINER)a2, s2 ) ));
        if ( bAtomFlagsSet( a1, ATOMPOSITIONKNOWN ) &&
            bAtomFlagsSet( a2, ATOMPOSITIONKNOWN ) ) {
	    dValue = dVectorAtomLength( &vAtomPosition(a1),
					    &vAtomPosition(a2) );
	    MESSAGE(( "Got bond length from externals\n" ));
	} else {
            dValue = dModelBondLength( a1, a2 );
	    MESSAGE(( "Got bond length from the model builder\n" ));
	}
        iInternalBond( a1, a2, dValue );
	MESSAGE(( "++++Bond INTERNAL: %lf  for %s - %s\n", 
		dValue,
		sContainerFullDescriptor((CONTAINER)a1,s1),
		sContainerFullDescriptor((CONTAINER)a2,s2) ));
    }

                /* Build chirality */
                /* Chirality is only defined for atoms with 3 or 4 bonds */
                /* Chirality is the only INTERNAL that does not */
                /* need to be specified */

    
    LOOPOVERALL( cCont, ATOMS, aTemp, ATOM, lTemp ) {
	dValue = dChiralityForAtom(aTemp);
	AtomSetFlags( aTemp, fSet );
	AtomResetFlags( aTemp, fReset );
	if ( dValue != 0.0 ) {
	    iInternalChirality( aTemp, dValue );
	    MESSAGE(( "Got chirality from external coordinates\n" ));
	    MESSAGE(( "++++Chirality INTERNAL: %lf  for %s\n", 
		dValue,
		sContainerFullDescriptor((CONTAINER)aTemp,s1) ));
	} else {
	    MESSAGE(( "Left chirality undefined\n" ));
	}
    }
}








/*
 *      BuildDestroyInternals
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy all of the INTERNALs on the ATOMs in the (lPLoop).
 *	(lPLoop) must be a LOOP over ATOMs.
 *
 */
void
BuildDestroyInternals( LOOP *lPLoop )
{
LOOP            lInternals;
INTERNAL        iInt;
ATOM            aAtom;
#ifdef	DEBUG
STRING		sTemp;
#endif

    MESSAGE(( "Destroying all INTERNALs in the LOOP\n" ));
    while ( ( aAtom=(ATOM)oNext(lPLoop) ) != NULL ) {
	MESSAGE(( "Destroying INTERNALs on: %s\n",
			sContainerFullDescriptor((CONTAINER)aAtom,sTemp) ));
        lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
#ifdef DEBUG
	while ( ( iInt=(INTERNAL)oNext(&lInternals) ) != NULL ) {
	    MESSAGE(( "    Internal type = %c  at: %lX\n", 
			iInternalType(iInt), iInt ));
	}
#endif
	lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
        while ( ( iInt=(INTERNAL)oNext(&lInternals) ) != NULL ) {
	    MESSAGE(( "    destroying    = %c  at: %lX\n", 
			iInternalType(iInt), iInt ));
            Destroy( (OBJEKT *)&iInt );
        }
    }
}





/*
 *      BuildFixInternals
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	When external coordinates are imposed on a structure
 *	that has internal coordinates defined, there will be
 *	inconsistancies in the torsions between the internals 
 *	and externals. This routine attempts to resolve these
 *	inconsistancies by changing the internal torsions to
 *	match the externals.
 *
 *      Loop through all bonds, looking for bonds that lie
 *      at the center of torsions.  Look at all of the torsions that
 *      share the bond as a center and calculate the actual torsion
 *      angle for each torsion.  If any of the calculated torsion angles
 *      is different than the value of the INTERNAL then choose one and
 *      calculate the difference between the INTERNAL and the calculated
 *      and add it to all of the INTERNAL torsions.
 */
void
BuildFixInternals( UNIT uUnit )
{
INTERNAL        iaTorsions[MAXTORSIONSAROUNDBOND], iInt;
int             iNextTorsion;
double          dTorsion;
BOOL            bFoundOne;
LOOP            lBonds;
ATOM            aAtom2, aAtom3;
int             i;
double          dSub, dNew;
#ifdef DEBUG
STRING		sAtom1, sAtom2, sAtom3, sAtom4;
#endif

    MESSAGE(( "In BuildFixInternals -----------\n" ));

                /* Loop over all bonds */
                /* looking for torsions */
                
    lBonds = lLoop( (OBJEKT)uUnit, BONDS );
    while ( oNext(&lBonds) != NULL ) {
        LoopGetBond( &lBonds, &aAtom2, &aAtom3 );
	MESSAGE(( "Looking at torsions around: %s - %s\n",
		sContainerFullDescriptor((CONTAINER)aAtom2,sAtom2),
		sContainerFullDescriptor((CONTAINER)aAtom3,sAtom3) ));
        
                /* Look at all the INTERNALs of aAtom2 */
                /* looking for torsions which have aAtom2 and */
                /* aAtom3 as the center of the torsion */

	iNextTorsion = iInternalFindAllTorsionInternalsAround( 
				aAtom2, aAtom3, iaTorsions );

                /* Look at these torsions and if all of the coordinates */
                /* are defined then calculate the angle.  Compare the */
                /* calculated with the IMPROPER and if they are different */
                /* calculate the difference */
                
        bFoundOne = FALSE;
        for ( i=0; i<iNextTorsion; i++ ) {
            iInt = iaTorsions[i];
	    MESSAGE(( "Measuring torsion of fixed atoms: %s - %s - %s - %s\n",
		sContainerFullDescriptor((CONTAINER)aInternalAtom1(iInt),sAtom1),
		sContainerFullDescriptor((CONTAINER)aInternalAtom2(iInt),sAtom2),
		sContainerFullDescriptor((CONTAINER)aInternalAtom3(iInt),sAtom3),
		sContainerFullDescriptor((CONTAINER)aInternalAtom4(iInt),sAtom4)));

            if ( !bAtomFlagsSet(aInternalAtom1(iInt),ATOMPOSITIONKNOWN) ) 
                continue;
            if ( !bAtomFlagsSet(aInternalAtom2(iInt),ATOMPOSITIONKNOWN) ) 
                continue;
            if ( !bAtomFlagsSet(aInternalAtom3(iInt),ATOMPOSITIONKNOWN) ) 
                continue;
            if ( !bAtomFlagsSet(aInternalAtom4(iInt),ATOMPOSITIONKNOWN) ) 
                continue;
            dTorsion = dVectorAtomTorsion( 
                &vAtomPosition(aInternalAtom1(iInt)),
                &vAtomPosition(aInternalAtom2(iInt)),
                &vAtomPosition(aInternalAtom3(iInt)),
                &vAtomPosition(aInternalAtom4(iInt)) );
            if ( dTorsion != dInternalValue(iInt) ) {
                dSub = dTorsion - dInternalValue(iInt);
                bFoundOne = TRUE;
            }
        }

                /* If a difference was found the modify ALL of the INTERNALS */

        if ( bFoundOne ) {
            MESSAGE(( "Twisting torsions centered on %s - %s by %lf degrees\n",
                        sContainerFullDescriptor((CONTAINER)aAtom2,sAtom2), 
                        sContainerFullDescriptor((CONTAINER)aAtom3,sAtom3), 
			dSub/DEGTORAD ));
            for ( i=0; i<iNextTorsion; i++ ) {
                dNew = dInternalValue(iaTorsions[i]) + dSub;
	        MESSAGE(( "Twisting torsion for atoms: %s-%s-%s-%s\n",
			sContainerName(aInternalAtom1(iaTorsions[i])),
			sContainerName(aInternalAtom2(iaTorsions[i])),
			sContainerName(aInternalAtom3(iaTorsions[i])),
			sContainerName(aInternalAtom4(iaTorsions[i])) ));
                MESSAGE(( "------- From %lf to %lf\n", 
                         dInternalValue(iaTorsions[i])/DEGTORAD, 
			 dNew/DEGTORAD ));
		
                InternalSetValue( iaTorsions[i], dNew );
            }
        }
    }
}







/*
 *      BuildHierarchy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	
 *      This function takes a CONTAINER and builds a hierarchy of:
 *      UNIT - MOLECULE - RESIDUE - ATOM under it.
 *      It will create RESIDUEs to contain ATOMs, and MOLECULEs to
 *      contain RESIDUEs.
 *
 *      Groupings of ATOMs that could be considered RESIDUEs are found
 *      by touching all ATOMs that are already in RESIDUEs, and then
 *      constructing spanning trees over all atoms that are not touched.
 *      The atoms within each spanning tree are then considered to be
 *      within separate RESIDUEs.
 *
 *      Groupings of ATOMs that could be considered MOLECULEs are found
 *      by touching all ATOMs that are already in MOLECULEs, and then
 *      constructing spanning trees over all atoms that are not touched.
 *      The atoms within each spanning tree are then considered to be
 *      within separate MOLECULEs, and the RESIDUEs that contain them
 *      are added to the MOLECULEs.  While the spanning trees are
 *      constructed, look for atoms which have been touched and store
 *      the MOLECULE that contains them, this is the MOLECULE that should
 *      be used to contain the other atoms, rather than creating a new
 *      one.
 */
void
BuildHierarchy( UNIT uUnit )
{
LOOP            lResidues;
CONTAINER       cCont;
RESIDUE         rRes;
LOOP            lAtoms, lSpan, lMolecules;
ATOM            aAtom, aAtom2, aTemp;
MOLECULE        mMol;

    ContainerResetAllAtomsFlags( (CONTAINER)uUnit, ATOMTOUCHED );

                /* Touch all ATOMs within RESIDUEs */

    lResidues = lLoop( (OBJEKT)uUnit, RESIDUES );
    while ( ( rRes = (RESIDUE)oNext(&lResidues) ) != NULL ) {
        lAtoms = lLoop( (OBJEKT)rRes, ATOMS );
        while ( ( aAtom=(ATOM)oNext(&lAtoms)) != NULL ) 
            AtomSetFlags( aAtom, ATOMTOUCHED );
    }

                /* Build spanning trees over the untouched ATOMs */

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    LoopDefineInvisibleAtoms( &lAtoms, ATOMTOUCHED );
    while ( (aAtom=(ATOM)oNext(&lAtoms)) != NULL ) {
        lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
        LoopDefineInvisibleAtoms( &lSpan, ATOMTOUCHED );

                /* Create a new RESIDUE for the atoms */

        rRes = (RESIDUE)oCreate(RESIDUEid);
        ContainerSetName( (CONTAINER)rRes, "???" );
        
                /* Remove all the the spanning tree ATOMs from the UNIT */
                /* and put them in the RESIDUE */

        cCont = (CONTAINER)uUnit;
        while ( (aAtom2=(ATOM)oNext(&lAtoms)) != NULL ) {
           bContainerRemove( (CONTAINER)uUnit, (OBJEKT)aAtom2 );
           ContainerAdd( (CONTAINER)rRes, (OBJEKT)aAtom2 );
           AtomSetFlags( aAtom2, ATOMTOUCHED );
           if ( iObjectType(cContainerWithin((CONTAINER)aAtom2)) == MOLECULEid )
                cCont = cContainerWithin((CONTAINER)aAtom2);
        }
        
                /* Now add the RESIDUE to the UNIT, or the MOLECULE that */
                /* contained one of the ATOMs */

                /* Figure out a way to set up the connections for the */
                /* RESIDUE */
                
        ContainerAdd( cCont, (OBJEKT)rRes );
    }

                /* untouch all ATOMs */
                
    ContainerResetAllAtomsFlags( (CONTAINER)uUnit, ATOMTOUCHED );

                /* Touch all ATOMs within MOLECULEs */

    lMolecules = lLoop( (OBJEKT)uUnit, MOLECULES );
    while ( ( mMol = (MOLECULE)oNext(&lMolecules) ) != NULL ) {
        lAtoms = lLoop( (OBJEKT)mMol, ATOMS );
        while ( ( aAtom=(ATOM)oNext(&lAtoms)) != NULL ) 
            AtomSetFlags( aAtom, ATOMTOUCHED );
    }

                /* Build spanning trees over the untouched ATOMs */

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    LoopDefineInvisibleAtoms( &lAtoms, ATOMTOUCHED );
    while ( (aAtom=(ATOM)oNext(&lAtoms)) != NULL ) {

                /* First run through the spanning tree looking for */
                /* collisions */

        lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
        LoopDefineInvisibleAtoms( &lSpan, ATOMTOUCHED );
        while ( oNext(&lSpan) != NULL ) /* Nothing */;
        if ( iLoopInvisibleCollisionCount(&lSpan) > 0 ) {
            aTemp = aLoopLastCollisionAtom(&lSpan);
            rRes = (RESIDUE)cContainerWithin(aTemp);
            mMol = (MOLECULE)cContainerWithin(rRes);
        } else {
                        /* Create a new MOLECULE */
            mMol = (MOLECULE)oCreate(MOLECULEid);
            ContainerSetName( (CONTAINER)mMol, "???" );
            ContainerAdd( (CONTAINER)uUnit, (OBJEKT)mMol );
        }

                /* Search for a collision with an invisible atom */
                /* which will be contained within a MOLECULE to which */
                /* all of the residues should be added */
                
                /* Remove all the the spanning tree ATOMs from the UNIT */
                /* and put them in the RESIDUE */

        lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
        LoopDefineInvisibleAtoms( &lSpan, ATOMTOUCHED );
        while ( (aAtom2=(ATOM)oNext(&lSpan)) != NULL ) {
           rRes = (RESIDUE)cContainerWithin((CONTAINER)aAtom2);
           if ( iObjectType(cContainerWithin((CONTAINER)rRes))==UNITid ) {
                bContainerRemove( (CONTAINER)uUnit, (OBJEKT)rRes );
                ContainerAdd( (CONTAINER)mMol, (OBJEKT)rRes );
           }
           AtomSetFlags( aAtom2, ATOMTOUCHED );
        }
    }
}







/*
 *	bBuildChangeInternalBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Change the bond INTERNAL for the specified atoms within
 *	the CONTAINER.
 *	This routine assumes that the atoms within the CONTAINER
 *	have all of their INTERNAL coordinates.
 *	This function searches through the CONTAINER for the
 *	atom with name sName, and then looks for the other atoms
 *	that are bonded in a chain.  If it doesn't find it, it
 *	returns FALSE, otherwise it finds the INTERNAL that
 *	contains all of the ATOMs and changes the value of
 *	the INTERNAL to what the caller specified.
 */
BOOL
bBuildChangeInternalBond( CONTAINER cCont, char *sAtom1, char *sAtom2, 
			double dValue )
{
ATOM		aAtom1, aAtom2;
INTERNAL	iInt;

		/* Find the first named ATOM within the container */

    aAtom1 = (ATOM)cContainerFindName( cCont, ATOMid, sAtom1 );
    if ( aAtom1 == NULL ) return(FALSE);

		/* Find the bonded ATOM that has the required name */

    aAtom2 = zaBuildFindNextAtomWithNameButNot( aAtom1, sAtom2, NULL );
    if ( aAtom2 == NULL ) return(FALSE);

		/* Find the INTERNAL for the ATOMs */

    iInt = iInternalFindBond( aAtom1, aAtom2 );

		/* Change it's value */

    InternalSetValue( iInt, dValue );

    return(TRUE);
}




/*
 *	bBuildChangeInternalAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Change the angle INTERNAL for the specified atoms within
 *	the CONTAINER.
 *	The value dValue is in RADIANS.
 *	This routine assumes that the atoms within the CONTAINER
 *	have all of their INTERNAL coordinates.
 *	This function searches through the CONTAINER for the
 *	atom with name sName, and then looks for the other atoms
 *	that are bonded in a chain.  If it doesn't find it, it
 *	returns FALSE, otherwise it finds the INTERNAL that
 *	contains all of the ATOMs and changes the value of
 *	the INTERNAL to what the caller specified.
 */
BOOL
bBuildChangeInternalAngle( CONTAINER cCont, 
	char *sAtom1, char *sAtom2, char *sAtom3, double dValue )
{
ATOM		aAtom1, aAtom2, aAtom3;
INTERNAL	iInt;

		/* Find the first named ATOM within the container */

    aAtom1 = (ATOM)cContainerFindName( cCont, ATOMid, sAtom1 );
    if ( aAtom1 == NULL ) return(FALSE);

		/* Find the bonded ATOM that has the required name */

    aAtom2 = zaBuildFindNextAtomWithNameButNot( aAtom1, sAtom2, NULL );
    if ( aAtom2 == NULL ) return(FALSE);
    aAtom3 = zaBuildFindNextAtomWithNameButNot( aAtom2, sAtom3, aAtom1 );
    if ( aAtom3 == NULL ) return(FALSE);

		/* Find the INTERNAL for the ATOMs */

    iInt = iInternalFindAngle( aAtom1, aAtom2, aAtom3 );

		/* Change it's value */

    InternalSetValue( iInt, dValue );

    return(TRUE);
}





/*
 *	bBuildChangeInternalTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Change the torsion INTERNAL for the specified atoms within
 *	the CONTAINER.
 *	The value dValue must be in RADIANS.
 *	This routine assumes that the atoms within the CONTAINER
 *	have all of their INTERNAL coordinates.
 *	This function searches through the CONTAINER for the
 *	atom with name sName, and then looks for the other atoms
 *	that are bonded in a chain.  If it doesn't find it, it
 *	returns FALSE, otherwise it finds the INTERNAL that
 *	contains all of the ATOMs and twists all torsions
 *	around the central pair of atoms to make the torsion
 *	the caller requested the correct value.
 */
BOOL
bBuildChangeInternalTorsion( CONTAINER cCont, char *sAtom1, char *sAtom2,
				char *sAtom3, char *sAtom4, double dValue )
{
ATOM		aAtom1, aAtom2, aAtom3, aAtom4;
INTERNAL	iInt, iaTorsions[MAXTORSIONSAROUNDBOND];
double		dAdd, dTorsion;
int		i, iTorsions;
#ifdef DEBUG
STRING		s1, s2, s3, s4;
#endif

		/* Find the first named ATOM within the container */

    aAtom1 = (ATOM)cContainerFindName( cCont, ATOMid, sAtom1 );
    if ( aAtom1 == NULL ) return(FALSE);

		/* Find the bonded ATOM that has the required name */

    aAtom2 = zaBuildFindNextAtomWithNameButNot( aAtom1, sAtom2, NULL );
    if ( aAtom2 == NULL ) return(FALSE);
    aAtom3 = zaBuildFindNextAtomWithNameButNot( aAtom2, sAtom3, aAtom1 );
    if ( aAtom3 == NULL ) return(FALSE);
    aAtom4 = zaBuildFindNextAtomWithNameButNot( aAtom3, sAtom4, aAtom2 );
    if ( aAtom4 == NULL ) return(FALSE);

		/* Find the INTERNAL for the ATOMs */

    iInt = iInternalFindTorsion( aAtom1, aAtom2, aAtom3, aAtom4 );
    MESSAGE(( "Currently the internal between: %s-%s-%s-%s is: %lf\n",
		sContainerFullDescriptor((CONTAINER)aAtom1,s1),
		sContainerFullDescriptor((CONTAINER)aAtom2,s2),
		sContainerFullDescriptor((CONTAINER)aAtom3,s3),
		sContainerFullDescriptor((CONTAINER)aAtom4,s4),
		dInternalValue(iInt)/DEGTORAD ));

		/* Find the amount we have to twist all the other torsions */

    dAdd = dValue - dInternalValue( iInt );
    MESSAGE(( "Going to twist all torsions around: %s - %s by %lf\n",
		sContainerName(aAtom2), sContainerName(aAtom3),
		dAdd/DEGTORAD ));

		/* Get all of the torsions around the center pair of atoms */

    iTorsions = iInternalFindAllTorsionInternalsAround( aAtom2, aAtom3, 
						iaTorsions );

		/* Change all of the torsions */

    for ( i=0; i<iTorsions; i++ ) {
	dTorsion = dAdd + dInternalValue( iaTorsions[i] );
	InternalSetValue( iaTorsions[i], dTorsion );
    }

    return(TRUE);
}





/*
 *	BuildRelaxInFramework
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Relax the ATOMs in the UNIT that have ATOMNEEDSMINIMIZER
 *	set.  Relax the ATOMs in the framework of the fixed ATOMs
 *	around them.
 *	Use the MINIMIZER mStrain.
 */
void
BuildRelaxInFramework( UNIT uUnit, MINIMIZER mStrain )
{
LOOP            lAtoms, lTemp;
ATOM            aAtom, aAtom1, aAtom2, aAtom3, aAtom4;
BOOL            bM1, bM2, bM3, bM4, bOneMinimizedAtom;
double          dKb, dR0, dKt, dT0, dTkub, dRkub, dKp, dP0;
double          dScEE, dScNB;
STRING		sAtom1, sAtom2, sAtom3, sAtom4, sDesc;
PARMSET		psTemp;
TORSION		tTorsion;
int		iTag, i, iIndex, iN, iDefaults;

                /* Now do the minimization */
                /* Add all atoms that */
                /* require the minimizer to them */
                /* Then add all the bond, angle, torsion interactions */
                /* that involve only atoms that need the minimizer and */
                /* all atoms whose coordinates are fixed that are involved */
                /* in the interactions. */
                /* Only the atoms whose coordinates are fixed need to be */
                /* added to the MINIMIZER object as the bonds, angles, */
                /* torsions are being looped through */ 

                /* Loop over all atoms, adding those that require */
                /* minimization to the MINIMIZER object */

    bOneMinimizedAtom = FALSE;
    MESSAGE(( "^^^Looping over atoms to add to MINIMIZER\n" ));
    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    LoopDefineVisibleAtoms( &lAtoms, ATOMNEEDSMINIMIZER );
    while ( ( aAtom=(ATOM)oNext(&lAtoms) ) != NULL ) {
        MinimizerAddAtom( mStrain, aAtom );
        bOneMinimizedAtom = TRUE;
    }

                /* If there are no atoms to minimize then return */

    if ( !bOneMinimizedAtom ) 
	return;

                /* Loop over all bonds, adding those that contain */
                /* at least one atom that needs the minimizer and */
                /* otherwise only atoms with fixed coordinates */
                /* Add atoms that have fixed coordinates to the */
                /* MINIMIZER object first */

    MESSAGE(( "^^^Looping over bonds to add to MINIMIZER\n" ));
    iDefaults = 0;
    lTemp = lLoop( (OBJEKT)uUnit, BONDS );
    while ( oNext(&lTemp) != NULL ) {
        LoopGetBond( &lTemp, &aAtom1, &aAtom2 );

        bM1 = bAtomFlagsSet( aAtom1, ATOMNEEDSMINIMIZER );
        bM2 = bAtomFlagsSet( aAtom2, ATOMNEEDSMINIMIZER );

                /* If neither of the atoms needs minimizing then continue */

        if ( ! (bM1 || bM2) ) 
		continue;

                /* Add atoms that may be missing to the minimizer */
        if ( !bM1 ) 
		MinimizerAddAtom( mStrain, aAtom1 );
        if ( !bM2 ) 
		MinimizerAddAtom( mStrain, aAtom2 );

		/* Try to find the bond parameter in the PARMLIB */
	iTag = PARM_NOT_FOUND;	
	PARMLIB_DEFAULT_LOOP( psTemp, 
		( iTag = iParmSetFindBond( psTemp,
					      sAtomType(aAtom1),
					      sAtomType(aAtom2) ) ) );
	if ( iTag != PARM_NOT_FOUND ) {
	    ParmSetBond( psTemp, iTag, sAtom1, sAtom2, &dKb, &dR0, sDesc );
	} else {

		/* Get a parameter from the model builder */

	    	ModelBondParm( aAtom1, aAtom2, &dKb, &dR0 );
		iDefaults++;
	}
        if ( !bMinimizerAddBond( mStrain, aAtom1, aAtom2, dKb, dR0 ) ) {
                DFATAL(( "Could not add bond to MINIMIZER" ));
        }
    }
    if ( iDefaults )
	VP0(( " (used %d default bond params)\n", iDefaults ));


                /* Loop over all angles, adding those that contain */
                /* at least one atom that needs the minimizer and */
                /* otherwise only atoms with fixed coordinates */
                /* Add atoms that have fixed coordinates to the */
                /* MINIMIZER object first */

    MESSAGE(( "^^^Looping over angles to add to MINIMIZER\n" ));
    iDefaults = 0;
    lTemp = lLoop( (OBJEKT)uUnit, ANGLES );
    while ( oNext(&lTemp) != NULL ) {
        LoopGetAngle( &lTemp, &aAtom1, &aAtom2, &aAtom3 );
	bM1 = bAtomFlagsSet( aAtom1, ATOMNEEDSMINIMIZER );
	bM2 = bAtomFlagsSet( aAtom2, ATOMNEEDSMINIMIZER );
	bM3 = bAtomFlagsSet( aAtom3, ATOMNEEDSMINIMIZER );

                /* If none of the atoms need minimization then continue */

        if ( !(bM1 || bM2 || bM3) ) 
		continue;

                /* Add any atoms that may not be in the Minimizer */

        if ( !bM1 ) 
		MinimizerAddAtom( mStrain, aAtom1 );
        if ( !bM2 ) 
		MinimizerAddAtom( mStrain, aAtom2 );
        if ( !bM3 ) 
		MinimizerAddAtom( mStrain, aAtom3 );
	
		/* Try to find the angle parameter in the PARMLIB */
		
	iTag = PARM_NOT_FOUND;	
	PARMLIB_DEFAULT_LOOP( psTemp, 
		( iTag = iParmSetFindAngle( psTemp,
					      sAtomType(aAtom1),
					      sAtomType(aAtom2),
					      sAtomType(aAtom3) ) ) );
	if ( iTag != PARM_NOT_FOUND ) {
	    ParmSetAngle( psTemp, iTag, 
	    			sAtom1, sAtom2, sAtom3,
				&dKt, &dT0, &dTkub, &dRkub, sDesc );
	} else {

		/* Get a parameter from the model builder */
		iDefaults++;
	    	ModelAngleParm( aAtom1, aAtom2, aAtom3, &dKt, &dT0 );
	}
	if (!bMinimizerAddAngle( mStrain, aAtom1, aAtom2, aAtom3, dKt, dT0 )) {
                DFATAL(( "Could not add angle to MINIMIZER" ));
        }
    }

    if ( iDefaults )
	VP0(( " (used %d default angle params)\n", iDefaults ));


                /* Loop over all torsions, adding those that contain */
                /* at least one atom that needs the minimizer and */
                /* otherwise only atoms with fixed coordinates */
                /* Add atoms that have fixed coordinates to the */
                /* MINIMIZER object first */

    MESSAGE(( "^^^Looping over torsions to add to MINIMIZER\n" ));
    iDefaults = 0;
    lTemp = lLoop( (OBJEKT)uUnit, PROPERS );
    while ( oNext(&lTemp) != NULL ) {
        LoopGetTorsion( &lTemp, &aAtom1, &aAtom2, &aAtom3, &aAtom4 );
	bM1 = bAtomFlagsSet( aAtom1, ATOMNEEDSMINIMIZER );
	bM2 = bAtomFlagsSet( aAtom2, ATOMNEEDSMINIMIZER );
	bM3 = bAtomFlagsSet( aAtom3, ATOMNEEDSMINIMIZER );
	bM4 = bAtomFlagsSet( aAtom4, ATOMNEEDSMINIMIZER );

                /* If none of the atoms need minimization then continue */

        if ( !(bM1 || bM2 || bM3 || bM4) ) 
		continue;

                /* Add any atoms that may not be in the Minimizer */

        if ( !bM1 ) 
		MinimizerAddAtom( mStrain, aAtom1 );
        if ( !bM2 ) 
		MinimizerAddAtom( mStrain, aAtom2 );
        if ( !bM3 ) 
		MinimizerAddAtom( mStrain, aAtom3 );
        if ( !bM4 ) 
		MinimizerAddAtom( mStrain, aAtom4 );

		/* Search for all of the torsion terms */

	tTorsion = tParmSetTORSIONCreate();
	PARMLIB_DEFAULT_LOOP_ALL( psTemp ) {
	    iParmSetFindProperTerms( psTemp, tTorsion, FALSE,
	    				sAtomType(aAtom1),
					sAtomType(aAtom2),
					sAtomType(aAtom3),
					sAtomType(aAtom4) );
	}
	if ( iParmSetTORSIONTermCount(tTorsion) == 0 ) {
		iDefaults++;
		ModelProperTerms( aAtom1, aAtom2, aAtom3, aAtom4, tTorsion );
	}

	for ( i=0; i<iParmSetTORSIONTermCount(tTorsion); i++ ) {
	    ParmSetTORSIONTerm( tTorsion, i,
	    		&iIndex,
			sAtom1, sAtom2, sAtom3, sAtom4,
			&iN, &dKp, &dP0, &dScEE, &dScNB,
			sDesc );
	    if ( !bMinimizerAddTorsion( mStrain, 
	    				aAtom1, aAtom2, aAtom3, aAtom4,
					(double)iN, dKp, dP0 )) {
		DFATAL(( "Could not add proper to MINIMIZER" ));
	    }
	}
        ParmSetTORSIONDestroy( &tTorsion );
    }

    if ( iDefaults )
	VP0(( " (used %d default torsion params)\n", iDefaults ));

    MinimizerMinimize( mStrain );
}





/*
 *	BuildRotateAroundBondFromTo
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Rotate ATOMs around a bond.  The bond is defined from
 *	aInv to aStart.  All of the ATOMs on the aStart side of
 *	the bond are rotated by the angle dRotate.  Special
 *	things are done if the variable bInRing is TRUE.
 */
void
BuildRotateAroundBondFromTo( CONTAINER cCont, ATOM aInv, ATOM aStart, 
		double dRotate, BOOL bInRing )
{
MATRIX		mTransform;
LOOP		lSpan;
ATOM		aAtom;
VECTOR		vPos, vNew;
int		j;

	/* Build the transformation matrix to rotate around the bond */

    MatrixRotateAround( mTransform, &vAtomPosition(aStart),
			&vAtomPosition(aInv), dRotate );

    MESSAGE(( "Rotating torsion around %s -> %s   inRing: %s\n",
		sAtomName(aInv), sAtomName(aStart), sBOOL(bInRing) ));

	/* If the bond to be rotated around is not within a ring */
	/* then just rotate, otherwise do something different */

    if ( !bInRing ) {
	lSpan = lLoop( (OBJEKT)aStart, SPANNINGTREE );
	LoopDefineInvisibleAtom( &lSpan, aInv );
	while ( (aAtom = (ATOM)oNext(&lSpan )) ) {
	    vPos = vAtomPosition(aAtom);
	    MatrixTimesVector( vNew, mTransform, vPos );
	    AtomSetPosition( aAtom, vNew );
	}
    } else {
	for ( j=0; j<iAtomCoordination(aStart); j++ ) {
	    if ( aAtomBondedNeighbor(aStart,j) != aInv ) {
		vPos = vAtomPosition(aAtomBondedNeighbor(aStart,j));
		MatrixTimesVector( vNew, mTransform, vPos );
		AtomSetPosition( aAtomBondedNeighbor(aStart,j), vNew );
	    }
	}
    }
}	    






/*
 *	bBuildFlipChiralityFor
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Flip the chirality for one ATOM.
 *	This is done by prioritizing the side chains coming
 *	off of the ATOM and moving the two chains that have
 *	the least number of ATOMs.  The moving is done by
 *	building the internal coordinates for the side chains and
 *	then swapping the values of the torsions that contain the
 *	first ATOMs in the moving side chains and the ATOM whos
 *	chirality is being flipped and the ATOMs bonded to the
 *	ATOM being flipped which are not moving.
 *	Once the torsions are flipped, the externals are build
 *	again.
 *
 *	Return TRUE if the chirality was flipped.
 *
 *TODO:	Make this routine handle rings properly.
 *
 */

typedef	struct	{
	ATOM	aAtom;
	int	iCount;
} COUNTATOMSt;

void
bBuildFlipChiralityFor( CONTAINER cContainer, ATOM aFlip )
{
COUNTATOMSt	caaCount[MAXBONDS];
ATOM		aAtom;
ATOM		aA, aB, aFixA, aFixB;
LOOP		lSpan, lInternals, lAtoms;
int		i, iCountRings, iDum;
INTERNAL	iRing;
INTERNAL	iFlipRing;
VECTOR		vNew, vA, vB;
VECTOR		vNormal;
MATRIX		mTransform;
BOOL		bPartOfRing;

	/* The ATOM must have at least 3 bonds and no more than 4 */

    if ( iAtomCoordination(aFlip) < 3 ||
	 iAtomCoordination(aFlip) > 4 ) return;

    MESSAGE(( "The coordination of the ATOM to flip is: %d\n",
			iAtomCoordination(aFlip) ));
    MESSAGE(( "About to flip: %s\n", sAtomName(aFlip) ));

	/* Find the two side chains with the least number of ATOMs */

    for ( i=0; i<iAtomCoordination(aFlip); i++ ) {
	aAtom = aAtomBondedNeighbor(aFlip,i);
	lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
	LoopDefineInvisibleAtom( &lSpan, aFlip );
	caaCount[i].aAtom = aAtom;
	caaCount[i].iCount = 0;
	while ( (aAtom = (ATOM)oNext(&lSpan)) ) caaCount[i].iCount++;
    }

	/* Sort the bonds by the number of ATOMs */

    SortByInteger( (GENP) caaCount, iAtomCoordination(aFlip),
                   sizeof(COUNTATOMSt), (GENP) &(caaCount[0].iCount), TRUE );

	/* Figure out which ATOMs are to be moved. */
	/* This depends on the number of atoms around the central atom */
	/* and whether they are in rings */

    GraphUtilFindAllSmallestRings( (UNIT) cContainer );

	/* Count the number of rings passing through the flip ATOM */
	/* If there is more than 1 then we cannot perform the flip */

    lInternals = lLoop( (OBJEKT)aFlip, INTERNALS );
    iCountRings = 0;
    iFlipRing = NULL;
    while ( (iRing = (INTERNAL)oNext(&lInternals)) ) {
	if ( iInternalType(iRing) == INTERNALRING ) {
	    iCountRings++;
	    iFlipRing = iRing;
	}
    }
    if ( iCountRings > 1 ) {
	lAtoms = lLoop( (OBJEKT)cContainer, ATOMS );
	BuildDestroyInternals( &lAtoms );
	MESSAGE(( "Cannot flip!  There is more than one ring\n" ));
	return;
    }

	/* If there is a ring then find the ATOMs that are */
	/* bonded but not in the ring and put them in aA, aB */
	/* And put the ATOMs that are in the ring in aFixA, aFixB */

    if ( iCountRings == 1 ) {
	MESSAGE(( "Atom to flip is in a ring\n" ));
	aA = NULL;
	aB = NULL;
	aFixA = NULL;
	aFixB = NULL;
	for ( i=0; i<iAtomCoordination(aFlip); i++ ) {
	    aAtom = aAtomBondedNeighbor( aFlip, i );
	    lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
	    bPartOfRing = FALSE;
	    while ( (iRing = (INTERNAL)oNext(&lInternals)) ) {
		if ( iInternalType(iRing) == INTERNALRING &&
		     iRing == iFlipRing ) {
		    if ( aFixA == NULL ) aFixA = aAtom;
		    else 		 aFixB = aAtom;
		    bPartOfRing = TRUE;
		} else break;
	    }
	    if ( !bPartOfRing ) {
		if ( aA == NULL )	aA = aAtom;
		else			aB = aAtom;
	    }
	}
    } else {

		/* If there are no rings then prioritize the ATOMs around */
		/* aFlip by the size of their side chains */

	if ( iAtomCoordination(aFlip) == 3 ) {
	    aA = caaCount[0].aAtom;
	    aB = NULL;
	    aFixA = caaCount[1].aAtom;
	    aFixB = caaCount[2].aAtom;
	} else {
	    aA = caaCount[0].aAtom;
	    aB = caaCount[1].aAtom;
	    aFixA = caaCount[2].aAtom;
	    aFixB = caaCount[3].aAtom;
	}
    }
    MESSAGE(( "Fixed atoms: %s and %s\n", 
		sAtomName(aFixA), sAtomName(aFixB) ));
    MESSAGE(( "Moving atom: %s\n", sAtomName(aA) ));
    if ( aB != NULL ) {
	MESSAGE(( "Moving atom: %s\n", sAtomName(aB) ));
    }

		/* Now destroy all of the ring INTERNALs */

    lAtoms = lLoop( (OBJEKT)cContainer, ATOMS );
    BuildDestroyInternals( &lAtoms );

		/* Setup the flags that will be required to properly */
		/* build external coordinates for only the two ATOMs */
		/* that will be moving and their side chains */

    ContainerResetAllAtomsFlags( cContainer, ATOMNEEDSBUILD );

    lSpan = lLoop( (OBJEKT)aA, SPANNINGTREE );
    LoopDefineInvisibleAtom( &lSpan, aFlip );
    while ( (aAtom = (ATOM)oNext(&lSpan)) ) {
	AtomSetFlags( aAtom, ATOMNEEDSBUILD );
    }
    if ( aB != NULL ) {
	lSpan = lLoop( (OBJEKT)aB, SPANNINGTREE );
	LoopDefineInvisibleAtom( &lSpan, aFlip );
	while ( (aAtom = (ATOM)oNext(&lSpan)) ) {
	    AtomSetFlags( aAtom, ATOMNEEDSBUILD );
	}
    }

		/* Build the INTERNAL coordinates */

    lSpan = lLoop( (OBJEKT)aA, SPANNINGTREE );
    LoopDefineInvisibleAtom( &lSpan, aFlip );
    BuildInternalsUsingFlags( &lSpan, 
				ATOMNEEDSBUILD, 
				0,
				0, 
				ATOMPOSITIONKNOWN );
    if ( aB != NULL ) {
	lSpan = lLoop( (OBJEKT)aB, SPANNINGTREE );
	LoopDefineInvisibleAtom( &lSpan, aFlip );
	BuildInternalsUsingFlags( &lSpan,
					ATOMNEEDSBUILD,
					0,
					0,
					ATOMPOSITIONKNOWN );
    }	

		/* Now build a MATRIX which will reflect the ATOMs */
		/* that are to be moved through the plane defined */
		/* by aFlip, aFixA, aFixB.  Only aA and aB(if != NULL) */
		/* will be flipped.  The side chains connected to */
		/* aA and aB will be rebuilt using their INTERNALs */

    vA = vVectorSub( &vAtomPosition(aFixA), &vAtomPosition(aFlip) );
    vB = vVectorSub( &vAtomPosition(aFixB), &vAtomPosition(aFlip) );
    vNormal = vVectorCross( &vA, &vB );

    MatrixReflectAcross( mTransform, &vAtomPosition(aFlip), &vNormal );

		/* Flip the ATOMs aA and aB(if != NULL) across the plane */

    MatrixTimesVector( vNew, mTransform, vAtomPosition(aA) );
    AtomSetPosition( aA, vNew );
    AtomResetFlags( aA, ATOMNEEDSBUILD );
    if ( aB != NULL ) {
	MatrixTimesVector( vNew, mTransform, vAtomPosition(aB) );
	AtomSetPosition( aB, vNew );
	AtomResetFlags( aB, ATOMNEEDSBUILD );
    }

		/* Build the rest of the coordinates for the side chains */

    lSpan = lLoop( (OBJEKT)aA, SPANNINGTREE );
    LoopDefineInvisibleAtom( &lSpan, aFlip );

			/* Skip the first ATOM because we flipped it */
			/* using the MATRIX */
    oNext(&lSpan);
    iDum = 0; /* for purify */
    BuildExternalsUsingFlags( &lSpan,
				ATOMNEEDSBUILD,
				0,
				ATOMPOSITIONKNOWN,
				ATOMNEEDSBUILD, 
				&iDum, &iDum, &iDum, TRUE );

    if ( aB != NULL ) {
	lSpan = lLoop( (OBJEKT)aB, SPANNINGTREE );
	LoopDefineInvisibleAtom( &lSpan, aFlip );

			/* Skip the first ATOM because we flipped it */
			/* using the MATRIX */
	oNext(&lSpan);
	BuildExternalsUsingFlags( &lSpan,
				ATOMNEEDSBUILD,
				0,
				ATOMPOSITIONKNOWN,
				ATOMNEEDSBUILD,
				&iDum, &iDum, &iDum, TRUE );
    }

		/* Destroy all of the INTERNALs */

    lSpan = lLoop( (OBJEKT)aA, SPANNINGTREE );
    LoopDefineInvisibleAtom( &lSpan, aFlip );
    BuildDestroyInternals( &lSpan );
    if ( aB != NULL ) {
	lSpan = lLoop( (OBJEKT)aB, SPANNINGTREE );
	LoopDefineInvisibleAtom( &lSpan, aFlip );
	BuildDestroyInternals( &lSpan );
    }

		/* Reset all of the flags that say positions fixed */

    ContainerResetAllAtomsFlags( cContainer, ATOMPOSITIONFIXED );
}





/*
 *	BuildInternalsBetweenUnitsUsingFlags
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build INTERNAL coordinates for the linkage between
 *	the two UNITs that are about to be sequenced.
 *	This is done by building a small spanning tree from the 
 *	CONNECT ATOMs between the UNITs.
 *
 *	Only build INTERNALs for ATOMs that have (fForSet) set
 *	and (fForReset) reset.
 */
void
BuildInternalsBetweenUnitsUsingFlags( UNIT uFirst, UNIT uSecond,
					FLAGS fForSet, FLAGS fForReset )
{
ATOM		aLast, aFirst;
LOOP		lSpan;

    MESSAGE(( "&   BuildInternalsBetweenUnitsUsingFlags\n" ));

    if ( uFirst != NULL && uSecond != NULL &&
	 bUnitHeadUsed(uSecond) &&
	 bUnitTailUsed(uFirst) ) {
	aFirst = aUnitHead(uSecond);
	aLast = aUnitTail(uFirst);

		/* Create a bond between two ATOMs in different */
		/* UNITs.  This is only temporary so that we can */
		/* build a spanning tree across a junction */

	if ( AtomTmpBondTo( aFirst, aLast ) == FALSE ) {
		VP0(("Skipping generating internals - effect unknown\n"));
		return;
	}

		/* Now create a spanning tree that is only 4 ATOMs */
		/* deep to create INTERNALs across the junction */
		/* Using 4 atoms will garentee that INTERNALs will */
		/* be built for ATOMs up to 3 ATOMs away from (aLast) */

	lSpan = lLoop( (OBJEKT)aFirst, SPANNINGTREE );
	LoopDefineMaxDistanceFromRoot( &lSpan, 4 );
	BuildInternalsUsingFlags( &lSpan, fForSet, fForReset, 0, 0 );

		/* Break the temporary bond */

	AtomRemoveBond( aFirst, aLast );

    }
}








/*
 *	BuildInternalsForSimpleRings
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	This facility allows the builder to handle simple
 *	six membered rings which will probably account for
 *	90% of all rings that are drawn within LEaP.
 *
 *	Find all of the smallest rings in the CONTAINER,
 *	and then search through the ring groups for simple rings:
 *
 *	BENZENE
 *	-	six membered rings where at most one atom is
 *			specified, and all ATOMs are SP2.
 *	CYCLOHEXANE
 *	-	six membered rings where at most one atom is
 *			specified, and all ATOMs are SP3.
 *
 *	For these rings, assign internal coordinates for
 *	the ring only.
 */
void
BuildInternalsForSimpleRings( CONTAINER cContainer )
{
VARARRAY	vaRingGroups;
int		i;
LIST		lRingGroup;

    GraphUtilFindAllSmallestRingsAndRingGroups( (UNIT) cContainer, &vaRingGroups );

		/* search through all of the RINGS */

    for ( i=0; i<iVarArrayElementCount(vaRingGroups); i++ ) {
	lRingGroup = *PVAI( vaRingGroups, LIST, i );
	zbBuildTorsionInternalsForRingGroupMaybe( lRingGroup );
    }

		/* Destroy the ring groups VARARRAY */

    GraphUtilDestroyRingGroupVarArray(&vaRingGroups);
}
