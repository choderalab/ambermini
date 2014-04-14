/*
 *      File: internal.c
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
 *              INTERNAL
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              An INTERNAL is an internal coordinate, containing
 *              pointers to the atoms that define the internal coordinate
 *              and the value of the internal coordinate.
 *              INTERNALs are contained within EVERY! ATOM that is
 *              referenced by the internal coordinate.
 *
 */




#include	"basics.h"

#include        "classes.h"

#ifndef	RING_H
#include	"ring.h"
#endif



/*
===================================================================

        Define private routines here.
*/


/*
*******************************************************************

        Define public message routines here.
*/



/*
 *      iInternalCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an INTERNAL, initialize it's contents.
 *      This is not the routine normally used to create an internal
 *      coordinate because it does not allow definition of the
 *      type of internal coordinate or the atoms/value of the
 *      INTERNAL.
 *      It is better to use iInternalBond, iInternalAngle etc.
 */
INTERNAL
iInternalCreate()
{
INTERNAL        iNew;

    MALLOC( iNew, INTERNAL, sizeof(INTERNALt) );
    
    iNew->cInternalType = INTERNALUNDEFINED;
    iNew->iType.niNormal.aAtom1 = NULL;
    iNew->iType.niNormal.aAtom2 = NULL;
    iNew->iType.niNormal.aAtom3 = NULL;
    iNew->iType.niNormal.aAtom4 = NULL;
    
    return(iNew);
}




/*
 *      InternalDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy an INTERNAL, first removing the INTERNAL from
 *      ALL of the atoms that it references.
 */
void    
InternalDestroy( INTERNAL *iPInt )
{
ATOM		aAtom;
#ifdef	DEBUG
STRING		s1, s2, s3, s4;
#endif

        /* Remove the IMPROPER from ALL containers that contain it */
        

    REF( *iPInt );  /* bContainerRemove() needs this */

    switch ( iInternalType(*iPInt) ) {
        case INTERNALIMPROPER:
        case INTERNALTORSION:
	    MESSAGE(( "Removing torsion INTERNAL: %s-%s-%s-%s\n",
		sContainerFullDescriptor((CONTAINER)(*iPInt)->iType.niNormal.aAtom1,s1),
		sContainerFullDescriptor((CONTAINER)(*iPInt)->iType.niNormal.aAtom2,s2),
		sContainerFullDescriptor((CONTAINER)(*iPInt)->iType.niNormal.aAtom3,s3),
		sContainerFullDescriptor((CONTAINER)(*iPInt)->iType.niNormal.aAtom4,s4) ));
            bContainerRemove( (CONTAINER)(*iPInt)->iType.niNormal.aAtom4, (OBJEKT)*iPInt );
	    bContainerRemove( (CONTAINER)(*iPInt)->iType.niNormal.aAtom3, (OBJEKT)*iPInt );
	    bContainerRemove( (CONTAINER)(*iPInt)->iType.niNormal.aAtom2, (OBJEKT)*iPInt );
	    if ( (*iPInt)->iType.niNormal.aAtom4 != 
			(*iPInt)->iType.niNormal.aAtom1 ) {
		bContainerRemove( (CONTAINER)(*iPInt)->iType.niNormal.aAtom1, (OBJEKT)*iPInt );
	    }
	    break;
        case INTERNALANGLE:
            bContainerRemove( (CONTAINER)(*iPInt)->iType.niNormal.aAtom3, (OBJEKT)*iPInt );
        case INTERNALBOND:
            bContainerRemove( (CONTAINER)(*iPInt)->iType.niNormal.aAtom2, (OBJEKT)*iPInt );
        case INTERNALCHIRALITY:
            bContainerRemove( (CONTAINER)(*iPInt)->iType.niNormal.aAtom1, (OBJEKT)*iPInt );
            break;
	case INTERNALRING:
	    InternalRingLoopAtoms(*iPInt);
	    while ( (aAtom = aInternalRingNextAtom(*iPInt)) ) {
		bContainerRemove( (CONTAINER)aAtom, (OBJEKT)*iPInt );
	    }
	    RingDestroy( &((*iPInt)->iType.riRing.rAtoms) );
	    break;
        default:
	    DFATAL(( "Illegal type of INTERNAL" ));
            break;
    }

        /* Now Free the INTERNAL */

    FREE( *iPInt );
    *iPInt = NULL;
}







/*
 *      InternalDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe the INTERNAL.
 *
 */
void
InternalDescribe( INTERNAL iInt )
{
    switch ( iInternalType(iInt) ) {
        case INTERNALCHIRALITY:
            VP0(( "Internal chirality.\n" ));
            VP0(( "Atom1=%s\n", sContainerName(iInt->iType.niNormal.aAtom1) ));
            break;
        case INTERNALBOND:
            VP0(( "Internal bond length.\n" ));
            VP0(( "Atom1=%s\n", sContainerName(iInt->iType.niNormal.aAtom1) ));
            VP0(( "Atom2=%s\n", sContainerName(iInt->iType.niNormal.aAtom2) ));
            break;
        case INTERNALANGLE:
            VP0(( "Internal bond angle.\n" ));
            VP0(( "Atom1=%s\n", sContainerName(iInt->iType.niNormal.aAtom1) ));
            VP0(( "Atom2=%s\n", sContainerName(iInt->iType.niNormal.aAtom2) ));
            VP0(( "Atom3=%s\n", sContainerName(iInt->iType.niNormal.aAtom3) ));
            break;
        case INTERNALTORSION:
            VP0(( "Internal torsion angle.\n" ));
            VP0(( "Atom1=%s\n", sContainerName(iInt->iType.niNormal.aAtom1) ));
            VP0(( "Atom2=%s\n", sContainerName(iInt->iType.niNormal.aAtom2) ));
            VP0(( "Atom3=%s\n", sContainerName(iInt->iType.niNormal.aAtom3) ));
            VP0(( "Atom4=%s\n", sContainerName(iInt->iType.niNormal.aAtom4) ));
            break;
        case INTERNALIMPROPER:
            VP0(( "Internal improper angle.\n" ));
            VP0(( "Atom1=%s\n", sContainerName(iInt->iType.niNormal.aAtom1) ));
            VP0(( "Atom2=%s\n", sContainerName(iInt->iType.niNormal.aAtom2) ));
            VP0(( "Atom3=%s\n", sContainerName(iInt->iType.niNormal.aAtom3) ));
            VP0(( "Atom4=%s\n", sContainerName(iInt->iType.niNormal.aAtom4) ));
            break;
        case INTERNALUNDEFINED:
            VP0(( "Internal undefined coordinate.\n" ));
            break;
	case INTERNALRING:
	    VP0(( "Ring with %d atoms\n", iInternalRingSize(iInt) ));
	    break;
        default:
            DFATAL(( "Unknown INTERNAL" ));
            break;
    }
    VP0(( "Internal value=%lf\n\n", iInt->iType.niNormal.dValue ));
}



/*
 *      iInternalChirality
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an INTERNAL and initialize it as a chirality.
 *      Then ADD the internal coordinate to all the atoms that it
 *      references.
 */
INTERNAL
iInternalChirality( ATOM aAtom1, double dValue )
{
INTERNAL        iInt;

    TESTMEMORY();
    iInt = (INTERNAL)oCreate(INTERNALid);
    iInt->cInternalType = INTERNALCHIRALITY;
    iInt->iType.niNormal.aAtom1 = aAtom1;
    iInt->iType.niNormal.dValue = dValue;
  
    TESTMEMORY(); 
    ContainerAdd( (CONTAINER)aAtom1, (OBJEKT)iInt );
    TESTMEMORY();
    return(iInt);
}





/*
 *      iInternalBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an INTERNAL and initialize it as a bond length.
 *      Then ADD the internal coordinate to all the atoms that it
 *      references.
 */
INTERNAL
iInternalBond( ATOM aAtom1, ATOM aAtom2, double dValue )
{
INTERNAL        iInt;

    iInt = (INTERNAL)oCreate(INTERNALid);
    iInt->cInternalType = INTERNALBOND;
    iInt->iType.niNormal.aAtom1 = aAtom1;
    iInt->iType.niNormal.aAtom2 = aAtom2;
    iInt->iType.niNormal.dValue = dValue;
    
    ContainerAdd( (CONTAINER)aAtom1, (OBJEKT)iInt );
    ContainerAdd( (CONTAINER)aAtom2, (OBJEKT)iInt );
    return(iInt);
}




/*
 *      iInternalAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an INTERNAL and initialize it as a bond angle.
 *      Then ADD the internal coordinate to all the atoms that it
 *      references.
 */
INTERNAL
iInternalAngle( ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, double dValue )
{
INTERNAL        iInt;

    iInt = (INTERNAL)oCreate(INTERNALid);
    iInt->cInternalType = INTERNALANGLE;
    iInt->iType.niNormal.aAtom1 = aAtom1;
    iInt->iType.niNormal.aAtom2 = aAtom2;
    iInt->iType.niNormal.aAtom3 = aAtom3;
    iInt->iType.niNormal.dValue = dValue;
   
    ContainerAdd( (CONTAINER)aAtom1, (OBJEKT)iInt );
    ContainerAdd( (CONTAINER)aAtom2, (OBJEKT)iInt );
    ContainerAdd( (CONTAINER)aAtom3, (OBJEKT)iInt );
    return(iInt);
}




/*
 *      iInternalTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an INTERNAL and initialize it as a torsion angle.
 *      Then ADD the internal coordinate to all the atoms that it
 *      references.
 */
INTERNAL
iInternalTorsion( ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, ATOM aAtom4, 
		double dValue )
{
INTERNAL        iInt;

    iInt = (INTERNAL)oCreate(INTERNALid);
    iInt->cInternalType = INTERNALTORSION;
    iInt->iType.niNormal.aAtom1 = aAtom1;
    iInt->iType.niNormal.aAtom2 = aAtom2;
    iInt->iType.niNormal.aAtom3 = aAtom3;
    iInt->iType.niNormal.aAtom4 = aAtom4;
    iInt->iType.niNormal.dValue = dValue;

    ContainerAdd( (CONTAINER)aAtom1, (OBJEKT)iInt );
    ContainerAdd( (CONTAINER)aAtom2, (OBJEKT)iInt );
    ContainerAdd( (CONTAINER)aAtom3, (OBJEKT)iInt );

		/* If the four atoms are not in a 3 member ring then */
		/* add the torsion to the fourth atom */

    if ( aAtom4 != aAtom1 ) {
    	ContainerAdd( (CONTAINER)aAtom4, (OBJEKT)iInt );
    }
    return(iInt);
}





/*
 *      iInternalFindChirality
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Find the INTERNAL chirality that contains the atom.
 */
INTERNAL
iInternalFindChirality( ATOM aAtom1 )
{
LOOP            lInternals;
INTERNAL        iInt;

    lInternals = lLoop( (OBJEKT)aAtom1, INTERNALS );
    while ( ( iInt = (INTERNAL)oNext(&lInternals) ) != NULL ) {
        if ( iInternalType(iInt) == INTERNALCHIRALITY ) return(iInt);
    }
    return(NULL);
}






/*
 *      iInternalFindTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Find the INTERNAL torsion that contains the four atoms:
 *      aAtom1, aAtom2, aAtom3, aAtom4;
 */
INTERNAL        
iInternalFindTorsion( ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, ATOM aAtom4 )
{
LOOP            lInternals;
INTERNAL        iInt;


    MESSAGE(( "Looking for torsion %s-%s-%s-%s\n",
		sContainerName(aAtom1),
		sContainerName(aAtom2),
		sContainerName(aAtom3),
		sContainerName(aAtom4) ));

    lInternals = lLoop( (OBJEKT)aAtom1, INTERNALS );
    while ( ( iInt = (INTERNAL)oNext(&lInternals) ) != NULL ) {
        if ( iInternalType(iInt) != INTERNALTORSION ) continue;
	MESSAGE(( "Looking at: %s-%s-%s-%s\n",
		sContainerName(aInternalAtom1(iInt)),
		sContainerName(aInternalAtom2(iInt)),
		sContainerName(aInternalAtom3(iInt)),
		sContainerName(aInternalAtom4(iInt)) ));
        if ( aAtom1 == aInternalAtom1(iInt) ) {
            if ( (aAtom2 == aInternalAtom2(iInt)) &&
                 (aAtom3 == aInternalAtom3(iInt)) &&
                 (aAtom4 == aInternalAtom4(iInt)) ) return(iInt);
        } else if ( aAtom1 == aInternalAtom4(iInt) ) {
            if ( (aAtom2 == aInternalAtom3(iInt)) &&
                 (aAtom3 == aInternalAtom2(iInt)) &&
                 (aAtom4 == aInternalAtom1(iInt)) ) return(iInt);
	}	    
    }
    MESSAGE(( "Found nothing\n" ));

    return(NULL);
}





/*
 *      iInternalFindAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Find the INTERNAL angle that contains the three atoms:
 *      aAtom1, aAtom2, aAtom3;
 */
INTERNAL
iInternalFindAngle( ATOM aAtom1, ATOM aAtom2, ATOM aAtom3 )
{
LOOP            lInternals;
INTERNAL        iInt;

    lInternals = lLoop( (OBJEKT)aAtom2, INTERNALS );
    while ( ( iInt = (INTERNAL)oNext(&lInternals) ) != NULL ) {
        if ( iInternalType(iInt) != INTERNALANGLE ) continue;
        if ( aAtom2 == aInternalAtom2(iInt) ) {
            if ( (aAtom1 == aInternalAtom1(iInt)) &&
                 (aAtom3 == aInternalAtom3(iInt)) ) return(iInt);
            if ( (aAtom1 == aInternalAtom3(iInt)) &&
                 (aAtom3 == aInternalAtom1(iInt)) ) return(iInt);
        }
    }
    return(NULL);
}






/*
 *      iInternalFindBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Find the INTERNAL bond that contains the two atoms:
 *      aAtom1, aAtom2;
 */
INTERNAL
iInternalFindBond( ATOM aAtom1, ATOM aAtom2 )
{
LOOP            lInternals;
INTERNAL        iInt;

    lInternals = lLoop( (OBJEKT)aAtom1, INTERNALS );
    while ( ( iInt = (INTERNAL)oNext(&lInternals) ) != NULL ) {
        if ( iInternalType(iInt) != INTERNALBOND ) continue;
        if ( aAtom2 == aInternalAtom2(iInt) ) return(iInt);
        if ( aAtom2 == aInternalAtom1(iInt) ) return(iInt);
    }
    return(NULL);
}







/*
 *      bInternalGoodTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the INTERNAL specified can be considered
 *      to be a good torsion.
 *      Goodness is measured in terms of whether or not the
 *      torsion can be used to construct the external coordinates
 *      of aAtom.
 *      The torsion is good if:
 *      1) aAtom is a terminal atom of the torsion.
 *      2) The other three atoms of the torsion have their coordinates
 *              specified.
 *      3) There exist an angle and a bond that lie on the torsion,
 *              and have aAtom as their terminal atoms.
 *
 *      If the torsion is good then return the bond, angle and torsion
 *      in iPBond, iPAngle, and iPTorsion.
 */
BOOL
bInternalGoodTorsion( INTERNAL iInt, ATOM aAtom, INTERNAL *iPTorsion, 
			INTERNAL *iPAngle, INTERNAL *iPBond )
{
BOOL		bFound;

    bFound = FALSE;
    if ( (aAtom == aInternalAtom1(iInt)) ) {
        if ( bAtomFlagsSet( aInternalAtom2(iInt), ATOMPOSITIONKNOWN ) &&
             bAtomFlagsSet( aInternalAtom3(iInt), ATOMPOSITIONKNOWN ) &&
             bAtomFlagsSet( aInternalAtom4(iInt), ATOMPOSITIONKNOWN ) ) {
             *iPTorsion = iInt;
             *iPAngle = iInternalFindAngle( aAtom, aInternalAtom2(iInt),
                                                aInternalAtom3(iInt) );
             *iPBond = iInternalFindBond( aAtom, aInternalAtom2(iInt) );
             bFound = TRUE;
        }
    } else if ( (aAtom == aInternalAtom4(iInt)) ) {
        if ( bAtomFlagsSet( aInternalAtom1(iInt), ATOMPOSITIONKNOWN ) &&
             bAtomFlagsSet( aInternalAtom2(iInt), ATOMPOSITIONKNOWN ) &&
             bAtomFlagsSet( aInternalAtom3(iInt), ATOMPOSITIONKNOWN ) ) {
             *iPTorsion = iInt;
             *iPAngle = iInternalFindAngle( aAtom, aInternalAtom3(iInt),
                                                aInternalAtom2(iInt) );
             *iPBond = iInternalFindBond( aAtom, aInternalAtom3(iInt) );
             bFound = TRUE;
        }
    }

#ifdef	DEBUG
    if ( bFound ) {
	ATOM	a1, a2, a3, a4;
	if ( aAtom == aInternalAtom1(iInt) ) {
	    a1 = aInternalAtom1(iInt);
	    a2 = aInternalAtom2(iInt);
	    a3 = aInternalAtom3(iInt);
	    a4 = aInternalAtom4(iInt);
	} else {
	    a1 = aInternalAtom4(iInt);
	    a2 = aInternalAtom3(iInt);
	    a3 = aInternalAtom2(iInt);
	    a4 = aInternalAtom1(iInt);
	}
	if ( !*iPBond || !*iPAngle || !*iPTorsion ) {
	    STRING	s1, s2, s3, s4;
	    VP0(( "From: %s to: %s - %s - %s\n",
			sContainerFullDescriptor((CONTAINER) a1, s1 ),
			sContainerFullDescriptor((CONTAINER) a2, s2 ),
			sContainerFullDescriptor((CONTAINER) a3, s3 ),
			sContainerFullDescriptor((CONTAINER) a4, s4 ) ));
	    DFATAL(( "Bad internal coordinate: iBond = %lX, iAngle = %lX\n",
			*iPBond, *iPAngle ));
	}
    }
#endif

    return(bFound);
}



        

/*
 *      bInternalGoodPairOfAngles
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the INTERNAL specified can be considered
 *      to be a good angle and another angle can be found that shares
 *      the terminal atom and the middle atom, and all atoms have
 *      coordinates specified.
 *      Goodness is measured in terms of whether or not the
 *      angle can be used to construct the external coordinates
 *      of aAtom.
 *      The angle is good if:
 *      1) aAtom is a terminal atom of the angle.
 *      2) The other two atoms of the angle have their coordinates
 *              specified.
 *      3) There exists another angle that has the atom as a terminal
 *              and shares the middle atom with the first angle and
 *              has all of its coordinates specified.
 *      4) There exists a bond that lies on the angle,
 *              and has aAtom as its terminal atom.
 *
 *      If the angle is good then return the angles and bond
 *      in iPAngle1, iPAngle2, and iPBond.
 */
BOOL
bInternalGoodPairOfAngles( INTERNAL iInt, ATOM aAtom, 
		INTERNAL *iPAngle1, INTERNAL *iPAngle2, INTERNAL *iPBond )
{
BOOL            bGotOne;
ATOM            aAtom2, aAtom3;
LOOP            lInternals;
INTERNAL        iNew;

    aAtom2 = NULL;
    bGotOne = FALSE;
    if ( (aAtom == aInternalAtom1(iInt)) ) {
        if ( bAtomFlagsSet( aInternalAtom2(iInt), ATOMPOSITIONKNOWN ) &&
             bAtomFlagsSet( aInternalAtom3(iInt), ATOMPOSITIONKNOWN ) ) {
             *iPAngle1 = iInt;
             *iPBond = iInternalFindBond( aAtom, aInternalAtom2(iInt) );
             aAtom2 = aInternalAtom2(iInt);
             aAtom3 = aInternalAtom3(iInt);
             bGotOne = TRUE;
        }
    }
    if ( (aAtom == aInternalAtom3(iInt)) ) {
        if ( bAtomFlagsSet( aInternalAtom1(iInt), ATOMPOSITIONKNOWN ) &&
             bAtomFlagsSet( aInternalAtom2(iInt), ATOMPOSITIONKNOWN )) {
             *iPAngle1 = iInt;
             *iPBond = iInternalFindBond( aAtom, aInternalAtom2(iInt) );
             aAtom2 = aInternalAtom2(iInt);
             aAtom3 = aInternalAtom1(iInt);
             bGotOne = TRUE;
        }
    }
    if ( !bGotOne ) return(FALSE);

                /* Loop over the internals of the first atom */
    lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
    while ( (iNew = (INTERNAL)oNext(&lInternals)) ) {
                /* Make sure that we don't get the same angle twice */
        if ( iNew == iInt ) continue;
        if ( iInternalType(iNew) == INTERNALANGLE ) {
                /* Check if the angle shares the center atom with the other */
                /* angle */
            if ( aInternalAtom2(iNew)==aAtom2 ) {
                if ( aAtom == aInternalAtom1(iNew) ) {
                    if ( bAtomFlagsSet( aInternalAtom3(iNew), 
                                        ATOMPOSITIONKNOWN ) ) {
                        *iPAngle2 = iNew;
                        return(TRUE);
                    }
                } else {
                    if ( bAtomFlagsSet( aInternalAtom1(iNew),
                                        ATOMPOSITIONKNOWN ) ) {
                        *iPAngle2 = iNew;
                        return(TRUE);
                    }
                }
            }
        }
    }
    return(FALSE);
}




/*
 *      bInternalGoodAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the INTERNAL specified can be considered
 *      to be a good angle.
 *      Goodness is measured in terms of whether or not the
 *      angle can be used to construct the external coordinates
 *      of aAtom.
 *      The angle is good if:
 *      1) aAtom is a terminal atom of the angle.
 *      2) The other two atoms of the angle have their coordinates
 *              specified.
 *      3) There exists a bond that lies on the angle,
 *              and have aAtom as its terminal atom.
 *
 *      If the angle is good then return the bond
 *      in iPAngle, and iPBond.
 */
BOOL    
bInternalGoodAngle( INTERNAL iInt, ATOM aAtom, INTERNAL *iPAngle, 
			INTERNAL *iPBond )
{

    if ( (aAtom == aInternalAtom1(iInt)) ) {
        if ( bAtomFlagsSet( aInternalAtom2(iInt), ATOMPOSITIONKNOWN ) &&
             bAtomFlagsSet( aInternalAtom3(iInt), ATOMPOSITIONKNOWN ) ) {
             *iPAngle = iInt;
             *iPBond = iInternalFindBond( aAtom, aInternalAtom2(iInt) );
             return(TRUE);
        }
    }
    if ( (aAtom == aInternalAtom3(iInt)) ) {
        if ( bAtomFlagsSet( aInternalAtom1(iInt), ATOMPOSITIONKNOWN ) &&
             bAtomFlagsSet( aInternalAtom2(iInt), ATOMPOSITIONKNOWN )) {
             *iPAngle = iInt;
             *iPBond = iInternalFindBond( aAtom, aInternalAtom2(iInt) );
             return(TRUE);
        }
    }
    return(FALSE);
}



        


/*
 *      bInternalGoodBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the INTERNAL specified can be considered
 *      to be a good bond.
 *      Goodness is measured in terms of whether or not the
 *      bond can be used to construct the external coordinates
 *      of aAtom.
 *      The torsion is good if:
 *      1) aAtom is a terminal atom of the bond.
 *      2) The other atom of the bond has its coordinates
 *              specified.
 *
 *      If the bond is good then return the bond in iPBond.
 */
BOOL
bInternalGoodBond( INTERNAL iInt, ATOM aAtom, INTERNAL *iPBond )
{
    if ( (aAtom == aInternalAtom1(iInt)) ) {
        if ( bAtomFlagsSet( aInternalAtom2(iInt), ATOMPOSITIONKNOWN )) {
             *iPBond = iInt;
             return(TRUE);
        }
    }
    if ( (aAtom == aInternalAtom2(iInt)) ) {
        if ( bAtomFlagsSet( aInternalAtom1(iInt), ATOMPOSITIONKNOWN )) {
             *iPBond = iInt;
             return(TRUE);
        }
    }
    return(FALSE);
}




/*
 *	iInternalFindAllTorsionInternalsAround
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Find all of the torsion INTERNALs around a bond.
 *	Put them into the iaTorsions array.
 *	Return the number found.
 *
 *	NOTE: Care must be taken to make sure that the array iaTorsions
 *		is large enough to accomidate all of the torsions.
 *		It should be dimensioned for MAXTORSIONSAROUNDBOND INTERNALs.
 */
int
iInternalFindAllTorsionInternalsAround( ATOM aAtom2, ATOM aAtom3, 
		INTERNAL iaTorsions[] )
{
int		iTorsions;
LOOP		lInternals;
INTERNAL	iInt;

    iTorsions = 0;
    lInternals = lLoop( (OBJEKT)aAtom2, INTERNALS );
    while ( (iInt = (INTERNAL)oNext(&lInternals)) ) {
	if ( iInternalType(iInt) == INTERNALTORSION ) {
	    if ( (aInternalAtom2(iInt) == aAtom2 &&
		  aInternalAtom3(iInt) == aAtom3 ) ||
		 (aInternalAtom2(iInt) == aAtom3 &&
		  aInternalAtom3(iInt) == aAtom2 ) ) {
		iaTorsions[iTorsions++] = iInt;
	    }
	}
    }
    return(iTorsions);
}





/*
 *	iInternalRing
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return an INTERNAL that is a ring.
 */
INTERNAL
iInternalRing()
{
INTERNAL	iInt;

    iInt = (INTERNAL)oCreate(INTERNALid);
    InternalSetType( iInt, INTERNALRING );
    iInt->iType.riRing.rAtoms = rRingCreate();
    return(iInt);
}





/*
 *	InternalRingAddAtomAfter
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add an ATOM to the ring INTERNAL and add the ring INTERNAL
 *	to that ATOM.  Add the ATOM after the ATOM (aPrev).
 */
void
InternalRingAddAtomAfter( INTERNAL iRing, ATOM aAtom, ATOM aPrev )
{
    RingAfterAdd( iRing->iType.riRing.rAtoms, (GENP)aPrev, (GENP)aAtom );
    ContainerAdd( (CONTAINER)aAtom, (OBJEKT)iRing );
}



/*
 *	InternalRingAddAtomBefore
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add an ATOM to the ring INTERNAL and add the ring INTERNAL
 *	to that ATOM.  Add the ATOM before the ATOM (aBefore).
 */
void
InternalRingAddNextAdjacentAtom( INTERNAL iInt, ATOM aAtom, ATOM aBefore )
{
    RingBeforeAdd( iInt->iType.riRing.rAtoms, (GENP)aBefore, (GENP)aAtom );
    ContainerAdd( (CONTAINER)aAtom, (OBJEKT)iInt );
}



/*
 *	InternalRingLoopAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Initialize the loop over the ring ATOMs.
 */
void
InternalRingLoopAtoms( INTERNAL iInt )
{
    rlRingLoop(iInt->iType.riRing.rAtoms, &iInt->iType.riRing.rlLoop);
}



/*
 *	aInternalRingNextAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the next ATOM in the ring, if there are no more,
 *	return NULL.
 */
ATOM
aInternalRingNextAtom( INTERNAL iInt )
{
    return((ATOM)PRingNext( &(iInt->iType.riRing.rlLoop) ));
}



/*
 *	bInternalRingRemoveAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Remove the ATOM from the ring INTERNAL, return TRUE if it
 *	was found, otherwise FALSE.
 */
BOOL
bInternalRingRemoveAtom( INTERNAL iInt, ATOM aAtom )
{
    if ( bRingRemove( iInt->iType.riRing.rAtoms, (GENP)aAtom ) ) {
	REF( iInt );  /* bContainerRemove() needs this */
	if ( !bContainerRemove( (CONTAINER)aAtom, (OBJEKT)iInt ) ) {
	    DFATAL(( "The ATOM did not contain the INTERNAL" ));
	}
	DEREF( iInt );  /* reset after bContainerRemove() */
	return(TRUE);
    }
    return(FALSE);
}









        
