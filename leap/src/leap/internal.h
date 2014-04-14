#ifndef INTERNAL_H
#define INTERNAL_H
/*
 *      File: internal.h
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
 *              An INTERNAL is an object representing an internal coordinate.
 *              Either a chirality, a bond length, a bond angle, 
 *              a torsion angle, or an improper torsion angle
 *              (not implemented yet).
 *              INTERNALs are contained within ATOMs.  Each atom
 *              that is referenced by the INTERNAL also contains the
 *              INTERNAL.  ATOMs contain internals as part of their
 *              CONTAINER contents.
 *		A special type of INTERNAL is the RING.  This allows
 *		the system to calculate and keep track of where RINGs
 *		are.
 *		
 *
 *
 */
 

#include	"ring.h"


/*
 *-----------------------------------------------------------------------
 *
 *      Define object typedefs here.
 *        
 *      Object typedef MUST include the superclass object as
 *      its first structure element.
*/


#define INTERNALUNDEFINED       '?'
#define INTERNALCHIRALITY       'C'
#define INTERNALBOND            'B'
#define INTERNALANGLE           'A'
#define INTERNALTORSION         'T'
#define INTERNALIMPROPER        'I'
#define	INTERNALRING		'R'

typedef	struct	{
	ATOM		aAtom1, aAtom2, aAtom3, aAtom4;
	double		dValue;
} NORMALINTERNALt;


/*
 *	New type used to store ring INTERNALs
 */


typedef	struct	{
	RING		rAtoms;		/* Store RING atoms here */
	RINGLOOP	rlLoop;		/* Used to loop over ATOMs */
	int		iRingTempInt;	/* Used to separate rings */
} RINGINTERNALt;

typedef struct  {
	OBJEKTt         oSuper;
	char            cInternalType;
	union	{
		NORMALINTERNALt		niNormal;
		RINGINTERNALt		riRing;
	} iType;
} INTERNALt;
                
typedef INTERNALt	*INTERNAL;


#define         MAXTORSIONSAROUNDBOND   12      /* Actually never more than 9 */




/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */

extern INTERNAL		iInternalCreate();
extern void		InternalDestroy(INTERNAL *iPInt);
extern void		InternalDescribe(INTERNAL iInt);

extern INTERNAL		iInternalChirality(ATOM aAtom1, double dValue);
extern INTERNAL		iInternalBond(ATOM aAtom1, ATOM aAtom2, double dValue);
extern INTERNAL		iInternalAngle(ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, 
				double dValue);
extern INTERNAL		iInternalTorsion(ATOM aAtom1, ATOM aAtom2, 
				ATOM aAtom3, ATOM aAtom4, double dValue);
extern INTERNAL		iInternalRing();

extern void		InternalRingAddAtomAfter(INTERNAL iRing, ATOM aAtom, 
				ATOM aPrev);
extern void		InternalRingAddAtomBefore(INTERNAL iInt, ATOM aAtom, 
				ATOM aBefore);
#define	iInternalRingSize(II) iRingSize(((INTERNAL)(II))->iType.riRing.rAtoms)
extern void		InternalRingLoopAtoms(INTERNAL iInt);
extern ATOM		aInternalRingNextAtom(INTERNAL iInt);
extern BOOL		bInternalRingRemoveAtom(INTERNAL iInt, ATOM aAtom);
#define	InternalRingSetTempInt(i,in)	\
		((INTERNAL)(i)->iType.riRing.iTempInt = in )
#define	iInternalRingTempInt(i)	((INTERNAL)(i)->iType.riRing.iTempInt)


extern INTERNAL		iInternalFindChirality(ATOM aAtom1);
extern INTERNAL		iInternalFindTorsion(ATOM aAtom1, ATOM aAtom2, 
				ATOM aAtom3, ATOM aAtom4);
extern INTERNAL		iInternalFindAngle(ATOM aAtom1, ATOM aAtom2, 
				ATOM aAtom3);
extern INTERNAL		iInternalFindBond(ATOM aAtom1, ATOM aAtom2);

extern BOOL		bInternalGoodTorsion(INTERNAL iInt, ATOM aAtom, 
				INTERNAL *iPTorsion, INTERNAL *iPAngle, 
				INTERNAL *iPBond);
extern BOOL		bInternalGoodAngle(INTERNAL iInt, ATOM aAtom, 
				INTERNAL *iPAngle, INTERNAL *iPBond);
extern BOOL		bInternalGoodPairOfAngles( INTERNAL iInt, ATOM aAtom,
				INTERNAL *iPAngle1, INTERNAL *iPAngle2, 
				INTERNAL *iPBond); 
extern BOOL		bInternalGoodBond(INTERNAL iInt, ATOM aAtom, 
				INTERNAL *iPBond);


#define InternalSetType(II,TT)  (((INTERNAL)II)->cInternalType=TT)
#define iInternalType(II)       ((int)(((INTERNAL)II)->cInternalType))
#define aInternalAtom1(II)      (((INTERNAL)(II))->iType.niNormal.aAtom1)
#define aInternalAtom2(II)      (((INTERNAL)(II))->iType.niNormal.aAtom2)
#define aInternalAtom3(II)      (((INTERNAL)(II))->iType.niNormal.aAtom3)
#define aInternalAtom4(II)      (((INTERNAL)(II))->iType.niNormal.aAtom4)
#define dInternalValue(II)      (((INTERNAL)(II))->iType.niNormal.dValue)
#define InternalSetValue(II,VV) (((INTERNAL)(II))->iType.niNormal.dValue=VV)


extern int	iInternalFindAllTorsionInternalsAround(ATOM aAtom2, 
			ATOM aAtom3, INTERNAL iaTorsions[]);

#endif/* INTERNAL_H */
