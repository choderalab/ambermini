/*
 *	File:	restraint.h
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
 *		RESTRAINTs are objects that maintain
 *		a bond/angle/torsion term that the user
 *		wants to add to a UNIT.
 *		These RESTRAINTs are added to the energy
 *		expression for the UNIT.
 *		RESTRAINTs all share the same data structure
 *		even though there are different types.
 *		The types are RESTRAINTBOND, RESTRAINTANGLE,
 *		and RESTRAINTTORSION.
 *		RESTRAINTs maintain a set of atoms across
 *		which the restraint is defined, and 
 *		several constants which define the energy
 *		function that makes up the RESTRAINT.
 *
 *		RESTRAINTs also have flags that determine when
 *		in a Free Energy Perturbation calculation they
 *		are active.  By DEFAULT they are active for
 *		both lambda=0 and lambda=1 (nonpert, pert)
 *
 */

#ifndef	RESTRAINT_H
#define RESTRAINT_H

#define	RESTRAINTNONE		0
#define	RESTRAINTBOND		1
#define	RESTRAINTANGLE		2
#define	RESTRAINTTORSION	3

		/* Flags that tell when the RESTRAINT is active */

#define	RESTRAINTNEVER		0x00000000
#define	RESTRAINTNONPERT	0x00000001
#define	RESTRAINTPERT		0x00000002
#define	RESTRAINTALWAYS		(RESTRAINTNONPERT | RESTRAINTPERT)


typedef	struct	{
	ATOM		aAtom1;
	ATOM		aAtom2;
	double		dKr;
	double		dR0;
} RESTRAINTBONDt;

typedef	struct	{
	ATOM		aAtom1;
	ATOM		aAtom2;
	ATOM		aAtom3;
	ATOM		aAtom4;
	double		dKt;
	double		dT0;
} RESTRAINTANGLEt;

typedef	struct	{
	ATOM		aAtom1;
	ATOM		aAtom2;
	ATOM		aAtom3;
	ATOM		aAtom4;
	double		dKp;
	double		dP0;
	double		dN;
} RESTRAINTTORSIONt;

typedef	struct	{
	int		iType;
	FLAGS		fFlags;
	union	{
		RESTRAINTBONDt		rbBond;
		RESTRAINTANGLEt		raAngle;
		RESTRAINTTORSIONt	rtTorsion;
	} rType;
} RESTRAINTt;

typedef	RESTRAINTt	*RESTRAINT;



extern RESTRAINT	rRestraintCreate();
extern void		RestraintDestroy(RESTRAINT *rPRes);

extern BOOL	bRestraintBondMatchAtoms(RESTRAINT rRes, 
			ATOM aAtom1, ATOM aAtom2);
extern BOOL	bRestraintAngleMatchAtoms(RESTRAINT rRes,
			ATOM aAtom1, ATOM aAtom2, ATOM aAtom3);
extern BOOL	bRestraintTorsionMatchAtoms(RESTRAINT rRes,
			ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, ATOM aAtom4);

extern BOOL	bRestraintContainsAtom(RESTRAINT rRes, ATOM aAtom);


#define	iRestraintType(r)		( (r)->iType )
#define	bRestraintFlagsSet(r,f)		( ((r)->fFlags | (f))!= 0 )
#define	RestraintDefineFlags(r,f)	( (r)->fFlags = (f) )
#define	fRestraintFlags(r)		( (r)->fFlags )

#define	RestraintBondSet(r,a1,a2,kr,r0) ( \
	(r)->iType = RESTRAINTBOND, \
	(r)->fFlags = (RESTRAINTALWAYS), \
	(r)->rType.rbBond.aAtom1 = a1, \
	(r)->rType.rbBond.aAtom2 = a2, \
	(r)->rType.rbBond.dKr = kr, \
	(r)->rType.rbBond.dR0 = r0 )

#define	RestraintAngleSet(r,a1,a2,a3,kt,t0) ( \
	(r)->iType = RESTRAINTANGLE, \
	(r)->fFlags = (RESTRAINTALWAYS), \
	(r)->rType.raAngle.aAtom1 = a1, \
	(r)->rType.raAngle.aAtom2 = a2, \
	(r)->rType.raAngle.aAtom3 = a3, \
	(r)->rType.raAngle.dKt = kt, \
	(r)->rType.raAngle.dT0 = t0 )

#define	RestraintTorsionSet(r,a1,a2,a3,a4,kp,p0,n) ( \
	(r)->iType = RESTRAINTTORSION, \
	(r)->fFlags = (RESTRAINTALWAYS), \
	(r)->rType.rtTorsion.aAtom1 = a1, \
	(r)->rType.rtTorsion.aAtom2 = a2, \
	(r)->rType.rtTorsion.aAtom3 = a3, \
	(r)->rType.rtTorsion.aAtom4 = a4, \
	(r)->rType.rtTorsion.dKp = kp, \
	(r)->rType.rtTorsion.dN = n, \
	(r)->rType.rtTorsion.dP0 = p0 )

#define	RestraintBondGet(r,aP1,aP2,dPkr,dPr0) ( \
	*aP1  = (r)->rType.rbBond.aAtom1, \
	*aP2  = (r)->rType.rbBond.aAtom2, \
	*dPkr = (r)->rType.rbBond.dKr, \
	*dPr0 = (r)->rType.rbBond.dR0 )

#define	RestraintAngleGet(r,aP1,aP2,aP3,dPkt,dPt0) ( \
	*aP1  = (r)->rType.raAngle.aAtom1, \
	*aP2  = (r)->rType.raAngle.aAtom2, \
	*aP3  = (r)->rType.raAngle.aAtom3, \
	*dPkt = (r)->rType.raAngle.dKt, \
	*dPt0 = (r)->rType.raAngle.dT0 )

#define	RestraintTorsionGet(r,aP1,aP2,aP3,aP4,dPkp,dPp0,dPn) ( \
	*aP1  = (r)->rType.rtTorsion.aAtom1, \
	*aP2  = (r)->rType.rtTorsion.aAtom2, \
	*aP3  = (r)->rType.rtTorsion.aAtom3, \
	*aP4  = (r)->rType.rtTorsion.aAtom4, \
	*dPkp = (r)->rType.rtTorsion.dKp, \
	*dPp0 = (r)->rType.rtTorsion.dP0, \
	*dPn  = (r)->rType.rtTorsion.dN )


#endif	/* RESTRAINT_H */





