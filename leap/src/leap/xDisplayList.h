/*
 *	File:	xDisplayList.h
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
 *		The display list maintains a list of objects that
 *		can be displayed within the viewing window.
 *		The display list is a linked list of records,
 *		each record contains an operation and 
 *		data for the operation.  Operations
 *		include those to draw bonds between atoms,
 *		to draw labels on atoms, to draw labels for
 *		residues etc.
 */




#define	DLBONDop	1
#define	DLTEXTop	2
#define	DLLINEop	3
#define	DLDISTANCEop	4
#define	DLRESIDUEop	5


typedef	struct	{
		ATOM		aAtom1;
		ATOM		aAtom2;
		BONDt		bOrder;
		} DLBONDt;

typedef	struct	{
		ATOM		aAtom;
		} DLATOMt;

typedef	struct	{
		GVECTOR		gvStart;
		GVECTOR		gvStop;
		int		iColor;
		} DLLINEt;

typedef	struct	{
		ATOM		aAtom1;
		ATOM		aAtom2;
		double		dDistance;
		} DLDISTANCEt;

typedef	struct	{
		RESIDUE		rRes;
		} DLRESIDUEt;


typedef	struct	DLENTRYs	{
				struct	DLENTRYs*	dleNext;
				int			iDLCommand;
				union	{
					DLBONDt		bond;
					DLATOMt		atom;
					DLLINEt		line;
					DLDISTANCEt	distance;
					DLRESIDUEt	residue;
					} data;
				} DLENTRYt;

typedef	DLENTRYt*	DLENTRY;


typedef	struct	{
		int		iSize;
		DLENTRY		dleFirst;
		DLENTRY		dleLoop;
		}
typedef	DLt*	DL;



/*
 *	Functions
 *

extern	DL	dlDLCreate();		/* () */
extern	void	DLDestroy();		/* ( DL ) */

extern	DLENTRY	dleDLAddEntry();	/* ( DL ) */
extern	void	DLEmpty();		/* ( DL ) */

extern	void	DLLoop();		/* ( DL ) */
extern	DLENTRY	dleDLNext();		/* ( DL ) */




/*
 *	Macros for accessing and altering display list elements
 */

#define	iDLCommand(dle)	(dle->iDLCommand)

#define	DLDefineBond( dle, a1, a2, bOrder ) (\
		dle->iDLCommand = DLBONDop,\
		dle->data.bond.aAtom1 = a1,\
		dle->data.bond.aAtom2 = a2,\
		dle->data.bond.bOrder = bOrder )
#define	aDLBondAtom1(dle)	(dle->data.bond.aAtom1)
#define	aDLBondAtom2(dle)	(dle->data.bond.aAtom2)
#define	bDLBondOrder(dle)	(dle->data.bond.bOrder)


#define	DLDefineAtom( dle, a1 ) (\
		dle->iDLCommand = DLATOMop,\
		dle->data.atom.aAtom = a1 )
#define	aDLAtomAtom( dle )	(dle->data.atom.aAtom)


#define	DLDefineLine( dle, gv1, gv2, iColor ) ( \
		dle->iDLCommand = DLLINEop,\
		dle->data.line.gvStart = gv1,\
		dle->data.line.gvStop = gv2,\
		dle->data.line.iColor = iColor )
#define	gvDLLineStart(dle)	(dle->data.line.gvStart)
#define	gvDLLineStop(dle)	(dle->data.line.gvStop)
#define	gvDLLineColor(dle)	(dle->data.line.iColor)


#define	DLDefineDistance(dle,a1,a2,dDist) (\
		dle->iDLCommand = DLLINEop,\
		dle->data.distance.aAtom1 = a1,\
		dle->data.distance.aAtom2 = a2,\
		dle->data.distance.dDistance = dDist )
#define	aDLDistanceAtom1(dle)	(dle->data.atom.aAtom1)
#define	aDLDistanceAtom2(dle)	(dle->data.atom.aAtom2)
#define	aDLDistanceAtom2(dle)	(dle->data.atom.dDistance)


#define	DLDefineResidue(dle,r) (\
		dle->iDLCommand = DLRESIDUEop,\
		dle->data.residue.rRes = r )
#define	rDLResidueResidue(dle)	(dle->data.residue.rRes)




















