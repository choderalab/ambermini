/*
 *	File:	ring.h
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
 *		Maintain a RING structure.  A doubly linked list of
 *		nodes where insertions can be made relative to other
 *		elements within the RING.  Each node of the RING
 *		has a pointer to an arbitrary data type.
 */


#ifndef	RING_H
#define	RING_H



typedef	struct	RINGNODEs {
	GENP			PData;
	struct RINGNODEs	*rnPPrev;
	struct RINGNODEs	*rnPNext;
} RINGNODEt;

typedef	struct	{
	RINGNODEt		*rnPFirst;
	int			iElements;
} RINGt;

typedef	struct	{
	RINGNODEt	*rnPNext;
	RINGNODEt	*rnPCur;
	int		iCount;
} RINGLOOP;

typedef	RINGt *		RING;



extern RING	rRingCreate(void);
extern void	RingDestroy(RING *rPRing);

#define	bRingEmpty(r)	(((RING)(r))->rnPFirst == NULL )
#define	iRingSize(r)	(((RING)(r))->iElements)

extern void	RingAfterAdd(RING rRing, GENP PAfter, GENP PData);
extern void	RingBeforeAdd(RING rRing, GENP PBefore, GENP PData);

extern void	rlRingLoop(RING rRing, RINGLOOP *rlP);
extern GENP	PRingNext(RINGLOOP *rlPLoop);
extern BOOL	bRingRemove(RING rRing, GENP PData);

#define	rnPRingNode(rl)		(rl.rnPCur)




#endif	/* RING_H */














