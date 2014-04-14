/*
 *      File: bag.h
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
 *      Define BAG class.
 *	A BAG maintains a collection of objects in a linked
 *	list.  The objects are arbitrary pointers.
 *	This code was taken from list.c/list.h and
 *	modified to handle general objects.
 *
 */

#ifndef	BAG_H
#define BAG_H 

typedef	struct	BNodeStruct	{
	struct BNodeStruct	*bnPNextNode;
	GENP			PObject;
} BNODE;

typedef BNODE    *BNODEP;

typedef	struct	{
	BNODE	*bnPCur;
	BNODE	*bnPNext;
} BAGLOOP;

typedef struct  {
	int		iElementCount;
	BNODEP		bnPFirstNode;
	BNODEP		bnPLastNode;
} BAGt;
                
typedef BAGt	*BAG;


#define		iBagSize(b)            ( b->iElementCount )


/*  bag.c  */

extern BAG		bBagCreate();
extern void 		BagDestroy( BAG *bPBag );
extern void		BagAdd( BAG bBag, GENP PObj );
extern void		BagAddToEnd( BAG bBag, GENP PObj );
extern void		BagDescribe( BAG bBag );
extern BOOL		bBagRemove( BAG bBag, GENP PPtr );
extern BOOL		bBagContains( BAG bBag, GENP PPtr );
extern void		BagUniqueAdd( BAG bBag, GENP PObj );
extern BAGLOOP		blBagLoop( BAG bBag );
extern GENP		PBagNext( BAGLOOP *blPBagLoop );


#endif	/* BAG_H */

