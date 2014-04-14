/*
 *      File: loop.h
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
 *              Loop object, allows looping over elements of any
 *              object which acts as a container/collection class.
 *
 *              LOOP is not a class like other MolecularObjects
 *              it does not share the OBJEKT header.
 *              It is only used to loop over objects.
 *              The main reason for this is because LOOPs do not need
 *              to be explicitly Destroyed, they do
 *              not need to have the same persistance that REAL objects
 *              have.  LOOPs are NOT stored on the heap! they are stored
 *              on the stack.  Not treating LOOPs as objects saves time and
 *              memory.
 *
 *		LOOPs can also be set up to record the OBJEKTs
 *		that they generate and then be told to replay
 *		their contents using (oNext).
 *
 *		LoopUseMemory(LOOP*) is used after (lLoop) to
 *			start recording OBJEKTs.
 *		LoopRewindMemory(LOOP*) is used after the LOOP
 *			has been run through once to rewind the
 *			pointers to the start of the memory.
 *		(oNext) is used to get each successive element of
 *			the memory.
 *		LoopDestroyMemory(LOOP*) is used to free the 
 *			space used by the memory.
 */
 


#ifndef	LOOP_H
#define	LOOP_H

        /* LOOP subgoals,  Complete goals are constructed by bit ORing */
        /* a loop goal with a loop type */
        /* If no loop type is specified then the default is taken */


#define GOALONLY                0x0000FFFF
#define MOLECULES               MOLECULEid
#define RESIDUES                RESIDUEid       
#define ATOMS                   ATOMid
#define CONTAINERS              CONTAINERid
#define INTERNALS               INTERNALid
#define DIRECTCONTENTS          0x00000100
#define DIRECTCONTENTSBYSEQNUM  0x00000200

#define NONCONTAINERLOOP        0x00001000
#define BONDS                   0x00001001
#define ANGLES                  0x00001002
#define PROPERS			0x00001003      
#define IMPROPERS               0x00001004

        /* LOOP ways */

		/* ALLOWDUPLICATES allows the loops over BONDS, */
		/* ANGLES, TORSIONS to return duplicates.  Eg:  */
		/* A bond between atoms A-B can be returned from */
		/* both A and B */

#define WAYONLY                 0xFFFF0000
#define DEFAULT                 0x00000000
#define SPANNINGTREE            0x00010000
#define	ALLOWDUPLICATES		0x00020000


#define MAXOBJ          4
#define MAXSUBLOOPS     5


		/* Flags that control what determines ATOM visibility */

#define	TEMPINTVISIBLE		0x00000001
#define	TEMPINTINVISIBLE	0x00000002
#define	TEMPINTUSED		0x00000004


typedef	struct	LOOPNODESTRUCT{
		struct	LOOPNODESTRUCT*	lnNext;
		OBJEKT			oData;
		} LOOPNODEt;

typedef	LOOPNODEt*	LOOPNODE;

typedef struct  LOOPSTRUCT{
	OBJEKT          oOver;
	BOOL		bInitialized;
	BOOL            bLoopDone;
	int             iGoal;
	int             iIndex0, iIndex1, iIndex2, iIndex3;
	OBJEKT          oaObj[MAXOBJ];
	int             iCurSubLoop;
	OBJEKT          oaSubLoopOver[MAXSUBLOOPS];
	LISTLOOP        PaSubLoopList[MAXSUBLOOPS];
    /* Spanning tree stuff */
	OBJEKT          aCurSpan;
	OBJEKT          aLastSpan;
	int             iSeenId;
	int             fVisibleFlagsOn;
	int             fVisibleFlagsOff;
	int		iMaxDistanceFromRoot;
	OBJEKT		aInvisibleAtom;
	int             iInvisibleCollisions;
	OBJEKT          aLastCollision;
	FLAGS		fVisibilityFlags;
	int		iTempInt;
    /* Loop memory fields */
	LOOPNODE	lnMemory;
	LOOPNODE	lnLast;
	BOOL		bUsingMemory;
	BOOL		bReplayingMemory;
} LOOP;




/*
======================================================================

        This is not a regular OBJEKT, there is only a constructor
        'lLoop' which returns a structure, and 'oNext' which returns
        the next element of the loop.  Because 'lLoop' returns a structure
        there does not need to be an explicite destructor.

*/


extern LOOP		lLoop(OBJEKT oOver, int iGoal);
extern OBJEKT		oNext(LOOP *lPLoop);
extern void		LoopDestroyMemory(LOOP *lPLoop);

#define oLoopOver( LPL )                ( (LPL)->oOver )
#define iLoopGoal( LPL )                ( (LPL)->iGoal )
#define LoopGetBond( xl,aP1,aP2 )\
    ((*(aP1))=(ATOM)(xl)->oaObj[0],(*(aP2))=(ATOM)(xl)->oaObj[1] )
    
#define LoopGetAngle( xl,aP1,aP2,aP3 )  \
    ((*(aP1))=(ATOM)(xl)->oaObj[0],(*(aP2))=(ATOM)(xl)->oaObj[1],\
    (*(aP3))=(ATOM)(xl)->oaObj[2])
    
#define LoopGetTorsion( xl,aP1,aP2,aP3,aP4 )    \
((*(aP1))=(ATOM)(xl)->oaObj[0],(*(aP2))=(ATOM)(xl)->oaObj[1],\
(*(aP3))=(ATOM)(xl)->oaObj[2],(*(aP4))=(ATOM)(xl)->oaObj[3])

                /* NOTE: In LoopGetImproper, the central atom is */
                /* returned in aP3 !!!!!!!!!!!!! */
		/* (this is the same as LoopGetTorsion for now, maybe forever */

#define LoopGetImproper( xl,aP1,aP2,aP3,aP4 )   \
((*(aP1))=(ATOM)(xl)->oaObj[0],(*(aP2))=(ATOM)(xl)->oaObj[1],\
(*(aP3))=(ATOM)(xl)->oaObj[2],(*(aP4))=(ATOM)(xl)->oaObj[3])

		/* LOOP tags are used by SPANNINGTREE to keep track */
		/* of what ATOMs have been touched by the SPANNINGTREE */
		/* Using the Tag macros the caller can fiddle with */
		/* the SeenId to build multiple SPANNINGTREE's that don't */
		/* overlap. */

#define	iLoopSeenId(l)			( l.iSeenId )
#define	LoopUseSeenId(l,i)		( l.iSeenId = i )

#define	iLoopNextSeenId()		( SiUniqueId )
#define	iLoopNextSeenIdInc()		( SiUniqueId++ )

#define	LoopDefineMaxDistanceFromRoot(LPL,FF) \
			( (LPL)->iMaxDistanceFromRoot = (FF) )
#define LoopDefineVisibleAtoms(LPL,FF)   ( (LPL)->fVisibleFlagsOn=(FF))
#define LoopDefineInvisibleAtoms(LPL,FF) ( (LPL)->fVisibleFlagsOff=(FF))
#define iLoopInvisibleCollisionCount(LPL) ( (LPL)->iInvisibleCollisions )
#define aLoopLastCollisionAtom(LPL)       (ATOM)( (LPL)->aLastCollision )
#define	LoopDefineInvisibleAtom(LPL,A)	( (LPL)->aInvisibleAtom = (OBJEKT)(A) )


		/* The following two macros can be used to define */
		/* whether the TempInt property of an ATOM is used */
		/* to determine visibility OR invisibility of a particular ATOM */

#define	LoopDefineOnlyVisibleTempInt(LPL,A) \
		( (LPL)->fVisibilityFlags |= (TEMPINTUSED|TEMPINTVISIBLE),\
		  (LPL)->iTempInt = (A) )
#define	LoopDefineOnlyInvisibleTempInt(LPL,A) \
		( (LPL)->fVisibilityFlags |= (TEMPINTUSED|TEMPINTINVISIBLE),\
		  (LPL)->iTempInt = (A) )


		/* LOOP memory macros. */
		/* LoopUseMemory starts recording, */
		/* LoopRewindMemory starts replaying */

#define	LoopUseMemory(lPL)	((lPL)->bUsingMemory = TRUE )
#define	LoopRewindMemory(lPL)	((lPL)->bReplayingMemory = TRUE,\
				 (lPL)->lnLast = (lPL)->lnMemory )

                /* Some macros that make LOOPs a little easier to use */

#define FOREACH(x,type,l)       while ( ( x = (type)oNext(&l) ) != NULL )
#define LOOPOVERALL(over,goal,element,type,loop) \
        loop = lLoop((OBJEKT)over,goal);\
        while ( ( element=(type)oNext(&loop) ) != NULL )


#endif/* LOOP_H */
