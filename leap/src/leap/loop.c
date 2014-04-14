/*
 *      File: loop.c
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
 *              LOOP object is a general object which allows looping
 *              over elements of objects that act as container or collection
 *              objects.  LOOPs have two internal variables that define
 *              what they are:  oOver and iGoal.
 *              oOver contains the object over which the loop is
 *              being made, and iGoal contains information on what
 *              to search for and how to search for it.
 *
 *              The LOOP structure contains several fields that
 *              can be used to index the loop.  These fields
 *              can be used any way the programmer sees fit.
 *              More can be created as required.
 *
 */




#include	"basics.h"

#include        "classes.h"




/*
-------------------------------------------------------------------

        Define static variables here.
*/


static  int     SiUniqueId = 1;


/*
===================================================================

        Define private routines here.
*/



#define cCurLoopOver(lPL)       ((lPL)->oaSubLoopOver[(lPL)->iCurSubLoop])
#define PCurLoopNext(lPl)       ((lPl)->PaSubLoopList[(lPl)->iCurSubLoop])
#define CurLoopSetNext(lPl,N)   ((lPl)->PaSubLoopList[(lPl)->iCurSubLoop] =\
 (LISTLOOP)N )
#define PushSubLoop(lPl,Ov)     {\
(lPl)->iCurSubLoop++;\
(lPl)->oaSubLoopOver[(lPl)->iCurSubLoop]=Ov;\
if ( iObjectType(Ov) == LISTid ) \
  (lPl)->PaSubLoopList[(lPl)->iCurSubLoop]=llListLoop((LIST)Ov);\
  else\
  (lPl)->PaSubLoopList[(lPl)->iCurSubLoop]=llListLoop((LIST)cContainerContents(Ov));\
}
#define PopSubLoop(lPl)         ((lPl)->iCurSubLoop--)
#define bNoSubLoops(lPl)        ((lPl)->iCurSubLoop<0)



/*
 *      bLoopAtomVisible
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the atom flags satisfy the
 *      Visible/Invisible requirements of the loop.
 *      This means that the atom has ALL of the same flags 
 *      set as set in fVisibleFlagsOn and
 *      ALL the same flags reset as SET in fVisibleFlagsOff.
 */
static BOOL
bLoopAtomVisible( LOOP *lPLoop, ATOM aAtom )
{
FLAGS           fFlags;

    fFlags = fAtomFlags(aAtom);
    if ( (fFlags&lPLoop->fVisibleFlagsOn)!=lPLoop->fVisibleFlagsOn )
        return(FALSE);

    if ( (fFlags|(~(lPLoop->fVisibleFlagsOff)))!=(~(lPLoop->fVisibleFlagsOff)) )
        return(FALSE);

			/* If the TempInt field is being used to determine */
			/* visibility, then check if the ATOM has the */
			/* proper TempInt field */

    if ( lPLoop->fVisibilityFlags & TEMPINTUSED ) {
	if ( lPLoop->fVisibilityFlags & TEMPINTINVISIBLE ) {
	    if ( iAtomTempInt(aAtom) == lPLoop->iTempInt ) {
		return(FALSE);
	    }
	} else {
	    if ( iAtomTempInt(aAtom) != lPLoop->iTempInt ) {
		return(FALSE);
	    }
	}
    }	

		/* Check if the ATOM is the one invisible ATOM for the loop */
		/* Used to construct spanning trees that do not pass a */
		/* certian ATOM, good for looping over side chains */

    if ( lPLoop->aInvisibleAtom == (OBJEKT)aAtom ) return(FALSE);

    return(TRUE);
}
 



/*
 *      bSpanAtomVisible
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the atom is visible to a SPANNINGTREE loop.
 *      This means that the atom has not been seen before ( Set
 *      bPSeenBefore to TRUE if it has ) and that the atom
 *      has ALL of the same flags set as set in fVisibleFlagsOn and
 *      ALL the same flags reset as SET in fVisibleFlagsOff.
 */
static BOOL
bSpanAtomVisible( LOOP *lPLoop, ATOM aAtom, BOOL *bPSeenBefore )
{
    *bPSeenBefore = FALSE;
    if ( iAtomSeenId(aAtom) == lPLoop->iSeenId ) {
        *bPSeenBefore = TRUE;
        return(FALSE);
    }

		/* If the ATOM in the spanning tree is too far */
		/* from the root ATOM then say it is invisible */

    if ( lPLoop->iMaxDistanceFromRoot >= 0 ) {
	if ( lPLoop->iMaxDistanceFromRoot <
		iAtomBackCount(lPLoop->aCurSpan)+1 ) {
	    return(FALSE);
	}
    }
    return(bLoopAtomVisible( lPLoop, aAtom ));
}
 


/*
 *      bLoopSatisfiedBy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the loop is satisfied by the object.
 */
static BOOL
bLoopSatisfiedBy( LOOP *lPLoop, OBJEKT oObject )
{
int             iGoal;

    iGoal = lPLoop->iGoal;
    switch ( iGoal ) {
        case MOLECULES:
            return(( iObjectType(oObject) == MOLECULEid ));
        case RESIDUES:
            return(( iObjectType(oObject) == RESIDUEid ));
        case ATOMS:
            return(( iObjectType(oObject) == ATOMid ));
        case INTERNALS:
            return(( iObjectType(oObject) == INTERNALid ));
        case CONTAINERS:
            return(( bObjectInClass( oObject, CONTAINERid ) ));
        case DIRECTCONTENTS:
            return(TRUE);
    }
    return(FALSE);
}






/*
 *      InitLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Initialize the loop counters.
 */
static void
InitLoop( LOOP *lPLoop )
{
int             i;

    lPLoop->iIndex0   = 0;
    lPLoop->iIndex1   = 0;
    lPLoop->iIndex2   = 0;
    lPLoop->iIndex3   = 0;
    for ( i=0; i<MAXOBJ; i++ ) lPLoop->oaObj[i] = NULL;
}          

/*
 *      bNextObjectInAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if there is another BOND, ANGLE, TORSION, IMPROPER
 *      connected to the atom, FALSE if there is not.
 *      If there is such a thing then place information in the loop
 *      so that it may be found again.
 */
static BOOL
bNextObjectInAtom( LOOP *lPLoop )
{
CONTAINER       cCont;
ATOM            aAtom1, aAtom2, aAtom3;
BOOL            bDone, bAllowDuplicates;

    cCont = (CONTAINER)cCurLoopOver(lPLoop);

		/* Check if we should allow duplicate bonds, angles, and */
		/* torsions */

    bAllowDuplicates = ( (lPLoop->iGoal&ALLOWDUPLICATES) != 0 );

    switch ( (lPLoop->iGoal&GOALONLY) ) {

                        /* When LOOPing over BONDS, use iIndex0 as      */
                        /* the bond index for the next ATOM to check    */
                        /* keep looking until a suitable bond is found  */
                        /* Eg: where id(ATOM) < id(BONDED_ATOM)         */
                        /* and place the two atoms in oaObj[0] and      */
                        /* oaObj[1] of the Top Loop                     */
                        /* then return the atom we are looping over     */
                        /* or NULL if there are no more bonds.          */
                        /* The atoms in the bonds can then be retrieved */
                        /* Using the LoopGetBond routine                */
        case BONDS:
                while ( lPLoop->iIndex0 < iAtomCoordination(cCont) ) {
                    if ( (iAtomId(cCont) < 
                iAtomId(aAtomBondedNeighbor(cCont,lPLoop->iIndex0))) ||
			  bAllowDuplicates ) {
                        lPLoop->oaObj[0] = (OBJEKT)cCont;
                        lPLoop->oaObj[1] = (OBJEKT)
                                aAtomBondedNeighbor(cCont,lPLoop->iIndex0);
                        lPLoop->iIndex0++;
                        return(TRUE);
                    }
                    lPLoop->iIndex0++;
                }
                break;

			/* When LOOPing over ANGLES, use iIndex0 and    */
                        /* iIndex1 as the bond angle indices.           */
                        /* An outer ATOM is the current ATOM.           */
                        /* iIndex0 ranges through the current atoms     */
                        /* coordination and iIndex1, ranges from 0	*/
                        /* to the current atoms coordination.		*/
                        /* When an angle is found place the three atoms */
                        /* in order in oaObj, and return the current    */
                        /* atom.  Otherwise return NULL                 */
                        /* The atoms of the angle are retrieved using   */
                        /*the LoopGetAngle routine                      */
        
	case ANGLES:
                bDone = FALSE;
                while ( ! bDone ) {
                                /* First test if all the indices are valid */
                    if ( lPLoop->iIndex0>=iAtomCoordination(cCont) ) goto ANONE;
                    aAtom1 = aAtomBondedNeighbor( cCont, lPLoop->iIndex0 );
                    if ( lPLoop->iIndex1>=iAtomCoordination(aAtom1) ) goto AINC1;
                    aAtom2 = aAtomBondedNeighbor( aAtom1, lPLoop->iIndex1 );

                                /* Test if the id's are in the right order */

                    if ( !(iAtomId(cCont) < iAtomId(aAtom2)) &&
			 !bAllowDuplicates ) goto AINC1;

                                /* Test if the atoms are unique */

                    if ( (ATOM)cCont ==aAtom2 ) goto AINC1;

                                /* If it passed all these tests then    */
                                /* it is a valid angle                 */

                    bDone = TRUE;

                                /* Increment everything */
AINC1:                   
                    lPLoop->iIndex1++;
                    if ( lPLoop->iIndex1 >= iAtomCoordination(aAtom1) ) {
                        lPLoop->iIndex0++;
                        lPLoop->iIndex1 = 0;
                    }
                }
                        /* An angle was found, place it in the oaObj array */
                lPLoop->oaObj[0] = (OBJEKT)cCont;
                lPLoop->oaObj[1] = (OBJEKT)aAtom1;
                lPLoop->oaObj[2] = (OBJEKT)aAtom2;
                return(TRUE);
ANONE:          break;


                        /* When LOOPing over TORSIONs, use iIndex0,     */
                        /* iIndex1, and iIndex2 as indices for the      */
                        /* next torsion.                                */
                        /* iIndex0 gives the bond index from the current*/
                        /* atom to the second atom in the torsion.      */
                        /* iIndex1 gives the bond index from the second */
                        /* to the third and iIndex2 is from the third   */
                        /* to the fourth                                */
                        /* Keep searching for torsions until there      */
                        /* are no more bonds to search on on the first  */
                        /* atom.  Return only torsions where            */
                        /* id(first) < id(last) and where every atom in */
                        /* the torsion is unique                        */
                        /* The search is done using a test/increment    */
                        /* sequence.                                    */
        case PROPERS:
                bDone = FALSE;
                while ( ! bDone ) {
                                /* First test if all the indices are valid */
                    if ( lPLoop->iIndex0>=iAtomCoordination(cCont) ) goto TNONE;
                    aAtom1 = aAtomBondedNeighbor( cCont, lPLoop->iIndex0 );
                    if ( lPLoop->iIndex1>=iAtomCoordination(aAtom1) ) goto TINC1;
                    aAtom2 = aAtomBondedNeighbor( aAtom1, lPLoop->iIndex1 );
                    if ( lPLoop->iIndex2>=iAtomCoordination(aAtom2) ) goto TINC2;
                    aAtom3 = aAtomBondedNeighbor( aAtom2, lPLoop->iIndex2 );
                                /* Test if the id's are in the right order */
                    if ( !(iAtomId(cCont) < iAtomId(aAtom3)) &&
				!bAllowDuplicates ) goto TINC2;
                                /* Test if the atoms are unique */
                    if ( (ATOM)cCont ==aAtom2 ) goto TINC1;
                    if ( (ATOM)cCont == aAtom3 ) goto TINC2;
                    if ( aAtom1 == aAtom3 ) goto TINC2;
                                /* If it passed all these tests then    */
                                /* it is a valid torsion                */
                    bDone = TRUE;
                                /* Increment everything */
TINC2:
                    lPLoop->iIndex2++;
                    if ( lPLoop->iIndex2 >= iAtomCoordination(aAtom2) ) {
TINC1:                   
                        lPLoop->iIndex1++;
                        lPLoop->iIndex2 = 0;
                        if ( lPLoop->iIndex1 >= iAtomCoordination(aAtom1) ) {
                            lPLoop->iIndex0++;
                            lPLoop->iIndex1 = 0;
                        }
                    }
                }
                        /* A torsion was found, place it in the oaObj array */
                lPLoop->oaObj[0] = (OBJEKT)cCont;
                lPLoop->oaObj[1] = (OBJEKT)aAtom1;
                lPLoop->oaObj[2] = (OBJEKT)aAtom2;
                lPLoop->oaObj[3] = (OBJEKT)aAtom3;
                return(TRUE);
TNONE:          break;

                        /* When LOOPing over IMPROPERS, use iIndex0 and */
                        /* iIndex1,iIndex2 as the improper indices.     */
                        /* The central atom of the improper is the current */
                        /* improper being looped over.                  */
                        /* iIndex0 ranges through the current atoms     */
                        /* coordination and iIndex1, ranges from iIndex0*/
                        /* +1 to the current atoms coordination.        */
                        /* iIndex2 ranges form iIndex1+1 to the current */
                        /* atoms coordination                           */
                        /* When an angle is found place the three atoms */
                        /* in order in oaObj, and return the current    */
                        /* atom.  Otherwise return NULL                 */
                        /* The atoms of the angle are retrieved using   */
                        /*the LoopGetImproper routine                   */
                        /* The central atom is stored in oaObj[2], the  */
                        /* third atom in the list of four.  This is to  */
                        /* maintain consistancy with the way impropers  */
                        /* are defined in PARM89A.DAT                   */
        case IMPROPERS:
                        /* Anticipate that the indexes all start at 0 */
                if ( lPLoop->iIndex1 <= lPLoop->iIndex0 ) 
                                lPLoop->iIndex1 = lPLoop->iIndex0+1;
                if ( lPLoop->iIndex2 <= lPLoop->iIndex1 ) 
                                lPLoop->iIndex2 = lPLoop->iIndex1+1;
                else lPLoop->iIndex2++;
                if ( lPLoop->iIndex2 >= iAtomCoordination(cCont) ) {
                    lPLoop->iIndex1++;
                    lPLoop->iIndex2 = lPLoop->iIndex1 + 1;
                    if ( lPLoop->iIndex2 >= iAtomCoordination(cCont) ) {
                        lPLoop->iIndex0++;
                        lPLoop->iIndex1 = lPLoop->iIndex0 + 1;
                        lPLoop->iIndex2 = lPLoop->iIndex1 + 1;
                    }
                }
                if ( lPLoop->iIndex2 < (iAtomCoordination(cCont)) ) {
                    lPLoop->oaObj[2] = (OBJEKT)cCont;
                    lPLoop->oaObj[0] = (OBJEKT)
                        aAtomBondedNeighbor( cCont, lPLoop->iIndex0 );
                    lPLoop->oaObj[1] = (OBJEKT)
                        aAtomBondedNeighbor( cCont, lPLoop->iIndex1 );
                    lPLoop->oaObj[3] = (OBJEKT)
                        aAtomBondedNeighbor( cCont, lPLoop->iIndex2 );
                    return(TRUE);
                }
                break;

        default:
                break;
    }
    
    InitLoop(lPLoop);
    return(FALSE);
}




/*
 *      iLoopContainerMatch
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Compare the two CONTAINERs pointed to
 *      and return -1, 0, 1 if the sequence number of cA,cB
 *      are cA<cB, cA==cB, cA>cB.
 */
static int
iLoopContainerMatch( CONTAINER *cPA, CONTAINER *cPB )
{
    if ( iContainerSequence(*cPA)<iContainerSequence(*cPB) ) 
	return(-1);
    if ( iContainerSequence(*cPA)==iContainerSequence(*cPB) ) 
	return(0);
    return(1);
}





/*
 *	zLoopAddToMemory
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add an OBJEKT to the LOOP memory.
 */
static void
zLoopAddToMemory( LOOP *lPLoop, OBJEKT oObject )
{
LOOPNODE	lnNew;

    MALLOC( lnNew, LOOPNODE, sizeof(LOOPNODEt) );
    if ( lPLoop->lnLast != NULL ) {
	lPLoop->lnLast->lnNext = lnNew;
    } else {
	lPLoop->lnMemory = lnNew;
    }
    lPLoop->lnLast = lnNew;
    lnNew->lnNext = NULL;
    lnNew->oData = oObject;
}




/*
*******************************************************************

        Define public message routines here.
*/



/*
 *      lLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an empty loop, over nothing.
 *      Further calls must be made to define what the loop is over
 *      and what it's goal is.
 *
 *      A loop over DIRECTCONTENTS returns only the OBJEKTs which
 *      are DIRECTLY contained within the cOver.
 *
 *      A LOOP over DIRECTCONTENTSBYSEQNUM returns only the OBJEKTS
 *      which are DIRECTLY contained within cOver, sorted by
 *      sequence number in accending order.
 */
LOOP
lLoop( OBJEKT oOver, int iGoal )
{
LOOP            lL, lTemp;
int             i;
VARARRAY        vaSeqNum;
CONTAINER       *cPCur;
CONTAINER	cTemp;

    if ( oOver == NULL ) {
         DFATAL( ("Attempt to define LOOP over NULL") );
    }
    memset(&lL, 0, sizeof(lL));		/* for Purify */

    MESSAGE(( "===Creating LOOP\n" ));
    lL.bInitialized 	  = FALSE;
    lL.bLoopDone          = FALSE;

		/* Initialize the LOOP memory fields to the default */
		/* of no memory. */

    lL.bUsingMemory	  = FALSE;
    lL.bReplayingMemory   = FALSE;
    lL.lnMemory  = NULL;
    lL.lnLast = NULL;

    lL.fVisibilityFlags = 0;

    lL.iCurSubLoop = -1;
    lL.oOver     = oOver;
    lL.iGoal     = iGoal;
    lL.fVisibleFlagsOn = 0;
    lL.fVisibleFlagsOff = 0;
    lL.aInvisibleAtom   = NULL;

    InitLoop(&lL);

                /* If the LOOP is by sequence number then first sort */
                /* the CONTAINERs by sequence number and build a */
                /* linked list */
    if ( iGoal == DIRECTCONTENTSBYSEQNUM ) {
        i = iCollectionSize( cContainerContents(oOver) );
	if ( i != 0 ) {
            vaSeqNum = vaVarArrayCreate( sizeof(CONTAINER) );
	    TESTMEMORY();
            VarArraySetSize( vaSeqNum, i );

                /* Fill the array with the CONTAINERs */
	    TESTMEMORY();
            cPCur = PVAI( vaSeqNum, CONTAINER, 0 );
            lTemp = lLoop( oOver, DIRECTCONTENTS );
            while ( (cTemp = (CONTAINER)oNext(&lTemp)) ) {
		*cPCur = cTemp;
		cPCur++;
	    }
	    TESTMEMORY();
            qsort( PVAI( vaSeqNum, CONTAINER, 0 ), i, sizeof(CONTAINER),
                (int (*) (const void *, const void *) )iLoopContainerMatch );
	    TESTMEMORY();
                /* Create the sorted linked list */
            cPCur = PVAI( vaSeqNum, CONTAINER, 0 );
            lL.oaSubLoopOver[0] = (OBJEKT)(*cPCur);
	    TESTMEMORY();
            for ( i=0; i< iVarArrayElementCount(vaSeqNum)-1; i++ ) {
                ContainerSetLoopNext( *cPCur, *(cPCur+1) );
                cPCur++;
            }
            ContainerSetLoopNext( *cPCur, NULL );
            VarArrayDestroy( &vaSeqNum );
	} else {
	    lL.oaSubLoopOver[0] = (OBJEKT)NULL;
	}
	TESTMEMORY();
    } else if ( (iGoal&WAYONLY)==SPANNINGTREE ) {

                /* Set up the LOOP for a spanning tree */

        VERIFYOBJEKT( oOver, ATOMid );
	lL.iMaxDistanceFromRoot = -1;
        lL.aCurSpan  = oOver;
        lL.aLastSpan = oOver;
        lL.iSeenId   = iLoopNextSeenIdInc();
        lL.iInvisibleCollisions = 0;
        lL.aLastCollision = NULL;
    } else {
	PushSubLoop( &lL, oOver );
    }
    return(lL);
}





/*
 *      oNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the next element in the LOOP that satisfies the goal.
 *
 */
OBJEKT
oNext( LOOP *lPLoop )
{
OBJEKT                  oObject;
int                     iGoal, i;
LISTLOOP                llPLoop;
BOOL                    bSeenBefore;
ATOM                    aPrev, aBond;

		/* At the end of the function, there is code */
		/* to get the next OBJEKT out of the memory list */


    if ( lPLoop->bReplayingMemory ) goto DONE;

                /* Handle the loop differently if it's goal is */
                /* DIRECTCONTENTSBYSEQNUM */

    if ( lPLoop->iGoal == DIRECTCONTENTSBYSEQNUM ) {
        oObject = lPLoop->oaSubLoopOver[0];
        if ( oObject == NULL ) goto DONE;
        lPLoop->oaSubLoopOver[0] = (OBJEKT)cContainerLoopNext(oObject);
        goto DONE;
    }

                /* Handle the loop differently if it is a SPANNINGTREE */

    if ( (lPLoop->iGoal&WAYONLY) == SPANNINGTREE ) {

			/* The first time the LOOP is entered, set up */
			/* some stuff on the first ATOM */

	if ( !lPLoop->bInitialized ) {
	    AtomSetBackSpan( lPLoop->oOver, NULL );
            AtomSetBackCount( lPLoop->oOver, 0 );
            AtomSetSeenId( lPLoop->oOver, lPLoop->iSeenId );
	}

        if ( lPLoop->aCurSpan==NULL ) {
	    oObject = NULL;
	    goto DONE;
	}

                /* Connect all visible atoms bonded to the current atom */
                /* into the spanning tree */

        aPrev = (ATOM)lPLoop->aLastSpan;
        AtomSetNextSpan( aPrev, NULL );
        for ( i=0; i<iAtomCoordination(lPLoop->aCurSpan); i++ ) {
            aBond = aAtomBondedNeighbor( lPLoop->aCurSpan, i );
            MESSAGE(( "---Looking at: %s\n", sContainerName(aBond) ));

                        /* If the atom is visible then add it */
            if ( bSpanAtomVisible( lPLoop, aBond, &bSeenBefore ) ) {
                MESSAGE(( "--- Visible\n" ));
                AtomSetSeenId( aBond, lPLoop->iSeenId );
                AtomSetBackCount( aBond, iAtomBackCount(lPLoop->aCurSpan)+1 );
                AtomSetBackSpan( aBond, (ATOM)lPLoop->aCurSpan );
                AtomSetNextSpan( aPrev, aBond );
                AtomSetNextSpan( aBond, NULL );
                aPrev = aBond;
            } else {
                MESSAGE(( "--- NOT Visible\n" ));

                        /* If the atom is invisible, but has not been seen */
                        /* before then it counts as an invisible collision */
                        /* Increment the collision count and save the atom */
                        /* which caused the collision                      */
                if ( !bSeenBefore ) {
                    MESSAGE(( "--- COLLISION!\n" ));
                    lPLoop->iInvisibleCollisions++; 
                    lPLoop->aLastCollision = (OBJEKT)aBond;
                }
            }
        }

                /* Advance to the next atom and return the current one */
        lPLoop->aLastSpan = (OBJEKT)aPrev;
        oObject = lPLoop->aCurSpan;
        lPLoop->aCurSpan = (OBJEKT)aAtomNextSpan(lPLoop->aCurSpan);

		/* Return with oObject */

        goto DONE;

    } else if ( (lPLoop->iGoal & GOALONLY) != iObjectType(lPLoop->oOver) ) {

        iGoal = ( lPLoop->iGoal & GOALONLY );

                /* If the object being looped over is an ATOM and the goal */
                /* is not another OBJEKT eg: the goal is a bond, angle,    */
                /* torsions, impropers then handle differently */

        while (1) {
                /* Get the next object in the list             */
                /* Advance the list pointer                    */
                /* Then check if the object satisfies the goal */
                /* if not then create a subloop over it        */
                /* If the object over which the loop is taken is */
                /* an atom then do some special searching       */

            if ( (iObjectType(cCurLoopOver(lPLoop))==ATOMid) &&
                        ((iGoal&NONCONTAINERLOOP)!=0) ) {
                if ( bNextObjectInAtom(lPLoop) ) {
                    oObject = lPLoop->oOver;
                    goto DONE;
                } else oObject = NULL;
            } else {
                llPLoop = (LISTLOOP)PCurLoopNext(lPLoop);
                oObject = oListNext( &llPLoop );
                CurLoopSetNext( lPLoop, llPLoop );
            }
    
                    /* If there are no more elements in the list then */
                    /* pop this subloop                               */
    
            if ( oObject==NULL ) {
                PopSubLoop(lPLoop);
                if ( bNoSubLoops(lPLoop) ) {
                    oObject = NULL;
                    goto DONE;
                }
                continue;
            }
    
            if ( bLoopSatisfiedBy( lPLoop, oObject ) ) {
            
                    /* If LOOPing over CONTAINERS then create a subloop */
                    /* over the CONTAINER along with returning the      */
                    /* CONTAINER                                        */

                if ( (iGoal&GOALONLY) == CONTAINERS ) 
                    if ( iObjectType(oObject) != ATOMid ) 
                                PushSubLoop( lPLoop, oObject );

                        /* If the object that satisfies the loop is an */
                        /* atom then make sure that it satisfies the   */
                        /* flags */

                if ( iObjectType(oObject) == ATOMid ) 
                    if ( !bLoopAtomVisible( lPLoop, (ATOM)oObject )) continue;

                goto DONE;
            }
    
                    /* The current object does not satisfy the LOOP */
                    /* So use it to construct a subloop with the    */
                    /* same goal                                    */
    
            if ( bObjectInClass( oObject, CONTAINERid ) )
                PushSubLoop( lPLoop, oObject );

        } 

    } else {

			/* If the goal of the loop is the same type */
			/* of OBJEKT as the loop is over then there is */
			/* only one OBJEKT that will satisfy the loop */
			/* The OBJECT (oOver) */

	if ( !lPLoop->bInitialized ) {
	    lPLoop->bInitialized = TRUE;
	    oObject = lPLoop->oOver;
	    goto DONE;
	}
	oObject = NULL;
	goto DONE;
    }



DONE:

		/* If replaying the memory, then get the next */
		/* OBJEKT out of the memory list and advance */
		/* the list */
    if ( lPLoop->bReplayingMemory ) {
	if ( lPLoop->lnLast != NULL ) {
	    oObject = lPLoop->lnLast->oData;
	    lPLoop->lnLast = lPLoop->lnLast->lnNext;
	} else oObject = NULL;
    } else if ( lPLoop->bUsingMemory ) {
	if ( oObject != NULL ) zLoopAddToMemory( lPLoop, oObject );
    }

    return(oObject);
}



/*
 *	LoopDestroyMemory
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the LOOP memory.
 */
void
LoopDestroyMemory( LOOP *lPLoop )
{
LOOPNODE	lnCur, lnFree;

    lnCur = lPLoop->lnMemory;
    while ( lnCur != NULL ) {
	lnFree = lnCur;
	lnCur = lnCur->lnNext;
	FREE(lnFree);
    }
    lPLoop->lnMemory = NULL;
    lPLoop->bUsingMemory = FALSE;
    lPLoop->bReplayingMemory = FALSE;
}

