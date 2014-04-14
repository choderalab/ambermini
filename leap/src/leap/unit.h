/*
 *      File: unit.h
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
 *              UNIT
 *      Superclass: 
 *              CONTAINER
 *
 *      Description:
 *
 *              UNITs can contain molecules residues and/or atoms
 *              also parameters.
 */
#ifndef UNIT_H
#define UNIT_H

/*
 *-----------------------------------------------------------------------
 *
 *       Define object typedefs here.
 *        
 *        Object typedef MUST include the superclass object as
 *        its first structure element.
 *        
 *        UNITs have a mode, in UNITNORMAL mode, the unit is defined
 *        by the contents of the container and the pointers between
 *        the contents.  In UNITNORMAL mode the coordinates of the atoms
 *        are EXTERNAL coordinates.
 *        In UNITTABLES mode the UNIT is represented by the tables
 *        that have been read in, or generated from the contents.
 *        In UNITINTERNAL the UNIT is represented by the contents and 
 *        pointers between the contents, the coordinates of the atoms
 *        are represented by INTERNAL coordinates attached to the
 *        atoms.
 */


#include        "bag.h"
#include        "restraint.h"
#include        "dictionary.h"
#include        "parmLib.h"


#define UNITNORMAL      1
#define UNITTABLES      2
#define UNITINTERNAL    3


        /* UNIT flags */

#define UNITALLFLAGS            0xFFFFFFFF
#define UNITUSEBOUNDINGBOX      0x00000001
#define UNITBEINGEDITED         0x00000002
#define UNITUSESOLVENTCAP       0x00000004
#define UNITBOXOCT              0x00000008




typedef struct  {
        CONTAINERt      cHeader;

/* TODO: Get rid of aaConnect in UNITs, instead use a first atom/last atom */
/* TODO: paradigm */
        OBJEKT          aHead;
        OBJEKT          aTail;
        PARMSET         psParameters;
        FLAGS           fFlags;
        int             iMode;
        double          dBeta;
        double          dXWidth;
        double          dYWidth;
        double          dZWidth;
        VECTOR          vCapOrigin;
        double          dCapRadius;

                        /* Stuff to maintain named groups of ATOMs */
                        
        DICTIONARY      dAtomGroups;

                        /* Stuff to maintain RESTRAINTs */

        BAGLOOP         blRestraintLoop;
        BAG             bRestraints;
                
                        /* The following are PRIVATE, they are used */
                        /* to store interactions and pointers to parameters */
                        /* for writing UNITs to DATABASEs and for writing */
                        /* UNITs with parameters for SPASMS */

        int             iCapTempInt;
        VARARRAY        vaAtoms;
        VARARRAY        vaBonds;
        VARARRAY        vaAngles;
        VARARRAY        vaTorsions;
        VARARRAY        vaConnectivity;
        VARARRAY        vaRestraints;
        VARARRAY        vaResidues;
        VARARRAY        vaMolecules;
        VARARRAY        vaHierarchy;
        VARARRAY        vaConnect;
        VARARRAY        vaGroupNames;
        VARARRAY        vaGroupAtoms;
} UNITt;

typedef UNITt   *UNIT;





/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primitive superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */

extern UNIT     uUnitCreate();
extern void     UnitDelete(UNIT *uPUnit);
extern void     UnitDescribe(UNIT uUnit);
extern UNIT     uUnitDuplicate(UNIT uOld);
extern UNIT     uUnitCopy(UNIT uOld);
extern void     UnitResetPointers(UNIT uUnit);

extern void     UnitJoin(UNIT uA, UNIT uB);
extern void     UnitSequence(UNIT uA, UNIT uB);


extern void     UnitSave(UNIT uUnit, DATABASE db, PARMLIB plParameters);
extern UNIT     uUnitLoad(DATABASE db);
extern void     UnitCheck(UNIT uUnit, int *iPErrors, int *iPWarnings);
extern void     UnitCheckForParms( UNIT uUnit, PARMLIB plParms, 
                        PARMSET psParmSet );

extern void     UnitSaveAmberParmFile(UNIT uUnit, FILE *fOut, char *crdName,
                        PARMLIB plParms, BOOL bPolar, BOOL bPert, BOOL bNetcdf);

extern void     UnitYouAreBeingRemoved(UNIT uUnit);
extern void     UnitIAmBeingRemoved(UNIT uUnit, CONTAINER cRemoved);

extern BOOL     bUnitCanBePerturbed(UNIT uUnit);

                /* Restraint add/remove/loop */

extern void     UnitAddRestraint(UNIT uUnit, RESTRAINT rRest );
extern BOOL     bUnitRemoveRestraint(UNIT uUnit, RESTRAINT rRest);
extern void     UnitLoopRestraints(UNIT uUnit);
extern RESTRAINT        rUnitNextRestraint(UNIT uUnit);

extern int      iUnitRestraintTypeCount(UNIT uUnit, int iType);

extern void     UnitSetAttribute(UNIT uUnit, STRING sAttr, OBJEKT oAttr);

extern BOOL     bUnitCapContainsAtom(UNIT uUnit, ATOM aAtom);
extern BOOL     bUnitCapContainsContainer(UNIT uUnit, CONTAINER cCont);

extern BOOL     bUnitGroupCreate(UNIT uUnit, char *cPName);
extern LIST     lUnitGroup(UNIT uUnit, char *sGroup);
extern BOOL     bUnitGroupAddAtom(UNIT uUnit, char *sGroup, ATOM aAtom);
extern BOOL     bUnitGroupFindAtom(UNIT uUnit, char *sGroup, ATOM aAtom, 
                        BOOL *bPFound);
extern BOOL     bUnitGroupRemoveAtom(UNIT uUnit, char *sGroup, ATOM aAtom);
extern BOOL     bUnitGroupDestroy(UNIT uUnit, char *sGroup);
extern void     UnitFindBoundingBox(UNIT uUnit, VECTOR *vPLower, 
                        VECTOR *vPUpper);

extern void     UnitSetUseBox(UNIT uUnit, BOOL b);
extern void     UnitSetBoxOct(UNIT uUnit, BOOL b);
extern void     UnitSetUseSolventCap(UNIT uUnit, BOOL b);
extern void     UnitDestroy( UNIT *uPUnit );
extern BOOL     zbUnitIgnoreHwHwOwAngle( STRING sA, STRING sB, STRING sC );
extern BOOL     zbUnitIgnoreAngle( STRING sA, STRING sB, STRING sC );

#define UnitUseParameters( u, p ) \
        ( ((UNIT)(u))->psParameters = (OBJEKT)(p),CDU(u) )
#define psUnitParameters( u )     ((PARMSET)( ((UNIT)(u))->psParameters ))

#define bUnitHeadUsed(u)                ( ((UNIT)(u))->aHead != NULL )
#define aUnitHead(u)                    ( (ATOM)(((UNIT)(u))->aHead) )
#define UnitSetHead( u, a )             ( ((UNIT)(u))->aHead = (OBJEKT)(a),CDU(u) )

#define bUnitTailUsed(u)        ( ((UNIT)(u))->aTail != NULL )
#define aUnitTail(u)            ( (ATOM)(((UNIT)(u))->aTail) )
#define UnitSetTail( u, a )     ( ((UNIT)(u))->aTail = (OBJEKT)(a),CDU(u) )

#define UnitSetMode( u, m )             (((UNIT)(u))->iMode = m,CDU(u) )
#define iUnitMode( u )                  (((UNIT)(u))->iMode)
#define dUnitAtomGroups(u)              (((UNIT)(u))->dAtomGroups)
#define UnitGetBox( u, xP, yP, zP ) (\
        *(xP) = ((UNIT)(u))->dXWidth,\
        *(yP) = ((UNIT)(u))->dYWidth,\
        *(zP) = ((UNIT)(u))->dZWidth )
#define UnitSetBox( u, x, y, z ) (\
        ((UNIT)(u))->dXWidth = (x),\
        ((UNIT)(u))->dYWidth = (y),\
        ((UNIT)(u))->dZWidth = (z),CDU(u) )
#define UnitGetSolventCap( u, xP, yP, zP, rP ) {\
        (*xP) = dVX(&(u->vCapOrigin));\
        (*yP) = dVY(&(u->vCapOrigin));\
        (*zP) = dVZ(&(u->vCapOrigin));\
        (*rP) = u->dCapRadius;}
#define UnitSetSolventCap( u, x, y, z, r ) {\
        VectorDef( &(u->vCapOrigin), x, y, z );\
        u->dCapRadius = r;CDU(u); }
#define UnitSetBeta( u, d )     ( ((UNIT)(u))->dBeta = d,CDU(u) )
#define dUnitBeta( u )          ( ((UNIT)(u))->dBeta )
#define UnitDefineFlags( u, f ) ( ((UNIT)(u))->fFlags = (f), CDU(u) )
#define UnitSetFlags( u, f )    { ((UNIT)(u))->fFlags |= (f); CDU(u); }
#define UnitResetFlags( u, f )  { ((UNIT)(u))->fFlags &= ~(f); CDU(u); }
#define bUnitFlagsSet( u, f )   ( (((UNIT)(u))->fFlags & (f)) == (f) )
#define bUnitUseBox(u)          (bUnitFlagsSet(u,UNITUSEBOUNDINGBOX))
#define bUnitBoxOct(u)          (bUnitFlagsSet(u,UNITBOXOCT))
#define bUnitUseSolventCap(u)           bUnitFlagsSet(u,UNITUSESOLVENTCAP)


#endif /* UNIT_H */
