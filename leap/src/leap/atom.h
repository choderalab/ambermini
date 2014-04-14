#ifndef ATOM_H
#define ATOM_H

/*
 *      File: atom.h
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
 *             David A. Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Class: 
 *              ATOM
 *      Superclass: 
 *              CONTAINER
 *
 *      Description:
 *
 *              Atoms can contain anything.
 */
 
#include        "vector.h"

#define ATOM_DEFAULT_RADIUS     1.5

                        /* 5 characters to represent the type */
#define ATOMTYPELEN     5
                        /* 10 characters to represent the name */
#define NAMELEN         10
                        /* maximum 8 bonds out of each atom */
#define MAXBONDS        8


/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/


#include        "elements.h"

                /* Hybridization type */

#define HUNDEFINED      -1
#define HUNKNOWN        0
#define HSP1            1
#define HSP2            2
#define HSP3            3


                /* Bond flags */
#define BONDORDERONLY   0x0000000F
#define BONDNONE        0x00000000
#define BONDSINGLE      0x00000001
#define BONDDOUBLE      0x00000002
#define BONDTRIPLE      0x00000003
#define BONDAROMATIC    0x00000004

#define BONDTEMPORARY   0x00000010


#define ATOM_SEGID_LEN  5               /* 4 char + \0 */
typedef char            SEGIDt[ATOM_SEGID_LEN];
typedef char            BONDt;
typedef char            ATOMTYPEt[ATOMTYPELEN];

typedef struct  ATOMSTRUCT {
        CONTAINERt              cHeader;
        int                     iUniqueId;
        int                     iAtomicNumber;
        int                     iPertAtomicNumber;
        CONTAINERNAMEt          sPertName;
        ATOMTYPEt               sType;
        ATOMTYPEt               sPertType;
        double                  dCharge;
        double                  dPertCharge;
        double                  dPolar;
        double                  dPertPolar;
	double			dScreenF;
        int                     iIndex;
        VECTOR                  vPosition;
        VECTOR                  vVelocity;
        FLAGS                   fFlags;         /* Atom flags */
        int                     iCoordination;
        struct ATOMSTRUCT       *aaBonds[MAXBONDS];
        FLAGS                   faBondFlags[MAXBONDS];
        double                  dTemp;
        GENP                    PTemp;
        SEGIDt                  siSegid;
                /* Spanning tree stuff */
        struct ATOMSTRUCT       *aNextSpan;
        struct ATOMSTRUCT       *aBackSpan;
        int                     iSeenId;
        int                     iBackCount;
                /* Graphics stuff */
        GENP                    PGraphicsData;
} ATOMt;

typedef ATOMt   *ATOM;


                /* Atom flags */

#define ATOMALLFLAGS            0xFFFFFFFF

#define ATOMPERMANENTFLAGS      0x0000FFFF

                                /* Differentiate between the atoms position */
                                /* being FIXED or being BUILT               */
                                /* Depending on whether the atoms position  */
                                /* was read in from PDB or built using the  */
                                /* builder, Also define whether the position*/
                                /* is known or not.  This makes visibility */
                                /* testing in spanning tree loops easy      */
#define ATOMSELECTED            0x00000001
                                /* This flag must be set if the ATOM is to */
                                /* be perturbed */

#define ATOMPERTURB             0x00000002
#define ATOMNOTDISPLAYED        0x00000004
#define RESIDUEIMAGEATOM        0x00000010      /* not used (so far) */

#define ATOMTEMPORARYFLAGS      0xFFFF0000
        
#define ATOMTOUCHED             0x00010000
#define ATOMPOSITIONKNOWN       0x00020000      /* ATOMPOSITIONEXTERNAL */
#define ATOMPOSITIONINTERNAL    0x00040000
#define ATOMNEEDSMINIMIZER      0x00080000
#define ATOMNEEDSBUILD          0x00100000
#define ATOMTOUCHED2            0x00200000
#define ATOMPOSITIONFIXED       0x00400000
#define ATOMPOSITIONBUILT       0x00800000
#define ATOMPOSITIONDRAWN       0x01000000
#define ATOMMOCKPOSITIONKNOWN   0x02000000
#define ATOMFLAGUNUSED0         0x04000000
#define ATOMFLAGUNUSED1         0x08000000
#define ATOMFLAG0               0x10000000
#define ATOMFLAG1               0x20000000
#define ATOMFLAG2               0x40000000
#define ATOMFLAG3               0x80000000




/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primitive superclass
        of this class that it will eventually make it into these routines.
*/



#define AtomDefineBondFlags(a,i,f)      (((ATOM)(a))->faBondFlags[i] = f,CDU(a) )
#define AtomSetBondFlags(a,i,f)         (((ATOM)(a))->faBondFlags[i] |= f,CDU(a) )
#define AtomResetBondFlags(a,i,f)       \
                (((ATOM)(a))->faBondFlags[i] &= ~(f),CDU(a) )
#define AtomSetBondOrder(a,i,f) ( AtomResetBondFlags(a,i,BONDORDERONLY),\
                                  AtomSetBondFlags(a,i,f),CDU(a) )
#define iAtomBondOrder(a,i)     ((((ATOM)(a))->faBondFlags[i])&(BONDORDERONLY))
#define fAtomBondFlags(a,i)     (((ATOM)(a))->faBondFlags[i])
#define bAtomBondFlagsSet(a,i,f)        ((((ATOM)(a))->faBondFlags & f)!=0)
#define aAtomBondedNeighbor(a,i)        (((ATOM)(a))->aaBonds[i])
#define AtomSetTempPtr( a, p )  (((ATOM)(a))->PTemp = (GENP)(p))
#define PAtomTempPtr(a)                 (((ATOM)(a))->PTemp)
#define AtomSetElement( a, n )          (((ATOM)(a))->iAtomicNumber = n,CDU(a))
#define iAtomElement(a)                 (((ATOM)(a))->iAtomicNumber)
#define AtomSetPertElement( a, n )      (((ATOM)(a))->iPertAtomicNumber = n,CDU(a))
#define iAtomPertElement(a)             (((ATOM)(a))->iPertAtomicNumber)
#define AtomDupPosition(a,vv)   {((ATOM)(a))->fFlags|=ATOMPOSITIONKNOWN;\
                                        VectorCopy(vv, vAtomPosition((ATOM)(a))); CDU(a);}
#define AtomSetPosition(a,vv)   (((ATOM)(a))->fFlags|=ATOMPOSITIONKNOWN,\
                                        vAtomPosition((ATOM)(a)) = vv, CDU(a))
#define AtomSetPositionNoFlags(a,vv)    {VectorCopy(vv, vAtomPosition((ATOM)(a))); CDU(a));}
#define vAtomPosition(a)                (((ATOM)(a))->vPosition)
#define AtomDupVelocity(a,vv)           {VectorCopy(vv, vAtomVelocity((ATOM)(a))); CDU(a);}
#define AtomSetVelocity(a,vv)           (vAtomVelocity((ATOM)(a)) = vv, CDU(a))
#define vAtomVelocity(a)                (((ATOM)(a))->vVelocity)
#define iAtomCoordination(a)            (((ATOM)(a))->iCoordination)
#define iAtomId(a)                      (((ATOM)(a))->iUniqueId)
#define AtomSetType( a, x )             (strcpy(((ATOM)(a))->sType,x),CDU(a))
#define sAtomType( a )                  (((ATOM)(a))->sType )
#define AtomSetName( a, s )             ContainerSetName( a, s )
#define sAtomName( a )                  sContainerName( a )
#define AtomSetPertName( a, s )         (strcpy(((ATOM)(a))->sPertName, s ),CDU(a))
#define sAtomPertName( a )              (((ATOM)(a))->sPertName)
#define AtomSetPertType( a, x )         (strcpy(((ATOM)(a))->sPertType,x),CDU(a))
#define sAtomPertType( a )              (((ATOM)(a))->sPertType )
#define bAtomPerturbed( a )             (bAtomFlagsSet(((ATOM)(a)),ATOMPERTURB))
#define AtomDefineFlags(a,f)            (((ATOM)(a))->fFlags = f,CDU(a) )
#define AtomSetFlags(a,f)               (((ATOM)(a))->fFlags |= f,CDU(a) )
#define AtomResetFlags(a,f)             (((ATOM)(a))->fFlags &= (~f),CDU(a) )
#define fAtomFlags(a)                   (((ATOM)(a))->fFlags)
#define bAtomFlagsSet(a,f)              ((((ATOM)(a))->fFlags&(f))==(f))
#define bAtomFlagsReset(a,f)     ((~(((ATOM)(a))->fFlags)&(f))==(f))
#define aAtomNextSpan(AA)               (ATOM)(((ATOM)(AA))->aNextSpan)
#define AtomSetNextSpan(AA,SS)          (((ATOM)(AA))->aNextSpan = SS )
#define aAtomBackSpan(AA)               (ATOM)(((ATOM)(AA))->aBackSpan) 
#define AtomSetBackSpan(AA,SS)          (((ATOM)(AA))->aBackSpan = SS )
#define iAtomBackCount(AA)              (((ATOM)(AA))->iBackCount)
#define AtomSetBackCount(AA,I)          (((ATOM)(AA))->iBackCount=(I))
#define iAtomSeenId(AA)                 (((ATOM)(AA))->iSeenId)
#define AtomSetSeenId(AA,II)            (((ATOM)(AA))->iSeenId =(II))
#define AtomSetTempInt(AA,II)           (ContainerSetTempInt(AA,II))
#define iAtomTempInt(AA)                (iContainerTempInt(AA))
#define AtomSetCharge( a, x )           (((ATOM)(a))->dCharge = x,CDU(a))
#define AtomSetPolar( a, x )            (((ATOM)(a))->dPolar = x,CDU(a))
#define dAtomCharge( a )                (((ATOM)(a))->dCharge)
#define dAtomPolar( a )                 (((ATOM)(a))->dPolar)
#define AtomSetPertCharge( a, x )       (((ATOM)(a))->dPertCharge = x,CDU(a))
#define AtomSetPertPolar( a, x )        (((ATOM)(a))->dPertPolar = x,CDU(a))
#define dAtomPertCharge( a )            (((ATOM)(a))->dPertCharge)
#define dAtomPertPolar( a )             (((ATOM)(a))->dPertPolar)
#define iAtomIndex( a )                 (((ATOM)(a))->iIndex)
#define AtomSetIndex( a, i )            (((ATOM)(a))->iIndex = i)
#define AtomSetGraphicsPointer(a,p) (((ATOM)(a))->PGraphicsData=((GENP)p))
#define PAtomGraphicsPointer(a)  (((ATOM)(a))->PGraphicsData)

#define AtomSetTempDouble(a,v)          (((ATOM)(a))->dTemp = v)
#define AtomTempDoubleIncrement(a,v)    (((ATOM)(a))->dTemp += v)
#define AtomTempDoubleSquare(a)         (((ATOM)(a))->dTemp *= \
                                                ((ATOM)(a))->dTemp)
#define AtomTempDoubleSquareRoot(a)     (((ATOM)(a))->dTemp = \
                                                sqrt(((ATOM)(a))->dTemp) )
#define dAtomTemp(a)                    (((ATOM)(a))->dTemp)

#define bAtomVisible(a) ((bAtomFlagsSet(a,ATOMPOSITIONKNOWN)||\
                          bAtomFlagsSet(a,ATOMPOSITIONDRAWN))&&\
                          !bAtomFlagsSet(a,ATOMNOTDISPLAYED))
#define bAtomHasPosition(a)     (bAtomFlagsSet(a,ATOMPOSITIONKNOWN)||\
                                bAtomFlagsSet(a,ATOMPOSITIONDRAWN))

#define sAtomSegid(a)           ((a)->siSegid)
#define AtomSetSegid(a,s)       (strncpy((a)->siSegid,(s),ATOM_SEGID_LEN-1), \
                                (a)->siSegid[ATOM_SEGID_LEN-1]='\0')


/*  atom.c  */

extern ATOM             aAtomCreate();
extern void             AtomDestroy( ATOM *aPAtom );
extern void             AtomDescribe( ATOM aAtom );
extern void             AtomDescStr( ATOM aA, BOOL bResNum, char *cPDesc );
extern BOOL             bAtomCoordinationSaturated( ATOM aAtom );
extern BOOL             AtomTmpBondTo( ATOM aAtom1, ATOM aAtom2 );
extern void             AtomBondToOrder( ATOM aAtom1, ATOM aAtom2, int iOrder );
extern void             AtomBondToFlags( ATOM aAtom1, ATOM aAtom2, 
                                FLAGS fFlags );
extern void             AtomRemoveBond( ATOM aAtom1, ATOM aAtom2 );
extern BOOL             bAtomBondedTo( ATOM aAtom1, ATOM aAtom2 );
extern int              iAtomFindBondOrder( ATOM aAtom, ATOM aNeighbor );
extern void             AtomFindSetBondOrder( ATOM aAtom, ATOM aNeighbor, 
                                int iOrder );
extern FLAGS            fAtomFindBondFlags( ATOM aAtom, ATOM aNeighbor );
extern ATOM             aAtomDuplicate( ATOM aAtom );
extern void             AtomResetPointers( ATOM aAtom );
extern void             AtomCheck( ATOM aAtom, int *iPErrors, int *iPWarnings );
extern void             AtomYouAreBeingRemoved( ATOM aAtom );
extern void             AtomIAmBeingRemoved( ATOM aAtom, CONTAINER cRemoved );
extern void             AtomSetAttribute( ATOM aAtom, STRING sAttr, 
                                OBJEKT oAttr );
extern int              iAtomHybridization( ATOM aAtom );
extern int              iAtomBondOrderFromName( char *sName );
extern BOOL             bAtomSpaceConflict( ATOM aAtom1, ATOM aAtom2 );
extern double           dAtomVanderWaals( ATOM aAtom );
extern int              iAtomSetTmpRadius( ATOM aAtom );

extern void             AtomBondTo( ATOM aAtom1, ATOM aAtom2 );

#endif /* ATOM_H */
