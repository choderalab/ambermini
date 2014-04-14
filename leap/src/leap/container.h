/*
 *      File: container.h
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
 *              CONTAINER
 *      Superclass: 
 *              OBJEKT, LOOP
 *      SubClasses:
 *              MOLSYSTEM, MOLECULE, RESIDUE, ATOM
 *
 *      Description:
 *
 *              CONTAINER is a superclass for all objects used to
 *              define a system of molecules.  The properties that
 *              it defines are:
 *              Every CONTAINER can be contained by another.
 *              Every CONTAINER can have contents.
 *              Every CONTAINER has a name.
 *              Looping over CONTAINER contents is handled by
 *              the LOOP class.
 */
 
#ifndef CONTAINER_H
#define CONTAINER_H


# include       "displayer.h"
# include       "matrix.h"
# include       "stringExtra.h"


/*
 *-----------------------------------------------------------------------
 *
 *        Define object typedefs here.
 *        
 *        Object typedef MUST include the superclass object as
 *        its first structure element.
 */



#define CONTAINERNAMELEN        32
typedef char    CONTAINERNAMEt[CONTAINERNAMELEN];

typedef struct  CONTAINERSTRUCT {
        OBJEKTt                 oHeader;
        struct CONTAINERSTRUCT  *cCopy;
        CONTAINERNAMEt          sName;
        struct CONTAINERSTRUCT  *cContainedBy;
        int                     iNextChildsSequence;
        DISPLAYER               dDisp;
        int                     iTempInt;
        int                     iSequence;
        struct CONTAINERSTRUCT  *cLoopNext;     /* Used for LOOPing */
        OBJEKT                  lContents;
} CONTAINERt;

typedef CONTAINERt      *CONTAINER;





/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primitive superclass
        of this class that it will eventually make it into these routines.
*/


#define CDU(c)  ContainerDisplayerUpdate((CONTAINER)c)

#define sContainerName( c )     ( ((CONTAINER)c)->sName )
#define ContainerSetName( c, n ) ( StringCopyMax( ((CONTAINER)c)->sName, n,\
                                        sizeof(CONTAINERNAMEt) ),CDU(c) )
#define cContainerWithin( c )   ( ((CONTAINER)c)->cContainedBy )
#define ContainerSetWithin( c, w ) (((CONTAINER)c)->cContainedBy = w,CDU(c) )
#define cContainerContents( c )    (((CONTAINER)c)->lContents)
#define iContainerNumberOfChildren(c) (iListSize(cContainerContents(c)))


                /* If the copy is NULL then return the original */
                /* otherwise return the copy                    */
#define cContainerCopyPointer( c ) \
( c == NULL ? NULL : ( ((((CONTAINER)c)->cCopy)!=NULL) ? \
(CONTAINER)(((CONTAINER)c)->cCopy) : ((CONTAINER)c) ))

#define ContainerSetCopyPointer( c, n )   (((CONTAINER)c)->cCopy=n)

#define ContainerSetTempInt(c,i)        (((CONTAINER)(c))->iTempInt = i )
#define iContainerTempInt(c)            (((CONTAINER)(c))->iTempInt)
#define iContainerSequence( c )         (((CONTAINER)c)->iSequence)
#define ContainerSetSequence( c,xn )    (((CONTAINER)c)->iSequence = xn,CDU(c))
#define iContainerNextChildsSequenceInc( c ) \
        (((CONTAINER)(c))->iNextChildsSequence++)
#define iContainerNextChildsSequence( c ) \
        (((CONTAINER)(c))->iNextChildsSequence)
#define ContainerSetNextChildsSequence( c, i ) \
        (((CONTAINER)(c))->iNextChildsSequence = i,CDU(c))


#define ContainerSetLoopNext( c, n )    ( ((CONTAINER)(c))->cLoopNext = n )
#define cContainerLoopNext(c)           ( ((CONTAINER)(c))->cLoopNext )

#define dContainerDisplayer(c)  (((CONTAINER)c)->dDisp)


extern CONTAINER        cContainerCreate( int iType );
extern void             ContainerDestroy( CONTAINER *cPContainer );
extern void             ContainerDescribe( CONTAINER cContainer  ) ;
extern void             ContainerAdd( CONTAINER cContainer, OBJEKT oObject );
extern BOOL             bContainerRemove( CONTAINER cContainer, 
                                OBJEKT oObject );
extern CONTAINER        cContainerDuplicate( CONTAINER cOld );
extern void             ContainerResetPointers( CONTAINER cContainer );
extern CONTAINER        cContainerFindSequence( CONTAINER cCont, 
                                int iContainerType, int iSeq );
extern CONTAINER        cContainerFindName( CONTAINER cCont, 
                                int iContainerType, char *sName );
extern VECTOR           vContainerGeometricCenter( CONTAINER cCont );
extern void             ContainerBoundingBox( CONTAINER cCont, 
                                VECTOR *vPMin, VECTOR *vPMax );
extern void             ContainerCenterAt( CONTAINER cCont, VECTOR vCenter );
extern void             ContainerTransformBy( CONTAINER cCont, 
                                MATRIX mTransform );
extern void             ContainerTranslateBy( CONTAINER cCont, VECTOR vOffset );
extern void             ContainerSetAllAtomsFlags( CONTAINER cCont, 
                                FLAGS fFlags );
extern void             ContainerResetAllAtomsFlags( CONTAINER cCont, 
                                FLAGS fFlags );
extern void             ContainerWithFlagsSetAtomFlags( CONTAINER cCont, 
                                FLAGS fNeed, FLAGS fFlags );
extern void             ContainerWithFlagsResetAtomFlags( CONTAINER cCont, 
                                FLAGS fNeed, FLAGS fFlags );
extern void             ContainerWithoutFlagsSetAtomFlags( CONTAINER cCont, 
                                FLAGS fNeed, FLAGS fFlags );
extern void             ContainerWithoutFlagsResetAtomFlags( CONTAINER cCont,
                                FLAGS fNeed, FLAGS fFlags );
extern char             *sContainerDescriptor( CONTAINER cCont, char *sDesc );
extern char             *sContainerFullDescriptor( CONTAINER cCont, 
                                char *sFullDesc );
extern void             ContainerCheck( CONTAINER cCont, int *iPErrors, 
                                int *iPWarnings );
extern BOOL             bContainerContainedBy( CONTAINER cIn, CONTAINER cOut );
extern void             ContainerYouAreBeingRemoved( CONTAINER cCont ) ;
extern void             ContainerIAmBeingRemoved( CONTAINER cCont, 
                                CONTAINER cRemoved );
extern void             ContainerSetAttribute( CONTAINER cCont, 
                                STRING sAttribute, OBJEKT oValue );
extern void             ContainerTotalCharge( CONTAINER cCont, 
                                double *dPCharge, double *dPPertCharge );
extern void             ContainerDisplayerUpdate( CONTAINER cCont );
extern void             ContainerTreeMakeInsensitive( CONTAINER cCont );
extern void             ContainerTreeMakeSensitive( CONTAINER cCont );
extern void             ContainerResetAllCopyPointers( CONTAINER cTop );
extern BOOL             bContainerSpaceConflict( CONTAINER cCont1, 
                                CONTAINER cCont2 );


#endif /* CONTAINER_H */
