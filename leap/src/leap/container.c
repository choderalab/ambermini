/*
 *      File: container.c
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
 *              CONTAINER
 *      Superclass:
 *              OBJEKT
 *
 *      Description:
 *              CONTAINER is a superclass for all objects used to
 *              define a system of molecules.  The properties that
 *              it defines are:
 *              Every CONTAINER can be contained within another CONTAINER.
 *              Every CONTAINER can contain other containers, or objects.
 *              Every CONTAINER has a name.
 *              Looping over CONTAINER contents is handled by
 *              the LOOP object.
 *              Containers have a sequence number that can be
 *              used to search for a particular type of
 *              container within another.
 *
 *              When a CONTAINER is being removed from another CONTAINER
 *              the CONTAINER that is being removed is informed using
 *              ContainerYouAreBeingRemoved and it then informs the
 *              CONTAINER that contains it that it is being removed
 *              using ContainerIAmBeingRemoved, a message that is passed
 *              all the way up the chain.
 *              This is to garentee that the CONTAINERs that contain
 *              the CONTAINER being removed can NULL any pointers
 *              to that CONTAINER, leaving no loose pointers.
 *
 *              CONTAINERs use DISPLAYERs to update any windows
 *              that may be displaying themselves and their contents.
 *
 *              Whenever a CONTAINER is modified in a way that would
 *              change the way it is being displayed call
 *              ContainerDisplayerUpdate or use the macro CDU
 *
 */
 



#include        "basics.h"

#include        "vector.h"

#include        "classes.h"

#include        "container.h"


#define NOSEQUENCE      -1



/*
 *      cContainerCreate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to create objects of
 *      a subclass of CONTAINER.
 *
 *      Return:
 *              Return a pointer to the object created.
 */
CONTAINER 
cContainerCreate( int iType )
{
CONTAINER       c;

    c = NULL;
    switch ( iType ) {
        case CONTAINERid:
            DFATAL( ("Illegal to create raw CONTAINERs") );
            break;
        case UNITid:
            c = (CONTAINER)uUnitCreate();
            break;
        case MOLECULEid:
            c = (CONTAINER)mMoleculeCreate();
            break;
        case RESIDUEid:
            c = (CONTAINER)rResidueCreate();
            break;
        case ATOMid:
            c = (CONTAINER)aAtomCreate();
            break;
        default:
            DFATAL( ("ERROR, can't create unknown CONTAINER id:%d", iType) );
            break;
    }

        /* Initialize the CONTAINER itself */
    c->lContents = oCreate(LISTid);
    strcpy( c->sName, "" );
    c->cContainedBy = NULL;
    c->dDisp = dDisplayerCreate((GENP)c);
    ContainerSetSequence(c,NOSEQUENCE);

                /* Set the next sequence number given to a Child */
                
    c->iNextChildsSequence = 1;
                                        
    
    return(c);
}







/*
 *      ContainerDestroy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to destroy CONTAINERS.
 *      After the call to Destroy, the object will be destroyed and
 *      the value pointed to by oObject will be undefined.
 *
 *      Arguments:
 *      oObject -       A pointer to the object.
 */
void    
ContainerDestroy( CONTAINER *cPContainer )
{
LISTLOOP        llLoop;
OBJEKT          oObj;

    if ( *cPContainer == NULL ) 
        return;
 
    if ((*cPContainer)->dDisp != NULL)  
        DisplayerDestroy(&((*cPContainer)->dDisp));
 
                /* DEREF the contents of the CONTAINER */
                
    llLoop = llListLoop((LIST)(*cPContainer)->lContents);
    while ( (oObj=oListNext(&llLoop))!= NULL ) DEREF(oObj);
    ListDestroy((LIST *)&((*cPContainer)->lContents));
    
    switch ( iObjectType(*cPContainer) ) {
        case CONTAINERid:
            DFATAL( ("Illegal to destroy raw CONTAINERs") );
            break;
        case UNITid:
            UnitDestroy((UNIT *)cPContainer);
            break;
        case MOLECULEid:
            MoleculeDestroy((MOLECULE *)cPContainer);
            break;
        case RESIDUEid:
            ResidueDestroy((RESIDUE *)cPContainer);
            break;
        case ATOMid:
            AtomDestroy((ATOM *)cPContainer);
            break;
        default:
            DFATAL( ("Attempting to destroy CONTAINER with id:%c", 
                        iObjectType(*cPContainer)) );
            break;
    }
}






/*
 *      ContainerDescribe
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to return a descriptor
 *      string for the CONTAINER object.
 *      It is the callers responsibility to ensure that there is room
 *      enough in the buffer to describe the object.
 *
 *      Arguments:
 *      oObject -       The object to describe.
 *
 */
void    
ContainerDescribe( CONTAINER cContainer  )
{
    switch ( iObjectType(cContainer) ) {
        case CONTAINERid:
            DFATAL( ("Attempt to describe a raw CONTAINER") );
            break;
        case UNITid:
            UnitDescribe( (UNIT)cContainer );
            break;
        case MOLECULEid:
            MoleculeDescribe( (MOLECULE)cContainer );
            break;
        case RESIDUEid:
            ResidueDescribe( (RESIDUE)cContainer );
            break;
        case ATOMid:
            AtomDescribe( (ATOM)cContainer );
            break;
        default:
            DFATAL( ("Attempting to describe CONTAINER with id:%d", 
                        iObjectType(cContainer)) );
            break;
    }
}




/*
 *      ContainerAdd
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add an object to this container.
 *      This routine has the side effect of updating the sequence
 *      number of oObject.
 *
 */
void    
ContainerAdd( CONTAINER cContainer, OBJEKT oObject )
{
    TESTMEMORY();

                /* Add the OBJEKT to the end of the list to perserve */
                /* the NATURAL ORDER OF THINGS */
    ListAddToEnd( (LIST)cContainer->lContents, oObject );
        
    TESTMEMORY();
                /* If the object being placed within this container */
                /* is itself a container then record within the object */
                /* which container it is being placed within. */
    if ( bObjectInClass( oObject, CONTAINERid ) ) {
        ContainerSetSequence( oObject, 
                iContainerNextChildsSequenceInc(cContainer) );
        ContainerSetWithin( oObject, cContainer );
    }
    
    TESTMEMORY();

    CDU(cContainer);
}







/*
 *      bContainerRemove
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove an object from this container or it's subcontainers.
 *      Return FALSE if the OBJEKT could not be found
 *      within the CONTAINER. 
 *
 *      UNTIL A BETTER FIX IS DEVISED: for safety, always REF() the
 *      thing being removed before calling this routine, since it
 *      assumes the thing survives the intitial list removal, which
 *      may not be the case, which leads to accessing of freed space
 *      and memory corruption.
 */
BOOL    
bContainerRemove( CONTAINER cContainer, OBJEKT oObject )
{
LISTLOOP        llList;
CONTAINER       cSub;
BOOL            bResult;

        /* First try and remove it directly from this container */
        
    if ( bListRemove( (LIST)(cContainer->lContents), (GENP)oObject )) {
        /*
         *  here's where oObject had better still exist!
         */
        if ( bObjectInClass( oObject, CONTAINERid ) ) {
            ContainerYouAreBeingRemoved( (CONTAINER)oObject );
            ContainerSetWithin( oObject, NULL );
        }
        bResult = TRUE;
    } else {
        /* Now try and remove it from a sub-container */
    
        llList = llListLoop( (LIST)cContainer->lContents );
        while ( (cSub = (CONTAINER)oListNext(&llList)) != NULL ) {
            if ( bContainerRemove( cSub, oObject ) ) return(TRUE);
        }
    
        /* None of the sub-containers contained the object */
        /* so return FALSE                                 */

        bResult = FALSE;
    }

    CDU(cContainer);
    return(bResult);
}







/*
 *      cContainerDuplicate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Send the message to the appropriate subclass to return
 *      a duplicate. The <subclasse>Duplicate() routines should
 *      only be called from here, so that the container-specific
 *      properties can be updated here. By the same token, this
 *      routine should only be called by oObjectDuplicate(), so
 *      that it can update the objekt-specific properties.
 */
CONTAINER       
cContainerDuplicate( CONTAINER cOld )
{
CONTAINER       cNew;
OBJEKT          oContents;
        
        /* Duplicate the specific container information */

    cNew = NULL;        
    switch ( iObjectType(cOld) ) {
        case UNITid:
            cNew = (CONTAINER)uUnitDuplicate((UNIT)cOld);
            break;
        case MOLECULEid:
            cNew = (CONTAINER)mMoleculeDuplicate((MOLECULE)cOld);
            break;
        case RESIDUEid:
            cNew = (CONTAINER)rResidueDuplicate((RESIDUE)cOld);
            break;
        case ATOMid:
            cNew = (CONTAINER)aAtomDuplicate((ATOM)cOld);
            break;
    }
   
        /* The individual duplicate routines should have duplicated     */
        /* the entire object                                            */
    
        /* Duplicate the contents */
    oContents = oObjectDuplicate( cContainerContents(cOld) );
    cNew->lContents = oContents;
    
        /* Set redirection pointer for fixing up internal pointers */
        /* later */
        
    ContainerSetCopyPointer( cNew, NULL );
    ContainerSetCopyPointer( cOld, cNew );
    
        /* Create a new displayer for the new container. */
        /* The displayer is not to be copied */

    cNew->dDisp = dDisplayerCreate((GENP)cNew);
 
    return(cNew);
}



/*
 *      ContainerResetPointers
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Reset all internal pointers.
 *      When copies are made of CONTAINERs, internal pointers, such
 *      as terminis, are left pointing to the old CONTAINERs.
 *      However, the CONTAINER class supports an extra pointer
 *      on the old CONTAINER which points to the new CONTAINER.
 *      This can be used to reset the internal pointers of the RESIDUE.
 *      The value of the pointer in new CONTAINERs is NULL.
 *
 *      EVERY internal pointer of the CONTAINER has to be resolved.
 */
void    
ContainerResetPointers( CONTAINER cContainer )
{
    cContainer->cContainedBy = 
        cContainerCopyPointer(cContainer->cContainedBy);

    switch ( iObjectType(cContainer) ) {
        case UNITid:
            UnitResetPointers((UNIT)cContainer);
            break;
        case MOLECULEid:
            MoleculeResetPointers((MOLECULE)cContainer);
            break;
        case RESIDUEid:
            ResidueResetPointers((RESIDUE)cContainer);
            break;
        case ATOMid:
            AtomResetPointers((ATOM)cContainer);
            break;
        default:
            DFATAL( ("Can not ResetPointers for unknown class id: %d",
                        iObjectType(cContainer) ) );
            break;
    }
}






/*
 *      cContainerFindSequence
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Search the CONTAINER for another CONTAINER of iContainerType
 *      with the sequence number iSeq.
 *      If it is found, return it, otherwise return NULL.
 */
CONTAINER       
cContainerFindSequence( CONTAINER cCont, int iContainerType, int iSeq )
{
LOOP    lContainers;
OBJEKT  oObj;

    if ( iSeq > iContainerNextChildsSequence(cCont) ) {
        return(NULL);
    }
    lContainers = lLoop( (OBJEKT)cCont, iContainerType );
    while ( (oObj = oNext(&lContainers) ) != NULL ) {
        if ( iContainerSequence(oObj) == iSeq ) break;
    }
    return((CONTAINER)oObj);
}




/*
 *      cContainerFindName
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Search the CONTAINER for another CONTAINER of iContainerType
 *      with the name sName.
 *      If it is found, return it, otherwise return NULL.
 */
CONTAINER       
cContainerFindName( CONTAINER cCont, int iContainerType, char *sName )
{
LOOP    lContainers;
OBJEKT  oObj;

    lContainers = lLoop( (OBJEKT)cCont, iContainerType );
    while ( (oObj = oNext(&lContainers) ) != NULL ) {
        if ( strcmp( sContainerName(oObj), sName) == 0 ) break;
    }
    return((CONTAINER)oObj);
}





/*
 *      vContainerGeometricCenter
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all the ATOMs of the CONTAINER and return the
 *      Geometric center of all the vectors.
 */
VECTOR  
vContainerGeometricCenter( CONTAINER cCont )
{
VECTOR  vSum;
int     iVectorCount;
LOOP    lAtoms;
OBJEKT  aAtom;

    iVectorCount = 0;
    VectorDef( &vSum, 0.0, 0.0, 0.0 );
    if ( iObjectType(cCont) == ATOMid ) {
        vSum = vAtomPosition(cCont);
    } else {
        lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
        while ( (aAtom = oNext(&lAtoms)) != NULL ) {
            if ( bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) ) {
                vSum = vVectorAdd( &vSum, &vAtomPosition(aAtom) );
                iVectorCount++;
            }
        }
        if ( iVectorCount != 0 ) {
            vSum = vVectorTimesScalar( &vSum, 
                                (double)(1.0/(double)iVectorCount) );
        }
    }
    return(vSum);
}





/*
 *      ContainerBoundingBox
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Modified: David A. Rivkin (1-9-93)
 *
 *      Calculate bounding box for the ATOMs within the CONTAINER.
 *      Return the min X,Y,Z in vPMin and max X,Y,Z in vPMax.
 */
void    
ContainerBoundingBox( CONTAINER cCont, VECTOR *vPMin, VECTOR *vPMax )
{
double          dXMin, dYMin, dZMin;
double          dXMax, dYMax, dZMax;
LOOP            lAtoms;
ATOM            aAtom;

    dXMin = dYMin = dZMin =  9.0E99;
    dXMax = dYMax = dZMax = -9.0E99;

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
            while ( (aAtom = (ATOM)oNext( &lAtoms )) ) {
                if ( dXMin > dVX( &(vAtomPosition( aAtom )))) {
                    dXMin = dVX( &(vAtomPosition( aAtom )));
                    MESSAGE(( "Min X Value:  %4.2lf\n", dXMin ));
                }
                if ( dYMin > dVY( &(vAtomPosition( aAtom )))) {
                    dYMin = dVY( &(vAtomPosition( aAtom )));
                    MESSAGE(( "Min Y Value:  %4.2lf\n", dYMin ));
                }
                if ( dZMin > dVZ( &(vAtomPosition( aAtom )))) {
                    dZMin = dVZ( &(vAtomPosition( aAtom )));
                    MESSAGE(( "Min Z Value:  %4.2lf\n", dZMin ));
                }
                if ( dXMax < dVX( &(vAtomPosition( aAtom )))) {
                    dXMax = dVX( &(vAtomPosition( aAtom )));
                    MESSAGE(( "Max X Value:  %4.2lf\n", dXMax ));
                }
                if ( dYMax < dVY( &(vAtomPosition( aAtom )))) {
                    dYMax = dVY( &(vAtomPosition( aAtom )));
                    MESSAGE(( "Max Y Value:  %4.2lf\n", dYMax ));
                }
                if ( dZMax < dVZ( &(vAtomPosition( aAtom )))) {
                    dZMax = dVZ( &(vAtomPosition( aAtom )));
                    MESSAGE(( "Max Z Value:  %4.2lf\n", dZMax ));
                }
            }
    VectorDef( vPMin, dXMin, dYMin, dZMin );
    VectorDef( vPMax, dXMax, dYMax, dZMax );

}



/*
 *      ContainerCenterAt
 *
 *      Author: Christian Schafmeister (1991)
 *      Modified: David A. Rivkin (1-15-93)  Fixed vectors.
 *
 *      Loop through all the ATOMs of the CONTAINER and translate 
 *      all of the atoms so that the geometric center of the UNIT
 *      lies at the position vCenter.
 */
void    
ContainerCenterAt( CONTAINER cCont, VECTOR vCenter )
{
VECTOR  vOffset, vPos;
VECTOR  vCur;
LOOP    lAtoms;
OBJEKT  aAtom;

    vCur = vContainerGeometricCenter( cCont ) ;
    vOffset = vVectorSub( &vCenter, &vCur );
    MESSAGE(( "Containers' Geometric Center: (%4.2lf, %4.2lf, %4.2lf)\n",
        dVX(&vCur), dVY(&vCur), dVZ(&vCur)));
        
    if ( iObjectType(cCont) == ATOMid ) {
        AtomSetPosition( cCont, vVectorAdd( &vOffset, &vAtomPosition(cCont) ));
    } else {
        lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
        while ( (aAtom = oNext(&lAtoms)) != NULL ) {
            if ( bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN )) {
                vPos = vVectorAdd( &vOffset, &vAtomPosition(aAtom));
                MESSAGE(( "Atoms old Z:  %4.2lf\n", dVX(&vAtomPosition(aAtom)) ));
                MESSAGE(( "Vectors' Z:  %4.2lf\n", dVX(&vPos) ));
                AtomSetPosition( aAtom, vPos);
                MESSAGE(( "Atoms new Z:  %4.2lf\n", dVX(&vAtomPosition(aAtom)) ));
            } else
                VP0(( "ContainerCenterAt(): %s: %s:%s\n",
                        "Skipping atom w/ unknown position",
                        sContainerName(aAtom), sAtomName(aAtom) ));
        }
    }
}






/*
 *      ContainerTransformBy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all the ATOMs of the CONTAINER and transform 
 *      all of the atoms by the matrix.
 */
void    
ContainerTransformBy( CONTAINER cCont, MATRIX mTransform )
{
VECTOR          vPos;
LOOP            lAtoms;
OBJEKT          aAtom;

    if ( iObjectType(cCont) == ATOMid ) {
        MatrixTimesVector( vPos, mTransform, vAtomPosition(cCont) );
        AtomSetPosition( cCont, vPos );
    } else {
        lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
        while ( (aAtom = oNext(&lAtoms)) != NULL ) {
            MatrixTimesVector( vPos, mTransform, vAtomPosition(aAtom) );
            AtomSetPosition( aAtom, vPos );
        }
    }
}



/*
 *      ContainerTranslateBy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all the ATOMs of the CONTAINER and translate 
 *      all of the atoms by the amount vOffset.
 */
void    
ContainerTranslateBy( CONTAINER cCont, VECTOR vOffset )
{
VECTOR          vPos;
LOOP            lAtoms;
OBJEKT          aAtom;

    if ( iObjectType(cCont) == ATOMid ) {
        AtomSetPosition( cCont, vVectorAdd( &vOffset, &vAtomPosition(cCont) ));
    } else {
        lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
        while ( (aAtom = oNext(&lAtoms)) != NULL ) {
            if ( bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) ) {
                vPos = vVectorAdd( &vOffset, &vAtomPosition(aAtom) );
                AtomSetPosition( aAtom, vPos );
            }
        }
    }
}





/*
 *      ContainerSetAllAtomsFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all ATOMs in the CONTAINER and set
 *      the flags desired.
 */
void    
ContainerSetAllAtomsFlags( CONTAINER cCont, FLAGS fFlags )
{
LOOP            lAtoms;
ATOM            aAtom;

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) != NULL ) 
                AtomSetFlags( aAtom, fFlags );
}






/*
 *      ContainerResetAllAtomsFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all ATOMs in the CONTAINER and reset
 *      the flags desired.
 */
void    
ContainerResetAllAtomsFlags( CONTAINER cCont, FLAGS fFlags )
{
LOOP            lAtoms;
ATOM            aAtom;

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) 
        AtomResetFlags( aAtom, fFlags );
}







/*
 *      ContainerWithFlagsSetAtomFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all ATOMs in the CONTAINER and set
 *      the flags desired for all ATOMs that have fNeed set.
 */
void    
ContainerWithFlagsSetAtomFlags( CONTAINER cCont, FLAGS fNeed, FLAGS fFlags )
{
LOOP            lAtoms;
ATOM            aAtom;

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) != NULL ) {
        if ( !bAtomFlagsSet( aAtom, fNeed ) ) continue;
        AtomSetFlags( aAtom, fFlags );
    }
}




/*
 *      ContainerWithFlagsResetAtomFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all ATOMs in the CONTAINER and set
 *      the flags desired for all ATOMs that have fNeed set.
 */
void    
ContainerWithFlagsResetAtomFlags( CONTAINER cCont, FLAGS fNeed, FLAGS fFlags )
{
LOOP            lAtoms;
ATOM            aAtom;

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) != NULL ) {
        if ( bAtomFlagsSet( aAtom, fNeed ) ) 
                AtomResetFlags( aAtom, fFlags );
    }
}




/*
 *      ContainerWithoutFlagsSetAtomFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all ATOMs in the CONTAINER and set
 *      the flags desired for all ATOMs that have fNeed reset.
 */
void    
ContainerWithoutFlagsSetAtomFlags( CONTAINER cCont, FLAGS fNeed, FLAGS fFlags )
{
LOOP            lAtoms;
ATOM            aAtom;

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) != NULL ) {
        if ( bAtomFlagsSet( aAtom, fNeed ) ) continue;
        AtomSetFlags( aAtom, fFlags );
    }
}




/*
 *      ContainerWithoutFlagsResetAtomFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Loop through all ATOMs in the CONTAINER and set
 *      the flags desired for all ATOMs that have fNeed reset.
 */
void    
ContainerWithoutFlagsResetAtomFlags( CONTAINER cCont, 
                FLAGS fNeed, FLAGS fFlags )
{
LOOP            lAtoms;
ATOM            aAtom;

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) != NULL ) {
        if ( bAtomFlagsSet( aAtom, fNeed ) ) continue;
        AtomResetFlags( aAtom, fFlags );
    }
}




/*
 *      sContainerDescriptor
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the descriptor for the container.
 *      The descriptor is the name : sequence number.
 *
 *      Arguments:
 *              cCont - The container whose full descriptor is required.
 *              sDesc - A string capable of containing the
 *                      descriptor.
 *      Returns:
 *              sDesc
 */
char *
sContainerDescriptor( CONTAINER cCont, char *sDesc )
{
    if ( iObjectType(cCont) == UNITid ) {
        strcpy( sDesc, "" );
    } else {
        sprintf( sDesc, "%c<%s %d>", iObjectType(cCont), 
                sContainerName(cCont), 
                iContainerSequence(cCont) );
    }
    return(sDesc);
}




/*
 *      sContainerFullDescriptor
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Print the FULL name of the container.
 *      The FULL name consists of the names of the containers 
 *      of this container from outer most to this one.
 *      The name of each container is the descriptor of
 *      the container, the name:sequence#.
 *
 *      Arguments:
 *              cCont - The container whose full descriptor is required.
 *              sFullDesc -     A string capable of containing the
 *                              full descriptor.
 *      Returns:
 *              sFullDesc
 */
char *
sContainerFullDescriptor( CONTAINER cCont, char *sFullDesc )
{
STRING          sDesc;
CONTAINER       cTemp;

    if ( cCont == NULL ) {
        strcpy( sFullDesc, "!NULL!" );
        return(sFullDesc);
    }
    cTemp = cCont;
    strcpy( sFullDesc, "" );
    do {
        sContainerDescriptor( cTemp, sDesc );
        strcat( sDesc, "." );
        strcat( sDesc, sFullDesc );
        strcpy( sFullDesc, sDesc );
        cTemp = cContainerWithin( cTemp );
    } while ( cTemp != NULL );
    sFullDesc[strlen(sFullDesc)-1] = '\0';
    return(sFullDesc);
}




/*
 *      ContainerCheck
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Run tests on the CONTAINER, printing warning messages 
 *      that tell whether the CONTAINER is ready to have calculations
 *      run on it.
 *      The CONTAINER then runs tests on its contents.
 *      The routine counts the number of errors and warnings.
 *
 *      Arguments:
 *              cCont - The CONTAINER to test.
 *              iPErrors - Add the number of errors found.
 *              iPWarnings - Add the number of warnings found.
 *
 *TODO: Check for unit charge, undefined residue types, etc.
 */
void    
ContainerCheck( CONTAINER cCont, int *iPErrors, int *iPWarnings )
{
    switch ( iObjectType(cCont) ) {
        case CONTAINERid:
            DFATAL( ("Attempt to test a raw CONTAINER") );
            break;
        case UNITid:
            UnitCheck( (UNIT)cCont, iPErrors, iPWarnings );
            break;
        case MOLECULEid:
            MoleculeCheck( (MOLECULE)cCont, iPErrors, iPWarnings );
            break;
        case RESIDUEid:
            ResidueCheck( (RESIDUE)cCont, iPErrors, iPWarnings );
            break;
        case ATOMid:
            AtomCheck( (ATOM)cCont, iPErrors, iPWarnings );
            break;
        default:
            DFATAL( ("Attempting to check CONTAINER with id:%c", 
                        iObjectType(cCont)) );
            break;
    }
}




/*
 *      bContainerContainedBy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return TRUE if the first CONTAINER ( cIn ) is contained by
 *      ( cOut ).  Check this by following up the linked list of
 *      CONTAINERs that contain cIn.
 */
BOOL    
bContainerContainedBy( CONTAINER cIn, CONTAINER cOut )
{
CONTAINER       cTemp;

    cTemp = cIn;
    while ( cTemp != NULL ) {
        cTemp = cContainerWithin( cTemp );
        if ( cTemp == cOut ) return(TRUE);
    }
    return(FALSE);
}




/*
 *      ContainerYouAreBeingRemoved
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This routine calls the correct routine for each
 *      type of CONTAINER to inform it that it is
 *      being removed from the CONTAINER that contains it.
 *      There is a check to see if the object being
 *      informed is an INTERNAL, which is not strictly a
 *      CONTAINER but is treated like one for some purposes.
 */
void    
ContainerYouAreBeingRemoved( CONTAINER cCont )
{
    switch ( iObjectType(cCont) ) {
        case CONTAINERid:
            DFATAL( ("Attempt to inform a raw CONTAINER") );
            break;
        case UNITid:
            UnitYouAreBeingRemoved( (UNIT)cCont );
            break;
        case MOLECULEid:
            MoleculeYouAreBeingRemoved( (MOLECULE)cCont );
            break;
        case RESIDUEid:
            ResidueYouAreBeingRemoved( (RESIDUE)cCont );
            break;
        case ATOMid:
            AtomYouAreBeingRemoved( (ATOM)cCont );
            break;
        case INTERNALid:
            break;
        default:
            DFATAL( ("Attempting to inform CONTAINER with id:%c", 
                        iObjectType(cCont)) );
            break;
    }
}





/*
 *      ContainerIAmBeingRemoved
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This routine calls the correct routine for each
 *      type of CONTAINER to inform it that a CONTAINER
 *      within it is being removed.
 *      There is a check to see if the object being
 *      informed is an INTERNAL, which is not strictly a
 *      CONTAINER but is treated like one for some purposes.
 */
void    
ContainerIAmBeingRemoved( CONTAINER cCont, CONTAINER cRemoved )
{

        /* If the CONTAINER being informed is NULL then we are */
        /* at the top of the hierarchy, just return */

    if ( cCont == NULL ) return;

    switch ( iObjectType(cCont) ) {
        case CONTAINERid:
            DFATAL( ("Attempt to inform a raw CONTAINER") );
            break;
        case UNITid:
            UnitIAmBeingRemoved( (UNIT)cCont, (CONTAINER)cRemoved );
            break;
        case MOLECULEid:
            MoleculeIAmBeingRemoved( (MOLECULE)cCont, (CONTAINER)cRemoved );
            break;
        case RESIDUEid:
            ResidueIAmBeingRemoved( (RESIDUE)cCont, (CONTAINER)cRemoved );
            break;
        case ATOMid:
            AtomIAmBeingRemoved( (ATOM)cCont, (CONTAINER)cRemoved );
            break;
        default:
            DFATAL( ("Attempting to inform CONTAINER with id:%c", 
                        iObjectType(cCont)) );
            break;
    }
}





/*
 *      ContainerSetAttribute
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This routine changes an attribute of the CONTAINER.
 */
void    
ContainerSetAttribute( CONTAINER cCont, STRING sAttribute, OBJEKT oValue )
{

        /* If the CONTAINER being informed is NULL then we are */
        /* at the top of the hierarchy, just return */

    if ( cCont == NULL ) return;

    switch ( iObjectType(cCont) ) {
        case CONTAINERid:
            DFATAL( ("Attempt to set an attribute of a raw CONTAINER") );
            break;
        case UNITid:
            UnitSetAttribute( (UNIT)cCont, sAttribute, oValue );
            break;
        case MOLECULEid:
            MoleculeSetAttribute( (MOLECULE)cCont, sAttribute, oValue );
            break;
        case RESIDUEid:
            ResidueSetAttribute( (RESIDUE)cCont, sAttribute, oValue );
            break;
        case ATOMid:
            AtomSetAttribute( (ATOM)cCont, sAttribute, oValue );
            break;
        default:
            DFATAL( ("Attempting to change attribute of unknown OBJEKT ID:%c", 
                        iObjectType(cCont)) );
            break;
    }
}




/*
 *      ContainerTotalCharge
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Total up the unperturbed and perturbed charge of
 *      all ATOMs within the CONTAINER.
 */
void    
ContainerTotalCharge( CONTAINER cCont, double *dPCharge, double *dPPertCharge )
{
LOOP            lAtoms;
ATOM            aAtom;

    *dPCharge = 0.0;
    *dPPertCharge = 0.0;
    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        (*dPCharge) += dAtomCharge(aAtom);
        (*dPPertCharge) += dAtomPertCharge(aAtom);
    }
}




/*
 *-----------------------------------------------------------
 *
 *      CONTAINER DISPLAYER routines.
 *
 *      These routines are used in Event Driven environments
 *      to update any windows that are displaying the contents
 *      of CONTAINERS.
 *
 *      Whenever a CONTAINER is modified, it calls
 *      ContainerDisplayerUpdate which checks if the
 *      CONTAINER is sensitive and updates its DISPLAYERs
 *      if it is, otherwise it will be marked as MODIFIED
 *      and when it becomes sensitive again its DISPLAYERs
 *      will be updated.
 *      Most of this functionality is handled in DISPLAYER routines.
 *
 *      When CONTAINERs are modified, they update their own DISPLAYER
 *      and then send the modified message up to their parent.
 *      This allows things like UNIT displayers to be notified
 *      when an ATOM within them has been modified.
 *
 */

/*
 *      ContainerDisplayerUpdate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Update the DISPLAYER for this container, and then
 *      send the message up to the CONTAINER that contains this one.
 */
void    
ContainerDisplayerUpdate( CONTAINER cCont )
{

    if ( dContainerDisplayer(cCont) == NULL ) 
        return;

    if ( bDisplayerSensitive(dContainerDisplayer(cCont)) ) {
        DisplayerUpdate(dContainerDisplayer(cCont));

                /* If (cCont) is within another CONTAINER */
                /* then send the ContainerDisplayerUpdate */
                /* message to it */

        if ( cContainerWithin(cCont) != NULL ) {
            ContainerDisplayerUpdate(cContainerWithin(cCont));
        }
    } else {
        DisplayerSetModified(dContainerDisplayer(cCont));
    }
}



/*
 *      ContainerTreeMakeInsensitive
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Make all of the CONTAINERs within this one, and
 *      including this one insensitive.
 *
 *      This is a performance enhancement that can be called
 *      to prevent a large number of updates to DISPLAYERs
 *      being made during a large number of modifications
 *      to a CONTAINER and its contents.
 */
void    
ContainerTreeMakeInsensitive( CONTAINER cCont )
{
LOOP            lContainer;
CONTAINER       cCur;

    lContainer = lLoop( (OBJEKT)cCont, CONTAINERS );
    while ( (cCur = (CONTAINER)oNext(&lContainer)) ) {
        DisplayerSetSensitive(dContainerDisplayer(cCur),FALSE);
    }
}




/*
 *      ContainerTreeMakeSensitive
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Make the CONTAINER tree under (cCont) sensitive and
 *      then make this CONTAINER sensitive.
 *      If this CONTAINER has been modified then update its
 *      DISPLAYER.
 *
 *      This is a performance enhancement that can be called
 *      to update all windows after a large number of modifications
 *      have been made.
 */
void    
ContainerTreeMakeSensitive( CONTAINER cCont )
{
LOOP            lContents;
CONTAINER       cCur;

        /* If the OBJEKT is an ATOM then it will have no  */
        /* CONTAINERs within it */

    if ( iObjectType(cCont) != ATOMid ) {
        lContents = lLoop( (OBJEKT)cCont, DIRECTCONTENTS );
        while ( (cCur = (CONTAINER)oNext(&lContents)) ) {
            ContainerTreeMakeSensitive(cCur);
        }
    }

                /* Make this CONTAINER sensitive and */
                /* if it has been modified then update it */

    DisplayerSetSensitive(dContainerDisplayer(cCont),TRUE);
    if ( bDisplayerModified(dContainerDisplayer(cCont)) ) {
        DisplayerUpdate(dContainerDisplayer(cCont));
    }
}



/*
 *      ContainerResetAllCopyPointers
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Reset all of the internal pointers within CONTAINERS
 *      contained by this CONTAINER.
 *      This is required after Duplication.
 */
void    
ContainerResetAllCopyPointers( CONTAINER cTop )
{
LOOP            lContainerLoop;
CONTAINER       cCont;

    lContainerLoop = lLoop( (OBJEKT)cTop, CONTAINERS );
    while ( ( cCont = (CONTAINER)oNext( &lContainerLoop ) ) != NULL ) {
        ContainerResetPointers(cCont);
    }

            /* Reset the pointers for the top CONTAINER itself */

    ContainerResetPointers( cTop );
}

/*
 *      bContainerSpaceConflict
 *      
 *      Author: David A. Rivkin (1992)
 *
 *      Determine if any of the atoms of either container are within
 *      their resective Van der Waals radii of each other.
 *
 */
BOOL 
bContainerSpaceConflict( CONTAINER cCont1, CONTAINER cCont2 )
{
LOOP    lAtoms1, lAtoms2;
ATOM    aAtom1, aAtom2;

    lAtoms1 = lLoop( (OBJEKT)cCont1, ATOMS );
    while (( aAtom1 = (ATOM)oNext( &lAtoms1 ))) {
        lAtoms2 = lLoop( (OBJEKT)cCont2, ATOMS );
        while (( aAtom2 = (ATOM)oNext( &lAtoms2 ))) {
            if ( bAtomSpaceConflict( aAtom1, aAtom2 )) {
                return( TRUE );
            }
        }
    }
    return(FALSE);
}



