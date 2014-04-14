/*
 *      File: xAction.c
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
 *              Handle the user interaction with the TANK.
 *
 *      WARNING:
 *              This file directly accesses the TANK internal variables.
 *              This is because this file is conceptually part of the
 *              xTank.c file but is kept seperate because it deals
 *              with the user input to the TANK.
 */


#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#include <X11/cursorfont.h>

#include        "basics.h"


#include        "build.h"
#include        "vector.h"
#include        "matrix.h"
#include        "threed.h"
#include        "classes.h"
#include        "select.h"
#include        "varArray.h"
#include        "xTank.h"
#include        "xAction.h"









/*
 *      ActionPickAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the atom that the user picked.
 *      Return the atom with the largest Z coordinate that was
 *      drawn near the mouse hot spot.
 *      Return NULL if there is no atom.
 */
static void
ActionPickAtom( TANK tTank, int iX, int iY, ATOM *aPAtom, int *iPX, int *iPY )
{
GRAPHATOMt      *gaPA;
ATOM            *aPA;
ATOM            aA, aAtom;
int             i, iMax, iDX, iDY;
GVECTOR         *gvPA;
NUMBR           nDist;

    *aPAtom = NULL;
    iMax = iVarArrayElementCount( tTank->tank.vaAtomPtrs );
    if ( !iMax ) {
        return;
    }

        /* Now quickly loop through all of the ATOMs */
        /* Looking for those that would have coordinates */
        /* near the mouse position in iX, iY */

    aAtom = NULL;
    nDist = -9999999.9;
    aPA = PVAI( tTank->tank.vaAtomPtrs, ATOM, 0 );
    for ( i=0; i<iMax; i++, aPA++ ) {
/*
fprintf(stderr, " atom %s type %c/%d\n", 
sContainerName( *aPA ),
iObjectType( *aPA ),
iObjectType( *aPA ));
*/
        aA = *aPA;
        gaPA = (GRAPHATOMt*) PAtomGraphicsPointer(aA);
        gvPA = &(gaPA->gvScreen);
        if ( !bAtomVisible(aA) ) continue;
        if ( !bGVectorMiddle( *gvPA ) ) continue;
        if ( (iDX = abs(((int)(nGVX(*gvPA))) - iX)) < PICKDIST &&
             (iDY = abs(((int)(nGVY(*gvPA))) - iY)) < PICKDIST ) {
            if ( nDist < nGVZ(*gvPA) ) {
                aAtom = aA;
                nDist = nGVZ(*gvPA);
                *iPX = nGVX(*gvPA);
                *iPY = nGVY(*gvPA);
            }
        }
    }
    *aPAtom = aAtom;
/*
fprintf(stderr, "done/ found one %s 0x%x type %d\n",  
sContainerName(aAtom), aAtom, iObjectType(aAtom));
*/
}






/*
 *      ActionPickAtomsInBox
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Call the function (fCallback) with each ATOM that
 *      lies within the box marked by (iXTopLeft), (iYTopLeft)
 *      to (iXBottomRight), (iYBottomRight).
 *
 *      Pass the data (Pdata) to the callback after the ATOM.
 *
 *      The callback has the structure:
 *              void callback( aAtom, PData )
 *              ATOM            aAtom;
 *              GENP            PData;
 *
 */
static void
ActionPickAtomsInBox( TANK tTank, int iXTopLeft, int iYTopLeft,
                int iXBottomRight, int iYBottomRight, 
                VFUNCTION fCallback, GENP PData )
{
ATOM            *aPA;
GRAPHATOMt      *gaPA;
ATOM            aA;
int             i, iMax;
GVECTOR         *gvPA;
int             iX, iY;

    MESSAGE(( "ActionPickAtomsInBox: %d, %d - %d, %d\n", 
                iXTopLeft, iYTopLeft, iXBottomRight, iYBottomRight ));

    iMax = iVarArrayElementCount( tTank->tank.vaAtomPtrs );
    if ( !iMax )
        return;

        /* Now quickly loop through all of the ATOMs */
        /* Looking for those that would have coordinates */
        /* near the mouse position in iX, iY */

    aPA = PVAI( tTank->tank.vaAtomPtrs, ATOM, 0 );
    for ( i=0; i<iMax; i++, aPA++ ) {
        aA = *aPA;
        gaPA = (GRAPHATOMt*)PAtomGraphicsPointer(aA);
        gvPA = &(gaPA->gvScreen);
        if ( !bAtomVisible(aA) ) continue;
        if ( !bGVectorMiddle( *gvPA ) ) continue;
        iX = (int)nGVX(*gvPA);
        iY = (int)nGVY(*gvPA);
        if ( ( iXTopLeft < iX && iX < iXBottomRight ) &&
             ( iYTopLeft < iY && iY < iYBottomRight ) ) {
            fCallback( aA, PData );
        }
    }

}






/*
 *      ActionPickBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the bond that the user picked.
 *      Return the first bond found that is close to the X,Y position.
 *      Return the two atoms of the bond.
 *      Return the atom that the pick was closest to in aPA and
 *      the other atom in aPB.  Return TRUE in *bPPickedSide if 
 *      the pick was DEFINITLY on one side of the bond. 
 *      Return FALSE if there is no bond.
 */
static void
ActionPickBond( TANK tTank, int iX, int iY, ATOM *aPAtomA, ATOM *aPAtomB, 
                BOOL *bPFound, BOOL *bPPickedSide )
{
ATOM            *aPA;
ATOM            aA, aB;
int             i, j, iAtoms;
int             iXA, iYA, iXB, iYB;
GVECTOR         *gvPA, *gvPB;
BOOL            bInside;
int             iTop, iBottom, iDX, iDY, iDistSquared;
double          dFraction = 0.0;
BOOL            bFoundOne;

                /* By default we found nothing */


    *bPFound = FALSE;
    bFoundOne = FALSE;

    iAtoms = iVarArrayElementCount( tTank->tank.vaAtomPtrs );
    if ( !iAtoms )
        return;

        /* Now quickly loop through all of the ATOMs */
        /* Looking for those that would have coordinates */
        /* near the mouse position in iX, iY */


    aPA = PVAI( tTank->tank.vaAtomPtrs, ATOM, 0 );
    for ( i=0; !bFoundOne && i<iAtoms; i++, aPA++ ) {
        aA = *aPA;
        if ( !bAtomVisible(aA) ) 
            continue;
        for ( j=0; j<iAtomCoordination(aA); j++ ) {
            aB = aAtomBondedNeighbor( aA, j );
            if ( iAtomId(aA) < iAtomId(aB) ) 
                continue;
            if ( !bAtomVisible(aB) ) 
                continue;
            gvPA = &(((GRAPHATOMt*)PAtomGraphicsPointer(aA))->gvScreen);
            gvPB = &(((GRAPHATOMt*)PAtomGraphicsPointer(aB))->gvScreen);
            if ( !bX3dVisibleLineClip( tTank->tank.x3dEngine, 
                                 &iXA, &iYA, &iXB, &iYB, 
                                 gvPA, gvPB ) ) 
                continue;

            MESSAGE(( "Looking at bond: %s - %s\n",
                        sContainerName(aA), sContainerName(aB) ));

                /* Check if the pick point is within the bounding box */
                /* of the line */

            bInside = FALSE;
            if ( iXA<iXB ) {
                if ( iXA-PICKDIST<iX && iX<iXB+PICKDIST ) {
                    if ( iYA<iYB ) {
                        if ( iYA-PICKDIST<iY && iY<iYB+PICKDIST ) 
                                bInside = TRUE;
                    } else {
                        if ( iYB-PICKDIST<iY && iY<iYA+PICKDIST ) 
                                bInside = TRUE;
                    }
                }
            } else {
                if ( iXB-PICKDIST<iX && iX<iXA+PICKDIST ) {
                    if ( iYA<iYB ) { 
                        if ( iYA-PICKDIST<iY && iY<iYB+PICKDIST ) 
                                bInside = TRUE;
                    } else {
                        if ( iYB-PICKDIST<iY && iY<iYA+PICKDIST ) 
                                bInside = TRUE;
                    }
                }
            }
            if ( !bInside ) 
                continue;

            MESSAGE(( "Inside bounding box\n" ));

                /* Find the fraction along the line where the */
                /* pick is closest to */

            iTop =    (iX -iXA)*(iXB-iXA)+(iY -iYA)*(iYB-iYA);
            iBottom = (iXB-iXA)*(iXB-iXA)+(iYB-iYA)*(iYB-iYA);
            dFraction = ((double)iTop)/((double)iBottom);
            iDX = iXA + ( iXB - iXA )*dFraction;
            iDY = iYA + ( iYB - iYA )*dFraction;
            iDistSquared = (iX-iDX)*(iX-iDX)+(iY-iDY)*(iY-iDY);

            MESSAGE(( "Distance^2 from line = %d\n", iDistSquared ));

            if ( iDistSquared < PICKLINE ) {
                bFoundOne = TRUE;
                break;
            }
        }
    }

    if ( dFraction > 0.5 ) {
        *aPAtomA = aB;
        *aPAtomB = aA;
        dFraction = 1.0 - dFraction;
    } else {
        *aPAtomA = aA;
        *aPAtomB = aB;
    }
    *bPFound = bFoundOne;

    *bPPickedSide = ( dFraction < PICKLINESIDE );
}






/*
 *      zActionRotateXYTransform
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Rotate the object in the TANK around the vector
 *      perpendicular to the vector
 *      (nX,nY,0) and the vector going into the screen.
 *      Rotate around the line defined by gvAround and the above
 *      vector.
 *
 *      If bToWorld is TRUE then first transform the rotation
 *      axis to World coordinates.
 */
static void
zActionRotateXYTransform( TANK tTank, MATRIX mTransform, VECTOR *vPAround, 
                NUMBR nX1, NUMBR nY1, NUMBR nX2, NUMBR nY2, BOOL bToWorld )
{
VECTOR          vDrag, vNormal, vLine, vA, vB;
GVECTOR         gvTemp, gvLine, gvOrigin;
double          dAngle;
NUMBR           nX, nY;

    nX = nX2 - nX1;
    nY = nY2 - nY1;
    VectorDef( &vDrag, nX, nY, 0.0 );
    VectorDef( &vNormal, 0.0, 0.0, -1.0 );
    MESSAGE(( "Drag vector   = %lf, %lf, %lf\n",
                dVX(&vDrag), dVY(&vDrag), dVZ(&vDrag) ));
    MESSAGE(( "Normal vector = %lf, %lf, %lf\n",
                dVX(&vNormal), dVY(&vNormal), dVZ(&vNormal) ));
    vLine = vVectorCross( &vDrag, &vNormal );
    if ( bToWorld ) {
        GVectorDef( gvTemp, dVX(&vLine), dVY(&vLine), dVZ(&vLine) );
        X3dScreenToWorld( tTank->tank.x3dEngine, &gvTemp, &gvLine );
        GVectorDef( gvTemp, 0.0, 0.0, 0.0 );
        X3dScreenToWorld( tTank->tank.x3dEngine, &gvTemp, &gvOrigin );
        GVectorSubtract( gvLine, gvOrigin, gvLine );
        VectorDef( &vLine, nGVX(gvLine), nGVY(gvLine), nGVZ(gvLine) );
    }
    vA = *vPAround;
    vB = vVectorAdd( &vA, &vLine );
    MESSAGE(( "vA            = %lf, %lf, %lf\n",
                dVX(&vA), dVY(&vA), dVZ(&vA) ));
    MESSAGE(( "vB            = %lf, %lf, %lf\n",
                dVX(&vB), dVY(&vB), dVZ(&vB) ));
    dAngle = sqrt(nX*nX + nY*nY)*DEGTORAD;
    MatrixRotateAround( mTransform, &vA, &vB, dAngle );
    MESSAGE(( "Transformation matrix: \n" ));
MESSAGEEXECUTE( {
    MatrixPrint( mTransform );
                } );

}





/*
 *      ActionRotateXY
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Rotate the object in the TANK around the vector
 *      perpendicular to the vector
 *      (nX,nY,0) and the vector going into the screen.
 *      Rotate around the line defined by gvAround and the above
 *      vector.
 */
static void
ActionRotateXY( TANK tTank, VECTOR *vPAround, 
                NUMBR nX1, NUMBR nY1, NUMBR nX2, NUMBR nY2 )
{
MATRIX          mTransform;
GMATRIX         gmTransform, gmY, gmRotate;

    zActionRotateXYTransform( tTank, mTransform, vPAround, 
                                nX1, nY1, nX2, nY2, FALSE );
    GMatrixFromMatrix( gmTransform, mTransform );

    X3dGetWorldRotate( tTank->tank.x3dEngine, gmRotate );
    GMatrixMultiply( gmY, gmTransform, gmRotate );
    X3dSetWorldRotate( tTank->tank.x3dEngine, gmY );

}


/*
 *      zActionRotateZTransform
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Rotate the object in the TANK around the vector
 *      perpendicular to the vector
 *      (nX,nY,0) and the vector going into the screen.
 *      Rotate around the line defined by gvAround and the above
 *      vector.
 *      If bToWorld is TRUE then first transform the rotation
 *      axis to World coordinates.
 */
static void
zActionRotateZTransform( TANK tTank, MATRIX mTransform, 
                VECTOR *vPAround, NUMBR nX1, NUMBR nY1, NUMBR nX2, NUMBR nY2, 
                BOOL bToWorld )
{
double          dAngle1, dAngle2, dAngle;
VECTOR          vLine, vA, vB;
GVECTOR         gvTemp, gvLine, gvOrigin;

    MESSAGE(( "nX1, nY1 = %lf, %lf\n", nX1, nY1 ));
    MESSAGE(( "nX2, nY2 = %lf, %lf\n", nX2, nY2 ));

    dAngle1 = myAcos( nX1/sqrt(nY1*nY1 + nX1*nX1) );
    if ( nY1 < 0.0 ) dAngle1 = -dAngle1;
    dAngle2 = myAcos( nX2/sqrt(nY2*nY2 + nX2*nX2) );
    if ( nY2 < 0.0 ) dAngle2 = -dAngle2;
    MESSAGE(( "Angle 1 = %lf    Angle 2 = %lf\n",
                dAngle1/DEGTORAD, dAngle2/DEGTORAD ));

    vA = *vPAround;
    VectorDef( &vLine, 0.0, 0.0, -1.0 );
     if ( bToWorld ) {
        GVectorDef( gvTemp, dVX(&vLine), dVY(&vLine), dVZ(&vLine) );
        X3dViewToWorld( tTank->tank.x3dEngine, &gvTemp, &gvLine );
        GVectorDef( gvTemp, 0.0, 0.0, 0.0 );
        X3dViewToWorld( tTank->tank.x3dEngine, &gvTemp, &gvOrigin );
        GVectorSubtract( gvLine, gvLine, gvOrigin );
        VectorDef( &vLine, nGVX(gvLine), nGVY(gvLine), nGVZ(gvLine) );
    }
    vB = vVectorAdd( &vA, &vLine );
    dAngle = dAngle2 - dAngle1;
    MESSAGE(( "vA            = %lf, %lf, %lf\n",
                dVX(&vA), dVY(&vA), dVZ(&vA) ));
    MESSAGE(( "vB            = %lf, %lf, %lf\n",
                dVX(&vB), dVY(&vB), dVZ(&vB) ));
    MatrixRotateAround( mTransform, &vA, &vB, dAngle );
    MESSAGE(( "Transformation matrix: \n" ));
MESSAGEEXECUTE( {
    MatrixPrint( mTransform );
                } );

}


/*
 *      ActionRotateZ
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set up the TANK to rotate the object in the TANK around 
 *      the Z axis by the angle that is the difference between the
 *      angles cut by nX1,nY1 and nX2,nY2.
 */
static void
ActionRotateZ( TANK tTank, VECTOR *vPAround, 
                NUMBR nX1, NUMBR nY1, NUMBR nX2, NUMBR nY2 )
{
MATRIX          mTransform;
GMATRIX         gmTransform, gmY, gmRotate;

    zActionRotateZTransform( tTank, mTransform, vPAround, 
                                nX1, nY1, nX2, nY2, FALSE );
    GMatrixFromMatrix( gmTransform, mTransform );

    X3dGetWorldRotate( tTank->tank.x3dEngine, gmRotate );
    GMatrixMultiply( gmY, gmTransform, gmRotate );
    X3dSetWorldRotate( tTank->tank.x3dEngine, gmY );
}






/*
 *      aActionCreateAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create an ATOM using the current settings.
 */
static ATOM
aActionCreateAtom( TANK tTank, int iNumber, int iX, int iY )
{
ATOM            a1;
STRING          sName, sElement;
GVECTOR         gvScreen, gvWorld;
VECTOR          vPos;

                /* Obtain the WORLD coordinates of the atom */

    nGVX(gvScreen) = iX;
    nGVY(gvScreen) = iY;
    nGVZ(gvScreen) = nX3dViewCenter( tTank->tank.x3dEngine );
    X3dScreenToWorld( tTank->tank.x3dEngine, &gvScreen, &gvWorld );
    VectorDef( &vPos, nGVX(gvWorld), nGVY(gvWorld), nGVZ(gvWorld) );

    MESSAGE(( "Creating an atom at: %f, %f, %f\n",
                nGVX(gvWorld), nGVY(gvWorld), nGVZ(gvWorld) ));

    a1 = (ATOM)oCreate(ATOMid);
    AtomSetElement( a1, tTank->tank.iCurrentDrawingElement );
    sprintf( sName, "%s%d", 
        sElementName( tTank->tank.iCurrentDrawingElement, sElement ),
        iNumber );
    ContainerSetName( a1, sName );
    AtomSetPosition( a1, vPos );

                /* Define the flags, to turn off ATOMPOSITIONKNOWN */
                /* instead tell the ATOM that its position has only */
                /* been drawn and should only be used to display */
    AtomDefineFlags( a1, ATOMNEEDSBUILD | ATOMPOSITIONDRAWN );
    return(a1);
}








/*
 *      ActionUseLineSegmentToDraw
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Figure out what to do with the line segment. 
 *      If it was shorter than PICKXDIST,PICKYDIST then
 *      it is not a line segment, it is an attempt to change the
 *      bond order of a bond.
 *
 *      If the user is attempting to create a bond then:
 *      search for the two atoms at the ends of the line segment.
 */
static void
ActionUseLineSegmentToDraw( TANK tTank )
{
int             iLen;
int             iXStart, iYStart, iXStop, iYStop;
ATOM            a1, a2;
BOOL            bGot1, bGot2;
int             iX, iY, iXTemp, iYTemp;
LOOP            lResidues;
RESIDUE         rRes;
BOOL            bFound, bPickedSide;
int             iOrder;
STRING          sNewName;
int             iNumber;


    DisplayerAccumulateUpdates();

    iXStart = tTank->tank.iXStart;
    iYStart = tTank->tank.iYStart;
    iXStop  = tTank->tank.iXStop;
    iYStop  = tTank->tank.iYStop;
    iLen = (iXStop-iXStart)*(iXStop-iXStart)+
                (iYStop-iYStart)*(iYStop-iYStart);
    MESSAGE(( "Length^2 of line = %d\n", iLen ));

                /* If the line segment is very short then */
                /* the user is not drawing a bond, but rather */
                /* trying to change the order of a bond, or trying */
                /* to change the element type of an atom. */
    if ( iLen < pow2(PICKDIST) ) {

        iX = ( iXStart + iXStop ) / 2;
        iY = ( iYStart + iYStop ) / 2;

        /* Search for an atom or bond to change type or */
        /* increment the bond order */

                /* First check if the user wants to change the element */
                /* type of an atom */

        ActionPickAtom( tTank, iX, iY, &a1, &iXTemp, &iYTemp );
        if ( a1 != NULL ) {
            AtomSetElement( a1, tTank->tank.iCurrentDrawingElement );
            sElementName( tTank->tank.iCurrentDrawingElement, sNewName );
            sprintf( sNewName, "%s%d", sNewName, iContainerSequence(a1) );
            AtomSetName( a1, sNewName );
            TankRedisplayUnit( tTank );
        } else {

            MESSAGE(( "The user did not pick an atom. Checking bond\n" ));
            
                /* The user did not pick an atom, so check if they */
                /* picked a bond.  If they did, increment the order */
                /* of that bond */

            ActionPickBond( tTank, iX, iY, &a1, &a2, &bFound, &bPickedSide );
            if ( bFound ) {
                iOrder = iAtomFindBondOrder( a1, a2 );
                MESSAGE(( "Found bond between: %s - %s  bond order = %d\n", 
                        sContainerName(a1), sContainerName(a2), iOrder ));

                        /* Increment the bond order */

                switch ( iOrder ) {
                    case BONDSINGLE:
                        iOrder = BONDDOUBLE;
                        break;
                    case BONDDOUBLE:
                        iOrder = BONDTRIPLE;
                        break;
                    case BONDTRIPLE:
                        iOrder = BONDAROMATIC;
                        break;
                    case BONDAROMATIC:
                        iOrder = BONDSINGLE;
                        break;
                    default:
                        DFATAL(( "Illegal bond order" ));
                        break;
                }

                        /* Change the bond order */

                AtomFindSetBondOrder( a1, a2, iOrder );
                TankRedisplayUnit( tTank );

                        /* Did not find a bond, this means the user */
                        /* wants to draw a new atom */
                        /* Create a new atom and add it to an arbitrary */
                        /* RESIDUE */

            } else {
                lResidues = lLoop( (OBJEKT)tTank->tank.uUnit, RESIDUES );
                rRes = (RESIDUE)oNext(&lResidues);
                iNumber = iContainerNextChildsSequence(rRes);
                a1 = aActionCreateAtom( tTank, iNumber, iXStart, iYStart );
                ContainerAdd( (CONTAINER)rRes, (OBJEKT)a1 );
                TankRedisplayUnit(tTank);
            }
        }       
    } else {

                /* Find the atoms at the ends of the line segment */

        bGot1 = TRUE;
        bGot2 = TRUE;
        ActionPickAtom( tTank, iXStart, iYStart, &a1, &iX, &iY );
        ActionPickAtom( tTank, iXStop, iYStop, &a2, &iX, &iY );
        if ( a1 == NULL ) bGot1 = FALSE;
        if ( a2 == NULL ) bGot2 = FALSE;

                /* If both atoms are new then add them to the first */
                /* residue in the UNIT */
                /* TODO: Figure out a better place to put new */
                /* TODO: atoms, rather than the first residue found */

        if ( !bGot1 && !bGot2 ) {
            lResidues = lLoop( (OBJEKT)tTank->tank.uUnit, RESIDUES );
            rRes = (RESIDUE)oNext( &lResidues );
            iNumber = iContainerNextChildsSequence(rRes);
            a1 = aActionCreateAtom( tTank, iNumber, iXStart, iYStart );
            ContainerAdd( (CONTAINER)rRes, (OBJEKT)a1 );
            iNumber = iContainerNextChildsSequence(rRes);
            a2 = aActionCreateAtom( tTank, iNumber, iXStop, iYStop );
            ContainerAdd( (CONTAINER)rRes, (OBJEKT)a2 );
        } else if ( !bGot1 ) {
            rRes = (RESIDUE)cContainerWithin(a2);
            iNumber = iContainerNextChildsSequence(rRes);
            a1 = aActionCreateAtom( tTank, iNumber, iXStart, iYStart );
            ContainerAdd( (CONTAINER)rRes, (OBJEKT)a1 );
        } else if ( !bGot2 ) {
            rRes = (RESIDUE)cContainerWithin(a1);
            iNumber = iContainerNextChildsSequence(rRes);
            a2 = aActionCreateAtom( tTank, iNumber, iXStop, iYStop );
            ContainerAdd( (CONTAINER)rRes, (OBJEKT)a2 );
        }
                /* Create a bond between the atoms */
                /* Only if there is not already a bond, and */
                /* the maximum number of bonds from each ATOM */
                /* is not exceeded */

        if ( !bAtomBondedTo( a1, a2 ) ) {
            if ( iAtomCoordination(a1) < MAXBONDS &&
                 iAtomCoordination(a2) < MAXBONDS ) {
                AtomBondToOrder( a1, a2, BONDSINGLE );
            } else {
                TankBeep(tTank);
            }
        } else {
            TankBeep(tTank);
        }
                /* Display the result */

        TankRedisplayUnit(tTank);
    }

    DisplayerReleaseUpdates();

    return;
}






/*
 *      ActionUseLineSegmentToErase
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Erase the atom or bond that the user picked.
 *
 */
static void
ActionUseLineSegmentToErase( TANK tTank )
{
ATOM            a1, a2;
int             iX, iY, iXTemp, iYTemp;
BOOL            bFound, bPickedSide;
CONTAINER       cContainer1, cContainer2;
RESIDUE         rRes;
LOOP            lResidues;
int             iDelPdbSeq, iPdbSeq;


    iX = tTank->tank.iXStart;
    iY = tTank->tank.iYStart;

    ActionPickAtom( tTank, iX, iY, &a1, &iXTemp, &iYTemp );
    if ( a1 != NULL ) {
        cContainer1 = cContainerWithin(a1);
        REF( a1 );  /* bContainerRemove() needs this */
        if ( !bContainerRemove( cContainer1, (OBJEKT)a1 ) ) {
            DFATAL(( "Atom was not within its container!!!!" ));
        }
        Destroy( (OBJEKT *)&a1 );  
        if ( iContainerNumberOfChildren( cContainer1 ) == 0 ) {
            /*
             *  residue gone; renumber all in unit
             */
            iDelPdbSeq = iResiduePdbSequence( cContainer1 );

            lResidues = lLoop( (OBJEKT)tTank->tank.uUnit, RESIDUES );
            while ( (rRes = (RESIDUE)oNext( &lResidues )) ) {
                iPdbSeq = iResiduePdbSequence( rRes );
                if ( iPdbSeq > iDelPdbSeq )
                        ResidueSetPdbSequence( rRes, iPdbSeq - 1 );
            }
            /*
             *  delete all empty things in the hierarchy, except
             *  the unit itself
             */
            while ( (iContainerNumberOfChildren( cContainer1 ) == 0)
                   && ( cContainer1 != (CONTAINER)tTank->tank.uUnit) ) {
                cContainer2 = cContainerWithin( cContainer1 );
                REF( cContainer1 );  /* bContainerRemove() needs this */
                if ( !bContainerRemove( cContainer2, (OBJEKT)cContainer1 ) )
                        DFATAL(( "containment problem" ));
                Destroy( (OBJEKT *)&cContainer1 );
                cContainer1 = cContainer2;
            }
        }
        TankRedisplayUnit( tTank );
    } else {

        MESSAGE(( "The user did not pick an atom. Checking bond\n" ));
            
                /* The user did not pick an atom, so check if they */
                /* picked a bond.  If they did, increment the order */
                /* of that bond */

        ActionPickBond( tTank, iX, iY, &a1, &a2, &bFound, &bPickedSide );
        if ( bFound ) {
            MESSAGE(( "Found bond between: %s - %s\n", 
                        sContainerName(a1), sContainerName(a2) ));

                        /* Remove the bond */

            AtomRemoveBond( a1, a2 );

            TankRedisplayUnit( tTank );
        }       
    }
}






/*
 *      zActionSelectHierarchy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Select an atom/residue/molecule, depending
 *      on the current select state.  The ATOM passed
 *      is an ATOM within the RESIDUE/molecule that
 *      the user wants to select.
 */
static void
zActionSelectHierarchy( TANK tTank, ATOM aAtom, BOOL bOn )
{
UNIT            uUnit;

    MESSAGE(( "Selecting\n" ));

    uUnit = tTank->tank.uUnit;

        /* Select the next level of organization: atom/ring/residue/mol */
        /* Loop until the level is selected.  Mainly loop if no ring */
        /* was found, then go straight to RESIDUE */

    switch ( tTank->tank.iTankSelectLastState ) {
        case TANKSELECTATOM:
            SelectAtom( aAtom, bOn );
            break;
        case TANKSELECTRESIDUE:
            SelectResidueWithAtom( uUnit, aAtom, bOn );
            break;
        case TANKSELECTMOLECULE:
            SelectMoleculeWithAtom( uUnit, aAtom, bOn );
            break;
        case TANKSELECTALL:
            SelectEverything( uUnit, bOn );
            break;
        default:
            /* Nothing */
            break;
    }
}





/*
 *      zActionSelectCallback
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Callback for picking ATOMs within a box.
 *      Call SelectAtom for each ATOM.
 */
static void
zActionSelectCallback( ATOM aAtom, BOOL *bPOn )
{
    SelectAtom( aAtom, *bPOn );
}



/*
 *      ActionSelectDeselectAtoms
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Select/DeSelect ATOMs within the TANK using data in the TANK.
 *      If iSelectionType is (ATOM_OR_CHAIN) then 
 *      use iXStart, iYStart, iXStop, iYStop to select
 *      an ATOM/RESIDUE/molecule/everything or a chain of ATOMs.
 *
 *      If iSelectionType is (TANK_SELECT_BOX) then
 *      use iXStart, iYStart, iXStop, iYStop to select all
 *      ATOMs within a box.
 *
 *      If (bOn) is TRUE then select the ATOMs, otherwise
 *      deselect them.
 *
 */
static void
ActionSelectDeselectAtoms( TANK tTank, BOOL bOn )
{
int             iLen, iXStart, iYStart, iXStop, iYStop;
int             iX, iY, iXTemp, iYTemp, iTemp;
ATOM            a1, a2;

    DisplayerAccumulateUpdates();

    MESSAGE(( "In ActionSelectDeselectAtoms\n" ));
    iXStart = tTank->tank.iXStart;
    iYStart = tTank->tank.iYStart;
    iXStop  = tTank->tank.iXStop;
    iYStop  = tTank->tank.iYStop;

                /* If the selection is of an ATOM or a CHAIN then */
                /* find the ATOM or the two ATOMs at the end of the */
                /* chain */

    if ( tTank->tank.iSelectionType == TANK_SELECT_ATOM_OR_CHAIN ) {
        iLen = (iXStop-iXStart)*(iXStop-iXStart)+
                (iYStop-iYStart)*(iYStop-iYStart);
        MESSAGE(( "Length^2 of line = %d\n", iLen ));

                /* If the line segment is very short then */
                /* the user is not drawing a bond, but rather */
                /* trying to select an atom/residue/everything */
        if ( iLen < pow2(PICKDIST) ) {

            iX = ( iXStart + iXStop ) / 2;
            iY = ( iYStart + iYStop ) / 2;

                /* Search for an atom or bond to change type or */
                /* increment the bond order */

                /* First check if the user wants to change the element */
                /* type of an atom */

            ActionPickAtom( tTank, iX, iY, &a1, &iXTemp, &iYTemp );
            if ( a1 != NULL ) {
                zActionSelectHierarchy( tTank, a1, bOn );
                TankRedisplayUnit( tTank );
            } else {
                TankBeep(tTank);
            }   
        } else {
            ActionPickAtom( tTank, iXStart, iYStart, &a1, &iX, &iY );
            ActionPickAtom( tTank, iXStop, iYStop, &a2, &iX, &iY );
            if ( a1 == NULL ) goto DONE;
            if ( a2 == NULL ) goto DONE;

            if ( !bSelectChainBetween( tTank->tank.uUnit, a1, a2, bOn ) ) {
                TankBeep(tTank);
            } else {
                TankRedisplayUnit(tTank);
            }
        }
    } else {

                /* If the user made a double click then select everything */

        if ( tTank->tank.tTimeSinceLastClick < SELECT_DELAY ) {
            MESSAGE(( "Selecting everything, selection = %s\n", sBOOL(bOn) ));
            SelectEverything( tTank->tank.uUnit, bOn );
            TankRedisplayUnit(tTank);
            goto DONE;
        }
                /* Otherwise the selection is of all ATOMs */
                /* within a box.  Create the box and */
                /* find all of the ATOMs within it and select */
                /* them */

        if ( iXStart > iXStop ) SWAP( iXStart, iXStop, iTemp );
        if ( iYStart > iYStop ) SWAP( iYStart, iYStop, iTemp );

        ActionPickAtomsInBox( tTank, iXStart, iYStart, iXStop, iYStop,
                                        zActionSelectCallback, (GENP)&bOn );
        TankRedisplayUnit(tTank);
    }

DONE:

    DisplayerReleaseUpdates();
    return;
}




/*
 *      ActionSelect
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Select ATOMs according to data in the TANK.
 */
static void
ActionSelect( TANK tTank )
{
    MESSAGE(( "In ActionSelect\n" ));
    ActionSelectDeselectAtoms( tTank, TRUE );
}




/*
 *      ActionDeSelect
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      DeSelect ATOMs according to data in the TANK.
 */
static void
ActionDeSelect( TANK tTank )
{
    MESSAGE(( "In ActionDeSelect\n" ));
    ActionSelectDeselectAtoms( tTank, FALSE );
}





/*
 *      zActionDrag
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Drag the selected ATOMs.
 */
static void
zActionDrag( TANK tTank )
{
VECTOR          vPos, vDir, vStart, vStop;
GVECTOR         gvScreen, gvWorld;
LOOP            lAtoms;
ATOM            aAtom;

    nGVX(gvScreen) = tTank->tank.iXStart;
    nGVY(gvScreen) = tTank->tank.iYStart;
    nGVZ(gvScreen) = nX3dViewCenter(tTank->tank.x3dEngine);
    X3dScreenToWorld( tTank->tank.x3dEngine, &gvScreen, &gvWorld );
    VectorDef( &vStart, nGVX(gvWorld), nGVY(gvWorld), nGVZ(gvWorld) );

    nGVX(gvScreen) = tTank->tank.iXStop;
    nGVY(gvScreen) = tTank->tank.iYStop;
    nGVZ(gvScreen) = nX3dViewCenter(tTank->tank.x3dEngine);
    X3dScreenToWorld( tTank->tank.x3dEngine, &gvScreen, &gvWorld );
    VectorDef( &vStop, nGVX(gvWorld), nGVY(gvWorld), nGVZ(gvWorld) );

    vDir = vVectorSub( &vStop, &vStart );

    MESSAGE(( "Moving by: %lf, %lf, %lf\n",
                dVX(&vDir), dVY(&vDir), dVZ(&vDir) ));

    lAtoms = lLoop( (OBJEKT)tTank->tank.uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
            vPos = vAtomPosition(aAtom);
            vPos = vVectorAdd( &vPos, &vDir );
            AtomSetPosition( aAtom, vPos );
        }
    }

    TankFastRedisplayUnit( tTank );

    tTank->tank.iXStart = tTank->tank.iXStop;
    tTank->tank.iYStart = tTank->tank.iYStop;

}





/*
 *      zActionTwistTorsions
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Twist around the selected torsions.
 *      Use the displacement of the mouse in the vertical direction
 *      to calculate the angle to rotate by.
 *      If the torsion is not in a ring then rotate all of the ATOMs
 *      to one side of the torsion.  If it is in a ring then
 *      rotate only the atoms bonded directly to the ATOM on
 *      one side of the torsion.
 */
static void
zActionTwistTorsions( TANK tTank )
{
double          dTwist;
int             i, iMax;
TWISTTORSIONt   *tT;

                /* Calculate angle to rotate by */


    iMax = iVarArrayElementCount(tTank->tank.vaTwistTorsions);
    if ( !iMax )
        return;

    dTwist = tTank->tank.iYStop - tTank->tank.iYStart;
    dTwist *= DEGTORAD;

    MESSAGE(( "Going to rotate torsions by: %lf\n",
                        dTwist/DEGTORAD ));

    tT = PVAI( tTank->tank.vaTwistTorsions, TWISTTORSIONt, 0 );
    for ( i=0; i<iMax; i++, tT++ ) {
        BuildRotateAroundBondFromTo( tTank->tank.uUnit,
                                        tT->aAtomInvisible, 
                                        tT->aAtomStart, 
                                        dTwist, 
                                        tT->bInRing );
    }

    tTank->tank.iXStart = tTank->tank.iXStop;
    tTank->tank.iYStart = tTank->tank.iYStop;

    TankFastRedisplayUnit( tTank );

}


/*
 *      zActionRotateAtoms
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Rotate the ATOMs in the vaRotateAtoms list around
 *      the point in aRotateCenter.
 */
static void
zActionRotateAtoms( TANK tTank )
{
ATOM            *aPCur, *aPLast;
MATRIX          mTransform;
NUMBR           nX1, nY1, nX2, nY2;
VECTOR          vNew;
int             iWidth, iHeight;

                /* Calculate the transform */

    iWidth = iX3dScreenWidth(tTank->tank.x3dEngine);
    iHeight = iX3dScreenHeight(tTank->tank.x3dEngine);
    nX1 = ( tTank->tank.iXStart - (iWidth>>1) );
    nY1 = -( tTank->tank.iYStart - (iHeight>>1));
    nX2 = ( tTank->tank.iXStop - (iWidth>>1));
    nY2 = -( tTank->tank.iYStop - (iHeight>>1));
    if ( bTankPointInCircle( tTank,
                tTank->tank.iXStop, tTank->tank.iYStop ) ) {
        zActionRotateXYTransform( tTank, mTransform, 
                &(vAtomPosition(tTank->tank.aRotateCenter)), 
                -nX1, nY1, -nX2, nY2, TRUE );
    } else {
        zActionRotateZTransform( tTank, mTransform,
                &(vAtomPosition(tTank->tank.aRotateCenter)), 
                nX1, nY1, nX2, nY2, TRUE );
    }

                /* Rotate the ATOMs */

    aPLast = PVAI( tTank->tank.vaRotateAtoms, ATOM, 
                iVarArrayElementCount(tTank->tank.vaRotateAtoms)-1 );
    aPCur = PVAI( tTank->tank.vaRotateAtoms, ATOM, 0 );
    while ( aPCur <= aPLast ) {
        MatrixTimesVector( vNew, mTransform, vAtomPosition(*aPCur) );
        AtomSetPosition( *aPCur, vNew );
        aPCur++;
    }

    tTank->tank.iXStart = tTank->tank.iXStop;
    tTank->tank.iYStart = tTank->tank.iYStop;

    TankFastRedisplayUnit( tTank );

}






/*
 *      zActionDrawRubberBandOrBox
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Draw the rubber band or rubber box depending
 *      on the value in (iSelectionType).
 */
static void
zActionDrawRubberBandOrBox( TANK tTank )
{
    switch( tTank->tank.iSelectionType ) {
        case TANK_SELECT_ATOM_OR_CHAIN:
            TankToggleDrawRubberBand(tTank);
            break;
        case TANK_SELECT_BOX:
            TankToggleDrawRubberBox(tTank);
            break;
        default:
            DFATAL(( "Illegal iSelectionType" ));
            break;
    }
}




/*
 *      zActionHandlePointerDown
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Handle pointer down event.
 */
static void
zActionHandlePointerDown( TANK tTank, int iX, int iY )
{
ATOM            aAtom;

    tTank->tank.bRubberOn = FALSE;
    tTank->tank.iXStart = iX;
    tTank->tank.iYStart = iY;

    switch ( tTank->tank.iTankState ) {
        case TANKDRAGROTATE:

            DisplayerAccumulateUpdates();

            ActionPickAtom( tTank, iX, iY, &aAtom, &iX, &iY );
            if ( aAtom != NULL ) {
                tTank->tank.aRotateCenter = aAtom;
            }
            TankSetFlags( tTank, TANKDRAWCIRCLE );
            break;
        case TANKTWIST:
            DisplayerAccumulateUpdates();
            break;
        default:
            break;
    }
}





/*
 *      zActionHandlePointerMove
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Handle pointer move event.
 */
static void
zActionHandlePointerMove( TANK tTank, int iX, int iY )
{

    if ( tTank->tank.bRubberOn ) 
        zActionDrawRubberBandOrBox(tTank);
    tTank->tank.iXStop = iX;
    tTank->tank.iYStop = iY;

    switch( tTank->tank.iTankState ) {
        case TANKDRAGROTATE:
            if ( tTank->tank.iButtonState == BPOINT ) {
                zActionDrag( tTank );
            } else if ( tTank->tank.iButtonState == BSHIFT_POINT ) {
                zActionRotateAtoms(tTank);
            }
            goto DONE;
            break;
        case TANKTWIST:
            zActionTwistTorsions( tTank );
            goto DONE;
            break;
    }
    zActionDrawRubberBandOrBox(tTank);

DONE:
    return;
}





/*
 *      zActionHandlePointerUp
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Handle pointer up event.
 */
static void
zActionHandlePointerUp( TANK tTank, int iX, int iY )
{

    if ( tTank->tank.bRubberOn ) 
        zActionDrawRubberBandOrBox(tTank);
    tTank->tank.iXStop = iX;
    tTank->tank.iYStop = iY;

        /* Now figure out what to do with the line segment that was drawn */

    switch ( tTank->tank.iTankState ) {
        case TANKDRAW:
            ActionUseLineSegmentToDraw( tTank );
            break;
        case TANKERASE:
            ActionUseLineSegmentToErase( tTank );
            break;
        case TANKSELECT:

            if ( tTank->tank.iButtonState == BSHIFT_POINT ) {
                ActionDeSelect( tTank );
            } else {
                ActionSelect( tTank );
            }
            break;
        case TANKTWIST:
            TankRedisplayUnit( tTank );
            DisplayerReleaseUpdates();
            break;
        case TANKDRAGROTATE:
            TankResetFlags( tTank, TANKDRAWCIRCLE );
            TankRedisplayUnit( tTank );
            DisplayerReleaseUpdates();
            break;
        default:
            DFATAL(( "Illegal TANK state" ));
            break;
    }

}




/*
 *      zActionSnapToAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the X,Y position of an ATOM if the mouse
 *      is near the ATOM, otherwise return the X,Y position
 *      of the mouse.
 */
static BOOL
zbActionSnapToAtom( TANK tTank, ATOM *aPAtom, int *iPX, int *iPY )
{
int             iAX, iAY;
BOOL            bGotOne;

                /* If the pointer is near an ATOM then return the */
                /* coordinates of that ATOM, otherwise return the */
                /* coordinates of the mouse.  This causes it to SNAP */
                /* to ATOMs */


    bGotOne = FALSE;
    *aPAtom = NULL;
    ActionPickAtom( tTank, *iPX, *iPY, aPAtom, &iAX, &iAY );
    if ( *aPAtom != NULL ) {
        bGotOne = TRUE;
        *iPX = iAX;
        *iPY = iAY;
    }
    return(bGotOne);
}




/*
 *      ziActionTranslateButtonsToState
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Translate the combination of buttons currently
 *      being held down to a button state.
 *      bDown is TRUE if the new button is being pressed down
 *      or FALSE if it is being released.  iButtonMask is
 *      the mask of the buttons before the event occured.
 *      iButton is the number of the button pressed, and iOldState
 *      is the old state of the TANK before the button was pressed.
 *
 *      Return the new TANK state.
 *
 */
static int
ziActionTranslateButtonsToState( BOOL bDown, int *iPButtonMask, 
                                        int iButton, int iOldState )
{
int             iNewMask;

    iNewMask = iButton;

    if ( bDown )        (*iPButtonMask) |= iNewMask;
    else                (*iPButtonMask) &= ~iNewMask;

    if ( *iPButtonMask == 0 ) return(BNONE);
    switch (*iPButtonMask) {
        case 0:
            return(BNONE);
        case ButtonA:
            return(BPOINT);
        case ButtonB:
            return(BROTATE);
        case ButtonC:
            return(BTRANSXY);
        case ButtonB | ButtonC:
            return(BSCALE);
        case ButtonA | ButtonD:
            return(BSHIFT_POINT);
        default:
            return(iOldState);
    }
    return(iOldState);          /* for lint */
}

            
            

/*
 *      ActionPointDown
 *
 *      Author: Christian Schafmeister (1991)
 *      ActionPointMotion
 *      ActionPointUp
 *
 *      These three routines translate mouse motions into
 *      pointing.  Pointing does different things depending
 *      on the state of the TANK.  The state of the
 *      TANK can be changed by calling other actions.
 */

static void
ActionPointDown( TANK tTank, int iX, int iY, Time tTime )
{
ATOM            aAtom;
BOOL            bOnAtom;
Time            tElapsed;

    if ( tTank->tank.uUnit == NULL ) 
        return;

                /* Check what the time is for double clicking */

    tTank->tank.tTimeSinceLastClick = tTime - tTank->tank.tLastClickMSec;
    tTank->tank.tLastClickMSec = tTime;

    MESSAGE(( "ActionPointDown\n" ));

    bOnAtom = FALSE;
    if ( !( tTank->tank.iTankState == TANKTWIST ||
            tTank->tank.iTankState == TANKDRAGROTATE ) ) {
        bOnAtom = zbActionSnapToAtom( tTank, &aAtom, &iX, &iY );
    }
    tTank->tank.iSelectionType = TANK_SELECT_ATOM_OR_CHAIN;
    if ( tTank->tank.iTankState == TANKSELECT ) {
        if ( bOnAtom ) {
            tTank->tank.iSelectionType = TANK_SELECT_ATOM_OR_CHAIN;
        } else {
            tTank->tank.iSelectionType = TANK_SELECT_BOX;
        }
    }

                /* Check if this is a multi-click */

    tElapsed = tTank->tank.tTimeSinceLastClick;
    if ( tElapsed > SELECT_DELAY ) {
        tTank->tank.iTankSelectLastState = TANKSELECTNOTHING;
        MESSAGE(( "New click\n" ));
    } else {
        MESSAGE(( "Multi click\n" ));
    }
    tTank->tank.iTankSelectLastState++;

    zActionHandlePointerDown( tTank, iX, iY );
}


static void
ActionPointMotion( TANK tTank, int iX, int iY )
{
ATOM            aAtom;

    if ( tTank->tank.uUnit == NULL ) 
        return;

    MESSAGE(( "ActionPointMotion\n" ));

    if ( !( tTank->tank.iTankState == TANKTWIST ||
            tTank->tank.iTankState == TANKDRAGROTATE ||
             ( tTank->tank.iTankState == TANKSELECT  &&
                tTank->tank.iSelectionType == TANK_SELECT_BOX ) ) ) {
        zbActionSnapToAtom( tTank, &aAtom, &iX, &iY );
    }

    zActionHandlePointerMove( tTank, iX, iY );

}





static void
ActionPointUp( TANK tTank, int iX, int iY )
{
ATOM            aAtom;

    if ( tTank->tank.uUnit == NULL ) 
        return;

    MESSAGE(( "ActionPointUp\n" ));

    if ( !( tTank->tank.iTankState == TANKTWIST ||
            tTank->tank.iTankState == TANKDRAGROTATE ||
             ( tTank->tank.iTankState == TANKSELECT &&
                tTank->tank.iSelectionType == TANK_SELECT_BOX ) ) ) {
        zbActionSnapToAtom( tTank, &aAtom, &iX, &iY );
    }

    zActionHandlePointerUp( tTank, iX, iY );

}




/*
 *      ActionRotateDown
 *      ActionRotateMotion
 *      ActionRotateUp
 *
 *      These three routines take mouse events and translate them into
 *      TANK rotation.
 */

static void
ActionRotateDown( TANK tTank, int iX, int iY )
{
    TankSetFlags( tTank, TANKDRAWCIRCLE );
    MESSAGE(( "ActionRotateDown\n" ));
    tTank->tank.iActionX = iX;
    tTank->tank.iActionY = iY;
}


static void
ActionRotateMotion( TANK tTank, int iX, int iY )
{
NUMBR           nX1, nY1, nX2, nY2;
VECTOR          vWorld;
int             iWidth, iHeight;


    TankSetFlags( tTank, TANKDRAWCIRCLE );
    bTankPointInCircle( tTank, iX, iY ); 
    MESSAGE(( "ActionRotateMotion\n" ));

    iWidth = iX3dScreenWidth(tTank->tank.x3dEngine);
    iHeight = iX3dScreenHeight(tTank->tank.x3dEngine);
    nX1 = ( tTank->tank.iActionX - (iWidth>>1) ); 
    nY1 = -( tTank->tank.iActionY - (iHeight>>1) );
    nX2 = ( iX - (iWidth>>1) );
    nY2 = -( iY - (iHeight>>1));

    if ( bTankPointInCircle( tTank, iX, iY ) ) {
        VectorDef( &vWorld, 0.0, 0.0, 0.0 );
        ActionRotateXY( tTank, &vWorld, nX1, nY1, nX2, nY2 );
    } else {
        VectorDef( &vWorld, 0.0, 0.0, 0.0 );
        ActionRotateZ( tTank, &vWorld, nX1, nY1, nX2, nY2 );
    }
    tTank->tank.iActionX = iX; 
    tTank->tank.iActionY = iY; 

    X3dBuildTransform( tTank->tank.x3dEngine );
    TankDisplay( tTank, TRUE );
    
}


static void
ActionRotateUp( TANK tTank)
{
    TankResetFlags( tTank, TANKDRAWCIRCLE );
    tTank->tank.bNeedToRedraw = TRUE;
    MESSAGE(( "ActionRotateUp\n" ));
    TankRedisplayUnit( tTank );
}





/*
 *      ActionScaleDown
 *      ActionScaleMotion
 *      ActionScaleUp
 *
 *      These three routines translate mouse motions into
 *      TANK scaling.
 */


static void
ActionScaleDown( TANK tTank, int iY )
{
    MESSAGE(( "ActionScaleDown\n" ));
    tTank->tank.iActionY = iY;
}


static void
ActionScaleMotion( TANK tTank, int iY )
{
NUMBR           nScale;

    MESSAGE(( "ActionScaleMotion\n" ));
    nScale = nX3dScale( tTank->tank.x3dEngine );
    nScale += ( tTank->tank.iActionY - iY )*SCALE_SPEED;
    if ( nScale < 0.2 ) 
        nScale = 0.2;
    X3dSetScale( tTank->tank.x3dEngine, nScale );
    
    MESSAGE(( "Old Y = %d  New Y = %d  scale = %f\n",
                tTank->tank.iActionY, iY, nScale ));
                
    tTank->tank.iActionY = iY; 

/* TODO:  Change this so scale does not become 0 */

    if ( nX3dScale(tTank->tank.x3dEngine) < 0.0 ) 
                X3dSetScale( tTank->tank.x3dEngine, 0.0 );

    X3dBuildTransform( tTank->tank.x3dEngine );
    TankDisplay( tTank, TRUE );
}


static void
ActionScaleUp( TANK tTank )
{
    MESSAGE(( "ActionScaleUp\n" ));
    TankRedisplayUnit( tTank );
}



/*
 *      ActionTransXYDown
 *      ActionTransXYMotion
 *      ActionTransXYUp
 *
 *      These three routines translate mouse motions into
 *      TANK scaling.
 */

static void
ActionTransXYDown( TANK tTank, int iX, int iY )
{
    MESSAGE(( "ActionTransXYDown\n" ));
    MESSAGE(( "TransXY iX,iY = %d, %d\n", iX, iY ));
    tTank->tank.iActionX = iX;
    tTank->tank.iActionY = iY;
}



static void
ActionTransXYMotion( TANK tTank, int iX, int iY )
{
GVECTOR         gvScreen, gvA, gvB, gvDiff, gvPos;

    MESSAGE(( "ActionTransXYMotion\n" ));

    MESSAGE(( "TransXY iX,iY = %d, %d\n", iX, iY ));

    nGVX(gvScreen) = iX;
    nGVY(gvScreen) = iY;
    nGVZ(gvScreen) = 0;         /* for purify & Dec */
    X3dScreenToView( tTank->tank.x3dEngine, &gvScreen, &gvA );

    nGVX(gvScreen) = tTank->tank.iActionX;
    nGVY(gvScreen) = tTank->tank.iActionY;
    X3dScreenToView( tTank->tank.x3dEngine, &gvScreen, &gvB );

    GVectorSubtract( gvDiff, gvB, gvA );
    GVectorTimesScalar( gvDiff, gvDiff, TRANSLATE_SPEED );

    MESSAGE(( "Translating by: %lf, %lf, %lf\n",
                nGVX(gvDiff), nGVY(gvDiff), nGVZ(gvDiff) ));

    gvPos = gvX3dCenterOfView( tTank->tank.x3dEngine );
    GVectorAdd( gvPos, gvPos, gvDiff );
    X3dSetCenterOfView( tTank->tank.x3dEngine, gvPos );

    tTank->tank.iActionX = iX;
    tTank->tank.iActionY = iY;
    X3dBuildTransform( tTank->tank.x3dEngine );
    TankDisplay( tTank, TRUE );
}


static void
ActionTransXYUp( TANK tTank )
{
    MESSAGE(( "ActionTransXYUp\n" ));
    TankRedisplayUnit( tTank );
}



/*
 *-------------------------------------------------------
 *
 *      Action procedures through which
 *      all xAction functions are entered.
 *
 */


/*
 *      ActionButtonDown
 *      ActionButtonMotion
 *      ActionButtonUp
 *
 *      These three actions translate button events into
 *      the rotate/scale/translate/point actions.
 */

static void
ActionButtonDown( Widget wWidget, int iButton, int iX, int iY, Time tTime )
{
TANK            tTank;
int             iNewState;


    tTank = (TANK)wWidget;
    if ( tTank->tank.uUnit == NULL ) 
        return;

    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);

    MESSAGE(( "ActionButtonDown\n" ));

    iNewState = ziActionTranslateButtonsToState( TRUE,
                        &(tTank->tank.iRawButtonState),
                        iButton, tTank->tank.iButtonState );

                /* Execute an Up action for the old state */

    switch( tTank->tank.iButtonState ) {
        case BNONE:
            break;
        case BPOINT:
        case BSHIFT_POINT:
            ActionPointUp( tTank, iX, iY );
            break;
        case BROTATE:
            ActionRotateUp( tTank );
            break;
        case BTRANSXY:
            ActionTransXYUp( tTank );
            break;
        case BSCALE:
            ActionScaleUp( tTank );
            break;
        default:
            DFATAL(( "Illegal button state for TANK: %d\n",
                        tTank->tank.iButtonState ));
            break;
    }

                /* Then send a Down Action for the new state */

    tTank->tank.iButtonState = iNewState;

    switch( tTank->tank.iButtonState ) {
        case BNONE:
            break;
        case BPOINT:
        case BSHIFT_POINT:
            ActionPointDown( tTank, iX, iY, tTime );
            break;
        case BROTATE:
            ActionRotateDown( tTank, iX, iY );
            break;
        case BTRANSXY:
            ActionTransXYDown( tTank, iX, iY );
            break;
        case BSCALE:
            ActionScaleDown( tTank, iY );
            break;
        default:
            DFATAL(( "Illegal button state for TANK: %d\n",
                        tTank->tank.iButtonState ));
            break;
    }
    PopCurrentPrintSink();
}


static void
ActionButtonMotion( Widget wWidget, int iX, int iY )
{
TANK            tTank;

    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);

    tTank = (TANK)wWidget;
    
    switch ( tTank->tank.iButtonState ) {
        case BNONE:
            break;
        case BPOINT:
        case BSHIFT_POINT:
            ActionPointMotion( tTank, iX, iY );
            break;
        case BROTATE:
            ActionRotateMotion( tTank, iX, iY );
            break;
        case BTRANSXY:
            ActionTransXYMotion( tTank, iX, iY );
            break;
        case BSCALE:
            ActionScaleMotion( tTank, iY );
            break;
        default:
            DFATAL(( "Illegal button state for TANK: %d\n",
                        tTank->tank.iButtonState ));
            break;
    }
    PopCurrentPrintSink();
}





static void
ActionButtonUp( Widget wWidget, int iButton, int iX, int iY, Time tTime )
{
TANK            tTank;
int             iNewState;


    tTank = (TANK)wWidget;

    if ( tTank->tank.uUnit == NULL ) 
        return;

    PushCurrentPrintSink(tTank->tank.iPrintSink);

                /* Send an Up Action for the old state */

    switch ( tTank->tank.iButtonState ) {
        case BNONE:
            break;
        case BPOINT:
        case BSHIFT_POINT:
            ActionPointUp( tTank, iX, iY );
            break;
        case BROTATE:
            ActionRotateUp( tTank );
            break;
        case BTRANSXY:
            ActionTransXYUp( tTank );
            break;
        case BSCALE:
            ActionScaleUp( tTank );
            break;
        default:
            DFATAL(( "Illegal button state for TANK: %d\n",
                        tTank->tank.iButtonState ));
            break;
    }

                /* Get the new state and send a Down Action for it */

    iNewState = ziActionTranslateButtonsToState( FALSE,
                        &(tTank->tank.iRawButtonState),
                        iButton,
                        tTank->tank.iButtonState );
    tTank->tank.iButtonState = iNewState;

    switch( tTank->tank.iButtonState ) {
        case BNONE:
            break;
        case BPOINT:
        case BSHIFT_POINT:
            ActionPointDown( tTank, iX, iY, tTime );
            break;
        case BROTATE:
            ActionRotateDown( tTank, iX, iY );
            break;
        case BTRANSXY:
            ActionTransXYDown( tTank, iX, iY );
            break;
        case BSCALE:
            ActionScaleDown( tTank, iY );
            break;
        default:
            DFATAL(( "Illegal button state for TANK: %d\n",
                        tTank->tank.iButtonState ));
            break;
    }
    PopCurrentPrintSink();
}



 
/*
 ***************************************************************
 ***************************************************************
 ***************************************************************
 *
 *      ACTION procedures
 *      
 *      Called from translation tables
 * 
 ***************************************************************
 ***************************************************************
 ***************************************************************
 */




XtActionProc
xtapActionCenterUnit( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);

    MESSAGE(( "ActionCenterUnit\n" ));
    X3dBuildCenteredScaledTransform(
                        ((TANK)wWidget)->tank.x3dEngine,
                        ((TANK)wWidget)->tank.uUnit );
    TankDisplay( (TANK)wWidget, FALSE );

    PopCurrentPrintSink();

    return NULL;
}




XtActionProc    
xtapActionSetCancelHit( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{

    MESSAGE(( "ActionSetCancelHit\n" ));
    BasicsSetInterrupt();

    return NULL;
}


/*
 ****************************************************
 *
 *      There are three sets of three button actions.
 *      These can be connected to buttons or keyboard keys.
 *
 *      ButtonA is set up to be used as the pointer button.
 *      ButtonB is set up to be used as the rotate button.
 *      ButtonC is set up to be used as the translate button.
 *      
 *      ButtonB and ButtonC together are set up to scale.
 */


XtActionProc    
xtapActionMouseMotion( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


        /* Change the event to look like a mouse event */

    ActionButtonMotion( wWidget, xePEvent->xbutton.x, xePEvent->xbutton.y );

    PopCurrentPrintSink();

    return NULL;
}


/*
 *      ButtonA action callbacks.
 */

XtActionProc    
xtapActionButtonADown( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);

    MESSAGE(( "ActionButtonADown\n" ));

        /* Change the event to look like a mouse event */

    xePEvent->xbutton.button = Button1;
    ActionButtonDown( wWidget, ButtonA,
                                xePEvent->xbutton.x, xePEvent->xbutton.y,
                                xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}



XtActionProc    
xtapActionButtonAUp( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


    MESSAGE(( "ActionButtonAUp\n" ));

        /* Change the event to look like a mouse event */

    ActionButtonUp( wWidget, ButtonA,
                        xePEvent->xbutton.x, xePEvent->xbutton.y,
                        xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}



/*
 *      ButtonB action callbacks.
 */


XtActionProc    
xtapActionButtonBDown( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


    MESSAGE(( "ActionButtonBDown\n" ));

        /* Change the event to look like a mouse event */

    ActionButtonDown( wWidget, ButtonB,
                                xePEvent->xbutton.x, xePEvent->xbutton.y,
                                xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}




XtActionProc    
xtapActionButtonBUp( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


    MESSAGE(( "ActionButtonBUp\n" ));

        /* Change the event to look like a mouse event */

    ActionButtonUp( wWidget, ButtonB,
                        xePEvent->xbutton.x, xePEvent->xbutton.y,
                        xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}



/*
 *      ButtonC action callbacks.
 */

XtActionProc    
xtapActionButtonCDown( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


    MESSAGE(( "ActionButtonCDown\n" ));

        /* Change the event to look like a mouse event */

    ActionButtonDown( wWidget, ButtonC,
                                xePEvent->xbutton.x, xePEvent->xbutton.y,
                                xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}



XtActionProc    
xtapActionButtonCUp( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);



        /* Change the event to look like a mouse event */

    ActionButtonUp( wWidget, ButtonC,
                        xePEvent->xbutton.x, xePEvent->xbutton.y,
                        xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}



/*
 *      ButtonD action callbacks.
 */

XtActionProc    
xtapActionButtonDDown( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


    MESSAGE(( "ActionButtonDDown\n" ));

        /* Change the event to look like a mouse event */

    ActionButtonDown( wWidget, ButtonD,
                                xePEvent->xbutton.x, xePEvent->xbutton.y,
                                xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}



XtActionProc    
xtapActionButtonDUp( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


    MESSAGE(( "ActionButtonDUp\n" ));

        /* Change the event to look like a mouse event */

    ActionButtonUp( wWidget, ButtonD,
                        xePEvent->xbutton.x, xePEvent->xbutton.y,
                        xePEvent->xbutton.time );
    PopCurrentPrintSink();
 
    return NULL;
}





/*
 *      ButtonE action callbacks.
 */

XtActionProc    
xtapActionButtonEDown( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


    MESSAGE(( "ActionButtonEDown\n" ));

        /* Change the event to look like a mouse event */

    ActionButtonDown( wWidget, ButtonE,
                                xePEvent->xbutton.x, xePEvent->xbutton.y,
                                xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}



XtActionProc    
xtapActionButtonEUp( Widget wWidget, XEvent *xePEvent, 
                String *sPArgs, Cardinal *cPNum )
{
    PushCurrentPrintSink(((TANK)(wWidget))->tank.iPrintSink);


    MESSAGE(( "ActionButtonEUp\n" ));

        /* Change the event to look like a mouse event */

    ActionButtonUp( wWidget, ButtonE,
                        xePEvent->xbutton.x, xePEvent->xbutton.y,
                        xePEvent->xbutton.time );

    PopCurrentPrintSink();

    return NULL;
}







/*
 *      ActionInit
 *
 *      Initialize the Action functions.
 */
void
ActionInit( )
{

}

