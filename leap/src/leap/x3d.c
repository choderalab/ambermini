/*
 *      File: x3d.c
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
 *		Code to do 3D transformations for X-Windows.
 *
 *		A word about coordinate systems:
 *
 *		World coordinates are the coordinate system of the atoms.
 *		View coordinates are the coordinate system of the viewer.
 *			The viewer is at the origin and is looking down the
 *			negative Z axis with up as the positive Y axis and
 *			right is the positive X axis.
 *		Screen coordinates are 2D coordinates of the screen.
 *			When points are converted from screen coordinates to
 *			view coordinates, their Z value is the average of
 *			the front and back clipping planes.
 *
 */

#define		DEBUG_TRANSFORM	1

#include	"basics.h"

#include	"vector.h"
#include        "threed.h"

#include	"varArray.h"

#include	"classes.h"

#include	"x3d.h"



/*
 *---------------------------------------------------
 *
 *	Private routines
 *
 */

#define	RESCALE_MULTIPLIER	3.0

void
GMatrixTimesVector( GVECTOR *result, GMATRIX m, GVECTOR *b )
{
	result->nX = m[0][0] * b->nX + 
		     m[1][0] * b->nY + 
		     m[2][0] * b->nZ + 
		     m[3][0];
	result->nY = m[0][1] * b->nX + 
		     m[1][1] * b->nY + 
		     m[2][1] * b->nZ + 
		     m[3][1];
	result->nZ = m[0][2] * b->nX + 
		     m[1][2] * b->nY + 
		     m[2][2] * b->nZ + 
		     m[3][2];
}

/*
 *      zX3dClipToPlane
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Clip the line defined by *gvPA, *gvPB  to the plane with
 *	the Z value (nClip).  Return the intersection in (*gvPClip).
 */
void
zX3dClipToPlane( GVECTOR *gvPClip, GVECTOR *gvPA, GVECTOR *gvPB, NUMBR nClip )
{
NUMBR           nT;
GVECTOR         gvC;

                /* Solve for the intersection of a line and a plane */

    nT = ( nClip - nGVZ(*gvPA) ) / ( nGVZ(*gvPB) - nGVZ(*gvPA) );
    GVectorDef( gvC,
                nT*(nGVX(*gvPB) - nGVX(*gvPA)) + nGVX(*gvPA),
                nT*(nGVY(*gvPB) - nGVY(*gvPA)) + nGVY(*gvPA),
                nClip );
    *gvPClip = gvC;
    GVectorSetClipStatus( *gvPClip, CLIPMIDDLE );
}




/*
 *------------------------------------------------------
 *
 *	Public routines
 *
 */


/*
 *	x3dX3dCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create an X3DENGINE.
 */
X3DENGINE
x3dX3dCreate()
{
X3DENGINE	x3dEngine;

    MALLOC( x3dEngine, X3DENGINE, sizeof(X3DENGINEt) );

    return(x3dEngine);
}




/*
 *	X3dDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the X3DENGINE.
 */
void
X3dDestroy( X3DENGINE *x3dPEngine )
{
    FREE( *x3dPEngine );
    *x3dPEngine = NULL;
}



/*
 *	X3dInit
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Initialize the matrices that will be used to
 *	build the Transformation matrix.
 */
void
X3dInit( X3DENGINE x3dEngine )
{
    GVectorDef( x3dEngine->gvCenterOfRotation, 0.0, 0.0, 0.0 );
    x3dEngine->nViewCenter = (DEFAULTFARCLIPPLANE+DEFAULTNEARCLIPPLANE)/2.0;
    GVectorDef( x3dEngine->gvCenterOfView, 
    	0.0, 0.0, x3dEngine->nViewCenter );
    GMatrixIdentity( x3dEngine->gmRotate );
    x3dEngine->nViewFar      = DEFAULTFARCLIPPLANE;
    x3dEngine->nViewNear     = DEFAULTNEARCLIPPLANE;
    x3dEngine->nViewCenter = (x3dEngine->nViewFar+x3dEngine->nViewNear)/2.0;
    GVectorDef( x3dEngine->gvCenterOfView, 
    	0.0, 0.0, x3dEngine->nViewCenter );
    x3dEngine->nViewWidth    = DEFAULTVIEWWIDTH;
    x3dEngine->nViewHeight   = DEFAULTVIEWHEIGHT;
    x3dEngine->nViewScreenDistance = DEFAULTSCREENDISTANCE;
    x3dEngine->nScale        = 1.0;
    x3dEngine->iScreenWidth = 500;
    x3dEngine->iScreenHeight = 500;
    x3dEngine->iUseWidth = 500;
    x3dEngine->iUseHeight = 500;

}



/*
 *      X3dBuildTransform
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the Transform matrix from the parameters
 *      defined for the X3DENGINE. 
 *
 */
void
X3dBuildTransform( X3DENGINE x3dEngine )
{
GMATRIX         gmT, gmU, gmV;

		/* Build the prescaling translation matrix */

    GMatrixTranslate( gmU, 
                        -nGVX(x3dEngine->gvCenterOfRotation),
                        -nGVY(x3dEngine->gvCenterOfRotation),
                        -nGVZ(x3dEngine->gvCenterOfRotation) );


                /* Multiply in the rotation matrix */

    GMatrixMultiply( gmT, x3dEngine->gmRotate, gmU );

                /* Build and multiply in the translation */
                /* matrix to move the system to the point */
                /* in View space that will coorespond to the */
                /* center of the window */

    GMatrixTranslate( gmU, nGVX(x3dEngine->gvCenterOfView),
    			   nGVY(x3dEngine->gvCenterOfView),
			   nGVZ(x3dEngine->gvCenterOfView) );
			   
    GMatrixMultiply( gmV, gmU, gmT );

		/* Save World to View matrix */

    GMatrixMultiply( x3dEngine->gmWorldToView, gmT, gmV );

                /* Build the matrix to convert from View to */
                /* unprojected Screen Coordinates */

    GMatrixDiagonal( gmT,
        ((NUMBR)(x3dEngine->iUseWidth))/x3dEngine->nViewWidth,
        ((NUMBR)(x3dEngine->iUseHeight))/x3dEngine->nViewHeight,
        1.0, 1.0 );

    GMatrixMultiply( x3dEngine->gmWorldToPartialScreen, gmT, gmV );


                /* Now make a note that the inverse of */
                /* this matrix has not been calculated */
                /* It will be calculated when it is needed */

    x3dEngine->bRecalculateForward = FALSE;

    x3dEngine->nScreenXOffset = ((NUMBR)(x3dEngine->iScreenWidth))/2.0;
    x3dEngine->nScreenYOffset = ((NUMBR)(x3dEngine->iScreenHeight))/2.0;

#ifdef	DEBUG_TRANSFORM
    MESSAGE(( "nViewWidth, nViewHeight = %f, %f\n",
		x3dEngine->nViewWidth,
		x3dEngine->nViewHeight ));
    MESSAGE(( "iUseWidth, iUseHeight = %d, %d\n", 
		x3dEngine->iUseWidth, x3dEngine->iUseHeight ));
    MESSAGE(( "To world origin: %lf, %lf, %lf\n",
    		-nGVX(x3dEngine->gvCenterOfRotation),
    		-nGVY(x3dEngine->gvCenterOfRotation),
    		-nGVZ(x3dEngine->gvCenterOfRotation) ));
    MESSAGE(( "Rotation matrix\n" ));
    GMatrixPrint(x3dEngine->gmRotate );
    MESSAGE(( "To View origin: %lf, %lf, %lf\n",
    		nGVX(x3dEngine->gvCenterOfView),
    		nGVY(x3dEngine->gvCenterOfView),
    		nGVZ(x3dEngine->gvCenterOfView) ));
    MESSAGE(( "WorldToPartialScreen matrix: \n" ));
    GMatrixPrint(x3dEngine->gmWorldToPartialScreen);
    MESSAGE(( "Scale = %lf\n", x3dEngine->nScale ));
#endif

}



/*
 *      X3dBuildInverseTransform
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the inverse transform to transform from
 *      Screen coordinates to World coordinates by
 *      applying the inverse of each individual matrix
 *      in 'X3dBuildTransform' in reverse order.
 */
void
X3dBuildInverseTransform( X3DENGINE x3dEngine )
{
GMATRIX         gmT, gmU, gmV;

        /* If it has already been calculated then return */


    x3dEngine->bRecalculateBackward = FALSE;


        /* Build the Inverse of the matrix for scaling the View */
        /* to the screen coordinates */

    GMatrixDiagonal( gmT,
        1.0/(((NUMBR)(x3dEngine->iUseWidth))/x3dEngine->nViewWidth),
        1.0/(((NUMBR)(x3dEngine->iUseHeight))/x3dEngine->nViewHeight),
        1.0, 1.0 );


	/* Translate from view origin to world origin */

    GMatrixTranslate( gmU, -nGVX(x3dEngine->gvCenterOfView),
    			   -nGVY(x3dEngine->gvCenterOfView),
			   -nGVZ(x3dEngine->gvCenterOfView) );

    GMatrixMultiply( gmV, gmU, gmT );

    GMatrixCopy( x3dEngine->gmPartialScreenToView, gmV );


	/* Get the inverse of the rotation matrix */

    GMatrixTranspose( gmT, x3dEngine->gmRotate );

	/* Inverse of translate to view origin matrix */

    GMatrixTranslate( gmV, 
                        nGVX(x3dEngine->gvCenterOfRotation),
                        nGVY(x3dEngine->gvCenterOfRotation),
                        nGVZ(x3dEngine->gvCenterOfRotation) );

    GMatrixMultiply( x3dEngine->gmViewToWorld, gmV, gmT );
    GMatrixMultiply( x3dEngine->gmPartialScreenToWorld,
				x3dEngine->gmViewToWorld,
				x3dEngine->gmPartialScreenToView );


        /* Note that the ToWorld coordinate transformation matrix */
        /* has now been built */

#ifdef	DEBUG_TRANSFORM
    MESSAGE(( "PartialScreenToView matrix\n" ));
    GMatrixPrint(x3dEngine->gmPartialScreenToView );
    MESSAGE(( "ViewToWorld matrix\n" ));
    GMatrixPrint(x3dEngine->gmViewToWorld );
    MESSAGE(( "PartialScreenToWorld matrix\n" ));
    GMatrixPrint(x3dEngine->gmPartialScreenToWorld );
#endif


}


/*
 *--------------------------------------------------------------
 *
 *      Three Dimension routines
 *
 */


/*
 *      X3dWorldToScreen
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Convert a point from world coordinates to screen.
 */
void
X3dWorldToScreen( X3DENGINE x3dEngine, GVECTOR *gvPPoint, GVECTOR *gvPNew )
{
    if ( x3dEngine->bRecalculateForward ) 
	X3dBuildTransform( x3dEngine );

    GMatrixTimesVector( gvPNew, x3dEngine->gmWorldToPartialScreen, gvPPoint );
    GVectorSetClipStatus( *gvPNew, CLIPMIDDLE );
    if ( nGVZ(*gvPNew) > x3dEngine->nViewNear ) 
        GVectorSetClipStatus( *gvPNew, CLIPFRONT );
    else if ( nGVZ(*gvPNew) < x3dEngine->nViewFar )
        GVectorSetClipStatus( *gvPNew, CLIPBACK );

                /* If PERSPECTIVE is on then do the perspective */
                /* calculation */

    if ( x3dEngine->bPerspectiveOn ) {
        nGVX(*gvPNew) = nGVX(*gvPNew)*(x3dEngine->nViewScreenDistance)/
                                nGVZ(*gvPNew);
        nGVY(*gvPNew) = nGVY(*gvPNew)*(x3dEngine->nViewScreenDistance)/
                                nGVZ(*gvPNew);
    }
    
    nGVX(*gvPNew) *= x3dEngine->nScale;
    nGVY(*gvPNew) *= x3dEngine->nScale;

                /* Translate the origin of the view space into the */
                /* center of the screen */

    nGVX(*gvPNew) += x3dEngine->nScreenXOffset;
    nGVY(*gvPNew) = x3dEngine->nScreenYOffset - nGVY(*gvPNew);

}



/*
 *      X3dScreenToWorld
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Convert a point from screen coordinates to world.
 */
void
X3dScreenToWorld( X3DENGINE x3dEngine, GVECTOR *gvPPoint, GVECTOR *gvPNew )
{
GVECTOR         gvTemp;

    if ( x3dEngine->bRecalculateBackward ) 
        X3dBuildInverseTransform( x3dEngine );

    nGVX(gvTemp) = nGVX(*gvPPoint) - x3dEngine->nScreenXOffset;
    nGVY(gvTemp) = x3dEngine->nScreenYOffset - nGVY(*gvPPoint) ;
    
    nGVX(gvTemp) /= x3dEngine->nScale;
    nGVY(gvTemp) /= x3dEngine->nScale;

    nGVZ(gvTemp) = x3dEngine->nViewCenter;

    if ( x3dEngine->bPerspectiveOn ) {
        nGVX(gvTemp) = nGVX(gvTemp)*nGVZ(gvTemp)/
			(x3dEngine->nViewScreenDistance);
	nGVY(gvTemp) = nGVY(gvTemp)*nGVZ(gvTemp)/
			(x3dEngine->nViewScreenDistance);
    }

    GMatrixTimesVector( gvPNew, x3dEngine->gmPartialScreenToWorld, &gvTemp );

}




/*
 *      X3dScreenToView
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Convert a point from screen coordinates to view coordinates.
 *	The Z coordinate of the view coordinate returned is halfway
 *	between the near and far clip planes.
 */
void
X3dScreenToView( X3DENGINE x3dEngine, GVECTOR *gvPPoint, GVECTOR *gvPNew )
{
GVECTOR         gvTemp;

    if ( x3dEngine->bRecalculateBackward ) 
        X3dBuildInverseTransform( x3dEngine );

    nGVX(gvTemp) = nGVX(*gvPPoint) - x3dEngine->nScreenXOffset;
    nGVY(gvTemp) = x3dEngine->nScreenYOffset - nGVY(*gvPPoint);
    
    nGVX(gvTemp) /= x3dEngine->nScale;
    nGVY(gvTemp) /= x3dEngine->nScale;
    
    nGVZ(gvTemp) = nGVZ(*gvPPoint) - x3dEngine->nViewScreenDistance;

    GMatrixTimesVector( gvPNew, x3dEngine->gmPartialScreenToView, &gvTemp );

}




/*
 *	X3dViewToWorld
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Convert view coordinates to world coordinates.
 */
void	
X3dViewToWorld( X3DENGINE x3dEngine, GVECTOR *gvPView, GVECTOR *gvPWorld )
{
    if ( !x3dEngine->bRecalculateBackward )
	X3dBuildInverseTransform( x3dEngine );

    GMatrixTimesVector( gvPWorld, x3dEngine->gmViewToWorld, gvPView );
}





/*
 *	X3dBuildCenteredScaledTransform
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build a transformation that centers and scales the molecule
 *	so that it fits the screen perfectly.
 */
void	
X3dBuildCenteredScaledTransform( X3DENGINE x3dEngine, UNIT uUnit )
{
VECTOR		vCenter;
LOOP		lAtoms;
double		dTopLeftX, dTopLeftY, dBottomRightX, dBottomRightY;
ATOM		aAtom;
GVECTOR		gvTemp, gvScreen;
double		dXWidth, dYWidth, dScale;
int		iCount;
BOOL		bPerspective;

		/* Save perspective mode and turn it off */


    if ( uUnit == NULL ) {
	DFATAL(( "NULL UNIT passed to X3dBuildCenteredScaledTransform" ));
    }

    bPerspective = x3dEngine->bPerspectiveOn;
    x3dEngine->bPerspectiveOn = FALSE;

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    iCount = 0;
    VectorDef( &vCenter, 0.0, 0.0, 0.0 );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
	if ( bAtomVisible(aAtom) ) {
	    vCenter = vVectorAdd( &vCenter, &vAtomPosition(aAtom) );
	    iCount++;
	}
    }
    if ( iCount != 0 ) {
    	vCenter = vVectorTimesScalar( &vCenter, 1.0/(double)iCount );
    }

    GVectorDef( x3dEngine->gvCenterOfRotation,
			dVX(&vCenter),
			dVY(&vCenter),
			dVZ(&vCenter) );

    MESSAGE(( "Using UNIT: %s\n", sContainerName(uUnit) ));
    MESSAGE(( "Geometric center @ %f, %f, %f\n",
		nGVX(x3dEngine->gvCenterOfRotation),
		nGVY(x3dEngine->gvCenterOfRotation),
		nGVZ(x3dEngine->gvCenterOfRotation) ));

    GVectorDef( x3dEngine->gvCenterOfView, 
    	0.0, 0.0, x3dEngine->nViewCenter );
    x3dEngine->nScale        = 1.0;

		/* Build the transform with scale of 1.0 to find */
		/* the bounding box of the UNIT */

    X3dBuildTransform( x3dEngine );
    MESSAGE(( "Screen offsets %f, %f\n", 
		x3dEngine->nScreenXOffset,
		x3dEngine->nScreenYOffset ));
    MESSAGE(( "Scale = %f\n", x3dEngine->nScale ));
    MESSAGE(( "ViewCenter = %f\n", x3dEngine->nViewCenter ));
    MESSAGE(( "Screen distance = %f\n", x3dEngine->nViewScreenDistance ));

		/* Find the bounding box of the visible ATOMs in the UNIT */

    dTopLeftX = 99999.0;
    dTopLeftY = 99999.0;
    dBottomRightX = -99999.0;
    dBottomRightY = -99999.0;
    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
	if ( bAtomVisible(aAtom) ) {
            GVectorDef( gvTemp, dVX(&vAtomPosition(aAtom)),
                      dVY(&vAtomPosition(aAtom)),
                      dVZ(&vAtomPosition(aAtom)) );
            X3dWorldToScreen( x3dEngine, &gvTemp, &gvScreen );
	    MESSAGE(( "Atom %s at view: %lf, %lf, %lf\n",
			sContainerName(aAtom),
			nGVX(gvScreen), nGVY(gvScreen) ));
	    if ( dTopLeftX > nGVX(gvScreen) ) 
			dTopLeftX = nGVX(gvScreen);
	    if ( dTopLeftY > nGVY(gvScreen) ) 
			dTopLeftY = nGVY(gvScreen);
	    if ( dBottomRightX < nGVX(gvScreen) ) 
			dBottomRightX = nGVX(gvScreen);
	    if ( dBottomRightY < nGVY(gvScreen) ) 
			dBottomRightY = nGVY(gvScreen);
	}
    }

    dXWidth = dBottomRightX - dTopLeftX;
    dYWidth = dBottomRightY - dTopLeftY;

    MESSAGE(( "Unscaled image size: %lf, %lf\n", dXWidth, dYWidth ));
    MESSAGE(( "Window size: %d, %d\n", 
		x3dEngine->iScreenWidth, x3dEngine->iScreenHeight ));

    MESSAGE(( "Checking if scaling in X direction fits Y direction\n" ));
    if ( dXWidth < VERYSMALL )
	dScale = 1;
    else
	dScale = x3dEngine->iScreenWidth / dXWidth;
	
    if ( dScale * dYWidth < x3dEngine->iScreenWidth ) {
	MESSAGE(( "Scaling X   scale = %lf\n", dScale ));
    } else if ( dYWidth > VERYSMALL ) {
	dScale = x3dEngine->iScreenHeight/dYWidth;
	MESSAGE(( "Scaling Y   scale = %lf\n", dScale ));
    }
    x3dEngine->nScale = dScale*RESCALE_MULTIPLIER;

		/* Remember whether perspective is being used */

    x3dEngine->bPerspectiveOn = bPerspective;
    X3dBuildTransform( x3dEngine );


}






/*
 *--------------------------------------------------------------------
 *
 *	Clip to back and front planes
 */



/*
 *      bX3dVisibleLineClip
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Take two points in Screen space and check if the
 *      line between them will be displayed in the window.
 *	If not then return FALSE, otherwise clip them and return
 *	the integer coordinates.
 */
BOOL
bX3dVisibleLineClip( X3DENGINE x3dEngine, int *iPX1, int *iPY1, 
	int *iPX2, int *iPY2, GVECTOR *gvP1, GVECTOR *gvP2 )
{
NUMBR           nPlane;
GVECTOR         gv1, gv2;

    gv1 = *gvP1;
    gv2 = *gvP2;

    if ( !bGVectorBothMiddle( gv1, gv2 ) ) {
        if ( bGVectorSameSide( gv1, gv2 ) ) return(FALSE);

        if ( bGVectorDifferentSides( gv1, gv2 ) ) {
            zX3dClipToPlane( &gv1, &gv1, &gv2, x3dEngine->nViewNear );
            zX3dClipToPlane( &gv2, &gv1, &gv2, x3dEngine->nViewFar );
        } else {
            if ( bGVectorMiddle( gv1 ) ) {
                if ( bGVectorFront(gv2) )   nPlane = x3dEngine->nViewNear;
                else                        nPlane = x3dEngine->nViewFar;
                zX3dClipToPlane( &gv2, &gv1, &gv2, nPlane );
            } else {
                if ( bGVectorFront(gv1) )   nPlane = x3dEngine->nViewNear;
                else                        nPlane = x3dEngine->nViewFar;
                zX3dClipToPlane( &gv1, &gv1, &gv2, nPlane );
            }
        }
    }
    *iPX1 = (int)nGVX(gv1);
    *iPY1 = (int)nGVY(gv1);
    *iPX2 = (int)nGVX(gv2);
    *iPY2 = (int)nGVY(gv2);
    return(TRUE);   
}



