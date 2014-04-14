#ifndef	X3D_H

#define	X3D_H
/*
 *      File: x3d.h
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



#include "threed.h"

/*
 *      Default Viewspace parameters
 */

#define DEFAULTFARCLIPPLANE     -210.0
#define DEFAULTNEARCLIPPLANE    -50.0
#define DEFAULTVIEWWIDTH        10.0
#define DEFAULTVIEWHEIGHT       10.0
#define	TRANSLATEMULTIPLIER	1.0
#define	CIRCLERADIUS		0.8

#define DEFAULTSCREENDISTANCE   9.0     /* This effects the perspective */




typedef	struct	{
	GVECTOR		gvCenterOfRotation;
	GVECTOR		gvCenterOfView;
	NUMBR           nScale;
	GMATRIX         gmRotate;
	NUMBR           nViewFar;
	NUMBR           nViewNear;
	NUMBR           nViewHeight;
	NUMBR           nViewWidth;
	NUMBR           nViewCenter;
	NUMBR           nViewScreenDistance;
	NUMBR           nScreenXOffset;
	NUMBR           nScreenYOffset;
	BOOL            bPerspectiveOn;
	int             iScreenHeight;
	int             iScreenWidth;
	int             iUseHeight;
	int             iUseWidth;
	BOOL		bRecalculateForward;
	BOOL		bRecalculateBackward;
	GMATRIX         gmWorldToView;
	GMATRIX		gmWorldToPartialScreen;
	GMATRIX		gmPartialScreenToView;
	GMATRIX		gmViewToWorld;
	GMATRIX		gmPartialScreenToWorld;
} X3DENGINEt;

typedef	X3DENGINEt	*X3DENGINE;

extern void	GMatrixTimesVector( GVECTOR *result, GMATRIX m, GVECTOR *b );
extern void	zX3dClipToPlane( GVECTOR *gvPClip, GVECTOR *gvPA, 
			GVECTOR *gvPB, NUMBR nClip );

extern X3DENGINE	x3dX3dCreate();
extern void		X3dInit(X3DENGINE x3dEngine);
extern void		X3dDestroy(X3DENGINE *x3dPEngine);


extern void	X3dBuildTransform(X3DENGINE x3dEngine);
extern void	X3dBuildInverseTransform(X3DENGINE x3dEngine);
extern void	X3dWorldToScreen(X3DENGINE x3dEngine, GVECTOR *gvPPoint, 
			GVECTOR *gvPNew);
extern void	X3dScreenToWorld(X3DENGINE x3dEngine, GVECTOR *gvPPoint, 
			GVECTOR *gvPNew);
extern void	X3dScreenToView(X3DENGINE x3dEngine, GVECTOR *gvPPoint, 
			GVECTOR *gvPNew);
extern void	X3dViewToWorld(X3DENGINE x3dEngine, GVECTOR *gvPView, 
			GVECTOR *gvPWorld);
extern void	X3dBuildCenteredScaledTransform(X3DENGINE x3dEngine, 
			UNIT uUnit);
extern BOOL	bX3dVisibleLineClip(X3DENGINE x3dEngine, int *iPX1, int *iPY1,
			int *iPX2, int *iPY2, GVECTOR *gvP1, GVECTOR *gvP2);


#define	zX3dSetRecalculateBoth(x3d)	\
	(x3d->bRecalculateForward=TRUE,x3d->bRecalculateBackward=TRUE)


#define	gvX3dCenterOfRotation(x3d)	(x3d->gvCenterOfRotation)
#define	X3dSetCenterOfRotation(x3d,d)	\
		(x3d->gvCenterOfRotation=d,zX3dSetRecalculateBoth(x3d))

#define	gvX3dCenterOfView(x3d)		(x3d->gvCenterOfView)
#define	X3dSetCenterOfView(x3d,d)	\
		(x3d->gvCenterOfView=d,zX3dSetRecalculateBoth(x3d))

#define	nX3dScale(x3d)			(x3d->nScale)
#define	X3dSetScale(x3d,d)		\
		(x3d->nScale=d,zX3dSetRecalculateBoth(x3d))

#define	nX3dFarViewPlane(x3d)		(x3d->nViewFar)
#define	X3dSetFarViewPlane(x3d,d)	\
		(x3d->nViewFar=d,zX3dSetRecalculateBoth(x3d))

#define	nX3dNearViewPlane(x3d)		(x3d->nViewNear)
#define	X3dSetNearViewPlane(x3d,d)	\
		(x3d->nViewNear=d,zX3dSetRecalculateBoth(x3d))

#define	X3dGetWorldRotate(x3d,d)	{GMatrixCopy(d,(x3d)->gmRotate);}
#define	X3dSetWorldRotate(x3d,d)	{GMatrixCopy(x3d->gmRotate,d); \
					 zX3dSetRecalculateBoth(x3d);};

#define	nX3dViewHeight(x3d)		(x3d->nViewHeight)
#define	X3dSetViewHeight(x3d,d)		\
		(x3d->nViewHeight=d,zX3dSetRecalculateBoth(x3d))

#define	nX3dViewWidth(x3d)		(x3d->nViewWidth)
#define	X3dSetViewWidth(x3d,d)		\
		(x3d->nViewWidth=d,zX3dSetRecalculateBoth(x3d))

#define	nX3dViewCenter(x3d)		(x3d->nViewCenter)
#define	X3dSetViewCenter(x3d,d)		\
		(x3d->nViewCenter=d,zX3dSetRecalculateBoth(x3d))

#define	nX3dViewScreenDistance(x3d)	(x3d->nViewScreenDistance)
#define	X3dSetViewScreenDistance(x3d,d)	\
		(x3d->nViewScreenDistance=d,zX3dSetRecalculateBoth(x3d))

#define	nX3dScreenXOffset(x3d)		(x3d->nScreenXOffset)
#define	X3dSetScreenXOffset(x3d,d)	\
		(x3d->nScreenXOffset=d,zX3dSetRecalculateBoth(x3d))

#define	nX3dScreenYOffset(x3d)		(x3d->nScreenYOffset)
#define	X3dSetScreenYOffset(x3d,d)	\
		(x3d->nScreenYOffset=d,zX3dSetRecalculateBoth(x3d))

#define	bX3dPerspectiveOn(x3d)		(x3d->bPerspectiveOn)
#define	X3dSetPerspectiveOn(x3d,d)	(x3d->bPerspectiveOn=d)

#define	iX3dScreenHeight(x3d)		(x3d->iScreenHeight)
#define	X3dSetScreenHeight(x3d,d)	(x3d->iScreenHeight=d)

#define	iX3dScreenWidth(x3d)		(x3d->iScreenWidth)
#define	X3dSetScreenWidth(x3d,d)	(x3d->iScreenWidth=d)

#define	iX3dScreenUseHeight(x3d)	(x3d->iUseHeight)
#define	X3dSetScreenUseHeight(x3d,d)	(x3d->iUseHeight=d)

#define	iX3dScreenUseWidth(x3d)		(x3d->iUseWidth)
#define	X3dSetScreenUseWidth(x3d,d)	(x3d->iUseWidth=d)


#endif	/* X3D_H */
