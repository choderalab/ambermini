/*
 *      File:   matrix.c
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
 *              Routines for handling 4x4 Homogenous coordinate system
 *              matrices.   The MATRIX elements can be accessed
 *              m[COLUMN][ROW]
 *
 *		Most of the routines for handling the matrices are 
 *		macros in the matrix.h file.
 */





#include	"basics.h"

#include	"vector.h"
#include	"matrix.h"




/*
 *	MatrixRotateAround
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build a rotation matrix around an arbitrary axis in 3-space.
 *	The axis is along the line connecting the points *vPA and *vPB
 *	and the angle to rotate by is in dRotate.
 */
void
MatrixRotateAround( MATRIX mTransform, VECTOR *vPA, VECTOR *vPB, double dRotate )
{
MATRIX		mTranslate, mIntoXZ, mIntoX;
MATRIX		mRotate;
MATRIX		mTranslateInv, mIntoXZInv, mIntoXInv;
VECTOR		vTemp;
VECTOR		vXAxis, vYAxis, vZAxis;
VECTOR		vDir;
double		dAngle;
MATRIX		mT, mU;

    vDir = vVectorSub( vPB, vPA );
    VectorDef( &vXAxis, 1.0, 0.0, 0.0 );
    VectorDef( &vYAxis, 0.0, 1.0, 0.0 );
    VectorDef( &vZAxis, 0.0, 0.0, 1.0 );

	/* Build the MATRIX to place the axis on the origin */

    MatrixTranslate( mTranslate, -dVX(vPA), -dVY(vPA), -dVZ(vPA) );
    MatrixTranslate( mTranslateInv, dVX(vPA), dVY(vPA), dVZ(vPA) );

	/* Build the MATRIX to rotate the axis onto the XZ plane */

    VectorDef( &vTemp, dVX(&vDir), dVY(&vDir), 0.0 );
    dAngle = dVectorAngle( &vTemp, &vXAxis );
    if ( dVY(&vDir) < 0.0 ) dAngle = -dAngle;
    MatrixZRotate( mIntoXZ, dAngle );
    MatrixZRotate( mIntoXZInv, -dAngle );

	/* Build the MATRIX to rotate onto the X axis */

    MatrixTimesVector( vTemp, mIntoXZ, vDir );
    dAngle = dVectorAngle( &vTemp, &vXAxis );
    if ( dVZ(&vTemp) > 0.0 ) dAngle = -dAngle;
    MatrixYRotate( mIntoX, dAngle );
    MatrixYRotate( mIntoXInv, -dAngle );

	/* Build the rotation matrix */

    MatrixXRotate( mRotate, dRotate );

	/* Build the transformation matrix */

    MatrixMultiply( mT, mIntoXZ, mTranslate );
    MatrixMultiply( mU, mIntoX, mT );
    MatrixMultiply( mT, mRotate, mU );
    MatrixMultiply( mU, mIntoXInv, mT );
    MatrixMultiply( mT, mIntoXZInv, mU );
    MatrixMultiply( mTransform, mTranslateInv, mT );

}



/*
 *	MatrixReflectAcross
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build a matrix which reflects points across a plane
 *	The plain lies on vPA and is normal to vPNormal.
 */
void
MatrixReflectAcross( MATRIX mTransform, VECTOR *vPA, VECTOR *vPNormal )
{
MATRIX		mTranslate, mIntoXZ, mIntoX;
MATRIX		mReflect;
MATRIX		mTranslateInv, mIntoXZInv, mIntoXInv;
VECTOR		vTemp;
VECTOR		vXAxis, vYAxis, vZAxis;
double		dAngle;
MATRIX		mT, mU;


    VectorDef( &vXAxis, 1.0, 0.0, 0.0 );
    VectorDef( &vYAxis, 0.0, 1.0, 0.0 );
    VectorDef( &vZAxis, 0.0, 0.0, 1.0 );

	/* Build the MATRIX to place the plane on the origin */

    MatrixTranslate( mTranslate, -dVX(vPA), -dVY(vPA), -dVZ(vPA) );
    MatrixTranslate( mTranslateInv, dVX(vPA), dVY(vPA), dVZ(vPA) );

	/* Build the MATRIX to rotate the normal onto the XZ plane */

    VectorDef( &vTemp, dVX(vPNormal), dVY(vPNormal), 0.0 );
    dAngle = dVectorAngle( &vTemp, &vXAxis );
    if ( dVY(vPNormal) < 0.0 ) dAngle = -dAngle;
    MatrixZRotate( mIntoXZ, dAngle );
    MatrixZRotate( mIntoXZInv, -dAngle );

	/* Build the MATRIX to rotate onto the X axis */

    MatrixTimesVector( vTemp, mIntoXZ, *vPNormal );
    dAngle = dVectorAngle( &vTemp, &vXAxis );
    if ( dVZ(&vTemp) > 0.0 ) dAngle = -dAngle;
    MatrixYRotate( mIntoX, dAngle );
    MatrixYRotate( mIntoXInv, -dAngle );

	/* Build the reflection matrix */

    MatrixDiagonal( mReflect, -1.0, 1.0, 1.0 );

	/* Build the transformation matrix */

    MatrixMultiply( mT, mIntoXZ, mTranslate );
    MatrixMultiply( mU, mIntoX, mT );
    MatrixMultiply( mT, mReflect, mU );
    MatrixMultiply( mU, mIntoXInv, mT );
    MatrixMultiply( mT, mIntoXZInv, mU );
    MatrixMultiply( mTransform, mTranslateInv, mT );

}

