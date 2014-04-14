/*
 *      File:   zMatrix.c
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
 *              Calculate external coordinates from zMatrix entries.
 */




#include	"basics.h"

#include        "vector.h"
#include        "matrix.h"
#include        "zMatrix.h"


#define         MAXNEWTONSTEPS  20

VECTOR	vXAxis, vYAxis, vZAxis;

void
zMatrixInit()
{
	vXAxis.dX = 1.0;
	vXAxis.dY = 0.0;
	vXAxis.dZ = 0.0;

	vYAxis.dX = 0.0;
	vYAxis.dY = 1.0;
	vYAxis.dZ = 0.0;

	vZAxis.dX = 0.0;
	vZAxis.dY = 0.0;
	vZAxis.dZ = 1.0;
}

/*
 *      ZMatrixBondAngleTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Calculate the vector for an atom that is
 *      specified by a bond length, a bond angle, and a torsion angle
 *      wrt three other vectors.
 */
void
ZMatrixBondAngleTorsion( VECTOR *vPPos, VECTOR *vPBond, VECTOR *vPAngle, VECTOR *vPTorsion,
	double dBond, double dAngle, double dTorsion )
{
MATRIX          mT;
double          dAngleX, dAngleY, dAngleZ;
VECTOR          vTrans, vTemp32, vTemp43, vTempXZ, vNew;

    MESSAGE(( "Bond: %lf   Angle: %lf    Torsion: %lf\n",
                dBond, dAngle, dTorsion ));
                   
                /* The procedure for finding the the coordinate is: */
                /* Translate vAtom2 to the origin -> 3'-2' */
                /* Find angle between PROJ((3'-2'),YZ plane) & Y axis */
                /* Rotate into XZ plane */
                /* Find angle between (3''-2'') and X axis */
                /* Rotate onto X axis */
                /* Find angle between PROJ((4'''-3'''),YZ plane) and  Z axis */
                /* Rotate onto XZ plane */
                /* Calculate coordinates in 3Space */
                /* Apply the reverse transformation to the new point */
                /* Actually, all that is done is the elements for the */
                /* forward transformations are calculated then used */
                /* to generate an inverse transform matrix */
    vTrans = *vPBond;

    	/* translate the two known bond vectors to the origin */
    vTemp32 = vVectorSub( vPAngle, vPBond );	/* 'bond' -> 'angle' atoms */
    vTemp43 = vVectorSub( vPTorsion, vPAngle );	/* 'angle' -> 'torsion' atoms */

    	/* project the 'bond->angle' vector on the XZ plane */
    vTempXZ = vTemp32;
    VectorSetY( &vTempXZ, 0.0 );

	/* measure the angle of the projection to the X axis 
	  around the Y axis in a given ('Abs') direction */
    if ( dVectorLen(&vTempXZ) != 0.0 ) {
        dAngleY = dVectorAbsAngle( &vTempXZ, &vXAxis, &vYAxis );
	/* rotate the origin-translated vectors around Y axis
	   so that the 'bond->angle' one is in the XY plane */
    	MatrixYRotate( mT, -dAngleY );
    	MatrixTimesVector( vTemp32, mT, vTemp32 );
    	MatrixTimesVector( vTemp43, mT, vTemp43 );
    } else 
	dAngleY = 0.0;
 
	/* measure the angle of the rotated 'bond->angle' vector to 
	   the X axis around the Z axis in a given direction */
    dAngleZ = dVectorAbsAngle( &vTemp32, &vXAxis, &vZAxis );

	/* rotate the 'angle->torsion' vector around the Z axis
	   so that the 'bond->angle' one would be on the X axis */
    MatrixZRotate( mT, -dAngleZ );
    MatrixTimesVector( vTemp43, mT, vTemp43 );

	/* project 'angle->torsion' vector on YZ plane and measure its angle
	   to the Z axis around the X axis in a given direction. */
    VectorSetX( &vTemp43, 0.0 );
    dAngleX = dVectorAbsAngle( &vTemp43, &vZAxis, &vXAxis );
 
	/* create the new point as a vector from the origin */
    VectorDef( &vNew, dBond*cos(dAngle),
                      dBond*sin(dAngle)*sin(dTorsion),
                      dBond*sin(dAngle)*cos(dTorsion) );

	/* rotate and translate new point into position in 
	   the original context by applying transformations 
	   in reverse */
    MatrixXRotate( mT, dAngleX );
    MatrixTimesVector( vNew, mT, vNew );
    MatrixZRotate( mT, dAngleZ );
    MatrixTimesVector( vNew, mT, vNew );
    MatrixYRotate( mT, dAngleY );
    MatrixTimesVector( vNew, mT, vNew );
    MatrixTranslate( mT, dVX(&vTrans), dVY(&vTrans), dVZ(&vTrans) );
    MatrixTimesVector( vNew, mT, vNew );
    *vPPos = vNew;
    MESSAGE(( "ZMatrixAll:  %lf,%lf,%lf\n", 
			dVX(vPPos), dVY(vPPos), dVZ(vPPos) ));
}


/*
 *      zvZMatrixCalculatePositionFromAngles
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Use NEWTON-RAPHSON method for finding the coordinate for the
 *      vector(vC) which is dAngleA from the vector vA (on the X axis)
 *	and dAngleB from a vector (vB) which lies in the XY plane 
 *	dAngleC from the X axis.
 *
 *      The point is dBond from the origin.
 *
 */
VECTOR
zvZMatrixCalculatePositionFromAngles( double dAngleA, double dAngleB, 
                                               double dAngleC, double dBond )
{
int             iCount;
double		dCosA, dSinA;
double		dCosB, dSinB;
double		dCosC, dSinC;
double		dCosX, dSinX;
double		dX, dXNew;
double		dF1, dF2;
VECTOR		vNew;

    dCosA = cos(dAngleA);
    dSinA = sin(dAngleA);
    dCosB = cos(dAngleB);
    dSinB = sin(dAngleB);
    dCosC = cos(dAngleC);
    dSinC = sin(dAngleC);

		/* The idea is to minimize the function: */
		/* E = ( DOT(vC,vB) - cos(dAngleB) )^2 */
		/* using NEWTONS method */
		/* The vector vC is constrained to make the angle */
		/* dAngleA with vA */
		/* The vector vC makes the angle (dX) with the XY plane */
		/* and the parameter that is optimized is dX */


		/* A reasonable starting point */

    dX = dAngleB;
    iCount = 0;
    while ( iCount <MAXNEWTONSTEPS ) {

	dCosX = cos(dX);
	dSinX = sin(dX);

	dF1 = -2*dSinA*dSinC*(-dCosB + dCosA*dCosC + 
		    dCosX*dSinA*dSinC)*dSinX;

        dF2 = -2.0*dCosX*dSinA*dSinC*
		   (-dCosB + dCosA*dCosC + dCosX*dSinA*dSinC) + 
		  2.0*pow2(dSinA)*pow2(dSinC)*pow2(dSinX);

        MESSAGE( ( "Iteration %d dF1=%lf  dF2=%lf  dB=%lf\n",
                         iCount, dF1, dF2, dX ) );
        if ( fabs(dF1) < VERYSMALL*10.0 ) break;
        if ( fabs(dF2) < VERYSMALL ) {
DFATAL(( "Could not optimize! dF1 = %lf, dF2 = %lfdX = %lf steps=%d", 
		dF1, dF2, dX, iCount ));
	}
        dXNew = dX - dF1/dF2;
        if ( fabs(dXNew - dX) < VERYSMALL ) break;
        dX = dXNew;
        iCount++;
    }
   
#ifdef  DEBUG 
    if ( iCount > MAXNEWTONSTEPS ) 
        DDEBUG( ("Exceeded maximum number of Newton Raphson steps: %d\n",
                        MAXNEWTONSTEPS) );
#endif

                /* Generate new coordinate */

    VectorDef( &vNew, dBond*cos(dAngleA), 
                      dBond*sin(dAngleA)*cos(dX),
                      dBond*sin(dAngleA)*sin(dX) );
    return(vNew);
}




/*
 *      ZMatrixBondTwoAnglesOrientation
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinate for the atom when
 *      the orientation, a bond length and two angles are supplied.
 *      The orientation is a positive or negative number which specifies
 *      the orientation of the new position.  It is calculated by:
 *              a=crossProduct( vPAtomA-vPCenter, vPAtomB-vPCenter );
 *              orientation = dotProduct( vPPos-vPCenter, a );
 *
 *      vPAtomC points to the position of the central atom.
 *
 */
void
ZMatrixBondTwoAnglesOrientation( VECTOR *vPPos, VECTOR *vPAtomC, VECTOR *vPAtomA, 
	VECTOR *vPAtomB, double dBond, double dAngleA, double dAngleB, double dOrient )
{
MATRIX          mT, mT1, mT2, mTX, mTY, mTZ, mTT;
double          dAngleX, dAngleY, dAngleZ;
double          dAngle;
VECTOR          vTrans, vTempAC, vTempBC, vTempXZ, vNew, vLab;

                /* The procedure for finding the the coordinate is: */
                /* Translate vAtomC to the origin -> A'-C' */
                /* Find angle between PROJ((A'-C'),YZ plane) & Y axis */
                /* Rotate into XZ plane */
                /* Find angle between (A''-C'') and X axis */
                /* Rotate onto X axis */
                /* Find angle between PROJ((B'''-C'''),YZ plane) and Y axis */
                /* Rotate onto XY plane */
                /* Calculate coordinates in 3Space */
                /* Apply the reverse transformation to the new point */
                /* Actually, all that is done is the elements for the */
                /* forward transformations are calculated then used */
                /* to generate an inverse transform matrix */

    vTrans = *vPAtomC;
    vTempAC = vVectorSub( vPAtomA, vPAtomC );
    vTempBC = vVectorSub( vPAtomB, vPAtomC );
MESSAGE(( "AC= %lf, %lf, %lf\n", 
        dVX(&vTempAC), dVY(&vTempAC), dVZ(&vTempAC) ));
MESSAGE(( "BC= %lf, %lf, %lf\n", 
        dVX(&vTempBC), dVY(&vTempBC), dVZ(&vTempBC) ));
    vTempXZ = vTempAC;
    VectorSetY( &vTempXZ, 0.0 );
    if ( dVectorLen(&vTempXZ) != 0.0 ) {
        dAngleY = dVectorAbsAngle( &vTempXZ, &vXAxis, &vYAxis );
    } else dAngleY = 0.0;
    
    MatrixYRotate( mT, -dAngleY );
    MatrixTimesVector( vTempAC, mT, vTempAC );
    MatrixTimesVector( vTempBC, mT, vTempBC );
MESSAGE(( "Rotated around Y\n" ));
MESSAGE(( "New AC= %lf, %lf, %lf\n", 
        dVX(&vTempAC), dVY(&vTempAC), dVZ(&vTempAC) ));
MESSAGE(( "New BC= %lf, %lf, %lf\n", 
        dVX(&vTempBC), dVY(&vTempBC), dVZ(&vTempBC) ));

    dAngleZ = dVectorAbsAngle( &vTempAC, &vXAxis, &vZAxis );
    MatrixZRotate( mT, -dAngleZ );
    MatrixTimesVector( vTempBC, mT, vTempBC );
#ifdef DEBUG
    MatrixTimesVector( vTempAC, mT, vTempAC );
#endif
MESSAGE(( "Rotated around Z\n" ));
MESSAGE(( "New AC= %lf, %lf, %lf\n", 
        dVX(&vTempAC), dVY(&vTempAC), dVZ(&vTempAC) ));
MESSAGE(( "New BC= %lf, %lf, %lf\n", 
        dVX(&vTempBC), dVY(&vTempBC), dVZ(&vTempBC) ));
        
    VectorSetX( &vTempBC, 0.0 );

    dAngleX = dVectorAbsAngle( &vTempBC, &vYAxis, &vXAxis );

                /* Build the transformation matrix to convert from */
                /* lab coordinates to molecule coordinates in mT*/

    MatrixXRotate( mTX, dAngleX );
    MatrixZRotate( mTZ, dAngleZ );
    MatrixYRotate( mTY, dAngleY );
    MatrixTranslate( mTT, dVX(&vTrans), dVY(&vTrans), dVZ(&vTrans) );
    MatrixMultiply( mT1, mTZ, mTX );
    MatrixMultiply( mT2, mTY, mT1 );
    MatrixMultiply( mT, mTT, mT2 );
    
                /* Calculate coordinates of new atom */
    dAngle = dVectorAtomAngle( vPAtomA, vPAtomC, vPAtomB );         
    vLab = zvZMatrixCalculatePositionFromAngles( dAngleA, dAngleB, 
						dAngle, dBond );

    if ( dOrient != 0.0 ) {
        VectorSetZ( &vLab, dOrient*dVZ(&vLab) );
    }

                /* If there is no chirality defined yet then just */
                /* leave it the way it is */
        
    MatrixTimesVector( vNew, mT, vLab );
    *vPPos = vNew;
    MESSAGE(( "ZMatrix2Angle:  %lf,%lf,%lf\n", 
			dVX(vPPos), dVY(vPPos), dVZ(vPPos) ));
}






/*
 *      ZMatrixBondAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinates for an atom using only
 *      a bond length and a bond angle.
 *
 *      Place the bond along the X axis and put the new vector in the
 *      X-Y plane, then rotate it back.
 */
void
ZMatrixBondAngle( VECTOR *vPPos, VECTOR *vPAtom2, VECTOR *vPAtom3, 
	double dBond, double dAngle )
{
MATRIX          mT;
double          dAngleY, dAngleZ;
VECTOR          vTrans, vTempX, vTempXZ, vNew;

                /* The procedure for finding the the coordinate is: */
                /* Translate vAtom2 to the origin -> 3'-2' */
                /* Find angle between PROJ((3'-2'),YZ plane) & Y axis */
                /* Rotate into XZ plane */
                /* Find angle between (3''-2'') and X axis */
                /* Rotate onto X axis */
                /* Calculate coordinates in XY plane */
                /* Apply the reverse transformation to the new point */
                /* Actually, all that is done is the elements for the */
                /* forward transformations are calculated then used */
                /* to generate an inverse transform matrix */

    vTrans = *vPAtom2;
    vTempX = vVectorSub( vPAtom3, vPAtom2 );
    vTempXZ = vTempX;
    VectorSetY( &vTempXZ, 0.0 );
    if ( dVectorLen(&vTempXZ) != 0.0 ) {
        dAngleY = dVectorAbsAngle( &vTempXZ, &vXAxis, &vYAxis );
    } else dAngleY = 0.0;
   
    MESSAGE(( "Angle around Y=%lf\n", dAngleY ));
 
    MatrixYRotate( mT, -dAngleY );
    MatrixTimesVector( vTempX, mT, vTempX );

MESSAGE(( "Rotated around Y = %lf, %lf, %lf\n",
                dVX(&vTempX), dVY(&vTempX), dVZ(&vTempX) ));
    
    dAngleZ = dVectorAbsAngle( &vTempX, &vXAxis, &vZAxis );

#ifdef  DEBUG  
MatrixZRotate( mT, -dAngleZ );
MatrixTimesVector( vTempX, mT, vTempX );
MESSAGE(( "Rotated around Z = %lf, %lf, %lf\n",
                dVX(&vTempX), dVY(&vTempX), dVZ(&vTempX) ));
#endif
 
    MESSAGE(( "Angle around Z=%lf\n", dAngleZ ));
 
    VectorDef( &vNew, dBond*cos(dAngle), dBond*sin(dAngle), 0.0 );
    MatrixZRotate( mT, dAngleZ );
    MatrixTimesVector( vNew, mT, vNew );
    MatrixYRotate( mT, dAngleY );
    MatrixTimesVector( vNew, mT, vNew );
    MatrixTranslate( mT, dVX(&vTrans), dVY(&vTrans), dVZ(&vTrans) );
    MatrixTimesVector( vNew, mT, vNew );
    *vPPos = vNew;
    MESSAGE(( "ZMatrixBondAngle:  %lf,%lf,%lf\n", 
			dVX(vPPos), dVY(vPPos), dVZ(vPPos) ));
}





/*
 *      ZMatrixBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinates for an atom with only a bond
 *      length.
 *
 *      Place the vector at a distance dBond along the x-axis.
 */
void
ZMatrixBond( VECTOR *vPPos, VECTOR *vPAtom2, double dBond )
{
VECTOR          vNew;

    VectorDef( &vNew, dBond, 0.0, 0.0 );
    vNew = vVectorAdd( &vNew, vPAtom2 );
    *vPPos = vNew;
    MESSAGE(( "ZMatrixBond:  %lf,%lf,%lf\n", 
			dVX(vPPos), dVY(vPPos), dVZ(vPPos) ));
}





/*
 *      ZMatrixNothing
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinates for an atom with no
 *      internal coordinates.
 *
 *      Place the vector at the origin.
 */
void
ZMatrixNothing( VECTOR *vPPos )
{
    VectorDef( vPPos, 0.0, 0.0, 0.0 );
    MESSAGE(( "ZMatrixNothing:  %lf,%lf,%lf\n", 
			dVX(vPPos), dVY(vPPos), dVZ(vPPos) ));
}






