/*
 *      File:   vector.c
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
 */



#include        "basics.h"

#include        "vector.h"


/*
 *
 *      Vector objects.
 *
 *      Perform vector arithmatic.
 */



/*
 *      vVectorAdd
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add two vectors.
 */
VECTOR
vVectorAdd( VECTOR *vPX, VECTOR *vPY )
{
VECTOR  vA;

    vA.dX = vPX->dX + vPY->dX;
    vA.dY = vPX->dY + vPY->dY;
    vA.dZ = vPX->dZ + vPY->dZ;

    return(vA);
}




/*
 *      vVectorSub
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Subtract two vectors.
 */
VECTOR
vVectorSub( VECTOR *vPX, VECTOR *vPY )
{
VECTOR  vA;

    vA.dX = vPX->dX - vPY->dX;
    vA.dY = vPX->dY - vPY->dY;
    vA.dZ = vPX->dZ - vPY->dZ;

    return(vA);
}





/*
 *      dVectorLen
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the length of the vector
 */
double
dVectorLen( VECTOR *vPX )
{
double  dLen;

    dLen = sqrt(  (vPX->dX)*(vPX->dX) 
                + (vPX->dY)*(vPX->dY) 
                + (vPX->dZ)*(vPX->dZ) );
    return(dLen);
} 



/*
 *      dVectorAngle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the angle between two vectors in Radians.
 */
double
dVectorAngle( VECTOR *vPX, VECTOR *vPY )
{
VECTOR  vT1, vT2;
double  dLen, dDot;

    VectorDef( &vT1, 0.0, 0.0, 0.0 );
    VectorDef( &vT2, 0.0, 0.0, 0.0 );
    dLen = dVectorLen(vPX);
    if ( dLen != 0.0 ) 
        vT1 = vVectorTimesScalar( vPX, 1/dLen );
    dLen = dVectorLen(vPY);
    if ( dLen != 0.0 ) 
        vT2 = vVectorTimesScalar( vPY, 1/dLen );
    dDot = dVectorDot( &vT1, &vT2 );
    return(myAcos(dDot));
}


/*
 *      dVectorAbsAngle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the absolute angle between two vectors in Radians
 *      the sign of the angle is determined by the direction of
 *      the cross product.
 */
double
dVectorAbsAngle( VECTOR *vPX, VECTOR *vPY, VECTOR *vPRef )
{
VECTOR  vT1, vT2, vT;
double  dLen, dAngle;

    VectorDef( &vT1, 0.0, 0.0, 0.0 );
    VectorDef( &vT2, 0.0, 0.0, 0.0 );
    dLen = dVectorLen(vPX);
    if ( dLen != 0.0 ) vT1 = vVectorTimesScalar( vPX, 1/dLen );
    dLen = dVectorLen(vPY);
    if ( dLen != 0.0 ) vT2 = vVectorTimesScalar( vPY, 1/dLen );

    vT = vVectorCross( &vT1, &vT2 );
    dAngle = myAcos(dVectorDot( &vT1, &vT2 ));
    if ( dVectorDot( &vT, vPRef ) < 0.0 ) dAngle = -dAngle;

    return(dAngle);
}

/*
 *      vVectorTimesScalar
 *
 *      Author: Christian Schafmeister (1991)
 */
VECTOR
vVectorTimesScalar( VECTOR *vPX, double dS )
{
VECTOR  V;

    V.dX = vPX->dX * dS;
    V.dY = vPX->dY * dS;
    V.dZ = vPX->dZ * dS;

    return(V);
}


/*
 *      dVectorDot
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the dot product
 */
double
dVectorDot( VECTOR *vPX, VECTOR *vPY )
{
    return(( vPX->dX*vPY->dX + vPX->dY*vPY->dY + vPX->dZ*vPY->dZ ));
}



/*
 *      vVectorCross
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the cross product of two vectors.
 */
VECTOR
vVectorCross( VECTOR *vPX, VECTOR *vPY )
{
VECTOR  vT;

    vT.dX =  vPX->dY*vPY->dZ - vPX->dZ*vPY->dY;
    vT.dY = -vPX->dX*vPY->dZ + vPX->dZ*vPY->dX;
    vT.dZ =  vPX->dX*vPY->dY - vPX->dY*vPY->dX;
    return(vT);
}



/*
----------------------------------------------------------------------

        The following routines are not PURE vector functions.
        They are used by some routines to perform more complex vector
        calculations.
*/


/*
 *      dVectorAtomChirality
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Calculate the chirality of the vectors A,B,C
 *      This is done by crossing A to B and then dotting the
 *      result with C.
 *      The chirality of the vectors is determined by the sign of
 *      the result, which is determined by whether or not C has
 *      a component in the direction AxB or in the opposite direction.
 */
double  
dVectorAtomChirality( VECTOR *vPCenter, VECTOR *vPA, VECTOR *vPB, VECTOR *vPC )
{
VECTOR          vA, vB, vC, vCross;
double          dDot;

    vA = vVectorSub( vPA, vPCenter );
    vB = vVectorSub( vPB, vPCenter );
    vC = vVectorSub( vPC, vPCenter );
    
    vCross = vVectorCross( &vA, &vB );
    dDot = dVectorDot( &vCross, &vC );
    if ( dDot == 0.0 ) return(0.0);
    if ( dDot < 0.0 ) return(-1.0);
    return(1.0);
}






/*
 *      dVectorAtomNormalizedChirality
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the chirality of the atom.
 *      Either +1.0, -1.0,   or 0.0 if not known or undefined.
 *      Currently the criteria for chirality is absolute orientation
 *      of the vectors joining this atom to its neighbors.
 *      The neighbors are passed as vPA, vPB, vPC, vPD and 
 *      bA, bB, bC, bD define whether or not the position is
 *      defined.
 *
 *      This routine calculates the chirality using the defined
 *      vectors and then flips the sign depending on which
 *      vectors were used to calculate the chirality.
 *      If the ATOMs have (fKnown) set then their coordinate
 *      is considered to be known.
 */
double
dVectorAtomNormalizedChirality( VECTOR *vPCenter, int iCenterCoord,
                                        VECTOR *vPA, BOOL bA,
                                        VECTOR *vPB, BOOL bB,
                                        VECTOR *vPC, BOOL bC,
                                        VECTOR *vPD, BOOL bD )
{
double          dChi;

    dChi = 0.0;
    
                /* Only atoms with 3,4 neighbors have chirality */
                
    if ( (iCenterCoord == 3) ||
         (iCenterCoord == 4) ) {

                        /* If A is not known then use B,C,D to calc chirality */
                        /* The chirality calculated will be negative wrt the */
                        /* correct chirality */
        if ( !bA ) {
            if ( !bB || !bC || !bD ) goto DONE;
            dChi = -dVectorAtomChirality( vPCenter, vPB, vPC, vPD );
            goto DONE;
        }

                        /* If B is not known then use A,C,D to calc chirality */
                        /* The chirality calculated will be correct */

        if ( !bB ) {
            if ( !bA || !bC || !bD ) goto DONE;
            dChi = dVectorAtomChirality( vPCenter, vPA, vPC, vPD );
            goto DONE;
        }

                        /* If C is not known then use A,B,D to calc chirality */
                        /* The chirality calculated will be negative wrt the */
                        /* correct chirality */
        if ( !bC ) {
            if ( !bA || !bB || !bD ) goto DONE;
            dChi = -dVectorAtomChirality( vPCenter, vPA, vPB, vPD );
            goto DONE;
        }

        if ( !bA || !bB || !bC ) goto DONE;
        dChi = dVectorAtomChirality( vPCenter, vPA, vPB, vPC );
    }
    
DONE:   
    return(dChi);
}




/*
 *      dVectorAtomLength
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the distance between the Vectors A,B
 *
 */
double
dVectorAtomLength( VECTOR *vPA, VECTOR *vPB )
{
VECTOR          vX;

    vX = vVectorSub( vPA, vPB );
    return(dVectorLen( &vX ));
}





/*
 *      dVectorAtomAngle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the angle between the Vectors A,B,C
 *
 */
double
dVectorAtomAngle( VECTOR *vPA, VECTOR *vPB, VECTOR *vPC )
{
VECTOR          vX, vY;

    vX = vVectorSub( vPA, vPB );
    vY = vVectorSub( vPC, vPB );
    return(dVectorAngle( &vX, &vY ));
}



/*
 *      dVectorAtomTorsion
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the torsion angle between the Vectors A,B,C,D
 *
 *      Of four atoms     a        d
 *                         \      /
 *                          b----c
 *
 *      The sign of the torsion is positive if looking down b-c
 *      going from front(b-a) to back(c-d) you go clock-wise.
 */
double
dVectorAtomTorsion( VECTOR *vPA, VECTOR *vPB, VECTOR *vPC, VECTOR *vPD )
{
VECTOR          vX, vY, vZ, vN1, vN2, vN;
double          dTorsion, dSign;

    vX = vVectorSub( vPA, vPB );
    vY = vVectorSub( vPC, vPB );
    vZ = vVectorSub( vPD, vPC );
    vN1 = vVectorCross( &vX, &vY );
    vN2 = vVectorCross( &vZ, &vY );
    dTorsion = dVectorAngle( &vN1, &vN2 );
    vN = vVectorCross( &vN1, &vN2 );
    dSign = dVectorDot( &vN, &vY );
    if ( dSign < 0.0 ) 
        dTorsion = -dTorsion;
    return(dTorsion);
}

