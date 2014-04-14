/*
 *      File:   mathop.c
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
 *              This file tries to contain many of the mathematical
 *              operations that are used in molecular modeling.
 *              Eg: converting nonbond parameters for two individual
 *                      atoms into a single A,C pair. etc.
 */





#include	"basics.h"

#include	"vector.h"
#include	"matrix.h"

#include	"classes.h"



/*
 *      MathOpConvertNonBondToAC
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Convert a pair of nonbond parameters dE1, dR1 and dE2, dR2 to
 *      a set of A and C for a non-bond interaction.
 */
void
MathOpConvertNonBondToAC( double dE1, double dR1, double dE2, double dR2, 
		double *dPA, double *dPC )
{
double  dE, dR, dR6, dER6;


    dR = dR1 + dR2;
    dE = sqrt( dE1 * dE2 );
    dR6 = pow6(dR);
    dER6 = dE*dR6;
    *dPC = 2.0*dER6;    /* C = 2*dE*dR^6 */
    *dPA = dER6*dR6;    /* A = dE*dR^12 */
}





/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 *
 *	Calculate principle axis of rotation.
 *
 */

		/* Currently only work with 3x3 matrices and 3vectors */
#define	VSIZE	3
/*
 *	MathOpDiagonalize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Diagonalize a VSIZExVSIZE matrix using the Jacobi method.
 *	Currently works only on the upper 3x3 matrix within a
 *	4x4 matrix.  Diagonalize the (mA)
 *	matrix, returning the Eigenvalues in (vD) and the Eigenvectors
 *	in (mV).  Return the number of rotations performed in (*iPNrot).
 */
void
MathOpDiagonalize( MATRIX mA, VECTORASPTR vPEigen, MATRIX mV, int *iPNrot )
{
VECTORASPTR	vD;
VECTORASARRAY	vB, vZ;
int		j, iq, ip, i;
double		thresh, theta, tau, t, sm, s, h, g, c;
#define	ROTATE(a,i,j,k,l)	g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
				a[k][l]=h+s*(g-h*tau);


    vD = vPEigen;
    for ( ip=0; ip<VSIZE; ip++ ) {
	for ( iq=0; iq<VSIZE; iq++ ) mV[ip][iq] = 0.0;
	mV[ip][ip] = 1.0;
    }
    for ( ip=0; ip<VSIZE; ip++ ) {
	vB[ip] = vD[ip] = mA[ip][ip];
	vZ[ip] = 0.0;
    }
    *iPNrot = 0;
    for ( i=1; i<=50; i++ ) {
	sm = 0.0;
	for ( ip=0; ip< VSIZE-1; ip++ ) {
	    for ( iq=ip+1; iq<VSIZE; iq++ )
		sm += fabs(mA[ip][iq]);
	}
	if ( sm == 0.0 ) return;
	if ( i<4 )	thresh = 0.2*sm/(VSIZE*VSIZE);
	else		thresh = 0.0;
	for ( ip=0; ip<VSIZE-1; ip++ ) {
	    for ( iq=ip+1; iq<VSIZE; iq++ ) {
		g = 100.0*fabs(mA[ip][iq]);
		if ( i>4 && fabs(vD[ip])+g == fabs(vD[ip])
			&& fabs(vD[iq])+g  == fabs(vD[iq]))
		    mA[ip][iq] = 0.0;
		else if ( fabs(mA[ip][iq]) > thresh ) {
		    h = vD[iq]-vD[ip];
		    if ( fabs(h)+g == fabs(h))
			t = (mA[ip][iq])/h;
		    else {
			theta = 0.5*h/(mA[ip][iq]);
			t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
			if ( theta < 0.0 ) t = -t;
		    }
		    c = 1.0/sqrt(1+t*t);
		    s = t*c;
		    tau = s/(1.0+c);
		    h = t*mA[ip][iq];
		    vZ[ip] -= h;
		    vZ[iq] += h;
		    vD[ip] -= h;
		    vD[iq] += h;
		    mA[ip][iq] = 0.0;
		    for ( j=0; j<=ip-1; j++ ) {
			ROTATE( mA, j, ip, j, iq );
		    }
		    for ( j=ip+1; j<=iq-1; j++ ) {
			ROTATE( mA, ip, j, j, iq );
		    }
		    for ( j=iq+1; j<VSIZE; j++ ) {
			ROTATE( mA, ip, j, iq, j );
		    }
		    for ( j=0; j<VSIZE; j++ ) {
			ROTATE( mV, j, ip, j, iq );
		    }
		    ++(*iPNrot);
		}
	    }
	}
	for ( ip=0; ip<VSIZE; ip++ ) {
	    vB[ip] += vZ[ip];
	    vD[ip] = vB[ip];
	    vZ[ip] = 0.0;
	}
    }
    DFATAL(( "Too many iterations" ));
}



/*
 *	MathOpMomentOfInertia
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Calculate the Moment of Inertia, using 1 as the
 *	mass of all atoms.  The moment of inertia
 *	is with respect to the coordinate axis passing
 *	through the origin.  Return the moment in the upper 3x3 
 *	matrix of a 4x4 matrix.
 */
void
MathOpMomentOfInertia( UNIT uUnit, MATRIX mMoment )
{
double		dSxx, dSyy, dSzz, dSxy, dSxz, dSyz;
LOOP		lAtoms;
ATOM		aAtom;
VECTOR		vPos;
double		dX2, dY2, dZ2;

    MatrixIdentity( mMoment );
    dSxx = 0.0;
    dSyy = 0.0;
    dSzz = 0.0;
    dSxy = 0.0;
    dSxz = 0.0;
    dSyz = 0.0;

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
	vPos = vAtomPosition(aAtom);
	dX2 = pow2(dVX(&vPos));
	dY2 = pow2(dVY(&vPos));
	dZ2 = pow2(dVZ(&vPos));
	dSxx += (dY2+dZ2);
	dSyy += (dX2+dZ2);
	dSzz += (dX2+dY2);
	dSxy += dVX(&vPos)*dVY(&vPos);
	dSxz += dVX(&vPos)*dVZ(&vPos);
	dSyz += dVY(&vPos)*dVZ(&vPos);
    }
    mMoment[0][0] =  dSxx; mMoment[0][1] = -dSxy; mMoment[0][2] = -dSxz;
    mMoment[1][0] = -dSxy; mMoment[1][1] =  dSyy; mMoment[1][2] = -dSyz;
    mMoment[2][0] = -dSxz; mMoment[2][1] = -dSyz; mMoment[2][2] =  dSzz;
}

    		







