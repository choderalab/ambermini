/*
 *      File:   matrix.h
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
 *              Macros for handling 4x4 Homogenous coordinate system
 *              matrices.   The MATRIX elements can be accessed
 *              m[COLUMN][ROW]
 */

#ifndef	MATRIX_H
#define	MATRIX_H

#include    <math.h>

#define  ROWS     4
#define  COLUMNS  4

typedef  double MATRIX[COLUMNS][ROWS];

typedef  struct   {
	MATRIX   m;
} MATRIXSTRUCT;

/*
 *	MATRIX routines
 */


extern void	MatrixRotateAround(MATRIX mTransform, VECTOR *vPA, VECTOR *vPB, 
			double dRotate);
extern void	MatrixReflectAcross(MATRIX mTransform, VECTOR *vPA, VECTOR *vPNormal);

                           /* Create an empty matrix */
#define  MatrixZero( m ) {\
m[0][0]=0.0;m[1][0]=0.0;m[2][0]=0.0;m[3][0]=0.0;\
m[0][1]=0.0;m[1][1]=0.0;m[2][1]=0.0;m[3][1]=0.0;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]=0.0;m[3][2]=0.0;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]=0.0;}

                           /* Create an identity matrix */
#define  MatrixIdentity( m ) {\
m[0][0]=1.0;m[1][0]=0.0;m[2][0]=0.0;m[3][0]=0.0;\
m[0][1]=0.0;m[1][1]=1.0;m[2][1]=0.0;m[3][1]=0.0;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]=1.0;m[3][2]=0.0;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]=1.0;}

                           /* Create an identity matrix */
#define  MatrixDef( m,a00,a10,a20,a30,a01,a11,a21,a31,a02,a12,a22,a32,a03,a13,a23,a33 ) {\
m[0][0]=(a00);m[1][0]=(a10);m[2][0]=(a20);m[3][0]=(a30);\
m[0][1]=(a01);m[1][1]=(a11);m[2][1]=(a21);m[3][1]=(a31);\
m[0][2]=(a02);m[1][2]=(a12);m[2][2]=(a22);m[3][2]=(a32);\
m[0][3]=(a03);m[1][3]=(a13);m[2][3]=(a23);m[3][3]=(a33);}

			/* Create a diagonal matrix */
#define  MatrixDiagonal( m, x, y, z ) {\
m[0][0]= x ;m[1][0]=0.0;m[2][0]=0.0;m[3][0]=0.0;\
m[0][1]=0.0;m[1][1]= y ;m[2][1]=0.0;m[3][1]=0.0;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]= z ;m[3][2]=0.0;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]=1.0;}

            /* Return a rotation matrix around the X axis */
#define  MatrixXRotate( m, a ) {\
m[0][0]=1.0    ;m[1][0]=0.0    ;m[2][0]=0.0    ;m[3][0]=0.0;\
m[0][1]=0.0    ;m[1][1]= cos(a);m[2][1]= sin(a);m[3][1]=0.0;\
m[0][2]=0.0    ;m[1][2]= -sin(a);m[2][2]= cos(a);m[3][2]=0.0;\
m[0][3]=0.0    ;m[1][3]=0.0    ;m[2][3]=0.0    ;m[3][3]=1.0;}

            /* Return a rotation matrix around the Y axis */
#define  MatrixYRotate( m, a ) {\
m[0][0]= cos(a);m[1][0]=0.0    ;m[2][0]= -sin(a);m[3][0]=0.0;\
m[0][1]=0.0    ;m[1][1]=1.0    ;m[2][1]=0.0    ;m[3][1]=0.0;\
m[0][2]= sin(a);m[1][2]=0.0    ;m[2][2]= cos(a);m[3][2]=0.0;\
m[0][3]=0.0    ;m[1][3]=0.0    ;m[2][3]=0.0    ;m[3][3]=1.0;}

            /* Return a rotation matrix around the Z axis */
#define  MatrixZRotate( m, a ) {\
m[0][0]= cos(a);m[1][0]= sin(a);m[2][0]=0.0    ;m[3][0]=0.0;\
m[0][1]= -sin(a);m[1][1]= cos(a);m[2][1]=0.0    ;m[3][1]=0.0;\
m[0][2]=0.0    ;m[1][2]=0.0    ;m[2][2]=1.0    ;m[3][2]=0.0;\
m[0][3]=0.0    ;m[1][3]=0.0    ;m[2][3]=0.0    ;m[3][3]=1.0;}

                     /* Return a scaling matrix */
#define  MatrixScale( m, s ) {\
m[0][0]= s ;m[1][0]=0.0;m[2][0]=0.0;m[3][0]=0.0;\
m[0][1]=0.0;m[1][1]= s ;m[2][1]=0.0;m[3][1]=0.0;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]= s ;m[3][2]=0.0;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]=1.0;}

                     /* Return a translation matrix */
#define  MatrixTranslate( m, x, y, z ) {\
m[0][0]=1.0;m[1][0]=0.0;m[2][0]=0.0;m[3][0]= x ;\
m[0][1]=0.0;m[1][1]=1.0;m[2][1]=0.0;m[3][1]= y ;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]=1.0;m[3][2]= z ;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]=1.0;}

                     /* Multiply a matrix by a scalar */
#define  MatrixTimesScalar( m, s ) {\
m[0][0]*=s;m[1][0]*=s;m[2][0]*=s;m[3][0]*=s;\
m[0][1]*=s;m[1][1]*=s;m[2][1]*=s;m[3][1]*=s;\
m[0][2]*=s;m[1][2]*=s;m[2][2]*=s;m[3][2]*=s;\
m[0][3]*=s;m[1][3]*=s;m[2][3]*=s;m[3][3]*=s;}

                     /* add two matrices */
#define  MatrixAdd( m, b ) {\
m[0][0]+=b[0][0];m[1][0]+=b[1][0];m[2][0]+=b[2][0];m[3][0]+=b[3][0];\
m[0][1]+=b[0][1];m[1][1]+=b[1][1];m[2][1]+=b[2][1];m[3][1]+=b[3][1];\
m[0][2]+=b[0][2];m[1][2]+=b[1][2];m[2][2]+=b[2][2];m[3][2]+=b[3][2];\
m[0][3]+=b[0][3];m[1][3]+=b[1][3];m[2][3]+=b[2][3];m[3][3]+=b[3][3];}

                     /* Subtract two matrices */
#define  MatrixSubtract( m, b ) {\
m[0][0]-=b[0][0];m[1][0]-=b[1][0];m[2][0]-=b[2][0];m[3][0]-=b[3][0];\
m[0][1]-=b[0][1];m[1][1]-=b[1][1];m[2][1]-=b[2][1];m[3][1]-=b[3][1];\
m[0][2]-=b[0][2];m[1][2]-=b[1][2];m[2][2]-=b[2][2];m[3][2]-=b[3][2];\
m[0][3]-=b[0][3];m[1][3]-=b[1][3];m[2][3]-=b[2][3];m[3][3]-=b[3][3];}

                     /* Create a copy of the matrix */
#define  MatrixCopy( m, b ) {\
m[0][0]=b[0][0];m[1][0]=b[1][0];m[2][0]=b[2][0];m[3][0]=b[3][0];\
m[0][1]=b[0][1];m[1][1]=b[1][1];m[2][1]=b[2][1];m[3][1]=b[3][1];\
m[0][2]=b[0][2];m[1][2]=b[1][2];m[2][2]=b[2][2];m[3][2]=b[3][2];\
m[0][3]=b[0][3];m[1][3]=b[1][3];m[2][3]=b[2][3];m[3][3]=b[3][3];}


                     /* Transpose the matrix */
#define  MatrixTranspose( m, b ) {\
m[0][0]=b[0][0]; m[0][1]=b[1][0]; m[0][2]=b[2][0]; m[0][3]=b[3][0];\
m[1][0]=b[0][1]; m[1][1]=b[1][1]; m[1][2]=b[2][1]; m[1][3]=b[3][1];\
m[2][0]=b[0][2]; m[2][1]=b[1][2]; m[2][2]=b[2][2]; m[2][3]=b[3][2];\
m[3][0]=b[0][3]; m[3][1]=b[1][3]; m[3][2]=b[2][3]; m[3][3]=b[3][3];}

                     /* Invert the upper 3x3 matrix */
#define MatrixInvert33( m, a ) {\
double   fDet;\
    fDet = ( -a[0][2]*a[1][1]*a[2][0] + a[0][1]*a[1][2]*a[2][0] \
           + a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1]\
           - a[0][1]*a[1][0]*a[2][2] + a[0][0]*a[1][1]*a[2][2]);\
    m[0][0] =( a[1][1]*a[2][2] - a[2][1]*a[1][2])/fDet;\
    m[1][0] =( a[1][2]*a[2][0] - a[1][0]*a[2][2])/fDet;\
    m[2][0] =(-a[1][1]*a[2][0] + a[1][0]*a[2][1])/fDet;\
    m[0][1] =( a[0][2]*a[2][1] - a[0][1]*a[2][2])/fDet; \
    m[1][1] =(-a[0][2]*a[2][0] + a[0][0]*a[2][2])/fDet;\
    m[2][1] =( a[0][1]*a[2][0] - a[0][0]*a[2][1])/fDet;\
    m[0][2] =(-a[0][2]*a[1][1] + a[0][1]*a[1][2])/fDet;\
    m[1][2] =( a[0][2]*a[1][0] - a[0][0]*a[1][2])/fDet;\
    m[2][2] =(-a[0][1]*a[1][0] + a[0][0]*a[1][1])/fDet;}


                     /* multiply two matrices, and place the result in r */
                     /* r must be different from m and b */
#define  MatrixMultiply( r, m, b ) {\
int      XiR, XiC, XiI;\
   for (XiR=0;XiR<ROWS;XiR++)\
      for (XiC=0;XiC<COLUMNS;XiC++) {\
         r[XiC][XiR] =0.0;\
         for (XiI=0;XiI<COLUMNS;XiI++)\
            r[XiC][XiR]+=m[XiI][XiR]*b[XiC][XiI]; }}


#define MatrixTimesVector( v, m, b ) {\
VECTOR t;\
t.dX = m[0][0]*(b).dX + m[1][0]*(b).dY + m[2][0]*(b).dZ + m[3][0];\
t.dY = m[0][1]*(b).dX + m[1][1]*(b).dY + m[2][1]*(b).dZ + m[3][1];\
t.dZ = m[0][2]*(b).dX + m[1][2]*(b).dY + m[2][2]*(b).dZ + m[3][2]; v = t;}



                     /* Print a graphics matrix */
#define  MatrixPrint( m ) {\
MESSAGE(( "%10.6f, %10.6f, %10.6f, %10.6f\n",m[0][0],m[1][0],m[2][0],m[3][0]));\
MESSAGE(( "%10.6f, %10.6f, %10.6f, %10.6f\n",m[0][1],m[1][1],m[2][1],m[3][1]));\
MESSAGE(( "%10.6f, %10.6f, %10.6f, %10.6f\n",m[0][2],m[1][2],m[2][2],m[3][2]));\
MESSAGE(( "%10.6f, %10.6f, %10.6f, %10.6f\n",m[0][3],m[1][3],m[2][3],m[3][3]));}

                     /* Print a FLOAT matrix */
#define  MATRIXPrint( m ) {\
MESSAGE(( "%10.6lf, %10.6lf, %10.6lf, %10.6lf\n",m[0][0],m[1][0],m[2][0],m[3][0]));\
MESSAGE(( "%10.6lf, %10.6lf, %10.6lf, %10.6lf\n",m[0][1],m[1][1],m[2][1],m[3][1]));\
MESSAGE(( "%10.6lf, %10.6lf, %10.6lf, %10.6lf\n",m[0][2],m[1][2],m[2][2],m[3][2]));\
MESSAGE(( "%10.6lf, %10.6lf, %10.6lf, %10.6lf\n",m[0][3],m[1][3],m[2][3],m[3][3]));}



#endif
