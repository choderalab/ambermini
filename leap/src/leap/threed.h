#ifndef	THREED_H
#define	THREED_H

/*
 *      File:   threed.h
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
 *              m[COLUMN][ROW].
 *              For speed these matrices are stored as NUMBR.
 *
 *		This stuff is used in the X-Windows version of LEaP
 */




typedef	float	NUMBR;


#define nabs(x) fabs(x)
  
#define ROWS            4
#define COLUMNS         4

typedef NUMBR GMATRIX[COLUMNS][ROWS];


#define CLIPMIDDLE      0x01
#define CLIPFRONT       0x02
#define CLIPBACK        0x04

typedef	struct  {
	NUMBR		nX;
	NUMBR		nY;
	NUMBR		nZ;
	unsigned char	cClipPlaneStatus;
} GVECTOR;

typedef	struct   {
	GMATRIX   m;
} GMATRIXSTRUCT;



/*
 *      GVECTOR routines
 */

#define GVectorDef( g, x, y, z ) ((g).nX=(x),(g).nY=(y),(g).nZ=(z))
#define cGVectorClipStatus( g )  ((g).cClipPlaneStatus)
#define GVectorSetClipStatus( g, s ) ((g).cClipPlaneStatus = (s) )
#define nGVX(g) (g).nX
#define nGVY(g) (g).nY
#define nGVZ(g) (g).nZ
#define	GVectorSubtract( a, b, c ) ( 	(a).nX = (b).nX-(c).nX,\
					(a).nY = (b).nY-(c).nY,\
					(a).nZ = (b).nZ-(c).nZ )
#define	GVectorAdd( a, b, c ) ( 	(a).nX = (b).nX+(c).nX,\
					(a).nY = (b).nY+(c).nY,\
					(a).nZ = (b).nZ+(c).nZ )
#define	GVectorTimesScalar( a, b, c ) ( (a).nX = (b).nX*(c),\
					(a).nY = (b).nY*(c),\
					(a).nZ = (b).nZ*(c) )
#define	GVectorCross( a, b, c ) (\
(a).nX = (b).nY*(c).nZ-(b).nZ*(c).nY,\
(a).nY = -(b).nX*(c).nZ+(b).nZ*(c).nX,\
(a).nZ = (b).nX*(c).nY-(b).nY*(c).nX)

#define bGVectorBothMiddle( a, b )     \
( ((a).cClipPlaneStatus|(b).cClipPlaneStatus) == CLIPMIDDLE )
#define bGVectorSameSide( a, b )       \
( ((a).cClipPlaneStatus&(b).cClipPlaneStatus) != 0 )
#define bGVectorMiddle( a )     ( (a).cClipPlaneStatus == CLIPMIDDLE )
#define bGVectorDifferentSides( a, b ) \
( ((a).cClipPlaneStatus|(b).cClipPlaneStatus) == (CLIPBACK|CLIPFRONT) )
#define bGVectorBack(a)                ( (a).cClipPlaneStatus == CLIPBACK )
#define bGVectorFront(a)               ( (a).cClipPlaneStatus == CLIPFRONT )

/*
   GMATRIX routines
*/


                           /* Create an identity matrix */
#define  GMatrixIdentity( m ) {\
m[0][0]=1.0;m[1][0]=0.0;m[2][0]=0.0;m[3][0]=0.0;\
m[0][1]=0.0;m[1][1]=1.0;m[2][1]=0.0;m[3][1]=0.0;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]=1.0;m[3][2]=0.0;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]=1.0;}

                           /* Create a diagonal matrix */
#define  GMatrixDiagonal( m, x, y, z, w ) {\
m[0][0]= x ;m[1][0]=0.0;m[2][0]=0.0;m[3][0]=0.0;\
m[0][1]=0.0;m[1][1]= y ;m[2][1]=0.0;m[3][1]=0.0;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]= z ;m[3][2]=0.0;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]= w ;}

            /* Return a rotation matrix around the X axis */
#define  GMatrixXRotate( m, a ) {\
m[0][0]=1.0    ;m[1][0]=0.0    ;m[2][0]=0.0    ;m[3][0]=0.0;\
m[0][1]=0.0    ;m[1][1]= cos(a);m[2][1]= sin(a);m[3][1]=0.0;\
m[0][2]=0.0    ;m[1][2]= -sin(a);m[2][2]= cos(a);m[3][2]=0.0;\
m[0][3]=0.0    ;m[1][3]=0.0    ;m[2][3]=0.0    ;m[3][3]=1.0;}

            /* Return a rotation matrix around the Y axis */
#define  GMatrixYRotate( m, a ) {\
m[0][0]= cos(a);m[1][0]=0.0    ;m[2][0]= -sin(a);m[3][0]=0.0;\
m[0][1]=0.0    ;m[1][1]=1.0    ;m[2][1]=0.0    ;m[3][1]=0.0;\
m[0][2]= sin(a);m[1][2]=0.0    ;m[2][2]= cos(a);m[3][2]=0.0;\
m[0][3]=0.0    ;m[1][3]=0.0    ;m[2][3]=0.0    ;m[3][3]=1.0;}

            /* Return a rotation matrix around the Z axis */
#define  GMatrixZRotate( m, a ) {\
m[0][0]= cos(a);m[1][0]= sin(a);m[2][0]=0.0    ;m[3][0]=0.0;\
m[0][1]= -sin(a);m[1][1]= cos(a);m[2][1]=0.0    ;m[3][1]=0.0;\
m[0][2]=0.0    ;m[1][2]=0.0    ;m[2][2]=1.0    ;m[3][2]=0.0;\
m[0][3]=0.0    ;m[1][3]=0.0    ;m[2][3]=0.0    ;m[3][3]=1.0;}

                     /* Return a scaling matrix */
#define  GMatrixScale( m, s ) {\
m[0][0]= s ;m[1][0]=0.0;m[2][0]=0.0;m[3][0]=0.0;\
m[0][1]=0.0;m[1][1]= s ;m[2][1]=0.0;m[3][1]=0.0;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]= s ;m[3][2]=0.0;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]=1.0;}

                     /* Return a translation matrix */
#define  GMatrixTranslate( m, x, y, z ) {\
m[0][0]=1.0;m[1][0]=0.0;m[2][0]=0.0;m[3][0]= x ;\
m[0][1]=0.0;m[1][1]=1.0;m[2][1]=0.0;m[3][1]= y ;\
m[0][2]=0.0;m[1][2]=0.0;m[2][2]=1.0;m[3][2]= z ;\
m[0][3]=0.0;m[1][3]=0.0;m[2][3]=0.0;m[3][3]=1.0;}

                     /* Multiply a matrix by a scalar */
#define  GMatrixTimesScalar( m, s ) {\
m[0][0]*=s;m[1][0]*=s;m[2][0]*=s;m[3][0]*=s;\
m[0][1]*=s;m[1][1]*=s;m[2][1]*=s;m[3][1]*=s;\
m[0][2]*=s;m[1][2]*=s;m[2][2]*=s;m[3][2]*=s;\
m[0][3]*=s;m[1][3]*=s;m[2][3]*=s;m[3][3]*=s;}

                     /* add two matrices */
#define  GMatrixAdd( m, b ) {\
m[0][0]+=b[0][0];m[1][0]+=b[1][0];m[2][0]+=b[2][0];m[3][0]+=b[3][0];\
m[0][1]+=b[0][1];m[1][1]+=b[1][1];m[2][1]+=b[2][1];m[3][1]+=b[3][1];\
m[0][2]+=b[0][2];m[1][2]+=b[1][2];m[2][2]+=b[2][2];m[3][2]+=b[3][2];\
m[0][3]+=b[0][3];m[1][3]+=b[1][3];m[2][3]+=b[2][3];m[3][3]+=b[3][3];}

                     /* Subtract two matrices */
#define  GMatrixSubtract( m, b ) {\
m[0][0]-=b[0][0];m[1][0]-=b[1][0];m[2][0]-=b[2][0];m[3][0]-=b[3][0];\
m[0][1]-=b[0][1];m[1][1]-=b[1][1];m[2][1]-=b[2][1];m[3][1]-=b[3][1];\
m[0][2]-=b[0][2];m[1][2]-=b[1][2];m[2][2]-=b[2][2];m[3][2]-=b[3][2];\
m[0][3]-=b[0][3];m[1][3]-=b[1][3];m[2][3]-=b[2][3];m[3][3]-=b[3][3];}

                     /* Create a copy of the matrix */
#define  GMatrixCopy( m, b ) {\
m[0][0]=b[0][0];m[1][0]=b[1][0];m[2][0]=b[2][0];m[3][0]=b[3][0];\
m[0][1]=b[0][1];m[1][1]=b[1][1];m[2][1]=b[2][1];m[3][1]=b[3][1];\
m[0][2]=b[0][2];m[1][2]=b[1][2];m[2][2]=b[2][2];m[3][2]=b[3][2];\
m[0][3]=b[0][3];m[1][3]=b[1][3];m[2][3]=b[2][3];m[3][3]=b[3][3];}

#define  GMatrixFromMatrix( m, b ) {\
m[0][0]=b[0][0];m[1][0]=b[1][0];m[2][0]=b[2][0];m[3][0]=b[3][0];\
m[0][1]=b[0][1];m[1][1]=b[1][1];m[2][1]=b[2][1];m[3][1]=b[3][1];\
m[0][2]=b[0][2];m[1][2]=b[1][2];m[2][2]=b[2][2];m[3][2]=b[3][2];\
m[0][3]=b[0][3];m[1][3]=b[1][3];m[2][3]=b[2][3];m[3][3]=b[3][3];}


                     /* Transpose the matrix */
#define  GMatrixTranspose( m, b ) {\
m[0][0]=b[0][0]; m[0][1]=b[1][0]; m[0][2]=b[2][0]; m[0][3]=b[3][0];\
m[1][0]=b[0][1]; m[1][1]=b[1][1]; m[1][2]=b[2][1]; m[1][3]=b[3][1];\
m[2][0]=b[0][2]; m[2][1]=b[1][2]; m[2][2]=b[2][2]; m[2][3]=b[3][2];\
m[3][0]=b[0][3]; m[3][1]=b[1][3]; m[3][2]=b[2][3]; m[3][3]=b[3][3];}

                     /* Invert the upper 3x3 matrix */
#define GMatrixInvert33( m, a ) {\
NUMBR   fDet;\
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
#define  GMatrixMultiply( r, m, b ) {\
int      XiR, XiC, XiI;\
   for (XiR=0;XiR<ROWS;XiR++)\
      for (XiC=0;XiC<COLUMNS;XiC++) {\
         r[XiC][XiR] =0.0;\
         for (XiI=0;XiI<COLUMNS;XiI++)\
            r[XiC][XiR]+=m[XiI][XiR]*b[XiC][XiI]; }}


                     /* Print a graphics matrix */
#define  GMatrixPrint( m ) {\
MESSAGE(( "%10.6f, %10.6f, %10.6f,%10.6f\n",m[0][0],m[1][0],m[2][0],m[3][0]));\
MESSAGE(( "%10.6f, %10.6f, %10.6f, %10.6f\n",m[0][1],m[1][1],m[2][1],m[3][1]));\
MESSAGE(( "%10.6f, %10.6f, %10.6f, %10.6f\n",m[0][2],m[1][2],m[2][2],m[3][2]));\
MESSAGE(( "%10.6f, %10.6f, %10.6f, %10.6f\n",m[0][3],m[1][3],m[2][3],m[3][3]));}



#endif
