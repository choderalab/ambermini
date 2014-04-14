/*
 *      File:   zMatrix.h
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



#ifndef	ZMATRIX_H
#define	ZMATRIX_H

#include	"vector.h"

extern void	zMatrixInit();

extern void	ZMatrixBondAngleTorsion( VECTOR *vPPos, VECTOR *vPBond, 
			VECTOR *vPAngle, VECTOR *vPTorsion,
			double dBond, double dAngle, double dTorsion );

extern VECTOR	zvZMatrixCalculatePositionFromAngles( double dAngleA,
			double dAngleB, double dAngleC, double dBond );

extern void	ZMatrixBondTwoAnglesOrientation( VECTOR *vPPos, VECTOR *vPAtomC, 
			VECTOR *vPAtomA, VECTOR *vPAtomB, 
			double dBond, double dAngleA, double dAngleB,
			double dOrient );

extern void	ZMatrixBondAngle( VECTOR *vPPos, VECTOR *vPAtom2, VECTOR *vPAtom3,
			double dBond, double dAngle );

extern void	ZMatrixBond( VECTOR *vPPos, VECTOR *vPAtom2, double dBond );

extern void	ZMatrixNothing( VECTOR *vPPos );

#endif
