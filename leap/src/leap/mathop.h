/*
 *      File:   mathop.h
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



#ifndef	MATHOP_H
#define	MATHOP_H

extern void	MathOpConvertNonBondToAC(double dE1, double dR1, 
			double dE2, double dR2, double *dPA, double *dPC);
extern void	MathOpDiagonalize( MATRIX mA, VECTORASPTR vPEigen, MATRIX mV, 
			int *iPNrot );
extern void	MathOpMomentOfInertia( UNIT uUnit, MATRIX mMoment );


#endif
