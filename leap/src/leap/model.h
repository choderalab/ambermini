/*
 *      File:   model.h
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
 *              The MODEL object is an object with a single instance
 *              that interacts with the BUILD object and generates
 *              model bond lengths, angles, torsions etc which
 *              are determined from atom types.
 */


#ifndef	MODEL_H
#define	MODEL_H

extern double	dModelBondLength(ATOM aAtom1, ATOM aAtom2);
extern double	dModelBondAngle(ATOM aAtom1, ATOM aAtom2, ATOM aAtom3 );
extern void	ModelAddHydrogens( UNIT uUnit );
extern void	ModelProperTerms(ATOM a1, ATOM a2, ATOM a3, ATOM a4, TORSION tTorsion);
extern void	ModelAngleParm(ATOM a1, ATOM a2, ATOM a3, double *dPK, double *dPE);
extern void	ModelBondParm(ATOM a1, ATOM a2, double *dPK, double *dPE);
extern void	ModelAssignTorsionsAround(ATOM aX, ATOM aY, ATOM aOnly);

#endif
