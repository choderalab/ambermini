/*
 *      File:   build.h
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
 */
 


#ifndef	BUILD_H
#define	BUILD_H

#include "classes.h"
#include "minimizer.h"
#include "loop.h"
/*  build.c  */

extern void		BuildInternalsUsingFlags( LOOP *lPAtoms, FLAGS fForSet, 
				FLAGS fForReset, FLAGS fSetFlags, 
				FLAGS fResetFlags );
extern void		BuildExternalsUsingFlags( LOOP *lPLoop, FLAGS fForSet, 
				FLAGS fForReset, FLAGS fSetFlags, 
				FLAGS fResetFlags,
				int *iPAddH, int *iPAddHeavy, int *iPAddUnk,
				BOOL bMsg );
extern void		BuildInternalsForContainer( CONTAINER cCont, 
				FLAGS fSet, FLAGS fReset );
extern void		BuildDestroyInternals( LOOP *lPLoop );
extern void		BuildFixInternals( UNIT uUnit );
extern void		BuildHierarchy( UNIT uUnit );
extern BOOL		bBuildChangeInternalBond( CONTAINER cCont, 
				char *sAtom1, char *sAtom2, double dValue );
extern BOOL		bBuildChangeInternalAngle( CONTAINER cCont, 
				char *sAtom1, char *sAtom2, char *sAtom3, 
				double dValue );
extern BOOL		bBuildChangeInternalTorsion( CONTAINER cCont, 
				char *sAtom1, char *sAtom2,
                                char *sAtom3, char *sAtom4, double dValue );
extern void		BuildRelaxInFramework( UNIT uUnit, MINIMIZER mStrain );
extern void		BuildRotateAroundBondFromTo( CONTAINER cCont, 
				ATOM aInv, ATOM aStart,
				double dRotate, BOOL bInRing );
extern void		bBuildFlipChiralityFor( CONTAINER cContainer, 
				ATOM aFlip );
extern void		BuildInternalsBetweenUnitsUsingFlags( UNIT uFirst, 
				UNIT uSecond, FLAGS fForSet, FLAGS fForReset );
extern void		BuildInternalsForSimpleRings( CONTAINER cContainer );



#endif

