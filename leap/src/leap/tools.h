/*
 *      File:   tools.h
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
 *              This file contains several routines that
 *              perform different functions.
 *              These functions are currently called from
 *              the command line interpretor.
 */



                /* Different functions that can be performed during */
                /* a distance search */

#define DISTANCE_SEARCH_CREATE_BONDS    1
#define DISTANCE_SEARCH_PRINT_WARNINGS  2

#define DEFAULT_DISTANCE_SEARCH         2.0


                /* NOSHELL is used to specify a ridiculously thick */
                /* solvent shell, so that NO solvent is removed */
                /* because of being too far from the solute */

#define NOSHELL         999999.0

extern BOOL     bToolGeometricCenter( OBJEKT oObjekt, VECTOR *vPCenter );
extern int      iToolDistanceSearch( CONTAINER cCont, double dCloseness, 
                        BOOL bAbsoluteDistance, int iOperation );
extern LIST     lToolListOfResidues( UNIT uUnit, LIST lResidues );
extern void     ToolCenterUnitByRadii( UNIT uUnit, BOOL bOrient );
extern void     ToolInitSolventPotential( UNIT uUnit, VARARRAY vaSolvent,
                        int *iPMinPotRes, int *iPMaxPotRes );
extern void     ToolOctBoxCheck( UNIT uSolute, double *dPBuf, BOOL bMsg );
extern void     ToolOrientPrincipleAxisAlongCoordinateAxis( UNIT uUnit );
extern void     ToolReplaceSolvent( UNIT uUnit, VARARRAY vaSolvent, int iSolv,
                        UNIT uIon, double dCharge, 
                        int *iPMinPotRes, int *iPMaxPotRes );
extern void     ToolSanityCheckBox( UNIT uUnit );
extern UNIT     zToolSetupSolvent( UNIT uSolvent );
extern void     zToolSolvateAndShell( UNIT uSolute, UNIT uSolvent,
                        double dXW, double dYW, double dZW, double dCloseness,
                        double dFarness, BOOL bShell, BOOL bClip, BOOL bOct,
                        BOOL bIsotropic );
extern void     zToolSolvateInSphere( UNIT uSolute, UNIT uSolvent, 
                        VECTOR *vPCenter, double dRadius, double dCloseness );
extern void     ToolSetUnitBoxByCenters( UNIT uUnit );

