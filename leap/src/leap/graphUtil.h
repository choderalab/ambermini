/*
 *	File:	graphUtil.h
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
 *	Description:
 *		This file contains functions that perform
 *		various operations on graphs.  Some of the
 *		functions are finding the shortest paths between
 *		two vertices, finding all of the smallest cycles,
 *		so on and so on.
 */





/*
 *	Functions
 */

extern BOOL	bGraphUtilFindShortestPath( UNIT uUnit, ATOM aStart, ATOM aStop,
                	void (*fPCallback)());
extern void	zGraphUtilDescribeRing( INTERNAL iRing );
extern BOOL	zbGraphUtilSeparate( INTERNAL iBig, int iBigIndex, 
			INTERNAL iSmall, int iSmallIndex );
extern int	ziGraphUtilJoinRingGroups( VARARRAY vaRingGroup, int iGroup1, 				int iGroup2 );
extern void	GraphUtilFindAllSmallestRingsAndRingGroups( UNIT uUnit, 
			VARARRAY *vaPRingGroups);
extern void	GraphUtilDestroyRingGroupVarArray( VARARRAY *vaPRingGroup );
extern void	GraphUtilFindAllSmallestRings( UNIT uUnit );



