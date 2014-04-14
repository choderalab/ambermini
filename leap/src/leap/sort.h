/*
 *      File: sort.h
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
 *              This is an interface to the general sorting routine
 *              'Sort'.  It allows for sorting by string, integer, and
 *              double keys.
 */

#ifndef	SORT_H
#define	SORT_H

#define	SORT_ASCENDING		TRUE
#define	SORT_DESCENDING		FALSE

extern void	SortByInteger(GENP PStart, int iElements, int iSize, 
			GENP PFirst, BOOL bAscending);
extern void	SortByDouble( GENP PStart, int iElements, int iSize, 
			GENP PFirst, BOOL bAscending);
extern void	SortByString(GENP PStart, int iElements, int iSize,
			GENP PFirst, BOOL bAscending);

typedef BOOL	(*SIFTFUNCTION)(); 

extern void	Sift(GENP PData, int iElementSize, int iElements, 
			SIFTFUNCTION bFCriteria, int *iPFirstFalse);
#endif
