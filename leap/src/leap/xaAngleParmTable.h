/*
 *      File: xaAngleParmTable.h
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
 *      Author: 	David Rivkin
 *      Date Created:	19 August 1992
 *	Dates Changed:
 *
 *      Description:
 *
 *              Bring up and manage a table based editor for 
 *		force field angle parameter set.
 *
 *
 */
 

#ifndef	HASH_H
# include	"hash.h"
#endif

#include	"basics.h"

/*  Exported functions */


extern void	XAVPTPopupTable(Widget wCreated, Widget wWidget, 
			PARMSET psParmSet, VFUNCTION fCallback);
extern void	XAVPTNewPopupTable( Widget wCreated, Widget wWidget, 
			PARMSET psParmSet, VFUNCTION fCallback ) ;
