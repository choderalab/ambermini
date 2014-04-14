/*
 *	File:	xaAtomTable.c
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
 *		Handle editing of ATOM properties in a table
 *		format.
 *
 *		The ATOMs are found using their ATOM ID's
 *		so that ATOMs can be Destroyed in the UnitEditor
 *		while their properties are being edited within
 *		the TABLE.
 *
 * TODO: Take care of the problem that if the CONTAINER for the TABLE
 * TODO: is Destroyed while it's ATOMs properties are being edited then
 * TODO: the TABLE will attempt to write the properties back to
 * TODO: nowhere.  Maybe fix this by REF/DEREF the container.
 */

#ifndef _xaAtomTable_h_
#define _xaAtomTable_h_

extern Widget	wXAATPopupTable( Widget wCreated, Widget wWidget, UNIT uUnit,
				VFUNCTION fCallback );

#endif /* _xaAtomTable_h_ */


