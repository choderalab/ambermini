/*
 *	File: 	action.h
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
 *		Handle the user interaction with the TANK.
 */

extern	XtActionProc	xtapActionSetCancelHit();
extern	XtActionProc	xtapActionCenterUnit();
extern	XtActionProc    xtapActionMouseMotion();
extern	XtActionProc    xtapActionButtonADown();
extern	XtActionProc    xtapActionButtonAUp();
extern	XtActionProc    xtapActionButtonBDown();
extern	XtActionProc    xtapActionButtonBUp();
extern	XtActionProc    xtapActionButtonCDown();
extern	XtActionProc    xtapActionButtonCUp();
extern	XtActionProc    xtapActionButtonDDown();
extern	XtActionProc    xtapActionButtonDUp();
extern	XtActionProc    xtapActionButtonEDown();
extern	XtActionProc    xtapActionButtonEUp();

extern	void		ActionInit();


/* Remember to add any new actions to SxtaaActions[] in tank.c */


/* TANK button states */

#define	ButtonA		0x0001
#define	ButtonB		0x0002
#define	ButtonC		0x0004
#define	ButtonD		0x0008
#define	ButtonE		0x0010

#define	SCALE_SPEED	0.1
#define	TRANSLATE_SPEED	10.0

#define	BSHIFT		ButtonD



