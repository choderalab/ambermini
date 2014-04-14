/*
 *	File:	xaLeap.h
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
 *		Provide prototypes for functions in xaLeap.c
 *		that are available to all X-Windows Athena Widget
 *		specific routines.
 *
 *		XALDirectOutput(wWidget) redirects all VPx screen output to
 *					 the Text Widget wWidget.
 *
 */

#ifndef	XALEAP_H
#define	XALEAP_H	

#include "XrawRegistr.h"

#ifndef _
#define	_(argsInParens)	argsInParens	/* ANSI --> int foo ( int, char* ); */
#endif


#define IS_GADGET(a)  (                                               \
		       ((Widget)(a)!=(Widget)NULL                     \
			&& XtIsWidget(a)                              \
			&& !((Widget)(a))->core.being_destroyed       \
			)                                             \
		       ||((RectObj)(a)!=(RectObj)NULL                 \
			  && XtIsRectObj(a)                           \
			  && !((RectObj)(a))->object.being_destroyed) \
		       ||((Object)(a)!=(Object)NULL                   \
			  && XtIsObject(a)                            \
			  && !((Object)(a))->object.being_destroyed)  \
		       )

#define GetWidgetByName(name) (Widget)WcFullNameToWidget(GwTopWidget,(name))
#define UnInt     unsigned int
#define UnCh      unsigned char


extern	Widget		GwTopWidget;
extern	Widget		GwCmd;
extern  XtAppContext	GxtacApp;
extern	int		GiCmdSink;
extern	Display		*GdDisplay;
extern	XrmDatabase	GXrmDbase, GXrmMyDbase;
extern	char		*GsResources;
extern	Screen		*GsScreen;
extern	Window		GwRootWindow;
extern  GC              GgcDark;
extern  GC              GgcLight;

extern	void XALPrintStringToWidget( char*, Widget ) ;
extern  Widget XAGetWidgetFromString(Widget, char** ) ;
#endif	/* XALEAP_H */









