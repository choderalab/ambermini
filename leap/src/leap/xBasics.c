/*
 *	File - xBasics.c
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
 *	Basic X Windows functions
 *
 */

#include <X11/X.h>
#include <X11/Intrinsic.h>
#include <X11/Xlib.h>

#include "basics.h"

/*
 *	RedrawWidget -
 *
 *	Tells a widget to redraw (expose) itself, and hopefully 
 *	all its children.
 *
 */

void
RedrawWidget( Widget wWidget )
{
XEvent	eeEvent;

    eeEvent.type = Expose;
    eeEvent.xexpose.display = XtDisplay( wWidget );
    eeEvent.xexpose.window = XtWindow( wWidget );
    
    XtDispatchEvent ( &eeEvent );
}
 

/*
 *	RaiseWidget( wWidget )
 *
 */

void
RaiseWidget( Widget wWidget )
{
    XRaiseWindow( XtDisplay( wWidget ), XtWindow( wWidget ));
}

	
