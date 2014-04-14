/*
* $XConsortium: xPaned.h,v 1.13 91/02/17 13:16:15 rws Exp $
*/


/***********************************************************
Copyright 1987, 1988 by Digital Equipment Corporation, Maynard, Massachusetts,
and the Massachusetts Institute of Technology, Cambridge, Massachusetts.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and its 
documentation for any purpose and without fee is hereby granted, 
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in 
supporting documentation, and that the names of Digital or MIT not be
used in advertising or publicity pertaining to distribution of the
software without specific, written prior permission.  

DIGITAL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL
DIGITAL BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
SOFTWARE.

******************************************************************/

/*
 * xPaned.h - xPaned Composite Widget's public header file.
 *
 * Updated and significantly modifided from the Athena VxPaned Widget.
 *
 * Date:    March 1, 1989
 *
 * By:      Chris D. Peterson
 *          MIT X Consortium
 *          kit@expo.lcs.mit.edu
 */

#ifndef _XxPaned_h
#define _XxPaned_h

#include <X11/Constraint.h>
#include <X11/Xmu/Converters.h>


/****************************************************************
 *
 * Vertical xPaned Widget (SubClass of CompositeClass)
 *
 ****************************************************************/

/* RESOURCES:

 Name		         Class		   RepType	    Default Value
 ----		         -----		   -------	    -------------
 background	         Background	   Pixel	    XtDefaultBackground
 betweenCursor	         Cursor	           Cursor	    **
 border		         BorderColor       Pixel	    XtDefaultForeground
 borderWidth	         BorderWidth       Dimension	    1
 cursor		         Cursor	           Cursor	    None
 destroyCallback         Callback	   Pointer	    NULL
 height		         Height	           Dimension	    0
 gripIndent	         GripIndent	   Position	    16
 gripCursor	         Cursor	           Cursor	    **
 horizontalGripCursol    Cursor	           Cursor	    sb_h_double_arrow
 horizontalBetweencursor Cursor	           Cursor	    sb_up_arrow
 internalBorderColor     BorderColor	   Pixel	    XtDefaultForeground
 internalBorderWidth     BorderWidth	   Position	    1
 leftCursor	         Cursor	           Cursor	    sb_left_arrow
 lowerCursor	         Cursor	           Cursor	    sb_down_arrow
 mappedWhenManaged       MappedWhenManaged Boolean	    True
 orientation             Orientation       XtOrientation    XtorientVertical
 refigureMode	         Boolean	   Boolean	    On
 rightCursor	         Cursor	           Cursor           sb_right_arrow
 sensitive	         Sensitive	   Boolean	    True
 upperCursor	         Cursor	           Cursor	    sb_up_arrow
 verticalBetweenCursor   Cursor	           Cursor           sb_left_arrow
 verticalGripCursor      Cursor	           Cursor           sb_v_double_arrow
 width		         Width	           Dimension	    0
 x		         Position	   Position	    0
 y		         Position	   Position    	    0



 --- V. Romanovski ---
 alongResizable          Boolean           Boolean          False
 acrossResizable         Boolean           Boolean          False
 --- V. Romanovski ---
 
** These resources now are set to the vertical or horizontal cursor
   depending upon orientation, by default.  If a value is specified here
   then that cursor will be used reguardless of orientation.


CONSTRAINT RESOURCES:

 Name		      Class		RepType		Default Value
 ----		      -----		-------		-------------
 allowResize	      Boolean	        Boolean         False
 max		      Max	        Dimension	unlimited
 min		      Min		Dimension	Grip Size
 preferredPaneSize    PerferredPaneSize Dimension	PANED_ASK_CHILD
 resizeToPreferred    Boolean		Boolean	 	False
 showGrip	      ShowGrip		Boolean		True
 skipAdjust	      Boolean	        Boolean         False

 --- V. Romanovski ---        
 order                Order             Int             unlimited
 handled              Boolean           Boolean         True
 --- V. Romanovski ---
  
*/



#define PANED_ASK_CHILD 0
#define PANED_GRIP_SIZE 0

/* New Fields */
#define XtNallowResize "allowResize"
#define XtNbetweenCursor "betweenCursor"
#define XtNverticalBetweenCursor "verticalBetweenCursor"
#define XtNhorizontalBetweenCursor "horizontalBetweenCursor"
#define XtNgripCursor "gripCursor"
#define XtNgripIndent "gripIndent"
#define XtNhorizontalGripCursor "horizontalGripCursor"
#define XtNinternalBorderColor "internalBorderColor"
#define XtNinternalBorderWidth "internalBorderWidth"
#define XtNleftCursor "leftCursor"
#define XtNlowerCursor "lowerCursor"
#define XtNrefigureMode "refigureMode"
#define XtNposition "position"
#define XtNmin "min"
#define XtNmax "max"
#define XtNpreferredPaneSize "preferredPaneSize"
#define XtNresizeToPreferred "resizeToPreferred"
#define XtNrightCursor "rightCursor"
#define XtNshowGrip "showGrip"
#define XtNskipAdjust "skipAdjust"
#define XtNupperCursor "upperCursor"
#define XtNverticalGripCursor "verticalGripCursor"

#define XtCGripIndent "GripIndent"
#define XtCMin "Min"
#define XtCMax "Max"
#define XtCPreferredPaneSize "PreferredPaneSize"
#define XtCShowGrip "ShowGrip"


/* V. Romanovski */
#define XtNorder           "order"
#define XtCOrder           "Order"
#define XtNhandled         "handled"
#define XtChandled         "Chandled"
#define XtNalongResizable  "alongResizable"  
#define XtNacrossResizable "acrossResizable"
#define XtCResizable       "Resizable"
#ifndef XtNrealizeCallback
#define XtNrealizeCallback "realizeCallback"
#endif
/* V. Romanovski */



/* Class record constant */
extern WidgetClass xpanedWidgetClass;

typedef struct _xPanedClassRec	*xPanedWidgetClass;
typedef struct _xPanedRec	*xPanedWidget;

/************************************************************
 *
 *  Public Procedures 
 *
 ************************************************************/

/* _XFUNCPROTOBEGIN */

/*	Function Name: XxPanedSetMinMax
 *	Description: Sets the min and max size for a pane.
 *	Arguments: widget - the widget that is a child of the xPaned widget.
 *                 min, max - the new min and max size for the pane.
 *	Returns: none.
 */

extern void XxPanedSetMinMax(
    Widget		/* w   */,
    int			/* min */,
    int			/* max */
);

/*	Function Name: XxPanedGetMinMax
 *	Description: Gets the min and max size for a pane.
 *	Arguments: widget - the widget that is a child of the xPaned widget.
 ** RETURNED **    min, max - the current min and max size for the pane.
 *	Returns: none.
 */

extern void XxPanedGetMinMax(
    Widget		/* w */,
    int *		/* min_return */,
    int *		/* max_return */
);

/*	Function Name: XxPanedSetRefigureMode
 *	Description: Allows a flag to be set the will inhibit 
 *                   the xpaned widgets relayout routine.
 *	Arguments: w - the xpaned widget.
 *                 mode - if FALSE then inhibit refigure.
 *	Returns: none.
 */

extern void XxPanedSetRefigureMode(
    Widget		/* w */,
    /* Boolean */ int	/* mode */
);

/*	Function Name: XxPanedGetNumSub
 *	Description: Returns the number of panes in the xpaned widget.
 *	Arguments: w - the xpaned widget.
 *	Returns: the number of panes in the xpaned widget.
 */

extern int XxPanedGetNumSub(
    Widget		/* w */
);

/*	Function Name: XxPanedAllowResize
 *	Description: Allows a flag to be set that determines if the xpaned
 *                   widget will allow geometry requests from this child
 *	Arguments: widget - a child of the xpaned widget.
 *	Returns: none.
 */

extern void XxPanedAllowResize(
    Widget		/* w */,
    /* Boolean */ int	/* allow_resize */
);

/* V. Romanovski */
extern void XxPanedChangePaneToOnes( Widget, Widget);
extern void XxPanedPutPaneBefore( Widget, Widget );
/* V. Romanovski */
			    
/* _XFUNCPROTOEND */

#endif /* _XxPaned_h */
