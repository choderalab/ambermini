/***********************************************************

  $XConsortium: xPanedP.h,v 1.5 91/05/09 20:58:23 gildea Exp $

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
 * xPanedP.h - xPaned Composite Widget's private header file.
 *
 * Updated and significantly modified from the Athena VxPaned Widget.
 *
 * Date:    March 1, 1989
 *
 * By:      Chris D. Peterson
 *          MIT X Consortium
 *          kit@expo.lcs.mit.edu
 */

#ifndef _XxPanedP_h
#define _XxPanedP_h

#include "xPaned.h"

/*********************************************************************
 *
 * xPaned Widget Private Data
 *
 *********************************************************************/

/* New fields for the xPaned widget class record */

typedef struct _xPanedClassPart {
    int foo;			/* keep compiler happy. */
} xPanedClassPart;

/* Full Class record declaration */
typedef struct _xPanedClassRec {
    CoreClassPart       core_class;
    CompositeClassPart  composite_class;
    ConstraintClassPart constraint_class;
    xPanedClassPart     xpaned_class;
} xPanedClassRec;

extern xPanedClassRec xpanedClassRec;

/* xPaned constraint record */
typedef struct _xPanedConstraintsPart {
  /* Resources. */
    Dimension	min;		/* Minimum height */
    Dimension	max;		/* Maximum height */
    Boolean	allow_resize;	/* TRUE iff child resize requests are ok */
    Boolean     show_grip;	/* TRUE iff child will have grip below it,
				   when it is not the bottom pane. */
    Boolean	skip_adjust;	/* TRUE iff child's height should not be */
				/* changed without explicit user action. */
    int		position;	/* position location in xPaned (relative to
				   other children) ** NIY ** */
    Dimension   preferred_size;	/* The Preferred size of the pane.
				   Iff this is zero then ask child for size.*/
    Boolean     resize_to_pref;	/* resize this pane to its preferred size
				   on a resize or change managed after 
				   realize. */
    /* V. Romanovski */
    int            order;
    Boolean        handled;
    /* V. Romanovski */
    
    
  /* Private state. */
    Position	delta;		/* Desired Location */
    Position	olddelta;	/* The last value of dy. */
    Boolean     xpaned_adjusted_me; /* Has the vxpaned adjusted this widget w/o
				     user interaction to make things fit? */
    Dimension	wp_size;	/* widget's preferred size */ 
    int         size;		/* the size the widget will actually get. */
    Widget	grip;		/* The grip for this child */

} xPanedConstraintsPart, *Pane;

typedef struct _xPanedConstraintsRec {
    xPanedConstraintsPart xpaned;
} xPanedConstraintsRec, *xPanedConstraints;

/*
 * The Pane Stack Structure.
 */

typedef struct _PaneStack {
    struct _PaneStack * next;	/* The next element on the stack. */
    Pane pane;			/* The pane in this element on the stack. */
    int start_size;		/* The size of this element when it was pushed
				   onto the stack. */
} PaneStack;

/* New Fields for the xPaned widget record */
typedef struct {
    /* resources */
    Position    grip_indent;               /* Location of grips (offset	
					      from right margin) */
    Boolean     refiguremode;              /* Whether to refigure changes 
					      right now */
    XtTranslations grip_translations;      /* grip translation table */
    Pixel       internal_bp;               /* color of internal borders. */
    Dimension   internal_bw;	           /* internal border width. */
    XtOrientation orientation;	           /* Orientation of xpaned widget. */

    Cursor	cursor;		           /* Cursor for xpaned window */
    Cursor	grip_cursor;               /* inactive grip cursor */
    Cursor	v_grip_cursor;             /* inactive vert grip cursor */
    Cursor	h_grip_cursor;             /* inactive horiz grip cursor */
    Cursor	adjust_this_cursor;        /* active grip cursor: T */
    Cursor	v_adjust_this_cursor;      /* active vert grip cursor: T */
    Cursor	h_adjust_this_cursor;      /* active horiz grip cursor: T */

				          /* vertical. */
    Cursor	adjust_upper_cursor;      /* active grip cursor: U */
    Cursor	adjust_lower_cursor;      /* active grip cursor: D */

				          /* horizontal. */
    Cursor	adjust_left_cursor;       /* active grip cursor: U */
    Cursor	adjust_right_cursor;      /* active grip cursor: D */

    /* V. Romanovski */
    
    Boolean     alongResizable;
    Boolean     acrossResizable;
    XtCallbackList realizeCallback;
    
    /* V. Romanovski */

    
    /* private */
    Boolean	recursively_called;        /* for ChangeManaged */
    Boolean	resize_children_to_pref;   /* override constrain resources
					      and resize all children to
					      preferred size. */
    int         start_loc;	           /* mouse origin when adjusting */
    Widget      whichadd;                  /* Which pane to add changes to */
    Widget      whichsub;                  /* Which pane to sub changes from */
    GC          normgc;                    /* GC to use when drawing borders */
    GC          invgc;                     /* GC to use when erasing borders */
    GC          flipgc;                    /* GC to use when animating
					      borders */
    int		num_panes;                 /* count of managed panes */
    PaneStack * stack;		           /* The pane stack for this widget.*/
    
} xPanedPart;

/**************************************************************************
 *
 * Full instance record declaration
 *
 **************************************************************************/

typedef struct _xPanedRec {
    CorePart       core;
    CompositePart  composite;
    ConstraintPart constraint;
    xPanedPart     xpaned;
} xPanedRec;

#endif /* _XxPanedP_h */

