/* $XConsortium: Viewport.c,v 1.68 91/07/24 18:56:11 converse Exp $ */

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

#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#include <X11/Composite.h>

#include "XawInit.h"
#include "../Xmu/Misc.h"
#include "Scrollbar.h"
#include "ViewportP.h"

#include "XrawDebug.h"


#define VW vw->viewport
#define DEFAULT_SCROLL_WIDTH 15
#define MIN_SIZE (DEFAULT_SCROLL_WIDTH + 2 * VW.distance)

#define X_MARGIN(vw) (VW.distance)
#define Y_MARGIN(vw) (VW.distance)

#define IS_MANAGED(w) (w && XtIsManaged(w))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

static void DefDistance();

#define Offset(field) XtOffsetOf(ViewportRec, viewport.field)
static XtResource resources[] = {
    {
      XtNframeWidth, XtCFrameWidth, XtRDimension, sizeof(Dimension),
      Offset(frame_width), XtRImmediate, (XtPointer)2
    },
    {
      XtNscrollbarWidth, XtCScrollbarWidth, XtRDimension, sizeof(Dimension),
      Offset(scrollbar_width), XtRImmediate, (XtPointer)DEFAULT_SCROLL_WIDTH
    },
    {
      XtNframeType, XtCFrameType, XtRFrameType, sizeof(XawFrameType),
      Offset(frame_type), XtRImmediate, (XtPointer) XawSUNKEN
    },
    {
      XtNforceBars, XtCBoolean, XtRBoolean, sizeof(Boolean),
      Offset(forcebars), XtRImmediate, (XtPointer)False
    },
    {
      XtNallowHoriz, XtCBoolean, XtRBoolean, sizeof(Boolean),
      Offset(allowhoriz), XtRImmediate, (XtPointer)False
    },
    {
      XtNallowVert, XtCBoolean, XtRBoolean, sizeof(Boolean),
      Offset(allowvert), XtRImmediate, (XtPointer)False
    },
    {
      XtNuseBottom, XtCBoolean, XtRBoolean, sizeof(Boolean),
      Offset(usebottom), XtRImmediate, (XtPointer)True
    },
    {
      XtNuseRight, XtCBoolean, XtRBoolean, sizeof(Boolean),
      Offset(useright), XtRImmediate, (XtPointer)True
    },
    {
      XtNreportCallback, XtCReportCallback, XtRCallback, sizeof(XtPointer),
      Offset(report_callbacks), XtRImmediate, (XtPointer) NULL
    },
    {
      XtNdistance, XtCDistance, XtRDimension, sizeof(Dimension),
      Offset(distance), XtRCallProc, (XtPointer)DefDistance
    },
};

/* ARGSUSED */
static void DefDistance( w, offset, value)
     Widget w;
     int offset;
     XrmValue *value;
{
  ViewportWidget vw = (ViewportWidget) w;
  static Dimension distance;

  distance = 2 * VW.frame_width + 4;

  value->addr = (caddr_t)&distance;
}

#undef Offset


static void ThumbProc();


static void Initialize();
static void ConstraintInitialize();
static void Realize();
static void Resize();
static void ChangeManaged();
static Boolean SetValues();
static void Redisplay();

static void Layout();
static XtGeometryResult GeometryManager(), PreferredGeometry();
//static XtGeometryResult QueryGeometry();
//static Boolean GetGeometry();

#define SuperClass	(&containerClassRec)
ViewportClassRec viewportClassRec = {
  { /* core_class fields */
    /* superclass	  */	(WidgetClass) SuperClass,
    /* class_name	  */	"Viewport",
    /* widget_size	  */	sizeof(ViewportRec),
    /* class_initialize	  */	XawInitializeWidgetSet,
    /* class_part_init    */    NULL,
    /* class_inited	  */	FALSE,
    /* initialize	  */	Initialize,
    /* initialize_hook    */    NULL,
    /* realize		  */	Realize,
    /* actions		  */	NULL,
    /* num_actions	  */	0,
    /* resources	  */	resources,
    /* num_resources	  */	XtNumber(resources),
    /* xrm_class	  */	NULLQUARK,
    /* compress_motion	  */	TRUE,
    /* compress_exposure  */	TRUE,
    /* compress_enterleave*/    TRUE,
    /* visible_interest	  */	FALSE,
    /* destroy		  */	NULL,
    /* resize		  */	Resize,
    /* expose		  */	Redisplay,
    /* set_values	  */	SetValues,
    /* set_values_hook    */    NULL,
    /* set_values_almost  */    XtInheritSetValuesAlmost,
    /* get_values_hook    */	NULL,
    /* accept_focus	  */	NULL,
    /* version            */    XtVersion,
    /* callback_private	  */	NULL,
    /* tm_table    	  */	NULL,
    /* query_geometry     */    PreferredGeometry,
    /* display_accelerator*/	XtInheritDisplayAccelerator,
    /* extension          */	NULL
  },
  { /* composite_class fields */
    /* geometry_manager	  */	GeometryManager,
    /* change_managed	  */	ChangeManaged,
    /* insert_child	  */	XtInheritInsertChild,
    /* delete_child	  */	XtInheritDeleteChild,
    /* extension          */	NULL
  },
  { /* constraint_class fields */
    /* subresourses	  */	NULL,
    /* subresource_count  */	0,
    /* constraint_size	  */	sizeof(ViewportConstraintsRec),
    /* initialize	  */	ConstraintInitialize,
    /* destroy		  */	NULL,
    /* set_values	  */	NULL,
    /* extension          */	NULL
  },
    /* Container class part */
  {
    /* unused            */     0,
  },
  { /* viewport_class fields */
    /* empty		  */	0
  }
};


WidgetClass viewportWidgetClass = (WidgetClass)&viewportClassRec;

static Widget CreateScrollbar (vw, vertical)
     ViewportWidget vw;
     Boolean vertical;
{
  Widget scroll;
  
  scroll = XtVaCreateWidget 
    ((vertical ? "vertical" : "horizontal"),
     scrollbarWidgetClass, (Widget)vw,
     XtNborderWidth, 0,
     XtNorientation, ( vertical ? XtorientVertical : XtorientHorizontal ),
     XtNhighlightThickness, 0,
     NULL);
  
  if (vertical) 
    VW.vert_bar = scroll;
  else 
    VW.horiz_bar = scroll;

  return scroll;
}

/* ARGSUSED */
static void Initialize(request,new,args,n_args)
     Widget   request; /* unused */
     Widget   new;     
     ArgList  args;    /* unused */
     Cardinal *n_args; /* unused */
{
  register ViewportWidget vw = (ViewportWidget)new;
  static Arg clip_args[5];
  Cardinal num_args;
  XtWidgetGeometry preferred;

  /* 
   * Initialize all widget pointers to NULL.
   */

  VW.not_manage = False;
  VW.child      = WNULL;
  VW.horiz_bar  = WNULL;
  VW.vert_bar   = WNULL;

  
  /* 
   * Create Clip Widget.
   */

  num_args = 0;
  XtSetArg(clip_args[num_args], XtNborderWidth,  0);    num_args++;
  XtSetArg(clip_args[num_args], XtNshadowWidth,  0);    num_args++;
  XtSetArg(clip_args[num_args], XtNwidth,        1);    num_args++;
  XtSetArg(clip_args[num_args], XtNheight,       1);    num_args++;

  VW.clip = XtCreateWidget("clip", coreWidgetClass, new, clip_args, num_args);


  /*
   * Create scrollbar which is needed.
   */

  if (VW.allowhoriz) 
    (void) CreateScrollbar(vw, False);

  if (VW.allowvert) 
    (void) CreateScrollbar(vw, True);
  

  /* 
   * Initialize geometry.
   */

  _XawQueryGeometry (new, &preferred);

  if (vw->core.width == 0)
    vw->core.width = preferred.width;
  
  if (vw->core.height == 0)
    vw->core.height = preferred.height;

}


static void SendReport (vw, changed)
    ViewportWidget vw;
    unsigned int changed;
{
  XawPannerReport rep;
  
  if (XtCallbackHasSome == XtHasCallbacks ((Widget)vw, XtNreportCallback))
  {
    register Widget child = VW.child;
    register Widget clip = VW.clip;
    
    rep.changed       = changed;
    rep.slider_x      = -child->core.x;	/* child is canvas */
    rep.slider_y      = -child->core.y;	/* clip is slider */
    rep.slider_width  = clip->core.width;
    rep.slider_height = clip->core.height;
    rep.canvas_width  = child->core.width;
    rep.canvas_height = child->core.height;

    XtCallCallbackList ((Widget) vw, VW.report_callbacks, (XtPointer) &rep);
  }
}

#define VERT_BAR_NEEDED  (1L<<0)
#define HORIZ_BAR_NEEDED (1L<<1)


static int long WhatBarsNeeded(vw)
     ViewportWidget vw;
{
  register Widget child = VW.child;
  Dimension scroll_cross;
  int long bars_needed;
  int stuff_w;
  int stuff_h;

  scroll_cross = VW.scrollbar_width + VW.distance;

  stuff_w = vw->core.width - (scroll_cross + 2 * X_MARGIN(vw));
  stuff_h = vw->core.height - (scroll_cross + 2 * Y_MARGIN(vw));



  /* 
   * Correlation full size of child widget and clip one 
   */
#define V_LESS     (child == WNULL ? True : (stuff_h >= FULL_HEIGHT(child)))
#define V_MORE     (child == WNULL ? False :  \
		    ((stuff_h + scroll_cross) < FULL_HEIGHT(child)))

#define H_LESS     (child == WNULL ? True : (stuff_w >= FULL_WIDTH(child)))
#define H_MORE     (child == WNULL ? False :  \
		    ((stuff_w + scroll_cross) < FULL_WIDTH(child)))

  /* 
   * Grades of necessity of scrollbars 
   */
#define V_NO_BAR   (!VW.allowvert || VW.vert_bar == WNULL )
#define V_NO_FORCE (VW.allowvert && (!VW.forcebars))
#define V_FORCE    (VW.allowvert && VW.forcebars)

#define V_MUST     !V_NO_BAR && (V_FORCE  || (V_NO_FORCE && V_MORE))
#define V_REJECT   (V_NO_BAR || (V_NO_FORCE && V_LESS))

#define H_NO_BAR   (!VW.allowhoriz || VW.horiz_bar == WNULL)
#define H_NO_FORCE (VW.allowhoriz && (!VW.forcebars))
#define H_FORCE    (VW.allowhoriz && VW.forcebars)

#define H_MUST     !H_NO_BAR && (H_FORCE  || (H_NO_FORCE && H_MORE))
#define H_REJECT   (H_NO_BAR || (H_NO_FORCE && H_LESS))


  bars_needed = 0;

  if (H_MUST)
  {
    bars_needed |= HORIZ_BAR_NEEDED;
    if (V_MUST)
    {
      bars_needed |= VERT_BAR_NEEDED;
    }
    else if (V_REJECT)
    {
      bars_needed &= ~VERT_BAR_NEEDED;
    }
    else
    {
      bars_needed |= VERT_BAR_NEEDED;
    }
  }
  else if (H_REJECT)
  {
    bars_needed &= ~HORIZ_BAR_NEEDED;
    if (V_MUST)
    {
      bars_needed |= VERT_BAR_NEEDED;
    }
    else if (V_REJECT)
    {
      bars_needed &= ~VERT_BAR_NEEDED;
    }
    else
    {
      bars_needed &= ~VERT_BAR_NEEDED;
    }
  }
  else
  {
    if (V_MUST)
    {
      bars_needed = HORIZ_BAR_NEEDED | VERT_BAR_NEEDED;
    }
  }

  return bars_needed;
}


static void Layout(vw)
     ViewportWidget vw;
{
  register Widget child = VW.child;
  register Widget clip = VW.clip;
  int long bars_needed;
  float top;
  float shown;
  Dimension scroll_cross;
  Position x;
  Position y;
  Position stuff_x;
  Position stuff_y;
  int      stuff_w;
  int      stuff_h;


  bars_needed = WhatBarsNeeded(vw);

  scroll_cross = VW.scrollbar_width + VW.distance;

  /* 
   * Default geometry of clip widget
   */
  stuff_x = X_MARGIN(vw);
  stuff_y = Y_MARGIN(vw);
  stuff_w = vw->core.width - 2 * stuff_x;
  stuff_h = vw->core.height -  2 * stuff_y;


  /* 
   * Vertical bar must be on screen
   */
  if (bars_needed & VERT_BAR_NEEDED)
  {
    int w;
    int h;

    x = vw->core.width - scroll_cross;
    y = Y_MARGIN(vw);
    w = VW.scrollbar_width;
    h = vw->core.height - 2 * y;


    if (bars_needed & HORIZ_BAR_NEEDED)
    {
      h -= scroll_cross;

      if (!VW.usebottom)
	y += scroll_cross;
    }

    if (!VW.useright) 
    {
      x = VW.distance;
      stuff_x = X_MARGIN(vw) + scroll_cross;
    }

    h = MAX (h, MIN_SIZE);

    stuff_y = y;
    stuff_h = h;
    stuff_w -= VW.scrollbar_width + VW.distance;

    XtConfigureWidget (VW.vert_bar, x, y, (Dimension)w, (Dimension)h, 0);
  }
  else if (bars_needed & HORIZ_BAR_NEEDED)
    stuff_h -= VW.scrollbar_width + VW.distance;


  /* 
   * Horizontal bar must be on screen
   */
  if (bars_needed & HORIZ_BAR_NEEDED)
  {
    int w;
    int h;

    x = X_MARGIN(vw);
    y = vw->core.height-scroll_cross;
    w = vw->core.width - 2 * x;
    h = VW.scrollbar_width;

    if (bars_needed & VERT_BAR_NEEDED)
    {
      w -= scroll_cross;

      if (!VW.useright)
	x += scroll_cross;
    }

    if (!VW.usebottom) 
    {
      y = VW.distance;
      stuff_y = X_MARGIN(vw) + scroll_cross;
    }

    w = MAX (w, MIN_SIZE);

    stuff_x = x;
    stuff_w = w;
    
    XtConfigureWidget (VW.horiz_bar, x, y, (Dimension)w, (Dimension)h, 0);
  }


  /* 
   * Do not allow to shrink clip widget less then MIN_SIZE
   */
  stuff_w = MAX (stuff_w, MIN_SIZE);
  stuff_h = MAX (stuff_h, MIN_SIZE);


  /* 
   * Replace clip widget
   */
  XtConfigureWidget (VW.clip, stuff_x, stuff_y, 
		     (Dimension)stuff_w, (Dimension)stuff_h, 0);


  /* 
   * Do not allow to shrink child widget less then clip one
   */
  stuff_w = MAX (stuff_w, VW.child->core.width);
  stuff_h = MAX (stuff_h, VW.child->core.height);

  if ( XtIsRealized(VW.child))
    XtResizeWidget (VW.child, (Dimension)stuff_w, (Dimension)stuff_h, 
		    VW.child->core.border_width);
  
  /* 
   * Replace child widget in vertical direction
   */
  top   = 0.0;
  shown = 0.1;
  
  if (IS_MANAGED(VW.child))
    shown = (float)VW.clip->core.height / (float)VW.child->core.height;
  
  if (shown >= 1.0) 
  {
    shown = 1.0;
    top   = 0.0;
  }
  else if (IS_MANAGED(VW.child)) 
  {
    top = (float)VW.child->core.y / (float)VW.child->core.height;
    top = -top;
    top = MIN (top, 1.0 - shown);
  }
  
  if (IS_MANAGED(VW.child)) 
  {
    y = (Position)((float)VW.child->core.height * top);
    
    /* make sure we never move past border */
    if (y + (int)clip->core.height > (int)child->core.height)
      y = child->core.height - clip->core.height;
    
    if (y < 0) y = 0;
    
    XtMoveWidget (VW.child, VW.child->core.x, -y);
  }

  /* 
   * Control for vertical bar
   */
  if (VW.vert_bar && XtIsRealized(VW.vert_bar))
  {
    if (bars_needed & VERT_BAR_NEEDED)
    {
      XtAddCallback(VW.vert_bar, XtNvalueChangedProc, ThumbProc, (Widget)vw);
      XtAddCallback(VW.vert_bar, XtNdragProc, ThumbProc, (Widget)vw);
  
      XMapWindow(XtDisplay(VW.vert_bar),XtWindow(VW.vert_bar));
      XawScrollbarSetThumb(VW.vert_bar, top, shown);
      SendReport (vw, XawPRSliderY);
    }
    else
    {
      XtRemoveCallback(VW.vert_bar, XtNvalueChangedProc, 
		       ThumbProc, (Widget)vw);
      XtRemoveCallback(VW.vert_bar, XtNdragProc, 
		       ThumbProc, (Widget)vw);

      XUnmapWindow(XtDisplay(VW.vert_bar),XtWindow(VW.vert_bar));
    }
  }

  /* 
   * Replace child widget in vertical direction
   */
  top = 0.0;
  shown = 0.1;
  
  if (IS_MANAGED(VW.child))
    shown = (float)VW.clip->core.width / (float)VW.child->core.width;
  
  if (shown >= 1.0) 
  {
    shown = 1.0;
    top   = 0.0;
  }
  else if (IS_MANAGED(VW.child)) 
  {
    top = (float)VW.child->core.x / (float)VW.child->core.width;
    top = -top;
    top = MIN (top, 1.0 - shown);
  }
  
  if (IS_MANAGED(VW.child)) 
  {
    x = (Position)((float)VW.child->core.width * top);
    
    /* make sure we never move past border */
    if (x + (int)clip->core.width > (int)child->core.width)
      x = child->core.width - clip->core.width;
    
    if (x < 0) x = 0;
    
    XtMoveWidget (VW.child, -x , VW.child->core.y );
  }


  /* 
   * Control for hirizontal bar 
   */
  if (VW.horiz_bar && XtIsRealized(VW.horiz_bar))
  {
    if (bars_needed & HORIZ_BAR_NEEDED)
    {
      XtAddCallback(VW.horiz_bar, XtNvalueChangedProc, ThumbProc, (Widget)vw);
      XtAddCallback(VW.horiz_bar, XtNdragProc, ThumbProc, (Widget)vw);
  
      XMapWindow(XtDisplay(VW.horiz_bar),XtWindow(VW.horiz_bar));
      XawScrollbarSetThumb(VW.horiz_bar, top, shown);
      SendReport (vw, XawPRSliderX);
    }
    else
    {
      XtRemoveCallback(VW.horiz_bar, XtNvalueChangedProc, 
		       ThumbProc, (Widget)vw);
      XtRemoveCallback(VW.horiz_bar, XtNdragProc, 
		       ThumbProc, (Widget)vw);

      XUnmapWindow(XtDisplay(VW.horiz_bar),XtWindow(VW.horiz_bar));
    }
  }

}



/* ARGSUSED */
static void ConstraintInitialize(request, new)
    Widget request, new;
{
  ((ViewportConstraints)new->core.constraints)->viewport.reparented = False;
}


static void Realize(widget, value_mask, attributes)
    Widget widget;
    XtValueMask *value_mask;
    XSetWindowAttributes *attributes;
{
  ViewportWidget vw = (ViewportWidget)widget;
  register Widget child = vw->viewport.child;
  register Widget clip = vw->viewport.clip;
  
  *value_mask |= CWBitGravity;
  attributes->bit_gravity = NorthWestGravity;
  (*SuperClass->core_class.realize)(widget, value_mask, attributes);
  
  XtManageChild (clip);
  XtRealizeWidget (clip);

  if (child != WNULL) 
  {
    XtMoveWidget (child, (Position)0, (Position)0 );
    XtRealizeWidget (child);
    XReparentWindow (XtDisplay(widget), XtWindow(child), XtWindow(clip),
		    (Position)0, (Position)0 );

    if (child->core.mapped_when_managed)
      XtMapWidget (child);
  }

  Layout(vw);
}

/* ARGSUSED */
static Boolean SetValues(current, request, new, args, num_args)
    Widget current, request, new;
    ArgList args;
    Cardinal *num_args;
{
    ViewportWidget vw = (ViewportWidget)new;
    ViewportWidget cw = (ViewportWidget)current;

#define NE(field) (vw->viewport.field != cw->viewport.field)
    if (NE(forcebars)  ||
	NE(allowvert)  ||
	NE(allowhoriz) ||
	NE(useright)   ||
	NE(usebottom)  ||
	NE(distance)) 
    {
      (*vw->core.widget_class->core_class.resize)(new);
    }

    return False;

#undef NE
}

static void ChangeManaged(widget)
    Widget widget;
{
  register ViewportWidget vw = (ViewportWidget)widget;
  register Widget  child, *childP;
  int i;
  int num_children = vw->composite.num_children;
  XtWidgetGeometry vw_request;

  if (VW.not_manage)
    return;

  VW.not_manage = True;

  child = WNULL;

  for (childP=vw->composite.children, i=0; i < num_children; childP++, i++) 
  {
    if (XtIsManaged(*childP)
	&& *childP != VW.clip
	&& *childP != VW.horiz_bar
	&& *childP != VW.vert_bar)
    {
      if (child == WNULL) 
	child = *childP;
      else
	XtUnmanageChild(*childP);
    }
  }
  

  if (child != VW.child) 
  {
    VW.child = child;
    
    if (child != WNULL) 
    {
      if (XtIsRealized(widget)) 
      {
	ViewportConstraints constraints =
	  (ViewportConstraints)child->core.constraints;
	
	if (!XtIsRealized(child)) 
        {
	  XtMoveWidget (child, (Position)0, (Position)0);

	  XtRealizeWidget (child);

	  XReparentWindow (XtDisplay(widget), XtWindow(child),
			  XtWindow(VW.clip), (Position)0, (Position)0);

	  if (child->core.mapped_when_managed)
	    XtMapWidget (child);
	  
	  constraints->viewport.reparented = True;
	}
	else if (!constraints->viewport.reparented) 
        {
	  XReparentWindow (XtDisplay(widget), XtWindow(child),
			  XtWindow(VW.clip), (Position)0, (Position)0 );

	  constraints->viewport.reparented = True;

	  if (child->core.mapped_when_managed)
	    XtMapWidget (child);
	}
      }
    }
  }

  
  _XawQueryGeometry((Widget)vw, &vw_request);

  while (XtGeometryAlmost ==
	 XtMakeGeometryRequest((Widget)vw, &vw_request, &vw_request))
    {/* EMPTY */}

  Layout(vw);
  VW.not_manage = False;;
}


static void SetBar(w, top, length, total)
    Widget w;
    Position top;
    Dimension length, total;
{
  float top_f;
  float shown;

  top_f = (float)top/(float)total;
  shown = (float)length/(float)total;

  XawScrollbarSetThumb(w, top_f, shown);
}

static void RedrawThumbs(vw)
  ViewportWidget vw;
{

  register Widget child = VW.child;
  register Widget clip = VW.clip;
  
  if (VW.horiz_bar != WNULL)
    SetBar( VW.horiz_bar, -(child->core.x),
	   clip->core.width, child->core.width );
  
  if (VW.vert_bar != WNULL)
    SetBar( VW.vert_bar, -(child->core.y),
	   clip->core.height, child->core.height );

}



static void MoveChild(vw, x, y, flag)
     ViewportWidget vw;
     Position x, y;
     Boolean  flag;
{
    register Widget child = VW.child;
    register Widget clip = VW.clip;

    /* make sure we never move past right/bottom borders */
    if (-x + (int)clip->core.width > (int)child->core.width)
	x = -(child->core.width - clip->core.width);

    if (-y + (int)clip->core.height > (int)child->core.height)
	y = -(child->core.height - clip->core.height);

    /* make sure we never move past left/top borders */
    if (x >= 0) x = 0;
    if (y >= 0) y = 0;

    XtMoveWidget(child, x, y);
    SendReport (vw, (XawPRSliderX | XawPRSliderY));

    if (flag)
      RedrawThumbs(vw);
}

static void Resize(widget)
    Widget widget;
{
  ViewportWidget vw = (ViewportWidget) widget;
  XtWidgetGeometry query;

  if (VW.child != WNULL)
  {
    _XawQueryGeometry (VW.child, &query);
    XtResizeWidget (VW.child, query.width, query.height, 
		    VW.child->core.border_width);
  }

  Layout ((ViewportWidget)widget);

  if (XtIsRealized(widget))
  {
    XClearWindow (XtDisplay(widget), XtWindow(widget));
    Redisplay (widget);
    XFlush (XtDisplay(widget));
  }
}

/* ARGSUSED */
static void Redisplay(gw,  event, region)
     Widget gw;
     XEvent * event;                 /* unused. */
     Region region;                  /* unused. */
{
  ViewportWidget vw = (ViewportWidget) gw;
  Position x;
  Position y;
  Dimension w;
  Dimension h;

  x = VW.clip->core.x - VW.frame_width;
  y = VW.clip->core.y - VW.frame_width;
  w = VW.clip->core.width + 2 * VW.frame_width;
  h = VW.clip->core.height + 2 * VW.frame_width;

  (*SuperClass->core_class.expose) (gw, event, region);
  
  XawDrawFrame (gw, x, y, w, h, 
		VW.frame_type,
		VW.frame_width,
		vw->container.top_shadow_GC,
		vw->container.bottom_shadow_GC);
}



static void ThumbProc(widget, client_data, call_data)
    Widget    widget;
    XtPointer client_data;
    XtPointer call_data;
{
  ViewportWidget vw = (ViewportWidget) client_data;
  XawScrollBarCallbackStruct *scrollBarStruct;
  register Widget child = VW.child;
  register Widget clip = VW.clip;
  float percent;
  Position x;
  Position y;
    

  scrollBarStruct = (XawScrollBarCallbackStruct*)call_data;
  percent = scrollBarStruct->top;

  if (widget == VW.horiz_bar)
  {
    if (IS_MANAGED(VW.child)) {
      x = (Position)((float)VW.child->core.width * percent);
      
      /* make sure we never move past border */
      if (x + (int)clip->core.width > (int)child->core.width)
	x = child->core.width - clip->core.width;

      if (x < 0) x = 0;

      XtMoveWidget (VW.child, -x , VW.child->core.y );
      SendReport (vw, XawPRSliderX);
    }
  }
  else if (widget == VW.vert_bar)
  {
    if (IS_MANAGED(VW.child)) {
      y = (Position)((float)VW.child->core.height * percent);

      /* make sure we never move past border */
      if (y + (int)clip->core.height > (int)child->core.height)
	y = child->core.height - clip->core.height;
      
      if (y < 0) y = 0;

      XtMoveWidget (VW.child, VW.child->core.x, -y);
      SendReport (vw, XawPRSliderY);
    }
  }

}


#if 0
static XtGeometryResult TestSmaller(vw, request, reply_return)
     ViewportWidget vw; 
     XtWidgetGeometry *request, *reply_return;
{
  if (request->width < vw->core.width || request->height < vw->core.height)
    return XtMakeGeometryRequest((Widget)vw, request, reply_return);
  else
    return XtGeometryYes;  
}

static XtGeometryResult
GeometryRequestPlusScrollbar(vw, horizontal, request, reply_return)
     Boolean horizontal;
     ViewportWidget vw; 
     XtWidgetGeometry *request, *reply_return;
{
  Widget sb;
  XtWidgetGeometry plusScrollbars;
  plusScrollbars = *request;
  if ((sb = VW.horiz_bar) == (Widget)NULL)
    sb = CreateScrollbar( vw, horizontal);
  request->width += sb->core.width;
  request->height += sb->core.height;
  XtDestroyWidget(sb);
  return XtMakeGeometryRequest((Widget) vw, &plusScrollbars, reply_return);
 }
#endif

#define WidthChange() (request->width != vw->core.width)
#define HeightChange() (request->height != vw->core.height)

#if 0
static XtGeometryResult QueryGeometry(vw, request, reply_return)
     ViewportWidget vw; 
     XtWidgetGeometry *request, *reply_return;
{	
  if (VW.allowhoriz && VW.allowvert) 
    return TestSmaller(vw, request, reply_return);

  else if (VW.allowhoriz && !VW.allowvert) {
    if (WidthChange() && !HeightChange())
      return TestSmaller(vw, request, reply_return);
    else if (!WidthChange() && HeightChange())
      return XtMakeGeometryRequest((Widget) vw, request, reply_return);
    else if (WidthChange() && HeightChange()) /* hard part */
      return GeometryRequestPlusScrollbar(vw, True, request, reply_return);
    else /* !WidthChange() && !HeightChange() */
      return XtGeometryYes;
  }
  else if (!VW.allowhoriz && VW.allowvert) {
    if (!WidthChange() && HeightChange())
      return TestSmaller(vw, request, reply_return);
    else if (WidthChange() && !HeightChange())
      return XtMakeGeometryRequest((Widget)vw, request, reply_return);
    else if (WidthChange() && HeightChange()) /* hard part */
      return GeometryRequestPlusScrollbar(vw, False, request, reply_return);
    else /* !WidthChange() && !HeightChange() */
      return XtGeometryYes;
  }      
  else /* (!VW.allowhoriz && !VW.allowvert) */
    return XtMakeGeometryRequest((Widget) vw, request, reply_return);
}
#endif

#undef WidthChange
#undef HeightChange



/*
 *      GetSurmiseViewportGeometry 
 *      
 * Calculate surmise size for Viewport widget 
 * base on child size.
 *
 */

static void GetSurmiseViewportGeometry(child, width, height)
     Widget     child;
     Dimension *width;      /* return */
     Dimension *height;     /* return */
{
  register ViewportWidget vw = (ViewportWidget) XtParent(child);
  Dimension child_width  = *width;
  Dimension child_height = *height;

  *width  += 2 * X_MARGIN(vw);
  *height += 2 * Y_MARGIN(vw);

  if (FULL_WIDTH(child) < child_width)
    *height += VW.distance + VW.scrollbar_width;    /* exactly height */

  if (FULL_HEIGHT(child) < child_height)
    *width  += VW.distance + VW.scrollbar_width;    /* exactly width */
}


#if 0
/*
 *      GetSurmiseChildGeometry 
 *      
 * Calculate surmise size for child widget 
 * base on size of Viewport widget.
 *
 */


static void GetSurmiseChildGeometry(child, view_width, view_height,
				    child_request_width, child_request_height)
     Widget     child;
     Dimension *view_width;      /* return */
     Dimension *view_height;     /* return */
     Dimension  child_request_width;
     Dimension  child_request_height;
{
  register ViewportWidget vw = (ViewportWidget) XtParent(child);
  int width  = *view_width;
  int height = *view_height;
  int long bars_needed;

  bars_needed = WhatBarsNeeded(vw);

  if (!(bars_needed & HORIZ_BAR_NEEDED))
  {
    width  -= 2 * X_MARGIN(vw);
    
    if (bars_needed & VERT_BAR_NEEDED)
      width  -= VW.distance + VW.scrollbar_width;    /* exactly width */
  }
  else
  {
    width = child_request_width;
  }

  if (!(bars_needed & VERT_BAR_NEEDED))
  {
    height -= 2 * Y_MARGIN(vw);
    
    if (bars_needed & HORIZ_BAR_NEEDED)
      height -= VW.distance + VW.scrollbar_width;    /* exactly height */
  }
  else
  {
    height = child_request_height;
  }

  *view_width  = width  > 0 ? width  : 1;
  *view_height = height > 0 ? height : 1;
}
#endif

#define SET_Width(request) (request->request_mode & CWWidth)
#define SET_Height(request) (request->request_mode & CWHeight)

/* ARGSUSED */
static XtGeometryResult GeometryManager(child, request, reply)
     Widget            child;
     XtWidgetGeometry *request;
     XtWidgetGeometry *reply;    /* unused */
{
  register ViewportWidget vw = (ViewportWidget) XtParent(child);
  Dimension width;
  Dimension height;
  Boolean relayout = False;
  Boolean ask_parent_for_width = False;
  Boolean ask_parent_for_height = False;
  XtGeometryResult result;

  if (child != VW.child)
    return XtGeometryNo;
  
  if (SET_Width(request) && 
      (request->width < VW.clip->core.width  || !VW.allowhoriz))
    ask_parent_for_width = XtIsRealized((Widget)vw) || vw->core.width == 0;

  if (SET_Height(request) &&
      (request->height < VW.clip->core.height || !VW.allowvert))
    ask_parent_for_height = XtIsRealized((Widget)vw) || vw->core.height == 0;
  
  /* 
   * Make geometry management with viewport's parent participation..
   */

  if (ask_parent_for_width || ask_parent_for_height)
  {
    XtWidgetGeometry intended;

    if (SET_Width(request))      width = request->width;
    else                         width = child->core.width;
    
    if (SET_Height(request))     height = request->height;
    else                         height = child->core.height;
    
    GetSurmiseViewportGeometry(child, &width, &height);

    intended.width  = width;
    intended.height = height;
    intended.request_mode = 0;

    if (ask_parent_for_width)
      intended.request_mode |= CWWidth;

    if (ask_parent_for_height)
      intended.request_mode |= CWHeight;

    if (request->request_mode & XtCWQueryOnly)
      intended.request_mode |= XtCWQueryOnly;
    
    do {
      result = XtMakeGeometryRequest((Widget)vw, &intended, &intended);
    } while(result == XtGeometryAlmost);
  }


  if (SET_Width(request))    width = request->width;
  else                       width = child->core.width;
  
  if (SET_Height(request))   height = request->height;
  else                       height = child->core.height;
  
  if (!(request->request_mode & XtCWQueryOnly))
    XtResizeWidget (child, width, height, child->core.border_width);

  if ((request->request_mode & CWWidth) || (request->request_mode & CWHeight))
    relayout = True;
  
  if (relayout && !(request->request_mode & XtCWQueryOnly))
  { 
    Layout (vw);
    
    if (XtIsRealized((Widget)vw))
    {
      XClearWindow (XtDisplay((Widget)vw), XtWindow((Widget)vw));
      Redisplay ((Widget)vw);
      XFlush (XtDisplay((Widget)vw));
    }
  }

  return XtGeometryYes;
}

#if 0
static XtGeometryResult _GeometryManager(child, request, reply)
    Widget child;
    XtWidgetGeometry *request, *reply;
{
    ViewportWidget vw = (ViewportWidget)child->core.parent;
    Boolean rWidth = (Boolean)(request->request_mode & CWWidth);
    Boolean rHeight = (Boolean)(request->request_mode & CWHeight);
    XtWidgetGeometry allowed;
    XtGeometryResult result;
    Boolean reconfigured;
    Boolean child_changed_size;
    Dimension height_remaining;

    if (request->request_mode & XtCWQueryOnly)
      return QueryGeometry(vw, request, reply);

    if (child != VW.child
        || request->request_mode & ~(CWWidth | CWHeight
				     | CWBorderWidth)
	|| ((request->request_mode & CWBorderWidth)
	    && request->border_width > 0))
	return XtGeometryNo;

    allowed = *request;

    reconfigured = GetGeometry( (Widget)vw,
			        (rWidth ? request->width : vw->core.width),
			        (rHeight ? request->height : vw->core.height)
			      );

    child_changed_size = ((rWidth && child->core.width != request->width) ||
			  (rHeight && child->core.height != request->height));

    height_remaining = vw->core.height;
    if (rWidth && vw->core.width != request->width) {
	if (VW.allowhoriz && request->width > vw->core.width) {
	    /* horizontal scrollbar will be needed so possibly reduce height */
	    Widget bar; 
	    if ((bar = VW.horiz_bar) == (Widget)NULL)
		bar = CreateScrollbar( vw, True );
	    height_remaining -= bar->core.height + bar->core.border_width;
	    reconfigured = True;
	}
	else {
	    allowed.width = vw->core.width;
	}
    }
    if (rHeight && height_remaining != request->height) {
	if (VW.allowvert && request->height > height_remaining) {
	    /* vertical scrollbar will be needed, so possibly reduce width */
	    if (!VW.allowhoriz || request->width < vw->core.width) {
		Widget bar;
		if ((bar = VW.vert_bar) == (Widget)NULL)
		    bar = CreateScrollbar( vw, False );
		if (!rWidth) {
		    allowed.width = vw->core.width;
		    allowed.request_mode |= CWWidth;
		}
		if ( (int)allowed.width >
		     (int)(bar->core.width + bar->core.border_width) )
		    allowed.width -= bar->core.width + bar->core.border_width;
		else
		    allowed.width = 1;
		reconfigured = True;
	    }
	}
	else {
	    allowed.height = height_remaining;
	}
    }

    if (allowed.width != request->width || allowed.height != request->height) {
	*reply = allowed;
	result = XtGeometryAlmost;
    }
    else {
	if (rWidth)  child->core.width = request->width;
	if (rHeight) child->core.height = request->height;
	result = XtGeometryYes;
    }

    if (reconfigured || child_changed_size)
	Layout(vw);

/*
    if (child_changed_size)
	ComputeLayout ((Widget)vw, True, True);
  */  
    return result;
}
#endif


#if 0
static Boolean GetGeometry(w, width, height)
    Widget w;
    Dimension width, height;
{
    XtWidgetGeometry geometry, return_geom;
    XtGeometryResult result;

    if (width == w->core.width && height == w->core.height)
	return False;

    geometry.request_mode = CWWidth | CWHeight;
    geometry.width = width;
    geometry.height = height;

    if (XtIsRealized(w)) {
	if (((ViewportWidget)w)->viewport.allowhoriz && width > w->core.width)
	    geometry.width = w->core.width;
	if (((ViewportWidget)w)->viewport.allowvert && height > w->core.height)
	    geometry.height = w->core.height;
    } else {
	/* This is the Realize call; we'll inherit a w&h iff none currently */
	if (w->core.width != 0) {
	    if (w->core.height != 0) return False;
	    geometry.width = w->core.width;
	}
	if (w->core.height != 0) geometry.height = w->core.height;
    }

    result = XtMakeGeometryRequest(w, &geometry, &return_geom);
    if (result == XtGeometryAlmost)
	result = XtMakeGeometryRequest(w, &return_geom, NULL);

    return (result == XtGeometryYes);
}
#endif


static XtGeometryResult PreferredGeometry(w, intended, preferred)
     Widget            w;
     XtWidgetGeometry *intended;
     XtWidgetGeometry *preferred;
{
  register ViewportWidget vw = (ViewportWidget) w;
  XtWidgetGeometry reply_return;

  reply_return.width  = 0;
  reply_return.height = 0;

  _XawQueryGeometry (VW.child, &reply_return);
  
  /*
   * Preferred size
   */

  preferred->request_mode = 0;

  if (XtIsRealized((Widget)vw) || vw->core.width == 0) {
    preferred->width        = reply_return.width + 2 * VW.distance;
    preferred->request_mode |= CWWidth;
  } else

  if (XtIsRealized((Widget)vw) || vw->core.height == 0) {
    preferred->height       = reply_return.height + 2 * VW.distance;
    preferred->request_mode |= CWHeight;
  }

/*  preferred->request_mode = CWWidth | CWHeight; */
  
  if (intended
      && ((intended->request_mode & (CWWidth | CWHeight)) ==
	  (CWWidth | CWHeight))
      && intended->width == preferred->width
      && intended->height == preferred->height)
  {  
    return XtGeometryYes;
  }
  else if (preferred->width == w->core.width
	   && preferred->height == w->core.height)
  {
    return XtGeometryNo;
  }
  else
  {
    return XtGeometryAlmost;
  }

}


void
XawViewportSetLocation (Widget gw,
			double xoff, double yoff)
{
  ViewportWidget vw = (ViewportWidget) gw;
  register Widget child = VW.child;
  Position x, y;
  
  
  if (!XtIsSubclass (gw, viewportWidgetClass)) {
    String subs[1];
    Cardinal num_subs = 1;
    subs[0] = gw->core.name;
    XtAppWarningMsg(XtWidgetToApplicationContext(gw),
		    "SetLocation", "SetLocation","XawToolkitError",
"routine: XawViewportSetLocation\n\
 Widget '%s' is not subclass of viewportWidgetClass",
		    subs, &num_subs);
    return;
  }
  
  
  
  vw = (ViewportWidget) gw;
  child = VW.child;
  
  if (xoff > 1.0)			/* scroll to right */
    x = child->core.width;
  else if (xoff < 0.0)		/* if the offset is < 0.0 nothing */ 
    x = child->core.x;
  else
    x = (Position) (((float) child->core.width) * xoff);
  
  if (yoff > 1.0) 
    y = child->core.height;
  else if (yoff < 0.0)
    y = child->core.y;
  else
    y = (Position) (((float) child->core.height) * yoff);
  
  MoveChild (vw, -x, -y, True);
}

void
XawViewportSetCoordinates (Widget gw,
			   int x, int y)
{
    ViewportWidget vw = (ViewportWidget) gw;
    Widget child = VW.child;

    if (x > (int)child->core.width) 
      x = child->core.width;
    else if (x < 0)
      x = child->core.x;

    if (y > (int)child->core.height)
      y = child->core.height;
    else if (y < 0)
      y = child->core.y;

    MoveChild (vw, -x, -y, True);
}

