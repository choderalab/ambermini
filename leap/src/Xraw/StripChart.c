/* $XConsortium: StripChart.c,v 1.20 91/05/24 17:20:42 converse Exp $ */

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

#include <stdio.h>
#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#include "XawInit.h"
#include "3d.h"
#include "StripCharP.h"

#include "XrawDebug.h"


#define MS_PER_SEC 1000

/* Private Data */

#define offset(field) XtOffsetOf(StripChartRec, field)

static XtResource resources[] = {
    {XtNwidth, XtCWidth, XtRDimension, sizeof(Dimension),
	offset(core.width), XtRImmediate, (XtPointer) 120},
    {XtNheight, XtCHeight, XtRDimension, sizeof(Dimension),
	offset(core.height), XtRImmediate, (XtPointer) 120},
    {XtNupdate, XtCInterval, XtRInt, sizeof(int),
        offset(strip_chart.update), XtRImmediate, (XtPointer) 10},
    {XtNminScale, XtCScale, XtRInt, sizeof(int),
        offset(strip_chart.min_scale), XtRImmediate, (XtPointer) 1},
    {XtNforeground, XtCForeground, XtRPixel, sizeof(Pixel),
        offset(strip_chart.fgpixel), XtRString, XtDefaultForeground},
    {XtNhighlight, XtCForeground, XtRPixel, sizeof(Pixel),
        offset(strip_chart.hipixel), XtRString, XtDefaultForeground},
    {XtNgetValue, XtCCallback, XtRCallback, sizeof(XtPointer),
        offset(strip_chart.get_value), XtRImmediate, (XtPointer) NULL},
    {XtNjumpScroll, XtCJumpScroll, XtRInt, sizeof(int),
        offset(strip_chart.jump_val), XtRImmediate, (XtPointer) DEFAULT_JUMP},
};

#undef offset

static void Initialize(), Destroy(), Redisplay(), MoveChart(), Resize();
static Boolean SetValues();
static int repaint_window();

#define SuperClass ((SimpleWidgetClass)&simpleClassRec)

StripChartClassRec stripChartClassRec = {
    { /* core fields */
    /* superclass		*/	(WidgetClass) SuperClass,
    /* class_name		*/	"StripChart",
    /* size			*/	sizeof(StripChartRec),
    /* class_initialize		*/	XawInitializeWidgetSet,
    /* class_part_initialize	*/	NULL,
    /* class_inited		*/	FALSE,
    /* initialize		*/	Initialize,
    /* initialize_hook		*/	NULL,
    /* realize			*/	XtInheritRealize,
    /* actions			*/	NULL,
    /* num_actions		*/	0,
    /* resources		*/	resources,
    /* num_resources		*/	XtNumber(resources),
    /* xrm_class		*/	NULLQUARK,
    /* compress_motion		*/	TRUE,
    /* compress_exposure	*/	XtExposeCompressMultiple |
	                                XtExposeGraphicsExposeMerged,
    /* compress_enterleave	*/	TRUE,
    /* visible_interest		*/	FALSE,
    /* destroy			*/	Destroy,
    /* resize			*/	Resize,
    /* expose			*/	Redisplay,
    /* set_values		*/	SetValues,
    /* set_values_hook		*/	NULL,
    /* set_values_almost	*/	NULL,
    /* get_values_hook		*/	NULL,
    /* accept_focus		*/	NULL,
    /* version			*/	XtVersion,
    /* callback_private		*/	NULL,
    /* tm_table			*/	NULL,
    /* query_geometry		*/	XtInheritQueryGeometry,
    /* display_accelerator	*/	XtInheritDisplayAccelerator,
    /* extension		*/	NULL
    },
    { /* Simple fields */
    /* change_sensitive         */      XtInheritChangeSensitive,
    /* display_rect             */      XtInheritDisplayRectProc,
    /* extension                */      NULL
    },
    { /* Stripchart fields */
    /* ignore			*/	0
    }
};

WidgetClass stripChartWidgetClass = (WidgetClass) &stripChartClassRec;

/****************************************************************
 *
 * Private Procedures
 *
 ****************************************************************/

static void draw_it();

/*	Function Name: CreateGC
 *	Description: Creates the GC's
 *	Arguments: w - the strip chart widget.
 *                 which - which GC's to create.
 *	Returns: none
 */

static void
CreateGC(w, which)
StripChartWidget w;
unsigned int which;
{
  XGCValues	myXGCV;

  if (which & FOREGROUND) {
    myXGCV.foreground = w->strip_chart.fgpixel;
    w->strip_chart.fgGC = XtGetGC((Widget) w, GCForeground, &myXGCV);
  }

  if (which & HIGHLIGHT) {
    myXGCV.foreground = w->strip_chart.hipixel;
    w->strip_chart.hiGC = XtGetGC((Widget) w, GCForeground, &myXGCV);
  }
}

/*	Function Name: DestroyGC
 *	Description: Destroys the GC's
 *	Arguments: w - the strip chart widget.
 *                 which - which GC's to destroy.
 *	Returns: none
 */

static void
DestroyGC(w, which)
StripChartWidget w;
unsigned int which;
{
  if (which & FOREGROUND) 
    XtReleaseGC((Widget) w, w->strip_chart.fgGC);

  if (which & HIGHLIGHT) 
    XtReleaseGC((Widget) w, w->strip_chart.hiGC);
}

/* ARGSUSED */
static void Initialize (greq, gnew)
    Widget greq, gnew;
{
    StripChartWidget w = (StripChartWidget)gnew;

    if (w->strip_chart.update > 0)
        w->strip_chart.interval_id = XtAppAddTimeOut( 
					XtWidgetToApplicationContext(gnew),
					w->strip_chart.update * MS_PER_SEC, 
					draw_it, (XtPointer) gnew);
    CreateGC(w, (unsigned int) ALL_GCS);

    w->strip_chart.scale = w->strip_chart.min_scale;
    w->strip_chart.interval = 0;
    w->strip_chart.max_value = 0.0;
    w->strip_chart.points = NULL;
    Resize(w);
}
 
static void Destroy (gw)
     Widget gw;
{
     StripChartWidget w = (StripChartWidget)gw;

     if (w->strip_chart.update > 0)
         XtRemoveTimeOut (w->strip_chart.interval_id);
     if (w->strip_chart.points)
	 XtFree((char *) w->strip_chart.points);
     DestroyGC(w, (unsigned int) ALL_GCS);
}

/*
 * NOTE: This function really needs to recieve graphics exposure 
 *       events, but since this is not easily supported until R4 I am
 *       going to hold off until then.
 */

/* ARGSUSED */
static void Redisplay(gw, event, region)
     Widget gw;
     XEvent *event;
     Region region;
{
    StripChartWidget w = (StripChartWidget) gw;
    Dimension	s = w->simple.shadow_thickness;

    (*SuperClass->core_class.expose) (gw, event, region);

    if (event->type == GraphicsExpose)
	(void) repaint_window (w, event->xgraphicsexpose.x,
			       event->xgraphicsexpose.width - 2 * s);
    else
	(void) repaint_window (w, event->xexpose.x,
			       event->xexpose.width - 2 * s);
}

/* ARGSUSED */
static void 
draw_it(client_data, id)
XtPointer client_data;
XtIntervalId *id;		/* unused */
{
    StripChartWidget w = (StripChartWidget)client_data;
    Dimension        s = w->simple.shadow_thickness;
    double      value;
   
    if (w->strip_chart.update > 0)
	w->strip_chart.interval_id =
	    XtAppAddTimeOut(XtWidgetToApplicationContext( (Widget) w),
		       w->strip_chart.update * MS_PER_SEC,draw_it,client_data);

    if (w->strip_chart.interval >= (int)(w->core.width - 2 * s))
	MoveChart( (StripChartWidget) w, TRUE);

    /* Get the value, stash the point and draw corresponding line. */

    if (w->strip_chart.get_value == NULL)
	return;

    XtCallCallbacks( (Widget)w, XtNgetValue, (XtPointer)&value );

    /* 
     * Keep w->strip_chart.max_value up to date, and if this data 
     * point is off the graph, change the scale to make it fit. 
     */
   
    if (value > w->strip_chart.max_value) {
	w->strip_chart.max_value = value;
	if (w->strip_chart.max_value > w->strip_chart.scale) {
	    XClearWindow( XtDisplay ((Widget) w), XtWindow ((Widget) w));
	    w->strip_chart.interval = 
		repaint_window(w, 0, (int) (w->core.width - 2 * s));
	    XawDrawFrame ((Widget)w,
			  w->simple.highlight_thickness,
			  w->simple.highlight_thickness,
		(Dimension)(w->core.width - 2*w->simple.highlight_thickness),
		(Dimension)(w->core.height - 2*w->simple.highlight_thickness),
			  XawRAISED,
			  w->simple.shadow_thickness,
			  w->simple.top_shadow_GC,
			  w->simple.bottom_shadow_GC);
        }
    }

    w->strip_chart.valuedata[w->strip_chart.interval] = value;
    if (XtIsRealized((Widget) w)) {
	int y = (int) ((w->core.height - 2 * s)
			- (int)((w->core.height - 2 * s) * value) 
			/ w->strip_chart.scale);

        XFillRectangle(XtDisplay((Widget) w), XtWindow((Widget) w), w->strip_chart.fgGC,
		      w->strip_chart.interval + s, y + s, 
		      (unsigned int) 1, (w->core.height - 2 * s) - y);
        /*
	 * Fill in the graph lines we just painted over.
	 */

        if (w->strip_chart.points != NULL) {
	    w->strip_chart.points[0].x = w->strip_chart.interval + s;
	    XDrawPoints(XtDisplay(w), XtWindow(w), w->strip_chart.hiGC,
		       w->strip_chart.points, w->strip_chart.scale - 1,
		       CoordModePrevious);
	}

	XFlush(XtDisplay(w));		    /* Flush output buffers */
    }
    w->strip_chart.interval++;		    /* Next point */
} /* draw_it */

/* Blts data according to current size, then redraws the stripChart window.
 * Next represents the number of valid points in data.  Returns the (possibly)
 * adjusted value of next.  If next is 0, this routine draws an empty window
 * (scale - 1 lines for graph).  If next is less than the current window width,
 * the returned value is identical to the initial value of next and data is
 * unchanged.  Otherwise keeps half a window's worth of data.  If data is
 * changed, then w->strip_chart.max_value is updated to reflect the
 * largest data point.
 */

static int 
repaint_window(w, left, width)
StripChartWidget w;
int left, width;
{
    Dimension s = w->simple.shadow_thickness;

    register int i, j;
    register int next = w->strip_chart.interval;
    int scale = w->strip_chart.scale;
    int scalewidth = 0;

    /* Compute the minimum scale required to graph the data, but don't go
       lower than min_scale. */
    if (w->strip_chart.interval != 0 || scale <= (int)w->strip_chart.max_value)
      scale = ((int) (w->strip_chart.max_value)) + 1;
    if (scale < w->strip_chart.min_scale)
      scale = w->strip_chart.min_scale;

    if (scale != w->strip_chart.scale) {
      w->strip_chart.scale = scale;
      left = 0;
      width = next;
      scalewidth = w->core.width - 2 * s;

      Resize(w);

      if (XtIsRealized ((Widget) w)) {
	XClearWindow (XtDisplay ((Widget) w), XtWindow ((Widget) w));
	    XawDrawFrame ((Widget)w,
			  w->simple.highlight_thickness,
			  w->simple.highlight_thickness,
		(Dimension)(w->core.width - 2*w->simple.highlight_thickness),
		(Dimension)(w->core.height - 2*w->simple.highlight_thickness),
			  XawRAISED,
			  w->simple.shadow_thickness,
			  w->simple.top_shadow_GC,
			  w->simple.bottom_shadow_GC);
      }

    }

    if (XtIsRealized((Widget)w)) {
	Display *dpy = XtDisplay((Widget) w);
	Window win = XtWindow((Widget) w);

	width += left - 1;
	if (!scalewidth) scalewidth = width;

	if (next < (++width - s)) width = next + s;

	/* Draw data point lines. */
	for (i = left; i < width; i++) {
	    int y = (int) ((w->core.height - 2 * s)
			- (int)((w->core.height - 2 * s)
			* w->strip_chart.valuedata[i]) / w->strip_chart.scale);

	    XFillRectangle(dpy, win, w->strip_chart.fgGC, 
			   i, y, (unsigned int) 1,
			   (unsigned int) (w->core.height -2 * s - y));
	}

	/* Draw graph reference lines */
	for (i = 1; i < w->strip_chart.scale; i++) {
	    j = i * ((int)(w->core.height - 2 * s) / w->strip_chart.scale);
	    XDrawLine(dpy, win, w->strip_chart.hiGC, 
			left + s, j + s, scalewidth + s, j + s);
	}
    }
    return(next);
}

/*	Function Name: MoveChart
 *	Description: moves the chart over when it would run off the end.
 *	Arguments: w - the load widget.
 *                 blit - blit the bits? (TRUE/FALSE).
 *	Returns: none.
 */

static void
MoveChart(w, blit)
StripChartWidget w;
Boolean blit;
{
    Dimension s = w->simple.shadow_thickness;
    double old_max;
    int left, i, j;
    register int next = w->strip_chart.interval;

    if (!XtIsRealized((Widget) w)) return;

    if (w->strip_chart.jump_val == DEFAULT_JUMP)
        j = (w->core.width - 2 * s) / 2; /* Half the window width. */
    else {
        j = (w->core.width - 2 * s) - w->strip_chart.jump_val;
	if (j < 0) j = 0;
    }

    bcopy((char *)(w->strip_chart.valuedata + next - j),
	  (char *)(w->strip_chart.valuedata), j * sizeof(double));
    next = w->strip_chart.interval = j;
	
    /*
     * Since we just lost some data, recompute the 
     * w->strip_chart.max_value. 
     */

    old_max = w->strip_chart.max_value;
    w->strip_chart.max_value = 0.0;
    for (i = 0; i < next; i++) {
      if (w->strip_chart.valuedata[i] > w->strip_chart.max_value) 
	w->strip_chart.max_value = w->strip_chart.valuedata[i];
    }

    if (!blit) return;		/* we are done... */

    if ( ((int) old_max) != ( (int) w->strip_chart.max_value) ) {
      XClearWindow(XtDisplay(w), XtWindow(w));
      repaint_window(w, 0, (int) (w->core.width - 2 * s));
      return;
    }

    XCopyArea(XtDisplay((Widget)w), XtWindow((Widget)w), XtWindow((Widget)w),
		w->strip_chart.hiGC, (int) w->core.width - j + s, s,
		(unsigned int) j, (unsigned int) (w->core.height - 2 * s),
		s, s);

    XClearArea(XtDisplay((Widget)w), XtWindow((Widget)w), 
		(int) j + s, s, 
		(unsigned int) w->core.width - j, 
		(unsigned int) (w->core.height - 2 * s),
		FALSE);

    /* Draw graph reference lines */
    left = j;
    for (i = 1; i < w->strip_chart.scale; i++) {
      j = i * ((int)(w->core.height - 2 * s) / w->strip_chart.scale);
      XDrawLine(XtDisplay((Widget) w), XtWindow( (Widget) w),
		w->strip_chart.hiGC, 
		left, j + s, (int)(w->core.width - s), j + s);
    }
    return;
}

/* ARGSUSED */
static Boolean SetValues (current, request, new, args, num_args)
    Widget current, request, new;
    ArgList args;
    Cardinal *num_args;
{
    StripChartWidget old = (StripChartWidget)current;
    StripChartWidget w = (StripChartWidget)new;
    Boolean ret_val = FALSE;
    unsigned int new_gc = NO_GCS;

    if (w->strip_chart.update != old->strip_chart.update) {
	if (old->strip_chart.update > 0)
	    XtRemoveTimeOut (old->strip_chart.interval_id);
	if (w->strip_chart.update > 0)
	    w->strip_chart.interval_id =
		XtAppAddTimeOut(XtWidgetToApplicationContext(new),
				w->strip_chart.update * MS_PER_SEC,
				draw_it, (XtPointer)w);
    }

    if ( w->strip_chart.min_scale > (int) ((w->strip_chart.max_value) + 1) )
      ret_val = TRUE;
     
    if ( w->strip_chart.fgpixel != old->strip_chart.fgpixel ) {
      new_gc |= FOREGROUND;
      ret_val = True;
    }
    
    if ( w->strip_chart.hipixel != old->strip_chart.hipixel ) {
      new_gc |= HIGHLIGHT;
      ret_val = True;
    }
    
    DestroyGC(old, new_gc);
    CreateGC(w, new_gc);

    return( ret_val );
}

/*	Function Name: Resize
 *	Description: Sets up the polypoint that will be used to draw in
 *                   the graph lines.
 *	Arguments: w - the StripChart widget.
 *	Returns: none.
 */

#define HEIGHT ( (unsigned int) w->core.height)

static void
Resize(w)
StripChartWidget w;
{
    Dimension s = w->simple.shadow_thickness;
    XPoint * points;
    Cardinal size;
    int i;

    if (w->strip_chart.scale <= 1) { /* no scale lines. */
	XtFree ((char *) w->strip_chart.points);
	w->strip_chart.points = NULL;
	return;
    }
    
    size = sizeof(XPoint) * (w->strip_chart.scale - 1);

    points = (XPoint *) XtRealloc( (XtPointer) w->strip_chart.points, size);
    w->strip_chart.points = points;

    /* Draw graph reference lines into clip mask */

    for (i = 1; i < w->strip_chart.scale; i++) {
	points[i - 1].x = 0;
	points[i - 1].y =  (short)(((double)w->core.height - 2.0 * (double) s) / (double) w->strip_chart.scale);
    }
    points[0].y += (short)s;
}
