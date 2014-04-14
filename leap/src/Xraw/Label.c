/* $XConsortium: Label.c,v 1.93 91/10/16 21:34:35 eswu Exp $ */

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
 * Label.c - Label widget
 *
 */
#include <stdio.h>
#include <ctype.h>

#ifdef SYSV
#include <string.h>
#else
#include <strings.h>
#endif

#include <X11/IntrinsicP.h>
#include <X11/RectObjP.h>
#include <X11/StringDefs.h>
#include <X11/Xos.h>

#include "../Xmu/Converters.h"
#include "../Xmu/Drawing.h"

#include "XawInit.h"
#include "LabelP.h"
#include "MenuButtoP.h"

#include "XrawDebug.h"


#define MULTI_LINE_LABEL 32767

#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#ifdef CRAY
#define WORD64
#endif

/****************************************************************
 *
 * Full class record constant
 *
 ****************************************************************/

/* Private Data */

static void Def_font();


#define offset(field) XtOffsetOf(LabelRec, field)
static XtResource resources[] = {
    {
      XtNforeground, XtCForeground, XtRPixel, sizeof(Pixel),
      XtOffsetOf(LabelRec,simple.foreground), XtRString, 
      (caddr_t) XtDefaultForeground
    },
    {
      XtNfont,  XtCFont, XtRFontStruct, sizeof(XFontStruct *),
      offset(label.font),XtRCallProc, (XtPointer)Def_font
	/*XtRString, XtDefaultFont*/ /* vtr */
    },
    {
      XtNlabel,  XtCLabel, XtRString, sizeof(String),
      offset(label.label), XtRString, NULL
    },
    {
      XtNencoding, XtCEncoding, XtRUnsignedChar, sizeof(unsigned char),
      offset(label.encoding), XtRImmediate, (XtPointer)XawTextEncoding8bit
    },
    {
      XtNjustify, XtCJustify, XtRJustify, sizeof(XtJustify),
      offset(label.justify), XtRImmediate, (XtPointer)XtJustifyCenter
    },
    {
      XtNinternalWidth, XtCWidth, XtRDimension,  sizeof(Dimension),
      offset(label.internal_width), XtRImmediate, (XtPointer)4
    },
    {
      XtNinternalHeight, XtCHeight, XtRDimension, sizeof(Dimension),
      offset(label.internal_height), XtRImmediate, (XtPointer)2
    },
    {
      XtNleftBitmap, XtCLeftBitmap, XtRPixmap, sizeof(Pixmap),
      offset(label.left_bitmap), XtRImmediate, (XtPointer) None
    },
    {
      XtNbitmap, XtCPixmap, XtRPixmap, sizeof(Pixmap),
      offset(label.pixmap), XtRImmediate, (XtPointer)None
    },
    {
      XtNresize, XtCResize, XtRBoolean, sizeof(Boolean),
      offset(label.resize), XtRImmediate, (XtPointer)True
    },
    {
      XtNshadowWidth, XtCShadowWidth, XtRDimension, sizeof(Dimension),
      offset(simple.shadow_thickness), XtRImmediate, (XtPointer) 0
    },
    {
      XtNborderWidth, XtCBorderWidth, XtRDimension, sizeof(Dimension),
      XtOffsetOf(RectObjRec,rectangle.border_width), XtRImmediate,
      (XtPointer)0
    }
};

/* ARGSUSED */
static void Def_font(w, off, value)
     Widget    w;
     int       off;
     XrmValue *value;
{
  static XFontStruct *font;

  font = XLoadQueryFont(XtDisplay(w),"-*-helvetica-bold-r-*-*-14-*-*-*-*-*-*-*");

  if (font == NULL)
  {
    XrmValue source, dest;

    source.size = strlen(XtDefaultFont)+1;
    source.addr = XtDefaultFont;
    dest.size = sizeof(XFontStruct *);
    dest.addr = (caddr_t)&font;
    
  (void) XtConvertAndStore(w, XtRString, &source, XtRFont, &dest);

  }
    
  value->addr = (caddr_t)&font;
}

#undef offset

static void Initialize();
static void ClassPartInitialize();
static void Resize();
static void Redisplay();
static Boolean SetValues();
static void ClassInitialize();
static void Destroy();
static XtGeometryResult QueryGeometry();
static void LabelPosition();

LabelClassRec labelClassRec = {
  {
/* core_class fields */	
    /* superclass	  	*/	(WidgetClass) &simpleClassRec,
    /* class_name	  	*/	"Label",
    /* widget_size	  	*/	sizeof(LabelRec),
    /* class_initialize   	*/	ClassInitialize,
    /* class_part_initialize	*/	ClassPartInitialize,
    /* class_inited       	*/	FALSE,
    /* initialize	  	*/	Initialize,
    /* initialize_hook		*/	NULL,
    /* realize		  	*/	XtInheritRealize,
    /* actions		  	*/	NULL,
    /* num_actions	  	*/	0,
    /* resources	  	*/	resources,
    /* num_resources	  	*/	XtNumber(resources),
    /* xrm_class	  	*/	NULLQUARK,
    /* compress_motion	  	*/	TRUE,
    /* compress_exposure  	*/	TRUE,
    /* compress_enterleave	*/	TRUE,
    /* visible_interest	  	*/	FALSE,
    /* destroy		  	*/	Destroy,
    /* resize		  	*/	Resize,
    /* expose		  	*/	Redisplay,
    /* set_values	  	*/	SetValues,
    /* set_values_hook		*/	NULL,
    /* set_values_almost	*/	XtInheritSetValuesAlmost,
    /* get_values_hook		*/	NULL,
    /* accept_focus	 	*/	NULL,
    /* version			*/	XtVersion,
    /* callback_private   	*/	NULL,
    /* tm_table		   	*/	NULL,
    /* query_geometry		*/	QueryGeometry,
    /* display_accelerator	*/	XtInheritDisplayAccelerator,
    /* extension		*/	NULL
  },
    /* Simple class fields initialization */
  {
    /* change_sensitive         */      XtInheritChangeSensitive,
    /* display_rect		*/	XtInheritDisplayRectProc,
    /* extension                */      NULL
  },
    /* Label class fields initialization */
  {
    /* label_position 		*/	LabelPosition
  }
};
WidgetClass labelWidgetClass = (WidgetClass)&labelClassRec;
/****************************************************************
 *
 * Private Procedures
 *
 ****************************************************************/

static void ClassInitialize()
{
    XawInitializeWidgetSet();
    XtAddConverter( XtRString, XtRJustify, XmuCvtStringToJustify, 
		    (XtConvertArgList)NULL, 0 );
}

static void ClassPartInitialize(class)
    WidgetClass class;
{
  register LabelWidgetClass c = (LabelWidgetClass)class;

  if (c->label_class.label_position == XtInheritLabelPosition)
    c->label_class.label_position = LabelPosition;
}

#ifndef WORD64

#define TXT16 XChar2b

#else

#define TXT16 char

static XChar2b *buf2b;
static int buf2blen = 0;

_XawLabelWidth16(fs, str, n)
    XFontStruct *fs;
    char *str;
    int	n;
{
    int i;
    XChar2b *ptr;

    if (n > buf2blen) {
	buf2b = (XChar2b *)XtRealloc((char *)buf2b, n * sizeof(XChar2b));
	buf2blen = n;
    }
    for (ptr = buf2b, i = n; --i >= 0; ptr++) {
	ptr->byte1 = *str++;
	ptr->byte2 = *str++;
    }
    return XTextWidth16(fs, buf2b, n);
}

_XawLabelDraw16(dpy, d, gc, x, y, str, n)
    Display *dpy;
    Drawable d;
    GC gc;
    int x, y;
    char *str;
    int n;
{
    int i;
    XChar2b *ptr;

    if (n > buf2blen) {
	buf2b = (XChar2b *)XtRealloc((char *)buf2b, n * sizeof(XChar2b));
	buf2blen = n;
    }
    for (ptr = buf2b, i = n; --i >= 0; ptr++) {
	ptr->byte1 = *str++;
	ptr->byte2 = *str++;
    }
    XDrawString16(dpy, d, gc, x, y, buf2b, n);
}

#define XTextWidth16 _XawLabelWidth16
#define XDrawString16 _XawLabelDraw16

#endif /* WORD64 */

/*
 * Calculate width and height of displayed text in pixels
 */

static void SetTextWidthAndHeight(w)
    LabelWidget w;
{
  register LabelWidget lw = w;
  register XFontStruct *fs = lw->label.font;
  char *nl;

  if (lw->label.pixmap != None) {
    Window root;
    int x, y;
    unsigned int width, height, bw, depth;
    if (XGetGeometry(XtDisplay(lw), lw->label.pixmap, &root, &x, &y,
		     &width, &height, &bw, &depth)) {
      lw->label.label_height = height;
      lw->label.label_width = width;
      lw->label.label_len = depth;
      return;
    }
  }
  
  lw->label.label_height   = fs->max_bounds.ascent + fs->max_bounds.descent;
  if (lw->label.label == NULL) {
    lw->label.label_len = 0;
    lw->label.label_width = 0;
  }
  else if ((nl = index(lw->label.label, '\n')) != NULL) {
    char *label;
    lw->label.label_len = MULTI_LINE_LABEL;
    lw->label.label_width = 0;
    for (label = lw->label.label; nl != NULL; nl = index(label, '\n')) {
      int width;
      
      if (lw->label.encoding)
	width = XTextWidth16(fs, (TXT16*)label, (int)(nl - label)/2);
      else
	width = XTextWidth(fs, label, (int)(nl - label));
      if (width > (int)lw->label.label_width)
	lw->label.label_width = width;
      label = nl + 1;
      if (*label)
	lw->label.label_height +=
	  fs->max_bounds.ascent + fs->max_bounds.descent;
    }
    if (*label) {
      int width;
      
      if (lw->label.encoding)
	width = XTextWidth16(fs, (TXT16*)label, (int)strlen(label)/2);
      else
	width = XTextWidth(fs, label, strlen(label));

      if (width > (int) lw->label.label_width)
	lw->label.label_width = width;
    }
  } else {
    lw->label.label_len = strlen(lw->label.label);
    if (lw->label.encoding)
      lw->label.label_width =
	XTextWidth16(fs, (TXT16*)lw->label.label,
		     (int) lw->label.label_len/2);
    else
      lw->label.label_width =
	XTextWidth(fs, lw->label.label, (int) lw->label.label_len);
  }

/*
  lw->label.label_height += 2 * SIMPLE_MARGIN(w);
  lw->label.label_width  += 2 * SIMPLE_MARGIN(w);
*/
}

static void GetnormalGC(lw)
    LabelWidget lw;
{
    XGCValues	values;

    values.foreground	= lw->simple.foreground;
    values.background	= lw->core.background_pixel;
    values.font		= lw->label.font->fid;
    values.graphics_exposures = False;

    lw->label.normal_GC = XtGetGC(
	(Widget)lw,
	(unsigned) GCForeground | GCBackground | GCFont | GCGraphicsExposures,
	&values);
}

static void GetgrayGC(lw)
    LabelWidget lw;
{
    XGCValues	values;

    values.foreground = lw->simple.foreground;
    values.background = lw->core.background_pixel;
    values.font	      = lw->label.font->fid;
    values.fill_style = FillTiled;
    values.tile       = XmuCreateStippledPixmap(XtScreen((Widget)lw),
						lw->simple.foreground, 
						lw->core.background_pixel,
						lw->core.depth);
    values.graphics_exposures = False;

    lw->label.stipple = values.tile;
    lw->label.gray_GC = XtGetGC((Widget)lw, 
				(unsigned) GCForeground | GCBackground |
					   GCFont | GCTile | GCFillStyle |
					   GCGraphicsExposures,
				&values);
}

static void compute_bitmap_offsets (lw)
    LabelWidget lw;
{
    /*
     * bitmap will be eventually be displayed at 
     * (internal_width, internal_height + lbm_y)
     */
    if (lw->label.lbm_height != 0) {
	lw->label.lbm_y = (lw->core.height -
			  (SIMPLE_MARGIN(lw) * 2 + 
			   lw->label.internal_height * 2 +
			   lw->label.lbm_height)) / 2;
    } else {
	lw->label.lbm_y = 0;
    }
}


static void set_bitmap_info (lw)
    LabelWidget lw;
{
    Window root;
    int x, y;
    unsigned int bw, depth;

    if (!(lw->label.left_bitmap &&
	  XGetGeometry (XtDisplay((Widget)lw), lw->label.left_bitmap, 
			&root, &x, &y,
			&lw->label.lbm_width, &lw->label.lbm_height,
			&bw, &depth))) {
	lw->label.lbm_width = lw->label.lbm_height = 0;
    }
    compute_bitmap_offsets (lw);
}



/* ARGSUSED */
static void Initialize(request, new, args, num_args)
    Widget request, new;
    ArgList args;
    Cardinal *num_args;
{
    LabelWidget lw = (LabelWidget) new;

    if (lw->label.label == NULL) 
        lw->label.label = XtNewString(lw->core.name);
    else 
        lw->label.label = XtNewString(lw->label.label);

    GetnormalGC(lw);
    GetgrayGC(lw);

    if (lw->core.height == 0 || lw->core.width == 0){
      XtWidgetGeometry preferred;
      
      (void)(*XtClass(new)->core_class.query_geometry)(new, NULL, &preferred);
      
      if (lw->core.width  == 0) lw->core.width = preferred.width;
      if (lw->core.height == 0) lw->core.height = preferred.height;
    }
    else
      SetTextWidthAndHeight(lw);

    set_bitmap_info(lw);
    
    lw->label.label_x = lw->label.label_y = 0;

    (*XtClass(new)->core_class.resize) (new);

} /* Initialize */

/*
 * Repaint the widget window
 */

/* ARGSUSED */
static void Redisplay(gw, event, region)
    Widget gw;
    XEvent *event;
    Region region;
{
  LabelWidget w = (LabelWidget) gw;
  CommandWidget cw = (CommandWidget) gw;  
  XRectangle rectangle[1];
  GC gc;

  if (!XtIsRealized(gw))
    return;
  
  if (w->simple.shadow_thickness && !XtIsSubclass(gw, menuButtonWidgetClass))
    (*simpleWidgetClass->core_class.expose) (gw, event, region);
  
  /*
   * now we'll see if we need to draw the rest of the label
   */
  if (region != NULL) {
    int x = w->label.label_x + SIMPLE_MARGIN(gw);
    unsigned int width = w->label.label_width;
    if (w->label.lbm_width) {
      if (w->label.label_x > (x = w->label.internal_width))
	width += w->label.label_x - x;
    }
    if (XRectInRegion(region, x, w->label.label_y,
		      width, w->label.label_height) == RectangleOut){
      return;
    }
  }
  
  gc = XtIsSensitive(gw) ? w->label.normal_GC : w->label.gray_GC;
  
  (*XtLabelClass(w)->simple_class.display_rect)(gw, rectangle);
  XSetClipRectangles (XtDisplay(gw), gc, 0, 0, rectangle, 1, YSorted);
  
  if (w->label.pixmap == None) {
    int len = w->label.label_len;
    char *label = w->label.label;
    Position y = w->label.label_y + w->label.font->max_bounds.ascent;

    y += SIMPLE_MARGIN(gw);
    
    /* 
     * display left bitmap 
     */
    if (w->label.left_bitmap && w->label.lbm_width != 0) 
    {
      Window root;
      int x, y;
      unsigned int width, height, bw, depth;
      
      (void) XGetGeometry (XtDisplay(gw), w->label.left_bitmap, &root,
			   &x, &y, &width, &height, &bw, &depth);
      
      if (depth == 1) 
      {
	GC bitmap_gc;
	int bitmap_x;
	int bitmap_y;
	
	bitmap_x = w->label.internal_width + SIMPLE_MARGIN(w);
	bitmap_y = w->label.internal_height + SIMPLE_MARGIN(w)+ w->label.lbm_y;
	
	if (XtIsSubclass(gw, commandWidgetClass) && cw->command.set) {
	  bitmap_gc = cw->command.shift_GC;
	  bitmap_x += 1;
	  bitmap_y += 1;
	} else {
	  bitmap_gc = gc;
	}
	
	XCopyPlane (XtDisplay(gw),
		    w->label.left_bitmap,
		    XtWindow(gw),
		    bitmap_gc,
		    0, 0, w->label.lbm_width, w->label.lbm_height,
		    bitmap_x, bitmap_y,
		    (unsigned long) 1L);
      } else {
	int bitmap_x;
	int bitmap_y;

	bitmap_x = w->label.internal_width + SIMPLE_MARGIN(w);
	bitmap_y = w->label.internal_height + SIMPLE_MARGIN(w)+ w->label.lbm_y;
	
	if (XtIsSubclass(gw, commandWidgetClass) && cw->command.set) {
	  bitmap_x += 1;
	  bitmap_y += 1;
	}
	
	XCopyArea(XtDisplay(gw),
		  w->label.left_bitmap,
		  XtWindow(gw),
		  gc,
		  0,
		  0,
		  (unsigned int)w->label.lbm_width,
		  (unsigned int)w->label.lbm_height,
		  bitmap_x,
		  bitmap_y);
      }
    }
    
    if (len == MULTI_LINE_LABEL) {
      char *nl;
      while ((nl = index(label, '\n')) != NULL) {
	if (w->label.encoding)
	  XDrawString16(XtDisplay(gw), XtWindow(gw), gc,
			w->label.label_x, 
			y - SIMPLE_MARGIN(gw),
			(TXT16*)label, (int)(nl - label)/2);
	else
	  XDrawString(XtDisplay(gw), XtWindow(gw), gc,
		      w->label.label_x, 
		      y - SIMPLE_MARGIN(gw), 
		      label, (int)(nl - label));
	y += w->label.font->max_bounds.ascent + 
	  w->label.font->max_bounds.descent;
	label = nl + 1;
      }
      len = strlen(label);
    }
    if (len) {
      if (w->label.encoding)
	XDrawString16(XtDisplay(gw), XtWindow(gw), gc,
		      w->label.label_x /*+ SIMPLE_MARGIN(gw)*/, 
		      y - SIMPLE_MARGIN(gw), 
		      (TXT16*)label, len/2);
      else
	XDrawString(XtDisplay(gw), XtWindow(gw), gc,
		    w->label.label_x /*+ SIMPLE_MARGIN(gw)*/, 
		    y - SIMPLE_MARGIN(gw), label, len);
    }
  } else {
    
    /* Vladimir Romanovski start adding*/
    
    Window root;
    int x, y;
    unsigned int width, height, bw, depth;
    
    (void)XGetGeometry(XtDisplay(gw), w->label.pixmap, &root,
		       &x, &y, &width, &height, &bw, &depth);
    
    if (depth == 1) 
      XCopyPlane (XtDisplay(gw), w->label.pixmap, XtWindow(gw), gc,
		  0, 0, w->label.label_width, w->label.label_height,
		  (int) w->label.label_x,
		  (int) w->label.label_y,
		  (unsigned long) 1L);
    else
      /* Vladimir Romanovski stop adding */
      
      XCopyArea(XtDisplay(gw), w->label.pixmap, XtWindow(gw), gc,
		0, 0, w->label.label_width, w->label.label_height,
		w->label.label_x, w->label.label_y);
  }
  
  XSetClipMask (XtDisplay(w), gc, None);
}

static void Resize(w)
    Widget w;
{
    LabelWidget lw = (LabelWidget)w;

    (*(XtLabelClass(lw)->label_class.label_position)) (w);
    compute_bitmap_offsets (lw);
}

/*
 * Set specified arguments into widget
 */

#define PIXMAP 0
#define WIDTH 1
#define HEIGHT 2
#define NUM_CHECKS 3

static Boolean SetValues(current, request, new, args, num_args)
    Widget current, request, new;
    ArgList args;
    Cardinal *num_args;
{
    register LabelWidget curlw = (LabelWidget) current;
    register LabelWidget reqlw = (LabelWidget) request;
    register LabelWidget newlw = (LabelWidget) new;
    int i;
    Boolean was_resized = False, redisplay = False, checks[NUM_CHECKS];

#define NE(field) (curlw->label.field !=  newlw->label.field)
    
    for (i = 0; i < NUM_CHECKS; i++)
	checks[i] = FALSE;

    for (i = 0; i < *num_args; i++) {
	if (streq(XtNbitmap, args[i].name))
	    checks[PIXMAP] = TRUE;
	if (streq(XtNwidth, args[i].name))
	    checks[WIDTH] = TRUE;
	if (streq(XtNheight, args[i].name))
	    checks[HEIGHT] = TRUE;
    }

    if (newlw->label.label == NULL) {
	newlw->label.label = newlw->core.name;
    }

    /*
     * resize on bitmap change
     */

    if (NE(left_bitmap) || NE(encoding)) 
      was_resized = True;

    if (NE(label))
    {
      if (streq(curlw->label.label, newlw->label.label))
      {
	/* It has not sense to redraw the same string */
	newlw->label.label = curlw->label.label;
      }
      else
      {
	if (curlw->label.label != curlw->core.name)
	  XtFree( (char *)curlw->label.label );
	
	if (newlw->label.label != newlw->core.name) {
	  newlw->label.label = XtNewString (newlw->label.label);
	}
	was_resized = True;
      }
    }

    if (was_resized || NE(font) || NE(justify) || checks[PIXMAP]) {

	SetTextWidthAndHeight(newlw);
	was_resized = True;
    }

    /*
     * recalculate the window size if something has changed.
     */
    if (newlw->label.resize && was_resized) {
      XtWidgetGeometry preferred;

      (void) (*XtClass(new)->core_class.query_geometry)(new, NULL, &preferred);
      
      if ((curlw->core.width == reqlw->core.width) && !checks[WIDTH])
	newlw->core.width = preferred.width;

      if ((curlw->core.height == reqlw->core.height) && !checks[HEIGHT]) 
	newlw->core.height = preferred.height;
    }

    if (curlw->simple.foreground != newlw->simple.foreground
	|| curlw->core.background_pixel != newlw->core.background_pixel
	|| curlw->label.font->fid != newlw->label.font->fid) {

	XtReleaseGC(new, curlw->label.normal_GC);
	XtReleaseGC(new, curlw->label.gray_GC);
	XmuReleaseStippledPixmap( XtScreen(current), curlw->label.stipple );
	GetnormalGC(newlw);
	GetgrayGC(newlw);
	redisplay = True;
    }

    return was_resized || redisplay ||
      XtIsSensitive(current) != XtIsSensitive(new);
}

static void Destroy(w)
    Widget w;
{
    LabelWidget lw = (LabelWidget)w;

    XtFree( lw->label.label );
    XtReleaseGC( w, lw->label.normal_GC );
    XtReleaseGC( w, lw->label.gray_GC);
    XmuReleaseStippledPixmap( XtScreen(w), lw->label.stipple );
}


static XtGeometryResult QueryGeometry(w, intended, preferred)
    Widget w;
    XtWidgetGeometry *intended, *preferred;
{
    register LabelWidget lw = (LabelWidget)w;


    preferred->request_mode = CWWidth | CWHeight;

    SetTextWidthAndHeight(lw);
    
    /*
     *  get geometry of left bitmap
     */
    {
      Window root;
      int x, y;
      unsigned int bw, depth;
      
      if (!(lw->label.left_bitmap &&
	    XGetGeometry (XtDisplay(lw), lw->label.left_bitmap, &root,
			  &x, &y, &lw->label.lbm_width, &lw->label.lbm_height,
			  &bw, &depth)))
      {
	lw->label.lbm_width = lw->label.lbm_height = 0;
      }
      
    }
    
    preferred->height = (MAX(lw->label.label_height,lw->label.lbm_height) +
			 2 * lw->label.internal_height +
			 2 * SIMPLE_MARGIN(lw) );
    
    
    /*
     * here it is used left bitmap width which was defined above
     */
    preferred->width = (lw->label.label_width +
			2 * lw->label.internal_width +
			2 * SIMPLE_MARGIN(lw) + 
			LEFT_OFFSET(lw));
#define WIDTH_HEIGHT (CWWidth | CWHeight)
    
    if (intended
	&& ((intended->request_mode & WIDTH_HEIGHT) == WIDTH_HEIGHT)
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


static void LabelPosition(w)
     Widget w;
{
  register LabelWidget lw = (LabelWidget)w;
  Position width, height;
  Position newPos;
  Position leftedge;
  

  if (!XtIsSubclass(w, labelWidgetClass))
    return ;

  width    = lw->core.width;
  height   = lw->core.height;
  leftedge = lw->label.internal_width + SIMPLE_MARGIN(lw) + LEFT_OFFSET(lw);


  switch (lw->label.justify) {

  case XtJustifyLeft   :	 newPos = leftedge;

    break;
  case XtJustifyRight  :	 newPos = width - (lw->label.label_width +
						   lw->label.internal_width +
						   SIMPLE_MARGIN(lw));
    break;
  case XtJustifyCenter :
  default              :         newPos = (width - lw->label.label_width) / 2;
    break;
  }

  if (newPos < leftedge)
    lw->label.label_x = leftedge;
  else
    lw->label.label_x = newPos;

  lw->label.label_y = (height - lw->label.label_height) / 2;
  

  return;
  
}

