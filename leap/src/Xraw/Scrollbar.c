/* $XConsortium: Scrollbar.c,v 1.69 91/05/04 23:07:32 keith Exp $ */

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

/* ScrollBar.c */
/* created by weissman, Mon Jul  7 13:20:03 1986 */
/* converted by swick, Thu Aug 27 1987 */

#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#include <stdint.h>

#include "../Xmu/Drawing.h"

#include "XawInit.h"
#include "3d.h"
#include "ScrollbarP.h"
#include "SimpleP.h"
#include "ContainerP.h"
#include "Arrow.h"
#include "Text.h"


#ifdef EBUG_XRAW_MALLOC
#include <dbmalloc/malloc.h>

#include "XrawDebug.h"

#endif


/***********************************************************************
 *
 *   Private definitions
 *
 ***********************************************************************/

static float floatZero = 0.0;

static void Def_pixel();
static void Def_highlight_thickness();

#define Offset(field) XtOffsetOf(ScrollbarRec, field)

static XtResource resources[] = {
  {
    XtNincrement, XtCIncrement, XtRDimension, sizeof(Dimension),
    Offset(scrollbar.increment), XtRImmediate, (XtPointer) 10
  },
  {
    XtNpageIncrement, XtCIncrement, XtRDimension, sizeof(Dimension),
    Offset(scrollbar.page_increment), XtRImmediate, (XtPointer) 20
  },
  {
    XtNrepeatDelay, XtCRepeatDelay, XtRInt, sizeof(int),
    Offset(scrollbar.delay), XtRImmediate, (XtPointer) 150
  },
  {
    XtNlength, XtCLength, XtRDimension, sizeof(Dimension),
    Offset(scrollbar.along), XtRImmediate, (XtPointer) 1},

  {
    XtNthickness, XtCThickness, XtRDimension, sizeof(Dimension),
    Offset(scrollbar.cross), XtRImmediate, (XtPointer) 10
  },
  {
    XtNorientation, XtCOrientation, XtROrientation, sizeof(XtOrientation),
    Offset(scrollbar.orientation), XtRImmediate, (XtPointer) XtorientVertical
  },
  {
    XtNstripeColor, XtCStripeColor, XtRPixel, sizeof(Pixel),
    Offset(core.background_pixel), XtRCallProc, (XtPointer) Def_pixel
  },
  {
    XtNthumbColor, XtCThumbColor, XtRPixel, sizeof(Pixel),
    Offset(scrollbar.foreground), XtRCallProc, (XtPointer) Def_pixel
  },
  {
    XtNshown, XtCShown, XtRFloat, sizeof(float),
    Offset(scrollbar.shown), XtRFloat, (XtPointer)&floatZero
  },
  {
    XtNtopOfThumb, XtCTopOfThumb, XtRFloat, sizeof(float),
    Offset(scrollbar.top), XtRFloat, (XtPointer)&floatZero
  },
  {
    XtNminimumThumb, XtCMinimumThumb, XtRDimension, sizeof(Dimension),
    Offset(scrollbar.min_thumb), XtRImmediate, (XtPointer) 7
  },
  {
    XtNshowArrows, XtCShowArrows, XtRBoolean, sizeof (Boolean),
    Offset(scrollbar.showArrows), XtRImmediate, (XtPointer) True
  },

      /*
       *                C A L L B A C K S
       *
       *
       */

  {
    XtNincrementProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.incrementProc), XtRCallback, (XtPointer)NULL
  },
  {
    XtNdecrementProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.decrementProc), XtRCallback, (XtPointer)NULL
  },
  {
    XtNpageIncrementProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.pageIncrementProc), XtRCallback, (XtPointer)NULL
  },
  {
    XtNpageDecrementProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.pageDecrementProc), XtRCallback, (XtPointer)NULL
  },
  {
    XtNvalueChangedProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.valueChangedProc), XtRCallback, (XtPointer)NULL
  },
  {
    XtNdragProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.dragProc), XtRCallback, (XtPointer)NULL
  },
  {
    XtNtoTopProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.toTopProc), XtRCallback, (XtPointer)NULL
  },
  {
    XtNtoBottomProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.toBottomProc), XtRCallback, (XtPointer)NULL
  },
  {
    XtNscrollProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.scrollProc), XtRCallback, NULL
  },
  {
    XtNthumbProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.thumbProc), XtRCallback, NULL
  },
  {
    XtNjumpProc, XtCCallback, XtRCallback, sizeof(XtPointer),
    Offset(scrollbar.jumpProc), XtRCallback, NULL
  },

      /*
       *
       *                   O V E R W R I T E
       *
       *
       *
       */

  {
    XtNhighlightThickness, XtCHighlightThickness, XtRDimension,
    sizeof(Dimension),
    XtOffsetOf(SimpleRec,simple.highlight_thickness), XtRCallProc,
    (XtPointer)Def_highlight_thickness 
  },
  {
    XtNshadowWidth, XtCShadowWidth, XtRDimension, sizeof(Dimension),
    XtOffsetOf(SimpleRec, simple.shadow_thickness), XtRImmediate,
    (XtPointer) 2
  },
  {
    XtNcursor, XtCCursor, XtRCursor, sizeof(Cursor),
    Offset(simple.cursor), XtRString, "left_ptr"
  }
};


/* ARGSUSED */
static void Def_highlight_thickness( w, offset, value)
     Widget w;
     int offset;
     XrmValue *value;
{
  Widget parent = XtParent(w);
  static Dimension high;
  
  if ( XtIsSubclass(parent, containerWidgetClass))
    high = 2;
  else
    high = 0;

  value->addr = (caddr_t)& high;
}

static void Def_pixel( w, offset, value)
     Widget w;
     int offset;
     XrmValue *value;
{
  Widget parent = XtParent(w);
  static Pixel pixel;

  if (offset ==  Offset(core.background_pixel))
    ArmedColor (w, parent->core.background_pixel, &pixel);
  else
    pixel = parent->core.background_pixel;
  
  value->addr = (caddr_t)&pixel;
}
#undef Offset

static void ClassInitialize();
static void Initialize();
static void Destroy();
static void Realize();
static void Resize();
static void Redisplay();
static Boolean SetValues();
static XtGeometryResult QueryGeometry();

//static void XawDrawFrameWindow();
static void Select();
static void Moved();
static void Jumped();
static void Release();
static void PageUpOrLeft();
static void PageDownOrRight();
static void IncrementUpOrLeft();
static void IncrementDownOrRight();
static void ToTop();
static void ToBottom();

#ifdef like_motif

static char defaultTranslations[] =
    "<Btn1Down>:     Select()                    \n\
     <Btn1Motion>:   Moved()                     \n\
     <BtnUp>:        Release()                   \n\
     <Btn2Down>:     Jumped()                    \n\
     Ctrl<Key>Right: PageDownOrRight(Right)      \n\
     Ctrl<Key>Left:  PageUpOrLeft(Left)          \n\
     Ctrl<Key>Down:  PageDownOrRight(Down)       \n\
     Ctrl<Key>Up:    PageUpOrLeft(Up)            \n\
     Meta<Key>Right: ToBottom(Right)             \n\
     Meta<Key>Left:  ToTop(Left)                 \n\
     Meta<Key>Down:  ToBottom(Down)              \n\
     Meta<Key>Up:    ToTop(Up)                   \n\
     <Key>Right:     IncrementDownOrRight(Right) \n\
     <Key>Left:      IncrementUpOrLeft(Left)     \n\
     <Key>Down:      IncrementDownOrRight(Down)  \n\
     <Key>Up:        IncrementUpOrLeft(Up)";

static XtActionsRec actions[] = {
	{"Select",                   Select},
	{"Moved",                    Moved},
	{"Release",                  Release},
	{"Jumped",                   Jumped},
        {"PageDownOrRight",          PageDownOrRight},
        {"PageUpOrLeft",             PageUpOrLeft},
        {"IncrementDownOrRight",     IncrementDownOrRight},
        {"IncrementUpOrLeft",        IncrementUpOrLeft},
	{"ToTop",                    ToTop},
	{"ToBottom",                 ToBottom},
};

#else

static char defaultTranslations[] =
    "<Btn1Down>:     StartScroll()                 \n\
     <Btn2Down>:     StartScroll() NotifyThumb()   \n\
     <Btn3Down>:     StartScroll()                 \n\
     <Btn1Motion>:   Moved()                       \n\
     <Btn2Motion>:   NotifyThumb()                 \n\
     <BtnUp>:        NotifyScroll() EndScroll()    \n\
     Ctrl<Key>Right: PageDownOrRight(Right)      \n\
     Ctrl<Key>Left:  PageUpOrLeft(Left)          \n\
     Ctrl<Key>Down:  PageDownOrRight(Down)       \n\
     Ctrl<Key>Up:    PageUpOrLeft(Up)            \n\
     Meta<Key>Right: ToBottom(Right)             \n\
     Meta<Key>Left:  ToTop(Left)                 \n\
     Meta<Key>Down:  ToBottom(Down)              \n\
     Meta<Key>Up:    ToTop(Up)                   \n\
     <Key>Right:     IncrementDownOrRight(Right) \n\
     <Key>Left:      IncrementUpOrLeft(Left)     \n\
     <Key>Down:      IncrementDownOrRight(Down)  \n\
     <Key>Up:        IncrementUpOrLeft(Up)";


static XtActionsRec actions[] = {
	{"StartScroll",              Select},
	{"EndScroll",                Release},
	{"NotifyScroll",             Release},
	{"NotifyThumb",              Jumped},
	{"MoveThumb",                Jumped},
	{"Moved",                    Moved},
        {"PageDownOrRight",          PageDownOrRight},
        {"PageUpOrLeft",             PageUpOrLeft},
        {"IncrementDownOrRight",     IncrementDownOrRight},
        {"IncrementUpOrLeft",        IncrementUpOrLeft},
	{"ToTop",                    ToTop},
	{"ToBottom",                 ToBottom},
};


#endif

#define SuperClass ((SimpleWidgetClass)&simpleClassRec)

ScrollbarClassRec scrollbarClassRec = {
  { /* core fields */
    /* superclass       */      (WidgetClass) SuperClass,
    /* class_name       */      "Scrollbar",
    /* size             */      sizeof(ScrollbarRec),
    /* class_initialize	*/	ClassInitialize,
    /* class_part_init  */	NULL,
    /* class_inited	*/	FALSE,
    /* initialize       */      Initialize,
    /* initialize_hook  */	NULL,
    /* realize          */      Realize,
    /* actions          */      actions,
    /* num_actions	*/	XtNumber(actions),
    /* resources        */      resources,
    /* num_resources    */      XtNumber(resources),
    /* xrm_class        */      NULLQUARK,
    /* compress_motion	*/	TRUE,
    /* compress_exposure*/	TRUE,
    /* compress_enterleave*/	TRUE,
    /* visible_interest */      TRUE,
    /* destroy          */      Destroy,
    /* resize           */      Resize,
    /* expose           */      Redisplay,
    /* set_values       */      SetValues,
    /* set_values_hook  */	NULL,
    /* set_values_almost */	XtInheritSetValuesAlmost,
    /* get_values_hook  */	NULL,
    /* accept_focus     */      NULL,
    /* version          */	XtVersion,
    /* callback_private */      NULL,
    /* tm_table         */      defaultTranslations,
    /* query_geometry	*/	QueryGeometry,
    /* display_accelerator*/	XtInheritDisplayAccelerator,
    /* extension        */	NULL
  },
  { /* simple fields */
    /* change_sensitive */      XtInheritChangeSensitive,
    /* display_rect     */      XtInheritDisplayRectProc,
    /* extension        */      NULL
  },
  { /* scrollbar fields */
    /* ignore           */      0
  }

};

WidgetClass scrollbarWidgetClass = (WidgetClass)&scrollbarClassRec;

/*
 *
 *  Scrollbar state
 *
 */
     
#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define FREE_SCROLL     (0)
#define BESY_SCROLL     (1)     
#define MOVE_SCROLL     (2)
     
#define NoButton        -1
#define LOCATION        True
#define PROPORTION      False
#define InRange(n,s,b)  MAX(s,MIN(n, b))
#define MARGIN          0
#define IS_VERTICAL(w)  if(BAR(w).orientation == XtorientVertical)
#define SHADOW_WIDTH(w) (((SimpleWidget)w)->simple.shadow_thickness + \
			 ((SimpleWidget)w)->simple.highlight_thickness)
#define ARROW_SIZE(w)   (BAR(w).showArrows ? BAR(w).cross + 2 * MARGIN : 0)
#define INNER_MARGIN(w) (SHADOW_WIDTH(w) + MARGIN)

#define LOW_LIMIT(w)     (int)(INNER_MARGIN(w) + ARROW_SIZE(w))
#define HIGH_LIMIT(w)    (int)(BAR(w).along + INNER_MARGIN(w) + ARROW_SIZE(w) \
			      - BAR(w).shownLength)
#define Xaw_UP          True
#define Xaw_DOWN        False
#define Xaw_BOTTOM      (1 << 1)
#define Xaw_TOP         (1 << 2)


#define THUMB_SHADOW_THICKNESS(w) 2

#define DO_CALLBACK(w,callback,data) \
       if (XtCallbackHasSome == XtHasCallbacks(w, callback)) \
       XtCallCallbacks ((Widget)w, callback, data) 

#define ADD_TIMEOUT(w,proc) \
  XtAppAddTimeOut (XtWidgetToApplicationContext ((Widget) w), \
		   (unsigned long) BAR(w).delay, proc, (XtPointer)w)

#define CLEAR_TIMEOUT(w) \
  if (BAR(w).timer) { \
      XtRemoveTimeOut (BAR(w).timer); \
      BAR(w).timer = 0; \
  }

static void PaintArrows();
static void TimerIncrementProc();
static void TimerDecrementProc();
static void FillXawScrollBarCallbackStruct();


static Dimension int_to_Dimension(i, n)
     register int i;
     register int n;
{
  i = MAX (i, n);
  return (Dimension) i;
}

#define INT_DIM_0(i) int_to_Dimension(i, 0)
#define INT_DIM_1(i) int_to_Dimension(i, 1)

static void ClassInitialize()
{
    XawInitializeWidgetSet();
    XtAddConverter( XtRString, XtROrientation, XmuCvtStringToOrientation,
		   (XtConvertArgList)NULL, (Cardinal)0 );
}


static float exectly(shown)
     float shown;
{
  float top;
  float tmp;
  
  top = 1.0 - shown;
  tmp = top + shown - 1.0;
  top = 1.0 - shown - tmp;

  if ((tmp = top + shown) > 1.0)
  {
    tmp = 2.0 - tmp;
    top = tmp - shown;
  }
  
  return top;
}


static int CalculateThumb( sw, x, y, w, h, flag)
     ScrollbarWidget sw;
     int            *x;
     int            *y;
     unsigned int   *w;
     unsigned int   *h;
     Boolean        flag;
{
  int cross, along;
  int inner_margin;
  int new_top;
  int thumb_size;
  
  inner_margin = INNER_MARGIN(sw);
  
  cross = BAR(sw).cross;
  along = BAR(sw).along;
  
  if (along < BAR(sw).min_thumb)
    return TRUE;

  thumb_size = (int)((double)along * (double)BAR(sw).shown);
  BAR(sw).shownLength = MAX (thumb_size, (int)BAR(sw).min_thumb);
  
  if (flag == LOCATION)
    new_top = BAR(sw).newLoc - (inner_margin + ARROW_SIZE(sw));
  else 
    new_top = (along * BAR(sw).top);

  if (flag == LOCATION) {
    double tmp = (double) (along + thumb_size - BAR(sw).shownLength);
    BAR(sw).top = (double)(new_top) / tmp;
  }
  
  new_top = InRange (new_top, 0, along - (int)BAR(sw).shownLength);

  
  if (BAR(sw).top + BAR(sw).shown > 1.0)
    BAR(sw).top = exectly(BAR(sw).shown);

  new_top += inner_margin + ARROW_SIZE(sw);
  
  IS_VERTICAL(sw)
  {
    *x = (int)inner_margin;
    *y = (int)new_top;
    *w = (unsigned int)cross;
    *h = (unsigned int)BAR(sw).shownLength;
  }
  else
  {
    *x = (int)new_top;
    *y = (int)inner_margin;
    *w = (unsigned int)BAR(sw).shownLength;
    *h = (unsigned int)cross;
  } 

  BAR(sw).newLoc      =
  BAR(sw).topLoc      = new_top;

  return FALSE;

}


/*
 *                  PaintThumb
 *
 *
 *  Paint the thumb in the area specified by w->top and w->shown.
 *  The old area is erased.  The painting and erasing is done cleverly
 *  so that no flickering will occur.
 *
 */

static void PaintThumb( gw ,region, flag)
     Widget gw;
     Region region;
     Boolean flag;
{
  ScrollbarWidget sbw = (ScrollbarWidget)gw;
  int          x, y;
  unsigned int w, h;
  Position oldtop, oldbot, newtop, newbot;
  
  oldbot = (oldtop = sbw->scrollbar.topLoc) + (Position)BAR(sbw).shownLength;
  
  if (CalculateThumb( sbw, &x, &y, &w, &h, flag))
    return;

  if (!XtIsRealized(gw) || !sbw->core.visible)
    return;

  IS_VERTICAL(sbw)
  {
    newtop = newbot = y;
    newbot += h;
  }
  else
  {
    newtop = newbot = x;
    newbot += w;
  }

  if (newtop > oldtop)
  {
    IS_VERTICAL(sbw)
    {
      XClearArea(XtDisplay(gw), XtWindow(gw),
		 INNER_MARGIN(sbw),
		 oldtop,
		 BAR(sbw).cross,
		 (newtop > oldbot) ? oldbot - oldtop  + 1: newtop - oldtop,
		 FALSE);
    }
    else
    {
      XClearArea(XtDisplay(gw), XtWindow(gw),
		 oldtop,
		 INNER_MARGIN(sbw),
		 MIN(newtop, oldbot) - oldtop,
		 BAR(sbw).cross,
		 FALSE);
    }
  }

  if ( newbot < oldbot)
  {
    IS_VERTICAL(sbw)
    {
      XClearArea(XtDisplay(gw),XtWindow(gw),
		 INNER_MARGIN(sbw),
		 MAX(newbot, oldtop),
		 BAR(sbw).cross,
		 (newbot > oldtop) ? oldbot - newbot: oldbot - oldtop + 1,
		 FALSE);
    }
    else
    {
      XClearArea(XtDisplay(gw),XtWindow(gw),
		 MAX(newbot, oldtop),
		 INNER_MARGIN(sbw),
		 oldbot - MAX(newbot, oldtop),
		 BAR(sbw).cross,
		 FALSE);
      
    }
  }

  if ( region == NULL || XRectInRegion(region, x, y, w, h) != RectangleOut )
  {
    XFillRectangle(XtDisplay(gw),
		   XtWindow(gw),
		   BAR(sbw).gc,
		   x + sbw->simple.shadow_thickness,
		   y + sbw->simple.shadow_thickness,
		   (unsigned int)(w - 2 * sbw->simple.shadow_thickness),
		   (unsigned int)(h - 2 * sbw->simple.shadow_thickness));
    
    XawDrawFrame( gw,
		 (Position)x,
		 (Position)y,
		 (Dimension)w,
		 (Dimension)h,
		 XawRAISED,
		 THUMB_SHADOW_THICKNESS(gw) /*SIMPLE(sbw).shadow_thickness*/,
		 SIMPLE(sbw).top_shadow_GC,
		 SIMPLE(sbw).bottom_shadow_GC);
  }
  
  XFlush(XtDisplay(gw));

}

static void Destroy(w)
     Widget w;
{
  XtReleaseGC(w, BAR(w).gc);
}

static void  CreateGC(w)
     Widget w;
{
  XGCValues gcValues;
  XtGCMask mask;
  
  mask = GCForeground ;
  gcValues.foreground = BAR(w).foreground;
  BAR(w).gc = XtGetGC( w, mask, &gcValues);
}

static void PositionChild(gw)
     Widget gw;
{
  int x, y;
  unsigned int width, height;
  unsigned int shadow_thickness, arrow_size;

  if (!BAR(gw).showArrows)
    return ;

  shadow_thickness  = SHADOW_WIDTH(gw);
  arrow_size    = ARROW_SIZE(gw);
  
  x      = shadow_thickness;
  y      = shadow_thickness;
  width  = arrow_size;
  height = arrow_size;
  
  XawMakeDrawArrowStruct( x, y, width, height, THUMB_SHADOW_THICKNESS(gw),
			 ( (BAR(gw).orientation == XtorientVertical)
			  ? XawTop : XawLeft),
			 (XawDrawArrowStruct*)&BAR(gw).top_arrow);
  
  x      = gw->core.width  - arrow_size - shadow_thickness;
  y      = gw->core.height - arrow_size - shadow_thickness;
/*width  = arrow_size; the same */
/*height = arrow_size; the same */
  
  XawMakeDrawArrowStruct( x, y, width, height, THUMB_SHADOW_THICKNESS(gw),
			 ( (BAR(gw).orientation == XtorientVertical)
			  ? XawBottom : XawRight),
			 (XawDrawArrowStruct*)&BAR(gw).bot_arrow);

}


/* ARGSUSED */
static void Initialize( request, new )
   Widget request;
   Widget new;
{
  ScrollbarWidget w = (ScrollbarWidget) new;
  Dimension       arrow_size;
  
  
  CreateGC(new);

  BAR(w).delay = MAX(BAR(w).delay, 75);
  BAR(w).timer = 0;

  /*
   * Recalculate BAR(w).cross if width and height of Scrollbar are set.
   */
  IS_VERTICAL(w)
  {
    if (w->core.width != 0)
      BAR(w).cross =
	INT_DIM_1 (w->core.width - 2 * INNER_MARGIN(w));
  }
  else
  {
    if (w->core.height != 0)
      BAR(w).cross =
	INT_DIM_1 (w->core.height - 2 * INNER_MARGIN(w)); 
  }

  BAR(w).along = MAX (BAR(w).along, BAR(w).min_thumb);
  /*
   * Calculate size of arrow in the dependence of BAR(w).cross .
   */
  arrow_size = BAR(w).cross + 2 * MARGIN;

  

  /*
   * Recalculate BAR(w).along if width and height of Scrollbar are set.
   */
  IS_VERTICAL(w)
  {
    if (w->core.height != 0)
      BAR(w).along =
	INT_DIM_1 (w->core.height - 2 * (INNER_MARGIN(w) + arrow_size)); 
  }
  else
  {
    if (w->core.width != 0)
      BAR(w).along =
	INT_DIM_1 (w->core.width - 2 * (INNER_MARGIN(w) + arrow_size));
  }

  BAR(w).along = MAX (BAR(w).along, BAR(w).min_thumb);

  

  /*
   * Recalculate width and height of Scrollbar if they are set to zero
   * by taking into acount BAR(w).cross, BAR(w).along, arrow_size .
   */
  if (w->core.width == 0) {
    IS_VERTICAL(w){
      w->core.width = BAR(w).cross + 2 * INNER_MARGIN(w);
    }else{
      w->core.width = BAR(w).along + 2 * (INNER_MARGIN(w) + arrow_size);
    }
  }
  
  if (w->core.height == 0) {
    IS_VERTICAL(w){
      w->core.height = BAR(w).along + 2 * (INNER_MARGIN(w) + arrow_size);
    }else{
      w->core.height = BAR(w).cross + 2 * INNER_MARGIN(w);
    }
  }

  arrow_size = BAR(w).cross + 2 * MARGIN;

  BAR(w).scrolling   = FREE_SCROLL;
  BAR(w).topLoc      = INNER_MARGIN(w) + arrow_size;
  BAR(w).shownLength = BAR(w).along;

  if (BAR(w).showArrows)
    PositionChild(new);

}

static void Realize( gw, valueMask, attributes )
   Widget gw;
   Mask *valueMask;
   XSetWindowAttributes *attributes;
{
  /* 
   * The Simple widget actually stuffs the value in the valuemask. 
   */
  
  (*SuperClass->core_class.realize) (gw, valueMask, attributes);

}

/* ARGSUSED */
static Boolean SetValues( current, request, desired, args, num_args )
     Widget    current, request, desired;
     ArgList   args;
     Cardinal *num_args;
{
  ScrollbarWidget old = (ScrollbarWidget) current;
  ScrollbarWidget new = (ScrollbarWidget) desired;
  Boolean         redraw = FALSE;
  
  /*
   * If these values are outside the acceptable range ignore them...
   */
  
  if (BAR(new).top < 0.0 || BAR(new).top > 1.0)
    BAR(new).top = BAR(old).top;
  
  if (BAR(new).shown < 0.0 || BAR(new).shown > 1.0)
    BAR(new).shown = BAR(old).shown;
  
  if(BAR(new).delay < 50)
    BAR(new).delay = BAR(old).delay;


#define NE(field) (BAR(old).field != BAR(new).field)    
  
  if (NE(cross))
  { 
    IS_VERTICAL(new)
    {
      new->core.width = BAR(new).cross + 2 * INNER_MARGIN(new);
    }
    else
    {
      new->core.height = BAR(new).cross + 2 * INNER_MARGIN(new);
    }
  }

  
  if (NE(along))
  { 
    IS_VERTICAL(new) {
      new->core.height =
	BAR(new).along + 2 * (INNER_MARGIN(new) + ARROW_SIZE(new));
    }
    else {
      new->core.width =
	BAR(new).along + 2 * (INNER_MARGIN(new) + ARROW_SIZE(new));
    }
  }

  if (new->core.width == 0) 
  {
    IS_VERTICAL(new) {
      new->core.width = BAR(new).cross + 2 * INNER_MARGIN(new);
    }
    else {
      new->core.width =
	BAR(new).along + 2 * (INNER_MARGIN(new) + ARROW_SIZE(new));
    }
  }
  
  if (new->core.height == 0) 
  {
    IS_VERTICAL(new) {
      new->core.height =
	BAR(new).along + 2 * (INNER_MARGIN(new) + ARROW_SIZE(new));
    }
    else {
      new->core.height = BAR(new).cross + 2 * INNER_MARGIN(new);
    }
  }

  if ( NE(showArrows) )
  {
    if (BAR(new).showArrows)
      PositionChild(current);

    redraw = TRUE;
  }

  if ( old->core.background_pixel != new->core.background_pixel)
    redraw = TRUE;
  
  if (NE(foreground))
  {
    XtReleaseGC (current, BAR(new).gc);
    CreateGC (current);
    redraw = TRUE;
  }
  
  if (NE(top) || NE(shown))
    redraw = TRUE;
  
  return (redraw);
}

static void Resize( gw )
   Widget gw;
{
  ScrollbarWidget w = (ScrollbarWidget) gw;

  IS_VERTICAL(w)
  {
    BAR(w).cross = INT_DIM_1(w->core.width - 2 * INNER_MARGIN(w));
    BAR(w).along =
      INT_DIM_1(w->core.height - 2 * (INNER_MARGIN(w) + ARROW_SIZE(w)));
  }
  else
  {
    BAR(w).cross = INT_DIM_1(w->core.height - 2 * INNER_MARGIN(w));
    BAR(w).along =
      INT_DIM_1(w->core.width - 2 * (INNER_MARGIN(w) + ARROW_SIZE(w)));
  }

  if (BAR(gw).showArrows)
    PositionChild(gw);

  if (XtIsRealized(gw) && gw->core.visible)
  {
    XClearArea(XtDisplay(gw),XtWindow(gw), 0, 0, 0, 0, FALSE);
    Redisplay (gw, (XEvent*)NULL, (Region)NULL);
  }

}


/* ARGSUSED */
static void Redisplay( gw, event, region )
   Widget gw;
   XEvent *event;
   Region region;
{

  if (XtIsRealized(gw) && gw->core.visible) {

    register SimpleWidget s = (SimpleWidget)gw;

     if ( XtIsSubclass(XtParent(gw), containerWidgetClass)){
      ContainerWidget cw = (ContainerWidget) XtParent(gw);
      XawDrawFrame (gw,
		    0,
		    0,
		    s->core.width,
		    s->core.height,
		    XawRAISED,
		    s->simple.highlight_thickness,
		    cw->container.background_GC,
		    cw->container.background_GC);
    } 
    
    XawDrawFrame (gw,
		  s->simple.highlight_thickness,
		  s->simple.highlight_thickness,
		  s->core.width - 2 * s->simple.highlight_thickness,
		  s->core.height - 2 * s->simple.highlight_thickness,
		  XawRAISED,
		  s->simple.shadow_thickness,
		  s->simple.bottom_shadow_GC,
		  s->simple.top_shadow_GC);
    
    PaintArrows(gw, Xaw_UP, (unsigned int)(Xaw_BOTTOM | Xaw_TOP));

    PaintThumb( gw, region, PROPORTION);

  }
}

/* ARGSUSED */
static void PaintArrows(gw, up_down, top_bottom)
     Widget       gw;
     Boolean      up_down;
     unsigned int top_bottom;
{
  if (!BAR(gw).showArrows || !XtIsRealized(gw) || !gw->core.visible)
    return ;

  if ( up_down == Xaw_UP )
  {
    if (top_bottom & Xaw_TOP)
      XawDrawArrow(gw,
		   BAR(gw).gc,
		   SIMPLE(gw).top_shadow_GC,
		   SIMPLE(gw).bottom_shadow_GC,
		   (XawDrawArrowStruct*)&BAR(gw).top_arrow);

    if (top_bottom & Xaw_BOTTOM)
      XawDrawArrow(gw,
		   BAR(gw).gc,
		   SIMPLE(gw).top_shadow_GC,
		   SIMPLE(gw).bottom_shadow_GC,
		   (XawDrawArrowStruct*)&BAR(gw).bot_arrow);
  }
  else
  {
    if (top_bottom & Xaw_TOP)
      XawDrawArrow(gw,
		   BAR(gw).gc,
		   SIMPLE(gw).bottom_shadow_GC,
		   SIMPLE(gw).top_shadow_GC,
		   (XawDrawArrowStruct*)&BAR(gw).top_arrow);
    
    if (top_bottom & Xaw_BOTTOM)
      XawDrawArrow(gw,
		   BAR(gw).gc,
		   SIMPLE(gw).bottom_shadow_GC,
		   SIMPLE(gw).top_shadow_GC,
		   (XawDrawArrowStruct*)&BAR(gw).bot_arrow);
  } 
  
}


static void ExtractPosition( event, x, y )
    XEvent *event;
    Position *x, *y;		/* RETURN */
{
    switch( event->type ) {
      case MotionNotify:
		*x = event->xmotion.x;	 *y = event->xmotion.y;
		break;
      case ButtonPress:
      case ButtonRelease:
		*x = event->xbutton.x;   *y = event->xbutton.y;
		break;
      case KeyPress:
      case KeyRelease:
		*x = event->xkey.x;      *y = event->xkey.y;
		break;
      case EnterNotify:
      case LeaveNotify:
		*x = event->xcrossing.x; *y = event->xcrossing.y;
		break;
      default:
		*x = 0; *y = 0;
    }
}

/* ARGSUSED */
static void Select( gw, event, params, num_params )
  Widget gw;
  XEvent *event;
  String *params;		/* direction: Back|Forward|Smooth */
  Cardinal *num_params;		/* we only support 1 */
{
  ScrollbarWidget w = (ScrollbarWidget) gw;
  Position  x,y;
  Position  n;
  unsigned int p;
  unsigned int tmp;
  XawScrollBarCallbackStruct bar_data;
  
  if (BAR(w).scrolling)
    return;                /* if we're already in progress */
  
  ExtractPosition( event, &x, &y );

  IS_VERTICAL(w)
  {
    n = y;
    p = w->core.height;
  }
  else
  {
    n = x;
    p = w->core.width;
  }

  tmp = (unsigned int)(ARROW_SIZE(w) + INNER_MARGIN(w));


  if(BAR(w).showArrows && (n > (Position)(p - tmp)))     
  {
    /***************************************************************
     *
     *
     *           Bottom arrow
     *
     *
     ***************************************************************/

    int newLoc;
    

    newLoc = BAR(w).topLoc + BAR(w).increment;

    newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));
    

    if (newLoc == (int)BAR(w).topLoc)
      return;

    BAR(w).newLoc = newLoc;
    
    PaintThumb(gw, (Region)NULL, LOCATION);
    PaintArrows(gw, Xaw_DOWN, (unsigned int)Xaw_BOTTOM);
    
    FillXawScrollBarCallbackStruct(gw, &bar_data, XawSB_INCREMENT, event);
    
    DO_CALLBACK (gw, XtNincrementProc, (XtPointer)&bar_data);
    else
    DO_CALLBACK (gw, XtNvalueChangedProc,(XtPointer)&bar_data);

    {
      float top = BAR(w).top;
      DO_CALLBACK (gw,XtNthumbProc, *(XtPointer*)&top);
      DO_CALLBACK (gw,XtNjumpProc, (XtPointer)&top);
    }
    
    XFlush(XtDisplay(gw));
    BAR(w).timer = ADD_TIMEOUT(w, TimerIncrementProc);

  }


  else if( n > BAR(w).topLoc + BAR(w).shownLength )
  {
    /***************************************************************
     *
     *
     *           Bottom zone ( under thumb )
     *
     *
     ***************************************************************/

    PageDownOrRight( gw, event, params, num_params );     
    {
      float top = BAR(w).top;
      DO_CALLBACK (gw,XtNthumbProc, *(XtPointer*)&top);
      DO_CALLBACK (gw,XtNjumpProc, (XtPointer)&top);
    }

  }


  else if( n > BAR(w).topLoc )
  {
    /***************************************************************
     *
     *
     *           Thumb
     *
     *
     ***************************************************************/

    BAR(w).scrolling = MOVE_SCROLL;                       
    BAR(w).shift = n - BAR(w).topLoc;
  }


  else if( BAR(w).showArrows && (n > (Position)tmp))
  {
    /***************************************************************
     *
     *
     *           Top zone ( above thumb )
     *
     *
     ***************************************************************/

    PageUpOrLeft( gw, event, params, num_params );        
    {
      float top = BAR(w).top;

      DO_CALLBACK (gw,XtNthumbProc, *(XtPointer*)&top);
      DO_CALLBACK (gw,XtNjumpProc, (XtPointer)&top);
    }
  }


  else
  {
    /***************************************************************
     *
     *
     *           Top arrow
     *
     *
     ***************************************************************/

    int newLoc;

    newLoc = BAR(w).topLoc - BAR(w).increment;

    newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));

    if (newLoc == BAR(w).topLoc)
      return;

    BAR(w).newLoc = newLoc;
    
    PaintThumb(gw, (Region)NULL, LOCATION);
    PaintArrows(gw, Xaw_DOWN, (unsigned int)Xaw_TOP);
    
    FillXawScrollBarCallbackStruct(gw, &bar_data, XawSB_DECREMENT, event);
    
    DO_CALLBACK (gw, XtNdecrementProc, (XtPointer)&bar_data);
    else
    DO_CALLBACK (gw, XtNvalueChangedProc,(XtPointer)&bar_data);

    {
      float top = BAR(w).top;
      DO_CALLBACK (gw,XtNthumbProc, *(XtPointer*)&top);
      DO_CALLBACK (gw,XtNjumpProc, (XtPointer)&top);
    }

    XFlush(XtDisplay(gw));
    BAR(w).timer = ADD_TIMEOUT(w, TimerDecrementProc);
  }
  
}


static Boolean AreDifferentEvents( oldEvent, newEvent )
    XEvent *oldEvent, *newEvent;
{
#define Check(field) if (newEvent->field != oldEvent->field) return True

  Check(xany.display);
  Check(xany.type);
  Check(xany.window);
  
  switch( newEvent->type ) {
  case MotionNotify:
                     Check(xmotion.state);
		     break;
  case ButtonPress:
  case ButtonRelease:
                     Check(xbutton.state);
		     Check(xbutton.button);
		     break;
  case KeyPress:
  case KeyRelease:
		     Check(xkey.state);
		     Check(xkey.keycode);
		     break;
  case EnterNotify:
  case LeaveNotify:
		     Check(xcrossing.mode);
		     Check(xcrossing.detail);
		     Check(xcrossing.state);
		     break;
  }

#undef Check

    return False;
}

struct EventData {
	XEvent *oldEvent;
	int count;
};

static Bool PeekNotifyEvent( dpy, event, args )
    Display *dpy;
    XEvent *event;
    char *args;
{
    struct EventData *eventData = (struct EventData*)args;

    return ((++eventData->count == QLength(dpy)) /* since PeekIf blocks */
	    || AreDifferentEvents(event, eventData->oldEvent));
}


static Boolean LookAhead( w, event )
    Widget w;
    XEvent *event;
{
    XEvent newEvent;
    struct EventData eventData;

    if (QLength(XtDisplay(w)) == 0) return False;

    eventData.count = 0;
    eventData.oldEvent = event;

    XPeekIfEvent(XtDisplay(w), &newEvent, PeekNotifyEvent, (char*)&eventData);

    if (AreDifferentEvents(event, &newEvent))
      return True;
    else
      return False;
}

#ifdef notdef

static void LookAhead( w, event )
    Widget w;
    XEvent *event;
{
    XEvent   ahead;
    Display *dpy = XtDisplay(w);
    
    while( XPending(dpy) > 0)
    {
      XPeekEvent(dpy, &ahead);
      if (AreDifferentEvents(event, &ahead))
	break;
      XNextEvent( dpy, event);
    }
}

#endif


/* ARGSUSED */
static void Release(w, event, params, num_params )
   Widget w;
   XEvent *event;		/* unused */
   String *params;		/* unused */
   Cardinal *num_params;	/* unused */
{
  int position;
  int scrolling = BAR(w).scrolling;
  
  BAR(w).scrolling = FREE_SCROLL;
  CLEAR_TIMEOUT(w);  
  PaintArrows(w, Xaw_UP, (unsigned int)(Xaw_BOTTOM | Xaw_TOP));

#if 0
  FillXawScrollBarCallbackStruct(w, &bar_data, XawSB_VALUE_CHANGED, event);

  DO_CALLBACK (w, XtNvalueChangedProc, (XtPointer)&bar_data);
#endif
 
  if (scrolling == BESY_SCROLL)
  {
    position  = (int)((float)(BAR(w).along) * BAR(w).top);
    position += 2 * (ARROW_SIZE(w) + INNER_MARGIN(w) + SIMPLE_MARGIN(w));

    DO_CALLBACK (w, XtNscrollProc, (XtPointer)(intptr_t)position);
  }
  XFlush(XtDisplay(w));

}

/* ARGSUSED */
static void ToTop( gw, event, params, num_params )
   Widget    gw;
   XEvent   *event;             /* unused */
   String   *params;
   Cardinal *num_params;	/* unused */
{
  ScrollbarWidget w = (ScrollbarWidget) gw;
  XawScrollBarCallbackStruct bar_data;
  
  if (BAR(w).scrolling)
    return;                               /* if we're already scrolling */

  if (*num_params > 0) {
    if ((*params[0] == 'L' || *params[0] == 'l') &&
	BAR(w).orientation == XtorientVertical )
      return;
    if ((*params[0] == 'U' || *params[0] == 'u') &&
	BAR(w).orientation != XtorientVertical )
      return;
  }
  
  BAR(w).top = 0.0;
  PaintThumb(gw, (Region)NULL, PROPORTION);

  FillXawScrollBarCallbackStruct(gw, (XawScrollBarCallbackStruct*)&bar_data,
				 XawSB_TOP, event);

  DO_CALLBACK (gw, XtNtoTopProc, (XtPointer)&bar_data);
  else
  DO_CALLBACK (gw, XtNvalueChangedProc, (XtPointer)&bar_data);

  XFlush(XtDisplay(gw));

}

/* ARGSUSED */
static void ToBottom( gw, event, params, num_params )
   Widget gw;
   XEvent *event;               /* unused */
   String *params;
   Cardinal *num_params;	/* unused */
{
  ScrollbarWidget w = (ScrollbarWidget) gw;
  XawScrollBarCallbackStruct bar_data;
  
  if (BAR(w).scrolling)
    return;                             /* if we're already scrolling */

  if (*num_params > 0) {
    if ((*params[0] == 'R' || *params[0] == 'r') &&
	BAR(w).orientation == XtorientVertical )
      return;
    if ((*params[0] == 'D' || *params[0] == 'd') &&
	BAR(w).orientation != XtorientVertical )
      return;
  }
  
  BAR(w).top = 1.0;
  PaintThumb(gw, (Region)NULL, PROPORTION);

  FillXawScrollBarCallbackStruct(gw, (XawScrollBarCallbackStruct*)&bar_data,
				 XawSB_BOTTOM, event);

  DO_CALLBACK (gw, XtNtoBottomProc, (XtPointer)&bar_data);
  else 
  DO_CALLBACK (gw, XtNvalueChangedProc, (XtPointer)&bar_data);

  XFlush(XtDisplay(gw));

}

/* ARGSUSED */
static void PageUpOrLeft( gw, event, params, num_params )
   Widget    gw;
   XEvent   *event;             /* unused */
   String   *params;		/* unused */
   Cardinal *num_params;	/* unused */
{
  ScrollbarWidget w = (ScrollbarWidget) gw;
  XawScrollBarCallbackStruct bar_data;
  int newLoc;
  
  if (BAR(w).scrolling)
    return;                               /* if we're already scrolling */

  if (*num_params > 0) {
    if ((*params[0] == 'L' || *params[0] == 'l') &&
	BAR(w).orientation == XtorientVertical )
      return;
    if ((*params[0] == 'U' || *params[0] == 'u') &&
	BAR(w).orientation != XtorientVertical )
      return;
  }
  
  newLoc = BAR(w).topLoc - BAR(w).page_increment;

  newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));
  
  if (newLoc == BAR(w).topLoc)
    return;

  BAR(w).newLoc = newLoc;
  
  PaintThumb(gw, (Region)NULL, LOCATION);

  FillXawScrollBarCallbackStruct(gw, &bar_data, XawSB_PAGE_DECREMENT, event);

  DO_CALLBACK (gw, XtNpageDecrementProc, (XtPointer)&bar_data);
  else
  DO_CALLBACK (gw, XtNvalueChangedProc, (XtPointer)&bar_data);

  XFlush(XtDisplay(gw));

}


/* ARGSUSED */
static void PageDownOrRight( gw, event, params, num_params )
   Widget    gw;
   XEvent   *event;             /* unused */
   String   *params;		/* unused */
   Cardinal *num_params;	/* unused */
{
  ScrollbarWidget w = (ScrollbarWidget) gw;
  XawScrollBarCallbackStruct bar_data;
  int newLoc;
  
  if (BAR(w).scrolling)
    return;                             /* if we're already scrolling */

  if (*num_params > 0)
  {
    if ((*params[0] == 'R' || *params[0] == 'r') &&
	BAR(w).orientation == XtorientVertical )
      return;
    if ((*params[0] == 'D' || *params[0] == 'd') &&
	BAR(w).orientation != XtorientVertical )
      return;
  }
  
  newLoc = BAR(w).topLoc + BAR(w).page_increment;
  
  newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));
  
  if (newLoc == BAR(w).topLoc)
    return;

  BAR(w).newLoc = newLoc;
  
  PaintThumb(gw, (Region)NULL, LOCATION);

  FillXawScrollBarCallbackStruct(gw, (XawScrollBarCallbackStruct*)&bar_data,
				 XawSB_PAGE_INCREMENT, event);

  DO_CALLBACK (gw, XtNpageIncrementProc, (XtPointer)&bar_data);
  else
  DO_CALLBACK (gw, XtNvalueChangedProc, (XtPointer)&bar_data);

  XFlush(XtDisplay(gw));
}

static void IncrementUpOrLeft( gw, event, params, num_params )
   Widget gw;
   XEvent *event;
   String *params;		/* unused */
   Cardinal *num_params;	/* unused */
{
  ScrollbarWidget w = (ScrollbarWidget) gw;
  XawScrollBarCallbackStruct bar_data;
  int newLoc;
  
  if (BAR(w).scrolling)
    return;                                 /* if we're already scrolling */

  if (*num_params > 0)
  {
    if ((*params[0] == 'L' || *params[0] == 'l') &&
	BAR(w).orientation == XtorientVertical )
      return;
    if ((*params[0] == 'U' || *params[0] == 'u') &&
	BAR(w).orientation != XtorientVertical )
      return;
  }
  
  newLoc = BAR(w).topLoc - BAR(w).increment;
  
  newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));
  
  if (newLoc == BAR(w).topLoc)
    return;
  
  BAR(w).newLoc = newLoc;
  
  PaintThumb(gw, (Region)NULL, LOCATION);
  
  FillXawScrollBarCallbackStruct(gw, &bar_data, XawSB_DECREMENT, event);

  DO_CALLBACK (gw, XtNdecrementProc, (XtPointer)&bar_data);
  else
  DO_CALLBACK (gw, XtNvalueChangedProc, (XtPointer)&bar_data);

  XFlush(XtDisplay(gw));
}


/* ARGSUSED */
static void IncrementDownOrRight( w, event, params, num_params )
   Widget    w;
   XEvent   *event;             /* unused */
   String   *params;
   Cardinal *num_params;	/* unused */
{
  XawScrollBarCallbackStruct bar_data;
  int newLoc;
  
  if (BAR(w).scrolling)
    return;                                    /* if we're already scrolling */

  if (*num_params > 0)
  {
    if ( (*params[0] == 'R' || *params[0] == 'r') &&
	BAR(w).orientation == XtorientVertical )
      return;
    if ( (*params[0] == 'D' || *params[0] == 'd') &&
	BAR(w).orientation != XtorientVertical )
      return;
  }

  newLoc = BAR(w).topLoc + BAR(w).increment;
  newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));

  if (newLoc == BAR(w).topLoc)
    return;

  BAR(w).newLoc = newLoc;

  PaintThumb(w, (Region)NULL, LOCATION);

  FillXawScrollBarCallbackStruct(w, &bar_data, XawSB_INCREMENT, event);

  DO_CALLBACK (w, XtNincrementProc, (XtPointer)&bar_data);
  else
  DO_CALLBACK (w, XtNvalueChangedProc, (XtPointer)&bar_data);

  XFlush(XtDisplay(w));
}



/* ARGSUSED */
static void Moved( w, event, params, num_params )
   Widget w;
   XEvent *event;
   String *params;		/* unused */
   Cardinal *num_params;	/* unused */
{
  Position                   x, y;
  int                   newLoc;
  XawScrollBarCallbackStruct bar_data;
  
  if (BAR(w).scrolling == FREE_SCROLL)
    return;                                  /* if no Select */
  
  if (LookAhead(w, event))
    return;
  
  if (!event->xmotion.same_screen)
    return;

  BAR(w).scrolling = MOVE_SCROLL;
  
  ExtractPosition( event, &x, &y );

  IS_VERTICAL(w) 
    newLoc = y - BAR(w).shift;
  else
    newLoc = x - BAR(w).shift;

  newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));

  if (newLoc == BAR(w).topLoc)
    return;

  BAR(w).newLoc = newLoc;
  
  PaintThumb(w, (Region)NULL, LOCATION);
  
  FillXawScrollBarCallbackStruct(w, &bar_data, XawSB_DRAG, event);
  
  DO_CALLBACK (w, XtNdragProc, (XtPointer)&bar_data);
  
  {
    float top = BAR(w).top;
    
    DO_CALLBACK (w,XtNthumbProc, *(XtPointer*)&top);
    DO_CALLBACK (w,XtNjumpProc, (XtPointer)&top);
  }

  XFlush(XtDisplay(w));
}

/* ARGSUSED */
static void Jumped( w, event, params, num_params )
   Widget    w;
   XEvent   *event;
   String   *params;		/* unused */
   Cardinal *num_params;	/* unused */
{
  Position                   x, y;
  Position                   newLoc;
  XawScrollBarCallbackStruct bar_data;
  
  if (BAR(w).scrolling)
    return;                                          /* if Select */
  
  if (LookAhead(w, event))
    return;
  
  if (!event->xmotion.same_screen)
    return;
  
  ExtractPosition( event, &x, &y );

  IS_VERTICAL(w) 
    newLoc = y - ((float)BAR(w).along * BAR(w).shown / 2.);
  else
    newLoc = x - ((float)BAR(w).along * BAR(w).shown / 2.);

  newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));
  
  if (newLoc == BAR(w).topLoc)
    return;

  BAR(w).newLoc = newLoc;

  PaintThumb(w, (Region)NULL, LOCATION);
  
  FillXawScrollBarCallbackStruct(w, &bar_data, XawSB_VALUE_CHANGED, event);
    
  DO_CALLBACK (w, XtNvalueChangedProc, (XtPointer)&bar_data);

  {
    float top = BAR(w).top;
    DO_CALLBACK (w, XtNthumbProc, *(XtPointer*)&top);
    DO_CALLBACK (w, XtNjumpProc, (XtPointer)&top);
  }

  XFlush(XtDisplay(w));
}

#if 0
static void
XawDrawFrameWindow(dpy, win, x, y, w, h, frame_type, t, lightgc,darkgc)
     Display      *dpy;
     Window        win;
     int           x, y;
     unsigned int  w, h; 
     XawFrameType  frame_type;
     int           t;
     GC            lightgc;
     GC            darkgc;
{
  XPoint topShadowPolygon[7];
  XPoint bottomShadowPolygon[7];
  
  if (t == 0) return;

#define topPolygon(i,xx,yy) \
  topShadowPolygon[i].x = (short)( xx); \
  topShadowPolygon[i].y = (short)( yy)

#define bottomPolygon(i,xx,yy) \
  bottomShadowPolygon[i].x = (short)( xx); \
  bottomShadowPolygon[i].y = (short)( yy)
  
    if( frame_type == XawRAISED || frame_type == XawSUNKEN ) {

      topPolygon(0,x    ,y    ); bottomPolygon(0,x+w  ,y+h  ); 
      topPolygon(1,x+w  ,y    ); bottomPolygon(1,x    ,y+h  );
      topPolygon(2,x+w-t,y+t  ); bottomPolygon(2,x+t  ,y+h-t);
      topPolygon(3,x+t  ,y+t  ); bottomPolygon(3,x+w-t,y+h-t);
      topPolygon(4,x+t  ,y+h-t); bottomPolygon(4,x+w-t,y+t  );
      topPolygon(5,x    ,y+h  ); bottomPolygon(5,x+w  ,y    );
      topPolygon(6,x    ,y    ); bottomPolygon(6,x+w  ,y+h  );

      if (frame_type == XawSUNKEN) {

	XFillPolygon(dpy, win, darkgc,
		     topShadowPolygon, 7, Nonconvex, CoordModeOrigin);

	XFillPolygon(dpy, win, lightgc,
		     bottomShadowPolygon, 7, Nonconvex, CoordModeOrigin);
      } else {

	XFillPolygon(dpy, win, lightgc,
		     topShadowPolygon, 7, Nonconvex, CoordModeOrigin);
	XFillPolygon(dpy, win,  darkgc,
		     bottomShadowPolygon, 7, Nonconvex, CoordModeOrigin);
      }

    } else if ( frame_type == XawLEDGED ) {

      XawDrawFrameWindow( dpy, win, x, y, w, h, XawRAISED, t/2,
			 lightgc, darkgc);

      XawDrawFrameWindow(dpy, win, (Position)(x+t/2), (Position)(y+t/2),
			 (Dimension)(w-2*(int)(t/2)),
			 (Dimension)(h-2*(int)(t/2)),
			 XawSUNKEN, t/2, lightgc, darkgc);

    } else if ( frame_type == XawCHISELED ) {

      XawDrawFrameWindow(dpy, win, x, y, w, h,
			 XawSUNKEN, t/2, lightgc, darkgc);

      XawDrawFrameWindow(dpy, win, (Position)(x+t/2),(Position)(y+t/2),
			 (Dimension)(w-2*(int)(t/2)),
			 (Dimension)(h-2*(int)(t/2)),
			 XawRAISED, t/2, lightgc, darkgc);
    }

#undef topPolygon
#undef bottomPolygon
}
#endif


/* ARGSUSED */
static void FillXawScrollBarCallbackStruct(w , s, reason, event)
     Widget w;
     XawScrollBarCallbackStruct *s;
     int reason;
     XEvent *event;
{
  s->reason      = reason;
  s->event       = event;
  s->top         = BAR(w).top;
  s->length      = BAR(w).along;
  s->shown       = BAR(w).shown;
  s->topLoc      = BAR(w).topLoc;
  s->shownLength = (Dimension) (BAR(w).along * BAR(w).shown);
  /* BAR(w).shownLength; */
}

/* ARGSUSED */
static void TimerIncrementProc (client_data, id)
    XtPointer client_data;
    XtIntervalId *id;      /* unused */
{
  Widget w = (Widget)client_data;
  XawScrollBarCallbackStruct bar_data;
  int newLoc;
  
  newLoc = BAR(w).topLoc + BAR(w).increment;
  
  newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));

  if (newLoc == BAR(w).topLoc)
    return;
    
  BAR(w).newLoc = newLoc;
    
  PaintThumb(w, (Region)NULL, LOCATION);
    
  PaintArrows(w, Xaw_DOWN, (unsigned int)Xaw_BOTTOM);
    
  FillXawScrollBarCallbackStruct(w, &bar_data, XawSB_INCREMENT, NULL);
    
  DO_CALLBACK (w, XtNincrementProc, (XtPointer)&bar_data);
  else
  DO_CALLBACK (w, XtNvalueChangedProc,(XtPointer)&bar_data);

  {
    float top = BAR(w).top;
    DO_CALLBACK (w,XtNthumbProc, *(XtPointer*)&top);
    DO_CALLBACK (w,XtNjumpProc, (XtPointer)&top);
  }
  
  XFlush(XtDisplay(w));
  BAR(w).timer = ADD_TIMEOUT(w, TimerIncrementProc);
}

/* ARGSUSED */
static void TimerDecrementProc (client_data, id)
         XtPointer client_data;
         XtIntervalId *id; /* unused */
{
  Widget w = (Widget)client_data;
  XawScrollBarCallbackStruct bar_data;
  int newLoc;
  
  newLoc = BAR(w).topLoc - BAR(w).increment;
  
  newLoc = InRange (newLoc, LOW_LIMIT(w), HIGH_LIMIT(w));
  
  if (newLoc == BAR(w).topLoc)
    return;

  BAR(w).newLoc = newLoc;
  
  PaintThumb(w, (Region)NULL, LOCATION);
  
  FillXawScrollBarCallbackStruct(w, &bar_data, XawSB_DECREMENT, (XEvent*)NULL);

  DO_CALLBACK (w, XtNdecrementProc, (XtPointer)&bar_data);
  else
  DO_CALLBACK (w, XtNvalueChangedProc, (XtPointer)&bar_data);
  
  {
    float top = BAR(w).top;
    DO_CALLBACK (w,XtNthumbProc, *(XtPointer*)&top);
    DO_CALLBACK (w,XtNjumpProc, (XtPointer)&top);
  }

  XFlush(XtDisplay(w));
  BAR(w).timer = ADD_TIMEOUT(w, TimerDecrementProc);

}

/************************************************************
 *
 *  Public routines. 
 *
 ************************************************************/

/*
 *      XawScrollbarSetThumb
 *
 * Set the scroll bar to the given location.
 */

void
  XawScrollbarSetThumb(Widget gw, double top, double shown)
{
  ScrollbarWidget w = (ScrollbarWidget)gw;

  if (BAR(w).scrolling != FREE_SCROLL)
    return;                   /* if we're already scrolling */

  top   = (top   > 1.0) ? 1.0 : (top   >= 0.0) ? top : BAR(w).top;
  shown = (shown > 1.0) ? 1.0 : (shown >= 0.0) ? shown : BAR(w).shown;

  if ( top != BAR(w).top || shown != BAR(w).shown)
  {
    BAR(w).top   = (float)top;
    BAR(w).shown = (float)shown;
    
    PaintThumb( gw, (Region)NULL, PROPORTION);
  }

}


static void PreferableSize (gw, width, height)
     Widget gw;
     Dimension *width;
     Dimension *height;
{
  ScrollbarWidget w = (ScrollbarWidget) gw;
  Dimension       arrow_size;
  int             cross = 0;
  int             along = 0;

  /*
   *    Recalculate cross if width and height of Scrollbar are set.
   */

  IS_VERTICAL(w)
  {
    if (w->core.width != 0)
      cross = w->core.width - 2 * INNER_MARGIN(w);
    else
      cross = BAR(w).cross;
  }
  else
  {
    if (w->core.height != 0)
      cross = w->core.height - 2 * INNER_MARGIN(w);
    else
      cross = BAR(w).cross;
  }

  /*
   *    Calculate size of arrow in the dependence of BAR(w).cross .
   */
  
  cross      = MAX (cross, 1);
  arrow_size = cross + 2 * MARGIN;

  /*
   *    Recalculate along if width and height of Scrollbar are set.
   */

  IS_VERTICAL(w)
  {
    if (w->core.height != 0)
      along = w->core.height - 2 * (INNER_MARGIN(w) + arrow_size); 
    else
      along = BAR(w).along;
  } else {
    if (w->core.width != 0)
      along = w->core.width - 2 * (INNER_MARGIN(w) + arrow_size);
    else
      along = BAR(w).along;
  }

  along = MAX (along, BAR(w).min_thumb);

  /*
   * Recalculate width and height of Scrollbar if they are set to zero
   * by taking into acount BAR(w).cross, BAR(w).along, arrow_size .
   */
   
  IS_VERTICAL(w) {
    *width = cross + 2 * INNER_MARGIN(w);
  }else{
    *width = along + 2 * (INNER_MARGIN(w) + arrow_size);
  }

  IS_VERTICAL(w){
    *height = along + 2 * (INNER_MARGIN(w) + arrow_size);
  }else{
    *height = cross + 2 * INNER_MARGIN(w);
  }
}

static XtGeometryResult QueryGeometry(w, intended, preferred)
    Widget w;
    XtWidgetGeometry *intended, *preferred;
{

    preferred->request_mode = CWWidth | CWHeight;

    PreferableSize (w, &preferred->width, &preferred->height);

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

