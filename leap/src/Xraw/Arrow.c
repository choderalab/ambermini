/* 
 * Xraw: Arrow.c
 *
 * If You would like to use accurate method to culculate math expressions,
 * use USE_MATH_LIB definition. 
 * #define USE_MATH_LIB 
 */

#include <math.h>
#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>

#include "XawInit.h"
#include "ArrowP.h"

#include "XrawDebug.h"
#include "../Xmu/CharSet.h"

#define CORE(w) (w)->core

#define MARGIN    1

#define InnerMargin(new) (SIMPLE_MARGIN(new) + MARGIN)

#define offset(field) XtOffsetOf(ArrowRec, arrow.field)

static XtResource resources[] = {
  {
    XtNdirection, XtCDirection, XtRDirection, sizeof (XawDirection),
    offset(outline.direction), XtRImmediate, (XtPointer) XawTop
  },
  {
    XtNforeground, XtCForeground, XtRPixel, sizeof (Pixel),
    XtOffsetOf(ArrowRec,simple.foreground), XtRString, 
    (XtPointer) XtDefaultBackground
  },
  {
    XtNarrowShadow, XtCArrowShadow, XtRDimension, sizeof (Dimension),
    offset(arrowShadow), XtRImmediate, (XtPointer) 2
  },
  {
    XtNshadowWidth, XtCShadowWidth, XtRDimension, sizeof(Dimension),
    XtOffsetOf(ArrowRec,simple.shadow_thickness), XtRImmediate, (XtPointer) 2
  }
};
#undef offset

static void push_up();
static void push_down();

static XtActionsRec actions[] =  {
  {"push_up",   push_up},
  {"push_down", push_down},
};

static char defaultTranslations[] =
  "<Btn1Down>: push_down() start() \n\
   <Btn1Up>:   push_up()   stop() ";

static void ClassInitialize();
static void Initialize ();
static Boolean SetValues ();
static void Redisplay ();
static void Resize ();

#define SuperClass (RepeaterWidgetClass)&repeaterClassRec
ArrowClassRec arrowClassRec ={
  { /* core fields */
    /* superclass		*/    (WidgetClass) SuperClass,
    /* class_name		*/    "Arrow",
    /* widget_size		*/    sizeof (ArrowRec),
    /* class_initialise         */    ClassInitialize,
    /* class_part_initialize	*/    NULL,
    /* class_inited		*/    FALSE,
    /* initialize		*/    Initialize,
    /* initialize_hook		*/    NULL,
    /* realize			*/    XtInheritRealize,
    /* actions			*/    actions,
    /* num_actions		*/    XtNumber(actions),
    /* resources		*/    resources,
    /* num_resources		*/    XtNumber(resources),
    /* xrm_class		*/    NULLQUARK,
    /* compress_motion		*/    TRUE,
    /* compress_exposure	*/    TRUE,
    /* compress_enterleave	*/    TRUE,
    /* visible_interest		*/    FALSE,
    /* destroy			*/    NULL,
    /* resize			*/    Resize,
    /* expose			*/    Redisplay,
    /* set_values		*/    SetValues,
    /* set_values_hook		*/    NULL,
    /* set_values_almost	*/    XtInheritSetValuesAlmost,
    /* get_values_hook		*/    NULL,
    /* accept_focus		*/    XtInheritAcceptFocus,
    /* version			*/    XtVersion,
    /* callback_private		*/    NULL,
    /* tm_table			*/    defaultTranslations,
    /* query_geometry		*/    XtInheritQueryGeometry,
    /* display_accelerator	*/    XtInheritDisplayAccelerator,
    /* extension		*/    NULL
  },
  { /* simple fields */
    /* change_sensitive		*/    XtInheritChangeSensitive,
    /* display_rect		*/    XtInheritDisplayRectProc,
    /* extension		*/    NULL				
  },
  { /* label fields */
    /* label_position		*/    XtInheritLabelPosition,
  },
  { /* command fields */
    /* ignore			*/    0
  },
  { /* repeater fields */
    /* ignore                   */    0
  },
  { /* arrow    fields */
    /* ignore                   */    0
  },
};

WidgetClass arrowWidgetClass = (WidgetClass) & arrowClassRec;

#define done(type, value)                    \
      if (to->addr != NULL) {                \
	  if (to->size < sizeof(type)) {     \
	      to->size = sizeof(type);       \
	      return False;                  \
	  }                                  \
	  *(type*)(to->addr) = (value);      \
      } else {                               \
	  static type static_val;            \
	  static_val = (value);              \
	  to->addr = (XtPointer)&static_val; \
      }                                      \
      to->size = sizeof(type);               \
      return True

/* ARGSUSED */
static  Boolean
cvtStringToDirection ( display, args, num_args,  from, to, converter_data)
     Display *display;
     XrmValuePtr args;
     Cardinal *num_args;
     XrmValuePtr from;
     XrmValuePtr to;
     XtPointer *converter_data;
{
  String s = (String) from->addr;
  
  if (*num_args != 0)
    XtAppErrorMsg(XtDisplayToApplicationContext(display),
		  "cvtStringToDirection", "wrongParameters",
		  "XtToolkitError",
		  "String to XawDirection conversion needs no arguments",
		  (String*) NULL, (Cardinal*) NULL);

#define IS_EQUAL(resource) (XmuCompareISOLatin1(s, resource) == 0)  

  if      (IS_EQUAL("top")   ) { done(XawDirection, XawTop);    }
  else if (IS_EQUAL("bottom")) { done(XawDirection, XawBottom); }
  else if (IS_EQUAL("left")  ) { done(XawDirection, XawLeft);   }
  else if (IS_EQUAL("right") ) { done(XawDirection, XawRight);  }

  XtDisplayStringConversionWarning(display, s, XtRDirection);
  done(XawDirection, XawTop);

#undef IS_EQUAL
}

static void ClassInitialize()
{
  XawInitializeWidgetSet();

  XtSetTypeConverter(XtRString, XtRDirection, cvtStringToDirection,
		     (XtConvertArgList)NULL, 0, XtCacheNone, NULL);
}

static void CreateArrowGC (new)
     Widget new;
{
  XtGCMask mask;
  XGCValues values;

  if (ARROW(new).arrowgc != NULL)
    XtReleaseGC (new, ARROW(new).arrowgc);
  mask = GCForeground;
  values.foreground = ((ArrowWidget)new)->simple.foreground;
  ARROW(new).arrowgc = XtGetGC (new, mask, &values);
}

/*ARGSUSED */
static void push_up (new, event, params, num_params)
     Widget new;
     XEvent *event;
     String *params;
     Cardinal *num_params;
{

  switch (ARROW(new).outline.direction) {
  case XawTop:
    XFillPolygon (XtDisplay (new), XtWindow (new),
		  ARROW(new).arrowgc, ARROW(new).outline.p1, 3, Convex,
		  CoordModeOrigin);

    if (ARROW(new).arrowShadow <= 0)
      break;

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p2, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p3, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p4, 
		  4, Convex, CoordModeOrigin);
    break;
  case XawLeft:

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  ARROW(new).arrowgc, ARROW(new).outline.p1, 3, Convex,
		  CoordModeOrigin);

    if (ARROW(new).arrowShadow <= 0)
      break;

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p3, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p4, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p2, 
		  4, Convex, CoordModeOrigin);

    break;
  case XawBottom:

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  ARROW(new).arrowgc, ARROW(new).outline.p1, 3, Convex,
		  CoordModeOrigin);

    if (ARROW(new).arrowShadow <= 0)
      break;

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p2, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p4, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p3, 
		  4, Convex, CoordModeOrigin);

    break;
  case XawRight:

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  ARROW(new).arrowgc, ARROW(new).outline.p1, 
		  3, Convex, CoordModeOrigin);

    if (ARROW(new).arrowShadow <= 0)
      break;

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p3, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p4, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p2, 
		  4, Convex, CoordModeOrigin);

    break;
  }
}

/*ARGSUSED */
static void push_down (new, event, params, num_params)
     Widget new;
     XEvent *event;
     String *params;
     Cardinal *num_params;
{

  switch (ARROW(new).outline.direction) {
  case XawTop:
    XFillPolygon (XtDisplay (new), XtWindow (new),
		  ARROW(new).arrowgc, ARROW(new).outline.p1, 3, Convex,
		  CoordModeOrigin);

    if (ARROW(new).arrowShadow <= 0)
      break;

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p2, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p3, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p4, 
		  4, Convex, CoordModeOrigin);
    break;
  case XawLeft:

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  ARROW(new).arrowgc, ARROW(new).outline.p1, 3, Convex,
		  CoordModeOrigin);

    if (ARROW(new).arrowShadow <= 0)
      break;

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p3, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p4, 4, Convex,
		  CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p2, 
		  4, Convex, CoordModeOrigin);

    break;
  case XawBottom:

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  ARROW(new).arrowgc, ARROW(new).outline.p1, 3, Convex,
		  CoordModeOrigin);

    if (ARROW(new).arrowShadow <= 0)
      break;

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p2, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p4, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p3, 4, Convex,
		  CoordModeOrigin);
    break;
  case XawRight:

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  ARROW(new).arrowgc, ARROW(new).outline.p1, 3, Convex,
		  CoordModeOrigin);

    if (ARROW(new).arrowShadow <= 0)
      break;

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p3, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).bottom_shadow_GC, ARROW(new).outline.p4, 
		  4, Convex, CoordModeOrigin);

    XFillPolygon (XtDisplay (new), XtWindow (new),
		  SIMPLE(new).top_shadow_GC, ARROW(new).outline.p2, 4, Convex,
		  CoordModeOrigin);

    break;
  }
}

static void Set_a2_a3 (new)
     Widget new;
{
  Dimension a = ARROW(new).arrowShadow;
  Dimension a2, a3;
  int width, height;

  width  = CORE(new).width  - 2 * InnerMargin(new);
  height = CORE(new).height - 2 * InnerMargin(new);

#define f(x) ((float)x)
#if USE_MATH_LIB  
  {
    float alfa;
    alfa = atan( (double)(f(f(width)/2. ) / f(height)));
    a2   = (Dimension)(f(a)/cos(alfa) + f(a)*tan(alfa));
    a3   = (Dimension)(f(a)/sin(alfa));
  }
#else
  a2   = (Dimension)((1.0 + 0.71 * f(width) / f(height)) * f(a));
  a3   = (Dimension)((1.0 + 0.83 * f(height) / f(width)) * f(a));
#endif  
#undef f

  switch (ARROW(new).outline.direction) {
  case XawTop:
  case XawBottom:
    ARROW(new).outline.a2 = a2;
    ARROW(new).outline.a3 = a3;
    break;
  case XawLeft:
  case XawRight:
    ARROW(new).outline.a2 = a3;
    ARROW(new).outline.a3 = a2;
    break;
  }

  ARROW(new).set_a2_a3 = TRUE;

}

/*ARGSUSED */ 
static void Initialize (request, new, args, num_args)
     Widget request;
     Widget new;
     ArgList args;
     Cardinal *num_args;
{
  if (ARROW(new).outline.direction != XawTop &&
      ARROW(new).outline.direction != XawLeft &&
      ARROW(new).outline.direction != XawRight &&
      ARROW(new).outline.direction != XawBottom) {
    XtWarning ("Direction of Arrow widget incorrect; So hi's set to ``top''");
    ARROW(new).outline.direction = XawTop;
  }

  ARROW(new).set_a2_a3 = FALSE;
  ARROW(new).arrowgc = NULL;
  CreateArrowGC (new);
}

/*ARGSUSED */ 
static Boolean SetValues (old, request, new, args, num_args)
     Widget old;
     Widget request;
     Widget new;
     ArgList args;
     Cardinal *num_args;
{
  Boolean need_redisplay = False;

  if (ARROW(new).outline.direction != XawTop &&
      ARROW(new).outline.direction != XawLeft &&
      ARROW(new).outline.direction != XawRight &&
      ARROW(new).outline.direction != XawBottom) {
    XtWarning ("Direction of Arrow widget incorrect; So hi has set to `top'");
    ARROW(new).outline.direction = XawTop;
  }
#define NE(field) (ARROW(new).field != ARROW(old).field)

  if (NE(outline.direction))
    need_redisplay = True;

  if (((ArrowWidget)new)->simple.foreground != 
      ((ArrowWidget)old)->simple.foreground)
  {
    CreateArrowGC (new);
    need_redisplay = True;
  }

  if (NE(arrowShadow))
  {
    Set_a2_a3(new);
    need_redisplay = True;
  }
      
  return need_redisplay;
#undef NE
}


static void Redisplay (new, event, region)
     Widget new;
     XEvent *event;
     Region region;
{
  if (!XtIsRealized (new))
    return;

  (*simpleWidgetClass->core_class.expose) (new, event, region);
  
  XawDrawArrow(new,
	       ARROW(new).arrowgc,
	       SIMPLE(new).top_shadow_GC,
	       SIMPLE(new).bottom_shadow_GC,
	       (XawDrawArrowStruct*)&ARROW(new).outline);
}


static void Resize( w )
     Widget w;
{
  Position x, y;
  Dimension width, height;

  Set_a2_a3(w);

  if (ARROW(w).outline.direction != XawTop &&
      ARROW(w).outline.direction != XawLeft &&
      ARROW(w).outline.direction != XawRight &&
      ARROW(w).outline.direction != XawBottom) 
  {
    XtWarning ("Direction of Arrow widget incorrect; So hi's set to ``top''");
    ARROW(w).outline.direction = XawTop;
  }

  x = y  = InnerMargin(w);
  width  = CORE(w).width  - 2 * InnerMargin(w);
  height = CORE(w).height - 2 * InnerMargin(w);

  XawMakeDrawArrowStruct( x, y, width, height, ARROW(w).arrowShadow,
			 ARROW(w).outline.direction,
			 (XawDrawArrowStruct*)&ARROW(w).outline);

  Redisplay (w, (XEvent*)NULL, (Region)NULL );
}


void
  XawMakeDrawArrowStruct(int x,
			 int y,
			 unsigned int w,
			 unsigned int h,
			 unsigned int thickness,
			 XawDirection direction,
			 XawDrawArrowStruct *result)
{
  int a, a2, a3;
  
  a = thickness;

  y --;  /* that is magic */
  h ++;  /* that is magic */
  
#define f(x) ((float)x)
  if ( direction == XawTop || direction == XawBottom) {
#if USE_MATH_LIB  
    float alfa;
    alfa = atan( (double)(f(f(w)/2. ) / f(h)));
    a2   = (int)(f(a)/cos(alfa) + f(a)*tan(alfa));
    a3   = (int)(f(a)/sin(alfa));
#else
    a2   = (int)((1.0 + 0.71 * f(w) / f(h)) * f(a));
    a3   = (int)((1.0 + 0.83 * f(h) / f(w)) * f(a));
#endif  
  } else {
#if USE_MATH_LIB  
    float alfa;
    alfa = atan( (double)(f(f(h)/2. ) / f(w)));
    a3   = (int)(f(a)/cos(alfa) + f(a)*tan(alfa));
    a2   = (int)(f(a)/sin(alfa));
#else
    a3   = (int)((1.0 + 0.71 * f(w) / f(h)) * f(a));
    a2   = (int)((1.0 + 0.83 * f(h) / f(w)) * f(a));
#endif  
  }
#undef f

  result->a         = a;
  result->a2        = a2;
  result->a3        = a3;
  result->direction = direction;
  
#define PUT_COORD(array,n,xx,yy) result->array[n].x=(short)(xx);\
                                 result->array[n].y=(short)(yy)

  switch (direction) {
  case XawTop:

    PUT_COORD(p1,0,x + w / 2 , y + a3);
    PUT_COORD(p1,1,x + a2    , y + h - a);
    PUT_COORD(p1,2,x + w - a2, y + h - a);

    if (a <= 0)
      break;

    PUT_COORD(p2,0,x + w / 2 , y);
    PUT_COORD(p2,1,x + w / 2 , y + a3);
    PUT_COORD(p2,2,x + w - a2, y + h - a);
    PUT_COORD(p2,3,x + w     , y + h);

    PUT_COORD(p3,0,x + a2    , y + h - a);
    PUT_COORD(p3,1,x         , y + h);
    PUT_COORD(p3,2,x + w     , y + h);
    PUT_COORD(p3,3,x + w - a2, y + h - a);

    PUT_COORD(p4,0,x + w / 2 , y);
    PUT_COORD(p4,1,x         , y + h);
    PUT_COORD(p4,2,x + a2    , y + h - a);
    PUT_COORD(p4,3,x + w / 2 , y + a3);

    break;
  case XawLeft:

    PUT_COORD(p1,0,x + a2        , y + h / 2);
    PUT_COORD(p1,1,x + w - a , y + a3);
    PUT_COORD(p1,2,x + w - a , y + h - a3);

    if (a <= 0)
      break;

    PUT_COORD(p2,0,x + w     , y);
    PUT_COORD(p2,1,x         , y + h / 2);
    PUT_COORD(p2,2,x + a2    , y + h / 2);
    PUT_COORD(p2,3,x + w - a , y + a3);

    PUT_COORD(p3,0,x         , y + h / 2);
    PUT_COORD(p3,1,x + w     , y + h);
    PUT_COORD(p3,2,x + w - a , y + h - a3);
    PUT_COORD(p3,3,x + a2        , y + h / 2);

    PUT_COORD(p4,0,x + w     , y);
    PUT_COORD(p4,1,x + w - a , y + a3);
    PUT_COORD(p4,2,x + w - a , y + h - a3);
    PUT_COORD(p4,3,x + w     , y + h);

    break;
  case XawBottom:

    y ++;  /* that is magic */
    h --;  /* that is magic */
  
    PUT_COORD(p1,0,x + w / 2 , y + h - a3);
    PUT_COORD(p1,1,x + a2    , y + a);
    PUT_COORD(p1,2,x + w - a2, y + a);

    if (a <= 0)
      break;

    PUT_COORD(p2,0,x         , y);
    PUT_COORD(p2,1,x + w / 2 , y + h);
    PUT_COORD(p2,2,x + w / 2 , y + h - a3);
    PUT_COORD(p2,3,x + a2    , y + a);

    PUT_COORD(p3,0,x + w     , y);
    PUT_COORD(p3,1,x + w - a2, y + a);
    PUT_COORD(p3,2,x + w / 2 , y + h - a3);
    PUT_COORD(p3,3,x + w / 2 , y + h);

    PUT_COORD(p4,0,x         , y);
    PUT_COORD(p4,1,x + a2    , y + a);
    PUT_COORD(p4,2,x + w - a2, y + a);
    PUT_COORD(p4,3,x + w     , y);

    break;
  case XawRight:

    PUT_COORD(p1,0,x + w - a2, y + h / 2);
    PUT_COORD(p1,1,x + a     , y + a3);
    PUT_COORD(p1,2,x + a     , y + h - a3);

    if (a <= 0)
      break;

    PUT_COORD(p2,0,x         , y + h);
    PUT_COORD(p2,1,x + w     , y + h / 2);
    PUT_COORD(p2,2,x + w - a2, y + h / 2);
    PUT_COORD(p2,3,x + a     , y + h - a3);

    PUT_COORD(p3,0,x         , y);
    PUT_COORD(p3,1,x + a     , y + a3);
    PUT_COORD(p3,2,x + w - a2, y + h / 2);
    PUT_COORD(p3,3,x + w     , y + h / 2);

    PUT_COORD(p4,0,x         , y);
    PUT_COORD(p4,1,x         , y + h);
    PUT_COORD(p4,2,x + a     , y + h - a3);
    PUT_COORD(p4,3,x + a     , y + a3);

    break;
  }
#undef PUT_COORD
}


void
XawDrawArrow (Widget w,
	      GC inner_gc,
	      GC top_gc,
	      GC bottom_gc,
	      XawDrawArrowStruct *draw_struct)
{
  if (!XtIsRealized (w))
    return;

  switch (draw_struct->direction) {
  case XawTop:
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  inner_gc, draw_struct->p1, 3, Convex,
		  CoordModeOrigin);

    if (draw_struct->a <= 0)
      break;

    XFillPolygon (XtDisplay (w), XtWindow (w),
		  bottom_gc, draw_struct->p2, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  bottom_gc, draw_struct->p3, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  top_gc, draw_struct->p4, 4, Convex,
		  CoordModeOrigin);
    break;

  case XawLeft:
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  inner_gc, draw_struct->p1, 3, Convex,
		  CoordModeOrigin);

    if (draw_struct->a <= 0)
      break;

    XFillPolygon (XtDisplay (w), XtWindow (w),
		  bottom_gc, draw_struct->p3, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  bottom_gc, draw_struct->p4, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  top_gc, draw_struct->p2, 4, Convex,
		  CoordModeOrigin);
    break;

  case XawBottom:
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  inner_gc, draw_struct->p1, 3, Convex,
		  CoordModeOrigin);

    if (draw_struct->a <= 0)
      break;

    XFillPolygon (XtDisplay (w), XtWindow (w),
		  top_gc, draw_struct->p2, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  top_gc, draw_struct->p4, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  bottom_gc, draw_struct->p3, 4, Convex,
		  CoordModeOrigin);
    break;

  case XawRight:
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  inner_gc, draw_struct->p1, 3, Convex,
		  CoordModeOrigin);

    if (draw_struct->a <= 0)
      break;

    XFillPolygon (XtDisplay (w), XtWindow (w),
		  top_gc, draw_struct->p3, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  top_gc, draw_struct->p4, 4, Convex,
		  CoordModeOrigin);
    XFillPolygon (XtDisplay (w), XtWindow (w),
		  bottom_gc, draw_struct->p2, 4, Convex,
		  CoordModeOrigin);
    break;
  }
}

