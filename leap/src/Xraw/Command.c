/*
 * Command.c - Command button widget
 */

#include <stdio.h>

#include <X11/IntrinsicP.h>
#include <X11/RectObjP.h>
#include <X11/StringDefs.h>

#include "../Xmu/Misc.h"
#include "../Xmu/Converters.h"

#include "XawInit.h"
#include "3d.h"
#include "CommandP.h"
#include "MenuButtoP.h"
#include "Container.h"

#include "XrawDebug.h"


#define DEFAULT_HIGHLIGHT_THICKNESS 2
#define DEFAULT_SHAPE_HIGHLIGHT 32767

/****************************************************************
 *
 * Full class record constant
 *
 ****************************************************************/

/* Private Data */

static char defaultTranslations[] =
#ifdef OLD_TRANSLATIONS
    "<EnterWindow>:	highlight()		\n\
     <LeaveWindow>:	reset()			\n\
     <Btn1Down>:	set()			\n\
     <Btn1Up>:		notify() unset()	";
#else
    "<Btn1Down>:	highlight() set()       \n\
     <Btn1Up>:		notify() unset() reset() ";
#endif

static void InsPixel();

#define offset(field) XtOffsetOf(CommandRec, field)
static XtResource resources[] = { 
    {
      XtNcallback, XtCCallback, XtRCallback, sizeof(XtPointer), 
      offset(command.callbacks), XtRCallback, (XtPointer)NULL
    },
    {
      XtNhighlightThickness, XtCThickness, XtRDimension, sizeof(Dimension),
      offset(command.highlight_thickness), XtRImmediate,
      (XtPointer) DEFAULT_SHAPE_HIGHLIGHT
    },
    {
      XtNshapeStyle, XtCShapeStyle, XtRShapeStyle, sizeof(int),
      offset(command.shape_style), XtRImmediate, (XtPointer)XawShapeRectangle
    },
    {
      XtNcornerRoundPercent, XtCCornerRoundPercent,
      XtRDimension, sizeof(Dimension),
      offset(command.corner_round), XtRImmediate, (XtPointer) 25
    },
    {
      XtNshadowWidth, XtCShadowWidth, XtRDimension, sizeof(Dimension),
      offset(simple.shadow_thickness), XtRImmediate, (XtPointer) 2
    },
    {
      XtNborderWidth, XtCBorderWidth, XtRDimension, sizeof(Dimension),
      XtOffsetOf(RectObjRec,rectangle.border_width), XtRImmediate,
      (XtPointer)0
    },
    {
      XtNarmedColor, XtCArmedColor, XtRPixel, sizeof(Pixel),
      offset(command.armed), XtRCallProc, (XtPointer)InsPixel
    }
};

/*ARGSUSED*/
static void InsPixel( w, off, value)
     Widget w;
     int off;
     XrmValue *value;
{
  static Pixel pixel;
  
  ArmedColor (w, w->core.background_pixel, &pixel);
  
  value->addr = (caddr_t)&pixel;
}

#undef offset

static void Initialize(), Redisplay(), Set(), Reset(), Notify(), Unset();
static void Highlight(), Unhighlight(), Destroy();
static void ClassInitialize();
static void Realize(), Resize();
static Boolean SetValues();
static Boolean ShapeButton();

extern int XShapeQueryExtension(Display*, int*, int* );
     
static XtActionsRec actionsList[] = {
  {"set",		Set},
  {"notify",		Notify},
  {"highlight",		Highlight},
  {"reset",		Reset},
  {"unset",		Unset},
  {"unhighlight",	Unhighlight}
};

XtActionList xaw_command_actions_list = actionsList;

#define SuperClass ((LabelWidgetClass)&labelClassRec)

CommandClassRec commandClassRec = {
    /* CoreClass fields initialization */
  {
    (WidgetClass) SuperClass,		/* superclass		  */	
    "Command",				/* class_name		  */
    sizeof(CommandRec),			/* size			  */
    ClassInitialize,			/* class_initialize	  */
    NULL,	                	/* class_part_initialize  */
    FALSE,				/* class_inited		  */
    Initialize,				/* initialize		  */
    NULL,				/* initialize_hook	  */
    Realize,				/* realize		  */
    actionsList,			/* actions		  */
    XtNumber(actionsList),		/* num_actions		  */
    resources,				/* resources		  */
    XtNumber(resources),		/* resource_count	  */
    NULLQUARK,				/* xrm_class		  */
    FALSE,				/* compress_motion	  */
    TRUE,				/* compress_exposure	  */
    TRUE,				/* compress_enterleave    */
    FALSE,				/* visible_interest	  */
    Destroy,				/* destroy		  */
    Resize,				/* resize		  */
    Redisplay,				/* expose		  */
    SetValues,				/* set_values		  */
    NULL,				/* set_values_hook	  */
    XtInheritSetValuesAlmost,		/* set_values_almost	  */
    NULL,				/* get_values_hook	  */
    NULL,				/* accept_focus		  */
    XtVersion,				/* version		  */
    NULL,				/* callback_private	  */
    defaultTranslations,		/* tm_table		  */
    XtInheritQueryGeometry,		/* query_geometry	  */
    XtInheritDisplayAccelerator,	/* display_accelerator	  */
    NULL				/* extension		  */
  },
    /* SimpleClass fields initialization */
  {
    XtInheritChangeSensitive,           /* change_sensitive       */
    XtInheritDisplayRectProc,		/* display_rect    	  */
    NULL                                /* extension              */
  },
    /* LabelClass fields initialization */  
  {
    XtInheritLabelPosition              /* label_position    */
  },
    /* CommandClass fields initialization */
  {
    0                                   /* field not used    */
  }
};

  /* for public consumption */
WidgetClass commandWidgetClass = (WidgetClass) &commandClassRec;

/****************************************************************
 *
 * Private Procedures
 *
 ****************************************************************/

static GC 
Get_GC(cbw, fg, bg)
CommandWidget cbw;
Pixel fg, bg;
{
  XGCValues	values;
  
  values.foreground   = fg;
  values.background	= bg;
  values.font		= cbw->label.font->fid;
  values.cap_style = CapProjecting;
  
  if (cbw->command.highlight_thickness > 1 )
    values.line_width   = cbw->command.highlight_thickness;
  else 
    values.line_width   = 0;
  
  return XtGetGC((Widget)cbw,
		 (GCForeground|GCBackground|GCFont|GCLineWidth|GCCapStyle),
		 &values);
}


/* ARGSUSED */
static void 
Initialize(request, new, args, num_args)
Widget request, new;
ArgList args;			/* unused */
Cardinal *num_args;		/* unused */
{
    CommandWidget cbw = (CommandWidget) new;
    int shape_event_base, shape_error_base;
    
    if (cbw->command.shape_style != XawShapeRectangle
	&& !XShapeQueryExtension(XtDisplay(new), &shape_event_base, 
			       &shape_error_base))
	cbw->command.shape_style = XawShapeRectangle;

    if (cbw->command.highlight_thickness == DEFAULT_SHAPE_HIGHLIGHT) {
	if (cbw->command.shape_style != XawShapeRectangle)
	    cbw->command.highlight_thickness = 0;
	else
	    cbw->command.highlight_thickness = DEFAULT_HIGHLIGHT_THICKNESS;
    }

    if (cbw->command.shape_style != XawShapeRectangle) {
	cbw->simple.shadow_thickness = 0;
	cbw->core.border_width = 1;
    }

    cbw->command.shift_GC = (GC)NULL;
    
    if (cbw->label.left_bitmap != None)
    {
      Window root;
      int x, y;
      unsigned int bw, depth;
      
      if (XGetGeometry (XtDisplay(new), cbw->label.left_bitmap, 
			&root, &x, &y,
			&cbw->label.lbm_width, &cbw->label.lbm_height,
			&bw, &depth) 
	  && depth == 1 )
	cbw->command.shift_GC = Get_GC(cbw, cbw->simple.foreground,
				       cbw->command.armed);
    }
    
    cbw->command.armed_GC = Get_GC(cbw, cbw->command.armed,
				   cbw->core.background_pixel);
    cbw->command.normal_GC = Get_GC(cbw, cbw->simple.foreground, 
				  cbw->core.background_pixel);
    cbw->command.inverse_GC = Get_GC(cbw, cbw->core.background_pixel, 
				   cbw->simple.foreground);
    XtReleaseGC(new, cbw->label.normal_GC);
    cbw->label.normal_GC = cbw->command.normal_GC;

    cbw->command.set               = FALSE;
    cbw->command.highlighted       = HighlightNone;
    cbw->command.changed_highlight = FALSE;
    cbw->command.changed_set       = FALSE;
}

/***************************
*
*  Action Procedures
*
***************************/

/* ARGSUSED */
static void 
Set(w,event,params,num_params)
Widget w;
XEvent *event;
String *params;		/* unused */
Cardinal *num_params;	/* unused */
{
    CommandWidget cbw = (CommandWidget)w;

    if (cbw->command.set)
	return;

    cbw->command.set= TRUE;
    cbw->command.changed_set = True;
    
    (*XtClass(w)->core_class.expose)(w, event, (Region) NULL); 
}

/* ARGSUSED */
static void Unset(w,event,params,num_params)
     Widget w;
     XEvent *event;
     String *params;		/* unused */
     Cardinal *num_params;
{
    CommandWidget cbw = (CommandWidget)w;

    if (!cbw->command.set)
	return;

    cbw->command.set = FALSE;
    cbw->command.changed_set = True;
    
    (*XtClass(w)->core_class.expose)(w, event, (Region) NULL); 
}

/* ARGSUSED */
static void Reset(w,event,params,num_params)
     Widget w;
     XEvent *event;
     String *params;		/* unused */
     Cardinal *num_params;   /* unused */
{
    CommandWidget cbw = (CommandWidget)w;

    if (cbw->command.set) {
      cbw->command.highlighted = HighlightNone;
	Unset(w, event, params, num_params);
    } else
	Unhighlight(w, event, params, num_params);
}

/* ARGSUSED */
static void Highlight(w,event,params,num_params)
     Widget w;
     XEvent *event;
     String *params;		
     Cardinal *num_params;	
{
  CommandWidget cbw = (CommandWidget)w;
  Boolean       highlighted = cbw->command.highlighted;
  
  
  if ( *num_params == (Cardinal) 0) 
    cbw->command.highlighted = HighlightWhenUnset;
  else {
    if ( *num_params != (Cardinal) 1) 
      XtWarning("Too many parameters passed to highlight action table.");
    switch (params[0][0]) {
    case 'A':
    case 'a':
      cbw->command.highlighted = HighlightAlways;
      break;
    default:
      cbw->command.highlighted = HighlightWhenUnset;
      break;
    }
  }
  
  cbw->command.changed_highlight = highlighted != cbw->command.highlighted;
  
  (*XtClass(w)->core_class.expose)(w, event, (Region) NULL);
  
}

/* ARGSUSED */
static void Unhighlight(w,event,params,num_params)
     Widget w;
     XEvent *event;
     String *params;		/* unused */
     Cardinal *num_params;	/* unused */
{
    CommandWidget cbw = (CommandWidget)w;

    cbw->command.changed_highlight = cbw->command.highlighted != HighlightNone;
    cbw->command.highlighted = HighlightNone;

    (*XtClass(w)->core_class.expose)(w, event, (Region) NULL);
    
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
static void 
Notify(w,event,params,num_params)
Widget w;
XEvent *event;
String *params;		/* unused */
Cardinal *num_params;	/* unused */
{
  CommandWidget cbw = (CommandWidget)w; 
  Position   x, y;
  XRectangle rectangle;
  Boolean    rectangle_exist;
  Boolean    point_is_in_box;
  
  ExtractPosition (event, &x, &y);
  
  rectangle_exist = (*SuperClass->simple_class.display_rect) (w, &rectangle);
  
  point_is_in_box  = x >= rectangle.x && x < (rectangle.x + rectangle.width);
  point_is_in_box &= y >= rectangle.y && y < (rectangle.y + rectangle.height);

  /* check to be sure state is still Set so that user can cancel
     the action (e.g. by moving outside the window, in the default
     bindings.
  */

  if (cbw->command.set && (!rectangle_exist || point_is_in_box))
    XtCallCallbackList(w, cbw->command.callbacks, NULL);
}

/*
 * Repaint the widget window
 */

/************************
*
*  REDISPLAY (DRAW)
*
************************/

/* ARGSUSED */
static void Redisplay(gw, event, region)
     Widget gw;
     XEvent *event;
     Region region;
{
  CommandWidget w = (CommandWidget) gw;
  XRectangle rectangle;
    
  if (!XtIsRealized(gw))
    return;

  (*XtLabelClass(w)->simple_class.display_rect)(gw, &rectangle);
  
  if ( !XtIsSubclass(gw, menuButtonWidgetClass) &&
      (w->command.changed_highlight || w->command.changed_set))
  {
    if (w->command.highlighted != HighlightNone)
    {
      XFillRectangle(XtDisplay(w), XtWindow(w), w->command.armed_GC,
		     rectangle.x , rectangle.y,
		     rectangle.width, rectangle.height);
    }
    else
    {
      XClearArea(XtDisplay(gw), XtWindow(gw),
		 rectangle.x , rectangle.y,
		 rectangle.width, rectangle.height, False);
    }
    w->command.changed_highlight = False;
    w->command.changed_set       = False;
  }
  
  if (w->command.set && region != (Region)NULL)
  {
    XFillRectangle(XtDisplay(w), XtWindow(w), w->command.armed_GC,
		   rectangle.x , rectangle.y,
		   rectangle.width, rectangle.height);
  }

#ifndef COMMAND_FACE_TACK  
  if (w->command.set) { 
    w->label.label_x += 1; w->label.label_y += 1; 
    
  }
#endif

  (*SuperClass->core_class.expose)(gw, event, region);

#ifndef COMMAND_FACE_TACK
  if (w->command.set) { w->label.label_x -= 1; w->label.label_y -= 1; }
#endif

#ifdef COMMAND_FACE_TACK  
  XawDrawFrame (gw,
		w->simple.highlight_thickness,
		w->simple.highlight_thickness,
		(Dimension)(w->core.width - 2*w->simple.highlight_thickness),
		(Dimension)(w->core.height - 2*w->simple.highlight_thickness),
		/*  XawSUNKEN */ XawTACK,
		w->simple.shadow_thickness,
		w->simple.top_shadow_GC,
		w->simple.bottom_shadow_GC);
#endif

  if (w->command.set)
  {
#ifdef COMMAND_FACE_TACK
    XDrawRectangle (XtDisplay(gw), XtWindow(gw), w->simple.bottom_shadow_GC,
		    w->simple.highlight_thickness,
		    w->simple.highlight_thickness,
		    (w->core.width - 2*w->simple.highlight_thickness) - 1,
		    (w->core.height - 2*w->simple.highlight_thickness) - 1 );
#else
    XawDrawFrame (gw,
		  (Position)w->simple.highlight_thickness,
		  (Position)w->simple.highlight_thickness,
	     (Dimension)(w->core.width - 2*w->simple.highlight_thickness),
	     (Dimension)(w->core.height - 2*w->simple.highlight_thickness),
		  XawSUNKEN,
		  w->simple.shadow_thickness,
		  w->simple.top_shadow_GC,
		  w->simple.bottom_shadow_GC);
#endif
  }
  
}

static void 
Destroy(w)
Widget w;
{
    CommandWidget cbw = (CommandWidget) w;

    /* so Label can release it */
    if (cbw->label.normal_GC == cbw->command.normal_GC)
	XtReleaseGC( w, cbw->command.inverse_GC );
    else
	XtReleaseGC( w, cbw->command.normal_GC );
}

/*
 * Set specified arguments into widget
 */

/* ARGSUSED */
static Boolean SetValues (current, request, new, args, num_args)
     Widget current, request, new;
     ArgList args;
     Cardinal *num_args;
{
  CommandWidget oldcbw = (CommandWidget) current;
  CommandWidget cbw = (CommandWidget) new;
  Boolean redisplay = False;

  if ( oldcbw->core.sensitive != cbw->core.sensitive && !cbw->core.sensitive)
  {
    /* about to become insensitive */
    cbw->command.set = FALSE;
    cbw->command.highlighted = HighlightNone;
    redisplay = TRUE;
  }

#define NE(field) (oldcbw->command.field != cbw->command.field)
  
  if (oldcbw->label.left_bitmap != cbw->label.left_bitmap) 
  {
    if (cbw->command.shift_GC != (GC)NULL) {
      XtReleaseGC (current, cbw->command.shift_GC);
      cbw->command.shift_GC = (GC)NULL;
    }
    
    if (cbw->label.left_bitmap != None)
    {
      Window root;
      int x, y;
      unsigned int bw, depth;
      
      if (XGetGeometry (XtDisplay(new), cbw->label.left_bitmap, 
			&root, &x, &y,
			&cbw->label.lbm_width, &cbw->label.lbm_height,
			&bw, &depth) 
	  && depth == 1 )
	cbw->command.shift_GC = 
	  Get_GC(cbw, cbw->simple.foreground, cbw->command.armed);
    }
  }
  else  if (cbw->command.shift_GC != (GC)NULL && 
      (NE(armed) || (oldcbw->simple.foreground != cbw->simple.foreground)))
  {
    XtReleaseGC (current, cbw->command.shift_GC);
    cbw->command.shift_GC = 
      Get_GC(cbw, cbw->simple.foreground, cbw->command.armed);
  }

  if (NE(armed)||(oldcbw->core.background_pixel != cbw->core.background_pixel))
  {
    XtReleaseGC (current, cbw->command.armed_GC);
    cbw->command.armed_GC = 
      Get_GC(cbw, cbw->command.armed, cbw->core.background_pixel);
  }

  if ( (oldcbw->simple.foreground != cbw->simple.foreground)           ||
       (oldcbw->core.background_pixel != cbw->core.background_pixel) ||
       (oldcbw->label.font != cbw->label.font) ) 
  {
    if (oldcbw->label.normal_GC == oldcbw->command.normal_GC)
	/* Label has release one of these */
      XtReleaseGC(new, cbw->command.inverse_GC);
    else
      XtReleaseGC(new, cbw->command.normal_GC);

    cbw->command.normal_GC = Get_GC(cbw, cbw->simple.foreground, 
				    cbw->core.background_pixel);
    cbw->command.inverse_GC = Get_GC(cbw, cbw->core.background_pixel, 
				     cbw->simple.foreground);
    XtReleaseGC(new, cbw->label.normal_GC);

    cbw->label.normal_GC = cbw->command.normal_GC;
/* 
    cbw->label.normal_GC = (cbw->command.set
			    ? cbw->command.inverse_GC
			    : cbw->command.normal_GC);
*/    
    redisplay = True;
  }

  if ( XtIsRealized(new)
       && oldcbw->command.shape_style != cbw->command.shape_style
       && !ShapeButton(cbw, TRUE))
  {
      cbw->command.shape_style = oldcbw->command.shape_style;
  }

  return (redisplay);
#undef NE
}

static void ClassInitialize()
{
  XawInitializeWidgetSet();
  XtSetTypeConverter( XtRString, XtRShapeStyle, XmuCvtStringToShapeStyle,
		     NULL, 0, XtCacheNone, NULL );

  commandWidgetClass->core_class.query_geometry =
    labelWidgetClass->core_class.query_geometry;
}


static Boolean
ShapeButton(cbw, checkRectangular)
CommandWidget cbw;
Boolean checkRectangular;
{
    Dimension corner_size;

    if ( (cbw->command.shape_style == XawShapeRoundedRectangle) ) {
	corner_size = (cbw->core.width < cbw->core.height) ? cbw->core.width 
	                                                   : cbw->core.height;
	corner_size = (int) (corner_size * cbw->command.corner_round) / 100;
    }

    if (checkRectangular || cbw->command.shape_style != XawShapeRectangle) {
	if (!XmuReshapeWidget((Widget) cbw, cbw->command.shape_style,
			      corner_size, corner_size)) {
	    cbw->command.shape_style = XawShapeRectangle;
	    return(False);
	}
    }
    return(TRUE);
}

static void Realize(w, valueMask, attributes)
    Widget w;
    Mask *valueMask;
    XSetWindowAttributes *attributes;
{
    (*SuperClass->core_class.realize) (w, valueMask, attributes);

    (void) ShapeButton( (CommandWidget) w, FALSE);
}

static void Resize(w)
    Widget w;
{
    if (XtIsRealized(w)) 
	(void) ShapeButton( (CommandWidget) w, FALSE);

    (*SuperClass->core_class.resize)(w);
}


void
XawSetCommandColours (Widget widget, Pixel background)
{     
  Pixel top_shadow;
  Pixel bottom_shadow;
  Pixel armed;
  Arg   args[4];
  Cardinal num_args = 0;

  /* 
   * Simple and Container Widget Class have shadow resources
   */
  if (!XtIsSubclass(widget, simpleWidgetClass) && 
      !XtIsSubclass(widget, containerWidgetClass))
    return;

  if (!TopShadowColor (widget, background, &top_shadow))
    top_shadow = background;

  if (!BottomShadowColor (widget, background, &bottom_shadow))
    bottom_shadow = background;

  XtSetArg(args[num_args], XtNtopShadowPixel,    top_shadow);    num_args++;
  XtSetArg(args[num_args], XtNbottomShadowPixel, bottom_shadow); num_args++;
  XtSetArg(args[num_args], XtNbackground,        background);    num_args++;

  if (XtIsSubclass(widget, commandWidgetClass))
  {
    if (!ArmedColor (widget, background, &armed))
      armed = background;

    XtSetArg(args[num_args], XtNarmedColor, armed); num_args++;
  }

  XtSetValues(widget, args, num_args);
}
