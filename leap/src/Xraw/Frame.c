/* 
 * Frame.c - Frame composite widget
 * 
 */

#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>

#include "../Xmu/Converters.h"
#include "../Xmu/CharSet.h"

#include "XawInit.h"
#include "Label.h"
#include "FrameP.h"

#include "XrawDebug.h"


#define SHADOW(w)    (CONTAINER(w).shadow_thickness)

#define X_MARGIN(w)  (CONTAINER(w).shadow_thickness + FRAME(w).h_space)
#define Y_MARGIN(w)  (CONTAINER(w).shadow_thickness + FRAME(w).v_space)

#define MIN(a,b) ((a)<(b)?(a):(b))

#define CAPTION_JUSTIFY_SIDE              10
#define CAPTION_JUSTIFY_OPPOSITE_SIDE     4

#define CAPTION_MARGIN(sw)                \
  (FRAME(sw).justify == XtJustifyCenter ? \
   2 * CAPTION_JUSTIFY_SIDE :             \
   CAPTION_JUSTIFY_SIDE + CAPTION_JUSTIFY_OPPOSITE_SIDE)

#define LINE_WIDTH 2

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifdef CRAY
#define WORD64
#endif

#ifndef WORD64

#define TXT16 XChar2b

#else

#define TXT16 char

static XChar2b *buf2b;
static int buf2blen = 0;

#define XTextWidth16 _XawLabelWidth16
#define XDrawString16 _XawLabelDraw16

#endif /* WORD64 */


/****************************************************************
 *
 * Frame Resources
 *
 ****************************************************************/

#define Offset(field) XtOffsetOf(FrameConstraintsRec, frame.field)

static XtResource frameConstraintResources[] = {
    {
      XtNtop, XtCFraction, XtRInt, sizeof(int),
      Offset(top), XtRImmediate, (XtPointer) 0
    },
    {
      XtNbottom, XtCFraction, XtRInt, sizeof(int),
      Offset(bottom), XtRImmediate, (XtPointer) 100
    },
    {
      XtNleft, XtCFraction, XtRInt, sizeof(int),
      Offset(left), XtRImmediate, (XtPointer) 0
    },
    {
      XtNright, XtCFraction, XtRInt, sizeof(int),
      Offset(right), XtRImmediate, (XtPointer) 100
    }
};
#undef Offset


#define offset(name) XtOffsetOf(FrameRec, frame.name)

static XtResource resources[] = {
    {
      XtNhSpace, XtCHSpace, XtRDimension, sizeof(Dimension),
      offset(h_space), XtRImmediate, (XtPointer)4
    },
    {
      XtNvSpace, XtCVSpace, XtRDimension, sizeof(Dimension),
      offset(v_space), XtRImmediate, (XtPointer)4
    },
    {
      XtNshadowWidth, XtCShadowWidth, XtRDimension, sizeof(Dimension),
      XtOffsetOf(FrameRec, container.shadow_thickness),
      XtRImmediate, (XtPointer) 2
    },
    {
      XtNframeType, XtCFrameType, XtRFrameType, sizeof(XawFrameType),
      offset(frame_type), XtRImmediate, (XtPointer) XawCHISELED
    },
    {
      XtNlayoutPolicy, XtCLayoutPolicy, XtRLayoutPolicy, 
      sizeof(XawLayoutPolicy),  
      offset(policy), XtRImmediate, (XtPointer) XawSINGLE
    }, 
    {
      XtNxFraction, XtCFraction, XtRInt, sizeof(int),
      offset(x_fraction), XtRImmediate, (XtPointer) 100
    },
    {
      XtNyFraction, XtCFraction, XtRInt, sizeof(int),
      offset(y_fraction), XtRImmediate, (XtPointer) 100
    },
    {
      XtNborderWidth, XtCBorderWidth, XtRDimension, sizeof(Dimension),
      XtOffsetOf(RectObjRec,rectangle.border_width), XtRImmediate,
      (XtPointer)0
    },

	/* C A P T I O N */
    {
      XtNcaptionOn, XtCCaptionOn, XtRBoolean, sizeof(Boolean),
      offset(caption), XtRImmediate, (XtPointer)False
    },
    {
      XtNcaptionLabel, XtCCaptionLabel, XtRString, sizeof(char*),
      offset(label), XtRImmediate, (XtPointer)NULL
    },
    {
      XtNfont,  XtCFont, XtRFontStruct, sizeof(XFontStruct *),
      offset(font),XtRString, (XtPointer) XtDefaultFont
    },
    {
      XtNencoding, XtCEncoding, XtRUnsignedChar, sizeof(unsigned char),
      offset(encoding), XtRImmediate, (XtPointer)XawTextEncoding8bit
    },
    {
      XtNjustify, XtCJustify, XtRJustify, sizeof(XtJustify),
      offset(justify), XtRImmediate, (XtPointer)XtJustifyLeft
    },
};
#undef offset

/***************************************************************************
 *
 * Frame  class record
 *
 ***************************************************************************/

static void ClassInitialize();
static void Resize();
static void Redisplay();
static void Initialize();
static void Destroy();
static Boolean SetValues();
static Boolean ConstraintSetValues();
static void ChangeManaged();
static XtGeometryResult GeometryManager();
static XtGeometryResult PreferredGeometry();
static void CalculateNewSize();

#define SuperClass (&containerClassRec)

FrameClassRec frameClassRec = {
  {
    /* superclass	  */	(WidgetClass)SuperClass,
    /* class_name	  */	"Frame",
    /* size		  */	sizeof(FrameRec),
    /* class_initialize	  */	ClassInitialize,
    /* class_part_initialize*/	NULL,
    /* Class init'ed ?	  */	FALSE,
    /* initialize	  */	Initialize,
    /* initialize_hook	  */	NULL,		
    /* realize		  */	XtInheritRealize,
    /* actions		  */	NULL,
    /* num_actions	  */	0,
    /* resources	  */	resources,
    /* resource_count	  */	XtNumber(resources),
    /* xrm_class	  */	NULLQUARK,
    /* compress_motion	  */	FALSE,
    /* compress_exposure  */	TRUE,
    /* compress_enterleave*/	FALSE,
    /* visible_interest	  */	FALSE,
    /* destroy		  */	Destroy,
    /* resize		  */	Resize,
    /* expose		  */	Redisplay,
    /* set_values	  */	SetValues,
    /* set_values_hook	  */	NULL,			
    /* set_values_almost  */	XtInheritSetValuesAlmost,  
    /* get_values_hook	  */	NULL,
    /* accept_focus	  */	NULL,
    /* intrinsics version */	XtVersion,
    /* callback offsets	  */	NULL,
    /* tm_table		  */	NULL,
    /* query_geometry	  */	PreferredGeometry,
    /* display_accelerator*/	XtInheritDisplayAccelerator,
    /* extension	  */	NULL
  },
  {
    /* geometry_manager	  */	GeometryManager,
    /* change_managed	  */	ChangeManaged,
    /* insert_child	  */	XtInheritInsertChild,
    /* delete_child	  */	XtInheritDeleteChild,
    /* extension	  */	NULL
  },
  { /* constraint_class fields */
    /* subresourses       */   frameConstraintResources,
    /* subresource_count  */   XtNumber(frameConstraintResources),
    /* constraint_size    */   sizeof(FrameConstraintsRec),
    /* initialize         */   NULL,
    /* destroy            */   NULL,
    /* set_values         */   ConstraintSetValues,
    /* extension          */   NULL
  },
    /* container class part */
  {
    /* unused             */   0,
  },
  { /* frame_class fields */
    /* extension          */   NULL
  }
};

WidgetClass frameWidgetClass =	(WidgetClass) (&frameClassRec);

/****************************************************************
 *
 * Private Routines
 *
 ****************************************************************/
#define done(type, value)  \
  { \
      if (to->addr != NULL) { \
	  if (to->size < sizeof(type)) { \
	      to->size = sizeof(type); \
	      return False; \
	  } \
	  *(type*)(to->addr) = (value); \
      } else { \
	  static type static_val; \
	  static_val = (value); \
	  to->addr = (XtPointer)&static_val; \
      } \
      to->size = sizeof(type); \
      return True; \
  }

/* ARGSUSED */
static  Boolean
cvtStringToFrameType ( display, args, num_args,  from, to, converter_data)
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
		  "cvtStringToFrameType", "wrongParameters",
		  "XtToolkitError",
		  "String to FrameType conversion needs no arguments",
		  (String*) NULL, (Cardinal*) NULL);
  
  if (XmuCompareISOLatin1(s, "raised")   == 0) done(XawFrameType, XawRAISED);
  if (XmuCompareISOLatin1(s, "sunken")   == 0) done(XawFrameType, XawSUNKEN);
  if (XmuCompareISOLatin1(s, "chiseled") == 0) done(XawFrameType, XawCHISELED);
  if (XmuCompareISOLatin1(s, "ledged")   == 0) done(XawFrameType, XawLEDGED);
  if (XmuCompareISOLatin1(s, "tack")     == 0) done(XawFrameType, XawTACK);

  XtDisplayStringConversionWarning(display, s, XtRFrameType);
  printf("Frame.c");

  done(XawFrameType, XawRAISED);
}

/* ARGSUSED */
static  Boolean
cvtStringToLayoutPolicy ( display, args, num_args,  from, to, converter_data)
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
		  "cvtStringToLayoutPolicy", "wrongParameters",
		  "XtToolkitError",
		  "String to frame type conversion needs no arguments",
		  (String*) NULL, (Cardinal*) NULL);
  
  if (XmuCompareISOLatin1(s,"single") == 0) done(XawLayoutPolicy, XawSINGLE);
  if (XmuCompareISOLatin1(s,"fraction")== 0)done(XawLayoutPolicy, XawFRACTION);
  if (XmuCompareISOLatin1(s,"center") == 0) done(XawLayoutPolicy, XawCENTER);
  XtDisplayStringConversionWarning(display, s, XtRLayoutPolicy);
  done(XawLayoutPolicy, XawSINGLE);
}

static void ClassInitialize()
{
  XawInitializeWidgetSet();

  XtSetTypeConverter(XtRString, XtRFrameType, cvtStringToFrameType,
		     (XtConvertArgList)NULL, 0, XtCacheNone, NULL);
  
  XtSetTypeConverter(XtRString, XtRLayoutPolicy, cvtStringToLayoutPolicy,
		     (XtConvertArgList)NULL, 0, XtCacheNone, NULL);
}

static void SetClip(w)
     Widget w;
{
  XRectangle rectangle[1];
  int label_width;
  int label_x;
  int edge = X_MARGIN(w) + CAPTION_JUSTIFY_SIDE + LINE_WIDTH;

  if (FRAME(w).gc == (GC)NULL)
    return;

  if (FRAME(w).encoding)
    label_width = XTextWidth16 (FRAME(w).font, (TXT16*)FRAME(w).label, 
				(int)strlen(FRAME(w).label)/2);
  else
    label_width = XTextWidth (FRAME(w).font, FRAME(w).label, 
			      (int)strlen(FRAME(w).label));

  switch (FRAME(w).justify) 
  {
  case  XtJustifyLeft :
                                label_x = edge;
    break;
  case XtJustifyCenter :
                                label_x = (CORE(w).width - label_width) / 2;
				label_x = MAX (label_x, edge);
    break;
  case XtJustifyRight :
	                        label_x = 
				  CORE(w).width - label_width - LINE_WIDTH -
				    CAPTION_JUSTIFY_SIDE - X_MARGIN(w);
				label_x = MAX (label_x, edge);
    break;
  }
      
  rectangle[0].x = (short)label_x;
  rectangle[0].y = (short)Y_MARGIN(w);
  rectangle[0].width = (unsigned short)
    MIN (label_width, CORE(w).width - CAPTION_MARGIN(w) - 
	 2 * (LINE_WIDTH + X_MARGIN(w)));
  rectangle[0].height = (unsigned short)
    (FRAME(w).font->max_bounds.ascent + FRAME(w).font->max_bounds.descent);
  
  XSetClipRectangles(XtDisplay(w), FRAME(w).gc, 0, 0, rectangle, 1, YSorted);
}


static void CreateGC(w)
     Widget w;
{
  XtGCMask mask;
  XGCValues values;
  Display *dpy = XtDisplay(w);
  Window   win = XtIsRealized(w) ? XtWindow(w) : XDefaultRootWindow(dpy);

  if (FRAME(w).gc != (GC)NULL)
    XFreeGC(XtDisplay(w), FRAME(w).gc);

  mask = GCForeground | GCFont;

  values.foreground = CONTAINER(w).foreground;
  values.font       = FRAME(w).font->fid;

  FRAME(w).gc = XCreateGC(dpy, win, mask, &values);
  
  SetClip(w); 
}

/* ARGSUSED */
static void Initialize(req, new, args, num_args)
     Widget req, new;
     ArgList args;
     Cardinal *num_args;
{
  FRAME(new).preferred_width  = 1;
  FRAME(new).preferred_height = 1;
  
  if (CORE(new).width == 0)
    CORE(new).width = FRAME(new).preferred_width;
  
  if (CORE(new).height == 0)
    CORE(new).height = FRAME(new).preferred_height;

  FRAME(new).gc = (GC)NULL;
  
  if (FRAME(new).caption && FRAME(new).label != NULL)
    CreateGC (new);
}

static void Destroy(w)
    Widget w;
{
  if (FRAME(w).gc != (GC)NULL)
    XFreeGC(XtDisplay(w), FRAME(w).gc);
}

/* ARGSUSED */
static Boolean SetValues(current, request, new, args, num_args)
    Widget current, request, new;
    ArgList args;
    Cardinal *num_args;
{
    Boolean redisplay = False;

#define NE(name) (FRAME(new).name != FRAME(current).name)

    FRAME(new).policy = FRAME(current).policy;
    
    if (CONTAINER(new).shadow_thickness != CONTAINER(current).shadow_thickness
	|| NE(h_space) || NE(v_space) || NE(caption) 
	|| (FRAME(new).caption && NE(label) && 
	   (FRAME(new).label == NULL || FRAME(current).label == NULL)) )
    {
      CalculateNewSize(new, &CORE(new).width, &CORE(new).height);
      redisplay =  True;
    }

    if (CONTAINER(new).foreground != CONTAINER(current).foreground || NE(font))
    {
      XGCValues values;
      XtGCMask mask;

      values.font       = FRAME(new).font->fid;
      values.foreground = CONTAINER(new).foreground;
      mask = GCForeground | GCFont;

      XChangeGC (XtDisplay(new), FRAME(new).gc, mask, &values);
      redisplay =  True;
    }

    if ((NE(justify) && FRAME(new).caption && FRAME(new).label != NULL) ||
	(NE(label) && FRAME(new).label != NULL ))
    {
      SetClip (new);
      redisplay =  True;
    }

    return redisplay;

#undef NE
}

static void InnerSize (sw, width, height)
     register Widget sw;
     int *width;
     int *height;
{

  *width = CORE(sw).width - 2 * X_MARGIN(sw);
  
  *height = CORE(sw).height - 2 * Y_MARGIN(sw);

  if (FRAME(sw).caption) 
  {
    if (FRAME(sw).label != NULL) {
      *height -= FRAME(sw).font->max_bounds.ascent;
    } 
    else {
      *height -= LINE_WIDTH;
    }
    *height -=      LINE_WIDTH + 2 * FRAME(sw).v_space;
    *width  -= 2 * (LINE_WIDTH + FRAME(sw).h_space);
  }

}

static void Resize(sw)
     register Widget sw;
{
  Widget *child = COMPOSITE(sw).children;
  int i;
  int width;
  int height;
  int x_base;
  int y_base;
  

  InnerSize (sw, &width, &height);

  if (FRAME(sw).caption) 
  {
    x_base =  SHADOW(sw) + LINE_WIDTH + 2 * FRAME(sw).h_space;

    if (FRAME(sw).label != NULL)
      y_base = SHADOW(sw) + 2 * FRAME(sw).v_space + 
	FRAME(sw).font->max_bounds.ascent;
    else
      y_base = SHADOW(sw) + LINE_WIDTH + 2 * FRAME(sw).v_space;
  } else {
    x_base = X_MARGIN(sw);
    y_base = Y_MARGIN(sw);
  }

  if (FRAME(sw).policy == XawCENTER) {
    for(i = 0; i < COMPOSITE(sw).num_children; i++, child++) {
      if (XtIsManaged(*child) && XtIsSubclass (*child, coreWidgetClass)) 
      {
	Position x;
	Position y;
	XtWidgetGeometry preferred;
	
	_XawQueryGeometry (*child, &preferred);

	preferred.width  -= 2 * (*child)->core.border_width;
	preferred.height -= 2 * (*child)->core.border_width;

	width  = MIN (width, preferred.width);
	height = MIN (height,preferred.height);

	x = (CORE(sw).width - width) / 2;

	if (FRAME(sw).caption && FRAME(sw).label != NULL)
	  y = (CORE(sw).height + FRAME(sw).font->max_bounds.ascent - height)/2;
	else
	  y = (CORE(sw).height - height) / 2;

	XtConfigureWidget (*child, x, y, (Dimension)width, (Dimension)height, 
			   (*child)->core.border_width);
			   
	break;
      }
    }
  }
  else if (FRAME(sw).policy == XawFRACTION)
  {
    float xratio = (float)width / (float)(FRAME(sw).x_fraction);
    float yratio = (float)height / (float)(FRAME(sw).y_fraction);

    for(i = 0; i < COMPOSITE(sw).num_children; i++, child++) {
      if (XtIsManaged(*child) && XtIsSubclass(*child,coreWidgetClass)) {
        FrameConstraints frame = (FrameConstraints)CORE(*child).constraints;
	Dimension b = (*child)->core.border_width;
	int x = x_base + xratio * frame->frame.left + b;
	int y = y_base + yratio * frame->frame.top  + b;
	int w = xratio * (frame->frame.right - frame->frame.left) - 2*b;
        int h = yratio * (frame->frame.bottom - frame->frame.top) - 2*b;

	if (w > 0 && h > 0) 
	  XtConfigureWidget( *child, x, y, (Dimension)w, (Dimension)h, b);
      }
    }
  }
  else
  {
    register Dimension child_width;
    register Dimension child_height;

    child_width  = CORE(sw).width - 2 * x_base;
    child_height = CORE(sw).height - y_base - Y_MARGIN(sw) - FRAME(sw).v_space;

    for(i = 0; i < COMPOSITE(sw).num_children; i++, child++) {
      if (XtIsManaged(*child) && XtIsSubclass(*child,coreWidgetClass)) {
	XtConfigureWidget( *child, x_base, y_base, child_width, child_height,
			  (*child)->core.border_width);
	break;
      }
      
    }
  }
  
  SetClip(sw);

}

/* ARGSUSED */
static void Redisplay(gw, event, region)
    Widget gw;
    XEvent *event;		/* unused */
    Region region;		/* unused */
{
  XRectangle rectangle[3];

  
  if (FRAME(gw).caption) 
  {
    if (FRAME(gw).label != NULL) 
    {
      int label_len = strlen(FRAME(gw).label);
      int label_x;
      int label_y;
      int label_width;
      int edge = X_MARGIN(gw) + CAPTION_JUSTIFY_SIDE + LINE_WIDTH;

      if (FRAME(gw).encoding)
	label_width = XTextWidth16 (FRAME(gw).font, (TXT16*)FRAME(gw).label, 
				   (int)strlen(FRAME(gw).label)/2);
      else
	label_width = XTextWidth (FRAME(gw).font, FRAME(gw).label, 
				  (int)strlen(FRAME(gw).label));

      switch (FRAME(gw).justify) 
      {
      case  XtJustifyLeft :
                                 label_x = edge;
	break;
      case XtJustifyCenter :
	                         label_x = (CORE(gw).width - label_width) / 2;
				 label_x = MAX (label_x, edge);
	break;
      case XtJustifyRight :
	                         label_x = 
				   CORE(gw).width - label_width - LINE_WIDTH -
				     CAPTION_JUSTIFY_SIDE - X_MARGIN(gw);
				 label_x = MAX (label_x, edge);
      break;
      }
      
      label_y = Y_MARGIN(gw) + FRAME(gw).font->max_bounds.ascent;
      
      if (FRAME(gw).encoding)
	XDrawString16(XtDisplay(gw), XtWindow(gw), FRAME(gw).gc,
		      label_x, label_y, (TXT16*)FRAME(gw).label,
		      label_len/2);
      else
	XDrawString(XtDisplay(gw), XtWindow(gw), FRAME(gw).gc,
		    label_x, label_y, FRAME(gw).label, 
		    label_len);

      /* Left rectangle */
      rectangle[0].x = (short)0;
      rectangle[0].y = (short)0;
      rectangle[0].width = (unsigned short) (label_x - LINE_WIDTH);
      rectangle[0].height = (unsigned short) CORE(gw).height;

      /* Right rectangle */
      label_width = MIN (label_width, CORE(gw).width - CAPTION_MARGIN(gw) 
			 - 2 * (LINE_WIDTH + X_MARGIN(gw)));

      rectangle[1].x = (short) (label_width + label_x + LINE_WIDTH); 
      rectangle[1].y = (short) 0;
      rectangle[1].width = (unsigned short) (CORE(gw).width - rectangle[1].x);
      rectangle[1].height = (unsigned short) CORE(gw).height;


      /* Bottom rectangle */
      rectangle[2].x = (short)0;
      rectangle[2].y = (short) (CORE(gw).height - Y_MARGIN(gw) - LINE_WIDTH);
      rectangle[2].width = (unsigned short) CORE(gw).width;
      rectangle[2].height = (unsigned short) LINE_WIDTH;

      XSetClipRectangles(XtDisplay(gw), CONTAINER(gw).top_shadow_GC, 
			 0, 0, rectangle, 3, YSorted);

      XSetClipRectangles(XtDisplay(gw), CONTAINER(gw).bottom_shadow_GC, 
			 0, 0, rectangle, 3, YSorted);
    }
/*     else  */
    {
      Position x;
      Position y;
      Dimension w;
      Dimension h;
      
      x = X_MARGIN(gw);
      y = Y_MARGIN(gw);

      if (FRAME(gw).label != NULL)
	y += FRAME(gw).font->max_bounds.ascent -  LINE_WIDTH;

      w = CORE(gw).width - 2 * x;
      h = CORE(gw).height - y - Y_MARGIN(gw);
      
      XawDrawFrame(gw, x, y, w, h, XawCHISELED, LINE_WIDTH,
		   CONTAINER(gw).top_shadow_GC,
		   CONTAINER(gw).bottom_shadow_GC);

      if (FRAME(gw).caption && FRAME(gw).label != NULL) 
      {
	XSetClipMask(XtDisplay(gw), CONTAINER(gw).top_shadow_GC, 
		     (Pixmap)None);
	XSetClipMask(XtDisplay(gw), CONTAINER(gw).bottom_shadow_GC, 
		     (Pixmap)None);
      }
    }
  }
    
  XawDrawFrame(gw,
	       (Position)0,
	       (Position)0,
	       (Dimension)CORE(gw).width,
	       (Dimension)CORE(gw).height,
	       FRAME(gw).frame_type,
	       CONTAINER(gw).shadow_thickness,
	       CONTAINER(gw).top_shadow_GC,
	       CONTAINER(gw).bottom_shadow_GC);

  
}



static void CalculateNewSize(w, width, height)
     Widget w;
     register Dimension *width;
     register Dimension *height;
{

  register Widget sw = w;    
  Widget *child = COMPOSITE(sw).children;
  XtWidgetGeometry reply_return;
  XtGeometryResult result;
  int i;


  *width  = 0;
  *height = 0;

  if (FRAME(sw).policy == XawFRACTION)
  {
    float xratio= -1.;
    float yratio= -1.;

    for(i = 0; i < COMPOSITE(sw).num_children; i++, child++)
    {
      if (XtIsManaged(*child) && XtIsSubclass(*child,coreWidgetClass))
      {
	FrameConstraints frame = (FrameConstraints)CORE(*child).constraints;
	float tmp;
	
	result = XtQueryGeometry (*child, NULL, &reply_return);

	if (result == XtGeometryYes || result == XtGeometryNo) {
	  reply_return.width  = CORE(*child).width;
	  reply_return.height = CORE(*child).height;
	}else {
	  if ((reply_return.request_mode & CWWidth))
	    reply_return.width  = CORE(*child).width;
	  if ((reply_return.request_mode & CWHeight))
	    reply_return.height = CORE(*child).height;
	}
	  
	
	tmp = (float)(reply_return.width + (CORE(*child).border_width << 1)) /
	  (float)(frame->frame.right - frame->frame.left);

	if( tmp > xratio)
	  xratio = tmp;

	tmp = (float)(reply_return.height + (CORE(*child).border_width << 1)) /
	  (float)(frame->frame.bottom - frame->frame.top);

	if( tmp > yratio)
	  yratio = tmp;
      }
    }      
    *width  = FRAME(sw).x_fraction * xratio;
    *height = FRAME(sw).y_fraction * yratio;
  }
  else
  {
    for(i = 0; i < COMPOSITE(sw).num_children; i++, child++)
    {
      if (XtIsManaged(*child) && XtIsSubclass(*child, coreWidgetClass)) 
      {
	result = XtQueryGeometry (*child, NULL, &reply_return);

	if (result == XtGeometryYes || result == XtGeometryNo) {
	  reply_return.width  = CORE(*child).width;
	  reply_return.height = CORE(*child).height;
	}else {
	  if ((reply_return.request_mode & CWWidth))
	    reply_return.width  = CORE(*child).width;
	  if ((reply_return.request_mode & CWHeight))
	    reply_return.height = CORE(*child).height;
	}
	  
	*width  = reply_return.width;
	*height = reply_return.height;
	
	break;
      }
    }
  }

  if (FRAME(sw).caption)
  {
    if (FRAME(sw).label != NULL)
    {
      int label_width;

      *height += FRAME(sw).font->max_bounds.ascent;

      if (FRAME(sw).encoding)
	label_width = XTextWidth16(FRAME(sw).font, (TXT16*)FRAME(sw).label, 
				   (int)strlen(FRAME(sw).label)/2);
      else
	label_width = XTextWidth(FRAME(sw).font, FRAME(sw).label, 
				 strlen(FRAME(sw).label));

      *width = MAX (*width, label_width + CAPTION_MARGIN(sw));
    } else {
      *height += LINE_WIDTH;
    }

    *height +=      LINE_WIDTH + 2 * FRAME(sw).v_space;
    *width  += 2 * (LINE_WIDTH + 2 * FRAME(sw).h_space);
  }

  *width  += 2 * X_MARGIN(sw);
  *height += 2 * Y_MARGIN(sw);

  if (*width == 0)  *width = 1;
  if (*height == 0) *height = 1;

}

/* ARGSUSED */
static void ChangeManaged(w)
    Widget w;
{
  FrameWidget fw = (FrameWidget) w;
  XtWidgetGeometry request;
  XtGeometryResult result;
  
  CalculateNewSize(w, &request.width, &request.height);
  
  if (request.width != fw->core.width || request.height != fw->core.height) {
    request.request_mode = (CWWidth | CWHeight);
    do {
      result = XtMakeGeometryRequest(w, &request, &request);
    } while (result == XtGeometryAlmost);
  }
  Resize(w);
}

/* ARGSUSED */
static XtGeometryResult GeometryManager(w, desired, allowed)
     Widget w;
     XtWidgetGeometry *desired;
     XtWidgetGeometry *allowed; /* unused */
{
  FrameWidget fw = (FrameWidget) XtParent(w);
  XtWidgetGeometry request;
  XtGeometryResult result;
  XtWidgetGeometry save;
  
#define Wants(flag) (desired->request_mode & flag)
  
  if (Wants(XtCWQueryOnly)) return XtGeometryYes;

#define SWAP(name) { save.name = w->core.name; w->core.name = desired->name;}
  if (Wants(CWX))           SWAP(x);
  if (Wants(CWY))           SWAP(y);
  if (Wants(CWWidth))       SWAP(width);
  if (Wants(CWHeight))      SWAP(height);
  if (Wants(CWBorderWidth)) SWAP(border_width);
#undef SWAP
    
  CalculateNewSize((Widget)fw, &request.width, &request.height);

#define SWAP(name) w->core.name = save.name  
  if (Wants(CWX))           SWAP(x);
  if (Wants(CWY))           SWAP(y);
  if (Wants(CWWidth))       SWAP(width);
  if (Wants(CWHeight))      SWAP(height);
  if (Wants(CWBorderWidth)) SWAP(border_width);
#undef SWAP

  if (request.width != fw->core.width || request.height != fw->core.height) {
    request.request_mode = (CWWidth | CWHeight);
    do {
      result = XtMakeGeometryRequest((Widget) fw, &request, &request);
    } while (result == XtGeometryAlmost);
  } 
  Resize((Widget)fw);
  return XtGeometryYes;
  
#undef Wants
}

static XtGeometryResult PreferredGeometry( widget, request, reply  )
    Widget widget;
    XtWidgetGeometry *request, *reply;
{
    FrameWidget w = (FrameWidget)widget;
    
    CalculateNewSize(widget,
		     &w->frame.preferred_width,
		     &w->frame.preferred_height);
    
    reply->width        = w->frame.preferred_width;
    reply->height       = w->frame.preferred_height;
    reply->request_mode = (CWWidth | CWHeight);

    if ((request->request_mode & (CWWidth | CWHeight)) == (CWWidth | CWHeight)
	&&
	request->width  == reply->width
	&&
	request->height == reply->height)
    {
      return XtGeometryYes;
    }
    else if (reply->width == w->core.width && reply->height == w->core.height)
    {	
      return XtGeometryNo;
    }
    else
    {
      return XtGeometryAlmost;
    }
}

/* ARGSUSED */
static Boolean ConstraintSetValues(current, request, new, args, num_args)
    Widget current, request, new;
    ArgList args;
    Cardinal *num_args;
{
  register FrameConstraints cfc = (FrameConstraints) current->core.constraints;
  register FrameConstraints nfc = (FrameConstraints) new->core.constraints;

#define NE(name) (cfc->frame.name != nfc->frame.name)
  return(NE(top) || NE(bottom) || NE(left) || NE(right));
#undef NE
}

