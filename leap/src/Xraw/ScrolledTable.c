#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>

#include "Table2.h"
#include "ScrolledTableP.h"
#include "Scrollbar.h"

#include "XrawDebug.h"


#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define IS_MANAGED(w) (w && XtIsManaged(w))
#define DEFAULT_SCROLL_WIDTH 15
#define STW stw->scrolled_table
#define MIN_SIZE (DEFAULT_SCROLL_WIDTH + 2 * STW.distance)
#define PRINTF_BOOL(a) ( (a) ? "TRUE":"FALSE" )
#define PRINTF_NULL(a) ( (a)==NULL ? "NULL":"not NULL" )

/****************************************************************
 *
 * Full class record constant
 *
 ****************************************************************/

/* Private Data */


static void DefScrollbarWidth(), DefDistance();

#define Offset(field) XtOffsetOf(ScrolledTableRec, scrolled_table.field)

static XtResource resources[] = {
    {
      XtNframeWidth, XtCFrameWidth, XtRDimension, sizeof(Dimension),
      Offset(frame_width), XtRImmediate, (XtPointer)2
    },
    {
      XtNscrollbarWidth, XtCScrollbarWidth, XtRDimension, sizeof(Dimension),
      Offset(scrollbar_width), XtRCallProc, (XtPointer)DefScrollbarWidth 
    },
    {
      XtNframeType, XtCFrameType, XtRFrameType, sizeof(XawFrameType),
      Offset(frame_type), XtRImmediate, (XtPointer) XawSUNKEN
    },
    {
      XtNsignWidget, XtCSignWidget, XtRWidget, sizeof(Widget),
	Offset(sign_widget), XtRImmediate, (XtPointer)NULL
    },
    {
      XtNrowWidget, XtCRowWidget, XtRWidget, sizeof(Widget),
	Offset(row_widget), XtRImmediate, (XtPointer)NULL
    },
    {
      XtNcolumnWidget, XtCColumnWidget, XtRWidget, sizeof(Widget),
	Offset(column_widget), XtRImmediate, (XtPointer)NULL
    },
    {
      XtNstuffWidget, XtCStuffWidget, XtRWidget, sizeof(Widget),
	Offset(stuff_widget), XtRImmediate, (XtPointer)NULL
    },
    {
      XtNallowHoriz, XtCBoolean, XtRBoolean, sizeof(Boolean),
      Offset(allow_horizontal_scrollbar), XtRImmediate, (XtPointer) True
    },
    {
      XtNallowVert, XtCBoolean, XtRBoolean, sizeof(Boolean),
      Offset(allow_vertical_scrollbar), XtRImmediate, (XtPointer) True
    },
    {
      XtNdistance, XtCDistance, XtRDimension, sizeof(Dimension),
      Offset(distance), XtRCallProc, (XtPointer)DefDistance
    },
    {
      XtNforceBars, XtCBoolean, XtCBoolean, sizeof(Boolean),
      Offset(force_scrollbar), XtRImmediate, (XtPointer)False
    },
};

/* ARGSUSED */
static void DefScrollbarWidth( w, offset, value)
     Widget w;
     int offset;
     XrmValue *value;
{
  static Dimension scrollbar_width = DEFAULT_SCROLL_WIDTH;

  value->addr = (caddr_t)&scrollbar_width;
}

/* ARGSUSED */
static void DefDistance( w, offset, value)
     Widget w;
     int offset;
     XrmValue *value;
{
  XrawScrolledTableWidget stw = (XrawScrolledTableWidget) w;
  static Dimension distance;

  distance = 2 * STW.frame_width + 4;

  value->addr = (caddr_t)&distance;
}

#undef offset

static void Initialize();
static void Realize();
static void Redisplay();
static void Destroy();
static void ChangeManaged();
static void DeleteChild();
static void Resize();
static void ThumbProc();
static void ReSetScrollbar();
static void CreateScrollbars ();
static void CalculateLayout();

static Boolean SetValues();

static XtGeometryResult QueryGeometry();
static XtGeometryResult GeometryManager();

static void MoveChild();

#define SuperClass ((ContainerWidgetClass)&containerClassRec)

ScrolledTableClassRec scrolledTableClassRec = {
  { /* core_class fields */	
    /* superclass	  	*/	(WidgetClass) SuperClass,
    /* class_name	  	*/	"ScrolledTable",
    /* widget_size	  	*/	sizeof(ScrolledTableRec),
    /* classInitialize   	*/	NULL,
    /* class_partInitialize	*/	NULL,
    /* class_inited       	*/	FALSE,
    /* initialize	  	*/	Initialize,
    /* initialize_hook		*/	NULL,
    /* realize		  	*/	Realize,
    /* actions		  	*/	NULL,
    /* num_actions	  	*/	0,
    /* resources	  	*/	resources,
    /* num_resources	  	*/	XtNumber(resources),
    /* xrm_class	  	*/	NULLQUARK,
    /* compress_motion	  	*/	TRUE,
    /* compress_exposure  	*/	TRUE,
    /* compress_enterleave	*/	TRUE,
    /* visible_interest	  	*/	TRUE,
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
    /* Composite class part */
  {
    /* geometry manager         */      GeometryManager,
    /* change_managed           */      ChangeManaged,
    /* insert_child             */      XtInheritInsertChild,
    /* delete_child             */      DeleteChild,
    /* extension                */      NULL
  },
    /* Constraint class part */
  {
    /* subresources             */      NULL,
    /* subresource_count        */      0,   
    /* constraint_size          */      0,   
    /* initialize               */      NULL,
    /* destroy                  */      NULL,
    /* set_values               */      NULL,
    /* extension                */      NULL
  },
    /* Container class part */
  {
    /* ignore                   */      0
  },
  { /* ScrolledTable class fields  */
    /* ignore 			*/	0
  }
};

WidgetClass scrolledTableWidgetClass = (WidgetClass)&scrolledTableClassRec;


     
static void SetNewSize(stw)
     XrawScrolledTableWidget stw;
{
  XtWidgetGeometry preferred;
  
  (void)QueryGeometry ((Widget)stw, NULL, &preferred);

  stw->core.width = preferred.width;
  
  stw->core.height = preferred.height;
}


/* ARGSUSED */
static void Initialize(request,new,args,n_args)
     Widget   request;
     Widget   new;
     ArgList  args;
     Cardinal *n_args;
{
  register XrawScrolledTableWidget stw = (XrawScrolledTableWidget)new;
  Arg clip_args[5];
  Cardinal   num_args;
  
    
  num_args = 0;
  XtSetArg(clip_args[num_args], XtNborderWidth,      0);    num_args++;
  XtSetArg(clip_args[num_args], XtNshadowWidth,      0);    num_args++;
  XtSetArg(clip_args[num_args], XtNwidth,            1);    num_args++;
  XtSetArg(clip_args[num_args], XtNheight,           1);    num_args++;

#define CLIP_WIDGET_CLASS  compositeWidgetClass

  /*------------------------ Column ---------------------------------*/
  STW.column_widget = NULL;
  STW.column_clip = 
    XtCreateWidget ("column_clip", CLIP_WIDGET_CLASS, new,
		    clip_args, num_args);
  
  /*------------------------ Row ------------------------------------*/
  STW.row_widget = NULL;
  STW.row_clip = 
    XtCreateWidget ("row_clip", CLIP_WIDGET_CLASS, new,
		    clip_args, num_args);
  
  /*------------------------ Stuff -----------------------------------*/
  STW.stuff_widget = NULL;
  STW.stuff_clip = 
    XtCreateWidget ("stuff_clip", CLIP_WIDGET_CLASS, new,
		    clip_args, num_args);

  /*------------------------ Scrollbars ------------------------------*/

  if (STW.allow_vertical_scrollbar)
    CreateScrollbars (stw, True);

  if (STW.allow_horizontal_scrollbar)
    CreateScrollbars (stw, False);

  /*-----------------------------------------------------------------*/

  if (stw->core.width == 0)
    stw->core.width = MIN_SIZE;
  
  if (stw->core.height == 0)
    stw->core.height = MIN_SIZE;

#undef CLIP_WIDGET_CLASS 
}


#define VERT_BAR_NEEDED (1L<<0)
#define HORIZ_BAR_NEEDED (1L<<1)


static void CalculateLayout(stw)
     XrawScrolledTableWidget stw;
{
  Dimension row_height   = 0;
  Dimension column_width = 0;
  Dimension vertical_scroll_offset;
  Dimension horizontal_scroll_offset;
  Dimension scroll_cross;
  int need = 0;

  int stuff_width;
  int stuff_height;

  if (IS_MANAGED(STW.row_widget))  row_height = FULL_HEIGHT(STW.row_widget);
  else
  if (IS_MANAGED(STW.sign_widget)) row_height = FULL_HEIGHT(STW.sign_widget);
    
  if (IS_MANAGED(STW.column_widget))
    column_width = FULL_WIDTH(STW.column_widget);
  else
  if (IS_MANAGED(STW.sign_widget))
    column_width = FULL_WIDTH(STW.sign_widget);


  STW.stuff_x = column_width != 0 ? column_width + 2 * STW.distance :
    STW.distance;
  STW.stuff_y =   row_height != 0 ?   row_height + 2 * STW.distance :
    STW.distance;

  scroll_cross = STW.scrollbar_width + STW.distance;

  stuff_height = stw->core.height - STW.stuff_y - scroll_cross - STW.distance;

  stuff_width = stw->core.width - STW.stuff_x - scroll_cross - STW.distance;

/* definitions describe correlation full size of children and clip */
#define V_LESS   (STW.stuff_widget == WNULL ? True : \
		  (stuff_height >= FULL_HEIGHT(STW.stuff_widget)))

#define V_MORE   (STW.stuff_widget == WNULL ? False :  \
		  ((stuff_height+scroll_cross)<FULL_HEIGHT(STW.stuff_widget)))


#define H_LESS   (STW.stuff_widget == WNULL ? True : \
		  (stuff_width >= FULL_WIDTH(STW.stuff_widget)))

#define H_MORE   (STW.stuff_widget == WNULL ? False :  \
		  ((stuff_width+scroll_cross) < FULL_WIDTH(STW.stuff_widget)))

/* definitions describe the grade of necessity of scrollbar */
#define V_NO_BAR   (!STW.allow_vertical_scrollbar || STW.v_scroll == WNULL)
#define V_NO_FORCE (STW.allow_vertical_scrollbar && (!STW.force_scrollbar))
#define V_FORCE    (STW.allow_vertical_scrollbar && STW.force_scrollbar)

#define V_MUST     !V_NO_BAR && (V_FORCE  || (V_NO_FORCE && V_MORE))
#define V_REJECT   (V_NO_BAR || (V_NO_FORCE && V_LESS))

#define H_NO_BAR   (!STW.allow_horizontal_scrollbar || STW.h_scroll == WNULL)
#define H_NO_FORCE (STW.allow_horizontal_scrollbar && (!STW.force_scrollbar))
#define H_FORCE    (STW.allow_horizontal_scrollbar && STW.force_scrollbar)

#define H_MUST     !H_NO_BAR && (H_FORCE  || (H_NO_FORCE && H_MORE))
#define H_REJECT   (H_NO_BAR || (H_NO_FORCE && H_LESS))

  if (H_MUST)
  {
    need |= HORIZ_BAR_NEEDED;
    if (V_MUST)
    {
      need |= VERT_BAR_NEEDED;
    }
    else if (V_REJECT)
    {
      need &= ~VERT_BAR_NEEDED;
    }
    else
    {
      need |= VERT_BAR_NEEDED;
    }
  }
  else if (H_REJECT)
  {
    need &= ~HORIZ_BAR_NEEDED;
    if (V_MUST)
    {
      need |= VERT_BAR_NEEDED;
    }
    else if (V_REJECT)
    {
      need &= ~VERT_BAR_NEEDED;
    }
    else
    {
      need &= ~VERT_BAR_NEEDED;
    }
  }
  else
  {
    if (V_MUST)
    {
      need = HORIZ_BAR_NEEDED | VERT_BAR_NEEDED;
    }
  }


  vertical_scroll_offset = (need & VERT_BAR_NEEDED ? 
	    STW.scrollbar_width + 2 * STW.distance : STW.distance);
 
  horizontal_scroll_offset = (need & HORIZ_BAR_NEEDED ? 
	    STW.scrollbar_width + 2 * STW.distance : STW.distance);


  stuff_height = stw->core.height - STW.stuff_y - horizontal_scroll_offset;

  stuff_width = stw->core.width - STW.stuff_x - vertical_scroll_offset;

  STW.stuff_width  = MAX (stuff_width,  MIN_SIZE);
  STW.stuff_height = MAX (stuff_height, MIN_SIZE);

  STW.stuff_width  = MIN (STW.stuff_width, FULL_WIDTH(STW.stuff_widget));
  STW.stuff_height = MIN (STW.stuff_height, FULL_HEIGHT(STW.stuff_widget));

  if (XtIsRealized(STW.v_scroll))
  {
    if (need & VERT_BAR_NEEDED)
      XMapWindow(XtDisplay(STW.v_scroll),XtWindow(STW.v_scroll));
    else
      XUnmapWindow(XtDisplay(STW.v_scroll),XtWindow(STW.v_scroll));
  }

  if (XtIsRealized(STW.h_scroll))
  {
    if (need & HORIZ_BAR_NEEDED)
      XMapWindow(XtDisplay(STW.h_scroll),XtWindow(STW.h_scroll));
    else
      XUnmapWindow(XtDisplay(STW.h_scroll),XtWindow(STW.h_scroll));
  }
}


#undef VERT_BAR_NEEDED
#undef HORIZ_BAR_NEEDED



static void ConfigureChildren(stw)
     XrawScrolledTableWidget stw;
{
  CalculateLayout(stw);
  
  if (IS_MANAGED(STW.sign_widget))
  {
    Dimension w;
    Dimension h;

    w = STW.stuff_x - 2 * (STW.distance + STW.sign_widget->core.border_width);
    h = STW.stuff_y - 2 * (STW.distance + STW.sign_widget->core.border_width);
    
    XtConfigureWidget(STW.sign_widget,
		      STW.distance - STW.sign_widget->core.border_width,
		      STW.distance - STW.sign_widget->core.border_width,
		      w,
		      h,
		      STW.sign_widget->core.border_width);
  }

  if (IS_MANAGED(STW.column_widget))
  {
    XtConfigureWidget(STW.column_clip,
		      STW.distance,
		      STW.stuff_y,
		      MAX (STW.stuff_x - 2 * STW.distance, 1),
		      MAX (STW.stuff_height, 1),
		      0);

    MoveChild (STW.column_widget, STW.column_clip, 
	       STW.column_widget->core.x, STW.column_widget->core.y); 
  }

  if (STW.allow_vertical_scrollbar)
    XtConfigureWidget(STW.v_scroll,
		      STW.stuff_x + STW.stuff_width + STW.distance,
		      STW.stuff_y,
		      STW.scrollbar_width,
		      STW.stuff_height,
		      0);

  if (IS_MANAGED(STW.row_widget))
  {
    XtConfigureWidget(STW.row_clip,
		      STW.stuff_x,
		      STW.distance,
		      MAX (STW.stuff_width, 1),
		      MAX (STW.stuff_y - 2 * STW.distance, 1),
		      0);
  
    MoveChild (STW.row_widget, STW.row_clip, 
	       STW.row_widget->core.x, STW.row_widget->core.y); 
  }
  
  if (STW.allow_horizontal_scrollbar)
    XtConfigureWidget(STW.h_scroll,
		      STW.stuff_x,
		      STW.stuff_y + STW.stuff_height + STW.distance,
		      STW.stuff_width,
		      STW.scrollbar_width,
		      0);

  if (IS_MANAGED(STW.stuff_widget))
  {
    XtConfigureWidget(STW.stuff_clip,
		      STW.stuff_x,
		      STW.stuff_y,
		      STW.stuff_width,
		      STW.stuff_height,
		      0);

    MoveChild (STW.stuff_widget, STW.stuff_clip, 
	       STW.stuff_widget->core.x, STW.stuff_widget->core.y); 
  }

  ReSetScrollbar(stw);
}

static void Destroy(w)
     Widget w;
{
  XrawScrolledTableWidget stw = (XrawScrolledTableWidget)w;

/* returning now since Purify gives bad access msgs on each of
   the following 'if's */
return;

  if (STW.column_clip != WNULL && !STW.column_clip->core.being_destroyed)
    XtDestroyWidget (STW.column_clip);


  if (STW.row_clip != WNULL && !STW.row_clip->core.being_destroyed)
    XtDestroyWidget (STW.row_clip);

  if (STW.stuff_clip != WNULL && !STW.stuff_clip->core.being_destroyed)
    XtDestroyWidget (STW.stuff_clip);

  if (STW.allow_vertical_scrollbar && !STW.v_scroll->core.being_destroyed)
    XtDestroyWidget (STW.v_scroll);

  if (STW.allow_horizontal_scrollbar && !STW.h_scroll->core.being_destroyed)
    XtDestroyWidget (STW.h_scroll);
}
     

static void DeleteChild(w)
     Widget w;
{
  XrawScrolledTableWidget stw =(XrawScrolledTableWidget) XtParent(w);
  /* Boolean need_configure = False; */
  

  if (w == STW.sign_widget)
  {
    STW.sign_widget = NULL;
    /* need_configure = True; */
  }
  else if (w == STW.stuff_widget)
  {
    STW.stuff_widget = NULL;
    /* need_configure = True; */
  }
  else if (w == STW.column_widget)
  {
    STW.column_widget = NULL;
    /* need_configure = True; */
  }
  else if (w == STW.row_widget)
  {
    STW.row_widget = NULL;
    /* need_configure = True; */
  }

  (*SuperClass->composite_class.delete_child) (w);

/*  
  if (need_configure)
    ConfigureChildren(stw);
*/  
}



static void Realize(w, valueMask, attributes)
     Widget w;
     Mask *valueMask;
     XSetWindowAttributes *attributes;
{
  XrawScrolledTableWidget stw = (XrawScrolledTableWidget) w;
    

  *valueMask |= CWBitGravity;
  attributes->bit_gravity = NorthWestGravity;

  (*SuperClass->core_class.realize)(w, valueMask, attributes);


  /*------------------------ Column ---------------------------------*/
  if (IS_MANAGED(STW.column_widget))
    {
      XtManageChild   (STW.column_clip);
      
      XtMoveWidget    (STW.column_widget, (Position)0, (Position)0 );
      XtRealizeWidget (STW.column_clip); 
      XtRealizeWidget (STW.column_widget); 
      
      
      XReparentWindow (XtDisplay(w), 
		       XtWindow(STW.column_widget),
		       XtWindow(STW.column_clip),
		       (Position)0, (Position)0 );

      if (STW.column_widget->core.mapped_when_managed)
	XtMapWidget     (STW.column_widget);
    }
  
  /*------------------------ Row ------------------------------------*/
  if (IS_MANAGED(STW.row_widget))
    {
      XtManageChild   (STW.row_clip);
      
      XtMoveWidget    (STW.row_widget, (Position)0, (Position)0 );
      XtRealizeWidget (STW.row_clip); 
      XtRealizeWidget (STW.row_widget); 
      
      
      XReparentWindow (XtDisplay(w), 
		       XtWindow(STW.row_widget),
		       XtWindow(STW.row_clip),
		       (Position)0, (Position)0 );
      
      if (STW.row_widget->core.mapped_when_managed)
	XtMapWidget     (STW.row_widget);
    }
  
  /*------------------------ Stuff ----------------------------------*/
  if (IS_MANAGED(STW.stuff_widget))
    {
      XtManageChild   (STW.stuff_clip);
      
      XtMoveWidget    (STW.stuff_widget, (Position)0, (Position)0 );
      XtRealizeWidget (STW.stuff_clip); 
      XtRealizeWidget (STW.stuff_widget); 
      
      
      XReparentWindow (XtDisplay(w), 
		       XtWindow(STW.stuff_widget),
		       XtWindow(STW.stuff_clip),
		       (Position)0, (Position)0 );
      
      if (STW.stuff_widget->core.mapped_when_managed)
	XtMapWidget     (STW.stuff_widget);
    }
  
}

static void ChangeManaged(w)
     Widget w;
{
  XrawScrolledTableWidget stw = (XrawScrolledTableWidget)w;
  XtWidgetGeometry stw_request;
  XtWidgetGeometry stw_return;
  XtGeometryResult result;

  _XawQueryGeometry((Widget)stw, &stw_request);

  result = XtMakeGeometryRequest((Widget)stw, &stw_request, &stw_return);
  
  if (result == XtGeometryAlmost)
  {
    if (stw_request.width == stw_return.width &&
	stw_request.height < stw_return.height)
    {
      stw_request.width += STW.scrollbar_width + STW.distance;
    } 
    else if (stw_request.height == stw_return.height &&
	     stw_request.width < stw_return.width)
    {
      stw_request.height += STW.scrollbar_width + STW.distance;
    }
    
    do{
      result = XtMakeGeometryRequest((Widget)stw, &stw_request, &stw_return);
      
    } while (result == XtGeometryAlmost);
  }

  ConfigureChildren (stw);
}

static void CreateScrollbars (stw, vertical)
     XrawScrolledTableWidget stw;
     Boolean vertical;
{

  if (vertical) 
  {
    STW.v_scroll =
      XtVaCreateWidget ("vertical",
			scrollbarWidgetClass, (Widget)stw,
			XtNborderWidth, 0,
			XtNorientation, XtorientVertical,
			XtNhighlightThickness, 0,
			NULL);

    XtAddCallback(STW.v_scroll, XtNvalueChangedProc, ThumbProc, (Widget)stw);
    XtAddCallback(STW.v_scroll, XtNdragProc, ThumbProc, (Widget)stw);
  } 
  else 
  {
    STW.h_scroll =
      XtVaCreateWidget ("horizontal",
			scrollbarWidgetClass, (Widget)stw,
			XtNborderWidth, 0, 
			XtNorientation, XtorientHorizontal,
			XtNhighlightThickness, 0,
			NULL);

    XtAddCallback(STW.h_scroll, XtNvalueChangedProc, ThumbProc, (Widget)stw);
    XtAddCallback(STW.h_scroll, XtNdragProc, ThumbProc, (Widget)stw);
  }
}

/* ARGSUSED */
static Boolean SetValues(old, request, new, args, num_args)
     Widget old, request, new;
     ArgList args;
     Cardinal *num_args;
{
  XrawScrolledTableWidget stw     = (XrawScrolledTableWidget) new;
  XrawScrolledTableWidget stw_old = (XrawScrolledTableWidget) old;
  Boolean redraw = True;
  
#define NE(field) (STW.field != stw_old->scrolled_table.field)
  
  if (NE(allow_vertical_scrollbar)) 
  {
    if (STW.allow_vertical_scrollbar) {
      CreateScrollbars (stw, True);
    }else{
      XtDestroyWidget (STW.v_scroll);
    }
  }
  
  if (NE(allow_horizontal_scrollbar)) 
  {
    if (STW.allow_horizontal_scrollbar) {
      CreateScrollbars (stw, False);
    }else{
      XtDestroyWidget (STW.h_scroll);
    }
  }
  

  if (NE(stuff_widget)             ||
      NE(row_widget)               ||
      NE(column_widget)            ||
      NE(sign_widget)              || 
      NE(allow_vertical_scrollbar) || 
      NE(allow_horizontal_scrollbar))
  {
    SetNewSize(stw);
  }

  if (NE(stuff_widget))
    XtVaSetValues(STW.stuff_widget, 
		  XtNvetricalScroll,   STW.v_scroll,
		  XtNhorizontalScroll, STW.h_scroll,
		  NULL);

  if (NE(row_widget))
    XtVaSetValues(STW.row_widget, 
		  XtNhorizontalScroll, STW.h_scroll,
		  NULL);

  if (NE(column_widget))
    XtVaSetValues(STW.column_widget, 
		  XtNvetricalScroll,   STW.v_scroll,
		  NULL);

  if (NE(frame_type)               ||
      NE(frame_width)              ||
      NE(allow_vertical_scrollbar) || 
      NE(allow_horizontal_scrollbar))
  {
    redraw = True;
  }

  return redraw;
}
#undef NE



static void Redisplay(w,  event, region)
     Widget w;
     XEvent * event;
     Region region;
{
  XrawScrolledTableWidget stw = (XrawScrolledTableWidget) w;

  (*SuperClass->core_class.expose) (w, event, region);
  
  /*------------------------ Column ---------------------------------*/
  if (IS_MANAGED(STW.column_widget))
    XawDrawFrame (w,
		  STW.distance - STW.frame_width,
		  STW.stuff_y - STW.frame_width,
		  STW.stuff_x - 2 * STW.distance + 2 * STW.frame_width,
		  STW.stuff_height + 2 * STW.frame_width,
		  STW.frame_type,
		  STW.frame_width,
		  stw->container.top_shadow_GC,
		  stw->container.bottom_shadow_GC);

  /*------------------------ Row ------------------------------------*/
  if (IS_MANAGED(STW.row_widget))
    XawDrawFrame (w,
		  STW.stuff_x - STW.frame_width,
		  STW.distance - STW.frame_width,
		  STW.stuff_width + 2 * STW.frame_width,
		  STW.stuff_y - 2 * STW.distance + 2 * STW.frame_width,
		  STW.frame_type,
		  STW.frame_width,
		  stw->container.top_shadow_GC,
		  stw->container.bottom_shadow_GC);

  /*------------------------ Stuff ----------------------------------*/
  if (IS_MANAGED(STW.stuff_widget))
    XawDrawFrame (w,
		  STW.stuff_x - STW.frame_width,
		  STW.stuff_y - STW.frame_width,
		  STW.stuff_width + 2 * STW.frame_width,
		  STW.stuff_height + 2 * STW.frame_width,
		  STW.frame_type,
		  STW.frame_width,
		  stw->container.top_shadow_GC,
		  stw->container.bottom_shadow_GC);

}

static void Resize(w)
     Widget w;
{
  XrawScrolledTableWidget stw = (XrawScrolledTableWidget) w;

  
  ConfigureChildren (stw);
  
  if (XtIsRealized(w))
  {
    XClearWindow(XtDisplay(w), XtWindow(w));
    Redisplay(w);
    XFlush(XtDisplay(w));
  }
}



static void ThumbProc(widget, client_data, call_data)
    Widget    widget;
    XtPointer client_data;
    XtPointer call_data;
{
  XrawScrolledTableWidget stw = (XrawScrolledTableWidget) client_data;
  XawScrollBarCallbackStruct *scrollBarStruct;
  float percent;

  scrollBarStruct = (XawScrollBarCallbackStruct*)call_data;
  percent = scrollBarStruct->top;

  if (widget == STW.h_scroll)
  {
    Position x;
    Position y;
    
    if (IS_MANAGED(STW.stuff_widget))
      x = (Position)((float)STW.stuff_widget->core.width * percent);
    else
    if (IS_MANAGED(STW.row_widget))
      x = (Position)((float)STW.row_widget->core.width * percent);
    
    if(IS_MANAGED(STW.row_widget))
    {
      XtVaGetValues(STW.row_widget, XtNy, &y, NULL);
      XtMoveWidget (STW.row_widget, -x , y );
    }
    
    if(IS_MANAGED(STW.stuff_widget))
    {
      XtVaGetValues(STW.stuff_widget, XtNy, &y, NULL);
      XtMoveWidget (STW.stuff_widget, -x , y );
    }
  }
  else if (widget == STW.v_scroll)
  {
    Position x;
    Position y;
    
    if (IS_MANAGED(STW.stuff_widget))
      y = (Position)((float)STW.stuff_widget->core.height * percent);
    else
    if (IS_MANAGED(STW.column_widget))
      y = (Position)((float)STW.column_widget->core.height * percent);
    
    if(IS_MANAGED(STW.column_widget))
    {
      XtVaGetValues(STW.column_widget, XtNx, &x, NULL);
      XtMoveWidget (STW.column_widget, x, -y);
    }
      
    if(IS_MANAGED(STW.stuff_widget))
    {
      XtVaGetValues(STW.stuff_widget, XtNx, &x, NULL);
      XtMoveWidget (STW.stuff_widget, x, -y);
    }
  }

}

static void ReSetScrollbar(stw)
     XrawScrolledTableWidget stw;
{
  Position x;
  Position y;
  float top;
  float shown;
  
  if (STW.h_scroll && XtIsRealized(STW.h_scroll))
  {
    top = 0.0;
    shown = 0.1;
    
    if (IS_MANAGED(STW.stuff_widget))
      shown = (float)STW.stuff_width / (float)STW.stuff_widget->core.width;
    else
      if (IS_MANAGED(STW.row_widget))
	shown = (float)STW.stuff_width / (float)STW.row_widget->core.width;
    
    shown = MIN (shown, 1.0);
    
    if (IS_MANAGED(STW.stuff_widget))
      top = (float)STW.stuff_widget->core.x / 
	(float)STW.stuff_widget->core.width;
    else
      if (IS_MANAGED(STW.row_widget))
	top = (float)STW.row_widget->core.x / 
	  (float)STW.row_widget->core.width;
    
    top = -top;
    
    top = MIN (top, 1.0 - shown);
    
    XawScrollbarSetThumb(STW.h_scroll, top, shown);
      
    if (IS_MANAGED(STW.stuff_widget))
      x = (Position)((float)STW.stuff_widget->core.width * top);
    else
      if (IS_MANAGED(STW.row_widget))
	x = (Position)((float)STW.row_widget->core.width * top);
    
    if(IS_MANAGED(STW.row_widget)) {
      XtVaGetValues(STW.row_widget, XtNy, &y, NULL);
      XtMoveWidget (STW.row_widget, -x , y );
    }
      
    if(IS_MANAGED(STW.stuff_widget)) {
      XtVaGetValues(STW.stuff_widget, XtNy, &y, NULL);
      XtMoveWidget (STW.stuff_widget, -x , y );
    }
  }
      
  if (STW.v_scroll && XtIsRealized(STW.v_scroll))
  {
    top   = 0.0;
    shown = 0.1;
    
    if (IS_MANAGED(STW.stuff_widget))
      shown = (float)STW.stuff_height / (float)STW.stuff_widget->core.height;
    else
      if (IS_MANAGED(STW.column_widget))
	shown = (float)STW.stuff_height / 
	  (float)STW.column_widget->core.height;
  
    shown = MIN (shown, 1.0);
    
    if (IS_MANAGED(STW.stuff_widget))
      top = (float)STW.stuff_widget->core.y / 
	(float)STW.stuff_widget->core.height;
    else
      if (IS_MANAGED(STW.column_widget))
	top = (float)STW.column_widget->core.y / 
	  (float)STW.column_widget->core.height;
    
    top = -top;
    
    top = MIN (top, 1.0 - shown);
    
    XawScrollbarSetThumb(STW.v_scroll, top, shown);
    
    if (IS_MANAGED(STW.stuff_widget))
      y = (Position)((float)STW.stuff_widget->core.height * top);
    else
      if (IS_MANAGED(STW.column_widget))
	y = (Position)((float)STW.column_widget->core.height * top);
    
    if(IS_MANAGED(STW.column_widget)) {
      XtVaGetValues(STW.column_widget, XtNx, &x, NULL);
      XtMoveWidget (STW.column_widget, x, -y);
    }
      
    if(IS_MANAGED(STW.stuff_widget)) {
      XtVaGetValues(STW.stuff_widget, XtNx, &x, NULL);
      XtMoveWidget (STW.stuff_widget, x, -y);
    }
  }
  
}


static XtGeometryResult QueryGeometry(w, intended, preferred)
     Widget            w;
     XtWidgetGeometry *intended;
     XtWidgetGeometry *preferred;
{
  register XrawScrolledTableWidget stw = (XrawScrolledTableWidget) w;

  XtWidgetGeometry sign;
  XtWidgetGeometry row;
  XtWidgetGeometry column;
  XtWidgetGeometry stuff;

  Dimension row_height   = 0;
  Dimension column_width = 0;

  Dimension stuff_width  = 0;
  Dimension stuff_height = 0;

  Dimension row_offset;
  Dimension column_offset;

  /*
   *  QueryGeometry requests to children
   */

  if (IS_MANAGED(STW.sign_widget))
    _XawQueryGeometry (STW.sign_widget, &sign);
  
  if (IS_MANAGED(STW.row_widget))
    _XawQueryGeometry (STW.row_widget, &row);
  
  if (IS_MANAGED(STW.column_widget))
    _XawQueryGeometry (STW.column_widget, &column);
  
  if (IS_MANAGED(STW.stuff_widget))
    _XawQueryGeometry (STW.stuff_widget, &stuff);


  /*
   * Prepare sub sizes
   */
    
  if (IS_MANAGED(STW.row_widget))    row_height = row.height;
  else
  if (IS_MANAGED(STW.sign_widget))   row_height = sign.height;
    
  if (IS_MANAGED(STW.column_widget)) column_width = column.width;
  else
  if (IS_MANAGED(STW.sign_widget))   column_width = sign.width;


  if (IS_MANAGED(STW.stuff_widget))  stuff_height = stuff.height;
  else
  if (IS_MANAGED(STW.column_widget)) stuff_height = column.height;
    
  if (IS_MANAGED(STW.stuff_widget))  stuff_width = stuff.width;
  else
  if (IS_MANAGED(STW.row_widget))    stuff_width = row.width;
    

  column_offset = (column_width != 0 ?
		   column_width + 2 * STW.distance : STW.distance);

  row_offset = (row_height != 0 ? row_height + 2 * STW.distance: STW.distance);

  /*
   * Preferred size
   */
  
  preferred->width        = stuff_width + column_offset + STW.distance;
  preferred->height       = stuff_height + row_offset + STW.distance;
  preferred->request_mode = CWWidth | CWHeight;

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
    if (intended && (intended->request_mode & CWWidth) &&
	intended->width != preferred->width)
      preferred->height += STW.scrollbar_width + STW.distance;

    if (intended && (intended->request_mode & CWHeight) &&
	intended->height != preferred->height)
      preferred->width += STW.scrollbar_width + STW.distance;
    
    return XtGeometryAlmost;
  }

}

/* ARGSUSED */
static XtGeometryResult GeometryManager(child, request, reply)
     Widget            child;
     XtWidgetGeometry *request;
     XtWidgetGeometry *reply;
{
  register XrawScrolledTableWidget stw;
  Dimension width  = child->core.width;
  Dimension height = child->core.height;

  stw = (XrawScrolledTableWidget) XtParent(child);

  if (child == STW.sign_widget &&
      (IS_MANAGED(STW.row_widget) || IS_MANAGED(STW.column_widget)))
    return XtGeometryNo;
  
  if (request->request_mode & CWWidth)
    width = request->width;
  
  if (request->request_mode & CWHeight)
    height = request->height;
  
  if ((request->request_mode & XtCWQueryOnly) != XtCWQueryOnly)
    XtResizeWidget (child, width, height, child->core.border_width);

  if (child == STW.stuff_widget)
    ConfigureChildren (stw);

  return XtGeometryYes;
}

static void MoveChild(child, clip, x, y)
     register Widget child;
     register Widget clip;
     Position x, y;
{
  
  /* make sure we never move past right/bottom borders */
  if (-x + (int)clip->core.width > (int)child->core.width)
    x = -(child->core.width - clip->core.width);
  
  if (-y + (int)clip->core.height > (int)child->core.height)
    y = -(child->core.height - clip->core.height);
  
  /* make sure we never move past left/top borders */
  if (x >= 0) x = 0;
  if (y >= 0) y = 0;
  
  XtMoveWidget(child, x, y);
}

void
XawScrolledTableSetLocation (Widget gw, double xoff, double yoff)
{
  XrawScrolledTableWidget stw = (XrawScrolledTableWidget) gw;
  register Widget child;
  Position x, y;
  
  
  if (!XtIsSubclass (gw, scrolledTableWidgetClass)) {
    String subs[1];
    Cardinal num_subs = 1;
    subs[0] = gw->core.name;
    XtAppWarningMsg(XtWidgetToApplicationContext(gw),
		    "SetLocation", "SetLocation","XawToolkitError",
"routine: XawScrolledTableSetLocation\n\
 Widget '%s' is not subclass of ScrolledTableWidgetClass",
		    subs, &num_subs);
    return;
  }
  
  
  
  stw = (XrawScrolledTableWidget) gw;

  if ((child = STW.stuff_widget) != NULL)
  {
    if (xoff > 1.0)			
      x = child->core.width;
    else if (xoff < 0.0)		
      x = child->core.x;
    else
      x = (Position) (((float) child->core.width) * xoff);
    
    if (yoff > 1.0) 
      y = child->core.height;
    else if (yoff < 0.0)
      y = child->core.y;
    else
      y = (Position) (((float) child->core.height) * yoff);
    
    MoveChild (child, STW.stuff_clip, -x, -y);

    if (STW.h_scroll) 
    {
      float top   = (float)(-(child->core.x))/(float)child->core.width;
      float shown = (float)STW.stuff_clip->core.width/(float)child->core.width;
    
      XawScrollbarSetThumb (STW.h_scroll, top, shown);
    }

    if (STW.v_scroll) 
    {
      float top   = (float)(-(child->core.y))/(float)child->core.height;
      float shown = (float)STW.stuff_clip->core.height/
	(float)child->core.height;
    
      XawScrollbarSetThumb (STW.v_scroll, top, shown);
    }
  }

}

