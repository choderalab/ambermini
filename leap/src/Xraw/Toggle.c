/* $XConsortium: Toggle.c,v 1.24 91/07/25 14:07:48 converse Exp $ */

/*
 * Copyright 1989 Massachusetts Institute of Technology
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of M.I.T. not be used in advertising or
 * publicity pertaining to distribution of the software without specific,
 * written prior permission.  M.I.T. makes no representations about the
 * suitability of this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 *
 * M.I.T. DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL M.I.T.
 * BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN 
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

/*
 * Toggle.c - Toggle button widget
 *
 * Author: Chris D. Peterson
 *         MIT X Consortium 
 *         kit@expo.lcs.mit.edu
 *  
 * Date:   January 12, 1989
 *
 */

#include <stdio.h>
#include <stdint.h>

#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>

#include "../Xmu/Converters.h"
#include "../Xmu/Misc.h"

#include "XawInit.h"
#include "3d.h"
#include "ToggleP.h"

#include "XrawDebug.h"

/****************************************************************
 *
 * Full class record constant
 *
 ****************************************************************/

#define TOGGLE_SING_SIZE 10

#ifndef MIN
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

/* 
 *Private Data 
 */

extern XtActionList xaw_command_actions_list;

/* 
 * The order of toggle and notify are important, as the state has
 * to be set when we call the notify proc.
 */

static char defaultTranslations[] =
#ifdef OLD_TRANSLATIONS
    "<EnterWindow>:	    highlight()	   \n\
     <LeaveWindow>:	    unhighlight()  \n\
     <Btn1Down>,<Btn1Up>:   toggle() notify()";
#else
     "<Btn1Down>,<Btn1Up>:   toggle() notify()";
#endif

#define offset(field) XtOffsetOf(ToggleRec, field)

static XtResource resources[] = { 
   {
     XtNjustify, XtCJustify, XtRJustify, sizeof(XtJustify),
     XtOffsetOf(ToggleRec,label.justify), XtRImmediate,
     (XtPointer)XtJustifyLeft
   },
   {
     XtNstate, XtCState, XtRBoolean, sizeof(Boolean), 
     offset(command.set), XtRString, "off"
   },
   {
     XtNradioGroup, XtCWidget, XtRWidget, sizeof(Widget), 
     offset(toggle.widget), XtRWidget, (XtPointer) NULL
   },
   {
     XtNradioData, XtCRadioData, XtRPointer, sizeof(XtPointer), 
     offset(toggle.radio_data), XtRPointer, (XtPointer) NULL
   },
   {
     XtNshadowWidth, XtCShadowWidth, XtRDimension, sizeof(Dimension),
     offset(simple.shadow_thickness), XtRImmediate, (XtPointer) 0
   }
};
#undef offset


static void ClassInit();
static void Initialize();
static void Redisplay();
static Boolean SetValues();
static XtGeometryResult QueryGeometry();

static void Toggle();
static void Notify();
static void ToggleSet();
static void ToggleDestroy();

/* Functions for handling the Radio Group. */

static RadioGroup * GetRadioGroup();
static void CreateRadioGroup();
static void AddToRadioGroup();
static void TurnOffRadioSiblings();
static void RemoveFromRadioGroup();
static void LabelPosition();

static XtActionsRec actionsList[] =
{
  {"toggle",	        Toggle},
  {"notify",	        Notify},
  {"set",	        ToggleSet}
};

#define SuperClass ((CommandWidgetClass)&commandClassRec)

ToggleClassRec toggleClassRec = {
  {
    (WidgetClass) SuperClass,		/* superclass		  */	
    "Toggle",				/* class_name		  */
    sizeof(ToggleRec),			/* size			  */
    ClassInit,				/* class_initialize	  */
    NULL,				/* class_part_initialize  */
    FALSE,				/* class_inited		  */
    Initialize,				/* initialize		  */
    NULL,				/* initialize_hook	  */
    XtInheritRealize,			/* realize		  */
    actionsList,			/* actions		  */
    XtNumber(actionsList),		/* num_actions		  */
    resources,				/* resources		  */
    XtNumber(resources),		/* resource_count	  */
    NULLQUARK,				/* xrm_class		  */
    FALSE,				/* compress_motion	  */
    TRUE,				/* compress_exposure	  */
    TRUE,				/* compress_enterleave    */
    FALSE,				/* visible_interest	  */
    NULL,         			/* destroy		  */
    XtInheritResize,			/* resize		  */
    Redisplay,  			/* expose		  */
    SetValues,				/* set_values		  */
    NULL,				/* set_values_hook	  */
    XtInheritSetValuesAlmost,		/* set_values_almost	  */
    NULL,				/* get_values_hook	  */
    NULL,				/* accept_focus		  */
    XtVersion,				/* version		  */
    NULL,				/* callback_private	  */
    defaultTranslations,		/* tm_table		  */
    QueryGeometry,		        /* query_geometry	  */
    XtInheritDisplayAccelerator,	/* display_accelerator	  */
    NULL				/* extension		  */
  },  
    /* Simple fields initialization */
  {
    XtInheritChangeSensitive,           /* change_sensitive       */
    XtInheritDisplayRectProc,		/* display_rect           */
    NULL                                /* extension              */
  }, 
    /* Label fields initialization */
  {
    LabelPosition,                      /* label_position    */
  }, 
    /* Commmand fields initialization */
  {
    0                                     /* field not used    */
  }, 
    /* Toggle fields initialization */
  {
      NULL,			        /* Set Procedure. */
      NULL,			        /* Unset Procedure. */
      NULL			        /* extension. */
  } 
};

  /* for public consumption */
WidgetClass toggleWidgetClass = (WidgetClass) &toggleClassRec;

/****************************************************************
 *
 * Private Procedures
 *
 ****************************************************************/

static void
ClassInit()
{
  XtActionList actions;
  Cardinal num_actions;
  Cardinal i;
  ToggleWidgetClass class = (ToggleWidgetClass) toggleWidgetClass;
  static XtConvertArgRec parentCvtArgs[] = {
      {XtBaseOffset, (XtPointer)XtOffsetOf(WidgetRec, core.parent),
	   sizeof(Widget)}
  };

  XawInitializeWidgetSet();

  XtAddConverter(XtRString, XtRWidget, XmuCvtStringToWidget,
		     parentCvtArgs, XtNumber(parentCvtArgs));
/* 
 * Find the set and unset actions in the command widget's action table. 
 */

#if defined(XtSpecificationRelease) && XtSpecificationRelease > 4

  XtGetActionList(commandWidgetClass, &actions, &num_actions);

#else

  actions = xaw_command_actions_list;
  num_actions = SuperClass->core_class.num_actions;

#endif
  
  for (i = 0 ; i < num_actions ; i++) {
    if (streq(actions[i].string, "set"))
	class->toggle_class.Set = actions[i].proc;
    if (streq(actions[i].string, "unset")) 
	class->toggle_class.Unset = actions[i].proc;

    if ( (class->toggle_class.Set != NULL) &&
	 (class->toggle_class.Unset != NULL) ) {
#if defined(XtSpecificationRelease) && XtSpecificationRelease > 4
	XtFree((char *) actions);
#endif


	return;
    }
  }  

/* We should never get here. */
  XtError("Aborting, due to errors resolving bindings in the Toggle widget.");
}

static void Initialize(request, new)
 Widget request, new;
{
    ToggleWidget tw = (ToggleWidget) new;
    ToggleWidget tw_req = (ToggleWidget) request;


    if (tw->label.left_bitmap == None && tw->label.need_calculate_width) {
      tw->core.width = tw->label.label_width +
	               3 * tw->label.internal_width +
		       2 * SIMPLE_MARGIN(tw) +
		       TOGGLE_SING_SIZE;
    }
      
    tw->toggle.radio_group = NULL;

    if (tw->toggle.radio_data == NULL) 
      tw->toggle.radio_data = (XtPointer) new->core.name;

    if (tw->toggle.widget != NULL) {
      if ( GetRadioGroup(tw->toggle.widget) == NULL) 
	CreateRadioGroup(new, tw->toggle.widget);
      else
	AddToRadioGroup( GetRadioGroup(tw->toggle.widget), new);
    }      
    XtAddCallback(new, XtNdestroyCallback, ToggleDestroy, NULL);

/*
 * Command widget assumes that the widget is unset, so we only 
 * have to handle the case where it needs to be set.
 *
 * If this widget is in a radio group then it may cause another
 * widget to be unset, thus calling the notify proceedure.
 *
 * I want to set the toggle if the user set the state to "On" in 
 * the resource group, reguardless of what my ancestors did.
 */

    if (tw_req->command.set)
      ToggleSet(new, NULL, NULL, 0);
}


static void Redisplay(gw, event, region)
     Widget gw;
     XEvent *event;
     Region region;
{
    register ToggleWidget w = (ToggleWidget) gw;
    Position x,y;
    XRectangle rectangle;
    short g = TOGGLE_SING_SIZE / 2 + 1;

    if (!XtIsRealized(gw))
      return;

    (*XtLabelClass(w)->simple_class.display_rect)(gw, &rectangle);
    
    if (w->command.changed_highlight) {
      if (w->command.highlighted != HighlightNone) {
	XFillRectangle(XtDisplay(w), XtWindow(w), w->command.armed_GC,
		       rectangle.x , rectangle.y,
		       rectangle.width, rectangle.height);
      } else {
	XClearArea(XtDisplay(gw), XtWindow(gw),
		   rectangle.x , rectangle.y,
		   rectangle.width, rectangle.height, False);
      }
      w->command.changed_highlight = False;
    } 
    else if (w->command.changed_set) 
    {
      XClearArea(XtDisplay(gw), XtWindow(gw),
		 rectangle.x , rectangle.y,
		 rectangle.width, rectangle.height, False);
      w->command.changed_set = False;
    }
  
    (*labelClassRec.core_class.expose)(gw, event, region);

      
    if (w->toggle.radio_group)
    {
      x = SIMPLE_MARGIN(w) + (w)->label.internal_width + g;
      y = w->core.height / 2;
      
      DrawRhombus(gw,
		  (short)x, (short)y, g, 2,
		  w->simple.top_shadow_GC,
		  w->command.armed_GC,
		  w->simple.bottom_shadow_GC,
		  w->command.set);
      
    } else {
      x = SIMPLE_MARGIN(w) + (w)->label.internal_width;
      y = (w->core.height - TOGGLE_SING_SIZE) / 2;
    
      if (w->command.set) {
	XFillRectangle(XtDisplay(w), XtWindow(w), w->command.armed_GC,
		       x, y, TOGGLE_SING_SIZE, TOGGLE_SING_SIZE);
      }
      
      XawDrawFrame (gw,
		    x,
		    y,
		    (Dimension)TOGGLE_SING_SIZE,
		    (Dimension)TOGGLE_SING_SIZE,
		    w->command.set ? XawSUNKEN : XawRAISED,
		    2,
		    w->simple.top_shadow_GC,
		    w->simple.bottom_shadow_GC);
    }
    
}



/************************************************************
 *
 *  Action Procedures
 *
 ************************************************************/

/* ARGSUSED */
static void 
ToggleSet(w,event,params,num_params)
Widget w;
XEvent *event;
String *params;		/* unused */
Cardinal *num_params;	/* unused */
{
    ToggleWidgetClass class = (ToggleWidgetClass) w->core.widget_class;

    TurnOffRadioSiblings(w);
    (*class->toggle_class.Set)(w, event, NULL, 0);
}

/* ARGSUSED */
static void Toggle(w,event,params,num_params)
     Widget w;
     XEvent *event;
     String *params;		/* unused */
     Cardinal *num_params;	/* unused */
{
  ToggleWidget tw = (ToggleWidget)w;
  ToggleWidgetClass class = (ToggleWidgetClass) w->core.widget_class;

  if (tw->command.set && tw->toggle.radio_group == NULL) /* vtr */
    class->toggle_class.Unset(w, event, NULL, 0);
  else 
    ToggleSet(w, event, params, num_params);
}

/* ARGSUSED */
static void Notify(w,event,params,num_params)
     Widget w;
     XEvent *event;
     String *params;		/* unused */
     Cardinal *num_params;	/* unused */
{
  ToggleWidget tw = (ToggleWidget) w;
  XtCallCallbacks(w, XtNcallback, (XtPointer)(intptr_t)tw->command.set);
}

/************************************************************
 *
 * Set specified arguments into widget
 *
 ***********************************************************/

/* ARGSUSED */
static Boolean SetValues (current, request, new, args, num_args)
     Widget current, request, new;
     ArgList args;
     Cardinal *num_args;
{
    ToggleWidget oldtw = (ToggleWidget) current;
    ToggleWidget tw = (ToggleWidget) new;
    ToggleWidget rtw = (ToggleWidget) request;

    if (oldtw->toggle.widget != tw->toggle.widget)
      XawToggleChangeRadioGroup(new, tw->toggle.widget);

    if (!tw->core.sensitive && oldtw->core.sensitive && rtw->command.set)
	tw->command.set = True;

    if (oldtw->command.set != tw->command.set) {
	tw->command.set = oldtw->command.set;
	Toggle(new, NULL, NULL, 0); /* Does a redisplay. */
	Notify(new, NULL, NULL, 0); /* Invokes callback V.T.R. */
    }
    return(FALSE);
}

/*	Function Name: ToggleDestroy
 *	Description: Destroy Callback for toggle widget.
 *	Arguments: w - the toggle widget that is being destroyed.
 *                 junk, grabage - not used.
 *	Returns: none.
 */

/* ARGSUSED */
static void
ToggleDestroy(w, junk, garbage)
Widget w;
XtPointer junk, garbage;
{
  RemoveFromRadioGroup(w);
}

/************************************************************
 *
 * Below are all the private proceedures that handle 
 * radio toggle buttons.
 *
 ************************************************************/

/*	Function Name: GetRadioGroup
 *	Description: Gets the radio group associated with a give toggle
 *                   widget.
 *	Arguments: w - the toggle widget who's radio group we are getting.
 *	Returns: the radio group associated with this toggle group.
 */

static RadioGroup *
GetRadioGroup(w)
Widget w;
{
  ToggleWidget tw = (ToggleWidget) w;

  if (tw == NULL) return(NULL);
  return( tw->toggle.radio_group );
}

/*	Function Name: CreateRadioGroup
 *	Description: Creates a radio group. give two widgets.
 *	Arguments: w1, w2 - the toggle widgets to add to the radio group.
 *	Returns: none.
 * 
 *      NOTE:  A pointer to the group is added to each widget's radio_group
 *             field.
 */

static void
CreateRadioGroup(w1, w2)
Widget w1, w2;
{
  char error_buf[BUFSIZ];
  ToggleWidget tw1 = (ToggleWidget) w1;
  ToggleWidget tw2 = (ToggleWidget) w2;

  if ( (tw1->toggle.radio_group != NULL) || (tw2->toggle.radio_group != NULL) ) {
    sprintf(error_buf, "%s %s", "Toggle Widget Error - Attempting",
	    "to create a new toggle group, when one already exists.");
    XtWarning(error_buf);
  }

  AddToRadioGroup( NULL, w1 );
  AddToRadioGroup( GetRadioGroup(w1), w2 );
}

/*	Function Name: AddToRadioGroup
 *	Description: Adds a toggle to the radio group.
 *	Arguments: group - any element of the radio group the we are adding to.
 *                 w - the new toggle widget to add to the group.
 *	Returns: none.
 */

static void
AddToRadioGroup(group, w)
RadioGroup * group;
Widget w;
{
  ToggleWidget tw = (ToggleWidget) w;
  RadioGroup * local;

  local = (RadioGroup *) XtMalloc( sizeof(RadioGroup) );
  local->widget = w;
  tw->toggle.radio_group = local;

  if (group == NULL) {		/* Creating new group. */
    group = local;
    group->next = NULL;
    group->prev = NULL;
    return;
  }
  local->prev = group;		/* Adding to previous group. */
  if ((local->next = group->next) != NULL)
      local->next->prev = local;
  group->next = local;
}

/*	Function Name: TurnOffRadioSiblings
 *	Description: Deactivates all radio siblings.
 *	Arguments: widget - a toggle widget.
 *	Returns: none.
 */

static void
TurnOffRadioSiblings(w)
Widget w;
{
  RadioGroup * group;
  ToggleWidgetClass class = (ToggleWidgetClass) w->core.widget_class;

  if ( (group = GetRadioGroup(w)) == NULL)  /* Punt if there is no group */
    return;

  /* Go to the top of the group. */

  for ( ; group->prev != NULL ; group = group->prev );

  while ( group != NULL ) {
    ToggleWidget local_tog = (ToggleWidget) group->widget;
    if ( local_tog->command.set ) {
      class->toggle_class.Unset(group->widget, NULL, NULL, 0);
      Notify( group->widget, NULL, NULL, 0);
    }
    group = group->next;
  }
}

/*	Function Name: RemoveFromRadioGroup
 *	Description: Removes a toggle from a RadioGroup.
 *	Arguments: w - the toggle widget to remove.
 *	Returns: none.
 */

static void
RemoveFromRadioGroup(w)
Widget w;
{
  RadioGroup * group = GetRadioGroup(w);
  if (group != NULL) {
    if (group->prev != NULL)
      (group->prev)->next = group->next;
    if (group->next != NULL)
      (group->next)->prev = group->prev;
    XtFree((char *) group);
  }
}

/************************************************************
 *
 * Public Routines
 *
 ************************************************************/
   
/*	Function Name: XawToggleChangeRadioGroup
 *	Description: Allows a toggle widget to change radio groups.
 *	Arguments: w - The toggle widget to change groups.
 *                 radio_group - any widget in the new group.
 *	Returns: none.
 */

void
XawToggleChangeRadioGroup(Widget w, Widget radio_group)
{
  ToggleWidget tw = (ToggleWidget) w;
  RadioGroup * group;

  RemoveFromRadioGroup(w);

/*
 * If the toggle that we are about to add is set then we will 
 * unset all toggles in the new radio group.
 */

  if ( tw->command.set && radio_group != NULL )
    XawToggleUnsetCurrent(radio_group);

  if (radio_group != NULL) {
      if ((group = GetRadioGroup(radio_group)) == NULL)
	  CreateRadioGroup(w, radio_group);
      else AddToRadioGroup(group, w);
  }
}

/*	Function Name: XawToggleGetCurrent
 *	Description: Returns the RadioData associated with the toggle
 *                   widget that is currently active in a toggle group.
 *	Arguments: w - any toggle widget in the toggle group.
 *	Returns: The XtNradioData associated with the toggle widget.
 */

XtPointer
XawToggleGetCurrent(Widget w)
{
  RadioGroup * group;

  if ( (group = GetRadioGroup(w)) == NULL) return(NULL);
  for ( ; group->prev != NULL ; group = group->prev);

  while ( group != NULL ) {
    ToggleWidget local_tog = (ToggleWidget) group->widget;
    if ( local_tog->command.set )
      return( local_tog->toggle.radio_data );
    group = group->next;
  }
  return(NULL);
}

/*	Function Name: XawToggleSetCurrent
 *	Description: Sets the Toggle widget associated with the
 *                   radio_data specified.
 *	Arguments: radio_group - any toggle widget in the toggle group.
 *                 radio_data - radio data of the toggle widget to set.
 *	Returns: none.
 */

void
XawToggleSetCurrent(Widget radio_group, XtPointer radio_data)
{
  RadioGroup * group;
  ToggleWidget local_tog; 

/* Special case case of no radio group. */

  if ( (group = GetRadioGroup(radio_group)) == NULL) {
    local_tog = (ToggleWidget) radio_group;
    if ( (local_tog->toggle.radio_data == radio_data) )     
      if (!local_tog->command.set) {
	ToggleSet((Widget) local_tog, NULL, NULL, 0);
	Notify((Widget) local_tog, NULL, NULL, 0);
      }
    return;
  }

/*
 * find top of radio_roup 
 */

  for ( ; group->prev != NULL ; group = group->prev);

/*
 * search for matching radio data.
 */

  while ( group != NULL ) {
    local_tog = (ToggleWidget) group->widget;
    if ( (local_tog->toggle.radio_data == radio_data) ) {
      if (!local_tog->command.set) { /* if not already set. */
	ToggleSet((Widget) local_tog, NULL, NULL, 0);
	Notify((Widget) local_tog, NULL, NULL, 0);
      }
      return;			/* found it, done */
    }
    group = group->next;
  }
}
 
/*	Function Name: XawToggleUnsetCurrent
 *	Description: Unsets all Toggles in the radio_group specified.
 *	Arguments: radio_group - any toggle widget in the toggle group.
 *	Returns: none.
 */

void
XawToggleUnsetCurrent(Widget radio_group)
{
  ToggleWidgetClass class;
  ToggleWidget local_tog = (ToggleWidget) radio_group;

  /* Special Case no radio group. */

  if (local_tog->command.set) {
    class = (ToggleWidgetClass) local_tog->core.widget_class;
    class->toggle_class.Unset(radio_group, NULL, NULL, 0);
    Notify(radio_group, NULL, NULL, 0);
  }
  if ( GetRadioGroup(radio_group) == NULL) return;
  TurnOffRadioSiblings(radio_group);
}


static void LabelPosition(w)
     Widget w;
{
  register ToggleWidget tw = (ToggleWidget)w;
  Position width, height;
  Position newPos;
  Position leftedge;
  

  if (!XtIsSubclass(w, toggleWidgetClass))
    return ;

  width    = tw->core.width;
  height   = tw->core.height;
  leftedge = tw->label.internal_width * 2 +
             SIMPLE_MARGIN(tw) + TOGGLE_SING_SIZE;


  switch (tw->label.justify) {

  case XtJustifyLeft   :	 newPos = leftedge;

    break;
  case XtJustifyRight  :	 newPos = width - (tw->label.label_width +
						   tw->label.internal_width +
						   SIMPLE_MARGIN(tw));
    break;
  case XtJustifyCenter :
  default              :         newPos = (width - tw->label.label_width) / 2;
    break;
  }

  if (newPos < leftedge)
    tw->label.label_x = leftedge;
  else
    tw->label.label_x = newPos;

  tw->label.label_y = (height - tw->label.label_height) / 2;
  

  return;
  
}

static XtGeometryResult QueryGeometry(w, intended, preferred)
    Widget w;
    XtWidgetGeometry *intended, *preferred;
{
    register ToggleWidget tw = (ToggleWidget)w;


    (void)(*labelWidgetClass->core_class.query_geometry) (w, NULL, preferred);
    
    if (tw->label.left_bitmap == None)
    {
      preferred->width = (tw->label.label_width +
			  3 * tw->label.internal_width +
			  2 * SIMPLE_MARGIN(tw) +
			  TOGGLE_SING_SIZE);

      preferred->height = (MAX(tw->label.label_height,TOGGLE_SING_SIZE) +
			   2 * tw->label.internal_height +
			   2 * SIMPLE_MARGIN(tw));
      
    }
    
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

