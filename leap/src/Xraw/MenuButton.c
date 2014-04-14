/* $XConsortium: MenuButton.c,v 1.18 91/06/22 18:03:46 rws Exp $ */

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


/***********************************************************************
 *
 * MenuButton Widget
 *
 ***********************************************************************/

/*
 * MenuButton.c - Source code for MenuButton widget.
 *
 * This is the source code for the Athena MenuButton widget.
 * It is intended to provide an easy method of activating pulldown menus.
 *
 * Date:    May 2, 1989
 *
 * By:      Chris D. Peterson
 *          MIT X Consortium 
 *          kit@expo.lcs.mit.edu
 *
 * --------------------------------
 *
 * Date:    February, 1996
 *
 * Changes: Vladimir T. Romanovski
 *          romsky@oea.ihep.su       // IHEP (Russia)
 *          romsky@munin.ucsf.edu    // University of California San Francisco
 *
 */


#include <stdio.h>
#include <signal.h>

#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#include <X11/cursorfont.h>


#include "XawInit.h"
#include "PanedP.h"
#include "MenuButtoP.h"
#include "SimpleMenP.h"

#include "XrawDebug.h"


#define CALLOC(num,type) (num * sizeof(type) <= 0 ? (type*)NULL : \
  (type*)XtMalloc((Cardinal)(num) * (Cardinal)sizeof(type)))


static char defaultTranslations[] = 
    "<EnterWindow>:     EnterWindow() highlight()   \n\
     <LeaveWindow>:     reset()                     \n\
     <BtnDown>:         PopupMenu()                 \n\
     <BtnUp>:           reset()                     \n";


/****************************************************************
 *
 * Full class record constant
 *
 ****************************************************************/

/* Private Data */

#define offset(field) XtOffsetOf(MenuButtonRec, field)
static XtResource resources[] = {
  {
    XtNmenuName, XtCMenuName, XtRString, sizeof(String), 
    offset(menu_button.menu_name), XtRString, (XtPointer)"menu"
  },
  {
    XtNhighlightThickness, XtCHighlightThickness, XtRDimension,
    sizeof(Dimension),
    offset(simple.highlight_thickness), XtRImmediate, (XtPointer) 2
  },
/*
  {
    XtNinternalHeight, XtCHeight, XtRDimension, sizeof(Dimension),
    offset(label.internal_height), XtRImmediate, (XtPointer) 6
  },
  {
    XtNinternalWidth, XtCWidth, XtRDimension, sizeof(Dimension),
    offset(label.internal_width), XtRImmediate, (XtPointer) 6
  }
*/
};
#undef offset

static void ClassInitialize();
static void PopupMenu();
static void EnterWindow();

static void Unhighlight();
static void Handler_ButtonRelease();

static XtActionsRec actionsList[] =
{
  {"PopupMenu",	  PopupMenu},
  {"EnterWindow", EnterWindow}
};

#define SuperClass ((CommandWidgetClass)&commandClassRec)

MenuButtonClassRec menuButtonClassRec = {
  {
    (WidgetClass) SuperClass,		/* superclass		  */	
    "MenuButton",			/* class_name		  */
    sizeof(MenuButtonRec),       	/* size			  */
    ClassInitialize,			/* class_initialize	  */
    NULL,				/* class_part_initialize  */
    FALSE,				/* class_inited		  */
    NULL,				/* initialize		  */
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
    NULL,				/* destroy		  */
    XtInheritResize,			/* resize		  */
    XtInheritExpose,			/* expose		  */
    NULL,				/* set_values		  */
    NULL,				/* set_values_hook	  */
    XtInheritSetValuesAlmost,		/* set_values_almost	  */
    NULL,				/* get_values_hook	  */
    NULL,				/* accept_focus		  */
    XtVersion,				/* version		  */
    NULL,				/* callback_private	  */
    defaultTranslations,               	/* tm_table		  */
    XtInheritQueryGeometry,		/* query_geometry	  */
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
    XtInheritLabelPosition              /* label_position    */
  }, 
  { /* Command fields initialization */
    0                                   /* field not used    */
  },
  { /* MenuButton fields initialization */
    0                                   /* field not used    */
  }  
};

  /* for public consumption */
WidgetClass menuButtonWidgetClass = (WidgetClass) &menuButtonClassRec;

/****************************************************************
 *
 * Private Data Types
 *
 ****************************************************************/

typedef struct _GrabDataRec{
  struct _GrabDataRec *pred;
  struct _GrabDataRec *next;
  Widget menu;
  Widget menubar;
  Widget menu_button;
}GrabDataRec, *GrabData;


static struct {
  GrabData       list;
} multiDisplayGrabData;


static GrabData CreateGrabData(menubar, button, menu)
     Widget menubar;
     Widget button;
     Widget menu;
{
  GrabDataRec* data = XtNew (GrabDataRec);

  data->menu         = menu;
  data->menubar      = menubar;
  data->menu_button  = button;

  data->pred = (GrabData)NULL;
  data->next = multiDisplayGrabData.list;

  if (data->next != (GrabData)NULL)
    data->next->pred = (GrabData)data;

  multiDisplayGrabData.list = (GrabData)data;
  return (GrabData)data;
}


static void FreeGrabData(data)
     register GrabData data;
{
  if (data->next != (GrabData)NULL) 
    data->next->pred = data->pred;

  if (data->pred != (GrabData)NULL)
    data->pred->next = data->next;
  else
    multiDisplayGrabData.list = data->next;

  XtFree((char*)data);

}

/****************************************************************
 *
 * Private Procedures
 *
 ****************************************************************/

static void ClassInitialize()
{
  XawInitializeWidgetSet();
  XtRegisterGrabAction(PopupMenu, True, ButtonPressMask | ButtonReleaseMask,
		       GrabModeAsync, GrabModeAsync);

  multiDisplayGrabData.list = NULL;
}

static Widget FindPopupMenu(w)
     register Widget w;
{
  register MenuButtonWidget mbw = (MenuButtonWidget) w;
  Widget                    menu = WNULL;

  while(w != NULL) 
  {
    if ((menu = XtNameToWidget(w, mbw->menu_button.menu_name)) == NULL) 
      w = XtParent(w);
    else
      break;
  }
  return menu;
}

static Widget GetPopupMenu(w)
     Widget w;
{
  MenuButtonWidget mbw = (MenuButtonWidget) w;
  Widget           menu;
  Position         menu_x;
  Position         menu_y;
  Dimension        menu_width;
  Dimension        menu_height;

  menu = FindPopupMenu(w);

  if (menu == NULL) {
    char error_buf[BUFSIZ];
    sprintf(error_buf, "MenuButton: %s %s.",
	    "Could not find menu widget named", mbw->menu_button.menu_name);
    XtAppWarning(XtWidgetToApplicationContext(w), error_buf);
    return WNULL;
  }
  
  if (!XtIsRealized(menu))
    XtRealizeWidget(menu);
  
  menu_width  = menu->core.width + 2 * menu->core.border_width;
  menu_height = menu->core.height + 2 * menu->core.border_width;

  XtTranslateCoords(w, 0, 0, &menu_x, &menu_y);

  menu_x += mbw->simple.highlight_thickness;
  menu_y += w->core.height -
            mbw->simple.highlight_thickness +
	    menu->core.border_width;

  if (menu_x >= 0)  {
    int scr_width = WidthOfScreen(XtScreen(menu));
    if (menu_x + menu_width > scr_width)
      menu_x = scr_width - menu_width;
  }
  if (menu_x < 0) 
    menu_x = 0;

  if (menu_y >= 0)  {
    int scr_height = HeightOfScreen(XtScreen(menu));
    if (menu_y + menu_height > scr_height)
      menu_y = scr_height - menu_height;
  }
  if (menu_y < 0)
    menu_y = 0;

  XtVaSetValues(menu, XtNx, menu_x, XtNy, menu_y, NULL);

  return menu;
}

static void Unhighlight (gw)
     Widget gw;
{
  Dimension thick;
  Dimension shadow;
  Dimension width;
  Dimension height;
  Display   *dpy = XtDisplay(gw);
  Window    win = XtWindow(gw);
  

  if ((shadow = ((MenuButtonWidget)gw)->simple.shadow_thickness) == 0)
    return;
  
  thick  = ((MenuButtonWidget)gw)->simple.highlight_thickness;
  width  = gw->core.width - 2 * thick;
  height = gw->core.height - 2 * thick;
  
    
  XClearArea(dpy, win, thick, thick, shadow, height , False);
  XClearArea(dpy, win, thick, gw->core.height - thick -shadow, 
	     width, shadow, False);

  XClearArea(dpy, win, thick, thick, width, shadow, False);
  XClearArea(dpy, win, gw->core.width - thick - shadow, thick,
	     shadow, height, False);
}

static GrabData GetGrabData(w)
     Widget w;
{
  register GrabData data;
 
  for (data = multiDisplayGrabData.list;
       data != (GrabData)NULL;
       data = data->next)
    {
      if (XtDisplay(w) == XtDisplay(data->menu))
        return data;
    }
  return (GrabData)NULL;
}

/****************************************************************
 *
 * Actions
 *
 ****************************************************************/

/* ARGSUSED */
static void EnterWindow(w, event, params, num_params)
     Widget     w;
     XEvent     *event;      /* unused */
     String     *params;     /* unused */
     Cardinal   *num_params; /* unused */
{
  MenuButtonWidget mbw = (MenuButtonWidget) w;
  register GrabData data;
  register Widget new_menu;

  if ((event->xcrossing.state & Button1Mask) != Button1Mask)
    return ;


  if (((data = GetGrabData (w)) != NULL)
      && (data->menubar == XtParent(w))
      && ((new_menu = GetPopupMenu(w)) != WNULL)
      && data->menu != new_menu)
  {
    XtPopdown (data->menu);
    XtRemoveEventHandler (data->menu,
			  ButtonReleaseMask,
			  False,
			  Handler_ButtonRelease,
			  (XtPointer)data);

    Unhighlight(data->menu_button);

    XawDrawFrame (w,
		  mbw->simple.highlight_thickness,
		  mbw->simple.highlight_thickness,
	     (Dimension)(w->core.width - 2*mbw->simple.highlight_thickness),
	     (Dimension)(w->core.height - 2*mbw->simple.highlight_thickness),
	          XawRAISED,
		  mbw->simple.shadow_thickness,
		  mbw->simple.top_shadow_GC,
		  mbw->simple.bottom_shadow_GC);

    data->menu_button = w;
    data->menu        = new_menu;

    XtInsertEventHandler (data->menu,
			  ButtonReleaseMask,
			  False,
			  Handler_ButtonRelease,
			  (XtPointer)data,
			  XtListHead);
    
    XtPopupSpringLoaded (data->menu);
    XtAddGrab (data->menubar, False, False);
  }
}


/* ARGSUSED */
static void PopupMenu(w, event, params, num_params)
     Widget     w;
     XEvent *   event;      /* unused */
     String *   params;     /* unused */
     Cardinal * num_params; /* unused */
{
  MenuButtonWidget mbw = (MenuButtonWidget) w;
  Widget           menu;
  GrabData         data;

  if (((data = GetGrabData (w)) != NULL) && (data->menu_button == w))
    return;

  if ((menu = GetPopupMenu(w)) != WNULL )
  {
    (void)CreateGrabData(XtParent(w), w, menu);

    XtInsertEventHandler(menu,
			 ButtonReleaseMask,
			 False,
			 Handler_ButtonRelease,
			 NULL,
			 XtListHead);

    XawDrawFrame (w,
		  mbw->simple.highlight_thickness,
		  mbw->simple.highlight_thickness,
	     (Dimension)(w->core.width - 2*mbw->simple.highlight_thickness),
	     (Dimension)(w->core.height - 2*mbw->simple.highlight_thickness),
	          XawRAISED,
		  mbw->simple.shadow_thickness,
		  mbw->simple.top_shadow_GC,
		  mbw->simple.bottom_shadow_GC);

    XtPopupSpringLoaded(menu);
    XtAddGrab(XtParent(w), False, False);
  }
}


/****************************************************************
 *
 * Handlers
 *
 ****************************************************************/

/* ARGSUSED */
static void Handler_ButtonRelease(w, client_data, event, need_continue)
     Widget    w;             /* menu popp-up shell */
     XtPointer client_data;   /* unused */
     XEvent*   event;         /* unused */
     Boolean*  need_continue; /* unused */
{
  GrabData data;

  if((data = GetGrabData(w)) != (GrabData)NULL)
  {
    Unhighlight(data->menu_button);
    FreeGrabData(data);
  }

  XtRemoveEventHandler(w,
		       ButtonReleaseMask,
		       False,
		       Handler_ButtonRelease,
		       NULL);
}


#if 0
static Boolean IsInRectangle(x, y, rectangle)
     int x;
     int y;
     register XRectangle* rectangle;
{
  return (x >= rectangle->x && x <= (rectangle->x + rectangle->width)
	  &&
	  y >= rectangle->y && y <= (rectangle->y + rectangle->height));
}
#endif
