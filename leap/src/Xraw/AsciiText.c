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

/***********************************************************************
 *
 * AsciiText Widget
 *
 ***********************************************************************/

/*
 * AsciiText.c - Source code for AsciiText Widget.
 *
 * This Widget is intended to be used as a simple front end to the 
 * text widget with an ascii source and ascii sink attached to it.
 *
 * Date:    June 29, 1989
 *
 * By:      Chris D. Peterson
 *          MIT X Consortium 
 *          kit@expo.lcs.mit.edu
 */

#include <stdio.h>
#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>

#include "XawInit.h"
#include "AsciiSrc.h"
#include "AsciiSink.h"
#include "AsciiTextP.h"
#include "Cardinals.h"

#include "XrawDebug.h"


#define TAB_COUNT 32

static void ClassInitialize(), Initialize(), CreateSourceSink(), Destroy();
static Boolean SetValues();

AsciiTextClassRec asciiTextClassRec = {
  { /* core fields */
    /* superclass       */      (WidgetClass) &textClassRec,
    /* class_name       */      "Text",
    /* widget_size      */      sizeof(AsciiRec),
    /* class_initialize */      ClassInitialize,
    /* class_part_init  */	NULL,
    /* class_inited     */      FALSE,
    /* initialize       */      Initialize,
    /* initialize_hook  */	CreateSourceSink,
    /* realize          */      XtInheritRealize,
    /* actions          */      textActionsTable,
    /* num_actions      */      0, /* let's look in ClassInitialize */
    /* resources        */      NULL,
    /* num_ resource    */      0,
    /* xrm_class        */      NULLQUARK,
    /* compress_motion  */      TRUE,
    /* compress_exposure*/      XtExposeGraphicsExpose,
    /* compress_enterleave*/	TRUE,
    /* visible_interest */      FALSE,
    /* destroy          */      Destroy,
    /* resize           */      XtInheritResize,
    /* expose           */      XtInheritExpose,
    /* set_values       */      SetValues,
    /* set_values_hook  */	NULL,
    /* set_values_almost*/	XtInheritSetValuesAlmost,
    /* get_values_hook  */	NULL,
    /* accept_focus     */      XtInheritAcceptFocus,
    /* version          */	XtVersion,
    /* callback_private */      NULL,
    /* tm_table         */      XtInheritTranslations,
    /* query_geometry	*/	XtInheritQueryGeometry
  },
  { /* Simple fields */
    /* change_sensitive	*/	XtInheritChangeSensitive,
    /* display_rect	*/	XtInheritDisplayRectProc,
    /* extension	*/	NULL				
  },
  { /* text fields */
    /* empty            */      0
  },
  { /* ascii fields */
    /* empty            */      0
  }
};

WidgetClass asciiTextWidgetClass = (WidgetClass)&asciiTextClassRec;

static void 
ClassInitialize()
{
  XawInitializeWidgetSet();

  asciiTextClassRec.core_class.num_actions = textActionsTableCount;

#ifdef notdef
#if defined(XtSpecificationRelease) && XtSpecificationRelease < 5
  asciiTextClassRec.core_class.num_actions = textActionsTableCount;
#else
  XtInitializeWidgetClass(textWidgetClass);
  asciiTextClassRec.core_class.num_actions =
    textClassRec.core_class.num_actions;
  asciiTextClassRec.core_class.actions = textClassRec.core_class.actions;
#endif
#endif

}

/* ARGSUSED */
static void
Initialize(request, new)
Widget request, new;
{
  /* superclass Initialize can't set the following,
   * as it didn't know the source or sink when it was called */
  if (request->core.height == DEFAULT_TEXT_HEIGHT)
    new->core.height = DEFAULT_TEXT_HEIGHT;

}

static void 
CreateSourceSink(widget, args, num_args)
Widget widget;
ArgList args;
Cardinal *num_args;
{
  AsciiWidget w = (AsciiWidget) widget;
  int i;
  int tabs[TAB_COUNT], tab;
  
  w->text.source = XtCreateWidget( "textSource", asciiSrcObjectClass,
				   widget, args, *num_args );
  w->text.sink = XtCreateWidget( "textSink", asciiSinkObjectClass,
				 widget, args, *num_args );

  if (w->core.height == DEFAULT_TEXT_HEIGHT)
    w->core.height = VMargins(w) + XawTextSinkMaxHeight(w->text.sink, 1);

  for (i=0, tab=0 ; i < TAB_COUNT ; i++) 
    tabs[i] = (tab += 8);
  
  XawTextSinkSetTabs(w->text.sink, TAB_COUNT, tabs);

  XawTextDisableRedisplay(widget);
  XawTextEnableRedisplay(widget);
}

static void 
Destroy(w)
Widget w;
{
  XtDestroyWidget( ((AsciiWidget)w)->text.source);
  XtDestroyWidget( ((AsciiWidget)w)->text.sink );
}

#ifdef ASCII_STRING

/************************************************************
 *
 * Ascii String Compatibility Code.
 *
 ************************************************************/

AsciiStringClassRec asciiStringClassRec = {
  { /* core fields */
    /* superclass       */      (WidgetClass) &asciiTextClassRec,
    /* class_name       */      "Text",
    /* widget_size      */      sizeof(AsciiStringRec),
    /* class_initialize */      NULL,
    /* class_part_init  */	NULL,
    /* class_inited     */      FALSE,
    /* initialize       */      NULL,
    /* initialize_hook  */	NULL,
    /* realize          */      XtInheritRealize,
    /* actions          */      textActionsTable,
    /* num_actions      */      0,
    /* resources        */      NULL,
    /* num_ resource    */      0,
    /* xrm_class        */      NULLQUARK,
    /* compress_motion  */      TRUE,
    /* compress_exposure*/      XtExposeGraphicsExpose,
    /* compress_enterleave*/	TRUE,
    /* visible_interest */      FALSE,
    /* destroy          */      NULL,
    /* resize           */      XtInheritResize,
    /* expose           */      XtInheritExpose,
    /* set_values       */      NULL,
    /* set_values_hook  */	NULL,
    /* set_values_almost*/	XtInheritSetValuesAlmost,
    /* get_values_hook  */	NULL,
    /* accept_focus     */      XtInheritAcceptFocus,
    /* version          */	XtVersion,
    /* callback_private */      NULL,
    /* tm_table         */      XtInheritTranslations,
    /* query_geometry	*/	XtInheritQueryGeometry
  },
  { /* Simple fields */
    /* change_sensitive	*/	XtInheritDisplayRectProc
  },
  { /* text fields */
    /* empty            */      0
  },
  { /* ascii fields */
    /* empty            */      0
  }
};

WidgetClass asciiStringWidgetClass = (WidgetClass)&asciiStringClassRec;

#endif /* ASCII_STRING */

#ifdef ASCII_DISK

/************************************************************
 *
 * Ascii Disk Compatibility Code.
 *
 ************************************************************/

AsciiDiskClassRec asciiDiskClassRec = {
  { /* core fields */
    /* superclass       */      (WidgetClass) &asciiTextClassRec,
    /* class_name       */      "Text",
    /* widget_size      */      sizeof(AsciiDiskRec),
    /* class_initialize */      NULL,
    /* class_part_init  */	NULL,
    /* class_inited     */      FALSE,
    /* initialize       */      NULL,
    /* initialize_hook  */	NULL,
    /* realize          */      XtInheritRealize,
    /* actions          */      textActionsTable,
    /* num_actions      */      0,
    /* resources        */      NULL,
    /* num_ resource    */      0,
    /* xrm_class        */      NULLQUARK,
    /* compress_motion  */      TRUE,
    /* compress_exposure*/      XtExposeGraphicsExpose,
    /* compress_enterleave*/	TRUE,
    /* visible_interest */      FALSE,
    /* destroy          */      NULL,
    /* resize           */      XtInheritResize,
    /* expose           */      XtInheritExpose,
    /* set_values       */      NULL,
    /* set_values_hook  */	NULL,
    /* set_values_almost*/	XtInheritSetValuesAlmost,
    /* get_values_hook  */	NULL,
    /* accept_focus     */      XtInheritAcceptFocus,
    /* version          */	XtVersion,
    /* callback_private */      NULL,
    /* tm_table         */      XtInheritTranslations,
    /* query_geometry	*/	XtInheritQueryGeometry
  },
  { /* Simple fields */
    /* change_sensitive	*/	XtInheritDisplayRectProc
  },
  { /* text fields */
    /* empty            */      0
  },
  { /* ascii fields */
    /* empty            */      0
  }
};

WidgetClass asciiDiskWidgetClass = (WidgetClass)&asciiDiskClassRec;

#endif /* ASCII_DISK */

/*ARGSUSED*/
static Boolean 
SetValues(current, request, new, args, num_args)
Widget current, request, new;
ArgList args;
Cardinal *num_args;
{
  AsciiWidget oldtw = (AsciiWidget) current;
  AsciiWidget newtw = (AsciiWidget) new;
  Boolean    redisplay = FALSE;
  
  if (oldtw->core.background_pixel != newtw->core.background_pixel) {
    XtVaSetValues(newtw->text.sink,
		  "background",newtw->core.background_pixel, NULL);
    redisplay = TRUE;
  }

  return redisplay;
  
}  
