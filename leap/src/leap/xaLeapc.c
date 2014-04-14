/*
 *	File:	xaLeapc.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *	Description:
 *		Front end for LEaP under X-Windows using the ATHENA WIDGET set.
 *		Redirect the output using PushCurrentPrintSink.
 *		to a routine that will append the output
 *		to the Text Widget.
 */

#include        <stdio.h>
#include        <stdlib.h>

#include 	<X11/IntrinsicP.h>
#ifdef sun
#include <X11/ObjectP.h>	/* why don't they just use X from mit!?! */
#include <X11/RectObjP.h>
#endif
#include	<X11/StringDefs.h>
#include 	<X11/Shell.h>
#include	<X11/keysym.h>
#include        <X11/cursorfont.h>


#include 	"../Wc/WcCreate.h"
#include	"../Xpm/xpm.h"

#ifdef FUNC
#undef FUNC
#endif

#include	"basics.h"
  
#include        "xaLeap.h"
#include	"classes.h"
#include	"xaTools.h"
#include	"xTank.h"
#include	"xaUnitEditor.h"
#include	"xaCommand.h"
#include	"xaTable.h"

#include	"parser.h"
#include        "block.h"
#include        "../Xraw/Command.h"
#include        "../Xraw/3d.h"
#define ICON_NAME ("xleap_icon")

extern VFUNCTION	GfAtomClassGraphicsCreator;
extern VFUNCTION	GfAtomClassGraphicsDestructor;
extern int isatty();

//static int error_handler();

void ParseInit( RESULTt *rPResult );
void ParseArguments( int argc, char *argv[] );
void ParseShutdown();
void XAPEInitialize( XtAppContext app );

/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 *
 *	Global variables
 *
 *	Maintain a pointer to the Widget which maintains the
 *	Text Widget which acts as the CommandLine interface.
 */

Widget		GwTopWidget;		/* Top widget */
Widget		GwCmd;			/* Command line interface */
XtAppContext	GxtacApp;
int		GiCmdSink;
Display		*GdDisplay;
XrmDatabase	GXrmDbase, GXrmMyDbase;
Screen		*GsScreen;
Window		GwRootWindow;
GC              GgcDark;
GC              GgcLight;
FILE           	*GfInput;

static
void SetCommandColoursCB( Widget w, XtPointer client_data, XtPointer call_data)
{
  char	*param = (char*)client_data;
  Pixel  background = 0;
  Widget to_set;
  STRING name;

  if (param != NULL && *param != '\0') {
    param = WcCleanName( param, name);
    to_set = (Widget)WcFullNameToWidget (w, name);

    if (param != NULL && *param != '\0' && to_set != (Widget)NULL) {
      param = WcCleanName( param, name);
      (void) FetchPixel (w, name, &background);
      XawSetCommandColours (to_set, background);
    }
  }
  
}

static
void SetSensitiveCB( Widget w, XtPointer client_data, XtPointer call_data)
{
  char		*param = (char*)client_data;
  Widget	to_set;
  STRING	name;

  while ( param != NULL && *param != '\0') {
    
    param = WcCleanName( param, name);
    to_set = GetWidgetByName( name);
    
    if ( IS_GADGET(to_set) )
      XtSetSensitive( to_set, True);
  }
  
}

static
void SetInsensitiveCB( Widget w, XtPointer client_data, XtPointer call_data)
{
  char		*param = (char*)client_data;
  Widget	to_set;
  STRING	name;

  while ( param != NULL && *param != '\0') {
    
    param = WcCleanName( param, name);
    to_set = GetWidgetByName( name);
    
    if ( IS_GADGET(to_set) )
      XtSetSensitive( to_set, False);
  }

}

static void 
SetParentShellAsGroup( Widget wWid, caddr_t client_data, caddr_t call_data )
{
  char		*param = (char*)client_data;
  STRING        name;
  Widget        wTop;
  Boolean       trans;
  Window        group;
  
  if ( param != NULL && *param != '\0') {

    param = WcCleanName( param, name);
    wTop = (Widget)WcFullNameToWidget( wWid, name); 
    
    if ( XtIsShell(wTop) ) {
      group = XtWindow(wTop);
      trans = True;
      
      XtVaSetValues( wWid,
		    XtNtransient, (XtArgVal) trans, 
		    XtNwindowGroup, (XtArgVal) group,
		    XtNtransientFor, (XtArgVal) wTop, 
		    NULL );
    }

  }

}

static void
SetIconForShell( Widget wWid, caddr_t client_data, caddr_t call_data )
{
  char		*param = (char*)client_data;
  STRING        name;
  
  if ( param != NULL && *param != '\0' && XtIsShell(wWid)) {

    XWMHints  xwmhHints;
    XrmValue source, dest;
    Pixmap   pixmap;
    
    param = WcCleanName( param, name);
    
    source.size = strlen( name) + 1;
    source.addr = name;
    
    dest.size = sizeof(Pixmap);
    dest.addr = (caddr_t)&pixmap;
    
    if (!XtConvertAndStore( wWid, XtRString, &source, XtRPixmap, &dest))
      PRINTF(( "Icon pixmap \"%s\" is not accessible \n", name ));
    
    xwmhHints.flags       = IconPixmapHint;
    xwmhHints.icon_pixmap = pixmap;

    if (XtWindow( wWid))
      XSetWMHints( GdDisplay, XtWindow( wWid), &xwmhHints );
    else {
      PRINTF(( "Icon pixmap \"%s\", window is absent \n", name ));
      XtVaSetValues( wWid, XtNiconPixmap, (XtArgVal) pixmap, NULL);
      
    }
  }
  
}

Widget
XAGetWidgetFromString( Widget wWid, char **string)
{
  STRING widget_name;

  if ( (*string) == (char*)NULL || (**string) == '\0' )
    return (Widget)NULL;
  
  (*string) = WcCleanName( (*string), widget_name);
  
  return WcFullNameToWidget( wWid, widget_name);
  
}


static void
XASetKeyboardFocus( Widget wWid, caddr_t client_data, caddr_t call_data )
{
  Widget shell = XAGetWidgetFromString (wWid, (char**)&client_data);

  XtSetKeyboardFocus (shell, wWid);
}

static void
XAClearText(Widget wWid, caddr_t client_data, caddr_t call_data )
{
  Widget widget = XAGetWidgetFromString (wWid, (char**)&client_data);

  XtVaSetValues (widget, XtNstring, "", NULL);
}

static void
XASetPopdown( Widget wWid, caddr_t client_data, caddr_t call_data )
{
  static XtPopdownIDRec popdown;
  char	 *string = (char*)client_data;
  Widget  wDismiss;
  
  popdown.shell_widget  = XAGetWidgetFromString( wWid, &string);
  popdown.enable_widget = XAGetWidgetFromString( wWid, &string);
  wDismiss              = XAGetWidgetFromString( wWid, &string);

  if ( popdown.shell_widget == NULL ||
      popdown.enable_widget == NULL ||
      wDismiss == NULL){
    PRINTF(("XASetPopdown requires 3 params (shell enable dismiss).\n")); 
    
  }

  XtAddCallback( wDismiss, XtNcallback, XtCallbackPopdown, (caddr_t) &popdown);
  
}


/*
 *	zXaLeapCallback
 *
 *	Process all pending events.
 */
void	zXaLeapCallback()
{
    XSynchronize( GdDisplay, True );
}


#define	done(type, value) \
	{							\
	    if (toVal->addr != NULL) {				\
		if (toVal->size < sizeof(type)) {		\
		    toVal->size = sizeof(type);			\
		    return False;				\
		}						\
		*(type*)(toVal->addr) = (value);		\
	    }							\
	    else {						\
		static type static_val;				\
		static_val = (value);				\
		toVal->addr = (XtPointer)&static_val;		\
	    }							\
	    toVal->size = sizeof(type);				\
	    return True;					\
	}


/*ARGSUSED*/
static Boolean
XACvtStringToBitmap(Display *dpy, XrmValuePtr args, Cardinal *num_args, 
	XrmValuePtr fromVal, XrmValuePtr toVal, XtPointer *data)
{
  static Pixmap pixmap;		/*** static for cvt magic ***/
  Window        drawable = XDefaultRootWindow(GdDisplay);
  char		*system_env;
  char         *name     = (char *)fromVal->addr;
  static char   sufix[] = ".xbm";
  char          sepr     = ':';
  char          file_name[200];
  unsigned int  w, h;
  int           x_hot,y_hot;

  if (*num_args != 0)
    XtErrorMsg("wrongParameters","cvtStringToBitmap","XtToolkitError",
	       "String to Pixmap conversion needs no any arguments",
	       (String *)NULL, (Cardinal *)NULL);
  
  system_env = (char *) getenv("PIXMAP_PATH");
  
  if (system_env != NULL ) {
    int   i;
    char  *ch;
    
    i = 0;
    for(ch=system_env;  *ch != '\0' ;ch++) {
      if ( *ch != sepr ){
	file_name[i++] = *ch;
      } else {
	file_name[i++] = '/';
	file_name[i=0] = '\0'; 
	strcat(file_name, name);	

	(void)XReadBitmapFile(GdDisplay, drawable, file_name, &w,&h,
			      &pixmap, &x_hot, &y_hot);
	
	if (pixmap != None)
	  break;

	strcat(file_name, sufix);

	(void)XReadBitmapFile(GdDisplay, drawable, file_name, &w,&h,
			      &pixmap, &x_hot, &y_hot);
	
	if (pixmap != None)
	  break;
	
      }
    }
    
    if ( *ch == '\0' ) {
      file_name[i++] = '/';
      file_name[i]   = '\0';
      strcat(file_name, name);

      (void)XReadBitmapFile(GdDisplay, drawable, file_name, &w,&h,
			      &pixmap, &x_hot, &y_hot);
      
      if ( pixmap == None ) {
	
	strcat(file_name, sufix);
	
	(void)XReadBitmapFile(GdDisplay, drawable, file_name, &w,&h,
			      &pixmap, &x_hot, &y_hot);

      }
      
    }
    
  } else {
    
    system_env = (char *) getenv("AMBERHOME");
    file_name[0] = '\0';
    strcat(file_name, system_env);
    strcat(file_name, "/dat/pixmap/");
    strcat(file_name, name);

    (void)XReadBitmapFile(GdDisplay, drawable, file_name, &w,&h,
			  &pixmap, &x_hot, &y_hot);

    if ( pixmap == None ) {
      
      strcat(file_name, sufix);
      
      (void)XReadBitmapFile(GdDisplay, drawable, file_name, &w,&h,
			    &pixmap, &x_hot, &y_hot);
      
    }

  }
  
  if (pixmap != None) {

    done(Pixmap, pixmap)
	    
  } else {

    XtStringConversionWarning (name, "Pixmap");
    return False;

  }

}
static Boolean  
XACvtStringToPixmap(Display *dpy, XrmValuePtr args, Cardinal *num_args, 
	XrmValuePtr fromVal, XrmValuePtr toVal, XtPointer *data)
{
  static Pixmap pixmap;		/*** static for cvt magic ***/
  static Pixmap shape_pixmap;
  Window        drawable = XDefaultRootWindow(GdDisplay);
  char		*system_env;
  char         *name     = (char *)fromVal->addr;
  char		*sufix   = ".xpm";
  char          sepr     = ':';
  char          file_name[260];


  if (*num_args != 0)
    XtErrorMsg("wrongParameters","cvtStringToPixmap","XtToolkitError",
	       "String to Pixmap conversion needs no any arguments",
	       (String *)NULL, (Cardinal *)NULL);
  
  system_env = (char *) getenv("PIXMAP_PATH");
  
  if (system_env != NULL ) {
    int   i;
    char  *ch;
    
    i = 0;
    for(ch=system_env;  *ch != '\0' ;ch++) {
      if ( *ch != sepr ){
	file_name[i++] = *ch;
      } else {
	file_name[i++] = '/';
	file_name[i]   = '\0'; i = 0;
	strcat(file_name, name);	

	(void)XpmReadFileToPixmap(GdDisplay, drawable,
				  file_name, &pixmap, &shape_pixmap, NULL);
	
	if (pixmap != None)
	  break;

	strcat(file_name, sufix);

	(void)XpmReadFileToPixmap(GdDisplay, drawable,
				  file_name, &pixmap, &shape_pixmap, NULL);
	
	if (pixmap != None)
	  break;
	
      }
    }
    
    if ( *ch == '\0' ) {
      file_name[i++] = '/';
      file_name[i]   = '\0';
      strcat(file_name, name);

      (void)XpmReadFileToPixmap(GdDisplay, drawable,
				file_name, &pixmap, &shape_pixmap, NULL);
      
      if ( pixmap == None ) {
	
	strcat(file_name, sufix);
	
	(void)XpmReadFileToPixmap(GdDisplay, drawable,
				  file_name, &pixmap, &shape_pixmap, NULL);
      }
      
    }
    
  } else {
    
    system_env = (char *) getenv("AMBERHOME");
    if (system_env == NULL) {
	fprintf(stderr, "`AMBERHOME' environment not set\n");
	exit(1);
    }
    file_name[0] = '\0';
    strcat(file_name, system_env);
    strcat(file_name, "/dat/pixmap/");
    strcat(file_name, name);
    
    (void)XpmReadFileToPixmap(GdDisplay, drawable,
			      file_name, &pixmap, &shape_pixmap, NULL);

    if ( pixmap == None ) {
      
      strcat(file_name, sufix);
      
      (void)XpmReadFileToPixmap(GdDisplay, drawable,
				file_name, &pixmap, &shape_pixmap, NULL);
    }

  }
  
  if (pixmap != None) 
  {
    done(Pixmap, pixmap)
  } 
  else 
  {
    return XACvtStringToBitmap(dpy,args, num_args, fromVal, toVal, data);
  }

}



static void
XAPixmapDestructor(XtAppContext app, XrmValuePtr to, XtPointer converter_data,
			       XrmValuePtr args, Cardinal *num_args)
{
  XFreePixmap(GdDisplay,*(Pixmap*)to->addr);

}
     

static void 
RegisterXACvtStringToPixmap()
{

  XtSetTypeConverter("String", "Pixmap", XACvtStringToPixmap,
                    (XtConvertArgList)NULL, 0,
		     XtCacheAll,XAPixmapDestructor );
  XtSetTypeConverter("String", "Bitmap", XACvtStringToBitmap,
                    (XtConvertArgList)NULL, 0,
		     XtCacheAll,XAPixmapDestructor );

}

static void
zXaLeapRegister()
{

#define RCP( name, class ) WcRegisterClassPtr  ( GxtacApp, name, class )
#define RCC( name, cback ) WcRegisterCallback  ( GxtacApp, name, (XtCallbackProc)cback , NULL)

  /*** Register the TANK class ***/
  RCP("Tank",                          tankWidgetClass );
  RCP("tankWidgetClass",               tankWidgetClass );
  
  /*** Register global (instead of Wcl local) sensitive callbacks  ***/

  RCC( "XASetPopdown",                 XASetPopdown                  );

  RCC( "SetSensitiveCB",               SetSensitiveCB                );
  RCC( "SetSensitive",                 SetSensitiveCB                );
  RCC( "SetInsensitiveCB",             SetInsensitiveCB              );
  RCC( "SetInsensitive",               SetInsensitiveCB              );
  RCC( "SetParentShellAsGroup",        SetParentShellAsGroup         );
  RCC( "SetIconForShell",              SetIconForShell               );
  RCC( "SetCommandColoursCB",          SetCommandColoursCB           );

  RCC( "XASetKeyboardFocus",           XASetKeyboardFocus            );
  RCC( "XAClearText",                  XAClearText                   );

#undef RCP  
#undef RCC
}

static void
XAInputAvailable( XtPointer xtPData, int *iPSource, XtInputId *xtiiPId )
{
STRING          sBuffer;
char		*cPRead;

  cPRead = fgets( sBuffer, sizeof(STRING), GfInput );
  if ( !cPRead ) {
    return;
  }

  if ( strlen(sBuffer) > 0 ) {
    PushCurrentPrintSink(GiCmdSink);

    for(; *cPRead; cPRead++) 
      zXACTextProcessOneChar( *cPRead);

    PopCurrentPrintSink();
  }

}


/******************************************************************************
 *   MAIN function
 *
 *	Author:	Christian Schafmeister (1991)
 ******************************************************************************/
extern int itest;
int 
main( int argc, char *argv[] )
{   
  STRING    sAppClass;
  Widget    wAppShell;
  int       i, iLastSlash;
  RESULTt   rResult;

IMem();

  BasicsInitialize();
  /*
   *  Set up the graphics-specific stuff
   */
  
  GbGraphicalEnvironment = TRUE;
  
  /*** Define the ATOM graphics data creator ***/
  /*** destructor before ANY ATOMs are created ***/
  
  GfAtomClassGraphicsCreator = TankAtomGraphicsDataCreator;
  GfAtomClassGraphicsDestructor = TankAtomGraphicsDataDestructor;
  
  /*** Initialize X-Windows stuff ***/
  
  iLastSlash = -1;

  for ( i=0; i<strlen(argv[0]); i++ ) 
    if ( argv[0][i] == '/' )
      iLastSlash = i;

  for ( i=iLastSlash+1; i<strlen(argv[0]); i++ ) 
    GsProgramName[i-iLastSlash-1] = argv[0][i];

  GsProgramName[strlen(argv[0])-iLastSlash-1] = '\0';
  
  sprintf( sAppClass, "%s", GsProgramName );
  
  /***
   * initialize first letter to make class, or first two if
   * first is already capitalized, or don't worry about it.
   ***/
  
    if (islower(sAppClass[0]))
      sAppClass[0] = cUpper(sAppClass[0]);
    else if (islower(sAppClass[1]))
      sAppClass[1] = cUpper(sAppClass[1]);
  
  /***   Initialize Toolkit creating the application shell ***/

  wAppShell = XtInitialize(GsProgramName, sAppClass,
			   (XrmOptionDescRec*)NULL, (Cardinal)0,
			   &argc, argv );

  GxtacApp     = XtWidgetToApplicationContext(wAppShell);
  GdDisplay    = XtDisplay(wAppShell);
/*
This is for doing away w/ XaLeap resources file if
it can ever be packed into a string that the compiler
will accept
  GXrmDbase    = XtDatabase(GdDisplay);
  GXrmMyDbase  = XrmGetStringDatabase(GsResources);
  XrmMergeDatabases( GXrmMyDbase, &GXrmDbase );
*/
  GsScreen     = DefaultScreenOfDisplay(GdDisplay);
  GwRootWindow = RootWindowOfScreen(GsScreen);
  GwTopWidget  = wAppShell;
  
#ifdef	DEBUG
    PRINTF(( "Program name= %s, application class= %s\n", GsProgramName,
		sAppClass ));
#endif

    
  /*** Parse command line arguments ***/
  
  ParseArguments( argc, argv );
  
  /*** Register new widgets and global callbacks ***/

  zXaLeapRegister();
  
  /***  Register converter for color pixmap ***/
    
  RegisterXACvtStringToPixmap();
  
  /*** Register Callbacks for xaCommand ***/
  /*** Register all of the callbacks for xaUnitEditor ***/
  /*** Register Callbacks for the xaTable ***/

  
  XATInitializeTools( GxtacApp);
  XACInitialize( GxtacApp );
  XAUEInitialize( GxtacApp );
  XATInitialize( GxtacApp,GwTopWidget );

  
  /***	Davids Changes ***/
  /*** initialize the parmeditor ***/

  XAPEInitialize( GxtacApp );

  /***	End of Davids Changes ***/
  
  
  /***   Register all Athena and Public widget classes ***/
  
  XrawRegisterAll( GxtacApp );
  
  /***   Create widget tree below toplevel shell using Xrm database ***/

  WcWidgetCreation( wAppShell );
  
  /***   Realize the widget tree and enter the main application loop ***/

  XtRealizeWidget( wAppShell );
  
  /***  Set up the window manager hints ***/

  {
    XWMHints  xwmhHints;
    XrmValue source, dest;
    static Pixmap   pixmap;

    source.size = strlen(ICON_NAME) + 1;
    source.addr = ICON_NAME;
    
    dest.size = sizeof(Pixmap);
    dest.addr = (caddr_t)&pixmap;

    if (!XtConvertAndStore( wAppShell, XtRString, &source, XtRPixmap, &dest)){
      PRINTF(( "Icon pixmap \"%s\" is not accessible \n", ICON_NAME ));
    } else {
    
      xwmhHints.flags       = InputHint | IconPixmapHint;
      xwmhHints.icon_pixmap = pixmap;
      xwmhHints.input       = TRUE;
      
      XSetWMHints( GdDisplay, XtWindow(wAppShell), &xwmhHints );
      XFlush( GdDisplay );
    }
  }
  
  /*** Create the print sink for the command line interface ***/
  
  GiCmdSink = iCreatePrintSink( XATPrintStringToWidget, "", (GENP)GwCmd );

  GgcLight = XCreateGC(GdDisplay, GwRootWindow, 0, 0);
  GgcDark  = XCreateGC(GdDisplay, GwRootWindow, 0, 0);
  
  PushCurrentPrintSink(GiCmdSink);
  
  /*** Check if Num Lock is on and print a message if it is ***/
  XModifierKeymap *map = XGetModifierMapping(GdDisplay);
  KeyCode code = XKeysymToKeycode(GdDisplay, XK_Num_Lock);
  Window w1; int i1;
  unsigned int iMask;
  int iKeyMask = 0;
  
  for (i=0; i<8; ++i)
    if (map->modifiermap[map->max_keypermod*i] == code)
      iKeyMask = 1 << i;
    XQueryPointer(GdDisplay, DefaultRootWindow(GdDisplay), &w1, &w1, &i1, &i1, &i1, &i1, &iMask);
  if (iMask & iKeyMask) VPDISPLAY( ("Please turn Num Lock off for the menus to function!\n"));
  XFreeModifiermap(map);
  
  /*** Initialize the parser after all the widgets ***/
  /*** have been created ***/

  ParseInit( &rResult );

  setbuf(stdout,NULL);

  if (!isatty(fileno(stdin))) {
    GfInput = stdin;
    (void) XtAppAddInput( GxtacApp, fileno(stdin),
		         (XtPointer)XtInputReadMask,
		         (XtInputCallbackProc)XAInputAvailable, 
			 (XtPointer)NULL );
  }

#if 0
  XSynchronize( XtDisplay(wAppShell), True );
  (void)XSetErrorHandler(error_handler);
#endif

  if ( rResult.iCommand != CQUIT ) {
    VPDISPLAY(( "> " ));	/*** Print the initial prompt ***/
    
    PopCurrentPrintSink();

    XtMainLoop ( ); 

  }
  
  exit(0);

}






#if 0
static int 
error_handler(Display *dpy, XErrorEvent *err)
{
  char *p = NULL;
  /* Segmentation fault is needed to get coredump */
  *p = ' ';
  return 1;
}
#endif
