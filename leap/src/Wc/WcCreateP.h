#ifndef _WcCreateP_h
#define _WcCreateP_h
#include "COPY.h"

/*
* SCCS_data: %Z% %M%	%I% %E% %U%
*
* Widget Creation Library - WcCreateP.h
*
* Private defines for the Widget Creation Library.  Clients generally
* should not need to include this file.
*
* Anything and everything in here may change dramatically.
*
*******************************************************************************
*/

#include <stdint.h>
#include <ctype.h>		/* isupper, tolower, atoi macros */
#include "WcCreate.h"
#include "MapAg.h"

BEGIN_NOT_Cxx

/*
*******************************************************************************
* Private_constant_declarations.
*******************************************************************************
*/

#undef  NUL
#define NUL '\0'
#define MAX_XRMSTRING   4096		/* max length of the Xrm DB string  */
#define MAX_ERRMSG      1024		/* max length of error message      */
#define MAX_CHILDREN    1024		/* max number of widget's children  */
#define MAX_PATHNAME    1024		/* max length of the pathname       */
#define INCR_ALLOC        32		/* initial incr of malloc'd arrays  */
#define MAX_CALLBACKS     64            /* max callbacks per Xrm resource   */
#define MAX_WIDGETS      512		/* max depth of a widget tree       */
#define MAX_ROOT_WIDGETS  32		/* max # separate widget trees	    */
#define MAX_RES_FILES    512		/* max # res file names per interf  */
#ifndef MAX_ARGS
#define MAX_ARGS	 128		/* max # args for callback / action */
#define NAME_RESOLUTION	 128		/* #chars resolving registered names*/
#endif

/*
*******************************************************************************
* Private_type_declarations.
*******************************************************************************
*/

/*  -- Wcl
*******************************************************************************
    This contains parameters which are specific to an applications use of the
    Widget Creation Library.  It is somewhat analogous to the XtAppContext.
    There is one global instance (per application which maps the library into
    its address space).
*/
typedef struct _WclRec
{
    /* Application Wide Resources - Fetched During Wcl Initialization
    */
    char*	initResFile;
    char*	resFiles;
    Boolean	traceResFiles;
    char*	errorDatabaseFile;
    char*	widgetResourceFiles;
    char*	templateFiles;
    Boolean	traceTemplateDef;
    Boolean	verboseWarnings;
    char*	dynamicLibs;		/* just to force conversion, not used */
    Boolean	slowTightNames;
    char*	appClassName;

} WclRec, *WclRecPtr, *Wcl;

extern Wcl wcl;


/*  -- QuarkRec
*******************************************************************************
    This is used by the string-to-quark converter for the wcCreate resource.
*/
typedef struct	_QuarkRec
{
    char*	string;		/* as seen in resource db		*/
    XrmQuark	quark;		/* made from lower case of string	*/
} QuarkRec;



typedef struct  _ResourceRec
{
    /* Pre-Creation Resources
    */
    Boolean         preCreateDump;	/* dump resources pre-create	*/
    XrmQuark	    template;		/* name of template		*/
    Boolean	    traceTemplate;	/* template trace required	*/
    Boolean         postTemplateDump;	/* dump resources post-template	*/
    QuarkRec	    create;		/* WcMapFind gives class/constr	*/
    QuarkRec	    className;		/* backward compatibility	*/
    QuarkRec	    class;		/* backward compatibility	*/
    QuarkRec	    constructor;	/* backward compatibility	*/

    /* Post-Creation Resources
    */
    Boolean	    trace;		/* creation trace required          */
    Boolean         postCreateDump;	/* dump resources post-create	    */
    XtCallbackList  callback;		/* invoked after creation	    */
    char*	    popups;		/* list of popup children to create */
    char*	    children;		/* list of children names to create */
    Boolean	    managed;		/* created  managed (default TRUE)  */
    XtCallbackList  afterPopups; 	/* invoked after popups created     */
    XtCallbackList  afterChildren; 	/* invoked after children created   */
    XtCallbackList  afterManageChildren;/* invoked after kids are managed   */

} ResourceRec, *ResourceRecPtr;



/* Used by WcLateBinderCB for late binding of callbacks.
*/
typedef struct _WcLateBind 
{
 XtAppContext	app;		/* Always seem to need an app!		*/
 Widget		widget;		/* widget invoking the callback		*/
 XtPointer	callData;	/* callback specific data from widget	*/

 /* These are quarkified versions of what we saw as the callback resource
  * value.  args is also specially derived from what we saw in the callback
  * resource value (careful! see WcxClosureFromSeg() for concerns about
  * being able to remove callbacks).
  */
 XrmQuark	libQ;		/* shared library abbreviation: -lXt	*/
 XrmQuark	classQ;		/* class, for methods			*/
 XrmQuark	nameQ;		/* case sensitive callback name		*/
 XrmQuark	nameq;		/* case insensitive callback name braindamage */

 /* These are needed for dynamic linking of shared libraries
 */
 char*		libFullPath;	/* if doing dynamic linking		*/
 void*		libHandle;	/* from dlopen()			*/

 /* The callback procedure address which gets invoked.
 */ 
 XtCallbackProc	Callback;	/* what finally gets invoked		*/
 char*		args;		/* careful! see WcxClosureFromSeg()	*/
 XtPointer	regClosure;	/* client data as registered		*/
 XtPointer	object;		/* instance of the class - this changes */

} WcLateBindRec, *WcLateBind;

/* This is a gross hack, simply to support WcOnceOnlyCB().  It may change...
*/
void WcLateBinder_RemoveCallback _(( char* callbackName ));

/*
*******************************************************************************
* Private_macro_definitions.
*******************************************************************************
    done() and new_done() macros used to return converted values from
    Xrm resource converters.
*/

#include "done.h"

/*  ONCE_PER_XtAppContext(app) should be invoked at the beginning of each 
    function which performs registration, like WcRegisterWcCallbacks.
    Note that this IS a macro: therefore, the return statement actually
    causes the return from the registration function.
*/
#define ONCE_PER_XtAppContext( app )	\
{					\
    static XtAppContext already[1024];	\
    static int numApps = 0;		\
    int i;				\
					\
    for (i = 0; i < numApps ; i++)	\
        if (app == already[i])		\
            return;			\
					\
    already[numApps++] = app;		\
}

/* For compatibility with old Xt libraries
*/
#ifndef XtIsWidget
#ifdef XtSpecificationRelease
#define XtIsWidget(obj) XtIsSubclass((obj),(WidgetClass)coreWidgetClass)
#else
#define XtIsWidget(obj) XtIsSubclass((obj),(WidgetClass)widgetClass)
#endif
#endif

/*
*******************************************************************************
* Private_function_declarations.
*******************************************************************************
    The following functions are generally private functions to the
    WcCreate routines, but they may be defined in different files from
    where they are used.  Client programs probably should not invoke
    these functions directly.
*/

#if NeedFunctionPrototypes
/****************************** ANSI FUNC DECLS ******************************/
#define CONVERTER(arg) XrmValue*, Cardinal*, XrmValue*, XrmValue*
#define CALLBACK(arg) Widget, XtPointer, XtPointer
#define ACTION(arg) Widget, XEvent*, char**, Cardinal*
#else
/****************************** K&R FUNC DECLS ******************************/
#define CONVERTER(arg) /**/
#define CALLBACK(arg) /**/
#define ACTION(arg) /**/
#endif

/*  -- Wcl Initialization Procedures
*/
extern void WcWarningInitialize		_(( XtAppContext, WclRecPtr ));
extern void WcWidgetResourcesInitialize	_(( XtAppContext, WclRecPtr ));
extern void WcTemplateInitialize	_(( XtAppContext, WclRecPtr ));
extern void WcRegisterIntrinsic		_(( XtAppContext ));
extern void WcRegisterWcCallbacks	_(( XtAppContext ));
extern void WcRegisterWcActions		_(( XtAppContext ));

extern int WcMoreResourceFilesToLoad _(( Widget, WclRecPtr ));

/*  -- Merge resources from file into XrmDatabase pointer to by argument,
 *     or create a new XrmDatabase and return its address in argument.
 *     Return true if database was aactually loaded.
 */
extern int WcLoadResourceFileIntoDatabase _(( Widget, char*, XrmDatabase* ));

/*  -- Wcl Templates
*/
extern int WcApplyTemplate _(( XrmQuark, Widget, char*, int ));

/*  -- Support for providing sub-part resources for widgets (safe WcSetValues)
*/
#define WcWidgetResourcesInitialize(a,w) /* not yet implemented */

/*  -- Find root widget of argument, remember if never seen before
*/
extern Widget WcRootWidget _(( Widget ));

/*  -- Use HOME etc to find user's home directory (static storage of rtn val)
*/
extern char* WcHomeDirectory _(( char* /*user*/ ));

/*  -- String to Widget Converter.  So libXmp and libXp can use easily
*/
extern int wcWidgetCvtArgsCount;
extern XtConvertArgRec wcWidgetCvtArgs[];
#ifdef XtSpecificationRelease
/* new style converter */
extern Boolean WcCvtStringToWidget _((Display*,
                                   XrmValue*, Cardinal*, XrmValue*, XrmValue*,
                                   XtPointer*));
#else
/* old style converter */
extern void WcCvtStringToWidget _((XrmValue*, Cardinal*, XrmValue*, XrmValue*));
#endif

/*  -- Callback List stuff for Add/Remove callbacks
*/
XtCallbackRec* WcStringToCallbackList _(( Widget, char* ));
void WcFreeCallbackList _(( XtCallbackRec* ));

/*  -- Perform late binding of callbacks and methods - NOT REGISTERED!
*/
extern void WcLateBinderCB _(( Widget, XtPointer, XtPointer ));

/*  -- Experimental External Interfaces to WcLateBinder
================================================================
    All callback invocations are now resolved at invocation-time by the
    late binder: no callbacks are fully resolved at widget creation time.

    The WcLateBinder gets invoked at callback invocation time for every
    callback which was added to a widget's callback list by Wcl: i.e., whenever
    Wcl "converts" a string into a callback list.  This happens when a widget
    is created and some of its callbacks are specified by name in the resource
    database, or when callbacks are added using WcSetValues, or when callbacks
    are added using WcAddCallbackCB or WcAddCallbackACT.

    The intention is to allow a client program to try to resolve callbacks.
    All the information which Wcl uses to resolve callbacks is passed to one or
    more client provided "Hook" function.  The client can provide any number of
    Hook functions.

    When a Hook function resolves a callback address sucessfully, it returns
    True, which causes Wcl to stop trying to resolve the address.

    The Hook functions get called before Wcl itself tries to resolve the
    address.  Therefore, if the Hook ALWAYS returns True, it is free to use
    the WcLateBind structure any way it wants, as the Wcl late binder will
    in that case never look at the structure contents.

    If the Hook can return False, then the Hook should treat WcLateBind
    structure fields as read only or write only:

      read only:  app, widget, callData, libQ, classQ, nameQ, nameq, args
      write only: libFullPath, libHandle, Callback, regClosure, object
 
    It is possible for the Hook to remember bindings by treating the 
    write only fields as read/write fields: the values are persistent.
    However, Wcl itself always tries to resolve the callback address and
    client data again, so the client can change the callback name to
    callback procedure mapping dynamically.  It is suggested that the
    Hooks also follow this philosophy.

    This interface will hopefully be useful to allow Wcl to be used as a
    front-end to interpretive languages and environments.  
*/
typedef Boolean (*WcLateBinderHook) _(( XtPointer hookData, WcLateBind lb ));

void WcAddLateBinderHook    _(( WcLateBinderHook Hook, XtPointer hookData ));
void WcRemoveLateBinderHook _(( WcLateBinderHook Hook, XtPointer hookData ));
 
/*  -- Similar to XtWarningMsg()
================================================================
    WcWARN* procedures use WcPrint() to actually display the messages.

    The Widget argument to identify the application context, so messages
    can be overridden from the appropriate error message database.

    1st str: name of procedure where problem occurred.
    2nd str: name of warning message.
    3rd str: default warning message.
    4th... : args to be put in for %s in warning message.

  TIP:	Wherever possible, use no more than one argument - makes
	translating or changing the messages ALOT easier.
*/
#ifndef WCL_ERRORDB
#define WCL_ERRORDB "/usr/lib/X11/WclErrorDB"
#endif

void WcWARN  _(( Widget, char*, char*, char* ));
void WcWARN1 _(( Widget, char*, char*, char*, char* ));
void WcWARN2 _(( Widget, char*, char*, char*, char*, char* ));
void WcWARN3 _(( Widget, char*, char*, char*, char*, char*, char* ));

/*  -- Override for Default Wcl Error Printing
================================================================
    By default, Wcl uses "fputs( msg, stderr )" to print all messages.
    However, many applications may want to get these messages and print
    them to a log file or display them in a Widget.

    The Wcl client software may want to assign a replacement to the WcPrint
    procedure pointer BEFORE ANY calls to Wcl are made in order to get ALL
    Wcl messages.  In this case, your prcedure should NOT assume any widgets
    have been created!

    The WcPrint procedure pointer is statically initialized, so the client
    can make a copy of this procedure pointer if the intent is to provide
    a wrapper around the default procedure.

    WcPrint() provides the end-of-message newline: i.e., Wcl warning messages
    themselves do NOT have a terminating newline.  Long messages may have
    embedded newlines.

    Note that the WcPrint function is a varargs function.  Each string is
    sent to the output device (stderr by default), finally followed by a
    newline.

    The procedure can be used like this:

	WcPrint( "Wcl Warning: WcErrorDatabaseText: ",
		 "strlen(func+msgName) > BUFSIZ ",
		 "for func = (", func, ") and msgName = (", msgName, ")",
		 NULL );
*/

typedef void (*WcPrintProc)  _(( String, ... ));

extern WcPrintProc  WcPrint;

/*  -- Buffer calls to WcPrint
================================================================
    When using Widgets to display messages (especially inefficient
    Widgets like XmText) it is useful to buffer messages.  To build
    up a single big message from lotsa little messages, make
    multiple calls to WcBuffer(), and then pass the resultant buffer
    to WcPrint(), like this:

	WcBuffer buf = WcBuffer_New();
	while ( moreToPrint ) {
	    WcBuffer_Append( buf, some, strings, to, append, NULL );
	}
	WcPrint( WcBuffer_String( buf ), NULL );
	WcBuffer_Free( buf );
*/

typedef struct _WcBuffer* WcBuffer;

extern WcBuffer WcBuffer_New();
extern void     WcBuffer_Free   _(( WcBuffer this ));
extern void     WcBuffer_Append _(( WcBuffer this, char* msg, ... ));
extern String   WcBuffer_String _(( WcBuffer this ));



/*  -- Get Error Messages from the Error Database
================================================================
    Wcl uses the Xt error database mechanism to allow overrides for
    all Wcl messages.  Clients may want to do the same.  Wcl uses a
    Widget to find the appropriate display/screen database, the
    name of the procedure or capability where the error was detected,
    and the name of the error.  For example:

	userSpecifiedMsg = WcErrorDatabaseText( widget,
						"SomeCapability",
						"someError" );
*/
    
char* WcErrorDatabaseText _(( Widget, char*, char* ));

/*  -- Mapping Agents used by Wcl
================================================================
    Need to have these visible only so macros below can be used.
*/
extern MapAg dlAgent, cbAgent, cdAgent, clAgent, conAgent, temAgent, rfAgent;

/*  -- Mapping Agent Access Macros
=====================================================================
    Used for consistent access to mapping agents.  If a single agent is
    used for more than one type of data, be certain that the arguments
    do NOT collide!
*/

/*============== dlAgent ==============*/
#define WcMapDynLib( app, quark, name ) \
	MapAg_Define( dlAgent, (app), (quark), 1, (name) )

#define WcMapDynLibFind( app, quark ) \
	(char*)MapAg_Find( dlAgent, (app), (quark), 1 )

/*============== cbAgent ==============*/
#define WcMapCallback( app, quark, cbRecPtr ) \
	MapAg_Define( cbAgent, (app), (quark), NULL, (cbRecPtr) )

#define WcMapCallbackFind( app, quark ) \
	(XtCallbackRec*)MapAg_Find( cbAgent, (app), (quark), NULL )

#define WcMapCallbackMethod( app, classQ, nameQ, cbRecPtr ) \
	MapAg_Define( cbAgent, (app), (classQ), (nameQ), (cbRecPtr) )

#define WcMapCallbackMethodFind( app, classQ, nameQ ) \
	(XtCallbackRec*)MapAg_Find( cbAgent, (app), (classQ), (nameQ) )

/*============== cdAgent ==============*/
#define WcMapClosure( quark, string ) \
	MapAg_Define( cdAgent, (quark), NULL, NULL, (string) )

#define WcMapClosureFind( quark ) \
	(char*)MapAg_Find( cdAgent, (quark), NULL, NULL )

#define WcMapObject( wid, classQ, object ) \
	MapAg_Define( cdAgent, (wid), (classQ), NULL, (object) )

#define WcMapObjectFind( wid, classQ ) \
	(XtPointer)MapAg_Find( cdAgent, (wid), (classQ), NULL )

#define WcMapObjectForget( wid, classQ ) \
	MapAg_Forget( cdAgent, (wid), (classQ), NULL )

#define WcMapObjectPurge( object ) \
	MapAg_Purge( cdAgent, object )

/*============== conAgent ==============*/
#define WcMapConstructor( app, quark, Constr ) \
	MapAg_Define( conAgent, (app), (quark), NULL, (Constr) )

#define WcMapConstructorFind( app, quark ) \
	(WcWidgetConstructor)MapAg_Find( conAgent, (app), (quark), NULL )

/*============== clAgent ==============*/
#define WcMapClass( app, quark, class ) \
	MapAg_Define( clAgent, (app), (quark), NULL, (class) )

#define WcMapClassFind( app, quark ) \
	(WidgetClass)MapAg_Find( clAgent, (app), (quark), NULL )

/*============== temAgent -- Map Template Names to Template XrmDatabase =====*/
#define WcMapTemplateFileLoaded( name ) \
	MapAg_Define( temAgent, NULL, name, NULL, 1 )

#define WcMapTemplateFileAlreadyLoaded( name ) \
	(intptr_t)MapAg_Find( temAgent, NULL, name, NULL )

#define WcMapTemplateNameToDatabaseDefine( name, db ) \
	MapAg_Define( temAgent, name, NULL, NULL, db )

#define WcMapTemplateNameToDatabase( name ) \
	(XrmDatabase)MapAg_Find( temAgent, name, NULL, NULL )

END_NOT_Cxx

#endif /* _WcCreateP_h */
