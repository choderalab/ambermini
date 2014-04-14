#include "COPY.h"

/*
* SCCS_data: %Z% %M% %I% %E% %U%
*
* Widget Creation Library - WcWarn.c
*
* Implements simple procedures for printing warning messages similar to
* XtWarningMsg().  The number or string substitutions are checked for
* run-time safety, and no fixed sized buffers are used (instead, fprintf()
* is called multiple times which may be inefficient, but who cares for
* warning messages!).
*
* Error messages can be changed on a system wide basis be changing the Wcl
* error database file.  Note the way messages can be changed:
*	<func>.<msgName>:	<msg>
*
*******************************************************************************
*/

#include <X11/Intrinsic.h>
#include "WcCreateP.h"

#ifdef _NO_VSPRINTF
static void vsprintf _(( String buf, String fmt, va_list argPtr ));

static void vsprintf( buf, fmt, argPtr )
    String	buf, fmt;
    va_list	argPtr;
{
    strcpy( buf, fmt );
}
#endif /*_NO_VSPRINTF*/


/* Wcl Print Error Default Procedure Pointer
**===========================================**
*/

static void WcxPrint  _(( String, ... ));

WcPrintProc  WcPrint  = WcxPrint;

/* Wcl Xrm Error Database
**========================**
   It is OK if wcxErrorDB is left NULL.
*/
static XrmDatabase wcxErrorDB;

/* save calls to scanf
**=====================**
*/
static String WcxNumToStr _((int));

static String WcxNumToStr( num )
    int num;
{
    switch ( num )
    {
	case 0:	 return "0";
	case 1:	 return "1";
	case 2:	 return "2";
	case 3:	 return "3";
	case 4:	 return "4";
	case 5:	 return "5";
	case 6:	 return "6";
	case 7:	 return "7";
	case 8:	 return "8";
	case 9:	 return "9";
	case 10: return "10";
	default: return "?";
    }
}



/* Initialize ErrorDB
**====================**
*/
/*ARGSUSED*/
void WcWarningInitialize( app, wcl )
    XtAppContext	app;
    WclRecPtr		wcl;
{
    if ( (String)0 != wcl->errorDatabaseFile )
    {
	wcxErrorDB = XrmGetFileDatabase( wcl->errorDatabaseFile );
    }
}

/* Get Error Message From ErrorDB
**================================**
*/
/*ARGSUSED*/
String WcErrorDatabaseText( w, func, msgName )
    Widget w;			/* we don't need this anymore! */
    String func, msgName;
{
    char	func_name[BUFSIZ];
    String	ignored_type;
    XrmValue	result;

    if ( (Widget)0 == w )
	w = WcRootWidget(0);

    if ( BUFSIZ < WcStrLen(func) + WcStrLen(msgName) + 2 )
    {
	WcPrint( "Wcl Warning: WcErrorDatabaseText: ",
		 "strlen(func+msgName) > BUFSIZ ",
		 "for func = (", func, ") and msgName = (", msgName, ")",
		 NULL );
	return (String)0;
    }

    sprintf( func_name, "%s.%s", func, msgName );

    if ( (XrmDatabase)0 != wcxErrorDB )
    {
	XrmGetResource( wcxErrorDB,
		func_name, func_name,
		&ignored_type, &result );
	return (String)result.addr;
    }
    else
	return (String)0;
}

/* Default message print procedure
**=================================**
   Invoked via pointer to procedure WcPrint.
*/
#if NeedFunctionPrototypes
static void WcxPrint( String msg, ... )
#else
static void WcxPrint( msg, va_alist )
    String	msg;
    va_dcl
#endif
{
    if ( WcNonNull( msg ) )
    {
	va_list argPtr;

	(void)fputs( msg, stderr );

	Va_start( argPtr, msg );		/* point at 1st unnamed arg */

	msg = va_arg( argPtr, String );		/* get unnamed string arg */ 

	while ( msg != NULL )
	{
	    (void)fputs( msg, stderr );
	    msg = va_arg( argPtr, String );	/* get next string arg	*/
	}

	va_end( argPtr );			/* clean up */
    }

    (void)fputs( "\n", stderr );
}


/* Use dbFmt if non-null, make sure arg count matches %s in format string
**========================================================================**
*/
static void WcxWarn _(( int, String, String, ... ));

#if NeedFunctionPrototypes
static void WcxWarn( int numStringArgs, String defFmt, String dbFmt, ... )
#else
static void WcxWarn( numStringArgs, defFmt, dbFmt, va_alist )
    int		numStringArgs;
    String	defFmt, dbFmt;
    va_dcl
#endif
{
    char	buf[10*BUFSIZ];		/* a HUGE buffer	*/
    int		numPercentS;
    va_list	argPtr;

    Va_start( argPtr, dbFmt );		/* argPtr points at 1st unnamed arg */

    if ( WcNonNull( dbFmt ) )
    {
	numPercentS = WcPrintfFormatStrings( dbFmt );

	if ( numStringArgs != numPercentS )
	{
	    WcPrint( "Wcl Warning: Wcl Error DB Problem: ",
		     "Wrong number of `%s' in message from ErrorDB ( expected ",
		     WcxNumToStr( numStringArgs ),
		     ", found ",
		     WcxNumToStr( numPercentS ),
		     ") - message from Error DB is: ",
		     dbFmt,
		     NULL );
	}
	else if ( numStringArgs == 0 )
	{
	    WcPrint( dbFmt, NULL );
	}
	else
	{
	    vsprintf( buf, dbFmt, argPtr );
	    WcPrint( buf, NULL );
	}
    }
    else if ( WcNonNull( defFmt ) )
    {
	numPercentS = WcPrintfFormatStrings( defFmt );

	if ( numStringArgs != numPercentS )
	{
	    WcPrint( "Wcl Warning: Wcl Programming Error: ",
		     "Wrong number of `%s' in Default Error Msg ( expected ",
		     WcxNumToStr( numStringArgs ),
		     ", found ",
		     WcxNumToStr( numPercentS ),
		     ") - Default Error Msg is: ",
		     defFmt,
		     NULL );
	}
	else if ( numStringArgs == 0 )
	{
	    WcPrint( defFmt, NULL );
	}
	else
	{
	    vsprintf( buf, defFmt, argPtr );
	    WcPrint( buf, NULL );
	}
    }

    va_end( argPtr );
}


typedef struct _WcBuffer {
    int allocated, inUse;
    char* buffer;
} WcBufferStruct;

static char WcBuffer_couldNotAllocate[] =
"Wcl Warning: Could not allocate %d bytes for WcBuffer object";
static char WcBuffer_couldNotAllocateBuffer[] =
"Wcl Warning: Could not allocate %d bytes for msg buffer in WcBuffer object";
static char WcBuffer_couldNotGrowBuffer[] =
"Wcl Warning: Could not allocate %d bytes to grow buffer in WcBuffer object";

WcBuffer WcBuffer_New()
{
    WcBuffer this = (WcBuffer)XtCalloc( 1, sizeof(WcBufferStruct) );

    if ( this == NULL )
    {
	char buf[sizeof(WcBuffer_couldNotAllocate)+10];
	sprintf( buf, WcBuffer_couldNotAllocate, sizeof(WcBufferStruct) );
	WcWARN( (Widget)0, "WcBuffer_New", "couldNotAllocate", buf );
    }
    else
    {
	/* Note: MUST Calloc this, so we can always use StrCat instead
	 * of StrCpy the first (and only the first) time.
	 */
	this->buffer = XtCalloc( BUFSIZ, sizeof(char) );

	if ( this->buffer == NULL )
	{
	    char buf[BUFSIZ * sizeof(char) + 10];
	    sprintf( buf,WcBuffer_couldNotAllocateBuffer,BUFSIZ * sizeof(char));
	    WcWARN( (Widget)0, "WcBuffer_New", "couldNotAllocateBuffer", buf );
	}
	else
	{
	    this->allocated = BUFSIZ;
	    this->inUse     = 0;
	}
    }

    return this;
}

void WcBuffer_Free( this )
    WcBuffer this;
{
    if ( this != NULL )
    {
	if ( this->buffer != NULL )
	{
	    XtFree( this->buffer );
	}
	XtFree( (char*)this );
    }
}

#if NeedFunctionPrototypes
void WcBuffer_Append( WcBuffer this, String msg, ... )
#else
void WcBuffer_Append( this, msg, va_alist )
    WcBuffer this;
    String   msg;
    va_dcl
#endif
{
    va_list argPtr;

    if ( this == NULL || this->buffer == NULL || msg == NULL )
    {
	/* XtCalloc failed, already gave warning, or nothing to append.
	*/
	return;
    }

    Va_start( argPtr, msg );		/* point at 1st unnamed arg */

    while ( msg != NULL )
    {
	int msgLen;

	msgLen = strlen( msg );

	if ( msgLen + this->inUse >= this->allocated )
	{
	    /* Need to grow buffer.  At least double the size of the buffer.
	     * If the new msg component plus the existing buffered message
	     * is larger than twice the existing buffer then make it that much
	     * bigger.
	     */
	    if ( this->allocated * 2 > msgLen + this->inUse )
		this->allocated *= 2;
	    else
		this->allocated = msgLen + this->inUse + 1;

	    this->buffer = XtRealloc( this->buffer, this->allocated );

	    if ( this->buffer == NULL )
	    {
		char buf[sizeof(WcBuffer_couldNotGrowBuffer)+10];
		sprintf( buf, WcBuffer_couldNotGrowBuffer, this->allocated );
		WcWARN((Widget)0, "WcBuffer_Append", "couldNotGrowBuffer", buf);
		return;
	    }
	}

	/* Append the msg to the buffer.
	*/
	WcStrCat( this->buffer, msg );
	this->inUse += msgLen;

	msg = va_arg( argPtr, String );		/* get next string arg */ 
    }

    va_end( argPtr );			/* clean up */
}

String WcBuffer_String( this )
    WcBuffer this;
{
    if ( this == NULL || this->buffer == NULL )
	return "";
    else
	return this->buffer;
}




void WcWARN( w, func, msgName, defMsg )
    Widget w;
    String func, msgName, defMsg;
{
    WcxWarn( 0, defMsg, WcErrorDatabaseText( w, func, msgName ) );
}

void WcWARN1( w, func, msgName, defMsg, a1 )
    Widget w;
    String func, msgName, defMsg, a1;
{
    WcxWarn( 1, defMsg, WcErrorDatabaseText( w, func, msgName ), a1 );
}

void WcWARN2( w, func, msgName, defMsg, a1, a2 )
    Widget w;
    String func, msgName, defMsg, a1, a2;
{
    WcxWarn( 2, defMsg, WcErrorDatabaseText( w, func, msgName ), a1, a2 );
}

void WcWARN3( w, func, msgName, defMsg, a1, a2, a3 )
    Widget w;
    String func, msgName, defMsg, a1, a2, a3;
{
    WcxWarn( 3, defMsg, WcErrorDatabaseText( w, func, msgName ), a1, a2, a3 );
}
