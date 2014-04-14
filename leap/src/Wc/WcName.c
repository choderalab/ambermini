#include <stdlib.h>
#include "COPY.h"

/*
* SCCS_data: %Z% %M% %I% %E% %U%
*
* Widget Creation Library - WcName.c
*
* Implements several name-to-widget and widget-to-name functions as well as
* other general purpose parsing functions.
*
*******************************************************************************
*/

#include <X11/IntrinsicP.h>

#ifdef sun
#include <X11/ObjectP.h>	/* why don't they just use X from mit!?! */
#include <X11/RectObjP.h>
#endif

#include <X11/StringDefs.h>

#include "WcCreateP.h"


/* Private Data involving the root widget list
*******************************************************************************
*/

static int    numRoots = 0;
static Widget rootWidgets[MAX_ROOT_WIDGETS];


/*  -- Skip leading whitespace
*******************************************************************************
*/

char* WcSkipWhitespace( cp )
    char* cp;
{
    while ( *cp && isspace(*cp) )
        cp++;
    return cp;
}

char* WcSkipWhitespace_Comma( cp )
    char* cp;
{
    if (*cp == NUL) return cp;

    while ( *cp && isspace(*cp) )	cp++;
    if    ( *cp == ',' )		cp++;
    while ( *cp && isspace(*cp) )	cp++;
    return cp;
}

/*  -- Extract "clean" name from input string
*******************************************************************************
    Copies clean name into output buffer which caller provides.  Returns
    pointer to whitespace or NUL following clean name.
*/

char* WcCleanName( in, out )
    char* in;
    char* out;
{
    if (*in == NUL) { *out = NUL ; return in; }

    while ( *in && isspace(*in) )	/* in = WcSkipWhitespace( in ); */
	in++;
    while ( *in && !isspace(*in) && *in != ',' && *in != '(' && *in != ')' )
        *out++ = *in++;
    *out = NUL;
    return in;  /* this points at 1st whitespace or comma following "out" */
}

/*
******************************************************************************
**  -- Full Name To Widget and Supporting Static Functions
******************************************************************************
*/

/*  -- Map Name To Widget starting from any root widget.
******************************************************************************
    Goes down list of all RootWidgets maintained by WcRootWidget, and
    tries each one in turn as the root used to resolve the name into
    a widget.  Therefore, the name must be a path name starting at
    a root: i.e., it can be `*foo' or `root*foo' or similar.
*/

static Widget WcxPathFromAnyRootToWidget( widgetName )
    char* widgetName;
{
    int	widgetNameLen = WcStrLen(widgetName);
    int	i;

    for ( i = 0 ; i < numRoots ; i++ )
    {
	Widget	root		   = rootWidgets[i];
	char*	rootName	   = XtName( root );
	int	rootNameLen	   = WcStrLen( rootName );
	int	startsWithRootName = WcStrEqN(widgetName,rootName,rootNameLen);

	if ( startsWithRootName && widgetName[rootNameLen] == '*' )
	{
	    /* the root widget name is followed by a `*' so strip the
	    ** root name, but keep the star as it implies loose binding.
	    */
	    char* stripped = &widgetName[rootNameLen];
	    return WcChildNameToWidget( root, stripped );
	}

	else if ( startsWithRootName && widgetName[rootNameLen] == '.' )
	{
	    /* the root widget name is followed by a `.' so strip the
	    ** root name and the period to imply tight binding.
	    */
	    char* stripped = &widgetName[++rootNameLen];
	    return WcChildNameToWidget( root, stripped );
	}

	else if ( startsWithRootName && (widgetNameLen == rootNameLen) )
	{
	    /* widgetName is the root widget. */
	    return root;
	}

	/* Did not find the root name.  Try the next.
	*/
    }
    /* Did not find the named widget under any of the known roots.
    */
    return (Widget)NULL;
}

/*  -- Find the named widget
*******************************************************************************

    WcFullNameToWidget examines the first few characters of the `name' argument
    in order to determine which widget is the reference widget where the search
    will begin.

    The possibilities are these:  The `name' can begin with one or more of the
    relative prefix characters: ^ or ~ or /, which means the reference widget 
    will be some relative node in the widget tree, such as the parent or the 
    shell ancestor or the WMShell (dialogShell, topLevelShell, applicationShell)
    ancestor.  Otherwise, the  root widget will be the starting point.

    Once the starting point is found, XtNameToWidget is invoked.
*/

Widget WcFullNameToWidget( refWidget, name )
    Widget refWidget;
    char*  name;
{
    Widget    retWidget;
/*
    static int       check = 0;
*/
    if ( WcNull( name ))
	return NULL;
/*
#define TESTING(string, num)  if (WcStrEq( string, name )) check = num

    TESTING(("this")       ,1);
    TESTING(("aFind")      ,2);
    TESTING(("*aFindDialog"),3);
    TESTING(("xaLeap")     ,4);
    TESTING(("aAgain")     ,5);
*/
    
    if (WcStrEq( "this", name ))
	return refWidget;

    /* OK, now it gets more compilicated.
    ** See if we have a relative name.
    */

    if (   name[0] == '^'
	|| name[0] == '~'
	|| name[0] == '/'
	|| name[0] == '.'
	|| (   WcStrEqN( "this", name, 4 )
	    && (   name[4] == '.'
		|| name[4] == '*'
		|| name[4] == '^'
		|| name[4] == '~'
		|| name[4] == '/'	) ) )
    {
	/* Relative naming: widget must be in the same widget tree, and
	** must be the relative-named widget or a child thereof.
	*/
	int i = 0;

	if (name[0] == 't')	/* skip the "this" - not needed */
	    i = 4;

	while (name[i] == '^' 	/* parent */
	    || name[i] == '~' 	/* shell ancestor */
	    || name[i] == '/' 	/* WM shell ancestor */
	    || name[i] == '.')	/* eaten and ignored */
	{
	    /* Careful! don't go above a root widget! 
	    */
	    if (!XtParent(refWidget)
	     && (name[i] == '^' || name[i] == '~'|| name[i] == '/'))
	    {
		return NULL;
	    }
	    else if (name[i] == '^')
	    {
		refWidget = XtParent( refWidget );
	    }
	    else if (name[i] == '~')
	    {
		/* Some IntrinsicP.h forget to parenthesize the argument in
		** the definition of XtIsShell(), so the seemingly redundant
		** parenthesis are NECESSARY on some machines.
		**
		** Go up widget tree until we find a Shell:
		*/
		while ( XtParent(refWidget) &&
			(! XtIsShell( (refWidget = XtParent(refWidget)) ) ))
		    /*EMPTY*/;
	    }
	    else if (name[i] == '/')
	    {
		/* Find shell the window manager plays with (not menu)
		*/
		while ( XtParent(refWidget) &&
			(! XtIsWMShell( (refWidget = XtParent(refWidget)) ) ))
		    /*EMPTY*/;
	    }
	    i++;	/* eat this character */
	}
	/* We have eaten "this" and all of the relative naming characters
	** leaving nothing, or a name, or a name starting with '*'
	*/
	if ( name[i] == NUL )
	{
	    return refWidget;
	}
	return WcChildNameToWidget( refWidget, (char*)&(name[i]) );
    }
    /* Name did not begin with a relative naming operator.
    ** See if named widget is descended from reference widget.
    */
    /*ASSIGN_IN_IF*/
    if ((retWidget = WcChildNameToWidget( refWidget, name )))
    {
	return retWidget;
    }

    /* The Xaw CvtStringToWidget converter walks up the widget tree,
    ** parent by parent, looking for a child to match.  Perhaps this
    ** is a good idea.  It is certainly slow.  By default, we don't do
    ** this.  A Wcl library-wide resource can turn this on.
    */
    if ( wcl->slowTightNames )
    {
	Widget parent = XtParent( refWidget );
/*	while ( refWidget != (Widget)0 ) */
	while ( parent != (Widget)NULL )           /* Vladimir Romanovski */
	{
	    /*ASSIGN_IN_IF*/
	    if ( (retWidget = WcChildNameToWidget( parent, name )) )
		return retWidget;
	    parent = XtParent( parent );
	}
    }

    /* Start at top of this tree, and try to find named widget
    */
    /*ASSIGN_IN_IF*/
    if ((retWidget = WcChildNameToWidget( WcRootWidget(refWidget), name )))
    {
	return retWidget;
    }

    /* name not in widget tree under this root, or it begins with a root name.
    */
    /*ASSIGN_IN_IF*/
    if ((retWidget = WcxPathFromAnyRootToWidget( name )))
    {
	return retWidget;
    }

    /* Completely unsucessful in parsing this name.
    */
    return NULL;
}

/*  -- WidgetToFullName
*******************************************************************************
    Traverse up the widget tree, sprintf each name right up to
    the root of the widget tree.  sprintf the names to buffer.  Use
    recursion so order of names comes out right.  Client MUST free
    the char string alloc's and returned by WcWidgetToFullName().
*/

static char* nextChar;

static int WcxFullNameLen( w )
    Widget w;
{
    int len = 1 + WcStrLen( XtName(w) );

    if (w->core.parent)
	len += WcxFullNameLen(w->core.parent);
    return len;
}

static void WcxWidgetToFullName( w )
    Widget w;
{
    char* cp;

    if (w->core.parent)
    {
        WcxWidgetToFullName( w->core.parent );	/* nextChar AFTER parent name */
	*nextChar++ = '.';			/* inter-name `dot' */
    }

    cp = XtName(w);

    while (*cp)
	*nextChar++ = *cp++;
}

char* WcWidgetToFullName( w )
    Widget w;
{
    char* buff = XtMalloc( WcxFullNameLen( w ) );

    nextChar = buff;

    WcxWidgetToFullName( w );
    *nextChar = NUL;

    return buff;
}

static char* WcxLastSlash( cp )
    char* cp;
{
    char* slash = (char*)0;

    while ( *cp )
    {
	if ( *cp == '/' )
	    slash = cp;

	++cp;
    }

    return slash;
}

#ifdef VMS
static char* WcxLastRbrak( cp )
    char* cp;
{
    char* rbrak = (char*)0;

    while ( *cp )
    {
	if ( *cp == ']' )
	    rbrak = cp;

	++cp;
    }

    return rbrak;
}

static char* WcxFirstDot( cp )
    char* cp;
{
    while ( *cp )
    {
	if ( *cp == '.' )
	    return cp;

	++cp;
    }

    return (char*)0;
}
#endif

/*  -- Create an Application Class Name from Application Name
*******************************************************************************
    Initialize first letter to make class, or first two if first is already
    capitalized, or don't worry about it.
*/

char* WcAppNameToAppClass( appName )
    char* appName;
{
#ifndef VMS
    char* lastSlash = WcxLastSlash( appName );
    char* appClass  = XtNewString( lastSlash ? lastSlash+1 : appName );
#else
    /* VMS names look like this:
     * <VMS file name> ::= <bracketed option><base name><extension option>
     * <bracketed option> ::= <empty>
     *                    |   [ <stuff> ]
     *                    |   <bracketed option> [ <stuff> ]
     * <base name>        ::= <character string>
     * <extension option> ::= <empty>
     *                    |   .<stuff>
     *
     * We need the <base name> part.  Find the last ']' and then the first '.'
     * after any such optional ']'.
     */
    char* lastRbrak = WcxLastRbrak( appName );
    char* firstDot  = WcxFirstDot( lastRbrak ? lastLbrak+1 : appName );
    char* appClass;

    if ( firstDot )
	*firstDot = '\0';	/* truncate appName at this firstDot */

    /* <base name> out of <VMS file name>
    */
    appClass = XtNewString( lastRbrak  ? lastRbrak+1  : appName );
#endif

    if (islower(appClass[0]))
	appClass[0] = toupper(appClass[0]);
    else if (islower(appClass[1]))
	appClass[1] = toupper(appClass[1]);

    return appClass;
}

/*  -- Find Application Name and Class Name in argc, argv
*******************************************************************************
    Looks for "-name" followed by application class.
    If there is no "-name", it uses argv[0].
*/

char* WcAppName( argc, argv )
    int argc;
    char** argv;
{
    int i;
    for ( i = 0 ; i < argc ; i++ )
    {
	if ( WcStrEq( argv[i], "-name" ) )
	{
	    if ( i < argc-1 )
		return XtNewString( argv[i+1] );
	    else
		return XtNewString( argv[0] );
	}
    }
    return XtNewString( argv[0] );
}

char* WcAppClass( argc, argv )
    int argc;
    char** argv;
{
    String appName   = WcAppName( argc, argv );
    String className = WcAppNameToAppClass( appName );
    XtFree( appName );
    return className;
}

/*  --  Quark from Lower Case String
*******************************************************************************
    This allows strings to be quarkified in a case-insensitve manner,
    as used by the registration and conversion routines.
*/
#define isname(c) (isalnum(c) || c == '_')

XrmQuark WcStringToQuark( cp )
    char* cp;
{
    char  lower[NAME_RESOLUTION];
    int   i;

    while ( *cp && !isname(*cp) )
	++cp;
    for ( i = 0  ;  isname(*cp) && i < NAME_RESOLUTION-1  ;  cp++, i++ )
	lower[i] = ( isupper(*cp) ? tolower(*cp) : *cp );
    lower[i] = NUL;

    return XrmStringToQuark(lower);
}
#undef isname

/*   --- Quark from segment of string
*******************************************************************************
    Two versions to maintain /&'$%/'& backward compatibility.  Yuck!
    The first is the old a stupid version, where for some goddamn reason
    I flatten the case of the callback proc names, so the resource files can
    be more case insensitive.  

    For dynamic binding, however, we MUST have a case sensitive quark, so
    we can go from quark to actual string.  We don't need registration,
    so the string in the resource file is the only thing we can go on,
    and the dynamic linking stuff is (of course) case sensitive.  Man, I
    am so pissed I ever took the suggestion to flatten the case.  Its been
    a pain in the ass for years.
*/
XrmQuark WcSubStringToQuark( cp, end )
    char* cp;	/* first character of substring */
    char* end;	/* last character of substring */
{
    char  lower[NAME_RESOLUTION];
    int   i;

    if ( cp == (char*)0 || end == (char*)0 )
	return NULLQUARK;

    for ( i = 0 ; cp <= end && i < NAME_RESOLUTION-1 ; cp++, i++ )
	lower[i] = ( isupper(*cp) ? tolower(*cp) : *cp );
    lower[i] = NUL;

    return XrmStringToQuark(lower);
}

XrmQuark WcCsSubStringToQuark( cp, end )
    char* cp;	/* first character of substring */
    char* end;	/* last character of substring */
{
    char  name[NAME_RESOLUTION];
    int   i;

    if ( cp == (char*)0 || end == (char*)0 )
	return NULLQUARK;

    for ( i = 0 ; cp <= end && i < NAME_RESOLUTION-1 ; cp++, i++ )
	name[i] = *cp;
    name[i] = NUL;

    return XrmStringToQuark(name);
}

/*  --  Recursive PrintTree
*******************************************************************************
    Prints the names and classes of all children of a given widget.
    Each nested level is indented 2 additional spaces.
*/

static void WcxPrintTree _(( WcBuffer buf, Widget w, int depth ));

static void WcxPrintTree( buf, w, depth )
    WcBuffer buf;
    Widget   w;
    int      depth;
{
    int		i;
    char*	class = WcWidgetClassName( w );		/* DONT FREE */
    char*	name  = WcWidgetToFullName( w );	/* MUST FREE */
    char	spaces[81];

    *spaces = NUL;

    for (i = 0 ; i < depth && i < 40 ; i++)
	WcStrCat( spaces, "  " );
    depth++;

    WcBuffer_Append( buf, spaces, name, " of class ", class, "\n", NULL );

    XtFree( name );

    if ( !XtIsWidget( w ) )
	return; 		/* Gadgets cannot have children of any kind */

    if ( 0 < w->core.num_popups )
    {
	WcBuffer_Append( buf, spaces, "Popups:\n", NULL );

	for (i = 0 ; i < w->core.num_popups ; i++ )
	    WcxPrintTree( buf, w->core.popup_list[i], depth );
    }

    if ( XtIsComposite( w ) )
    {
	CompositeWidget cw = (CompositeWidget) w;
	if ( 0 == cw->composite.num_children )
	{
	    WcBuffer_Append( buf, spaces, "No Children.\n", NULL );
	}
	else
	{
	    WcBuffer_Append( buf, spaces, "Children:\n", NULL );

	    for (i = 0 ; i < cw->composite.num_children ; i++ )
	        WcxPrintTree( buf, cw->composite.children[i], depth );
	}
    }
}

void WcPrintTree( w )
    Widget w;
{
    WcBuffer buf = WcBuffer_New();

    WcxPrintTree( buf, w, 0 );

    WcPrint( WcBuffer_String( buf ), NULL );

    WcBuffer_Free( buf );
}

/*
*******************************************************************************
* Private Data involving the root widget list, declared at top of this file
*	static int numRoots = 0;
*	static Widget rootWidgets[MAX_ROOT_WIDGETS];
*******************************************************************************
*/


/*  -- Forget about a root widget
*******************************************************************************
    When a root widget gets destroyed, we need to take that widget out
    of our list of root widgets.  This is a destroy callback routine
    which is added to a root widget's destroy callback list by WcRootWidget.
*/
/*ARGSUSED*/
static void ForgetRootCB( w, ignored, unused )
    Widget	w;
    XtPointer	ignored;
    XtPointer	unused;
{
    int i;
    for (i = 0 ; i < numRoots ; i++ )
    {
        if ( w == rootWidgets[i] )
	{
	    /* move all following widgets up to close the gap */
	    for ( ; i < numRoots ; i++ )
	    {
		rootWidgets[i] = rootWidgets[i+1];
	    }
	    numRoots-- ;
	    return;
	}
    }
    /* should never get here */
}

/*  -- Find root widget
*******************************************************************************
    If a widget is passed, then find the root of that widget.  See if
    it is one of the root widgets we already know about.  Add to list
    if not.  Return the root widget.

    If no widget is passed, then return the first root widget we know
    about.  If we know of no root widgets, then we will return a NULL
    since the rootWidgets[] array starts out filled with nulls, and
    gets re-filled as roots are destroyed.
*/

Widget WcRootWidget( w )
    Widget w;
{
    int i;

    if (w)
    {
	while ( XtParent(w) )
	    w = XtParent(w);

	for (i = 0 ; i < numRoots ; i++)
	{
	    if ( w == rootWidgets[i] )
		return w;
	}

	rootWidgets[i] = w;
	numRoots++;
	XtAddCallback( w, XtNdestroyCallback, ForgetRootCB, NULL );
	return w;
    }
    else
    {
	return rootWidgets[0];
    }
}

/* -- Equivalent to ANSI C library function strstr()
*******************************************************************************
*/

char* WcStrStr( s1, s2 )
    char* s1;
    char* s2;
{
    while (*s1)
    {
	if (*s1 == *s2)
	{
	    char* start = s1;
	    char* c = s2;
	    while (*++s1 & *++c && *s1 == *c)
		/*EMPTY*/;
	    if (*c == '\0')
		return start;
	    else
		s1 = ++start;
	}
	else
	{
	    s1++ ;
	}
    }
    return (char*)0;
}

/* Safe forms of strlen, strcpy, strcat, strcmp, and strncmp (also streq)
** Some are macros in WcCreate.h
*/

char* WcStrCpy( str1, str2 )
    char* str1;
    char* str2;
{
    if ( str1 == NULL )
    {
	WcWARN( (Widget)0, "WcStrCpy", "nullTarget",
		"Wcl Warning: WcStrCpy() called with NULL target!");
	abort();
    }
    if ( str2 == NULL )
	*str1 = NUL;
    else
	str1 = strcpy( str1, str2 );

    return str1;
}

char* WcStrCat( str1, str2 )
    char* str1;
    char* str2;
{
    if ( str1 == NULL )
    {
	WcWARN( (Widget)0, "WcStrCat", "nullTarget",
		"Wcl Warning: WcStrCat() called with NULL target!");
	abort();
    }
    if ( str2 == NULL )
	return strcat( str1, "" );
    else
	return strcat( str1, str2 );
}

int WcStrCmp( str1, str2 )
    char* str1;
    char* str2;
{
    if ( WcNonNull(str1) && WcNonNull(str2) )
	return strcmp( str1, str2 );	/* str1 and str2 are non-null */
    if ( WcNull(str1) && WcNull(str2) )
	return 0;			/* both null, I guess that's equal */
    if ( WcNull(str1) )
	return -1;			/* a NULL is less than anything */
    /* if ( WcNull(str2) ) */
	return 1;			/* anything is greater than NULL */
}

int WcStrCmpN( str1, str2, num )
    char* str1;
    char* str2;
    int   num;
{
    if ( num > 0 && (WcNonNull(str1) && WcNonNull(str2)) )
	return strncmp( str1, str2, num);
    if ( num == 0 || (WcNull(str1) && WcNull(str2)) )
	return 0;			/* equal, I guess */
    if ( WcNull(str1) )
	return -1;			/* NULL is less than anything */
    /* if ( WcNull(str2) ) */
	return 1;			/* anything is greater than NULL */
}

/**************************************************************************
These functions make it easy to determine if a format string
can be used by printf, fprintf, or sprintf to print the appropriate
number of strings.  This allows error messages to be re-configured 
but remain safe at execution time.  This makes it much safer to change
error messages used by the WARN procedures.

It also allows textual values to be set on arbitrary widgets like this:

void RegisterNameDisplayCB( w, sprintFmt, ignored )
    Widget w; XtPointer sprintFmt, ignored;
{
    format = (char*)sprintFmt;
    if ( 2 != WcPrintfFormatStrings( format ) )
    {
	WARN( w, "RegisterNameDisplay", "needTwoStringFormats",
		"Sorry, you blew it!" )
	format = "this.labelString: Send %s To %s";
    }
}

void DisplaySendMessage( w, note, name )
    Widget w; String note, name;
{
    char buf[MAX_XRMSTRING];
    sprintf( buf, format, note, name );
    WcSetValues( w, buf );
}

See? If the user interface instead uses Athena widgets, the format string
passed from the RegisterNameDisplay() in the resource file passes something
like "this.label: Send %s to %s" or if the UI really wants the name to show
up in the shell, pass "~.title: Send %s to %s".
**************************************************************************/

/* Number Of Strings Which May Be Printf'd
**=========================================**
   Returns 0 if something other than strings are in the format
   specification.  Returns 0 if a string format specification contains a
   '*' which means the width or precision is parametric.  Therefore, a
   positive return value means printf() can be safely called with the
   returned number of strings as arguments following the format.
*/
int WcPrintfFormatStrings( cp )
    char* cp;
{
    int strings = 0;
#ifdef SUNOS_PRINTF
    int digitDollarSpec, field[10], i;
    for ( i = 1 ; i <= 9 ; i++ )	/* NB: 1 <= field <= 9 */
	field[i] = 0;
#endif

    while ( cp && *cp )
    {
	if ( *cp != '%' )
	    ++cp;
	else if ( *++cp == '%' )	/* look at character after %	*/
	    ++cp;			/* also a % - not a conversion	*/
	else
	{
	    /* conversion specification: already ate the %
	    */

#ifdef SUNOS_PRINTF
	    /* WARNING! Most machines do NOT support this %digit$ concept.
	    */
	    if ( '1' <= *cp && *cp <= '9' && *(cp+1) == '$' )
	    {
		digitDollarSpec = 1;
		field[*cp - '0']++;
		cp = cp+2;
	    }
	    else
		digitDollarSpec = 0;
#endif

	    /* optional flags:
	    */
	    while (*cp == '+' || *cp == '-' || *cp == ' ' || *cp == '#' )
		cp++;

	    /* optional field width:
	    */
	    if ( *cp == '*' )
		return 0;		/* No parametric field widths! */
	    while (*cp && '0' <= *cp && *cp <= '9')
		++cp;

	    /* optional precision:
	    */
	    if (*cp == '.')
	    {
		++cp;
		if ( *cp == '*' )
		    return 0;		/* No parametric precision! */
		while (*cp && '0' <= *cp && *cp <= '9')
		    ++cp;
	    }

	    /* optional long conversion
	    */
	    if (*cp == 'l')
		++cp;

	    /* type of conversion: anything but string returns zero
	    */
	    if (*cp != 's')
		return 0;

#ifdef SUNOS_PRINTF
	    if ( !digitDollarSpec )
		strings++;
#else
	    strings++;
#endif
	}
    }
#ifdef SUNOS_PRINTF
    /* If format spec named specific fields with %digit$ specification,
    ** then the number of strings is at least that number of fields.
    */
    for ( i = 1 ; i <= 9 ; i++ )
	if ( field[i] && strings < i )
	    strings = i;
#endif 
    return strings;
}

/* Break String into Lines
**********************************************************************
*/

String WcBreakIntoLines( string, columnsPerLine )
    String	string;
    int		columnsPerLine;
{
    int		len, charPos;
    String	cp, retVal, lastWhitespace;

    /* Allocate enough space for worst case, where subject is continuous
     * text without any whitespace.
     */
    len = WcStrLen( string );

    cp = retVal = XtMalloc( len + len/columnsPerLine + 1 );
    if ( cp == NULL )
    {
	/* malloc failed */
	return NULL;
    }

    for ( lastWhitespace = (char*)0, charPos = 0  ;  *string  ;  )
    {
	if ( *string == '\n' )
	{
	    *cp++ = *string++;
	    lastWhitespace = (char*)0, charPos = 0;
	    continue;
	}
	if ( charPos < columnsPerLine )
	{
	    if ( isspace( *string ) )
		lastWhitespace = cp;

	    if ( ispunct( *string )
	      && ( string[1] && string[1] == ' ' )
	      && ( string[2] && string[2] == ' ' )
	      && ( string[3] && string[3] != ' ' ) )
	    {
		/* An American!
		*/
		*cp++ = *string++;		/* the punctuation mark */
		string++;			/* skip one blank	*/
	    }
	    else
		*cp++ = *string++;

	    ++charPos;
	}
	else
	{
	    /* At number of columns
	    */
	    if ( lastWhitespace == (char*)0 )
	    {
		/* No whitespace in line - just break line here.  We ADD
		 * a character to the output buffer.
		 */
		*cp++ = '\n';
		*cp++ = *string++;
		charPos = 1;
	    }
	    else if ( ispunct( *string ) )
	    {
		/* Can end line on punctuation mark.  Drop multiple spaces!
		 * This ADDs a character to output buff if next char != ' '
		 */
		*cp++ = *string++;
		*cp++ = '\n';
		charPos = 0;

		if ( *string == '\n' )
		    string++;
		else
		    while ( *string == ' ' )
			string++;
	    }
	    else
	    {
		/* Overwrite last whitespace with newline
		*/
		*lastWhitespace = '\n';
		*cp++ = *string++;
		charPos = cp - lastWhitespace;
	    }
	    lastWhitespace = (char*)0;
	}
    }
    *cp = '\0';

    return retVal;
}

String WcSpliceLines( cp, columnsPerLine )
    String	cp;
    int		columnsPerLine;
{
    String	rp, retVal;
    int		col;

    rp = retVal = (String)XtMalloc( WcStrLen( cp ) + 1 );

    col = 1;
    while ( *cp )
    {
	if ( *cp == '\n' && col <= columnsPerLine )
	{
	    /* Maybe we will remove or replace (with space) this newline.
	     * Look ahead beyond the newline for the first whitespace.
	     */
	    int ahead = 1;

	    while ( cp[ahead] && !isspace( cp[ahead] ) )
		++ahead;

	    if ( ahead == 1 )
	    {
		/* The next character after this newline is whitespace
		 * (or the end of the string), therefore this newline
		 * is intentional - keep it.
		 */
		*rp++ = *cp++;
	    }
	    else if ( col + 1 + ahead < columnsPerLine )
	    {
		/* The next word could easily fit on the line: therefore,
		 * assume this newline is intentional, and keep it.
		 */
		*rp++ = *cp++;
	    }
	    else if ( ahead < columnsPerLine && columnsPerLine < col+1+ahead )
	    {
		/* The next character after the newline is NOT whitespace.
		 * The next "word" plus a whitespace (instead of the newline)
		 * would have made a line over 74 characters long - hence,
		 * this line was (probably) broken to make the lines less
		 * than 74 characters.  Therefore, this is a newline we can
		 * replace with a blank.
		 */
		++cp;		/* skip this newline			*/
		*rp++ = ' ';	/* put a space in the return value	*/
	    }
	    else
	    {
		/* We looked columnsPerLine or more characters ahead before we
		 * found another whitespace (or the EOS).  Therefore, we can
		 * assume this is a long word broken in the middle.  We drop
		 * the newline, the return value is now SHORTER.
		 */
		++cp;		/* skip this newline 			*/
		*rp++ = *cp++;	/* copy NEXT char to return value	*/
	    }
	    /* Reset column count on every newline in source string.
	    */
	    col = 1;
	}
	else if ( *cp == '\n'   /* && col > columnsPerLine */ )
	{
	    /* This newline seems to be intentional (actually, this entire
	     * cover letter appears that it has NOT been broken into 74
	     * character lines - therefore, ANY newline must be intentional).
	     */
	    *rp++ = *cp++;
	    col = 1;
	}
	else
	{
	    /* Pass character through
	    */
	    *rp++ = *cp++;
	    ++col;
	}
    }
    *rp = *cp;		/* null terminator */

    return retVal;
}

static char* XEvent_names[] = {
"",
"",
"KeyPress",             /* 2 */
"KeyRelease",           /* 3 */
"ButtonPress",          /* 4 */
"ButtonRelease",        /* 5 */
"MotionNotify",         /* 6 */
"EnterNotify",          /* 7 */
"LeaveNotify",          /* 8 */
"FocusIn",              /* 9 */
"FocusOut",             /* 10 */
"KeymapNotify",         /* 11 */
"Expose",               /* 12 */
"GraphicsExpose",       /* 13 */
"NoExpose",             /* 14 */
"VisibilityNotify",     /* 15 */
"CreateNotify",         /* 16 */
"DestroyNotify",        /* 17 */
"UnmapNotify",          /* 18 */
"MapNotify",            /* 19 */
"MapRequest",           /* 20 */
"ReparentNotify",       /* 21 */
"ConfigureNotify",      /* 22 */
"ConfigureRequest",     /* 23 */
"GravityNotify",        /* 24 */
"ResizeRequest",        /* 25 */
"CirculateNotify",      /* 26 */
"CirculateRequest",     /* 27 */
"PropertyNotify",       /* 28 */
"SelectionClear",       /* 29 */
"SelectionRequest",     /* 30 */
"SelectionNotify",      /* 31 */
"ColormapNotify",       /* 32 */
"ClientMessage",        /* 33 */
"MappingNotify",        /* 34 */
"LASTEvent",            /* 35 */        /* must be bigger than any event # */
};
int XEvent_numNames = 35;

extern char* WcXEventName( event )
    XEvent* event;
{
    return XEvent_names[ event->type ];
}

extern char* WcWidgetClassName( widget )
    Widget widget;
{
    return widget->core.widget_class->core_class.class_name;
}

