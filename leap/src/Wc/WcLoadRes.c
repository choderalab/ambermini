#include "COPY.h"

/*
* SCCS_data: %Z% %M% %I% %E% %U%
*
* Widget Creation Library - WcLoadRes.c
*
* WcLoadResourceFile remembers which resource files have been loaded (by
* resource file name, not complete path name).  It handles absolute pathnames
* starting at the current working dorectory, root or tilda, or uses the same
* search path algorithm used by Xt: it tries to find the file in the system
* apps-default directory/directories, and then in the users app-defaults
* directory/directories.   It uses XFILESEARCHPATH, XUSERFILESEARCHPATH,
* XAPPLRESDIR, and HOME in the same way R3, R4, or R5 does (as controlled
* by XtSpecificationRelease from Intrinsic.h).
*
* This code tries to mimic GetAppUserDefaults() in mit/lib/Xt/Initialize.c
*
*******************************************************************************
*/

#include <X11/IntrinsicP.h>

#ifdef sun
#include <X11/ObjectP.h>	/* why don't they just use X from mit!?! */
#include <X11/RectObjP.h>
#endif

#include "WcCreateP.h"
#include <stdlib.h>		/* getenv, getuid	*/
#include <unistd.h>
#ifndef VMS
#include <pwd.h>		/* getpwuid, getpwnam	*/
#endif


/*
*******************************************************************************
* Private_function_declarations.
*******************************************************************************
*/

static int WcxLoadFileFromCurrentDirectory _(( Widget, char*, XrmDatabase* ));
static int WcxLoadFileUsingTildaPathname _(( Widget, char*, XrmDatabase* ));
static int WcxLoadSystemFile _(( Widget, char*, XrmDatabase* ));
static int WcxLoadUserFile _(( Widget, char*, XrmDatabase* ));

static void WcxLoadTrace _(( char* filename ));

static void WcxLoadWarnNoFileName _(( Widget ));
static void WcxLoadWarnNameTooLong _(( Widget, char* name ));
static void WcxLoadWarnSlashTerminated _(( Widget, char* name ));
static void WcxLoadWarnUnreadable _(( Widget, char* name, char* filename ));
static void WcxLoadWarnTooLong _(( Widget, char* name, char* path, char* file));
static void WcxLoadWarnMalloc _(( Widget, char* name ));

/* Find Home Directory
*******************************************************************************
    Used by WcLoadResourceFile and WcUserDefined.
    Argument can be THIS_USER (a NULL) or can be a user's name.
    The character string returned is in static storage, as returned
    by getpwnam().
*/

#define THIS_USER ((char*)0)

char* WcHomeDirectory( user )
    char* user;
{
    struct passwd *pw;
    char* homeDir = "";

    if ( WcNull(user) )
    {
	homeDir = getenv("HOME");

#ifndef VMS
        if (WcNull(homeDir))
        {
            char* thisUser = getenv("USER");

            if (WcNull(thisUser))
                pw = getpwuid( getuid() );
            else
                pw = getpwnam( thisUser );
            if (pw)
                homeDir = pw->pw_dir;
        }
    }
    else
    {
	/* some other user */
	pw = getpwnam( user );
	if (pw)
	    homeDir = pw->pw_dir;
#endif /*!VMS*/
    }
    return homeDir;
}

/*  Load Resource File - return TRUE if file was actually loaded
*******************************************************************************
    This function loads specified resource file into application resource 
    database.  The directory search for the file is (or tries hard to be)
    the same as for application class resource file.
*/

int WcLoadResourceFile ( w,  name )
    Widget w;		/* Wcl 1.x had this declared as a Display*	*/
    char*  name;	/* X resources file name			*/
{
    XrmQuark	filenameQuark;
    XrmDatabase	dbUsedByXt;

#if defined(XtSpecificationRelease) && XtSpecificationRelease > 4
    Screen*	correlation;
#else
    Display*	correlation;
#endif

    /*
     * Hack to prevent core-dumps if someone tries to pass us a display
     * as the first parameter (i.e., they think they are using Wcl 1.x)
     * When they re-compile, instead of just re-link, they will figure it
     * out (or on a Sun, when they run lint or Saber - whoops - Centerline),
     * so no need to give warning message.  The problem is that per-screen
     * resources are not *exactly* correct in the rare case of someone
     * having 2 kinds of screens on one display.
     */
    if ( w == (Widget)0 || w != w->core.self )
    {
	w = WcRootWidget(0);
    }

#if defined(XtSpecificationRelease) && XtSpecificationRelease > 4
    correlation = XtScreenOfObject( w );
    dbUsedByXt  = XtScreenDatabase( correlation );
#else
    correlation = XtDisplay( w );
    dbUsedByXt  = correlation->db;
#endif

    /* Check for repeated load, remember the first time.
    */
    filenameQuark = XrmStringToQuark( name );
    if ( MapAg_Find( rfAgent, correlation, (intptr_t)filenameQuark, NULL ) )
	return 0;	/* already loaded this file - OK, not an error */

    /* Remember we have loaded (or tried to load) this file.
    */
    MapAg_Define( rfAgent, correlation, (intptr_t)filenameQuark, NULL, 1 );

    return WcLoadResourceFileIntoDatabase( w, name, &dbUsedByXt );
}


/*  Load Resource File into Database - return TRUE if file was actually loaded
*******************************************************************************
    This function loads specified resource file into the resource database
    pointed to by an argument.  If the database is NULL, a new database is
    created and returned.  If the argument is non-NULL, then the file
    is merged with the existing database.

    The directory search for the file is (or tries hard to be) the same as
    for application class resource file.

    The return value is true if a resource file was found and loaded.
*/

int WcLoadResourceFileIntoDatabase( w, fileName, dbPointer )
    Widget       w;		/* for warning messages only		      */
    char*        fileName;	/* file name to load	     		      */
    XrmDatabase* dbPointer;	/* existing db or NULL: see XrmMergeDatabases */
{
    int loaded;

    /* Check for garbled pathname (missing or too long)
    */
    if ( WcNull( fileName ) )
    {
	WcxLoadWarnNoFileName( w );
        return 0;
    }
    if ( WcStrLen( fileName ) >= (int)MAX_PATHNAME )
    {
	WcxLoadWarnNameTooLong( w, fileName );
	return 0;
    }
    if ( fileName[ WcStrLen( fileName ) - 1 ] == '/' )
    {
	WcxLoadWarnSlashTerminated( w, fileName );
	return 0;
    }

    /* See if filename is a pathname from tilda
    */
    if ( '~' == fileName[0] )
    {
	return WcxLoadFileUsingTildaPathname( w, fileName, dbPointer );
    }

    /* Try to load file starting from current working directory.  This
     * works for "Hello", "./Hello", "../Mri/Hello", "subdir/Hello",
     * "/home/david/Hello", "../../MX400/app-defaults/Hello" etc ...
     */
    loaded = WcxLoadFileFromCurrentDirectory( w, fileName, dbPointer );

    if ( loaded )
    {
	/* File was found, so do NOT try to find using Xt-like file search.
	 * This means resource files in the current working directory 
	 * override any system or user resource files in other places.
	 */
	return 1;
    }
#ifndef VMS
    else if ( fileName[0] == '/' || fileName[0] == '.' )
#else
    else if ( fileName[0] == '[' )
#endif
    {
	/* Specific absolute paths must work, or give a warning message.
	*/
	WcxLoadWarnUnreadable( w, fileName, fileName );
	return 0;
    }

    /* Try loading like Xt: first system app-defaults, then user files.
    */
    loaded  = WcxLoadSystemFile( w, fileName, dbPointer );
    loaded += WcxLoadUserFile(   w, fileName, dbPointer );

    if ( loaded >= 1 )
	return 1;
    else
    {
	WcxLoadWarnUnreadable( w, fileName, fileName );
	return 0;
    }
}

/*ARGSUSED*/
static int WcxLoadFileFromCurrentDirectory( w, name, dbPointer )
    Widget	 w;		/* for calling consistency only */
    char*        name;
    XrmDatabase* dbPointer;
{
    XrmDatabase rdb = XrmGetFileDatabase( name );

    if ( (XrmDatabase)0 != rdb )
    {
	XrmMergeDatabases( rdb, dbPointer );
	WcxLoadTrace( name );
	return 1;
    }
    return 0;
}


static int WcxLoadFileUsingTildaPathname( w, name, dbPointer )
    Widget	 w;
    char*        name;
    XrmDatabase* dbPointer;
{
    XrmDatabase rdb;
    char*	homeDir;
    char	path[ MAX_PATHNAME+1 ];
    char	user[ MAX_PATHNAME+1 ];
    char*	from = &name[1];		/* skip the tilda */
    char*	to   = &user[0];

    /* Copy user's name into user[].  I.e., all characters between
    ** the `~' and the first `/`.  This may mean no characters (a NUL)
    ** which means we will need to look up the user's name.  We will
    ** not overflow the user[] array because we already checked the 
    ** length of name at the beginning.  `from' will point into `name'
    ** after the `~user', right at the `/`
    */
    while (*from && *from != '/')
	*to++ = *from++;
    *to = '\0';

    homeDir = WcHomeDirectory( user );

    if ( WcNull(homeDir) )
    {
	WcxLoadWarnUnreadable( w, name, name );
	return 0;
    }

    if( WcStrLen(homeDir) + WcStrLen(from) >= MAX_PATHNAME )
    {
	WcxLoadWarnTooLong( w, name, homeDir, from );
	return 0; 
    }

    WcStrCpy( path, homeDir );
    WcStrCat( path, from );

    if ( (XrmDatabase)0 != (rdb = XrmGetFileDatabase( path )) )
    {
	XrmMergeDatabases( rdb, dbPointer );
	WcxLoadTrace( path );
	return 1;
    }
    else
    {
	WcxLoadWarnUnreadable( w, name, path );
	return 0;
    }
}

static int WcxLoadSystemFile( w, name, dbPointer )
    Widget	 w;
    char*        name;
    XrmDatabase* dbPointer;
{
#ifdef XtSpecificationRelease
    /*
     * XtRelease 4 has XtResolvePathname() to do all the hunting around
     * in different directories using XFILESEARCHPATH or the compiled in
     * constant path XAPPLOADDIR.
     */
    char* filename = XtResolvePathname( XtDisplay(w), "app-defaults",
					name, NULL, NULL,
					NULL, 0, NULL );
    if ( WcNonNull( filename ) )
    {
	XrmDatabase rdb = XrmGetFileDatabase( filename );
	XrmMergeDatabases( rdb, dbPointer );
	WcxLoadTrace( filename );
	XtFree( filename );
	return 1;
    }
    return 0;

#else
    /*
     * For pre XtR4 (Motif 1.0), look in the site defined directory for
     * application defaults (commonly /usr/lib/X11/app-defaults).
     */
#ifndef XAPPDIR
#define XAPPDIR "/usr/lib/X11/app-defaults"
#endif

    XrmDatabase rdb;
    char        filename[ MAX_PATHNAME + 1 ];
    char*       sysResDir = XAPPDIR;

    if ( WcNull(sysResDir) )
	return 0;		/* already tried to load cwd files */

    if( ( WcStrLen(sysResDir) + 1 + WcStrLen(name) >= MAX_PATHNAME ) )
    {
	WcxLoadWarnTooLong( w, name, sysResDir, name );
	return 0;
    }

    WcStrCpy( filename, sysResDir );
    WcStrCat( filename, "/" );
    WcStrCat( filename, name );

    if ( (XrmDatabase)0 != (rdb = XrmGetFileDatabase(filename)) )
    {
	XrmMergeDatabases( rdb, dbPointer );
	WcxLoadTrace( filename );
	return 1;
    }
    return 0;

#endif
}

static int WcxLoadUserFile( w, name, dbPointer )
    Widget	 w;
    char*        name;
    XrmDatabase* dbPointer;
{
#ifdef XtSpecificationRelease
    /*
     * XtRelease 4 has XtResolvePathname() to do all the hunting around,
     * but we must build the search path like Xt does when looking for user
     * preferences: use XUSERFILESEARCHPATH and XAPPLRESDIR.
     */
    char* searchPath;
    int   deallocate;
    char* filename;
    char* homeDir             = WcHomeDirectory( (char*)0 );
#ifndef VMS
    char* xUserFileSearchPath = getenv("XUSERFILESEARCHPATH");
    char* xApplResDir         = getenv("XAPPLRESDIR");
#else
    char* xUserFileSearchPath = "DECW$USER_DEFAULTS%:%N";
    char* xApplResDir         = "";  /* What is it in VMS??? */
#endif

#ifdef VMS
    /*
     * This is an almost hopelessly incomplete attempt to make this
     * work with VMS, but it does at least make -resourceFile work...
     */
    if ( WcNull( getenv("DECW$USER_DEFAULTS") ) )
	xUserFileSearchPath = "DECW$SYSTEM_DEFAULTS%:%N";
#endif

    if ( WcNonNull( xUserFileSearchPath ) )
    {
	searchPath = xUserFileSearchPath;
	deallocate = 0;
    }
    else if ( WcNonNull( xApplResDir ) )
    {
	char* pattern =
#if XtSpecificationRelease > 4
"%s/%%L/%%N%%C:%s/%%l/%%N%%C:%s/%%N%%C:%s/%%N%%C:%s/%%L/%%N:%s/%%l/%%N:%s/%%N:%s/%%N";
#else
"%s/%%L/%%N%:%s/%%l/%%N%:%s/%%N%:%s/%%N%:%s/%%L/%%N:%s/%%l/%%N:%s/%%N:%s/%%N";
#endif
	searchPath = XtMalloc( 6*WcStrLen(xApplResDir) +
				     2*WcStrLen(homeDir) +
				     WcStrLen(pattern) );

	if ( (char*)0 == searchPath )
	{
	    WcxLoadWarnMalloc( w, name );
	    return 0;
	}
	deallocate = 1;

	sprintf( searchPath, pattern,
			xApplResDir, xApplResDir, xApplResDir,
			homeDir,
			xApplResDir, xApplResDir, xApplResDir,
			homeDir );
    }
    else
    {
	char* pattern =
#if XtSpecificationRelease > 4
"%s/%%L/%%N%%C:%s/%%l/%%N%%C:%s/%%N%%C:%s/%%L/%%N:%s/%%l/%%N:%s/%%N";
#else
"%s/%%L/%%N%:%s/%%l/%%N%:%s/%%N%:%s/%%L/%%N:%s/%%l/%%N:%s/%%N";
#endif
	searchPath = XtMalloc( 6*WcStrLen(homeDir) +
				     WcStrLen(pattern) );

	if ( (char*)0 == searchPath )
	{
	    WcxLoadWarnMalloc( w, name );
	    return 0;
	}
	deallocate = 1;

	sprintf( searchPath, pattern,
		 homeDir, homeDir, homeDir, homeDir, homeDir, homeDir );

    }

    filename = XtResolvePathname(XtDisplay(w), 0, name, 0, searchPath, 0, 0, 0);

    if ( WcNonNull(filename) )
    {
	XrmDatabase rdb = XrmGetFileDatabase( filename );
	XrmMergeDatabases( rdb, dbPointer );
	WcxLoadTrace( filename );
	XtFree( filename );
	if (deallocate)
	    XtFree( searchPath );
	return 1;
    }
    if (deallocate)
	XtFree( searchPath );
    return 0;

#else
    /*
     * For pre XtR4 (Motif 1.0), look for the file in XAPPLRESDIR and HOME
     */
    char* xApplResDir = getenv("XAPPLRESDIR");
    char* homeDir     = WcHomeDirectory( (char*)0 );
    char* filename;

    if ( WcNonNull( xApplResDir ) )
    {
	filename = XtMalloc( WcStrLen(xApplResDir)+WcStrLen(name)+1 );

	if ( (char*)0 == filename )
	{
	    WcxLoadWarnMalloc( w, name );
	    return 0;
	}
	WcStrCpy( filename, xApplResDir );
	WcStrCat( filename, name );

	if ( (XrmDatabase)0 != ( rdb = XrmGetFileDatabase( filename ) ) )
	{
	    XrmMergeDatabases( rdb, dbPointer );
	    WcxLoadTrace( filename );
	    XtFree( filename );
	    return 1;
	}
	XtFree( filename );
    }

    filename = XtMalloc( WcStrLen(homeDir)+WcStrLen(name)+2 );

    if ( (char*)0 == filename )
    {
	WcxLoadWarnMalloc( w, name );
	return 0;
    }
    WcStrCpy( filename, homeDir );
    WcStrCat( filename, "/" );
    WcStrCat( filename, name );

    if ( (XrmDatabase)0 != ( rdb = XrmGetFileDatabase( filename ) ) )

    {
	XrmMergeDatabases( rdb, dbPointer );
	WcxLoadTrace( filename );
	XtFree( filename );
	return 1;
    }
    XtFree( filename );
    return 0;

#endif
}

/*  -- Load resource files specified in Wcl record
*******************************************************************************
    Called by WcInitialize(), this loads all the resource files specified
    by the wclResFiles resource returning True if any file was loaded.
*/
int WcMoreResourceFilesToLoad( root, wcl )
    Widget	root;
    WclRecPtr	wcl;
{
    char* wcl_resFiles_copy;
    char* next;
    char  cleanName[MAX_XRMSTRING];
    int   loadedSomething = 0;

    /* Make copy of wcl->resFiles, as WcLoadResourceFile() may well 
    ** overwrite (and XtFree()) wcl->resFiles.
    */
    wcl_resFiles_copy = XtNewString( wcl->resFiles );
    next = WcCleanName( wcl_resFiles_copy, cleanName );

    while ( cleanName[0] != '\0' )
    {
	/* We may not need to load any given file
	*/
	loadedSomething += WcLoadResourceFile( root, cleanName );

	next = WcSkipWhitespace_Comma( next );
	next = WcCleanName( next, cleanName );
    }

    XtFree( wcl_resFiles_copy );
    return loadedSomething;
}

static void WcxLoadTrace( filename )
    String filename;
{
    if ( wcl->traceResFiles )
    {
	WcPrint( "Wcl loaded: ", filename, NULL );
    }
}

static void WcxLoadWarnNoFileName( w )
    Widget w;
{
WcWARN( w, "WcLoadResourceFile", "noFileName",
"Wcl Warning: WcLoadResourceFile() - No file name provided." );
}

static void WcxLoadWarnNameTooLong( w, fileName )
    Widget	w;
    String	fileName;
{
WcWARN1( w, "WcLoadResourceFile", "tooLong",
"Wcl Warning: WcLoadResourceFile( %s ) - File name too long.",
fileName );
}

static void WcxLoadWarnSlashTerminated( w, fileName )
    Widget	w;
    String	fileName;
{
WcWARN1( w, "WcLoadResourceFile", "slashTerminated",
"Wcl Warning: WcLoadResourceFile( %s ) - File must not end with slash ('/')",
fileName );
}

static void WcxLoadWarnUnreadable( w, name, filename )
    Widget	w;
    String	name, filename;
{
WcWARN2( w, "WcLoadResourceFile", "notReadable",
"Wcl Warning: WcLoadResourceFile( %s )\n\t- File <%s> not found, not readable, or wrong type.",
name, filename );
}

static void WcxLoadWarnTooLong( w, name, path, file )
    Widget	w;
    String	name, path, file;
{
WcWARN3( w, "WcLoadResourceFile", "tooLongExpanded",
"Wcl Warning: WcLoadResourceFile( %s )\n\t- Expanded name too long: %s%s", 
name, path, file );
}

static void WcxLoadWarnMalloc( w, name )
    Widget w;
    String name;
{
WcWARN1( w,  "WcLoadResourceFile", "mallocFailed",
"Wcl Warning: WcLoadResourceFile( %s ) - malloc failed",
name );
}
