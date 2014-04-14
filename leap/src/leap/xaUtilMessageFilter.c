/*
 *	File:	xaUtilMessFilter.c
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
 *		Read from the standard input.
 */




#include	"basics.h"

#include 	<X11/IntrinsicP.h>
#include	<X11/StringDefs.h>


#include	<X11/keysym.h>

#include 	"../Wc/WcCreate.h"

#include	"varArray.h"



#include	"function.h"
#include	"mess.h"


#include	"X11/Xaw/List.h"
#include	"X11/Xaw/Toggle.h"
#include	"X11/Xaw/AsciiText.h"
#include	"X11/Xaw/Label.h"
#include	<X11/Shell.h>




/*
 *-------------------------------------------------------------
 *
 *	Global and Static variables
 */


#define	TOGGLE_WIDGET_NAME	"function"
#define	FILE_WIDGET_NAME	"file"
#define	FUNCTION_POPUP_NAME	"funcShell"

#define	XUMF_SET		1
#define	XUMF_RESET		2
#define	XUMF_TOGGLE		3


static	int	SiCurrentFunction = 0;
static	int	SiCurrentFile = 0;



#define	BUFFER_INC	16*1024

Widget		GwText;
Widget		GwTextSource = NULL;
Widget		GwTopWidget;
VARARRAY	GvaInitialOnFunctions = NULL;
BOOL		GbInputFile = FALSE;
STRING		GsInitialOnFunctionFilename;
STRING		GsInputFilename;
XtAppContext	GxtacApp;
char*		GcPMessBuffer = NULL;
char*		GcPMessTail;
int		GiSize = 0;
int		GiInputDescriptor;
XtInputId	GxtiiId;
FILE*		GfInput;
BOOL		GbTextNeedsUpdate = FALSE;


/*
 *-------------------------------------------------------------
 *
 *	Private functions
 *
 */




/*
 *	XUMFUpdateTextWidget
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Tell the Text Widget that it has new text.
 */
static void
XUMFUpdateTextWidget()
{

    if ( GwTextSource == NULL ) {
	GwTextSource = XawTextGetSource(GwText);
    }
    if ( GcPMessBuffer ) {
	XtVaSetValues( GwTextSource,
			XtNstring, (XtArgVal) GcPMessBuffer,
			NULL );
   } else {
	XtVaSetValues( GwTextSource,
			XtNstring, (XtArgVal) "---nothing---",
			NULL );
   }

}


/*
 *	XUMFAppendToText
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Loop through all the messages, and for all those
 *	that have bPrint on, write them into the character
 *	buffer in (*cPPBuffer).
 */
static void
XUMFAppendToText( char *cPData )
{
char	*cPTemp;
int	iLen;

    iLen = strlen(cPData);
    if ( GcPMessBuffer == NULL ) {
	MALLOC( GcPMessBuffer, char*, BUFFER_INC+10 );
	GiSize = BUFFER_INC;
	GcPMessTail = GcPMessBuffer;
    } else {
	if ( GiSize < ( (GcPMessTail+iLen) - GcPMessBuffer ) ) {
	    REALLOC( cPTemp, char*, GcPMessBuffer, 
		    GiSize+BUFFER_INC+10 );
	    GcPMessTail = cPTemp + ( GcPMessTail - GcPMessBuffer );
	    GcPMessBuffer = cPTemp;
	    GiSize += BUFFER_INC;
	}
    }
    strcpy( GcPMessTail, cPData );
    GcPMessTail += iLen;

}





/*
 *	XUMFAppendMessToText
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Append the message to the Text buffer.
 */
static void
XUMFAppendMessToText( MESS mMess )
{
STRING		sLine;

    if ( iMessFunction(mMess) == NO_FUNCTION ) {
	sprintf( sLine, "  %s\n", sMessText(mMess) );
    } else {
	sprintf( sLine, "%s|%s\n",
		sFunctionFunction(iMessFunction(mMess)),
		sMessText(mMess) );
    }
    XUMFAppendToText(sLine);
}






/*
 *	XUMFInputAvailable
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback to handle data from input file.
 */
static void
XUMFInputAvailable( XtPointer xtPData, int *iPSource, XtInputId *xtiiPId )
{
STRING		sLine;
char		*cPRead;
MESS		mMess;
int		iLen;

    cPRead = fgets( sLine, sizeof(STRING), GfInput );
    if ( !cPRead ) {
	return;
    }
    iLen = strlen(sLine);
    if ( iLen > 0 ) {
	if ( sLine[iLen-1] == '\n' ) sLine[iLen-1] = '\0';
    }
    mMess = mMessAdd( sLine );
    if ( bMessPrint(mMess) ) {
	XUMFAppendMessToText(mMess);

		/* Do not update the display, the user */
		/* must select 'Rebuild Mess List' */
    }

}

    

	


/*
 *	XUMFReadLineFile
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Read a file that describes where all of the functions
 *	in a source file can be found.
 */
static void
XUMFReadLineFile( char *sFilename )
{
FILE		*fIn;
STRING		sLine, sFile, sFunction;
int		iStart, iStop, iRead;
BOOL		bGotFile;
char		*cPRead;

    fIn = fopen( sFilename, "r" );
    if ( fIn == NULL ) {
	DFATAL(( "Could not read file: %s\n", sFilename ));
    }
    bGotFile = FALSE;
    while ( !feof(fIn) ) {
	strcpy( sLine, "" );
	cPRead = fgets( sLine, sizeof(STRING), fIn );
	if ( cPRead == NULL ) break;
	sscanf( sLine, "%s %d %d %s", sFile, &iStart, &iStop, sFunction );
	FunctionAdd( sFile, iStart, iStop, sFunction );
	if ( !bGotFile ) {
	    bGotFile = TRUE;
	    FunctionFileAdd( sFile );
	}
    }
    fclose(fIn);
}



/*
 *	XUMFReadInitialOnFile
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Read a file that describes which functions are
 *	initially on.  Read it into the global VARARRAY
 *	(GvaInitialOnFunctions).
 */
static void
XUMFReadInitialOnFile( char *sFilename )
{
FILE		*fIn;
STRING		sLine;

    GvaInitialOnFunctions = vaVarArrayCreate(sizeof(STRING));

    fIn = fopen( sFilename, "r" );
    if ( fIn == NULL ) {
	DFATAL(( "Could not read file: %s\n", sFilename ));
    }
    while ( !feof(fIn) ) {
	strcpy( sLine, "" );
	fgets( sLine, sizeof(STRING), fIn );
	if ( sLine[0]=='\n' ) break;
	if ( strlen(sLine)>0 ) {
	    sLine[strlen(sLine)-1] = '\0';
	}
	VarArrayAdd( GvaInitialOnFunctions, (GENP)sLine );
    }
    fclose(fIn);
}




/*
 *	XUMFChangeMessPrint
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Change the print status for the messages
 *	that belong to (iFunction).
 */
static void
XUMFChangeMessPrint( int iFunction, BOOL bPrint )
{
MESS		mMess;


		/* Search through all the messages for */
		/* those from this function */
    mMess = mMessLoop();
    while ( mMessNext(&mMess) ) {
	if ( iMessFunction(mMess) == iFunction ) {
	    MessSetPrint( mMess, bPrint );
	}
    }

}


	

/*
 *	XUMFChangePrintForMesssFrom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Change the Print status for messages
 *	from the function in (sFunction).
 */
static void
XUMFChangePrintForMesssFrom( char *sFunction, BOOL bPrint )
{
int		iFunction;

    iFunction = iFunctionFindWithFunction(sFunction);
    if ( iFunction == NO_FUNCTION ) {
	DFATAL(( "Unknown function: %s\n", sFunction ));
    }
    XUMFChangeMessPrint( iFunction, bPrint );
}


	

/*
 *	XUMFChangeToggleState
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Change the state of the Toggle that is associated
 *	with a function, also change the state of all
 *	the messages associated with that function, and
 *	the state of the function.
 *
 *	(iAction) can be one of XUMF_SET, XUMF_RESET, XUMF_TOGGLE
 */
static void
XUMFChangeToggleState( Widget wToggle, int iAction )
{
char		*cPFunction;
int		iFunction;
BOOL		bOldState, bNewState;

    XtVaGetValues( wToggle,
			XtNlabel, &cPFunction,
			XtNstate, &bOldState,
			NULL );

    iFunction = iFunctionFindWithFunction( cPFunction );
    if ( iFunction == NO_FUNCTION ) {
	DFATAL(( "Unknown function name: %s\n", cPFunction ));
    }

    switch ( iAction ) {
	case XUMF_SET:
	    bNewState = TRUE;
	    break;
	case XUMF_RESET:
	    bNewState = FALSE;
	    break;
	case XUMF_TOGGLE:
	    bNewState = !bOldState;
	    break;
	default:
	    DFATAL(( "Illegal action \n" ));
	    break;
    }

			/* Change the Widget */

    XtVaSetValues( wToggle,
			XtNstate, (XtArgVal) bNewState,
			NULL );

			/* Change the FUNCTION */

    FunctionSetPrint( iFunction, bNewState );

			/* Change the messages */

    XUMFChangeMessPrint( iFunction, bNewState );

}


    
    




/*
 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
 *
 *	Callbacks
 *
 */




/*
 *	XUMFToggleNotify
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Save the Widget which maintains the command line interface.
 *	The argument to this function must be the
 *	name of the Action that is used by the TextWidget
 *	to insert strings into the Text Widget.
 */
XtCallbackProc
XUMFToggleNotify( Widget wNew, char *sActionName, caddr_t PData2 )
{
char		*cPFunction;
int		iFunction;
BOOL		bPrint;

    XtVaGetValues( wNew,
			XtNlabel, &cPFunction,
			NULL );

    iFunction = iFunctionFindWithFunction(cPFunction);

    bPrint = !bFunctionPrint(iFunction);
    FunctionSetPrint( iFunction, bPrint );
    XUMFChangeMessPrint( iFunction, bPrint );

		/* Insert code to possibly update the TEXT Widget */

}







/*
 *	XUMFSetAllWithin
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Flip the state of all children of sContainerName.
 */
XtCallbackProc
XUMFSetAllWithin( Widget wNew, char *sContainerName, caddr_t PData2 )
{
Widget		wCont, wCur;
Widget		*wPChildren;
int		iChildren, i;
BOOL		bState, bNewState;
char		*cPName;

    wCont = WcFullNameToWidget( wNew, sContainerName );

    XtVaGetValues( wCont,
			XtNnumChildren, &iChildren,
			XtNchildren, &wPChildren,
			NULL );

    for ( i=0; i<iChildren; i++ ) {
	wCur = wPChildren[i];
	XUMFChangeToggleState( wCur, XUMF_SET );
    }

}



/*
 *	XUMFResetAllWithin
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Flip the state of all children of sContainerName.
 */
XtCallbackProc
XUMFResetAllWithin( Widget wNew, char *sContainerName, caddr_t PData2 )
{
Widget		wCont, wCur;
Widget		*wPChildren;
int		iChildren, i;
BOOL		bState, bNewState;
char		*cPName;

    wCont = WcFullNameToWidget( wNew, sContainerName );

    XtVaGetValues( wCont,
			XtNnumChildren, &iChildren,
			XtNchildren, &wPChildren,
			NULL );

    for ( i=0; i<iChildren; i++ ) {
	wCur = wPChildren[i];
	XUMFChangeToggleState( wCur, XUMF_RESET );
    }

}



/*
 *	XUMFFilePopup
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Popup the popup of (wNew).
 */
XtCallbackProc
XUMFFilePopup( Widget wNew, char *sActionName, caddr_t PData2 )
{
Widget		wPopup;

    if ( !XtIsSubclass(wNew,coreWidgetClass) ) {
	DFATAL(( "Widget is not a subclass of Core\n" ));
    }
    if ( wNew->core.num_popups != 1 ) {
	DFATAL(( "Wrong number of popups of: %s\n", XtName(wNew) ));
    }

    wPopup = wNew->core.popup_list[0];
    XtPopup( wPopup, FALSE );

}



/*
 *	XUMFRebuildMessBuffer
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Loop through all the messages, and for all those
 *	that have bPrint on, write them into the character
 *	buffer in (*cPPBuffer).
 */
XtCallbackProc
XUMFRebuildMessBuffer( Widget wNew, char *sActionName, caddr_t PData2 )
{
MESS		mCur;
STRING		sLine;
char		*cPTail, *cPBuffer, *cPTemp;
int		iLen;

    if ( GcPMessBuffer != NULL ) FREE( GcPMessBuffer );
    GcPMessBuffer = NULL;
    mCur = mMessLoop();
    while ( mMessNext(&mCur) ) {
	if ( bMessPrint(mCur) ) {
	    XUMFAppendMessToText(mCur);
	}
    }
    XUMFUpdateTextWidget();

}



/*
 *	XUMFTextWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Save the Widget which maintains the command line interface.
 *	The argument to this function must be the
 *	name of the Action that is used by the TextWidget
 *	to insert strings into the Text Widget.
 */
XtCallbackProc
XUMFTextWidgetRegister( Widget wNew, caddr_t PData1, caddr_t PData2 )
{
    GwText = wNew;
}




/*
 *	XUMFFileContainerWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 *
 */
XtCallbackProc
XUMFFileContainerWidgetRegister( Widget wNew, char *sActionName, caddr_t PData2)
{
    for ( SiCurrentFile=0; SiCurrentFile<iFunctionFileCount(); 
			SiCurrentFile++ ) {
	WcCreateNamedChildren( wNew, FILE_WIDGET_NAME );
    }
}





/*
 *	XUMFFileWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Save the Widget which maintains the command line interface.
 *	The argument to this function must be the
 *	name of the Action that is used by the TextWidget
 *	to insert strings into the Text Widget.
 */
XtCallbackProc
XUMFFileWidgetRegister( Widget wNew, char *sActionName, caddr_t PData2 )
{
STRING		sName;

		/* Define the name of the TOGGLE Widget to be the */
		/* function name */

    XtVaSetValues( wNew, XtNlabel, 
		(XtArgVal) sFunctionFileFilename(SiCurrentFile), NULL );

    WcCreateNamedPopups( wNew, FUNCTION_POPUP_NAME );

}



/*
 *	XUMFFunctionShellWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 *
 */
XtCallbackProc
XUMFFunctionShellWidgetRegister( Widget wNew, char *sActionName, caddr_t PData2)
{
    XtVaSetValues( wNew,
			XtNtitle, 
			(XtArgVal) sFunctionFileFilename(SiCurrentFile),
			NULL );
}




/*
 *	XUMFToggleContainerWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Save the Widget which maintains the command line interface.
 *	The argument to this function must be the
 *	name of the Action that is used by the TextWidget
 *	to insert strings into the Text Widget.
 */
XtCallbackProc
XUMFToggleContainerWidgetRegister( Widget wNew, char *sActionName, caddr_t PData2 )
{
    for ( SiCurrentFunction=0; SiCurrentFunction<iFunctionCount(); 
			SiCurrentFunction++ ) {
	if ( strcmp( sFunctionFilename(SiCurrentFunction),
			sFunctionFileFilename(SiCurrentFile) ) == 0 ) {
	    WcCreateNamedChildren( wNew, TOGGLE_WIDGET_NAME );
	}
    }

}





/*
 *	XUMFToggleWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Save the Widget which maintains the command line interface.
 *	The argument to this function must be the
 *	name of the Action that is used by the TextWidget
 *	to insert strings into the Text Widget.
 */
XtCallbackProc
XUMFToggleWidgetRegister( Widget wNew, char *sActionName, caddr_t PData2 )
{
STRING		sName;

		/* Define the name of the TOGGLE Widget to be the */
		/* function name */

    XtVaSetValues( wNew,
		XtNlabel, (XtArgVal) sFunctionFunction(SiCurrentFunction),
		XtNstate, (XtArgVal) bFunctionPrint(SiCurrentFunction),
		NULL );


}




/******************************************************************************
*   MAIN function
 *
 *	Author:	Christian Schafmeister (1991)
******************************************************************************/

void
main( int argc, char *argv[] )
{   
STRING		sAppClass;
Widget		appShell;
XWMHints	xwmhHints;
BOOL		bUseStartup;
char		c;
extern	char*	optarg;
extern	int	optind;
int		i;



    BasicsInitialize();

		/* Initialize X-Windows stuff */

    strcpy( sAppClass, argv[0] );
    if (islower(sAppClass[0])) sAppClass[0] = cUpper(sAppClass[0]);
    else if (islower(sAppClass[1])) sAppClass[1] = cUpper(sAppClass[1]);

	/*
	 * initialize first letter to make class, or first two if
	 * first is already capitalized, or don't worry about it.
	 */

    
    	/*  -- Intialize Toolkit creating the application shell */
    appShell = XtInitialize( 
	argv[0], sAppClass,		/* app name and class */
	NULL, 0, 			/* description of cmd line options */
	&argc, argv 
    );


		/* Parse command line arguments */

		/* -s {file_of_on_functions} */
		/* -i {input_file} if not included, get from stdin */
		/* All files after first arguments contain function names */

    GbInputFile = FALSE;
    while ( (c = getopt( argc, argv, "m:s:" )) != (char)EOF ) {
	switch (c) {
	    case 's':
		strcpy( GsInitialOnFunctionFilename, optarg );
		break;
	    case 'm':
		strcpy( GsInputFilename, optarg );
		printf( "Reading input from: %s\n", GsInputFilename );
		GbInputFile = TRUE;
		break;
	    default:
		printf( "Usage: %s [-m input] [-s on_file] line_files\n",
			argv[0] );
		exit(1);
	}
    }
    for ( i=optind; i<argc; i++ ) {
	XUMFReadLineFile(argv[i]);
    }


    GwTopWidget = appShell;
    GxtacApp = XtWidgetToApplicationContext(appShell);

		/* Get the descriptor for the input file */
		/* If the filename was not supplied then use STDIN=0 */

    if ( GbInputFile ) {
	GiInputDescriptor = open( GsInputFilename, O_RDONLY );
    } else GiInputDescriptor = 0;
    GfInput = fdopen( GiInputDescriptor, "r" );


    GxtiiId = XtAppAddInput( GxtacApp, GiInputDescriptor,
			XtInputReadMask,
			XUMFInputAvailable, 0 );


#define	REGCB(n,r) WcRegisterCallback( GxtacApp, n, r, NULL );
    REGCB( "XUMFTextWidgetRegister", XUMFTextWidgetRegister );
    REGCB( "XUMFToggleWidgetRegister", XUMFToggleWidgetRegister );
    REGCB( "XUMFToggleContainerWidgetRegister", 
		XUMFToggleContainerWidgetRegister );
    REGCB( "XUMFFileWidgetRegister", XUMFFileWidgetRegister );
    REGCB( "XUMFFileContainerWidgetRegister", 
		XUMFFileContainerWidgetRegister );
    REGCB( "XUMFToggleNotify", XUMFToggleNotify );
    REGCB( "XUMFFilePopup", XUMFFilePopup );
    REGCB( "XUMFFunctionShellWidgetRegister",
		XUMFFunctionShellWidgetRegister );
    REGCB( "XUMFSetAllWithin", XUMFSetAllWithin );
    REGCB( "XUMFResetAllWithin", XUMFResetAllWithin );
    REGCB( "XUMFRebuildMessBuffer",
		XUMFRebuildMessBuffer );

	/*  -- Register all Athena and Public widget classes */

    XpRegisterAll( GxtacApp );


	/*  -- Create widget tree below toplevel shell using Xrm database */
    WcWidgetCreation( appShell );

	/*  -- Realize the widget tree and enter the main application loop */

    XtRealizeWidget( appShell );

	/* Set up the window manager hints */

    xwmhHints.flags = InputHint;
    xwmhHints.input = TRUE;
    XSetWMHints( XtDisplay(appShell),
		 XtWindow(appShell),
		 &xwmhHints );
    XFlush( XtDisplay(appShell) );


    XtMainLoop ( );

}
