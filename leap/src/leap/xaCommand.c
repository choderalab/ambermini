/*
 *	File:	xaCommand.c
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
 *		Callbacks and routines for handling the command
 *		line interface using ATHENA Widgets.
 *
 */


#include 	<X11/IntrinsicP.h>
#include	<X11/StringDefs.h>
#include 	<X11/Xatom.h>
#include 	"../Xmu/Atoms.h"
#include 	"../Xmu/Misc.h"
#include	"../Xraw/List.h"
#include	"../Xraw/AsciiText.h"
#include    "../Xraw/AsciiSrc.h"
#include	"../Xraw/Label.h"
#include	"../Xraw/Toggle.h"
#include	<X11/keysym.h>
#include	<X11/Shell.h>

#include	"basics.h"

#include	"classes.h"

#include	"dictionary.h"
#include	"leap.h"

#include	"block.h"

#include	"parser.h"
#include	"xaLeap.h"
#include	"xaTools.h"

#include	"xaCommand.h"

#include        "../Wc/WcCreate.h"

#define	EMPTY_LIST	"      "
#define	XACEMPTYSTRING	"";

extern BLOCK	GbCommand;

extern void ParseBlock( BLOCK, RESULTt* );
extern void XAPEPopupParmEditor( Widget, char *, char *, PARMSET );
extern void XAUEPopupUnitEditor( Widget, char *, char *, UNIT );
	    
/*
 *	zXACDisplayChildren
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Display the children of this widget.
 *	This is used for debugging.
 */
static void
zXACDisplayChildren( Widget w, int iIndent )
{
Widget		*wPChildren;
int		iChildren, i, j;

    iIndent += 3;
    wPChildren = NULL;
    iChildren = 0;
    XtVaGetValues( w,
			"children", &wPChildren,
			"numChildren", &iChildren,
			NULL );
    if ( wPChildren != NULL ) {
	for ( i=0; i<iChildren; i++ ) {
	    PRINTF(( "%lX  ", *wPChildren ));
	    for ( j=0; j<iIndent; j++ ) { 
		if ( j%3 == 0 ) {
		    PRINTF_NO_PREFIX(( "|" ));
		} else {
		    PRINTF_NO_PREFIX(( " " )); 
		}
	    }
	    PRINTF_NO_PREFIX(( "Child:%-20s\n", XtName(*wPChildren) ));
	    zXACDisplayChildren( *wPChildren, iIndent );
	    wPChildren++;
	}
    }
    if ( XtIsSubclass(w,coreWidgetClass) ) {
    	i = w->core.num_popups;
    	for (wPChildren = w->core.popup_list; i; i--, wPChildren++) {
	    PRINTF(( "%lX  ", *wPChildren ));
	    for ( j=0; j<iIndent; j++ ) { 
		if ( j%3 == 0 ) {
		    PRINTF_NO_PREFIX(( "*" ));
		} else {
		    PRINTF_NO_PREFIX(( " " ));
		}
	    }
	    PRINTF_NO_PREFIX(( "Popup:%-20s\n", XtName(*wPChildren) ));
	    zXACDisplayChildren( *wPChildren, iIndent );
	}
    }
}




/*
 *	zXACDisplayWidgetTree
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Display the tree of Widgets whos root is w.
 *	Used for debugging.
 */
static void
zXACDisplayWidgetTree( Widget w )
{
    PRINTF(( "\n" ));
    PRINTF(( "%lX  %-20s\n", w, XtName(w) ));
    zXACDisplayChildren( w, 0 );
}





/*
 *	zXACSetWidgetValue
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Set the value of the Widget.
 *
 *	If the Widget is a Text widget, then set the XtNstring resource
 *	of its source Widget.
 *	If the Widget is a Label widget then set the XtNlabel resource.
 */
static void
zXACSetWidgetValue( Widget wWid, char *sValue )
{
Widget		wSource;

	if (XtClass(wWid) == asciiTextWidgetClass) {
		wSource = NULL;
		XtVaGetValues( wWid,
			XtNtextSource, &wSource,
			NULL );
		if ( wSource == NULL )
	    		DFATAL(( "Could not get Text source\n" ));
		XtVaSetValues( wSource,
			XtNstring, (XtArgVal) sValue,
			NULL );
	} else if (XtClass(wWid) == labelWidgetClass) {
		XtVaSetValues( wWid,
			XtNlabel, (XtArgVal) sValue,
			NULL );
	} else {
      		PRINTF (("zXACSetWidgetValue: I don't work with %s class\n",
				XtClass(wWid)->core_class.class_name ));
	}
}



/*
 *	zXACGetWidgetValue
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Get the value of the Widget.
 *	Return a pointer to the string that represents the Widgets
 *	current value.
 *
 *	If the Widget is a List then return the currently selected
 *	element.
 *	If the Widget is a AsciiTextWidget return the string.
 *	If the Widget is a Label then return the XtNlabel resource.
 *
 */
static void
zXACGetWidgetValue( Widget wWid, char **cPPValue )
{
XawListReturnStruct	*xawlrsPElement;
Widget			wSource;

    *cPPValue = XACEMPTYSTRING;
    if ( XtClass(wWid) == listWidgetClass ) {
	xawlrsPElement = XawListShowCurrent(wWid);
	if ( xawlrsPElement->list_index != XAW_LIST_NONE ) {
	    *cPPValue = xawlrsPElement->string;
	}
    } else if ( XtClass(wWid) == asciiTextWidgetClass ) {
	wSource = NULL;
	XtVaGetValues( wWid,
			XtNtextSource, &wSource,
			NULL );
	if ( wSource == NULL ) {
	    DFATAL(( "Could not get Text source\n" ));
	}
	XtVaGetValues( wSource,
			XtNstring, cPPValue,
			NULL );
    } else if ( XtClass(wWid) == labelWidgetClass ) {
	XtVaGetValues( wWid,
			XtNlabel, cPPValue,
			NULL );
    } else if ( XtClass(wWid) == toggleWidgetClass ) {
        *cPPValue = XawToggleGetCurrent(wWid);
      } else {
      PRINTF (("zXACGetWidgetValue: I don't work with %s class\n",
		XtClass(wWid)->core_class.class_name ));
    }

}



/*
 *	zXACCopyAndExpandWidgetNames
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Copy the string from sSource to sDest, while expanding
 *	Widget names within parenthesis to the Widget values.
 *
 */
static void
zXACCopyAndExpandWidgetNames( Widget wRef, char *sSource, char *sDest )
{
char	*cPSource, *cPDest, *cPVal, *cPWidgetName;
STRING	sWidgetName;

    cPSource = sSource;
    cPDest = sDest;

    while ( *cPSource ) {

			/* Look for parenthesis containing Widget */
			/* names, copy out the Widget names and */
			/* find those Widgets values */
	if ( *cPSource == '(' ) {
	    cPSource++;
	    cPWidgetName = sWidgetName;
	    while ( *cPSource && *cPSource != ')' ) {
		if ( *cPSource == ' ' ) {
		    cPSource++;
		    continue;
		}
		*cPWidgetName = *cPSource;
		cPWidgetName++;	
		cPSource++;
	    }
	    if ( !*cPSource ) {
		DFATAL(( "Missing close parenthesis in callback for: %s\n",
				XtName(wRef) ));
	    }
	    *cPWidgetName = '\0';
	    cPSource++;

			/* Get the Widget value */
	    XACWidgetNameToValue( wRef, sWidgetName, &cPVal );
	    MESSAGE(( "Converted Widget name: %s to value: %s\n",
			    sWidgetName, cPVal ));

			/* Copy the Widget value into cPDest */

	    while ( *cPVal ) {
		*cPDest = *cPVal;
		cPDest++;
		cPVal++;
	    }
	} else {

			/* Otherwise just copy the character into cPDest */
	    *cPDest = *cPSource;
	    cPDest++;
	    cPSource++;
	}
    }
    *cPDest = '\0';
}


/*
 *	ziXACCompareTwoDirectoryEntries
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return -1 if sA is less than sB, 0 if they are the same,
 *	and +1 if sA is larger than sB.
 *
 *	Order is determined by:
 *		a is less than b if a is not a directory and b is.
 *		a is less than b if a proceeds b in an alphabetical
 *			list.
 */
static int
ziXACCompareTwoDirectoryEntries( char *sA, char *sB )
{
char	cA, cB;

    if ( strlen(sA)==0 || strlen(sB)==0 ) {
	DFATAL(( "Zero length directory name\n" ));
    }

		/* First check if one is a directory and the other isnt */

    cA = sA[strlen(sA)-1];
    cB = sB[strlen(sB)-1];
    if ( cA == '/' && cB != '/' ) return(-1);
    if ( cA != '/' && cB == '/' ) return(1);

		/* Make sure that "." comes before ".." */

    if ( sA[0] == '.' && sB[0] == '.' ) {
	if ( strlen(sA) == 1 )  {
	    return(1);
	} else {
	    return(-1);
	}
    }
    return(strcmp( sA, sB ));
}






/*
 *	zbXACTextProcessOneChar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a single character to the BLOCK.
 *	If the block was ended then return TRUE.
 */
void
zXACTextProcessOneChar( char c )
{
RESULTt		rResult;

    VPDISPLAY(( "%c", c ));
    if ( bBlockAddChar( GbCommand, c ) ) {
	ParseBlock( GbCommand, &rResult );
	if ( rResult.iCommand == CQUIT ) {
	    exit(0);

		/* If the USER want's to EDIT a UNIT then */
		/* popup an editor */
	} else if ( rResult.iCommand == CEDIT ) {
		/************************************
		* Davids Changes 		    *
		* If PARMSET then open a parmEditor *
		* else open a unitEditor 	    *
		*************************************/
	    if(iObjectType(rResult.oObject) == PARMSETid) {
	      XAPEPopupParmEditor( GwTopWidget, 
				  "parmEditor",
				  rResult.sVariable,
				  rResult.oObject );
	    } else {
	      XAUEPopupUnitEditor( GwTopWidget,
				  "unitEditor", 
				  rResult.sVariable, 
				  rResult.oObject );
	    }
	} else if ( rResult.iCommand == CVERBOSITY ) {
	  Widget wVerbo = (Widget)NULL;

	  if (       rResult.sVariable[0] == (char)0){
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l0");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) True, NULL);

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l1");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) False, NULL);
			
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l2");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) False, NULL);
	  }else if (rResult.sVariable[0] == (char)1){
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l0");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) False, NULL);

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l1");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) True, NULL);
			
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l2");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) False, NULL);
	  }else if (rResult.sVariable[0] == (char)2){
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l0");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) False, NULL);

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l1");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) False, NULL);
			
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l2");
	    XtVaSetValues( wVerbo, XtNstate, (XtArgVal) True, NULL);
	  }
	}
	BlockEmpty( GbCommand );
	VPDISPLAY(( "> " ));
	return;
    }
    if ( c=='\n' ) 
	VPDISPLAY(( "> " ));
}





/*
 *---------------------------------------------------------------
 *
 *	Callbacks
 *
 */


/*
 *	xtcpXACCmdWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Save the Widget which maintains the command line interface.
 *	The argument to this function must be the
 *	name of the Action that is used by the TextWidget
 *	to insert strings into the Text Widget.
 */
XtCallbackProc
xtcpXACCmdWidgetRegister( Widget wNew, char *sActionName, caddr_t PData2 )
{
    MESSAGE(( "Registered the Cmd with the application\n" ));
    GwCmd = wNew;
    return(NULL);
}


XtCallbackProc	
xtcpXACCmdAbout( Widget wNew, caddr_t PData1, caddr_t PData2 )
{
  Widget text;

  text = WcFullNameToWidget(wNew, "*value");
  XtVaSetValues(text, XtNtype, (XtArgVal) XawAsciiFile, XtNstring, 
	(XtArgVal) "XaLeap", NULL);
  return(NULL);
}

XtCallbackProc
xtcpXACInitToggle( Widget wNew, caddr_t sArgs, caddr_t PData2 )
{
  char	*ss;
  
  ss = (char*)XtNewString((char*)sArgs);
  XtVaSetValues(wNew, XtNradioData, (XtArgVal) ss, NULL);
  return(NULL);
}

XtCallbackProc
xtcpXACReplaceShell( Widget wNew, caddr_t sArgs, caddr_t PData2 )
{
  Widget wTopShell;
  Position x,y;

  wTopShell = (Widget) WcFullNameToWidget( wNew, "/");
  XtVaGetValues( wTopShell, XtNx, &x, XtNy, &y, NULL);
  x += 300;
  y += 100;
  XtVaSetValues(wNew, XtNx, (XtArgVal) x, XtNy, (XtArgVal) y, NULL);
  return(NULL);
}




/*
 *	xtcpXACDisplayTree
 *
 *	Author:	Christian Schafmeister (1991)
 *
 */
XtCallbackProc
xtcpXACDisplayTree( Widget wNew, char *sRootName, caddr_t PData2 )
{
Widget		wRoot;

		/* Direct all VPx output to the GwCmd Widget */
    PushCurrentPrintSink(GiCmdSink);

    wRoot = WcFullNameToWidget( wNew, sRootName );
    zXACDisplayWidgetTree(wRoot);

    PopCurrentPrintSink();
    return(NULL);
}



/*
 *	xtcpXACInitializeList
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build a list of variable names and place it in the
 *	Widgets value field.
 *
 *	The variable names that are chosen are those that
 *	coorespond to OBJEKTs that are of one of the types in
 *	sObjektTypes.  sObjektTypes is a string of characters,
 *	each character MUST coorespond to an OBJEKT id.
 *	The OBJEKT ids are listed in 'objekt.h'
 *	An example of an sObjektTypes string is: 'UMRA' which
 *	will coorespond to all UNITs, MOLECULEs, RESIDUEs, and ATOMs.
 *	'O' will coorespond to all variables.
 */
XtCallbackProc
xtcpXACInitializeList( Widget wNew, char *sArgs, caddr_t PData2 )
{
String		*saList;
int		i, iNumArgs, iStrings;
DICTIONARY	dVariables;
DICTLOOP	dlVar;
OBJEKT		oVar;
char		*cPType, *sObjektTypes;
Widget		wWid;
ARGt		aArgs;

    PushCurrentPrintSink(GiCmdSink);

		/* Split up the arguments into seperate strings */

    XATParseCallbackArgs( sArgs, aArgs, &iNumArgs );
    if ( iNumArgs != 2 ) {
	DFATAL(( "Need 2 arguments for xtcpXACInitializeList widget: %s\n", 
			XtName(wNew) ));
    }

    wWid = WcFullNameToWidget( wNew, aArgs[0] );

    if ( !wWid ) {
	DFATAL(( "Callback:xtcpXACInitializeList arg#1 expected ListWidget\n"));
    }

    sObjektTypes = aArgs[1];

		/* Count the number of Variables the coorespond */
		/* to the types in sObjektTypes */
    iStrings = 0;
    dVariables = dVariablesDictionary();
    dlVar = ydlDictionaryLoop(dVariables);
    while ( (oVar = (OBJEKT)yPDictionaryNext( dVariables, &dlVar )) ) {
	for ( cPType=sObjektTypes;*cPType;cPType++ ) {
	    if ( bObjectInClass( oVar, *cPType ) ) 
		iStrings++;
	}
    }

		/* Allocate memory for the array of pointers to Strings */
		/* If there is nothing in the list then create a single */
		/* element filled with spaces */

    if ( iStrings == 0 ) {
	MALLOC( saList, String*, sizeof(String) );
	saList[0] = EMPTY_LIST;
	iStrings = 1;
    } else {
    	MALLOC( saList, String*, sizeof(String)*iStrings );

		/* Now initialize the pointers to the Strings */

    	i = 0;
    	dlVar = ydlDictionaryLoop(dVariables);
    	while ( (oVar = (OBJEKT)yPDictionaryNext( dVariables, &dlVar )) ) {
	    for ( cPType=sObjektTypes;*cPType;cPType++ ) {
	    	if ( bObjectInClass( oVar, *cPType ) ) {
		    saList[i++] = sDictLoopKey(dlVar);
		}
	    }
	}
    }

		/* Put it into the List Widget */

    XawListChange( wWid, saList, iStrings, 0, TRUE );

    PopCurrentPrintSink();
    return(NULL);
}



/*
 *	xtcpXACDestroyList
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Deallocate the memory used to store the List elements.
 *	This must be called for every XALInitializeListOfxxxx
 */
XtCallbackProc
xtcpXACDestroyList( Widget wNew, char *sArgs, caddr_t PData2 )
{
String		*saList;
STRING		sWidget;
Widget		wWid;

    PushCurrentPrintSink(GiCmdSink);

		/* Get the Widget with the list to destroy  */

    (void) WcCleanName( sArgs, sWidget );
    wWid = WcFullNameToWidget( wNew, sWidget );

		/* Get the pointer to the start of the list */

    saList = NULL;

    XtVaGetValues( wWid,
			XtNlist, &saList,
			NULL );

		/* Inform the List that its value is being changed */

    XawListUnhighlight( wWid );
    XtVaSetValues( wWid,
			XtNlist, (XtArgVal) NULL,
			XtNnumberStrings, (XtArgVal) 0,
			NULL );

    if ( saList != NULL ) 
	FREE(saList);

    PopCurrentPrintSink();
    return(NULL);
}



/*
 *	xtcpXACParseCommand
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback Arguments:
 *		#1 -	Command to parse.
 *
 *	Causes LEaP to parse the text passed.
 */
XtCallbackProc
xtcpXACParseCommand( Widget wCur, caddr_t PArg, caddr_t PData2 )
{
char		*cPArg;
STRING		sArg, sCommand;

		/* Direct all VPx output to the GwCmd Widget */

    PushCurrentPrintSink(GiCmdSink);

    MESSAGE(( "About to parse: %s\n", PArg ));

    strcpy( sArg, PArg );
    cPArg = sArg;
    while ( *cPArg != '\0' ) {
	if ( *cPArg == ';' ) *cPArg = '\n';
        cPArg++;
    }
    *cPArg = '\0';
    cPArg = sArg;

		/* Expand Widget names to Widget values */

    zXACCopyAndExpandWidgetNames( wCur, cPArg, sCommand );

    cPArg = sCommand;
    while ( *cPArg != '\0' ) {
	zXACTextProcessOneChar( *cPArg );
	cPArg++;
    }

    PopCurrentPrintSink();
    return(NULL);
}



/*
 *	xtcpXACCopyValueTo
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback Arguments:
 *		#1 -	Widget to copy value of current Widget into.
 *
 *	Copy the value of the current Widget into the Widget
 *	with the name in sWidget.
 */
XtCallbackProc
xtcpXACCopyValueTo( Widget wCur, char *sWidget, caddr_t PData2 )
{
Widget		wWid;
char		*cPValue;

    PushCurrentPrintSink(GiCmdSink);

    wWid = WcFullNameToWidget( wCur, sWidget );
    zXACGetWidgetValue( wCur, &cPValue );
    zXACSetWidgetValue( wWid, cPValue );

    PopCurrentPrintSink();
    return(NULL);
}





/*
 *	xtcpXACInitializeWorkingDirectoryName
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback Arguments:
 *		#1 -	Name of Widget to copy the current path to.
 *
 *	Copy the current working directory path name
 *	into the widget in (sWidget).
 */
XtCallbackProc
xtcpXACInitializeWorkingDirectoryName( Widget wCur, char *sWidget, 
	caddr_t PData2 )
{
Widget		wWid;
STRING		sPath;

    PushCurrentPrintSink(GiCmdSink);

    wWid = WcFullNameToWidget( wCur, sWidget );

		/* Get the current working pathname */

    BasicsCurrentWorkingDirectory( sPath );

    zXACSetWidgetValue( wWid, sPath );

    PopCurrentPrintSink();
    return(NULL);
}


/*
 *	xtcpXACInitializeDirectoryList
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback Arguments:
 *		#1 -	Name of List Widget to initialize
 *		#2 -	Name of Widget whose value contains the path
 *			to list.
 *
 *	Create a sorted list of directory entry names, sorting by
 *	alphabetical order, and sifting directory names to the top.
 *	Suffix directory names with '/'.
 *
 */
XtCallbackProc
xtcpXACInitializeDirectoryList( Widget wCur, char *sArgs, caddr_t PData2 )
{
Widget		wList, wValue;
char		*cPPath;
int		iNumber, i, iNumArgs, iSize;
FILESTATUSt	fsStatus;
STRING		*saPNames;
char		**cPaList;
char		*cPName;
ARGt		aArgs;
STRING		sName;

    PushCurrentPrintSink(GiCmdSink);


		/* Get the arguments to the callback */


    XATParseCallbackArgs( sArgs, aArgs, &iNumArgs );

    wList = WcFullNameToWidget( wCur, aArgs[0] );
    wValue = WcFullNameToWidget( wCur, aArgs[1] );

    if ( !wList ) {
	DFATAL(( "Callback:%s, arg#1 need: ListWidget\n",
			"xtcpXACInitializeDirectoryList" ));
    }
    if ( !wValue ) {
	DFATAL(( "Callback:%s, arg#2 need: PathWidget\n",
			"xtcpXACInitializeDirectoryList" ));
    }

    zXACGetWidgetValue( wValue, &cPPath );

    /* 
     *  List the directory
     */
    SysdependDirectoryList (cPPath, &saPNames, &iNumber);

    for ( i=0; i<iNumber; i++ ) {
        /*
         *  if it's a dir, append '/'
         */
        strcpy( sName, cPPath );
        strcat( sName, saPNames[i]);
        fsStatus = fsBasicsFileStatus( sName );
        if ( fsStatus.fMode & FILEDIRECTORY ) {
            strcat( saPNames[i], "/" );
        }
    }

    /* 
     *  Sort the list into alphabetical order, w/ directories
     *	at the top of the list
     */
    qsort( saPNames, iNumber, sizeof(STRING), 
		(int (*) (const void *, const void *) )ziXACCompareTwoDirectoryEntries );

    /*
     *  Reformulate as an array of pointers to char
     *	TODO - do this right off, & adjust ziXACCompareTwoDirectoryEntries
     */
    MESSAGE(( "Directory list\n" ));
    MALLOC( cPaList, char**, sizeof(char*) * iNumber );
    for ( i=0; i<iNumber; i++ ) {
	iSize = strlen(saPNames[i])+2;
	MALLOC( cPName, char*, iSize );
	strcpy( cPName, saPNames[i] );
	MESSAGE(( "List[%d](size=%4d)=|%s|\n", i, iSize, cPName ));
	cPaList[i] = cPName;
    }
    FREE( saPNames );
    MESSAGE(( "---------------\n" ));

		/* Set the value of the List Widget */

    XawListChange( wList, cPaList, iNumber, 0, TRUE );
    /* TODO - XawListChange() doesn't free prev list */

    PopCurrentPrintSink();
    return(NULL);
}





/*
 *	xtcpXACDestroyDirectoryList
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback Arguments:
 *		#1 -	Name of List Widget which contains the list
 *			to destroy.
 *
 *	Destroy the Directory name list in the Widget.
 *
 */
XtCallbackProc
xtcpXACDestroyDirectoryList( Widget wCur, char *sArgs, caddr_t PData2 )
{
char		**cPaList;
int		iNumber, i;
Widget		wList;

    PushCurrentPrintSink(GiCmdSink);

    wList = WcFullNameToWidget( wCur, sArgs );
    if ( !wList ) {
	DFATAL(( "Callback:%s arg#1 expected: ListWidget\n",
			"xtcpXACDestroyDirectoryList" ));
    }

		/* Get the List of the current Widget */

    XtVaGetValues( wList,
			XtNlist, &cPaList,
			XtNnumberStrings, &iNumber,
			NULL );

		/* Say that the list is being changed to nothing */

    XawListUnhighlight(wList);
    XtVaSetValues( wList,
			XtNlist, (XtArgVal) NULL,
			XtNnumberStrings, (XtArgVal) 0,
			NULL );

		/* FREE each string */

    for ( i=0; i<iNumber; i++ ) {
	FREE( cPaList[i] );
    }
    FREE( cPaList );

    PopCurrentPrintSink();
    return(NULL);
}





/*
 *	xtcpXACDirectoryPick
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback Arguments:
 *		#1 -	Text Widget for file name
 *		#2 -	Label Widget for path name
 *
 *	When the user selects an item out of a Directory list,
 *	get the item, and if it is a directory, clear the
 *	Text Widget, destroy the current list, move into 
 *	the new one directory, and regenerate
 *	a directory list, then update the current path name.
 *	If the item picked is a file name then copy it into the Text
 *	Widget.
 */
XtCallbackProc
xtcpXACDirectoryPick( Widget wCur, char *sArgs, caddr_t PData2 )
{
ARGt		aArgs;
int		iNumArgs, i;
char		*cPPick;
char		*cPCwd;
Widget		wText, wPath;
STRING		sPath, sTempArgs, sPick;
char		**cPaOldList;
int		iNumberInOldList;

    PushCurrentPrintSink(GiCmdSink);

		/* Get the arguments */

    XATParseCallbackArgs( sArgs, aArgs, &iNumArgs );
    wText = WcFullNameToWidget( wCur, aArgs[0] );
    wPath = WcFullNameToWidget( wCur, aArgs[1] );

    if ( !wText ) {
	DFATAL(( "Callback:xtcpXACDirectoryPick, arg#1 need: TextWidget\n" ));
    }
    if ( !wPath ) {
	DFATAL(( "Callback:xtcpXACDirectoryPick, arg#2 need: PathWidget\n" ));
    }
		/* Get the item picked */

    zXACGetWidgetValue( wCur, &cPPick );
    strcpy( sPick, cPPick );

		/* Check if the pick is a directory */

    if ( cPPick[strlen(sPick)-1] == '/' ) {

			/* Clear the Text Widget */

	zXACSetWidgetValue( wText, "" );

			/* Get the current List */

	XtVaGetValues( wCur,
			XtNlist, &cPaOldList,
			XtNnumberStrings, &iNumberInOldList,
			NULL );

			/* Change the path Label Widget to reflect the */
			/* change in path */

	zXACGetWidgetValue( wPath, &cPCwd );
	strcpy( sPath, cPCwd );
	if ( strcmp( sPick, "./" ) != 0 ) {
	    if ( strcmp( sPick, "../" ) == 0 ) {

			/* Move up one level */
		if ( strlen(sPath) > 0 ) {
		    for ( i=strlen(sPath)-2; i>=0; i-- ) {
			if ( sPath[i] == '/' ) {
			    sPath[i+1] = '\0';
			    break;
			}
		    }
		}
	    } else {
		strcat( sPath, sPick );
	    }
	}
	zXACSetWidgetValue( wPath, sPath );

			/* Now generate a new List */

	strcpy( sTempArgs, "*" );
	strcat( sTempArgs, XtName(wCur) );
	strcat( sTempArgs, "," );
	strcat( sTempArgs, aArgs[1] );

			/* Execute the callback to generate the new list */

	xtcpXACInitializeDirectoryList( wCur, sTempArgs, NULL );

			/* Destroy the old list */

	for ( i=0; i<iNumberInOldList; i++ ) {
	     FREE( cPaOldList[i] );
	}
	FREE( cPaOldList );

    } else {

		/* Just copy the value of the pick into the Text Widget */

	zXACSetWidgetValue( wText, sPick );
    }
    PopCurrentPrintSink();
    return(NULL);
}





/*
 **********************************************************************
 *
 *	Action procedures.
 *
 */






/*
 *	xtapXACTextProcessRubout
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Action procedure to process backspace and delete.
 */
XtActionProc
xtapXACTextProcessRubout( Widget wText, XEvent *xePEvent, String *sPArgs, 
	Cardinal *cPNum )
{

		/* Direct all VPx output to the GwCmd Widget */

    PushCurrentPrintSink(GiCmdSink);

    if ( cBlockLastWrittenChar(GbCommand) &&
		cBlockLastWrittenChar(GbCommand) != '\n' ) {
	bBlockRemoveChar(GbCommand);
	XtCallActionProc( wText, "delete-previous-character",
				xePEvent, NULL, 0 );
    }

    PopCurrentPrintSink();
    return(NULL);
}








/*
 *	xtapXACTextProcessRuboutLine
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Action procedure to process backspace and delete.
 */
XtActionProc
xtapXACTextProcessRuboutLine( Widget wText, XEvent *xePEvent, String *sPArgs, 
	Cardinal *cPNum )
{

		/* Direct all VPx output to the GwCmd Widget */

    PushCurrentPrintSink(GiCmdSink);

    while ( cBlockLastWrittenChar(GbCommand) &&
		cBlockLastWrittenChar(GbCommand) != '\n' ) {
	if ( bBlockRemoveChar(GbCommand) ) {
	    XtCallActionProc( wText, "delete-previous-character",
				xePEvent, NULL, 0 );
	}
    }

    PopCurrentPrintSink();
    return(NULL);
}







/*
 *	xtapXACTextProcessKeys
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Action procedure to process key strokes for
 *	the command line interface.
 */
XtActionProc
xtapXACTextProcessKeys( Widget wText, XEvent *xePEvent, String *sPArgs, 
	Cardinal *cPNum )
{
#define	BUFSIZE		100
char			sBuffer[BUFSIZE];
KeySym			ksSym;
int			i, iNum;

    if (xePEvent->type != KeyPress)
	return NULL;

		/* Direct all VPx output to the GwCmd Widget */

    PushCurrentPrintSink(GiCmdSink);

	/* Read the keypress event and print it to the display */

    iNum = XLookupString( &(xePEvent->xkey), sBuffer,
			 BUFSIZE-1, &ksSym, NULL );
    sBuffer[iNum] = '\0';
    if ( (ksSym >= XK_space && ksSym <= XK_asciitilde) ) {
	for ( i=0; i<iNum; i++ ) 
		zXACTextProcessOneChar( sBuffer[i] );
    } else if ( ksSym == XK_Return ) {
	zXACTextProcessOneChar( '\n' );
    }

    PopCurrentPrintSink();
    return(NULL);
}






/*
 *	xtapXACTextProcessCancel
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Action procedure to process CNTRL-C which
 *	causes the current block to be abandoned.
 */
XtActionProc
xtapXACTextProcessCancel( Widget wText, XEvent *xePEvent, String *sPArgs, 
	Cardinal *cPNum )
{
		/* Direct all VPx output to the GwCmd Widget */

    PushCurrentPrintSink(GiCmdSink);

    XtCallActionProc( wText, "insert-char", xePEvent, NULL, 0 );
    XtCallActionProc( wText, "newline", xePEvent, NULL, 0 );
    BlockEmpty( GbCommand );
    VPDISPLAY(( "> " ));

		/* Tell the system that a CONTROL-C has been hit */

    BasicsSetInterrupt();

    PopCurrentPrintSink();
    return(NULL);
}


/*
 *	xtapXACTextProcessInsertSelection
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Action procedure to process CNTRL-C which
 *	causes the current block to be abandoned.
 */
XtActionProc
xtapXACTextProcessInsertSelection( Widget wText, XEvent *xePEvent, 
			String *sPArgs, Cardinal *cPNum )
{
Atom		aSelection;
int		i, iBuffer, iBytes;
char		*cPData;

		/* Direct all VPx output to the GwCmd Widget */

    PushCurrentPrintSink(GiCmdSink);

    MESSAGE(( "About to insert selection\n" ));

    XmuInternStrings( XtDisplay(wText), sPArgs, (Cardinal)1, &aSelection );

    switch (aSelection) {
      case XA_CUT_BUFFER0: iBuffer = 0; break;
      case XA_CUT_BUFFER1: iBuffer = 1; break;
      case XA_CUT_BUFFER2: iBuffer = 2; break;
      case XA_CUT_BUFFER3: iBuffer = 3; break;
      case XA_CUT_BUFFER4: iBuffer = 4; break;
      case XA_CUT_BUFFER5: iBuffer = 5; break;
      case XA_CUT_BUFFER6: iBuffer = 6; break;
      case XA_CUT_BUFFER7: iBuffer = 7; break;
      default:	       iBuffer = -1;
    }
    MESSAGE(( "Getting selection from cut buffer: %d\n", iBuffer ));
    if ( iBuffer >= 0) {
	cPData = XFetchBuffer( XtDisplay(wText), &iBytes, iBuffer );
	MESSAGE(( "Got %d bytes\n", iBytes ));
	for ( i=0; i<iBytes; i++ ) {
	    zXACTextProcessOneChar( *cPData );
	    cPData++;
	}

	if ( iBytes > 0 )
	  XFree((char*)(cPData - iBytes));

      } else {
	DFATAL(( "Illegal cut buffer: %s\n", sPArgs[0] ));
    }

    PopCurrentPrintSink();
    return(NULL);
}

/*
 *	xManageDelete
 *	Author: David A. Rivkin
 *	Catches WM_PROTCOLS messages and closes windows as needed
 *
 */
Atom wmDeleteWindow = 0L;

void
xManageDelete(Widget w, XEvent *event, String *params, Cardinal *num_params)
{
    if (event->type == ClientMessage) {
        if (event->xclient.data.l[0] != wmDeleteWindow) {
            XBell(XtDisplay(w), 0);
        } else {
            /* handle WM_DELETE_WINDOW... */
            if ( XtIsApplicationShell( w )) 
		exit(0);
            else
		XtDestroyWidget( w );
        }
    }
    
}





static	XtActionsRec	SxtarCommandActions[] = {
	{ "xtapXACTextProcessKeys", (VFUNCTION)xtapXACTextProcessKeys },
	{ "xtapXACTextProcessRubout", (VFUNCTION)xtapXACTextProcessRubout },
	{ "xtapXACTextProcessRuboutLine", 
				(VFUNCTION)xtapXACTextProcessRuboutLine },
	{ "xtapXACTextProcessCancel", (VFUNCTION)xtapXACTextProcessCancel },
	{ "xtapXACTextProcessInsertSelection", 
				(VFUNCTION)xtapXACTextProcessInsertSelection},
	{ "xManageDelete", (VFUNCTION)xManageDelete }

};







/*
 *=================================================================
 *
 *	Public routines
 *
 */


/*
 *	XACWidgetNameToValue
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Get the value of the Widget with the name sName.
 *
 */
void
XACWidgetNameToValue( Widget wRef, char *sName, char **cPPValue )
{
Widget			wWid, WcFullNameToWidget();

    wWid = WcFullNameToWidget( wRef, sName );
    zXACGetWidgetValue( wWid, cPPValue );
}

	


/*
 *	XACInitialize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Register all of the callbacks that
 *	this file makes available.
 */
void
XACInitialize( XtAppContext xtacApp )
{
#define	REGCB(n,r) WcRegisterCallback( xtacApp, n, (XtCallbackProc)r, NULL );


		/* Create the block that will contain the commands */

    GbCommand = bBlockCreate();
    
    REGCB( "xtcpXACDisplayTree", xtcpXACDisplayTree );
    REGCB( "xtcpXACInitializeList", xtcpXACInitializeList );
    REGCB( "xtcpXACDestroyList", xtcpXACDestroyList );
    REGCB( "xtcpXACParseCommand", xtcpXACParseCommand );
    REGCB( "xtcpXACCopyValueTo", xtcpXACCopyValueTo );
    REGCB( "xtcpXACInitializeWorkingDirectoryName", 
		xtcpXACInitializeWorkingDirectoryName );
    REGCB( "xtcpXACInitializeDirectoryList",
		xtcpXACInitializeDirectoryList );
    REGCB( "xtcpXACDestroyDirectoryList",
		xtcpXACDestroyDirectoryList );
    REGCB( "xtcpXACDirectoryPick", xtcpXACDirectoryPick );
    REGCB( "xtcpXACCmdWidgetRegister", xtcpXACCmdWidgetRegister );

    REGCB( "xtcpXACCmdAbout", xtcpXACCmdAbout ); /* V. Romanovski */
    REGCB( "xtcpXACInitToggle", xtcpXACInitToggle ); /* V. Romanovski */

    REGCB( "xtcpXACReplaceShell", xtcpXACReplaceShell ); /* V. Romanovski */


    XtAppAddActions( xtacApp, SxtarCommandActions,
    			XtNumber(SxtarCommandActions) );
			
}

#ifdef notdef
	  
#define BACKGROUND "tan"     /* look  xaCommand.rm4 SET_BUTTON */
#define SET        "maroon"  /* look  xaCommand.rm4 SET_BUTTON */

	  if (       rResult.sVariable[0] == (char)0){
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l0");

	    WcSetValueFromString( wVerbo, "background", SET);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "-30");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "-20");

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l1");
	    WcSetValueFromString( wVerbo, "background", BACKGROUND);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "20");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "40");

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l2");
	    WcSetValueFromString( wVerbo, "background", BACKGROUND);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "20");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "40");

	  }else if (rResult.sVariable[0] == (char)1){
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l1");
	    WcSetValueFromString( wVerbo, "background", SET);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "-30");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "-20");

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l0");
	    WcSetValueFromString( wVerbo, "background", BACKGROUND);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "20");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "40");

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l2");
	    WcSetValueFromString( wVerbo, "background", BACKGROUND);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "20");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "40");

	  }else if (rResult.sVariable[0] == (char)2){
	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l2");
	    WcSetValueFromString( wVerbo, "background", SET);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "-30");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "-20");

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l0");
	    WcSetValueFromString( wVerbo, "background", BACKGROUND);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "20");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "40");

	    wVerbo = WcFullNameToWidget ( GwTopWidget, "*l1");
	    WcSetValueFromString( wVerbo, "background", BACKGROUND);
	    WcSetValueFromString( wVerbo, "topShadowContrast", "20");
	    WcSetValueFromString( wVerbo, "bottomShadowContrast", "40");
	  }
#endif
