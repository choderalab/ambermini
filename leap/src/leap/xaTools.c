/*
 *	File:	xaTools.c
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
 *		Contains various routines for dealing with ATHENA Widgets.
 *		Define all bitmaps used by xaLeap within this file
 *		and initialize the widgets that use them using
 *		xtcpXATDefineBitmap.
 *
 */


#include 	<X11/Intrinsic.h>
#include	"basics.h"

#include	"dictionary.h"

#include	"xaLeap.h"
#include	"xaTools.h"
#include        "../Wc/WcCreate.h"

/*
 *------------------------------------------------------------
 *
 *	Global variables
 *
 */

static	DICTIONARY	SdBitmaps = NULL;


/*
 *------------------------------------------------------------
 *
 *	Handle printing to different Text Widgets
 *
 */


/*
 *      XATPrintStringToWidget
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Print the string to the Command line interface Text Widget.
 */
void
XATPrintStringToWidget( char *sOutput, Widget wWid )
{
XEvent		xeEvent;

    /* this gets called for 'nonevents' like pressing shift key */
    if (*sOutput == '\0')
	return;
    xeEvent.type = 0;
    XtCallActionProc( wWid, "insert-string", &xeEvent, &sOutput, 1 );
    XFlush( XtDisplay(wWid) );
}




/*
 *------------------------------------------------------------
 *
 *	Handle parsing callback arguments.
 *
 */
 
/*
 *	XATParseCallbackArgs
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Split up a string containing several arguments seperated by comments.
 *	Also remove leading spaces and tabs from each argument.
 */
void
XATParseCallbackArgs( char *sArgs, ARGt aArgs, int *iPNumArgs )
{
char	*cPArgs, *cPNew;
int	i;
BOOL	bLeadingSpace;

    MESSAGE(( "Splitting: %s\n", sArgs ));
    i = 0;
    cPArgs = sArgs;
    while ( *cPArgs ) {
	bLeadingSpace = TRUE;
	cPNew = aArgs[i];
	while ( *cPArgs && *cPArgs != ',' ) {
	    if ( *cPArgs == ' ' && bLeadingSpace ) {
		cPArgs++;
		continue;
	    }
	    bLeadingSpace = FALSE;
	    *cPNew = *cPArgs;
	    cPArgs++;
	    cPNew++;
	}
	*cPNew = '\0';
	MESSAGE(( "Split: |%s|\n", aArgs[i] ));
	i++;
	if ( !*cPArgs ) break;
	cPArgs++;
    }
    *iPNumArgs = i;
    MESSAGE(( "Split %d arguments\n", i ));
}



/*
 ***********************************************
 *
 *	Bitmap routines.
 *
 */

/*
 *	zXATAddBitmap
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a named bitmap to the DICTIONARY of bitmaps.
 */
static void
zXATAddBitmap( char *sName, char *cPData, int iWidth, int iHeight )
{
Pixmap	*pPBitmap;

		/* Allocate memory for the Pixmap and create it */
		/* from the data */
		
    MALLOC( pPBitmap, Pixmap*, sizeof(Pixmap));
    *pPBitmap = XCreateBitmapFromData( GdDisplay,
    					GwRootWindow,
					cPData, iWidth, iHeight );

		/* Now put the Pixmap into the DICTIONARY */

    if ( !SdBitmaps ) {
    	SdBitmaps = dDictionaryCreate();
    }
    
    DictionaryAdd( SdBitmaps, sName, (GENP)pPBitmap );
}




/*
 *==============================================================
 *
 *	Callbacks
 *
 */
 

/*
 *	xtcpXATDefineBitmap
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback that defines the bitmap of the Widget for which
 *	this callback was invoked.
 *	The callback takes two arguments:
 *		xtcpXATDefineBitmap( sResource, sBitmapName )
 *		sResource - name of the resource to set.
 *		sBitmapName - name of the bitmap to use.
 */
XtCallbackProc
xtcpXATDefineBitmap( Widget wWidget, char *cPArguments, GENP PData )
{
ARGt	aArguments;
int	iNum;
Pixmap	*pPBitmap;


    XATParseCallbackArgs( cPArguments, aArguments, &iNum );
    if ( iNum != 2 ) {
	DFATAL(( 
	  "xtcpXATDefineBitmap called with %d arguments, need (sRes,sBitmap)",
	  iNum ));
    }

    pPBitmap = (Pixmap*)yPDictionaryFind( SdBitmaps, aArguments[1] );
    if ( pPBitmap == NULL ) {
    	DFATAL(( "xtcpXATDefineBitmap called with unknown bitmap name: %s",
		aArguments[1] ));
    }
    XtVaSetValues( wWidget,
    			aArguments[0], (XtArgVal) *pPBitmap,
			NULL );
    return(NULL);
}




#include	"xaLogo.xbm"
#include        "exclaim.xbm"

/*
 *	XATInitializeTools
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Initialize the Tools routines.
 */
void
XATInitializeTools( XtAppContext xtacApp )
{
#define	REGCB(n,r) WcRegisterCallback( xtacApp, n, (XtCallbackProc)r, NULL );

		/* Define all of the bitmaps */
		
    zXATAddBitmap( "xaLogo", xaLogo_bits, xaLogo_width, xaLogo_height );
    zXATAddBitmap( "exclaim", exclaim_bits, exclaim_width, exclaim_height );


    REGCB( "xtcpXATDefineBitmap", xtcpXATDefineBitmap );
}

    
