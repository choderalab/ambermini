/*
 *	File:	xaParmEditor.c
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
 *
 *		This file manages the AMBER parameter files viewing,
 *		ediiting and MOD file creation.
 *
 *		This Edit Parameter code uses the ATHENA set of
 *		Widgets.
 *
 *		It maintains all of the callback functions that 
 *		are called when the user clicks buttons etc.
 *
 *
 * 	NOTE:	In order for these routines to work the programmer
 *		must make the Widget that is the root of the UNIT
 *		editor have the name in the macro UNITEDITORWIDGETNAME!!!!!
 */


#define	PARMEDITORCLASS		transientShellWidgetClass


#include 	<X11/IntrinsicP.h>
#include 	<X11/StringDefs.h>
#include 	<X11/cursorfont.h>
#include 	<X11/Shell.h>
#include	<X11/keysym.h>


#include	"basics.h"

#include        "classes.h"

#include	"varArray.h"

#include	"xaTools.h"
#include	"xaLeap.h"
#include	"xaCommand.h"
#include	"xaAtomParmTable.h"
#include	"xaBondParmTable.h"
#include	"xaAngleParmTable.h"
#include	"xaTorsionParmTable.h"
#include	"xaHBondParmTable.h"
#include        "../Wc/WcCreate.h"

typedef	struct	{
	Widget		wTop;
	PARMSET		psParmSet;
} PARMEDITORt;


/**************************************************************************
	Static Stuff

*/

static	STRING		SsParmSet;
static	PARMSET		SpsParmSet;
static	XtAppContext	SxtacApp = NULL;
static  VARARRAY 	SvaParmEditors = NULL;


extern void XAIPTNewPopupTable( Widget, Widget, PARMSET, VFUNCTION );
extern void XATPTNewPopupTable( Widget, Widget, PARMSET, VFUNCTION );
extern void XATPTPopupTable( Widget, Widget, PARMSET, VFUNCTION );
extern void XAIPTPopupTable( Widget, Widget, PARMSET, VFUNCTION );
/*=========================================================================
	Private Functions


*/

/*
 *	zXAPEAddParmEditor
 *
 *	Add a new ParmEditor to the collection.
 *	The Shell which contains the UnitEditor is used to identify
 *	the ParmEditor from all the others.
 */
static void
zXAPEAddParmEditor( Widget wShell )
{
int		i;
PARMEDITORt	peNew;
PARMEDITORt	*pePNew;

		/* If there are no ParmEditors yet then create the */
		/* VARARRAY to store them */

    if ( SvaParmEditors == NULL ) {
	SvaParmEditors = vaVarArrayCreate(sizeof(PARMEDITORt));
    }

		/* Now try to find an empty place to put it */

    for ( i=0; i<iVarArrayElementCount(SvaParmEditors); i++ ) {
	if ( PVAI(SvaParmEditors,PARMEDITORt,i)->wTop == NULL ) break;
    }

		/* If no empty space was found then add one */

    if ( iVarArrayElementCount(SvaParmEditors) >= i ) {
	VarArrayAdd( SvaParmEditors, (GENP)&peNew );
	pePNew = PVAI(SvaParmEditors,PARMEDITORt,
			iVarArrayElementCount(SvaParmEditors)-1);
    } else {
	pePNew = PVAI(SvaParmEditors,PARMEDITORt,i);
    }

		/* Set up the new ParmEditor */

    pePNew->wTop = wShell;
    pePNew->psParmSet = SpsParmSet;
    
}



#define	zwXAPETop(w) (zuePXAPEParmEditorForSubWidget(w)->wTop)
#define zpsXAPEParmSet(w) (zuePXAPEParmEditorForSubWidget(w)->psParmSet)


/* 
 *	zuePXAPEParmEditorForSubWidget - Find the TOP Widget for the dialog.
 * 
 *	Do this by climbing up the Widget hierarchy until we find 
 *	a Widget of the proper class 
 *
 */
 
static PARMEDITORt *
zuePXAPEParmEditorForSubWidget(Widget wSub)
{
Widget 		wTop;
int		i;
PARMEDITORt	*pePCur;

		/* First find the TOP Widget for the UNITEDITOR */
		/* Do this by climbing up the Widget hierarchy */
		/* until we find a Widget of the proper class */

    wTop = wSub;
    while ( wTop != NULL ) {
	if ( XtClass(wTop) == PARMEDITORCLASS ) break;
	wTop = XtParent(wTop);
    }
    		/* Now find the PARMEDITOR that is associated with */
		/* the TOP Widget */

    pePCur = PVAI(SvaParmEditors,PARMEDITORt,0);
    for ( i=0; i<iVarArrayElementCount(SvaParmEditors); i++ ) {
	if ( pePCur->wTop == wTop ) return(pePCur);
	pePCur++;
    }

    DFATAL(( "Could not find a ParmEditor" ));
    return(NULL);	/* for lint */
}


/*
 *	zXAPEEditDone
 *	Does nothing at this point
 */
static void
zXAPEEditDone(Widget wTop)
{
} 

/*-------------------------------------------------------------------------
	Public Functions -  Callbacks


*/


/*
 *	XAPERaiseTable
 *
 *	Raises the table the is to be edited (get around the lack of
 *	control offered by the window manager
 *
 */
 
static XtCallbackProc
XAPERaiseTable( Widget wWid, caddr_t PAppData, caddr_t PArg )
{
    XRaiseWindow( XtDisplay( wWid ), XtWindow( wWid ));
    return(NULL);
}


/*
 *	XAPEEditAtomParameters
 *
 *	Create a TABLE to edit the Atom Parameters.
 */
static XtCallbackProc
XAPEEditAtomParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{

  /* 
   *  create the table if there are parameters 
   *	or if the user chooses to create it for adding parms 
   */
  if ( !iParmSetTotalAtomParms( zpsXAPEParmSet( wCur ))) {
    WcCreateNamedChildren( XtParent( wCur ), "aParmNoParmsShell" );
  } else {
    MESSAGE(( "Creating Atom Parameter Table..." ));		
    XAAPTPopupTable( wCur, zwXAPETop(wCur), zpsXAPEParmSet( wCur ) ,
	  zXAPEEditDone );
  }
  return(NULL);
}

/*
 *	XAPENewtAtomParameters
 *
 *	Create a TABLE to edit the Atom Parameters.
 */
XtCallbackProc
XAPENewAtomParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget		wParmEdit;
  
  MESSAGE(( "Creating Atom Parameter Table..." ));		
  wParmEdit = XtParent( XtParent( XtParent( wCur )));		
  XAAPTNewPopupTable( wCur, zwXAPETop( wParmEdit ), 
		     zpsXAPEParmSet( wParmEdit ), zXAPEEditDone );
  return(NULL);
}


/*
 *	XAPEEditBondParameters
 *
 *	Create a TABLE to edit the Bond Parameters.
 */
XtCallbackProc
XAPEEditBondParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  
  /* create the table if there are parameters or*/
  /* if the user chooses to create it for adding parms */

  if ( !iParmSetTotalBondParms( zpsXAPEParmSet( wCur ))) {
    WcCreateNamedChildren( XtParent( wCur ), "bParmNoParmsShell" );
  } else {
    MESSAGE(( "Creating Bond Parameter Table..." ));		
    XABPTPopupTable( wCur, zwXAPETop(wCur), zpsXAPEParmSet( wCur ),
		    zXAPEEditDone );
  }
  return(NULL);
}

/*
 *	XAPENewBondParameters
 *
 *	Create a TABLE to edit the Bond Parameters.
 */
XtCallbackProc
XAPENewBondParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget	wParmEdit;
  
  MESSAGE(( "Creating Bond Parameter Table..." ));
  wParmEdit = XtParent( XtParent( XtParent( wCur )));		
  XABPTNewPopupTable( wCur, zwXAPETop( wParmEdit ), 
		     zpsXAPEParmSet( wParmEdit ), zXAPEEditDone );
  return(NULL);
}

/*
 *	XAPEEditAngleParameters
 *
 *	Create a TABLE to edit the Angle Parameters.
 */
XtCallbackProc
XAPEEditAngleParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  
  /*
   * create the table if there are parameters or
   * if the user chooses to create it for adding parms
   */
  if ( !iParmSetTotalAngleParms( zpsXAPEParmSet( wCur ))) {
    WcCreateNamedChildren( XtParent( wCur ), "vParmNoParmsShell" );
  } else {
    MESSAGE(( "Creating Angle Parameter Table..." ));		
    XAVPTPopupTable( wCur, zwXAPETop(wCur), zpsXAPEParmSet( wCur ),
		    zXAPEEditDone );
  }
  return(NULL);
}

/*
 *	XAPENewAngleParameters
 *
 *	Create a TABLE to edit the Angle Parameters.
 */
XtCallbackProc
XAPENewAngleParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget		wParmEdit;
  
  MESSAGE(( "Creating Angle Parameter Table..." ));		
  wParmEdit = XtParent( XtParent( XtParent( wCur )));		
  XAVPTNewPopupTable( wCur, zwXAPETop( wParmEdit ), 
		     zpsXAPEParmSet( wParmEdit ), zXAPEEditDone );
  return(NULL);
}

/*
 *	XAPEEditTorsionParameters
 *
 *	Create a TABLE to edit the Torsion parameters.
 */
XtCallbackProc
XAPEEditTorsionParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{

  /*
   * create the table if there are parameters or
   *  if the user chooses to create it for adding parms
   */
  if ( !iParmSetProperCount( zpsXAPEParmSet( wCur ) )) {
    WcCreateNamedChildren( XtParent( wCur ), "tParmNoParmsShell" );
  } else {
    MESSAGE(( "Creating Proper Torsion Parameter Table..." ));		
    XATPTPopupTable( wCur, zwXAPETop(wCur), zpsXAPEParmSet( wCur ),
		    zXAPEEditDone );
  }
  return(NULL);
}

/*
 *	XAPENewTorsionParameters
 *
 *	Create a TABLE to edit the Torsion parameters.
 */
XtCallbackProc
XAPENewTorsionParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget		wParmEdit;

  MESSAGE(( "Creating Proper Torsion Parameter Table..." ));		
  wParmEdit = XtParent( XtParent( XtParent( wCur )));		
  XATPTNewPopupTable( wCur, zwXAPETop( wParmEdit ), 
		     zpsXAPEParmSet( wParmEdit ), zXAPEEditDone );
  return(NULL);
}


/*
 *	XAPEEditImproperParameters
 *
 *	Create a TABLE to edit the Improper Torsion parameters.
 */
XtCallbackProc
XAPEEditImproperParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  /*
   * create the table if there are parameters or
   * if the user chooses to create it for adding parms
   */
  if ( !iParmSetImproperCount( zpsXAPEParmSet( wCur ) )) {
    WcCreateNamedChildren( XtParent( wCur ), "iParmNoParmsShell" );
  } else {
    MESSAGE(( "Creating Improper Torsion Parameter Table..." )); 
    XAIPTPopupTable( wCur, zwXAPETop(wCur), zpsXAPEParmSet( wCur ),
		    zXAPEEditDone );
  }
  return(NULL);
}

/*
 *	XAPENewImproperParameters
 *
 *	Create a TABLE to edit the Improper Torsion parameters.
 */
XtCallbackProc
XAPENewImproperParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget 		wParmEdit;
  
  MESSAGE(( "Creating Improper Torsion Parameter Table..." ));		
  wParmEdit = XtParent( XtParent( XtParent( wCur )));		
  XAIPTNewPopupTable( wCur, zwXAPETop( wParmEdit ), 
		     zpsXAPEParmSet( wParmEdit ), zXAPEEditDone );
  return(NULL);
}

/*
 *	XAPEEditHBondParameters
 *
 *	Create a TABLE to edit the Torsion parameters.
 */
XtCallbackProc
XAPEEditHBondParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  /*
   * create the table if there are parameters or
   * if the user chooses to create it for adding parms
   */
  if ( !iParmSetTotalHBondParms( zpsXAPEParmSet( wCur ))) {
    WcCreateNamedChildren( XtParent( wCur ), "hParmNoParmsShell" );
  } else {
    MESSAGE(( "Creating Hydrogen Bond Parameter Table..." ));		
    XAHPTPopupTable( wCur, zwXAPETop(wCur), zpsXAPEParmSet( wCur ),
		    zXAPEEditDone );
  }
  return(NULL);
}

/*
 *	XAPENewHBondParameters
 *
 *	Create a TABLE to edit the Torsion parameters.
 */
XtCallbackProc
XAPENewHBondParameters( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget		wParmEdit;
  
  MESSAGE(( "Creating Hydrogen Bond Parameter Table..." ));		
  wParmEdit = XtParent( XtParent( XtParent( wCur )));		
  XAHPTNewPopupTable( wCur, zwXAPETop( wParmEdit ), 
		     zpsXAPEParmSet( wParmEdit ), zXAPEEditDone );
  return(NULL);
}

/*------------- Initialization functions */

/*
 *	XAPEPopupParmEditor
 *
 *	Popup a Parameter Set editor.
 */
void
XAPEPopupParmEditor( Widget wOrig, char *sName, char *sParmSet, 
	PARMSET psParmSet )
{
  /* if the parameter set is already being edited,
   * then do not open a new on to be edited
   * if it is not being edited,
   * then make it so that it will not be edited again
   */
  if( bParmSetBeingEdited(psParmSet ) == TRUE ) {
    MESSAGE(( "Second request to edit a parameter set.\n" ));
    VPDISPLAY(( "The parameter set %s is already being edited!\n", sParmSet ));
    return;
  }
  
  ParmSetEditing(psParmSet, TRUE);

  strcpy( SsParmSet, sParmSet );
  SpsParmSet = psParmSet;
  
  MESSAGE(( "Creating PARMEDITOR parent: %s\n", XtName(wOrig) ));
  MESSAGE(( "Creating widget with name: %s\n", sName ));
  WcCreateNamedChildren( wOrig, sName );
}



/*
 *	XAPERegisterParmEditor
 *
 *	Callback to register the ParmEditor.
 */
XtCallbackProc
XAPERegisterParmEditor( Widget wShell, caddr_t PAppData, caddr_t PArg )
{
STRING		sTitle;

		/* Define the title of the UnitEditor */

    sprintf( sTitle, "XLEaP: ParmSet Editor: %s", SsParmSet );
    XtVaSetValues( wShell,
			XtNtitle, (XtArgVal) sTitle,
			XtNiconName, (XtArgVal) SsParmSet,
			NULL );


		/* Create the new PARMEDITORt */

    zXAPEAddParmEditor( wShell );
    return(NULL);


}


/*
 *	XAPEEditDone
 *
 *	Once editing of the parameters is done
 *	this routine is called.
 *	Turns off the ParmSet BeingEdited flag
 */
XtCallbackProc
XAPEEditDone( Widget wWidget, caddr_t PAppData, caddr_t PArg)
{
    ParmSetEditing(SpsParmSet, FALSE);
    return(NULL);
}

/*
 *	XAPESetAtomLabel
 *
 *	Set the label to tell how many Atom Parms there are.
 */
XtCallbackProc
XAPESetAtomLabel( Widget wWidget, caddr_t PAppData, caddr_t PArg)
{
STRING	sLabel;
   
    sprintf( sLabel, "Atom Parameters (%d)", 
			iParmSetTotalAtomParms( SpsParmSet ));
    XtVaSetValues( wWidget, "label", (XtArgVal) sLabel, NULL );
    return(NULL);
}

/*
 *	XAPESetBondLabel
 *
 *	Set the label to tell how many Bond Parms there are.
 */
XtCallbackProc
XAPESetBondLabel( Widget wWidget, caddr_t PAppData, caddr_t PArg)
{
STRING	sLabel;
   
    sprintf( sLabel, "Bond Parameters (%d)", 
		iParmSetTotalBondParms( SpsParmSet ));
    XtVaSetValues( wWidget, "label", (XtArgVal) sLabel, NULL );
    
    return(NULL);
}


/*
 *	XAPESetAngleLabel
 *
 *	Set the label to tell how many Angle Parms there are.
 */
XtCallbackProc
XAPESetAngleLabel( Widget wWidget, caddr_t PAppData, caddr_t PArg)
{
STRING	sLabel;
   
    sprintf( sLabel, "Angle Parameters (%d)", 
		iParmSetTotalAngleParms( SpsParmSet ));
    XtVaSetValues( wWidget, "label", (XtArgVal) sLabel, NULL );
    return(NULL);
}


/*
 *	XAPESetTorsionLabel
 *
 *	Set the label to tell how many Torsion Parms there are.
 */
XtCallbackProc
XAPESetTorsionLabel( Widget wWidget, caddr_t PAppData, caddr_t PArg)
{
STRING	sLabel;
   
    sprintf( sLabel, "Proper Torsion Parameters (%d)", 
		iParmSetProperCount( SpsParmSet ));
    XtVaSetValues( wWidget, "label", (XtArgVal) sLabel, NULL );
    return(NULL);
}

/*
 *	XAPESetAtomLabel
 *
 *	Set the label to tell how many Improper Torsion Parms there are.
 */
XtCallbackProc
XAPESetImproperLabel( Widget wWidget, caddr_t PAppData, caddr_t PArg)
{
STRING	sLabel;

    sprintf( sLabel, "Improper Torsion Parameters (%d)", 
		iParmSetImproperCount( SpsParmSet ));
    XtVaSetValues( wWidget, "label", (XtArgVal) sLabel, NULL );
    return(NULL);
}


/*
 *	XAPESetHBondLabel
 *
 *	Set the label to tell how many Hydrogen Bond Parms there are.
 */
XtCallbackProc
XAPESetHBondLabel( Widget wWidget, caddr_t PAppData, caddr_t PArg)
{
STRING	sLabel;
   
    sprintf( sLabel, "Atom Parameters (%d)", 
		iParmSetTotalHBondParms( SpsParmSet ));
    XtVaSetValues( wWidget, "label", (XtArgVal) sLabel, NULL );
    return(NULL);
}



/*
 *	XAPEInitialize
 *
 *	Register all of the callback functions with Wcl.
 */
void
XAPEInitialize( XtAppContext app )
{
#define	WCLREG( n, f )	WcRegisterCallback( app, n, (XtCallbackProc)f, NULL )


    SxtacApp = app;

    WCLREG( "XAPERegisterParmEditor", XAPERegisterParmEditor );
    WCLREG( "XAPEEditAtomParameters", XAPEEditAtomParameters );
    WCLREG( "XAPEEditBondParameters", XAPEEditBondParameters );
    WCLREG( "XAPEEditAngleParameters", XAPEEditAngleParameters );
    WCLREG( "XAPEEditTorsionParameters", XAPEEditTorsionParameters );
    WCLREG( "XAPEEditImproperParameters", XAPEEditImproperParameters );
    WCLREG( "XAPEEditHBondParameters", XAPEEditHBondParameters );
    WCLREG( "XAPENewAtomParameters", XAPENewAtomParameters );
    WCLREG( "XAPENewBondParameters", XAPENewBondParameters );
    WCLREG( "XAPENewAngleParameters", XAPENewAngleParameters );
    WCLREG( "XAPENewTorsionParameters", XAPENewTorsionParameters );
    WCLREG( "XAPENewImproperParameters", XAPENewImproperParameters );
    WCLREG( "XAPENewHBondParameters", XAPENewHBondParameters );
    WCLREG( "XAPESetAtomLabel", XAPESetAtomLabel );
    WCLREG( "XAPESetBondLabel", XAPESetBondLabel );
    WCLREG( "XAPESetAngleLabel", XAPESetAngleLabel );
    WCLREG( "XAPESetTorsionLabel", XAPESetTorsionLabel );
    WCLREG( "XAPESetImproperLabel", XAPESetImproperLabel );
    WCLREG( "XAPESetHBondLabel", XAPESetHBondLabel );
    WCLREG( "XAPERaiseTable", XAPERaiseTable );
    WCLREG( "XAPEEditDone", XAPEEditDone );
    
}



