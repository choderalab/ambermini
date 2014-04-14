/*
 *	File:	xaUnitEditor.c
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
 *		Manage a collection of UnitEditors which allow 
 *		the USER to edit a UNIT in a graphical window.
 *		This file maintains a VARARRAY of structures, one
 *		for each UnitEditor that is currently active.
 *
 *		This UnitEditor uses the ATHENA set of
 *		Widgets.
 *
 *		It maintains all of the callback functions that 
 *		are called when the user clicks buttons etc.
 *
 *		One of the functions of the Unit Editor is to invoke
 *		the MINIMIZER on a selected portion of the UNIT.
 *		If the user hits Ctrl-C within the TANK UnitEditor while
 *		the MINIMIZER is active the MINIMIZER will stop.
 *
 * 	NOTE:	In order for these routines to work the programmer
 *		must make the Widget that is the root of the UNIT
 *		editor have the name in the macro UNITEDITORWIDGETNAME!!!!!
 */

#define TRANSIENT_SHELL transientShellWidgetClass
#define TOPLEVEL_SHELL  topLevelShellWidgetClass
#define	UNITEDITORCLASS TOPLEVEL_SHELL


#include 	<X11/IntrinsicP.h>
#include 	<X11/StringDefs.h>
#include 	<X11/cursorfont.h>
#include 	<X11/Shell.h>
#include	<X11/keysym.h>
#include    "../Xraw/SmeBSB.h"

#include	"basics.h"
#include        "classes.h"
#include        "commands.h"
#include	"varArray.h"
#include	"minimizer.h"
#include 	"xTank.h"
#include        "select.h"
#include        "build.h"
#include        "model.h"
#include        "graphUtil.h"
#include        "leap.h"
#include	"xaTools.h"
#include	"xaLeap.h"
#include	"xaCommand.h"
#include	"xaUnitEditor.h"
#include	"xaAtomTable.h"
#include        "../Wc/WcCreate.h"
extern int		GiUnitEditors;

/*
 *----------------------------------------------------------------------
 *
 *	Private routines.
 *
 *		The following routines maintain the VARARRAY
 *		of currently active UNITEDITORt.
 *
 *		The caller can use ANY Widget that is in the subtree
 *		that defines the UNITEDITORt to access the UNITEDITORt
 *		fields for that UNITEDITORt.
 *
 */
		/* Store the Application context which is the same */
		/* over ALL Widgets in this application */

static	XtAppContext	SxtacApp            = NULL;

typedef	struct	{
 	Widget		wTop;
	Widget		wTank;
	Widget		wStatus;
	STRING		sVariable;
	UNIT		uUnit;
	MINIMIZER	mStrain;
	int		iPrintSink;
	XrmQuark        qName;              /* Vladimir Romanovski */
	BOOL            bTableExist;        /* Vladimir Romanovski */
} UNITEDITORt;


static	VARARRAY	SvaEditors = NULL;


		/* Temporary storage for data between */
		/* creation of UnitEditor and registration */
		/* of the shell Widget */
		/* UNIT to edit is VERY temporarily stored here between */
		/* the call to XAUEPopupUnitEditor and the call to */
		/* XAEUTankWidgetRegister when the TANK Widget is */
		/* registered */

static	STRING		SsUnit;
static	UNIT		SuUnit = NULL;



#define	UNITEDITORPREFIX	"Unit Editor: "

#define SELECTION_EDIT          ("*selEdit")



#define	zwXAUETop(w) (zuePXAUEUnitEditorForSubWidget(w)->wTop)
#define	zwXAUETank(w)	(zuePXAUEUnitEditorForSubWidget(w)->wTank)
#define	zwXAUEStatus(w)	(zuePXAUEUnitEditorForSubWidget(w)->wStatus)
#define	zsXAUEVariable(w) (zuePXAUEUnitEditorForSubWidget(w)->sVariable)
#define	zuXAUEUnit(w)	(zuePXAUEUnitEditorForSubWidget(w)->uUnit)
#define	zmXAUEMinimizer(w)	(zuePXAUEUnitEditorForSubWidget(w)->mStrain)
#define	ziXAUEPrintSink(w) (zuePXAUEUnitEditorForSubWidget(w)->iPrintSink)
#define	zXAUESetTank(w,d)	(zuePXAUEUnitEditorForSubWidget(w)->wTank=d)
#define	zXAUESetStatus(w,d)	(zuePXAUEUnitEditorForSubWidget(w)->wStatus=d)
#define	zXAUESetVariable(w,d) \
		strcpy(zuePXAUEUnitEditorForSubWidget(w)->sVariable,d)
#define	zXAUESetUnit(w,d)	(zuePXAUEUnitEditorForSubWidget(w)->uUnit=d )
#define	zXAUESetMinimizer(w,d) \
		(zuePXAUEUnitEditorForSubWidget(w)->mStrain=d )
#define	zXAUESetAcceptCallback(w,d) \
		(zuePXAUEUnitEditorForSubWidget(w)->fAcceptCallback=d )
#define	zXAUESetRejectCallback(w,d) \
		(zuePXAUEUnitEditorForSubWidget(w)->fRejectCallback=d )
#define	zXAUESetPrintSink(w,d) \
		(zuePXAUEUnitEditorForSubWidget(w)->iPrintSink=d )
#define	zXAUETableExist(w,d) \
		(zuePXAUEUnitEditorForSubWidget(w)->bTableExist=d )


/*
 *	zuePXAUEUnitEditorForSubWidget
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a pointer to the UnitEditor field in the
 *	collection of UnitEditors for the particular
 *	Widget.
 *
 *	The caller can provide ANY 
 */
static UNITEDITORt *
zuePXAUEUnitEditorForSubWidget( Widget wSub )
{
int		i;
UNITEDITORt	*uePCur;
Widget		wTop;

		/* First find the TOP Widget for the UNITEDITOR */
		/* Do this by climbing up the Widget hierarchy */
		/* until we find a Widget of the proper class */


    wTop = wSub;
    while ( wTop != NULL ) {
	if ( XtClass(wTop) == UNITEDITORCLASS ) break;
	wTop = XtParent(wTop);
    }

		/* Now find the UNITEDITOR that is associated with */
		/* the TOP Widget */

    uePCur = PVAI(SvaEditors,UNITEDITORt,0);
    for ( i=0; i<iVarArrayElementCount(SvaEditors); i++, uePCur++ ) {
	if ( uePCur->wTop == wTop ) 
	    return(uePCur);
    }

    DFATAL(( "Could not find a UnitEditor" ));
    return(NULL);	/* for lint */
}






/*
 *	zXAUEAddUnitEditor
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a new UnitEditor to the collection.
 *	The Shell which contains the UnitEditor is used to identify
 *	the UnitEditor from all the others.
 */
static void
zXAUEAddUnitEditor( Widget wShell )
{
int		i;
UNITEDITORt	ueNew;
UNITEDITORt	*uePNew;

		/* If there are no UnitEditors yet then create the */
		/* VARARRAY to store them */

    if ( SvaEditors == NULL ) {
	SvaEditors = vaVarArrayCreate(sizeof(UNITEDITORt));
    }

		/* Now try to find an empty place to put it */

    for ( i=0; i<iVarArrayElementCount(SvaEditors); i++ ) {
	if ( PVAI(SvaEditors,UNITEDITORt,i)->wTop == NULL ) break;
    }

		/* If no empty space was found then add one */

    if ( iVarArrayElementCount(SvaEditors) >= i ) {
	VarArrayAdd( SvaEditors, (GENP)&ueNew );
	uePNew = PVAI(SvaEditors,UNITEDITORt,
			iVarArrayElementCount(SvaEditors)-1);
    } else {
	uePNew = PVAI(SvaEditors,UNITEDITORt,i);
    }

		/* Set up the new UnitEditor */

    uePNew->wTop        = wShell;
    uePNew->qName       = XrmStringToQuark(SsUnit);
    uePNew->bTableExist = FALSE;

                /* To avoid a dependence on ordering in what 
                   tank and status widgets would be created   V.T.R. 1995*/ 

    uePNew->wTank       = NULL;
    uePNew->wStatus     = NULL;

}


/*
 *      XAUECheckUnitEdit
 *
 *      Author: Vladimir Romanovski (1994)
 *
 *      Check for existing table for units with the same name.
 *
 */

XtCallbackProc
XAUECheckUnitEdit( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget       wTop;
  XrmQuark     qTest;
  int          iCount;
  char		*sUnitName;
  UNITEDITORt  *unit;
  register int i;

  /* wCur is nemu button in Unit editor */

  wTop = (Widget)WcFullNameToWidget( wCur,"/");
  
  XtVaGetValues( wTop, XtNiconName, &sUnitName, NULL);
  
  qTest = XrmStringToQuark(sUnitName);
  
  iCount = iVarArrayElementCount(SvaEditors);
  
  for ( i = 0; i < iCount ; i++){
      
    unit = PVAI(SvaEditors,UNITEDITORt,i);

    if (unit->wTop  != NULL  &&	unit->qName == qTest &&	unit->bTableExist)
      {
	XtSetSensitive( wCur, False);
	break;
      }
  }
  return NULL;
}


  



/*
 *	zXAUERemoveUnitEditor
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy a UNITEDITOR referenced by a Widget that is
 *	in its hierarchy.
 */
static void
zXAUERemoveUnitEditor( Widget wSub )
{
UNITEDITORt	*uePEditor;

    uePEditor = zuePXAUEUnitEditorForSubWidget(wSub);

		/* Destroy the print sink for the Unit Editor */

    DestroyPrintSink(uePEditor->iPrintSink);

		/* Just set the Top Widget to NULL */

    uePEditor->wTop = NULL;

}




/*
 *	zXAUEEditDone
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Once editing of the properties of the ATOMs
 *	this routine is called to redisplay the UNIT.
 */
static void
zXAUEEditDone( Widget wTop )
{
Widget		wTank;

    wTank = zwXAUETank(wTop);

    TankRedisplayUnit((TANK)wTank);
}






/*
 *	bXAUERelaxCallback
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback to redisplay the UNIT every step of a minimization.
 *	Return TRUE to indicate that the minimization should continue.
 */
static BOOL
bXAUERelaxCallback( Widget wTank )
{
XEvent		xeEvent;
BOOL		bContinue;

		/* Redisplay the UNIT in the TANK */


    TankFastRedisplayUnit((TANK) wTank );

		/* If there is an event pending then handle it */

    bContinue = TRUE;
    while ( XtAppPending( SxtacApp ) ) {
	MESSAGE(( "There is an event pending\n" ));
	XtAppNextEvent( SxtacApp, &xeEvent );
	if ( bBasicsInterrupt() ) {
	    bContinue = FALSE;
	    BasicsResetInterrupt();
	}
	XtDispatchEvent( &xeEvent );
    }

    return(bContinue);
}




/*
 *-------------------------------------------
 *
 *	Callback routines called by DISPLAYERs
 */


/*
 *	XAUEDisplayerUpdate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Update the unit editor when this routine
 *	is called by a DISPLAYER.
 */
static void
XAUEDisplayerUpdate( DISPLAYER dDisp, UNIT uUnit, Widget wWidget )
{
Widget		wTank;

    wTank = zwXAUETank(wWidget);
    TankRedisplayUnit((TANK)wTank);
}


/*--------------------------------------------------*/

static void
CheckUnCheck( Widget wCur, BOOL bCurrent )       /* Vladimir Romanovski */
{
  XtVaSetValues( wCur, XtNstate, (XtArgVal) !bCurrent, NULL);
}



/*
 *-------------------------------------------------------------------
 *
 *	UnitEditor callbacks.
 *
 */




/*
 *	XAUEToggleShowNames
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEToggleShowNames( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
BOOL		bCurrent;
Widget		wTank;
  
    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

    wTank = zwXAUETank( wCur );

    bCurrent = bTankFlagsSet((TANK) wTank, TANKSHOWNAMES );
    CheckUnCheck(wCur,bCurrent);
    if ( bCurrent ) 
	TankResetFlags((TANK) wTank, TANKSHOWNAMES );
    else
	TankSetFlags((TANK) wTank, TANKSHOWNAMES );
    TankRedisplayUnit((TANK)wTank);

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEToggleShowPertNames
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEToggleShowPertNames( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  BOOL		bCurrent;
  Widget		wTank;
  
    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

  wTank = zwXAUETank( wCur );
  
  
  bCurrent = bTankFlagsSet((TANK) wTank, TANKSHOWPERTNAMES );
  CheckUnCheck(wCur,bCurrent);
  if ( bCurrent ) 
	TankResetFlags((TANK) wTank, TANKSHOWPERTNAMES );
  else
	TankSetFlags((TANK) wTank, TANKSHOWPERTNAMES );
  TankRedisplayUnit((TANK)wTank);

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEToggleShowResidues
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEToggleShowResidues( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
BOOL		bCurrent;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

    wTank = zwXAUETank( wCur );

    bCurrent = bTankFlagsSet((TANK) wTank, TANKSHOWRESIDUES );
    CheckUnCheck(wCur,bCurrent);
    if ( bCurrent ) 
	TankResetFlags((TANK) wTank, TANKSHOWRESIDUES );
    else
	TankSetFlags((TANK) wTank, TANKSHOWRESIDUES );
    TankRedisplayUnit((TANK)wTank);

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}
   



/*
 *	XAUEToggleShowTypes
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEToggleShowTypes( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
BOOL		bCurrent;
FLAGS		fFlag;
Widget		wTank;


    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

    wTank = zwXAUETank( wCur );


    fFlag = TANKSHOWTYPES;

    bCurrent = bTankFlagsSet((TANK) wTank, fFlag );
    CheckUnCheck(wCur,bCurrent);
    if ( bCurrent ) 
	TankResetFlags((TANK) wTank, fFlag );
    else
	TankSetFlags((TANK) wTank, fFlag );
    TankRedisplayUnit((TANK)wTank);

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEToggleShowPertTypes
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEToggleShowPertTypes( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
BOOL		bCurrent;
FLAGS		fFlag;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

    wTank = zwXAUETank( wCur );


    fFlag = TANKSHOWPERTTYPES;

    bCurrent = bTankFlagsSet((TANK) wTank, fFlag );
    CheckUnCheck(wCur,bCurrent);
    if ( bCurrent ) 
	TankResetFlags((TANK) wTank, fFlag );
    else
	TankSetFlags((TANK) wTank, fFlag );
    TankRedisplayUnit((TANK)wTank);

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEToggleShowCharges
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEToggleShowCharges( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
BOOL		bCurrent;
FLAGS		fFlag;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    fFlag = TANKSHOWCHARGES;

    bCurrent = bTankFlagsSet((TANK) wTank, fFlag );
    CheckUnCheck(wCur,bCurrent);
    if ( bCurrent ) 
	TankResetFlags((TANK) wTank, fFlag );
    else
	TankSetFlags((TANK) wTank, fFlag );
    TankRedisplayUnit((TANK)wTank);

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEToggleShowAxis
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEToggleShowAxis( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
BOOL		bCurrent;
FLAGS		fFlag;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    fFlag = TANKSHOWAXIS;

    bCurrent = bTankFlagsSet((TANK) wTank, fFlag );
    CheckUnCheck(wCur,bCurrent);
    if ( bCurrent ) 
	TankResetFlags((TANK) wTank, fFlag );
    else
	TankSetFlags((TANK) wTank, fFlag );
    TankRedisplayUnit((TANK)wTank);

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}





/*
 *	XAUEToggleShowBox
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEToggleShowBox( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
BOOL		bCurrent;
FLAGS		fFlag;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    fFlag = TANKSHOWBOX;

    bCurrent = bTankFlagsSet((TANK) wTank, fFlag );
    CheckUnCheck(wCur,bCurrent);
    if ( bCurrent ) 
	TankResetFlags((TANK) wTank, fFlag );
    else
	TankSetFlags((TANK) wTank, fFlag );
    TankRedisplayUnit((TANK)wTank);

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}




/*
 *	XAUESetDrawingElement
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUESetDrawingElement( Widget wCur, caddr_t PArg, caddr_t PAppData )
{
int		iElement;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

    wTank = zwXAUETank( wCur );


    iElement = iElementNumber(PArg);
    if ( iElement == NOELEMENT ) {
	DFATAL(( "ACTION: Illegal element in XAUESetDrawingElement: %s", 
								PArg ));
    } else {
	TankSetDrawingElement((TANK) wTank, iElement );
	MESSAGE(( "Set current drawing element to: %s\n", PArg ));
    }
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}










/*
 *	XAUEBuildExternals
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build the external coordinates for the UNIT.
 */
XtCallbackProc
XAUEBuildExternals( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
LOOP		lAtom, lSpan;
ATOM		aAtom, aStart;
UNIT		uUnit;
Widget		wTank;
int		iDum;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );

    uUnit = uTankUnit((TANK)wTank);

		/* Turn off ATOMNEEDSBUILD flags, they are only
			for internal use */

    ContainerResetAllAtomsFlags((CONTAINER) uUnit, ATOMNEEDSBUILD );

    if ( uUnit == NULL ) return NULL;

    MESSAGE(( "Building external coordinates for the UNIT\n" ));

	/* Try to build geometries for simple rings */

    BuildInternalsForSimpleRings( (CONTAINER)uUnit );


	/* Assign internal coordinates for all internals that */
	/* include atoms that need building */

    lAtom = lLoop((OBJEKT) uUnit, ATOMS );
    BuildInternalsUsingFlags( &lAtom, ATOMPOSITIONDRAWN, 0,
				ATOMNEEDSBUILD,
				ATOMPOSITIONDRAWN);

	/* Build spanning trees for all atoms that need building */
	/* and build external coordinates for those atoms */

    lAtom = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtom)) ) {
	if ( bAtomFlagsSet( aAtom, ATOMNEEDSBUILD ) ) {

			/* Look for a collision with an ATOM that has */
			/* already been built */
			/* Then start building from there */

	    lSpan = lLoop((OBJEKT) aAtom, SPANNINGTREE );
	    LoopDefineVisibleAtoms( &lSpan, ATOMNEEDSBUILD );
	    while ( oNext(&lSpan) );
	    if ( iLoopInvisibleCollisionCount(&lSpan) > 0 ) {
		aStart = aLoopLastCollisionAtom(&lSpan);
		lSpan = lLoop((OBJEKT) aStart, SPANNINGTREE );
	    } else {
		lSpan = lLoop((OBJEKT) aAtom, SPANNINGTREE );
	    }
	    LoopDefineVisibleAtoms( &lSpan, ATOMNEEDSBUILD );
	    iDum = 0;	/* for purify */
	    BuildExternalsUsingFlags( &lSpan, ATOMNEEDSBUILD, 0,
					ATOMPOSITIONKNOWN,
					ATOMNEEDSBUILD,
					&iDum, &iDum, &iDum, TRUE );
	}
    }

		/* Destroy all of the INTERNALs */

    lAtom = lLoop((OBJEKT) uUnit, ATOMS );
    BuildDestroyInternals( &lAtom );

    TankRedisplayUnit((TANK) wTank );
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}





/*
 *	XAUEAddHydrogensBuildExternals
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add missing hydrogens and build the external coordinates
 *	for the UNIT.
 */
XtCallbackProc
XAUEAddHydrogensAndBuildExternals( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
LOOP		lAtom, lSpan;
ATOM		aAtom, aStart;
UNIT		uUnit;
Widget		wTank;
int		iDum;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    uUnit = uTankUnit((TANK)wTank);
    if ( uUnit == NULL ) return NULL;

    MESSAGE(( "Building external coordinates for the UNIT\n" ));

	/* Add hydrogens */

    ModelAddHydrogens( uUnit );

	/* Try to build geometries for simple rings */

    BuildInternalsForSimpleRings( uUnit );

	/* Assign internal coordinates for all internals that */
	/* include atoms that need building */

    lAtom = lLoop((OBJEKT) uUnit, ATOMS );
    BuildInternalsUsingFlags( &lAtom, ATOMPOSITIONDRAWN, 0,
				ATOMNEEDSBUILD,
				ATOMPOSITIONDRAWN );

	/* Build spanning trees for all atoms that need building */
	/* and build external coordinates for those atoms */

    lAtom = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtom)) ) {
	if ( bAtomFlagsSet( aAtom, ATOMNEEDSBUILD ) ) {

			/* Look for a collision with an ATOM that has */
			/* already been built */
			/* Then start building from there */

	    lSpan = lLoop((OBJEKT) aAtom, SPANNINGTREE );
	    LoopDefineVisibleAtoms( &lSpan, ATOMNEEDSBUILD );
	    while ( oNext(&lSpan) );
	    if ( iLoopInvisibleCollisionCount(&lSpan) > 0 ) {
		aStart = aLoopLastCollisionAtom(&lSpan);
		lSpan = lLoop((OBJEKT) aStart, SPANNINGTREE );
	    } else {
		lSpan = lLoop((OBJEKT) aAtom, SPANNINGTREE );
	    }
	    LoopDefineVisibleAtoms( &lSpan, ATOMNEEDSBUILD );
	    iDum = 0;	/* for purify */
	    BuildExternalsUsingFlags( &lSpan, ATOMNEEDSBUILD, 0,
					ATOMPOSITIONKNOWN,
					ATOMNEEDSBUILD,
					&iDum, &iDum, &iDum, TRUE );
	}
    }

		/* Destroy all of the INTERNALs */

    lAtom = lLoop((OBJEKT) uUnit, ATOMS );
    BuildDestroyInternals( &lAtom );

    TankRedisplayUnit((TANK) wTank );
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}






/*
 *	XAUEStateDraw
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEStateDraw( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
Widget		wTank;

    wTank = zwXAUETank( wCur );

   TankSetState((TANK) wTank, TANKDRAW );
   return NULL;
}



/*
 *	XAUEStateErase
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEStateErase( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
Widget		wTank;

    wTank = zwXAUETank( wCur );

    TankSetState((TANK) wTank, TANKERASE );
    return NULL;
}



/*
 *	XAUEStateSelect
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEStateSelect( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
Widget		wTank;

    wTank = zwXAUETank( wCur );

    TankSetState((TANK) wTank, TANKSELECT );
    return NULL;
}



/*
 *	XAUEStateSelectDragRotate
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEStateSelectDragRotate( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
Widget		wTank;

    wTank = zwXAUETank( wCur );

    TankSetState((TANK) wTank, TANKDRAGROTATE );
    return NULL;
}



/*
 *	XAUEStateSelectTwist
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEStateSelectTwist( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
Widget		wTank;

    wTank = zwXAUETank( wCur );

    TankSetState((TANK) wTank, TANKTWIST );
    return NULL;
}



/*
 *	XAUESelectRings
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Select all of the rings that have ATOMs selected.
 */
XtCallbackProc
XAUESelectRings( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
LOOP		lInternals, lAtoms;
int		i;
INTERNAL	iInt;
UNIT		uUnit;
VARARRAY	vaRings;
ATOM		aFirst, aAtom;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    uUnit = uTankUnit((TANK)wTank);
    vaRings = vaVarArrayCreate(sizeof(INTERNAL));

	/* Find all of the rings in the UNIT and search through the */
	/* INTERNALs of the atom that the user selected for the INTERNAL */
	/* that represents the smallest ring the ATOM is in */

    GraphUtilFindAllSmallestRings(uUnit);

    lAtoms = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aFirst = (ATOM)oNext(&lAtoms)) ) {
	if ( !bAtomFlagsSet( aFirst, ATOMSELECTED ) ) 
	    continue;
	lInternals = lLoop((OBJEKT) aFirst, INTERNALS );
	while ( (iInt = (INTERNAL)oNext(&lInternals)) ) {
	    if ( iInternalType(iInt) == INTERNALRING ) {
		VarArrayAdd( vaRings, (GENP)&iInt );
	    }
	}
    }

    for ( i=0; i<iVarArrayElementCount(vaRings); i++ ) {
	iInt = *PVAI( vaRings, INTERNAL, i );
	InternalRingLoopAtoms(iInt);
	MESSAGE(( "Ring size: %d\n", iInternalRingSize(iInt) ));
	while ( (aAtom = aInternalRingNextAtom(iInt)) ) {
	    MESSAGE(( "Selecting atom: %s\n", sAtomName(aAtom) ));
	    AtomSetFlags( aAtom, ATOMSELECTED );
	}
    }

    lAtoms = lLoop((OBJEKT) uUnit, ATOMS );
    BuildDestroyInternals(&lAtoms);
    VarArrayDestroy( &vaRings );

    TankRedisplayUnit((TANK)wTank);
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUESelectResidues
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Select all of the RESIDUES that have ATOMs selected.
 */
XtCallbackProc
XAUESelectResidues( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
UNIT		uUnit;
LOOP		lAtoms;
ATOM		aFirst;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    uUnit = uTankUnit((TANK)wTank);
    if ( uUnit == NULL ) return NULL;

		/* First turn off the ATOMTOUCHED flags for all ATOMs */

    ContainerResetAllAtomsFlags((CONTAINER) uUnit, ATOMTOUCHED );

		/* Look over all the ATOMs that are selected and select */
		/* all of the ATOMs within the RESIDUEs that they are */
		/* part of.  Use the ATOMTOUCHED to prevent duplicating */
		/* work. */

    lAtoms = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aFirst = (ATOM)oNext(&lAtoms)) ) {
	if ( bAtomFlagsSet( aFirst, ATOMSELECTED ) &&
	     !bAtomFlagsSet( aFirst, ATOMTOUCHED ) ) {
	    ContainerSetAllAtomsFlags( 
		    cContainerWithin(aFirst), ATOMSELECTED|ATOMTOUCHED );
	}
    }
    TankRedisplayUnit((TANK) wTank );
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
} 


/*
 *	XAUESelectMolecules
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Select all of the MOLECULEs that have ATOMs selected.
 */
XtCallbackProc
XAUESelectMolecules( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
UNIT		uUnit;
LOOP		lAtoms, lSpan;
ATOM		aFirst, aAtom;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );

    uUnit = uTankUnit((TANK)wTank);
    if ( uUnit == NULL ) return NULL;

		/* First turn off the ATOMTOUCHED flags for all ATOMs */

    ContainerResetAllAtomsFlags((CONTAINER) uUnit, ATOMTOUCHED );

		/* Look over all the ATOMs that are selected and select */
		/* all of the ATOMs within the RESIDUEs that they are */
		/* part of.  Use the ATOMTOUCHED to prevent duplicating */
		/* work. */

    lAtoms = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aFirst = (ATOM)oNext(&lAtoms)) ) {
	if ( bAtomFlagsSet( aFirst, ATOMSELECTED ) &&
		!bAtomFlagsSet( aFirst, ATOMTOUCHED ) ) {

	    lSpan = lLoop((OBJEKT) aFirst, SPANNINGTREE );
	    while ( (aAtom = (ATOM)oNext(&lSpan)) ) {
		AtomSetFlags( aAtom, ATOMSELECTED|ATOMTOUCHED );
	    }
	}
    }
    TankRedisplayUnit((TANK) wTank );
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}






/*
 *	XAUEHideSelection
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEHideSelection( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
LOOP		lAtoms;
ATOM		aAtom;
UNIT		uUnit;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    uUnit = uTankUnit((TANK)wTank);

    if ( uUnit == NULL ) return NULL;
    lAtoms = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
	if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
	    AtomSetFlags( aAtom, ATOMNOTDISPLAYED );
	}
    }
    TankRedisplayUnit((TANK) wTank );
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEHideAllButSelection
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEHideAllButSelection( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
LOOP		lAtoms;
ATOM		aAtom;
UNIT		uUnit;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    uUnit = uTankUnit((TANK)wTank);

    if ( uUnit == NULL ) return NULL;
    lAtoms = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
	if ( !bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
	    AtomSetFlags( aAtom, ATOMNOTDISPLAYED );
	}
    }
    TankRedisplayUnit((TANK) wTank );
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}




/*
 *	XAUEUnhideAll
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEUnhideAll( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
UNIT		uUnit;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

    wTank = zwXAUETank( wCur );


    uUnit = uTankUnit((TANK)wTank);
    if ( uUnit == NULL ) 
	return NULL;
    ContainerResetAllAtomsFlags((CONTAINER) uUnit, ATOMNOTDISPLAYED );
    TankRedisplayUnit((TANK) wTank );
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUESelectedAtomsFlipChirality
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUESelectedAtomsFlipChirality( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
UNIT		uUnit;
LOOP		lAtoms;
ATOM		aAtom;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );


    uUnit = uTankUnit((TANK)wTank);
    if ( uUnit == NULL ) return(NULL);

    lAtoms = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
	if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
	    bBuildFlipChiralityFor( uUnit, aAtom );
	}
    }
    TankRedisplayUnit((TANK) wTank );
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}







/*
 *	XAUERelaxSelectionInFramework
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Relax the selected ATOMs within the framework
 *	of the unselected ATOMs.
 */
XtCallbackProc
XAUERelaxSelectionInFramework( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
MINIMIZER       mStrain;
Widget		wTank;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );

    MESSAGE(( "XAUERelaxSelectionInFramework\n" ));

		/* Setup a MINIMIZER and give it a callback to use */
		/* to update the display every step of the minimization */

    mStrain = mMinimizerCreate();
    zXAUESetMinimizer( wCur, mStrain );

		/* Set up the MINIMIZER to use, and turn off any */
		/* control-c that may have been hit before */

    BasicsResetInterrupt();
    MinimizerSetCallback( mStrain, bXAUERelaxCallback, wTank );

    SelectRelaxInFramework( uTankUnit((TANK)wTank), mStrain );

    MinimizerDestroy( &mStrain );
    zXAUESetMinimizer( wCur, NULL );

		/* Redisplay the UNIT */

    TankRedisplayUnit((TANK) wTank );

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}

/*
 *      XAUESensitiveUnit
 *
 *      Author: Vladimir Romanovski (1994)
 *
 *      Make sensetive EditSelectedAtoms menu button for Units
 *         with the same name.
 */

XtCallbackProc
XAUESensitiveUnit( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget	wTop;
  XrmQuark	qTest;
  int		iCount;
  char		*sUnitName;
  UNITEDITORt	*unit;
  register int	i;
  
  /* wCur is shell widget for table */

  wTop = (Widget)WcFullNameToWidget( wCur,"/");
  
  XtVaGetValues( wTop, XtNiconName, &sUnitName, NULL);
  
  qTest = XrmStringToQuark(sUnitName);
  
  iCount = iVarArrayElementCount(SvaEditors);
  
  for ( i = 0; i < iCount ; i++){
      
    unit = PVAI(SvaEditors,UNITEDITORt,i);

    if ( unit->wTop != NULL && unit->qName == qTest )
      {
	XtSetSensitive( (Widget)WcFullNameToWidget( unit->wTop,
						   SELECTION_EDIT), True);
	unit->bTableExist = FALSE;
      }
  }
    return NULL;
}

/*
 *	XAUEEditSelectedAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create a TABLE to edit the selected ATOMs.
 */
XtCallbackProc
XAUEEditSelectedAtoms( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
  Widget	wTank;
  Widget	wTableShell;
  XrmQuark	qTest;
  int		iCount;
  char		*sUnitName;
  UNITEDITORt	*unit;
  register int	i;
  STRING	stTitle;
  char		*s;
  
    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

  wTank = zwXAUETank( wCur );
  
  wTableShell = wXAATPopupTable ( wCur, zwXAUETop(wCur),
					 uTankUnit((TANK)wTank), zXAUEEditDone);

  if ( wTableShell != (Widget)NULL)
    {
      XtVaGetValues( zwXAUETop ( wCur), XtNiconName, &sUnitName, NULL);

      XtVaGetValues( wTableShell, XtNtitle, &s, NULL);

      stTitle[0] = '\0';
      (void)strcat( stTitle, s);
      (void)strcat( stTitle, ": ");
      (void)strcat( stTitle, sUnitName);

      XtVaSetValues(  wTableShell,
		    XtNiconName, (XtArgVal) stTitle,
		    XtNtitle, (XtArgVal) stTitle,
		    NULL);

      qTest = XrmStringToQuark(sUnitName);

      iCount = iVarArrayElementCount(SvaEditors);
      
      for ( i = 0; i < iCount ; i++){
	
	unit = PVAI(SvaEditors,UNITEDITORt,i);

	if (unit->wTop != NULL && unit->qName == qTest )
	  {
	    XtSetSensitive( (Widget)WcFullNameToWidget( unit->wTop,
						       SELECTION_EDIT), False);
	    
	  }
      }

      zXAUETableExist( wCur, TRUE);

    }
  
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEMarkUnBuilt
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEMarkUnBuilt( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
UNIT		uTemp;
LOOP		lAtoms;
Widget		wTank;
ATOM		aAtom;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );
    uTemp = uTankUnit((TANK)wTank);
    lAtoms = lLoop((OBJEKT) uTemp, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
	    AtomSetFlags( aAtom, ATOMPOSITIONDRAWN|ATOMNEEDSBUILD);
	    AtomResetFlags( aAtom, ATOMPOSITIONKNOWN );
	}
    }
    TankRedisplayUnit((TANK)wTank);
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEMarkBuilt
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEMarkBuilt( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
UNIT		uTemp;
LOOP		lAtoms;
Widget		wTank;
ATOM		aAtom;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


    wTank = zwXAUETank( wCur );
    uTemp = uTankUnit((TANK)wTank);
    lAtoms = lLoop((OBJEKT) uTemp, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
	    AtomResetFlags( aAtom, ATOMPOSITIONDRAWN|ATOMNEEDSBUILD );
	    AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
	}
    }
    TankRedisplayUnit((TANK)wTank);
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}





/*
 *	XAUEPrintAllInternals
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	List all of the bonds, angles, and torsions that are on the
 *	selected ATOMs.
 */
XtCallbackProc
XAUEPrintAllInternals( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
UNIT		uUnit;
ATOM		aAtom;
LOOP		lSelected;
LOOP		lTemp;
ATOM		a1, a2, a3, a4;
STRING		s1, s2, s3, s4;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

		/* Define where the output should go */


		/* Get the UNIT being edited and put it */
		/* in an ASSOC to pass it to the command */
		/* handler */

    VP0(( "--- Bonds, angles, torsions of selected ATOMs\n" ));
    uUnit = zuXAUEUnit( wCur );
    lSelected = lLoop((OBJEKT) uUnit, ATOMS );
    LoopDefineVisibleAtoms( &lSelected, ATOMSELECTED );
    while ( (aAtom = (ATOM)oNext(&lSelected)) ) {
        VP0(( "Looking at atom: %s\n", 
			sContainerFullDescriptor((CONTAINER)aAtom,s1) ));

	LOOPOVERALL( aAtom, BONDS, a1, ATOM, lTemp ) {
	    LoopGetBond( &lTemp, &a1, &a2 );
	    VP0(( "BOND: %s - %s\n",
			sContainerFullDescriptor((CONTAINER)a1,s1),
			sContainerFullDescriptor((CONTAINER)a2,s2) ));
	}
	LOOPOVERALL( aAtom, ANGLES, a1, ATOM, lTemp ) {
	    LoopGetAngle( &lTemp, &a1, &a2, &a3 );
	    VP0(( "ANGLE: %s - %s - %s\n",
			sContainerFullDescriptor((CONTAINER)a1,s1),
			sContainerFullDescriptor((CONTAINER)a2,s2),
			sContainerFullDescriptor((CONTAINER)a3,s3) ));
	}
	LOOPOVERALL( aAtom, PROPERS, a1, ATOM, lTemp ) {
	    LoopGetTorsion( &lTemp, &a1, &a2, &a3, &a4 );
	    VP0(( "TORSION: %s - %s - %s - %s\n",
			sContainerFullDescriptor((CONTAINER)a1,s1),
			sContainerFullDescriptor((CONTAINER)a2,s2),
			sContainerFullDescriptor((CONTAINER)a3,s3),
			sContainerFullDescriptor((CONTAINER)a4,s4) ));
	}
	VP0(( "\n" ));
    }
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}





/*
 *	XAUECheck
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Execute the 'check' command on the UNIT and print
 *	the result in the Status window.
 */
XtCallbackProc
XAUECheck( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
ASSOC		aAssoc;
UNIT		uUnit;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

		/* Define where the output should go */


		/* Get the UNIT being edited and put it */
		/* in an ASSOC to pass it to the command */
		/* handler */

    uUnit = zuXAUEUnit( wCur );

    aAssoc = (ASSOC)oCreate(ASSOCid);
    AssocSetObject( aAssoc, uUnit );
    AssocSetName( aAssoc, zsXAUEVariable(wCur) );

    VP0(( "> check %s\n", zsXAUEVariable(wCur) ));

    oCmd_check( 1, &aAssoc );

    Destroy( (OBJEKT *)&aAssoc );

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}




/*
 *	XAUECharge
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Execute the 'charge' command on the UNIT and print
 *	the result in the Status window.
 */
XtCallbackProc
XAUECharge( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
ASSOC		aAssoc;
UNIT		uUnit;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

		/* Define where the output should go */


		/* Get the UNIT being edited and put it */
		/* in an ASSOC to pass it to the command */
		/* handler */

    uUnit = zuXAUEUnit( wCur );

    aAssoc = (ASSOC)oCreate(ASSOCid);
    AssocSetObject( aAssoc, uUnit );

    VP0(( "> charge %s\n", zsXAUEVariable(wCur) ));
    oCmd_charge( 1, &aAssoc );

    Destroy( (OBJEKT *)&aAssoc );

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL; 
}




/*
 *	XAUEImport
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Combine the UNIT in this unit editor with the UNIT
 *	whose variable name is in the widget whose name is in sOtherUnit.
 */
XtCallbackProc
XAUEImport( Widget wCur, char *sWidget, caddr_t PArg )
{
Widget		wTank;
char		*cPVar;
UNIT		uUnit, uOtherUnit, uTempUnit;

		/* Define where the output should go */

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();


		/* Get the UNIT being edited and put it */
		/* in an ASSOC to pass it to the command */
		/* handler */
		/* Get the variable name associated with the UNIT */
		/* from the Widget whos name is in sWidget */

    uUnit = zuXAUEUnit( wCur );
    XACWidgetNameToValue( wCur, sWidget, &cPVar );
    uOtherUnit = (UNIT)oVariable(cPVar);

    if ( uOtherUnit == NULL ) {
	VP0(( "Import: no unit selected.\n" ));
    	DisplayerReleaseUpdates(); 
    	PopCurrentPrintSink();
	return NULL;
    }
    if ( iObjectType(uOtherUnit) != UNITid ) {
	VP0(( "Only UNITs can be imported.\n" ));
    	DisplayerReleaseUpdates(); 
    	PopCurrentPrintSink();
	return NULL;
    }

		/* Make a copy of the other UNIT and join it to the first */

    uTempUnit = (UNIT)oCopy((OBJEKT)uOtherUnit);
    UnitJoin( uUnit, uTempUnit );

    wTank = zwXAUETank( wCur );
    TankRedisplayUnit((TANK) wTank );

    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}



/*
 *	XAUEClose
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEClose( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
Widget		wTop;
UNIT		uUnit;

    PushCurrentPrintSink(ziXAUEPrintSink(wCur));
    DisplayerAccumulateUpdates();

		/* If a MINIMIZER is defined then do not exit */


    if ( zmXAUEMinimizer(wCur) != NULL ) return NULL;

    MESSAGE(( "Discarding changes\n" ));

    wTop = zwXAUETop(wCur);
    uUnit = zuXAUEUnit(wCur);

    bDisplayerRemove( dContainerDisplayer(uUnit),
			XAUEDisplayerUpdate,
			(GENP)wTop );

    zXAUERemoveUnitEditor(wCur);
    XtDestroyWidget( wTop );
    GiUnitEditors--;
    DisplayerReleaseUpdates(); 
    PopCurrentPrintSink();
    return NULL;
}





/*
 *------------------------------------------------------------------
 *
 *	Callbacks to register Widgets.
 *
 */




/*
 *	XAUETankWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUETankWidgetRegister( Widget wCur, caddr_t PAppData, caddr_t PArg )
{

    MESSAGE(( "Registering the TANK Widget\n" ));

		/* Tell the TANK which UNIT it is working on */

    TankUseUnit((TANK) wCur, zuXAUEUnit(wCur) );

                /* To avoid a dependence on ordering in what 
                   tank and status widgets would be created   V.T.R. 1995*/ 

    if (zwXAUEStatus(wCur) != (Widget)NULL)
      TankDefinePrintSink((TANK) wCur, ziXAUEPrintSink(wCur));

		/* Tell the UNITEDITORt which TANK it contains */

    zXAUESetTank( wCur, wCur );

    return NULL;
}


/*
 *	XAUEStatusWidgetRegister
 *
 *	Author:	Christian Schafmeister (1991)
 */
XtCallbackProc
XAUEStatusWidgetRegister( Widget wCur, caddr_t PAppData, caddr_t PArg )
{
int		iSink;


    MESSAGE(( "Registering the Text Widget\n" ));

    zXAUESetStatus( wCur, wCur );

		/* Define the Print Sink for the Unit Editor */

    iSink = iCreatePrintSink( XATPrintStringToWidget, UNITEDITORPREFIX, 
		(GENP)wCur );
    zXAUESetPrintSink( wCur, iSink );

                /* To avoid a dependence on ordering in what 
                   tank and status widgets would be created   V.T.R. 1995*/ 

    if (zwXAUETank(wCur) != (Widget)NULL)
      TankDefinePrintSink((TANK) zwXAUETank(wCur), iSink);

    return NULL;
}




/*
 *	XAUERegisterUnitEditor
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Callback to register the UnitEditor.
 */
XtCallbackProc
XAUERegisterUnitEditor( Widget wShell, caddr_t PAppData, caddr_t PArg )
{
  STRING sTitle;
  
  /* Define the title of the UnitEditor */

  sprintf( sTitle, "XLEaP: Unit editor: %s", SsUnit );
  XtVaSetValues( wShell,
		XtNtitle, (XtArgVal) sTitle,
		XtNiconName, (XtArgVal) SsUnit,
		NULL );
  
  
  /* Create the new UNITEDITORt */
  
  zXAUEAddUnitEditor( wShell );
  zXAUESetUnit( wShell, SuUnit );
  zXAUESetVariable( wShell, SsUnit );
  zXAUESetMinimizer( wShell, NULL );
  
  DisplayerAdd( dContainerDisplayer(SuUnit),
	       XAUEDisplayerUpdate,
	       (GENP)wShell );
  
  GiUnitEditors++;
  return NULL;
}









/*
 *-------------------------------------------
 *-------------------------------------------
 *-------------------------------------------
 *-------------------------------------------
 *-------------------------------------------
 *
 *	Public routines
 *
 */






/*
 *	XAUEPopupUnitEditor
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Popup a UNIT editor.
 */
void
XAUEPopupUnitEditor( Widget wOrig, char *sName, char *sUnit, UNIT uUnit )
{
  strcpy( SsUnit, sUnit );
  SuUnit = uUnit;
  
  MESSAGE(( "Creating UNITEDITOR parent: %s\n", XtName(wOrig) ));
  WcCreateNamedChildren( wOrig, sName );
  
}





/*
 *	XAUEInitialize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Register all of the callback functions with Wcl.
 */
void
XAUEInitialize( XtAppContext app )
{
#define	WCLREG( n, f )	WcRegisterCallback( app, n, (XtCallbackProc)f, NULL )


    SxtacApp = app;

    WCLREG( "XAUEToggleShowNames", XAUEToggleShowNames );
    WCLREG( "XAUEToggleShowPertNames", XAUEToggleShowPertNames );
    WCLREG( "XAUEToggleShowResidues", XAUEToggleShowResidues );
    WCLREG( "XAUEToggleShowTypes", XAUEToggleShowTypes );
    WCLREG( "XAUEToggleShowPertTypes", XAUEToggleShowPertTypes );
    WCLREG( "XAUEToggleShowCharges", XAUEToggleShowCharges );
    WCLREG( "XAUEToggleShowAxis", XAUEToggleShowAxis );
    WCLREG( "XAUEToggleShowBox", XAUEToggleShowBox );
    WCLREG( "XAUESetDrawingElement", XAUESetDrawingElement );
    WCLREG( "XAUEBuildExternals", XAUEBuildExternals );
    WCLREG( "XAUEAddHydrogensAndBuildExternals", 
			XAUEAddHydrogensAndBuildExternals );
    WCLREG( "XAUEStateDraw", XAUEStateDraw );
    WCLREG( "XAUEStateErase", XAUEStateErase );
    WCLREG( "XAUEStateSelect", XAUEStateSelect );
    WCLREG( "XAUESelectResidues", XAUESelectResidues );
    WCLREG( "XAUESelectMolecules", XAUESelectMolecules );
    WCLREG( "XAUEHideSelection", XAUEHideSelection );
    WCLREG( "XAUEHideAllButSelection", XAUEHideAllButSelection );
    WCLREG( "XAUEUnhideAll", XAUEUnhideAll );
    WCLREG( "XAUESelectRings", XAUESelectRings );
    WCLREG( "XAUERelaxSelectionInFramework", 
			XAUERelaxSelectionInFramework );
    WCLREG( "XAUEMarkUnBuilt", XAUEMarkUnBuilt );
    WCLREG( "XAUEMarkBuilt", XAUEMarkBuilt );
    WCLREG( "XAUEStateSelectDragRotate", XAUEStateSelectDragRotate );
    WCLREG( "XAUEStateSelectTwist", XAUEStateSelectTwist );
    WCLREG( "XAUESelectedAtomsFlipChirality",
		XAUESelectedAtomsFlipChirality );
    WCLREG( "XAUEEditSelectedAtoms", XAUEEditSelectedAtoms );
    WCLREG( "XAUETankWidgetRegister", XAUETankWidgetRegister );
    WCLREG( "XAUEStatusWidgetRegister", XAUEStatusWidgetRegister );
    WCLREG( "XAUECheck", XAUECheck );
    WCLREG( "XAUECharge", XAUECharge );
    WCLREG( "XAUEImport", XAUEImport );
    WCLREG( "XAUEClose", XAUEClose );
    WCLREG( "XAUEPrintAllInternals", XAUEPrintAllInternals );
    WCLREG( "XAUERegisterUnitEditor", XAUERegisterUnitEditor );

    WCLREG( "XAUESensitiveUnit",    XAUESensitiveUnit    );
    WCLREG( "XAUECheckUnitEdit",    XAUECheckUnitEdit    );       
}
