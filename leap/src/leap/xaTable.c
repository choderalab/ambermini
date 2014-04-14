/*
 *    File:   xaTable.c
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
 *
 *      Description:
 *              Handle entry of data in a table format.
 *              The bulk of the table is defined in a resource file
 *              making appropriate callbacks into this code to register
 *              the individual ATHENA Text widgets that will be used
 *              to actually enter the data.
 */

/*
 * Davids Changes 28-August-1992 
 */


#define TABLESHELLCLASS      topLevelShellWidgetClass
#define PARMTABLECLASS       transientShellWidgetClass

#include        <stdlib.h>

#include        <X11/IntrinsicP.h>

#ifdef sun
#include <X11/ObjectP.h>     /* why don't they just use X from mit ?!! */
#include <X11/RectObjP.h>
#endif

#include        <X11/StringDefs.h>
#include        <X11/Shell.h>
#include        <X11/cursorfont.h>
#include        <X11/keysym.h>
#include        "../Xraw/ScrolledTable.h"
#include        "../Xraw/Scrollbar.h"
#include        "../Xraw/Table2.h"
#include        "../Wc/WcCreate.h"
#include        "../Xraw/3d.h"

#include        "basics.h"
#include        "varArray.h"

#include        "xaTable.h"
#include        "xaTools.h"
#include        "xaLeap.h"


extern Widget WcFullNameToWidget();
#define UPGRADE_STRING(a,b)                               \
    if ( (a) != (String)NULL ) XtFree(a);                   \
    if ( (b) != (String)NULL ) (a) = XtNewString((String)(b))

#define TABLEPREFIX  "Table Editor: "   /* For now leave blank */
#define TABLE_VP0(a) {tTable->bCouldNotFind = FALSE; VP0(a);}
void VarArrayDeleteMore( VARARRAY, int, int );

typedef struct {
    int row;
    int column;
    String sMessage;
} ErrorStruct, *ErrorStructPtr;


/*
 *----------------------------------------------------------------------
 *
 *      Actions and Translations
 *
 */

static void XATFindReturn ();
static XtActionsRec XATActions[] =
{
    {"find_return", XATFindReturn}
};


/*
 *----------------------------------------------------------------------
 *
 *      Private routines
 *
 */

static XtAppContext SxtacApp = NULL;
static VARARRAY SvaTables = NULL;

/*
 *    zXATAddTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add a new TABLEt to the collection.
 *      The Shell which contains the TABLEt is used to identify
 *      the TABLEt from all the others.
 */

static
TABLEt *
ztXATAddTable(Widget wShell)
{
    int i;
    TABLEt tNew;
    TABLEt *tPNew;

    /*
   * If there are no UnitEditors yet then we ought to create the 
   * VARARRAY to store them.
   */

    if (SvaTables == NULL)
    SvaTables = vaVarArrayCreate (sizeof (TABLEt));

    /*
   * Now try to find an empty place to put it.
   */

    for (i = 0; i < iVarArrayElementCount (SvaTables); i++)
    if (PVAI (SvaTables, TABLEt, i)->wTop == NULL)
        break;


    /*
   * If no empty space was found then add one. 
   */

    if (iVarArrayElementCount (SvaTables) >= i) {
        VarArrayAdd (SvaTables, (GENP)&tNew);
        tPNew = PVAI (SvaTables, TABLEt, iVarArrayElementCount (SvaTables) - 1);
    } else {
        tPNew = PVAI (SvaTables, TABLEt, i);
    }

    /*
   * Set up the new TABLEt.
   */

    tPNew->wTop = wShell;
    return ((TABLE) (tPNew));

}


/*
 *    WidgetToTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a pointer to the TABLEt field in the
 *      collection of TABLEt for the particular
 *      Widget.
 *
 *      The caller can provide ANY Widget that is a sub-Widget
 *      of the ShellWidget that contains the TABLEt.
 *
 *      Also register the TABLEs Status Widget as the current
 *      output Widget.
 */

static TABLEt *
WidgetToTable(Widget wSub)
{
    int i, iMax;
    TABLEt *tPCur;
    Widget wTop;

    /*
   * First find the TOP Widget for the TABLE.
   * Do this by climbing up the Widget hierarchy 
   * until we find a Widget of the proper class.
   */

    for (wTop = wSub; wTop != NULL; wTop = XtParent (wTop)) {
        if (XtClass (wTop) == TABLESHELLCLASS ||
            XtClass (wTop) == PARMTABLECLASS)
            break;
    }
    if (wTop == NULL)
        DFATAL (("Could not find parent TABLESHELLCLASS or PARMTABLECLASS"));

    /*
   * Now find the TABLE that is associated with the TOP Widget
   * OBSOLETE Then register the Status Widget as the Widget 
   * to which output should be sent from the VPx macros
   */

    if ((iMax = iVarArrayElementCount (SvaTables))) {
        tPCur = PVAI (SvaTables, TABLEt, 0);
    for (i = 0; i < iMax; i++, tPCur++)
        if (tPCur->wTop == wTop)
            return (tPCur);
    }

    DFATAL (("Could not find a TABLE"));

    return ((TABLEt*)NULL);

}




/*
 *    zXATRemoveTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Destroy a TABLE referenced by a Widget that is
 *      in its hierarchy.
 */

static void
zXATRemoveTable(TABLEt *tTable)
{

    /*
   * Just set the Top Widget to NULL.  This renders the table 
   * a dummy in the SvaTables vararray. The space for the TABLEt 
   * can be reused later.
   */

    tTable->wTop = NULL;
}

/*
 *    XATRemoveVarArrayOfError 
 *
 *      Author: Vladimir Romanovski (1994)
 *
 *      Destroy an array of errors in table.
 *      Each element in the array has string, so they must be free
 *      separately.
 */

static void
XATRemoveVarArrayOfError(TABLEt *tTable)
{
    register int i;
    int iCount = iVarArrayElementCount (tTable->vaErrors);

    for (i = 0; i < iCount; i++) {
        ErrorStructPtr erError = PVAI (tTable->vaErrors, ErrorStruct, i);

        XawTableSetCellDefaultColours(tTable->wTable, 
                                      erError->row, erError->column);

        if (erError->sMessage)
            XtFree (erError->sMessage);
    }

    tTable->iError = 0;
    VarArraySetSize (tTable->vaErrors, 0);

}


/*
 *    XATDestroy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Destroy the table, freeing any space allocated for it.
 */

static void
XATDestroy(TABLEt *tTable)
{

    if (!IS_GADGET (tTable->wWidgetCreatedTable))
        fprintf(stderr,"XATDestroy: There are problems with WidgetCreatedTable\n");
    else
        XtSetSensitive (tTable->wWidgetCreatedTable, True);

    /*
   * Call the clients Table Destroy callback.
   */

    (tTable->fDestroyTable) (tTable);

    /*
   * Destroy the print sink.
   */

    DestroyPrintSink (tTable->iPrintSink);

    /*
   * Destroy the TABLEt components.
   */

    if (tTable->sFindString != (String) NULL)
        XtFree (tTable->sFindString);

    XATRemoveVarArrayOfError (tTable);
    VarArrayDestroy (&tTable->vaErrors);

    XtDestroyWidget (tTable->wTop);
    zXATRemoveTable (tTable);

}



/*
 *      XATMakeRowVisible
 *
 *      Author: Vladimir Romanovski (1994)
 *
 *      Make iCurent row visible on ScrolledTable.
 *
 */
static void
XATMakeRowVisible(TABLEt *tTable, int iCurent)
{
    float fCurent;
    float fNext;
    float fTop, fShown;

    PushCurrentPrintSink( tTable->iPrintSink );

    fCurent = ((float)iCurent)       / ((float)tTable->iRows);
    fNext   = ((float)(iCurent + 1)) / ((float)tTable->iRows);

    XtVaGetValues (tTable->wScroll,
                 XtNtopOfThumb, &fTop,
                 XtNshown, &fShown,
                 NULL);

    if (fCurent < fTop) {
        XawScrolledTableSetLocation (tTable->wView, 0.0, (double) fCurent);
    } else if (fNext > (fTop + fShown)) {
        XawScrolledTableSetLocation (tTable->wView, 0.0, (double) (fNext - fShown));
    }

    PopCurrentPrintSink();
}


/********************************************************************
 *
 *      XATFindElement
 *
 *      Original : Christian Schafmeister (1991)
 *      Rewrite  : Vladimir Romanovski    (1995) 
 *
 *      Find an element that matches the string given and
 *      then jumps to it.
 *
 ********************************************************************/
static void
XATFindElement(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    char *param = (char *) client_data;
    TABLEt *tTable;
    String sText;
    Widget wAgain;
    Widget wText;
  
    tTable = WidgetToTable (XAGetWidgetFromString (wWid, (char**)&param));

    PushCurrentPrintSink( tTable->iPrintSink );

    wAgain = XAGetWidgetFromString (wWid, (char**)&param);
    wText  = XAGetWidgetFromString (wWid, (char**)&param);
  
    if (tTable->iRows == 0) {
        XtSetSensitive (wAgain, False);
        PopCurrentPrintSink();
        return;
    }

    /* Clear previous failed search process. */
    if (tTable->bCouldNotFind) {
        int i;

        for (i = 0; i < 5; i++)
            XtCallActionProc (tTable->wStatus, "delete-previous-word",
                              NULL, NULL, 0);

        tTable->bCouldNotFind = FALSE;
    }

    /* Get the string from the dialog box. */
  
    XtVaGetValues (wText, XtNstring, &sText, NULL);

    UPGRADE_STRING (tTable->sFindString, sText);

    /* Glance over the table collating cell's labels with the string. */
  
    tTable->find_row    = 0;      /* set original cell to begin a search */
    tTable->find_column = 0;                

    if ( XawTableSearchLabel( tTable->wTable, tTable->sFindString,
                   &tTable->find_row, &tTable->find_column)) {
        /*
         * If the string is find out...
         */

        /* Scroll table to make the cell as visible. */
        XATMakeRowVisible (tTable, tTable->find_row);

        /* Raise edit cell in fetching cell. */
        XawTableSetEdit( tTable->wTable, tTable->find_row, tTable->find_column);
    
        /* Permit ``Find again...'' operation. */
        XtSetSensitive (wAgain, True);

    } else {

        /*
         * If the string is NOT find out...
         */

        /* Ban ``Find again...'' operation. */
        XtSetSensitive (wAgain, False);

        /* Make a message. */
        TABLE_VP0 (("Could not find '%s'.\n", tTable->sFindString));

        /* Make remark about failed search to clear this next find process. */
        tTable->bCouldNotFind = TRUE;
    }

    /* Permit ``Find...'' operation. */
    XtSetSensitive (tTable->wFind, True);

    PopCurrentPrintSink();
}

/********************************************************************
 *
 *      XATFindElementAgain
 *
 *      Original : Christian Schafmeister (1991)
 *      Rewrite  : Vladimir Romanovski    (1995) 
 *
 *      Find again the previous element...
 *
 ********************************************************************/
static void
XATFindElementAgain(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);

    PushCurrentPrintSink( tTable->iPrintSink );

    if (tTable->iRows == 0) {
        PopCurrentPrintSink();
        return;
    }

    if ( tTable->find_column != tTable->iColumns-1){
        tTable->find_column++;
    } else if ( tTable->find_row != tTable->iRows-1){
        tTable->find_column = 0;
        tTable->find_row++;
    }else{
        tTable->find_column = 0;
        tTable->find_row    = 0;
    }

    if ( XawTableSearchLabel( tTable->wTable, tTable->sFindString,
                             &tTable->find_row, &tTable->find_column)) {

        XATMakeRowVisible (tTable, tTable->find_row);

        XawTableSetEdit( tTable->wTable, tTable->find_row, tTable->find_column);
    
    } else
        TABLE_VP0 (("Could not find '%s'.\n", tTable->sFindString));

    PopCurrentPrintSink();
}


/*/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/
  
                                ACTIONS
  
/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/*/
  

static void
XATFindReturn (Widget wWid, XEvent *xPEvent, String saParms[], int *iPParms)
{
    KeySym ksSymbol;

    if (xPEvent->type == KeyPress) {

        ksSymbol = XKeycodeToKeysym (GdDisplay, xPEvent->xkey.keycode, 0);

        if (ksSymbol == XK_Return) {
            Widget wButtonOK = WcFullNameToWidget (wWid, "~*ok");
            Widget wShell = WcFullNameToWidget (wWid, "~");
            char sParam[200];
            int i;

            sParam[0] = '\0';
            for (i = 0; i < (*iPParms); i++) {
                (void) strcat (sParam, saParms[i]);
                (void) strcat (sParam, " ");
            }

            XATFindElement (wButtonOK, (caddr_t)sParam, (caddr_t)NULL);
            XtPopdown (wShell);

        }
    }
}


/*/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/
  
                                CALLBACKS
              
/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/*/
  
/********************************************************************
 *      XATSetWidgetFind 
 *
 *      Author: Vladimir Romanovski (1994)
 *
 *      Mark Menu button ``Find...'' to make it sensitive after
 *      processing find.
 *
 ********************************************************************/
/* ARGSUSED */
XtCallbackProc
XATSetWidgetFind(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);

    tTable->wFind = wWid;

    return NULL;
}

/********************************************************************
 *      XATVerifyChangedCell
 *
 *      Author: Vladimir Romanovski (1995)
 *
 *      A cell was changed so we need to make notice about it.
 *      processing find.
 *
 ********************************************************************/
/* ARGSUSED */
static XtCallbackProc
XATVerifyChangedCell(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);

    tTable->bTableChanged = TRUE;
/* %%% TODO - call VerifyCell */

    return NULL;
}

/********************************************************************
 *      XATVerifyAddRow        really unused in this version
 *
 *      Author: Vladimir Romanovski (1995)
 *
 *      New row was added so we need to update locations of errors
 *      and start cell for ``Find again''.
 *
 ********************************************************************/
/* ARGSUSED */
static XtCallbackProc
XATVerifyAddRow(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    int iCount = iVarArrayElementCount (tTable->vaErrors);
    register int i;
    XawTableCallbackStruct        *srt = (XawTableCallbackStruct*)call_data;
  
    PushCurrentPrintSink( tTable->iPrintSink );

    for (i = 0; i < iCount; i++) {
        ErrorStructPtr erError = PVAI (tTable->vaErrors, ErrorStruct, i);

    if ( erError->row >= srt->row)
        erError->row++;
    }

    if ( tTable->find_row >= srt->row)
        tTable->find_row++;

    PopCurrentPrintSink();

    return NULL;
}

/********************************************************************
 *      XATVerifyDeleteRow
 *
 *      Author: Vladimir Romanovski (1995)
 *
 *      A row was deleted so we need to update locations of errors
 *      and start cell for ``Find again''.
 *
 ********************************************************************/
/* ARGSUSED */
static XtCallbackProc
XATVerifyDeleteRow(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    int iCount = iVarArrayElementCount (tTable->vaErrors);
    XawTableCallbackStruct* srt = (XawTableCallbackStruct*)call_data;
  
    PushCurrentPrintSink( tTable->iPrintSink );

    for (iCount-- ; iCount >= 0; iCount--) {
        ErrorStructPtr erError = PVAI (tTable->vaErrors, ErrorStruct, iCount);

    if ( erError->row > srt->row)
        erError->row--;
    else if ( erError->row == srt->row)
        VarArrayDelete(tTable->vaErrors, iCount);
    }

    if ( tTable->find_row > srt->row)
        tTable->find_row--;
    else if ( tTable->find_row == srt->row)
        tTable->find_column = tTable->find_row = 0;

    PopCurrentPrintSink();

    return NULL;
}

/********************************************************************
 *      XATAddRow
 *
 *      Original : Christian Schafmeister (1991)
 *      Rewrite  : Vladimir Romanovski    (1995)
 *
 *      Add a row to the top of the table.
 *
 ********************************************************************/
/* ARGSUSED */
static XtCallbackProc
XATAddRow(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);

    PushCurrentPrintSink( tTable->iPrintSink );

    XawTableDoLayout (tTable->wTable, False);

    XawTableUnsetEdit (tTable->wTable);
  
    XawTableAppendRow (tTable->wTable);

    XawScrolledTableSetLocation (tTable->wView, 0.0, 1.0);

    XawTableDoLayout (tTable->wTable, True);

    TABLE_VP0 (("Added row at end.\n"));

    tTable->iRows++;
    tTable->bTableChanged = TRUE;

    PopCurrentPrintSink();

    return NULL;
}

/********************************************************************
 *      XATInsertRow
 *
 *      Original : Christian Schafmeister (1991)
 *      Rewrite  : Vladimir Romanovski    (1995)
 *
 *      Insert a row in the table before cell which is editable now.
 *
 ********************************************************************/
/* ARGSUSED */
static XtCallbackProc
XATInsertRow(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    int row, column;
  
    PushCurrentPrintSink( tTable->iPrintSink );

    if ((tTable->iRows == 0) || !XawTableIsEditManaged (tTable->wTable)) {
        TABLE_VP0 (("Insert row: must select a row.\n"));
        PopCurrentPrintSink();
        return NULL;
    }

    XawTableGetEditPosition (tTable->wTable, &row, &column);

    XawTableInsertRow (tTable->wTable, row);

    XATMakeRowVisible (tTable, row);

    TABLE_VP0 (("A new row was inserted to the table.\n"));

    tTable->iRows++;
    tTable->bTableChanged = TRUE;

    PopCurrentPrintSink();

    return NULL;
}

/********************************************************************
 *      XATDeleteRow
 *
 *      Original : Christian Schafmeister (1991)
 *      Rewrite  : Vladimir Romanovski    (1995)
 *
 *      Delete the row where is editable cell.
 *
 ********************************************************************/
/* ARGSUSED */
static XtCallbackProc
XATDeleteRow(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    int row, column;

    PushCurrentPrintSink( tTable->iPrintSink );

    if ((tTable->iRows == 0) || !XawTableIsEditManaged (tTable->wTable)) {
        TABLE_VP0 (("Must select row to be deleted.\n"));
        PopCurrentPrintSink();
        return NULL;
    }

    XawTableGetEditPosition (tTable->wTable, &row, &column);
    XawTableDeleteRow (tTable->wTable, row);

    tTable->iRows--;
    tTable->bTableChanged = TRUE;

    PopCurrentPrintSink();

    return NULL;
}

/********************************************************************
 *      XATGetLabelFromCell
 *
 *      Author: Vladimir Romanovski    (1995)
 *
 *      Get label from table's cell and copy out it in array 
 *
 ********************************************************************/
/* ARGSUSED */
static Boolean
XATGetLabelFromCell( Widget wWid, int iRow, int iCol, XawTableCell cell, 
        caddr_t client_data)
{
    STRING *col = (STRING*)client_data;

    strcpy (col[iCol], XawTableGetLabelByCell(cell));

    return False;
}

/********************************************************************
 *      XATAcceptRow
 *
 *      Authors: Vladimir Romanovski & Bill Ross    (1995)
 *
 ********************************************************************/
/* ARGSUSED */
static Boolean
XATAcceptRow( Widget wWid, int iRow, int iCol, XawTableCell cell, 
        caddr_t client_data)
{
    TABLEt *tTable = (TABLEt*)client_data;
    STRING col[20];
    int iRowReturn;
    int iColReturn;
  
    PushCurrentPrintSink( tTable->iPrintSink );

    if ( iCol != 0 )
        DFATAL(( "XATAcceptRow iCol (%d) expected 0\n", iCol ));
  
    /*
   *  get rest of row 
   */
    (void) XawTableWalk (wWid, (XawTableProc) XATGetLabelFromCell,
                       iRow, iRow, 0, tTable->iColumns - 1,
                       XawTABLE_RIGHT_DOWN,
                       &iRowReturn, &iColReturn, (XtPointer)col);

    (*tTable->fAcceptRow) (tTable, iRow, col);
 
    PopCurrentPrintSink();
 
    return False;
}

/********************************************************************
 *      XATReallySaveTable
 *
 *      Author: Vladimir Romanovski (1995)
 *
 *      Save the Table contents on the disk.
 *
 ********************************************************************/
/* ARGSUSED */
static void
XATReallySaveTable(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    char *param = (char *) client_data;
    TABLEt *tTable;
    int iRowReturn;
    int iColReturn;
  

    if (client_data != (caddr_t)NULL)
        tTable = WidgetToTable (XAGetWidgetFromString (wWid, (char**)&param));
    else
        tTable = WidgetToTable (wWid);

    PushCurrentPrintSink( tTable->iPrintSink );

    XawTableUnsetEdit(tTable->wTable);

    if (tTable->bTableChanged) {

        /*
         * If some changes was made...
         */

        if (tTable->fBeginAcceptTable)
            tTable->fBeginAcceptTable (tTable);

        /* Glance over the Table with accepting cell contents. */
        (void) XawTableWalk (tTable->wTable, (XawTableProc)XATAcceptRow,
                             0, tTable->iRows - 1, 0, 0,
                             XawTABLE_DOWN_RIGHT,
                             &iRowReturn, &iColReturn,(caddr_t)tTable);

        if (tTable->fEndAcceptTable)
            (*tTable->fEndAcceptTable) (tTable);

        tTable->bTableChanged = FALSE;
    }

    PopCurrentPrintSink();
}


/********************************************************************
 *      XATVerifyRow
 *
 *      Authors: Vladimir Romanovski & Bill Ross    (1995)
 *
 ********************************************************************/
/* ARGSUSED */
static Boolean
XATVerifyRow( Widget wWid, int iRow, int iCol, XawTableCell cell, 
        caddr_t client_data)
{
    TABLEt  *tTable = (TABLEt*)client_data;
    STRING  col[20];
    char    *cPError;
    Pixel   er_fore = 0;
    Pixel   er_back = 1;
    int     iRowReturn;
    int     iColReturn;
    int     iErrCol;
  
    PushCurrentPrintSink( tTable->iPrintSink );

    if ( iCol != 0 )
        DFATAL(( "XATAcceptRow iCol (%d) expected 0\n", iCol ));

    /*
   *  get columns
   */
    (void) XawTableWalk (wWid, (XawTableProc) XATGetLabelFromCell,
                         iRow, iRow, 0, tTable->iColumns - 1,
                         XawTABLE_RIGHT_DOWN,
                         &iRowReturn, &iColReturn, (XtPointer)col);
  
    /*
   *  check row w/ table-specific call
   */
    cPError = (tTable->fVerifyRow) (wWid, tTable, iRow, col, &iErrCol );
  
    /*
   *  handle err (in 1st column found; others may exist)
   */

    if (cPError != (char*)NULL) {
        ErrorStruct error;
    
        error.sMessage = XtNewString (cPError);
        error.row = iRow;
        error.column = iErrCol;
        VarArrayAdd (tTable->vaErrors, (GENP)&error);
        (void) FetchPixel (wWid, "yellow", &er_fore);
        (void) FetchPixel (wWid, "brown1", &er_back);
/*  fprintf(stderr, "SetLabelColours r/c %d %d\n", iRow,iCol);  */
        XawTableSetCellColours (wWid, iRow, iErrCol, er_fore, er_back);
    }

    PopCurrentPrintSink();

    return False;
}

/********************************************************************
 *      XATVerifyTable
 *
 *      Author: Vladimir Romanovski    (1995)
 *
 *      The base procedure for checking the table.
 *      Return a count of error.
 *
 ********************************************************************/
static int
XATVerifyTable(TABLEt *tTable) 
{
    int iErrorNumber = 0;
    int iX, iY;

    PushCurrentPrintSink( tTable->iPrintSink );

    /*
   * Release the current active cell.
   */

    XawTableUnsetEdit(tTable->wTable);

    /*
   * Clear varArray which keeps previous errors and reset warning counter.
   */

    XATRemoveVarArrayOfError (tTable);
    tTable->iWarningCount = 0;

    /*
   * First validate the data in the TABLE.
   * Call the validation callback for each row.
   */

    TABLE_VP0 (("Check table...\n"));

    iX = iY = 0;
    (void) XawTableWalk (tTable->wTable, (XawTableProc)XATVerifyRow,
                         0, tTable->iRows - 1, 0, 0,
                         XawTABLE_DOWN_RIGHT,
                         &iX, &iY, (caddr_t)tTable);
  
    iErrorNumber = iVarArrayElementCount (tTable->vaErrors);
   
    if (iErrorNumber == 0) {

        TABLE_VP0 (("The table has no errors and %d warnings.\n",
                    tTable->iWarningCount));

    } else {

        ErrorStructPtr erError;

        /*
         * Make the first none correct cell as current active cell.
         */

        erError = PVAI (tTable->vaErrors, ErrorStruct, 0);

        if (iErrorNumber == 1) {
            TABLE_VP0 (("The table has one error and %d warnings.\n",
                        tTable->iWarningCount));
        } else {
            TABLE_VP0 (("The table has %d errors and %d warnings.\n",
                        iErrorNumber, tTable->iWarningCount));
        }

        TABLE_VP0 (("> %s\n", erError->sMessage));

        XawTableSetEdit(tTable->wTable, erError->row, erError->column);

        XATMakeRowVisible (tTable, erError->row);
    }

    PopCurrentPrintSink();

    return (iErrorNumber);
}



/*
 *    XATSaveTable 
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Accept the changes made.
 *      This means cycling through all the Widgets
 *      getting the values of the Widgets and feeding
 *      them back to the caller through the AcceptRowCallback.
 */

static void
XATSaveTable(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    Widget widget = (Widget) NULL;
    STRING widget_name;

    PushCurrentPrintSink( tTable->iPrintSink );

    if (client_data != (caddr_t)NULL) {
        (void) WcCleanName ((char *) client_data, widget_name);
        widget = WcFullNameToWidget (wWid, widget_name);
    }
    if (XATVerifyTable (tTable) == 0) {

        XATReallySaveTable (wWid, (caddr_t)NULL, (caddr_t)NULL);

        TABLE_VP0 (("Table has been 'saved' back to the program.\n"));

        if (widget != (Widget) NULL)
            XtSetSensitive (widget, False);

    } else {

        if (widget != (Widget) NULL)
            XtSetSensitive (widget, True);

        widget = WcFullNameToWidget (wWid, "/*verify_save_none_corect_table");

        XtPopup (widget, XtGrabNonexclusive);

    }

    PopCurrentPrintSink();
}

static void XATReject ();

/*
 *    XATSaveQuitTable 
 *
 *      Author: 
 *
 *      Accept the changes made.
 *      This means cycling through all the Widgets
 *      getting the values of the Widgets and feeding
 *      them back to the caller through the AcceptElementCallback.
 */

static void
XATSaveQuitTable(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    Widget widget = (Widget) NULL;
    STRING widget_name;

    PushCurrentPrintSink( tTable->iPrintSink );

    if (client_data != (caddr_t)NULL) {
        (void) WcCleanName ((char *) client_data, widget_name);
        widget = WcFullNameToWidget (wWid, widget_name);
    }
    if (XATVerifyTable (tTable) == 0) {

        XATReallySaveTable (wWid, (caddr_t)NULL, (caddr_t)NULL);

        TABLE_VP0 (("Table has been 'saved' back to the program.\n"));

        if (widget != (Widget) NULL)
          XtSetSensitive (widget, False);

        XATReject (wWid, (caddr_t)NULL, (caddr_t)NULL);

    } else {

        if (widget != (Widget) NULL)
          XtSetSensitive (widget, True);

        widget = WcFullNameToWidget (wWid, 
                        "/*verify_save_none_corect_table_and_quit");

        XtPopup (widget, XtGrabNonexclusive);

    }

    PopCurrentPrintSink();
}

/*
 *    XATCloseTable 
 *
 *      Author: Vladimir Romanovski (1994)
 *
 *      
 *
 *
 *
 */
static void
XATCloseTable(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    Widget wDialog;

    PushCurrentPrintSink( tTable->iPrintSink );

    /*
   * Release the current active cell.
   */

    XawTableUnsetEdit(tTable->wTable);

    if (tTable->bTableChanged) {
    
        wDialog = WcFullNameToWidget (wWid, "/*verify_close_modified");

        XtPopup (wDialog, XtGrabNonexclusive);

    } else {
        XATReject (wWid, (caddr_t)NULL, (caddr_t)NULL);
    }

    PopCurrentPrintSink();
}



/*
 *    XATCheckTable 
 *
 *      Author: Vladimir Romanovski (1994)
 *
 *      Check for errors in the Table.
 *
 */
static void
XATCheckTable(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    Widget widget = (Widget) NULL;
    STRING widget_name;

    PushCurrentPrintSink( tTable->iPrintSink );

    if (client_data != (caddr_t)NULL) {
        (void) WcCleanName ((char *) client_data, widget_name);
        widget = WcFullNameToWidget (wWid, widget_name);
    }
    if (XATVerifyTable (tTable) == 0) {

    if (widget != (Widget) NULL)
      XtSetSensitive (widget, False);

    } else {

        if (widget != (Widget) NULL)
          XtSetSensitive (widget, True);
    }

    PopCurrentPrintSink();
}

/*
 *    XATGoToNextError 
 *
 *      Author: Vladimir Romanovski (1994)
 *
 *      Go to the next error in the Table.
 *
 */
static void
XATGoToNextError (Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable = WidgetToTable (wWid);
    ErrorStructPtr erError;
    int iErrorCount;

    PushCurrentPrintSink( tTable->iPrintSink );

    iErrorCount = iVarArrayElementCount (tTable->vaErrors);

    if ( iErrorCount > 1) {
        XawTableUnsetEdit(tTable->wTable);
    
        tTable->iError = ++tTable->iError % iErrorCount;
    
        erError = PVAI (tTable->vaErrors, ErrorStruct, tTable->iError);
    
        TABLE_VP0 ((" %s\n", erError->sMessage));
    
        XawTableSetEdit(tTable->wTable, erError->row, erError->column);
    
        XATMakeRowVisible (tTable, erError->row);
    }
  
    PopCurrentPrintSink();
}


/*
 *    XATReject
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Reject the changes made.
 */

static void
XATReject(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable;
    char *param = (char *) client_data;

    if (client_data != (caddr_t)NULL)
        tTable = WidgetToTable (XAGetWidgetFromString (wWid, (char**)&param));
    else
        tTable = WidgetToTable (wWid);

    XATDestroy (tTable);

}

static void
XATSetColumnsWidth(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    char *param = (char*)client_data;
    int i;
    int width;
    STRING label;

    if ( (param == (char*)NULL) || ((*param) == '\0') )
        return;
  
    i = 0;
    do{
        param = WcCleanName( param, label);
        width = atoi (label);
        (void) XawTableSetColumnWidth (wWid, i++, width);
    }while ( (param != (char*)NULL) && ((*param) != '\0') );

}

/*
 *-----------------------------------------------
 *
 *      Register Widgets
 *
 */

/*
 *    XATRegisterStatusWidget
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Register the AsciiText widget used to display status messages.
 */

static void
XATRegisterStatusWidget(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    TABLEt *tTable;

    tTable = WidgetToTable (wWid);

    tTable->wStatus = wWid;
    tTable->iPrintSink =
    iCreatePrintSink (XATPrintStringToWidget, TABLEPREFIX, (GENP)wWid);

}

static void
XATRegisterTitle (Widget wWid, caddr_t client_data, caddr_t call_data)
{
    char *param = (char*)client_data;
    int i;
    STRING label;
    char *c;

    if ( (param == (char*)NULL) || ((*param) == '\0') )
        return;
  
    i = 0;
    do {
        param = WcCleanName( param, label);

        for (c = label; (c != (char*)NULL) && ((*c) != '\0'); c++)
            if (*c == '_') *c = ' ';

        (void) XawTableSetLabel (wWid, 0, i++, label);
    }while ( (param != (char*)NULL) && ((*param) != '\0') );

}

/*
 *    XATRegisterTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Construct the Row Container Widget using the Widget
 *      name passed by the Constructor.
 *      Then create all of the children of the row container Widget
 *      which are the row Widgets.
 */

static void
XATRegisterTable(Widget wWid, caddr_t client_data, caddr_t call_data)
{
    int i, j;
    TABLEt *tTable;
    Widget wView;
  
    for (wView = wWid;
         XtClass (wView) != scrolledTableWidgetClass;
         wView = XtParent (wView))
      /* EMPTY */;

    tTable          = WidgetToTable (wWid);
    tTable->wTable  = wWid;
    tTable->wView   = wView;
    tTable->wScroll = WcFullNameToWidget(wView, "vertical");
  
    MESSAGE (("Creating Container for Table: %5d rows\n", tTable->iRows));

    XtVaSetValues (wWid,
                   XtNcolumns, (XtArgVal) tTable->iColumns,
                   XtNrows, (XtArgVal) tTable->iRows,
                   NULL);

    for (i = 0; i < tTable->iRows; i++) {
        for (j = 0; j < tTable->iColumns; j++) {
            (void) XawTableSetLabel (wWid, i, j, 
                                     (*tTable->fGetElement) (tTable, j, i));
        }
    }

    tTable->bTableChanged = FALSE;
  
}



static Widget SwWidgetCreatedTable;     /* Vladimir Romanovski */
static STRING SsTopWidgetName;
static STRING SsRowWidgetName;
static STRING SsAddedRowWidgetName;
static int SiColumns;
static int SiRows;
static SFUNCTION SfGetElement;
static SFUNCTION SfVerifyElement;
static SFUNCTION SfVerifyRow;
static VFUNCTION SfAcceptRow;
static VFUNCTION SfBeginAcceptTable;
static VFUNCTION SfEndAcceptTable;
static VFUNCTION SfDestroyTable;
static int SiClientInt;
static GENP SPClientPointer;



/*
 *    XATRegisterShell
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Callback to register the Table.
 */

static void
XATRegisterShell(Widget wShell, caddr_t PAppData, caddr_t PArg)
{
    TABLEt *tNew;

    tNew = ztXATAddTable (wShell);

    strcpy (tNew->sRowWidget, SsRowWidgetName);
    strcpy (tNew->sAddedRowWidget, SsAddedRowWidgetName);

    tNew->vaErrors = vaVarArrayCreate (sizeof (ErrorStruct));
    tNew->wWidgetCreatedTable = SwWidgetCreatedTable;
    tNew->fGetElement = SfGetElement;
    tNew->fVerifyElement = SfVerifyElement;
    tNew->fVerifyRow = SfVerifyRow;
    tNew->fAcceptRow = SfAcceptRow;
    tNew->fDestroyTable = SfDestroyTable;
    tNew->fBeginAcceptTable = SfBeginAcceptTable;
    tNew->fEndAcceptTable = SfEndAcceptTable;
    tNew->iColumns = SiColumns;
    tNew->iRows = SiRows;
    tNew->wStatus = (Widget) NULL;
    tNew->wFind = (Widget) NULL;
    tNew->sFindString = (String) NULL;
    tNew->iError = 0;
    tNew->iWarningCount = 0;
    tNew->bCouldNotFind = FALSE;
    tNew->bTableChanged = FALSE;
    tNew->sSaveBeforeEdit = (String) NULL;
    tNew->find_row = 0;
    tNew->find_column = 0;
    tNew->wTable = (Widget) NULL;
  
    XATSetClientInt (tNew, SiClientInt);
    XATSetClientPointer (tNew, SPClientPointer);

    if (!IS_GADGET (tNew->wWidgetCreatedTable))
      fprintf(stderr,"XATRegisterShell: There are problems with WidgetCreatedTable\n");

}


/*******************************************************************

                           Public routines
 
 *******************************************************************/

/*
 *    XATPopupTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a Table.
 *      The name of the Shell Widget that contains everything
 *      is passed in sTopWidgetName.
 *      The initial number of columns and rows is in iColumns and
 *      iRows.
 */

void
XATPopupTable( Widget wWidgetCreatedTable, Widget wOrig,
               char *sTopWidgetName, int iColumns, int iRows,
               char *sRowWidgetName, char *sAddedRowWidgetName,
               SFUNCTION fGetElement, SFUNCTION fVerifyElement,
               SFUNCTION fVerifyRow, VFUNCTION fAcceptRow,
               VFUNCTION fBeginAcceptTable, VFUNCTION fEndAcceptTable,
               VFUNCTION fDestroyTable,
               int iClientInt, GENP PClientPointer)
{


/*
 * Copy the arguments into static variables 
 * so that when the wcCallback is called 
 * everything can be initialized 
 */


    strcpy (SsTopWidgetName, sTopWidgetName);
    strcpy (SsRowWidgetName, sRowWidgetName);
    strcpy (SsAddedRowWidgetName, sAddedRowWidgetName);


    SwWidgetCreatedTable = wWidgetCreatedTable;
    SiColumns = iColumns;
    SiRows = iRows;
    SfGetElement = fGetElement;
    SfVerifyElement = fVerifyElement;
    SfVerifyRow = fVerifyRow;
    SfAcceptRow = fAcceptRow;
    SfBeginAcceptTable = fBeginAcceptTable;
    SfEndAcceptTable = fEndAcceptTable;
    SfDestroyTable = fDestroyTable;
    SiClientInt = iClientInt;
    SPClientPointer = PClientPointer;


    WcCreateNamedChildren (wOrig, sTopWidgetName);

}

void
XATSetClientInt(TABLEt *t, int i)
{
    t->iTableInt = i;
}


int
iXATClientInt(TABLEt *t)
{
    return t->iTableInt;
}


int
iXATRowsInt(TABLEt *t)
{
    return t->iRows;
}

void
XATSetClientPointer(TABLEt *t, GENP p)
{
    t->PTablePointer = p;
}

GENP
PXATClientPointer(TABLEt *t)
{
    return t->PTablePointer;
}

/*
 *    XATInitialize 
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Register all of the callback functions.
 */

void
XATInitialize(XtAppContext xtacApp, Widget topShell)
{
#define WCLREG(name,func) WcRegisterCallback( xtacApp, name, (XtCallbackProc)func, NULL)

    SxtacApp = xtacApp;

    WCLREG ("XATSetWidgetFind", XATSetWidgetFind);
    WCLREG ("XATReallySaveTable", XATReallySaveTable);
    WCLREG ("XATCheckTable", XATCheckTable);
    WCLREG ("XATCloseTable", XATCloseTable);
    WCLREG ("XATSaveQuitTable", XATSaveQuitTable);
    WCLREG ("XATGoToNextError", XATGoToNextError);
    WCLREG ("XATVerifyAddRow", XATVerifyAddRow);
    WCLREG ("XATVerifyDeleteRow", XATVerifyDeleteRow);
    WCLREG ("XATVerifyChangedCell", XATVerifyChangedCell);
    WCLREG ("XATRegisterTitle", XATRegisterTitle);

    /*
   *Davids Changes
   */
    WCLREG ("XATAddRow", XATAddRow);
    WCLREG ("XATInsertRow", XATInsertRow);
    WCLREG ("XATDeleteRow", XATDeleteRow);
    WCLREG ("XATFindElement", XATFindElement);
    WCLREG ("XATFindElementAgain", XATFindElementAgain);

    /*
   * End of Davids Changes
   */

    WCLREG ("XATSaveTable", XATSaveTable);
    WCLREG ("XATReject", XATReject);
    WCLREG ("XATRegisterStatusWidget", XATRegisterStatusWidget);
    WCLREG ("XATRegisterTable", XATRegisterTable);
    WCLREG ("XATRegisterShell", XATRegisterShell);

    WCLREG ("XATSetColumnsWidth", XATSetColumnsWidth);

    XtAppAddActions (xtacApp, XATActions, XtNumber (XATActions));

}

