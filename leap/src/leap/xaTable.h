/*
 *      File:   xaTable.h
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
 *      Description:
 *              This file contains routines to interface with WCL
 *              to create a Table which the user can use to edit
 *              and input large quantities of data.
 *
 *              The caller is able to store an integer and a pointer
 *              to some arbitrary data in the client fields of
 *              the TABLE. 
 *              
 */


#ifndef XATABLE_H
#define XATABLE_H

#include "varArray.h"
#include "xaLeap.h"
#include "varArray.h"


typedef char*   (*SFUNCTION)();
typedef struct  {

  /* Private fields */
  
  String         sFindString;         /* Vladimir Romanovski */
  VARARRAY       vaErrors;            /* Vladimir Romanovski */
  int            iError;              /* Vladimir Romanovski */
  int            iWarningCount;       /* Scott Brozell */
                 /* Count of Warnings from Table Editor's Check Table */

  Widget         wWidgetCreatedTable; /* Vladimir Romanovski */
  Widget         wFind;               /* Vladimir Romanovski */
  BOOL           bCouldNotFind;       /* Vladimir Romanovski */
  BOOL           bTableChanged;       /* Vladimir Romanovski */
  String         sSaveBeforeEdit;     /* Vladimir Romanovski */
  Widget         wView;               /* Vladimir Romanovski */
  Widget         wScroll;             /* Vladimir Romanovski */
  int            find_row;
  int            find_column;
  Widget         wTop;
  Widget         wTable;              /* Vladimir Romanovski */
  Widget         wTitle;              /* Vladimir Romanovski */
  int            iPrintSink;
  STRING         sRowWidget;
  STRING         sAddedRowWidget;     /* Davids Changes */
  int            iColumns;
  int            iRows;
  BOOL           bClearInputOnFirstCharacter;
  Widget         wStatus;
  SFUNCTION      fGetElement;
  SFUNCTION      fVerifyElement;
  SFUNCTION      fVerifyRow;
  VFUNCTION      fAcceptRow;
  VFUNCTION      fBeginAcceptTable;
  VFUNCTION      fEndAcceptTable;
  VFUNCTION      fDestroyTable;
  
  /* Public fields */
  
  int            iTableInt;
  GENP           PTablePointer;
  
} TABLEt, *TABLE;

/*
 *      Functions
 */

extern void XATSetClientInt(TABLEt *  , int );
extern int  iXATClientInt(TABLEt * );
extern int  iXATRowsInt(TABLEt * );
extern void XATSetClientPointer(TABLEt *  , GENP );
extern GENP PXATClientPointer(TABLEt * );
extern void XATInitialize(XtAppContext  , Widget );
extern void XATPopupTable(Widget  , Widget  , char *  , int  , int  ,
                                char *  , char *  , 
                                SFUNCTION  , SFUNCTION  , SFUNCTION  ,
                                VFUNCTION  , VFUNCTION  , VFUNCTION  ,
                                VFUNCTION  , int  , GENP );


#endif /* XATABLE_H */


