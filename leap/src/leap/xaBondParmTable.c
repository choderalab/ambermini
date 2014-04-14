/*
 *	File:		xaBondParmTable.c
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
 *	Author:		David A. Rivkin
 *	Date Created:	(Completed) 17 August 1992
 *	This file is based on xaBondParmTable.c
 *
 *	Description:
 *		Handle editing of parameters in a table format.
 */



#include	<X11/IntrinsicP.h>
#include	<X11/StringDefs.h>


#include	"basics.h"

#include	"varArray.h"

#include	"classes.h"

#include	"xaTable.h"

#include        "../Xraw/table.h"

#define TYPE1C		0
#define	TYPE2C		1
#define	KBC		2
#define	R0C		3
#define	DESCC		4

#define MAXTYPELEN	5
#define DESCLEN		32

typedef	struct	{
		PARMSET		psParmSet;
		VFUNCTION	fDestroyCallback;
		Widget		wTop;
		} BONDPARMTABLEt;

/*
 *----------------------------------------------------------------
 *
 *	Static variables
 *
 */

static	STRING	SsBuffer;
static	STRING	SsError;


/*
 *	zXABPTDestroyTable
 *
 *	Callback that is called when the TABLE is destroyed.
 *	This routine cleans up the extra data that was MALLOC'd
 *	for the BOND properties table.
 */
static void
zXABPTDestroyTable( TABLE tTable )
{
BONDPARMTABLEt	*bptPTemp;

    bptPTemp = (BONDPARMTABLEt*)PXATClientPointer(tTable);
	/* 
	 *  Check to see if there are any rows left.  
	 *	Set the parmSet to having no Bond parameters 
	 */
    if ( !iXATRowsInt( tTable) ) {
    	ParmSetNewBonds( bptPTemp->psParmSet, 0 );
    }
    bptPTemp->fDestroyCallback(bptPTemp->wTop);
        /* Destroy the VARARRAY that holds the table->ParmSet references */
    FREE( bptPTemp );
 
}


/*
 *	zcPXABPTGetElement
 *
 *	Get the values for the elements of the TABLE from the
 *	particular Bond Parmeter Entry.
 */
static char *
zcPXABPTGetElement( TABLE tTable, int iX, int iY )
{
#define	DBLFMT	"%1.4lf"
BONDPARMTABLEt*	bptPCur;
char          	sType1[MAXTYPELEN];
char         	sType2[MAXTYPELEN];
double         	dKb, dR0;
char		sDesc[DESCLEN];


    bptPCur = (BONDPARMTABLEt*)PXATClientPointer(tTable);

    ParmSetBond( bptPCur->psParmSet, iY , sType1, sType2, &dKb, &dR0, sDesc);

    switch ( iX ) {
	    case TYPE1C:
	    	strcpy( SsBuffer, sType1 );
		return( SsBuffer );
		break;
	    case TYPE2C:
	    	strcpy( SsBuffer, sType2 );
		return( SsBuffer );
		break;
	    case KBC:
		sprintf( SsBuffer, DBLFMT, dKb);
		return(SsBuffer);
		break;
	    case R0C:
	    	sprintf( SsBuffer, DBLFMT, dR0);
		return(SsBuffer);
		break;
	    case DESCC:
	    	strcpy( SsBuffer, sDesc );
		return( SsBuffer );
		break;
	    default:
		DFATAL(("Unexpected iX in zcPXABPTGetElement: %d", iX));
    }
    return(NULL);	/* for lint */
}


/*
 *	zcPXABPTVerifyElement
 *
 *	Verify the element to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the error.
 */
static char *
zcPXABPTVerifyElement( TABLE tTable, int iCol, int iRow, char *cPData )
{
double		dValue;

    switch ( iCol ) {
	    case TYPE1C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type1 must have name.");
		} else if ( !isalpha(*cPData) ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		    sprintf( SsError, 
			"Type1 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case TYPE2C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type2 must have name.");
		} else if ( !isalpha(*cPData) ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		    sprintf( SsError, 
			"Type2 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case KBC:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in Kb field.");
		}
		break;
	    case R0C:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in R0 field.");
		}
		break;
	    case DESCC:
		if ( strlen(cPData)>DESCLEN-1 ) {
		    sprintf( SsError, "%s %d characters.",
			"Parameter Description cannot be longer than",
			DESCLEN-1 );
		    return(SsError);
		}
		break;
	    default:
		DFATAL(("Unexpected column in zcPXABPTVerifyElement: %d", 
				iCol ));
    }
    return(NULL);
}


/*
 *	zcPXABPTVerifyRow
 *
 *	Verify the row to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the error.
 */
static char *
zcPXABPTVerifyRow(Widget w, TABLE tTable, int iRow, STRING col[], int *iPErrCol)
{
double		dValue;
char		*cPData;

    	/*
	 *  Type1
	 */
	*iPErrCol = 0;
	cPData = col[0];
	if ( strlen( cPData ) == 0 ) {
        	return("Type1 must have name.");
	} else if ( !isalpha(*cPData) ) {
		return("Type1: first character must be alphabetic.");
	} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		sprintf( SsError, 
			"Type1 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		return(SsError);
	}

	/*
	 *  Type2
	 */
	*iPErrCol = 1;
	cPData = col[1];
	if ( strlen( cPData ) == 0 ) {
	        return("Type2 must have name.");
	} else if ( !isalpha(*cPData) ) {
		return("Type2: first character must be alphabetic.");
	} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		sprintf( SsError, 
			"Type2 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		return(SsError);
	}

	/*
	 *  check/fix/update type order
	 */
	if ( strcmp( col[0], col[1] ) > 0 ) {
		STRING	sTemp;
		SWAP_STRINGS( col[0], col[1], sTemp );
		(void) XawTableSetLabel (w, iRow, 0, col[0] );
		(void) XawTableSetLabel (w, iRow, 1, col[1] );
	}

	/*
	 *  Kb
	 */
	*iPErrCol = 2;
	cPData = col[2];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in Kb field.");
	} else if ( dValue < 0.0 ) {
		return("Kb cannot be negative.");
	}

	/*
	 *  R0
	 */
	*iPErrCol = 3;
	cPData = col[3];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in R0 field.");
	} else if ( dValue < 0.0 ) {
		return("R0 cannot be negative.");
	}

	/*
	 *  Desc
	 */
	*iPErrCol = 4;
	cPData = col[4];
	if ( strlen(cPData) > DESCLEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Parameter Description cannot be longer than",
			DESCLEN-1 );
		return(SsError);
	}

	return(NULL);
}


/*
 *	zXABPTAcceptRow
 * 
 *	Accept/update table. Caller must check row 1st.
 */
static void
zXABPTAcceptRow( TABLE tTable, int iRow, STRING col[])
{
double		dKb, dR0;
PARMSET		psParmSet;
BONDPARMTABLEt	*bptPCur;
  	
    bptPCur = (BONDPARMTABLEt*)PXATClientPointer(tTable);
    psParmSet = bptPCur->psParmSet;

    /*
     *  convert the numeric fields
     */
    bStringToDouble( col[2], &dKb);
    bStringToDouble( col[3], &dR0);

    /*
     *  add or update row
     */
    if ( iRow > iParmSetTotalBondParms( psParmSet ) ) 
	DFATAL(( "programming err 1 in zXABPTAcceptRow\n" ));

    if ( iRow == iParmSetTotalBondParms( psParmSet ) ) {
	/*
	 *  need to add row to table
	 */
    	VP0(( "Adding bond parameter %s-%s.\n", col[0], col[1] ));
        if ( iRow != iParmSetAddBond( psParmSet, col[0], col[1], 
				dKb, dR0, col[4] ) ) 
		DFATAL(( "programming err 2 in zXABPTAcceptRow\n" ));
    } else {
	/*
	 *  update row in place
	 */
	ParmSetUpdateBond( psParmSet, iRow, col[0], col[1], &dKb, &dR0, col[4]);
    }

}




/*
 *------------------------------------------------------------------
 *
 *	Public routines
 *
 */



/*
 *	bXABPTPopupTable (X, Athena, Bond Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XABPTPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
BONDPARMTABLEt	*bptPNew;
int		iCount;


	/* First create a place to put the PARM TABLE data */

    MALLOC( bptPNew, BONDPARMTABLEt*, sizeof(BONDPARMTABLEt) );

	/* Set the variables that need to be passed on */
	
    bptPNew->fDestroyCallback = fCallback;
    bptPNew->wTop = wWidget;
    bptPNew->psParmSet = psParmSet;

	/* Determine how may Bond Parameters there are so that the nubmer of
		rows can be set */
		
    iCount = iParmSetTotalBondParms( psParmSet );
    
	/* Create the table for the BOND PARM TABLE */

    XATPopupTable(
		  wCreated,
		  wWidget,
		  "bParmTableShell",
		  5, iCount,
		  "bParmRow",
		  "bParmAddRow",
		  zcPXABPTGetElement,
		  zcPXABPTVerifyElement,
		  zcPXABPTVerifyRow,
		  zXABPTAcceptRow,
		  NULL,
		  NULL,
		  zXABPTDestroyTable,
		  0,
		  (GENP)bptPNew );
    
}

/*
 *	bXABPTNewPopupTable (X, Athena, Bond Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XABPTNewPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
BONDPARMTABLEt	*bptPNew;


		/* First create a place to put the PARM TABLE data */

    MALLOC( bptPNew, BONDPARMTABLEt*, sizeof(BONDPARMTABLEt) );

	/* Set the variables that need to be passed on */
	
    bptPNew->fDestroyCallback = fCallback;
    bptPNew->wTop = wWidget;
    bptPNew->psParmSet = psParmSet;
    
	/* Create the table for the BOND PARM TABLE */

    XATPopupTable(
		  wCreated,
		  wWidget,
		  "bParmTableShell",
		  5, 1,
		  "bParmAddRow",
		  "bParmAddRow",
		  zcPXABPTGetElement,
		  zcPXABPTVerifyElement,
		  zcPXABPTVerifyRow,
		  zXABPTAcceptRow,
		  NULL,
		  NULL,
		  zXABPTDestroyTable,
		  0,
		  (GENP)bptPNew );

}


