/*
 *	File:		xaHBondParmTable.c
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
 *	Date Created:	11 August 1992
 *	Dates Changed:	
 *	
 *			This file is based on xaAtomParmTable.c
 *
 *	Description:
 *			Handle editing of parameters in a table format.
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
#define	AC		2
#define	BC		3
#define	DESCC		4

#define MAXTYPELEN	5
#define DESCLEN		32

typedef	struct	{
		PARMSET		psParmSet;
		VFUNCTION	fDestroyCallback;
		Widget		wTop;
		} HBONDPARMTABLEt;


/*
 *----------------------------------------------------------------
 *
 *	Static variables
 *
 */

static	STRING	SsBuffer;
static	STRING	SsError;


/*
 *	zXAHPTDestroyTable
 *
 *	Callback that is called when the TABLE is destroyed.
 *	This routine cleans up the extra data that was MALLOC'd
 *	for the HBond properties table.
 */
static void
zXAHPTDestroyTable( TABLE tTable )
{
HBONDPARMTABLEt	*hptPTemp;

    hptPTemp = (HBONDPARMTABLEt*)PXATClientPointer(tTable);
	/* 
	 *  Check to see if there are any rows left.  
	 *	Set the parmSet to having no Hydrogen Bond parameters 
	 */
    if ( !iXATRowsInt(tTable) ) {
    	ParmSetNewHBonds( hptPTemp->psParmSet, 0 );
    }
    hptPTemp->fDestroyCallback(hptPTemp->wTop);
    	/* Destroy the VARARRAY that holds the table->ParmSet references */
    FREE( hptPTemp );
        
}



/*
 *	zcPXAHPTGetElement
 *
 *	Get the values for the elements of the TABLE from the
 *	particular Atom Parmeter Entry.
 */
static char *
zcPXAHPTGetElement( TABLE tTable, int iX, int iY )
{
#define	DBLFMT	"%1.4lf"
HBONDPARMTABLEt		*hptPCur;
PARMSET			psParmSet;
char			sType1[MAXTYPELEN], sType2[MAXTYPELEN];
double			dA, dB;
char			sDesc[DESCLEN];

    hptPCur = (HBONDPARMTABLEt*)PXATClientPointer(tTable);
    psParmSet = hptPCur->psParmSet;

    ParmSetHBond( psParmSet, iY, sType1, sType2, &dA, &dB, sDesc );

    switch ( iX ) {
	    case TYPE1C:
	    	strcpy( SsBuffer , sType1 );
		break;
	    case TYPE2C:
	    	strcpy( SsBuffer , sType2 );
		break;
	    case AC:
		sprintf( SsBuffer, DBLFMT, dA);
		break;
	    case BC:
	    	sprintf( SsBuffer, DBLFMT, dB);
		break;
	    case DESCC:
	    	strcpy( SsBuffer , sDesc );
		break;
	    default:
		DFATAL(("Unexpected iX in zcPXAHPTGetElement: %d", iX));
    }
    return( SsBuffer );
}


/*
 *	zcPXAHPTVerifyElement
 *
 *	Verify the element to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the error.
 */
static char *
zcPXAHPTVerifyElement( TABLE tTable, int iX, int iY, char *cPData )
{
double		dValue;

    switch ( iX ) {
	    case TYPE1C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type1 must have name.");
		} else if ( !isalpha(*cPData) ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		    sprintf( SsError, 
			"Type2 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	   case TYPE2C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type 2 must have name.");
		} else if ( !isalpha(*cPData) ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		    sprintf( SsError, 
			"Type 2 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case AC:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in A field.");
		}
		break;
	    case BC:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in B field.");
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
		DFATAL(("Unexpected iX in zcPXAHPTVerifyElement: %d", iX));
    }
    return(NULL);
}


/*
 *	zcPXAHPTVerifyRow
 *
 *	Verify the row to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the 1st column error.
 */
static char *
zcPXAHPTVerifyRow( Widget w, TABLE tTable, int iRow, STRING *col, int *iPErrCol)
{
char	*cPData;
double	dValue;

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
	 *  A
	 */
	*iPErrCol = 2;
	cPData = col[2];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in A field.");
	} else if ( dValue < 0.0 ) {
		return("A cannot be negative.");
	}

	/*
	 *  B
	 */
	*iPErrCol = 3;
	cPData = col[3];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in B field.");
	} else if ( dValue < 0.0 ) {
		return("B cannot be negative.");
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
 *	zXAHPTAcceptRow
 * 
 */
static void
zXAHPTAcceptRow( TABLE tTable, int iRow, STRING col[] )
{
double			dA, dB;
HBONDPARMTABLEt		*hptPCur;
PARMSET			psParmSet;

	hptPCur = (HBONDPARMTABLEt*)PXATClientPointer(tTable);
	psParmSet = hptPCur->psParmSet;
    
	/*
	 *  convert the numeric fields
	 */
	bStringToDouble( col[2], &dA );
	bStringToDouble( col[3], &dB );

	/*
	 *  add or update row
	 */
	if ( iRow > iParmSetTotalHBondParms( psParmSet ) )
		DFATAL(( "programming err 1 in zXAHPTAcceptRow\n" ));

	if ( iRow == iParmSetTotalHBondParms( psParmSet ) ) {
		/*
		 *  need to add row to table
		 */
		MESSAGE(( "Adding hbond parameter %s-%s.\n", col[0], col[1] ));
        	if ( iRow != iParmSetAddHBond( psParmSet, col[0], col[1],
				dA, dB, col[4] ) ) 
			DFATAL(( "programming err 2 in zXABHPTAcceptRow\n" ));
	} else {
		/*
		 *  update row in place
		 */
		ParmSetUpdateHBond( psParmSet, iRow, col[0], col[1],
				&dA, &dB, col[4] );
	}

}




/*
 *------------------------------------------------------------------
 *
 *	Public routines
 *
 */



/*
 *	XAHPTPopupTable (X, Athena, HydrogenBond Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XAHPTPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
HBONDPARMTABLEt		*hptPNew;
int			iCount;

		/* First create a place to put the PARM TABLE data */

    MALLOC( hptPNew, HBONDPARMTABLEt*, sizeof(HBONDPARMTABLEt) );

    hptPNew->fDestroyCallback = fCallback;
    hptPNew->wTop = wWidget;
    hptPNew->psParmSet = psParmSet;

		/* First create the table for the PARM TABLE */

    iCount = iParmSetTotalHBondParms( psParmSet );
 
    XATPopupTable(
		  wCreated,
		  wWidget,
		  "hParmTableShell",
		  5, iCount,
		  "hParmRow",
		  "hParmAddRow",
		  zcPXAHPTGetElement,
		  zcPXAHPTVerifyElement,
		  zcPXAHPTVerifyRow,
		  zXAHPTAcceptRow,
		  NULL,
		  NULL,
		  zXAHPTDestroyTable,
		  0,
		  (GENP)hptPNew );
    if ( !iCount ) {
	VP0(( "There are no Hydrogen Bond parameters in this PARMSET.\n"));
	VP0(("You may add some using the menu if you wish.\n" ));
    }

}

/*
 *	XAHPTNewPopupTable (X, Athena, HydrogenBond Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XAHPTNewPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
HBONDPARMTABLEt	*hptPNew;
	

		/* First create a place to put the PARM TABLE data */

    MALLOC( hptPNew, HBONDPARMTABLEt*, sizeof(HBONDPARMTABLEt) );

    hptPNew->fDestroyCallback = fCallback;
    hptPNew->wTop = wWidget;
    hptPNew->psParmSet = psParmSet;

		/* First create the table for the PARM TABLE */
 
    XATPopupTable(
		  wCreated,
		  wWidget,
		  "hParmTableShell",
		  5, 1,
		  "hParmAddRow",
		  "hParmAddRow",
		  zcPXAHPTGetElement,
		  zcPXAHPTVerifyElement,
		  zcPXAHPTVerifyRow,
		  zXAHPTAcceptRow,
		  NULL,
		  NULL,
		  zXAHPTDestroyTable,
		  0,
		  (GENP)hptPNew );

}

