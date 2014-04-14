/*
 *	File:		xaAtomParmTable.c
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
 *	This file is based on xaAtomTable.c
 *
 *	Description:
 *		Handle editing of parameters in a table format.
 */


#include	<X11/IntrinsicP.h>
#include	<X11/StringDefs.h>

#include	"basics.h"

#include	"varArray.h"

#include	"classes.h"
#include 	"amber.h"

#include	"elements.h"
#include	"hybrid.h"

#include	"xaTable.h"


#define TYPEC		0
#define	MASSC		1
#define POLARC		2
#define	EPSILONC	3
#define	RC		4
#define	ELEMENTC	5
#define	HYBRIDC		6
#define	DESCC		7


#define MAXTYPELEN	5
#define	DESCLEN		32

typedef	struct	{
		PARMSET		psParmSet;
		VFUNCTION	fDestroyCallback;
		Widget		wTop;
		} ATOMPARMTABLEt;


/*
 *----------------------------------------------------------------
 *
 *	Static variables
 *
 */

static	STRING	SsBuffer;
static	STRING	SsError;


/*
 *	zXAAPTDestroyTable
 *
 *	Callback that is called when the TABLE is destroyed.
 *	This routine cleans up the extra data that was MALLOC'd
 *	for the ATOM properties table.
 */
static void
zXAAPTDestroyTable( TABLE tTable )
{
ATOMPARMTABLEt	*aptPTemp;


    aptPTemp = (ATOMPARMTABLEt*)PXATClientPointer(tTable);
	/* Check to see if there are any rows left.  Set the parmSet to having
		no Atom parameters */
/*
    if ( !iXATRowsInt(tTable) ) {
    	ParmSetNewAtoms( aptPTemp->psParmSet, 0 );
    }
*/
    aptPTemp->fDestroyCallback(aptPTemp->wTop);
        /* Destroy the VARARRAY that holds the table->ParmSet references */
    FREE( aptPTemp );
        
}


/*
 *	zcPXAAPTGetElement
 *
 *	Get the values for the elements of the TABLE from the
 *	particular Atom Parameter Entry.
 */
 
static char *
zcPXAAPTGetElement( TABLE tTable, int iCol, int iRow )
{
#define	DBLFMT	"%1.4lf"
ATOMPARMTABLEt	*aptPCur;
PARMSET		psParmSet;
char		sType[MAXTYPELEN], sDesc[DESCLEN];
double		dMass, dPolar, dEpsilon, dR, dEpsilon14, dR14;
double		dScreenF;
int		iElement, iHybrid;


    aptPCur = (ATOMPARMTABLEt*)PXATClientPointer(tTable);
    psParmSet = aptPCur->psParmSet;
    
    ParmSetAtom( psParmSet, iRow, sType, &dMass, &dPolar, &dEpsilon, &dR, 
		 &dEpsilon14, &dR14, &dScreenF, &iElement, &iHybrid, sDesc );

    switch ( iCol ) {
	case TYPEC:
		strcpy( SsBuffer, sType );
		return( SsBuffer );
		break;
	case MASSC:
	    	sprintf( SsBuffer, DBLFMT, dMass );
		return( SsBuffer );
		break;
	case POLARC:
		if ( dPolar == -1.0 )
			sprintf( SsBuffer, "?" );
		else
			sprintf( SsBuffer, DBLFMT, dPolar );
		return( SsBuffer );
		break;
	case EPSILONC:
		sprintf( SsBuffer, DBLFMT, dEpsilon );
		return( SsBuffer );
		break;
	case RC:
	    	sprintf( SsBuffer, DBLFMT, dR );
		return( SsBuffer );
		break;
	case ELEMENTC:
		sElementName( iElement, SsBuffer );
		return( SsBuffer );
		break;
	case HYBRIDC:
		sHybridName( iHybrid, SsBuffer );
		return( SsBuffer );
		break;
	case DESCC:
		strcpy( SsBuffer, sDesc );
		return( SsBuffer );
		break;
	default:
		DFATAL(("Unexpected column in zcPXAAPTGetElement: %d", iCol ));
    }
    return( NULL );	/* for lint */
}


/*
 *	zcPXAAPTVerifyElement
 *
 *	Verify the element to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the error.
 */

static char *
zcPXAAPTVerifyElement( TABLE tTable, int iCol, int iRow, char *cPData )
{
double		dValue;
int		iLen;

    switch ( iCol ) {
	    case TYPEC:
		iLen = strlen( cPData );
		if ( iLen == 0 )
		    return("Atom Type must have name.");
		if ( !isalpha( *cPData ) )
		    return("First character must be alphabetic.");
		if ( iLen > MAXTYPELEN-1 ) {
		    sprintf( SsError, 
			"Atom Type name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case MASSC:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in Mass field.");
		}
		break;
	    case POLARC:
		if ( strcmp( cPData, "?" ) == 0 )
			break;
		if ( bStringToDouble( cPData, &dValue ) )
			break;
		return("Invalid character in Polar field.");
		break;
	    case EPSILONC:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in Epsilon field.");
		}
		break;
	    case RC:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in R field.");
		}
		break;
	    case ELEMENTC:
		if ( strcmp( cPData, "?" ) && !isalpha(*cPData) ) {
		    sprintf( SsError, "%s %s",
		        "Element names (Atomic Symbol) must begin with",
			"an alphabetic character.");
		    return(SsError);
		} else if ( strlen(cPData)>ELEMENTNAMELEN-1 ) {
		    sprintf( SsError, "%s %d characters.",
			"Element name (Atomic Symbol) cannot be longer than",
			ELEMENTNAMELEN-1 );
		    return(SsError);
		}
		break;
	    case HYBRIDC:
		if ( strlen(cPData)>HYBRIDNAMELEN-1 ) {
		    sprintf( SsError, "%s %d characters.",
			"Atom Hybridization type cannot be longer than",
			HYBRIDNAMELEN-1 );
		    return(SsError);
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
		DFATAL(("unexpected column in zcPXAAPTVerifyElement: %d", 
			iCol));
    }
    return(NULL);	/* for lint */
}


/*
 *	zcPXAAPTVerifyRow
 *
 *	Verify the row to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the 1st column error encountered.
 */

static char *
zcPXAAPTVerifyRow( Widget w, TABLE tTable, int iRow, STRING *col, int *iPErrCol)
{
char		*cPData;
double		dValue;

	/*
	 *  Type
	 */
	*iPErrCol = TYPEC;
	cPData = col[TYPEC];
	if ( strlen( cPData ) == 0 )
		return("Type must have name.");
	if ( !isalpha(*cPData) )
		return("Type: first character must be alphabetic.");
	if ( strlen(cPData) > MAXTYPELEN-1 ) {
		sprintf( SsError, 
			"Type: name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		return(SsError);
	}

	/*
	 *  Mass
	 */
	*iPErrCol = MASSC;
	cPData = col[MASSC];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in Mass field.");
	} 
	if ( dValue < 0.0 ) {
		return("Mass cannot be negative.");
	}

	/*
	 *  Polarizability
	 */
	*iPErrCol = POLARC;
	cPData = col[POLARC];
	if ( strcmp( cPData, "?" ) != 0 ) {
		if ( !bStringToDouble( cPData, &dValue ) ) {
			return("Invalid character in Polariz field.");
		}
		if ( dValue < 0.0 ) {
			return("Polarizability cannot be negative.");
		}
		if ( dValue > 15.0 )
			VP0(( "Note: %s polarizability %lf beyond normal range\n",
					col[TYPEC], dValue ));
	}

	/*
	 *  Epsilon
	 */
	*iPErrCol = EPSILONC;
	cPData = col[EPSILONC];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in Epsilon field.");
	} 
	if ( dValue < 0.0 ) {
		return("Epsilon cannot be negative.");
	}

	/*
	 *  R*
	 */
	*iPErrCol = RC;
	cPData = col[RC];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in R* field.");
	} else if ( dValue < 0.0 ) {
		return("R* cannot be negative.");
	}

	/*
	 *  Element
	 */
	*iPErrCol = ELEMENTC;
	cPData = col[ELEMENTC];
	if ( strcmp( cPData, "?" ) && !isalpha(*cPData) ) {
		sprintf( SsError, "%s %s",
		        "Element names (Atomic Symbol) must begin with",
			"an alphabetic character.");
		return(SsError);
	} 
	if ( strlen(cPData) > ELEMENTNAMELEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Element name (Atomic Symbol) cannot be longer than",
			ELEMENTNAMELEN-1 );
		return(SsError);
	}

	/*
	 *  Hybrid
	 */
	*iPErrCol = HYBRIDC;
	cPData = col[HYBRIDC];
	if ( strlen(cPData) > HYBRIDNAMELEN-1 ) {
		    sprintf( SsError, "%s %d characters.",
			"Atom Hybridization type cannot be longer than",
			HYBRIDNAMELEN-1 );
		return(SsError);
	}

	/*
	 *  Desc
	 */
	*iPErrCol = DESCC;
	cPData = col[DESCC];
	if ( strlen(cPData) > DESCLEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Parameter Description cannot be longer than",
			DESCLEN-1 );
		return(SsError);
	}

	return(NULL);
}


/*
 *	zXAAPTAccepRow
 * 
 *	Accept/update table blindly - caller must check row 1st.
 */
 
static void
zXAAPTAcceptRow( TABLE tTable, int iRow, STRING col[] )
{
ATOMPARMTABLEt	*aptPCur;
PARMSET		psParmSet;
double		dMass, dPolar, dEps, dR, dEps14 = 0.0, dR14 = 0.0;
double		dScreenF = 0.0;
int		iElement, iHybrid;

	aptPCur = (ATOMPARMTABLEt*)PXATClientPointer(tTable);
	psParmSet = aptPCur->psParmSet;
    
	/*
	 *  convert the numeric fields
	 */

	if ( bStringToDouble( col[MASSC], &dMass ) == FALSE )
		VP0((" mass conversion failed (%s)\n", col[MASSC] ));

	if ( strcmp( col[POLARC], "?" ) == 0 )
		dPolar = -1.0;
	else if ( bStringToDouble( col[POLARC], &dPolar ) == FALSE )
		VP0((" polar conversion failed (%s)\n", col[POLARC] ));

	if ( bStringToDouble( col[EPSILONC], &dEps ) == FALSE )
		VP0((" epsilon conversion failed (%s)\n", col[EPSILONC] ));

	if ( bStringToDouble( col[RC], &dR ) == FALSE )
		VP0((" R* conversion failed (%s)\n", col[RC] ));

	iElement = iElementNumber(col[ELEMENTC]);
	iHybrid = iHybridNumber(col[HYBRIDC]);

	/*
	 *  add or update row
	 */
	if ( iRow > iParmSetTotalAtomParms( psParmSet ) ) {
		DFATAL(( "programming err 1 in zXAAPTAcceptRow\n" ));
	} else if ( iRow == iParmSetTotalAtomParms( psParmSet ) ) {
		/*
		 *  need to add row to table
		 */
		MESSAGE(( "Adding atom parameter %s.\n", col[0] ));
        	if ( iRow != iParmSetAddAtom( psParmSet, col[0], 
				dMass, dPolar, dEps, dR, dEps14, dR14, 
				dScreenF,
				iElement, iHybrid, col[DESCC] ) ) 
			DFATAL(( "programming err 2 in zXABPTAcceptRow\n" ));
	} else {
		/*
		 *  update row in place
		 */
		ParmSetUpdateAtom( psParmSet, iRow, col[TYPEC],
				&dMass, &dPolar, &dEps, &dR, 
				&dScreenF,
				&iElement, &iHybrid, col[DESCC] );

	}
}




/*
 *------------------------------------------------------------------
 *
 *	Public routines
 *
 */

/*
 *	bXAAPTPopupTable (X, Athena, Atom Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XAAPTPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
ATOMPARMTABLEt	*aptPNew;
int		iCount;


	/* First create a place to put the PARM TABLE data */

    MALLOC( aptPNew, ATOMPARMTABLEt*, sizeof(ATOMPARMTABLEt) );

	/* Set the variables that need to be passed on */
	
    aptPNew->fDestroyCallback = fCallback;
    aptPNew->wTop = wWidget;
    aptPNew->psParmSet = psParmSet;
       	
	/* Create the table for the ATOM PARM TABLE */
    iCount = iParmSetTotalAtomParms( psParmSet );

    XATPopupTable(
		  wCreated,
		  wWidget,
		  "aParmTableShell",
		  8, iCount, 	/* 8 columns */
		  "aParmRow",
		  "aParmAddRow",
		  zcPXAAPTGetElement,	
		  zcPXAAPTVerifyElement,	
		  zcPXAAPTVerifyRow,	
		  zXAAPTAcceptRow,
		  NULL,
		  NULL,
		  zXAAPTDestroyTable,
		  0,
		  (GENP)aptPNew );


}


/*
 *	bXAAPTNewPopupTable (X, Athena, Atom Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XAAPTNewPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
ATOMPARMTABLEt	*aptPNew;

	/* First create a place to put the PARM TABLE data */

    MALLOC( aptPNew, ATOMPARMTABLEt*, sizeof(ATOMPARMTABLEt) );

	/* Set the variables that need to be passed on */
	
    aptPNew->fDestroyCallback = fCallback;
    aptPNew->wTop = wWidget;
    aptPNew->psParmSet = psParmSet;
       	
	/* Create the table for the ATOM PARM TABLE */

    XATPopupTable(
		  wCreated,
		  wWidget,
		  "aParmTableShell",
		  8, 1,
		  "aParmAddRow",
		  "aParmAddRow",
		  zcPXAAPTGetElement,
		  zcPXAAPTVerifyElement,
		  zcPXAAPTVerifyRow,
		  zXAAPTAcceptRow,
		  NULL,
		  NULL,
		  zXAAPTDestroyTable,
		  0,
		  (GENP)aptPNew );


}


