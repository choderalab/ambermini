/*
 *	File:	xaAngleParmTable.c
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
 *	Author:	David A. Rivkin
 *	This file is based on xaAngleParmTable.c
 *
 *	Date Created:	11 August 1992
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

#include        "../Wc/WcCreate.h"

#define TYPE1C		0
#define	TYPE2C		1
#define	TYPE3C		2
#define	KTC		3
#define	T0C		4
#define	DESCC		5

#define MAXTYPELEN	5
#define DESCLEN		32

typedef	struct	{
		PARMSET		psParmSet;
		VFUNCTION	fDestroyCallback;
		Widget		wTop;
		} ANGLEPARMTABLEt;


/*
 *----------------------------------------------------------------
 *
 *	Static variables
 *
 */

static	STRING	SsBuffer;
static	STRING	SsError;


/*
 *	zXAVPTDestroyTable
 *
 *	Callback that is called when the TABLE is destroyed.
 *	This routine cleans up the extra data that was MALLOC'd
 *	for the ANGLE properties table.
 */
static void
zXAVPTDestroyTable( TABLE tTable )
{
ANGLEPARMTABLEt	*vptPTemp; /* v and V for angle */


    vptPTemp = (ANGLEPARMTABLEt*)PXATClientPointer(tTable);

	/* 
	 *  Check to see if there are any rows left.  
	 *	Set the parmSet to having no Angle parameters 
	 */
    if ( !iXATRowsInt(tTable) ) {
    	ParmSetNewAngles( vptPTemp->psParmSet, 0 );
    }
    vptPTemp->fDestroyCallback(vptPTemp->wTop);
    	/* Destroy the VARARRAY that holds the table->ParmSet references */
    FREE( vptPTemp );
       
}

/*
 *	zcPXAVPTGetElement
 *
 *	Get the values for the elements of the TABLE from the
 *	particular Angle Parameter Entry.
 */
static char *
zcPXAVPTGetElement( TABLE tTable, int iCol, int iRow )
{
#define	DBLFMT	"%1.4lf"

ANGLEPARMTABLEt		*vptPCur;
PARMSET			psParmSet;
char			sType1[MAXTYPELEN];
char			sType2[MAXTYPELEN];
char			sType3[MAXTYPELEN];
double			dKt, dT0, dTkub, dRkub;
char			sDesc[DESCLEN];


    vptPCur = (ANGLEPARMTABLEt*)PXATClientPointer(tTable);
    psParmSet = vptPCur->psParmSet;
 
    ParmSetAngle(psParmSet, iRow, sType1, sType2, sType3, &dKt, &dT0, 
		&dTkub, &dRkub, sDesc);

    switch ( iCol ) {
	    case TYPE1C:
	    	strcpy( SsBuffer, sType1 );
	    	return( SsBuffer );
		break;
	    case TYPE2C:
	    	strcpy( SsBuffer, sType2 );
	    	return( SsBuffer );
		break;
	    case TYPE3C:
	    	strcpy( SsBuffer, sType3 );
	    	return( SsBuffer );
		break;
	    case KTC:
	    	sprintf( SsBuffer, DBLFMT, dKt);
		return(SsBuffer);
		break;
	    case T0C:
		sprintf( SsBuffer, DBLFMT, dT0/DEGTORAD);
		return(SsBuffer);
		break;
	    case DESCC:
	    	strcpy( SsBuffer, sDesc );
	    	return( SsBuffer );
		break;
	    default:
		DFATAL(("Unexpected column in zcPXAVPTGetElement: %d", iCol ));
    }
    return(NULL);	/* for lint */
}


/*
 *	zcPXAVPTVerifyElement
 *
 *	Verify the element to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the error.
 */
static char *
zcPXAVPTVerifyElement( TABLE tTable, int iCol, int iRow, char *cPData )
{
double		dValue;

    switch ( iRow ) {
	    case TYPE1C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type1 must have name.");
		} else if ( *cPData >= '0' && *cPData <= '9' ) {
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
		} else if ( *cPData >= '0' && *cPData <= '9' ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		    sprintf( SsError, 
			"Type2 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case TYPE3C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type3 must have name.");
		} else if ( *cPData >= '0' && *cPData <= '9' ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>CONTAINERNAMELEN-1 ) {
		    sprintf( SsError, 
			"Type3 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case KTC:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in Kt field.");
		}
		break;
	    case T0C:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in T0 field.");
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
		DFATAL(("Unexpected column in zcPXAVPTVerifyElement: %d", 
						iCol));
    }
    return(NULL);
}


/*
 *	zcPXAVPTVerifyRow
 *
 *	Verify the row to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the first column error.
 */
static char *
zcPXAVPTVerifyRow(Widget w, TABLE tTable, int iRow, STRING *col, int *iPErrCol )
{
char		*cPData;
double		dValue;

	/*
	 *  Type1
	 */
	*iPErrCol = 0;
	cPData = col[0];
	if ( strlen( cPData ) == 0 ) {
		return("Type1 must have name.");
	} else if ( *cPData >= '0' && *cPData <= '9' ) {
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
	} else if ( *cPData >= '0' && *cPData <= '9' ) {
		return("Type2: first character must be alphabetic.");
	} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		sprintf( SsError, 
			"Type2 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		return(SsError);
	}

	/*
	 *  Type3
	 */
	*iPErrCol = 2;
	cPData = col[2];
	if ( strlen( cPData ) == 0 ) {
		return("Type3 must have name.");
	} else if ( *cPData >= '0' && *cPData <= '9' ) {
		return("Type3: first character must be alphabetic.");
	} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		sprintf( SsError, 
			"Type3 name cannot be longer than %d characters.",
			MAXTYPELEN-1 );
		return(SsError);
	}

	/*
	 *  check/fix/update type order
	 */
	if ( strcmp( col[0], col[2] ) > 0 ) {
		STRING	sTemp;
		SWAP_STRINGS( col[0], col[2], sTemp );
		(void) XawTableSetLabel (w, iRow, 0, col[0] );
		(void) XawTableSetLabel (w, iRow, 2, col[2] );
	}

	/*
	 *  Kt
	 */
	*iPErrCol = 3;
	cPData = col[3];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in Kt field.");
	} else if ( dValue < 0.0 ) {
		return("Kt cannot be negative.");
	}

	/*
	 *  T0
	 */
	*iPErrCol = 4;
	cPData = col[4];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in T0 field.");
	} else if ( dValue < 0.0 ) {
		return("T0 cannot be negative.");
	}

	/*
	 *  Desc
	 */
	*iPErrCol = 5;
	cPData = col[5];
	if ( strlen(cPData) > DESCLEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Parameter Description cannot be longer than",
			DESCLEN-1 );
		return(SsError);
	}

	return(NULL);
}


/*
 *	zXAVPTAcceptRow
 * 
 *	Accept/enter the row blindly - caller must check it 1st.
 */
static void
zXAVPTAcceptRow( TABLE tTable, int iRow, STRING col[] )
{
double		dKt, dT0, zero;
PARMSET		psParmSet;
ANGLEPARMTABLEt	*vptPCur;

	zero = 0.0;
	vptPCur = (ANGLEPARMTABLEt*)PXATClientPointer(tTable);
	psParmSet = vptPCur->psParmSet;
    
	/*
	 *  convert the numeric fields
	 */
	bStringToDouble( col[3], &dKt);
	bStringToDouble( col[4], &dT0);
	dT0 *= DEGTORAD;

	/*
	 *  add or update row
	 */
	if ( iRow > iParmSetTotalAngleParms( psParmSet ) ) {
		DFATAL(( "programming err 1 in zXABPTAcceptAngle\n" ));
	} else if ( iRow == iParmSetTotalAngleParms( psParmSet ) ) {
		/*
		 *  need to add row to table
		 */
		MESSAGE(( "Adding angle parameter %s-%s-%s.\n", 
				col[0], col[1], col[2] ));
		if ( iRow != iParmSetAddAngle( psParmSet, 
					col[0], col[1], col[2],
					dKt, dT0, zero,zero, col[5] ) ) 
			DFATAL(( "programming err 2 in zXABPTAcceptRow\n" ));
	} else {
		/*
		 *  update row in place
		 */
		ParmSetUpdateAngle( psParmSet, iRow, 
					col[0], col[1], col[2],
					&dKt, &dT0, col[5]);
	}
}




/*
 *------------------------------------------------------------------
 *
 *	Public routines
 *
 */



/*
 *	bXAVPTPopupTable (X, Athena, Angle Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XAVPTPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
		VFUNCTION fCallback )
{
ANGLEPARMTABLEt	*vptPNew;
int		iCount;

		/* First create a place to put the PARM TABLE data */

    MALLOC( vptPNew, ANGLEPARMTABLEt*, sizeof(ANGLEPARMTABLEt) );

    vptPNew->fDestroyCallback = fCallback;
    vptPNew->wTop = wWidget;
    vptPNew->psParmSet = psParmSet;

	/* Determine how may Bond Parameters there are so that the nubmer of
		rows can be set */
		
    iCount = iParmSetTotalAngleParms( psParmSet );

		/* First create the table for the PARM TABLE */

    XATPopupTable(
		  wCreated,
		  wWidget,
		  "vParmTableShell",
		  6, iCount,
		  "vParmRow",
		  "vParmAddRow",
		  zcPXAVPTGetElement,
		  zcPXAVPTVerifyElement,
		  zcPXAVPTVerifyRow,
		  zXAVPTAcceptRow,
		  NULL,
		  NULL,
		  zXAVPTDestroyTable,
		  0,
		  (GENP)vptPNew );

    if ( iCount == 0) {
/*    	WcCreateNamedChildren( "*vParmTableShell", "vParmNoParmsShell" );*/
      WcCreateNamedChildren( wWidget, "vParmNoParmsShell" ); /* V.Romanovski */
	/* VP0(( "There are no Angle parameters in this PARMSET.\n"));
	   VP0(("You may add some using the menu if you wish.\n" ));
    	*/
    }

}


/*
 *	bXAVPTNewPopupTable (X, Athena, Angle Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XAVPTNewPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
  ANGLEPARMTABLEt*	vptPNew;

		/* First create a place to put the PARM TABLE data */

    MALLOC( vptPNew, ANGLEPARMTABLEt*, sizeof(ANGLEPARMTABLEt) );

    vptPNew->fDestroyCallback = fCallback;
    vptPNew->wTop = wWidget;
    vptPNew->psParmSet = psParmSet;

		/* First create the table for the PARM TABLE */

    XATPopupTable(
		  wCreated,
		  wWidget,
		  "vParmTableShell",
		  6, 1,
		  "vParmAddRow",
		  "vParmAddRow",
		  zcPXAVPTGetElement,
		  zcPXAVPTVerifyElement,
		  zcPXAVPTVerifyRow,
		  zXAVPTAcceptRow,
		  NULL,
		  NULL,
		  zXAVPTDestroyTable,
		  0,
		  (GENP)vptPNew );
}

