/*
 *	File:	xaImproperParmTable.c
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

#include	"xaTable.h"


#define	TYPE1C		0
#define	TYPE2C		1
#define	TYPE3C		2
#define	TYPE4C		3
#define	NC		4
#define	KPC		5
#define P0C		6
#define DESCC		7

#define	MAXTYPELEN	5
#define DESCLEN		32


typedef	struct	{
	PARMSET		psParmSet;
	VFUNCTION	fDestroyCallback;
	Widget		wTop;
} IMPROPERPARMTABLEt;


/*
 *----------------------------------------------------------------
 *
 *	Static variables
 *
 */

static	STRING	SsBuffer;
static	STRING	SsError;


/*
 *	zXAIPTDestroyTable
 *
 *	Callback that is called when the TABLE is destroyed.
 *	This routine cleans up the extra data that was MALLOC'd
 *	for the Improper Torsion properties table.
 */
static void
zXAIPTDestroyTable( TABLE tTable )
{
IMPROPERPARMTABLEt	*iptPTemp;


    iptPTemp = (IMPROPERPARMTABLEt*)PXATClientPointer(tTable);
	/* 
	 *  Check to see if there are any rows left.  
	 *	Set the parmSet to having no Improper parameters
	 */
    if ( !iXATRowsInt(tTable) ) {
    	ParmSetNewImpropers( iptPTemp->psParmSet, 0 );
    }
    iptPTemp->fDestroyCallback(iptPTemp->wTop);
    	/* Destroy the VARARRAY that holds the table->ParmSet references */
    FREE( iptPTemp );
       
}


/*
 *	zcPXATPTGetElement
 *
 *	Get the values for the elements of the TABLE from the
 *	particular Torsion Parameter Entry.
 */
static char *
zcPXAIPTGetElement( TABLE tTable, int iCol, int iRow )
{
#define	DBLFMT	"%1.4lf"
IMPROPERPARMTABLEt	*iptPCur;
char			sType1[MAXTYPELEN];
char			sType2[MAXTYPELEN];
char			sType3[MAXTYPELEN];
char			sType4[MAXTYPELEN];
int			iN;
double			dKp;
double			dP0;
char			sDesc[DESCLEN];


    iptPCur = (IMPROPERPARMTABLEt*)PXATClientPointer(tTable);
    
    ParmSetImproper( iptPCur->psParmSet, iRow, 
	    sType1, sType2, sType3, sType4, &iN, &dKp, &dP0, sDesc );

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
	    case TYPE4C:
	    	strcpy( SsBuffer, sType4 );
	    	MESSAGE(( "Type4:  %s, %s\n", SsBuffer, sType4 ));
		return( SsBuffer );
		break;
	    case NC:
	    	sprintf( SsBuffer, "%i", iN );
		return(SsBuffer);
		break;
	    case KPC:
		sprintf( SsBuffer, DBLFMT, dKp );
		return(SsBuffer);
		break;
	    case P0C:
	    	if ( dP0 ) {
	    	    sprintf( SsBuffer, "Pi" );
	    	} else {
	    	    sprintf( SsBuffer, "0" );
	    	}
		return(SsBuffer);
		break;
	    case DESCC:
	    	strcpy( SsBuffer, sDesc );
		return( SsBuffer );
		break;
	    default:
		DFATAL(("Unexpected column in zcPXAIPTGetElement: %d", iCol));
    }
    return(NULL);	/* for lint */
}


/*
 *	zcPXAIPTVerifyElement
 *
 *	Verify the element to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the error.
 */
static char *
zcPXAIPTVerifyElement( TABLE tTable, int iCol, int iRow, char *cPData )
{
double		dValue;

    switch ( iCol ) {
	    case TYPE1C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type1 must have name.");
		} else if ( !isalpha(*cPData) ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		    sprintf( SsError, "%s %d characters.",
			"Type1 name cannot be longer than",
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
		    sprintf( SsError, "%s %d characters.",
			"Type2 name cannot be longer than",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case TYPE3C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type3 must have name.");
		} else if ( !isalpha(*cPData) ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		    sprintf( SsError, "%s %d characters.",
			"Type3 name cannot be longer than",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case TYPE4C:
		if ( strlen( cPData ) == 0 ) {
		    return("Type4 must have name.");
		} else if ( !isalpha(*cPData) ) {
		    return("First character must be alphabetic.");
		} else if ( strlen(cPData)>MAXTYPELEN-1 ) {
		    sprintf( SsError, "%s %d characters.",
			"Type4 name cannot be longer than",
			MAXTYPELEN-1 );
		    return(SsError);
		}
		break;
	    case KPC:
		if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in Kp field.");
		}
		break;
	    case NC:
	    	if ( !bStringToDouble( cPData, &dValue ) ) {
		    return("Invalid character in P0 field.");
		}
		break;
	    case P0C:
	        StringLower( cPData );
		if ( !strcmp( cPData, "pi" )) break;
		if ( !strcmp( cPData, "0" ))  break;
		return("P0 field must be '0' or 'Pi'.");
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
		DFATAL(("Unexpected column in zcPXAIPTVerifyElement: %d", 
					iCol));
    }
    return(NULL);
}


/*
 *	zcPXAIPTVerifyRow
 *
 *	Verify the row to make sure that it is acceptable.
 *	If it is then return NULL, otherwise return a message
 *	describing the 1st column error.
 */
static char *
zcPXAIPTVerifyRow( Widget w, TABLE tTable, int iRow, STRING *col, int *iPErrCol)
{
char	*cPData;
double	dValue;
int	iValue;


	/*
	 *  Type1
	 */
	*iPErrCol = 0;
	cPData = col[0];
	if ( strlen( cPData ) == 0 ) {
		return("Type1 must have name.");
	} else if ( !isalpha(*cPData) ) {
		return("Type1: first character must be alphabetic.");
	} else if ( strlen(cPData) > MAXTYPELEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Type1 name cannot be longer than",
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
	} else if ( strlen(cPData) > MAXTYPELEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Type2 name cannot be longer than",
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
	} else if ( !isalpha(*cPData) ) {
		return("Type3: first character must be alphabetic.");
	} else if ( strlen(cPData) > MAXTYPELEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Type3 name cannot be longer than",
			MAXTYPELEN-1 );
		return(SsError);
	}

	/*
	 *  Type4
	 */
	*iPErrCol = 3;
	cPData = col[3];
	if ( strlen( cPData ) == 0 ) {
		return("Type4 must have name.");
	} else if ( !isalpha(*cPData) ) {
		return("Type4: first character must be alphabetic.");
	} else if ( strlen(cPData) > MAXTYPELEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Type4 name cannot be longer than",
			MAXTYPELEN-1 );
		return(SsError);
	}

	/*
	 *  Kp
	 */
	*iPErrCol = 4;
	cPData = col[4];
	if ( !bStringToDouble( cPData, &dValue ) ) {
		return("Invalid character in Kp field.");
	} else if ( dValue < 0.0 ) {
		return("Kp cannot be negative.");
	}

	/*
	 *  N
	 */
	*iPErrCol = 5;
	cPData = col[5];
	if ( !bStringToInt( cPData, &iValue ) ) {
		return("Invalid character in N field.");
	} 

	/*
	 *  P0
	 */
	*iPErrCol = 6;
	cPData = col[6];
	StringLower( cPData );
	if ( strcmp( cPData, "pi" )  &&  strcmp( cPData, "0" )) {
		return("P0 field must be '0' or 'Pi'.");
	}

	/*
	 *  Desc
	 */
	*iPErrCol = 7;
	cPData = col[7];
	if ( strlen(cPData) > DESCLEN-1 ) {
		sprintf( SsError, "%s %d characters.",
			"Parameter Description cannot be longer than",
			DESCLEN-1 );
		return(SsError);
	}

	return(NULL);
}


/*
 *	zXATPTAcceptRow
 * 
 *	Accept/update. Caller must do ckecking.
 */
static void
zXAIPTAcceptRow( TABLE tTable, int iRow, STRING col[] )
{
IMPROPERPARMTABLEt	*iptPCur;
PARMSET			psParmSet;
int			iN;
double			dKp;
double			dP0;
double                  dScEE=0.0;
double                  dScNB=0.0;


	iptPCur = (IMPROPERPARMTABLEt*) PXATClientPointer(tTable);
	psParmSet = iptPCur->psParmSet;
    
	/*
	 *  convert the numeric fields
	 */
	bStringToInt( col[4], &iN );
	bStringToDouble( col[5], &dKp );
	StringLower(col[6]);
	if ( !strcmp( col[6], "pi" )  ) {
		dP0 = PI;
	} else {
		dP0 = 0;
	}

	/*
	 *  add or update row
	 */
	if ( iRow > iParmSetTotalImproperParms( psParmSet ) )
		DFATAL(( "programming err 1 in zXAIPTAcceptRow\n" ));

	if ( iRow == iParmSetTotalImproperParms( psParmSet ) ) {
		/*
		 *  need to add row to table
		 */
		MESSAGE(( "Adding improper parameter %s-%s-%s-%s.\n", 
					col[0], col[1], col[2], col[3] ));
        	if ( iRow != iParmSetAddImproperTerm( psParmSet,
					col[0], col[1], col[2], col[3],
					iN, dKp, dP0, dScEE, dScNB, col[7] ) )
			DFATAL(( "programming err 2 in zXAIPTAcceptRow\n" ));
	} else {
		/*
		 *  update row in place
		 */
		ParmSetUpdateImproper( psParmSet, iRow, 
				col[0], col[1], col[2], col[3],
				&iN, &dKp, &dP0, &dScEE, &dScNB, col[7] );
	}
}




/*
 *------------------------------------------------------------------
 *
 *	Public routines
 *
 */



/*
 *	XAIPTPopupTable (X, Athena, Torsion Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XAIPTPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
IMPROPERPARMTABLEt	*iptPNew;
int			iCount;
	

		/* First create a place to put the PARM TABLE data */

    MALLOC( iptPNew, IMPROPERPARMTABLEt*, sizeof(IMPROPERPARMTABLEt) );

    iptPNew->fDestroyCallback = fCallback;
    iptPNew->wTop = wWidget;
    iptPNew->psParmSet = psParmSet;

	/* First create the table for the PARM TABLE */
    iCount = iParmSetTotalImproperParms( psParmSet );
    
    XATPopupTable(
		  wCreated,
		  wWidget,
		  "iParmTableShell",
		  8, iCount,
		  "iParmRow",
		  "iParmAddRow",
		  zcPXAIPTGetElement,
		  zcPXAIPTVerifyElement,
		  zcPXAIPTVerifyRow,
		  zXAIPTAcceptRow,
		  NULL, 
		  NULL,
		  zXAIPTDestroyTable,
		  0,
		  (GENP)iptPNew );
    if ( !iCount ) {
	VP0(( "There are no Improper Torsion parameters in this PARMSET.\n"));
	VP0(( "You may add some using the menu if you wish.\n" ));
    }

}

/*
 *	XAIPTNewPopupTable (X, Athena, Torsion Parameter Table)
 *
 *	Popup a table for editing the parameter files.
 *	Register the callback to be called when the TABLE
 *	is destroyed.
 */
void
XAIPTNewPopupTable( Widget wCreated, Widget wWidget, PARMSET psParmSet, 
	VFUNCTION fCallback )
{
IMPROPERPARMTABLEt	*iptPNew;
	

		/* First create a place to put the PARM TABLE data */

    MALLOC( iptPNew, IMPROPERPARMTABLEt*, sizeof(IMPROPERPARMTABLEt) );

    iptPNew->fDestroyCallback = fCallback;
    iptPNew->wTop = wWidget;
    iptPNew->psParmSet = psParmSet;
   		
    XATPopupTable(
		  wCreated,
		  wWidget,
		  "iParmTableShell",
		  8, 1,
		  "iParmAddRow",
		  "iParmAddRow",
		  zcPXAIPTGetElement,
		  zcPXAIPTVerifyElement,
		  zcPXAIPTVerifyRow,
		  zXAIPTAcceptRow,
		  NULL,
		  NULL,
		  zXAIPTDestroyTable,
		  0,
		  (GENP)iptPNew );

}

