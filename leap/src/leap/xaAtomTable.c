/*
 *      File:   xaAtomTable.c
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
 *              Handle editing of ATOM properties in a table format.
 *
 *              The ATOMs are found using their ATOM ID's
 *              so that ATOMs can be Destroyed in the UnitEditor
 *              while their properties are being edited within
 *              the TABLE.
 *
 * TODO: Take care of the problem that if the CONTAINER for the TABLE
 * TODO: is Destroyed while it's ATOMs properties are being edited then
 * TODO: the TABLE will attempt to write the properties back to
 * TODO: nowhere.  Maybe fix this by REF/DEREF the container.
 */

#include        <X11/IntrinsicP.h>
#include        <X11/StringDefs.h>


#include        "basics.h"
#include        "varArray.h"
#include        "classes.h"
#include        "defaults.h"

#include        "xaTable.h"
#include        "xaAtomTable.h"
#include        "../Wc/WcCreate.h"

#define NAMEC           0
#define TYPEC           1
#define CHARGEC         2
#define ELEMENTC        3
#define PERTURB_ATOM_C  4
#define PERTNAMEC       5
#define PERTTYPEC       6
#define PERTCHARGEC     7
/*
element (i.e. mass) is not perturbed
#define PERTELEMENTC    8
*/

typedef struct  {
                CONTAINER       cContainer;
                VARARRAY        vaAtomIds;
                VFUNCTION       fDestroyCallback;
                Widget          wTop;
                } ATOMTABLEt;


/*
 *----------------------------------------------------------------
 *
 *      Static variables
 *
 */

static  STRING  SsBuffer;
static  STRING  SsError;


/*
 *      zXAATDestroyTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Callback that is called when the TABLE is destroyed.
 *      This routine cleans up the extra data that was MALLOC'd
 *      for the ATOM properties table.
 */
static void
zXAATDestroyTable( TABLE tTable )
{
ATOMTABLEt      *atPTemp;

    atPTemp = (ATOMTABLEt*)PXATClientPointer(tTable);
    atPTemp->fDestroyCallback(atPTemp->wTop);
    VarArrayDestroy(&(atPTemp->vaAtomIds) );
    FREE( atPTemp );
}



/*
 *      zaXAATGetAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the ATOM that the caller is specifying using
 *      the TABLE and the iY value.
 */
static ATOM
zaXAATGetAtom( TABLE tTable, int iRow )
{
ATOMTABLEt      *atPCur;
CONTAINER       cCont;
int             iId;
ATOM            aAtom, aCur;
LOOP            lAtoms;


    atPCur = (ATOMTABLEt*)PXATClientPointer(tTable);
    cCont = atPCur->cContainer;

    iId = *PVAI( atPCur->vaAtomIds, int, iRow );

                /* Find the ATOM with the AtomId */

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    aAtom = NULL;
    while ( (aCur = (ATOM)oNext(&lAtoms)) ) {
        if ( iAtomId(aCur) == iId ) aAtom = aCur;
    }

    return(aAtom);
}




/*
 *      zcPXAATGetElement
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Get the values for the elements of the TABLE from the
 *      particular ATOMs.
 */
static char *
zcPXAATGetElement( TABLE tTable, int iCol, int iRow )
{
#define DBLFMT  "%10.6lf"
ATOM            aAtom;

        aAtom = zaXAATGetAtom( tTable, iRow );

        if ( aAtom == NULL ) 
                return("Atom not found");

        switch ( iCol ) {
            case NAMEC:
                return(sContainerName(aAtom));
                break;
            case TYPEC:
                return(sAtomType(aAtom));
                break;
            case CHARGEC:
                sprintf( SsBuffer, DBLFMT, dAtomCharge(aAtom) );
                return(SsBuffer);
                break;
            case ELEMENTC:
                return(sElementName( iAtomElement(aAtom), SsBuffer ));
                break;
            case PERTURB_ATOM_C:
                if ( bAtomFlagsSet(aAtom,ATOMPERTURB) ) {
                    return("true");
                } else {
                    return("");
                }
                break;
            case PERTNAMEC:
                return(sAtomPertName(aAtom));
                break;
            case PERTTYPEC:
                return(sAtomPertType(aAtom));
                break;
            case PERTCHARGEC:
                sprintf( SsBuffer, DBLFMT, dAtomPertCharge(aAtom) );
                return(SsBuffer);
                break;
/*
            case PERTELEMENTC:
                return(sElementName( iAtomPertElement(aAtom), SsBuffer ));
                break;
*/
        }
        return("unexpected column");
}


/*
 *      zcPXAATVerifyElement
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Verify the element to make sure that it is acceptable.
 *      If it is then return NULL, otherwise return a message
 *      describing the error.
 */
static char *
zcPXAATVerifyElement( TABLE tTable, int iCol, int iRow, char *cPData )
{
ATOM            aAtom;
double          dValue;
STRING          sTemp;

    aAtom = zaXAATGetAtom( tTable, iRow );

    if ( aAtom == NULL ) 
        return("Unexpected error: atom is NULL");

    switch ( iCol ) {
            case NAMEC:
                if ( strlen( cPData ) == 0 )
                    return("Atom must have name.");
                if ( !isalpha(*cPData) )
                    return("First character must be alphabetic.");
                if ( strlen(cPData)>CONTAINERNAMELEN-1 ) {
                    sprintf( SsError, 
                        "Atom name cannot be longer than %d characters.",
                        CONTAINERNAMELEN-1 );
                    return(SsError);
                }
                break;
            case TYPEC:
                if ( strlen(cPData)>ATOMTYPELEN-1 ) {
                    sprintf( SsError, 
                        "Atom type cannot be longer than %d characters.",
                        ATOMTYPELEN-1 );
                    return(SsError);
                }
                break;
            case CHARGEC:
                if ( !bStringToDouble( cPData, &dValue ) ) {
                    return("Invalid character in charge field.");
                }
                break;
            case ELEMENTC:
                /* allow no element in case it's a dummy in pert */
                if ( !strcmp(cPData, NOELEMENTNAME) )
                    break;
                if ( iElementNumber(cPData) == NOELEMENT ) {
                    return("Unknown element name");
                }
                break;
            case PERTURB_ATOM_C:
                if ( strlen(cPData) != 0 ) {
                    strcpy( sTemp, cPData );
                    StringLower( sTemp );
                    if ( !(strcmp( sTemp, "true" ) == 0 || 
                         strcmp( sTemp, "false" ) == 0 )) {
                        return(
                          "Enter 'true', 'false', or nothing (same as false).");
                    }
                }
                break;
            case PERTNAMEC:
                if ( !isalpha(*cPData) )
                    return("First character must be alphabetic.");
                if ( strlen(cPData)>CONTAINERNAMELEN-1 ) {
                    sprintf( SsError, 
                     "Atom perturbed name cannot be longer than %d characters.",
                        CONTAINERNAMELEN-1 );
                    return(SsError);
                }
                break;
            case PERTTYPEC:
                if ( strlen(cPData)>ATOMTYPELEN-1 ) {
                    sprintf( SsError, 
                "Atom perturbed type cannot be longer than %d characters.",
                        ATOMTYPELEN-1 );
                    return(SsError);
                }
                break;
            case PERTCHARGEC:
                if ( !bStringToDouble( cPData, &dValue ) ) {
                    return("Invalid character in perturbed charge field.");
                }
                break;
    }

    return(NULL);
}

/*
 *      zcPXAATVerifyRow
 *
 *      Author: Bill Ross (1995)
 *
 *      Verify the row to make sure that it is acceptable.
 *      If it is then return NULL, otherwise return a message
 *      describing the 1st column error.
 */
static char *
zcPXAATVerifyRow( Widget w, TABLE tTable, int iRow, STRING col[], int *iPErrCol)
{
char            *cPData;
ATOM            aAtom;
double          dValue;
STRING          sTemp;
int             iLen, iPert;

        aAtom = zaXAATGetAtom( tTable, iRow );

        *iPErrCol = 0;
        if ( aAtom == NULL ) 
                return("Unexpected error: atom is NULL");

        /*
         *  Name
         */
        *iPErrCol = 0;
        cPData = col[0];
        iLen = strlen( cPData );

        if ( iLen == 0 )
                return("Atom must have name.");

        if ( !isalpha(*cPData) )
                return("First character of name  must be alphabetic.");

        if ( iLen > CONTAINERNAMELEN-1 ) {
                sprintf( SsError, 
                        "Atom name cannot be longer than %d characters.",
                        CONTAINERNAMELEN-1 );
                return(SsError);
        } 
        if ( iLen > 4 ) {
                VP0((
                 "Warning: atom name length exceeds PDB standard of 4: %s\n", 
                 cPData ));
                ++ tTable->iWarningCount; 
        } else if ( iLen == 4  &&  *cPData != 'H' ) {
                VP0((
                 "Warning: 4 char atom name must be H for PDB standard: %s\n", 
                 cPData ));
                ++ tTable->iWarningCount; 
        }

        /*
         *  Type
         */
        *iPErrCol = 1;
        cPData = col[1];
        if ( strlen(cPData) > ATOMTYPELEN-1 ) {
                sprintf( SsError, 
                        "Atom type cannot be longer than %d characters.",
                        ATOMTYPELEN-1 );
                return(SsError);
        }

        /*
         *  Charge
         */
        *iPErrCol = 2;
        cPData = col[2];
        if ( !bStringToDouble( cPData, &dValue ) ) {
                return("Invalid character in charge field.");
        }

        /*
         *  Element
         */
        *iPErrCol = 3;
        cPData = col[3];
        /* allow no element in case it's a dummy in pert */
        if ( strcmp(cPData, NOELEMENTNAME) ) {
                if ( iElementNumber(cPData) == NOELEMENT ) {
                    return("Unknown element name");
                }
        }

        /*
         *  pert flag
         */
        iPert = 0;
        *iPErrCol = 4;
        cPData = col[4];
        if ( strlen(cPData) != 0 ) {
                strcpy( sTemp, cPData );
                StringLower( sTemp );
                if ( !(strcmp( sTemp, "true" ) == 0 || 
                         strcmp( sTemp, "false" ) == 0 )) {
                        return(
                          "Enter 'true', 'false', or nothing (same as false).");
                }
                iPert = !strcmp( sTemp, "true" );
        }
        if( !GDefaults.iGibbs ) iPert = !strcmp( col[1], col[6] );

        /*
         *  pert name
         */
        *iPErrCol = 5;
        cPData = col[5];
        if ( iPert  &&  !isalpha(*cPData) )
                return("Perturbed name first character must be alphabetic.");
        if ( strlen(cPData) > CONTAINERNAMELEN-1 ) {
                sprintf( SsError, 
                        "Perturbed name cannot be longer than %d characters.",
                        CONTAINERNAMELEN-1 );
                return(SsError);
        }

        /*
         *  pert type
         */
        *iPErrCol = 6;
        cPData = col[6];
        if ( strlen(cPData) > ATOMTYPELEN-1 ) {
                sprintf( SsError, 
                        "Perturbed type cannot be longer than %d characters.",
                        ATOMTYPELEN-1 );
                return(SsError);
        }

        /*
         *  pert charge
         */
        *iPErrCol = 7;
        cPData = col[7];
        if ( !bStringToDouble( cPData, &dValue ) ) {
                return("Invalid character in perturbed charge field.");
        }

        return(NULL);
}

/*
 *      zXAATAcceptRow
 *
 */
static void
zXAATAcceptRow( TABLE tTable, int iRow, STRING col[] )
{
ATOM            aAtom;
double          dValue;
STRING          sTemp;

        aAtom = zaXAATGetAtom( tTable, iRow );

        if ( aAtom == NULL )
                DFATAL(( "atom %d is null (?)", iRow ));

        ContainerSetName( aAtom, col[0] );
        AtomSetType( aAtom, col[1] );
        bStringToDouble( col[2], &dValue );
        AtomSetCharge( aAtom, dValue );
        AtomSetElement( aAtom, iElementNumber(col[3]) );

        if ( strlen(col[4]) != 0 ) {
                strcpy( sTemp, col[4] );
                StringLower( sTemp );
                if ( strcmp( sTemp, "true" ) == 0 ) {
                        AtomSetFlags( aAtom, ATOMPERTURB );
                } else {
                        AtomResetFlags( aAtom, ATOMPERTURB );
                }
        } else {
                AtomResetFlags( aAtom, ATOMPERTURB );
        }
        if( !GDefaults.iGibbs ){
                if( strcmp( col[1], col[6] ) != 0 ) {
                        AtomSetFlags( aAtom, ATOMPERTURB );
                } else {
                        AtomResetFlags( aAtom, ATOMPERTURB );
                }
        }

        AtomSetPertName( aAtom, col[5] );
        AtomSetPertType( aAtom, col[6] );
        bStringToDouble( col[7], &dValue );
        AtomSetPertCharge( aAtom, dValue );
/*
        AtomSetPertElement( aAtom, iElementNumber(col[8]) );
*/
}





/*
 *      zXAATBeginAcceptTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Turn off DISPLAYER updates.
 */
static char *
zXAATBeginAcceptTable( TABLE tTable )
{
    DisplayerAccumulateUpdates();
    return NULL;
}



/*
 *      zXAATEndAcceptTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Turn on DISPLAYER updates.
 */
static char *
zXAATEndAcceptTable( TABLE tTable )
{
    DisplayerReleaseUpdates();
    return NULL;
}





/*
 *------------------------------------------------------------------
 *
 *      Public routines
 *
 */



/*
 *      wXAATPopupTable
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Popup a table for editing the properties of
 *      the ATOMs that are selected within the CONTAINER.
 *      Register the callback to be called when the TABLE
 *      is destroyed.
 */
Widget  
wXAATPopupTable( Widget wCreated, Widget wWidget, UNIT uUnit, 
        VFUNCTION fCallback )
{
  ATOMTABLEt    *atPNew;
  int           iCount;
  LOOP          lRes, lAtoms;
  ATOM          aAtom;
  int           iCur;
  RESIDUE       rRes;
  
                /* count how many ATOMs are selected */

    if ( iObjectType(uUnit) != UNITid ) {
      DFATAL(( "Call to wXAATPopupTable must pass a UNIT.\n" ));
    }
  
  iCount = 0;
  lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
  while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
    if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) 
      iCount++;
  }
  
  if ( iCount == 0 ) {
    VP0(("No atoms selected\n"));
    return(NULL);
  }

                /* create a place to put the ATOM TABLE data */

  MALLOC( atPNew, ATOMTABLEt*, sizeof(ATOMTABLEt) );
  
        /* 
         *  Create the VARARRAY that will be used to store the
         *      ATOM IDs associated with the TABLE rows
         *
         *  Order the atoms by residue sequence number/atom sequence num
         */

  atPNew->vaAtomIds = vaVarArrayCreate(sizeof(int));
  VarArraySetSize( (atPNew->vaAtomIds), iCount );
  iCur = 0;
  lRes = lLoop( (OBJEKT)uUnit, DIRECTCONTENTSBYSEQNUM );
  while ( (rRes = (RESIDUE)oNext(&lRes)) ) {
    lAtoms = lLoop( (OBJEKT)rRes, DIRECTCONTENTSBYSEQNUM );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
      if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
        *PVAI( atPNew->vaAtomIds, int, iCur ) = iAtomId(aAtom);
        iCur++;
      }
    }
  }
  
                /* Store the CONTAINER in the ATOM TABLE data */

  atPNew->cContainer = (CONTAINER) uUnit;
  atPNew->fDestroyCallback = fCallback;
  atPNew->wTop = wWidget;
  
                /* First create the table for the ATOM TABLE */

  XATPopupTable(
                wCreated,
                wWidget,
                "aTableShell",
                8,
                iCount,
                "aRow",
                "",
                zcPXAATGetElement,
                zcPXAATVerifyElement,
                zcPXAATVerifyRow,
                (VFUNCTION)zXAATAcceptRow,
                (VFUNCTION)zXAATBeginAcceptTable,
                (VFUNCTION)zXAATEndAcceptTable,
                zXAATDestroyTable,
                0,
                (GENP)atPNew
                );
  
  return( (Widget)WcFullNameToWidget( wWidget, ".aTableShell"));
  
}

