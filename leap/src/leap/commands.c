/*
 *      File:   commands.c
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
 *             David A. Rivkin                                          *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      All commands typed into the parser are
 *      executed here.
 *
 *      Each function name is prefixed with Cmd
 *      and has two arguments, an integer with
 *      the number of arguments, and an array
 *      of objects which are the arguments.
 *
 */

/*
 *TODO: Add the following commands to the documentation.
 *      New command line options.
 *      addPath
 *      measure
 *      select   - change this one
 *      deSelect
 */

/*        Modifications 
 *        Christine Cezard  (2007)
 *        Universite de Picardie - Jules Verne, Amiens
 *        http://q4md-forcefieldtools.org
 *
 *        Added SaveMol2
 *              Flip
 *           &  Relax
 *              commands
 *        Added case i (integer) in function bCmdMatchTypes
 * 
 * 
 *        Robin Betz (2011, 2013)
 *        UC San Diego
 * 
 *        Added addIonsRand
 *        Added saveAmberParmNetCDF
 * 
*/

/*        Modifications 
 *        Mason Louchart  (2011)
 *        Universite de Picardie - Jules Verne, Amiens
 *        http://q4md-forcefieldtools.org
 *
 *        Added SaveMol3
 *              LoadMol3
 */

#include        <stdio.h>
#include        <stdlib.h>
#include        <float.h>
  
#include        "basics.h"
#include        "vector.h"
#include        "matrix.h"
#include        "classes.h"
#include        "dictionary.h"
#include        "database.h"
#include        "library.h"
#include        "parmLib.h"
#include        "pdbFile.h"
#include        "help.h"
#include        "parser.h"
#include        "tools.h"
#include        "amber.h"
#include        "commands.h"
#include        "defaults.h"
#include        "leap.h"
#include        "octree.h"
#include        "tripos.h"
#include        "build.h"
#include        "zMatrix.h"
#include        "unitio.h"
#include        "mol2File.h"
#include        "mol3File.h"
#include        "minimizer.h"

int     iMemDebug = 0;

extern DICTIONARY       GdVariables;
extern BOOL             bCmdDeleteObj;


#define ATOMSINBOND     2
#define ATOMSINANGLE    3
#define ATOMSINTORSION  4

void VarArrayDeleteMore( VARARRAY header, int pos, int num);
void SelectRelaxInFramework(UNIT uUnit, MINIMIZER mMinimizer);

/* 
 *  COMMANDt    cCommands[]
 *
 *      table mapping command oCmd_ routines to strings
 *      - moved to end of this file.
 */




/*
 *---------------------------------------------------------------------
 *
 *        Private routines
 *
 */


/*
 *      bCmdMatchTypes
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return TRUE if the type of the OBJEKT matches one
 *      of the types in the string up to the first space or
 *      end of string.
 *
 *      (*sPNeedType) will return the name of the type required.
 *
 *      See 'bCmdGoodArguments' for a description of the type characters.
 */
static  BOOL    
bCmdMatchTypes( OBJEKT oObj, char **sPTypes, char *sNeedType, int *iPCount )
{
char    *sTemp;

    *iPCount = 0;
    strcpy( sNeedType, "" );

    sTemp = *sPTypes;
    while ( (**sPTypes!=' ') && (**sPTypes!='\0') ) {
        (*sPTypes)++;
    }
    if ( **sPTypes==' ' ) 
        (*sPTypes)++;
    do {
        switch ( *sTemp ) {
            case '*': 
                return(TRUE);
            case 'u': 
                if ( iObjectType(oObj) == UNITid ) 
                        return(TRUE);
                strcat( sNeedType, "unit" );
                break;
            case 'm': 
                if ( iObjectType(oObj) == MOLECULEid ) 
                        return(TRUE); 
                strcat( sNeedType, "molecule" );
                break;
            case 'r':
                if ( iObjectType(oObj) == RESIDUEid ) 
                        return(TRUE);
                strcat( sNeedType, "residue" );
                break;
            case 'a': 
                if ( iObjectType(oObj) == ATOMid ) 
                        return(TRUE);
                strcat( sNeedType, "atom" );
                break;
            case 'l':
                if ( iObjectType(oObj) == LISTid ) 
                        return(TRUE);
                strcat( sNeedType, "list" );
                break;
            case 'n':
                if ( iObjectType(oObj) == ODOUBLEid ) 
                        return(TRUE);
                strcat( sNeedType, "number" );
                break;
            case 'i':
                if ( iObjectType(oObj) == OINTEGERid ) 
                        return(TRUE);
                strcat( sNeedType, "integer" );
                break;
            case 's':
                if ( iObjectType(oObj) == OSTRINGid ) 
                        return(TRUE);
                strcat( sNeedType, "string" );
                break;
            case 'p':
                if ( iObjectType(oObj) == PARMSETid ) 
                        return(TRUE);
                strcat( sNeedType, "parameter_set" );
                break;
            case 'z':
                if ( oObj == NULL ) 
                        return(TRUE);
                strcat( sNeedType, "null" );
                break;
            case '\0':
            case ' ': 
                return(FALSE);
            default:
                DFATAL(( "ILLEGAL type character in bCmdMatchTypes" ));
        }
        sTemp++;
        if ( *sTemp != ' ' && *sTemp ) {
            strcat( sNeedType, " " );
            (*iPCount)++;
        }
    } while ( *sTemp != '\0' );
    return(FALSE);
}





/*
 *      bCmdGoodArguments
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return TRUE if iArgCount equals iNeeded, otherwise return
 *      FALSE and print an error message.
 *      Also check the types of the arguments against the type string
 *      A string of characters seperated by spaces.
 *      The characters corresponding to:
 *
 *      *       - Matches all types.
 *      u       - Matches a UNIT.
 *      m       - Matches a MOLECULE.
 *      r       - Matches a RESIDUE.
 *      a       - Matches an ATOM.
 *      l       - Matches a LIST.
 *      n       - Matches a NUMBER (ODOUBLE).
 *      i       - Matches an INTEGER
 *      s       - Matches an OSTRING.
 *      p       - Matches a PARMSET.
 *      z       - Matches a NULL.
 *
 *      Eg:  "* umra l s" will match anything in [0],
 *              a UNIT/MOLECULE/RESIDUE/ATOM in [1],
 *              a LIST in [2], and an OSTRING in [3].
 *      
 *      The arguments in the oaArgs array are either OBJEKTs or
 *      they are ASSOCs.  If they are ASSOCs then the OBJEKT within
 *      the ASSOC is tested.
 */
static  BOOL
bCmdGoodArguments( char *sCmd, int iArgCount, ASSOC aaArgs[], char *sTypes )
{
int             i;
OBJEKT          oObj;
int             iNeeded;
STRING          sNeed;
int             iNeedCount;

        /* Count how many arguments are required */

    if ( strlen(sTypes) == 0 ) {
        iNeeded = 0;
    } else {
        iNeeded = 1;
        for ( i=0; i<strlen(sTypes); i++ ) {
            if ( sTypes[i] == ' ' ) 
                iNeeded++;
        }
    }

                /* Check the number of arguments */

    if ( iArgCount != iNeeded ) {
        VP0(( "%s: Improper number of arguments!\n", sCmd ));
        return(FALSE);
    }

                /* Check the types of arguments */
                /* If the OBJEKT is an ASSOC then get the OBJEKT */
                /* within the ASSOC */

    for ( i=0; i<iArgCount; i++ ) {
        oObj = (OBJEKT)aaArgs[i];
        if ( iObjectType(oObj) == ASSOCid ) 
                oObj = oAssocObject(oObj);
        if ( !bCmdMatchTypes( oObj, &sTypes, sNeed, &iNeedCount ) ) {
            VP0(( "%s: Argument #%d is type %s must be of type: [%s]\n", 
                sCmd, i+1, sObjectType(oObj), sNeed ));
            return(FALSE);
        }
    }
    return(TRUE);
}


#if 0
/*
 *      vaCmdListToVarArray
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Copy the OBJEKTs from the LIST into a VARARRAY and
 *      return the VARARRAY.
 */
static VARARRAY
vaCmdListToVarArray( LIST lList )
{
VARARRAY        vaList;
LISTLOOP        llElements;
int             i, iSize;
OBJEKT          oObj;

    iSize = iListSize(lList);
    vaList = vaVarArrayCreate( sizeof(OBJEKT) );
    VarArraySetSize( vaList, iSize );

    i = 0;
    llElements = llListLoop(lList);
    while ( oObj = oListNext(&llElements) ) {
        *PVAI( vaList, OBJEKT, i ) = oObj;
    }
    return(vaList);
}
#endif


        
    


/*
 *=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 *
 *      Actual commands
 * 
 */


    
/*
 *      oCmd_quit
 *
 *      Author: Christian Schafmeister (1991)
 *
 */
OBJEKT
oCmd_quit( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "quit", iArgCount, aaArgs, (char *)"" ) ) {
        VP0(( "usage:  quit\n" ));
        return(NULL);
    }

    GrMainResult.iCommand = CQUIT;
    VP0(( "\tQuit\n" ));
 
    return(NULL);
}



/*
 *      oCmd_describe
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Describe the first OBJEKT in the array.
 */
OBJEKT
oCmd_describe( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "desc", iArgCount, aaArgs, "*" ) ) {
        VP0(( "usage:  desc <variable>\n" ));
        return(NULL);
    }

    Describe( oAssocObject(aaArgs[0]) );

    return(NULL);
}



/*
 *      oCmd_debugOn
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Turn on debugging.
 */
OBJEKT
oCmd_debugOn( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "debugOn", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  debugOn <filename>\n" ));
        return(NULL);
    }

    MessageAddFile( sOString(oAssocObject(aaArgs[0])) );
    VP0(( "Messages will be displayed from the files:\n" ));
    MessageFileList();
    return(NULL);
}


/*
 *      oCmd_debugOff
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Turn off debugging.
 */
OBJEKT
oCmd_debugOff( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "debugOff", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  debugOff\n" ));
        return(NULL);
    }

    MessageRemoveFile( sOString(oAssocObject(aaArgs[0])) );
    VP0(( "Messages will be displayed from the files:\n" ));
    MessageFileList();
    return(NULL);
}




/*
 *      oCmd_debugStatus
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Display status of LEaP, like memory usage etc.
 *
 *      Arguments:
 *      option  [0] - String
 *
 *      If [0] is a string then it can have the following values.
 *              testMemoryOn    - Turn memory testing on.
 *              testMemoryOff   - Turn memory testing off.
 *
 */
OBJEKT
oCmd_debugStatus( int iArgCount, ASSOC aaArgs[] )
{
STRING          sOption;
OBJEKT          oObj;
char            *sCmd = "debugStatus";
char            *sUsage = "usage:  debugStatus [testMemoryOff|testMemoryOn]\n";

    if ( iArgCount == 0 ) {
        if ( !bCmdGoodArguments( sCmd,  iArgCount, aaArgs, "" ) ) 
                return(NULL);
        VP0(( "Current memory usage: %ld bytes\n", GiMemoryAllocated ));
        VP0(( "Memory testing on = %s\n", sBOOL(bTEST_MEMORY_ON()) ));
        return(NULL);
    }
    if ( !bCmdGoodArguments( "debugStatus", iArgCount, aaArgs, "s" ) ) {
        VP0(( sUsage ));
        return(NULL);
    }

    oObj = oAssocObject(aaArgs[0]);
    strcpy( sOption, sOString(oObj) );
    if ( strcmp( sOption, "testMemoryOn" ) == 0 ) {
        TEST_MEMORY_ON( TRUE );
    } else if ( strcmp( sOption, "testMemoryOff" ) == 0 ) {
        TEST_MEMORY_ON( FALSE );
    } else {
        VP0(( "%s: Illegal option: %s\n", sCmd, sOption ));
        VP0(( sUsage ));
    }

    return(NULL);
}



                        


/*
 *      oCmd_help
 *
 *      Author: Christian Schafmeister (1991)
 *      Revised by:  David A. Rivkin (1992)
 *
 *      Offer help.
 */
OBJEKT
oCmd_help( int iArgCount, ASSOC aaArgs[] )
{
int             i = 1;
int             iColumns;
HELP            hTemp;

    iColumns = 4;

    if ( iArgCount == 0 ) {
        VP0(( "Help is available on the following subjects:\n\n" ));
        HelpLoop();
        while ( (hTemp = hHelpNext()) ) {
            if ( i % iColumns ) { 
                /* columns 0-(end-1) */
                VP0(( "%-20s", sHelpUpSubject(hTemp) ));
            } else {
                /* last column in line */
                VP0(( "%s\n", sHelpUpSubject(hTemp) ));
            }
            i++;
        }
        if ( (i-1) % iColumns )
            VP0(( "\n" ));      /* close final part line */

        i = 0;
        VP0(( "\nFor a list of the current aliases, type \"alias\".\n" ));
    } else if ( !bCmdGoodArguments( "help", iArgCount, aaArgs, "s" )) {
            VP0(( "usage:  help <command>\n" ));
    } else {
        hTemp = hHelp( sOString(oAssocObject(aaArgs[0])));
        if ( hTemp == NULL ) {
                VP0(( "No help available.\n" ));
        } else {
                VP0(( "\n" ));
                VP0(( "%s\n", sHelpText(hTemp) ));
        }
    }
    return(NULL);
}





/*
 *      oCmd_list
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      List all variables declared.
 */
OBJEKT
oCmd_list( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "list", iArgCount, aaArgs, "" ) ) 
        VP0(( "usage:  list\n" ));
    else
        VariablesList();
    return(NULL);
}





/*
 *      oCmd_loadOff
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load an OFF file and define all of the variables.
 */
OBJEKT
oCmd_loadOff( int iArgCount, ASSOC aaArgs[] )
{
STRING          sFilename;
char            *sNext;
OBJEKT          oObj;
LIBRARY         ul;

    if ( !bCmdGoodArguments( "loadOff", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  loadOff <filename>\n" ));
        return(NULL);
    }
    strcpy( sFilename, sOString(oAssocObject(aaArgs[0])) );

    ul = lLibraryOpen( sFilename, OPENREADONLY ); 
    if ( ul == NULL ) 
        return(NULL);

    VP0(( "Loading library: %s\n", GsBasicsFullName )); 

    LibraryLoop( ul );
    do {
        sNext = sLibraryNext( ul );
        if ( sNext != NULL ) {
            VP1(( "Loading: %s\n", sNext ));
            oObj = oLibraryLoad( ul, sNext );   /* comes w/ 1 REF */
            VariableSet( sNext, oObj );         /* adds 1 REF */
            DEREF( oObj );                      /* balance REF */

                        /* If the object loaded is a PARMSET then */
                        /* Add it to the PARMLIBRARY */
                        /* And make it the default PARMLIBRARY */

            if ( iObjectType(oObj)==PARMSETid ) {
                PARMSET psTemp = (PARMSET) oObj;

                strcpy( sParmName(psTemp), sFilename );
                ParmLibAddParmSet( GplAllParameters, (PARMSET)oObj );
                ParmLibDefineDefault( GplAllParameters );
            }
        } 
    } while ( sNext != NULL );

    LibraryClose( &ul );
    
    return(NULL);
}
        




/*
 *      oCmd_sequence
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a molecule from a sequence of residues.
 *
 *      Arguments:
 *              [0]     A LIST of units.
 */
OBJEKT
oCmd_sequence( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llElements;
UNIT            uFirst, uSecond, uUnit;
ASSOC           aAss;
LOOP            lTemp, lAtoms;
ATOM            aConnect;
char            *sCmd = "sequence";
RESIDUE         rRes;
int             iDum;

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "l" ) ) {
        VP0(( "usage:  sequence <LIST>\n" ));
        return(NULL);
    }
    llElements = llListLoop( (LIST)oAssocObject(aaArgs[0]) );

                /* Get the first element from the list */

    aAss = (ASSOC)oListNext(&llElements);
    MESSAGE(( "Copying the first UNIT\n" ));
    if ( iObjectType(oAssocObject(aAss)) != UNITid ) {
        VP0(( "%s: Illegal UNIT at position #1.\n", sCmd ));
        return(NULL);
    }

                /* Get the first UNIT and build INTERNALs */

    uFirst = (UNIT)oCopy(oAssocObject(aAss));
    VP1(( "Sequence: %s\n", sContainerName((CONTAINER) uFirst) ));
    while ( (aAss = (ASSOC)oListNext(&llElements)) ) {
        uUnit = (UNIT)oAssocObject(aAss);
        if ( uUnit == NULL ) {
                VP0(( "Unknown UNIT: %s\n", sAssocName(aAss) ));
        } else {

            MESSAGE(( "Copying a subsequent UNIT\n" ));

                /* If the object is not a UNIT then destroy what we have */
                /* up to this point and return */

            if ( iObjectType(uUnit) != UNITid ) {
                VP0(( "%s: Illegal UNIT named: %s\n", sCmd, sAssocName(aAss) ));
                Destroy((OBJEKT *)&uFirst);
                return(NULL);
            }

                        /* Copy the next UNIT */

            uSecond = (UNIT)oCopy( (OBJEKT)uUnit );
            VP1(( "Sequence: %s\n", sContainerName((CONTAINER) uSecond) ));

                        /* Build INTERNALs for the next UNIT */

            MESSAGE(( "Building internals for subsequent UNIT\n" ));
            BuildInternalsForContainer( (CONTAINER) uSecond, 
                        ATOMNEEDSBUILD, ATOMPOSITIONKNOWN );

            aConnect = aUnitHead( uSecond );

            if ( aConnect != NULL ) {
                BuildInternalsBetweenUnitsUsingFlags( uFirst, uSecond,
                                        ATOMNEEDSBUILD,
                                        0 );
            }
            MESSAGE(( "Joining two UNITS deleting the second\n" ));
            UnitSequence( uFirst, uSecond );

            if ( aConnect != NULL ) {
                lAtoms = lLoop( (OBJEKT)aConnect, SPANNINGTREE );
                iDum = 0;       /* for purify */
                BuildExternalsUsingFlags( &lAtoms, ATOMNEEDSBUILD, 0,
                                        ATOMPOSITIONKNOWN,
                                        ATOMNEEDSBUILD|ATOMPOSITIONFIXED,
                                        &iDum, &iDum, &iDum, FALSE );
            }
        }
    }

                /* Destroy INTERNALs to clean up */

    lAtoms = lLoop( (OBJEKT)uFirst, ATOMS );
    BuildDestroyInternals( &lAtoms );

                /* Define PDB sequence */

    lTemp = lLoop( (OBJEKT)uFirst, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lTemp)) ) {
        ResidueSetPdbSequence( rRes, iContainerSequence((CONTAINER) rRes) );
    }

    return((OBJEKT)uFirst);
}



/*
 *      oCmd_loadMol2
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load a UNIT from a Mol2 file.
 *
 *      Arguments:
 *              [0]     - OSTRING, filename.
 */
OBJEKT
oCmd_loadMol2( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fMol2;
UNIT            uUnit;

    if ( !bCmdGoodArguments( "loadMol2", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  <variable> = loadMol2 <filename>\n" ));
        return(NULL);
    }

    fMol2 = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );
    if ( fMol2 == NULL ) return(NULL);

    VP0(( "Loading Mol2 file: %s\n", GsBasicsFullName )); 
    uUnit = uTriposReadUnit( fMol2 );
    fclose(fMol2);

    return((OBJEKT)uUnit);
}


/*___ oCmd_loadMol3 ___________________________________________________________
|                                                                              |
|       Author: Mason Louchart (2011)                                          |
|       http://q4md-forcefieldtools.org                                        |
|       Universite de Picardie - Jules Verne, Amiens                           |
|                                                                              |
|       Tutorial available at                                                  |
|       http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php                 |
|                                                                              |
|       Load a UNIT from a Mol3 file.                                          |
|                                                                              |
|       Arguments:                                                             |
|               [0]     - OSTRING, filename.                                   |
|_____________________________________________________________________________*/
OBJEKT
oCmd_loadMol3( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fMol3;
UNIT            uUnit;

    if ( !bCmdGoodArguments( "loadMol3", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  <variable> = loadMol3 <filename>\n" ));
        return(NULL);
    }

    fMol3 = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );
    if ( fMol3 == NULL ) return(NULL);

    VP0(( "Loading Mol3 file: %s\n", GsBasicsFullName ));
    uUnit = uTriposReadUnit( fMol3 );

    fclose(fMol3);
    return((OBJEKT)uUnit);
}


/*
 *      oCmd_loadPdb
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load a UNIT from a PDB file.
 *
 *      Arguments:
 *              [0]     - OSTRING, filename.
 */
OBJEKT
oCmd_loadPdb( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fPdb;
UNIT            uUnit;

    if ( !bCmdGoodArguments( "loadPdb", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  <variable> = loadPdb <filename>\n" ));
        return(NULL);
    }

    fPdb = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );
    if ( fPdb == NULL ) return(NULL);

    VP0(( "Loading PDB file: %s\n", GsBasicsFullName )); 
    uUnit = uPdbRead( fPdb, NULL );
    fclose(fPdb);

    return((OBJEKT)uUnit);
}



/*
 *      oCmd_loadPdbUsingSeq
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load a UNIT from a PDB file.
 *
 *      Arguments:
 *              [0]     - OSTRING, filename.
 *              [1]     - LIST, list of UNITs to use.
 */
OBJEKT
oCmd_loadPdbUsingSeq( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fPdb;
UNIT            uUnit;
VARARRAY        vaUnits;
LIST            lUnits;
int             i, iErr;
ASSOC           aObj;
OBJEKT          oObj;
LISTLOOP        llLoop;
char            *sCmd = "loadPdbUsingSeq";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s l" ) ) {
        VP0(( "usage:  <variable> = loadPdbUsingSeq <filename> <unitLIST>\n" ));
        return(NULL);
    }

    fPdb = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );
    if ( fPdb == NULL ) return(NULL);

                /* Copy the list of UNITs into a VARARRAY */

    lUnits = (LIST)oAssocObject(aaArgs[1]);
    vaUnits = vaVarArrayCreate( sizeof(UNIT) );
    VarArraySetSize( vaUnits, iCollectionSize(lUnits) );
    llLoop = llListLoop(lUnits);
    i = 0;
    iErr = 0;
    while ( (aObj = (ASSOC)oListNext(&llLoop)) ) {
        oObj = oAssocObject(aObj);
        *PVAI( vaUnits, UNIT, i ) = (UNIT)oObj;
        if ( iObjectType(oObj) != UNITid ) {
            VP0(( "%s: %s is not a unit!\n", sCmd, sAssocName(aObj) ));
            iErr++;
        }
        i++;
    }
    if ( iErr ) {
        VarArrayDestroy( &vaUnits );
        fclose(fPdb);
        VP0(( "Not loaded\n" ));
        return(NULL);
    }


    VP0(( "Loading PDB file: %s using sequence %s\n", 
                GsBasicsFullName, sAssocName(aaArgs[1]) )); 
    uUnit = uPdbRead( fPdb, vaUnits );

    VarArrayDestroy( &vaUnits );

    fclose(fPdb);

    return((OBJEKT)uUnit);
}






/*
 *      oCmd_saveOff
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save a UNIT/PARMSET or a LIST of UNITs/PARMSETs to a UNITLIBRARY.
 *
 *      Arguments:
 *              [0]     - UNIT or LIST of UNITs to save.
 *              [1]     - OSTRING filename.
 */
OBJEKT
oCmd_saveOff( int iArgCount, ASSOC aaArgs[] )
{
LIBRARY         ul;
LISTLOOP        llUnits;
OBJEKT          oObj;
ASSOC           aObj;
char            *sCmd = "saveOff";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "upl s" ) ) {
        VP0(( "usage:  saveOff <object> <filename>\n" ));
        return(NULL);
    }

    ul = lLibraryOpen( sOString(oAssocObject(aaArgs[1])), OPENREADWRITE );
    if ( ul==NULL ) return(NULL);

    DisplayerAccumulateUpdates();
    if ( iObjectType( oAssocObject(aaArgs[0]) ) == UNITid ||
         iObjectType( oAssocObject(aaArgs[0]) ) == PARMSETid ) {
        VP1(( "Saving %s.\n", sAssocName(aaArgs[0]) ));
        LibrarySave( ul, sAssocName(aaArgs[0]), 
                        oAssocObject(aaArgs[0]), NULL );
    } else {
        oObj = oAssocObject(aaArgs[0]);
        llUnits = llListLoop( (LIST)oObj );
        while ( (aObj = (ASSOC)oListNext(&llUnits) ) != NULL ) {
            oObj = oAssocObject(aObj);
            if ( iObjectType(oObj) != UNITid &&
                 iObjectType(oObj) != PARMSETid ) {
                VP0(( "%s: Cannot save %s - type %s (ignoring).\n", 
                        sCmd, sAssocName(aObj), sObjectType(oObj) ));
            } else {
                VP1(( "Saving %s.\n", sAssocName(aObj) ));
                LibrarySave( ul, sAssocName(aObj), 
                                oAssocObject(aObj), NULL );
            }
        }
    }
    DisplayerReleaseUpdates();

    LibraryClose( &ul );
    return(NULL);
}


/*
 *      Davids Changes
 *
 */
 
/*
 *      oCmd_createParmset
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created PARMSET.
 *
 *      Arguments:
 *              [0]     - OSTRING PARMSET name.
 */
OBJEKT 
oCmd_createParmset( int iArgCount, ASSOC aaArgs[] )
{
PARMSET         psParmSet;

    if ( !bCmdGoodArguments( "createParmset",  iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  <variable> = createParmset <name>\n" ));
        return(NULL);
    }

    psParmSet = (PARMSET)oCreate(PARMSETid);
    AssocSetObject( aaArgs[0], (OBJEKT)psParmSet );
    return((OBJEKT)psParmSet);
}




/*
 *      oCmd_createUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created UNIT.
 *
 *      Arguments:
 *              [0]     - OSTRING UNIT name.
 */
OBJEKT
oCmd_createUnit( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;

    if ( !bCmdGoodArguments( "createUnit", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  <variable> = createUnit <name>\n" ));
        return(NULL);
    }

    uUnit = (UNIT)oCreate(UNITid);
    ContainerSetName( (CONTAINER) uUnit, sOString(oAssocObject(aaArgs[0])) );
    return((OBJEKT)uUnit);
}



/*
 *      oCmd_moleculeCreate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created MOLECULE.
 *
 *      Arguments:
 *              [0]     - OSTRING MOLECULE name.
 */
OBJEKT
oCmd_moleculeCreate( int iArgCount, ASSOC aaArgs[] )
{
MOLECULE        mMol;

/*
TODO - no help for this.. write & rename to createMolecule
or delete this cmd
*/
    if ( !bCmdGoodArguments( "moleculeCreate", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  <variable> = moleculecreate <name>\n" ));
        return(NULL);
    }

    mMol = (MOLECULE)oCreate(MOLECULEid);
    ContainerSetName( (CONTAINER) mMol, sOString(oAssocObject(aaArgs[0])) );
    return((OBJEKT)mMol);
}



/*
 *      oCmd_createResidue
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created RESIDUE.
 *
 *      Arguments:
 *              [0]     - OSTRING RESIDUE name.
 */
OBJEKT
oCmd_createResidue( int iArgCount, ASSOC aaArgs[] )
{
RESIDUE         rRes;

    if ( !bCmdGoodArguments( "createResidue", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  <variable> = createResidue <name>\n" ));
        return(NULL);
    }

    rRes = (RESIDUE)oCreate(RESIDUEid);
    ContainerSetName( (CONTAINER) rRes, sOString(oAssocObject(aaArgs[0])) );
    return((OBJEKT)rRes);
}






/*
 *      oCmd_createAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created ATOM.
 *
 *      Arguments:
 *              [0]     - OSTRING ATOM name.
 *              [1]     - OSTRING ATOM type.
 *              [2]     - ODOUBLE ATOM charge.
 */
OBJEKT
oCmd_createAtom( int iArgCount, ASSOC aaArgs[] )
{
ATOM            aAtom;
char            *sCmd = "createAtom";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s s n" ) ) {
        VP0(( "usage:  <variable> = createAtom <name> <type> <charge>\n" ));
        return(NULL);
    }

    aAtom = (ATOM)oCreate(ATOMid);
    ContainerSetName( (CONTAINER) aAtom, sOString(oAssocObject(aaArgs[0])) );
    AtomSetType( aAtom, sOString(oAssocObject(aaArgs[1])) );
    AtomSetCharge( aAtom, dODouble(oAssocObject(aaArgs[2])) );

    return((OBJEKT)aAtom);
}





/*
 *      oCmd_add
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add the oB to oA.
 *      This can only be done if oB is not contained by anything else.
 *
 *      Arguments:
 *              [0]     - OBJEKT oA
 *              [1]     - OBJEKT oB
 */
OBJEKT
oCmd_add( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oA, oB;
int             iA, iB;
STRING          sTemp;
char            *sCmd = "add";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umr mra" ) ) {
        VP0(( "usage:  add <unit/residue/atom> <unit/residue/atom>\n" ));
        return(NULL);
    }

    oA = oAssocObject(aaArgs[0]);
    oB = oAssocObject(aaArgs[1]);
    iA = iObjectType(oA);
    iB = iObjectType(oB);
    if ( (iA==UNITid) && 
        !((iB==MOLECULEid)||(iB==RESIDUEid)||(iB==ATOMid)) ) {
        VP0(( "%s: UNITs cannot contain UNITs.\n", sCmd ));
        return(NULL);
    }
    if ( (iA==MOLECULEid) && 
        !((iB==RESIDUEid)||(iB==ATOMid)) ) {
        VP0(( "%s: MOLECULEs can only contain RESIDUEs and ATOMs.\n", sCmd ));
        return(NULL);
    }
    if ( (iA==RESIDUEid) && (iB!=ATOMid) ) {
        VP0(( "%s: Residues can only contain ATOMs.\n", sCmd ));
        return(NULL);
    }
    if ( cContainerWithin((CONTAINER) oB) != NULL ) {
        VP0(( "%s: The object %s is already contained within %s\n",
                sCmd, sAssocName(aaArgs[1]), 
                sContainerFullDescriptor( cContainerWithin((CONTAINER) oB), sTemp ) ));
        return(NULL);
    }
    ContainerAdd( (CONTAINER)oA, oB );
    return(NULL);
}




/*
 *      oCmd_remove
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove oB from oA if it is contained by oA.
 *
 *      Arguments:
 *              [0]     - OBJEKT oA
 *              [1]     - OBJEKT oB
 */
OBJEKT
oCmd_remove( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oA, oB;
char            *sCmd = "remove";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umr mra" ) ) {
        VP0(( "usage:  remove <unit/residue/atom> <unit/residue/atom>\n" ));
        return(NULL);
    }

    oA = oAssocObject(aaArgs[0]);
    oB = oAssocObject(aaArgs[1]);

    REF( oB );  /* bContainerRemove() needs this */
    if ( !bContainerRemove( (CONTAINER)oA, oB ) ) {
        VP0(( "%s: Could not find %s within %s.\n",
                sCmd, sAssocName(aaArgs[1]), sAssocName(aaArgs[0]) ));
    }
    DEREF( oB ); /* reset count after bContainerRemove */

    return(NULL);
}



/*
 *      oCmd_bond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a bond between two atoms of the appropriate order.
 *      If the bond order is not given then use a single bond.
 *
 *      Arguments:
 *              [0]     - ATOM aA
 *              [1]     - ATOM aB
 *      option  [2]     - OSTRING, bond order.
 */
OBJEKT
oCmd_bond( int iArgCount, ASSOC aaArgs[] )
{
ATOM            aA, aB;
int             iOrder;
char            cOrder;
char            *sCmd = "bond";

   if ( iArgCount == 2 ) {
        cOrder = 'S';
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a" ) ) {
            VP0(( "usage:  bond <atom1> <atom2> [order]\n" ));
            return(NULL);
        }
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a s" ) ) {
            VP0(( "usage:  bond <atom1> <atom2> [order]\n" ));
            return(NULL);
        }
        cOrder = sOString(oAssocObject(aaArgs[2]))[0];
    }
    
    DisplayerAccumulateUpdates();
    
    aA = (ATOM)oAssocObject(aaArgs[0]);
    aB = (ATOM)oAssocObject(aaArgs[1]);
  
    switch ( cOrder ) {
        case 'S': 
            iOrder = BONDSINGLE;
            break;
        case 'D': 
            iOrder = BONDDOUBLE;
            break;
        case 'T': 
            iOrder = BONDTRIPLE;
            break;
        case 'A': 
            iOrder = BONDAROMATIC;
            break;
        default:
            VP0(( "%s: Unknown bond order, no bond made.\n Valid bond orders are: S, D, T, A.\n", sCmd ));
            goto DONE;
    }
    AtomBondToOrder( aA, aB, iOrder );
DONE:
    DisplayerReleaseUpdates();
    return(NULL);
}





/*
 *      oCmd_deleteBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove a bond between two atoms.
 *      Print a message if the bond does not exist.
 *
 *      Arguments:
 *              [0]     - ATOM aA
 *              [1]     - ATOM aB
 */
OBJEKT
oCmd_deleteBond( int iArgCount, ASSOC aaArgs[] )
{
ATOM            aA, aB;
char            *sCmd = "deleteBond";

                /* Check the arguments */

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a" ) ) {
        VP0(( "usage:  deleteBond <atom1> <atom2>\n" ));
        return(NULL);
    }

    DisplayerAccumulateUpdates();

    aA = (ATOM)oAssocObject(aaArgs[0]);
    aB = (ATOM)oAssocObject(aaArgs[1]);

    if ( bAtomBondedTo( aA, aB ) ) {
        AtomRemoveBond( aA, aB );
    } else 
        VP0(( "%s: Atoms are not bonded.\n", sCmd ));
    
    DisplayerReleaseUpdates();
    
    return(NULL);
}



/*
 *      oCmd_zMatrix
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Construct the external coordinates for the atoms.
 *
 *      Arguments:
 *              [0]     - Container that contains the atoms.
 *              [1]     - A list of atoms and internal coordinates.
 *
 *      The entries in the list of atoms and internal coordinates can
 *      look like:
 *
 *      a1 a2 b12
 *      a1 a2 a3 b12 t123
 *      a1 a2 a3 a4 b12 t123 p1234
 *      a1 a2 a3 a4 b12 t123 t124 orientation
 *
 *      Where a1,a2,a3,a4 can be an atom or an atom name which exists
 *      in the container.
 */
OBJEKT
oCmd_zMatrix( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llLines, llElements;
LIST            lLine;
int             iCount;
OBJEKT          oaElements[10];
VECTOR          vNew, vAtom2, vAtom3, vAtom4;
CONTAINER       cCont;
OBJEKT          oObj;
ASSOC           aAss, aAss2;
char            *sCmd = "zMatrix";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umr l" ) ) {
        VP0(( "usage:  zMatrix <unit/residue> <LIST>\n" ));
        return(NULL);
    }

    DisplayerAccumulateUpdates();
    
    cCont = (CONTAINER)oAssocObject(aaArgs[0]);
    llLines = llListLoop( (LIST)oAssocObject(aaArgs[1]) );    
    while ( (aAss = (ASSOC)oListNext(&llLines) ) != NULL ) {
        lLine = (LIST)oAssocObject(aAss);
        if ( iObjectType(lLine) != LISTid ) {
            VP0(( "%s: Illegal object in zMatrix list was ignored.\n", sCmd ));
            continue;
        }
        llElements = llListLoop( lLine );
        iCount = 0;
        while ( ( aAss2 = (ASSOC)oListNext(&llElements)) != NULL ) {
            oObj = oAssocObject(aAss2);
            if ( iObjectType(oObj) == OSTRINGid ) {
                oObj = (OBJEKT)cContainerFindName( cCont, 
                                ATOMid, sOString(oObj) );
            } 
            oaElements[iCount] = oObj;
            iCount++;
        }
        
        switch ( iCount ) {
            case 3:
                if ( !bCmdGoodArguments( sCmd, iCount, (ASSOC *)oaElements, 
                                                        "a a n" ) ) 
                    goto BAD;
                ZMatrixNothing( &vAtom2 );
                AtomSetPosition( oaElements[1], vAtom2 );
                ZMatrixBond( &vNew, &vAtom2, dODouble(oaElements[2]) );
                AtomSetPosition( oaElements[0], vNew );
                break;
            case 5:
                if ( !bCmdGoodArguments( sCmd, iCount, (ASSOC *)oaElements, 
                                                        "a a a n n"))
                    goto BAD;
                vAtom2 = vAtomPosition( oaElements[1] );
                vAtom3 = vAtomPosition( oaElements[2] );
                ZMatrixBondAngle( &vNew, &vAtom2, &vAtom3, 
                                dODouble(oaElements[3]),
                                dODouble(oaElements[4])*DEGTORAD );
                AtomSetPosition( oaElements[0], vNew );
                break;
            case 7:
                if ( !bCmdGoodArguments( sCmd, iCount, (ASSOC *)oaElements, 
                                                "a a a a n n n" ) )
                    goto BAD;
                vAtom2 = vAtomPosition( oaElements[1] );
                vAtom3 = vAtomPosition( oaElements[2] );
                vAtom4 = vAtomPosition( oaElements[3] );
                ZMatrixBondAngleTorsion( &vNew, &vAtom2, &vAtom3, &vAtom4,
                                dODouble(oaElements[4]),
                                dODouble(oaElements[5])*DEGTORAD,
                                dODouble(oaElements[6])*DEGTORAD );
                AtomSetPosition( oaElements[0], vNew );  
                break;
            case 8:
                if (!bCmdGoodArguments( sCmd, iCount, (ASSOC *)oaElements, 
                                                "a a a a n n n n" ))
                    goto BAD;
                vAtom2 = vAtomPosition( oaElements[1] );
                vAtom3 = vAtomPosition( oaElements[2] );
                vAtom4 = vAtomPosition( oaElements[3] );
                ZMatrixBondTwoAnglesOrientation( &vNew, 
                                &vAtom2, &vAtom3, &vAtom4,
                                dODouble(oaElements[4]),
                                dODouble(oaElements[5])*DEGTORAD,
                                dODouble(oaElements[6])*DEGTORAD,
                                dODouble(oaElements[7]) );
                AtomSetPosition( oaElements[0], vNew );  
                break;
            default:
                VP0(( "%s: Illegal zMatrix entry was ignored.\n", sCmd ));
                goto DONE;
        }
        continue;
BAD:    
        VP0(( "%s: Illegal object in zMatrix entry.  Entry was ignored\n",
                                                sCmd ));
        continue;
    }

DONE:
    DisplayerReleaseUpdates();
    return(NULL);
}



/*
 *      oCmd_saveOffParm
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save a UNIT or a LIST of UNITs to a UNITLIBRARY.
 *      Save it WITH parameters!!!!!!!!!!
 *
 *      Arguments:
 *              [0]     - UNIT or LIST of UNITs to save.
 *              [1]     - OSTRING filename.
 */
OBJEKT
oCmd_saveOffParm( int iArgCount, ASSOC aaArgs[] )
{
LIBRARY         ul;
LISTLOOP        llUnits;
UNIT            uUnit;
char            *sCmd = "saveOffParm";

    VP0(("saveOffParm: command deactivated\n" ));
    return(NULL);
    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "ul s" ) ) {
        VP0(( "usage:  saveOffParm <object> <filename>\n" ));
        return(NULL);
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
        return(NULL);
    }
    
    ul = lLibraryOpen( sOString( oAssocObject(aaArgs[1]) ), OPENREADWRITE );

    if ( iObjectType( oAssocObject(aaArgs[0]) ) == UNITid ) {
        LibrarySave( ul, sContainerName((CONTAINER) oAssocObject(aaArgs[0])), 
                        oAssocObject(aaArgs[0]), 
                        GplAllParameters );
    } else {
        llUnits = llListLoop( (LIST)oAssocObject(aaArgs[0]) );
        while ( (uUnit = (UNIT)oListNext(&llUnits) ) != NULL ) {
            if ( iObjectType(uUnit) != UNITid ) {
                VP0(( "%s: Illegal UNIT in list was ignored.\n", sCmd ));
                continue;
            }
            LibrarySave( ul, sContainerName((CONTAINER) uUnit), (OBJEKT) uUnit, GplAllParameters );
        }
    }

    LibraryClose( &ul );
    return(NULL);
}




/*
 *      oCmd_loadAmberPrep
 *
 *      Author: Christian Schafmeister (1991)
 *      Modified: David A. Rivkin ( 1-15-93 ) Fixed the Dummy atom reading
 *
 *      Load an AMBER PREP file, add all of the UNITs into
 *      the systems Variable list.
 *
 *      Arguments:
 *              [0] -   OSTRING, filename.
 *      option  [1] -   OSTRING, string to prefix to each unit name.
 *
 *      The option [1] is provided to change the names of Old AMBER
 *      types to distinguish All Atom residues from United Atom residues
 *      and to distinguish Terminating residues from main chain 
 *      residues.
 */
OBJEKT
oCmd_loadAmberPrep( int iArgCount, ASSOC aaArgs[] )
{
DICTIONARY      dUnits;
DICTLOOP        dlLoop;
UNIT            uUnit;

STRING          sName, sPrefix;
LOOP            lResidues;
RESIDUE         rRes;
char            *sCmd = "loadAmberPrep";

    if ( iArgCount == 1 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s" ) ) {
            VP0(( "usage:  loadAmberPrep <filename> [prefix]\n" ));
            return(NULL);
        }
        strcpy( sPrefix, "" );
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s s" ) ) {
            VP0(( "usage:  loadAmberPrep <filename> [prefix]\n" ));
            return(NULL);
        }
        strcpy( sPrefix, sOString(oAssocObject(aaArgs[1])) );
    }
    dUnits = dAmberReadPrepFile( sOString(oAssocObject(aaArgs[0])) );
    if ( dUnits != NULL ) {
        dlLoop = ydlDictionaryLoop( dUnits );
        while ( (uUnit=(UNIT)yPDictionaryNext( dUnits, &dlLoop )) != NULL ) {
            strcpy( sName, sPrefix ); 
            strcat( sName, sContainerName((CONTAINER) uUnit) );
            VP1(( "Loaded UNIT: %s\n", sName ));

                /* Set the name for the UNIT */

            ContainerSetName( (CONTAINER) uUnit, sName );

                /* Set the name for the only residue in the unit */

            lResidues = lLoop( (OBJEKT)uUnit, RESIDUES );
            rRes = (RESIDUE)oNext(&lResidues);
            ContainerSetName( (CONTAINER) rRes, sName );    
            VariableSet( sName, (OBJEKT) uUnit );       /* adds 1 REF */
        }
    }
    return(NULL);
}

/*
 *      oCmd_loadAmberParams
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load an AMBER parameter file, return the PARMSET 
 *      and add the PARMSET to the systems PARMLIBRARY.
 *
 *      Arguments:
 *              [0] -   OSTRING, filename.
 *
 */
OBJEKT
oCmd_loadAmberParams( int iArgCount, ASSOC aaArgs[] )
{
PARMSET         psParms;
STRING          sFile;
FILE            *fIn;
char            *sUsage = 
                  "usage:  <variable> = loadAmberParams <filename> \n";

    if ( !bCmdGoodArguments( "loadAmberParams", iArgCount, aaArgs, "s" ) ) {
        VP0(( sUsage ));
        return(NULL);
    }
    strcpy( sFile, sOString(oAssocObject(aaArgs[0])) );

    fIn = FOPENCOMPLAIN( sFile, "r" );
    if ( fIn == NULL ) 
        return(NULL);

    VP0(("Loading parameters: %s\n", GsBasicsFullName ));

    psParms = psAmberReadParmSet( fIn, sFile );

    if ( psParms != NULL ) {
        ParmLibAddParmSet( GplAllParameters, psParms );
        ParmLibDefineDefault( GplAllParameters );
    } else
        VP0(( "-- no parameters loaded" ));

    return((OBJEKT)psParms);
}



/*
 *      oCmd_savePdb
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save the UNIT as a PDB file.
 *      Create a duplicate of the UNIT and impose a hierarchy on
 *      it.
 *
 *      Arguments:
 *              [0] -   Unit to save.
 *              [1] -   Name of file to save to.
 */
OBJEKT
oCmd_savePdb( int iArgCount, ASSOC aaArgs[] )
{
char            *sString;
FILE            *fOut;

    if ( !bCmdGoodArguments( "savePdb", iArgCount, aaArgs, "u s" ) ) {
        VP0(( "usage:  savePdb <object> <filename>\n" ));
        return(NULL);
    }
    sString = sOString( oAssocObject(aaArgs[1]) );
    fOut = FOPENCOMPLAIN( sString, "w" );
    if ( fOut == NULL ) 
        return(NULL);
    VP0(("Writing pdb file: %s\n", sString));
    PdbWrite( fOut, (UNIT)oAssocObject(aaArgs[0]) );
    fclose( fOut );
    return(NULL);
}


/*
 *      oCmd_saveMol2
 *      Based on savepdb
 *      Author: Christine Cezard (2007)
 *      Universite de Picardie - Jules Verne, Amiens
 *      http://q4md-forcefieldtools.org
 *
 *      Tutorial available at
 *      http://q4md-forcefieldtools.org/Tutorial/leap.php
 *
 *      Save the UNIT as a Mol2 file.
 *      Create a duplicate of the UNIT and impose a hierarchy on it.
 *
 *      Arguments:
 *              [0] -   Unit to save.
 *              [1] -   Name of file to save to.
 *              [2] -   Option for column 6 (0 = Default, 1 = Amber Atom type) 
 */
OBJEKT
oCmd_saveMol2( int iArgCount, ASSOC aaArgs[] )
{
char            *sString;
FILE            *fOut;
double          choice;

    if ( !bCmdGoodArguments( "saveMol2", iArgCount, aaArgs, "u s n" ) ) {
        VP0(( "usage:  saveMol2 <object> <filename> <option>\n" ));
        VP0(( "<option> = 0 for Default atom types \n" ));
        VP0(( "<option> = 1 for AMBER atom types \n" ));
        return(NULL);
    }
    sString = sOString( oAssocObject(aaArgs[1]) );
    choice =  (int) dODouble( oAssocObject(aaArgs[2]) ) ;
    fOut = FOPENCOMPLAIN( sString, "w" );
    if ( fOut == NULL ) 
        return(NULL);
    VP0(("Writing mol2 file: %s\n", sString));
    Mol2Write( fOut, (UNIT)oAssocObject(aaArgs[0]), choice);
    fclose( fOut );
    return(NULL);
}


/*_____ oCmd_saveMol3 _________________________________________________________
|                                                                              |
|       Based on saveMol2                                                      |
|       Author: Mason Louchart (2011)                                          |
|       http://q4md-forcefieldtools.org                                        |
|       Universite de Picardie - Jules Verne, Amiens                           |
|                                                                              |
|       Tutorial available at                                                  |
|       http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php                 |
|                                                                              |
|       Save the UNIT as a Mol3 file.                                          |
|                                                                              |
|       Arguments:                                                             |
|               [0] -  Unit to save.                                           |
|               [1] -  Name of file to save to.                                |
|               [2] -  Option for column 6 (0 = Default, 1 = Amber Atom type)  |
|_____________________________________________________________________________*/
OBJEKT
oCmd_saveMol3( int iArgCount, ASSOC aaArgs[] )
{
char            *sString;
FILE            *fOut;
double          choice;

    if ( !bCmdGoodArguments( "saveMol3", iArgCount, aaArgs, "u s n" ) ) {
        VP0(( "usage:  saveMol3 <object> <filename> <option>\n" ));
        VP0(( "<option> = 0 for Default atom types \n" ));
        VP0(( "<option> = 1 for AMBER atom types \n" ));
        return(NULL);
    }
    sString = sOString( oAssocObject(aaArgs[1]) );
    choice =  (int) dODouble( oAssocObject(aaArgs[2]) );
    fOut = FOPENCOMPLAIN( sString, "w" );
    if ( fOut == NULL ) return(NULL);
        
    VP0(("Writing mol3 file: %s\n", sString));
    
    Mol3Write( fOut, (UNIT)oAssocObject(aaArgs[0]), choice);
    fclose( fOut );
    return(NULL);
}


/*
 *      oCmd_solvateBox
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Solvate the UNIT [0] within copies of a box of SOLVENT.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *              [2] -   LIST with ( x, y, z ) distances or
 *                      xyz, the closest the wall of the solvent box
 *                      can come to the smallest box which contains
 *                      the entire unit that is centered on the
 *                      origin.  This is called the BufferZone.
 *                      Or a single ODOUBLE when x,y,z distances are the
 *                      same.
 *      Option  [3] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateBox( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0, daBuffer[3];
ASSOC           aAss;
LISTLOOP        llNumbers;
int             i, iInitialSize, iFinalSize;
BOOL            bIsotropic = FALSE;
char            *sCmd = "solvateBox";
char            *sUsage =
                 "usage:  solvateBox <solute> <solvent> <buffer> [iso] [closeness]\n";

    /*
     *  check args - always need 2 units & cutoff
     */
    if ( !bCmdGoodArguments( sCmd, 3, aaArgs, "u u ln" ) ) {
        VP0(( sUsage ));
        return(NULL);
    }
    if ( iArgCount > 5 ) {
        VP0(( sUsage ));
        return(NULL);
    }
    if ( iArgCount != 3 ) {
        OBJEKT  oObj;

        /*
         *  handle possible 'iso' &/or closeness
         */
        oObj = (OBJEKT)aaArgs[3];
        if ( iObjectType(oObj) != ASSOCid )
                DFATAL(( "unexpected objtype %s\n", sObjectType(oObj) ));

        if ( strcmp("iso", sAssocName(oObj) ) ) {
                oObj = oAssocObject(oObj);
        } else {
                bIsotropic = TRUE;
                if ( iArgCount == 5 ) {
                        oObj = (OBJEKT)aaArgs[4];
                        oObj = oAssocObject(oObj);
                } else
                        oObj = NULL;
        }
        if ( oObj != NULL ) {
                if ( iObjectType(oObj) != ODOUBLEid ) {
                        VP0(( sUsage ));
                        return(NULL);
                }
                dCloseness = dODouble( oObj );
                if ( dCloseness <= 0.0 ) {
                        VP0(( sUsage ));
                        return(NULL);
                }
        }
    }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
        return(NULL);
    }

    uSolute = (UNIT)oAssocObject(aaArgs[0]);
    uSolvent = (UNIT)oAssocObject(aaArgs[1]);

    if ( iObjectType( oAssocObject(aaArgs[2]) ) == LISTid ) {
        /*
         *  x,y,z box clearances
         */
        if ( iListSize(oAssocObject(aaArgs[2])) != 3 ) {
            VP0(( "%s: Argument #3 must have three values: { x y z } or one.\n",
                                                sCmd ));
            VP0(( sUsage ));
            return(NULL);
        }
        llNumbers = llListLoop((LIST)oAssocObject(aaArgs[2]));
        for ( i=0; i<3; i++ ) {
            aAss = (ASSOC)oListNext(&llNumbers);
            if ( iObjectType(oAssocObject(aAss)) != ODOUBLEid ) {
                VP0(( "%s: Bad value #%d in the third argument.\n", sCmd, i ));
                VP0(( sUsage ));
                return(NULL);
            } 
            daBuffer[i] = dODouble(oAssocObject(aAss));
            if ( daBuffer[i] < 0.0 ) {
                VP0(( "%s: Bad value #%d in the third argument.\n", sCmd, i+1));
                VP0(( sUsage ));
                return(NULL);
            }
        }
    } else {
        daBuffer[0] = dODouble(oAssocObject(aaArgs[2]));
        if ( daBuffer[0] < 0.0 ) {
                VP0(( "%s: Bad box clearance.\n", sCmd ));
                VP0(( sUsage ));
                return(NULL);
        }
        daBuffer[2] = daBuffer[1] = daBuffer[0];
    }

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

    /*
     *  set up solute - centered, & if iso, w/ principal axes aligned
     */
    TurnOffDisplayerUpdates();

    if ( bIsotropic )   /* orient principal axes */
        ToolCenterUnitByRadii( uSolute, TRUE );
    else
        ToolCenterUnitByRadii( uSolute, FALSE );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( (CONTAINER) uSolute );

    /*
     *  do the solvation
     */
    TurnOffDisplayerUpdates();

    iInitialSize = iContainerNumberOfChildren( (CONTAINER) uSolute );

    zToolSolvateAndShell( uSolute, uSolvent, 
                daBuffer[0], daBuffer[1], daBuffer[2], dCloseness,
                NOSHELL, FALSE, TRUE, FALSE, bIsotropic );

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( (CONTAINER) uSolute );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate((CONTAINER) uSolute);


    VP0(( "  Added %d residues.\n", iFinalSize - iInitialSize ));

    return(NULL);
}


/*
 *      oCmd_solvateOct
 *
 *      Author: Bill Ross (1998)
 *
 *      Solvate the UNIT with copies of a box of SOLVENT in a
 *      box with shaved corners.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *              [2] -   LIST with ( x, y, z, d ) distances or
 *                      Or a single ODOUBLE when x,y,z,d distances are the
 *                      same.
 *      Option  [3] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateOct( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oObj;
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0, daBuffer[4];
ASSOC           aAss;
LISTLOOP        llNumbers;
int             i, iInitialSize, iFinalSize;
BOOL            bIsotropic = TRUE;
char            *sCmd = "solvateOct";
char            *sUsage =
                 "usage:  solvateOct <solute> <solvent> <buffer> [aniso] [closeness]\n";

        /*
         *  check args - always need 2 units & cutoff
         */
        if ( !bCmdGoodArguments( sCmd, 3, aaArgs, "u u ln" ) ) {
                VP0(( sUsage ));
                return(NULL);
        }
        if ( iArgCount > 5 ) {
                VP0(( sUsage ));
                return(NULL);
        }
        if ( iArgCount != 3 ) {

                /*
                *  handle possible 'aniso' &/or closeness
                */

                oObj = (OBJEKT)aaArgs[3];
                if ( iObjectType(oObj) != ASSOCid )
                DFATAL(( "unexpected objtype %s\n", sObjectType(oObj) ));

                if ( strcmp("aniso", sAssocName(oObj) ) ) {

                        /* aniso not found;
                           check for old-fashioned iso keyword, and ignore it:  */
                        if ( strcmp("iso", sAssocName(oObj) ) ) {
                                oObj = oAssocObject(oObj);
                        } else {  
                                if ( iArgCount == 5 ) {
                                        oObj = (OBJEKT)aaArgs[4];
                                        oObj = oAssocObject(oObj);
                                } else {
                                        oObj = NULL;
                                }
                        }

                } else {  

                        /* found aniso keyword: */
                        bIsotropic = FALSE;
                        if ( iArgCount == 5 ) {
                                oObj = (OBJEKT)aaArgs[4];
                                oObj = oAssocObject(oObj);
                        } else {
                                oObj = NULL;
                        }
                }

                if ( oObj != NULL ) {
                        if ( iObjectType(oObj) != ODOUBLEid ) {
                                VP0(( sUsage ));
                                return(NULL);
                        }
                        dCloseness = dODouble( oObj );
                        if ( dCloseness <= 0.0 ) {
                                VP0(( sUsage ));
                                return(NULL);
                        }
                }
        }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
        return(NULL);
    }

    uSolute = (UNIT)oAssocObject(aaArgs[0]);
    uSolvent = (UNIT)oAssocObject(aaArgs[1]);

    if ( iObjectType( oAssocObject(aaArgs[2]) ) == LISTid ) {

        if ( bIsotropic ) {
                VP0(( "%s: 'iso' requires a single clearance value\n", sCmd ));
                return(NULL);
        }

        /*
         *  x,y,z,d box clearances
         */
        if ( iListSize(oAssocObject(aaArgs[2])) != 4 ) {
            VP0(( 
                "%s: Argument #3 must have four values: { x y z d } or one.\n",
                                                sCmd ));
            VP0(( sUsage ));
            return(NULL);
        }
        llNumbers = llListLoop((LIST)oAssocObject(aaArgs[2]));
        for ( i=0; i<4; i++ ) {
            aAss = (ASSOC)oListNext(&llNumbers);
            if ( iObjectType(oAssocObject(aAss)) != ODOUBLEid ) {
                VP0(( "%s: Bad value #%d in the third argument.\n", sCmd, i ));
                VP0(( sUsage ));
                return(NULL);
            } 
            daBuffer[i] = dODouble(oAssocObject(aAss));
            if ( daBuffer[i] < 0.0 ) {
                VP0(( "%s: Bad value #%d in the third argument.\n", sCmd, i+1));
                VP0(( sUsage ));
                return(NULL);
            }
        }
    } else {
        daBuffer[0] = dODouble(oAssocObject(aaArgs[2]));
        if ( daBuffer[0] < 0.0 ) {
                VP0(( "%s: Bad box clearance.\n", sCmd ));
                VP0(( sUsage ));
                return(NULL);
        }
        daBuffer[3] = daBuffer[2] = daBuffer[1] = daBuffer[0];
    }

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

    /*
     *  center solute by vdw
     */
    TurnOffDisplayerUpdates();
    if ( bIsotropic )   /* orient principal axes */
        ToolCenterUnitByRadii( uSolute, TRUE );
    else
        ToolCenterUnitByRadii( uSolute, FALSE );
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( (CONTAINER) uSolute );

                /* adjust box for diagonal clearance if necc */

    ToolOctBoxCheck( uSolute, daBuffer, TRUE );

                /* Solvate the solute */

    TurnOffDisplayerUpdates();

    iInitialSize = iContainerNumberOfChildren( (CONTAINER) uSolute );
    zToolSolvateAndShell( uSolute, uSolvent, 
                daBuffer[0], daBuffer[1], daBuffer[2], 
                dCloseness, NOSHELL, FALSE, TRUE, TRUE, bIsotropic );

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( (CONTAINER) uSolute );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate((CONTAINER) uSolute);


    VP0(( "  Added %d residues.\n", iFinalSize - iInitialSize ));

    return(NULL);
}


/*
 *      oCmd_solvateDontClip
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Solvate the UNIT [0] within copies of a box of SOLVENT.
 *      Do not clip the large SOLVENT box that is created around
 *      the solute by duplicating the SOLVENT box.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *              [2] -   LIST with ( x, y, z ) distances or
 *                      xyz, the closest the wall of the solvent box
 *                      can come to the smallest box which contains
 *                      the entire unit that is centered on the
 *                      origin.  This is called the BufferZone.
 *                      Or a single ODOUBLE when x,y,z distances are the
 *                      same.
 *      Option  [3] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateDontClip( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0, daBuffer[3];
ASSOC           aAss;
LISTLOOP        llNumbers;
int             i, iInitialSize, iFinalSize;
BOOL            bIsotropic = FALSE;
char            *sCmd = "solvateDontClip";

    if ( iArgCount == 3 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u ln" ) ) {
            VP0(( 
         "usage:  solvateDontClip <solute> <solvent> <buffer> [closeness]\n" ));
            return(NULL);
        }
        dCloseness = 1.0;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u ln n" ) ) {
            VP0(( 
         "usage:  solvateDontClip <solute> <solvent> <buffer> [closeness]\n" ));
            return(NULL);
        }
        dCloseness = dODouble( oAssocObject(aaArgs[3]) );
    }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
        return(NULL);
    }

    uSolute = (UNIT)oAssocObject(aaArgs[0]);
    uSolvent = (UNIT)oAssocObject(aaArgs[1]);

    if ( iObjectType( oAssocObject(aaArgs[2]) ) == LISTid ) {
        if ( iListSize(oAssocObject(aaArgs[2])) != 3 ) {
            VP0(( "%s: Argument #3 must have three values: { x y z } or one.\n",
                                        sCmd ));
            VP0(( 
         "usage:  solvateDontClip <solute> <solvent> <buffer> [closeness]\n" ));
            return(NULL);
        }
        llNumbers = llListLoop((LIST)oAssocObject(aaArgs[2]));
        for ( i=0; i<3; i++ ) {
            aAss = (ASSOC)oListNext(&llNumbers);
            if ( iObjectType(oAssocObject(aAss)) != ODOUBLEid ) {
                VP0(( "%s: Bad value #%d in the third argument.\n", sCmd, i ));
                VP0(( 
         "usage:  solvateDontClip <solute> <solvent> <buffer> [closeness]\n" ));
                return(NULL);
            } else daBuffer[i] = dODouble(oAssocObject(aAss));
        }
    } else {
        daBuffer[0] = dODouble(oAssocObject(aaArgs[2]));
        daBuffer[2] = daBuffer[1] = daBuffer[0];
    }

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

        /* Orient the principle axis of the solute along the */
        /* coordinate axis, setting a box  */

    TurnOffDisplayerUpdates();
    ToolCenterUnitByRadii( uSolute, FALSE );
    TurnOnDisplayerUpdates();

    ContainerDisplayerUpdate((CONTAINER) uSolute);

    TurnOffDisplayerUpdates();
    iInitialSize = iContainerNumberOfChildren( (CONTAINER) uSolute );
    zToolSolvateAndShell( uSolute, uSolvent, 
                daBuffer[0], daBuffer[1], daBuffer[2], 
                dCloseness, NOSHELL, FALSE, FALSE, FALSE, bIsotropic );

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( (CONTAINER) uSolute );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate((CONTAINER) uSolute);

    VP0(( "  Added %d residues.\n", iFinalSize - iInitialSize ));
    
    return(NULL);
}





/*
 *      oCmd_solvateShell
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Solvate the UNIT [0] within copies of a box of SOLVENT.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *              [2] -   ODOUBLE with thickness parameter.
 *      Option  [3] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateShell( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0, dThickness;
int             iInitialSize, iFinalSize;
char            *sCmd = "solvateShell";
char            *sUsage =
        "usage:  solvateShell <solute> <solvent> <buffer> [closeness]\n";

    if ( iArgCount == 3 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" ) ) {
            VP0(( sUsage ));
            return(NULL);
        }
        dCloseness = 1.0;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n n" ) ) {
            VP0(( sUsage ));
            return(NULL);
        }
        dCloseness = dODouble( oAssocObject(aaArgs[3]) );
    }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
        return(NULL);
    }

    uSolute = (UNIT)oAssocObject(aaArgs[0]);
    uSolvent = (UNIT)oAssocObject(aaArgs[1]);

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

        /* Orient the principle axis of the solute along the */
        /* coordinate axis, setting a box */

    TurnOffDisplayerUpdates();
    ToolCenterUnitByRadii( uSolute, FALSE );
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate((CONTAINER) uSolute);

    dThickness = dODouble(oAssocObject(aaArgs[2]));
    if ( dThickness <= 0.0 ) {
        VP0(( "%s: bad thickness\n", sCmd ));
        VP0(( sUsage ));
        return( NULL );
    }

    TurnOffDisplayerUpdates();
    iInitialSize = iContainerNumberOfChildren( (CONTAINER) uSolute );
    zToolSolvateAndShell( uSolute, uSolvent, 
                dThickness, dThickness, dThickness,
                dCloseness, dThickness, TRUE, TRUE, FALSE, FALSE );

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( (CONTAINER) uSolute );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate((CONTAINER) uSolute);

    VP0(( " Added %d residues.\n", iFinalSize - iInitialSize ));

    
    return(NULL);
}








/*
 *      oCmd_verbosity
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the verbosity to the value passed.
 *
 *      Arguments:
 *              [0] -   Verbosity level 0-2.
 */
OBJEKT
oCmd_verbosity( int iArgCount, ASSOC aaArgs[] )
{
int     iVerb;

    if ( !bCmdGoodArguments( "verbosity", iArgCount, aaArgs, "n" ) ) {
        VP0(( "usage:  verbosity <level>\n" ));
        return(NULL);
    }

    iVerb = (int)dODouble(oAssocObject(aaArgs[0]));
    VerbositySet(iVerb);
    GrMainResult.iCommand     = CVERBOSITY;
    GrMainResult.sVariable[0] = (char)iVerb;

    VP2(( "Verbosity level: %d\n", iVerb ));
    return(NULL);
}




/*
 *      oCmd_logFile
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Define the LogFile.
 *
 *      Arguments:
 *              [0] -   OSTRING with Log File name.
 */
OBJEKT
oCmd_logFile( int iArgCount, ASSOC aaArgs[] )
{
STRING          sFilename;
FILE            *fLog;

    if ( !bCmdGoodArguments( "logFile", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  logFile <filename>\n" ));
        return(NULL);
    }

    strcpy( sFilename, sOString(oAssocObject(aaArgs[0])) );

                /* Close the old log file if it exists */
                
    if ( fVerbosityLogFile() != NULL ) fclose( fVerbosityLogFile() );
    VerbositySetLogFile( NULL );

                /* Open the new log file */

    fLog = FOPENCOMPLAIN( sFilename, "a" );
    if ( fLog != NULL ) {
        time_t  t = time( (time_t *)0 );

        fprintf(fLog, "log started: %s\n", ctime(&t) );
        VerbositySetLogFile( fLog );
        VP0(( "Log file: %s\n", GsBasicsFullName ));
    }
    return(NULL);
}




/*
 *      oCmd_combine
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Combine the contents of several UNITs.
 *      Return a UNIT the combined contents.
 *
 *      Arguments:
 *              [0]     A LIST of units.
 */
OBJEKT
oCmd_combine( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llElements;
UNIT            uCombined, uCurrent;
ASSOC           aAss;
LOOP            lTemp;
RESIDUE         rRes;
OBJEKT          oObj;
char            *sCmd = "combine";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "l" ) ) {
        VP0(( "usage:  <variable> = combine <LIST>\n" ));
        return(NULL);
    }
    llElements = llListLoop( (LIST)oAssocObject(aaArgs[0]) );

                /* Get the first element from the list */

    uCombined = NULL;
    while ( (aAss = (ASSOC)oListNext(&llElements)) ) {
        oObj = oAssocObject(aAss);
        if ( iObjectType( oObj ) != UNITid ) {
                VP0(( "%s: %s is type %s\n", sCmd, 
                                sAssocName(aAss), sObjectType(oObj) ));
                VP0(( "  expected UNIT\n" ));
                continue;
        }
        /*
        **      objekt is a unit, so make a copy
        */
        uCurrent = (UNIT)oCopy( oObj );
        VP1(( "  Sequence: %s\n", sContainerName((CONTAINER) uCurrent) ));
        if ( uCombined == NULL ) {
                MESSAGE(( "Copying the first UNIT\n" ));
                uCombined = uCurrent;
        } else {
                MESSAGE(( "Copied a subsequent UNIT\n" ));
                MESSAGE(( "Joining two UNITS deleting the second\n" ));
                UnitJoin( uCombined, uCurrent );
        }
    }
    if ( uCombined == NULL ) {
        VP0(( "No UNITS, so no combine performed\n" ));
        return(NULL);
    }

                /* Define PDB sequence */

    lTemp = lLoop( (OBJEKT)uCombined, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lTemp)) ) {
        ResidueSetPdbSequence( rRes, iContainerSequence((CONTAINER) rRes) );
    }

    return((OBJEKT)uCombined);
}




/*
 *      oCmd_source
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Execute a stream of commands from a file.
 *
 *      NOTE: This function MODIFIES the global variable GfCurrentInput.
 *
 *      Arguments:
 *              [0]     OSTRING containing the file.
 */
OBJEKT
oCmd_source( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fCmds;

    if ( !bCmdGoodArguments( "source", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  source <filename>\n" ));
        return(NULL);
    }

    fCmds = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );

    if ( fCmds != NULL ) {
        if ( bINPUTMAXDEPTH() ) {
            VP0(( "Source commands are nested too deep!\n" ));

        } else {
                /* Push the file onto the input file stack.  The main */
                /* parsing routine will continue parsing until there are */
                /* no more files on the input file stack. */
                VP0(( "----- Source: %s\n", GsBasicsFullName ));
                INPUTPUSHFILE( fCmds );
                VP0(( "----- Source of %s done\n", GsBasicsFullName ));
        }
    }

    return(NULL);
}



        
/*
 *      oCmd_check
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Check the OBJEKT to see if it is ready to have calculations
 *      run on it.  Various tests are run on the OBJEKT like,
 *      do the atoms have types defined, are all the RESIDUES known, etc.
 *
 *      Arguments:
 *              [0]     UNIT, MOLECULE, RESIDUE, ATOM to check.
 *      Davids Changes
 *      Option  [1]     optional PARMSET to place missing parms into.
 *
 */
OBJEKT
oCmd_check( int iArgCount, ASSOC aaArgs[] )
{
int             iErrors, iWarnings;
UNIT            uUnit;
PARMSET         psParmSet;
char            *sCmd = "check";

    switch ( iArgCount ) {
        case 1: 
            if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umra" )) {
                VP0(( "usage:  check <unit> [parmset]\n" ));
                return(NULL);
            }
            break;
        case 2: 
            if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u ps" )) {
                VP0(( "usage:  check <unit> [parmset]\n" ));
                return(NULL);
            }
            break;
        default: 
            VP0(( "%s: A maximum of two parameters are acceptable:\n", sCmd ));
            VP0(( "     a UNIT and an (optional) PARMSET.\n" ));
            return( NULL );
            break;
    } 
                /* Run the checks on the objects */

    iErrors = 0;
    iWarnings = 0;
    VP0(( "Checking '%s'....\n", sAssocName( aaArgs[0] )));
    ContainerCheck( (CONTAINER) oAssocObject( aaArgs[0] ), &iErrors, &iWarnings );

                        /* Look for close contacts */
    iWarnings += iToolDistanceSearch( (CONTAINER)oAssocObject(aaArgs[0]), 1.5,
                                        TRUE, DISTANCE_SEARCH_PRINT_WARNINGS );

    if ( iObjectType( oAssocObject( aaArgs[0] )) == UNITid ) {
        VP0(( "Checking parameters for unit '%s'.\n", sAssocName( aaArgs[0] ))); 
        uUnit = (UNIT)oAssocObject( aaArgs[0] );
        if ( iArgCount == 2 ) {
            if ( iObjectType(oAssocObject( aaArgs[1] )) == OSTRINGid ) {
                VP0(( "Creating empty parmset %s\n", sAssocName( aaArgs[1] ) ));
                psParmSet = (PARMSET)oCreate(PARMSETid);
            } else {
                psParmSet = (PARMSET)oAssocObject( aaArgs[1] );
            }
            UnitCheckForParms( uUnit, GplAllParameters, psParmSet );
            if ( iObjectType(oAssocObject( aaArgs[1] )) == OSTRINGid )
                Destroy((OBJEKT *) &psParmSet );
        } else {
            UnitCheckForParms( uUnit, GplAllParameters, NULL );
        }       
    }
    if ( iErrors || iWarnings )
        VP0(( "%s:  ", sCmd ));
    if ( iErrors )
        VP0(( "Errors:  %d   ", iErrors ));
    if ( iWarnings )
        VP0(( "Warnings: %d\n", iWarnings ));
    
    if ( iErrors == 0 ) 
        VP0(( "%s is OK.\n", sObjectType( oAssocObject( aaArgs[0] ))));
    
    return(NULL);
}



/*
 *      oCmd_set
 *
 *      Author: Christian Schafmeister (1991)
 *      Modified:  David A. Rivkin (1-22-93)
 *              Added support for the "Defaults" table.
 *
 *      Set attributes of an object or set a default value.
 *
 *      Arguments:
 *              [0] -  OBJEKT, LIST or STRING whose attribute is to be modified.
 *              [1] -  OSTRING with attribute to modify.
 *              [2] -  OBJEKT new value of attribute.
 */


static void
setUsage()
{
        VP0(( "usage:  set <container> <parameter> <object>\n" ));
        VP0(( "   or:  set default <parameter> <value>\n" ));
}

OBJEKT
oCmd_set( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llElements;
ASSOC           aAssoc;
char            *sString;
char             sStrings[100];
char            *sCmd = "set";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "lumras s *" ) ) {
        setUsage();
        return(NULL);
    }

    if ( iObjectType(oAssocObject(aaArgs[0])) == OSTRINGid ) {
        /*
         *  setting an 'environmental' default;
         *      catch user attempt to reference non-existent unit
         */
        sString = sOString(oAssocObject(aaArgs[0]));
        if (strchr(sString, '.')) {
                VP0(( "%s: not a container (e.g. residue)\n", sString ));
                setUsage();
                return(NULL);
        }
        /*
         *  we've got to have 3 args due to hardwired arg checker, so
         *      can't just 'set var value': 'set default var value' where
         *      'default' just occupies space
         */
        StringLower( sString );
        if ( strcmp( sString, "default" )) {
                VP0(( "%s: expected 'default'\n", sString ));
                setUsage();
                return(NULL);
        }

        /*
        **  handle the default 
        */
        if ( iObjectType(oAssocObject(aaArgs[1])) != OSTRINGid ) {
                VP0(( "%s: expected 'default <type> <value>'\n", sString ));
                setUsage();
                return(NULL);
        }

        sString = sOString(oAssocObject(aaArgs[1]));
        StringLower( sString );
        strcpy( sStrings, sString );
        sStrings[4] = 0;

        if ( !strcmp( sString, "pdbwritecharges" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.pdbwritecharges = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.pdbwritecharges = 0;
                } else {
                        VP0(("Set pdbwritecharges: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "nocenter" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.nocenter = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.nocenter = 0;
                } else {
                        VP0(("Set nocenter: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        }
        if ( !strcmp( sString, "reorder_residues" )) {
          
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.reorder_residues = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.reorder_residues = 0;
                } else {
                        VP0(("Set reorder_residues: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        }  
        if ( !strcmp( sString, "oldprmtopformat" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.iOldPrmtopFormat = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.iOldPrmtopFormat = 0;
                } else {
                        VP0(("Set OldPrmtopFormat: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "gibbs" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.iGibbs = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.iGibbs = 0;
                } else {
                        VP0(("Set Gibbs: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "useresids" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.iUseResIds = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.iUseResIds = 0;
                } else {
                        VP0(("Set UseResIds: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "charmm" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.iCharmm = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.iCharmm = 0;
                } else {
                        VP0(("Set Charmm: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "flexiblewater" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.iFlexibleWater = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.iFlexibleWater = 0;
                } else {
                        VP0(("Set FlexibleWater: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "deleteextrapointangles" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "on" ) ) {
                        GDefaults.iDeleteExtraPointAngles = 1;
                } else if ( !strcmp( sString, "off" ) ) {
                        GDefaults.iDeleteExtraPointAngles = 0;
                } else {
                        VP0(("Set DeleteExtraPointAngles: must be 'on' or 'off'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "pbradii" )) {
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );
                if ( !strcmp( sString, "bondi" ) ) {
                        GDefaults.iGBparm = 0;
                        VP0(("Using Bondi radii\n"));
                } else if ( !strcmp( sString, "amber6" ) ) {
                        GDefaults.iGBparm = 1;
                        VP0(("Using amber6 modified Bondi radii\n"));
                } else if ( !strcmp( sString, "mbondi" ) ) {
                        GDefaults.iGBparm = 2;
                        VP0(("Using modified Bondi radii\n"));
#if 0
                } else if ( !strcmp( sString, "pbamber" ) ) {
                        GDefaults.iGBparm = 3;
                        VP0(("Using radii optimized for Amber charges by Huo and Kollman\n"));
#endif
                } else if ( !strcmp( sString, "mbondi2" ) ) {
                        GDefaults.iGBparm = 6;
                        VP0(("Using H(N)-modified Bondi radii\n"));
                } else if ( !strcmp( sString, "parse")){
                       GDefaults.iGBparm = 7;
                       VP0(("Using PARSE radii\n"));
                } else if ( !strcmp( sString, "mbondi3")){
                        GDefaults.iGBparm = 8;
                        VP0(("Using ArgH and AspGluO modified Bondi2 radii\n"));
                } else {
                        VP0(("pbradii option must be 'bondi', 'amber6', 'mbondi', 'mbondi2', 'mbondi3', 'pbamber', or 'parse'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "searchdistance" )) {
                GDefaults.dDSearchDistance = dODouble(oAssocObject(aaArgs[2]));
                return(NULL);
        } 
        if ( !strcmp( sString, "gridspace" )) {
                GDefaults.dGridSpace = dODouble(oAssocObject(aaArgs[2]));
                if (GDefaults.dGridSpace < 0 ) {
                    GDefaults.dGridSpace = 1.0;
                    VP0(( "Must be greater than 0; resetting to %5.2f\n",
                                GDefaults.dGridSpace));
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "shellextent" )) {
                GDefaults.dShellExtent = dODouble(oAssocObject(aaArgs[2]));
                if (GDefaults.dShellExtent < 0 ) {
                    GDefaults.dShellExtent = 4.0;
                    VP0(( "%s: Shell extent must be greater than 0;\n", sCmd));
                    VP0(( "     resetting to %5.2f\n", GDefaults.dShellExtent));
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "dipole_damp_factor" )) {
            GDefaults.dDipoleDampFactor =  dODouble(oAssocObject(aaArgs[2]));
            // Add some verification here.
            if (GDefaults.dDipoleDampFactor < 0.0) {
                    VP0(( "Illegal value: Dipole damping factor can not be less than 0.0.\n"));
                    VP0(( "Set dipole_damp_factor to 0.0.\n"));
                    GDefaults.dDipoleDampFactor = 0.0;
                 }
            return(NULL);
        }
        if ( !strcmp( sString, "sceescalefactor")  || !strcmp( sStrings, "scee" )) {
            GDefaults.dSceeScaleFactor =  dODouble(oAssocObject(aaArgs[2]));
            if (GDefaults.dSceeScaleFactor < 0.0) {
                    VP0(( "Illegal value: 1-4 SCEE factor can not be less than 0.0.\n"));
                    VP0(( "Set sceescalefactor to 1.2.\n"));
                    GDefaults.dSceeScaleFactor = 1.2;
                 }
            return(NULL);
        }
        if ( !strcmp( sString, "scnbscalefactor" ) || !strcmp( sStrings, "scnb" )) {
            GDefaults.dScnbScaleFactor =  dODouble(oAssocObject(aaArgs[2]));
            if (GDefaults.dScnbScaleFactor < 0.0) {
                    VP0(( "Illegal value: 1-4 SCNB factor can not be less than 0.0.\n"));
                    VP0(( "Set scnbscalefactor to 2.0.\n"));
                    GDefaults.dScnbScaleFactor = 2.0;
                 }
            return(NULL);
        }
        if ( !strcmp( sString, "cmap" )) {
            sString = sOString(oAssocObject(aaArgs[2]));
            StringLower( sString );
            if ( !strcmp( sString, "on" ) ) {
                GDefaults.iCMAP = 1;
            } else if ( !strcmp( sString, "off" ) ) {
                GDefaults.iCMAP = 0;
            } else {
                VP0(("set default CMAP: must be 'on' or 'off'\n"));
                setUsage();
            }
            return(NULL);
        }
        if ( !strcmp( sString, "ipol" )) {
            int myinteger = dODouble(oAssocObject(aaArgs[2]));
        if ( myinteger < 0 || myinteger > 4 ) {
            VP0(("Only IPOL = 0 to 4 is supported. Fall back to IPOL = 0.\n"));
                myinteger = 0;
        } 
            if ( GDefaults.iIPOLset > 0 ) {
            VP0(("IPOL has already been set to %i in frcmod/parm.dat.\n", GDefaults.iIPOL));
            VP0(("Please change the setting in frcmod/parm.dat.\n"));
            } else {
                GDefaults.iIPOL = myinteger;
            }
            return(NULL);
    }
        if ( !strcmp( sString, "pdbloadsequential" )) {
                if ( iObjectType(oAssocObject(aaArgs[2])) != OSTRINGid ) {
                        VP0(("expected 'true' or 'false'\n"));
                        setUsage();
                        return(NULL);
                }
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );

                if (!strcmp(sString, "true"))
                        GDefaults.iPdbLoadSequential = 1;
                else if (!strcmp(sString, "false"))
                        GDefaults.iPdbLoadSequential = 0;
                else {
                        VP0(("expected 'true' or 'false'\n"));
                        setUsage();
                }
                return(NULL);
        }
        if ( !strcmp( sString, "dielectric" )) {
                if ( iObjectType(oAssocObject(aaArgs[2])) == OSTRINGid ) {
                    sString = sOString(oAssocObject(aaArgs[2]));
                    StringLower( sString );
                    if ( !strcmp( sString, "constant" ) ) {
                        GDefaults.iDielectricFlag = DIEL_R;
                    } else if ( !strcmp( sString, "distance" ) ) {
                        GDefaults.iDielectricFlag = DIEL_R2;
                    } else {
                        VP0(("%s: expected 'distance' or 'constant'\n", sCmd ));
                    }
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "residueimpropers" )) {

                if ( iObjectType(oAssocObject(aaArgs[2])) != OSTRINGid ) {
                        VP0(("expected 'true' or 'false'\n"));
                        setUsage();
                        return(NULL);
                }
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );

                if (!strcmp(sString, "true"))
                        GDefaults.iResidueImpropers = 1;
                else if (!strcmp(sString, "false"))
                        GDefaults.iResidueImpropers = 0;
                else {
                        VP0(("expected 'true' or 'false'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        if ( !strcmp( sString, "deleteextrapointangles" )) {

                if ( iObjectType(oAssocObject(aaArgs[2])) != OSTRINGid ) {
                        VP0(("expected 'true' or 'false'\n"));
                        setUsage();
                        return(NULL);
                }
                sString = sOString(oAssocObject(aaArgs[2]));
                StringLower( sString );

                if (!strcmp(sString, "true"))
                        GDefaults.iDeleteExtraPointAngles = 1;
                else if (!strcmp(sString, "false"))
                        GDefaults.iDeleteExtraPointAngles = 0;
                else {
                        VP0(("expected 'true' or 'false'\n"));
                        setUsage();
                }
                return(NULL);
        } 
        VP0(( "can't parse %s\n", sString ));
        setUsage();
        return(NULL);

    } 
    if ( iObjectType(oAssocObject(aaArgs[0])) == LISTid ) {
        DisplayerAccumulateUpdates();
        llElements = llListLoop((LIST)oAssocObject(aaArgs[0]));
        while ( (aAssoc = (ASSOC)oListNext(&llElements)) ) {
            MESSAGE(( "Setting attribute for: %s\n", sAssocName(aAssoc) ));
            if ( bObjectInClass( oAssocObject(aAssoc), CONTAINERid ) ) {
                ContainerSetAttribute( (CONTAINER) oAssocObject(aAssoc),
                          sOString(oAssocObject(aaArgs[1])),
                          oAssocObject(aaArgs[2]) );
            } else {
                VP0(( "%s: Cannot set attribute for %s - not a 'container'\n",
                        sCmd, sAssocName(aAssoc) ));
                VP0(( "\ttype %s\n", sObjectType(oAssocObject(aAssoc)) ));
            }
        }
        DisplayerReleaseUpdates();
        return(NULL);
    } 
    DisplayerAccumulateUpdates();
    ContainerSetAttribute( (CONTAINER) oAssocObject(aaArgs[0]),
                        sOString(oAssocObject(aaArgs[1])),
                        oAssocObject(aaArgs[2]) );
    
    DisplayerReleaseUpdates();
    
    return(NULL);
}

/*
 *      oCmd_setBox
 *
 *      Author: Bill Ross (1996)
 *
 *      Arguments:
 *              [0] -   OBJEKT to center/set vdw box.
 *              [1] -   'vdw' or 'centers'
 *              [2] -   optional offset (1 number or list of 3)
 */
OBJEKT
oCmd_setBox( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
double  dX, dY, dZ, daBuffer[3];
int     i;
BOOL    bUsage;
char    *sCmd = "setBox";
char    *sUsage = "usage:  setBox <unit> vdw|centers [ clearance | <clearance_xyz_list> ]\n";
char    *sOpt;
BOOL    bVdw;

    if ( iArgCount == 2 ) 
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s" );
    else 
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s ln" );

    if ( bUsage ) {
        VP0(( sUsage ));
        return(NULL);
    }

    uUnit = (UNIT)oAssocObject(aaArgs[0]);
    sOpt = (char *)sOString(oAssocObject(aaArgs[1]));
    if ( !strcmp(sOpt, "vdw")  ) {
        bVdw = TRUE;
    } else if ( !strcmp(sOpt, "centers")  ) {
        bVdw = FALSE;
    } else {
        VP0(( "%s: Expected 'vdw' or 'centers' for second argument.\n", sCmd ));
        VP0(( sUsage ));
        return(NULL);
    }

    if ( bUnitBoxOct(uUnit) ) {
        VP0(( " removing previous octbox..\n" ));
        UnitResetFlags(uUnit, UNITBOXOCT);
        UnitResetFlags(uUnit, UNITUSEBOUNDINGBOX);
    } else if ( bUnitUseBox( uUnit ) ) {
        VP0(( " removing previous box..\n" ));
        UnitResetFlags(uUnit, UNITUSEBOUNDINGBOX);
    }

    /*
     *  get any box offset
     */
    if ( iArgCount == 3 ) {
        if ( iObjectType( oAssocObject(aaArgs[2]) ) == LISTid ) {
            LISTLOOP        llNumbers;

            /*
             *  x,y,z box clearances
             */
            if ( iListSize(oAssocObject(aaArgs[2])) != 3 ) {
                VP0(( "%s: Expected 3 non-negative floating point numbers"
                        " { x y z } for third argument.\n",
                        sCmd ));
                VP0(( sUsage ));
                return(NULL);
            }
            llNumbers = llListLoop((LIST)oAssocObject(aaArgs[2]));
            for ( i=0; i<3; i++ ) {
                ASSOC   aAss = (ASSOC)oListNext(&llNumbers);
                if ( iObjectType(oAssocObject(aAss)) != ODOUBLEid ) {
                    VP0(( "%s: Bad value #%d in third argument.\n", sCmd, i ));
                    VP0(( "%s: Expected 3 non-negative floating point numbers"
                            " { x y z } for third argument.\n",
                            sCmd ));
                    VP0(( sUsage ));
                    return(NULL);
                } 
                daBuffer[i] = 2.0 * dODouble(oAssocObject(aAss));
                if ( daBuffer[i] < 0.0 ) {
                    VP0(( "%s: Bad value #%d in third argument.\n", sCmd, i ));
                    VP0(( "%s: Expected 3 non-negative floating point numbers"
                            " { x y z } for third argument.\n",
                            sCmd ));
                    VP0(( sUsage ));
                    return(NULL);
                }
            }
        } else {
            daBuffer[0] = 2.0 * dODouble(oAssocObject(aaArgs[2]));
            if ( daBuffer[0] < 0.0 ) {
                VP0(( "%s: Bad box clearance.\n", sCmd ));
                VP0(( "%s: Expected a non-negative floating point number"
                        " for third argument.\n",
                        sCmd ));
                VP0(( sUsage ));
                return(NULL);
            }
            daBuffer[2] = daBuffer[1] = daBuffer[0];
        }
    }

    DisplayerAccumulateUpdates();
    if( GDefaults.nocenter == 0 ){
        if (bVdw) {
            ToolCenterUnitByRadii( uUnit, FALSE );
        } else {
            VECTOR vOrigin;
            VectorDef( &vOrigin, 0.0, 0.0, 0.0 );
            ContainerCenterAt( (CONTAINER)uUnit, vOrigin );
            ToolSetUnitBoxByCenters( uUnit );
        }
    }
    DisplayerReleaseUpdates();

    UnitSetUseBox( uUnit, TRUE );
    UnitSetBeta( uUnit, 90.0*DEGTORAD );
    UnitGetBox( uUnit, &dX, &dY, &dZ );
    if ( iArgCount == 3 ) {
        dX += daBuffer[0];
        dY += daBuffer[1];
        dZ += daBuffer[2];
        UnitSetBox( uUnit, dX, dY, dZ );
    }
    VP0(( "Box dimensions:  %f %f %f\n", dX, dY, dZ ));
    return(NULL);
}
int
PrntOnOff(char* s1, char* s2, char* info, int iswitch, int ALL)
{
int l;
        l = strlen(s1);
        if ( !strncmp( s1, s2, l) || ALL) {
                VP0(("%25s :    ", info));
                if ( iswitch )
                        VP0(("on/true\n"));
                else
                        VP0(("off/false\n"));
return(1); } else return (0);

}
int
PrntOpt(char* s1, char* s2, char* info, char* opt[], int ind, int ALL)
{
int l;
        l = strlen(s1);
        if ( !strncmp( s1, s2, l) || ALL) {
                VP0(("%25s :    %s\n", info,opt[ind]));
return(1); } else return (0);
}

int
PrntInt(char* s1, char* s2, char* info, int i, int ALL)
{
int l;
        l = strlen(s1);
        if ( !strncmp( s1, s2, l) || ALL) {
                VP0(("%25s :    %i\n", info,i));
return(1); } else return (0);
}

int
PrntReal(char* s1, char* s2, char* info, double f, int ALL)
{
int l;
        l = strlen(s1);
        if ( !strncmp( s1, s2, l) || ALL) {
                VP0(("%25s :    %lf\n", info,f));
return(1); } else return (0);
}

/*
 * Show default variables
 *      Author: Yong Duan 
 *      (adapted from set default)
 */
OBJEKT
oCmd_showDefault( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llElements;
ASSOC           aAssoc;
char            *sString;
char            sStrings[]="all";
int             ALL, ind, found;

        ALL = 0;
        found = 0;
        if ( !iArgCount ) sString = sStrings;
        else sString = sOString(oAssocObject(aaArgs[0]));
        StringLower( sString );

        if ( !strcmp(sString, "all") || !strcmp(sString, "*" )) ALL = 1;

        if ( iArgCount) sString = sOString(oAssocObject(aaArgs[0]));
        StringLower( sString );

        found += PrntOnOff(sString, "pdbwritecharges", "PdbWriteCharges", GDefaults.pdbwritecharges, ALL);
        found += PrntOnOff(sString, "oldprmtopformat", "OldPrmtopFormat", GDefaults.iOldPrmtopFormat, ALL);
        found += PrntOnOff(sString, "gibbs", "Gibbs", GDefaults.iGibbs, ALL);
        found += PrntOnOff(sString, "useresids", "UseResIds", GDefaults.iUseResIds, ALL);
        found += PrntOnOff(sString, "charmm", "Charmm", GDefaults.iCharmm, ALL);
        found += PrntOnOff(sString, "flexiblewater", "FlexibleWater", GDefaults.iFlexibleWater, ALL);
        found += PrntOnOff(sString, "deleteextrapointangles", "DeleteExtraPointAngles", GDefaults.iDeleteExtraPointAngles, ALL);
        char* s[7]={"Bondi radii","Amber6 modified Bondi radii","Modified Bondi radii","","","","H(N)-modified Bondi radii"};
        found += PrntOpt(sString, "pbradii", "PB Radii", s,GDefaults.iGBparm, ALL);
        found += PrntReal(sString, "searchdistance", "SearchDistance", GDefaults.dDSearchDistance, ALL);
        found += PrntReal(sString, "gridspace", "GridSpace", GDefaults.dGridSpace, ALL);
        found += PrntReal(sString, "shellextent", "ShellExtent", GDefaults.dShellExtent, ALL);
        found += PrntReal(sString, "dipole_damp_factor", "Dipole Damping Factor", GDefaults.dDipoleDampFactor, ALL);
        found += PrntReal(sString, "sceescalefactor", "SCEE 1-4 Scale Factor", GDefaults.dSceeScaleFactor, ALL);
        found += PrntReal(sString, "scnbscalefactor", "SCNB 1-4 Scale Factor", GDefaults.dScnbScaleFactor, ALL);
        found += PrntOnOff(sString, "cmap", "CMAP", GDefaults.iCMAP, ALL);
        found += PrntInt(sString, "ipol", "IPOL", GDefaults.iIPOL, ALL);
        found += PrntOnOff(sString, "pdbloadsequential", "PdbLoadSequential", GDefaults.iPdbLoadSequential, ALL);
        char* s1[3]={"Constant","Distance","Undefined"};
        if (GDefaults.iDielectricFlag == DIEL_R ) ind = 0;
        else if (GDefaults.iDielectricFlag == DIEL_R2 ) ind = 1;
        else ind = 2;
        found += PrntOpt(sString, "dielectric", "Dielectric", s1,ind, ALL);
        found += PrntOnOff(sString, "residueimpropers", "ResidueImpropers", GDefaults.iResidueImpropers, ALL);
        found += PrntOnOff(sString, "deleteextrapointangles", "DeleteExtraPointAngles", GDefaults.iDeleteExtraPointAngles, ALL);
        if (!found) VP0(( "can't parse %s\n", sString ));

    return(NULL);
}

/*
 *      oCmd_charge
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   uNIT whose charge should be calculated.
 *
 */
OBJEKT
oCmd_charge( int iArgCount, ASSOC aaArgs[] )
{
CONTAINER               cCont;
double                  dCharge, dPertCharge;

    if ( !bCmdGoodArguments( "charge", iArgCount, aaArgs, "umral" ) ) {
        VP0(( "usage:  charge <object>\n" ));
        return(NULL);
    }

    cCont = (CONTAINER)oAssocObject(aaArgs[0]);

    ContainerTotalCharge( cCont, &dCharge, &dPertCharge );
    VP0(( "Total unperturbed charge: %10.6lf\n", dCharge ));
    /*dPertCharge is now the delta of the charge and not the actual perturbed charge*/
    VP0(( "Total perturbed charge:   %10.6lf\n", (dCharge+dPertCharge) ));

    return(NULL);
}



/*
 *      oCmd_saveAmberParm
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
void
sapUsage()
{
  VP0(("usage:  saveAmberParm <unit> <topologyfile> <coordfile> \n" ));
}
OBJEKT
oCmd_saveAmberParm( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
FILE    *fOut;
char    *sCmd = "saveAmberParm";

        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
                sapUsage();
                return(NULL);
        }

        if ( iParmLibSize(GplAllParameters) == 0 ) {
                VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
                return(NULL);
        }
        uUnit = (UNIT)oAssocObject(aaArgs[0]);
        fOut = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[1])), "w" );
        if ( fOut == NULL ) {
                VP0(( "%s: Could not open file: %s\n", 
                        sCmd, sOString(oAssocObject(aaArgs[1]))));
                return(NULL);
        }
        TurnOffDisplayerUpdates();
        UnitSaveAmberParmFile( uUnit, fOut, sOString(oAssocObject(aaArgs[2])), GplAllParameters, FALSE,
                FALSE, FALSE );
        TurnOnDisplayerUpdates();
    
        fclose( fOut );
        return(NULL);
}

/*
 *      oCmd_saveAmberParmNetCDF
 *
 *      Author: Robin Betz (2013)
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParmNetCDF( int iArgCount, ASSOC aaArgs[] )
{
  UNIT    uUnit;
  FILE    *fOut;
  char    *sCmd = "saveAmberParmNetCDF";
  
  if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
    VP0(("usage: saveAmberParmNetCDF <unit> <topologyfile> <coordfile> \n" ));
    return(NULL);
}

if ( iParmLibSize(GplAllParameters) == 0 ) {
  VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
  return(NULL);
}
uUnit = (UNIT)oAssocObject(aaArgs[0]);
fOut = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[1])), "w" );
if ( fOut == NULL ) {
  VP0(( "%s: Could not open file: %s\n", 
        sCmd, sOString(oAssocObject(aaArgs[1]))));
  return(NULL);
}
TurnOffDisplayerUpdates();
UnitSaveAmberParmFile( uUnit, fOut, sOString(oAssocObject(aaArgs[2])), GplAllParameters, FALSE,
                       FALSE, TRUE );
// UnitSaveAmberNetcdf( uUnit, sOString(oAssocObject(aaArgs[2])) );
TurnOnDisplayerUpdates();

fclose( fOut );
return(NULL);
}


/*
 *      oCmd_saveAmberParmPol
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParmPol( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
FILE    *fOut;
char    *sCmd = "saveAmberParmPol";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
        VP0(( "usage:  %s <unit> <topologyfile> <coordfile>\n", sCmd ));
        return(NULL);
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
        return(NULL);
    }
    uUnit = (UNIT)oAssocObject(aaArgs[0]);
    fOut = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[1])), "w" );
    if ( fOut == NULL ) {
        VP0(( "%s: Could not open file: %s\n", 
                sCmd, sOString(oAssocObject(aaArgs[1]))));
        return(NULL);
    }
    
    TurnOffDisplayerUpdates();
    UnitSaveAmberParmFile( uUnit, fOut, sOString(oAssocObject(aaArgs[2])), GplAllParameters, TRUE, FALSE, FALSE );
    TurnOnDisplayerUpdates();
    
    fclose( fOut );
    return(NULL);
}


/*
 *      oCmd_saveAmberParmPert
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParmPert( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
FILE    *fOut;
char    *sCmd = "saveAmberParmPert";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
        VP0(( 
          "usage:  %s <unit> <topologyfile> <coordfile>\n", sCmd ));
        return(NULL);
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
        return(NULL);
    }
    uUnit = (UNIT)oAssocObject(aaArgs[0]);
    fOut = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[1])), "w" );
    if ( fOut == NULL ) {
        VP0(( "%s: Could not open file: %s\n", 
                        sCmd, sOString(oAssocObject(aaArgs[1]))));
        return(NULL);
    }
    
    TurnOffDisplayerUpdates();
    UnitSaveAmberParmFile( uUnit, fOut, sOString(oAssocObject(aaArgs[1])), GplAllParameters, FALSE, TRUE, FALSE );
    TurnOnDisplayerUpdates();
    
    fclose( fOut );
    return(NULL);
}


/*
 *      oCmd_saveAmberParmPolPert
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParmPolPert( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
FILE    *fOut;
char    *sCmd = "saveAmberParmPert";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
        VP0(( 
          "usage:  %s <unit> <topologyfile> <coordfile>\n", sCmd ));
        return(NULL);
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VP0(( "%s: There are no parameter sets loaded\n", sCmd ));
        return(NULL);
    }
    uUnit = (UNIT)oAssocObject(aaArgs[0]);
    fOut = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[1])), "w" );
    if ( fOut == NULL ) {
        VP0(( "%s: Could not open file: %s\n", 
                        sCmd, sOString(oAssocObject(aaArgs[1]))));
        return(NULL);
    }
    
    TurnOffDisplayerUpdates();
    UnitSaveAmberParmFile( uUnit, fOut, sOString(oAssocObject(aaArgs[2])), GplAllParameters, TRUE, TRUE, FALSE );
    TurnOnDisplayerUpdates();
    
    fclose( fOut );
    return(NULL);
}



OBJEKT
oCmd_saveAmberPrep( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;
FILE            *fOut;
char            *sCmd = "saveAmberPrep";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s" ) ) {
        VP0(( "usage:  saveAmberPrep <unit> <file>\n" ));
        return(NULL);
    }

    uUnit = (UNIT)oAssocObject(aaArgs[0]);

    fOut = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[1])), "w" );
    if ( fOut == NULL ) {
        VP0(( "%s: Could not open file: %s\n", 
                sCmd, sOString(oAssocObject(aaArgs[1]))));
        return(NULL);
    }
    
    TurnOffDisplayerUpdates();
    UnitIOSaveAmberPrep( uUnit, fOut );
    TurnOnDisplayerUpdates();
    
    fclose( fOut );
    VP0(( "  -- Remember to delete unwanted IMPROPER terms!\n" ));
    return(NULL);
}


/*
 *      oCmd_clearVariables
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Erase all variables.
 */
OBJEKT
oCmd_clearVariables( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llVariables;
ASSOC           aAssoc;

    if ( iArgCount == 0 ) {
        VP0(( "Clearing all variables\n" ));
        VariablesDestroy();
        VariablesInit();
    } else {
        if ( !bCmdGoodArguments( "clearVariables", iArgCount, aaArgs, "l" ) ) {
            VP0(( "usage:  clearVariables [LIST]\n" ));
            return(NULL);
        }
        llVariables = llListLoop( (LIST)oAssocObject(aaArgs[0]) );
        while ( (aAssoc = (ASSOC)oListNext(&llVariables)) ) {
            VariableRemove( sAssocName(aAssoc) );
        }
    }
    return(NULL);
}






/*
 *      oCmd_matchVariables
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a list of all variables that match the name.
 *
 *      Arguments:
 *              [0] -   OSTRING, name with '*' and '?' wildcards to match.
 *
 */
OBJEKT
oCmd_matchVariables( int iArgCount, ASSOC aaArgs[] )
{
LIST            lVars;
DICTIONARY      dVariables;
DICTLOOP        dlEntries;
ASSOC           aAssoc;
char            *sCmd = "matchVariables";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  <variable> = matchVariables <string>\n" ));
        return(NULL);
    }
    dVariables = dVariablesDictionary();
    lVars = (LIST)oCreate(LISTid);

    dlEntries = ydlDictionaryLoop(dVariables);
    while ( yPDictionaryNext( dVariables, &dlEntries ) ) {
        if ( bStringMatchPattern( sDictLoopKey(dlEntries),
                                  sOString(oAssocObject(aaArgs[0])) ) ) {
            aAssoc = (ASSOC)oCreate(ASSOCid);
            AssocSetName( aAssoc, sDictLoopKey(dlEntries) );
            AssocSetObject( aAssoc, PDictLoopData(dlEntries) );
            ListAddToEnd( lVars, (OBJEKT)aAssoc );
        }
    }

    return((OBJEKT)lVars);
}




/*
 *      oCmd_edit
 *
 *      Author: Christian Schafmeister (1991)
 *      Modified: David A. Rivkin (1992)
 *
 *      If in a Graphical environment setup the GrMainResult structure
 *      to return the UNIT that the USER wishes to edit.  Otherwise
 *      print a disclaimer that 'edit' does not work in a non-graphical
 *      environment.
 *
 *      This command causes a COPY of the OBJEKT to be edited.
 *      If the OBJEKT that the user wants to edit is NULL then
 *      this command assumes that they want to edit a new
 *      UNIT with a single RESIDUE within it.
 *
 *      Arguments:
 *              [0] -   OBJEKT to edit.
 *      Davids Changes - changed "uz" to "ups" in bCmdGoodArguments call.
 *              Moved and copied the "GrMainResult.oObject = uUnit" to the
 *              first two if statements and added a check and 
 *              GrMainResult.oObject = psParmSet added PARMSET psParmSet
 *              to the local variables.
 */
OBJEKT
oCmd_edit( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;
PARMSET         psParmSet;
RESIDUE         rRes;

    if ( !GbGraphicalEnvironment ) {
        VP0(( "The edit command only works in a graphical environment.\n" ));
        GrMainResult.iCommand = CNONE;
        return(NULL);
    }
    if ( !bCmdGoodArguments( "edit", iArgCount, aaArgs, "ups" ) ) {
        VP0(( "usage:  edit <unit/parmset>\n" ));
        return(NULL);
    }

    GrMainResult.iCommand = CEDIT;
    strcpy( GrMainResult.sVariable, sAssocName(aaArgs[0]) );

                /* Check if the UNIT exists, if it doesn't then */
                /* create a new one with a single RESIDUE in it */

    if ( oAssocObject(aaArgs[0]) == NULL  || 
                        iObjectType(oAssocObject(aaArgs[0])) == OSTRINGid ) {
            VP0(( "Creating a new, empty UNIT \"%s\"\n", 
                        sOString( oAssocObject(aaArgs[0]))));
            uUnit = (UNIT)oCreate(UNITid);
            ContainerSetName( (CONTAINER) uUnit, 
                                sOString( oAssocObject(aaArgs[0])));
            rRes  = (RESIDUE)oCreate(RESIDUEid);
            ContainerSetName( (CONTAINER) rRes, 
                                sOString(oAssocObject(aaArgs[0])));
            ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rRes );
            VariableSet( sOString( oAssocObject(aaArgs[0])), (OBJEKT)uUnit );   /* adds 1 REF */
            GrMainResult.oObject = (OBJEKT) uUnit;

    } else if ( iObjectType(oAssocObject(aaArgs[0])) == UNITid ) {
            uUnit = (UNIT)(oAssocObject(aaArgs[0]));
            GrMainResult.oObject = (OBJEKT) uUnit;

    } else if ( iObjectType(oAssocObject(aaArgs[0])) == PARMSETid) {
            psParmSet = (PARMSET)(oAssocObject(aaArgs[0]));    
            GrMainResult.oObject = (OBJEKT) psParmSet;

    } else {
            VP0(( "Can't edit type %s\n", 
                                sObjectType(oAssocObject(aaArgs[0])) ));
            GrMainResult.iCommand = CNONE;
    }

    return(NULL);
}

        


/*
 *      oCmd_alignAxes
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Align the principle axis of the molecule along the
 *      coordinate axis.
 *
 *      Arguments:
 *              [0] -   UNIT to align.
 *
 */
OBJEKT
oCmd_alignAxes( int iArgCount, ASSOC aaArgs[] )
{

    if ( !bCmdGoodArguments( "alignAxes", iArgCount, aaArgs, "u" ) ) {
        VP0(( "usage:  alignAxes <unit>\n" ));
        return(NULL);
    }

    DisplayerAccumulateUpdates();
    
    ToolOrientPrincipleAxisAlongCoordinateAxis( (UNIT)oAssocObject(aaArgs[0]) );

                /* Instruct the graphics system to reset the viewing */
                /* matrices of the UNIT */

    DisplayerReleaseUpdates();
    
    return(NULL);
}





/*
 *      oCmd_select
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This command sets the ATOMSELECTED flag of atoms 
 *      within the OBJEKT
 *
 *      Arguments:
 *              [0]     - CONTAINER/LIST.
 */

OBJEKT
oCmd_select( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oOver;
LOOP            lAtoms;
ATOM            aAtom;

    if ( !bCmdGoodArguments( "select", iArgCount, aaArgs, "umral" ) ) {
        VP0(( "usage:  select <object>\n" ));
        return(NULL);
    }

    DisplayerAccumulateUpdates();
    
    oOver = oAssocObject(aaArgs[0]);
    lAtoms = lLoop( oOver, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        AtomSetFlags( aAtom, ATOMSELECTED );
    }
    DisplayerReleaseUpdates();
    return(NULL);
}

/*
 *      oCmd_scaleCharges
 *
 *      scale charges; useful for polar setup
 */
OBJEKT
oCmd_scaleCharges( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oOver;
LOOP            lAtoms;
ATOM            aAtom;
double          dScale;

    if ( !bCmdGoodArguments( "scaleCharges", iArgCount, aaArgs, "umral n" ) ) {
        VP0(( "usage:  scaleCharges <object> <scale_factor>\n" ));
        return(NULL);
    }

    dScale = dODouble(oAssocObject(aaArgs[1]));
    if ( dScale <= 0.0 ) {
        VP0(( "scaleCharges: scale_factor must be > 0\n" ));
        VP0(( "usage:  scaleCharges <object> <scale_factor>\n" ));
        return(NULL);
    }

    oOver = oAssocObject(aaArgs[0]);
    lAtoms = lLoop( oOver, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        /* HACK - resetting charge without marking for update */
        dAtomCharge( aAtom ) *= dScale;
        dAtomPertCharge( aAtom ) *= dScale;
    }
    return(NULL);
}


/*
 *      oCmd_deSelect
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This command resets the ATOMSELECTED flag of atoms 
 *      within the OBJEKT
 *
 *      Arguments:
 *              [0]     - CONTAINER/LIST.
 */

OBJEKT
oCmd_deSelect( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oOver;
LOOP            lAtoms;
ATOM            aAtom;

    if ( !bCmdGoodArguments( "deSelect", iArgCount, aaArgs, "umral" ) ) {
        VP0(( "usage:  deSelect <object>\n" ));
        return(NULL);
    }

    DisplayerAccumulateUpdates();
    
    oOver = oAssocObject(aaArgs[0]);
    lAtoms = lLoop( oOver, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        AtomResetFlags( aAtom, ATOMSELECTED );
    }
    
    DisplayerReleaseUpdates();
    
    return(NULL);
}




/*
 *      oCmd_restrainBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add a RESTRAINT bond to the UNIT, make sure that the atoms
 *      are contained within the UNIT.
 *
 *      Arguments:
 *              [0]     - The UNIT to add the RESTRAINT bond to.
 *              [1]     - The first ATOM of the RESTRAINT bond.
 *              [2]     - The second ATOM of the RESTRAINT bond.
 *              [3]     - ODOUBLE: The value of Kr.
 *              [4]     - ODOUBLE: The value of R0.
 */
OBJEKT
oCmd_restrainBond( int iArgCount, ASSOC aaArgs[] )
{
int             i;
UNIT            uUnit;
ATOM            aaAtoms[ATOMSINBOND];
double          dKr, dR0;
RESTRAINT       rRest;
char            *sCmd = "restrainBond";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a n n" ) ) {
        VP0(( "usage:  restrainBond <unit> <a> <b> <force> <length>\n" ));
        return(NULL);
    }

    uUnit  = (UNIT)oAssocObject(aaArgs[0]);
    for ( i=0; i<ATOMSINBOND; i++ )
        aaAtoms[i] = (ATOM)oAssocObject(aaArgs[1+i]);
    dKr    = dODouble(oAssocObject(aaArgs[3]));
    dR0    = dODouble(oAssocObject(aaArgs[4]));

    for ( i=0; i<ATOMSINBOND; i++ ) {
        if ( !bContainerContainedBy( (CONTAINER)aaAtoms[i], (CONTAINER)uUnit ) ) {
            VP0(( "%s: Atom#%d is not part of the UNIT\n", sCmd, i+1 ));
            return(NULL);
        }
    }

                /* Create the RESTRAINT and add it to the UNIT */

    rRest = rRestraintCreate();
    RestraintBondSet( rRest, aaAtoms[0], aaAtoms[1], dKr, dR0 );
    UnitAddRestraint( uUnit, rRest );

    return(NULL);
}





/*
 *      oCmd_restrainAngle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add a RESTRAINT angle to the UNIT, make sure that the atoms
 *      are contained within the UNIT.
 *
 *      Arguments:
 *              [0]     - The UNIT to add the RESTRAINT to.
 *              [1]     - The first ATOM of the RESTRAINT.
 *              [2]     - The second ATOM of the RESTRAINT.
 *              [3]     - The third ATOM of the RESTRAINT.
 *              [4]     - ODOUBLE: The value of Kt.
 *              [5]     - ODOUBLE: The value of T0.
 */
OBJEKT
oCmd_restrainAngle( int iArgCount, ASSOC aaArgs[] )
{
int             i;
UNIT            uUnit;
ATOM            aaAtoms[ATOMSINANGLE];
double          dKt, dT0;
RESTRAINT       rRest;
char            *sCmd = "restrainAngle";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a a n n" ) ) {
        VP0(( "usage:  restrainAngle <unit> <a> <b> <c> <force> <length>\n" ));
        return(NULL);
    }

    uUnit  = (UNIT)oAssocObject(aaArgs[0]);
    for ( i=0; i<ATOMSINANGLE; i++ )
        aaAtoms[i] = (ATOM)oAssocObject(aaArgs[1+i]);
    dKt    = dODouble(oAssocObject(aaArgs[4]));
    dT0    = DEGTORAD * dODouble(oAssocObject(aaArgs[5]));

    for ( i=0; i<ATOMSINANGLE; i++ ) {
        if ( !bContainerContainedBy( (CONTAINER)aaAtoms[i], (CONTAINER)uUnit ) ) {
            VP0(( "%s: Atom#%d is not part of the UNIT\n", sCmd, i+1 ));
            return(NULL);
        }
    }

                /* Create the RESTRAINT and add it to the UNIT */

    rRest = rRestraintCreate();
    RestraintAngleSet( rRest, aaAtoms[0], aaAtoms[1], aaAtoms[2], dKt, dT0 );
    UnitAddRestraint( uUnit, rRest );

    return(NULL);
}





/*
 *      oCmd_restrainTorsion
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add a RESTRAINT torsion to the UNIT, make sure that the atoms
 *      are contained within the UNIT.
 *
 *      Arguments:
 *              [0]     - The UNIT to add the RESTRAINT bond to.
 *              [1]     - The first ATOM of the RESTRAINT bond.
 *              [2]     - The second ATOM of the RESTRAINT bond.
 *              [3]     - The third ATOM of the RESTRAINT bond.
 *              [4]     - The fourth ATOM of the RESTRAINT bond.
 *              [5]     - ODOUBLE: The value of Kp.
 *              [6]     - ODOUBLE: The value of P0.
 *              [7]     - ODOUBLE: The value of N.
 */
OBJEKT
oCmd_restrainTorsion( int iArgCount, ASSOC aaArgs[] )
{
int             i;
UNIT            uUnit;
ATOM            aaAtoms[ATOMSINTORSION];
double          dKp, dP0, dN;
RESTRAINT       rRest;
char            *sCmd = "restrainTorsion";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a a a n n n" ) ) {
        VP0(( 
         "usage:  restrainTorsion <unit> <a> <b> <c> <d> <force> <length>\n" ));
        return(NULL);
    }

    uUnit  = (UNIT)oAssocObject(aaArgs[0]);
    for ( i=0; i<ATOMSINTORSION; i++ )
        aaAtoms[i] = (ATOM)oAssocObject(aaArgs[1+i]);
    dKp    = dODouble(oAssocObject(aaArgs[5]));
    dP0    = DEGTORAD * dODouble(oAssocObject(aaArgs[6]));
    dN     = dODouble(oAssocObject(aaArgs[7]));

    for ( i=0; i<ATOMSINTORSION; i++ ) {
        if ( !bContainerContainedBy( (CONTAINER)aaAtoms[i], (CONTAINER)uUnit ) ) {
            VP0(( "%s: Atom#%d is not part of the UNIT\n", sCmd, i+1 ));
            return(NULL);
        }
    }

                /* Create the RESTRAINT and add it to the UNIT */

    rRest = rRestraintCreate();
    RestraintTorsionSet( rRest, aaAtoms[0], aaAtoms[1], aaAtoms[2], aaAtoms[3],
                         dKp, dP0, dN );
    UnitAddRestraint( uUnit, rRest );

    return(NULL);
}



/*
 *      oCmd_deleteRestraint
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove a RESTRAINT bond, angle, or torsion depending
 *      on the number of atoms passed by the caller.
 */
OBJEKT
oCmd_deleteRestraint( int iArgCount, ASSOC aaArgs[] )
{
int             i, iAtomCount;
UNIT            uUnit;
ATOM            aaAtoms[4];
RESTRAINT       rRest;
BOOL            bGotOne;
char            *sCmd = "deleteRestraint";

    if ( iArgCount == 3 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a" ) ) {
          VP0(( "usage:  deleteRestraint <unit> <a> <b> [<c> <d>]\n" ));
          return(NULL);
        }
        iAtomCount = 2;
    } else if ( iArgCount == 4 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a a" ) ) {
          VP0(( "usage:  deleteRestraint <unit> <a> <b> [<c> <d>]\n" ));
          return(NULL);
        }
        iAtomCount = 3;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a a a" ) ) {
          VP0(( "usage:  deleteRestraint <unit> <a> <b> [<c> <d>]\n" ));
          return(NULL);
        }
        iAtomCount = 4;
    }

        /* Get the UNIT and the ATOMs that form the RESTRAINT */
    uUnit  = (UNIT)oAssocObject(aaArgs[0]);
    for ( i=0; i<iArgCount-1; i++ )
        aaAtoms[i] = (ATOM)oAssocObject(aaArgs[1+i]);

                /* Now loop through all the RESTRAINTs and find */
                /* the one the user specified */
    bGotOne = FALSE;
    UnitLoopRestraints( uUnit );
    while ( (rRest = rUnitNextRestraint(uUnit)) ) {
        switch ( iRestraintType(rRest) ) {
            case RESTRAINTBOND:
                if ( iAtomCount == 2 ) {
                    if ( bRestraintBondMatchAtoms( rRest, aaAtoms[0],
                                                   aaAtoms[1] ) ) {
                        bGotOne = TRUE;
                        goto DONE;
                    }
                }
                break;
            case RESTRAINTANGLE:
                if ( iAtomCount == 3 ) {
                    if ( bRestraintAngleMatchAtoms( rRest, aaAtoms[0],
                                                    aaAtoms[1], aaAtoms[2] )){
                        bGotOne = TRUE;
                        goto DONE;
                    }
                }
                break;
            case RESTRAINTTORSION:
                if ( iAtomCount == 4 ) {
                    if ( bRestraintTorsionMatchAtoms( rRest, aaAtoms[0],
                                                aaAtoms[1], aaAtoms[2],
                                                aaAtoms[3] ) ) {
                        bGotOne = TRUE;
                        goto DONE;
                    }
                }
                break;
            default:
                DFATAL(( "%s: Illegal RESTRAINT!\n", sCmd ));
        }
    }

DONE:
    if ( !bGotOne ) {
        VP0(( "%s: No such restraint could be found.\n", sCmd ));
    } else {
        VP1(( "Removing restraint.\n" ));
        bUnitRemoveRestraint( uUnit, rRest );
    }
    return(NULL);
}





/*
 *      oCmd_impose
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Impose the internal coordinates the user specified onto
 *      the RESIDUEs of the UNIT.
 *
 *      Arguments:
 *              [0]     - The UNIT to add change INTERNALs for.
 *              [1]     - The list of RESIDUEs to change.
 *              [2]     - The second ATOM of the RESTRAINT bond.
 *              [3]     - The third ATOM of the RESTRAINT bond.
 *              [4]     - The fourth ATOM of the RESTRAINT bond.
 *              [5]     - ODOUBLE: The value of Kp.
 *              [6]     - ODOUBLE: The value of P0.
 *              [7]     - ODOUBLE: The value of N.
 */
OBJEKT
oCmd_impose( int iArgCount, ASSOC aaArgs[] )
{
LIST            lResidues, lInternals;
OBJEKT          oOne;
LISTLOOP        llResidues, llOne;
RESIDUE         rRes;
#define         MAXOBJEKTSININTERNALLIST 101
OBJEKT          oaIntObjekts[MAXOBJEKTSININTERNALLIST];
int             iObjekts, iStart, iDum;
RESIDUE         rResModify;
LOOP            lAtoms, lSpanning;
UNIT            uUnit;
ASSOC           aInternal, aOne;
LISTLOOP        llInternals;
BOOL            bSkipSubList;
LIST            lOne;
char            *sCmd = "impose";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u l l" ) ) {
        VP0(( "usage:  %s <unit> <residueseqlist> <internalslistlist>\n",
                        sCmd ));
        return(NULL);
    }

    DisplayerAccumulateUpdates();
 
        /* Build the INTERNALs */

    uUnit = (UNIT)oAssocObject( aaArgs[0] );
    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    BuildInternalsUsingFlags( &lAtoms, 0, 0, 0, ATOMPOSITIONKNOWN );

        /* First get a list of RESIDUEs using oCmd_Residues */

    lResidues = lToolListOfResidues( uUnit, (LIST)oAssocObject(aaArgs[1]) );

        /* Now apply the INTERNALs to each of the RESIDUEs in turn */

    llResidues = llListLoop(lResidues);
    while ( (rRes = (RESIDUE)oListNext(&llResidues)) ) {

                /* Loop over all the sub-lists */

        lInternals = (LIST)oAssocObject(aaArgs[2]);
        llInternals = llListLoop(lInternals);
        while ( (aInternal = (ASSOC)oListNext(&llInternals)) ) {

            if ( iObjectType(oAssocObject(aInternal)) != LISTid ) {
                VP0(( "%s: Invalid internal list in internalslistlist !\n"
                      "        Note that argument #3 is a list of lists.\n"
                      "        Here is an example:\n"
                      "        impose x {167 168} { { \"CH3\" $C $N $CA 180.0"
                      " } }\n",
                      sCmd ));
                continue;
            }

                        /* Copy what is in the sub-list into an array */

            bSkipSubList = FALSE;
            iObjekts = 0;
            lOne = (LIST)oAssocObject(aInternal);
            llOne = llListLoop(lOne);
            while( (aOne = (ASSOC)oListNext(&llOne)) ) {
                oOne = (OBJEKT)oAssocObject(aOne);
                oaIntObjekts[iObjekts++] = oOne;
                if ( iObjekts > MAXOBJEKTSININTERNALLIST ) {
                    VP0(( "%s: Too many lists in argument #3\n", sCmd ));
                    bSkipSubList = TRUE;
                    break;
                }
            }

            if ( bSkipSubList ) continue;

                        /* Now apply the INTERNAL onto the RESIDUE */
            iStart = 0;
            if ( iObjectType(oaIntObjekts[0]) == ODOUBLEid ) iStart = 1;
        
            switch ( iObjekts - iStart ) {
                        /* Do a bond INTERNAL */
                case 3:
                    MESSAGE(( "Imposing a bond INTERNAL\n" ));
                    if ( iObjectType(oaIntObjekts[iStart+0]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+1]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+2]) != ODOUBLEid ) {
                        VP0(( "%s: Illegal bond internal definition\n", sCmd ));
                    } else {
                        if ( iStart == 1 ) {
                           rResModify = rResidueConnected( rRes, 
                                        (int)dODouble(oaIntObjekts[0]) );
                        } else rResModify = rRes;
                        MESSAGE(( "Imposing on Residue: %s:%d\n",
                                sContainerName((CONTAINER) rResModify),
                                iContainerSequence((CONTAINER) rResModify) ));
                        bBuildChangeInternalBond( (CONTAINER) rResModify,
                                sOString(oaIntObjekts[iStart+0]),
                                sOString(oaIntObjekts[iStart+1]),
                                dODouble(oaIntObjekts[iStart+2]) );
                    }
                    break;

                        /* Do an angle INTERNAL */
                case 4:
                    MESSAGE(( "Imposing an angle INTERNAL\n" ));
                    if ( iObjectType(oaIntObjekts[iStart+0]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+1]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+2]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+3]) != ODOUBLEid ) {
                        VP0(( "%s: Illegal angle internal definition\n", sCmd));
                    } else {
                        if ( iStart == 1 ) {
                           rResModify = rResidueConnected( rRes, 
                                        (int)dODouble(oaIntObjekts[0]) );
                        } else rResModify = rRes;
                        MESSAGE(( "Imposing on Residue: %s:%d\n",
                                sContainerName((CONTAINER) rResModify),
                                iContainerSequence((CONTAINER) rResModify) ));
                        bBuildChangeInternalAngle( (CONTAINER) rResModify,
                                sOString(oaIntObjekts[iStart+0]),
                                sOString(oaIntObjekts[iStart+1]),
                                sOString(oaIntObjekts[iStart+2]),
                                dODouble(oaIntObjekts[iStart+3])*DEGTORAD );
                    }
                    break;

                        /* Do a torsion INTERNAL */
                case 5:
                    MESSAGE(( "Imposing a torsion INTERNAL\n" ));
                    if ( iObjectType(oaIntObjekts[iStart+0]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+1]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+2]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+3]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+4]) != ODOUBLEid ) {
                        VP0(( "%s: Illegal angle internal definition\n", sCmd));
                    } else {
                        if ( iStart == 1 ) {
                           rResModify = rResidueConnected( rRes, 
                                        (int)dODouble(oaIntObjekts[0]) );
                        } else rResModify = rRes;
                        MESSAGE(( "Imposing on Residue: %s:%d\n",
                                sContainerName((CONTAINER) rResModify),
                                iContainerSequence((CONTAINER) rResModify) ));
                        bBuildChangeInternalTorsion( (CONTAINER) rResModify,
                                sOString(oaIntObjekts[iStart+0]),
                                sOString(oaIntObjekts[iStart+1]),
                                sOString(oaIntObjekts[iStart+2]),
                                sOString(oaIntObjekts[iStart+3]),
                                dODouble(oaIntObjekts[iStart+4])*DEGTORAD );
                    }
                    break;
                default:
                    VP0(( "%s: Improper internal definition\n", sCmd ));
                    break;
            }
        }
    }

                /* Now build the externals again */

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    lSpanning = lLoop( oNext(&lAtoms), SPANNINGTREE );
    iDum = 0;   /* for purify */
    BuildExternalsUsingFlags( &lSpanning, 0, 0,
                                ATOMPOSITIONKNOWN,
                                0,
                                &iDum, &iDum, &iDum, FALSE );

                /* Destroy all of the INTERNALs */

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    BuildDestroyInternals( &lAtoms );

    DisplayerReleaseUpdates();
 
    return(NULL);
}




/*
 *      oCmd_translate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Translate all of the ATOMs within the CONTAINER by
 *      the vector in [1].
 *
 *      Arguments:
 *              [0]     - The CONTAINER whose atoms to translate.
 *              [1]     - A list with the three coordinates of the vector.
 */
OBJEKT
oCmd_translate( int iArgCount, ASSOC aaArgs[] )
{
LIST            lVector;
LISTLOOP        llElements;
CONTAINER       cCont;
double          daVector[3];
int             i;
VECTOR          vOffset;
ASSOC           aAssoc;
char            *sCmd = "translate";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umral l" ) ) {
        VP0(( "usage:  %s <unit/residue/atom> <directionlist>\n", sCmd ));
        return(NULL);
    }

    DisplayerAccumulateUpdates();
    
                /* Get the CONTAINER to translate */

    cCont = (CONTAINER)oAssocObject( aaArgs[0] );

    lVector = (LIST)oAssocObject( aaArgs[1] );
    i = 0;
    llElements = llListLoop(lVector);
    while ( (aAssoc = (ASSOC)oListNext(&llElements)) ) {
        if ( iObjectType(oAssocObject(aAssoc)) == ODOUBLEid ) {
            if ( i<3 ) {
                daVector[i] = dODouble(oAssocObject(aAssoc));
                i++;
            } else {
                VP0(( "%s: Illegal vector\n", sCmd ));
                break;
            }
        } else {
            VP0(( "%s: Illegal vector\n", sCmd ));
            break;
        }
    }

    if ( i==3 ) {
        VectorDef( &vOffset, daVector[0], daVector[1], daVector[2] );
        ContainerTranslateBy( cCont, vOffset );
    }

    DisplayerReleaseUpdates();

    return(NULL);
}






/*
 *      oCmd_center
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Display the geometric center of the ATOMs in the CONTAINER.
 *
 *      Arguments:
 *              [0]     - The CONTAINER to use to find the geometric
 *                        center of ATOMs.
 */
OBJEKT
oCmd_center( int iArgCount, ASSOC aaArgs[] )
{
CONTAINER       cCont;
VECTOR          vOffset;

    if ( !bCmdGoodArguments( "center", iArgCount, aaArgs, "umral" ) ) {
          VP0(( "usage:  center <unit/residue/atom>\n" ));
          return(NULL);
    }

                /* Get the CONTAINER to translate */

    cCont = (CONTAINER)oAssocObject( aaArgs[0] );

    vOffset = vContainerGeometricCenter(cCont);

    VP0(( "The center is at: %4.2lf, %4.2lf, %4.2lf\n",
                dVX(&vOffset), dVY(&vOffset), dVZ(&vOffset) ));

    return(NULL);
}



/*
 *      oCmd_solvateCap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   UNIT to add cap of solvent to.
 *              [1] -   UNIT of solvent.
 *              [2] -   RESIDUE, MOLECULE, ATOM, LIST of atoms, or
 *                              LIST of 3 doubles to use to calculate
 *                              the point around which to place the cap.
 *              [3] -   ODOUBLE, radius of cap.
 *      Option  [4] -   ODOUBLE with closeness parameter.
 *
 */
OBJEKT
oCmd_solvateCap( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolute;
UNIT            uSolvent;
OBJEKT          oPosition;
VECTOR          vPos;
double          dRadius, dCloseness = 1.0;
int             iInitialSize, iFinalSize;
char            *sCmd = "solvateCap";

    if ( iArgCount == 4 ) {
      if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u mral n" ) ) {
        VP0(( 
"usage:  solvateCap <solute> <solvent> <position> <radius> <closeness>\n" ));
        return(NULL);
      }
      dCloseness = 1.0;
    } else {
      if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u mral n n" ) ) {
        VP0(( 
"usage:  solvateCap <solute> <solvent> <position> <radius> <closeness>\n" ));
        return(NULL);
      }
      dCloseness = dODouble(oAssocObject(aaArgs[4]));
    }

    uSolute = (UNIT)oAssocObject(aaArgs[0]);
    uSolvent = (UNIT)oAssocObject(aaArgs[1]);
    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );
    oPosition = oAssocObject(aaArgs[2]);

                /* Check if the oPosition argument is an array of three */
                /* doubles, if it is then treat it like a vector */

    if ( !bToolGeometricCenter( oPosition, &vPos ) ) {
        return(NULL);
    }

    
    dRadius = dODouble(oAssocObject(aaArgs[3]));
    if ( dRadius < 0.0 ) {
        VP0(( "radius (%f) must be > 0\n", dRadius ));
        return(NULL);
    }

    TurnOffDisplayerUpdates();
    iInitialSize = iContainerNumberOfChildren( (CONTAINER) uSolute );
    zToolSolvateInSphere( uSolute, uSolvent, &vPos, dRadius, dCloseness );

    iFinalSize = iContainerNumberOfChildren( (CONTAINER) uSolute );

    VP0(( "Added %d residues.\n", iFinalSize - iInitialSize ));

                /* Define the solvent cap within the UNIT */

    UnitSetSolventCap( uSolute, dVX(&vPos), dVY(&vPos), dVZ(&vPos), dRadius );
    UnitSetUseSolventCap( uSolute, TRUE );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( (CONTAINER) uSolute );
    
    return(NULL);
}





/*
 *      oCmd_addPath
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   OSTRING, directory to add to file search 
 *                      path list.
 *
 */
OBJEKT
oCmd_addPath( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "addPath", iArgCount, aaArgs, "s" ) ) {
        VP0(( "usage:  addPath <path>\n" ));
        return(NULL);
    }

    if ( BasicsAddDirectory( sOString(oAssocObject(aaArgs[0])), 0 ) )
        VP0(( "%s added to file search path.\n", 
                        sOString(oAssocObject(aaArgs[0])) ));

    return(NULL);
}




/*
 *      oCmd_crossLink
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a crosslink between two RESIDUEs
 *      according to the named connect ATOMs.
 *
 *      Arguments:
 *              [0] -   RESIDUE
 *              [1] -   OSTRING, name of connect ATOM on [0]
 *              [2] -   RESIDUE
 *              [3] -   OSTRING, name of connect ATOM on [2]
 *      option  [4] -   OSTRING, bond order, see BOND command.
 */
OBJEKT
oCmd_crossLink( int iArgCount, ASSOC aaArgs[] )
{
int             iConnectA, iConnectB;
RESIDUE         rA, rB;
OBJEKT          oObj;
int             iOrder;
char            *sCmd = "crossLink";
char            *usage =
        "usage:  crossLink <res1> <connect> <res2> <connect> [bondorder]\n";

    if ( iArgCount == 5 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "r s r s s" ) ) {
          VP0(( usage ));
          VP0(("(For multiple molecules, note that residue numbering jumps\n"));
          VP0((" by 1001 for each new molecule)\n"));
          return(NULL);
        }
        oObj = oAssocObject(aaArgs[4]);
        iOrder = iAtomBondOrderFromName(sOString(oObj));
        if ( iOrder == BONDNONE ) {
            VP0(( "%s: Illegal bond order\n", sCmd ));
            return(NULL);
        }
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "r s r s" ) ) {
          VP0(( 
"usage:  crossLink <res1> <connect> <res2> <connect> [bondorder]\n" ));
          return(NULL);
        }
        iOrder = BONDSINGLE;
    }

    oObj = oAssocObject(aaArgs[1]);
    iConnectA = iResidueConnectFromName(sOString(oObj));
    if ( iConnectA == NOEND ) {
        VP0(( "%s: Invalid connect atom: %s\n", sCmd, sOString(oObj) ));
        return(NULL);
    }
    oObj = oAssocObject(aaArgs[3]);
    iConnectB = iResidueConnectFromName(sOString(oObj));
    if ( iConnectB == NOEND ) {
        VP0(( "%s: Invalid connect atom: %s\n", sCmd, sOString(oObj) ));
        return(NULL);
    }

    rA = (RESIDUE)oAssocObject(aaArgs[0]);
    rB = (RESIDUE)oAssocObject(aaArgs[2]);

    DisplayerAccumulateUpdates();
    
    if ( !bResidueCrossLink( rA, iConnectA, rB, iConnectB, iOrder ) ) {
        VP0(( "%s: Could not form cross link, invalid connection atom\n", 
                                sCmd ));
    }

    DisplayerReleaseUpdates();
    
    return(NULL);
}


/*
 *      oCmd_addPdbResMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Append.
 *
 *      Arguments:
 *              [0] -   LIST of entries to add.
 *
 *      The Name Map is used to map residue names read from
 *      PDB files to variable names within LEaP.
 *      The LIST is a LIST of LISTs each sub-list containing
 *      two or three entries.  Each entry has the form:
 *
 *      { ODOUBLE OSTRING OSTRING }
 *      { OSTRING OSTRING }
 *
 *      ODOUBLE can be 0 or 1 
 *      The first OSTRING is the name within the PDB file.
 *      The second OSTRING is the variable name to map to.
 */
OBJEKT
oCmd_addPdbResMap( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "addPdbResMap", iArgCount, aaArgs, "l" ) ) {
         VP0(( "usage:  addPdbResMap <list_of_lists>\n" ));
         return(NULL);
    }
    PdbAppendToResMap( (LIST)oAssocObject(aaArgs[0]) );
    return(NULL);
}

/*
 *      oCmd_addPdbAtomMap
 *
 *      Author: Bill Ross
 */
OBJEKT
oCmd_addPdbAtomMap( int iArgCount, ASSOC aaArgs[] )
{

    if ( !bCmdGoodArguments( "addPdbMap", iArgCount, aaArgs, "l" ) ) {
         VP0(( "usage:  addPdbAtomMap <list_of_lists>\n" ));
         return(NULL);
    }

    PdbAppendToAtomMap( (LIST)oAssocObject(aaArgs[0]) );

    return(NULL);
}

/*
 *      oCmd_addAtomTypes
 *
 *      Author: Bill Ross (1996)
 */
OBJEKT
oCmd_addAtomTypes( int iArgCount, ASSOC aaArgs[] )
{
LIST            lList;
        if ( !bCmdGoodArguments( "addAtomTypes", iArgCount, aaArgs, "l" ) ) {
                VP0(( "usage:  addAtomTypes <list_of_lists>\n" ));
                return(NULL);
        }

        lList =  (LIST)oAssocObject(aaArgs[0]);

        AmberAddAtomTypes( lList );
        return(NULL);
}



/*
 *      oCmd_displayPdbResMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Display the Name Map.
 */
OBJEKT
oCmd_displayPdbResMap( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "displayPdbResMap", iArgCount, aaArgs, "" ) ) {
         VP0(( "usage:  displayPdbResMap\n" ));
         return(NULL);
    }
    PdbDisplayResMap();
    return(NULL);
}

OBJEKT
oCmd_displayPdbAtomMap( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "displayPdbAtomMap", iArgCount, aaArgs, "" ) ) {
         VP0(( "usage:  displayPdbAtomMap\n" ));
         return(NULL);
    }
    PdbDisplayAtomMap();
    return(NULL);
}





/*
 *      oCmd_clearPdbResMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Clear the Name Map.
 *
 *
 */
OBJEKT
oCmd_clearPdbResMap( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "clearPdbResMap", iArgCount, aaArgs, "" ) ) {
         VP0(( "usage:  clearPdbResMap\n" ));
         return(NULL);
    }
    PdbClearResMap();
    return(NULL);
}

OBJEKT
oCmd_clearPdbAtomMap( int iArgCount, ASSOC aaArgs[] )
{

    if ( !bCmdGoodArguments( "clearPdbAtomMap", iArgCount, aaArgs, "" ) ) {
         VP0(( "usage:  clearPdbAtomMap\n" ));
         return(NULL);
    }
    PdbClearAtomMap();
    return(NULL);
}





/*
 *      oCmd_measureGeom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *              Measure the distance, angle, torsion between
 *              two, three, or four ATOMs.
 *
 *      Arguments:
 *              [0] -   ATOM
 *              [1] -   ATOM
 *      option  [2] -   ATOM
 *      option  [3] -   ATOM
 *
 *      Return:
 *              ODOUBLE, return the internal coordinates value.
 */
OBJEKT
oCmd_measureGeom( int iArgCount, ASSOC aaArgs[] )
{
ATOM            aA, aB, aC, aD;
ODOUBLE         odVal;
double          dVal;
char            *sCmd = "measureGeom";

    if ( iArgCount == 4 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a a a" ) ) {
          VP0(( "usage:  measureGeom <atom> <atom> [atom atom]\n" ));
          return(NULL);
        }
        aA = (ATOM)oAssocObject(aaArgs[0]);
        aB = (ATOM)oAssocObject(aaArgs[1]);
        aC = (ATOM)oAssocObject(aaArgs[2]);
        aD = (ATOM)oAssocObject(aaArgs[3]);
        dVal = dVectorAtomTorsion(
                &vAtomPosition(aA), &vAtomPosition(aB),
                &vAtomPosition(aC), &vAtomPosition(aD) )/DEGTORAD;
        VP0(( "Torsion angle: %4.2lf degrees\n", dVal ));
    } else if ( iArgCount == 3 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a a" ) ) {
          VP0(( "usage:  measureGeom <atom> <atom> [atom atom]\n" ));
          return(NULL);
        }
        aA = (ATOM)oAssocObject(aaArgs[0]);
        aB = (ATOM)oAssocObject(aaArgs[1]);
        aC = (ATOM)oAssocObject(aaArgs[2]);
        dVal = dVectorAtomAngle(
                &vAtomPosition(aA), &vAtomPosition(aB),
                &vAtomPosition(aC) )/DEGTORAD;
        VP0(( "Angle: %4.2lf degrees\n", dVal ));
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a" ) ) {
          VP0(( "usage:  measureGeom <atom> <atom> [atom atom]\n" ));
          return(NULL);
        }
        aA = (ATOM)oAssocObject(aaArgs[0]);
        aB = (ATOM)oAssocObject(aaArgs[1]);
        dVal = dVectorAtomLength(
                &vAtomPosition(aA), &vAtomPosition(aB) );
        VP0(( "Distance: %4.2lf angstroms\n", dVal ));
    }
    odVal = (ODOUBLE)oCreate(ODOUBLEid);
    ODoubleSet( odVal, dVal );
    return((OBJEKT)odVal);
}


OBJEKT
oCmd_memDebug( int iArgCount, ASSOC aaArgs[] )
{
        iMemDebug = 1;
        return(NULL);
}

/*
 *      oCmd_bondByDistance
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Perform an N^2 search of all ATOMs within the
 *      container and create single bonds between all ATOMs
 *      that are within a certain distance of each other.
 *
 *      Arguments:
 *              [0] -   UNIT/MOLECULE/RESIDUE/LIST
 *      option  [1] -   ODOUBLE  maximum bonding distance
 */
OBJEKT
oCmd_bondByDistance( int iArgCount, ASSOC aaArgs[] )
{
double          dDist;
int             iBonds;
char            *sCmd = "bondByDistance";
char            *sUsage = "usage:  bondByDistance <unit> [maxdistance]\n";

    if ( iArgCount == 1 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umrl" ) ) {
            VP0(( sUsage ));
            return(NULL);
        }
        if (!(dDist = GDefaults.dDSearchDistance)) 
            dDist = DEFAULT_DISTANCE_SEARCH;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umrl n" ) ) {
            VP0(( sUsage ));
            return(NULL);
        }
        dDist = dODouble(oAssocObject(aaArgs[1]));
    }

    TurnOffDisplayerUpdates();

    iBonds = iToolDistanceSearch( (CONTAINER)oAssocObject(aaArgs[0]), dDist,
                                        TRUE, DISTANCE_SEARCH_CREATE_BONDS );

    VP0(( "Created %d bonds.\n", iBonds ));
    
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( (CONTAINER) oAssocObject(aaArgs[0]) );
    
    return(NULL);
}



/*
 *      oCmd_groupSelectedAtoms
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a group within the UNIT using all
 *      of the ATOMs that have ATOMSELECTED flag set.
 *
 *              [0] -   UNIT
 *              [1] -   OSTRING  name of group
 */
OBJEKT
oCmd_groupSelectedAtoms( int iArgCount, ASSOC aaArgs[] )
{
STRING          sGroup;
UNIT            uUnit;
LOOP            lAtoms;
ATOM            aAtom;
int             iAtoms;
char            *sCmd = "groupSelectedAtoms";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s" ) ) {
          VP0(( "usage:  groupSelectedAtoms <unit> <groupname>\n" ));
          return(NULL);
    }

    DisplayerAccumulateUpdates();

    uUnit = (UNIT)oAssocObject(aaArgs[0]);
    strcpy( sGroup, sOString(oAssocObject(aaArgs[1])));

    if ( lUnitGroup( uUnit, sGroup ) ) {
        bUnitGroupDestroy( uUnit, sGroup );
    }

    bUnitGroupCreate( uUnit, sGroup );

    iAtoms = 0;
    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
            bUnitGroupAddAtom( uUnit, sGroup, aAtom );
            iAtoms++;
        }
    }
    VP0(( "Added %d atoms.\n", iAtoms ));
    
    DisplayerReleaseUpdates();
    
    return(NULL);
}







/*
 *      oCmd_transform
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Translate all of the ATOMs within the CONTAINER by
 *      the vector in [1].
 *
 *      Arguments:
 *              [0]     - The CONTAINER whose atoms to transform.
 *              [1]     - A list of lists representing a 4x4 matrix.
 */
OBJEKT
oCmd_transform( int iArgCount, ASSOC aaArgs[] )
{
LIST            lVectorX;
LIST            lVectorY;
LISTLOOP        llElementsX;
LISTLOOP        llElementsY;
CONTAINER       cCont;
MATRIX          mTransform;
int             iX, iY;
ASSOC           aAssocX;
ASSOC           aAssocY;
char            *sCmd = "transform";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umral l" ) ) {
          VP0(( "usage:  transform <atoms> <matrixlist>\n" ));
          return(NULL);
    }

                /* Get the CONTAINER to translate */

    cCont = (CONTAINER)oAssocObject( aaArgs[0] );

    lVectorY = (LIST)oAssocObject( aaArgs[1] );
    iX = 0;
    iY = 0;
    MatrixIdentity( mTransform );
    llElementsY = llListLoop(lVectorY);
    while ( (aAssocY = (ASSOC)oListNext(&llElementsY)) ) {
        if ( iObjectType(oAssocObject(aAssocY)) == LISTid ) {
            if ( iY<4 ) {
                lVectorX = (LIST)oAssocObject(aAssocY);
                llElementsX = llListLoop(lVectorX);
                while ( (aAssocX = (ASSOC)oListNext(&llElementsX)) ) {
                    if ( iObjectType(oAssocObject(aAssocX)) == ODOUBLEid ) {
                        if ( iX<4 ) {
                            mTransform[iX][iY] = dODouble(oAssocObject(aAssocX));
                            iX++;
                        } else { goto ERROR; }
                    } else { goto ERROR; }
                }
                iX = 0;
                iY++;
            } else { goto ERROR; }
        } else { goto ERROR; }
    }

    DisplayerAccumulateUpdates();

    ContainerTransformBy( cCont, mTransform );

    DisplayerReleaseUpdates();
   
    return(NULL);
 
ERROR:

    VP0(( "%s: Illegal matrix\n", sCmd ));
    return(NULL);
}



/*
 *      oCmd_copy
 *
 *      Author: Christian Schafmeister (1991)
 *
 */
OBJEKT
oCmd_copy( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oNew;

    if ( !bCmdGoodArguments( "copy", iArgCount, aaArgs, "umranp" ) ) {
          VP0(( "usage:  <newvariable> = copy <variable>\n" ));
          return(NULL);
    }

    oNew = oCopy(oAssocObject(aaArgs[0]));

    return(oNew);
}



/*
 *      oCmd_listOff
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      List the contents of a library.
 */
OBJEKT
oCmd_listOff( int iArgCount, ASSOC aaArgs[] )
{
STRING          sFilename;
LIBRARY         lLib;
char            *cPNext;

    if ( !bCmdGoodArguments( "listOff", iArgCount, aaArgs, "s" ) ) {
          VP0(( "usage:  listOff <filename>\n" ));
          return(NULL);
    }
    strcpy( sFilename, sOString(oAssocObject(aaArgs[0])) );
   
    lLib = lLibraryOpen( sFilename, OPENREADONLY ); 
    if ( lLib == NULL ) return(NULL);

    VP0(( "Index of library: %s\n", sFilename )); 

    LibraryLoop( lLib );
    while ( (cPNext = sLibraryNext(lLib)) ) {
        VP0(( "%s\n", cPNext ));
    }

    LibraryClose( &lLib );
    
    return(NULL);
}








/*
 *      oCmd_deleteOffLibEntry
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove an entry from a library.
 */
OBJEKT
oCmd_deleteOffLibEntry( int iArgCount, ASSOC aaArgs[] )
{
STRING          sFilename;
LIBRARY         lLib;
BOOL            bRemoved;
STRING          sEntry;
char            *sCmd = "deleteOffLibEntry";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s s" ) ) {
          VP0(( "usage:  deleteOffLibEntry <filename> <entry>\n" ));
          return(NULL);
    }
    strcpy( sFilename, sOString(oAssocObject(aaArgs[0])) );
    strcpy( sEntry, sOString(oAssocObject(aaArgs[1])) );
   
    lLib = lLibraryOpen( sFilename, OPENREADWRITE ); 
    if ( lLib == NULL ) return(NULL);

    bRemoved = bLibraryRemove( lLib, sEntry );

    if ( !bRemoved ) {
        VP0(( "%s: %s was not found.\n", sCmd, sEntry ));
    } else {
        VP0(( "%s was removed.\n", sEntry ));
    }

    LibraryClose(&lLib);
    
    return(NULL);
}




/*
 *      oCmd_mutate
 *
 *      Author: Christian Schafmeister (1991)
 *
 */
OBJEKT
oCmd_mutate( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;
ODOUBLE         odSeqNum;
RESIDUE         rNew;
RESIDUE         rOld;
RESIDUE         rCopy;
int             iSeqNum;
LOOP            lAtoms;
ATOM            aAtom;
ATOM            aNeighbor;
int             i, iNext;
char            *sCmd = "mutate";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u n r" ) ) {
          VP0(( "usage:  mutate <unit> <number> <residue>\n" ));
          return(NULL);
    }
    uUnit = (UNIT)oAssocObject( aaArgs[0] );
    odSeqNum = (ODOUBLE)oAssocObject( aaArgs[1] );
    rNew = (RESIDUE)oAssocObject( aaArgs[2] );

    DisplayerAccumulateUpdates();

    iSeqNum = (int)dODouble(odSeqNum);

    rOld = (RESIDUE)cContainerFindSequence( (CONTAINER) uUnit, 
                                                RESIDUEid, iSeqNum );
    if ( rOld == NULL ) {
        VP0(( "%s: Could not find residue with sequence number: %d\n", 
                                        sCmd, iSeqNum ));
        goto FAIL;
    }

                /* Make sure there are no bonds out of the new residue */

    lAtoms = lLoop( (OBJEKT)rNew, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
            aNeighbor = aAtomBondedNeighbor( aAtom, i );
            if ( rNew != (RESIDUE)cContainerWithin((CONTAINER) aNeighbor) ) {
                VP0(( "%s: The mutant residue cannot be bonded to anything.\n",
                                                sCmd ));
                goto FAIL;
            }
        }
    }

    rCopy = (RESIDUE)oCopy( (OBJEKT)rNew );

                /* Perform the mutation */

    ResidueMutate( rCopy, rOld );

                /* Remove the old RESIDUE from the UNIT */
                /* And add the new RESIDUE, without incrementing */
                /* the UNITs next child sequence number */

    REF( rOld );  /* bContainerRemove() needs this */
    bContainerRemove( (CONTAINER)uUnit, (OBJEKT)rOld );
    iNext = iContainerNextChildsSequence( (CONTAINER) uUnit );
    ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rCopy );
    ContainerSetSequence( (CONTAINER) rCopy, 
                        iContainerSequence((CONTAINER) rOld) );
    ContainerSetNextChildsSequence( (CONTAINER) uUnit, iNext );

                /* DEREF the old RESIDUE */

    DEREF( rOld );
    goto RET;

FAIL:
    VP0(( "Mutation failed.\n" ));

RET:
    DisplayerReleaseUpdates();
    return(NULL);
}


/*
 *      oCmd_addIons
 *
 *      Author: Bill Ross (1993)
 *
 *      addIons unit, ion1, #ion1, [ion2, #ion2]
 *
 *              [0] -   UNIT
 *              [1] -   UNIT    ion to add
 *              [2] -   NUMBER  number of ions to add (0 = neutralize)
 *      Optional[3] -   UNIT    second ion to add
 *      Optional[4] -   NUMBER  number of the second ion to add
 *
 *      Adds Counter Ions.
 *      A Coulomb's Law ESP grid is used to calculate
 *      the appropriate location for the ion.
 *
 *      How it works:
 *
 *      1.  If [2] = 0, [0] must be charged and [1] must be opposite in charge;
 *              [0] is neutralized with [1].
 *      2.  If [2] != 0, [1] (and optionally [3]) are added (in alternation).
 *
 */

/*
 *  vaSolventResidues() - make array of residue pointers
 *      TODO - maybe restrict it to dShellExtent?
 */
VARARRAY        vaSolventResidues( uUnit )
UNIT            uUnit;
{
VARARRAY        vaSolvent;
LOOP            lRes;
RESIDUE         rRes;

        vaSolvent = vaVarArrayCreate( sizeof(RESIDUE) );
        lRes = lLoop( (OBJEKT)uUnit, RESIDUES );
        while ( (rRes = (RESIDUE)oNext( &lRes )) )
                if ( cResidueType( rRes ) == RESTYPESOLVENT ) {
                        VarArrayAdd( vaSolvent, (GENP)&rRes );
                }
        if ( !iVarArrayElementCount( vaSolvent ) )
                VarArrayDestroy( &vaSolvent );
        return( vaSolvent );
}

static void
CheckSolvent( UNIT uUnit, VARARRAY vaSolvent, UNIT uIon, VECTOR *PvIon )
{
RESIDUE         *PrRes, *PrClosest;
LOOP            lAtoms;
ATOM            aAtom;
VECTOR          vClosest;
double          d2, dmin2, x, y, z;
int             i, iCount;

        /*
         *  PrRes (pointer to residue pointer) is used so that
         *      the residue pointer can be set null if the residue
         *      is deleted
         */
        dmin2 = FLT_MAX;
        iCount = iVarArrayElementCount( vaSolvent );
        PrRes = PVAI( vaSolvent, RESIDUE, 0 );
        for (i=0; i<iCount; i++, PrRes++) {
                if ( *PrRes == NULL ) /* already deleted */
                        continue;
                lAtoms = lLoop( (OBJEKT)*PrRes, ATOMS );
                aAtom = (ATOM) oNext( &lAtoms );
                x = PvIon->dX - vAtomPosition( aAtom ).dX;
                x = x * x;
                y = PvIon->dY - vAtomPosition( aAtom ).dY;
                y = y * y;
                z = PvIon->dZ - vAtomPosition( aAtom ).dZ;
                z = z * z;
                d2 = x + y + z;
                if ( d2 < dmin2 ) {
                        PrClosest = PrRes;
                        vClosest = vAtomPosition( aAtom );
                        dmin2 = d2;
                }
        }
        if ( dmin2 < 9 ) {  /* HACK test */
                VP0(("(Replacing solvent molecule)\n"));
                REF( *PrClosest );  /* bContainerRemove() needs this */
                if ( bContainerRemove( (CONTAINER)uUnit, (OBJEKT)*PrClosest )) {
                        ContainerDestroy( (CONTAINER *) PrClosest );
                        *PrClosest = NULL;
                        *PvIon = vClosest;
                } else
                        VP0(( "solvent removal failed\n" ));
                DEREF( *PrClosest );  /* after bContainerRemove()/Add */
        } else
                VP0(( "(No solvent overlap)\n"));
        return;
}

OBJEKT
oCmd_addIons( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit=NULL, uIon1=NULL, uIon2=NULL, uPlace=NULL;
int             iIon1=0, iIon2=0;
double          dCharge, dPertCharge, dICharge1, dICharge2;
double          dIonSize1, dIonSize2, dMinSize;
int             i, iUnknown, ierr;
VECTOR          vNewPoint, vMaxPot, vMinPot;
HELP            hTemp;
LOOP            lAtoms;
ATOM            aAtom;
OCTREE          octTreeSolute; 
VARARRAY        vaSolvent = NULL;
char            *sCmd = "addIons";

    BasicsResetInterrupt();
    
    /*
     *  Test args
     */
    ierr = 0;
    switch( iArgCount ) {
        case 3: 
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" ))
                ierr++;
          break;
        case 5: 
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n" )) {
                ierr++;
                break;
          }
          /*
           *  Translate the 2 extra args
           */
          uIon2 = (UNIT)oAssocObject( aaArgs[3] );
          iIon2 = (int)dODouble( oAssocObject( aaArgs[4] )); 
          if ( uIon2  &&  iIon2 == 0 ) {
              VP0(( "%s: '0' is not allowed as the value for the second ion.\n",
                                                sCmd ));
              ierr++;
          }
          break;
        default:
          ierr++;
          break;
    } /* end of switch */
    
    if ( ierr ) {
          hTemp = hHelp( "addions" );
          if ( hTemp == NULL ) {
                VP0(( "No help available on addIons\n" ));
          } else {
                VP0(( "\n%s\n", sHelpText(hTemp) ));
          }
          return(NULL);
    }

    /*
     *  Translate args common to both cases
     */
    uUnit = (UNIT)oAssocObject( aaArgs[0] );
    uIon1 = (UNIT)oAssocObject( aaArgs[1] );
    iIon1 = (int)dODouble( oAssocObject( aaArgs[2] ));

    /*
     *  Consider target unit's charge
     */
    ContainerTotalCharge( (CONTAINER) uUnit, &dCharge, &dPertCharge );
    if ( !dCharge ) {
        VP0(( "%s has a charge of 0.\n", sAssocName( aaArgs[1] )));
        if ( iIon1 == 0 ) {
                VP0(( "%s: Can't neutralize.\n", sCmd ));
                return(NULL);
        }
        VP0(( "Adding the ions anyway.\n"));
    } else
        MESSAGE(( "dCharge:  %4.2lf\n", dCharge ));
    
    /*
     *  Consider ion(s) charge
     */
    ContainerTotalCharge((CONTAINER)uIon1, &dICharge1, &dPertCharge );
    if ( !dICharge1 ) {
        VP0(( "%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sAssocName( aaArgs[1] )));
        return(NULL);
    }
    if ( uIon2 ) {
        ContainerTotalCharge((CONTAINER)uIon2, &dICharge2, &dPertCharge );
        if ( !dICharge2 ) {
           VP0(( "%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sAssocName( aaArgs[3] )));
           return(NULL);
        }
    }

    /*
     *  Consider neutralization
     */
    if ( iIon1 == 0 ) {
        if ( (dICharge1 < 0  &&  dCharge < 0) ||
             (dICharge1 > 0  &&  dCharge > 0)) {
                VP0(( "%s: 1st Ion & target are the same charge:\n", sCmd ));
                VP0(( "     can't neutralize.\n" ));
                return(NULL);
        }
        /*
         *  Get the nearest integer number of ions that
         *      we need to add to get as close as possible
         *      to neutral
         */
        iIon1 = (int)lrint( fabs(dCharge) / fabs(dICharge1) );
        if ( iIon1 == 0 )
            VP0(( " %f %d %d %d\n", fabs( dCharge),
                (int)fabs( dCharge), (int)fabs( dICharge1 ),
                (int)(fabs( dCharge) / fabs( dICharge1 )) ));
        if ( uIon2 ) {
                VP0(( "%s: Neutralization - can't do 2nd ion.\n", sCmd ));
                return(NULL);
        }
        VP0(( "%d %s ion%s required to neutralize.\n", iIon1,
                sAssocName( aaArgs[1] ), (iIon1 > 1 ? "s" : "") ));
    } 

    /*
     *  Consider ion sizes and positions.
     */
    dIonSize1 = 0.0;
    iUnknown = 0;
    lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
    for(i=0; (aAtom = (ATOM)oNext(&lAtoms)); i++) {
        if ( iAtomSetTmpRadius( aAtom ) )
                VP0(( "Using default radius %5.2f for ion %s\n",
                        ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[1] ) ));
        dIonSize1 = MAX( dIonSize1, dAtomTemp( aAtom ) );
        if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                iUnknown++;
    }
    if ( i > 1 ) {
        if ( iUnknown ) {
            VP0(( "Ion %s is polyatomic and has %d atoms w/ no position\n",
                sAssocName( aaArgs[1] ), iUnknown ));
            return( NULL );
        }
        VP0(( "Ion %s is polyatomic; multiplying max radius %5.2f by # atoms\n",
                sAssocName( aaArgs[1] ), dIonSize1 ));
        dIonSize1 *= (double) i;
    } else if ( iUnknown ) {
        VECTOR          vPos;

        VectorDef( &vPos, 0.0, 0.0, 0.0 );
        lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
        aAtom = (ATOM)oNext(&lAtoms);
        AtomSetPosition( aAtom, vPos );
        AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
    }
    dMinSize = dIonSize1;
    dIonSize2 = 0.0;
    if ( uIon2 ) {
        iUnknown = 0;
        lAtoms = lLoop( (OBJEKT)uIon2, ATOMS );
        for(i=0; (aAtom = (ATOM)oNext(&lAtoms)); i++) {
                if ( iAtomSetTmpRadius( aAtom ) )
                        VP0(( "Using default radius %5.2f for ion %s\n",
                                ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[3] ) ));
                dIonSize2 = MAX( dIonSize2, dAtomTemp( aAtom ) );
                if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                    iUnknown++;
        }
        if ( i > 1 ) {
            if ( iUnknown ) {
                VP0(( "Ion %s is polyatomic and has %d atoms w/ no position\n",
                    sAssocName( aaArgs[1] ), iUnknown ));
                return( NULL );
            }
            VP0(( 
                "Ion %s is polyatomic; multiplying max radius %5.2f by # atoms",
                sAssocName( aaArgs[3] ), dIonSize2 ));
            dIonSize2 *= (double) i;
        } else if ( iUnknown ) {
            VECTOR              vPos;

            VectorDef( &vPos, 0.0, 0.0, 0.0 );
            lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
            aAtom = (ATOM)oNext(&lAtoms);
            AtomSetPosition( aAtom, vPos );
            AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
        }
        dMinSize = MIN( dIonSize1, dIonSize2 );
    }

    VP0(( "Adding %d counter ions to \"%s\" using 1A grid\n", 
                iIon1 + iIon2, sAssocName( aaArgs[0] )));

    if (iIon1 + iIon2 > 5) {
        const double xx = (double) (iIon1 + iIon2);
        const double ff = exp(log(xx + 1.0)/3.0);
        dMinSize = (dIonSize1 > dIonSize2 ? dIonSize1 : dIonSize2);
        dMinSize *= (ff > 1.0 ? ff : 1.0);
    }

    if ( iIon1 + iIon2 == 0 )
        return(NULL);

    /*
     *  Build grid and calc potential on it.
     */
    octTreeSolute = octOctTreeCreate( uUnit, OCT_SHELL, 
                                        GDefaults.dGridSpace, dMinSize, GDefaults.dShellExtent, 0 );
    if ( !octTreeSolute ) {
        VP0(( "%s: No solute to add ions to\n", sCmd ));
        return(NULL);
    }
    vaSolvent = vaSolventResidues( uUnit );
    
    if ( vaSolvent ) {
        VP0(( "Solvent present: replacing closest with ion\n" ));
        VP0(( "\t when steric overlaps occur\n" ));
    } else
        VP0((" (no solvent present)\n" ));

    TurnOffDisplayerUpdates();
    
    OctTreeInitCharges( octTreeSolute, AT_OCTREE, GDefaults.iDielectricFlag, 
                                        dIonSize1, &vMinPot, &vMaxPot );
/*
OctTreePrintGrid( octTreeSolute, "Charge", COLOR_RANGE );
*/

    while ( iIon1 || iIon2 ) {
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iIon1 ) {
                if ( dICharge1 < 0 )
                        vNewPoint = vMaxPot;
                else
                        vNewPoint = vMinPot;
                if ( vaSolvent )
                        CheckSolvent( uUnit, vaSolvent, uIon1, &vNewPoint );

                /*
                 *  Make a copy of ion unit and give it new point.
                 */
                uPlace = (UNIT) oCopy( (OBJEKT)uIon1 );
                ContainerCenterAt( (CONTAINER) uPlace, vNewPoint );

                /*
                 *  Add ion to solute.
                 */
                UnitJoin( uUnit, uPlace );
                VP0(( "Placed %s in %s at (%4.2lf, %4.2lf, %4.2lf).\n", 
                        sAssocName( aaArgs[1] ), sAssocName( aaArgs[0] ),
                        dVX(&(vNewPoint)), 
                        dVY(&(vNewPoint)), 
                        dVZ(&(vNewPoint))));
                /*
                 *  Delete ion from grid (allowing clearance to most likely
                 *      future adjacent ion) and update esp.
                 */
                OctTreeDeleteSphere( octTreeSolute, &vNewPoint, 
                                dIonSize1 + (iIon2 ? dIonSize2 : dIonSize1) );
                OctTreeUpdateCharge( octTreeSolute, &vNewPoint, 
                        (float)dICharge1, (iIon2 ? dIonSize2 : dIonSize1),
                        &vMaxPot, &vMinPot );
                iIon1--;
        }
        if ( iIon2 ) {
                if ( dICharge2 < 0 )
                        vNewPoint = vMaxPot;
                else
                        vNewPoint = vMinPot;
                if ( vaSolvent ) 
                        CheckSolvent( uUnit, vaSolvent, uIon2, &vNewPoint );

                /*
                 *  Make a copy of ion unit and give it new point.
                 */
                uPlace = (UNIT) oCopy( (OBJEKT)uIon2 );
                ContainerCenterAt( (CONTAINER) uPlace, vNewPoint );

                /*
                 *  Add ion to solute.
                 */
                UnitJoin( uUnit, uPlace );
                VP0(( "Placed %s in %s at (%4.2lf, %4.2lf, %4.2lf).\n", 
                        sAssocName( aaArgs[3] ), sAssocName( aaArgs[0] ),
                        dVX(&(vNewPoint)), 
                        dVY(&(vNewPoint)), 
                        dVZ(&(vNewPoint))));
                /*
                 *  Delete ion from grid (allowing clearance to most likely
                 *      future adjacent ion) and update esp.
                 */
                OctTreeDeleteSphere( octTreeSolute, &vNewPoint, 
                                dIonSize2 + (iIon1 ? dIonSize1 : dIonSize2) );
                OctTreeUpdateCharge( octTreeSolute, &vNewPoint, 
                        (float)dICharge2, (iIon1 ? dIonSize1 : dIonSize2),
                        &vMaxPot, &vMinPot );
                iIon2--;
        }
    }
/*
OctTreePrintGrid( octTreeSolute, "Charge2", COLOR_RANGE );
*/
    VP0(( "\nDone adding ions.\n" ));
    OctTreeDestroy( &octTreeSolute );
    if ( vaSolvent )
        VarArrayDestroy( &vaSolvent );
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( (CONTAINER) uUnit );
    return(NULL);

CANCEL:

    VP0(( "\n%s: Interrupted.\n", sCmd ));
    BasicsResetInterrupt();
    if ( octTreeSolute )
        OctTreeDestroy( &octTreeSolute );    
    if ( vaSolvent )
        VarArrayDestroy( &vaSolvent );
    DisplayerReleaseUpdates();
    return(NULL);
}

OBJEKT
oCmd_addIons2( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit=NULL, uIon1=NULL, uIon2=NULL, uPlace=NULL;
int             iIon1=0, iIon2=0;
double          dCharge, dPertCharge, dICharge1, dICharge2;
double          dIonSize1, dIonSize2, dMinSize;
int             i, iUnknown, ierr;
VECTOR          vNewPoint, vMaxPot, vMinPot;
HELP            hTemp;
LOOP            lAtoms;
ATOM            aAtom;
OCTREE          octTreeSolute; 
char            *sCmd = "addIons";

    BasicsResetInterrupt();
    
    /*
     *  Test args
     */
    ierr = 0;
    switch( iArgCount ) {
        case 3: 
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" ))
                ierr++;
          break;
        case 5: 
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n" )) {
                ierr++;
                break;
          }
          /*
           *  Translate the 2 extra args
           */
          uIon2 = (UNIT)oAssocObject( aaArgs[3] );
          iIon2 = (int)dODouble( oAssocObject( aaArgs[4] )); 
          if ( uIon2  &&  iIon2 == 0 ) {
              VP0(( "%s: '0' is not allowed as the value for the second ion.\n",
                                                sCmd ));
              ierr++;
          }
          break;
        default:
          ierr++;
          break;
    } /* end of switch */
    
    if ( ierr ) {
          hTemp = hHelp( "addions" );
          if ( hTemp == NULL ) {
                VP0(( "No help available on addIons\n" ));
          } else {
                VP0(( "\n%s\n", sHelpText(hTemp) ));
          }
          return(NULL);
    }

    /*
     *  Translate args common to both cases
     */
    uUnit = (UNIT)oAssocObject( aaArgs[0] );
    uIon1 = (UNIT)oAssocObject( aaArgs[1] );
    iIon1 = (int)dODouble( oAssocObject( aaArgs[2] ));

    /*
     *  Consider target unit's charge
     */
    ContainerTotalCharge( (CONTAINER) uUnit, &dCharge, &dPertCharge );
    if ( !dCharge ) {
        VP0(( "%s has a charge of 0.\n", sAssocName( aaArgs[1] )));
        if ( iIon1 == 0 ) {
                VP0(( "%s: Can't neutralize.\n", sCmd ));
                return(NULL);
        }
        VP0(( "Adding the ions anyway.\n"));
    } else
        MESSAGE(( "dCharge:  %4.2lf\n", dCharge ));
    
    /*
     *  Consider ion(s) charge
     */
    ContainerTotalCharge((CONTAINER)uIon1, &dICharge1, &dPertCharge );
    if ( !dICharge1 ) {
        VP0(( "%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sAssocName( aaArgs[1] )));
        return(NULL);
    }
    if ( uIon2 ) {
        ContainerTotalCharge((CONTAINER)uIon2, &dICharge2, &dPertCharge );
        if ( !dICharge2 ) {
           VP0(( "%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sAssocName( aaArgs[3] )));
           return(NULL);
        }
    }

    /*
     *  Consider neutralization
     */
    if ( iIon1 == 0 ) {
        if ( (dICharge1 < 0  &&  dCharge < 0) ||
             (dICharge1 > 0  &&  dCharge > 0)) {
                VP0(( "%s: 1st Ion & target are the same charge:\n", sCmd ));
                VP0(( "     can't neutralize.\n" ));
                return(NULL);
        }
        /*
         *  Get the nearest integer number of ions that
         *      we need to add to get as close as possible
         *      to neutral
         */
        iIon1 = (int)lrint( fabs(dCharge) / fabs(dICharge1) );
        if ( iIon1 == 0 )
            VP0(( " %f %d %d %d\n", fabs( dCharge),
                (int)fabs( dCharge), (int)fabs( dICharge1 ),
                (int)(fabs( dCharge) / fabs( dICharge1 )) ));
        if ( uIon2 ) {
                VP0(( "%s: Neutralization - can't do 2nd ion.\n", sCmd ));
                return(NULL);
        }
        VP0(( "%d %s ion%s required to neutralize.\n", iIon1,
                sAssocName( aaArgs[1] ), (iIon1 > 1 ? "s" : "") ));
    } 

    /*
     *  Consider ion sizes and positions.
     */
    dIonSize1 = 0.0;
    iUnknown = 0;
    lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
    for(i=0; (aAtom = (ATOM)oNext(&lAtoms)); i++) {
        if ( iAtomSetTmpRadius( aAtom ) )
                VP0(( "Using default radius %5.2f for ion %s\n",
                        ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[1] ) ));
        dIonSize1 = MAX( dIonSize1, dAtomTemp( aAtom ) );
        if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                iUnknown++;
    }
    if ( i > 1 ) {
        if ( iUnknown ) {
            VP0(( "Ion %s is polyatomic and has %d atoms w/ no position\n",
                sAssocName( aaArgs[1] ), iUnknown ));
            return( NULL );
        }
        VP0(( "Ion %s is polyatomic; multiplying max radius %5.2f by # atoms\n",
                sAssocName( aaArgs[1] ), dIonSize1 ));
        dIonSize1 *= (double) i;
    } else if ( iUnknown ) {
        VECTOR          vPos;

        VectorDef( &vPos, 0.0, 0.0, 0.0 );
        lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
        aAtom = (ATOM)oNext(&lAtoms);
        AtomSetPosition( aAtom, vPos );
        AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
    }
    dMinSize = dIonSize1;
    if ( uIon2 ) {
        dIonSize2 = 0.0;
        iUnknown = 0;
        lAtoms = lLoop( (OBJEKT)uIon2, ATOMS );
        for(i=0; (aAtom = (ATOM)oNext(&lAtoms)); i++) {
                if ( iAtomSetTmpRadius( aAtom ) )
                        VP0(( "Using default radius %5.2f for ion %s\n",
                                ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[3] ) ));
                dIonSize2 = MAX( dIonSize2, dAtomTemp( aAtom ) );
                if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                    iUnknown++;
        }
        if ( i > 1 ) {
            if ( iUnknown ) {
                VP0(( "Ion %s is polyatomic and has %d atoms w/ no position\n",
                    sAssocName( aaArgs[1] ), iUnknown ));
                return( NULL );
            }
            VP0(( 
                "Ion %s is polyatomic; multiplying max radius %5.2f by # atoms",
                sAssocName( aaArgs[3] ), dIonSize2 ));
            dIonSize2 *= (double) i;
        } else if ( iUnknown ) {
            VECTOR              vPos;

            VectorDef( &vPos, 0.0, 0.0, 0.0 );
            lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
            aAtom = (ATOM)oNext(&lAtoms);
            AtomSetPosition( aAtom, vPos );
            AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
        }
        dMinSize = MIN( dIonSize1, dIonSize2 );
    }

    VP0(( "Adding %d counter ions to \"%s\" using 1A grid\n", 
                iIon1 + iIon2, sAssocName( aaArgs[0] )));

    if ( iIon1 + iIon2 == 0 )
        return(NULL);

    /*
     *  Build grid and calc potential on it.
     */
    octTreeSolute = octOctTreeCreate( uUnit, OCT_SHELL, 
                                        GDefaults.dGridSpace, dMinSize, GDefaults.dShellExtent, 1 );
    if ( !octTreeSolute ) {
        VP0(( "%s: No atoms to add ions to\n", sCmd ));
        return(NULL);
    }
    

    TurnOffDisplayerUpdates();
    
    OctTreeInitCharges( octTreeSolute, AT_OCTREE, GDefaults.iDielectricFlag, 
                                        dIonSize1, &vMinPot, &vMaxPot );
/*
OctTreePrintGrid( octTreeSolute, "Charge", COLOR_RANGE );
*/

    while ( iIon1 || iIon2 ) {
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iIon1 ) {
                if ( dICharge1 < 0 )
                        vNewPoint = vMaxPot;
                else
                        vNewPoint = vMinPot;

                /*
                 *  Make a copy of ion unit and give it new point.
                 */
                uPlace = (UNIT) oCopy( (OBJEKT)uIon1 );
                ContainerCenterAt( (CONTAINER) uPlace, vNewPoint );

                /*
                 *  Add ion to solute.
                 */
                UnitJoin( uUnit, uPlace );
                VP0(( "Placed %s in %s at (%4.2lf, %4.2lf, %4.2lf).\n", 
                        sAssocName( aaArgs[1] ), sAssocName( aaArgs[0] ),
                        dVX(&(vNewPoint)), 
                        dVY(&(vNewPoint)), 
                        dVZ(&(vNewPoint))));
                /*
                 *  Delete ion from grid (allowing clearance to most likely
                 *      future adjacent ion) and update esp.
                 */
                OctTreeDeleteSphere( octTreeSolute, &vNewPoint, 
                                dIonSize1 + (iIon2 ? dIonSize2 : dIonSize1) );
                OctTreeUpdateCharge( octTreeSolute, &vNewPoint, 
                        (float)dICharge1, (iIon2 ? dIonSize2 : dIonSize1),
                        &vMaxPot, &vMinPot );
                iIon1--;
        }
        if ( iIon2 ) {
                if ( dICharge2 < 0 )
                        vNewPoint = vMaxPot;
                else
                        vNewPoint = vMinPot;

                /*
                 *  Make a copy of ion unit and give it new point.
                 */
                uPlace = (UNIT) oCopy( (OBJEKT)uIon2 );
                ContainerCenterAt( (CONTAINER) uPlace, vNewPoint );

                /*
                 *  Add ion to solute.
                 */
                UnitJoin( uUnit, uPlace );
                VP0(( "Placed %s in %s at (%4.2lf, %4.2lf, %4.2lf).\n", 
                        sAssocName( aaArgs[3] ), sAssocName( aaArgs[0] ),
                        dVX(&(vNewPoint)), 
                        dVY(&(vNewPoint)), 
                        dVZ(&(vNewPoint))));
                /*
                 *  Delete ion from grid (allowing clearance to most likely
                 *      future adjacent ion) and update esp.
                 */
                OctTreeDeleteSphere( octTreeSolute, &vNewPoint, 
                                dIonSize2 + (iIon1 ? dIonSize1 : dIonSize2) );
                OctTreeUpdateCharge( octTreeSolute, &vNewPoint, 
                        (float)dICharge2, (iIon1 ? dIonSize1 : dIonSize2),
                        &vMaxPot, &vMinPot );
                iIon2--;
        }
    }
/*
OctTreePrintGrid( octTreeSolute, "Charge2", COLOR_RANGE );
*/
    VP0(( "\nDone adding ions.\n" ));
    OctTreeDestroy( &octTreeSolute );
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( (CONTAINER) uUnit );
    return(NULL);

CANCEL:

    VP0(( "\n%s: Interrupted.\n", sCmd ));
    BasicsResetInterrupt();
    if ( octTreeSolute )
        OctTreeDestroy( &octTreeSolute );    
    DisplayerReleaseUpdates();
    return(NULL);
}

OBJEKT
oCmd_addIonSolv( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit=NULL, uIon1=NULL, uIon2=NULL;
int             iIon1=0, iIon2=0;
double          dCharge, dPertCharge, dICharge1, dICharge2;
double          dIonSize1, dIonSize2;
int             i, iUnknown, ierr, iMinPotRes, iMaxPotRes, iReplace;
HELP            hTemp;
LOOP            lAtoms;
ATOM            aAtom;
VARARRAY        vaSolvent;
char            *sCmd = "addIonSolv";

    BasicsResetInterrupt();
    
    /*
     *  Test args
     */
    ierr = 0;
    switch( iArgCount ) {
        case 3: 
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" ))
                ierr++;
          break;
        case 5: 
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n" )) {
                ierr++;
                break;
          }
          /*
           *  Translate the 2 extra args
           */
          uIon2 = (UNIT)oAssocObject( aaArgs[3] );
          iIon2 = (int)dODouble( oAssocObject( aaArgs[4] )); 
          if ( uIon2  &&  iIon2 == 0 ) {
              VP0(( "%s: '0' is not allowed as the value for the second ion.\n",
                                                sCmd ));
              ierr++;
          }
          break;
        default:
          ierr++;
          break;
    } /* end of switch */
    
    if ( ierr ) {
          hTemp = hHelp( "addionsolv" );
          if ( hTemp == NULL ) {
                VP0(( "No help available on addIonSolv\n" ));
          } else {
                VP0(( "\n%s\n", sHelpText(hTemp) ));
          }
          return(NULL);
    }

    /*
     *  Translate args common to both cases
     */
    uUnit = (UNIT)oAssocObject( aaArgs[0] );
    uIon1 = (UNIT)oAssocObject( aaArgs[1] );
    iIon1 = (int)dODouble( oAssocObject( aaArgs[2] ));

    /*
     *  make sure unit has solvent
     */
    vaSolvent = vaSolventResidues( uUnit );
    if ( vaSolvent == NULL ) {
        VP0(( "No solvent present: solvate 1st or use addIons\n" ));
        return(NULL);
    }

    /*
     *  Consider target unit's charge
     */
    ContainerTotalCharge( (CONTAINER) uUnit, &dCharge, &dPertCharge );
    if ( !dCharge ) {
        VP0(( "%s has a charge of 0.\n", sAssocName( aaArgs[1] )));
        if ( iIon1 == 0 ) {
                VP0(( "%s: Can't neutralize.\n", sCmd ));
                return(NULL);
        }
        VP0(( "Adding the ions anyway.\n"));
    } else
        MESSAGE(( "dCharge:  %4.2lf\n", dCharge ));
    
    /*
     *  Consider ion(s) charge
     */
    ContainerTotalCharge((CONTAINER)uIon1, &dICharge1, &dPertCharge );
    if ( !dICharge1 ) {
        VP0(( "%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sAssocName( aaArgs[1] )));
        return(NULL);
    }
    if ( uIon2 ) {
        ContainerTotalCharge((CONTAINER)uIon2, &dICharge2, &dPertCharge );
        if ( !dICharge2 ) {
           VP0(( "%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sAssocName( aaArgs[3] )));
           return(NULL);
        }
    }

    /*
     *  Consider neutralization
     */
    if ( iIon1 == 0 ) {
        if ( (dICharge1 < 0  &&  dCharge < 0) ||
             (dICharge1 > 0  &&  dCharge > 0)) {
                VP0(( "%s: 1st Ion & target are the same charge:\n", sCmd ));
                VP0(( "     can't neutralize.\n" ));
                return(NULL);
        }
        /*
         *  Get the nearest integer number of ions that
         *      we need to add to get as close as possible
         *      to neutral
         */
        iIon1 = (int)lrint( fabs(dCharge) / fabs(dICharge1) );
        if ( iIon1 == 0 )
            VP0(( " %f %d %d %d\n", fabs( dCharge),
                (int)fabs( dCharge), (int)fabs( dICharge1 ),
                (int)(fabs( dCharge) / fabs( dICharge1 )) ));
        if ( uIon2 ) {
                VP0(( "%s: Neutralization - can't do 2nd ion.\n", sCmd ));
                return(NULL);
        }
        VP0(( "%d %s ion%s required to neutralize.\n", iIon1,
                sAssocName( aaArgs[1] ), (iIon1 > 1 ? "s" : "") ));
    } 

    /*
     *  Consider ion sizes and positions.
     */
    dIonSize1 = 0.0;
    iUnknown = 0;
    lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
    for(i=0; (aAtom = (ATOM)oNext(&lAtoms)); i++) {
        if ( iAtomSetTmpRadius( aAtom ) )
                VP0(( "Using default radius %5.2f for ion %s\n",
                        ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[1] ) ));
        dIonSize1 = MAX( dIonSize1, dAtomTemp( aAtom ) );
        if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                iUnknown++;
    }
    if ( i > 1 ) {
        if ( iUnknown ) {
            VP0(( "Ion %s is polyatomic and has %d atoms w/ no position\n",
                sAssocName( aaArgs[1] ), iUnknown ));
            return( NULL );
        }
        VP0(( "Ion %s is polyatomic; multiplying max radius %5.2f by # atoms\n",
                sAssocName( aaArgs[1] ), dIonSize1 ));
        dIonSize1 *= (double) i;
    } else if ( iUnknown ) {
        VECTOR          vPos;

        VectorDef( &vPos, 0.0, 0.0, 0.0 );
        lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
        aAtom = (ATOM)oNext(&lAtoms);
        AtomSetPosition( aAtom, vPos );
        AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
    }
    if ( uIon2 ) {
        dIonSize2 = 0.0;
        iUnknown = 0;
        lAtoms = lLoop( (OBJEKT)uIon2, ATOMS );
        for(i=0; (aAtom = (ATOM)oNext(&lAtoms)); i++) {
                if ( iAtomSetTmpRadius( aAtom ) )
                        VP0(( "Using default radius %5.2f for ion %s\n",
                                ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[3] ) ));
                dIonSize2 = MAX( dIonSize2, dAtomTemp( aAtom ) );
                if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                    iUnknown++;
        }
        if ( i > 1 ) {
            if ( iUnknown ) {
                VP0(( "Ion %s is polyatomic and has %d atoms w/ no position\n",
                    sAssocName( aaArgs[1] ), iUnknown ));
                return( NULL );
            }
            VP0(( 
                "Ion %s is polyatomic; multiplying max radius %5.2f by # atoms",
                sAssocName( aaArgs[3] ), dIonSize2 ));
            dIonSize2 *= (double) i;
        } else if ( iUnknown ) {
            VECTOR              vPos;

            VectorDef( &vPos, 0.0, 0.0, 0.0 );
            lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
            aAtom = (ATOM)oNext(&lAtoms);
            AtomSetPosition( aAtom, vPos );
            AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
        }
    }

    VP0(( "Adding %d counter ions to \"%s\", substituting solvent\n", 
                        iIon1 + iIon2, sAssocName( aaArgs[0] )));

    if ( iIon1 + iIon2 == 0 ) {
        VarArrayDestroy( &vaSolvent );
        return(NULL);
    }

    if ( iIon1 + iIon2 > iVarArrayElementCount( vaSolvent ) ) {
        VP0(( "Can't do it - more ions than solvent\n" ));
        VarArrayDestroy( &vaSolvent );
        return(NULL);
    }

    /*
     *  calc potential on solvent centers.
     */
    VP0(("calculating initial potential at 1st atom in each solvent res..\n"));
    ToolInitSolventPotential( uUnit, vaSolvent, &iMinPotRes, &iMaxPotRes );

    VP0(( "placing ions..\n" ));
    TurnOffDisplayerUpdates();
    while ( iIon1 || iIon2 ) {
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iIon1 ) {
                if ( dICharge1 < 0 )
                        iReplace = iMaxPotRes;
                else
                        iReplace = iMinPotRes;
                ToolReplaceSolvent( uUnit, vaSolvent, 
                                        iReplace, uIon1, dICharge1,
                                        &iMinPotRes, &iMaxPotRes );
                VP0(( "Placed %s in %s.\n", 
                        sAssocName( aaArgs[1] ), sAssocName( aaArgs[0] ) ));

                iIon1--;
        }
        if ( iIon2 ) {
                if ( dICharge2 < 0 )
                        iReplace = iMaxPotRes;
                else
                        iReplace = iMinPotRes;

                ToolReplaceSolvent( uUnit, vaSolvent, 
                                        iReplace, uIon2, dICharge2,
                                        &iMinPotRes, &iMaxPotRes );
                VP0(( "Placed %s in %s.\n", 
                        sAssocName( aaArgs[3] ), sAssocName( aaArgs[0] ) ));
                iIon2--;
        }
    }
    VP0(( "\nDone adding ions.\n" ));
    VarArrayDestroy( &vaSolvent );
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( (CONTAINER) uUnit );
    return(NULL);

CANCEL:

    VP0(( "\n%s: Interrupted.\n", sCmd ));
    BasicsResetInterrupt();
    VarArrayDestroy( &vaSolvent );
    DisplayerReleaseUpdates();
    return(NULL);
}

/*
 *      oCmd_addIonsNear
 *
 */
OBJEKT
oCmd_addIonsNear( int iArgCount, ASSOC aaArgs[] )
{
        VP0(( "Not implemented\n"));
        return(NULL);
}

/*
 *     oCmd_addIonsRand
 * 
 *     Robin Betz (2011)
 */
OBJEKT
oCmd_addIonsRand( int iArgCount, ASSOC aaArgs[] )
{
  UNIT            uUnit=NULL, uIon1=NULL, uIon2=NULL, uPlace=NULL;
  int             iIon1=0, iIon2=0;
  double          dCharge, dPertCharge, dICharge1, dICharge2;
  double          dMinSeparation = 0.0;
  int             i, iUnknown, ierr, random;
  VECTOR          vNewPoint;
  HELP            hTemp;
  RESIDUE         *rPRes;
  LOOP            lAtoms;
  ATOM            aAtom;
  VARARRAY        vaSolvent = NULL;
  ATOM*           aIons = NULL;
  BOOL            bPlaceIon;
  int             counter = 0;
  int             iFailCounter = 0;
  double dIonDist;
  
  VECTOR vPoint;
  char            *sCmd = "addIonsRand";
  
  // Setup
  BasicsResetInterrupt();
  srand(time(NULL));
  
  // Test arguments
  ierr = 0;
  switch( iArgCount )
  {
    case 3:  // One ion and desired number / charge
      if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" ))
        ++ierr;
        break;
    case 4:  // One ion and desired number / charge and minimum separation
      if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n n" ))
      {
        ++ierr;
        break;
      }
      // Get the minimum separation
      dMinSeparation = dODouble( oAssocObject( aaArgs[3] ));
      if (dMinSeparation < 0.0)
      {
        VP0(( "%s: %d is not a valid minimum distance between ions.\n",
              sCmd, dMinSeparation ));
        ++ierr;
      }
        break;  
    case 5: // Two ions and number of each of them
      if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n" ))
      {
        ++ierr;
        break;
      }
      // Get the arguments for the second ion
      uIon2 = (UNIT)oAssocObject( aaArgs[3] );
      iIon2 = (int)dODouble( oAssocObject( aaArgs[4] )); 
      if ( uIon2  &&  iIon2 == 0 )
      {
        VP0(( "%s: '0' is not allowed as the value for the second ion.\n",
              sCmd ));
        ++ierr;
      }
      break;
      
    case 6: // Two ions and number of each of them and minimum separation
      if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n n" ))
      {
        ++ierr;
        break;
      }
      // Get the arguments for the second ion
      uIon2 = (UNIT)oAssocObject( aaArgs[3] );
      iIon2 = (int)dODouble( oAssocObject( aaArgs[4] )); 
      if ( uIon2  &&  iIon2 == 0 )
      {
        VP0(( "%s: '0' is not allowed as the value for the second ion.\n",
              sCmd ));
        ++ierr;
      }
      // Get the minimum separation
      dMinSeparation = dODouble( oAssocObject( aaArgs[3] ));
      if (dMinSeparation < 0.0)
      {
        VP0(( "%s: %d is not a valid minimum distance between ions.\n",
              sCmd, dMinSeparation ));
        ++ierr;
      }
      break;
      
    default:
      ++ierr;
      break;
  }
  
  // Display help if command is malformed
  if ( ierr )
  {
    hTemp = hHelp( "addionsrand" );
    if ( hTemp == NULL )
      VP0(( "No help available on addIons\n" ));
    else
      VP0(( "\n%s\n", sHelpText(hTemp) ));
    return(NULL);
  }
  
  // Translate unit, ion, and charge arguments
  uUnit = (UNIT)oAssocObject( aaArgs[0] );
  uIon1 = (UNIT)oAssocObject( aaArgs[1] );
  iIon1 = (int)dODouble( oAssocObject( aaArgs[2] ));
  
  // Check the unit's validity
  ContainerTotalCharge( (CONTAINER) uUnit, &dCharge, &dPertCharge );
  if ( !dCharge ) 
  {
    VP0(( "%s has a charge of 0.\n", sAssocName( aaArgs[1] )));
    if ( iIon1 == 0 )
    {
      VP0(( "%s: Can't neutralize.\n", sCmd ));
      return(NULL);
    }
    VP0(( "Adding the ions anyway.\n"));
  } 
  else
    MESSAGE(( "dCharge:  %4.2lf\n", dCharge ));
  
  // Make sure the ions are actually ions
  ContainerTotalCharge((CONTAINER)uIon1, &dICharge1, &dPertCharge );
  if ( !dICharge1 ) 
  {
    VP0(( "%s: %s is not an ion and is not appropriate for placement.\n",
          sCmd, sAssocName( aaArgs[1] )));
    return(NULL);
  }
  if ( uIon2 ) 
  {
    ContainerTotalCharge((CONTAINER)uIon2, &dICharge2, &dPertCharge );
    if ( !dICharge2 ) 
    {
      VP0(( "%s: %s is not an ion and is not appropriate for placement.\n",
            sCmd, sAssocName( aaArgs[3] )));
      return(NULL);
    }
  }
  
  // Check validity of neutralization
  if ( iIon1 == 0 ) 
  {
    if ( (dICharge1 < 0  &&  dCharge < 0) ||
      (dICharge1 > 0  &&  dCharge > 0)) {
      VP0(( "%s: 1st Ion & target are the same charge:\n", sCmd ));
      VP0(( "     can't neutralize.\n" ));
      return(NULL);
    }
    /*
     *  Get the nearest integer number of ions that
     *      we need to add to get as close as possible
     *      to neutral
     */
    iIon1 = (int)lrint( fabs(dCharge) / fabs(dICharge1) );
    if ( iIon1 == 0 ) {
      VP0(( " %f %d %d %d\n", fabs( dCharge),
            (int)fabs( dCharge), (int)fabs( dICharge1 ),
            (int)(fabs( dCharge) / fabs( dICharge1 )) ));
    }
    if ( uIon2 ) {
        VP0(( "%s: Neutralization - can't do 2nd ion.\n", sCmd ));
        return(NULL);
    }
    VP0(( "%d %s ion%s required to neutralize.\n", iIon1,
          sAssocName( aaArgs[1] ), (iIon1 > 1 ? "s" : "") ));
  } 
  
  // Check ion size and position
  iUnknown = 0;
  lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
  for(i=0; aAtom = (ATOM)oNext(&lAtoms); i++) {
    if ( iAtomSetTmpRadius( aAtom ) )
      VP0(( "Using default radius %5.2f for ion %s\n",
            ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[1] ) ));
    if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
      iUnknown++;
  }
  if ( i > 1 ) {
    if ( iUnknown ) {
      VP0(( "Ion %s is polyatomic and has %d atoms w/ no position\n",
            sAssocName( aaArgs[1] ), iUnknown ));
      return( NULL );
    }
  } else if ( iUnknown ) {
    VECTOR          vPos;
    
    VectorDef( &vPos, 0.0, 0.0, 0.0 );
    lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
    aAtom = (ATOM)oNext(&lAtoms);
    AtomSetPosition( aAtom, vPos );
    AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
  }
  if ( uIon2 ) {
    iUnknown = 0;
    lAtoms = lLoop( (OBJEKT)uIon2, ATOMS );
    for(i=0; aAtom = (ATOM)oNext(&lAtoms); i++) {
      if ( iAtomSetTmpRadius( aAtom ) )
        VP0(( "Using default radius %5.2f for ion %s\n",
              ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[3] ) ));
      if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
        iUnknown++;
    }
    if ( i > 1 ) {
      if ( iUnknown ) {
        VP0(( "Ion %s is polyatomic and has %d atoms w/ no position\n",
              sAssocName( aaArgs[1] ), iUnknown ));
        return( NULL );
      }
    } else if ( iUnknown ) {
      VECTOR              vPos;
      VectorDef( &vPos, 0.0, 0.0, 0.0 );
      lAtoms = lLoop( (OBJEKT)uIon1, ATOMS );
      aAtom = (ATOM)oNext(&lAtoms);
      AtomSetPosition( aAtom, vPos );
      AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
    }
  }
  
  vaSolvent = vaSolventResidues( uUnit );
  if ( !vaSolvent )
  {
    VP0(( "No solvent present. Add solvent first.\n"));
    return(NULL);
  }
  if ( iIon1 + iIon2 == 0 )
    return(NULL);
  if ( iVarArrayElementCount(vaSolvent)-iIon1-iIon2 <= 0)
  {
    VP0(( "Too few solvent molecules to add ions.\n" ));
    return(NULL);
  }
  VP0(( "Adding %d counter ions to \"%s\". %d solvent molecules will remain.\n", 
        iIon1 + iIon2, sAssocName( aaArgs[0] ), iVarArrayElementCount(vaSolvent)-iIon1-iIon2));
  
  TurnOffDisplayerUpdates();
  if( dMinSeparation )
    MALLOC(aIons, ATOM*, (iIon1+iIon2)*sizeof(ATOM));
  
  // now actually add the ions
  while ( iIon1 || iIon2 ) 
  {
    if ( bBasicsInterrupt() ) goto CANCEL;
    if ( iIon1 ) 
    {
      // Pick random solvent molecule to replace
      random = rand() % iVarArrayElementCount( vaSolvent );
      
      // Get position of solvent residue atom
      rPRes = (RESIDUE*)PVarArrayIndex ( vaSolvent, random );
      lAtoms = lLoop( (OBJEKT)*rPRes, ATOMS);
      aAtom = (ATOM)oNext(&lAtoms);
      vNewPoint = vAtomPosition( aAtom );
      
      // Check that new point isn't too close to other ions
      bPlaceIon = TRUE;
      for (i=0; i<counter; ++i)
      {
        vPoint = vAtomPosition( aIons[i] );
        vPoint = vVectorSub(&vPoint, &vNewPoint);
        dIonDist = dVectorLen(&vPoint);
        if (dIonDist < dMinSeparation)
        {
          ++iFailCounter;
          bPlaceIon = FALSE;
          break;
        }
      }
      
      if ( bPlaceIon )
      {
        VP0(( "%d: Placed %s in %s at (%4.2lf, %4.2lf, %4.2lf).\n", counter,
              sAssocName( aaArgs[1] ), sAssocName( aaArgs[0] ),
              dVX(&(vNewPoint)), 
              dVY(&(vNewPoint)), 
              dVZ(&(vNewPoint))));
        
        // Save this ion's position if desired
        uPlace = (UNIT) oCopy( (OBJEKT)uIon1 );
        ContainerCenterAt((CONTAINER) uPlace, vNewPoint );
        if ( dMinSeparation ) {
          lAtoms = lLoop( (OBJEKT)(uPlace), ATOMS);
          aIons[counter] = (ATOM)oNext(&lAtoms);
          ++counter;
        }
        // Copy ion unit, position, and add it to the unit
        UnitJoin( uUnit, uPlace );
        
        // Delete the solvent residue that was replaced
        REF( *rPRes );  /* bContainerRemove() needs this */
        ResidueYouAreBeingRemoved( *rPRes );
        if ( bContainerRemove( (CONTAINER)uUnit, (OBJEKT)*rPRes ) == FALSE)
          DFATAL(( "rmv solv %d failed\n", random ));
        ContainerDestroy((CONTAINER *) rPRes );
        rPRes = NULL;
        VarArrayDelete(vaSolvent, random);
        
        --iIon1;
      }
      if ( iFailCounter > 100 )
      {
        VarArrayDelete(vaSolvent, random);
        FREE( aIons );
        DFATAL(( "Impossible to place %d ions with minimum separation of %f.\n",
                 iIon1 + iIon2, dMinSeparation ));
      }
    }
    if ( iIon2 ) 
    {
      // Pick random solvent molecule to replace
      random = rand() % iVarArrayElementCount( vaSolvent );
      
      // Get position of solvent residue atom
      rPRes = PVAI( vaSolvent, RESIDUE, random );
      lAtoms = lLoop( (OBJEKT)*rPRes, ATOMS);
      aAtom = (ATOM)oNext(&lAtoms);
      vNewPoint = vAtomPosition( aAtom );
      
      // Check that new point isn't too close to other ions
      bPlaceIon = TRUE;
      for (i=0; i<counter; ++i)
      {
        vPoint = vAtomPosition( aIons[i] );
        vPoint = vVectorSub(&vPoint, &vNewPoint);
        dIonDist = dVectorLen(&vPoint);
        if (dIonDist < dMinSeparation)
        {
          ++iFailCounter;
          bPlaceIon = FALSE;
          break;
        }
      }
      
      if ( bPlaceIon )
      {
        VP0(( "Placed %s in %s at (%4.2lf, %4.2lf, %4.2lf).\n", 
              sAssocName( aaArgs[3] ), sAssocName( aaArgs[0] ),
              dVX(&(vNewPoint)), 
              dVY(&(vNewPoint)), 
              dVZ(&(vNewPoint))));
        
        // Save this ion's position
        uPlace = (UNIT) oCopy( (OBJEKT)uIon2 );
        ContainerCenterAt((CONTAINER) uPlace, vNewPoint );
        if (dMinSeparation) {
        lAtoms = lLoop( (OBJEKT)(uPlace), ATOMS);
        aIons[counter] = (ATOM)oNext(&lAtoms);
        ++counter;
        }
        // Copy ion unit, position, and add it to the unit
        UnitJoin( uUnit, uPlace );
        
        // Delete the solvent residue that was replaced
        REF( *rPRes );  /* bContainerRemove() needs this */
        ResidueYouAreBeingRemoved( *rPRes );
        if ( bContainerRemove( (CONTAINER)uUnit, (OBJEKT)*rPRes ) == FALSE)
          DFATAL(( "rmv solv %d failed\n", random ));
        ContainerDestroy((CONTAINER *) rPRes );
        rPRes = NULL;
        VarArrayDelete(vaSolvent, random);
        
        iIon2--;
      }
      if ( iFailCounter > 100 )
      {
        VarArrayDelete(vaSolvent, random);
        FREE( aIons );
        DFATAL(( "Impossible to place %d ions with minimum separation of %f.\n",
                 iIon1 + iIon2, dMinSeparation ));
      }
    }
  }
  
  // cleanup
  if ( vaSolvent )
    VarArrayDestroy( &vaSolvent );
  if ( aIons )
    FREE( aIons );
  TurnOnDisplayerUpdates();
  ContainerDisplayerUpdate( (CONTAINER) uUnit );
  return(NULL);
  
  // Error handling
  CANCEL:
  
  VP0(( "\n%s: Interrupted.\n", sCmd ));
  BasicsResetInterrupt();  
  if ( vaSolvent )
    VarArrayDestroy( &vaSolvent );
  DisplayerReleaseUpdates();
  return(NULL);
}


/*
 *      oCmd_alias
 *
 *      Author:  David A. Rivkin (1992)
 *
 *      Add or remove an entry to or list entries in the Alias table.
 *      If both strings are present, then add the alias to the table.
 *      If only the first one string is there then remove the alias.
 *      If no arguments are given, then list the current aliases.
 *
 *      alias [alias[, string]]
 *
 *      Arguments:
 *      optional[0] - OSTRING alias is the alias to use
 *      optional[1] - OSTRING string is the original command word
 *
 */
 
VARARRAY GvaAlias;


OBJEKT
oCmd_alias( int iArgCount, ASSOC aaArgs[] )
{
STRING  sCommand;
STRING  sAlias;
ALIASt  aAlias, *PaAlias;
BOOL    bOK;
int     iAliases, i;
STRING  sPossible;
HELP    hTemp;
char    *sCmd = "alias";

    switch ( iArgCount) {
        case 2 : /* Adds an alias */
                MESSAGE(( "Test:  %x, %s\n", 
                        oAssocObject(aaArgs[0]), 
                        sObjectType( oAssocObject(aaArgs[1]))));
                if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s s" )) {
                    strcpy( sAlias, sOString( oAssocObject( aaArgs[0] )));
                    strcpy( sCommand, sOString( oAssocObject( aaArgs[1] )));
                } else {
                    hTemp = hHelp( "alias" );
                    if ( hTemp == NULL ) {
                        VP0(( "No help available on \"alias\".\n" ));
                    } else {
                        VP0(( "\n" ));
                        VP0(( "%s\n", sHelpText(hTemp) ));
                    }
                    return(NULL);
                }
                StringLower( sAlias );
                StringLower( sCommand );
                bOK = FALSE;

                /* Make sure that sCommand is a command */
                for ( i=0; strlen(cCommands[i].sName) != 0; i++ ) {
                    strcpy( sPossible, cCommands[i].sName );
                    StringLower( sPossible );
                        if ( strcmp( sCommand, sPossible ) == 0 ) {
                            bOK = TRUE;
                            break;
                        }
                }
                if ( bOK == FALSE ) {
                    VP0(( "%s: '%s' is not a command.\n", sCmd, sCommand ));
                    VP0(( "Please check the spelling and try again.\n" ));
                    return( NULL );
                }
        
                /* Make sure that the alias is not an existing command */
                /* Make sure that sAlias is not a command */
                for ( i=0; strlen(cCommands[i].sName) != 0; i++ ) {
                    strcpy( sPossible, cCommands[i].sName );
                    StringLower( sPossible );
                    if ( strcmp( sAlias, sPossible ) == 0 ) {
                        bOK = FALSE;
                        break;
                    }
                }
                if ( bOK == FALSE ) {
                    VP0(( "%s: '%s' is already one of the commands.\n", 
                                                sCmd, sAlias ));
                    VP0(( "Please try something different.\n" ));
                    return( NULL );
                }       

                if ( GvaAlias == 0 ) {
                    MESSAGE(( "Creating a new Global Alias Structure.\n" ));
                    GvaAlias = vaVarArrayCreate( sizeof( ALIASt ));
                }
                /* 
                 *  Make sure that the alias does not already exist, 
                 *      if so, replace it 
                 */
                bOK = FALSE;
                iAliases = iVarArrayElementCount( GvaAlias );
                if ( iAliases ) {
                    PaAlias = PVAI( GvaAlias, ALIASt, 0 );
                    for ( i = 0; i < iAliases; i++, PaAlias++ ) {
                        if ( strcmp( sAlias, PaAlias->sName ) == 0 ) {
                            bOK = TRUE;
                            strcpy( PaAlias->sCommand, sCommand );
                            break;
                        }
                    }
                }
                if ( bOK == FALSE ) {
                    memset(&aAlias, 0, sizeof(aAlias)); /* for Purify */
                    strcpy( aAlias.sName, sAlias );
                    strcpy( aAlias.sCommand, sCommand );
                    VarArrayAdd( GvaAlias, (GENP)&aAlias );
                }
                break;
                
        case 1 : /* Remove an alias from the list */
                if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s" )) {
                    strcpy( sAlias, sOString( oAssocObject( aaArgs[0] )));
                } else if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "z" )) {
                    strcpy( sAlias, sAssocName( aaArgs[0] ));
                } else {
                    hTemp = hHelp( "alias" );
                    if ( hTemp == NULL ) {
                        VP0(( "No help available on \"alias\".\n" ));
                    } else {
                        VP0(( "\n" ));
                        VP0(( "%s\n", sHelpText(hTemp) ));
                    }
                    return(NULL);
                }
                StringLower( sAlias );
                iAliases = iVarArrayElementCount( GvaAlias );
                if ( !iAliases ) {
                    VP0(( "%s: There are no aliases loaded.\n", sCmd ));
                    return(NULL);
                }
                PaAlias = PVAI( GvaAlias, ALIASt, 0 );
                for ( i = 0; i < iAliases; i++, PaAlias++ ) {
                    if ( strcmp( sAlias, PaAlias->sName ) == 0 ) {
                        VarArrayDelete( GvaAlias, i );
                        break;
                    }
                }
                break;
        case 0 : /*  List all the aliases */
                iAliases = iVarArrayElementCount( GvaAlias );
                if ( !iAliases ) {
                    VP0(( "There are no aliases loaded.\n" ));
                    return(NULL);
                }
                VP0(( "Current Aliases  [alias....command]\n" ));
                for ( i = 0; i < iAliases; i++ ) {
                    aAlias = *PVAI( GvaAlias, ALIASt, i );
                    if (( i % 2 ) == 0 ) { /* An odd entry */
                        VP0(( "%-10s..%-24s", 
                                aAlias.sName, aAlias.sCommand ));
                    } else {
                        VP0(( "%-10s..%-10s\n", 
                                aAlias.sName, aAlias.sCommand ));
                    }
                }
                if (( i % 2 ) == 1 ) { 
                    /* Left over entry (odd number of entries) */
                    VP0(( "\n" ));
                }
                break;
                
        default: hTemp = hHelp( "alias" );
                if ( hTemp == NULL ) {
                    VP0(( "No help available on alias.\n" ));
                } else {
                    VP0(( "\n" ));
                    VP0(( "%s\n", sHelpText(hTemp) ));
                }
                break;
                
   } /* end of switch */
   return(NULL);
}



/*
 *      oCmd_update
 *
 *      Author: David A. Rivkin (1993)
 *
 *      Scans a UNIT for any residue that is different from 
 *      the residues in the UNIT and changes them to that residue.
 *      OR
 *      Checks a particular residue in a unit to see if the template
 *      residues have changed and modifies it as needed.
 *
 *      update unit/residue
 *              ARG[0] - UNIT/RESIDUE
 *
 */
OBJEKT
oCmd_update( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;
RESIDUE         rNew, rOld, rCopy, rTemp;
int             iResidues;
LOOP            lRes, lResidue;
int             i, iNext;
HELP            hTemp;
STRING          sResName, sTemp;
DICTLOOP        dlLoop;
char            *sCmd = "update";

    if ( iArgCount != 1 ) {
        hTemp = hHelp( "update" );
        if ( hTemp == NULL ) {
            VP0(( "No help available on update.\n" ));
        } else {
            VP0(( "\n" ));
            VP0(( "%s\n", sHelpText(hTemp) ));
        }
        return(NULL);
    }
    if (! bCmdGoodArguments( sCmd, iArgCount, aaArgs, "ru" )) {
        hTemp = hHelp( "update" );
        if ( hTemp == NULL ) {
            VP0(( "No help available on update.\n" ));
        } else {
            VP0(( "\n" ));
            VP0(( "%s\n", sHelpText(hTemp) ));
        }
        return(NULL);
    }
    switch ( iObjectType(oAssocObject( aaArgs[0] ))) {
        case    UNITid: 
            DisplayerAccumulateUpdates();
            uUnit = (UNIT)oAssocObject( aaArgs[0] );
            lResidue = lLoop( (OBJEKT)uUnit, RESIDUES );
            while ( (rOld = (RESIDUE)oNext( &lResidue )) ) {
                strcpy( sTemp, sContainerName((CONTAINER)rOld));
                for ( i = 0; i < strlen(sTemp); i++){
                    if ( sTemp[i] == ' ') break;
                    sResName[i] = sTemp[i];
                }
                sResName[i] = '\0';
                dlLoop = ydlDictionaryLoop(GdVariables);
                while ( yPDictionaryNext( GdVariables, &dlLoop )) {
                    if (!strcmp(sResName, sDictLoopKey(dlLoop))) {
                        iResidues = 0;
                        lRes = lLoop( (OBJEKT)PDictLoopData( dlLoop ), RESIDUES);
                        while ((rNew = (RESIDUE)oNext( &lRes ))) {
                            iResidues++;
                                /* Make sure that there are not more than 
                                        1 residue in the unit */
                            if ( iResidues > 1 ) {
                                MESSAGE(( 
                             "Unit name matches but has more than 1 residue.\n"
                                ));
                                goto NEXTRES1;
                            }
                            rCopy = (RESIDUE)oCopy((OBJEKT)rNew );
                        }
                        ResidueMutate( rCopy, rOld );

                        /* Remove the old RESIDUE from the UNIT */
                        /* And add the new RESIDUE, without incrementing */
                        /* the UNITs next child sequence number */

                        REF( rOld );  /* bContainerRemove() needs this */
                        bContainerRemove( (CONTAINER) uUnit, (OBJEKT)rOld );
                        iNext = iContainerNextChildsSequence( (CONTAINER) uUnit );
                        ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rCopy );
                        ContainerSetSequence( (CONTAINER) rCopy, iContainerSequence((CONTAINER) rOld) );
                        ContainerSetNextChildsSequence( (CONTAINER) uUnit, iNext );

                        DEREF( rOld );

                        VP0(( "Updating residue %s.\n", sResName ));
                        goto NEXTSEQ;
                    }
NEXTRES1:       ;
                }
NEXTSEQ:        ;
            }
            DisplayerReleaseUpdates();
            break;
            
        case    RESIDUEid: 
            DisplayerAccumulateUpdates();
            rOld = (RESIDUE)oAssocObject( aaArgs[0] );
            uUnit = (UNIT)cContainerWithin((CONTAINER)rOld );
            strcpy( sTemp, sContainerName((CONTAINER)rOld )); 
            for ( i = 0; i < strlen(sTemp); i++){
                if ( sTemp[i] == ' ') break;
                    sResName[i] = sTemp[i];
            }
            sResName[i] = '\0';                 
            
            dlLoop = ydlDictionaryLoop(GdVariables);
            while ( yPDictionaryNext( GdVariables, &dlLoop )) {
                if (!strcmp(sResName, sDictLoopKey(dlLoop))) {
                    iResidues = 0;
                    lRes = lLoop( (OBJEKT)PDictLoopData( dlLoop ), RESIDUES);
                    while ( (rTemp = (RESIDUE)oNext( &lRes)) ) {
                        rNew = rTemp;
                        iResidues++;
                    }
                        /* Make sure that there are not more than 
                                1 residue in the unit */
                    if ( iResidues > 1 ) {
                        MESSAGE(( 
                          "Unit name matches but has more than 1 residue.\n" ));
                        goto NEXTRES2;
                    }

                    rCopy = (RESIDUE)oCopy( (OBJEKT)rNew );
                    ResidueMutate( rCopy, rOld );

                        /* Remove the old RESIDUE from the UNIT */
                        /* And add the new RESIDUE, without incrementing */
                        /* the UNITs next child sequence number */

                    REF( rOld );  /* bContainerRemove() needs this */
                    bContainerRemove( (CONTAINER) uUnit, (OBJEKT)rOld );
                    iNext = iContainerNextChildsSequence( (CONTAINER) uUnit );
                    ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rCopy );
                    ContainerSetSequence( (CONTAINER) rCopy, iContainerSequence((CONTAINER) rOld) );
                    ContainerSetNextChildsSequence( (CONTAINER) uUnit, iNext );

                    DEREF( rOld );

                    VP0(( "Updating residue %s.\n", sResName ));
                    goto NEXTRES2;
                }
NEXTRES2:       ;
            }
            DisplayerReleaseUpdates();
            break;
        default:        
            hTemp = hHelp( "update" );
            if ( hTemp == NULL ) {
                VP0(( "No help available on update.\n" ));
            } else {
                VP0(( "\n" ));
                VP0(( "%s\n", sHelpText(hTemp) ));
            }
            break;
    }
    return(NULL);
}


/*
 *      oCmd_flip
 *      Based on XAUESelectedAtomsFlipChirality in xaUnitEditor.c (xleap source code)
 *      Author: Christine cezard (2007)
 *      Universite de Picardie - Jules Verne, Amiens 
 *      http://q4md-forcefieldtools.org
 *
 */
OBJEKT
oCmd_flip(int iArgCount, ASSOC aaArgs[])
{
UNIT        uUnit;
LOOP        lAtoms;
ATOM        aAtom;

    uUnit = (UNIT)oAssocObject(aaArgs[0]);
    if ( uUnit == NULL ) return(NULL);

    lAtoms = lLoop((OBJEKT) uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
            bBuildFlipChiralityFor((CONTAINER) uUnit, aAtom );
        }
    }
  

    return(NULL);
}

/*
 *      oCmd_relax
 *      Based on Based on XAUERelaxSelectionInFramework in xaUnitEditor.c (xleap source code)
 *      Author: Christine Cezard (2007)
 *      Universite de Picardie - Jules Verne, Amiens
 *      http://q4md-forcefieldtools.org
 *
 */
OBJEKT
oCmd_relax(int iArgCount, ASSOC aaArgs[])
{
MINIMIZER       mStrain;
UNIT            uUnit;
    
    uUnit = (UNIT)oAssocObject(aaArgs[0]);
    if ( uUnit == NULL ) return(NULL);
    
        /* Setup a MINIMIZER and give it a callback to use */
        /* to update the display every step of the minimization */

    mStrain = mMinimizerCreate();
        /* Set up the MINIMIZER to use, and turn off any */
        /* control-c that may have been hit before */

    BasicsResetInterrupt();
    SelectRelaxInFramework( (UNIT)uUnit, mStrain );

    return(NULL);
}

/*
 *      oCmd_addH
 *      Based on Based on XAUEAddHydrogensBuildExternals in xaUnitEditor.c (xleap source code)
 *      Author: D. Roe (2011)
 *      Rutgers University 
 *
 */
OBJEKT
oCmd_addH(int iArgCount, ASSOC aaArgs[])
{
LOOP            lAtom, lSpan;
ATOM            aAtom, aStart;
UNIT            uUnit;
int             iDum;

    DisplayerAccumulateUpdates();
    uUnit = (UNIT)oAssocObject(aaArgs[0]);
    if ( uUnit == NULL ) return(NULL);

        /* Add hydrogens */

    ModelAddHydrogens( uUnit );

        /* Try to build geometries for simple rings */

    BuildInternalsForSimpleRings( (CONTAINER)uUnit );

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
            iDum = 0;   /* for purify */
            BuildExternalsUsingFlags( &lSpan, ATOMNEEDSBUILD, 0,
                                        ATOMPOSITIONKNOWN,
                                        ATOMNEEDSBUILD,
                                        &iDum, &iDum, &iDum, TRUE );
        }
    }

                /* Destroy all of the INTERNALs */

    lAtom = lLoop((OBJEKT) uUnit, ATOMS );
    BuildDestroyInternals( &lAtom );

    DisplayerReleaseUpdates();

    return(NULL);
}


COMMANDt        cCommands[] = {

        { "add",                oCmd_add },
        { "addH",               oCmd_addH },
        { "addIons",            oCmd_addIons },
        { "addIons2",           oCmd_addIons2 },
        { "addIonSolv",         oCmd_addIonSolv },
        { "addIonsRand",        oCmd_addIonsRand },
        { "addPath",            oCmd_addPath },
        { "addPdbAtomMap",      oCmd_addPdbAtomMap },
        { "addPdbResMap",       oCmd_addPdbResMap },
        { "addAtomTypes",       oCmd_addAtomTypes },
        { "alias",              oCmd_alias },
        { "alignAxes",          oCmd_alignAxes },
        { "bond",               oCmd_bond },
        { "bondByDistance",     oCmd_bondByDistance },
        { "center",             oCmd_center },
        { "charge",             oCmd_charge },
        { "check",              oCmd_check },
        { "clearPdbAtomMap",    oCmd_clearPdbAtomMap },
        { "clearPdbResMap",     oCmd_clearPdbResMap },
        { "clearVariables",     oCmd_clearVariables },
        { "combine",            oCmd_combine },
        { "copy",               oCmd_copy },
        { "createAtom",         oCmd_createAtom },
        { "createParmset",      oCmd_createParmset },
        { "createResidue",      oCmd_createResidue },
        { "createUnit",         oCmd_createUnit },
        { "crossLink",          oCmd_crossLink },
        { "debugOff",           oCmd_debugOff },
        { "debugOn",            oCmd_debugOn },
        { "debugStatus",        oCmd_debugStatus },
        { "deleteBond",         oCmd_deleteBond },
        { "deleteOffLibEntry",  oCmd_deleteOffLibEntry },
        { "deleteRestraint",    oCmd_deleteRestraint },
        { "desc",               oCmd_describe },
        { "deSelect",           oCmd_deSelect },
        { "displayPdbAtomMap",  oCmd_displayPdbAtomMap },
        { "displayPdbResMap",   oCmd_displayPdbResMap },
        { "edit",               oCmd_edit },
        { "flip",               oCmd_flip },
        { "groupSelectedAtoms", oCmd_groupSelectedAtoms },
        { "help",               oCmd_help },
        { "impose",             oCmd_impose },
        { "list",               oCmd_list },
        { "listOff",            oCmd_listOff },
        { "loadAmberParams",    oCmd_loadAmberParams },
        { "loadAmberPrep",      oCmd_loadAmberPrep },
        { "loadOff",            oCmd_loadOff },
        { "loadMol2",           oCmd_loadMol2 },
        { "loadMol3",           oCmd_loadMol3 },
        { "loadPdb",            oCmd_loadPdb },
        { "loadPdbUsingSeq",    oCmd_loadPdbUsingSeq },
        { "logFile",            oCmd_logFile },
        { "matchVariables",     oCmd_matchVariables },
        { "measureGeom",        oCmd_measureGeom },
        { "memDebug",           oCmd_memDebug},
        { "quit",               oCmd_quit },
        { "relax",              oCmd_relax },
        { "remove",             oCmd_remove },
        { "restrainAngle",      oCmd_restrainAngle },
        { "restrainBond",       oCmd_restrainBond },
        { "restrainTorsion",    oCmd_restrainTorsion },
        { "saveAmberParm",      oCmd_saveAmberParm },
        { "saveAmberParmNetCDF",oCmd_saveAmberParmNetCDF },
        { "saveAmberParmPert",  oCmd_saveAmberParmPert },
        { "saveAmberParmPol",   oCmd_saveAmberParmPol },
        { "saveAmberParmPolPert", oCmd_saveAmberParmPolPert },
        { "saveAmberPrep",      oCmd_saveAmberPrep },
        { "saveMol2",           oCmd_saveMol2 },
        { "saveMol3",           oCmd_saveMol3 },
        { "saveOff",            oCmd_saveOff},
        { "saveOffParm",        oCmd_saveOffParm },
        { "savePdb",            oCmd_savePdb },
        { "scaleCharges",       oCmd_scaleCharges },
        { "select",             oCmd_select },
        { "sequence",           oCmd_sequence },
        { "set",                oCmd_set },
        { "setBox",             oCmd_setBox},
        { "showdefault",        oCmd_showDefault},
        { "solvateBox",         oCmd_solvateBox },
        { "solvateCap",         oCmd_solvateCap },
        { "solvateDontClip",    oCmd_solvateDontClip },
        { "solvateOct",         oCmd_solvateOct },
        { "solvateShell",       oCmd_solvateShell },
        { "source",             oCmd_source },
        { "translate",          oCmd_translate },
        { "transform",          oCmd_transform },
/*      { "update",             oCmd_update },          */
        { "verbosity",          oCmd_verbosity },
        { "zMatrix",            oCmd_zMatrix },
/* The last command must be blank */
        { "", NULL } 
};
