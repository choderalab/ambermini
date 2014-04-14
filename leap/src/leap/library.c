/*
 *      File:   library.c
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
 *              This object does very little except manage the
 *              reading and writing of several units into a 
 *              single DATABASE file.  It maintains an open
 *              DATABASE and a DICTIONARY of UNIT names.
 *              The object attatched to each UNIT name in the
 *              DICTIONARY is a dummy.
 */




#include	"basics.h"

#include        "classes.h"
#include        "dictionary.h"
#include        "database.h"
#include        "parmLib.h"

#include        "library.h"


static  int     SiDummy;




/*
 *	zLibraryRewriteIndex
 *
 *	Rewrite the index.
 */
static void
zLibraryRewriteIndex( LIBRARY ul )
{
VARARRAY	vaIndex;
DICTLOOP	dlLoop;
int		i;

	    /* Rewrite the INDEX! */

    vaIndex = vaVarArrayCreate( sizeof(STRING) );
    VarArraySetSize( vaIndex, iDictionaryElementCount(ul->dLibrary) );
    dlLoop = ydlDictionaryLoop( ul->dLibrary );
    for ( i=0; i<iDictionaryElementCount(ul->dLibrary); i++ ) {
	yPDictionaryNext( ul->dLibrary, &dlLoop );
	strcpy( (char*)(PVarArrayIndex( vaIndex, i )),
		    sDictLoopKey( dlLoop ) );
    }
    DBPutValue( ul->dbLibrary, "!index", ENTRYARRAY | ENTRYSTRING,
		    iDictionaryElementCount(ul->dLibrary),
		    PVarArrayIndex( vaIndex, 0 ), sizeof(STRING) );
    VarArrayDestroy( &vaIndex );

}




/*
 *      lLibraryOpen
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Open the LIBRARY.
 */
LIBRARY
lLibraryOpen( char *sName, int iOpenMode )
{
LIBRARY         ul;
int             iType, iLines;
VARARRAY        vaIndex;
int             i, iLen;
DATABASE	dbDatabase;

    dbDatabase = dbDBRndOpen( sName, iOpenMode );
    if ( !dbDatabase ) {
	if ( iDBLastError() == DB_ERROR_INVALID_FILE ) {
	    VP0(( "Could not open database: %s\n", sName ));
	} else if ( iDBLastError() == DB_ERROR_INVALID_DATABASE ) {
	    VP0(( "File: %s is not a valid database.\n", GsBasicsFullName ));
	}
	return(NULL);
    }
   
    MALLOC( ul, LIBRARY, sizeof(LIBRARYt) );
    
    ul->dLibrary = dDictionaryCreate();
    ul->dbLibrary = dbDatabase;
    
                /* Load in the index */

    if ( bDBGetType( ul->dbLibrary, "!index", &iType, &iLines ) ) {
        vaIndex = vaVarArrayCreate( sizeof(STRING) );
        VarArraySetSize( vaIndex, iLines );
        bDBGetValue( ul->dbLibrary, "!index", &iLen,
			 PVarArrayIndex( vaIndex, 0 ),
                        sizeof(STRING) );

                /* Copy the index into the DICTIONARY */

        for ( i=0; i<iLines; i++ ) {
            DictionaryAdd( ul->dLibrary, PVarArrayIndex( vaIndex, i ),
                                (GENP)&SiDummy );
        }
        VarArrayDestroy( &vaIndex );
    }
    
    return(ul);
}





/*
 *      oLibraryLoad
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Load an OBJEKT from the LIBRARY.
 *      If the unit is not in the DICTIONARY then return NULL.
 */
OBJEKT  
oLibraryLoad( LIBRARY ul, char *sName )
{
STRING          sPrefix;
OBJEKT          oObj;

                /* Check if the UNIT is even in the library */

    if ( yPDictionaryFind( ul->dLibrary, sName ) == NULL ) 
	return(NULL);

    strcpy( sPrefix, "entry." );
    strcat( sPrefix, sName );
    strcat( sPrefix, "." );
    DBPushZeroPrefix( ul->dbLibrary, sPrefix );
    
                /* Load the OBJEKT */

    oObj = (OBJEKT)uUnitLoad( ul->dbLibrary );
    if ( oObj != NULL ) 
	goto DONE;
    oObj = (OBJEKT)psParmSetLoad( ul->dbLibrary );
    if ( oObj == NULL)
	VP0(( " WARNING: %s: nothing loaded\n", sName ));

DONE:
    DBPopPrefix( ul->dbLibrary );
    return(oObj);
}




/*
 *      LibrarySave
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Save an OBJEKT into the LIBRARY.
 *      If the unit is not in the DICTIONARY then return NULL.
 */
void    
LibrarySave( LIBRARY ul, char *sName, OBJEKT oObj, PARMLIB plLibrary )
{
STRING          sPrefix;

                /* If the UNIT is not already in the LIBRARY then */
                /* add it */

    if ( yPDictionaryFind( ul->dLibrary, sName ) == NULL ) {
        DictionaryAdd( ul->dLibrary, sName, (GENP)&SiDummy );

	zLibraryRewriteIndex(ul);
    }

    strcpy( sPrefix, "entry." );
    strcat( sPrefix, sName );
    strcat( sPrefix, "." );
    DBPushZeroPrefix( ul->dbLibrary, sPrefix );
    
                /* Save the UNIT */

    if ( iObjectType(oObj)==UNITid ) 
	UnitSave( (UNIT)oObj, ul->dbLibrary, plLibrary );
    else if ( iObjectType(oObj)==PARMSETid ) 
	ParmSetSave( (PARMSET)oObj, ul->dbLibrary );
    
    DBPopPrefix( ul->dbLibrary );
}




/*
 *      bLibraryRemove
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Remove an entry from the library.
 */
BOOL	
bLibraryRemove( LIBRARY ul, char *sName )
{
STRING		sPrefix, sEntry;

    if ( yPDictionaryFind( ul->dLibrary, sName ) == NULL ) {
	return(FALSE);
    }

    strcpy( sPrefix, "entry." );
    strcat( sPrefix, sName );
    strcat( sPrefix, "." );
    DBZeroPrefix( ul->dbLibrary );

		/* Remove the entry from the dictionary */

    yPDictionaryDelete( ul->dLibrary, sName );
    zLibraryRewriteIndex(ul);

		/* List all of the entries with that prefix */

    DBRndLoopEntryWithPrefix( ul->dbLibrary, sPrefix );
    while ( bDBRndNextEntryWithPrefix( ul->dbLibrary, sEntry ) ) {
	bDBRndDeleteEntry( ul->dbLibrary, sEntry );
    }

    return(TRUE);
}





/*
 *      LibraryLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Initialize the UNITLIBRARY so that subsiquent calls to
 *      sULNext will return UNIT names.
 */
void
LibraryLoop( LIBRARY ul )
{
    ul->dlContents = ydlDictionaryLoop( ul->dLibrary );
}




/*
 *      sLibraryNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return UNIT names that are defined in the UNITLIBRARY
 */
char *
sLibraryNext( LIBRARY ul )
{
    yPDictionaryNext( ul->dLibrary, &(ul->dlContents) );
    if ( ul->dlContents == NULL ) return(NULL);
    return(sDictLoopKey( ul->dlContents ));
}







/*
 *      LibraryClose
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Close the UNITLIBRARY and release all memory used by it.
 */
void
LibraryClose( LIBRARY *lPLibrary )
{

    DBClose( &((*lPLibrary)->dbLibrary) );
    DictionaryDestroy( &((*lPLibrary)->dLibrary) );
    FREE( *lPLibrary );
    *lPLibrary = NULL;
}




/*
 *	LibraryDBPrefix
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the prefix used to access the 
 *	object whose name is (sStr).
 *	Return TRUE if the object is found, otherwise FALSE.
 */
BOOL
bLibraryDBPrefix( LIBRARY lLib, char *cPResult, char *sStr )
{
STRING          sPrefix;

                /* Check if the UNIT is even in the library */

    if ( yPDictionaryFind( lLib->dLibrary, sStr ) == NULL ) 
	return(FALSE);

    strcpy( sPrefix, "entry." );
    strcat( sPrefix, sStr );
    strcat( sPrefix, "." );
    strcpy( cPResult, sPrefix );
    return(TRUE);
}

    



