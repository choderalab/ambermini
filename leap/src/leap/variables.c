/*
 *      File:   variables.c
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
 *      Maintain the variables dictionary.
 */




#include	"basics.h"

#include        "classes.h"
#include        "dictionary.h"

#include        "leap.h"


DICTIONARY      GdVariables;



/*
 *      VariablesInit
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Define the Variables dictionary.
 */
void
VariablesInit()
{
    GdVariables = dDictionaryCreate();
}




/*
 *      VariablesList
 *
 *	Author:	Christian Schafmeister (1991)
 *	Modified: David A. Rivkin (1992)
 *		Improved Layout to 2 coloumns.
 *
 *      List the variables to stdout
 */
void
VariablesList()
{
DICTLOOP        dlLoop;
int             i = 1;
int             iColumns;

    /* TODO: figure columns for actual window size & max string size */
    if ( GbGraphicalEnvironment == TRUE )
        iColumns = 8;
    else
        iColumns = 8;

    BasicsResetInterrupt();
    dlLoop = ydlDictionaryLoop(GdVariables);
    while ( yPDictionaryNext( GdVariables, &dlLoop ) ) {
	if ( bBasicsInterrupt() ) goto CANCEL;
        if ( i % iColumns ) { 
	    /* columns 0-(end-1) */
            VP0(( "%-10s", sDictLoopKey(dlLoop) ));
        } else {
	   /* last column in line */
            VP0(( "%-10s\n", sDictLoopKey(dlLoop) ));
        }
        i++;
    }
    if ( (i-1) % iColumns )
        VP0(( "\n" ));	/* close final part line */
    return;

CANCEL:
    VP0(( "Interrupt.\n" ));
    BasicsResetInterrupt();
    return;
}



/*
 *      VariableSet
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Define a variable, making sure to REF and DEREF OBJEKTs.
 *
 */
void    VariableSet( char *sName, OBJEKT oObj )
{
OBJEKT  oOld;

    oOld = (OBJEKT)yPDictionaryDelete( GdVariables, sName );
    if ( oOld != NULL ) DEREF(oOld);
    if ( oObj == NULL ) return;
    DictionaryAdd( GdVariables, sName, (GENP)oObj );
    REF( oObj );
}




/*
 *      VariableRemove
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Remove a variable, DEREF'ing the obj.
 */
void
VariableRemove( char *sName )
{
OBJEKT  oOld;

    oOld = (OBJEKT)yPDictionaryDelete( GdVariables, sName );
    if ( oOld != NULL ) 
	DEREF(oOld);
}



/*
 *      oVariable
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the value of a variable, NULL if it does not exist.
 */
OBJEKT
oVariable( char *sName )
{
OBJEKT  oVar;

    oVar = (OBJEKT)yPDictionaryFind( GdVariables, sName );
    return(oVar);
}






/*
 *      VariablesDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy all of the variables.
 */
void
VariablesDestroy()
{
DICTLOOP        dlLoop;

                /* First DEREF the contents of the DICTIONARY */

    dlLoop = ydlDictionaryLoop(GdVariables);
    while ( yPDictionaryNext( GdVariables, &dlLoop ) ) {
        DEREF( PDictLoopData(dlLoop) );
    }

                /* Then Destroy the DICTIONARY */

    DictionaryDestroy( &GdVariables );
}



/*
 *      dVariablesDictionary
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the dictionary containing ALL variables.
 */
DICTIONARY
dVariablesDictionary()
{
    return(GdVariables);
}
