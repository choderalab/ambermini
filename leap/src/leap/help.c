/*
 *      File:   help.c
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
 *              This object manages the database of help information.
 *              Help is stored as a DICTIONARY of keywords which
 *              are associated with very large strings which contain
 *              the help information for that keyword.
 *              There is only one HELP object per system.
 */








#include	"basics.h"

#include        "dictionary.h"

#include        "help.h"




/*
 *-------------------------------------------------------------------------
 *
 *      Global variables
 *
 *
 */

static  DICTIONARY      SdHelp;         /* Contains the help information */
static  DICTLOOP        SdlLoop;        /* Used to loop through help keywords*/




/*
 *      HelpInitialize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Initialize the HELP database.
 */
void
HelpInitialize()
{

    SdHelp = dDictionaryCreate();

        /* Set up the HELP database */

    HTInit();
}



/*
 *	HelpShutdown
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Free all memory taken up by the help database.
 */
void
HelpShutdown()
{
DICTLOOP	dlEntries;
HELP		hCur;

    dlEntries = ydlDictionaryLoop( SdHelp );
    while ( (hCur = (HELP)yPDictionaryNext( SdHelp, &dlEntries )) != NULL ) {
	FREE( sHelpSubject(hCur) );
	FREE(hCur);
    }

    DictionaryDestroy( &SdHelp );
}





/*
 *      HelpAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a single definition to the HELP database.
 */
void
HelpAdd( char *sUpSubject, char *sText )
{
HELP    hNew;

    MALLOC( hNew, HELP, sizeof(HELPt) );
    sHelpUpSubject( hNew ) = sUpSubject;

    MALLOC( sHelpSubject( hNew ), char *, strlen(sUpSubject) + 1 );
    strcpy( sHelpSubject( hNew ), sUpSubject );

    StringLower( sHelpSubject( hNew ) );
    hNew->sText = sText;
    DictionaryAdd( SdHelp, sHelpSubject( hNew ), (GENP)hNew );
}

/*
 *      hHelp
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return HELP on the desired subject.
 *      Return NULL if there is no help.
 */
HELP
hHelp( char *sSubject )
{
HELP    hTemp;
STRING	sTmp;

    /*
     *  sometimes sSubject may be a hard-coded string, and in some
     *	Unixes these cannot be modified, so make a copy here for
     *	StringLower()
     */
    strcpy( sTmp, sSubject );
    StringLower( sTmp );
    hTemp = (HELP)yPDictionaryFind( SdHelp, sTmp );
    return(hTemp);
}







/*
 *      HelpLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Initialize the loop for all keywords.
 *      Subsequent calls to hHelpNextKeyword will return every HELP
 *      entry in the database.
 */
void
HelpLoop()
{
    SdlLoop = ydlDictionaryLoop(SdHelp);
}




/*
 *      hHelpNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return subsequent help entries.
 */
HELP
hHelpNext()
{
HELP    hNext;

    hNext = (HELP)yPDictionaryNext( SdHelp, &SdlLoop );
    return(hNext);
}

