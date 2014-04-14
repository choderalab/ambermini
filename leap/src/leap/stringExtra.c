/*
 *      File:   stringExtra.c
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
 *              Some extra routines to manipulate strings and characters.
 */


#include	"basics.h" 
 


 
/*
 *      sRemoveSpaces
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Remove all spaces from sIn and place the results in sOut.
 *
 */
char *
sRemoveSpaces( char *sIn, char *sOut )
{
char	*sWrite;

                /* Copy everything except spaces */

    sWrite = sOut;
    while ( (*sIn)!= '\0' ) {
        if ( (*sIn)!=' ' ) (*sWrite++)=(*sIn);
        sIn++;
    }
    (*sWrite) = '\0';
    return(sOut);
}




/*
 *      sRemoveControlAndPadding
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Remove all control characters and all Padding spaces,
 *      spaces at the start and end of the string.
 *
 *TODO This might not be machine independant
 */
char *
sRemoveControlAndPadding( char *sRaw, char *sResult )
{
char	*sCur, *sResultCur;

                /* First skip over the initial spaces and control characters */

    sCur = sRaw;
    while ( (*sCur!='\0') && (*sCur<=' ') ) sCur++;
    if ( *sCur == '\0' ) {
        *sResult = '\0';
        goto DONE;
    }
    
                /* Now copy the rest, minus control characters into sResult */
    sResultCur = sResult;
    while ( *sCur!='\0' ) {
        if ( *sCur >= ' ' ) {
            *sResultCur = *sCur;
            sResultCur++;
        }
        sCur++;
    }
    *sResultCur = '\0';
    
                /* Now remove the trailing spaces */
    if ( strlen(sResult) > 0 ) {
        sResultCur--;
        while ( *sResultCur == ' ' ) sResultCur--;
        sResultCur++;
        *sResultCur = '\0';
    }
    
DONE:
    return(sResult);
}




/*
 *      sRemoveLeadingSpaces
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Just like it says, return a string that has the first
 *      spaces removed.
 */
char *
sRemoveLeadingSpaces( char *sLine )
{
char	*sTemp;
int      mylength;
    sTemp = sLine;
    while ( (*sTemp==' ') && ( *sTemp!='\0' )) sTemp++;
    mylength=strlen(sTemp);
    //strcpy( sLine, sTemp );
    memmove(sLine, sTemp, mylength);
    sLine[mylength]='\0';
    return(sLine);
}




/*
 *      sRemoveFirstString
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Copy everything in the string up until the first space into sHead.
 *      Remove everything (including the space) from the string sLine.
 */
void
sRemoveFirstString( char *sLine, char *sHead )
{
char	*sTemp;
int      mylength;

    sTemp = sLine;
    while ( (*sTemp!=' ') && ( *sTemp!='\0' )) sTemp++;
    if ( *sTemp == '\0' ) {
        strcpy( sHead, sLine );
        *sLine = '\0';
        return;
    }
    *sTemp = '\0';
    strcpy( sHead, sLine );
    sTemp++;
    mylength=strlen(sTemp);
    //strcpy( sLine, sTemp );
    memmove(sLine, sTemp, mylength);
    sLine[mylength]='\0';
    
}










/*
 *      StringCopyMax
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Perform a string copy, but if the source string is longer
 *      than iMax then only copy iMax-2 characters, and place a '\0' 
 *      at the iMax-1'th position in sDest.
 */
void
StringCopyMax( char *sDest, char *sSource, int iMax )
{
    strncpy( sDest, sSource, iMax );
    sDest[iMax-1] = '\0';
}





/*
 *------------------------------------------------------------------
 *
 *      Pattern matching
 *
 *
 */





/*TODO
 *      bStringMatchPattern
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the string matches the pattern in
 *      sPattern.  The string sPattern can contain wildcards
 *      like '*' and '?' along with regular text.
 *      The wildcard '*' matches to anything, while '?'
 *      matches to a single character.
 *      This means you cannot explicitly search for '*' and '?'
 *      characters.
 *
 *      This pattern matcher is very poor, it does not try different
 *      ways of matching the pattern.  It must be rewritten!!!!!!
 *
 *
 */
BOOL
bStringMatchPattern( char *sString, char *sPattern )
{

    while ( 1 ) {
        if ( *sPattern == '\0' && *sString == '\0' ) return(TRUE);
        if ( *sString == '\0' || *sPattern == '\0' ) return(FALSE);
        if ( *sString == *sPattern ) {
            sPattern++;
            sString++;
            continue;
        }
        if ( *sPattern == '?' ) {
            sPattern++;
            sString++;
            continue;
        }
        if ( *sPattern == '*' ) {
            while ( 1 ) {
                if ( bStringMatchPattern( sString, sPattern+1 ) ) {
		    return(TRUE);
		} else {
                    sString++;
                    if ( *sString == '\0' ) {
                        if ( *(sPattern+1)=='\0' ) {
			    return(TRUE);
			} else {
			    return(FALSE);
			}
                    }
                }
            }
        }
	break;
    }
    return(FALSE);
}
