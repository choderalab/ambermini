/*
 *	File:	block.c
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
 *	Description:
 *		A BLOCK is an object that maintains an arbitrarily
 *		large string that contains a single command for the 
 *		parser.  The string is \0 terminated and
 *		can contains several lines seperated by \n.
 *		The purpose of breaking up the command over several
 *		lines is so that when the command is written to the journal
 *		they are not more than 80 characters.
 *
 *		BLOCKs will continue to accept characters until 
 *		it decides that it has recieved an entire command.
 *		The criteria that a BLOCK uses to determine whether
 *		or not a command is been completed are the following.
 *
 *		Commands are terminated by '\n' UNLESS
 *		1)	The last character on the line is a 
 *			continuation character '\'.
 *		2)	A LIST is not yet finished, eg: for every '{' seen
 *			there is not a '}'.
 *
 *		The BLOCK removes continuation characters from the input
 *		and expects the parser to put them back (even for LISTs)
 *		when they are written out to the LogFile.
 *
 */





#include	"basics.h"

#include	"block.h"


		/* By default grow the BLOCK text by 100 characters */
		/* each time */

#define	BLOCKGROW	100


/*
 *-----------------------------------------------------------------
 *
 *	Private routines
 *
 */




/*
 *	cBlockCurChar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the last character added to the BLOCK.
 */
static char	
cBlockLastChar( BLOCK bBlock )
{

    if ( bBlock->iTextNext <= 0 ) 
	return('\0');
    return(bBlock->cPText[bBlock->iTextNext-1]);
}





/*
 *	zBlockAddOneCharacter
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add one character to the BLOCK.  If there is not enough room
 *	in the BLOCK then grow it by BLOCKGROW characters.
 */
static void
zBlockAddOneCharacter( BLOCK bBlock, char c )
{
char	*cPText;
int	iSize;

		/* If there isn't enough space then grow the BLOCK */

    iSize = bBlock->iTextSize + BLOCKGROW;

    if ( (bBlock->iTextSize-1) <= bBlock->iTextNext ) {
	if ( bBlock->cPText == NULL ) {
	    MALLOC( cPText, char*, iSize );
	    MESSAGE(( "MALLOCed the BLOCK text\n" ));
        } else {
	    REALLOC( cPText, char*, bBlock->cPText, iSize );
	    MESSAGE(( "REALLOCed the BLOCK text\n" ));
	}
	bBlock->cPText = cPText;
        bBlock->iTextSize = iSize;
    }

		/* Now put in the character */

    bBlock->cPText[bBlock->iTextNext++] = c;
    bBlock->cPText[bBlock->iTextNext] = '\0';
}




/*
 *======================================================================
 *
 *	Public routines
 */






/*
 *	bBlockCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create an empty BLOCK and return it.
 */
BLOCK	
bBlockCreate()
{
BLOCK	bNew;

    MALLOC( bNew, BLOCK, sizeof(BLOCKt) );

    bNew->cPText = NULL;
    bNew->iTextSize = 0;
    bNew->iTextNext = 0;
    bNew->iReadPos = 0;
    bNew->iUnclosedLists = 0;
    return(bNew);
}







/*
 *	bBlockAddChar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add characters to the command buffer.
 *	If the character to add ends the block then return TRUE.
 */
BOOL
bBlockAddChar( BLOCK bBlock, char c )
{
BOOL	bEndOfBlock;


    bEndOfBlock = FALSE;
    if ( c == '\\' ) {
        MESSAGE(( "Got continuation\n" ));
    } else if ( c == '{' ) {
	bBlock->iUnclosedLists++;
    } else if ( c == '}' ) {
	bBlock->iUnclosedLists--;
 	if ( bBlock->iUnclosedLists < 0 ) bBlock->iUnclosedLists = 0;
    } else if ( c == '\n' ) {
	if ( cBlockLastChar(bBlock) != '\\' )
	    if ( bBlock->iUnclosedLists == 0 ) bEndOfBlock = TRUE;
    }

    zBlockAddOneCharacter( bBlock, c );

#ifdef	DEBUG
if ( bEndOfBlock ) {
    MESSAGE(( "Ending block.  Text=%s\n", sBlockText(bBlock) ));
}
#endif

    return(bEndOfBlock);
}




/*
 *	bBlockRemoveChar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Remove the last character from the BLOCK.
 *	If there are no characters to remove then return FALSE.
 */
BOOL
bBlockRemoveChar( BLOCK bBlock )
{
    if ( bBlock->iTextNext == 0 ) 
	return(FALSE);
    if ( cBlockLastChar(bBlock) == '{' ) 
	bBlock->iUnclosedLists--;
    else if ( cBlockLastChar(bBlock) == '}' ) 
	bBlock->iUnclosedLists++;
    bBlock->iTextNext--;
    return(TRUE);
}




/*
 *	BlockEmpty
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Empty the BLOCK.
 */
void	
BlockEmpty( BLOCK bBlock )
{
    bBlock->iTextNext = 0;
    bBlock->iReadPos = 0;
}





/*
 *	BlockDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the BLOCK.
 */
void	
BlockDestroy( BLOCK *bPBlock )
{
	/* If there is memory allocated for the BLOCK then free it */

    if ( (*bPBlock)->cPText != NULL ) FREE( (*bPBlock)->cPText );

	/* Now free the BLOCK itself */

    FREE( *bPBlock );
    *bPBlock = NULL;
}

  
    



/*
 *	bBlockReadLine
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Modified 10-Nov-92, David A. Rivkin
 *		Added test for missing close quotes in block.
 *			(Bug Fix)
 *
 *	Read a line from the block into the line buffer.
 *	Return TRUE if the next character in the BLOCK is a '\0'
 *	meaning that the line is the last line in the block.
 */
BOOL
bBlockReadLine( BLOCK bBlock, char *sLine )
{
int		i;
BOOL		bQuote = FALSE;

    for ( i=0; cBlockPeek( bBlock ) != '\n'; i++) {
	if (((sLine[i] = cBlockRead( bBlock ))) == '"' ) {
	    if ( bQuote == TRUE )
	    	bQuote = FALSE;
	    else 
	    	bQuote = TRUE;
	}
    }
    sLine[i] = cBlockRead(bBlock);
    if ( bQuote == TRUE ) {
	VP0(( "A close quote was missing.\n" ));
	VP0(( "I have added quotes to the end of the line and will try it.\n" ));
	sLine[i++] = '"';
	sLine[i] = '\n';
	sLine[i+1] = '\0';
	VP0((sLine));
    }
    sLine[++i] = '\0';
    return(( cBlockPeek(bBlock) == '\0' ));
}
