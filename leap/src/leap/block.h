/*
 *	File:	block.h
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



#ifndef	BLOCK_H
# define BLOCK_H

typedef	struct	{
	char		*cPText;
	int		iTextSize;
	int		iTextNext;
	int		iUnclosedLists;
	int		iReadPos;
} BLOCKt;

typedef	BLOCKt	*BLOCK;




#define	sBlockText(b)		( (b)->cPText )

#define	BlockResetRead(b)	( (b)->iReadPos = 0 )
#define	cBlockRead(b)		( (b)->cPText[(b)->iReadPos++] )
#define	cBlockPeek(b)		( (b)->cPText[(b)->iReadPos] )

#define	bBlockEndOfRead(b)	( (b)->iReadPos >= (b)->iTextNext )
#define	cBlockLastWrittenChar(b)	( \
	(b)->iTextNext>0 ? (b)->cPText[(b)->iTextNext-1] : '\0' )


/*  block.c  */

extern BLOCK		bBlockCreate();
extern BOOL		bBlockRemoveChar( BLOCK bBlock );
extern void		BlockEmpty( BLOCK bBlock );
extern void		BlockDestroy( BLOCK *bPBlock );
extern BOOL		bBlockReadLine( BLOCK bBlock, char *sLine );
extern BOOL             bBlockAddChar( BLOCK bBlock, char c );

#endif
