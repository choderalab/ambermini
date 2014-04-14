/*
 *      File:   amber.h
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
 *              Read AMBER prep and parameter database files into
 *              UNITs and PARMSETs.
 *
 *              This module was written to maintain compatability with
 *              AMBER 3a.
 *
 */



#ifndef	AMBER_H
#define	AMBER_H

#ifndef	DICTIONARY_H
#include	"dictionary.h"
#endif

#define	AMBER_WILD_CARD	"X"


/*  amber.c  */

extern void		AmberAddAtomTypes( LIST lEntries );
extern DICTIONARY	dAmberReadPrepFile( char *sFilename );
extern PARMSET		psAmberReadParmSet( FILE *fIn, char *sFilename );
extern int		iAtomTypeHybridization( char *sType );

#endif
