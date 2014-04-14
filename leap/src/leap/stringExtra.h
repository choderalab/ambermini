/*
 *      File:   stringExtra.h
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
 *              Some extra routines to manipulate strings.
 */
 
#ifndef	STRINGEXTRA_H
#define	STRINGEXTRA_H



extern char	*sRemoveSpaces(char *sIn, char *sOut);
extern char	*sRemoveControlAndPadding(char *sRaw, char *sResult);
extern char	*sRemoveLeadingSpaces(char *sLine);
extern char	*sRemoveFirstString(char *sLine, char *sHead);
extern BOOL	bStringMatchPattern(char *sString, char *sPattern);
extern void	StringCopyMax(char *sDest, char *sSource, int iMax);

#endif
