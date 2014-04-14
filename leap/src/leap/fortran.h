/*
 *      File:   fortran.h
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
 *
 *              This file contains routines for performing
 *              FORTRAN style formatted output.
 *              It contains a routine that allows the caller
 *              to specify the number of entries that should
 *              appear on each line and the format of each
 *              entry on the line.
 *              Then there is a routine that the user uses
 *              to write each subsequent piece of data.
 */

#ifndef	FORTRAN_H
# define FORTRAN_H

extern void	FortranFile(FILE *fOut);
extern void	FortranFormat(int iPerLine, char *sFormat);

extern void	FortranWriteInt(int iVal);
extern void	FortranWriteDouble(double dVal);
extern void	FortranWriteString(char *sVal);

extern void	FortranEndLine();

extern char	*sFortranReadString( char *sString );
extern char	*sFortranReadLabel( char *sString );

extern int	iFortranReadInt();
extern double	dFortranReadDouble();
extern void	FortranSkipLine();

extern void	FortranDebugOn();
extern void	FortranDebug(char *sStr);	

extern STRING	GsFortranDebug;


#endif
