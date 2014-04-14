/*
 *      File:   library.h
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



#ifndef	LIBRARY_H
#define	LIBRARY_H

#ifndef	DICTIONARY_H
# include "dictionary.h"
#endif

#ifndef	DATABASE_H
# include "database.h"
#endif


typedef struct  {
	DICTIONARY	dLibrary;
	DATABASE	dbLibrary;
	DICTLOOP	dlContents;
} LIBRARYt;

typedef LIBRARYt	*LIBRARY;




extern LIBRARY	lLibraryOpen(char *sName, int iOpenMode);
extern OBJEKT	oLibraryLoad(LIBRARY ul, char *sName);
extern void	LibrarySave(LIBRARY ul, char *sName, OBJEKT oObj, PARMLIB plLibrary);
extern BOOL	bLibraryRemove(LIBRARY ul, char *sName);

extern void	LibraryLoop(LIBRARY ul);
extern char	*sLibraryNext(LIBRARY ul);
extern void	LibraryClose(LIBRARY *lPLibrary);

extern BOOL	bLibraryDBPrefix(LIBRARY lLib, char *cPResult, char *sStr );
#define	dbLibraryDatabase(l)	(l->dbLibrary)


#endif

