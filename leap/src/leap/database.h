/*
 *      File:   database.h
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
 *              Routines to manage an object oriented database.
 *              The database is an ASCII file with
 *              very little internal structure, making
 *              it human readible, editable, and portable.
 */

#ifndef DATABASE_H
#define DATABASE_H

#ifndef	DICTIONARY_H
#include	"dictionary.h"
#endif


#define MAXDATALINELEN  1000
#define MAXPREFIXSTACK  10

typedef struct  {
	long            lFileOffset;
	STRING          sName;
	int             iType;
	int             iRows;
} ENTRYt;

typedef ENTRYt	*ENTRY;

#define	DB_READ		1
#define	DB_WRITE	2

#define	DB_RANDOM_ACCESS	1
#define	DB_SEQUENTIAL_ACCESS	2

typedef struct  {
	int		iAccessMode;
	FILE		*fDataBase;
	STRING          sFileName;
	int		iOpenMode;
	int             iPrefix;
	STRING          saPrefixStack[MAXPREFIXSTACK];
	BOOL            bCompactFileAtClose;
	DICTIONARY      dEntries;
	int             iCurrentLine;
	char            sLookAhead[MAXDATALINELEN];
	int		iLastSequentialOperation;

			/* Loop over entries with a certian prefix */
	STRING		sLoopPrefix;
	DICTLOOP	dlEntryLoop;

} DATABASEt;

typedef DATABASEt	*DATABASE;

#define	LENGTH_NOT_KNOWN	-1

#define ENTRYTYPE       0x0000000F
#define ENTRYINTEGER    0x00000001
#define ENTRYDOUBLE     0x00000002
#define ENTRYSTRING     0x00000003

#define ENTRYMODIFIER   0x000000F0
#define ENTRYSINGLE     0x00000010
#define ENTRYARRAY      0x00000020
#define ENTRYTABLE      0x00000040


/*
 *	DATABASE open modes.
 *
 *	File can be opened for read-only or read-write
 */

#define	OPENREADONLY	1
#define	OPENREADWRITE	2


/*
 *	DATABASE errors.
 */

extern	int	GiDBLastError;
#define	DB_ERROR_NONE			0
#define	DB_ERROR_INVALID_FILE		1
#define	DB_ERROR_INVALID_DATABASE	2

#define iDBLastError()			(GiDBLastError)

/*
 *--------------------------------------------------------------------
 *
 *	Routines for random access files
 */

extern DATABASE		dbDBRndOpen( char *sFileName, int iOpenMode );
extern BOOL		bDBRndDeleteEntry( DATABASE db, char *sOrgEntry );
extern void		DBRndLoopEntryWithPrefix( DATABASE db, 
				char *sOrgEntry );
extern BOOL		bDBRndNextEntryWithPrefix( DATABASE db, char *sEntry );

/*
 *--------------------------------------------------------------------
 *
 *	Routines for accessing DB_SEQUENTIAL_ACCESS files
 */

extern DATABASE		dbDBSeqOpen( char *sFileName, int iOpenMode );
extern void		DBSeqRewind( DATABASE db );
extern void		DBSeqSkipData( DATABASE db );
extern long		lDBSeqCurPos( DATABASE db );
extern void		DBSeqGoto( DATABASE db, long lPos );

/*
 *--------------------------------------------------------------------
 *
 *	Routines used on both DB_RANDOM_ACCESS and
 *	DB_SEQUENTIAL_ACCESS files.
 */


#define	sDBName(d)	(d->sName)

extern BOOL		bDBGetType( DATABASE db, char *sOrgEntry, 
				int *iPType, int *iPLength );
extern BOOL		bDBGetValue( DATABASE dbData, char *sOrgEntry, 
				int *iPLength, GENP PBuffer, int iBufferInc );
extern void		DBPutValue( DATABASE db, char *sOrgEntry, int iType, 
				int iCount, GENP PData, int iDataInc );
extern BOOL		bDBGetTableType( DATABASE db, char *sOrgEntry, 
				int *iPType, int *iPLength,
				int *iPInt1Column, char *sInt1Name,
				int *iPInt2Column, char *sInt2Name,
				int *iPInt3Column, char *sInt3Name,
				int *iPInt4Column, char *sInt4Name,
				int *iPInt5Column, char *sInt5Name,
				int *iPInt6Column, char *sInt6Name,
				int *iPInt7Column, char *sInt7Name,
				int *iPInt8Column, char *sInt8Name,
				int *iPDouble1Column, char *sDouble1Name,
				int *iPDouble2Column, char *sDouble2Name,
				int *iPDouble3Column, char *sDouble3Name,
				int *iPDouble4Column, char *sDouble4Name,
				int *iPString1Column, char *sString1Name,
				int *iPString2Column, char *sString2Name,
				int *iPString3Column, char *sString3Name,
				int *iPString4Column, char *sString4Name,
				int *iPString5Column, char *sString5Name );

extern BOOL	bDBGetTable( DATABASE db, char *sOrgEntry, int *iPLength,
			int iInt1Column, char *PInt1, int iInt1Skip,
			int iInt2Column, char *PInt2, int iInt2Skip,
			int iInt3Column, char *PInt3, int iInt3Skip,
			int iInt4Column, char *PInt4, int iInt4Skip,
			int iInt5Column, char *PInt5, int iInt5Skip,
			int iInt6Column, char *PInt6, int iInt6Skip,
			int iInt7Column, char *PInt7, int iInt7Skip,
			int iInt8Column, char *PInt8, int iInt8Skip,
			int iDouble1Column, char *PDouble1, int iDouble1Skip,
			int iDouble2Column, char *PDouble2, int iDouble2Skip,
			int iDouble3Column, char *PDouble3, int iDouble3Skip,
			int iDouble4Column, char *PDouble4, int iDouble4Skip,
			int iString1Column, char *PString1, int iString1Skip,
			int iString2Column, char *PString2, int iString2Skip,
			int iString3Column, char *PString3, int iString3Skip,
			int iString4Column, char *PString4, int iString4Skip,
			int iString5Column, char *PString5, int iString5Skip );

extern void	DBPutTable( DATABASE db, char *sOrgEntry, int iLines,
	int iInt1Column, char *sInt1, char *PInt1, int iInt1Skip,
	int iInt2Column, char *sInt2, char *PInt2, int iInt2Skip,
	int iInt3Column, char *sInt3, char *PInt3, int iInt3Skip,
	int iInt4Column, char *sInt4, char *PInt4, int iInt4Skip,
	int iInt5Column, char *sInt5, char *PInt5, int iInt5Skip,
	int iInt6Column, char *sInt6, char *PInt6, int iInt6Skip,
	int iInt7Column, char *sInt7, char *PInt7, int iInt7Skip,
	int iInt8Column, char *sInt8, char *PInt8, int iInt8Skip,
	int iDouble1Column, char *sDouble1, char *PDouble1, int iDouble1Skip,
	int iDouble2Column, char *sDouble2, char *PDouble2, int iDouble2Skip,
	int iDouble3Column, char *sDouble3, char *PDouble3, int iDouble3Skip,
	int iDouble4Column, char *sDouble4, char *PDouble4, int iDouble4Skip,
	int iString1Column, char *sString1, char *PString1, int iString1Skip,
	int iString2Column, char *sString2, char *PString2, int iString2Skip,
	int iString3Column, char *sString3, char *PString3, int iString3Skip,
	int iString4Column, char *sString4, char *PString4, int iString4Skip,
	int iString5Column, char *sString5, char *PString5, int iString5Skip );

extern void	DBClose( DATABASE *Pdb );
extern void	DBPushPrefix( DATABASE db, char *sStr );
extern void	DBPopPrefix( DATABASE db );
extern void	DBZeroPrefix( DATABASE db );
extern void	DBPushZeroPrefix( DATABASE db, char *sStr );


#endif  /* ifdef DATABASE_H */
