/*
 *	File:	test.c
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
 *	Author:	Christian Schafmeister (1991)
 *
 *	Description:
 *		Test the Serial Object File Format. (SOFFt).
 */



#include	"basics.h"

#include	"database.h"


FUNC
main( argc, argv )
int		argc;
char*		argv[];
{
STRING		sFile;
DATABASE	db;
STRING		sName;
int		iType, iLength, i;
int		iInt1Column;
int		iInt2Column;
int		iInt3Column;
int		iInt4Column;
int		iInt5Column;
int		iInt6Column;
int		iInt7Column;
int		iInt8Column;
int		iDouble1Column;
int		iDouble2Column;
int		iDouble3Column;
int		iDouble4Column;
int		iString1Column;
int		iString2Column;
int		iString3Column;
int		iString4Column;
int		iString5Column;
STRING		sInt1Name;
STRING		sInt2Name;
STRING		sInt3Name;
STRING		sInt4Name;
STRING		sInt5Name;
STRING		sInt6Name;
STRING		sInt7Name;
STRING		sInt8Name;
STRING		sDouble1Name;
STRING		sDouble2Name;
STRING		sDouble3Name;
STRING		sDouble4Name;
STRING		sString1Name;
STRING		sString2Name;
STRING		sString3Name;
STRING		sString4Name;
STRING		sString5Name;
int		iaData[100];
double		daData[100];
STRING		saData[100];

int		iaData1[100];
int		iaData2[100];
int		iaData3[100];
int		iaData4[100];
int		iaData5[100];
int		iaData6[100];
int		iaData7[100];
int		iaData8[100];

double		daData1[100];
double		daData2[100];
double		daData3[100];
double		daData4[100];

STRING		saData1[100];
STRING		saData2[100];
STRING		saData3[100];
STRING		saData4[100];
STRING		saData5[100];

int		iLines;

int		iVal;
double		dVal;
STRING		sVal;

BEGIN


    BasicsInitialize();

    if ( argc != 2 ) {
	printf( "Usage %s {OFF file}\n", argv[1] );
	exit(1);
    }

		/* Open the DATABASE */

    strcpy( sFile, argv[1] );
    db = dbDBSeqOpen( sFile, OPENREADWRITE );

    while ( bDBGetType( db, sName, &iType, &iLength ) ) {
	printf( "\nRead entry: %s\n", sName );
	switch ( iType & ENTRYMODIFIER ) {
	    case ENTRYSINGLE:
		printf( "++Single value\n" );
		break;
	    case ENTRYARRAY:
		printf( "++Array\n" );
		break;
	    case ENTRYTABLE:
		printf( "++Table\n" );
		break;
	    default:
		printf( "++    Illegal modifier type\n" );
		break;
	}

	if ( (iType&ENTRYMODIFIER) == ENTRYTABLE ) {

	    bDBGetTableType( db, sName, &iType, &iLength,
                        	&iInt1Column, sInt1Name,
       	               	 	&iInt2Column, sInt2Name,
                        	&iInt3Column, sInt3Name,
                        	&iInt4Column, sInt4Name,
                        	&iInt5Column, sInt5Name,
                        	&iInt6Column, sInt6Name,
                        	&iInt7Column, sInt7Name,
                        	&iInt8Column, sInt8Name,
                        	&iDouble1Column, sDouble1Name,
                        	&iDouble2Column, sDouble2Name,
                        	&iDouble3Column, sDouble3Name,
                        	&iDouble4Column, sDouble4Name,
                        	&iString1Column, sString1Name,
                        	&iString2Column, sString2Name,
                        	&iString3Column, sString3Name,
                        	&iString4Column, sString4Name,
                        	&iString5Column, sString5Name );

	    bDBGetTable( db, sName, &iLength,
				iInt1Column, iaData1, sizeof(int),
				iInt2Column, iaData2, sizeof(int),
				iInt3Column, iaData3, sizeof(int),
				iInt4Column, iaData4, sizeof(int),
				iInt5Column, iaData5, sizeof(int),
				iInt6Column, iaData6, sizeof(int),
				iInt7Column, iaData7, sizeof(int),
				iInt8Column, iaData8, sizeof(int),
				iDouble1Column, daData1, sizeof(double),
				iDouble2Column, daData2, sizeof(double),
				iDouble3Column, daData3, sizeof(double),
				iDouble4Column, daData4, sizeof(double),
				iString1Column, saData1, sizeof(STRING),
				iString2Column, saData2, sizeof(STRING),
				iString3Column, saData3, sizeof(STRING),
				iString4Column, saData4, sizeof(STRING),
				iString5Column, saData5, sizeof(STRING) );

#define	PPI(a)	for( i=0; i<iLength; i++ ) \
			printf( "Int[%2d] = %4d\n", i, a[i] );

	    if ( iInt1Column ) {
		printf( "Int col=%2d  name=%s\n", iInt1Column, sInt1Name );
		PPI(iaData1);
	    }
	    if ( iInt2Column ) {
		printf( "Int col=%2d  name=%s\n", iInt2Column, sInt2Name );
		PPI(iaData2);
	    }
	    if ( iInt3Column ) {
		printf( "Int col=%2d  name=%s\n", iInt3Column, sInt3Name );
		PPI(iaData3);
	    }
	    if ( iInt4Column ) {
		printf( "Int col=%2d  name=%s\n", iInt4Column, sInt4Name );
		PPI(iaData4);
	    }
	    if ( iInt5Column ) {
		printf( "Int col=%2d  name=%s\n", iInt5Column, sInt5Name );
		PPI(iaData5);
	    }
	    if ( iInt6Column ) {
		printf( "Int col=%2d  name=%s\n", iInt6Column, sInt6Name );
		PPI(iaData6);
	    }
	    if ( iInt7Column ) {
		printf( "Int col=%2d  name=%s\n", iInt7Column, sInt7Name );
		PPI(iaData7);
	    }
	    if ( iInt8Column ) {
		printf( "Int col=%2d  name=%s\n", iInt8Column, sInt8Name );
		PPI(iaData8);
	    }

#define	PD(a)	for( i=0; i<iLength; i++ ) \
			printf( "Double[%2d] = %lf\n", i, a[i] );
	    if ( iDouble1Column ) {
		printf( "Dbl col=%2d name=%s\n", iDouble1Column, sDouble1Name );
		PD(daData1);
	    }
	    if ( iDouble2Column ) {
		printf( "Dbl col=%2d name=%s\n", iDouble2Column, sDouble2Name );
		PD(daData2);
	    }
	    if ( iDouble3Column ) {
		printf( "Dbl col=%2d name=%s\n", iDouble3Column, sDouble3Name );
		PD(daData3);
	    }
	    if ( iDouble4Column ) {
		printf( "Dbl col=%2d name=%s\n", iDouble4Column, sDouble4Name );
		PD(daData4);
	    }

#define	PS(a)	for( i=0; i<iLength; i++ ) \
			printf( "String[%2d] = %s\n", i, a[i] );
	    if ( iString1Column ) {
		printf( "Str col=%2d name=%s\n", iString1Column, sString1Name );
		PS(saData1);
	    }
	    if ( iString2Column ) {
		printf( "Str col=%2d name=%s\n", iString2Column, sString2Name );
		PS(saData2);
	    }
	    if ( iString3Column ) {
		printf( "Str col=%2d name=%s\n", iString3Column, sString3Name );
		PS(saData3);
	    }
	    if ( iString4Column ) {
		printf( "Str col=%2d name=%s\n", iString4Column, sString4Name );
		PS(saData4);
	    }
	    if ( iString5Column ) {
		printf( "Str col=%2d name=%s\n", iString5Column, sString5Name );
		PS(saData5);
	    }

	} else {
	    switch ( iType & ENTRYTYPE ) {
	    	case ENTRYINTEGER:
		    printf( "--Integer\n" );
		    bDBGetValue( db, sName, &iLines, iaData, sizeof(int) );
		    printf( "Entry %s has %d lines\n", sName, iLines );
		    for ( i=0; i<iLines; i++ ) {
			printf( "Value[%2d] = %5d\n", i, iaData[i] );
		    }
		    break;
	    	case ENTRYDOUBLE:
		    printf( "--Double\n" );
		    bDBGetValue( db, sName, &iLines, daData, sizeof(double) );
		    printf( "Entry %s has %d lines\n", sName, iLines );
		    for ( i=0; i<iLines; i++ ) {
			printf( "Value[%2d] = %lf\n", i, daData[i] );
		    }
		    break;
	    	case ENTRYSTRING:
		    printf( "--String\n" );
		    bDBGetValue( db, sName, &iLines, saData, sizeof(STRING) );
		    printf( "Entry %s has %d lines\n", sName, iLines );
		    for ( i=0; i<iLines; i++ ) {
			printf( "Value[%2d] = %s\n", i, saData[i] );
		    }
		    break;
	    	default:
		    printf( "--    illegal type\n" );
		    break;
	    }
	}
    }

    iVal = 1234;
    DBPutValue( db, "new.single.integer", ENTRYSINGLE|ENTRYINTEGER,
			1, &iVal, sizeof(int) );
    DBPutValue( db, "new.array.integer", ENTRYARRAY|ENTRYINTEGER,
			10, &iVal, 0 );
    dVal = 1234.5678;
    DBPutValue( db, "new.single.double", ENTRYSINGLE|ENTRYDOUBLE,
			1, &dVal, sizeof(double) );
    DBPutValue( db, "new.array.double", ENTRYARRAY|ENTRYDOUBLE,
			10, &dVal, 0 );

    strcpy( sVal, "Hello there" );
    DBPutValue( db, "new.single.string", ENTRYSINGLE|ENTRYSTRING,
			1, sVal, sizeof(STRING) );
    DBPutValue( db, "new.array.string", ENTRYARRAY|ENTRYSTRING,
			10, sVal, 0 );


    DBPutTable( db, "test.sequential.table", 10,
		1, "an_int", &iVal, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,

		2, "a_double", &dVal, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,

		3, "a_string", sVal, 0,
		4, "a_string", sVal, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0,
		0, NULL, NULL, 0 );
	
    DBClose( &db );

END
}


