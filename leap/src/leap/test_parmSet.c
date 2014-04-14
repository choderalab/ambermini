/*
 *	File:	test_parmSet.c
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
 *		Test the parmSet object.
 */


#include	"basics.h"

#include	"classes.h"

#include	"database.h"
#include	"amber.h"



_FUNC	void
main( argc, argv )
int		argc;
char*		argv[];

_VAR
PARMSET	psNew;
DATABASE	dbTemp;
STRING		s1, s2, s3, s4;
TORSION		tTorsion;
int		i, iN;
double		dKp, dP0;
PARMSET		psCur;
int		iIndex;

_BEGIN

    BasicsInitialize();

    if ( argc != 2 ) {
    	PRINTF(( "usage: %s {OFF_parameter_file}\n", argv[0] ));
	exit(1);
    }
    PRINTF(( "Search parameter set: %s\n", argv[1] ));

    dbTemp = dbDBRndOpen( argv[1], OPENREADWRITE );
    DBPushPrefix( dbTemp, "entry.parameters." );
    psCur = psParmSetLoad( dbTemp );
    if ( psCur == NULL ) {
    	PRINTF(( "Could not read the PARMSET\n" ));
	exit(1);
    }
    
    while ( 1==1 ) {
        PRINTF(( "Enter torsion types: " ));
	scanf( "%s %s %s %s", s1, s2, s3, s4 );
	PRINTF(( "Looking for: %s-%s-%s-%s\n", s1, s2, s3, s4 ));
	tTorsion = tParmSetTORSIONCreate();
        iParmSetFindImproperTerms( psCur, tTorsion, TRUE,
				    s1, s2, s3, s4 );
	PRINTF(( "Found %d terms\n", iParmSetTORSIONTermCount(tTorsion) ));
	for ( i=0; i<iParmSetTORSIONTermCount(tTorsion); i++ ) {
	    ParmSetTORSIONTerm( tTorsion, i, &iIndex,
	    				s1, s2, s3, s4,
					&iN, &dKp, &dP0 );
	    PRINTF(( "%d: %s-%s-%s-%s   %lf  %d  %lf\n",
	    		iIndex, s1, s2, s3, s4, dKp, iN, dP0 ));
	}
	PRINTF(( "\n" ));
	ParmSetTORSIONDestroy(&tTorsion);
    }
_END
    
    


