/* 
 *	File:  		hybrid.c
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
 *	Author:		David A. Rivkin
 *	Date Created:	17 August 1992
 *	Dates Changed:
 *	Description:
 * 	Simple functions for translating Hybridization information.
 * 	Converts the interger representation to a string and vice versa.
 *	Due to the limited number of possible hybridizations, simple switch/case, if/then
 *	statements do the conversion.  This will be the fastest running as well.
 *
 */

#include	"basics.h"

	/* Hybridization type */

#define	HUNDEFINED	-1
#define	HUNKNOWN	0
#define	HSP1		1
#define	HSP2		2
#define	HSP3		3

/*
 *      iHybridNumber
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the element number for the element name.
 */
int
iHybridNumber( char *sName )
{
    StringLower(sName);

    if(!strcmp(sName, "undefined"))	return(-1);
    if(!strcmp(sName, "unknown"))	return(0);
    if(!strcmp(sName, "sp1"))		return(1);
    if(!strcmp(sName, "sp2"))		return(2);
    if(!strcmp(sName, "sp3"))		return(3);
    return(-1);	/* undefined */
}


/*
 *      sHybridName
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the element name for the element number.
 */
char *
sHybridName( int iHybrid, char *sName )
{
    switch(iHybrid){
	case -1: strcpy( sName, "undefined" );
		break;
	case 1:	strcpy( sName, "sp1" );
		break;
	case 2:	strcpy( sName, "sp2" );
		break;
	case 3:	strcpy( sName, "sp3" );
		break;
	case 0:	
	default: strcpy( sName, "unknown" );
		break;
    }
    return( sName );		
}

