#ifndef ELEMENTS_H
# define ELEMENTS_H

/*
 *      File:   elements.h
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
 *      Description
 *              Contains all atomic numbers and element names.
 */

#define MAXELEMENTS     	104

/*	Davids Changes */
#define ELEMENTNAMELEN		3
/*	End of Davids Changes */

#define NOELEMENT		-1

#define	NOELEMENTNAME	"?"

#define LONEPAIR        0
#define HYDROGEN        1
#define HELIUM          2
#define LITHIUM         3
#define BERILIUM        4
#define BORON           5
#define CARBON          6
#define NITROGEN        7
#define OXYGEN          8
#define FLOURINE        9
#define NEON            10
#define SODIUM          11
#define MAGNESIUM       12
#define ALUMINUM        13
#define SILICON         14
#define PHOSPHORUS      15
#define SULFUR          16
#define CHLORINE        17
#define ARGON           18
#define	POTASSIUM	19
#define	CALCIUM		20
#define SCANDIUM	21
#define TITANIUM	22
#define VANADIUM	23
#define CHROMIUM	24
#define MANGANESE	25
#define	IRON		26
#define	COBALT		27
#define	NICKLE		28
#define	COPPER		29
#define	ZINC		30
#define	GALLIUM		31
#define	GERMANIUM	32
#define	ARSENIC		33
#define	SELENIUM	34
#define	BROMINE		35
#define	KRYPTON		36

#define RUBIDIUM	37
#define STRONTIUM	38
#define YTTRIUM		39
#define ZIRCONIUM	40
#define NIOBIUM		41
#define MOLYBDENUM	42
#define TECHNETIUM	43
#define RUTHENIUM	44
#define RHODIUM		45
#define PALLADIUM	46
#define SILVER		47
#define CADMIUM		48
#define INDIUM		49
#define TIN		50
#define ANTIMONY	51
#define TELLURIUM	52
#define	IODINE		53
#define	XENON		54

#define CESIUM		55
#define BARIUM		56
#define LANTHANUM	57
#define GADOLINIUM	64
#define HAFNIUM		72
#define TANTALUM	73
#define TUNGSTEN	74
#define RHENIUM		75
#define OSMIUM		76
#define IRIDIUM		77
#define PLATINUM	78
#define GOLD		79
#define MERCURY		80
#define THALLIUM	81
#define LEAD		82
#define BISMUTH		83
#define POLONIUM	84
#define ASTATINE	85
#define RADON		86




typedef struct  {
	char    sName[4];
	int     iNumber;
} ELEMENTt;


extern  ELEMENTt        GeaElements[];



extern int	iElementNumber(char *sName);
extern char	*sElementName(int iNumber, char *sName);
extern int	iElementNumberFromAmber(char *sName);

#define	bElementLegalNumber(n)	((n>=0) && (n<=MAXELEMENTS))

#endif





