/*
 *      File:   elements.c
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
 *              Manages the atomic numbers and element names.
 */


#include	"basics.h"

#include        "elements.h"

/*	1st 7 rows of periodic table	*/

ELEMENTt        GeaElements[] = {

	{       "Lp",   LONEPAIR        },

	{       "H",    HYDROGEN        },
	{       "He",   HELIUM          },
	
	{       "Li",   LITHIUM         },
	{       "Be",   BERYLLIUM       },
	{       "B",    BORON           },
	{       "C",    CARBON          },
	{       "N",    NITROGEN        },
	{       "O",    OXYGEN          },
	{       "F",    FLUORINE        },
	{       "Ne",   NEON            },
	
	{       "Na",   SODIUM          },
	{       "Mg",   MAGNESIUM       },
	{       "Al",   ALUMINUM        },
	{       "Si",   SILICON         },
	{       "P",    PHOSPHORUS      },
	{       "S",    SULFUR          },
	{       "Cl",   CHLORINE        },
	{       "Ar",   ARGON           },
	
	{	"K",	POTASSIUM	},
	{	"Ca",	CALCIUM		},
	{	"Sc",	SCANDIUM	},
	{	"Ti",	TITANIUM	},
	{	"V",	VANADIUM	},
	{	"Cr",	CHROMIUM	},
	{	"Mn",	MANGANESE	},
	{	"Fe",	IRON		},
	{	"Co",	COBALT		},
	{	"Ni",	NICKLE		},
	{	"Cu",	COPPER		},
	{	"Zn",	ZINC		},
	{	"Ga",	GALLIUM		},
	{	"Ge",	GERMANIUM	},
	{	"As",	ARSENIC		},
	{	"Se",	SELENIUM	},
	{	"Br",	BROMINE		},
	{	"Kr",	KRYPTON		},
	
	{	"Rb",	RUBIDIUM	},
	{	"Sr",	STRONTIUM	},
	{	"Y",	YTTRIUM		},
	{	"Zr",	ZIRCONIUM	},
	{	"Nb",	NIOBIUM		},
	{	"Mo",	MOLYBDENUM	},
	{	"Tc",	TECHNETIUM	},
	{	"Ru",	RUTHENIUM	},
	{	"Rh",	RHODIUM		},
	{	"Pd",	PALLADIUM	},
	{	"Ag",	SILVER		},
	{	"Cd",	CADMIUM		},
	{	"In",	INDIUM		},
	{	"Sn",	TIN		},
	{	"Sb",	ANTIMONY	},
	{	"Te",	TELLURIUM	},
	{	"I",	IODINE		},
	{	"Xe",	XENON		},
	
	{	"Cs",	CESIUM		},
	{	"Ba",	BARIUM		},
	{	"La",	LANTHANUM	},
	{	"Ce",	CERIUM		},
	{	"Pr",	PRASEODYMIUM	},
	{	"Nd",	NEODYMIUM	},
	{	"Pm",	PROMETHIUM	},
	{	"Sm",	SAMARIUM	},
	{	"Eu",	EUROPIUM	},
	{	"Gd",	GADOLINIUM	},
	{	"Tb",	TERBIUM		},
	{	"Dy",	DYSPROSIUM	},
	{	"Ho",	HOLMIUM		},
	{	"Er",	ERBIUM		},
	{	"Tm",	THULIUM		},
	{	"Yb",	YTTERBIUM	},
	{	"Lu",	LUTETIUM	},
	{	"Hf",	HAFNIUM		},
	{	"Ta",	TANTALUM	},
	{	"W",	TUNGSTEN	},
	{	"Re",	RHENIUM		},
	{	"Os",	OSMIUM		},
	{	"Ir",	IRIDIUM		},
	{	"Pt",	PLATINUM	},
	{	"Au",	GOLD		},
	{	"Hg",	MERCURY		},
	{	"Tl",	THALLIUM	},
	{	"Pb",	LEAD		},
	{	"Bi",	BISMUTH		},
	{	"Po",	POLONIUM	},
	{	"At",	ASTATINE	},
	{	"Rn",	RADON		},

	{	"Fr",	FRANCIUM	},
	{	"Ra",	RADIUM		},
	{	"Ac",	ACTINIUM	},
	{	"Th",	THORIUM		},
	{	"Pa",	PROTACTINIUM	},
	{	"U",	URANIUM		},
	{	"Np",	NEPTUNIUM	},
	{	"Pu",	PLUTONIUM	},
	{	"Am",	AMERICIUM	},
	{	"Cm",	CURIUM		},
	{	"Bk",	BERKELIUM	},
	{	"Cf",	CALIFORNIUM	},
	{	"Es",	EINSTEINIUM	},
	{	"Fm",	FERMIUM		},
	{	"Md",	MENDELEVIUM	},
	{	"No",	NOBELIUM	},
	{	"Lr",	LAWRENCIUM	},
	{	"Rf",	RUTHERFORDIUM	},
	{	"Db",	DUBNIUM		},
	{	"Sg",	SEABORGIUM	},
	{	"Bh",	BOHRIUM		},
	{	"Hs",	HASSIUM		},
	{	"Mt",	MEITNERIUM	},
	{	"Ds",	DARMSTADTIUM	},
	{	"Rg",	ROENTGENIUM	},
	{	"Cn",	COPERNICIUM	},
	{	"Uut",  UNUNTRIUM	},
	{	"Fl",	FLEROVIUM	},
	{	"Uup",	UNUNPENTIUM	},
	{	"Lv",	LIVERMORIUM	},
	{	"Uus",	UNUNSEPTIUM	},
	{	"Uuo",	UNUNOCTIUM	},

	/* Last one must be empty */
	{       "",     NOELEMENT       }
};



/*
 *      iElementNumber
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the element number for the element name.
 */
int     
iElementNumber( char *sName )
{
int             i;

    for (i=0; GeaElements[i].iNumber != NOELEMENT; i++)
        if ( strcmp( GeaElements[i].sName, sName ) == 0 ) 
            return(GeaElements[i].iNumber);

    return(NOELEMENT);
}





/*
 *      sElementName
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the element name for the element number.
 */
char *
sElementName( int iNumber, char *sName )
{
int             i;


    for (i=0; GeaElements[i].iNumber != NOELEMENT; i++) {
        if ( GeaElements[i].iNumber == iNumber ) {
            strcpy( sName, GeaElements[i].sName );
            return(sName);
        }
    }
    strcpy( sName, NOELEMENTNAME );
    return(sName);
}




/*
 *      iElementNumberFromAmber
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the element number for the element name.
 *	Compare to uppercase names, and lower case names
 *	which means
 *	first compare the two character element names,
 *	then the one character element names.
 *	Only the first characters of sName are compared.
 */
int
iElementNumberFromAmber( char *sName )
{
int             i;
STRING		sTemp;

    /*
     *  strip off any spaces on general principles and
     *	numbers in case sName is from a pdb file
     */

    while ( *sName == ' ' ) sName++;
    while ( isdigit( *sName ) ) sName++;

		/* Compare the two-character element names */
    for (i=0; GeaElements[i].iNumber != NOELEMENT; i++) {

	if (strlen(GeaElements[i].sName) == 2 ) {
	    strcpy( sTemp, GeaElements[i].sName );
	    sTemp[1] = cUpper(sTemp[1]);
            if ( strncmp( sTemp, sName, 2 ) == 0 ||
	         strncmp( GeaElements[i].sName, sName, 2 ) == 0 ) 
	        return(GeaElements[i].iNumber);
	}
    }
		/* Compare the one-character element names */
    for (i=0; GeaElements[i].iNumber != NOELEMENT; i++) {

	if (strlen(GeaElements[i].sName) == 1 ) {
	    if ( GeaElements[i].sName[0] == sName[0] )
		return(GeaElements[i].iNumber);
	}
    }

    return(NOELEMENT);
}




