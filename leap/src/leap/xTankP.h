#ifndef _TankP_h
#define _TankP_h

/*
 *	File:	xTankP.h
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
 */

/* include superclass private header file */
#include <X11/CoreP.h>

#include       "elements.h"
# include       "varArray.h"
# include        "threed.h"

# include	"x3d.h"




/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 *
 *	Define the CLASS structure.
 */


typedef struct	{
	int		empty;
} TankClassPart;



typedef struct _TankClassRec {
	CoreClassPart	core_class;
	TankClassPart	tank_class;
} TankClassRec;

extern TankClassRec tankClassRec;


/*
 ********************************************************************
 *
 *	Define TANK Widget structure
 *
 */

		/* Button states */
#define	BNONE		0
#define	BPOINT		1
#define	BROTATE		2
#define	BTRANSXY	3
#define	BSCALE		4
#define	BSHIFT_POINT	5


		/* Selections can be specified different ways */
		/* depending on how the user starts the selection */

#define	TANK_SELECT_ATOM_OR_CHAIN	0
#define	TANK_SELECT_BOX			1




typedef	struct	{
	ATOM	aAtomStart;
	ATOM	aAtomInvisible;
	BOOL	bInRing;
} TWISTTORSIONt;



          
typedef struct {
/* Resources */
	String		cPColor0;
	String		cPColor1;
	String		cPColor2;
	String		cPColor3;
	String		cPColor4;
	String		cPColor5;
	String		cPColor6;
	String		cPColor7;
	String		cPMonoForeground;
	String		cPMonoBackground;
	int		iSelectDelay;
	String		cPRenderMethod;
/* Fields that are independant of display technique */
	BOOL            bGotInitialExpose;
	BOOL		bNeedToRedraw;
	X3DENGINE	x3dEngine;
	FLAGS		fFlags;		/* Control what to show */
	int		iPrintSink;
	Widget          wWidget;
	UNIT            uUnit;
	VARARRAY	vaAtomPtrs;
/* Maintain the state of the pointer TANKDRAW, TANKERASE, TANKSELECT */
	int		iTankState;
/* iButtonState is used to keep track of the current button mode */
	int		iButtonState;
	int		iRawButtonState;		/* Mask of buttons pressed */
	GVECTOR		gvTransXY;
	int		iActionX;
	int		iActionY;
/* The following fields are used for drawing */
	BOOL		bRubberOn;
	int		iSelectionType;
	int		iXStart;
	int		iYStart;
	int		iXStop;
	int		iYStop;
	int		iCurrentDrawingElement;
	int		iNextAtomNumber;
/* The following fields are for changing what happens when the pointer */
/* is clicked */
	int		iTankSelectLastState;
	Time		tLastClickMSec;
	Time		tTimeSinceLastClick;
/* The following fields are used to store data for different TANK states */
	VARARRAY	vaTwistTorsions;
	VARARRAY	vaRotateAtoms;
	ATOM		aRotateCenter;
/* Store special information for the different rendering techniques */
	Pixmap		pOffScreen;
	BOOL		bPageFlipFront;
	BOOL		bFastPageFlipLastTime;
} TankPart;


typedef struct _TankRec {
	CorePart	core;
	TankPart	tank;
} TankRec;

typedef	TankRec	*TANK;


typedef	struct	{
	GVECTOR		gvScreen;
} GRAPHATOMt;




		/* Some default values */

#define	TBLACK		0
#define	TRED		1
#define	TGREEN		2
#define	TBLUE		3
#define	TPURPLE		4
#define	TCYAN		5
#define	TYELLOW		6
#define	TWHITE		7


#define	DEFAULTSELECTTHICKNESS	4		/* Make monochrome selected lines 4 pixels thick */




                /* Mouse must be within 2 pixels of the atom */
                /* in order to pick it */
#define PICKDIST		8

		/* Distance squared in pixels for line/bond picks */
#define	PICKLINE		49

		/* The fraction of the bond that determines if the */
		/* user has definitely picked one side of the bond */
#define	PICKLINESIDE		0.3

		/* The length of the coordinate axis in ANGSTROMS */
#define	AXISLENGTH	3.0


		/* Number of milli seconds between mouse clicks */
		/* for multi clicking selections */

#define	SELECT_DELAY	500


extern void	TankToggleDrawRubberBand(TANK tTank);


#endif /* _TankP_h */







