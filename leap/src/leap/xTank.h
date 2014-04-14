#ifndef _XTank_h
#define _XTank_h

/*
 *	File:	xTank.h
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


#include	"x3d.h"

#include	"xTankP.h"


/****************************************************************
 *
 * Tank widget
 *
 ****************************************************************/

/* Resources:

 Name		     Class		RepType		Default Value
 ----		     -----		-------		-------------
 background	     Background		Pixel		XtDefaultBackground
 border		     BorderColor	Pixel		XtDefaultForeground
 borderWidth	     BorderWidth	Dimension	1
 destroyCallback     Callback		Pointer		NULL
 height		     Height		Dimension	0
 mappedWhenManaged   MappedWhenManaged	Boolean		True
 sensitive	     Sensitive		Boolean		True
 width		     Width		Dimension	0
 x		     Position		Position	0
 y		     Position		Position	0
*/






/* declare specific TemplateWidget class and instance datatypes */

typedef struct _TankClassRec	*TankWidgetClass;
typedef struct _TankRec		*TankWidget;

/* declare the class constant */

extern WidgetClass tankWidgetClass;


		/* Used to define the state of the TANK, what happens */
		/* when the user clicks the pointer button */

#define	TANKNOSTATE		0
#define	TANKDRAW		1
#define	TANKERASE		2
#define	TANKSELECT		3
#define	TANKDRAGROTATE		4
#define	TANKTWIST		5

#define	TANKSELECTNOTHING	0
#define	TANKSELECTATOM		1
#define	TANKSELECTRESIDUE	2
#define	TANKSELECTMOLECULE	3
#define	TANKSELECTALL		4


		/* Control various attributes of the TANK */

#define	TANKDEFAULTFLAGS	0x00000000
#define	TANKSHOWRESIDUES	0x00000001
#define	TANKSHOWNAMES		0x00000002
#define	TANKSHOWTYPES		0x00000004
#define	TANKSHOWCHARGES		0x00000008
#define	TANKSHOWAXIS		0x00000010
#define	TANKSHOWBOX		0x00000020
#define	TANKSHOWELEMENTS	0x00000040
#define	TANKSHOWPERTNAMES	0x00000080
#define	TANKSHOWPERTTYPES	0x00000100
#define	TANKDRAWCIRCLE		0x10000000



/*
 *-------------------------------------------------------------
 *
 *	TANK functions
 *
 */


extern void	TankUseUnit(TANK tTank, UNIT uUnit);
extern void	TankFastRedisplayUnit(TANK tTank);
extern void	TankRedisplayUnit(TANK tTank);
extern void	TankRedisplayRecenteredRescaledUnit(TANK tTank);
extern UNIT	uTankUnit(TANK tTank);

extern void	TankBeep(TANK tTank);

extern void	TankSetState(TANK tTank, int iNewState);
extern int	iTankState(TANK tTank);

extern void	TankSetFlags(TANK tTank, FLAGS fFlags);
extern void	TankResetFlags(TANK tTank, FLAGS fFlags);
extern BOOL	bTankFlagsSet(TANK tTank, FLAGS fFlags);
extern X3DENGINE	x3dTankX3dEngine(TANK tTank);

extern void	TankSetDrawingElement(TANK tTank, int iElement);
extern int	iTankDrawingElement(TANK tTank);

extern void	TankDefinePrintSink(TANK tTank, int iSink);


		/* The following are not to be called by programmers */
		/* they should be passed to  */
		/* AtomClassDefineGraphicsCreator and */
		/* AtomClassDefineGraphicsDestructor */
 
extern void	TankAtomGraphicsDataCreator(ATOM aAtom);
extern void	TankAtomGraphicsDataDestructor(ATOM aAtom);
extern BOOL	bTankPointInCircle( TANK tTank, int iX, int iY );
extern void	TankDisplay( TANK tTank, BOOL bFast );
extern void	TankToggleDrawRubberBox( TANK tTank );

#endif /* _XTank_h */
