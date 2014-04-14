/*
 *      File:   minimizer.h
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
 *              The MINIMIZER object contains atoms, bonds, angles, and
 *              torsions.  Its main function is to calculate
 *              the minimizer of the bonds, angles, and torsions 
 *              and minimize it, removing strain from the
 *              molecular model.
 */

#ifndef MINIMIZER_H
#define MINIMIZER_H

#ifndef VARARRAY_H
#include        "varArray.h"
#endif
#include        "atom.h"

typedef struct  {
	double          dMinRmsGradientSquared;
	VARARRAY        vaAtoms;
	VARARRAY        vaBonds;
	VARARRAY        vaAngles;
	VARARRAY        vaTorsions;
    /* The following fields are used to monitor the minimizer */
	BOOL		(*bFCallback)();
	GENP		PCallbackData;
	double		dEnergy;
	double		dRmsGradient;
	BOOL		bMinimizing;
} MINIMIZERt;
                
typedef MINIMIZERt	*MINIMIZER;



/*
----------------------------------------------------------------------

        Messages
*/

extern MINIMIZER	mMinimizerCreate();
extern void		MinimizerDestroy(MINIMIZER *mPMinimizer);
extern void		MinimizerAddAtom(MINIMIZER mMinimizer, ATOM aAtom);
extern BOOL		bMinimizerAddBond(MINIMIZER mMinimizer, 
				ATOM aAtom1, ATOM aAtom2, 
				double dKb, double dR0);
extern BOOL		bMinimizerAddAngle(MINIMIZER mMinimizer, 
				ATOM aAtom1, ATOM aAtom2, ATOM aAtom3,
				double dKt, double dT0 );
extern BOOL		bMinimizerAddTorsion(MINIMIZER mMinimizer,
				ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, ATOM aAtom4,
				double dN, double dKp, double dP0);
extern void		MinimizerMinimize(MINIMIZER mMinimizer);

#define MinimizerSetMinRms( m, r )	( m->dMinRmsGradientSquared = r*r )

	/* If the Callback function is not NULL then every iteration */
	/* of the minimization routine will result in a call to */
	/* the Callback routine, allowing the caller to display the */
	/* progress of the minimization.  If the Callback routine returns */
	/* TRUE then the minimization will continue, otherwise it will  */
	/* terminate.   Call the callback with the data the caller provides. */

#define	MinimizerSetCallback( m, c, d )	( m->bFCallback = c, m->PCallbackData = (GENP)(d) )
#define dMinimizerEnergy( m )		( m->dEnergy )
#define	dMinimizerCurrentRms(m)		( m->dRmsGradient )

#define	bMinimizerMinimizing(m)		( m->bMinimizing )




#endif  /* ifdef MINIMIZER_H */
