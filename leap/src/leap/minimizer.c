/*
 *      File:   minimizer.c
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
 *              torsions.  It's main function is to minimize strain
 *              in the molecular model using a quick and dirty minimizer
 *              minimizer.
 *
 *              The minimizer function has the form:
 *
 *              E = ALL_BOND_SUM( Kb (r-r0)^2 ) +
 *                  ALL_ANGLE_SUM( Kt ( t-t0 )^2 ) +
 *                  ALL_TORSION_SUM( Kp ( COS( n*p ) )
 *
 *              The caller must specify the atoms involved
 *              in the bond, angle, torsion and the force constants,
 *              equilibrium values, and multiplicity.
 *
 *              The minimizer modifies the actual coordinates of the
 *              atoms directly, except for atoms that have the
 *              ATOMFIXED flag set.  For these atoms, the forces are
 *              zeroed.
 *
 *              When adding bonds, angles, torsions to the minimizer,
 *              the atoms for those objects must already have been
 *              added.
 *
 *		MINIMIZERs have a special flag that are TRUE when
 *		they are actively minimizing a structure.
 *		This allows code in a graphical user environment to
 *		check whether a MINIMIZER is active and only
 *		suspended while the callback is being envoked,
 *		this is to provide a hook where the program
 *		can prevent a MINIMIZER from being shut down while
 *		it is active
 */

/*
#define		DEBUGINNER	1

#define		DEBUGTORSION	1
#define		DEBUGANGLE	1
#define		DEBUGBOND	1

#define		DEBUGTORSION	1

*/





#include	"basics.h"

#include        "varArray.h"
#include        "nVector.h"
#include        "classes.h"

#include        "minimizer.h"


#define NOTFOUND                        -1
//#define MAXSTEEPESTDESCENTSTEPS         10000           /* Should be 200 */
#define MAXCONJUGATEGRADIENTSTEPS       900             /* Should be 200 */
#define MAXLINESEARCHSTEPS              10
#define STARTSTEPSIZE                   0.0001
#define MINRMSSQUARED                   0.025

typedef struct  {
                ATOM            aAtom;
                } EATOMt;
                
typedef struct  {
                int             iAtom1;
                int             iAtom2;
                double          dKb;
                double          dR0;
                } EBONDt;
                
typedef struct  {
                int             iAtom1;
                int             iAtom2;
                int             iAtom3;
                double          dKt;
                double          dT0;
                } EANGLEt;
                
typedef struct  {
                int             iAtom1;
                int             iAtom2;
                int             iAtom3;
                int             iAtom4;
                double          dKp;
                double          dN;
                } ETORSIONt;



/*
-------------------------------------------------------------------------

        Private routines
        
*/



/*
 *      iFindAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the index of the atom in the MINIMIZER object.
 *      Return NOTFOUND if it is not found.
 */
static int
iFindAtom( MINIMIZER mMinimizer, ATOM aAtom )
{
int     i, iMax;
EATOMt	*eaPAtom;

    iMax = iVarArrayElementCount(mMinimizer->vaAtoms);
    if ( !iMax )
	return(NOTFOUND);
    eaPAtom = PVAI( mMinimizer->vaAtoms, EATOMt, 0 );
    for ( i=0; i<iMax; i++, eaPAtom++ ) {
        if ( aAtom == eaPAtom->aAtom ) 
		return(i);
    }
    return(NOTFOUND);
}





/*
 *      nvGetPositionVector
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a NVECTOR that has the position coordinates within it.
 */
static NVECTOR 
nvGetPositionVector( MINIMIZER mMinimizer )
{
EATOMt		*eaPAtom;
int             i, iMax;
NVECTOR         nvPosition;
VECTOR          vPos;

    MESSAGE(( "|||Before defining nvPosition\n" ));
    TESTMEMORY();
    MESSAGE(( "|||Passed\n" ));
    
    iMax = iVarArrayElementCount(mMinimizer->vaAtoms);
    nvPosition = nvNVectorCreate( iMax * 3 );
    MESSAGE(( "|||After defining nvPosition\n" ));
    TESTMEMORY();
    MESSAGE(( "|||Passed\n" ));
    
                /* Loop over all atoms and obtain the positions */
    if ( iMax ) {
        eaPAtom = PVAI( mMinimizer->vaAtoms, EATOMt, 0 );
    	for ( i=0; i<iMax; i++, eaPAtom++ ) {

            vPos = vAtomPosition( eaPAtom->aAtom );
            NVectorSetElement( nvPosition, i*3 + 0, dVX(&vPos) );
            NVectorSetElement( nvPosition, i*3 + 1, dVY(&vPos) );
            NVectorSetElement( nvPosition, i*3 + 2, dVZ(&vPos) );
        
    	}
    }
    MESSAGE(( "|||After filling nvPosition\n" ));
    TESTMEMORY();
    MESSAGE(( "|||Passed\n" ));
    return(nvPosition);
}



/*
 *	SetAtomCoordinates
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Copy the contents of the position vector into the atoms
 *      coordinates.
 */
static void
SetAtomCoordinates( MINIMIZER mMinimizer, NVECTOR nvPosition )
{
EATOMt		*eaPAtom;
int             i, iMax;
VECTOR          vPos;

                /* Loop over all atoms and set their positions */

    iMax = iVarArrayElementCount( mMinimizer->vaAtoms );
    if ( iMax ) {
    	eaPAtom = PVAI(mMinimizer->vaAtoms, EATOMt, 0 );
    	for ( i=0; i<iMax; i++, eaPAtom++ ) {

            VectorDef( &vPos, dNVectorElement( nvPosition, i*3 + 0 ), 
                              dNVectorElement( nvPosition, i*3 + 1 ), 
                              dNVectorElement( nvPosition, i*3 + 2 ) );
            AtomSetPosition( eaPAtom->aAtom, vPos );
    	}
    }
}


    
    
    

/*
 *      dBondMinimizer
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Evaluate the bond minimizer and the first derivative
 *      Evalueate the minimizer at the position nvPosition, and
 *      add the first derivatives to nvDerivative
 */
static double
dBondMinimizer( MINIMIZER mMinimizer, NVECTOR nvPosition, NVECTOR nvDerivative )
{
double		vv0, vv1, vv2, vv3, vv4, vv5, vv6, vv7, vv8, vv9;
double		vv10, vv11, vv12, vv13, vv14;
int             iAtom1, iAtom2;
EATOMt		*eaPAtom1, *eaPAtom2;
ATOM            aAtom1, aAtom2;
BOOL            bMove1, bMove2;
EBONDt		*ebPBond;
int             i, iBonds, iAtom1Index, iAtom2Index;
double          x1,y1,z1, x2,y2,z2;
double          dx1,dy1,dz1, dx2,dy2,dz2;
double          e, Kb, R0;

    e = 0.0;
                /* Loop over all bond interactions */

    iBonds = iVarArrayElementCount(mMinimizer->vaBonds);
    if ( !iBonds )
	return(e);

    ebPBond = PVAI(mMinimizer->vaBonds, EBONDt, 0);
    for ( i=0; i<iBonds; i++, ebPBond++ ) {

                /* Obtain all the parameters necessary to calculate */
                /* the minimizer and forces */

        Kb = ebPBond->dKb;
        R0 = ebPBond->dR0;
        iAtom1 = ebPBond->iAtom1;
        iAtom1Index = iAtom1*3;
        iAtom2 = ebPBond->iAtom2;
        iAtom2Index = iAtom2*3;

        eaPAtom1 = PVAI( mMinimizer->vaAtoms, EATOMt, iAtom1 );
        eaPAtom2 = PVAI( mMinimizer->vaAtoms, EATOMt, iAtom2 );
        aAtom1 = eaPAtom1->aAtom;
        aAtom2 = eaPAtom2->aAtom;

#ifdef	DEBUGBOND
        MESSAGE(( "Bond interaction %s - %s  Kb=%lf  R0=%lf\n",
                        sContainerName(aAtom1), sContainerName(aAtom2),
                        Kb, R0 ));
#endif

	bMove1 = bAtomFlagsSet( aAtom1, ATOMNEEDSMINIMIZER );
	bMove2 = bAtomFlagsSet( aAtom2, ATOMNEEDSMINIMIZER );

        x1 = dNVectorElement( nvPosition, iAtom1Index + 0 );
        y1 = dNVectorElement( nvPosition, iAtom1Index + 1 );
        z1 = dNVectorElement( nvPosition, iAtom1Index + 2 );
        x2 = dNVectorElement( nvPosition, iAtom2Index + 0 );
        y2 = dNVectorElement( nvPosition, iAtom2Index + 1 );
        z2 = dNVectorElement( nvPosition, iAtom2Index + 2 );

        vv0 = -R0 ;
        vv1 = -x2 ;
        vv2 = -y2 ;
        vv3 = -z2 ;
        vv4 = vv1 + x1 ;
        vv5 = vv2 + y1 ;
        vv6 = vv3 + z1 ;
        vv7 = vv4*vv4;
        vv8 = vv5*vv5 ;
        vv9 = vv6*vv6 ;
        vv10 = vv7 + vv8 + vv9 ;
        vv12 = sqrt(vv10) ;
        vv11 = 1/vv12 ;
        vv13 = vv0 + vv12 ;
        vv14 = 2*Kb*vv11*vv13;

                /* Calculate the minimizer and the forces */
        e += Kb*myPow(vv13,2.0) ;
        dx1 = vv14*vv4 ;
        dy1 = vv14*vv5 ;
        dz1 = vv14*vv6 ;
        dx2 = -vv14*vv4 ;
        dy2 = -vv14*vv5 ;
        dz2 = -vv14*vv6 ;

#ifdef	DEBUGBOND
        MESSAGE(( "x1 = %lf\n", x1 ));
        MESSAGE(( "y1 = %lf\n", y1 ));
        MESSAGE(( "z1 = %lf\n", z1 ));
        MESSAGE(( "x2 = %lf\n", x2 ));
        MESSAGE(( "y2 = %lf\n", y2 ));
        MESSAGE(( "z2 = %lf\n", z2 ));

        MESSAGE(( "Energy = %lf\n", e ));
        MESSAGE(( "dx1 = %lf\n", dx1 ));
        MESSAGE(( "dy1 = %lf\n", dy1 ));
        MESSAGE(( "dz1 = %lf\n", dz1 ));
        MESSAGE(( "dx2 = %lf\n", dx2 ));
        MESSAGE(( "dy2 = %lf\n", dy2 ));
        MESSAGE(( "dz2 = %lf\n", dz2 ));
#endif
                /* Add the forces */

        if ( bMove1 ) {
            NVectorSetElement( nvDerivative, iAtom1Index+0, 
                        dx1 + dNVectorElement( nvDerivative, iAtom1Index+0 ));
            NVectorSetElement( nvDerivative, iAtom1Index+1,
                        dy1 + dNVectorElement( nvDerivative, iAtom1Index+1 ));
            NVectorSetElement( nvDerivative, iAtom1Index+2,
                        dz1 + dNVectorElement( nvDerivative, iAtom1Index+2 ));
        }
        if ( bMove2 ) {
            NVectorSetElement( nvDerivative, iAtom2Index+0, 
                        dx2 + dNVectorElement( nvDerivative, iAtom2Index+0 ));
            NVectorSetElement( nvDerivative, iAtom2Index+1,
                        dy2 + dNVectorElement( nvDerivative, iAtom2Index+1 ));
            NVectorSetElement( nvDerivative, iAtom2Index+2,
                        dz2 + dNVectorElement( nvDerivative, iAtom2Index+2 ));
        }
    }

                /* return the minimizer */
    return(e);
}


            
/*
 *      dAngleMinimizer
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Evaluate the angle minimizer and the first derivative
 *      Evalueate the minimizer at the position nvPosition, and
 *      add the first derivatives to nvDerivative
 */
static double
dAngleMinimizer( MINIMIZER mMinimizer, NVECTOR nvPosition, NVECTOR nvDerivative )
{
double  	vv0, vv1, vv2, vv3, vv4, vv5, vv6, vv7, vv8, vv9 ;
double  	vv10, vv11, vv12, vv13, vv14, vv15, vv16, vv17, vv18, vv19 ;
double  	vv20, vv21, vv22, vv23, vv24, vv25, vv26, vv27, vv28, vv29 ;
double  	vv30, vv31, vv32, vv33, vv34, vv35 ;
double          Kt, T0;
int             iAtom1, iAtom2, iAtom3;
EATOMt		*eaPAtom1, *eaPAtom2, *eaPAtom3;
ATOM            aAtom1, aAtom2, aAtom3;
BOOL            bMove1, bMove2, bMove3;
EANGLEt		*eaPAngle;
int             i, iAngles, iAtom1Index, iAtom2Index, iAtom3Index;
double          x1,y1,z1, x2,y2,z2, x3,y3,z3;
double          dx1,dy1,dz1, dx2,dy2,dz2, dx3,dy3,dz3;
double          e;

    e = 0.0;
                /* Loop over all angle interactions */

    iAngles = iVarArrayElementCount(mMinimizer->vaAngles);
    if ( !iAngles )
	return(e);

    eaPAngle = PVAI(mMinimizer->vaAngles, EANGLEt, 0);
    for ( i=0; i<iAngles; i++, eaPAngle++ ) {

                /* Obtain all the parameters necessary to calculate */
                /* the minimizer and forces */

        Kt = eaPAngle->dKt;
        T0 = eaPAngle->dT0;
        iAtom1 = eaPAngle->iAtom1;
        iAtom2 = eaPAngle->iAtom2;
        iAtom3 = eaPAngle->iAtom3;
        iAtom1Index = iAtom1*3;
        iAtom2Index = iAtom2*3;
        iAtom3Index = iAtom3*3;
        eaPAtom1 = PVAI(mMinimizer->vaAtoms, EATOMt, iAtom1);
        eaPAtom2 = PVAI(mMinimizer->vaAtoms, EATOMt, iAtom2);
        eaPAtom3 = PVAI(mMinimizer->vaAtoms, EATOMt, iAtom3);
        aAtom1 = eaPAtom1->aAtom;
        aAtom2 = eaPAtom2->aAtom;
        aAtom3 = eaPAtom3->aAtom;
        
#ifdef	DEBUGANGLE
        MESSAGE(( "Angle interaction %s - %s - %s Kt=%lf  T0=%lf\n",
                        sContainerName(aAtom1), sContainerName(aAtom2),
                        sContainerName(aAtom3),
                        Kt, T0 ));
#endif
        
        bMove1 = bAtomFlagsSet( aAtom1, ATOMNEEDSMINIMIZER );
        bMove2 = bAtomFlagsSet( aAtom2, ATOMNEEDSMINIMIZER );
        bMove3 = bAtomFlagsSet( aAtom3, ATOMNEEDSMINIMIZER );

        if ( !bMove1 && !bMove2 && !bMove3 ) continue;

        x1 = dNVectorElement( nvPosition, iAtom1Index + 0 );
        y1 = dNVectorElement( nvPosition, iAtom1Index + 1 );
        z1 = dNVectorElement( nvPosition, iAtom1Index + 2 );

        x2 = dNVectorElement( nvPosition, iAtom2Index + 0 );
        y2 = dNVectorElement( nvPosition, iAtom2Index + 1 );
        z2 = dNVectorElement( nvPosition, iAtom2Index + 2 );

        x3 = dNVectorElement( nvPosition, iAtom3Index + 0 );
        y3 = dNVectorElement( nvPosition, iAtom3Index + 1 );
        z3 = dNVectorElement( nvPosition, iAtom3Index + 2 );

        vv0 = -T0 ;
        vv1 = -x2 ;
        vv2 = -y2 ;
        vv3 = -z2 ;
        vv4 = vv1 + x1 ;
        vv5 = vv1 + x3 ;
        vv6 = vv2 + y1 ;
        vv7 = vv2 + y3 ;
        vv8 = vv3 + z1 ;
        vv9 = vv3 + z3 ;
        vv10 = vv4*vv4; 
        vv11 = vv4*vv5 ;
        vv12 = vv5*vv5;
        vv13 = vv6*vv6;
        vv14 = vv6*vv7 ;
        vv15 = vv7*vv7;
        vv16 = vv8*vv8;
        vv17 = vv8*vv9 ;
        vv18 = vv9*vv9;
        vv19 = vv10 + vv13 + vv16 ;
        vv20 = vv11 + vv14 + vv17 ;
        vv21 = vv12 + vv15 + vv18 ;
        vv22 = myPow(vv19,(-3)/2.0) ;
        vv23 = 1/vv19 ;
        vv24 = 1/sqrt(vv19) ;
        vv25 = myPow(vv20,2.0) ;
        vv26 = myPow(vv21,(-3)/2.0) ;
        vv27 = 1/vv21 ;
        vv28 = 1/sqrt(vv21) ;
        vv29 = -(vv23*vv25*vv27) ;
        vv30 = vv20*vv24*vv28 ;
        vv31 = 1 + vv29 ;
        vv32 = myAcos(vv30) ;
        vv33 = 1/sqrt(vv31) ;
        vv34 = vv0 + vv32 ;
        vv35 = 2*Kt*vv33*vv34;
        e += Kt*myPow(vv34,2.0) ;
        dx1 = -vv35*(-(vv20*vv22*vv28*vv4) + vv24*vv28*vv5) ;
        dy1 = -vv35*(-(vv20*vv22*vv28*vv6) + vv24*vv28*vv7) ;
        dz1 = -vv35*(-(vv20*vv22*vv28*vv8) + vv24*vv28*vv9) ;
        dx2 = -vv35*(vv20*vv22*vv28*vv4 + vv24*vv28*(-vv4 - vv5) +  
               vv20*vv24*vv26*vv5) ;
        dy2 = -vv35*(vv20*vv22*vv28*vv6 + vv24*vv28*(-vv6 - vv7) +  
               vv20*vv24*vv26*vv7) ;
        dz2 = -vv35*(vv20*vv22*vv28*vv8 + vv24*vv28*(-vv8 - vv9) +  
               vv20*vv24*vv26*vv9) ;
        dx3 = -vv35*(vv24*vv28*vv4 - vv20*vv24*vv26*vv5) ;
        dy3 = -vv35*(vv24*vv28*vv6 - vv20*vv24*vv26*vv7) ;
        dz3 = -vv35*(vv24*vv28*vv8 - vv20*vv24*vv26*vv9) ;

#ifdef	DEBUGANGLE
        MESSAGE(( "x1 = %lf\n", x1 ));
        MESSAGE(( "y1 = %lf\n", y1 ));
        MESSAGE(( "z1 = %lf\n", z1 ));
        MESSAGE(( "x2 = %lf\n", x2 ));
        MESSAGE(( "y2 = %lf\n", y2 ));
        MESSAGE(( "z2 = %lf\n", z2 ));
        MESSAGE(( "x3 = %lf\n", x3 ));
        MESSAGE(( "y3 = %lf\n", y3 ));
        MESSAGE(( "z3 = %lf\n", z3 ));

        MESSAGE(( "Energy = %lf\n", e ));
        MESSAGE(( "dx1 = %lf\n", dx1 ));
        MESSAGE(( "dy1 = %lf\n", dy1 ));
        MESSAGE(( "dz1 = %lf\n", dz1 ));
        MESSAGE(( "dx2 = %lf\n", dx2 ));
        MESSAGE(( "dy2 = %lf\n", dy2 ));
        MESSAGE(( "dz2 = %lf\n", dz2 ));
        MESSAGE(( "dx3 = %lf\n", dx3 ));
        MESSAGE(( "dy3 = %lf\n", dy3 ));
        MESSAGE(( "dz3 = %lf\n", dz3 ));
#endif
                /* Add the forces */

        if ( bMove1 ) {
            NVectorSetElement( nvDerivative, iAtom1Index+0, 
                        dx1 + dNVectorElement( nvDerivative, iAtom1Index+0 ));
            NVectorSetElement( nvDerivative, iAtom1Index+1,
                        dy1 + dNVectorElement( nvDerivative, iAtom1Index+1 ));
            NVectorSetElement( nvDerivative, iAtom1Index+2,
                        dz1 + dNVectorElement( nvDerivative, iAtom1Index+2 ));
        }
        if ( bMove2 ) {
            NVectorSetElement( nvDerivative, iAtom2Index+0, 
                        dx2 + dNVectorElement( nvDerivative, iAtom2Index+0 ));
            NVectorSetElement( nvDerivative, iAtom2Index+1,
                        dy2 + dNVectorElement( nvDerivative, iAtom2Index+1 ));
            NVectorSetElement( nvDerivative, iAtom2Index+2,
                        dz2 + dNVectorElement( nvDerivative, iAtom2Index+2 ));
        }
        if ( bMove3 ) {
            NVectorSetElement( nvDerivative, iAtom3Index+0, 
                        dx3 + dNVectorElement( nvDerivative, iAtom3Index+0 ));
            NVectorSetElement( nvDerivative, iAtom3Index+1,
                        dy3 + dNVectorElement( nvDerivative, iAtom3Index+1 ));
            NVectorSetElement( nvDerivative, iAtom3Index+2,
                        dz3 + dNVectorElement( nvDerivative, iAtom3Index+2 ));
        }
    }
                /* return the minimizer */
    return(e);
}






/*
 *      dTorsionMinimizer
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Evaluate the torsion minimizer and the first derivative
 *      Evalueate the minimizer at the position nvPosition, and
 *      add the first derivatives to nvDerivative
 */
static double
dTorsionMinimizer( MINIMIZER mMinimizer, NVECTOR nvPosition, NVECTOR nvDerivative )
{
double  	vv0, vv1, vv2, vv3, vv4, vv5, vv6, vv7, vv8, vv9 ;
double  	vv10, vv11, vv12, vv13, vv14, vv15, vv16, vv17, vv18, vv19 ;
double  	vv20, vv21, vv22, vv23, vv24, vv25, vv26, vv27, vv28, vv29 ;
double  	vv30, vv31, vv32, vv33, vv34, vv35, vv36, vv37, vv38, vv39 ;
double  	vv40, vv41, vv42, vv43, vv44, vv45, vv46, vv47, vv48, vv49 ;
double  	vv50, vv51, vv52, vv53, vv54, vv55, vv56, vv57, vv58, vv59 ;
double  	vv60, vv61, vv62, vv63, vv64, vv65, vv66, vv67, vv68, vv69 ;
double  	vv70, vv71, vv72, vv73, vv74, vv75, vv76, vv77, vv78, vv79 ;
double  	vv80, vv81, vv82, vv83, vv84, vv85 ;
double          e, Kp;
int             iAtom1, iAtom2, iAtom3, iAtom4;
int             iAtom1Index, iAtom2Index, iAtom3Index, iAtom4Index;
EATOMt		*eaPAtom1, *eaPAtom2, *eaPAtom3, *eaPAtom4;
ATOM            aAtom1, aAtom2, aAtom3, aAtom4;
BOOL            bMove1, bMove2, bMove3, bMove4;
ETORSIONt	*etPTor;
int             i, iTorsions;
double          N;
double          x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4;
double          dx1,dy1,dz1, dx2,dy2,dz2, dx3,dy3,dz3, dx4,dy4,dz4;


    e = 0.0;

    iTorsions = iVarArrayElementCount(mMinimizer->vaTorsions);
    if ( !iTorsions )
	return(e);

                /* Loop over all torsion interactions */

    etPTor = PVAI( mMinimizer->vaTorsions, ETORSIONt, 0 );
    for ( i=0; i<iTorsions; i++, etPTor++ ) {

                /* Obtain all the parameters necessary to calculate */
                /* the minimizer and forces */

        Kp = etPTor->dKp;
        N  = etPTor->dN;
        iAtom1 = etPTor->iAtom1;
        iAtom2 = etPTor->iAtom2;
        iAtom3 = etPTor->iAtom3;
        iAtom4 = etPTor->iAtom4;
        iAtom1Index = iAtom1*3;
        iAtom2Index = iAtom2*3;
        iAtom3Index = iAtom3*3;
        iAtom4Index = iAtom4*3;
        eaPAtom1 = PVAI( mMinimizer->vaAtoms, EATOMt, iAtom1 );
        eaPAtom2 = PVAI( mMinimizer->vaAtoms, EATOMt, iAtom2 );
        eaPAtom3 = PVAI( mMinimizer->vaAtoms, EATOMt, iAtom3 );
        eaPAtom4 = PVAI( mMinimizer->vaAtoms, EATOMt, iAtom4 );
        aAtom1 = eaPAtom1->aAtom;
        aAtom2 = eaPAtom2->aAtom;
        aAtom3 = eaPAtom3->aAtom;
        aAtom4 = eaPAtom4->aAtom;

#ifdef	DEBUGTORSION
        MESSAGE(( "Torsion interaction %s - %s - %s - %s Kp=%lf  N=%lf\n",
                        sContainerName(aAtom1), sContainerName(aAtom2),
                        sContainerName(aAtom3), sContainerName(aAtom4),
                        Kp, N ));
#endif
        
        bMove1 = bAtomFlagsSet( aAtom1, ATOMNEEDSMINIMIZER );
        bMove2 = bAtomFlagsSet( aAtom2, ATOMNEEDSMINIMIZER );
        bMove3 = bAtomFlagsSet( aAtom3, ATOMNEEDSMINIMIZER );
        bMove4 = bAtomFlagsSet( aAtom4, ATOMNEEDSMINIMIZER );

        if ( !bMove1 && !bMove2 && !bMove3 && !bMove4 ) continue;

        x1 = dNVectorElement( nvPosition, iAtom1Index + 0 );
        y1 = dNVectorElement( nvPosition, iAtom1Index + 1 );
        z1 = dNVectorElement( nvPosition, iAtom1Index + 2 );

        x2 = dNVectorElement( nvPosition, iAtom2Index + 0 );
        y2 = dNVectorElement( nvPosition, iAtom2Index + 1 );
        z2 = dNVectorElement( nvPosition, iAtom2Index + 2 );

        x3 = dNVectorElement( nvPosition, iAtom3Index + 0 );
        y3 = dNVectorElement( nvPosition, iAtom3Index + 1 );
        z3 = dNVectorElement( nvPosition, iAtom3Index + 2 );

        x4 = dNVectorElement( nvPosition, iAtom4Index + 0 );
        y4 = dNVectorElement( nvPosition, iAtom4Index + 1 );
        z4 = dNVectorElement( nvPosition, iAtom4Index + 2 );

        dx1 = 0.0; dy1 = 0.0; dz1 = 0.0;
        dx2 = 0.0; dy2 = 0.0; dz2 = 0.0;
        dx3 = 0.0; dy3 = 0.0; dz3 = 0.0;
        dx4 = 0.0; dy4 = 0.0; dz4 = 0.0;

        vv0 = -x2 ;
        vv1 = -x3 ;
        vv2 = -y2 ;
        vv3 = -y3 ;
        vv4 = -z2 ;
        vv5 = -z3 ;
        vv6 = vv0 + x1 ;
        vv7 = vv1 + x2 ;
        vv8 = vv0 + x3 ;
        vv9 = vv1 + x4 ;
        vv10 = vv2 + y1 ;
        vv11 = vv3 + y2 ;
        vv12 = vv2 + y3 ;
        vv13 = vv3 + y4 ;
        vv14 = vv4 + z1 ;
        vv15 = vv5 + z2 ;
        vv16 = vv4 + z3 ;
        vv17 = vv5 + z4 ;
        vv18 = -vv10 ;
        vv19 = -vv11 ;
        vv20 = -vv12 ;
        vv21 = -vv13 ;
        vv22 = -vv14 ;
        vv23 = -vv15 ;
        vv24 = -vv16 ;
        vv25 = vv10*vv16 ;
        vv26 = -vv17 ;
        vv27 = vv11*vv17 ;
        vv28 = -vv6 ;
        vv29 = vv12*vv6 ;
        vv30 = -vv7 ;
        vv31 = vv13*vv7 ;
        vv32 = -vv8 ;
        vv33 = vv14*vv8 ;
        vv34 = -vv9 ;
        vv35 = vv15*vv9 ;
        vv36 = vv12 + vv18 ;
        vv37 = vv13 + vv19 ;
        vv38 = vv14*vv20 ;
        vv39 = vv10 + vv20 ;
        vv40 = vv15*vv21 ;
        vv41 = vv11 + vv21 ;
        vv42 = vv16 + vv22 ;
        vv43 = vv17 + vv23 ;
        vv44 = vv14 + vv24 ;
        vv45 = vv15 + vv26 ;
        vv46 = vv24*vv6 ;
        vv47 = vv32 + vv6 ;
        vv48 = vv26*vv7 ;
        vv49 = vv34 + vv7 ;
        vv50 = vv18*vv8 ;
        vv51 = vv28 + vv8 ;
        vv52 = vv19*vv9 ;
        vv53 = vv30 + vv9 ;
        vv54 = vv25 + vv38 ;
        vv55 = vv27 + vv40 ;
        vv56 = vv33 + vv46 ;
        vv57 = vv35 + vv48 ;
        vv58 = vv29 + vv50 ;
        vv59 = vv31 + vv52 ;
        vv60 = myPow(vv54,2.0) ;
        vv61 = vv54*vv55 ;
        vv62 = myPow(vv55,2.0) ;
        vv63 = myPow(vv56,2.0) ;
        vv64 = vv56*vv57 ;
        vv65 = myPow(vv57,2.0) ;
        vv66 = myPow(vv58,2.0) ;
        vv67 = vv58*vv59 ;
        vv68 = myPow(vv59,2.0) ;
        vv69 = vv60 + vv63 + vv66 ;
        vv70 = vv61 + vv64 + vv67 ;
        vv71 = vv62 + vv65 + vv68 ;
        vv72 = myPow(vv69,(-3)/2.0) ;
        vv73 = 1/vv69 ;
        vv74 = 1/sqrt(vv69) ;
        vv75 = myPow(vv70,2.0) ;
        vv76 = myPow(vv71,(-3)/2.0) ;
        vv77 = 1/vv71 ;
        vv78 = 1/sqrt(vv71) ;
        vv79 = -(vv73*vv75*vv77) ;
        vv80 = vv70*vv74*vv78 ;
        vv81 = 1 + vv79 ;
        vv82 = myAcos(vv80) ;           /* vv82 = torsion angle */
        vv84 = N*vv82 ;
        e += Kp*cos(vv84) ;

                /* If vv81 is too close to zero then the derivatives */
                /* cannot be calculated */

        if ( -VERYSMALL < vv81 && VERYSMALL > vv81 ) goto ZERODERIV;
        vv83 = 1/sqrt(vv81) ;
        vv85 = sin(vv84) ;
        dx1 = Kp*N*(-((-2*vv16*vv56 + 2*vv12*vv58)*vv70*vv72*vv78)/2 +
                (vv24*vv57 + vv12*vv59)*vv74*vv78)*vv83*vv85 ;
        dy1 = Kp*N*((vv16*vv55 + vv32*vv59)*vv74*vv78 -
                 vv70*vv72*vv78*(2*vv16*vv54 - 2*vv58*vv8)/2.0)*vv83*vv85 ;
        dz1 = Kp*N*(-(vv70*vv72*vv78*(-2*vv12*vv54 + 2*vv56*vv8))/2 +
                 vv74*vv78*(vv20*vv55 + vv57*vv8))*vv83*vv85 ;
        dx2 = Kp*N*(-((-2*vv17*vv57 + 2*vv13*vv59)*vv70*vv74*vv76)/2 -
                 (2*vv42*vv56 + 2*vv39*vv58)*vv70*vv72*vv78/2 +
             (vv26*vv56 + vv42*vv57 + vv13*vv58 + vv39*vv59)
                *vv74*vv78)*vv83*vv85 ;
        dy2 = Kp*N*vv83*vv85*(-((2*vv44*vv54 + 2*vv51*vv58)*vv70*vv72*vv78)/2 +
              (vv17*vv54 + vv44*vv55 + vv34*vv58 + vv51*vv59)*vv74*vv78 -
               vv70*vv74*vv76*(2*vv17*vv55 - 2*vv59*vv9)/2.0) ;
        dz2 = Kp*N*vv83*vv85*(-((2*vv36*vv54 + 2*vv47*vv56)*vv70*vv72*vv78)/2 +
                vv74*vv78*(vv21*vv54 + vv36*vv55 + vv47*vv57 + vv56*vv9) -
                vv70*vv74*vv76*(-2*vv13*vv55 + 2*vv57*vv9)/2.0) ;
        dx3 = Kp*N*(-((2*vv43*vv57 + 2*vv41*vv59)*vv70*vv74*vv76)/2 -
              (2*vv14*vv56 - 2*vv10*vv58)*vv70*vv72*vv78/2 +
              (vv43*vv56 + vv14*vv57 + vv41*vv58 + vv18*vv59)
                *vv74*vv78)*vv83*vv85 ;
        dy3 = Kp*N*(-((2*vv45*vv55 + 2*vv53*vv59)*vv70*vv74*vv76)/2 -
              (-2*vv14*vv54 + 2*vv58*vv6)*vv70*vv72*vv78/2 +
              (vv45*vv54 + vv22*vv55 + vv53*vv58 + vv59*vv6)
                *vv74*vv78)*vv83*vv85 ;
        dz3 = Kp*N*(-((2*vv37*vv55 + 2*vv49*vv57)*vv70*vv74*vv76)/2 -
              (2*vv10*vv54 - 2*vv56*vv6)*vv70*vv72*vv78/2 +
              (vv37*vv54 + vv10*vv55 + vv49*vv56 + vv28*vv57)
                *vv74*vv78)*vv83*vv85 ;
        dx4 = Kp*N*(-((2*vv15*vv57 - 2*vv11*vv59)*vv70*vv74*vv76)/2 +
                 (vv15*vv56 + vv19*vv58)*vv74*vv78)*vv83*vv85 ;
        dy4 = Kp*N*(-((-2*vv15*vv55 + 2*vv59*vv7)*vv70*vv74*vv76)/2 + 
                 (vv23*vv54 + vv58*vv7)*vv74*vv78)*vv83*vv85 ;
        dz4 = Kp*N*(-((2*vv11*vv55 - 2*vv57*vv7)*vv70*vv74*vv76)/2 +
                 (vv11*vv54 + vv30*vv56)*vv74*vv78)*vv83*vv85 ;

ZERODERIV:

#ifdef	DEBUGTORSION

        MESSAGE(( "x1 = %lf\n", x1 ));
        MESSAGE(( "y1 = %lf\n", y1 ));
        MESSAGE(( "z1 = %lf\n", z1 ));
        MESSAGE(( "x2 = %lf\n", x2 ));
        MESSAGE(( "y2 = %lf\n", y2 ));
        MESSAGE(( "z2 = %lf\n", z2 ));
        MESSAGE(( "x3 = %lf\n", x3 ));
        MESSAGE(( "y3 = %lf\n", y3 ));
        MESSAGE(( "z3 = %lf\n", z3 ));
        MESSAGE(( "x4 = %lf\n", x4 ));
        MESSAGE(( "y4 = %lf\n", y4 ));
        MESSAGE(( "z4 = %lf\n", z4 ));

        MESSAGE(( "Energy = %lf\n", e ));
        MESSAGE(( "dx1 = %lf\n", dx1 ));
        MESSAGE(( "dy1 = %lf\n", dy1 ));
        MESSAGE(( "dz1 = %lf\n", dz1 ));
        MESSAGE(( "dx2 = %lf\n", dx2 ));
        MESSAGE(( "dy2 = %lf\n", dy2 ));
        MESSAGE(( "dz2 = %lf\n", dz2 ));
        MESSAGE(( "dx3 = %lf\n", dx3 ));
        MESSAGE(( "dy3 = %lf\n", dy3 ));
        MESSAGE(( "dz3 = %lf\n", dz3 ));
        MESSAGE(( "dx4 = %lf\n", dx4 ));
        MESSAGE(( "dy4 = %lf\n", dy4 ));
        MESSAGE(( "dz4 = %lf\n", dz4 ));
#endif
                /* Add the forces */

        if ( bMove1 ) {
            NVectorSetElement( nvDerivative, iAtom1Index+0, 
                        dx1 + dNVectorElement( nvDerivative, iAtom1Index+0 ));
            NVectorSetElement( nvDerivative, iAtom1Index+1,
                        dy1 + dNVectorElement( nvDerivative, iAtom1Index+1 ));
            NVectorSetElement( nvDerivative, iAtom1Index+2,
                        dz1 + dNVectorElement( nvDerivative, iAtom1Index+2 ));
        }
        if ( bMove2 ) {
            NVectorSetElement( nvDerivative, iAtom2Index+0, 
                        dx2 + dNVectorElement( nvDerivative, iAtom2Index+0 ));
            NVectorSetElement( nvDerivative, iAtom2Index+1,
                        dy2 + dNVectorElement( nvDerivative, iAtom2Index+1 ));
            NVectorSetElement( nvDerivative, iAtom2Index+2,
                        dz2 + dNVectorElement( nvDerivative, iAtom2Index+2 ));
        }
        if ( bMove3 ) {
            NVectorSetElement( nvDerivative, iAtom3Index+0, 
                        dx3 + dNVectorElement( nvDerivative, iAtom3Index+0 ));
            NVectorSetElement( nvDerivative, iAtom3Index+1,
                        dy3 + dNVectorElement( nvDerivative, iAtom3Index+1 ));
            NVectorSetElement( nvDerivative, iAtom3Index+2,
                        dz3 + dNVectorElement( nvDerivative, iAtom3Index+2 ));
        }
        if ( bMove4 ) {
            NVectorSetElement( nvDerivative, iAtom4Index+0, 
                        dx4 + dNVectorElement( nvDerivative, iAtom4Index+0 ));
            NVectorSetElement( nvDerivative, iAtom4Index+1,
                        dy4 + dNVectorElement( nvDerivative, iAtom4Index+1 ));
            NVectorSetElement( nvDerivative, iAtom4Index+2,
                        dz4 + dNVectorElement( nvDerivative, iAtom4Index+2 ));
        }
    }
                /* return the minimizer */
    return(e);
}



/*
 *      dTotalEnergy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Evaluate the total minimizer.
 *      Evaluate the minimizer at nvPos and return the derivative
 *      in nvDerivative (which must be an initialized vector,
 *      but does not need to be filled with zeros).
 */
static double
dTotalEnergy( MINIMIZER mMinimizer, NVECTOR nvPos, NVECTOR nvDerivative )
{
double          dBond, dAngle, dTorsion;
double          dEnergy;

                /* Zero the derivative */

    NVectorZero( nvDerivative );
    dBond    = dBondMinimizer( mMinimizer, nvPos, nvDerivative );
    dAngle   = dAngleMinimizer( mMinimizer, nvPos, nvDerivative );
    dTorsion = dTorsionMinimizer( mMinimizer, nvPos, nvDerivative );
MESSAGE(( "dBond    = %lf\n", dBond ));
MESSAGE(( "dAngle   = %lf\n", dAngle ));
MESSAGE(( "dTorsion = %lf\n", dTorsion ));

    dEnergy = dBond + dAngle + dTorsion;
    return(dEnergy);
}

    




/*
 *      dLineSearch
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Perform a line search of the minimizer function.
 */
static double
dLineSearch( MINIMIZER mMinimizer, double *dPF0, double *dPD0, 
		NVECTOR nvOrigin, NVECTOR nvDir, 
		NVECTOR nvNewPos, NVECTOR nvNewDeriv, 
		NVECTOR nvTempPos, NVECTOR nvTempDeriv )
{
#define EPS                     0.1             /* Convergence parameter */
#define FUNCTIONEVALUATIONMAX   20              /* Maximum number of times */
                                                /* to evaluate function */
                                                
double          dAcur, dFcur, dDcur;
double          dApre, dFpre, dDpre;
double          dAmin, dFmin, dDmin;
double          dTest;
double          dF0, dD0;
int             iFunctionEvaluationCount;
double          dStart;

    iFunctionEvaluationCount = 0;

    dF0   = *dPF0;
    dD0   = *dPD0;
    dStart = STARTSTEPSIZE;
    
    dApre = 0.0;
    dFpre = dF0;
    dDpre = dD0;

    dAmin = 0.0;
    dFmin = dF0;
    dDmin = dD0;
    
    dAcur = dStart*fabs(dD0);

	/* If the derivative is zero then just return */
	
    if ( fabs(dD0) < VERYSMALL ) return( 0.0 );
    
    while (1) {

        /* ---------------------------------------------------- */

                /* Evaluate the function at the test point */
                /* Evaluate at dAcur along the line */
                /* Place the function value in dFcur and the */
                /* component of the derivative along the line */
                /* in dDcur  */
        NVectorTimesScalar( nvNewPos, nvDir, dAcur/dD0 );
        NVectorAdd( nvNewPos, nvNewPos, nvOrigin );
        
        dFcur = dTotalEnergy( mMinimizer, nvNewPos, nvNewDeriv );


                /* Obtain the component of the derivative in the direction */
                /* of the line */

        dDcur = dNVectorDot( nvNewDeriv, nvDir )/dD0;

MESSAGE(( "    ZERO       @%lf  F= %lf D= %lf\n", 0.0, dF0, dD0 ));
MESSAGE(( "    PREVIOUS   @%lf  F= %lf D= %lf\n", dApre, dFpre, dDpre ));
MESSAGE(( "    Evaluation @%lf  F= %lf D= %lf\n", dAcur, dFcur, dDcur ));

#if 0
                /* Keep track of the minimum value seen */
        if ( dFmin > dFcur ) {
            dAmin = dAcur;
            dFmin = dFcur;
            dDmin = dDcur;
            NVectorCopy( nvTempPos, nvNewPos );
            NVectorCopy( nvTempDeriv, nvNewDeriv );
        }
#endif

        /* -------------------------------------------------- */

                /* Check if the maximum number of function evaluations */
                /* has been reached */

        iFunctionEvaluationCount++;
        if ( iFunctionEvaluationCount > FUNCTIONEVALUATIONMAX ) {
            MESSAGE(( "WARNING: Evaluated function MAX number of times!\n" ));
            goto DONE;
        }


        /* -------------------------------------------------- */

                /* Test if the new point has a negative slope but a higher */
                /* function value than dAcur = 0.0.  If this is the case,  */
                /* the search has passed through a local max and is heading */
                /* for a distant local minimum */

        if ( (dFcur > dF0) && (dDcur < 0.0) ) {

MESSAGE(( "Passed through a local maximum, moving back\n" ));
            if ( dTest > 0.0 ) {
                        /* Reduce dAcur and restart */

                dTest = dAcur / 3.0;
                goto CONTINUE;
            }
        }


                /* Test whether the steplength criteria have been met */
                /* Goldstein test followed by Wolfe test */

        if ( dFcur > ( dF0 + 0.0001*dAcur*dD0 ) ) {
MESSAGE(( "Failed value too large test!\n" ));
            goto INTERPOLATE;
        }
        if ( fabs(dDcur/dD0) > 0.9 ) {
MESSAGE(( "Failed slope too steep test!\n" ));
MESSAGE(( "dDcur = %lf dD0 = %lf\n", dDcur, dD0 ));
            goto INTERPOLATE;
        }

MESSAGE(( "Passed steplength test\n" ));

                /* If the criteria has been satisfied then check */
                /* if the true line minimum has been found */

        if ( fabs(dDcur/dD0) > EPS ) goto INTERPOLATE;
        else                         goto DONE;

                /* Use CUBIC interpolation on the current values */
                /* and the previous values */
INTERPOLATE:
                /* CUBIC interpolation */

MESSAGE(( "Performing CUBIC interpolation\n" ));
#if 0
        u1 = dDpre + dDcur - 3*( dFpre - dFcur )/( dApre - dAcur );
        u2 = sqrt( u1*u1 - dDpre*dDcur );
        dTest = dAcur - ( dAcur - dApre ) * ( dDcur + u2 - u1 ) /
                         ( dDcur - dDpre + 2*u2 );
#endif
        dTest = dAcur - dDcur * ( dApre - dAcur ) / ( dDpre - dDcur );

MESSAGE(( "    dTest= %lf\n", dTest ));
                /* Test if the minimum has been bracketed */
                /* If the product of the current and previous derivatives */
                /* is < 0 then it has */

        if (  ( dDcur / dDpre ) < 0.0 ) {

MESSAGE(( "The minimum is bracketed\n" ));

                /* Test if the point is sufficiently within the interval */

            if ( (dTest < ( 1.01 * MIN( dAcur, dApre ))) ||
                 (dTest > ( 0.99 * MAX( dAcur, dApre ))) ) {
MESSAGE(( "Moving to halfway into the bracket\n" ));
                dTest = ( dAcur + dApre ) / 2.0;
            }

            goto CONTINUE;

        }
                /* Test if both points are greator than the minimum */
                /* and the trial point is sufficiently smaller than */
                /* either */

        if ( (dDcur > 0.0) &&
             (0.0 < dTest) &&
             ( dTest < ( 0.9999 * MIN(dAcur, dApre) ) ) ) {
MESSAGE(( "Both points are greator than the minimum\n" ));
                    goto CONTINUE;
        }

                /* Test if both points are less than the minimum and */
                /* the trial point is sufficiently large */

        if ( (dDcur < 0.0) &&
             ( dTest > ( 1.0001 * MAX(dApre,dAcur) ) ) ) {
MESSAGE(( "Both points are less than the minimum and dTest is large enough\n" ));
            goto CONTINUE;
        }

                /* If the trial point is too small , double */
                /* the largest prior point */

        if ( dDcur < 0.0 ) {
MESSAGE(( "Doubling the test point\n" ));
            dTest = 2.0 * MAX(dApre,dAcur);
        }

                /* If the trial point is too large, halve */
                /* the smallest prior point */

        if ( dDcur > 0.0 ) {
MESSAGE(( "Halving the test point\n" ));
            dTest = MIN(dApre,dAcur) / 2.0;
        }

        goto CONTINUE;

                /* Set dApre = dAcur, dAcur = dTest and continue */
CONTINUE:
MESSAGE(( "Continuing!\n" ));
            dApre = dAcur;
            dFpre = dFcur;
            dDpre = dDcur;
            dAcur = dTest;

    }

DONE:
    *dPF0 = dFcur;
    *dPD0 = dDcur;
    return(dAcur);

#if 0
    *dPF0 = dFmin;
    *dPD0 = dDmin;

        /* If the minimum value isn't the current one then copy the minimum */
        /* derivative */

    if ( dAmin != dAcur ) {
        NVectorCopy( nvNewPos, nvTempPos );
        NVectorCopy( nvNewDeriv, nvTempDeriv );
    }

MESSAGE(( "^^^^^^ Returning dA= %lf   dF= %lf   dD= %lf\n",
                dAmin, dFmin, dDmin ));

    return(dAmin);
#endif
}






/*
 *      ConjugateGradient
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Minimize the energy using ConjugateGradients
 *      The starting and ending postion are stored in nvPos.
 */
static void
ConjugateGradient( MINIMIZER mMinimizer, NVECTOR nvPos )
{
#define MINSLOPE        0.000001
#define MINCHANGE       0.01

NVECTOR         nvDir, nvNewPos, nvNewDeriv, nvTempPos, nvTempDeriv;
double          dX, dLastEnergy, dBeta, dDcur, dDdotpre, dDdotcur, dEnergy;
int             iK, iRestartSteps, iCount;

    iCount = 0;

                /* Calculate how many conjugate gradient steps can be */
                /* taken before a restart must be done */

    iRestartSteps = iNVectorSize(nvPos);
    nvDir      = nvNVectorCreate( iRestartSteps );
    nvNewPos   = nvNVectorCreate( iRestartSteps );
    nvNewDeriv = nvNVectorCreate( iRestartSteps );
    nvTempPos  = nvNVectorCreate( iRestartSteps );
    nvTempDeriv = nvNVectorCreate( iRestartSteps );

                /* Conjugate gradient minimization, evaluate the function */
                /* at the current position and use the gradient as the */
                /* direction along which to search. */
                /* keep the square of the length of the derivative in */
                /* dDdotcur */

    dEnergy = dTotalEnergy( mMinimizer, nvPos, nvNewDeriv );
    dDdotcur = dNVectorDot( nvNewDeriv, nvNewDeriv );
    do {

                /* Obtain the next direction vector, use the gradient */
                /* to restart the conjugate gradients method */
                /* nvDir/dDcur is a unit vector in the search direction */

        NVectorCopy( nvDir, nvNewDeriv );

                /* iK keeps track of number of steps until */
                /* conjugate gradients has to be restarted */

        for ( iK=0; iK<iRestartSteps; iK++ ) {

                /* Find the length of the current direction vector */

            dDcur = -dNVectorLen(nvDir);
	    MESSAGE(( "The starting derivative is %lf\n", dDcur ));

                /* Evaluate the energy and derivative at nvPos */

            dLastEnergy = dEnergy;
            dX = dLineSearch( mMinimizer, &dEnergy, &dDcur, 
                                nvPos, nvDir, 
                                nvNewPos, nvNewDeriv, 
                                nvTempPos, nvTempDeriv );

                        /* Test to see if done */
        
	    MESSAGE( ("---LineSearch step=%lf\n", dX ) );

                /* If the line search did not work then stop */

            if ( fabs(dX) < VERYSMALL ) {
	        MESSAGE(( "ERROR, the linesearch produced nothing\n" ));
                goto DONE;
            } else {
	    	NVectorCopy( nvPos, nvNewPos );
	    }


                /* Find the next search direction */        

            dDdotpre = dDdotcur;
            dDdotcur = dNVectorDot( nvNewDeriv, nvNewDeriv );

            if ( dDdotcur/(double)iRestartSteps < 
                        mMinimizer->dMinRmsGradientSquared ) {
VP1(( "Search terminated when RMS gradient became small enough!\n" ));
                goto DONE;
            }


	    if ( mMinimizer->bFCallback != NULL ) {
		mMinimizer->dEnergy = dEnergy;
		mMinimizer->dRmsGradient =
			sqrt(dDdotcur/(double)iRestartSteps);
		SetAtomCoordinates( mMinimizer, nvPos );
MESSAGE(( "Displaying UNIT\n" ));
		if ( !mMinimizer->bFCallback(mMinimizer->PCallbackData) ) {
		    VP0(( "Interrupted!\n" ));
		    goto DONE;
		}
	    }



                /* Make the direction vector a unit vector */

            if ( iCount % 10 == 0 ) 
                VP1(( "Gradient RMS=%lf\n", dDdotcur/(double)iRestartSteps ));

                /* If we dont have to restart this step calculate the new */
                /* conjugate gradient direction */

            if ( iK != ( iRestartSteps - 1 ) ) {
                dBeta = dDdotcur/dDdotpre;
                NVectorTimesScalar( nvDir, nvDir, dBeta );
                NVectorAdd( nvDir, nvDir, nvNewDeriv );
            }
            iCount++;
            if ( iCount >= MAXCONJUGATEGRADIENTSTEPS ) goto DONE;
        } 
    } while ( iCount < MAXCONJUGATEGRADIENTSTEPS );

DONE:
    if ( iCount >= MAXCONJUGATEGRADIENTSTEPS ) {
        VP1(( "Search terminated when max number of steps were taken!\n" ));
    }

    VP1(( "Number of conjugate gradient steps taken: %d\n", iCount ));

    NVectorDestroy( &nvDir );
    NVectorDestroy( &nvNewPos );
    NVectorDestroy( &nvNewDeriv );
    NVectorDestroy( &nvTempPos );
    NVectorDestroy( &nvTempDeriv );
}






/*
---------------------------------------------------------------------------
---------------------------------------------------------------------------
---------------------------------------------------------------------------

*/


#if 0
/*
 *      SteepestDescent
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Perform steepest descent minimization of the minimizer.
 *      The starting and ending postion are stored in nvPos.
 */
static void
SteepestDescent( MINIMIZER mMinimizer, NVECTOR nvPos )
{
#define MINSLOPE        0.000001
#define MINCHANGE       0.01

NVECTOR         nvDir, nvNewPos, nvNewDeriv, nvTempPos, nvTempDeriv;
double          dEnergy, dDeriv, dX, dLastEnergy, dRms;
int             iCount;

                /* Create the Direction vector */

    nvDir       = nvNVectorCreate( iNVectorSize(nvPos) );
    nvNewPos    = nvNVectorCreate( iNVectorSize(nvPos) );
    nvNewDeriv  = nvNVectorCreate( iNVectorSize(nvPos) );
    nvTempPos    = nvNVectorCreate( iNVectorSize(nvPos) );
    nvTempDeriv  = nvNVectorCreate( iNVectorSize(nvPos) );

    iCount = 0; 
                /* Steepest descent minimization, evaluate the function */
                /* at the current position and use the gradient as the */
                /* direction along which to search. */

    dEnergy = dTotalEnergy( mMinimizer, nvPos, nvNewDeriv );
    do {

                /* If the derivative is very small then stop */

        VP1(( "Energy = %12.5lf\n", dEnergy ));
        VP1(( "   Gradient RMS^2= %12.5lf    Goal: %12.5lf\n",
                dRms, mMinimizer->dMinRmsGradientSquared ));
        dDeriv = -dNVectorLen( nvNewDeriv );
        dRms = dDeriv*dDeriv/iNVectorSize(nvPos);
        if ( dRms < mMinimizer->dMinRmsGradientSquared ) {
            VP1(( "RMS gradient is small enough!\n" ));
            goto DONE;
        }

	if ( mMinimizer->bFCallback != NULL ) {
	    mMinimizer->dEnergy = dEnergy;
	    mMinimizer->dRmsGradient = sqrt(dRms);
	    SetAtomCoordinates( mMinimizer, nvPos );
MESSAGE(( "Displaying UNIT\n" ));
	    if ( !mMinimizer->bFCallback(mMinimizer->PCallbackData) ) {
		VP0(( "Interrupted!\n" ));
		goto DONE;
	    }
	}

                /* Make the direction vector a unit vector */

MESSAGE(( "The starting derivative is %lf\n", dDeriv ));
        dLastEnergy = dEnergy;
        NVectorCopy( nvDir, nvNewDeriv );
        dX = dLineSearch( mMinimizer, &dEnergy, &dDeriv, nvPos, nvDir, 
                                nvNewPos, nvNewDeriv, 
                                nvTempPos, nvTempDeriv );
        
        MESSAGE( ("---LineSearch step=%lf\n", dX ) );

                /* If the line search did not work then stop */

        if ( fabs(dX) < VERYSMALL ) {
            VP1(( "Linesearch step was too small\n" ));
            break;
        }

                /* Move to the new point */

        NVectorCopy( nvPos, nvNewPos );

        iCount++; 

    } while ( iCount < MAXSTEEPESTDESCENTSTEPS );

    VP1(( "Exceeded maximum number of steps\n" ));

DONE:
    NVectorDestroy( &nvDir );
    NVectorDestroy( &nvNewPos );
    NVectorDestroy( &nvNewDeriv );
    NVectorDestroy( &nvTempPos );
    NVectorDestroy( &nvTempDeriv );
}
#endif
                




/*
==========================================================================

        Public routines
        
*/

/*
 *      mMinimizerCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a new minimizer minimizer.
 *      Create the containers for atoms, bonds, angles, and torsions.
 */
MINIMIZER
mMinimizerCreate()
{
MINIMIZER       eNew;

    MALLOC( eNew, MINIMIZER, sizeof(MINIMIZERt) );
    
    eNew->vaAtoms      = vaVarArrayCreate( sizeof(EATOMt) );
    eNew->vaBonds      = vaVarArrayCreate( sizeof(EBONDt) );
    eNew->vaAngles     = vaVarArrayCreate( sizeof(EANGLEt) );
    eNew->vaTorsions   = vaVarArrayCreate( sizeof(ETORSIONt) );
    eNew->bFCallback   = NULL;
    eNew->dEnergy      = 0.0;
    eNew->dRmsGradient = 0.0;
    eNew->bMinimizing  = FALSE;		/* Flag to figure out if the */
					/* minimizer is minimizing */
    
    eNew->dMinRmsGradientSquared = MINRMSSQUARED;

    return(eNew);
}



/*
 *      MinimizerDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the minimizer minimizer and release the memory it uses.
 *      After this, the MINIMIZER object will be undefined.
 */
void
MinimizerDestroy( MINIMIZER *mPMinimizer )
{
                /* Destroy the VARARRAYs that contain the atoms, bonds, etc */

    VarArrayDestroy( &((*mPMinimizer)->vaAtoms) );
    VarArrayDestroy( &((*mPMinimizer)->vaBonds) );
    VarArrayDestroy( &((*mPMinimizer)->vaAngles) );
    VarArrayDestroy( &((*mPMinimizer)->vaTorsions) );
    
    FREE( *mPMinimizer );
    *mPMinimizer = NULL;
}





/*
 *      MinimizerAddAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add an atom to the minimizer's list.
 */
void
MinimizerAddAtom( MINIMIZER mMinimizer, ATOM aAtom )
{
EATOMt          eaAtom;
int             iAtom;

                /* First check if the atom has already been added */ 

    iAtom = iFindAtom( mMinimizer, aAtom );
    if ( iAtom != NOTFOUND ) 
	return;

                /* Add the atom */

    MESSAGE(( "Adding atom to MINIMIZER: %s\n", sContainerName(aAtom) ));
    
    eaAtom.aAtom = aAtom;
    VarArrayAdd( mMinimizer->vaAtoms, (GENP)&eaAtom );
}

 

/*
 *      MinimizerAddBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add bond interaction to the minimizer minimizer.
 *      Return FALSE if one of the atoms was not defined within the
 *      MINIMIZER object.
 */
BOOL
bMinimizerAddBond( MINIMIZER mMinimizer, 
		ATOM aAtom1, ATOM aAtom2, 
		double dKb, double dR0)
{
EBONDt          ebBond;
int             iAtom;

                /* Look up the indices of the atoms in the MINIMIZER object */

    iAtom = iFindAtom( mMinimizer, aAtom1 );
    if ( iAtom == NOTFOUND ) return(FALSE);
    ebBond.iAtom1 = iAtom;

    iAtom = iFindAtom( mMinimizer, aAtom2 );
    if ( iAtom == NOTFOUND ) return(FALSE);
    ebBond.iAtom2 = iAtom;

    ebBond.dKb = dKb;
    ebBond.dR0 = dR0;
    
    MESSAGE(( "Adding bond to MINIMIZER: %s - %s  Kb=%lf  R0=%lf\n",
                sContainerName(aAtom1), sContainerName(aAtom2),
                        dKb, dR0 ));
    
        /* Add the bond */
    VarArrayAdd( (mMinimizer->vaBonds), (GENP)&ebBond );
    
    return(TRUE);
}

 
 

/*
 *      MinimizerAddAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add angle interaction to the minimizer minimizer.
 *      Return FALSE if one of the atoms was not defined within the
 *      MINIMIZER object.
 */
BOOL
bMinimizerAddAngle( MINIMIZER mMinimizer, ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, 
		double dKt, double dT0 )
{
EANGLEt         eaAngle;
int             iAtom;

                /* Look up the indices of the atoms in the MINIMIZER object */

    iAtom = iFindAtom( mMinimizer, aAtom1 );
    if ( iAtom == NOTFOUND ) 
	return(FALSE);
    eaAngle.iAtom1 = iAtom;
    iAtom = iFindAtom( mMinimizer, aAtom2 );
    if ( iAtom == NOTFOUND ) 
	return(FALSE);
    eaAngle.iAtom2 = iAtom;
    iAtom = iFindAtom( mMinimizer, aAtom3 );
    if ( iAtom == NOTFOUND ) 
	return(FALSE);
    eaAngle.iAtom3 = iAtom;

    eaAngle.dKt = dKt;
    eaAngle.dT0 = dT0;
    
    MESSAGE(( "Adding angle to MINIMIZER: %s - %s - %s Kt=%lf  T0=%lf\n",
                sContainerName(aAtom1), sContainerName(aAtom2),
                sContainerName(aAtom3), dKt, dT0/DEGTORAD ));
    
        /* Add the angle */
    VarArrayAdd( (mMinimizer->vaAngles), (GENP)&eaAngle );

    return(TRUE);
}

 
 

/*
 *      MinimizerAddTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add torsion interaction to the minimizer minimizer.
 *      Return FALSE if one of the atoms was not defined within the
 *      MINIMIZER object.
 */
BOOL
bMinimizerAddTorsion( MINIMIZER mMinimizer, 
		ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, ATOM aAtom4,
		double dN, double dKp, double dP0 )
{
ETORSIONt       etTorsion;
int             iAtom;

                /* Look up the indices of the atoms in the MINIMIZER object */

    iAtom = iFindAtom( mMinimizer, aAtom1 );
    if ( iAtom == NOTFOUND ) return(FALSE);
    etTorsion.iAtom1 = iAtom;
    iAtom = iFindAtom( mMinimizer, aAtom2 );
    if ( iAtom == NOTFOUND ) return(FALSE);
    etTorsion.iAtom2 = iAtom;
    iAtom = iFindAtom( mMinimizer, aAtom3 );
    if ( iAtom == NOTFOUND ) return(FALSE);
    etTorsion.iAtom3 = iAtom;
    iAtom = iFindAtom( mMinimizer, aAtom4 );
    if ( iAtom == NOTFOUND ) return(FALSE);
    etTorsion.iAtom4 = iAtom;

    etTorsion.dKp = dKp;
    etTorsion.dN  = dN;

			/* If the phase angle is 180.0 degrees then */
			/* flip the sign of the force constant */

    if ( iDoubleCompare( dP0, 0.0 ) != 0 ) {
	etTorsion.dKp = -etTorsion.dKp;
    }

    
    MESSAGE(( "Adding torsion to MINIMIZER: %s - %s - %s - %s Kp=%lf N=%lf\n",
                sContainerName(aAtom1), sContainerName(aAtom2),
                sContainerName(aAtom3), sContainerName(aAtom4), dKp, dN ));
    
        /* Add the torsion */
    VarArrayAdd( (mMinimizer->vaTorsions), (GENP)&etTorsion );

    return(TRUE);
}


/*
 *      MinimizerMinimize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Minimize the minimizer.  Then write the new coordinates directly into
 *      the atoms.
 */
void
MinimizerMinimize( MINIMIZER mMinimizer )
{
NVECTOR nvPosition;

		/* Set the flag saying that the MINIMIZER is active */


    mMinimizer->bMinimizing = TRUE;

    nvPosition = nvGetPositionVector(mMinimizer);

    ConjugateGradient( mMinimizer, nvPosition );   

/*
    SteepestDescent( mMinimizer, nvPosition );
*/
    SetAtomCoordinates( mMinimizer, nvPosition );
    NVectorDestroy( &nvPosition );

		/* Reset the flag saying that the MINIMIZER is passive */

    mMinimizer->bMinimizing = FALSE;


}

