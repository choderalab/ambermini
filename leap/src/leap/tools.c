/*
 *      File:   tools.c
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
 *              This file contains several routines that
 *              perform different functions.
 *              These functions are currently called from
 *              the command line interpretor.
 */





#include        <float.h>
#include        "basics.h"

#include        "vector.h"
#include        "matrix.h"

#include        "classes.h"

#include        "dictionary.h"
#include        "database.h"

#include        "parser.h"

#include        "leap.h"

#include        "tools.h"

#include        "mathop.h"

#include        "sort.h"
/*
 *      dToolAtomR
 *
 *      Author: Bill Ross (1994)
 *
 *      Get atom size, negative if unknown. 
 *
 *      Note that the R* value in the force field is only a 1st 
 *      approximation of size, since actual atom-atom minimum 
 *      energy distances depend on the atom potential energy 
 *      depth (epsilon) as well.
 */

static double
dToolAtomR( ATOM aAtom )
{
PARMSET         psTemp;
int             iTag;
double          dMass, dPolar, dE, dR, dE14, dR14;
double		dScreenF;
int             iElement, iHybrid;
STRING          sType, sDesc;

        iTag = PARM_NOT_FOUND;
        PARMLIB_DEFAULT_LOOP( psTemp,
                ( iTag = iParmSetFindAtom( psTemp, sAtomType(aAtom) ) ));
        if ( iTag == PARM_NOT_FOUND )
                return( -1.0 );
        ParmSetAtom( psTemp, iTag, sType,
                   &dMass, &dPolar, &dE, &dR, &dE14, &dR14, &dScreenF,
                                   &iElement, &iHybrid, sDesc );
        return( dR );
}

/*
 *      dToolFindBiggestAtom
 *
 *      Author: Bill Ross (1994)
 *
 *      Find the biggest atom in the UNIT.
 */
static double
dToolFindBiggestAtom( UNIT uUnit )
{
LOOP            lTemp;
ATOM            aAtom;
double          dMaxR, dR;

    dMaxR = 0.0;

    lTemp = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
        dR = dToolAtomR( aAtom );
        if ( dR < 0.1 ) {
            /* R unknown or 0 */
            if ( iAtomElement( aAtom ) == HYDROGEN )
                dR = 1.0;
            else
                dR = ATOM_DEFAULT_RADIUS;
        }
        if ( dR > dMaxR )
                dMaxR = dR;
    }
    return( dMaxR );
}

/*
 *      dToolUnitMass
 *
 *      Author: Bill Ross (1996)
 *
 *      Get unit's mass, negative if any part unknown. 
 */

static double
dToolUnitMass( UNIT uUnit )
{
LOOP            lAtoms;
ATOM            aAtom;
PARMSET         psTemp;
int             iTag, iElement, iHybrid, iIncomplete;
double          dMass, dPolar, dE, dR, dE14, dR14, dTotalMass;
double		dScreenF;
STRING          sType, sDesc;

        dTotalMass = 0.0;
        iIncomplete = 0;
        lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
        while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
                /*
                 *  try ff=parmlibs 1st
                 */
                iTag = PARM_NOT_FOUND;
                PARMLIB_DEFAULT_LOOP( psTemp,
                    ( iTag = iParmSetFindAtom( psTemp, sAtomType(aAtom) ) ));
                if ( iTag != PARM_NOT_FOUND ) {
                        ParmSetAtom( psTemp, iTag, sType,
                    &dMass, &dPolar, &dE, &dR, &dE14, &dR14, &dScreenF,
                                        &iElement, &iHybrid, sDesc );
                        dTotalMass += dMass;
                        continue;
                }
                /*
                 *  TODO try element
                iElement = iAtomElement( aAtom );
                if ( iElement != NOELEMENT )
                 */
                iIncomplete = 1;
        }
        if ( iIncomplete )
                dTotalMass *= -1;
        return( dTotalMass );
}

/*
 * --------------------------------------------------------
 *
 *        Solvation routines
 *        
 */


#define DEFAULTR        1.5

                /* This parameter controls how close solvent can */
                /* get to solute before they are considered to be */
                /* overlapping.  */
                /* TODO: Play with the CLOSENESSMODIFIER */
                /* TODO: to get good solvent/solute boundary */

#define CLOSENESSMODIFIER       1.0

                /* Different criteria to use in determining whether a */
                /* solvent RESIDUE is discarded */

#define TOOLSOLUTECOLLISION     0x00000001
#define TOOLOUTSIDESHELL        0x00000002
#define TOOLOUTSIDESPHERE       0x00000004
#define TOOLOUTSIDEOFBOX        0x00000008
#define TOOLOUTSIDEOFOCTBOX     0x00000010

typedef struct  SPHEREtS {
        VECTOR                  vPos;
        double                  dR;
        struct SPHEREtS         *sPNextInteresting;
} SPHEREt;


typedef struct  {
        VECTOR          vCenter;
        double          dRadiusSqd;
        double          dX;
        double          dY;
        double          dZ;
} CRITERIAt;


void
ToolSanityCheckBox( UNIT uUnit )
{
        double  dX, dY, dZ;
        /*
         *  do a minimal sanity check just to be sure..
         */
        UnitGetBox( uUnit, &dX, &dY, &dZ);
        if ( dX < 0.1  ||  dY < 0.1  ||  dZ < 0.1  ||  
                        dUnitBeta( uUnit ) < 0.05 ) {
                /*
                 *  failed sanity check - clean it up 
                 */
                if ( bUnitUseBox( uUnit ) ) {
                        VP0(( " (turning off box flag on %s - bad param(s)\n", 
                                sContainerName((CONTAINER) uUnit ) ));
                        VP0(( "   beta %e  XYZ %e %e %e)\n",
                                dUnitBeta( uUnit ), dX, dY, dZ ));
                }
                UnitSetBox( uUnit, 0.0, 0.0, 0.0 );
                UnitSetBeta( uUnit, 0.0 );
                UnitSetUseBox( uUnit, FALSE );
        }
}

static void
zToolSetTempRadii( UNIT uUnit )
{
LOOP    lTemp;
ATOM    aAtom;
double  dR;

    lTemp = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
        dR = dToolAtomR( aAtom );
        if ( dR < 0.1 ) {
                /* R unknown or 0 */
                if ( iAtomElement( aAtom ) == HYDROGEN )
                        dR = 1.0;
                else
                        dR = ATOM_DEFAULT_RADIUS;
        }
        aAtom->dTemp = dR;
    }
}


/*
 *      ToolCenterUnitByRadii
 *
 *      Author: Bill Ross (1994)
 *
 *      Center unit according to a bounding box enclosing 
 *      atoms' vdw radii. Put the radii in each atom's dTemp 
 *      field. If requested, align unit on principal axes
 *      before centering.
 */
void
ToolCenterUnitByRadii( UNIT uUnit, BOOL bOrient )
{
LOOP            lTemp;
ATOM            aAtom;
double          dR, dXmin, dYmin, dZmin, dXmax, dYmax, dZmax;
double          dX, dY, dZ;
VECTOR          vPos;
int             iFirst = 1;

    ToolSanityCheckBox( uUnit );

    /*
     *  if unit already has a box, it is assumed to be properly
     *  set up already (e.g. orienting WATBOX216 here will
     *  put its diagonals on the axes, and forcing the box
     *  to be at the furthest vdw boundaries is not appropriate
     *  for such a solvent because it is equilibrated in a
     *  periodic system so that at the boundary, the vdw
     *  spheres mesh)
     */
    if ( bUnitUseBox( uUnit ) )
        return;

    if ( bOrient == TRUE )
        ToolOrientPrincipleAxisAlongCoordinateAxis( uUnit );

    lTemp = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
        dR = dToolAtomR( aAtom );
        if ( dR < 0.1 ) {
            /* R unknown or 0 */
            if ( iAtomElement( aAtom ) == HYDROGEN )
                dR = 1.0;
            else {
                dR = ATOM_DEFAULT_RADIUS;
                VP0(( " (using default radius %f for %s)\n", 
                                ATOM_DEFAULT_RADIUS, sAtomName(aAtom) ));
            }
        }
        aAtom->dTemp = dR;
        vPos = vAtomPosition(aAtom);

        dX = dVX(&vPos) + dR;
        dY = dVY(&vPos) + dR;
        dZ = dVZ(&vPos) + dR;

        if ( iFirst ) {
                dXmax = dX;
                dYmax = dY;
                dZmax = dZ;
        } else {
                if ( dX > dXmax ) 
                        dXmax = dX;
                if ( dY > dYmax ) 
                        dYmax = dY;
                if ( dZ > dZmax ) 
                        dZmax = dZ;
        }

        dX = dVX(&vPos) - dR;
        dY = dVY(&vPos) - dR;
        dZ = dVZ(&vPos) - dR;

        if ( iFirst ) {
                dXmin = dX;
                dYmin = dY;
                dZmin = dZ;
                iFirst = 0;
        } else {
                if ( dX < dXmin )
                        dXmin = dX;
                if ( dY < dYmin )
                        dYmin = dY;
                if ( dZ < dZmin )
                        dZmin = dZ;
        }
    }

    /*
     *  define center of bounding box
     */
    dX = dXmin + 0.5 * (dXmax - dXmin);
    dY = dYmin + 0.5 * (dYmax - dYmin);
    dZ = dZmin + 0.5 * (dZmax - dZmin);

    /*
     *  translate center to origin
     */
    VectorDef( &vPos, -dX, -dY, -dZ );
    ContainerTranslateBy((CONTAINER) uUnit, vPos );
    
    UnitSetBox( uUnit, dXmax-dXmin, dYmax-dYmin, dZmax-dZmin );
}

void
ToolSetUnitBoxByCenters( UNIT uUnit )
{
LOOP            lTemp;
ATOM            aAtom;
double          dXmin, dYmin, dZmin, dXmax, dYmax, dZmax;
double          dX, dY, dZ;
VECTOR          vPos;
int             iFirst = 1;

    lTemp = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
        vPos = vAtomPosition(aAtom);
        dX = dVX(&vPos);
        dY = dVY(&vPos);
        dZ = dVZ(&vPos);


        if ( iFirst ) {
                dXmax = dXmin = dX;
                dYmax = dYmin = dY;
                dZmax = dZmin = dZ;
                iFirst = 0;
                continue;
        }
        if ( dX > dXmax ) 
                dXmax = dX;
        else if ( dX < dXmin )
                dXmin = dX;

        if ( dY > dYmax ) 
                dYmax = dY;
        else if ( dY < dYmin )
                dYmin = dY;

        if ( dZ > dZmax ) 
                dZmax = dZ;
        else if ( dZ < dZmin )
                dZmin = dZ;
    }

    UnitSetBox( uUnit, dXmax-dXmin, dYmax-dYmin, dZmax-dZmin );
}

/*
 *      zToolSetupSolvent
 *
 *      Author: Bill Ross
 *
 *      Make a disposable copy of unit with all solvent attributes set.
 *
 */
UNIT
zToolSetupSolvent( UNIT uSolvent )
{
LOOP            lResidues;
RESIDUE         rRes;

    uSolvent = (UNIT)oCopy( (OBJEKT)uSolvent );

    if ( !bUnitUseBox( uSolvent ) ) {
        /*
         *  Center the unit in terms of its bounding box for periodic
         *      packing purposes.  
         */
        VP0(("Solvent has no box, so preparing by making box including vdw\n"));
        VP0(("(Use 'setBox centers' first if box was pre-equilibrated)\n"));
        ToolCenterUnitByRadii( uSolvent, FALSE );

        UnitSetUseBox( uSolvent, TRUE );
        UnitSetBeta( uSolvent, 90.0*DEGTORAD );
    }

    /* 
     *  Make sure that all solvent residues are marked as solvent
     *  and stuff the atom radii for solvent screening
     */
    lResidues = lLoop( (OBJEKT)uSolvent, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lResidues)) )
        ResidueSetType( rRes, RESTYPESOLVENT);

    /*
     *  set solvent atoms' dTemp = R*
     */
    zToolSetTempRadii( uSolvent );

    return( uSolvent );
}


/*
 *      zToolBuildSoluteArray
 *
 *      Author: Christian Schafmeister (1991)
 *      Amended: Bill Ross (1994)
 *
 *      Build a VARARRAY that contains the solute atoms
 *      so that they can be scanned over quickly.
 */
static void
zToolBuildSoluteArray( UNIT uSolute, double dCloseness, VARARRAY *vaPSolute )
{
int             i, iAtoms;
LOOP            lTemp;
VARARRAY        vaSolute;
ATOM            aAtom;
SPHEREt         *sPTmp;

                /* Collect all of the atoms in the solute UNIT */
                /* into an array of positions and RADII for faster */
                /* checking of close contacts */


    iAtoms = 0;
    lTemp = lLoop( (OBJEKT)uSolute, ATOMS );
    while ( oNext(&lTemp) != NULL ) 
        iAtoms++;
    
    vaSolute = vaVarArrayCreate( sizeof(SPHEREt) );
    VarArraySetSize( vaSolute, iAtoms );
    *vaPSolute = vaSolute;
 
    if ( iAtoms == 0 ) 
        return;

    sPTmp = PVAI( vaSolute, SPHEREt, 0 );
    lTemp = lLoop( (OBJEKT)uSolute, ATOMS );
    for (i = 0; i < iAtoms; i++, sPTmp++ ) {

        aAtom = (ATOM) oNext(&lTemp);
        sPTmp->dR = dToolAtomR( aAtom );
        if ( sPTmp->dR < 0.1 )  {
            /* R unknown or 0 */
            if ( iAtomElement( aAtom ) == HYDROGEN )
                sPTmp->dR = 1.0;
            else
                sPTmp->dR = ATOM_DEFAULT_RADIUS;
        }
        sPTmp->dR *= dCloseness;
        sPTmp->vPos = vAtomPosition( aAtom );
    }
}


/*
 *      zbToolAtomFailsCriteria
 *
 *      Author: Christian Schafmeister (1991)
 *      Revised: Bill Ross (1994)
 *
 *      Return TRUE if the atom fails any criteria determined 
 *      by iCriteria or if there are any contacts with the
 *      close spheres in the sPFirstInteresting list.
 */
static BOOL
zbToolAtomFailsCriteria( ATOM aAtom, int iCriteria, CRITERIAt *cPCriteria,
                                        SPHEREt *sPFirstInteresting, 
                                        double dCloseness, double dFarness2)
{
double          dR, dDist2, dRadii, dX, dY, dZ, dDist, dClosest;
VECTOR          vPos;
SPHEREt         *sPCurInteresting;

    vPos = vAtomPosition(aAtom);

    /*
     *  first check 'outside' clip criteria since this is faster
     */
    if ( iCriteria & (TOOLOUTSIDEOFBOX | TOOLOUTSIDEOFOCTBOX) ) {
        double  dXabs, dYabs, dZabs;

        dXabs = fabs(dVX(&vPos));
        dYabs = fabs(dVY(&vPos));
        dZabs = fabs(dVZ(&vPos));

        if ( iCriteria & TOOLOUTSIDEOFBOX ) {
                /*
                 *  check if atom falls outside of clipped rectangular box
                 */
                if ( dXabs >= cPCriteria->dX ) return(TRUE);
                if ( dYabs >= cPCriteria->dY ) return(TRUE);
                if ( dZabs >= cPCriteria->dZ ) return(TRUE);
        }
        if ( iCriteria & TOOLOUTSIDEOFOCTBOX ) {
                /*
                 *  check if atom falls outside of oct/diagonal clip 
                 */
                dX = 0.5 * (cPCriteria->dX - dXabs) / cPCriteria->dX;
                dX = fabs( dX - 0.5 );
                dY = 0.5 * (cPCriteria->dY - dYabs) / cPCriteria->dY;
                dY = fabs( dY - 0.5 );
                dZ = 0.5 * (cPCriteria->dZ - dZabs) / cPCriteria->dZ;
                dZ = fabs( dZ - 0.5 );

                if ( ( dX + dY + dZ )  >  0.75 )
                        return(TRUE);
        }
    }

    if ( iCriteria & TOOLOUTSIDESPHERE ) {
        /*
         *  check if atom falls outside of clipped sphere
         */
        dX = dVX(&vPos) - dVX(&(cPCriteria->vCenter));
        dY = dVY(&vPos) - dVY(&(cPCriteria->vCenter));
        dZ = dVZ(&vPos) - dVZ(&(cPCriteria->vCenter));
        dDist2 = dX * dX  +  dY * dY  +  dZ * dZ;
        if ( dDist2 > cPCriteria->dRadiusSqd ) return(TRUE);
    }

    /*
     *  consider overlap w/ solute atoms - if there are no candidates,
     *  there is no problem unless a shell is desired (in which case
     *  the caller is responsible for making sure all solute atoms within
     *  shell range of the solvent box are in the list).
     */
    sPCurInteresting = sPFirstInteresting;
    if ( sPCurInteresting == NULL )
        return((iCriteria&TOOLOUTSIDESHELL)!=0);

    dClosest = NOSHELL;
    dR = aAtom->dTemp;  /* radius */
    if ( dR > 0.0 )
        dR *= dCloseness * CLOSENESSMODIFIER;
    else
        dR = DEFAULTR * dCloseness * CLOSENESSMODIFIER;
 
    do {
                /* Check the distance-squared  between the centers */
        
        dX = dVX(&vPos) - dVX(&(sPCurInteresting->vPos));
        dY = dVY(&vPos) - dVY(&(sPCurInteresting->vPos));
        dZ = dVZ(&vPos) - dVZ(&(sPCurInteresting->vPos));

        dDist = dX*dX + dY*dY + dZ*dZ;

        dRadii = dR + sPCurInteresting->dR;
        dRadii *= dRadii;

        if ( dDist < dClosest ) 
            dClosest = dDist;

        if ( dDist < dRadii && (iCriteria & TOOLSOLUTECOLLISION) )
            return(TRUE);
        
                /* Now move the pointer to the next interesting sphere */

        sPCurInteresting = sPCurInteresting->sPNextInteresting;
    } while ( sPCurInteresting != NULL );

                /* Check if the closest solute atom to the solvent */
                /* is farther than the dFarness parameter.  If it is */
                /* then it will be outside the solvent shell and has */
                /* to be eliminated */

    if ( dClosest>=dFarness2 && (iCriteria & TOOLOUTSIDESHELL) )
        return(TRUE);

    return(FALSE);
}




/*
 *      zToolFindCloseSoluteAtoms
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Find the solute atoms that have a possibility of 
 *      interacting with the atoms within the solvent box
 *      using a simple 'square bound' test. The list is
 *      created by setting pointers in the existing solute
 *      scratch array.
 *
 *      The caller can specify a buffer distance by which
 *      the box should be expanded to make sure that all 
 *      possible interesting spheres are included. This
 *      is necessary for shell setup.
 */
static void
zToolFindCloseSoluteAtoms( VARARRAY vaPSolute, SPHEREt **sPPFirstInteresting,
                double dSoluteMaxR, double dCloseness, VECTOR *vPCenter, 
                double dXWidth, double dYWidth, double dZWidth,
                double dBufferZone )
{
SPHEREt         *sPLastInteresting, *sPCur;
int             i, iMax;
double          dZmax, dZmin, dXmin, dXmax, dYmin, dYmax;
double          dTemp;

    MESSAGE(( "Searching for close solute atoms\n" ));
    *sPPFirstInteresting = NULL;
    sPLastInteresting = NULL;

    iMax = iVarArrayElementCount( vaPSolute );
    if ( !iMax ) {
        MESSAGE(( "No solute atoms..\n" ));
        return;
    }

    /*
     *  Determine clearance from the box for testing whether 
     *  a solute atom might contact an atom in the box. Assumes
     *  solvent box includes vdw.
     */

    dTemp = dSoluteMaxR * 2.5 * dCloseness;

    /*
     *  when building a shell, solute atoms within the
     *  shell distance are also of interest
     */
    dTemp += dBufferZone;

    dXmin = dVX(vPCenter) - dXWidth/2.0 - dTemp;
    dXmax = dVX(vPCenter) + dXWidth/2.0 + dTemp;
    dYmin = dVY(vPCenter) - dYWidth/2.0 - dTemp;
    dYmax = dVY(vPCenter) + dYWidth/2.0 + dTemp;
    dZmin = dVZ(vPCenter) - dZWidth/2.0 - dTemp;
    dZmax = dVZ(vPCenter) + dZWidth/2.0 + dTemp;

    /*
     *  loop over the solute array
     */
    sPCur = PVAI( vaPSolute, SPHEREt, 0 );
    for ( i=0; i<iMax; i++, sPCur++ ) {
        /*    
         *  if atom.coord outside solvent limit, skip atom
         */
        if ( dVX(&(sPCur->vPos)) < dXmin ) continue;
        if ( dVX(&(sPCur->vPos)) > dXmax ) continue;
        if ( dVY(&(sPCur->vPos)) < dYmin ) continue;
        if ( dVY(&(sPCur->vPos)) > dYmax ) continue;
        if ( dVZ(&(sPCur->vPos)) < dZmin ) continue;
        if ( dVZ(&(sPCur->vPos)) > dZmax ) continue;

        /*    
         *  all atom.coords inside solvent limit, so add to list
         */

        MESSAGE(( "Found an interesting sphere\n" ));
        if ( sPLastInteresting == NULL )
                *sPPFirstInteresting = sPCur;
        else
                sPLastInteresting->sPNextInteresting = sPCur;

        sPLastInteresting = sPCur;
        sPCur->sPNextInteresting = NULL;
    }
}





/*
 *      zToolDiscardInvalidSolvent
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Check the solvent UNIT for solvent RESIDUES 
 *      that do not pass the tests specified by the iCriteria
 *      flags.
 *      This means checking for things like if the
 *      ATOMs of the solvent UNIT are colliding with
 *      solute, whether they fall outside of a shell around
 *      the solute, or if they fall outside of a solvent sphere.
 */
static void
zToolDiscardInvalidSolvent( UNIT uSolvent, int iCriteria, CRITERIAt *cPCriteria,
        SPHEREt *sPFirstInteresting, double dCloseness, double dFarness2 )
{
LOOP            lRes, lTemp;
RESIDUE         rRes;
ATOM            aAtom;

                /* Check each atom in each molecule of the */
                /* solvent unit.  If there is a collision with */
                /* an atom in the solute then discard the entire */
                /* molecule of solvent */


    lRes = lLoop( (OBJEKT)uSolvent, RESIDUES );
    if ( sPFirstInteresting == NULL  &&  iCriteria & TOOLOUTSIDESHELL ) {
        /* 
         *  remove all solvent (unlikely)
         */
        while ( ( rRes = (RESIDUE)oNext(&lRes) ) != NULL ) {
            REF( rRes );  /* bContainerRemove() needs this */
            bContainerRemove( (CONTAINER)uSolvent, (OBJEKT)rRes );
            ContainerDestroy((CONTAINER *) &rRes );
        }
        return;
    }
    while ( ( rRes = (RESIDUE)oNext(&lRes) ) != NULL ) {
            /*
             *  check each atom in the solvent residue against the
             *  'interesting' solute prospects
             */
            lTemp = lLoop( (OBJEKT)rRes, ATOMS );
            while ( ( aAtom = (ATOM)oNext(&lTemp)) != NULL ) {
                if ( zbToolAtomFailsCriteria( aAtom, iCriteria, cPCriteria,
                                sPFirstInteresting, 
                                dCloseness, dFarness2 ) ) {

                        /* Remove the RESIDUE from the solvent unit */
                    MESSAGE(( "Removing a residue\n" ));                

                    REF( rRes );  /* bContainerRemove() needs this */
                    bContainerRemove( (CONTAINER)uSolvent, (OBJEKT)rRes );
                    ContainerDestroy((CONTAINER *) &rRes );
                    break;
                }
            }
    }

}




/*
 *      zToolFindBoundingBoxRadii
 *
 *      Author: Bill Ross (1994)
 *
 *      Find the dimensions of the box which contains
 *      the UNIT, including atom r*'s. Assumes box is
 *      centered; this is the max of the absolute values,
 *      so is 1/2 of each edge.
 */
static void
zToolFindBoundingBoxRadii( UNIT uUnit, double *dPX, double *dPY, double *dPZ )
{
LOOP            lTemp;
ATOM            aAtom;
double          dR, dX, dY, dZ;
VECTOR          vPos;

    *dPX = 0.0;
    *dPY = 0.0;
    *dPZ = 0.0;

    lTemp = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
        dR = dToolAtomR( aAtom );
        if ( dR < 0.1 ) {
            /* R unknown or 0 */
            if ( iAtomElement( aAtom ) == HYDROGEN )
                dR = 1.0;
            else
                dR = ATOM_DEFAULT_RADIUS;
        }
        vPos = vAtomPosition(aAtom);
        dX = fabs(dVX(&vPos)) + dR;
        dY = fabs(dVY(&vPos)) + dR;
        dZ = fabs(dVZ(&vPos)) + dR;
        if ( dX > *dPX ) *dPX = dX;
        if ( dY > *dPY ) *dPY = dY;
        if ( dZ > *dPZ ) *dPZ = dZ;
    }
}

/*
 *      zToolAddAllBoxes
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Overlay many solvent boxes on the system.
 */
static void
zToolAddAllBoxes( UNIT uSolute, VARARRAY vaSolute, UNIT uSolvent, 
                double dCloseness,
                int iX, int iY, int iZ,
                double dXStart, double dYStart, double dZStart,
                double dXSolvent, double dYSolvent, double dZSolvent,
                double dBuffer, int iCriteria, CRITERIAt *cPCriteria )
{
double          dX, dY, dZ, dSoluteMaxR;
int             ix, iy, iz;
VECTOR          vPos;
SPHEREt         *sPFirstInteresting;
UNIT            uSolventDup;
double          dFarness2 = dBuffer * dBuffer;
                        

    dSoluteMaxR = dToolFindBiggestAtom( uSolute );

                /* Loop over all solvent boxes to add */

    VP2(( "The number of boxes:  x=%2d  y=%2d  z=%2d\n", iX, iY, iZ ));

    dX = dXStart;
    for ( ix=0; ix<iX; ix++, dX -= dXSolvent ) {
        dY = dYStart;
        for ( iy=0; iy<iY; iy++, dY -= dYSolvent ) {
            dZ = dZStart;
            for ( iz=0; iz<iZ; iz++, dZ -= dZSolvent ) {
                VP2(( "Adding box at: x=%d  y=%d  z=%d\n", ix, iy, iz));
                VectorDef( &vPos, dX, dY, dZ );

                /*
                 *  make a list of solute atoms that might contact 
                 *      solvent or (if dBuffer) be within shell distance
                 */
                zToolFindCloseSoluteAtoms( vaSolute, &sPFirstInteresting,
                                        dSoluteMaxR, dCloseness, &vPos,
                                        dXSolvent, dYSolvent, dZSolvent, 
                                        dBuffer );
                /*
                 *  copy the solvent box & place it
                 */
                uSolventDup = (UNIT)oCopy( (OBJEKT)uSolvent );
                ContainerCenterAt((CONTAINER) uSolventDup, vPos );
                VP1(( "Center of solvent box is: %lf, %lf, %lf\n",
                                dX, dY, dZ ));

                /*
                 *  discard unwanted residues from solvent box
                 *      TODO - for multi-residue solute molecules,
                 *      probably need to improve
                 */
                zToolDiscardInvalidSolvent( uSolventDup, iCriteria, cPCriteria,
                                        sPFirstInteresting, 
                                        dCloseness, dFarness2 );
                UnitJoin( uSolute, uSolventDup );

            }
        }
    }
}

static void
EwaldRotate( UNIT uUnit, double *dPAngle )
{
LOOP    lAtoms;
ATOM    aAtom;
double  tetra_angl, pi, phi, cos1, sin1, cos2, sin2;
double  t11, t12, t13, t21, t22, t23, t31, t32, t33;

/* for isotropic octbox only */

  tetra_angl=2*acos(1./sqrt(3.));
  pi=3.1415927;
  phi=pi/4.;
  cos1=cos(phi);
  sin1=sin(phi);
  phi=pi/2.-tetra_angl/2.;
  cos2=sqrt(2.)/sqrt(3.);
  sin2=1./sqrt(3.); 

  /********************************************************/
  /*       45 around z axis, (90-tetra/2) around y axis, */
  /*       90 around x axis                               */
  /*                                                      */
  /*   (1  0  0)    (cos2  0 -sin2)    (cos1 -sin1  0)    */
  /*   (0  0 -1)    (   0  1     0)    (sin1  cos1  0)    */
  /*   (0  1  0)    (sin2  0  cos2)    (   0     0  1)    */
  /*                                                      */
  /*   cntr-clk       clock              clock            */
  /*   Looking down + axis of rotation toward origin      */
  /********************************************************/

  t11= cos2*cos1;
  t12=-cos2*sin1;
  t13=-sin2;
  t21=-sin2*cos1;
  t22= sin2*sin1;
  t23=-cos2;
  t31= sin1;
  t32= cos1;
  t33=0;

  lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
  while ( (aAtom = (ATOM)oNext(&lAtoms)) != NULL ) {
        double  dX, dY, dZ;

        dX = vAtomPosition(aAtom).dX;
        dY = vAtomPosition(aAtom).dY;
        dZ = vAtomPosition(aAtom).dZ;

        vAtomPosition(aAtom).dX = t11*dX + t12*dY + t13*dZ;
        vAtomPosition(aAtom).dY = t21*dX + t22*dY + t23*dZ;
        vAtomPosition(aAtom).dZ = t31*dX + t32*dY + t33*dZ;
  }
  *dPAngle = tetra_angl*180./pi;
}

/*
 *      zToolSolvateAndShell
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This routine solvates the ENTIRE solute UNIT.
 *      It does this by determining the minimum number
 *      of solvent UNITS that will solvate the solute with
 *      a caller specified buffer zone on all sides.
 *      This routine makes repeated copies of the solvent UNIT.
 *      This routine also creates solvent SHELLS.  This is done
 *      by specifying a reasonable dFarness parameter which defines
 *      the maximum distance a solvent atom may be from the nearest
 *      solute atom.
 *      Define the bounding box of the solute/solvent system only if
 *      !bShell.
 *      If bClip is TRUE then solvent that falls out of the box centered on
 *      0,0,0 and with a width of dXW, dYW, dZW in the three
 *      coordinate directions will be thrown out.
 */
void
zToolSolvateAndShell( UNIT uSolute, UNIT uSolvent, 
                double dXW, double dYW, double dZW, double dCloseness,
                double dFarness, BOOL bShell, BOOL bClip, BOOL bOct,
                BOOL bIsotropic )
{
double          dXBox, dYBox, dZBox;
double          dXWidth, dYWidth, dZWidth;
int             iX, iY, iZ;
double          dXSolvent, dYSolvent, dZSolvent;
double          dXStart, dYStart, dZStart;
VARARRAY        vaSolute;
double          dBuffer;
int             iCriteria;
CRITERIAt       cCriteria;


    /*
     *  make temporary array of solute atoms for faster checking
     */
    zToolBuildSoluteArray( uSolute, dCloseness, &vaSolute );

                /* Setup the flags that control what criteria are used */
                /* to reject solvent molecules */

    iCriteria = TOOLSOLUTECOLLISION;
    if ( bShell ) 
        iCriteria |= TOOLOUTSIDESHELL;

    /*
     *  get bounding box of solute and set bound on solvated system
     *  - box set already by ToolCenterUnitByRadii()
     */
    UnitGetBox( uSolute, &dXBox, &dYBox, &dZBox );
    VP0(( "  Solute vdw bounding box:              %-5.3lf %-5.3lf %-5.3lf\n", 
                                dXBox, dYBox, dZBox ));

    dXWidth = dXBox + dXW * 2;
    dYWidth = dYBox + dYW * 2;
    dZWidth = dZBox + dZW * 2;

    if ( bIsotropic ) {
        double  dTemp, dMax;

        dTemp = dXWidth * dYWidth * dZWidth;

        dMax = MAX(dXWidth, dYWidth);
        dMax = MAX(dMax, dZWidth);
        dXWidth = dYWidth = dZWidth = dMax;

        dTemp = (dMax * dMax * dMax - dTemp ) / dTemp;

        VP0(( "  Total bounding box for atom centers:  %5.3lf %5.3lf %5.3lf\n", 
                        dXWidth, dYWidth, dZWidth ));
        VP0(( "      (box expansion for 'iso' is %5.1lf%%)\n", dTemp * 100.0 ));

        /*
         *  to make the actual clip right, 'iso' the solute box
         */
        dTemp = MAX(dXBox, dYBox);
        dTemp = MAX(dTemp, dZBox);
        dXBox = dYBox = dZBox = dTemp;
    } else
        VP0(( "  Total bounding box for atom centers:  %5.3lf %5.3lf %5.3lf\n", 
                        dXWidth, dYWidth, dZWidth ));

    if ( bClip ) {
        /*
         *  If the solvated system should be clipped to the exact
         *      size the user specified then note the criterion
         *      & dimensions (for 0,0,0-centered system)
         */
        iCriteria |= TOOLOUTSIDEOFBOX;
        cCriteria.dX = 0.5 * dXBox + dXW;
        cCriteria.dY = 0.5 * dYBox + dYW;
        cCriteria.dZ = 0.5 * dZBox + dZW;
    }
    if ( bOct ) {
        /* maybe allow oct clip on integer boxes someday.. but for now: */
        if ( !bClip )
                DFATAL(( "oct but no clip\n" ));
        iCriteria |= TOOLOUTSIDEOFOCTBOX;
    }
    /*
     *  see how many solvent boxes req'd in each dimension
     */
    UnitGetBox( uSolvent, &dXSolvent, &dYSolvent, &dZSolvent );
    VP0(( "  Solvent unit box:                     %5.3lf %5.3lf %5.3lf\n", 
                                        dXSolvent, dYSolvent, dZSolvent ));

    iX = ( dXWidth / dXSolvent ) + 1;
    iY = ( dYWidth / dYSolvent ) + 1;
    iZ = ( dZWidth / dZSolvent ) + 1;

    /* 
     *  Calculate the center of the first solvent box 
     *  (the one that goes in the max XYZ corner), given
     *  that the solute is centered at 0,0,0
     */
    dXStart = 0.5 * dXSolvent * (double) (iX-1);
    dYStart = 0.5 * dYSolvent * (double) (iY-1);
    dZStart = 0.5 * dZSolvent * (double) (iZ-1);

                /* If the caller wants a solvent shell then */
                /* make sure that the box used to find interesting solute */
                /* spheres takes into account the dFarness parameter */
                /* so that there are at least some solute spheres in */
                /* the interesting list to check against solvent */

    if ( bShell ) 
        dBuffer = dFarness;
    else 
        dBuffer = 0.0;


                /* Add all of the boxes of solvent */


    zToolAddAllBoxes( uSolute, vaSolute, uSolvent, dCloseness, 
                        iX, iY, iZ,
                        dXStart, dYStart, dZStart,
                        dXSolvent, dYSolvent, dZSolvent,
                        dBuffer, iCriteria, &cCriteria );

    VarArrayDestroy( &vaSolute );

                /* Define the size of the new solute/solvent system */

    if ( !bShell ) {
        double  dMaxX, dMaxY, dMaxZ, dVolume, dMass;

        ToolCenterUnitByRadii( uSolute, FALSE );
        UnitSetUseBox( uSolute, TRUE );
        if ( bOct ) {
                UnitSetBoxOct( uSolute, TRUE );
                if ( bIsotropic ) {
                        double  dAngle;

                        zToolFindBoundingBoxRadii( uSolute, &dMaxX, &dMaxY, 
                                                                &dMaxZ );
                        EwaldRotate( uSolute, &dAngle );
                        UnitSetBeta( uSolute, dAngle*DEGTORAD );
                        dAngle = ( dMaxX + dMaxY + dMaxZ ) / 3.0;
                        /* MFC just add an angstrom to the desired box size
                           rather than using the bounding box size */
                        dAngle = cCriteria.dX + .5;
                        dMaxX = dMaxY = dMaxZ = dAngle * sqrt(3.0) * 0.5;
                } else
                        zToolFindBoundingBoxRadii( uSolute, &dMaxX, &dMaxY, 
                                                                &dMaxZ );
        } else {
                zToolFindBoundingBoxRadii( uSolute, &dMaxX, &dMaxY, &dMaxZ );
                UnitSetBeta( uSolute, 90.0*DEGTORAD );
        }
        dMaxX *= 2;
        dMaxY *= 2;
        dMaxZ *= 2;
        UnitSetBox( uSolute, dMaxX, dMaxY, dMaxZ );

        if ( !bIsotropic )
                VP0(( 
                    "  Total vdw box size:%s%5.3lf %5.3lf %5.3lf angstroms.\n", 
                    "                   ", dMaxX, dMaxY, dMaxZ ));
        dVolume = dMaxX * dMaxY * dMaxZ;
        if ( bOct )
                dVolume *= 0.7698004;   /* = sqrt(1 - 3*cost^2 - 2*cost^3),
                                    where cost = -1/3 = cos(109.471)   */

        VP0(( "  Volume: %5.3lf A^3 %s\n", dVolume, (bOct ? "(oct)" : "") ));

        dMass = dToolUnitMass( uSolute );

        if ( dMass > 0.0 ) {
                VP0(( "  Total mass %5.3f amu,  Density %5.3lf g/cc\n", 
                        dMass, dMass / ( dVolume * 0.602204 ) )); 
        } else if ( dMass < 0.0 ) {
                VP0(( "  Mass > %5.3f amu,  Density > %5.3lf g/cc\n",
                        fabs(dMass), fabs(dMass) / ( dVolume * 0.602204 ) ));
                VP0(( "      (type - hence mass - of one or more atoms could not be found)\n" ));
        } else {
                VP0(( "  Mass could not be determined, so density unknown\n" ));
                VP0(( "      (i.e. type of all atoms could not be found)\n" ));
        }
    }
}





/*
 *      zToolSolvateInSphere
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This routine solvates a portion of the UNIT within a 
 *      sphere of solvent.
 *      It does this by determining the minimum number of
 *      solvent boxes that are required to contain the sphere
 *      specified in the arguments.  The routine
 *      makes repeated copies of the solvent UNIT.
 */
void
zToolSolvateInSphere( UNIT uSolute, UNIT uSolvent, VECTOR *vPCenter,
                                double dRadius, double dCloseness )
{
int             iX, iY, iZ;
double          dXSolvent, dYSolvent, dZSolvent;
double          dXStart, dYStart, dZStart;
VARARRAY        vaSolute;
double          dDiam, dBuffer;
int             iCriteria;
CRITERIAt       cCriteria;

    zToolBuildSoluteArray( uSolute, dCloseness, &vaSolute );

    UnitGetBox( uSolvent, &dXSolvent, &dYSolvent, &dZSolvent );

                /* Calculate how many solvent boxes are required */
                /* to enclose the sphere in each direction */

    dDiam = 2.0 * dRadius;
    iX = dDiam / dXSolvent + 1;
    iY = dDiam / dYSolvent + 1;
    iZ = dDiam / dZSolvent + 1;

                /* Calculate the position of the center of */
                /* the first solvent box */

    dXStart = 0.5 * dXSolvent * (double)(iX-1) + dVX(vPCenter);
    dYStart = 0.5 * dYSolvent * (double)(iY-1) + dVY(vPCenter);
    dZStart = 0.5 * dZSolvent * (double)(iZ-1) + dVZ(vPCenter);

    dBuffer = 0.0;
    iCriteria = TOOLSOLUTECOLLISION | TOOLOUTSIDESPHERE;
    cCriteria.dRadiusSqd = dRadius * dRadius;
    cCriteria.vCenter = *vPCenter;
    
                /* Add all of the boxes of solvent */

    zToolAddAllBoxes( uSolute, vaSolute, uSolvent, dCloseness, 
                        iX, iY, iZ,
                        dXStart, dYStart, dZStart,
                        dXSolvent, dYSolvent, dZSolvent,
                        dBuffer, iCriteria, &cCriteria );

    VarArrayDestroy( &vaSolute );
}




/*
        add ions at residue 1st-atom positions
*/

/*
 *  ToolInitSolventPotential() - for each solvent residue,
 *      calc potential at xyz of 1st atom.
 *  TODO - maybe restrict solv to those w/in some distance
 *      of solvent, building vaSolvent here.
 *
 *      Bill Ross, UCSF Jan 1995.
 */
void
ToolInitSolventPotential( UNIT uUnit, VARARRAY vaSolvent, 
        int *iPMinPotRes, int *iPMaxPotRes )
{
int             i, iSolvRes;
RESIDUE         rRes, *PrSolvRes;
LOOP            lResidues, lAtoms;
double          dMinPot, dMaxPot;

TMem();

        iSolvRes = iVarArrayElementCount( vaSolvent );
        PrSolvRes = PVAI( vaSolvent, RESIDUE, 0);
        for ( i=0; i<iSolvRes; i++, PrSolvRes++ ) {
                ATOM    aSolvAtom = aResidueImagingAtom( *PrSolvRes );
                double  dXSolv, dYSolv, dZSolv, dPot;
/*
VP1(("solv res %d 0x%lx\n", i, *PrSolvRes ));
*/
if ( aSolvAtom == NULL )
 DFATAL((" solvatom %d null\n", i ));
                dXSolv = vAtomPosition( aSolvAtom ).dX;
                dYSolv = vAtomPosition( aSolvAtom ).dY;
                dZSolv = vAtomPosition( aSolvAtom ).dZ;

                dPot = 0.0;
                lResidues = lLoop( (OBJEKT)uUnit, RESIDUES );
                while ( (rRes = (RESIDUE)oNext( &lResidues )) != NULL) {
                        ATOM    aSoluteAtom;

                        /* TODO - maybe possible to break here? */
                        if ( cResidueType( rRes ) == RESTYPESOLVENT )
                                continue;

                        /*
                         *  add this solute residue's atoms'
                         *      contributions to solvent residue's 
                         *      potential
                         */
                        lAtoms = lLoop( (OBJEKT)rRes, ATOMS );
                        while ( (aSoluteAtom = (ATOM)oNext(&lAtoms)) != NULL ) {
                                double  dX, dY, dZ, dR;

                                dX = dXSolv - vAtomPosition(aSoluteAtom).dX;
                                dY = dYSolv - vAtomPosition(aSoluteAtom).dY;
                                dZ = dZSolv - vAtomPosition(aSoluteAtom).dZ;

                                dR = sqrt( dX * dX  + dY * dY  +  dZ * dZ );
                                dPot += dAtomCharge( aSoluteAtom ) / dR;
                        }
                }
                /*
                 *  set the solvent res = the accumulated potential
                 */
                (*PrSolvRes)->dTemp = dPot;
/*
*/
                /*
                 *  keep track of max/min solvent res potentials
                 */
                if ( i == 0 ) {
                        dMinPot = dMaxPot = dPot;
                        *iPMinPotRes = *iPMaxPotRes = 0;
                } else if ( dPot < dMinPot ) {
                        *iPMinPotRes = i;
                        dMinPot = dPot;
                } else if ( dPot > dMaxPot ) {
                        *iPMaxPotRes = dPot;
                        dMaxPot = dPot;
                }
        }
}

void
ToolReplaceSolvent( UNIT uUnit, VARARRAY vaSolvent, int iSolv, 
        UNIT uIon, double dCharge, int *iPMinPotRes, int *iPMaxPotRes )
{
RESIDUE         *rPRes;
ATOM            aSolvAtom;
VECTOR          vPos;
UNIT            uPlace;
int             i, iSolvRes, iSeen;
double          dMinPot, dMaxPot;

        /*
         *  get position of 1st atom in solvent
         */
        rPRes = PVAI( vaSolvent, RESIDUE, iSolv );
        aSolvAtom = aResidueImagingAtom( *rPRes );
        vPos = vAtomPosition( aSolvAtom );

        /*
         *  Make a copy of ion unit, position it, & add to unit.
         */
        uPlace = (UNIT) oCopy( (OBJEKT)uIon );
        ContainerCenterAt((CONTAINER) uPlace, vPos );
        UnitJoin( uUnit, uPlace );

        /*
         *  delete the solvent residue & null the pointer in the array
         */
        REF( *rPRes );  /* bContainerRemove() needs this */
        if ( bContainerRemove( (CONTAINER)uUnit, (OBJEKT)*rPRes ) == FALSE)
                DFATAL(( "rmv solv %d failed\n", iSolv ));
        ContainerDestroy((CONTAINER *) rPRes );
        *rPRes = NULL;

        /*
         *  update the charges according to the new ion
         */
        dMinPot = FLT_MAX;
        dMaxPot = FLT_MIN;
        iSeen = 0;
        rPRes = PVAI( vaSolvent, RESIDUE, 0);
        iSolvRes = iVarArrayElementCount( vaSolvent );
        for ( i=0; i<iSolvRes; i++, rPRes++ ) {
                double  dX, dY, dZ, dR;

                if ( *rPRes == NULL )   /* already deleted */
                        continue;

                aSolvAtom = aResidueImagingAtom( *rPRes );
                dX = vPos.dX - vAtomPosition( aSolvAtom ).dX;
                dY = vPos.dY - vAtomPosition( aSolvAtom ).dY;
                dZ = vPos.dZ - vAtomPosition( aSolvAtom ).dZ;
                dR = sqrt( dX * dX  +  dY * dY  +  dZ * dZ );

                (*rPRes)->dTemp += dCharge / dR;

                if ( iSeen++ == 0 ) {
                        dMinPot = dMaxPot = (*rPRes)->dTemp;
                        *iPMinPotRes = *iPMaxPotRes = i;
                } else if ( (*rPRes)->dTemp < dMinPot ) {
                        *iPMinPotRes = i;
                        dMinPot = (*rPRes)->dTemp;
                } else if ( (*rPRes)->dTemp > dMaxPot ) {
                        *iPMaxPotRes = i;
                        dMaxPot = (*rPRes)->dTemp;
                }
        }
}



/*
--------------------------------------------------------

        Create bonds in a distance search
        
*/




        
/*
 *      iToolDistanceSearch
 *
 *      Author: Christian Schafmeister (1991)
 *      Amended: Bill Ross (add shortest bonds first) 1993
 *
 *      Do an N^2 search of all atoms, searching for
 *      pairs that are close enough to create a bond between
 *      them.  The argument bAbsoluteDistance is used to determine
 *      the criteria to use to create bonds.   If bAbsoluteDistance
 *      is TRUE then dCloseness is used as an absolute distance
 *      where atoms closer than that are bonded.  If FALSE then
 *      dCloseness is a multiplier of the average of the Van der Waals
 *      distances.
 *
 *      iOperation defines whether bonds should be created or
 *      if close contacts should be just reported.
 *
 *      Return the number of bonds created or close contacts found.
 */

typedef struct  {
        VECTOR          vPos;
        double          dR;
        ATOM            aAtom;
} BONDSEARCHt;

typedef struct {
        ATOM            aAtom1, aAtom2;
        double          dDist;
} PAIRt;

int
iToolDistanceSearch( CONTAINER cCont, double dCloseness, BOOL bAbsoluteDistance,
                                int iOperation )
{
LOOP                    lTemp;
VARARRAY                vaAtoms;
int                     iAtoms, i, j;
ATOM                    aAtom;
BONDSEARCHt*            bsPAtom1;
BONDSEARCHt*            bsPAtom2;
VECTOR                  vPos, vMin, vMax;
double                  dR, dX, dY, dZ;
double                  dDist;
STRING                  sTemp, sTemp1, sTemp2;
int                     iCount;
VARARRAY                vaPairs;
PAIRt                   *Ppair;
int                     iPairs;


                /* Collect all of the atoms in the CONTAINER cCont */
                /* into an array of positions and Rs for faster */
                /* checking of close contacts */

    iAtoms = 0;
    lTemp = lLoop( (OBJEKT)cCont, ATOMS );
    while ( oNext(&lTemp) != NULL ) iAtoms++;
    
    vaAtoms = vaVarArrayCreate( sizeof(BONDSEARCHt) );
    VarArraySetSize( vaAtoms, iAtoms );
 
    if ( iAtoms == 0 ) {
        VP0(( "No atoms to bond\n" ));
        return( 0 );
    }
    if ( iAtoms < 2 ) {
        VP0(( "Only one atom\n" ));
        return( 0 );
    }

    lTemp = lLoop( (OBJEKT)cCont, ATOMS );
    bsPAtom1= PVAI( vaAtoms, BONDSEARCHt, 0 );
    for (i=0; i<iAtoms; i++, bsPAtom1++ ) {
        aAtom = (ATOM) oNext(&lTemp);
        bsPAtom1->vPos = vAtomPosition(aAtom);
        bsPAtom1->aAtom = aAtom;
        bsPAtom1->dR = dToolAtomR( aAtom );
        if ( bsPAtom1->dR < 0.1 )       {
            /* R unknown or 0 */
            if ( iAtomElement( aAtom ) == HYDROGEN )
                bsPAtom1->dR = 1.0;
            else
                bsPAtom1->dR = ATOM_DEFAULT_RADIUS;
        }
        bsPAtom1->dR *= dCloseness;
    }

    if ( iOperation == DISTANCE_SEARCH_CREATE_BONDS ) {
        /*
         *  make an array for sorting distances; since
         *      a simple sizing of N^2/2 can exhaust
         *      memory for large #s of atoms, do a
         *      quick grid-based count to establish a 
         *      reasonable upper bound.
         */
        iCount = 0;
        bsPAtom1 = PVAI( vaAtoms, BONDSEARCHt, 0 );
        for ( i=0; i<iAtoms; i++, bsPAtom1++ ) {
                vPos = bsPAtom1->vPos;
                bsPAtom2 = bsPAtom1 + 1;
                for ( j=i+1; j<iAtoms; j++, bsPAtom2++ ) {
                        VECTOR  vPos2;

                        if ( bAtomBondedTo( bsPAtom1->aAtom, bsPAtom2->aAtom ) )
                                continue;
                        if ( iAtomElement(bsPAtom1->aAtom) == HYDROGEN &&
                             iAtomElement(bsPAtom2->aAtom) == HYDROGEN )
                                continue;

                        vPos2 = bsPAtom2->vPos;

                        if ( bAbsoluteDistance ) 
                                dR = dCloseness;
                        else
                                dR = 0.5 * (bsPAtom1->dR + bsPAtom2->dR);

                        VectorDef( &vMin, 
                                        dVX(&vPos2) - dR, 
                                        dVY(&vPos2) - dR, 
                                        dVZ(&vPos2) - dR );
                        VectorDef( &vMax, 
                                        dVX(&vPos2) + dR, 
                                        dVY(&vPos2) + dR, 
                                        dVZ(&vPos2) + dR );

                        if ( (dVX(&vPos) > dVX(&vMax))  ||
                             (dVY(&vPos) > dVY(&vMax))  ||
                             (dVZ(&vPos) > dVZ(&vMax)) ) 
                                continue;
                        if ( (dVX(&vPos) < dVX(&vMin))  ||
                             (dVY(&vPos) < dVY(&vMin))  ||
                             (dVZ(&vPos) < dVZ(&vMin)) )
                                continue;
                        iCount++;
                }
        }
        vaPairs = vaVarArrayCreate( sizeof(PAIRt) );
        VarArraySetSize( vaPairs, iCount );
        Ppair = PVAI( vaPairs, PAIRt, 0 );
        iPairs = 0;
    }

                /* Compare the position of every atom with every other */
                /* atom, looking for possible bonding contacts.  If one */
                /* is found and a bond can be created then create a bond */

    iCount = 0;
    bsPAtom1 = PVAI( vaAtoms, BONDSEARCHt, 0 );
    for ( i=0; i<iAtoms; i++, bsPAtom1++ ) {
        vPos = bsPAtom1->vPos;
        bsPAtom2 = bsPAtom1 + 1;
        for ( j=i+1; j<iAtoms; j++, bsPAtom2++ ) {
            VECTOR      vPos2;

            if ( bAtomBondedTo( bsPAtom1->aAtom, bsPAtom2->aAtom ) )
                continue;

            vPos2 = bsPAtom2->vPos;

            if ( bAbsoluteDistance ) 
                dR = dCloseness;
            else
                dR = 0.5 * (bsPAtom1->dR + bsPAtom2->dR);

           /*
            *  quick grid check to save cpu
            */
            VectorDef( &vMin, 
                        dVX(&vPos2) - dR, dVY(&vPos2) - dR, dVZ(&vPos2) - dR );
            VectorDef( &vMax, 
                        dVX(&vPos2) + dR, dVY(&vPos2) + dR, dVZ(&vPos2) + dR );

            if (  (dVX(&vPos) > dVX(&vMax)) ||
                  (dVY(&vPos) > dVY(&vMax)) ||
                  (dVZ(&vPos) > dVZ(&vMax)) || 
                  (dVX(&vPos) < dVX(&vMin)) ||
                  (dVY(&vPos) < dVY(&vMin)) ||
                  (dVZ(&vPos) < dVZ(&vMin)) )
                continue;

                /* Now check the actual distance between the centers */
        
            dX = dVX(&vPos) - dVX(&(bsPAtom2->vPos));
            dY = dVY(&vPos) - dVY(&(bsPAtom2->vPos));
            dZ = dVZ(&vPos) - dVZ(&(bsPAtom2->vPos));

            dR *= dR;
            dDist = dX*dX + dY*dY + dZ*dZ;

            if ( dDist < dR ) {

                        /* Check if there is already a bond between them */
                
                    iCount++;
                    switch ( iOperation ) {
                        case DISTANCE_SEARCH_CREATE_BONDS:

                            /*
                             *  create list of bondable atoms & the distances
                             */

                            if ( iAtomElement(bsPAtom1->aAtom) == HYDROGEN &&
                                 iAtomElement(bsPAtom2->aAtom) == HYDROGEN ) {
                                VP2(( 
                                "Hydrogens %s and %s within bonding distance\n",
                                   sContainerFullDescriptor( 
                                                (CONTAINER)bsPAtom1->aAtom, sTemp1 ),
                                   sContainerFullDescriptor( 
                                                (CONTAINER)bsPAtom2->aAtom, sTemp2 ) ));
                                continue;
                            }
                            Ppair->aAtom1 = bsPAtom1->aAtom;
                            Ppair->aAtom2 = bsPAtom2->aAtom;
                            Ppair->dDist = dDist;
                            Ppair++;
                            iPairs++;
                            break;
                        case DISTANCE_SEARCH_PRINT_WARNINGS:
			  if ( strcmp( sAtomName( bsPAtom1->aAtom ), "EPW" ) == 0 ||
			       strcmp( sAtomName( bsPAtom2->aAtom ), "EPW" ) == 0 ) {
			    /* Lone pair, ignore it  */
			  } else {
                            VP0(( "Warning: Close contact of %lf angstroms between %s and %s\n",
                                sqrt( dDist ),
                                sContainerFullDescriptor( (CONTAINER)bsPAtom1->aAtom, sTemp1 ),
                                sContainerFullDescriptor( (CONTAINER)bsPAtom2->aAtom, sTemp2 ) ));
			  }
			  break;
                        default:
                            DFATAL(( "Illegal iToolDistanceSearch operation: %d\n",
                                iOperation ));
                            break;
                    }
            }
        }                   
    }
    if ( iOperation == DISTANCE_SEARCH_CREATE_BONDS ) {

        /*
         *  assign bonds starting with the shortest distances
         */

        VP0(( "%d pairs of atoms within potential bonding distance\n", iCount));
        VP0(( "%d are not H-H pairs\n", iPairs ));

        Ppair = PVAI( vaPairs, PAIRt, 0 );
        SortByDouble( Ppair, iPairs, sizeof(PAIRt), &Ppair->dDist, TRUE );
        iCount = 0;
        for (i=0; i<iPairs; i++, Ppair++ ) {
                
                if ( bAtomCoordinationSaturated( Ppair->aAtom1 ) ) {
                        VP0(( "Atom %s has maximum number of bonds.\n",
                            sContainerFullDescriptor( (CONTAINER)Ppair->aAtom2, sTemp ) ));
                        continue;
                }
                if ( bAtomCoordinationSaturated( Ppair->aAtom2 ) ) {
                        VP0(( "Atom %s has maximum number of bonds.\n",
                            sContainerFullDescriptor( (CONTAINER)Ppair->aAtom2, sTemp ) ));
                        continue;
                }
                AtomBondToFlags( Ppair->aAtom1, Ppair->aAtom2, 
                                BONDTEMPORARY | BONDSINGLE );
                iCount++;
        }
        VarArrayDestroy( &vaPairs );
    }

    return(iCount);
}





/*
 *+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
 *
 *      Work with principle axis of the molecule.
 */


/*
 *      ToolOrientPrincipleAxisAlongCoordinateAxis
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Translate the GEOMETRIC center of the UNIT to the origin
 *      and calculate the principle axis, and then rotate the
 *      UNIT so that it's principle axis' are co-linear with
 *      the coordinate axis'.
 */
void
ToolOrientPrincipleAxisAlongCoordinateAxis( UNIT uUnit )
{
VECTOR          vTemp, vEigenvalues;
MATRIX          mMoment, mDiagonalize;
int             iSteps;
LOOP            lAtoms;
ATOM            aAtom;
VECTOR          vPos;
double          dDot;

        /* Translate the GEOMETRIC center of the UNIT to the origin */

    VectorDef( &vTemp, 0.0, 0.0, 0.0 );
    ContainerCenterAt((CONTAINER) uUnit, vTemp );

        /* Now calculate the moment of inertia in terms of GEOMETRY */
        /* meaning assign a mass of 1 to ALL ATOMS */

    MathOpMomentOfInertia( uUnit, mMoment );

        /* Diagonalize the moment of inertia matrix */

    VectorDef( &vEigenvalues, 1.0, 2.0, 3.0 );
    MatrixIdentity( mDiagonalize );
    iSteps = 1234;
    MathOpDiagonalize( mMoment, (double *)&vEigenvalues, mDiagonalize, &iSteps );
    // You should never treat a struct as an array, since the behavior
    // is undefined. EXTREMELY BAD BAD practice here, if you lie to the
    // compilers, eventually they will get their revenge. -- Mengjuei Hsieh
    // However I can't change that.....

        /* Check the handedness of the diagonalized matrix */
        /* If it is the wrong hand then it will transform the */
        /* UNIT into a mirror image */

    vPos = vVectorCross( (VECTOR *)mDiagonalize[0], (VECTOR *)mDiagonalize[1] );
    dDot = dVectorDot( &vPos, (VECTOR *)mDiagonalize[2] );
    MESSAGE(( "The handedness of the transformation is (+1=Right): %lf\n",
                dDot ));

        /* If the handedness of the matrix is wrong then change it */
    if ( dDot < 0.0 ) {
        mDiagonalize[2][0] *= -1.0;
        mDiagonalize[2][1] *= -1.0;
        mDiagonalize[2][2] *= -1.0;
    }

        /* Apply the matrix to the UNIT */

    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        vPos = vAtomPosition(aAtom);
        MatrixTimesVector( vPos, mDiagonalize, vPos );
        AtomSetPosition( aAtom, vPos );
    }
}





/*
 *----------------------------------------------------------------
 *
 */


/*
 *      bToolGeometricCenter
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a VECTOR that is obtained from the oObjekt
 *      argument.
 *
 *      If oObject is a CONTAINER then the geometric center
 *      of the ATOMs within the CONTAINER is returned.
 *
 *      If oObject is a LIST of CONTAINERs then the geometric
 *      center of all ATOMs withing all of the CONTAINERs is returned.
 *
 *      If oObject is a LIST of 3 ODOUBLEs then the vector
 *      defined by the three doubles is returned.
 *
 *      Return TRUE if the geometric center was defined.
 *
 */
BOOL
bToolGeometricCenter( OBJEKT oObjekt, VECTOR *vPCenter )
{
LISTLOOP        llElements;
OBJEKT          oElement, oNum;
CONTAINER       cCont;
VECTOR          vPos;
int             i, iCount;
LOOP            lAtoms;
ATOM            aAtom;
BOOL            bVector;
double          daElements[3];

    if ( bObjectInClass( oObjekt, CONTAINERid ) ) {
        *vPCenter = vContainerGeometricCenter((CONTAINER)oObjekt);
    } else if ( iObjectType(oObjekt) == LISTid ) {
        bVector = TRUE;
        llElements = llListLoop((LIST)oObjekt);
        i = 0;
        while ( (oElement = oListNext(&llElements)) ) {
            if ( iObjectType(oElement) == ASSOCid ) {
                oNum = oAssocObject(oElement);
                if ( iObjectType(oNum) == ODOUBLEid ) {
                    if ( i<3 ) {
                        daElements[i] = dODouble(oNum);
                        i++;
                    } else {
                        bVector = FALSE;
                        break;
                    }
                } else {
                    bVector = FALSE;
                    break;
                }
            } else {
                bVector = FALSE;
                break;
            }
        }
        if ( bVector ) {
            VectorDef( vPCenter, 
                        daElements[0], 
                        daElements[1], 
                        daElements[2] );
        } else {
            VectorDef( &vPos, 0.0, 0.0, 0.0 );
            iCount = 0;
            llElements = llListLoop( (LIST)oObjekt );
            while ( (oElement = oListNext(&llElements)) ) {
                if ( iObjectType(oElement) == ASSOCid ) {
                    cCont = (CONTAINER)oAssocObject(oElement);
                } else 
                    cCont = (CONTAINER)oElement;
                if ( !bObjectInClass( (OBJEKT)cCont, CONTAINERid ) ) 
                        return(FALSE);
                lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
                while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
                    vPos = vVectorAdd( &vPos, &vAtomPosition(aAtom) );
                    iCount++;
                }
            }
            if ( iCount > 0 ) {
                *vPCenter = vVectorTimesScalar( &vPos, 1/iCount );
            }
        }
    }
    return(TRUE);
}


/*
 *      lToolListOfResidues
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a list of residues by finding the residues within
 *      the UNIT with the correct sequence numbers.
 *
 *      Arguments:
 *              [0] -   Unit with residues.
 *              [1] -   List of residue numbers or ranges.
 *
 *      Ranges are expressed as two element lists within the main
 *      list.
 *
 */
LIST
lToolListOfResidues( UNIT uUnit, LIST lResidues )
{
LIST            lList;
LISTLOOP        llNumbers, llTemp;
CONTAINER       cAdd;
int             iSeq, iSeqStart, iSeqStop, i;
ASSOC           aAss;
OBJEKT          oObj;
                
    lList = (LIST)oCreate(LISTid);
    llNumbers = llListLoop( lResidues );
    while ( ( aAss = (ASSOC)oListNext(&llNumbers) ) != NULL ) {
        oObj = oAssocObject(aAss);
        if ( iObjectType(oObj) == ODOUBLEid ) {
            iSeq = (int)(dODouble(oObj));
            cAdd = cContainerFindSequence((CONTAINER) uUnit, RESIDUEid, iSeq );
            if ( cAdd != NULL ) {
                ListAdd( lList, (OBJEKT)cAdd );
                MESSAGE(( "Adding residue %s\n", sContainerName(cAdd) ));
            }
            
                /* If the object is a list then get the two values */
                /* that should be in it and search for the range */
                /* of RESIDUEs that they describe */
        } else {
            llTemp = llListLoop( (LIST)oObj );
            aAss = (ASSOC)oListNext(&llTemp);
            if ( aAss == NULL ) {
                VP0(( "Illegal residue range specified.\n" ));
                break;
            }
            iSeqStart = (int)(dODouble(oAssocObject(aAss)));
            aAss = (ASSOC)oListNext(&llTemp);
            if ( aAss == NULL ) {
                VP0(( "Illegal residue range specified.\n" ));
                break;
            }
            iSeqStop  = (int)(dODouble(oAssocObject(aAss)));
            for ( i=iSeqStart; i<=iSeqStop; i++ ) {
                cAdd = cContainerFindSequence((CONTAINER) uUnit, RESIDUEid, i );
                if ( cAdd != NULL ) {
                    ListAdd( lList, (OBJEKT)cAdd );
                    MESSAGE(( "Adding residue %s\n", sContainerName(cAdd) ));
                }
            }
        }
    }
    return(lList);
}


void
ToolOctBoxCheck( UNIT uSolute, double *dPBuf, BOOL bMsg )
{
LOOP            lTemp;
ATOM            aAtom;
VECTOR          vPos;
double          dXmin, dYmin, dZmin, dXmax, dYmax, dZmax;
double          dXhalf, dYhalf, dZhalf, dXunit, dYunit, dZunit;
double          dX, dY, dZ, dTmp, dMax, dBmax; 
int             iFirst = 1;


        /*
         *  get atom centers' x,y,z max,min
         */
        lTemp = lLoop( (OBJEKT)uSolute, ATOMS );
        while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
                vPos = vAtomPosition(aAtom);
                dX = dVX(&vPos);
                dY = dVY(&vPos);
                dZ = dVZ(&vPos);

                if ( iFirst ) {
                        dXmax = dXmin = dX;
                        dYmax = dYmin = dY;
                        dZmax = dZmin = dZ;
                        iFirst = 0;
                        continue;
                }
                if ( dX > dXmax ) 
                        dXmax = dX;
                else if ( dX < dXmin )
                        dXmin = dX;
                if ( dY > dYmax ) 
                        dYmax = dY;
                else if ( dY < dYmin )
                        dYmin = dY;
                if ( dZ > dZmax ) 
                        dZmax = dZ;
                else if ( dZ < dZmin )
                        dZmin = dZ;
        }

        /*
         *  calc halfbox on centers
         */
        dXhalf = 0.5 * (dXmax - dXmin) + dPBuf[0];
        dYhalf = 0.5 * (dYmax - dYmin) + dPBuf[1];
        dZhalf = 0.5 * (dZmax - dZmin) + dPBuf[2];

        /*
         *  find unit vector of diagonal
         */
        dX = dYhalf * dZhalf;
        dY = dXhalf * dZhalf;
        dZ = dXhalf * dYhalf;

        dTmp = 1.0 / sqrt( dX*dX + dY*dY + dZ*dZ );

        dXunit = dX * dTmp;
        dYunit = dY * dTmp;
        dZunit = dZ * dTmp;
        
        /*
         *  find max atom distance from origin along the diagonal
         */
        dMax = 0.0;
        lTemp = lLoop( (OBJEKT)uSolute, ATOMS );
        while ( (aAtom = (ATOM)oNext(&lTemp)) ) {
                vPos = vAtomPosition(aAtom);

                dX = dVX(&vPos);
                dY = dVY(&vPos);
                dZ = dVZ(&vPos);

                dTmp =  fabs(dX) * dXunit + 
                        fabs(dY) * dYunit + 
                        fabs(dZ) * dZunit;

                if ( dTmp > dMax )
                        dMax = dTmp;
        }

        /*
         *  calc distance of diagonal face from origin
         */
        dBmax = 0.5 * sqrt( dXhalf*dXhalf + dYhalf*dYhalf + dZhalf*dZhalf );
        
        /*
         *  see if diagonal clearance is satisfied
         */
        dTmp = dMax + dPBuf[3];
        if ( dTmp <= dBmax ) {
                if ( dPBuf[3] == 0.0 )
                        VP0(( "(Diagonal clearance is %f)\n", dMax ));
                return;
        }

        /*
         *  not satisfied: scale up box
         */
        dTmp /= dBmax;
        if ( bMsg )
          VP0(( 
            "Scaling up box by a factor of %f to meet diagonal cut criterion\n",
            dTmp ));

        dPBuf[0] *= dTmp;
        dPBuf[1] *= dTmp;
        dPBuf[2] *= dTmp;

}
