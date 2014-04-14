/*
 *      Class:  UNIT
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
 *      Superclass:
 *              CONTAINER, LOOP
 *
 *      Description:
 *              A UNIT is a subclass of CONTAINER.
 *              UNITS can contain molecules, residues, and atoms.
 *              UNITS are used to contain the entire molecular
 *              system plus other information about the system.
 *
 *      NOTE:   OFF files are used to write UNITs to files.
 *              In OFF files there is no implicit ordering
 *              of ANY DATA WHATSOEVER.  The PROGRAMMER is
 *              to assume that there is no order regardless
 *              of the ordering that the code in this file
 *              generates.
 *
 *      NOTE2:  When setting up tables for writing to OFF files,
 *              Indices are FORTRAN indices, where the first
 *              element has index=1.
 *
 */
 
#include        "basics.h"
#include        "vector.h"
#include        "classes.h"
#include        "restraint.h"
#include        "bag.h"
#include        "dictionary.h"
#include        "database.h"
#include        "parmLib.h"
#include        "unitio.h"
#include        "defaults.h"

/*
--------------------------------------------------------------------

        Type definitions, static variables etc.
*/


/*
 *----------------------------------------------------------------------

        Private routines
*/

/*
 *      zbUnitCheckBondParameters
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      For all of the bonds, search for their parameters
 *      within the passed PARMLIB.  If they are
 *      found then continue, otherwise
 *      add the bond parameter to the PARMSET.
 *      If either of the ATOMs in the bond are to be perturbed then
 *      do the same with the perturbation parameter.
 *
 *      Return TRUE if there was a problem generating parameters.
 */
static BOOL
zbUnitCheckBondParameters( PARMLIB plLib, UNIT uUnit)
{
LOOP            lTemp;
ATOM            aAtom1, aAtom2;
BOOL            bFailedGeneratingParameters;
STRING          sAtom1, sAtom2;
PARMSET         psTemp;
int             iIndex;

    bFailedGeneratingParameters = FALSE;
    
    VP0(( "Checking for bond parameters.\n" ));

    lTemp = lLoop( (OBJEKT)uUnit, BONDS );
    while ( oNext(&lTemp) != NULL ) {
        LoopGetBond( &lTemp, &aAtom1, &aAtom2 );
        strcpy( sAtom1, sAtomType(aAtom1) );
        strcpy( sAtom2, sAtomType(aAtom2) );
        PARMLIB_LOOP( plLib, psTemp,
                ( iIndex = iParmSetFindBond( psTemp, sAtom1, sAtom2 )));
        if ( iIndex == PARM_NOT_FOUND ) {
                bFailedGeneratingParameters = TRUE;
                VP0(( "Could not find bond parameter for: %s - %s\n", 
                    sAtom1, sAtom2 ));
        }

        if ( bAtomFlagsSet( aAtom1, ATOMPERTURB ) ||
                bAtomFlagsSet( aAtom2, ATOMPERTURB ) ) {

            /* Note that the bond is perturbed and whether or */
            /* not it is on the boundary between perturbed and */
            /* non-perturbed */
            if ( bAtomFlagsSet( aAtom1, ATOMPERTURB ))
                    strcpy( sAtom1, sAtomPertType( aAtom1 ));
            if ( bAtomFlagsSet( aAtom2, ATOMPERTURB ))
                    strcpy( sAtom2, sAtomPertType( aAtom2 ));
            PARMLIB_LOOP( plLib, psTemp,
                    ( iIndex = iParmSetFindBond( psTemp, sAtom1, sAtom2 ))) ;
            if ( iIndex == PARM_NOT_FOUND ) {
                    bFailedGeneratingParameters = TRUE;
                    VP0(( "No bond parameter for: %s - %s\n", sAtom1, sAtom2 ));
            }
        }
    }

    return(bFailedGeneratingParameters);
}


/*
 *      zbUnitCheckAngleParameters
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      For all of the angles, search for their parameters
 *      within the UNITs PARMSET.  If they are
 *      found then set the index to the entry, otherwise
 *      add the angle parameter to the PARMSET and set the index.
 *      If any of the ATOMs in the angle  are to be perturbed then
 *      do the same with the perturbation parameter.
 *
 *      Return TRUE if there was a problem generating parameters.
 */
static BOOL
zbUnitCheckAngleParameters( PARMLIB plLib, UNIT uUnit)
{
LOOP            lTemp;
ATOM            aAtom1, aAtom2, aAtom3;
BOOL            bFailedGeneratingParameters;
STRING          sAtom1, sAtom2, sAtom3;
PARMSET         psTemp;
int             iTemp = PARM_NOT_FOUND;

    bFailedGeneratingParameters = FALSE;

            /* Now generate the ANGLE table */
    VP0(( "Checking for angle parameters.\n" ));

    lTemp = lLoop( (OBJEKT)uUnit, ANGLES );
    while ( oNext(&lTemp) != NULL ) {
        LoopGetAngle( &lTemp, &aAtom1, &aAtom2, &aAtom3 );
        strcpy( sAtom1, sAtomType( aAtom1 ));
        strcpy( sAtom2, sAtomType( aAtom2 ));
        strcpy( sAtom3, sAtomType( aAtom3 ));

/* TODO:Fix this UNGODLY HACK, the PARMSET type has to be modified to */
/* TODO:Allow the user to specify interactions that should be ignored */
/* TODO:Like HW-HW-OW angles for TIP3 waters */

        if ( zbUnitIgnoreAngle( sAtom1, sAtom2, sAtom3 )) goto IGNORE1;

        PARMLIB_LOOP( plLib, psTemp, 
                ( iTemp = iParmSetFindAngle( psTemp, sAtom1, 
                                                sAtom2, sAtom3 ) ) );
        if ( iTemp == PARM_NOT_FOUND ) {
                bFailedGeneratingParameters = TRUE;
                VP0(( "Could not find angle parameter: %s - %s - %s\n", 
                            sAtom1, sAtom2, sAtom3 ));
        }

IGNORE1:
        ;

DONTIGNORE1:
        if ( bAtomFlagsSet( aAtom1, ATOMPERTURB ) ||
                bAtomFlagsSet( aAtom2, ATOMPERTURB ) ||
                bAtomFlagsSet( aAtom3, ATOMPERTURB ) ) {

                    /* Note that the angle is perturbed and whether or */
                    /* not it is on the boundary between perturbed and */
                    /* non-perturbed */
            if ( bAtomFlagsSet( aAtom1, ATOMPERTURB ))
                strcpy( sAtom1, sAtomPertType( aAtom1 ));
            if ( bAtomFlagsSet( aAtom2, ATOMPERTURB ))
                strcpy( sAtom2, sAtomPertType( aAtom2 ));
            if ( bAtomFlagsSet( aAtom3, ATOMPERTURB ))
                strcpy( sAtom3, sAtomPertType( aAtom3 ));

/* TODO:Fix this UNGODLY HACK, the PARMSET type has to be modified to */
/* TODO:Allow the user to specify interactions that should be ignored */
/* TODO:Like HW-HW-OW angles for TIP3 waters */

            if ( zbUnitIgnoreAngle( sAtom1, sAtom2, sAtom3 )) 
                        goto IGNORE2;

            PARMLIB_LOOP( plLib, psTemp, 
                        ( iTemp = iParmSetFindAngle( psTemp, sAtom1,
                                        sAtom2, sAtom3 )));
            if ( iTemp == PARM_NOT_FOUND ){
                    bFailedGeneratingParameters = TRUE;
                    VP0(( "Can't find angle parameter: %s - %s - %s\n", 
                            sAtom1, sAtom2, sAtom3 ));
            }
        }
IGNORE2:
        ;
    }
    return(bFailedGeneratingParameters);
}

/*
 *      zbUnitCheckTorsionParameters
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      For all of the angles, search for their parameters
 *      within the UNITs PARMSET.  If they are
 *      found then set the index to the entry, otherwise
 *      add the torsion parameter to the PARMSET and set the index.
 *      If any of the ATOMs in the torsion are to be perturbed then
 *      do the same with the perturbation parameter.
 *
 *      Return TRUE if there was a problem generating parameters.
 */
static BOOL
zbUnitCheckTorsionParameters( PARMLIB plLib, UNIT uUnit)
{
LOOP            lTemp;
ATOM            aAtom1, aAtom2, aAtom3, aAtom4;
STRING          sAtom1, sAtom2, sAtom3, sAtom4;
STRING          sPert1, sPert2, sPert3, sPert4;
TORSION         tTorsion, tPertTorsion;
BOOL            bPerturbTorsion;
PARMSET         psTemp;

#define         MAX_N           9999


    /* NOTE: In order for the 1-4 interactions to be calculated */
            /* properly, add constraint bonds and angles AFTER */
            /* the torsions are added.  This allows the 1-4 interaction */
            /* checker to use the bond lists and angle lists to check */
            /* connectivity */

    lTemp = lLoop( (OBJEKT)uUnit, PROPERS );
    while ( oNext(&lTemp) != NULL ) {
        LoopGetTorsion( &lTemp, &aAtom1, &aAtom2, &aAtom3, &aAtom4 );

                /* Define the names of the unperturbed ATOMs */         
        strcpy( sAtom1, sAtomType(aAtom1) );
        strcpy( sAtom2, sAtomType(aAtom2) );
        strcpy( sAtom3, sAtomType(aAtom3) );
        strcpy( sAtom4, sAtomType(aAtom4) );

                /* Don't include torsions relating to extra points:  */
        if( GDefaults.iDeleteExtraPointAngles ){
            if( strcmp( sAtom1, "EP" ) == 0 || strcmp( sAtom4, "EP" ) == 0 )
            continue;
        }

                /* Check if the torsion is to be perturbed, if it */
                /* is then set flags saying so, and create a TORSION for */
                /* the perturbation */
        bPerturbTorsion = FALSE;
        if ( bAtomFlagsSet( aAtom1, ATOMPERTURB ) ||
                bAtomFlagsSet( aAtom2, ATOMPERTURB ) ||
                bAtomFlagsSet( aAtom3, ATOMPERTURB ) ||
                bAtomFlagsSet( aAtom4, ATOMPERTURB )) {
                            
                /* Define the names of the perturbed atoms */
            if ( bAtomFlagsSet( aAtom1, ATOMPERTURB ) )
                 strcpy( sPert1, sAtomPertType(aAtom1) );
            else strcpy( sPert1, sAtom1 );
            if ( bAtomFlagsSet( aAtom2, ATOMPERTURB ) )
                 strcpy( sPert2, sAtomPertType(aAtom2) );
            else strcpy( sPert2, sAtom2 );
            if ( bAtomFlagsSet( aAtom3, ATOMPERTURB ) )
                 strcpy( sPert3, sAtomPertType(aAtom3) );
            else strcpy( sPert3, sAtom3 );
            if ( bAtomFlagsSet( aAtom4, ATOMPERTURB ) )
                 strcpy( sPert4, sAtomPertType(aAtom4) );
            else strcpy( sPert4, sAtom4 );
        }
        
                /* First search through the PARMSET for the */
                /* torsion parameters that we need. */
        tTorsion = tParmSetTORSIONCreate();
        if ( bPerturbTorsion ) {
            tPertTorsion = tParmSetTORSIONCreate();
        }

                /* Now search through the PARMLIB for all other */
                /* torsion terms that might not already be in */
                /* the PARMSET */
                        
        PARMLIB_LOOP_ALL( plLib, psTemp ) {
            iParmSetFindProperTerms( psTemp,
                                        tTorsion, FALSE,
                                        sAtom1, sAtom2,
                                        sAtom3, sAtom4 );
            if ( bPerturbTorsion ) {
                    iParmSetFindProperTerms( psTemp,
                                            tPertTorsion, FALSE,
                                            sPert1, sPert2,
                                            sPert3, sPert4 );
            }
        }
        ParmSetTORSIONDestroy( &tTorsion );
        if ( bPerturbTorsion ) ParmSetTORSIONDestroy( &tPertTorsion );
    }
    return(FALSE);
}


/*
 *      zbUnitParmsMissing
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Author: David Rivkin
 *
 *      Scan the unit and find any missing bond, 
 *      angle or torsion (Proper and Improper)
 *      parameters in the parmLib or the units' parameter set.
 *      If the parameter set is given, put the missing parms in that parmSet.
 *
 */
 
static BOOL
zbUnitParmsMissing( UNIT uUnit, PARMLIB plParameters)
{
BOOL    bMissing = FALSE;

    if ( plParameters == NULL ) {
        return(TRUE);
    }
    
    if ( !uUnit ) {
        DFATAL(( "Unit is NULL!\n" ));
    }

    bMissing |= 
        zbUnitCheckBondParameters( plParameters, uUnit);
                
        /* Now check the ANGLE parms */

    bMissing |=
        zbUnitCheckAngleParameters( plParameters, uUnit);

        /* Now check the TORSION parms */

    bMissing |=
        zbUnitCheckTorsionParameters( plParameters, uUnit);

    return( bMissing );
}



/*
 *======================================================================
 *
 *        Public routines
 *
 */



/*
 *      zbUnitIgnoreAngle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This is an ungodly hack that checks if the types are defining
 *      an angle interaction between HW-HW-OW.  If they are then return TRUE
 *      this will cause the code to ignore this interaction.
 *      This will allow LEAP to handle TIP3 waters.
 */
BOOL
zbUnitIgnoreAngle( STRING sA, STRING sB, STRING sC )
{

    if ( strcmp( sA, "OW" ) == 0 ) {
        if ( strcmp( sB, "HW" ) == 0 &&
             strcmp( sC, "HW" ) == 0 ) return(TRUE);
    }
	if( ! GDefaults.iFlexibleWater ){
		if ( strcmp( sB, "OW" ) == 0 ) {
			if ( strcmp( sA, "HW" ) == 0 &&
				 strcmp( sC, "HW" ) == 0 ) return(TRUE);
		}
	}
    if ( strcmp( sC, "OW" ) == 0 ) {
        if ( strcmp( sB, "HW" ) == 0 &&
             strcmp( sA, "HW" ) == 0 ) return(TRUE);
    }

/*  delete all angles related to extra points    */
        if( GDefaults.iDeleteExtraPointAngles ){
        if ( strcmp( sA, "EP" ) == 0 || strcmp( sC, "EP" ) == 0 ) return(TRUE);
        }

    return(FALSE);
}
 
/*
 *      uUnitCreate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a UNIT.
 *
 *      Return:
 *              Return a pointer to the UNIT created.
 */
UNIT
uUnitCreate()
{
UNIT    m;

    MALLOC( m, UNIT, sizeof(UNITt) );

        /* Allocate BAGs for RESTRAINTs */

    m->bRestraints = bBagCreate();

        /* Zero the VARARRAYs used to store the table representation of */
        /* the UNIT used when writing the UNIT to a file */
    m->cHeader.dDisp = NULL;
    m->psParameters = NULL;
    m->vaAtoms = NULL;
    m->vaBonds = NULL;
    m->vaAngles = NULL;
    m->vaTorsions = NULL;
    m->vaConnectivity = NULL;
    m->vaRestraints = NULL;
    m->vaResidues = NULL;
    m->vaMolecules = NULL;
    m->vaHierarchy = NULL;
    m->vaConnect = NULL;
    m->vaGroupNames = NULL;
    m->vaGroupAtoms = NULL;
    UnitSetMode( m, UNITNORMAL );
    UnitDefineFlags( m, 0 );
    UnitSetUseBox( m, FALSE );
    UnitSetBoxOct( m, FALSE );
    UnitSetBox( m, 0.0, 0.0, 0.0 );
    UnitSetBeta( m, 90.0*DEGTORAD );
    UnitSetHead( m, NULL );
    UnitSetTail( m, NULL );
    UnitSetSolventCap( m, 0.0, 0.0, 0.0, 0.0 );
    UnitSetUseSolventCap( m, FALSE );
    m->dAtomGroups = dDictionaryCreate();
    return(m);
}

/*
 *      UnitDestroy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Destroy the UNIT and all of the residues/atoms within it.
 *      Should only be called from ContainerDestroy() so that
 *      container contents can be zapped.
 *
 *      Arguments:
 *      oObject -       A pointer to the object.
 *
 */
void
UnitDestroy( UNIT *uPUnit )
{
BAGLOOP         blRestraints;
RESTRAINT       rRest;
LIST            lGroup;
DICTLOOP        dlGroups;

                /* Destroy everything in the UNIT */

    if ( (*uPUnit)->psParameters != NULL )
        ParmSetDestroy( &((*uPUnit)->psParameters) );

                /* Destroy the RESTRAINTs */

    blRestraints = blBagLoop((*uPUnit)->bRestraints);
    while ( (rRest = (RESTRAINT)PBagNext(&blRestraints)) ) {
        RestraintDestroy( &rRest );
    }
    BagDestroy( &((*uPUnit)->bRestraints) );

                /* Destroy the ATOM groups DICTIONARY */
                /* TODO: Check if this is slow because ATOMs */
                /* TODO: were DEREF'd before groups were destroyed */

    dlGroups = ydlDictionaryLoop( (*uPUnit)->dAtomGroups );
    while ( yPDictionaryNext( (*uPUnit)->dAtomGroups, &dlGroups ) ) {
        lGroup = (LIST)PDictLoopData(dlGroups);
        ListDestroy( &lGroup );
    }

    DictionaryDestroy( &((*uPUnit)->dAtomGroups) );

    if ( (*uPUnit)->vaAtoms )
        VarArrayDestroy( &(*uPUnit)->vaAtoms );
    if ( (*uPUnit)->vaBonds )
        VarArrayDestroy( &(*uPUnit)->vaBonds );
    if ( (*uPUnit)->vaAngles )
        VarArrayDestroy( &(*uPUnit)->vaAngles );
    if ( (*uPUnit)->vaTorsions )
        VarArrayDestroy( &(*uPUnit)->vaTorsions );
    if ( (*uPUnit)->vaConnectivity )
        VarArrayDestroy( &(*uPUnit)->vaConnectivity );
    if ( (*uPUnit)->vaRestraints )
        VarArrayDestroy( &(*uPUnit)->vaRestraints );
    if ( (*uPUnit)->vaResidues )
        VarArrayDestroy( &(*uPUnit)->vaResidues );
    if ( (*uPUnit)->vaMolecules )
        VarArrayDestroy( &(*uPUnit)->vaMolecules );
    if ( (*uPUnit)->vaHierarchy )
        VarArrayDestroy( &(*uPUnit)->vaHierarchy );
    if ( (*uPUnit)->vaConnect )
        VarArrayDestroy( &(*uPUnit)->vaConnect );
    if ( (*uPUnit)->vaGroupNames )
        VarArrayDestroy( &(*uPUnit)->vaGroupNames );
    if ( (*uPUnit)->vaGroupAtoms )
        VarArrayDestroy( &(*uPUnit)->vaGroupAtoms );

                /* Destroy the UNIT */

    FREE( *uPUnit );
    *uPUnit = NULL;
}






/*
 *      UnitDescribe
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to return a descriptor
 *      string for the UNIT object.
 *      It is the callers responsibility to ensure that there is room
 *      enough in the buffer to describe the object.
 *
 *      Arguments:
 *      oObject -       The object to describe.
 *
 */
void
UnitDescribe( UNIT uUnit )
{
LOOP            lContents;
OBJEKT          oObj;
STRING          sTemp;
double          dX, dY, dZ, dR;
BAGLOOP         blRestraints;
ATOM            aAtom1, aAtom2, aAtom3, aAtom4;
STRING          sAtom1, sAtom2, sAtom3, sAtom4;
double          dKr, dR0, dKt, dT0, dKp, dP0, dN;
RESTRAINT       rRest;
DICTLOOP        dlGroups;
LIST            lAtoms;
LISTLOOP        llAtoms;
ATOM            aAtom;
STRING          sAtom;


    VP0(( "UNIT name: %s\n", uUnit->cHeader.sName ));
    if ( bUnitHeadUsed(uUnit) ) 
        strcpy( sAtom1, sContainerFullDescriptor( (CONTAINER)aUnitHead(uUnit), sTemp ) );
    else strcpy( sAtom1, "null" );
    VP0(( "Head atom: %s\n", sAtom1 ));
    if ( bUnitTailUsed(uUnit) ) 
        strcpy( sAtom1, sContainerFullDescriptor( (CONTAINER)aUnitTail(uUnit), sTemp ) );
    else strcpy( sAtom1, "null" );
    VP0(( "Tail atom: %s\n", sAtom1 ));
    
    if ( bUnitUseBox(uUnit) ) {
        UnitGetBox( uUnit, &dX, &dY, &dZ );
        VP0(( "Periodic box: %10.5lf, %10.5lf, %10.5lf\n",
                dX, dY, dZ ));
    }
    if ( bUnitUseSolventCap(uUnit) ) {
        UnitGetSolventCap( uUnit, &dX, &dY, &dZ, &dR );
        VP0(( "Solvent cap origin:%10.5lf, %10.5lf, %10.5lf  radius:%10.5lf\n",
                dX, dY, dZ, dR ));
    }
    if ( iBagSize(uUnit->bRestraints) > 0 ) {
        blRestraints = blBagLoop(uUnit->bRestraints);
        while ( (rRest = (RESTRAINT)PBagNext(&blRestraints)) ) {
            if ( iRestraintType(rRest) == RESTRAINTBOND ) {
                RestraintBondGet( rRest, &aAtom1, &aAtom2, &dKr, &dR0 );
                VP0(( "Restraint BOND: %s - %s   Kr=%lf  R0=%lf\n",
                        sContainerFullDescriptor( (CONTAINER)aAtom1, sAtom1 ),
                        sContainerFullDescriptor( (CONTAINER)aAtom2, sAtom2 ),
                        dKr, dR0 ));
            }
        }
        blRestraints = blBagLoop(uUnit->bRestraints);
        while ( (rRest = (RESTRAINT)PBagNext(&blRestraints)) ) {
            if ( iRestraintType(rRest) == RESTRAINTANGLE ) {
                RestraintAngleGet( rRest, &aAtom1, &aAtom2, &aAtom3,
                                  &dKt, &dT0 );
                VP0(( "Restraint ANGLE: %s - %s - %s  Kt=%lf  T0=%lf\n",
                        sContainerFullDescriptor( (CONTAINER)aAtom1, sAtom1 ),
                        sContainerFullDescriptor( (CONTAINER)aAtom2, sAtom2 ),
                        sContainerFullDescriptor( (CONTAINER)aAtom3, sAtom3 ),
                        dKt, dT0 ));
            }
        }
        blRestraints = blBagLoop(uUnit->bRestraints);
        while ( (rRest = (RESTRAINT)PBagNext(&blRestraints)) ) {
            if ( iRestraintType(rRest) == RESTRAINTTORSION ) {
                RestraintTorsionGet( rRest, &aAtom1, &aAtom2, &aAtom3, &aAtom4,
                                  &dKp, &dP0, &dN );
        VP0(( "Restraint TORSION: %s - %s - %s - %s  Kt=%lf  T0=%lf  N=%lf\n",
                        sContainerFullDescriptor( (CONTAINER)aAtom1, sAtom1 ),
                        sContainerFullDescriptor( (CONTAINER)aAtom2, sAtom2 ),
                        sContainerFullDescriptor( (CONTAINER)aAtom3, sAtom3 ),
                        sContainerFullDescriptor( (CONTAINER)aAtom4, sAtom4 ),
                        dKp, dP0, dN ));
            }
        }
    }
    BasicsResetInterrupt();
    VP0(( "Contents: \n" ));
    lContents = lLoop( (OBJEKT)uUnit, DIRECTCONTENTSBYSEQNUM );
    while ( (oObj = oNext(&lContents )) ) {
        VP0(( "%s\n", sContainerDescriptor( (CONTAINER)oObj, sTemp ) ));
        if ( bBasicsInterrupt() ) {
            BasicsResetInterrupt();
            VP0(( "Interrupted\n" ));
            break;
        }
    }
    if ( iDictionaryElementCount(uUnit->dAtomGroups) != 0 ) {
        int     i;
        VP0(( "Atom groups:\n" ));
        dlGroups = ydlDictionaryLoop( uUnit->dAtomGroups );
        while ( yPDictionaryNext( uUnit->dAtomGroups, &dlGroups ) ) {
            VP0(( "  Group: %s\n", sDictLoopKey(dlGroups) ));
            lAtoms = (LIST)PDictLoopData(dlGroups);
            llAtoms = llListLoop(lAtoms);
            i = 0;
            while ( (aAtom = (ATOM)oListNext(&llAtoms)) ) {
                VP0(( "    %s\n", sContainerFullDescriptor((CONTAINER)aAtom,sAtom) ));
                i++;
            }
            if ( !i )
                VP0(( "    (empty group)\n" ));
        }
    }
    
}





/*
 *      uUnitDuplicate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a duplicate of the molecular system.
 *      Should only be called from cContainerDuplicate()
 *      which is expected to duplicate the CONTAINER variables.
 *
 *TODO: Duplicate the RESTRAINTS, and the groups.
 */
UNIT
uUnitDuplicate( UNIT uOld )
{
UNIT            uNew;

    MALLOC( uNew, UNIT, sizeof(UNITt) );
    memcpy( uNew, uOld, sizeof(UNITt) );

    uNew->bRestraints = bBagCreate();
    uNew->dAtomGroups = dDictionaryCreate();

    return(uNew);
}



void
UnitSetUseBox( UNIT uUnit, BOOL b )
{
        if ( b ) {
            UnitSetFlags(uUnit, UNITUSEBOUNDINGBOX);
        } else {
            UnitResetFlags(uUnit, UNITUSEBOUNDINGBOX);
        }
}

void
UnitSetBoxOct( UNIT uUnit, BOOL b )
{
        if ( b ) {
            UnitSetFlags(uUnit, UNITBOXOCT);
        } else {
            UnitResetFlags(uUnit, UNITBOXOCT);
        }
}

void
UnitSetUseSolventCap( UNIT uUnit, BOOL b ) 
{
        if ( b ) {
            UnitSetFlags(uUnit, UNITUSESOLVENTCAP);
        } else {
            UnitResetFlags(uUnit, UNITUSESOLVENTCAP);
        }
}




/*
 *      UnitResetPointers
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Reset all internal pointers.
 *      When copies are made of CONTAINERs, internal pointers, such
 *      as terminis, are left pointing to the old CONTAINERs.
 *      However, the CONTAINER class supports an extra pointer
 *      on the old CONTAINER which points to the new CONTAINER.
 *      This can be used to reset the internal pointers of the RESIDUE.
 *      The value of the pointer in new CONTAINERs is NULL.
 *
 *      EVERY internal pointer of the CONTAINER has to be resolved.
 *
 *TODO: Reset pointers for RESTRAINTs
 */
void
UnitResetPointers( UNIT uUnit )
{
                /* Reset the pointers for the connect pointers */

    UnitSetHead( uUnit, cContainerCopyPointer(aUnitHead(uUnit)) );
    UnitSetTail( uUnit, cContainerCopyPointer(aUnitTail(uUnit)) );
}




#if 0
/*
 *      uUnitCopy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a copy of the UNIT.
 *      This is the routine to use to copy the UNIT and have
 *      ALL internal pointers reset correctly.
 *
 *TODO: Copy RESTRAINTs
 */
UNIT
uUnitCopy( UNIT uUnit )
{
UNIT    uNewUnit;
LOOP    lContainerLoop;
OBJEKT  cCont;

    uNewUnit = (UNIT)oObjectDuplicate(uUnit);
    
    lContainerLoop = lLoop( (OBJEKT)uNewUnit, CONTAINERS );
    while ( ( cCont = oNext( &lContainerLoop ) ) != NULL ) {
        ContainerResetPointers(cCont);
    }

                /* Reset the pointers for the UNIT itself */

    UnitResetPointers( uNewUnit );

    return(uNewUnit);
}
#endif



/*
 *      UnitJoinTailHead
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Join the two UNITS together and create a bond
 *      between the Tail of the first and the Head of the last.
 *      The result of the joining is placed in uA, and uB is
 *      Destroyed.
 *
 *TODO: Figure out what you want to do with RESTRAINTs
 */
static void
UnitJoinTailHead( UNIT uA, UNIT uB, BONDt bBondOrder )
{
LOOP            lContents;
OBJEKT          oObj;
ATOM            aA, aB, aNext;
RESIDUE         rRes;
int             iPdbSeq;

    if ( !bUnitTailUsed(uA) ) { 
        VP0(( "UNIT: %s does not have a tail atom.\n",
                sContainerName(uA) ));
        goto DONE;
    } else aA = aUnitTail(uA);
    if ( !bUnitHeadUsed(uB) ) {
        VP0(( "UNIT: %s does not have a head atom.\n",
                sContainerName(uB) ));
        goto DONE;
    } else aB = aUnitHead(uB);

                /* Get the tail atom on uB */
    if ( !bUnitTailUsed(uB) ) aNext = NULL;
    else aNext = aUnitTail(uB);

                /* Get the highest PDB sequence number */
    iPdbSeq = 0;
    lContents = lLoop( (OBJEKT)uA, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lContents)) ) {
        if ( iPdbSeq < iResiduePdbSequence(rRes) ) {
            iPdbSeq = iResiduePdbSequence(rRes);
        }
    }
                /* Remove all of the objects from uB and put them in uA */
                /* This COMPLETELY empties the second UNIT, making it */
                /* have NULL connect atoms */

    lContents = lLoop( (OBJEKT)uB, DIRECTCONTENTSBYSEQNUM );
    while ( (oObj = oNext(&lContents)) ) {
        REF( oObj );    /* bContainerRemove() does a DEREF */
        bContainerRemove( (CONTAINER)uB, oObj );

                /* If the object being added is a RESIDUE then update */
                /* it's PDB sequence number */

        if ( iObjectType(oObj) == RESIDUEid ) {
            iPdbSeq++;
            ResidueSetPdbSequence( oObj, iPdbSeq );
        }
        ContainerAdd( (CONTAINER) uA, oObj );
        DEREF( oObj );  /* ContainerAdd() does a REF */
    }

                /* Now connect the atoms together */

    AtomBondToOrder( aA, aB, bBondOrder );
    UnitSetTail( uA, aNext );

DONE:
                /* Destroy UNIT uB */

    DisplayerUpdate(dContainerDisplayer(uA));
    Destroy((OBJEKT *) &uB );
    CDU(uA);
}
 
    

/*
 *      UnitJoin
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Join two units together without connecting atoms.
 *      The result of the joining is placed in uA, and uB is
 *      Destroyed.
 *TODO: Figure out what you want to do with RESTRAINTs
 */
void
UnitJoin( UNIT uA, UNIT uB )
{
LOOP            lContents;
OBJEKT          oObj;
ATOM            aB;
RESIDUE         rRes;
int             iPdbSeq;

        /* find max res seq of 1st unit */

    iPdbSeq = 0;
    lContents = lLoop( (OBJEKT)uA, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lContents)) ) {
        if ( iPdbSeq < iResiduePdbSequence(rRes) )
            iPdbSeq = iResiduePdbSequence(rRes);
    }

    aB = aUnitTail(uB);

                /* Remove all of the objects from uB and put them in uA */

    lContents = lLoop( (OBJEKT)uB, DIRECTCONTENTSBYSEQNUM );
    while ( (oObj = oNext(&lContents)) ) {
        REF( oObj );    /* bContainerRemove() does a DEREF */
        bContainerRemove( (CONTAINER)uB, oObj );

                /* If the object being added is a RESIDUE then update */
                /* its PDB sequence number */

        if ( iObjectType(oObj) == RESIDUEid ) {
            iPdbSeq++;
            ResidueSetPdbSequence( oObj, iPdbSeq );
        }
        ContainerAdd( (CONTAINER) uA, oObj );
        DEREF( oObj );  /* ContainerAdd() does a REF */
    }

                /* If the connect1 atom of the UNIT uB is defined then */
                /* set up the new UNIT uA to have the same connect1 atom */

    UnitSetTail( uA, aB );

    DisplayerUpdate(dContainerDisplayer(uA));
    Destroy( (OBJEKT *)&uB );
    CDU(uA);
}




/*
 *      UnitSequence
 *
 *      Author: Christian Schafmeister (1991)
 * 
 *      Join the two UNITs together into a sequence, meaning
 *      if they both have connect atoms defined for LASTEND and FIRSTEND
 *      then connect them.  If only one has its connect atom defined
 *      then print a warning and join them, 
 *      and if neither has its connect atom defined
 *      then simply join the UNITs. 
 *TODO: Figure out what you want to do with RESTRAINTs
 */
void
UnitSequence( UNIT uA, UNIT uB )
{
BOOL            bA, bB;

    bA = bUnitTailUsed(uA);
    bB = bUnitHeadUsed(uB);

    if ( bA && bB ) {
        VP2(( "Joining %s - %s\n", 
                sContainerName(cContainerWithin(aUnitTail(uA))),
                sContainerName(cContainerWithin(aUnitHead(uB)))
                ));
        UnitJoinTailHead( uA, uB, BONDSINGLE );
    } else if ( !bA && !bB ) {
        VP1(( "Starting new chain with %s\n", sContainerName(uB ) ));
        UnitJoin( uA, uB );
    } else if ( !bA ) {
        VP0(( "One sided connection. Residue: %s missing connect%d atom.\n",
                sContainerName(uA), LASTEND ));
        UnitJoin( uA, uB );
    } else {
        VP0(( "One sided connection. Residue: %s missing connect%d atom.\n",
                sContainerName(uB), FIRSTEND ));
/* %%% */
        UnitJoin( uA, uB );
    }
}
         


/*
 *      UnitSave
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save a description of the UNIT to a DATABASE.
 */
void
UnitSave( UNIT uUnit, DATABASE db, PARMLIB plParameters )
{
BOOL            bGeneratedParameters;

    zUnitIOBuildTables( uUnit, plParameters, &bGeneratedParameters, 
                        TRUE, FALSE );
    zUnitIOSaveTables( uUnit, db );
    zUnitIODestroyTables( uUnit );
}





/*
 *      uUnitLoad
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load a description of the UNIT from the DATABASE.
 */
UNIT
uUnitLoad( DATABASE db )
{
UNIT            uNew;

    uNew = (UNIT)oCreate(UNITid);
    if ( !zbUnitIOLoadTables( uNew, db ) ) {
        Destroy((OBJEKT *)&uNew);
        uNew = NULL;
        goto DONE;
    }
    zUnitIOBuildFromTables( uNew );
    zUnitIODestroyTables( uNew );
DONE:
    return(uNew);
}



/*
 *      UnitSaveAmberParmFile
 *
 *      Author: Christian Schafmeister (1991), Robin Betz (2013)
 *
 *      Write the UNIT to an AMBER parmfile.
 *
 *      Arguments:
 *              uUnit   - UNIT to save.
 *              fOut    - FILE to write to.
 *              fCrd    - FILE to write coordinates to.
 *              plParms - Parameter library from which to obtain parameters.
 *              bPert   - TRUE means write a perturbation file to be run by GIBBS.
 *              bNcdf   - TRUE means write the coordinates in NetCDF format
 */
void
UnitSaveAmberParmFile( UNIT uUnit, FILE *fOut, char *crdName, 
        PARMLIB plParms, BOOL bPolar, BOOL bPert, BOOL bNetcdf )
{
        BOOL            bGeneratedParameters;

        zUnitIOBuildTables( uUnit, plParms, &bGeneratedParameters, bPert, TRUE );
        if ( bGeneratedParameters == TRUE ) {
                if( GDefaults.iOldPrmtopFormat ) 
                        zUnitIOSaveAmberParmFormat_old( uUnit, fOut, crdName, 
                                                                bPolar, bPert );
                else
                        zUnitIOSaveAmberParmFormat( uUnit, fOut, crdName, 
                                                       bPolar, bPert, bNetcdf );
        } else {
                VP0(( "Parameter file was not saved.\n" ));
        }
        zUnitIODestroyTables( uUnit );
}

/*
 *      UnitYouAreBeingRemoved
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      UNITs should never be removed from anything.
 */
void
UnitYouAreBeingRemoved( UNIT uUnit )
{
    DFATAL(( "UNITs should never be removed from anything!" ));
}



/*
 *      UnitIAmBeingRemoved
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      If the thing being removed is an ATOM then check to
 *      make sure that it is not one of the UNITs connection
 *      ATOMs, or that it is not in a restraint bond, angle,
 *      or torsion.
 *      Or that it is not in an ATOM group.
 *      It is not interesting if RESIDUEs or MOLECULEs
 *      are being removed.
 */
void
UnitIAmBeingRemoved( UNIT uUnit, CONTAINER cRemoved )
{
BAGLOOP         blRestraints;
RESTRAINT       rRest;
BOOL            bFoundOne;
DICTLOOP        dlGroups;
LIST            lAtoms;

        /* If it is an ATOM being removed then... */

    if ( iObjectType(cRemoved) == ATOMid ) {

                /* Check the connection atoms */

        if ( uUnit->aHead == (OBJEKT)cRemoved ) {
            uUnit->aHead = NULL;
        }
        if ( uUnit->aTail == (OBJEKT)cRemoved ) {
            uUnit->aTail = NULL;
        }

                /* Check the restraint bonds/angles/torsions */
        bFoundOne = FALSE;
        blRestraints = blBagLoop(uUnit->bRestraints);
        while ( (rRest = (RESTRAINT)PBagNext(&blRestraints)) ) {
            if ( bRestraintContainsAtom( rRest, (ATOM)cRemoved ) ) {
                bFoundOne = TRUE;
                bBagRemove( uUnit->bRestraints, (GENP)rRest );
            }
        }
        if ( bFoundOne ) {
            VP1(( "Removed all restraints that contained the atom.\n" ));
        }

                /* Check the ATOM groups, remove (cRemoved) from */
                /* any ATOM groups that it may be in */
                
        dlGroups = ydlDictionaryLoop( uUnit->dAtomGroups );
        while ( (lAtoms = (LIST)yPDictionaryNext(uUnit->dAtomGroups,&dlGroups))){
            bListRemove( lAtoms, (GENP)cRemoved );
        }
    }
}




/*
 *      bUnitCanBePerturbed
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Check all the ATOMs within the UNIT to see if
 *      any of them can be perturbed.  If any of them can be
 *      then the UNIT can be perturbed.
 */
BOOL
bUnitCanBePerturbed( UNIT uUnit )
{
LOOP            lAtoms;
ATOM            aAtom;
BOOL            bCanBePerturbed;

    bCanBePerturbed = FALSE;
    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        bCanBePerturbed |= bAtomPerturbed(aAtom);
    }
    return(bCanBePerturbed);
}

        



/*
 *================================================================
 *
 *      RESTRAINT code
 *
 */


/*
 *      UnitAddRestraint
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add a restraint to the UNIT.
 *
 *TODO: Remember to redirect pointers and zero pointers when 
 *TODO: UNITs are copied, or atoms are removed.
 *
 */
void
UnitAddRestraint( UNIT uUnit, RESTRAINT rRest )
{
    BagAdd( uUnit->bRestraints, (GENP)rRest );
    CDU(uUnit);
}




/*
 *      bUnitRemoveRestraint
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove the RESTRAINT, return TRUE if it is found.
 */
BOOL
bUnitRemoveRestraint( UNIT uUnit, RESTRAINT rRest )
{
BOOL            bReturn;

    bReturn = bBagRemove( uUnit->bRestraints, (GENP)rRest );
    CDU(uUnit);    
    return(bReturn);
}





/*
 *      UnitLoopRestraints
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Initialize the one loop that is allowed at a time
 *      over the RESTRAINTs in the UNIT.
 */
void
UnitLoopRestraints( UNIT uUnit )
{
    uUnit->blRestraintLoop = blBagLoop(uUnit->bRestraints);
}



/*
 *      rUnitNextRestraint
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the next RESTRAINT in the RESTRAINT loop.
 */
RESTRAINT
rUnitNextRestraint( UNIT uUnit )
{
    return((RESTRAINT)PBagNext(&(uUnit->blRestraintLoop)));
}


/*
 *      iUnitRestraintTypeCount
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Count the number of restraints of a certian type.
 */
int
iUnitRestraintTypeCount( UNIT uUnit, int iType )
{
BAGLOOP         blRestraints;
RESTRAINT       rRest;
int             iCount;

    iCount = 0;
    blRestraints = blBagLoop(uUnit->bRestraints);
    while ( (rRest = (RESTRAINT)PBagNext(&blRestraints)) ) {
        if ( iRestraintType(rRest) == iType ) iCount++;
    }
    return(iCount);
}


/*
 *      bUnitCapContainsAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return TRUE if the UNITs solvent cap contains
 *      the ATOM.
 */
BOOL
bUnitCapContainsAtom( UNIT uUnit, ATOM aAtom )
{
VECTOR          vDiff;
double          dDist2;

        /* If there is no cap then return FALSE always */

    if ( !bUnitUseSolventCap(uUnit) ) return(FALSE);
    vDiff = vVectorSub( &vAtomPosition(aAtom), &(uUnit->vCapOrigin) );

    dDist2 = dVX(&vDiff)*dVX(&vDiff) +
                dVY(&vDiff)*dVY(&vDiff) +
                dVZ(&vDiff)*dVZ(&vDiff);

    if ( dDist2 < (uUnit->dCapRadius)*(uUnit->dCapRadius) ) return(TRUE);
    return(FALSE);
}

    

/*
 *      bUnitCapContainsContainer
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return TRUE if the UNITs solvent cap contains
 *      the at least one ATOM within the CONTAINER
 */
BOOL
bUnitCapContainsContainer( UNIT uUnit, CONTAINER cCont )
{
LOOP            lAtoms;
ATOM            aAtom;

        /* If there is no cap then return FALSE always */


    if ( !bUnitUseSolventCap(uUnit) ) return(FALSE);

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( bUnitCapContainsAtom( uUnit, aAtom ) ) return(TRUE);
    }
    return(FALSE);
}






/*
 *=====================================================================
 *
 *      Change UNIT attributes.
 *
 */


/*
 *      zUnitSetHead
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Make sure that the ATOM is within the UNIT and
 *      then set the connect ATOM, otherwise print an error.
 */
static void
zUnitSetHead( UNIT uUnit, ATOM aAtom )
{
STRING  sTemp1;

    if ( !((aAtom == NULL) || bObjektWarnType( (OBJEKT)aAtom, ATOMid )) ) 
        return;
    if ( aAtom != NULL ) {
        if ( !bContainerContainedBy( (CONTAINER)aAtom, (CONTAINER)uUnit ) ) {
            VP0(( "The UNIT must contain %s.\n",
                sContainerFullDescriptor( (CONTAINER)aAtom, sTemp1 ) ));
            return;
        }
    }
    UnitSetHead( uUnit, aAtom );
}


/*
 *      zUnitSetTail
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Make sure that the ATOM is within the UNIT and
 *      then set the connect ATOM, otherwise print an error.
 */
static void
zUnitSetTail( UNIT uUnit, ATOM aAtom )
{
STRING  sTemp1;

    if ( !((aAtom == NULL) || bObjektWarnType( (OBJEKT)aAtom, ATOMid )) ) 
        return;
    if ( aAtom != NULL ) {
        if ( !bContainerContainedBy( (CONTAINER)aAtom, (CONTAINER)uUnit ) ) {
            VP0(( "The UNIT must contain %s.\n",
                sContainerFullDescriptor( (CONTAINER)aAtom, sTemp1 ) ));
            return;
        }
    }
    UnitSetTail( uUnit, aAtom );
}



/*
 *      zUnitSetResidueType
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the RESIDUE type of all RESIDUEs within the UNIT.
 */
static void
zUnitSetResidueType( UNIT uUnit, OBJEKT oObj )
{
LOOP            lResidues;
RESIDUE         rRes;

    lResidues = lLoop( (OBJEKT)uUnit, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lResidues)) ) {
        ContainerSetAttribute( (CONTAINER)rRes, "restype", oObj );
    }
}




/*
 *      zUnitSetBox
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the bounding box of the UNIT to the
 *      dimensions passed in oAttr.
 *      If oAttr is a single double then use that for all three dimensions.
 *      If oAttr use the first, second, and third values as x,y,z.
 *      If oAttr is NULL then don't use the bounding box.
 */
static void
zUnitSetBox( UNIT uUnit, OBJEKT oAttr )
{
double          dX, daVector[3];
int             i;
LISTLOOP        llElements;
ASSOC           aAssoc;


    switch ( iObjectType(oAttr) ) {
        case ODOUBLEid:
            dX = dODouble(oAttr);
            UnitSetBox( uUnit, dX, dX, dX );
            UnitSetBeta( uUnit, 90.0*DEGTORAD );
            UnitSetUseBox( uUnit, TRUE );
            break;
        case LISTid:
            i = 0;
            llElements = llListLoop((LIST)oAttr);
            while ( (aAssoc = (ASSOC)oListNext(&llElements)) ) {
                if ( iObjectType(oAssocObject(aAssoc)) == ODOUBLEid ) {
                    if ( i<3 ) {
                        daVector[i++] = dODouble(oAssocObject(aAssoc));
                    } else {
                        VP0(( "Illegal vector\n" ));
                        break;
                    }
                } else {
                    VP0(( "Illegal vector\n" ));
                    break;
                }
            }
            if ( i==3 ) {
                UnitSetBox( uUnit, daVector[0], daVector[1], daVector[2] );
                UnitSetUseBox( uUnit, TRUE );
                UnitSetBeta( uUnit, 90.0*DEGTORAD );
            }
            break;
        case NULLid:
            UnitSetUseBox( uUnit, FALSE );
            break;
        default:
            VP0(( "Illegal box definition.\n" ));
            break;
    }
}


/*
 *      zUnitSetCap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the cap of the UNIT to the
 *      dimensions passed in oAttr.
 *      oAttr must be a list with four values, x,y,z of the origin, and
 *      the radius, or it can be NULL if the cap is to be turned off.
 */
static void
zUnitSetCap( UNIT uUnit, OBJEKT oAttr )
{
double          daVector[4];
int             i;
LISTLOOP        llElements;
ASSOC           aAssoc;

    switch ( iObjectType(oAttr) ) {
        case LISTid:
            i = 0;
            llElements = llListLoop((LIST)oAttr);
            while ( (aAssoc = (ASSOC)oListNext(&llElements)) ) {
                if ( iObjectType(oAssocObject(aAssoc)) == ODOUBLEid ) {
                    if ( i<4 ) {
                        daVector[i++] = dODouble(oAssocObject(aAssoc));
                    } else {
                        VP0(( "Illegal vector\n" ));
                        break;
                    }
                } else {
                    VP0(( "Illegal vector\n" ));
                    break;
                }
            }
            if ( i==4 ) {
                UnitSetSolventCap( uUnit, 
                                daVector[0], daVector[1], daVector[2],
                                daVector[3] );
                UnitSetUseSolventCap( uUnit, TRUE );
            }
            break;
        case NULLid:
            UnitSetUseSolventCap( uUnit, FALSE );
            break;
        default:
            VP0(( "Illegal box definition.\n" ));
            break;
    }
}







/*
 *      UnitSetAttribute
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Change the attribute of the UNIT.
 *      This function is used mainly by the command line interface
 *      to allow the user to change attributes of UNITs by
 *      name.
 */
void
UnitSetAttribute( UNIT uUnit, STRING sAttr, OBJEKT oAttr )
{

    if ( strcmp( sAttr, "head" ) == 0 ) {
        zUnitSetHead( uUnit, (ATOM)oAttr );
        goto DONE;
    } else if ( strcmp( sAttr, "tail" ) == 0 ) {
        zUnitSetTail( uUnit, (ATOM)oAttr );
        goto DONE;
    } else if ( strcmp( sAttr, "restype" ) == 0 ) {
        zUnitSetResidueType( uUnit, oAttr );
        goto DONE;
    } else if ( strcmp( sAttr, "box" ) == 0 ) {
        zUnitSetBox( uUnit, oAttr );
        goto DONE;
    } else if ( strcmp( sAttr, "cap" ) == 0 ) {
        zUnitSetCap( uUnit, oAttr );
        goto DONE;
    } else if ( strcmp( sAttr, "name" ) == 0 ) {
        if ( !bObjektWarnType( oAttr, OSTRINGid ) ) return;
        ContainerSetName( uUnit, sOString(oAttr) );
        goto DONE;
    }

    VP0(( "%s: %s is a non-existent attribute for a unit.\n", 
                                sContainerName(uUnit), sAttr ));
    VP0(( "\tUnit attributes: head, tail, restype, box, cap, name\n" ));
    if ( strncmp( sAttr, "connect", 7 ) == 0 )
        VP0(( "\t-- 'connectX' is used for residues, e.g. %s.1\n",
                                                sContainerName(uUnit) ));
    return;

DONE:
    CDU(uUnit);

}





/*
 *      bUnitGroupCreate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a new ATOM group within the UNIT.
 *      Return TRUE if the group was created.
 */
BOOL
bUnitGroupCreate( UNIT uUnit, char *cPName )
{
LIST    lAtoms;

    if ( yPDictionaryFind( uUnit->dAtomGroups, cPName ) ) {
        return(FALSE);
    }
    lAtoms = (LIST)oCreate(LISTid);
    DictionaryAdd( uUnit->dAtomGroups, cPName, (GENP)lAtoms );
    return(TRUE);
}



/*
 *      lUnitGroup
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the LIST for the group.
 */
LIST
lUnitGroup( UNIT uUnit, char *sGroup )
{
LIST    lGroup;

    lGroup = (LIST)yPDictionaryFind( uUnit->dAtomGroups, sGroup );
    return(lGroup);
}




/*
 *      bUnitGroupAddAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add an ATOM to the group.
 *      Return FALSE if the group does not exist.
 */
BOOL
bUnitGroupAddAtom( UNIT uUnit, char *sGroup, ATOM aAtom )
{
LIST    lGroup;

    lGroup = (LIST)yPDictionaryFind( uUnit->dAtomGroups, sGroup );
    if ( lGroup == NULL ) return(FALSE);
    ListAddUnique( lGroup, (GENP)aAtom );
    return(TRUE);
}




/*
 *      bUnitGroupFindAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set bFound to TRUE if the ATOM was found, otherwise FALSE.
 *      Return TRUE if the group exists.
 */
BOOL
bUnitGroupFindAtom( UNIT uUnit, char *sGroup, ATOM aAtom, BOOL *bPFound )
{
LIST    lAtoms;

    lAtoms = (LIST)yPDictionaryFind( uUnit->dAtomGroups, sGroup );
    if ( lAtoms == NULL ) return(FALSE);

    if ( bListContains( lAtoms, (GENP)aAtom ) ) 
        *bPFound = TRUE;
    else                                  *bPFound = FALSE;
    
    return(TRUE);
}



/*
 *      bUnitGroupRemoveAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove the ATOM from the group.
 *      Return FALSE if the group does not exist.
 */
BOOL
bUnitGroupRemoveAtom( UNIT uUnit, char *sGroup, ATOM aAtom )
{
LIST    lAtoms;

    lAtoms = (LIST)yPDictionaryFind( uUnit->dAtomGroups, sGroup );
    if ( lAtoms == NULL ) return(FALSE);

    bListRemove( lAtoms, (GENP)aAtom );
    return(TRUE);
}




/*
 *      bUnitGroupDestroy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Destroy the group.
 *      Return TRUE if the group existed, otherwise FALSE.
 */
BOOL
bUnitGroupDestroy( UNIT uUnit, char *sGroup )
{
LIST    lAtoms;

    lAtoms = (LIST)yPDictionaryDelete( uUnit->dAtomGroups, sGroup );
    if ( lAtoms == NULL ) return(FALSE);
    Destroy((OBJEKT *) &lAtoms );
    return(TRUE);
}
    
 
/*
 *      UnitCheck
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Run checks on the UNIT and its contents to see if the are
 *      ready to have calculations run on them.
 *
 *TODO: Figure out what you want to do with RESTRAINTs
 *
 *      Arguments:
 *              uUnit - The UNIT to check.
 *              iPErrors -      Add the number of errors found.
 *              iPWarnings -    Add the number of warnings found.
 */
void
UnitCheck( UNIT uUnit, int *iPErrors, int *iPWarnings )
{
LOOP            lContents;
CONTAINER       cCont;
LOOP            lBonds, lAtoms;
ATOM            aAtom1, aAtom2;
double          dLen;
STRING          sTemp1, sTemp2;
VECTOR          vVector;
double          dCharge, dPertCharge, dFrac, dAbs;
BOOL            bPert;


        /* Check to make sure that all bond lengths are 2.0 angstroms */
        /* or less, flag those that are as warnings */

    lBonds = lLoop( (OBJEKT)uUnit, BONDS );
    while ( oNext(&lBonds) ) {
        LoopGetBond( &lBonds, &aAtom1, &aAtom2 );
        vVector = vVectorSub( &vAtomPosition(aAtom1),
                              &vAtomPosition(aAtom2) );
        dLen = dVectorLen(&vVector);
        if ( dLen > MAXBONDLEN || dLen < MINBONDLEN ) {
            (*iPWarnings)++; 
	    if ( strcmp( sAtomName( aAtom1 ), "EPW" ) == 0 ||
		 strcmp( sAtomName( aAtom2 ), "EPW" ) == 0 ) {
	      /* Lone pair, ignore it  */
	    } else {
            VP0(( "WARNING: There is a bond of %lf angstroms between: \n",
                        dLen ));
            VP0(( "-------  %s and %s\n",
                        sContainerFullDescriptor( (CONTAINER)aAtom1, sTemp1 ),
                        sContainerFullDescriptor( (CONTAINER)aAtom2, sTemp2 ) ));
	    }
        }
    }

    /*
     *  note whether any perturbation is involved
     */

    bPert = FALSE;
    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    while ( (aAtom1 = (ATOM)oNext(&lAtoms)) ) {
        if ( bAtomPerturbed( aAtom1 ) ) {
            bPert = TRUE;
            break;
        }
    }
                /* Now check the charge of the UNIT */
                /* Flag an error if the charge is not integral */
                /* Flag a warning if the charge is not zero */

    ContainerTotalCharge( (CONTAINER)uUnit, &dCharge, &dPertCharge );

    dAbs = fabs(dCharge);
    dFrac = fabs( dAbs - (double)(int)(dAbs+0.5) );
    if ( dFrac > 0.01 ) {
        VP0(( 
          "ERROR: The unperturbed charge of the unit: %lf is not integral.\n",
          dCharge ));
        (*iPErrors)++;
    }
    if ( fabs(dCharge) > 0.01 ) {
        VP0(( "WARNING: The unperturbed charge of the unit: %lf is not zero.\n",
                dCharge ));
        (*iPWarnings)++;
    }
    if ( bPert == TRUE ) {
        dAbs = fabs(dCharge + dPertCharge);
        dFrac = fabs( dAbs - (double)(int)(dAbs+0.5) );
        if ( dFrac > 0.01 ) {
            VP0(( 
                "ERROR: The perturbed charge: %lf is not integral.\n",
                (dCharge+dPertCharge) ));
            (*iPErrors)++;
        }

        if ( fabs(dCharge+dPertCharge) > 0.01 ) {
            VP0(( 
                "WARNING: The perturbed charge: %lf is not zero.\n",
                (dCharge+dPertCharge) ));
            (*iPWarnings)++;
        }
    }
                /* Check the contents */

    lContents = lLoop( (OBJEKT)uUnit, DIRECTCONTENTS );
    while ( (cCont = (CONTAINER)oNext(&lContents)) ) {
        ContainerCheck( cCont, iPErrors, iPWarnings );
    }
}



/*
 *      UnitCheckForParms
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Check the unit or having all necessary parameters in the 
 *      parmLib given or in the units parmSet.
 *      If not, copy the parameters into the parmSet given if one is given.
 *
 */

void
UnitCheckForParms( UNIT uUnit, PARMLIB plParms, PARMSET psParmSet )
{

    if ( zbUnitParmsMissing( uUnit, plParms) == TRUE ) {
    
        VP0(( "There are missing parameters.\n" ));
        if ( psParmSet ) {
            VP0(( "Missing parameters have been added to the PARMSET.\n" ));
        }
    }
}
           



/*
 *      UnitFindBoundingBox
 *
 *      Author: David A. Rivkin (1992)
 *
 *      Determines the extent of the unit.      
 */
void
UnitFindBoundingBox( UNIT uUnit, VECTOR *vPLower, VECTOR *vPUpper )
{
LOOP    lAtoms;
ATOM    aAtom;
double  dXMin, dYMin, dZMin;
double  dXMax, dYMax, dZMax;


    lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
    dXMin = dYMin = dZMin = 0;
    dXMax = dYMax = dZMax = 0;
    while ( (aAtom = (ATOM)oNext( &lAtoms ))) {
        if ( dXMin > dVX( &(vAtomPosition( aAtom )))) {
            dXMin = dVX( &(vAtomPosition( aAtom )));
            MESSAGE(( "Min X Value:  %4.2lf\n", dXMin ));
        }
        if ( dYMin > dVY( &(vAtomPosition( aAtom )))) {
            dYMin = dVY( &(vAtomPosition( aAtom )));
            MESSAGE(( "Min Y Value:  %4.2lf\n", dYMin ));
        }
        if ( dZMin > dVZ( &(vAtomPosition( aAtom )))) {
            dZMin = dVZ( &(vAtomPosition( aAtom )));
            MESSAGE(( "Min Z Value:  %4.2lf\n", dZMin ));
        }
        if ( dXMax < dVX( &(vAtomPosition( aAtom )))) {
           dXMax = dVX( &(vAtomPosition( aAtom )));
           MESSAGE(( "Max X Value:  %4.2lf\n", dXMax ));
        }
        if ( dYMax < dVY( &(vAtomPosition( aAtom )))) {
            dYMax = dVY( &(vAtomPosition( aAtom )));
            MESSAGE(( "Max Y Value:  %4.2lf\n", dYMax ));
        }
        if ( dZMax < dVZ( &(vAtomPosition( aAtom )))) {
            dZMax = dVZ( &(vAtomPosition( aAtom )));
            MESSAGE(( "Max Z Value:  %4.2lf\n", dZMax ));
        }
    }
    VectorDef( vPLower, dXMin, dYMin, dZMin );
    VectorDef( vPUpper, dXMax, dYMax, dZMax );

}


