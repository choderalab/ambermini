/*
 *      File: atom.c
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
 *      Class:
 *              ATOM
 *      Superclass:
 *              CONTAINER, LOOP
 *
 *      Description:
 *              An ATOM is a subclass of CONTAINER.
 *              ATOMS do not contain any other container objects.
 *              Atoms can have connections to other atoms.
 *
 *              ATOMs also have a pointer that can be used
 *              to point to a structure that stores graphics
 *              information in 3D environments like X-Windows.
 *              This can be initialised whenever an ATOM
 *              is created by defining a creator and destructor
 *              using AtomClassDefineGraphicsCreator and
 *              AtomClassDefineGraphicsDestructor.
 *
 */
 



#include        "basics.h"
#include        "defaults.h"

#include        "classes.h"

#include        "elements.h"

extern int      iFatal;

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        Static variables
*/
static  int     SiUniqueId = 0;

typedef struct  {
                STRING          sName;
                int             iOrder;
                } ORDERt;

static  ORDERt  SoaOrders[] = {
        { "-",          BONDSINGLE },
        { "single",     BONDSINGLE },
        { "=",          BONDDOUBLE },
        { "double",     BONDDOUBLE },
        { "#",          BONDTRIPLE },
        { "triple",     BONDTRIPLE },
        { ":",          BONDAROMATIC },
        { "aromatic",   BONDAROMATIC } };

/*
 *  The Graphics function pointers are set in xaLeap.c main(),
 *      otherwise they remain null for the non-graphical tleap
 *      in which case the graphics routines aren't linked.
 */

VFUNCTION       GfAtomClassGraphicsCreator = NULL;
VFUNCTION       GfAtomClassGraphicsDestructor = NULL;

VARARRAY        GvaVDWTypes, GvaVDWValues;

/*
======================================================================

        Private messages
*/

/*
 *      zAtomRemoveBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Find and remove the bond from aAtom1 to aAtom2
 *      This routine must then be called to remove the bond from aAtom2
 *      to aAtom1.
 */
void
zAtomRemoveBond( ATOM aAtom1, ATOM aAtom2 )
{
int     i, j;

    for ( i=0; i<iAtomCoordination(aAtom1); i++ )
        if ( aAtom1->aaBonds[i] == aAtom2 ) break;
    if ( i == iAtomCoordination(aAtom1) ) 
        DFATAL( ("The atom does not contain that bond") );
    for ( j=i; j<iAtomCoordination(aAtom1)-1; j++ ) {
        aAtom1->aaBonds[j] = aAtom1->aaBonds[j+1];
        aAtom1->faBondFlags[j] = aAtom1->faBondFlags[j+1];
    }
    aAtom1->aaBonds[iAtomCoordination(aAtom1)-1] = NULL;
    aAtom1->iCoordination--;
}




/*
----------------------------------------------------------------------

        PUBLIC messages
        
*/


/*
 *      aAtomCreate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a molecule.
 *
 *      Return:
 *              Return a pointer to the ATOM created.
 */
ATOM 
aAtomCreate()
{
ATOM    a;
int     i;

    MALLOC( a, ATOM, sizeof(ATOMt) );
    a->iCoordination = 0;
    for ( i=0; i<MAXBONDS; i++ ) 
        a->aaBonds[i] = NULL;
    strcpy( a->sType, "" );
    strcpy( a->sPertType, "" );
    strcpy( a->sPertName, "" );
    strcpy( a->siSegid, "" );
    a->iUniqueId = SiUniqueId++;
    a->fFlags = 0;
    a->dCharge = 0.0;
    a->dPertCharge = 0.0;
    a->dPolar = 0.0;
    a->dPertPolar = 0.0;
    a->dScreenF = 0.0;
    a->iAtomicNumber = NOELEMENT;
    a->iPertAtomicNumber = NOELEMENT;
    a->iSeenId = 0;
    VectorDef( &(a->vPosition), 0.0, 0.0, 0.0 );
    VectorDef( &(a->vVelocity), 0.0, 0.0, 0.0 );
    a->PGraphicsData = NULL;

    if ( GfAtomClassGraphicsCreator != NULL ) {
        GfAtomClassGraphicsCreator(a);
    }

    return(a);
}







/*
 *      AtomDestroy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Destroy the ATOM and all of the residues/atoms within it.
 *
 *      Arguments:
 *      oObject -       A pointer to the object.
 */
void
AtomDestroy( ATOM *aPAtom )
{
    if ( GfAtomClassGraphicsDestructor != NULL ) {
        GfAtomClassGraphicsDestructor(*aPAtom);
    }
    FREE( *aPAtom );
    *aPAtom = NULL;
}






/*
 *      AtomDescribe
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Print a description of the atom to file fOut.
 *
 *      Arguments:
 *      oObject -       The object to describe.
 *
 */
void
AtomDescribe( ATOM aAtom )
{
int             i;
STRING          sTemp, sTemp1;

    if ( strlen(sAtomPertType(aAtom)) != 0 ||
         strlen(sAtomPertName(aAtom)) != 0 ) {
        VP0(( "ATOM\n" ));
        VP0(( "             Normal      Perturbed\n" ));
        VP0(( "Name:         %-5s         %-5s\n",
                sContainerName(aAtom), sAtomPertName(aAtom) ));
        VP0(( "Type:         %-5s         %-5s\n",
                sAtomType(aAtom), sAtomPertType(aAtom) ));
        VP0(( "Charge:       %6.4f        %5.3f\n",
                dAtomCharge(aAtom), dAtomPertCharge(aAtom) ));
        VP0(( "Polarization: %6.4f        %5.3f\n",
                dAtomPolar(aAtom), dAtomPertPolar(aAtom) ));
        VP0(( "Element: %-5s         (not affected by pert)\n",
                sElementName( iAtomElement(aAtom), sTemp1 ) ));
    } else {    
        VP0(( "ATOM\n" ));
        VP0(( "Name:         %-5s\n", sContainerName(aAtom) ));
        VP0(( "Type:         %-5s\n", sAtomType(aAtom) ));
        VP0(( "Charge:       %6.4f\n", dAtomCharge(aAtom) ));
        VP0(( "Polarization: %6.4f\n", dAtomPolar(aAtom) ));
        VP0(( "Element:      %-5s\n",
                sElementName( iAtomElement(aAtom), sTemp1 ) ));
    }

    VP0(( "Atom flags: (decimal %ld hex 0x%lx)\n", 
                                fAtomFlags(aAtom), fAtomFlags(aAtom) ));
    VP0(( "\tposfxd %c  posblt %c  posdrwn %c  selected %c\n",
        ( bAtomFlagsSet( aAtom, ATOMPOSITIONFIXED ) ? 'Y' : 'n' ),
        ( bAtomFlagsSet( aAtom, ATOMPOSITIONBUILT ) ? 'Y' : 'n' ),
        ( bAtomFlagsSet( aAtom, ATOMPOSITIONDRAWN ) ? 'Y' : 'n' ),
        ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ? 'Y' : 'n' ) ));
    VP0(( "\tpert %c  notdisp %c  touched %c  posknwn %c\n",
        ( bAtomFlagsSet( aAtom, ATOMPERTURB ) ? 'Y' : 'n' ),
        ( bAtomFlagsSet( aAtom, ATOMNOTDISPLAYED ) ? 'Y' : 'n' ),
        ( bAtomFlagsSet( aAtom, ATOMTOUCHED ) ? 'Y' : 'n' ),
        ( bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) ? 'Y' : 'n' ) ));
    VP0(( "\tinternal %c  needsmin %c  needsbuild %c\n",
        ( bAtomFlagsSet( aAtom, ATOMPOSITIONINTERNAL ) ? 'Y' : 'n' ),
        ( bAtomFlagsSet( aAtom, ATOMNEEDSMINIMIZER ) ? 'Y' : 'n' ),
        ( bAtomFlagsSet( aAtom, ATOMNEEDSBUILD ) ? 'Y' : 'n' ) ));

    VP0(( "Atom position: %lf, %lf, %lf\n", 
                dVX(&vAtomPosition(aAtom)),
                dVY(&vAtomPosition(aAtom)),
                dVZ(&vAtomPosition(aAtom)) ));
    VP0(( "Atom velocity: %lf, %lf, %lf\n", 
                dVX(&vAtomVelocity(aAtom)),
                dVY(&vAtomVelocity(aAtom)),
                dVZ(&vAtomVelocity(aAtom)) ));
    if ( iAtomCoordination(aAtom) == 0 ) 
        VP0(( "  NO BONDS\n" ));
    for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
        VP0(( "  Bonded to %s by a", 
             sContainerFullDescriptor((CONTAINER)aAtomBondedNeighbor(aAtom,i), sTemp) ));
        switch ( iAtomBondOrder(aAtom,i) ) {
            case BONDSINGLE:
                VP0(( " single" ));
                break;
            case BONDDOUBLE:
                VP0(( " double" ));
                break;
            case BONDTRIPLE:
                VP0(( " triple" ));
                break;
            case BONDAROMATIC:
                VP0(( "n aromatic" ));
                break;
            default:
                VP0(( " ????" ));
                break;
        }
        VP0(( " bond.\n" ));
    }
    if ( iListSize(aAtom->cHeader.lContents) != 0 ) {
        VP0(( "Number of internals = %d\n",
                iListSize(aAtom->cHeader.lContents) ));
    }
}


void
AtomDescStr( ATOM aA, BOOL bResNum, char *cPDesc )
{
CONTAINER       cTemp;
STRING          sTemp;
char            *cP;

        cTemp = cContainerWithin( aA );
        if ( cTemp == NULL )
                DFATAL(( " atom container is null\n" ));
        sContainerDescriptor( cTemp, sTemp );
        cP = strchr( sTemp+2, ' ' );
        if ( bResNum ) {
                *cP = '_';
                cP = strchr( cP, ' ' );
        }
        *cP = '\0';
        sprintf( cPDesc, "%s<%s>", sAtomName( aA ), sTemp+2 );
}


/*
 *-------------------------------------------------------------
 */
 
/*
 *      bAtomCoordinationSaturated
 *
 *      Author: Bill Ross (1993)
 *
 *      Simple routine to screen out the more egregious over-orders
 */

BOOL
bAtomCoordinationSaturated( ATOM aAtom )
{
        VERIFYOBJEKT( aAtom, ATOMid );

        /*
         *  perturbed atoms can have extra 'bonds' so allow the max
         */
        
        if ( bAtomFlagsSet( aAtom, ATOMPERTURB ) ) {
                if ( iAtomCoordination( aAtom ) >= MAXBONDS ) {
                        return(TRUE);
                } else {
                        return(FALSE);
                }
        }

        /*
         *  normal atoms have normal limits
         */

        switch ( iAtomElement( aAtom ) ) {
            case HYDROGEN:
                if ( iAtomCoordination( aAtom ) ) {
                        if ( !strncmp( aAtom->sType, "HW", ATOMTYPELEN)
                                        && iAtomCoordination( aAtom ) == 1 )
                                return(FALSE);
                        return(TRUE);
                }
                return(FALSE);
                break;
            case OXYGEN:
                if ( iAtomCoordination( aAtom ) >= 4)
                        return(TRUE);
                return(FALSE);
                break;
            default:
                if ( iAtomCoordination( aAtom ) >= MAXBONDS )
                        return(TRUE);
                return(FALSE);
                break;
        }

        return(FALSE); /* for lint */
}

/*
 *      bBondAtomProblem
 *
 *      Author: Bill Ross (1994)
 *
 *      Make sure both atoms are ok to bond
 */
BOOL
bBondAtomProblem( ATOM aAtom1, ATOM aAtom2 )
{
int     problem = 0;
STRING  sTemp;

    if ( bAtomCoordinationSaturated( aAtom1 ) ) {
        VP0(( "Bond: maximum coordination exceeded on %s\n",
                sContainerFullDescriptor( (CONTAINER)aAtom1, sTemp )) );
        problem = 1;
    }
    if ( bAtomCoordinationSaturated( aAtom2 ) ) {
        VP0(( "Bond: Maximum coordination exceeded on %s\n",
                sContainerFullDescriptor( (CONTAINER)aAtom2, sTemp )) );
        problem = 1;
    }
    if ( problem ) {
        VP0(( "      -- setting atoms pert=true overrides default limits\n" ));
        return( TRUE );
    }
    return( FALSE );
}

/*
 *      AtomBondTo
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a bond between two atoms, with bond order BONDSINGLE.
 *      The user is responsible for checking if there are already
 *      the maximum number of bonds on the atom.
 */
void
AtomBondTo( ATOM aAtom1, ATOM aAtom2 )
{

    VERIFYOBJEKT( aAtom1, ATOMid );
    VERIFYOBJEKT( aAtom2, ATOMid );

    if ( bBondAtomProblem( aAtom1, aAtom2 ) ) {
        STRING  sAtom1, sAtom2;
        sContainerFullDescriptor( (CONTAINER)aAtom1, sAtom1 );
        sContainerFullDescriptor( (CONTAINER)aAtom2, sAtom2 );
        VP0(( "ATOMS NOT BONDED: %s %s\n", sAtom1, sAtom2 ));
        DFATAL(("bondAtomProblem found\n"));
    }

    /*
     *  look up residues

    rR1 = cContainerWithin( aAtom1 );
    rR2 = cContainerWithin( aAtom2 );
    if ( rR1 != rR2 ) {
        int     i1, i2;
        STRING  sDesc1, sDesc2;

         *  see if the change affects residue numbering

        sContainerDescriptor( rR1, sDesc1 );
        sContainerDescriptor( rR2, sDesc2 );
        VP0(( "Connecting residues %s -- %s\n", sDesc1+1, sDesc2+1 ));
        VP0(( "Warning - leap can have problems with residue order %s\n",
                                        "and parm topology when this happens"));

TODO - if leap rewritten so that residue parents are
molecules, then it will be easy to check if molecules
are being joined and residues need to be renumbered
so the 'earlier' mol1 residue #s are < mol2 residues

        for ( i1=0; i1<MAXCONNECT; i1++ )
            if ( aResidueConnectAtom( rR1, i1 ) == aAtom1 )
                break;
        for ( i2=0; i2<MAXCONNECT; i2++ )
            if ( aResidueConnectAtom( rR2, i2 ) == aAtom2 )
                break;
        if ( i1 < MAXCONNECT  &&  i2 < MAXCONNECT ) {
             *  predefined connect atoms used, so make sure
             *  residue #s reflect connect order
            if ( i1 < i2 )
                OrderResidues( rR1, rR2 );
            else if ( i2 < i1 )
                OrderResidues( rR2, rR1 );
            else
                VP0(( "Connected same connect atoms (%d), %s\n",
                        i1, "residue order not checked" ));
        }
    }
     */

    /*
     *  make the bond
     */
    aAtom1->aaBonds[iAtomCoordination(aAtom1)] = aAtom2;
    AtomDefineBondFlags( aAtom1, iAtomCoordination(aAtom1), 0 );
    AtomSetBondOrder( aAtom1, iAtomCoordination(aAtom1), BONDSINGLE );
    aAtom1->iCoordination++;
    aAtom2->aaBonds[iAtomCoordination(aAtom2)] = aAtom1;
    AtomDefineBondFlags( aAtom2, iAtomCoordination(aAtom2), 0 );
    AtomSetBondOrder( aAtom2, iAtomCoordination(aAtom2), BONDSINGLE );
    aAtom2->iCoordination++;
    CDU(aAtom1);
    CDU(aAtom2);
}

/*
 *      AtomTmpBondTo
 *
 *      Create a bond between two atoms, with bond order BONDSINGLE.
 *      The user is responsible for checking if there are already
 *      the maximum number of bonds on the atom. Allow unreasonable
 *      bonds for hack purposes.
 */
BOOL
AtomTmpBondTo( ATOM aAtom1, ATOM aAtom2 )
{
STRING  sTemp;
int     ierr = 0;

    VERIFYOBJEKT( aAtom1, ATOMid );
    VERIFYOBJEKT( aAtom2, ATOMid );

    if ( iAtomCoordination( aAtom1 ) >= MAXBONDS ) {
        VP0(("Can't make temp bond - out of bonds for %s\n",
                        sContainerFullDescriptor( (CONTAINER)aAtom1, sTemp ) ));
        ierr = 1;
    }
    if ( iAtomCoordination( aAtom2 ) >= MAXBONDS ) {
        VP0(("Can't make temp bond - out of bonds for %s\n",
                        sContainerFullDescriptor( (CONTAINER)aAtom2, sTemp ) ));
        ierr = 1;
    }
    if ( ierr )
        return(FALSE);
        
    /*
     *  make the bond
     */
    aAtom1->aaBonds[iAtomCoordination(aAtom1)] = aAtom2;
    AtomDefineBondFlags( aAtom1, iAtomCoordination(aAtom1), 0 );
    AtomSetBondOrder( aAtom1, iAtomCoordination(aAtom1), BONDSINGLE );
    aAtom1->iCoordination++;
    aAtom2->aaBonds[iAtomCoordination(aAtom2)] = aAtom1;
    AtomDefineBondFlags( aAtom2, iAtomCoordination(aAtom2), 0 );
    AtomSetBondOrder( aAtom2, iAtomCoordination(aAtom2), BONDSINGLE );
    aAtom2->iCoordination++;

    return(TRUE);
}


/*
 *      AtomBondToOrder
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a bond between two atoms, with the specified ORDER
 *      The user is responsible for checking if there are already
 *      the maximum number of bonds on the atom.
 */
void
AtomBondToOrder( ATOM aAtom1, ATOM aAtom2, int iOrder )
{

    VERIFYOBJEKT( aAtom1, ATOMid );
    VERIFYOBJEKT( aAtom2, ATOMid );
    
    if ( bBondAtomProblem( aAtom1, aAtom2 ) )
        return;
 
    aAtom1->aaBonds[iAtomCoordination(aAtom1)] = aAtom2;
    AtomDefineBondFlags( aAtom1, iAtomCoordination(aAtom1), 0 );
    AtomSetBondOrder( aAtom1, iAtomCoordination(aAtom1), iOrder );
    aAtom1->iCoordination++;
    aAtom2->aaBonds[iAtomCoordination(aAtom2)] = aAtom1;
    AtomDefineBondFlags( aAtom2, iAtomCoordination(aAtom2), 0 );
    AtomSetBondOrder( aAtom2, iAtomCoordination(aAtom2), iOrder );
    aAtom2->iCoordination++;
    CDU(aAtom1);
    CDU(aAtom2);
}




/*
 *      AtomBondToFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a bond between two atoms, with the specified FLAGS 
 *      The user is responsible for checking if there are already
 *      the maximum number of bonds on the atom.
 */
void
AtomBondToFlags( ATOM aAtom1, ATOM aAtom2, FLAGS fFlags )
{

    VERIFYOBJEKT( aAtom1, ATOMid );
    VERIFYOBJEKT( aAtom2, ATOMid );
    
    if ( bBondAtomProblem( aAtom1, aAtom2 ) )
        return;

    aAtom1->aaBonds[iAtomCoordination(aAtom1)] = aAtom2;
    AtomDefineBondFlags( aAtom1, iAtomCoordination(aAtom1), fFlags );
    aAtom1->iCoordination++;
    aAtom2->aaBonds[iAtomCoordination(aAtom2)] = aAtom1;
    AtomDefineBondFlags( aAtom2, iAtomCoordination(aAtom2), fFlags );
    aAtom2->iCoordination++;
    CDU(aAtom1);
    CDU(aAtom2);
}



/*
 *      AtomRemoveBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove a bond between two atoms.
 *      The user is responsible for checking that a bond exists between
 *      the atoms.
 */
void
AtomRemoveBond( ATOM aAtom1, ATOM aAtom2 )
{

    VERIFYOBJEKT( aAtom1, ATOMid );
    VERIFYOBJEKT( aAtom2, ATOMid );
    
    if ( iAtomCoordination(aAtom1) <= 0 ) 
        DFATAL( ("There is no bond!") );
    if ( iAtomCoordination(aAtom2) <= 0 ) 
        DFATAL( ("There is no bond!") );
    
    zAtomRemoveBond( aAtom1, aAtom2 );
    zAtomRemoveBond( aAtom2, aAtom1 );
    CDU(aAtom1);
    CDU(aAtom2);
}




/*
 *      bAtomBondedTo
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return TRUE if aAtom1 is bonded to aAtom2.
 */
BOOL
bAtomBondedTo( ATOM aAtom1, ATOM aAtom2 )
{
int             i;

    VERIFYOBJEKT( aAtom1, ATOMid );
    VERIFYOBJEKT( aAtom2, ATOMid );
    
    if ( iAtomCoordination(aAtom1) <= 0 ) return(FALSE);
    if ( iAtomCoordination(aAtom2) <= 0 ) return(FALSE);
    
    for ( i=0; i<iAtomCoordination(aAtom1); i++ )
        if ( aAtom1->aaBonds[i] == aAtom2 ) break;
    if ( i == iAtomCoordination(aAtom1) ) return(FALSE);
    
    return(TRUE);
}




/*
 *      iAtomFindBondOrder
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the bond type between aAtom and aNeighbor.
 */
int
iAtomFindBondOrder( ATOM aAtom, ATOM aNeighbor )
{
int             i;

    for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
        if ( aAtomBondedNeighbor( aAtom, i ) == aNeighbor )
                return(iAtomBondOrder( aAtom, i ));
    }
    return(BONDNONE);
}




/*
 *      AtomFindSetBondOrder
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the bond type between aAtom and aNeighbor.
 */
void
AtomFindSetBondOrder( ATOM aAtom, ATOM aNeighbor, int iOrder )
{
int             i;

                /* Change the bond order for both the atom and its */
                /* neighbor */

    for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
        if ( aAtomBondedNeighbor( aAtom, i ) == aNeighbor )
                AtomSetBondOrder( aAtom, i, iOrder );
    }
    for ( i=0; i<iAtomCoordination(aNeighbor); i++ ) {
        if ( aAtomBondedNeighbor( aNeighbor, i ) == aAtom )
                AtomSetBondOrder( aNeighbor, i, iOrder );
    }
    CDU(aAtom);
    CDU(aNeighbor);
}




/*
 *      fAtomFindBondFlags
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the bond flags between aAtom and aNeighbor.
 */
FLAGS
fAtomFindBondFlags( ATOM aAtom, ATOM aNeighbor )
{
int             i;

    for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
        if ( aAtomBondedNeighbor( aAtom, i ) == aNeighbor )
                return(fAtomBondFlags( aAtom, i ));
    }
    return(BONDNONE);
}





/*
 *      aAtomDuplicate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Duplicate an atom and return it.
 *      Should only be called from cContainerDuplicate()
 *      which is expected to duplicate the CONTAINER variables.
 */
ATOM
aAtomDuplicate( ATOM aAtom )
{
ATOM    aNewAtom;

    MALLOC( aNewAtom, ATOM, sizeof(ATOMt) );
    memcpy( aNewAtom, aAtom, sizeof(ATOMt) );
    aNewAtom->iUniqueId = SiUniqueId++;
    if ( GfAtomClassGraphicsCreator != NULL ) {
        GfAtomClassGraphicsCreator(aNewAtom);
    }
    return(aNewAtom);   
}





/*
 *      AtomResetPointers
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Reset all internal pointers.
 *      When copies are made of ATOMs, internal pointers, such
 *      as bonds, are left pointing to the old ATOMs.
 *      However, the CONTAINER class supports an extra pointer
 *      on the old CONTAINER which points to the new CONTAINER.
 *      This can be used to reset the internal pointers of the ATOM.
 *      The value of the pointer in new CONTAINERs is NULL.
 *
 *      EVERY internal pointer of the ATOM has to be resolved.
 */
void
AtomResetPointers( ATOM aAtom )
{
int     iBond;

    for ( iBond=0; iBond<iAtomCoordination(aAtom); iBond++ ) 
        aAtom->aaBonds[iBond] = 
                (ATOM)cContainerCopyPointer(aAtom->aaBonds[iBond]);
}











/*
 *      AtomCheck
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Check the ATOM and its contents as the whether or
 *      not calculations can be run on it.
 *
 *      Arguments:
 *              aAtom -         Atom to check.
 *              iPErrors -      Add the number of errors found.
 *              iPWarnings -    Add the number of warnings found.
 */
void
AtomCheck( ATOM aAtom, int *iPErrors, int *iPWarnings )
{
STRING          sTemp;

        /* FATAL error if the atom doesn't have a type */

    if ( strlen(sAtomType(aAtom)) == 0 ) {
        VP0(( "FATAL:  Atom %s does not have a type.\n",
                sContainerFullDescriptor( (CONTAINER)aAtom, sTemp ) ));
        (*iPErrors)++;
        iFatal++;
    }

    if ( bAtomPerturbed( aAtom ) ) {

        /* WARNING if the atom has a perturbation type but no pert name */
        
        if ( strlen(sAtomPertName(aAtom)) == 0 &&
                strlen(sAtomPertType(aAtom)) != 0 ) {
                VP0(( 
                    "WARNING:  Atom %s has a perturbation type, but no name.\n",
                        sContainerFullDescriptor( (CONTAINER)aAtom, sTemp ) ));
                (*iPWarnings)++;
        }

        /* WARNING if the atom has a pert name but no pert type */

        if ( strlen(sAtomPertName(aAtom)) != 0 &&
                strlen(sAtomPertType(aAtom)) == 0 ) {
                VP0(( "WARNING:  Atom %s has a perturbation name,%s",
                        sContainerFullDescriptor( (CONTAINER)aAtom, sTemp ),
                                "\nbut no perturbation type.\n" ));
                (*iPWarnings)++;
        }
    }

        /* Warning for each bond that was made through a distance search */

        /* Warning if the atoms type was calculated using CALCTYPES */
}




/*
 *      AtomYouAreBeingRemoved
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      The ATOM is being told that it is being removed
 *      from some other CONTAINER.  This instructs it to 
 *      break ALL bonds with other ATOMs.
 *      It also removes all INTERNALs that it contains.
 *      Inform the USER about what is going on.
 */
void
AtomYouAreBeingRemoved( ATOM aAtom )
{
int             i, iNumberOfBonds;
ATOM            aBond;
STRING          sAtom1, sAtom2;
LOOP            lInternals;
INTERNAL        iInt;
ATOM            aaBond[MAXCONNECT];

        /* First break all bonds that this ATOM is associated with */

    /*
    **  get info about this atom
    */
    iNumberOfBonds = iAtomCoordination(aAtom);
    sContainerFullDescriptor( (CONTAINER)aAtom, sAtom1 );

    /*
    **  init bonded atom tracking array &
    **  make list of atoms bonded to
    */
    for ( i = 0; i < MAXCONNECT; i++ ) 
        aaBond[i] = NULL;
    for ( i = 0; i < iNumberOfBonds; i++ ) 
        aaBond[i] = aAtomBondedNeighbor( aAtom, i );

    /*
    **  delete the bond
    */
    for ( i = 0; i < iNumberOfBonds; i++ ) {
        aBond = aaBond[i];
        VP1(( "Breaking bond: %s - %s\n",
                sAtom1,
                sContainerFullDescriptor( (CONTAINER)aBond, sAtom2 ) ));
        AtomRemoveBond( aAtom, aBond );
    }

        /* Then remove all INTERNALs that this ATOM is associated with */

    lInternals = lLoop( (OBJEKT)aAtom, INTERNALS );
    while ( (iInt = (INTERNAL)oNext(&lInternals)) ) {
        InternalDestroy(&iInt);
    }

        /* Now pass the message up the hierarchy that this ATOM */
        /* is being removed */

    ContainerIAmBeingRemoved( cContainerWithin(aAtom), (CONTAINER)aAtom );
}




/*
 *      AtomIAmBeingRemoved
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This is a dummy procedure, nothing should
 *      every be contained within the ATOM to be removed.
 *      INTERNALs are not handled the same way as real
 *      CONTAINERs.
 */
void
AtomIAmBeingRemoved( ATOM aAtom, CONTAINER cRemoved )
{
    DFATAL(( "AtomIAmBeingRemoved should NEVER be called!\n" ));
}





/*
 *      AtomSetAttribute
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Change the named attribute of the ATOM.
 *      This routine is called mainly from the command
 *      line interface to allow users to change the
 *      attributes of ATOMs by name.
 */
void
AtomSetAttribute( ATOM aAtom, STRING sAttr, OBJEKT oAttr )
{
LISTLOOP        llElements;
ASSOC           aXCoord, aYCoord, aZCoord;
VECTOR          vPos;
int             iElement;
STRING          sTemp;

    if ( strcmp( sAttr, "charge" ) == 0 ) {
        if ( !bObjektWarnType( oAttr, ODOUBLEid ) ) return;
        AtomSetCharge( aAtom, dODouble(oAttr) );
        return;
    }
    if ( strcmp( sAttr, "pert" ) == 0 ) {
        if ( !bObjektWarnType( oAttr, OSTRINGid ) ) return;
        strcpy( sTemp, sOString(oAttr) );
        StringLower(sTemp);
        if ( strcmp( sTemp, "true" ) == 0 ) {
            AtomSetFlags( aAtom, ATOMPERTURB );
        } else if ( strcmp( sTemp, "false") == 0) {
            AtomResetFlags( aAtom, ATOMPERTURB );
        } else
            VP0(( "pert: expected 'true' or 'false'\n" ));
        return;
    } 
    if ( strcmp( sAttr, "pertCharge" ) == 0 ) {
        if ( !bObjektWarnType( oAttr, ODOUBLEid ) ) return;
        AtomSetPertCharge( aAtom, dODouble(oAttr) );
        return;
    }
    if ( strcmp( sAttr, "name" ) == 0 ) {
        if ( !bObjektWarnType( oAttr, OSTRINGid ) ) return;
        ContainerSetName( aAtom, sOString(oAttr) );
        return;
    }
    if ( strcmp( sAttr, "pertName" ) == 0 ) {
        if ( iObjectType(oAttr) == NULLid ) {
            AtomSetPertName( aAtom, "" );
        } else if ( !bObjektWarnType( oAttr, OSTRINGid ) ) return;
        AtomSetPertName( aAtom, sOString(oAttr) );
        return;
    }
    if ( strcmp( sAttr, "type" ) == 0 ) {
        if ( iObjectType(oAttr) == NULLid ) {
            AtomSetType( aAtom, "" );
        } else if ( !bObjektWarnType( oAttr, OSTRINGid ) ) 
            return;
        AtomSetType( aAtom, sOString(oAttr) );
        return;
    }
    if ( strcmp( sAttr, "element" ) == 0 ) {
        if ( iObjectType(oAttr) == NULLid ) {
            AtomSetElement( aAtom, NOELEMENT );
        } else if ( !bObjektWarnType( oAttr, OSTRINGid ) ) 
            return;
        iElement = iElementNumber(sOString(oAttr));
        if ( iElement == NOELEMENT ) {
            VP0(( "Unknown element: %s\n", sOString(oAttr) ));
        } else {
            AtomSetElement( aAtom, iElement );
        }
        return;
    }
    if ( strcmp( sAttr, "pertType" ) == 0 ) {
        if ( iObjectType(oAttr) == NULLid ) {
            AtomSetPertType( aAtom, "" );
        } else if ( !bObjektWarnType( oAttr, OSTRINGid ) ) return;
        AtomSetPertType( aAtom, sOString(oAttr) );
        if( !GDefaults.iGibbs ){
                if( strcmp( sAtomType( aAtom ), sAtomPertType( aAtom ))){
                        AtomSetFlags( aAtom, ATOMPERTURB );
                } else {
                        AtomResetFlags( aAtom, ATOMPERTURB );
                }
        }
        return;
    }
    if ( strcmp( sAttr, "position" ) == 0 ) {
        if ( !bObjektWarnType( oAttr, LISTid ) ) return;
        if ( iListSize( oAttr ) != 3 ) {
            VP0(( "Illegal vector\n" ));
        } else {
            llElements = llListLoop( (LIST)oAttr );
            aXCoord = (ASSOC)oListNext(&llElements);
            aYCoord = (ASSOC)oListNext(&llElements);
            aZCoord = (ASSOC)oListNext(&llElements);
            if ( !bObjektWarnType( oAssocObject(aXCoord), ODOUBLEid ) ||
                 !bObjektWarnType( oAssocObject(aYCoord), ODOUBLEid ) ||
                 !bObjektWarnType( oAssocObject(aZCoord), ODOUBLEid ) ) return;
            VectorDef( &vPos,   dODouble(oAssocObject(aXCoord)),
                                dODouble(oAssocObject(aYCoord)),
                                dODouble(oAssocObject(aZCoord)) );
            AtomSetPosition( aAtom, vPos );
        }
        return;
    }
    VP0(( "%s: non-existent attribute for an atom.\n", sAttr ));
}



/*
 *      iAtomHybridization
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the hybridization state of the atom.
 *      This is calculated by counting the
 *      number and types of bonds going into the 
 *      atom.  This routine will return a hybridization
 *      for all kinds of weird bonding arrangements.
 */
int
iAtomHybridization( ATOM aAtom )
{
int             iSingle, iDouble, iTriple, iAromatic;
int             i, iH;
STRING          sTemp, sDesc;
double          dMass, dPolar, dEpsilon, dRStar, dEpsilon14, dRStar14;
double		dScreenF;
int             iElement, iTag;
PARMSET         psSet;

        /* Check if the hybridization can be calculated from */
        /* the atom type */

                /* First check if we can find it in the */
                /* PARMLIBRARY */

    if ( bParmLibDefaultExists() ) {            
        PARMLIB_DEFAULT_LOOP( psSet, 
            (iTag = iParmSetFindAtom( psSet, sAtomType(aAtom) )) );
        if ( iTag != PARM_NOT_FOUND ) {
            ParmSetAtom( psSet, iTag, sTemp, 
                        &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14, &dRStar14,
			&dScreenF,
                        &iElement, &iH, sDesc );
            return(iH);
        }
    }

    iSingle = iDouble = iTriple = iAromatic = 0;
    for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
        switch ( iAtomBondOrder( aAtom, i ) ) {
            case BONDSINGLE:
                iSingle++;
                break;
            case BONDDOUBLE:
                iDouble++;
                break;
            case BONDTRIPLE:
                iTriple++;
                break;
            case BONDAROMATIC:
                iAromatic++;
                break;
            default:
                DFATAL(( "There is an illegal bond type" ));
                break;
        }
    }
                /* one or more aromatic bonds makes the atom SP2 */
    if ( iAromatic != 0 )       iH = HSP2;

                /* one or more triple bonds makes the atom SP1 */
    else if ( iTriple != 0 )    iH = HSP1;

                /* Two or more double bonds makes the atom linear, SP1 */
    else if ( iDouble >= 2 )    iH = HSP1;

                /* One double bond makes the atom SP2 */
    else if ( iDouble != 0 )    iH = HSP2;

                /* Otherwise the atom is SP3 */
    else                        iH = HSP3;

    return(iH);
}



/*
 *      iAtomBondOrderFromName
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the bond order associated with the name.
 */
int
iAtomBondOrderFromName( char *sName )
{
int             i, iNum, iOrder;

    iNum = sizeof(SoaOrders)/sizeof(SoaOrders[0]);
    for ( i=0; i<iNum; i++ ) {
        if ( strcmp( SoaOrders[i].sName, sName ) == 0 ) {
            iOrder = SoaOrders[i].iOrder;
            return(iOrder);
        }
    }
    return(BONDNONE);
}

/*
 *      bAtomSpaceConflict
 *
 *      Author: David A. Rivkin (1992)
 *
 *      Determines if the two atoms are have Van der Waals radii
 *      that are overlapping.
 *      NOTE:  Currently there is no test for TOTALLY missing parameters
 *              or if there are no types assigned to the atoms.
 *
 */
BOOL
bAtomSpaceConflict( ATOM aAtom1, ATOM aAtom2 )
{
double  dR1, dR2, dDist, dX, dY, dZ;

    /* Determine the distance that is needed */
    if ( (dR1 = dAtomVanderWaals( aAtom1 )) < 0.0 )    
        dR1 = 2.0;
    if ( (dR2 = dAtomVanderWaals( aAtom2 )) < 0.0 )    
        dR2 = 2.0;
    dDist = dR1 + dR2;

    /* Determine if atom1 is within the box that the Van Der Waals
       radii sum around atom2 */

    dX = dVX(&vAtomPosition( aAtom1 )) - dVX(&vAtomPosition( aAtom2 )) ;
    dY = dVY(&vAtomPosition( aAtom1 )) - dVY(&vAtomPosition( aAtom2 )) ;
    dZ = dVZ(&vAtomPosition( aAtom1 )) - dVZ(&vAtomPosition( aAtom2 )) ;
    if ( abs(dX) < dDist  &&  abs(dY) < dDist  &&  abs(dZ) < dDist ) {
        /*
         *  check actual distance
         */
        if ( dDist * dDist > dX * dX + dY * dY + dZ * dZ )
            return( TRUE );
    }
    return( FALSE );
}

/*
 *      dAtomVanderWaals
 *
 *      Author: David A. Rivkin (1992)
 *
 *      Determines the Van der Waals Radius of the atom.
 *
 */
double
dAtomVanderWaals( ATOM aAtom )
{
PARMSET psTemp;
int     iTemp, iElement, iHybrid;
STRING  sType, sDesc;
double  dMass, dPolar, dEpsilon, dR, dEpsilon14, dR14;
double  dScreenF;

    /* Check to see if it the Van der Waals for this type are already in the 
        Global VarArray */
    if ( ! GvaVDWTypes ) {
        GvaVDWTypes = vaVarArrayCreate( sizeof( char ) * ATOMTYPELEN );
        GvaVDWValues = vaVarArrayCreate( sizeof( double ));
    }
    
    for ( iTemp = 0; iTemp < iVarArrayElementCount( GvaVDWTypes ); iTemp++) {
        if ( !strcmp( *PVAI( GvaVDWTypes, char *, iTemp ), 
                sAtomType( aAtom ))) {
            return( *PVAI( GvaVDWValues, int, iTemp ) );
        }       
    }
    /*
     *  not found; look in parm libs & add to array if found
     *  HACK, since this array can't be superseded by changes
     *  to parmlibs
     */
    strcpy( sType, sAtomType( aAtom ));
    PARMLIB_LOOP( GplDefaultParmLib, psTemp, 
                ( iTemp = iParmSetFindAtom( psTemp, sType)));
    if( iTemp == PARM_NOT_FOUND) {
        return( -1.0 );
    }
    ParmSetAtom( psTemp, iTemp, sType, &dMass, &dPolar, &dEpsilon, &dR, 
                        &dEpsilon14, &dR14, &dScreenF, &iElement, &iHybrid,
		       	sDesc );
    MESSAGE(( "Type:  %s  r*:  %f\n", sType, dR ));
    VarArrayAdd( GvaVDWTypes, (GENP)sType );
    VarArrayAdd( GvaVDWValues, (GENP)&dR );

    return( dR );
}



int
iAtomSetTmpRadius( ATOM aAtom )
{
PARMSET psTemp;
int     iTag, iElement, iHybrid;
double  dMass, dPolar, dE, dR, dE14, dR14;
double  dScreenF;
STRING  sType, sDesc;
int     iDefaultedRadius = 0;

        iTag = PARM_NOT_FOUND;
        PARMLIB_DEFAULT_LOOP( psTemp,
                ( iTag = iParmSetFindAtom( psTemp, sAtomType(aAtom) ) ));
        if ( iTag != PARM_NOT_FOUND ) {
            ParmSetAtom( psTemp, iTag, sType,
                   &dMass, &dPolar, &dE, &dR, &dE14, &dR14, 
		   &dScreenF,
                                   &iElement, &iHybrid, sDesc );
        } else {
                iDefaultedRadius = 1;
                dR = ATOM_DEFAULT_RADIUS;
        }
        AtomSetTempDouble( aAtom, dR );
        return(iDefaultedRadius);
}
