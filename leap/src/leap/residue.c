/*
 *      File: residue.c
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
 *              RESIDUE
 *      Superclass:
 *              CONTAINER, LOOP
 *
 *      Description:
 *              A RESIDUE is a subclass of CONTAINER.
 *              RESIDUES can contain atoms only.
 *              The property that destinguishes RESIDUES
 *              from all other CONTAINER subclasses is that
 *              the graph defined by the bonds of the atoms
 *              within the RESIDUE is a SINGLE CLOSED or UNCLOSED graph.
 *              There MUST be a connection between any two atoms
 *              in the residue.
 *              Functions must be called to ensure that this
 *              property is maintained.
 *
 */
 



#include	"basics.h"

#include        "classes.h"

#include        "build.h"

/*
 *----------------------------------------------------------
 *
 *	Private variables
 *
 */


typedef	struct	{
		STRING		sName;
		int		iConnect;
		} CONNECTt;


static	CONNECTt	ScaConnectNames[] = {
	{ "connect0", 	CONNECT0 },
	{ "firstend", 	FIRSTEND },
	{ "nend", 	NEND },
	{ "connect1", 	CONNECT1 },
	{ "lastend", 	LASTEND },
	{ "cend", 	CEND },
	{ "connect2",	CONNECT2 },
	{ "send", 	SEND },
	{ "disulphide",	SEND },
	{ "connect3",	CONNECT3 },
	{ "connect4",	CONNECT4 },
	{ "connect5",	CONNECT5 } };

/*
 *	iResidueConnectFromName
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the connection ATOM index from the
 *	name.  Search for the name in the ScaConnectNames
 *	array.
 *	Return NOEND if no match could be found.
 */
int
iResidueConnectFromName( char *sName )
{
int		i, iNum;

    iNum = sizeof(ScaConnectNames)/sizeof(ScaConnectNames[0]);
    for ( i=0; i<iNum; i++ ) {
	if ( strcmp( sName, ScaConnectNames[i].sName ) == 0 ) {
	    return(ScaConnectNames[i].iConnect);
	}
    }
    return(NOEND);
}





/*
 *      rResidueCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a RESIDUE.
 *
 *      Return:
 *              Return a pointer to the RESIDUE created.
 */
RESIDUE
rResidueCreate()
{
RESIDUE m;
int     i;

    MALLOC( m, RESIDUE, sizeof(RESIDUEt) );
    m->cHeader.dDisp = NULL;
    m->vaImpropers = NULL;
    for ( i=0; i<MAXCONNECT; i++ ) {
        m->aaConnect[i] = NULL;
    }
    ResidueDefineFlags( m, 0 );
    ResidueSetType( m, RESTYPEUNDEFINED );
    ResidueSetImagingAtom( m, NULL );
    ResidueSetPdbSequence( m, 0 );
    return(m);
}







/*
 *      ResidueDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the RESIDUE and all of the residues/atoms within it.
 *
 *      Arguments:
 *      oObject -       A pointer to the object.
 */
void
ResidueDestroy( RESIDUE *rPResidue )
{
    if ( (*rPResidue)->vaImpropers )
	VarArrayDestroy( &(*rPResidue)->vaImpropers );

    FREE( *rPResidue );
    *rPResidue = NULL; 
}






/*
 *      ResidueDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to print a description 
 *      of the RESIDUE to the file fOut.
 *
 *      Arguments:
 *      oObject -       The object to describe.
 *
 */
void
ResidueDescribe( RESIDUE rResidue )
{
int             i;
LOOP            lContents;
OBJEKT          oObj;
STRING          sTemp;

    VP0(( "RESIDUE name: %s\n", rResidue->cHeader.sName ));
    if ( bResidueFlagsSet( rResidue, RESIDUEUNKNOWN ) ) 
        VP0(( "!!!This is an unknown residue!\n" ));
    VP0(( "RESIDUE sequence number: %d\n", 
                        iContainerSequence(rResidue) ));
    VP0(( "RESIDUE PDB sequence number: %d\n",
			iResiduePdbSequence(rResidue) ));
    VP0(( "Type: %s\n", 
	  sResidueTypeNameFromChar(cResidueType(rResidue)) ));
    if ( cResidueType(rResidue) == RESTYPESOLVENT ) {
	VP0(( "Solvent imaging atom: " ));
	if ( aResidueImagingAtom(rResidue) == NULL ) {
	    VP0(( "None.\n" ));
	} else {
	    VP0(( "%s\n", 
		sContainerDescriptor((CONTAINER)aResidueImagingAtom(rResidue),sTemp) ));
	}
    }	
    VP0(( "Connection atoms:\n" ));
     for ( i=0; i<MAXCONNECT; i++ ) {
        if ( bResidueConnectUsed(rResidue,i) ) {
            VP0(( " Connect atom %d: %s\n", i, 
                sContainerDescriptor( (CONTAINER)rResidue->aaConnect[i], sTemp ) ));
        }
    }
    if ( rResidue->vaImpropers ) {
	VP0(( "No improper torsions\n" ));
    } else {
	VP0(( "Improper torsions:\n" ));
    }
    BasicsResetInterrupt();
    VP0(( "Contents: \n" ));
    lContents = lLoop( (OBJEKT)rResidue, DIRECTCONTENTS );
    while ( (oObj = oNext(&lContents )) ) {
        VP0(( "%s\n", sContainerDescriptor( (CONTAINER)oObj, sTemp ) ));
	if ( bBasicsInterrupt() ) {
	    BasicsResetInterrupt();
	    VP0(( "Interrupted\n" ));
	    break;
	}
    }
}









/*
 *      rResidueDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a duplicate of the residue.
 *	Should only be called from cContainerDuplicate()
 *      which is expected to duplicate the CONTAINER variables.
 */
RESIDUE
rResidueDuplicate( RESIDUE rOld )
{
RESIDUE	rNew;

    MALLOC(rNew, RESIDUE, sizeof(RESIDUEt) );
    memcpy( rNew, rOld, sizeof(RESIDUEt) );

    if ( rOld->vaImpropers ) {
	rNew->vaImpropers = vaVarArrayCopy( rOld->vaImpropers );
    }
    return(rNew);
}




/*
 *      ResidueResetPointers
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
 */
void
ResidueResetPointers( RESIDUE rRes )
{
int     i;

    for ( i=0; i<MAXCONNECT; i++ ) {
        if ( bResidueConnectUsed(rRes,i) ) {
            ResidueSetConnectAtom( rRes, i, 
                 (OBJEKT)cContainerCopyPointer(rRes->aaConnect[i]));
        }
    }

	/* Reset the solvent imaging ATOM */

    if ( aResidueImagingAtom(rRes) != NULL ) {
        ResidueSetImagingAtom(rRes, 
		cContainerCopyPointer(aResidueImagingAtom(rRes)) );
    }
}





/*
 *      ResidueConnect
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Connect two residues together at the connection points
 *      specified
 *
 *      Arguments:
 *      rResA           - First Residue
 *      iConnectA       - The connection point on the first Residue
 *      rResB           - Second Residue
 *      iConnectB       - The connection point on the second Residue
 */
void
ResidueConnect( RESIDUE rResA, int iConnectA, RESIDUE rResB, int iConnectB )
{

    if ( !bResidueConnectUsed(rResA,iConnectA) )
        DFATAL( ("Connection %d is not used", iConnectA ) );
    if ( !bResidueConnectUsed(rResB,iConnectB) )
        DFATAL( ("Connection %d is not used", iConnectB ) );
        
    AtomBondTo( rResA->aaConnect[iConnectA],
                rResB->aaConnect[iConnectB] );
}




/*
 *      ResidueCheck
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Check the RESIDUE and its contents as the whether or
 *      not calculations can be run on it.
 *
 *      Arguments:
 *              rRes -          Residue to check.
 *              iPErrors -      Add the number of errors found.
 *              iPWarnings -    Add the number of warnings found.
 */
void
ResidueCheck( RESIDUE rRes, int *iPErrors, int *iPWarnings )
{
LOOP            lContents;
CONTAINER       cCont;
STRING          sTemp;

        /* Print a warning if the residue is an unknown one */

    if ( bResidueFlagsSet( rRes, RESIDUEUNKNOWN ) ) {
        (*iPWarnings)++;
        VP0(( "Warning: Unknown residue: %s\n", 
                sContainerFullDescriptor( (CONTAINER)rRes, sTemp ) ));
    }

    lContents = lLoop( (OBJEKT)rRes, DIRECTCONTENTS );
    while ( (cCont = (CONTAINER)oNext(&lContents)) ) {
        ContainerCheck( cCont, iPErrors, iPWarnings );
    }
}




/*
 *      ResidueMutate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Mutate the RESIDUE (rOld) into the RESIDUE (rNew).
 *      Do this by superimposing the coordinates from rOld onto
 *      rNew, breaking all the bonds from rOld to other residues,
 *      and rejoining them to identically named atoms in rNew,
 *      then building the coordinates for the rest of the atoms in rNew.
 *      The new RESIDUE rNew should not be bonded to anything else, it
 *      should also have EXTERNAL coordinates defined, but no INTERNALS.
 *
 */
void
ResidueMutate( RESIDUE rNew, RESIDUE rOld )
{
LOOP            lAtoms, lSpan;
ATOM            aNew, aNeighbor, aOld, aAtom, aSpan, aTemp;
#ifdef  DEBUG
STRING          sTemp, sSpan;
#endif
int             i, iDum;
FLAGS		fBondFlags;

    MESSAGE(( "Mutating: %s to: %s\n", sContainerName( rOld ),
                sContainerName(rNew) ));

                /* Build internal coordinates for the new RESIDUE */

    lAtoms = lLoop( (OBJEKT)rNew, ATOMS );
    BuildInternalsUsingFlags( &lAtoms, 
			ATOMPOSITIONKNOWN, 
			0,
			0, 
			ATOMPOSITIONKNOWN
			);

                /* Define the coordinates, and the flags */
                /* if there are bonds out of the old RESIDUE, break them */
                /* and rejoin them to identically named atoms in the new */
                /* RESIDUE */

    lAtoms = lLoop( (OBJEKT)rOld, ATOMS );
    FOREACH( aOld, ATOM, lAtoms ) {
        MESSAGE(( "Searching for atom in new residue with name: %s\n",
                        sContainerName(aOld) ));
        aNew = (ATOM)cContainerFindName( (CONTAINER)rNew, ATOMid,
                                             sContainerName(aOld) );
                                             
                /* If there is a cooresponding ATOM with the same name */
                /* then define its flags and coordinates */
                
        if ( aNew != NULL ) {
            MESSAGE(( "--- Found one\n" ));
            AtomSetPosition( aNew, vAtomPosition(aOld) );
            AtomDefineFlags( aNew, fAtomFlags(aOld) );
        } else {
            MESSAGE(( "--- No atom found\n" ));
        }

                /* Search for bonds out of the old RESIDUE */
        
        for ( i=0; i<iAtomCoordination(aOld); i++ ) {
            aNeighbor = aAtomBondedNeighbor( aOld, i );
	    MESSAGE(( "--- Looking at neighbor: %s\n",
			sContainerFullDescriptor( (CONTAINER)aNeighbor, sTemp ) ));
            if ( rOld != (RESIDUE)cContainerWithin(aNeighbor) ) {
                fBondFlags = fAtomBondFlags( aOld, i );
                AtomRemoveBond( aOld, aNeighbor );
                MESSAGE(( "Removing a bond to: %s\n",
                        sContainerFullDescriptor( (CONTAINER)aNeighbor, sTemp ) ));
                if ( aNew != NULL ) {
                    MESSAGE(( "--- And rejoining it to: %s\n",
                                sContainerFullDescriptor( (CONTAINER)aNew, sTemp ) ));
                    AtomBondToFlags( aNew, aNeighbor, fBondFlags );
                } else {
                    MESSAGE(( "--- Not rejoining it to anything.\n" ));
        VP1(( "There is no atom in residue: %s with the name: %s.\n" ));
        VP1(( "--- No bond could be made to the missing atom.\n" ));
                }
            }
        }
    }

        /* Build the coordinates for ATOMs that do not have them */

    BuildFixInternals( rNew );

                /* Loop through all ATOMs looking for those that */
                /* do not have positions known and build externals */
                /* for them and neighbors that are bonded to them */

    lAtoms = lLoop( (OBJEKT)rNew, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) ) {
            lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
            LoopDefineInvisibleAtoms( &lSpan, ATOMPOSITIONKNOWN );

                        /* Look for a collision with an ATOM whos */
                        /* ATOMPOSITIONKNOWN flag is set */
            aSpan = NULL;
            while ( (aTemp = (ATOM)oNext(&lSpan)) ) {
                if ( aSpan == NULL ) aSpan = aTemp;
                if ( aLoopLastCollisionAtom(&lSpan) != NULL ) {
                    aSpan = aLoopLastCollisionAtom(&lSpan);
                    break;
                }
            }

            MESSAGE(( "Building externals from: %s\n",
                        sContainerFullDescriptor( (CONTAINER)aSpan, sSpan ) ));

                        /* Build external coordinates */
            lSpan = lLoop( (OBJEKT)aSpan, SPANNINGTREE );
            LoopDefineInvisibleAtoms( &lSpan, ATOMPOSITIONKNOWN );
	    iDum = 0;	/* for purify */
            BuildExternalsUsingFlags( &lSpan,
                                        0, ATOMPOSITIONKNOWN,
                                        ATOMPOSITIONKNOWN, 0,
					&iDum, &iDum, &iDum, TRUE );
	}
    }

        /* Destroy the INTERNALs */

    lAtoms = lLoop( (OBJEKT)rNew, ATOMS );
    BuildDestroyInternals( &lAtoms );
    
}


/*
 *	rResidueConnected
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the RESIDUE that is connected to the main
 *	RESIDUE through connect ATOM iConnect.
 */
RESIDUE
rResidueConnected( RESIDUE rRes, int iConnect )
{
ATOM		aAtom;
int		i;

    if ( ( aAtom = aResidueConnectAtom( rRes, iConnect ) ) == NULL )
	return(NULL);

    for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
	if ( (RESIDUE)cContainerWithin(aAtomBondedNeighbor(aAtom,i)) != rRes )
	    return((RESIDUE)cContainerWithin(aAtomBondedNeighbor(aAtom,i)));
    }
    return(NULL);
}




/*
 *	sResidueTypeNameFromChar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a pointer to a string that defines
 *	the name of the RESIDUE type represented by the
 *	type character.
 */
char *
sResidueTypeNameFromChar( char c )
{
    switch ( c ) {
	case RESTYPEUNDEFINED	: return("undefined");
	case RESTYPESOLVENT	: return("solvent");
	case RESTYPEPROTEIN	: return("protein");
	case RESTYPENUCLEIC	: return("nucleic");
	case RESTYPESACCHARIDE	: return("saccharide");
	default:
	    return("UNKNOWN RESIDUE TYPE NAME!");
    }
    return("UNKNOWN RESIDUE TYPE NAME!");	/* for lint */
}



/*
 *	ResidueYouAreBeingRemoved
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	This routine handles the message to the RESIDUE
 *	that it is being removed from the CONTAINER that
 *	contains it.
 *	Currently there are no internal pointers between 
 *	RESIDUES to reset, so just pass the message up
 *	the hierarchy that this RESIDUE is being
 *	removed from something above it.
 */
void
ResidueYouAreBeingRemoved( RESIDUE rRes )
{
CONTAINER	cCont;
LOOP		lContents;

	/* Tell everything above that I am being removed */


    ContainerIAmBeingRemoved( cContainerWithin(rRes), (CONTAINER)rRes );

	/* Tell everything above that everything within me is */
	/* being removed */

    lContents = lLoop( (OBJEKT)rRes, CONTAINERS );
    while ( (cCont = (CONTAINER)oNext(&lContents)) )
	ContainerIAmBeingRemoved( cContainerWithin(rRes), cCont );

}



/*
 *	ResidueIAmBeingRemoved
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	If the OBJEKT being removed from this RESIDUE is an
 *	ATOM then check if it is one of the RESIDUEs connect
 *	atoms and set the connect pointer to NULL if this is
 *	true.  Likewise for the Imaging ATOM.
 *	Then pass the IAmBeingRemoved message up the
 *	hierarchy.
 */
void
ResidueIAmBeingRemoved( RESIDUE rRes, CONTAINER cRemoved )
{
STRING		sTemp;
int		i;

    if ( iObjectType(cRemoved) == ATOMid ) {
	for ( i=0; i<MAXCONNECT; i++ ) {
	    if ( rRes->aaConnect[i] == (OBJEKT)cRemoved ) {
		rRes->aaConnect[i] = NULL;
		VP1(( "Reseting connect#%d in %s\n", i,
			sContainerFullDescriptor( (CONTAINER)rRes, sTemp ) ));
	    }
	}
	if ( aResidueImagingAtom(rRes) == (ATOM)cRemoved ) {
	    ResidueSetImagingAtom( rRes, NULL );
	}
    } else {
	DFATAL(( "The only thing that can be removed is ATOMs" ));
    }

	/* Pass the message up the hierarchy */

    ContainerIAmBeingRemoved( cContainerWithin(rRes), cRemoved );
}




/*
 *	ResidueCrossLink
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create a bond between the two connection ATOMs.
 */
BOOL
bResidueCrossLink( RESIDUE rA, int iConnectA, 
	RESIDUE rB, int iConnectB, int iOrder )
{
ATOM		aA, aB;

    aA = aResidueConnectAtom( rA, iConnectA );
    if ( aA == NULL )
	VP0(( " %s: no such connect atom\n", rA->sDescription ));
    aB = aResidueConnectAtom( rB, iConnectB );
    if ( aB == NULL )
	VP0(( " %s: no such connect atom\n", rB->sDescription ));
    if ( aA == NULL || aB == NULL ) {
	return(FALSE);
    }
    AtomBondToOrder( aA, aB, iOrder );
    return(TRUE);
}



/*
 *===================================================================
 *
 *	Set RESIDUE attributes.
 */

/*
 *      zSetResidueConnect
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Make sure that the ATOM is within the RESIDUE and
 *      then set the connect ATOM, otherwise print an error.
 */
static void
zSetResidueConnect( RESIDUE rRes, int iConnect, ATOM aAtom )
{
STRING          sTemp, sTemp1;

    if (iObjectType(rRes) != RESIDUEid) {
	VP0(("SetResidueConnect: not a residue\n"));
	return;
    }
    if ( aAtom != NULL ) {
        if (iObjectType(aAtom) != ATOMid) {
	    VP0(("SetResidueConnect: not an atom\n"));
	    return;
	}
        if ( !bContainerContainedBy( (CONTAINER)aAtom, (CONTAINER)rRes ) ) {
            VP0(( "%s must contain %s.\n",
                sContainerFullDescriptor( (CONTAINER)rRes, sTemp ),
                sContainerFullDescriptor( (CONTAINER)aAtom, sTemp1 ) ));
            return;
        }
    }
    ResidueSetConnectAtom( rRes, iConnect, aAtom );
}





/*
 *	zSetResidueType
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Set the type of residue.
 */
static void
zSetResidueType( RESIDUE rRes, OBJEKT oObj )
{
STRING		sType;

    if ( iObjectType(oObj) != OSTRINGid ) {
	VP0(( "The residue type must be a string.\n" ));
	return;
    }
    strcpy( sType, sOString(oObj) );
    if ( strcmp( sType, sResidueTypeNameFromChar(RESTYPEUNDEFINED) ) == 0 ) {
	ResidueSetType( rRes, RESTYPEUNDEFINED );
    } else if ( strcmp( sType, sResidueTypeNameFromChar(RESTYPESOLVENT) )==0 ){
	ResidueSetType( rRes, RESTYPESOLVENT );
    } else if ( strcmp(sType,sResidueTypeNameFromChar(RESTYPENUCLEIC) )==0 ) {
	ResidueSetType( rRes, RESTYPENUCLEIC );
    } else if ( strcmp(sType,sResidueTypeNameFromChar(RESTYPEPROTEIN) )==0 ) {
	ResidueSetType( rRes, RESTYPEPROTEIN );
    } else if ( strcmp(sType,sResidueTypeNameFromChar(RESTYPESACCHARIDE) )==0){
	ResidueSetType( rRes, RESTYPESACCHARIDE );
    } else {
	VP0(( "Illegal residue type (%s).\n", sType ));
    }
}





/*
 *	ResidueSetAttribute
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Set the named attribute of the RESIDUE.
 *	This routine is used mainly by the command line interface
 *	to allow the user to change the attributes of RESIDUEs by
 *	name.
 */
void
ResidueSetAttribute( RESIDUE rRes, STRING sAttr, OBJEKT oAttr )
{
int		iConnect;

    if ( (iConnect = iResidueConnectFromName(sAttr)) != NOEND ) {
	zSetResidueConnect( rRes, iConnect, (ATOM)oAttr );
	goto DONE;
    } else if ( strcmp( sAttr, "restype" ) == 0 ) {
	zSetResidueType( rRes, oAttr );
	goto DONE;
    } else if ( strcmp( sAttr, "name" ) == 0 ) {
	if ( !bObjektWarnType( oAttr, OSTRINGid ) ) return;
	ContainerSetName( rRes, sOString(oAttr) );
	goto DONE;
    } else if ( strcmp( sAttr, "imagingAtom" ) == 0 ) {
	if ( oAttr != NULL ) {
	    if ( !bObjektWarnType( oAttr, ATOMid ) ) return;
	}
	ResidueSetImagingAtom( rRes, oAttr );
	goto DONE;
    }
    VP0(( "%s: non-existent attribute for a residue.\n", sAttr ));
    VP0(( "\tResidue attributes: restype, name, imagingAtom\n" ));
    return;

DONE:
    CDU(rRes);

}

/*
 *  OrderResidues() - make sure that if residues are in separate
 *	units, the orders are ok.
TODO - maybe leave this until leap is rewritten to use
hierarchical structure for molecules so that residue parents
can be checked to see if residues are already connected
void	OrderResidues( rRes1, rRes2 )
RESIDUE		rRes1, rRes2;
{


}
 */
