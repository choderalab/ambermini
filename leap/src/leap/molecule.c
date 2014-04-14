/*
 *      File: molecule.c
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
 *              MOLECULE
 *      Superclass:
 *              CONTAINER, LOOP
 *
 *      Description:
 *              A MOLECULE is a subclass of CONTAINER.
 *              MOLECULES can contain residues, and atoms.
 *              The property that destinguishes MOLECULES
 *              from all other CONTAINER subclasses is that
 *              the graph defined by the bonds of the atoms
 *              within the MOLECULE is a SINGLE CLOSED graph.
 *              Functions must be called to ensure that this
 *              property is maintained.
 *
 */
 



#include	"basics.h"

#include        "classes.h"



/*
 *      mMoleculeCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a molecule.
 *
 *      Return:
 *              Return a pointer to the MOLECULE created.
 */
MOLECULE
mMoleculeCreate()
{
MOLECULE        m;

    MALLOC( m, MOLECULE, sizeof(MOLECULEt) );
    return(m);
}







/*
 *      MoleculeDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the MOLECULE and all of the residues/atoms within it.
 *
 *      Arguments:
 *      oObject -       A pointer to the object.
 */
void    
MoleculeDestroy( MOLECULE *mPMolecule )
{
    FREE( *mPMolecule );
    *mPMolecule = NULL;
}






/*
 *      MoleculeDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to return a descriptor
 *      string for the MOLECULE object.
 *      It is the callers responsibility to ensure that there is room
 *      enough in the buffer to describe the object.
 *
 *      Arguments:
 *      oObject -       The object to describe.
 *
 */
void    
MoleculeDescribe( MOLECULE mMolecule )
{
LOOP    lContents;
OBJEKT  oObj;
STRING  sTemp;

    VP0(( "MOLECULE name: %s\n", mMolecule->cHeader.sName ));
    VP0(( "MOLECULE sequence number: %d\n", 
                        iContainerSequence(mMolecule) ));
    VP0(( "Contents: \n" ));

    BasicsResetInterrupt();
    lContents = lLoop( (OBJEKT)mMolecule, DIRECTCONTENTS );
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
 *      mMoleculeDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a duplicate of the molecule.
 *	Should only be called from cContainerDuplicate()
 *      which is expected to duplicate the CONTAINER variables.
 */
MOLECULE
mMoleculeDuplicate( MOLECULE mOld )
{
MOLECULE	mNew;

    MALLOC( mNew, MOLECULE, sizeof(MOLECULEt) );
    memcpy( mNew, mOld, sizeof(MOLECULEt) );
    return(mNew);
}





/*
 *      MoleculeResetPointers
 *
 *	Author:	Christian Schafmeister (1991)
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
MoleculeResetPointers( MOLECULE mMol )
{
        /* Currently there are no internal pointers */
}


/*
 *      MoleculeCheck
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Check the MOLECULE and its contents for whether or
 *      not calculations can be run on it.
 *
 *      Arguments:
 *              mMol -          Molecule to check.
 *              iPErrors -      Add the number of errors found.
 *              iPWarnings -    Add the number of warnings found.
 */
void    
MoleculeCheck( MOLECULE mMol, int *iPErrors, int *iPWarnings )
{
LOOP            lContents;
CONTAINER       cCont;

        /* I can't think of anything to check right now */

    lContents = lLoop( (OBJEKT)mMol, DIRECTCONTENTS );
    while ( (cCont = (CONTAINER)oNext(&lContents)) ) {
        ContainerCheck( cCont, iPErrors, iPWarnings );
    }
}




/*
 *	MoleculeYouAreBeingRemoved
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	There are no pointers between MOLECULEs that have
 *	to be maintained, so just pass the message up the
 *	hierarchy.
 */
void	
MoleculeYouAreBeingRemoved( MOLECULE mMol )
{
LOOP		lContents;
CONTAINER	cCont;

	/* Pass the message on up */

    ContainerIAmBeingRemoved( cContainerWithin(mMol), (CONTAINER)mMol );

	/* Tell everything above me that everything below is */
	/* also being removed */
    lContents = lLoop( (OBJEKT)mMol, CONTAINERS );
    while ( (cCont = (CONTAINER)oNext(&lContents)) )
	ContainerIAmBeingRemoved( cContainerWithin(mMol), cCont );

}




/*
 *	MoleculeIAmBeingRemoved
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	There are no pointers from MOLECULEs to things within
 *	it, so just pass the IAmBeingRemoved message up
 *	the hierarchy.
 */
void	
MoleculeIAmBeingRemoved( MOLECULE mMol, CONTAINER cRemoved )
{
	/* Pass the message on up */

    ContainerIAmBeingRemoved( cContainerWithin(mMol), cRemoved );
}





/*
 *	MoleculeSetAttribute
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Change the attribute of the MOLECULE.
 */
void
MoleculeSetAttribute( MOLECULE mMol, STRING sAttr, OBJEKT oAttr )
{
    if ( strcmp( sAttr, "name" ) == 0 ) {
	if ( !bObjektWarnType( oAttr, OSTRINGid ) ) return;
	ContainerSetName( mMol, sOString(oAttr) );
	return;
    }
    VP0(( "%s: non-existent attribute for a molecule.\n", sAttr ));
    VP0(( "\tMolecule attributes: name\n" ));
    
}


