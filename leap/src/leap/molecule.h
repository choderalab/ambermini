/*
 *      File: molecule.h
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
 *              CONTAINER
 *
 *      Description:
 *
 *              Molecules can contain residues and/or atom.s
 *              Molecules contain a single closed graph of
 *              atoms.
 */
 
#ifndef MOLECULE_H
#define MOLECULE_H

/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/


typedef struct  {
	CONTAINERt      cHeader;
} MOLECULEt;

typedef MOLECULEt	*MOLECULE;





/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */

extern MOLECULE		mMoleculeCreate();
extern void		MoleculeDelete(MOLECULE *mPMolecule);
extern void		MoleculeDescribe(MOLECULE mMolecule);
extern void             MoleculeDestroy(MOLECULE *mPMolecule);
extern MOLECULE		mMoleculeDuplicate(MOLECULE mOld);
extern void		MoleculeResetPointers(MOLECULE mMol);
extern void		MoleculeCheck(MOLECULE mMol, int *iPErrors, int *iPWarnings);

extern void		MoleculeYouAreBeingRemoved(MOLECULE mMol);
extern void		MoleculeIAmBeingRemoved(MOLECULE mMol, CONTAINER cRemoved);


extern void		MoleculeSetAttribute(MOLECULE mMol, STRING sAttr, OBJEKT oAttr);

#endif  /* MOLECULE_H */
