/*
 *      File: residue.h
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
 *              CONTAINER
 *
 *      Description:
 *
 *              Residues can contain atoms
 */
 
#ifndef RESIDUE_H
#define RESIDUE_H

/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/


                /* Some containers, UNITs and RESIDUEs have connect */
                /* pointers which point to atoms to which connections */
                /* can be made easily.  Other atoms can be connected to */
                /* but they have to be refered to directly */
#define MAXCONNECT      6

#define NOEND           -1

#define	CONNECT0	0
#define	CONNECT1	1
#define	CONNECT2	2
#define	CONNECT3	3
#define	CONNECT4	4
#define	CONNECT5	5


#define FIRSTEND        CONNECT0	/* This is the end that is bare in */
                                        /* the first RESIDUE */
#define LASTEND         CONNECT1	/* This is the end that is bare in */
                                        /* the last RESIDUE */

#define NEND            FIRSTEND        /* N terminus */
#define CEND            LASTEND         /* C terminus */
#define SEND            CONNECT2	/* Disulphide bridges */

		/* Residue types */

#define	RESTYPEUNDEFINED	'?'
#define	RESTYPESOLVENT		'w'
#define	RESTYPEPROTEIN		'p'
#define	RESTYPENUCLEIC		'n'
#define	RESTYPESACCHARIDE	's'

#define	SEGIDLEN		(6+1)
typedef	char		resSEGIDt[SEGIDLEN];

                /* Residue flags */

#define	RESIDUEPERM	0x0000FFFF
#define RESIDUEUNKNOWN  0x00000001

#define	RESIDUETEMP	0xFFFF0000
#define	RESIDUEINCAP	0x00010000

typedef struct {
	char	sName1[NAMELEN];
	char	sName2[NAMELEN];
	char	sName3[NAMELEN];
	char	sName4[NAMELEN];
} IMPROPERt;

typedef struct  {
	CONTAINERt      cHeader;
	FLAGS           fFlags;
	char		cResType;
	STRING          sDescription;
	OBJEKT          aaConnect[MAXCONNECT];
	OBJEKT		aSolventImagingAtom;	/* ATOM */
	int		iPdbSequenceNumber;
	char		cPdbChain;
	resSEGIDt	siPdbSegId;
	VARARRAY	vaImpropers;
	double		dTemp;
	int		iTemp;
} RESIDUEt;

typedef RESIDUEt	*RESIDUE;


typedef struct {
    STRING sName;
    int iSequenceNumber;
    int iaConnectIndex[MAXCONNECT];
    int iNextChildSequence;
    int iPdbSequenceNumber;
    int iAtomStartIndex;
    int iImagingAtomIndex;
    char sResidueType[3];
    RESIDUE rResidue;
} SAVERESIDUEt;



/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */

extern RESIDUE		rResidueCreate();
extern void		ResidueDelete(RESIDUE *rPResidue);
extern void		ResidueDescribe(RESIDUE rResidue);
extern void             ResidueDestroy(RESIDUE *rPResidue);
extern RESIDUE		rResidueDuplicate(RESIDUE rOld);
extern void		ResidueResetPointers(RESIDUE rRes);
extern void		ResidueCheck(RESIDUE rRes, 
				int *iPErrors, int *iPWarnings);
extern void		ResidueMutate(RESIDUE rNew, RESIDUE rOld);
extern char		*sResidueTypeNameFromChar(char c);

extern RESIDUE		rResidueConnected(RESIDUE rRes, int iConnect);
extern int		iResidueConnectFromName(char *sName);
extern void		ResidueSetAttribute( RESIDUE rRes, 
				STRING sAttr, OBJEKT oAttr );

extern BOOL	bResidueCrossLink(RESIDUE rA, int iConnectA,
			RESIDUE rB, int iConnectB, int iOrder);
extern void	ResidueYouAreBeingRemoved(RESIDUE rRes);
extern void	ResidueIAmBeingRemoved(RESIDUE rRes, CONTAINER cRemoved);


#define bResidueConnectUsed(r,c)     (((RESIDUE)r)->aaConnect[c]!=NULL)
#define ResidueSetConnectAtom(r,c,a) (((RESIDUE)r)->aaConnect[c]=(OBJEKT)(a),\
					CDU(r))
#define aResidueConnectAtom(r,c)        (ATOM)(((RESIDUE)r)->aaConnect[c])
#define ResidueSetDescription(r,s) (strcpy( ((RESIDUE)(r))->sDescription,s),\
					CDU(r))
#define sResidueDescription(r)          (((RESIDUE)(r))->sDescription)
#define bResidueFlagsSet(r,f)           ((((RESIDUE)(r))->fFlags & f)!= 0)
#define ResidueSetFlags(r,f)            (((RESIDUE)(r))->fFlags |= f,CDU(r) )
#define ResidueDefineFlags(r,f)  (((RESIDUE)(r))->fFlags = f,CDU(r))
#define ResidueResetFlags(r,f)    (((RESIDUE)(r))->fFlags &= ~f,CDU(r))
#define	ResidueSetType(r,c)	(((RESIDUE)(r))->cResType = (c),CDU(r))
#define	cResidueType(r)		(((RESIDUE)(r))->cResType)
#define ResidueSetImagingAtom(r,a) (((RESIDUE)r)->aSolventImagingAtom=(OBJEKT)(a),CDU(r))
#define aResidueImagingAtom(r)     (ATOM)(((RESIDUE)r)->aSolventImagingAtom)
#define	iResiduePdbSequence(r)	(((RESIDUE)r)->iPdbSequenceNumber)
#define	ResidueSetPdbSequence(r,i)	(((RESIDUE)r)->iPdbSequenceNumber=(i))



#endif /* RESIDUE_H */
