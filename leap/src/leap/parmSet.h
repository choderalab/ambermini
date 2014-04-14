/*
 *      File: parmSet.h
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
 *              PARMSET
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              A PARMSET is a repository of force field parameters.
 *
 *
 */
 
#ifndef PARMSET_H
#define PARMSET_H

# include	"hash.h"
# include	"database.h"


/*
 *-----------------------------------------------------------------------
 *
 *        Define object typedefs here.
 *        
 *        Object typedef MUST include the superclass object as
 *        its first structure element.
 */

typedef struct {
	OBJEKTt         oSuper;
	STRING		sFname;
	VARARRAY        vaAtoms;
	VARARRAY        vaBonds;
	VARARRAY        vaAngles;
	VARARRAY        vaTorsions;
	VARARRAY        vaImpropers;
	VARARRAY        vaHBonds;
        VARARRAY        vaNBEdits;
	BOOL		bBeingEdited;
} PARMSETt;

typedef PARMSETt	*PARMSET;


/*
 *  this one param type is moved here from parmSet.c
 *	to allow easy use in unitio.c
 */
#define MAXTYPELEN      5
typedef char    typeStr[MAXTYPELEN];

#define	MAXORDERLEN	5
typedef	char	orderStr[MAXORDERLEN];

#define PARMDESCRIPTIONLEN      32
typedef char            DESCRIPTION[PARMDESCRIPTIONLEN];

typedef struct  {
	typeStr         sType1;
	typeStr         sType2;
	typeStr         sType3;
	typeStr         sType4;
	int		iType;
	int		iN;
	double          dKp;
	double          dP0;
	double          dScEE;
	double          dScNB;
	orderStr	sOrder;		/* for impropers */
	DESCRIPTION     sDesc;
} TORSIONPARMt;


typedef	VARARRAY	TORSION;


#define	WILD_CARD_TYPE		"?"
#define	WILD_CARD_TYPE_CHAR	'?'



/*
 *======================================================================
 *
 *        Define object messages here.
 *        
 *        There must be at least a Create, Destroy, and Describe message.
 *        Hook into the messages of the superclasses so that
 *        when the message is sent to the most primative superclass
 *        of this class that it will eventually make it into these routines.
 */

#define	PARM_IGNORE		-2
#define	PARM_STOP_LOOKING	-3
#define	PARM_NOT_INDEX		-4
#define	PARM_FOUND_TERMS	-5
#define PARM_NOT_FOUND		-6
#define PARM_FOUND_WILD		-7
#define PARM_FOUND_EXACT	-8

/*
 *	The following routines are used to create/destroy/describe
 *	PARMSETs
 */
 
extern PARMSET	psParmSetCreate();
extern PARMSET	psParmSetDuplicate(PARMSET psOld);
extern void	ParmSetDestroy(PARMSET *psPLib);
extern void	ParmSetDescribe(PARMSET psLib);


/*
 *	The following routines add new parameters to the
 *	PARMSET.
 *
 *	Returns the index where the parameter was added.
 */

extern int	iParmSetAddAtom(PARMSET psLib, char *sType, 
			double dMass, double dPolar, 
			double dEpsilon, double dR, 
			double dEpsilon14, double dR14, 
			double dScreenF,
			int iElement, int iHybridization, char *sDesc);
extern int	iParmSetAddBond(PARMSET psLib, char *sType1, char *sType2,
			double dKb, double dR0, char *sDesc);
extern int	iParmSetAddAngle(PARMSET psLib, 
			char *sType1, char *sType2, char *sType3,
			double dKt, double dT0, double dTkub, double dRkub, char *sDesc);
extern int	iParmSetAddProperTerm(PARMSET psLib,
			char *sType1, char *sType2, char *sType3, char *sType4,
			int iN, double dKp, double dP0, double dScEE, double dScNB, char *sDesc);
extern int	iParmSetAddImproperTerm(PARMSET psLib,
			char *sType1, char *sType2, char *sType3, char *sType4,
			int iN, double dKp, double dP0, double dScEE, double dScNB, char *sDesc);
extern int	iParmSetAddHBond(PARMSET psLib, char *sType1, char *sType2,
			double dA, double dB, char *sDesc);
extern int      iParmSetAddNBEdit( PARMSET psLib, char *sType1, char *sType2,
				   double dEI, double dEJ, double dRI,
                                   double dRJ, char *sDesc );



/*
 *	Find the parameters within a PARMSET and return the index
 *	into the PARMSET for that parameter.
 *	For torsion parameters 
 *	a TORSION value is returned that is used with
 *	the ParmSetTorsion routines to obtain the values
 *	of the Torsion parameters.
 *	The routine ParmSetTorsion returns the index, types, and constants
 *	for the (i)th torsion parameter in the TORSION type.
 *	The index refers into the ParmSet from which the torsion
 *	was obtained.  Once the caller merges two TORSIONs from two
 *	different PARMSETs then the index is meaningless.
 */
 
extern int	iParmSetFindAtom(PARMSET psLib, char *sType);
extern int	iParmSetFindBond(PARMSET psLib, char *sType1, char *sType2);
extern int	iParmSetFindAngle(PARMSET psLib, 
			char *sType1, char *sType2, char *sType3);
extern int	iParmSetFindProperTerms(PARMSET psLib, TORSION tTorsion, 
			BOOL bUseIndex, 
			char *sType1, char *sType2, char *sType3, char *sType4);
extern int	iParmSetFindImproperTerms(PARMSET psLib, TORSION tTorsion, 
			BOOL bUseIndex,
			char *sType1, char *sType2, char *sType3, char *sType4);
extern int	iParmSetFindHBond(PARMSET psLib, char *sType1, char *sType2);
extern int	iParmSetFindNBEdit(PARMSET psLib, char *sType1, char *sType2);

/*
 *	TORSION routines
 *
 *	A TORSION is an object that maintains a list of Fourier components
 *	that make up a TORSION.  Each Fourier component is called a Term.
 *
 *	PROPER torsions are 'normal' interatomic potentials along bonds;
 *	IMPROPER ones are somewhat ad hoc contrivances to enforce planarity
 */

#define PROPER		0
#define IMPROPER	1

extern TORSION	tParmSetTORSIONCreate();

#define iParmSetTORSIONTermCount( tTorsion ) \
		( iVarArrayElementCount( (tTorsion) ) )
#define ParmSetTORSIONDestroy( tPT ) \
		{ VarArrayDestroy( (tPT) );}

extern void	ParmSetTORSIONTerm(TORSION tTorsion, int iTorsionIndex, 
			int *iPParmSetIndex,
			char *cPTyp1, char *cPTyp2, char *cPTyp3, char *cPTyp4,
			int *iPN, double *dPKp, double *dPP0, double *dPScEE,
			double *dPScNB, char *sDesc );
extern BOOL	bParmSetTORSIONAddProperTerm(TORSION tTorsion,
			char *cPType1, char *cPType2, 
			char *cPType3, char *cPType4,
			int iN, double dKp, double dP0, double dScEE,
			double dScNB, char *sDesc);
extern BOOL	bParmSetTORSIONAddImproperTerm(TORSION tTorsion,
			char *cPType1, char *cPType2, 
			char *cPType3, char *cPType4,
			int iN, double dKp, double dP0, double dScEE,
                        double dScNB, char *sDesc);
extern void	ParmSetTORSIONOrderAtoms();	
extern void	ParmSetImproperOrderAtoms( TORSION tTorsion, int iTorsionIndex,
			char *cPaTypes[4], int iaIndexes[4] );
extern BOOL	bParmSetCapableOfHBonding( PARMSET psParms, char *sType );

/*
 *	PARMSET information routines.
 */

#define sParmName( psParmSet )	(psParmSet)->sFname 

#define iParmSetTotalAtomParms( psParmSet ) \
			iVarArrayElementCount( (psParmSet)->vaAtoms )
#define iParmSetTotalBondParms( psParmSet ) \
			iVarArrayElementCount( (psParmSet)->vaBonds )
#define iParmSetTotalAngleParms( psParmSet ) \
			iVarArrayElementCount( (psParmSet)->vaAngles )
#define iParmSetTotalTorsionParms( psParmSet ) \
			iVarArrayElementCount( (psParmSet)->vaTorsions )
#define iParmSetTotalImproperParms( psParmSet ) \
			iVarArrayElementCount( (psParmSet)->vaImpropers )
#define iParmSetTotalHBondParms( psParmSet ) \
			iVarArrayElementCount( (psParmSet)->vaHBonds )
#define iParmSetTotalNBEdits( psParmSet ) \
			iVarArrayElementCount( (psParmSet)->vaNBEdits )

extern  BOOL    bParmSetCapableofHBonding();    /* ( PARMSET, char* ) */


/*
 *	Routines used to actually obtain PARMSET parameters
 *	Use the index obtained from iParmSetFindxxxxx
 */
extern void    	ParmSetAtom(PARMSET psLib, int i, char *sType, double *dPMass,
			double *dPPolar, double *dPEpsilon, double *dPR,
            double *dPEpsilon14, double *dPR14, double *dPScreenF,
			int *iPElement, int *iPHybridization, char *sDesc);
extern void	ParmSetBond(PARMSET psLib, int i, char *sType1, char *sType2,
			double *dPKb, double *dPR0, char *sDesc);
extern void	ParmSetAngle(PARMSET psLib, int i, 
			char *sType1, char *sType2, char *sType3,
			double *dPKt, double *dPT0, double *dPTkub, double *dPRkub,
			char *sDesc);
extern void	ParmSetTorsion(PARMSET psLib, int i,
			char *sType1, char *sType2, char *sType3, char *sType4,
			int *iPN, double *dPKp, double *dPP0, double *dPScEE,
			double *dPScNB, char *sDesc);
extern void	ParmSetImproper(PARMSET psLib, int i,
			char *sType1, char *sType2, char *sType3, char *sType4,
			int *iPN, double *dPKp, double *dPP0, char *sDesc);
extern void	ParmSetHBond(PARMSET psLib, int i, char *sType1, char *sType2,
			double *dPA, double *dPB, char *sDesc);
extern void     ParmSetNBEdit(PARMSET psLib, int i, char *sType1, char *sType2,
			      double *dEI, double *dEJ, double *dRI,
			      double *dRJ, char *sDesc);


/*	Davids Changes */
/*
 *	Routines used to update PARMSET parameters
 *	Use the idex obtained from iParmSetFindxxxxx
 *	Passing NULL ( 0 ) for any of the parameter variables
 *	will prevent that variable from being changed.
 */

extern void	ParmSetUpdateAtom(PARMSET psLib, int i, char *sType,
			double *dPMass, double *dPPolar, 
			double *dPEpsilon, double *dPR, double *dPScreenF,
			int *iPElement, int *iPHybrid, char *sDescription);
extern void	ParmSetUpdateBond( PARMSET psLib, int i, 
			char *sType1, char *sType2,
			double *dPKb, double *dPR0, char *sDescription);
extern void	ParmSetUpdateAngle(PARMSET psLib, int i,
			char *sType1, char *sType2, char *sType3,
			double *dPKt, double *dPT0, char *sDescription);
extern void	ParmSetUpdateTorsion(PARMSET psLib, int i,
			char *sType1, char *sType2, char *sType3, char *sType4,
			int *iPN, double *dPKp, double *dPP0, double *dPScEE,
			double *dPScNB, char *sDescription);
extern void	ParmSetUpdateImproper(PARMSET psLib, int i,
			char *sType1, char *sType2, char *sType3, char *sType4,
			int *iPN, double *dPKp, double *dPP0, double *dScEE, double *dScNB,
			char *sDescription);
extern void	ParmSetUpdateHBond(PARMSET psLib, int i, 
			char *sType1, char *sType2,
			double *dPA, double *dPB, char *sDescription);

extern void	ParmSetNewAtoms(PARMSET psParmSet, int iCount);
extern void	ParmSetNewBonds(PARMSET psParmSet, int iCount);
extern void	ParmSetNewAngles(PARMSET psParmSet, int iCount);
extern void	ParmSetNewTorsions(PARMSET psParmSet, int iCount);
extern void	ParmSetNewImpropers(PARMSET psParmSet, int iCount);
extern void	ParmSetNewHBonds(PARMSET psParmSet, int iCount);
extern int	iParmSetProperCount( PARMSET psParmSet );
extern int	iParmSetImproperCount( PARMSET psParmSet );


/*	End of Davids Changes */

/*
 *	Read and write PARMSETs from/to DATABASEs
 */
extern PARMSET	psParmSetLoad(DATABASE db);
extern void	ParmSetSave(PARMSET psLib, DATABASE db);

        /* iParmSetTorsionGenerality returns a number from 0 to 3 which */
        /* represents the generality of the torsion parameters */
        /* Proper torsions can be either 0 or 2 */
        /* Improper torsions can range from 0 to 3 */
        /* The third atom in an Improper torsion is the central atom */
        /* it should NEVER be WILD_CARD_TYPE */

#define iParmSetTorsionGenerality( a1,a2,a3,a4 ) ( \
	( ( a1[0] == WILD_CARD_TYPE_CHAR ) ? 1 : 0 ) + \
	( ( a2[0] == WILD_CARD_TYPE_CHAR ) ? 1 : 0 ) + \
	( ( a4[0] == WILD_CARD_TYPE_CHAR ) ? 1 : 0 ) )

/*	Davids Changes */

#define ParmSetEditing( parmSet, state ) ( parmSet->bBeingEdited = state )
#define bParmSetBeingEdited( parmSet )	( parmSet->bBeingEdited )
					
/*	End of Davids Changes */

#endif /* PARMSET_H */
