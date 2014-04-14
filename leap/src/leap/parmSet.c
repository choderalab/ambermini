/*
 *      File: parmSet.c
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
 *              Torsional parameters include both the proper and improper
 *              torsions, with a flag that indicates which are which.
 *
 *		A TORSION is an object that contains several torsional
 *		terms each of which cooresponds to a Fourier term
 *		of the torsion.
 *		A TORSION contains the terms in ascending order by
 *		multiplicity.
 *
 *		Entries in the bond, angle, proper/improper, and hbond
 *		VARARRAYs have their atom types sorted for quicker
 *		lookup.  THIS MUST BE ADHERED TOO.  All of the
 *		search routines REQUIRE pre-ordering.
 */



#include	"basics.h"

#include        "classes.h"

#include        "dictionary.h"
#include        "database.h"


#include	"amber.h"



typedef struct  {
                typeStr         sType;
                double          dMass;
		double		dPolar;
                double          dEpsilon;
                double          dR;
                double          dEpsilon14;
                double          dR14;
		double		dScreenF;
		int		iElement;
		int		iHybridization;
                DESCRIPTION     sDesc;
} ATOMPARMt;

typedef struct  {
                typeStr         sType1;
                typeStr         sType2;
                double          dKb;
                double          dR0;
                DESCRIPTION     sDesc;
} BONDPARMt;

typedef struct  {
                typeStr         sType1;
                typeStr         sType2;
                typeStr         sType3;
                double          dKt;
                double          dT0;
                double          dTkub;
                double          dRkub;
                DESCRIPTION     sDesc;
} ANGLEPARMt;


typedef struct  {
                typeStr         sType1;
                typeStr         sType2;
                double          dA;
                double          dB;
                DESCRIPTION     sDesc;
} HBONDPARMt;


typedef struct  {
                typeStr         sType1;
                typeStr         sType2;
                double          dEI;
                double          dEJ;
                double          dRI;
                double          dRJ;
                double          dA;
                double          dC;
                DESCRIPTION     sDesc;
} NBEDITt;



	/* TORSION_MATCHt are stored inside TORSIONs */
typedef	struct	{
		int		iIndex;
		TORSIONPARMt	tpTorsion;
} TORSION_MATCHt;


typedef struct CMNT_t {
	char *record;
	struct CMNT_t *next;
} CMNTt;

typedef struct {
                STRING          title;
		char            *reslist[8];
		int             nres;
		CMNTt           *cmnt;
		int		resolution;
		double          *map;
} CMAPt;

/*
 *-------------------------------------------------------------------
 *
 *        Define static variables here.
 */




/*
 *===================================================================
 *
 *       Define private routines here.
 */



/*
 *	zParmSetOrderBondAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Order the bond atom names into alphabetical order
 *	to speed up searches.
 */
static void
zParmSetOrderBondAtoms( char *sAtom1, char *sAtom2 )
{
STRING		sTemp;

    if ( strcmp( sAtom1, sAtom2 ) > 0 ) {
	SWAP_STRINGS( sAtom1, sAtom2, sTemp );
    }
}




/*
 *	zParmSetOrderAngleAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Order the angle atom types into alphabetical order
 *	to speed up searches.
 *
 *	The proper order is for the first atom type to
 *	be less than the third atom type.
 */
static void
zParmSetOrderAngleAtoms( char *sAtom1, char *sAtom2, char *sAtom3 )
{
STRING		sTemp;

    if ( strcmp( sAtom1, sAtom3 ) > 0 ) {
	SWAP_STRINGS( sAtom1, sAtom3, sTemp );
    }
}






/*
 *	zParmSetOrderTorsionAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Order the torsion atom types into alphabetical order
 *	to speed up searches.
 *
 *	The proper order is for the first type to be less than
 *	the fourth type.
 *	If the first and fourth are the same then the
 *	second type must be less than the third type.
 */
void
zParmSetOrderTorsionAtoms( char *sAtom1, char *sAtom2, 
		char *sAtom3, char *sAtom4 )
 {
STRING		sTemp;
int		iCmp14;

    iCmp14 = strcmp( sAtom1, sAtom4 );
    if ( iCmp14 == 0 ) {
	/* 
	 *  end types are the same, so can order inner pair
	 */
        zParmSetOrderBondAtoms( sAtom2, sAtom3 );
    } else if ( iCmp14 > 0 ) {
	/*
	 *  reorder end-to-end: ABCD->DCBA
	 */
	SWAP_STRINGS( sAtom1, sAtom4, sTemp );
	SWAP_STRINGS( sAtom2, sAtom3, sTemp );
    }
}



/*
 *	zParmSetOrderImproperAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Order the improper atom _types_ into alphabetical order
 *	to speed up searches.
 *
 *	The proper order is for the wild cards to be first
 *	and the rest of the types to be in alphabetical order.
 *	Only the first, second, and fourth atoms are affected,
 *	the third atom is ALWAYS the central atom.
 */
static void
zParmSetOrderImproperAtoms( char *sAtom1, char *sAtom2, 
	char *sAtom3, char *sAtom4, char *sOrder )
{
STRING		sTemp;
STRING		saAtoms[3];
char		cTemp;
#define	zzSwapC(a,b)	{(cTemp)=(a);(a)=(b);(b)=(cTemp);}

    strcpy( saAtoms[0], sAtom1 != (char*)NULL ? sAtom1 : " ");
    strcpy( saAtoms[1], sAtom2 != (char*)NULL ? sAtom2 : " ");
    strcpy( saAtoms[2], sAtom4 != (char*)NULL ? sAtom4 : " ");

		/* Change wild card types to spaces so that */
		/* they are alphabetically before all other types */
    if ( strcmp( saAtoms[0], WILD_CARD_TYPE ) == 0 ) strcpy( saAtoms[0], " " );
    if ( strcmp( saAtoms[1], WILD_CARD_TYPE ) == 0 ) strcpy( saAtoms[1], " " );
    if ( strcmp( saAtoms[2], WILD_CARD_TYPE ) == 0 ) strcpy( saAtoms[2], " " );

#ifndef BILL_NEW
    if ( strcmp( saAtoms[0], saAtoms[1] ) > 0 ) {
	SWAP_STRINGS(saAtoms[0], saAtoms[1], sTemp );
	zzSwapC(sOrder[0],sOrder[1]);
    }
    if ( strcmp( saAtoms[1], saAtoms[2] ) > 0 ) {
	SWAP_STRINGS(saAtoms[1], saAtoms[2], sTemp );
	zzSwapC(sOrder[1],sOrder[3]);
    }
    if ( strcmp( saAtoms[0], saAtoms[1] ) > 0 ) {
	SWAP_STRINGS(saAtoms[0], saAtoms[1], sTemp );
	zzSwapC(sOrder[0],sOrder[1]);
    }
    if ( strcmp( saAtoms[1], saAtoms[2] ) > 0 ) {
	SWAP_STRINGS(saAtoms[1], saAtoms[2], sTemp );
	zzSwapC(sOrder[1],sOrder[3]);
    }
#else
    if ( strcmp( saAtoms[0], saAtoms[1] ) > 0 ) {
	SWAP_STRINGS(saAtoms[0], saAtoms[1], sTemp );
	zzSwapC(sOrder[0],sOrder[1]);
    }
    if ( strcmp( saAtoms[1], saAtoms[3] ) > 0 ) {
	SWAP_STRINGS(saAtoms[1], saAtoms[3], sTemp );
	zzSwapC(sOrder[1],sOrder[3]);
    }
    if ( strcmp( saAtoms[0], saAtoms[1] ) > 0 ) {
	SWAP_STRINGS(saAtoms[0], saAtoms[1], sTemp );
	zzSwapC(sOrder[0],sOrder[1]);
    }
#endif
		/* Change wild card characters back */

    if ( strcmp( saAtoms[0], " " ) == 0 ) strcpy( saAtoms[0], WILD_CARD_TYPE );
    if ( strcmp( saAtoms[1], " " ) == 0 ) strcpy( saAtoms[1], WILD_CARD_TYPE );
    if ( strcmp( saAtoms[2], " " ) == 0 ) strcpy( saAtoms[2], WILD_CARD_TYPE );

    strcpy( sAtom1, saAtoms[0] );
    strcpy( sAtom2, saAtoms[1] );
    strcpy( sAtom4, saAtoms[2] );
}



/*
 *---------------------------------------------------------------------
 *
 */



/*
 *      zbParmSetMatchTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the four types match the torsion.
 */
static BOOL
zbParmSetMatchTorsion( TORSIONPARMt *tpPTorsion, char *s1, char *s2, 
	char *s3, char *s4 )
{
BOOL            bFoundOne;

    bFoundOne = FALSE;
   
                /* Check the torsion only one way, */
		/* this relies on the types being ordered properly */
		/* by zParmSetOrderTorsionAtoms */

    if ( strcmp( tpPTorsion->sType1, WILD_CARD_TYPE )!=0 ) {
        if ( strcmp( tpPTorsion->sType1, s1 ) == 0 &&
	     strcmp( tpPTorsion->sType2, s2 ) == 0 &&
	     strcmp( tpPTorsion->sType3, s3 ) == 0 &&
	     strcmp( tpPTorsion->sType4, s4 ) == 0 ) bFoundOne = TRUE;
    } else {
        if ( strcmp( s2, s3 ) <= 0 ) {
	    if ( strcmp( tpPTorsion->sType2, s2 ) == 0 &&
	    	 strcmp( tpPTorsion->sType3, s3 ) == 0 ) bFoundOne = TRUE;
	} else {
	    if ( strcmp( tpPTorsion->sType3, s2 ) == 0 &&
	    	 strcmp( tpPTorsion->sType2, s3 ) == 0 ) bFoundOne = TRUE;
	}
    }
#ifdef	DEBUG
    if ( bFoundOne ) {
        MESSAGE(( "Matched torsion %s-%s-%s-%s to: %s-%s-%s-%s\n",
		    tpPTorsion->sType1,
		    tpPTorsion->sType2,
		    tpPTorsion->sType3,
		    tpPTorsion->sType4,
		    s1, s2, s3, s4 ));
    }
#endif
		    
    return(bFoundOne);
}



/*
 *	zParmSetAddToTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a torsion to the TORSION.
 *	Check if the proper torsion is already in the TORSION,
 *	if it is, but is a general torsion rather than a specific
 *	one then replace it.
 *
 *	If (bUseIndex) is TRUE then write the index of the term
 *	within the PARMSET into the TORSION.iIndex, otherwise write
 *	PARM_NOT_FOUND.
 *
 *	Return TRUE if the proper term was actually added.
 *
 *TODO: Fix this routine, it will not properly compare proper torsions
 *TODO: with wild cards.
 */
static int
zbParmSetAddToTorsion( TORSION tTorsion, int iIndex, TORSIONPARMt *tpPTorsion, 
	BOOL bUseIndex )
{
TORSION_MATCHt	*tmPCur;
int		i;
TORSION_MATCHt	tmNew;

    memset( &tmNew, 0, sizeof(tmNew) );		/* for Purify */
    if ( !iVarArrayElementCount( tTorsion ) ) {

	/*
	 *  just add
	 */
        if ( bUseIndex ) 
	    	tmNew.iIndex = iIndex;
	else
		tmNew.iIndex = PARM_NOT_FOUND;

	tmNew.tpTorsion = *tpPTorsion;
	VarArrayAdd( tTorsion, (GENP)&tmNew );
	if ( !strcmp( tmNew.tpTorsion.sType1, WILD_CARD_TYPE ) )
		return(PARM_FOUND_WILD);
	return(PARM_FOUND_EXACT);
    }

    /*
     *  if the previous torsion is general and this one is
     *	specific, delete the previous one(s) and add the
     *	new one
     */
    tmPCur = PVAI(tTorsion, TORSION_MATCHt, 0);
    if ( strcmp( tmPCur->tpTorsion.sType1, WILD_CARD_TYPE ) == 0   &&
         strcmp( tpPTorsion->sType1, WILD_CARD_TYPE ) != 0 ) {
	VarArraySetSize( tTorsion, 0 );
        if ( bUseIndex ) 
	    	tmNew.iIndex = iIndex;
	else
		tmNew.iIndex = PARM_NOT_FOUND;
	tmNew.tpTorsion = *tpPTorsion;
	VarArrayAdd( tTorsion, (GENP)&tmNew );
	return(PARM_FOUND_EXACT);
    } 
		/* First check if it is already in the TORSION */
		/* VARARRAY */

    for ( tmPCur = PVAI(tTorsion, TORSION_MATCHt, 0),i=0; 
    		i<iVarArrayElementCount(tTorsion); 
		i++, tmPCur++ ) {
	    if ( tmPCur->tpTorsion.iN == tpPTorsion->iN ) {
		    return(PARM_NOT_FOUND);
	    }
    }

		/* If no match was found then simply add the */
		/* proper to the end of the TORSION */
    if ( bUseIndex ) 
    	tmNew.iIndex = iIndex;
    else
	tmNew.iIndex = PARM_NOT_FOUND;
    tmNew.tpTorsion = *tpPTorsion;

		/* Find the proper place to put the term */

    for ( tmPCur = PVAI(tTorsion, TORSION_MATCHt, 0),i=0; 
    		i<iVarArrayElementCount(tTorsion); 
		i++, tmPCur++ ) {
	    if ( tmNew.tpTorsion.iN < tmPCur->tpTorsion.iN ) 
		break;
    }
    VarArrayInsertBefore( tTorsion, i, (GENP)&tmNew );
    if ( !strcmp( tmNew.tpTorsion.sType1, WILD_CARD_TYPE ) )
	    return(PARM_FOUND_WILD);
    return(PARM_FOUND_EXACT);
}


/*
 *	zParmSetAddToImproper
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	If (bUseIndex) is TRUE then write the index of the term
 *	within the PARMSET into the TORSION.iIndex, otherwise write
 *	PARM_NOT_FOUND.
 *
 */
static int
zbParmSetAddToImproper( TORSION tTorsion, int iIndex, TORSIONPARMt *tpPTorsion,
	BOOL bUseIndex )
{
TORSION_MATCHt	*tmPCur;
TORSION_MATCHt	tmNew;
int		i;

    memset( &tmNew, 0, sizeof(tmNew) );		/* for Purify */
    if ( !iVarArrayElementCount( tTorsion ) ) {

	/*
	 *  just add
	 */
        if ( bUseIndex ) 
	    	tmNew.iIndex = iIndex;
	else
		tmNew.iIndex = PARM_NOT_FOUND;
	tmNew.tpTorsion = *tpPTorsion;
	VarArrayAdd( tTorsion, (GENP)&tmNew );
	if ( !strcmp( tmNew.tpTorsion.sType1, WILD_CARD_TYPE ) )
		return(PARM_FOUND_WILD);
	return(PARM_FOUND_EXACT);
    }

    /*
     *  if the previous torsion is general and this one is
     *	specific, delete the previous one(s) and add the
     *	new one
     */
    tmPCur = PVAI(tTorsion, TORSION_MATCHt, 0);
    if ( strcmp( tmPCur->tpTorsion.sType1, WILD_CARD_TYPE ) == 0   &&
         strcmp( tpPTorsion->sType1, WILD_CARD_TYPE ) != 0 ) {
	VarArraySetSize( tTorsion, 0 );
        if ( bUseIndex ) 
	    	tmNew.iIndex = iIndex;
	else
		tmNew.iIndex = PARM_NOT_FOUND;
	tmNew.tpTorsion = *tpPTorsion;
	VarArrayAdd( tTorsion, (GENP)&tmNew );
	return(PARM_FOUND_EXACT);
    } 
		/* First check if it is already in the (Improper) TORSION */
		/* VARARRAY */

    for ( tmPCur = PVAI(tTorsion, TORSION_MATCHt, 0),i=0; 
    		i<iVarArrayElementCount(tTorsion); 
		i++, tmPCur++ ) {
	    if ( tmPCur->tpTorsion.iN == tpPTorsion->iN ) {

		    /* If the improper that we are adding is more specific */
		    /* than the one that is there then replace the more */
		    /* general one */

		    int	iNew, iOld;

	            iOld = iParmSetTorsionGenerality( tmPCur->tpTorsion.sType1,
	    					 tmPCur->tpTorsion.sType2,
						 tmPCur->tpTorsion.sType3,
						 tmPCur->tpTorsion.sType4 );
	            iNew = iParmSetTorsionGenerality( tpPTorsion->sType1,
	    					 tpPTorsion->sType2,
						 tpPTorsion->sType3,
						 tpPTorsion->sType4 );
	            if ( iNew < iOld ) {  
	                if ( bUseIndex ) 
			    tmPCur->iIndex = iIndex;
		        else
			    tmPCur->iIndex = PARM_NOT_FOUND;
		        tmPCur->tpTorsion = *tpPTorsion;
		        return(PARM_FOUND_EXACT);
		    }
		    return(PARM_NOT_FOUND);
	    }
    }

		/* If no match was found then simply add the */
		/* proper to the end of the TORSION */
    if ( bUseIndex ) 
    	tmNew.iIndex = iIndex;
    else
	tmNew.iIndex = PARM_NOT_FOUND;
    tmNew.tpTorsion = *tpPTorsion;

		/* Find the proper place to put the term */

    for ( tmPCur = PVAI(tTorsion, TORSION_MATCHt, 0),i=0; 
    		i<iVarArrayElementCount(tTorsion); 
		i++, tmPCur++ ) {
	    if ( tmNew.tpTorsion.iN < tmPCur->tpTorsion.iN ) 
		break;
    }
    VarArrayInsertBefore( tTorsion, i, (GENP)&tmNew );
    if ( !strcmp( tmNew.tpTorsion.sType1, WILD_CARD_TYPE ) )
	    return(PARM_FOUND_WILD);
    return(PARM_FOUND_EXACT);
}



/*
 *	zbParmSetBuildTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build the TORSION by searching through the
 *	PARMSET for all proper torsions that match the
 *	atom types.
 *
 *	The atom types must be in canonical order.
 *
 *	Return PARM_* status according to terms added.
 */
static int
zbParmSetBuildTorsion( PARMSET psParmSet, char *s1, char *s2, 
	char *s3, char *s4, TORSION tTorsion, BOOL bUseIndex )
{
int		i;
TORSIONPARMt	*tpPCur;
int		iMax;
int		iRet = PARM_NOT_FOUND;
    
    if ( (iMax = iVarArrayElementCount( psParmSet->vaTorsions )) ) {
    	tpPCur = PVAI( psParmSet->vaTorsions, TORSIONPARMt, 0 );
    	for ( i=0; i<iMax; i++, tpPCur++ ) {
		if ( zbParmSetMatchTorsion( tpPCur, s1, s2, s3, s4 ) ) {
			int kk = 
	    		    zbParmSetAddToTorsion( tTorsion, i, tpPCur, 
						bUseIndex );
			if ( iRet == PARM_NOT_FOUND  ||
			     ( iRet == PARM_FOUND_WILD && 
				  kk ==  PARM_FOUND_EXACT ) )
				iRet = kk;
		}
	}
    } 

    return(iRet);
}



 
/*
 *      zbParmSetMatchImproper
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the four types match the torsion.
 *      If on entry bExact is true then only exact parameters will
 *      match, otherwise wildcards are OK.
 *      The CENTRAL atom MUST be the THIRD ATOM!!!!!!!!!!!!!!!!!!!
 *
 *      On entry, iWild equals the number of wildcards allowed.
 */
static BOOL
zbParmSetMatchImproper( TORSIONPARMt *tpPTorsion, char *s1, char *s2, 
	char *s3, char *s4 )
{
int             iWild;
BOOL            bFoundOne;

    bFoundOne = FALSE;

                /* Check the central atom */
    if ( strcmp( tpPTorsion->sType3, s3 ) != 0 ) return(FALSE);

    iWild = iParmSetTorsionGenerality( 
    			tpPTorsion->sType1,
			tpPTorsion->sType2,
			tpPTorsion->sType3,
			tpPTorsion->sType4 );

    switch ( iWild ) {
        case 0:
	    if ( strcmp( tpPTorsion->sType1, s1 ) == 0 &&
	    	 strcmp( tpPTorsion->sType2, s2 ) == 0 &&
		 strcmp( tpPTorsion->sType4, s4 ) == 0 ) bFoundOne = TRUE;
	    break;
	case 1:
	    if ( strcmp( tpPTorsion->sType2, s1 ) == 0 ) {
	        if ( strcmp( tpPTorsion->sType4, s2 ) == 0 ||
		     strcmp( tpPTorsion->sType4, s4 ) == 0 ) bFoundOne = TRUE;
	    } else if ( strcmp( tpPTorsion->sType2, s2 ) == 0 &&
	    	 strcmp( tpPTorsion->sType4, s4 ) == 0 ) bFoundOne = TRUE;
	    break;
	case 2:
	    if ( strcmp( tpPTorsion->sType4, s1 ) == 0 ||
	         strcmp( tpPTorsion->sType4, s2 ) == 0 ||
		 strcmp( tpPTorsion->sType4, s4 ) == 0 ) bFoundOne = TRUE;
	    break;
	case 3:
	    bFoundOne = TRUE;
	    break;
	default:
	    DFATAL(( "Illegal number of wildcards (0-3) got: %d\n",
	    		iWild ));
	    break;
    }
    if ( bFoundOne ) {
        MESSAGE(( "Matched improper: %s-%s-%s-%s\n",
                tpPTorsion->sType1,
                tpPTorsion->sType2,
                tpPTorsion->sType3,
                tpPTorsion->sType4 ));
    }

    return(bFoundOne);
}





/*
 *	zParmSetBuildImproperTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Build the tTorsion by searching through the
 *	PARMSET for all improper torsions that match the
 *	atom types.
 *
 *	If (bUseIndex) is TRUE then write the index of the term within
 *	the PARMSET into the tTorsion.iIndex field, otherwise write
 *	PARM_NOT_FOUND.
 *
 *	The atom types must be in canonical order.
 *	Return TRUE if the tTorsion was actually changed.
 */
static BOOL
zbParmSetBuildImproperTorsion( PARMSET psParmSet, char *s1, char *s2, 
	char *s3, char *s4, TORSION tTorsion, BOOL bUseIndex )
{
int		i;
TORSIONPARMt	*tpPCur;
BOOL		bAddedOne;
int		iMax;

    bAddedOne = FALSE;
    if ( (iMax = iVarArrayElementCount( psParmSet->vaImpropers )) ) {
        tpPCur = PVAI( psParmSet->vaImpropers, TORSIONPARMt, 0 );
        for ( i=0; i<iMax; i++, tpPCur++ ) {
	    if ( zbParmSetMatchImproper( tpPCur, s1, s2, s3, s4 ) ) {
	        if ( zbParmSetAddToImproper( tTorsion, i, tpPCur, bUseIndex )
			!= PARM_NOT_FOUND )
		    bAddedOne = TRUE;
	    }
	}
    }
    return(bAddedOne);
}




/*
 *----------------------------------------------------------------
 *
 *	Public functions
 *
 */




/*
 *      psParmSetCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Allocate memory for a PARMSET and initialize
 *      its contents.
 */
PARMSET
psParmSetCreate()
{
PARMSET psNew;

    MALLOC( psNew, PARMSET, sizeof(PARMSETt) );

    strcpy( psNew->sFname, "no name assigned" );

    psNew->vaAtoms       = vaVarArrayCreate( sizeof(ATOMPARMt) );
    psNew->vaBonds       = vaVarArrayCreate( sizeof(BONDPARMt) );
    psNew->vaAngles      = vaVarArrayCreate( sizeof(ANGLEPARMt) );
    psNew->vaTorsions    = vaVarArrayCreate( sizeof(TORSIONPARMt) );
    psNew->vaImpropers   = vaVarArrayCreate( sizeof(TORSIONPARMt) );
    psNew->vaHBonds      = vaVarArrayCreate( sizeof(HBONDPARMt) );
    psNew->vaNBEdits     = vaVarArrayCreate( sizeof(NBEDITt) );

    psNew->bBeingEdited  = FALSE; /* V. Romanovski */

    return(psNew);
}

/*
 *      psParmSetDuplicate
 *
 *	Author:	Bill Ross (1993)
 *
 *      Duplicate the contents of the PARMSET. 
 *	Should only be called by oObjectDuplicate(),
 *	which sets objekt attributes.
 */
PARMSET 
psParmSetDuplicate( PARMSET psOld )
{
PARMSET	psNew;

    MALLOC( psNew, PARMSET, sizeof(PARMSETt) );
    memcpy( psNew, psOld, sizeof(PARMSETt) );
    psNew->bBeingEdited = FALSE;

    psNew->vaAtoms = vaVarArrayCopy( psOld->vaAtoms );
    psNew->vaBonds = vaVarArrayCopy( psOld->vaBonds );
    psNew->vaAngles = vaVarArrayCopy( psOld->vaAngles );
    psNew->vaTorsions = vaVarArrayCopy( psOld->vaTorsions );
    psNew->vaImpropers = vaVarArrayCopy( psOld->vaImpropers );
    psNew->vaHBonds = vaVarArrayCopy( psOld->vaHBonds );
    psNew->vaNBEdits = vaVarArrayCopy( psOld->vaNBEdits );

    return(psNew);
}



/*
 *      ParmSetDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the contents of the PARMSET.
 */
void
ParmSetDestroy( PARMSET *psPLib )
{
    VarArrayDestroy( &((*psPLib)->vaAtoms) );
    VarArrayDestroy( &((*psPLib)->vaBonds) );
    VarArrayDestroy( &((*psPLib)->vaAngles) );
    VarArrayDestroy( &((*psPLib)->vaTorsions) );
    VarArrayDestroy( &((*psPLib)->vaImpropers) );
    VarArrayDestroy( &((*psPLib)->vaHBonds) );
    VarArrayDestroy( &((*psPLib)->vaNBEdits) );
    FREE( *psPLib );
    *psPLib = NULL;
}



/* Ignore this stuff if LINT is being used - seems to make lint coredump */

#ifndef LINT

/*
 *      ParmSetSave
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Save the PARMSET into a DATABASE
 */
void
ParmSetSave( PARMSET psLib, DATABASE db )
{
                /* If the PARMSET is NULL then delete all of the */
                /* parameter tables from the database, if they exist at all */

    if ( psLib == NULL ) {
	VP0(( "PARMSET is NULL => deleting all parameter tables from dbase\n"));
        bDBRndDeleteEntry( db, "parm.atoms" );
        bDBRndDeleteEntry( db, "parm.bonds" );
        bDBRndDeleteEntry( db, "parm.angles" );
	/*
	 *  (torsions & impropers are folded together for
	 *	backwards dbase compatibility)
	 */
        bDBRndDeleteEntry( db, "parm.torsions" );
        bDBRndDeleteEntry( db, "parm.hbonds" );
        return;
    }
    
    /*
     *  atoms
     */
    if ( iVarArrayElementCount( psLib->vaAtoms ) ) {
        DBPutTable( db, "parm.atoms", iVarArrayElementCount(psLib->vaAtoms),
    
                5, "element",
		    (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->iElement),
		    iVarArrayElementSize(psLib->vaAtoms),
                6, "hybrid",
		    (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->iHybridization),
		    iVarArrayElementSize(psLib->vaAtoms),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                2, "mass",
		    (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->dMass),
                    iVarArrayElementSize(psLib->vaAtoms),
		8, "polar",
		    (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->dPolar),
		    iVarArrayElementSize(psLib->vaAtoms),
                3, "e",
		   (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->dEpsilon),
                   iVarArrayElementSize(psLib->vaAtoms),
                4, "r",
		    (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->dR),
                    iVarArrayElementSize(psLib->vaAtoms),
                1, "type",
		    (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->sType),
                    iVarArrayElementSize(psLib->vaAtoms),
                7, "desc",
		    (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->sDesc),
                    iVarArrayElementSize(psLib->vaAtoms), 
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0 );
    }

    /*
     *  bonds
     */
    if ( iVarArrayElementCount( psLib->vaBonds ) ) {
        DBPutTable( db, "parm.bonds", iVarArrayElementCount(psLib->vaBonds),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                3, "kb", 
		    (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->dKb),
                    iVarArrayElementSize(psLib->vaBonds),
                4, "r0", 
		    (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->dR0),
                    iVarArrayElementSize(psLib->vaBonds),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                1, "type1",
		    (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->sType1),
                    iVarArrayElementSize(psLib->vaBonds),
                2, "type2",
		    (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->sType2),
                    iVarArrayElementSize(psLib->vaBonds),
                5, "desc",
		    (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->sDesc),
                    iVarArrayElementSize(psLib->vaBonds),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0 );
    }

    /*
     *  angles
     */
    if ( iVarArrayElementCount( psLib->vaAngles ) ) {
        DBPutTable( db, "parm.angles", iVarArrayElementCount(psLib->vaAngles),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                4, "kt", 
		    (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->dKt),
                    iVarArrayElementSize(psLib->vaAngles),
                5, "t0", 
		    (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->dT0),
                    iVarArrayElementSize(psLib->vaAngles),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                1, "type1",
		    (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->sType1),
                    iVarArrayElementSize(psLib->vaAngles),
                2, "type2",
		    (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->sType2),
                    iVarArrayElementSize(psLib->vaAngles),
                3, "type3",
		    (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->sType3),
                    iVarArrayElementSize(psLib->vaAngles),
                6, "desc",
		    (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->sDesc),
                    iVarArrayElementSize(psLib->vaAngles),
                0, NULL, NULL, 0 );
    }

    /*
     *  torsions & impropers - folded together here for
     *	backward dbase compatibility
     */
    if ( iVarArrayElementCount( psLib->vaTorsions )  ||
         iVarArrayElementCount( psLib->vaImpropers ) ) {
	VARARRAY	vaTorsTypes; 
	TORSIONPARMt	*tP;
	int		i;

	/*
	 *  copy into 1 vararray & set iType
	 */
	vaTorsTypes = vaVarArrayCopy2( psLib->vaTorsions, psLib->vaImpropers );
	tP = PVAI(vaTorsTypes,TORSIONPARMt,0);
	for (i=0; i<iVarArrayElementCount(psLib->vaTorsions); i++, tP++)
		tP->iType = PROPER;
	for (i=0; i<iVarArrayElementCount(psLib->vaImpropers); i++, tP++)
		tP->iType = IMPROPER;

        DBPutTable( db, "parm.torsions", 
		iVarArrayElementCount(vaTorsTypes),
                5, "type",
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->iType),
                    iVarArrayElementSize(vaTorsTypes),
                7, "n",
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->iN),
		    iVarArrayElementSize(vaTorsTypes),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                6, "kp", 
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->dKp),
                    iVarArrayElementSize(vaTorsTypes),
                8, "p0", 
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->dP0),
                    iVarArrayElementSize(vaTorsTypes),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                1, "type1",
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sType1),
                    iVarArrayElementSize(vaTorsTypes),
                2, "type2",
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sType2),
                    iVarArrayElementSize(vaTorsTypes),
                3, "type3",
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sType3),
                    iVarArrayElementSize(vaTorsTypes),
                4, "type4",
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sType4),
                    iVarArrayElementSize(vaTorsTypes),
                9, "desc",
		    (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sDesc),
                    iVarArrayElementSize(vaTorsTypes)
                );
        DBPutValue( db, "parm.torsionOrders", ENTRYSTRING|ENTRYARRAY, 
			iVarArrayElementCount(vaTorsTypes),
			PVAI(vaTorsTypes,TORSIONPARMt,0)->sOrder,
			iVarArrayElementSize(vaTorsTypes) );
	VarArrayDestroy( &vaTorsTypes );
    }

    /*
     *  hbonds
     */
    if ( iVarArrayElementCount( psLib->vaHBonds ) ) {
        DBPutTable( db, "parm.hbonds", 
			iVarArrayElementCount(psLib->vaHBonds),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                3, "a", 
		    (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->dA),
                    iVarArrayElementSize(psLib->vaHBonds),
                4, "b", 
		    (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->dB),
                    iVarArrayElementSize(psLib->vaHBonds),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0,

                1, "type1",
		    (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->sType1),
                    iVarArrayElementSize(psLib->vaHBonds),
                2, "type2",
		    (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->sType2),
                    iVarArrayElementSize(psLib->vaHBonds),
                5, "desc",
		    (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->sDesc),
                    iVarArrayElementSize(psLib->vaHBonds),
                0, NULL, NULL, 0,
                0, NULL, NULL, 0 );
    }
}






/*
 *      psParmSetLoad
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Load the PARMSET from a DATABASE.
 *      If there is no parameter set then return NULL.
 */
PARMSET
psParmSetLoad( DATABASE db )
{
PARMSET		psLib;
int             iType, iLines;
int		i;

                /* Load the atom type and non-bond parameters */

    if ( !bDBGetType( db, "parm.atoms", &iType, &iLines ) &&
         !bDBGetType( db, "parm.bonds", &iType, &iLines ) &&
         !bDBGetType( db, "parm.angles", &iType, &iLines ) &&
         !bDBGetType( db, "parm.torsions", &iType, &iLines ) &&
         !bDBGetType( db, "parm.hbonds", &iType, &iLines ) )
	return(NULL);
    
    psLib = (PARMSET)oCreate(PARMSETid);

    /*
     *  atoms
     */
    VarArraySetSize( (psLib->vaAtoms), iLines );
    if ( iLines ) 
    	bDBGetTable( db, "parm.atoms", &iLines,
                5, (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->iElement),
		    iVarArrayElementSize(psLib->vaAtoms),
                6, (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->iHybridization),
		    iVarArrayElementSize(psLib->vaAtoms),
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,

                2, (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->dMass),
                   iVarArrayElementSize(psLib->vaAtoms),
		8, (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->dPolar),
		   iVarArrayElementSize(psLib->vaAtoms),
                3, (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->dEpsilon),
                   iVarArrayElementSize(psLib->vaAtoms),
                4, (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->dR),
                   iVarArrayElementSize(psLib->vaAtoms),

                1, (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->sType),
                   iVarArrayElementSize(psLib->vaAtoms),
                7, (char *)&(PVAI(psLib->vaAtoms,ATOMPARMt,0)->sDesc),
                   iVarArrayElementSize(psLib->vaAtoms),
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0 );


    /*
     *  bonds
     */
    bDBGetType( db, "parm.bonds", &iType, &iLines );
    VarArraySetSize( (psLib->vaBonds), iLines );
    if ( iLines ) 
        bDBGetTable( db, "parm.bonds", &iLines,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,

                3, (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->dKb),
                   iVarArrayElementSize(psLib->vaBonds),
                4, (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->dR0),
                   iVarArrayElementSize(psLib->vaBonds),
                0, NULL, 0,
                0, NULL, 0,

                1, (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->sType1),
                   iVarArrayElementSize(psLib->vaBonds),
                2, (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->sType2),
                   iVarArrayElementSize(psLib->vaBonds),
                5, (char *)&(PVAI(psLib->vaBonds,BONDPARMt,0)->sDesc),
                   iVarArrayElementSize(psLib->vaBonds),
                0, NULL, 0,
                0, NULL, 0 );


    /*
     *  angles
     */
    bDBGetType( db, "parm.angles", &iType, &iLines );
    VarArraySetSize( (psLib->vaAngles), iLines );
    if ( iLines ) 
    	bDBGetTable( db, "parm.angles", &iLines,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,

                4, (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->dKt),
                   iVarArrayElementSize(psLib->vaAngles),
                5, (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->dT0),
                   iVarArrayElementSize(psLib->vaAngles),
                0, NULL, 0,
                0, NULL, 0,

                1, (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->sType1),
                   iVarArrayElementSize(psLib->vaAngles),
                2, (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->sType2),
                   iVarArrayElementSize(psLib->vaAngles),
                3, (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->sType3),
                   iVarArrayElementSize(psLib->vaAngles),
                6, (char *)&(PVAI(psLib->vaAngles,ANGLEPARMt,0)->sDesc),
                   iVarArrayElementSize(psLib->vaAngles),
                0, NULL, 0 );


    /*
     *  torsions & impropers are merged for backward dbase
     *	compatibility - disentangle here
     */
    bDBGetType( db, "parm.torsions", &iType, &iLines );
    MESSAGE(( "There are %d torsion+improper parameters.\n" ));
    if ( iLines ) {
	VARARRAY	vaTorsTypes; 
	TORSIONPARMt	*tP, *tP2;
	int		iCount = 0;

	vaTorsTypes = vaVarArrayCreate( sizeof(TORSIONPARMt) );
    	VarArraySetSize( vaTorsTypes, iLines );
    	bDBGetTable( db, "parm.torsions", &iLines,
                5, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->iType),
                 iVarArrayElementSize(vaTorsTypes),
                7, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->iN),
		   iVarArrayElementSize(vaTorsTypes),
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,

                6, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->dKp),
                   iVarArrayElementSize(vaTorsTypes),
                8, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->dP0),
                   iVarArrayElementSize(vaTorsTypes),
                0, NULL, 0,
                0, NULL, 0,

                1, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sType1),
                   iVarArrayElementSize(vaTorsTypes),
                2, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sType2),
                   iVarArrayElementSize(vaTorsTypes),
                3, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sType3),
                   iVarArrayElementSize(vaTorsTypes),
                4, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sType4),
                   iVarArrayElementSize(vaTorsTypes),
                9, (char *)&(PVAI(vaTorsTypes,TORSIONPARMt,0)->sDesc),
                   iVarArrayElementSize(vaTorsTypes)
                );
        if ( bDBGetType( db, "parm.torsionOrders", &iType, &iLines ) ) {
            bDBGetValue( db, "parm.torsionOrders", &iLines,
			PVAI(vaTorsTypes,TORSIONPARMt,0)->sOrder,
			iVarArrayElementSize(vaTorsTypes) );
	} else {
	    for ( i=0; i<iVarArrayElementCount(vaTorsTypes); i++ ) {
		strcpy( PVAI(vaTorsTypes,TORSIONPARMt,i)->sOrder,"0123");
	    }
	}

	/*
	 *  count propers, set size, & copy
	 */
	tP = PVAI(vaTorsTypes,TORSIONPARMt,0);
	for (i=0; i<iVarArrayElementCount(vaTorsTypes); i++, tP++)
		if ( tP->iType == PROPER )
			iCount++;
    	VarArraySetSize( psLib->vaTorsions, iCount );
	tP2 = PVAI(psLib->vaTorsions,TORSIONPARMt,0);
	tP = PVAI(vaTorsTypes,TORSIONPARMt,0);
	for (i=0; i<iVarArrayElementCount(vaTorsTypes); i++, tP++) {
		if ( tP->iType == PROPER ) {
			*tP2 = *tP;
			tP2++;
		}
	}

	/*
	 *  set improper size, copy
	 */
	iCount = iVarArrayElementCount(vaTorsTypes) - iCount;
    	VarArraySetSize( psLib->vaImpropers, iCount );
	tP2 = PVAI(psLib->vaImpropers,TORSIONPARMt,0);
	tP = PVAI(vaTorsTypes,TORSIONPARMt,0);
	for (i=0; i<iVarArrayElementCount(vaTorsTypes); i++, tP++) {
		if ( tP->iType == IMPROPER ) {
			*tP2 = *tP;
			tP2++;
		}
	}

	VarArrayDestroy( &vaTorsTypes );
    }

    /*
     *  hbonds
     */
    bDBGetType( db, "parm.hbonds", &iType, &iLines );
    VarArraySetSize( (psLib->vaHBonds), iLines );
    if ( iLines ) 
    	bDBGetTable( db, "parm.hbonds", &iLines,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,

                3, (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->dA),
                   iVarArrayElementSize(psLib->vaHBonds),
                4, (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->dB),
                   iVarArrayElementSize(psLib->vaHBonds),
                0, NULL, 0,
                0, NULL, 0,

                1, (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->sType1),
                   iVarArrayElementSize(psLib->vaHBonds),
                2, (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->sType2),
                   iVarArrayElementSize(psLib->vaHBonds),
                5, (char *)&(PVAI(psLib->vaHBonds,HBONDPARMt,0)->sDesc),
                   iVarArrayElementSize(psLib->vaHBonds),
                0, NULL, 0,
                0, NULL, 0 );

      return(psLib);
}


#endif          /* ifndef LINT */




/*
 *      ParmSetDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe the PARMSET
 */
void
ParmSetDescribe( PARMSET psLib )
{
ATOMPARMt      *apPAtom;
BONDPARMt      *bpPBond;
ANGLEPARMt     *apPAngle;
TORSIONPARMt   *tpPTorsion;
HBONDPARMt     *hpPHBond;
NBEDITt        *hpPNBEdit;
int             i, iMax;
STRING		sElement;

    BasicsResetInterrupt();
    VP0(( "PARMSET\n" ));

                /* Dump Atoms */

    VP0(( "--Atoms\n" ));
    iMax = iVarArrayElementCount( psLib->vaAtoms );

    if ( !iMax ) {
	VP0(( "  --None\n" ));
    } else {
	apPAtom = PVAI( psLib->vaAtoms, ATOMPARMt, 0 );
	for ( i=0; i<iMax; apPAtom++, i++ ) {
	    STRING	s1;
	    if ( apPAtom->dPolar == -1 )
		sprintf(s1, "def=0");
	    else
		sprintf(s1, "%8.2lf", apPAtom->dPolar);
	    VP0(( "    %4s  Mass=%8.2lf  Polar=%s  E = %8.2lf  R=%8.2lf\n",
			apPAtom->sType, apPAtom->dMass, s1,
			apPAtom->dEpsilon, apPAtom->dR ));
	    VP0(( "           Element=%s  Hybrid= Sp%d  Desc:%s\n",
			sElementName( apPAtom->iElement, sElement ), 
			apPAtom->iHybridization, apPAtom->sDesc ));
	    if ( bBasicsInterrupt() ) 
		goto QUIT;
	}
    }

                /* Dump bonds */

    VP0(( "--Bonds\n" ));
    iMax = iVarArrayElementCount( psLib->vaBonds );
    if ( !iMax ) {
	VP0(( "  --None\n" ));
    } else {
	bpPBond = PVAI( psLib->vaBonds, BONDPARMt, 0 );
	for ( i=0; i<iMax; bpPBond++, i++ ) {
	    VP0(( "    %4s - %4s   Kb=%8.2lf   R0=%8.2lf   Desc:%s\n",
			bpPBond->sType1, bpPBond->sType2,
			bpPBond->dKb, bpPBond->dR0, bpPBond->sDesc ));
	    if ( bBasicsInterrupt() ) 
		goto QUIT;
	}
    }

                /* Dump angles */

    VP0(( "--Angles\n" ));
    iMax = iVarArrayElementCount( psLib->vaAngles );
    if ( !iMax ) {
	VP0(( "  --None\n" ));
    } else {
	apPAngle = PVAI( psLib->vaAngles, ANGLEPARMt, 0 );
	for ( i=0; i<iMax; apPAngle++, i++ ) {
	    VP0(( "    %4s - %4s - %4s   Kt=%8.2lf   T0=%8.2lf   Desc:%s\n",
			apPAngle->sType1, apPAngle->sType2, apPAngle->sType3,
			apPAngle->dKt, apPAngle->dT0/DEGTORAD, 
			apPAngle->sDesc ));
	    if ( bBasicsInterrupt() ) 
		goto QUIT;
	}
    }
     
                /* Dump torsions */

    VP0(( "--Torsions\n" ));
    iMax = iVarArrayElementCount( psLib->vaTorsions );
    if ( !iMax ) {
	VP0(( "  --None\n" ));
    } else {
	tpPTorsion = PVAI( psLib->vaTorsions, TORSIONPARMt, 0 );
	for ( i=0; i<iMax; tpPTorsion++, i++ ) {
	    VP0(( "  %4s - %4s - %4s - %4s\n",
			tpPTorsion->sType1, tpPTorsion->sType2, 
			tpPTorsion->sType3, tpPTorsion->sType4 ));
	    VP0(( "        Kp=%8.2lf   N=%d   P0=%8.2lf   Order: %s  Desc:%s\n",
			tpPTorsion->dKp, tpPTorsion->iN,
			tpPTorsion->dP0/DEGTORAD, tpPTorsion->sOrder, 
			tpPTorsion->sDesc ));
	    if ( bBasicsInterrupt() ) 
		goto QUIT;
	}
    }

                /* Dump impropers */

    VP0(( "--Impropers\n" ));
    iMax = iVarArrayElementCount( psLib->vaImpropers );
    if ( !iMax ) {
	VP0(( "  --None\n" ));
    } else {
	tpPTorsion = PVAI( psLib->vaImpropers, TORSIONPARMt, 0 );
	for ( i=0; i<iMax; tpPTorsion++, i++ ) {
	    VP0(( "  %4s - %4s - %4s - %4s\n",
			tpPTorsion->sType1, tpPTorsion->sType2, 
			tpPTorsion->sType3, tpPTorsion->sType4 ));
	    VP0(( "        Kp=%8.2lf   N=%d   P0=%8.2lf   Order: %s  Desc:%s\n",
			tpPTorsion->dKp, tpPTorsion->iN,
			tpPTorsion->dP0/DEGTORAD, tpPTorsion->sOrder, 
			tpPTorsion->sDesc ));
	    if ( bBasicsInterrupt() ) 
		goto QUIT;
	}
    }

                /* Dump Hbonds */

    VP0(( "--HBonds\n" ));
    iMax = iVarArrayElementCount( psLib->vaHBonds );
    if ( !iMax ) {
	VP0(( "  --None\n" ));
    } else {
	hpPHBond = PVAI( psLib->vaHBonds, HBONDPARMt, 0 );
	for ( i=0; i<iMax; hpPHBond++, i++ ) {
	    VP0(( "    %4s - %4s   A=%8.2lf   B=%8.2lf   Desc:%s\n",
			hpPHBond->sType1, hpPHBond->sType2,
			hpPHBond->dA, hpPHBond->dB, hpPHBond->sDesc ));
	    if ( bBasicsInterrupt() ) 
		goto QUIT;
	}
    }

                /* Dump NB Edits */

    VP0(( "--Nonbonded Edits\n" ));
    iMax = iVarArrayElementCount( psLib->vaNBEdits );
    if ( !iMax ) {
      VP0(( "  --None\n" ));
    } else {
      hpPNBEdit = PVAI( psLib->vaNBEdits, NBEDITt, 0 );
      for ( i=0; i<iMax; hpPNBEdit++, i++ ) {
	VP0(( "    %4s - %4s   EI=%8.5lf  EJ=%8.5lf  RI=%8.5lf  RJ=%8.5lf   "
	      "Desc:%s\n", hpPNBEdit->sType1, hpPNBEdit->sType2,
	      hpPNBEdit->dEI, hpPNBEdit->dEJ, hpPNBEdit->dRI, hpPNBEdit->dRJ,
	      hpPNBEdit->sDesc ));
	if ( bBasicsInterrupt() )
	  goto QUIT;
      }
    }

    VP0(( "\n" ));
    return;

QUIT:
    VP0(( "Interrupted\n" ));
    BasicsResetInterrupt();
}






/*
 *      iParmSetAddAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add an atom parameter to the PARMSET.
 *      Return the index.
 */
int
iParmSetAddAtom( PARMSET psLib, char *sType, double dMass, double dPolar, 
	double dEpsilon, double dR, double dEpsilon14, double dR14,
	double dScreenF, int iElement, int iHybridization, char *sDesc )
{
ATOMPARMt	apAtom;

    memset ( &apAtom, 0, sizeof(apAtom) );	/* for Purify */
    strcpy( apAtom.sType, sType );
    apAtom.dMass 	= dMass;
    apAtom.dPolar	= dPolar;
    apAtom.dEpsilon 	= dEpsilon;
    apAtom.dR 		= dR;
    apAtom.dEpsilon14 	= dEpsilon14;
    apAtom.dR14 	= dR14;
    apAtom.dScreenF     = dScreenF;
    apAtom.iElement	= iElement;
    apAtom.iHybridization= iHybridization;
    if ( sDesc != NULL )
    	strcpy( apAtom.sDesc, sDesc );
    else
    	strcpy( apAtom.sDesc, "" );

    VarArrayAdd( (psLib->vaAtoms), (GENP)&apAtom );

    return(iVarArrayElementCount( psLib->vaAtoms )-1);
}





/*
 *      iParmSetAddBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a bond parameter to the PARMSET.
 *      Return the index.
 */
int
iParmSetAddBond( PARMSET psLib, char *sType1, char *sType2, 
	double dKb, double dR0, char *sDesc )
{
BONDPARMt       bpBond;

    memset( &bpBond, 0, sizeof(bpBond) );	/* for Purify */
    strcpy( bpBond.sType1, sType1 );
    strcpy( bpBond.sType2, sType2 );
    zParmSetOrderBondAtoms( bpBond.sType1, bpBond.sType2 );
    bpBond.dKb = dKb;
    bpBond.dR0 = dR0;
    if ( sDesc != NULL )
    	strcpy( bpBond.sDesc, sDesc);
    else
	strcpy( bpBond.sDesc, "" );

    VarArrayAdd( (psLib->vaBonds), (GENP)&bpBond ); 

    return(iVarArrayElementCount( psLib->vaBonds )-1);
}


/*
 *      iParmSetAddAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add an angle parameter to the PARMSET.
 *      Return the index.
 */
int
iParmSetAddAngle( PARMSET psLib, char *sType1, char *sType2, char *sType3, 
	double dKt, double dT0, double dTkub, double dRkub, char *sDesc )
{
ANGLEPARMt      apAngle;

    memset( &apAngle, 0, sizeof(apAngle) );	/* for Purify */
    strcpy( apAngle.sType1, sType1 );
    strcpy( apAngle.sType2, sType2 );
    strcpy( apAngle.sType3, sType3 );
    zParmSetOrderAngleAtoms( apAngle.sType1, apAngle.sType2, apAngle.sType3 );
    apAngle.dKt = dKt;
    apAngle.dT0 = dT0;
    apAngle.dTkub = dTkub;
    apAngle.dRkub = dRkub;
    if ( sDesc != NULL )
	strcpy( apAngle.sDesc, sDesc);
    else
	strcpy( apAngle.sDesc, "" );

    VarArrayAdd( (psLib->vaAngles), (GENP)&apAngle ); 

    return(iVarArrayElementCount( psLib->vaAngles )-1);
}



/*
 *      iParmSetAddProperTerm
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a torsion parameter to the PARMSET.
 *      Return the index. 
 */
int
iParmSetAddProperTerm( PARMSET psLib, 
	char *sType1, char *sType2, char *sType3, char *sType4, 
	int iN, double dKp, double dP0, double dScEE, double dScNB,
	char *sDesc )
{
TORSIONPARMt    tpTorsion;

    memset( &tpTorsion, 0, sizeof(tpTorsion) );		/* for Purify */
    strcpy( tpTorsion.sType1, sType1 );
    strcpy( tpTorsion.sType2, sType2 );
    strcpy( tpTorsion.sType3, sType3 );
    strcpy( tpTorsion.sType4, sType4 );
    zParmSetOrderTorsionAtoms( tpTorsion.sType1,
			  tpTorsion.sType2,
			  tpTorsion.sType3,
			  tpTorsion.sType4 );
    tpTorsion.dKp = dKp;
    tpTorsion.iN  = iN;
    tpTorsion.dP0 = dP0;
    tpTorsion.dScEE = dScEE;
    tpTorsion.dScNB = dScNB;
    strcpy( tpTorsion.sOrder, "0123" );
    if ( sDesc != NULL )
    	strcpy( tpTorsion.sDesc, sDesc );
    else
	strcpy( tpTorsion.sDesc, "" );

    VarArrayAdd( psLib->vaTorsions, (GENP)&tpTorsion ); 

    return(iVarArrayElementCount( psLib->vaTorsions )-1);
}






/*
 *      iParmSetAddImproperTerm
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add an improper torsion parameter to the PARMSET.
 *      Return the index.
 *      The THIRD atom MUST be the central atom!!!!!!!!!!
 */
int
iParmSetAddImproperTerm( PARMSET psLib, 
	char *sType1, char *sType2, char *sType3, char *sType4, 
	int iN, double dKp, double dP0, double dScEE, double dScNB, char *sDesc )
{
TORSIONPARMt    tpImproper;
orderStr	sOrder;

    memset( &tpImproper, 0, sizeof(tpImproper) );	/* for Purify */
    strcpy( sOrder, "0123" );
    strcpy( tpImproper.sType1, sType1 );
    strcpy( tpImproper.sType2, sType2 );
    strcpy( tpImproper.sType3, sType3 );
    strcpy( tpImproper.sType4, sType4 );

    zParmSetOrderImproperAtoms( tpImproper.sType1,
			   tpImproper.sType2,
			   tpImproper.sType3,
			   tpImproper.sType4, sOrder );
    tpImproper.dKp = dKp;
    tpImproper.iN  = iN;
    tpImproper.dP0 = dP0;
    tpImproper.dScEE = 0.;
    tpImproper.dScNB = 0.;

    if ( sDesc != NULL )
    	strcpy( tpImproper.sDesc, sDesc );
    else
	strcpy( tpImproper.sDesc, "" );

    strcpy( tpImproper.sOrder, sOrder );

    VarArrayAdd( psLib->vaImpropers, (GENP)&tpImproper ); 

    return(iVarArrayElementCount( psLib->vaImpropers )-1);
}






/*
 *      iParmSetAddHBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a bond parameter to the PARMSET.
 *      Return the index.
 */
int
iParmSetAddHBond( PARMSET psLib, char *sType1, char *sType2, 
	double dA, double dB, char *sDesc )
{
HBONDPARMt      hpHBond;

    memset( &hpHBond, 0, sizeof(hpHBond) );	/* for Purify */
    strcpy( hpHBond.sType1, sType1 );
    strcpy( hpHBond.sType2, sType2 );
    zParmSetOrderBondAtoms( hpHBond.sType1, hpHBond.sType2 );
    hpHBond.dA = dA;
    hpHBond.dB = dB;
    if ( sDesc != NULL )
    	strcpy( hpHBond.sDesc, sDesc );
    else
	strcpy( hpHBond.sDesc, "" );

    VarArrayAdd( (psLib->vaHBonds), (GENP)&hpHBond ); 

    return(iVarArrayElementCount( psLib->vaHBonds )-1);
}


/*
 *      iParmSetAddNBEdit
 *
 *      Author: David S. Cerutti (2013)
 *
 *      Adjust pairwise Lennard-Jones mixed parameters within the PARMSET.
 *      Return the index.
 */
int
iParmSetAddNBEdit( PARMSET psLib, char *sType1, char *sType2, double dEI,
		   double dEJ, double dRI, double dRJ, char *sDesc )
{
  int i;
  double       dA, dC;
  NBEDITt      hpNBEdit;

  memset( &hpNBEdit, 0, sizeof(hpNBEdit) );     /* for Purify */
  strcpy( hpNBEdit.sType1, sType1 );
  strcpy( hpNBEdit.sType2, sType2 );
  zParmSetOrderBondAtoms( hpNBEdit.sType1, hpNBEdit.sType2 );
  hpNBEdit.dEI = dEI;
  hpNBEdit.dEJ = dEJ;
  hpNBEdit.dRI = dRI;
  hpNBEdit.dRJ = dRJ;

  /* Compute the pair sigma and epsilon that WOULD */
  /* exist under the standard combining rules.     */
  dEI = 0.0;
  dEJ = 0.0;
  dRI = 0.0;
  dRJ = 0.0;
  for (i = 0; i < iVarArrayElementCount( psLib->vaAtoms ); i++) {
    if (strcmp(PVAI(psLib->vaAtoms, ATOMPARMt, i)->sType, sType1) == 0) {
      dEI = PVAI(psLib->vaAtoms, ATOMPARMt, i)->dEpsilon;
      dRI = PVAI(psLib->vaAtoms, ATOMPARMt, i)->dR;
    }
    if (strcmp(PVAI(psLib->vaAtoms, ATOMPARMt, i)->sType, sType2) == 0) {
      dEJ = PVAI(psLib->vaAtoms, ATOMPARMt, i)->dEpsilon;
      dRJ = PVAI(psLib->vaAtoms, ATOMPARMt, i)->dR;
    }
  }
  MathOpConvertNonBondToAC(dEI, dRI, dEJ, dRJ, &dA, &dC);
  hpNBEdit.dA = dA;
  hpNBEdit.dC = dC;
  if ( sDesc != NULL )
    strcpy( hpNBEdit.sDesc, sDesc );
  else
    strcpy( hpNBEdit.sDesc, "" );

  VarArrayAdd( (psLib->vaNBEdits), (GENP)&hpNBEdit );

  return(iVarArrayElementCount( psLib->vaNBEdits )-1);
}


// iParmSetAddCMAP
/*
int iParmSetCMAP(PARMSET psLib){
    return 0;
}
 */








/*
 *      iParmSetFindAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Search for an atom parameter in the ParmSet
 *      and return the index if it is found
 *      otherwise return PARM_NOT_FOUND.
 */
int
iParmSetFindAtom( PARMSET psLib, char *sType )
{
ATOMPARMt	*apPAtom;
int             i, iMax;
BOOL            bFoundOne;

    iMax = iVarArrayElementCount( psLib->vaAtoms );

    if ( !iMax )
	return(PARM_NOT_FOUND);

    bFoundOne = FALSE;
    apPAtom = PVAI( psLib->vaAtoms, ATOMPARMt, 0 );
    for ( i=0; i<iMax; apPAtom++, i++ ) {
        if ( strcmp( apPAtom->sType, sType ) == 0 ) {
                    bFoundOne = TRUE;
                    break;
        }
    }
    if ( bFoundOne )
	return(i);
    return(PARM_NOT_FOUND);
}






/*
 *      iParmSetFindBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Search for a bond parameter in the ParmSet
 *      and return the index if it is found
 *      otherwise return PARM_NOT_FOUND
 */
int
iParmSetFindBond( PARMSET psLib, char *sType1, char *sType2 )
{
BONDPARMt	*bpPBond;
int             i, iMax;
BOOL            bFoundOne;
STRING		s1, s2;

    iMax = iVarArrayElementCount( psLib->vaBonds );
    if ( !iMax )
	return(PARM_NOT_FOUND);

    strcpy( s1, sType1 );
    strcpy( s2, sType2 );
    zParmSetOrderBondAtoms( s1, s2 );

    bFoundOne = FALSE;
    bpPBond = PVAI( psLib->vaBonds, BONDPARMt, 0 );
    for ( i=0; i<iMax; bpPBond++, i++ ) {
        if ( strcmp( bpPBond->sType1, s1 ) == 0 ) {
            if ( strcmp( bpPBond->sType2, s2 ) == 0 ) {
                    bFoundOne = TRUE;
                    break;
            }
        }
    }
    if ( bFoundOne ) {
        MESSAGE(( "-Bond Parameter %s - %s\n", 
		sType1, sType2 ));
	return(i);
    }
    return(PARM_NOT_FOUND);
}





/*
 *      iParmSetFindAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Search for a angle parameter in the ParmSet
 *      and return the index if it is found
 *      otherwise return PARM_NOT_FOUND.
 */
int
iParmSetFindAngle( PARMSET psLib, char *sType1, char *sType2, char *sType3 )
{
ANGLEPARMt	*apPAngle;
int             i, iMax;
BOOL            bFoundOne;
STRING		s1, s2, s3;

    iMax = iVarArrayElementCount( psLib->vaAngles );
    if ( !iMax )
	return(PARM_NOT_FOUND);

    strcpy( s1, sType1 );
    strcpy( s2, sType2 );
    strcpy( s3, sType3 );
    zParmSetOrderAngleAtoms( s1, s2, s3 );

    bFoundOne = FALSE;
    apPAngle = PVAI( psLib->vaAngles, ANGLEPARMt, 0 );
    for ( i=0; i<iMax; apPAngle++, i++ ) {
        if ( strcmp( apPAngle->sType1, s1 ) == 0 ) {
            if ( (strcmp( apPAngle->sType2, s2 ) == 0) &&
                 (strcmp( apPAngle->sType3, s3 ) == 0) ) {
                    bFoundOne = TRUE;
                    break;
            }
        }
    }
    if ( bFoundOne ) {
	MESSAGE(( "-Angle Parameter %s - %s - %s\n", 
                sType1, sType2, sType3));
    }
    if ( bFoundOne )
	return(i);
    return(PARM_NOT_FOUND);
}





/*
 *      iParmSetFindProperTerms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Search for all torsion terms in the PARMSET
 *	that match the ATOM types passed.
 *	The terms are added to the TORSION.
 *	The terms are only added to the TORSION if there
 *	is not a term for that multiplicity or, the term in the
 *	TORSION is less specific than the term found within the PARMSET.
 *
 *	If bUseIndex is TRUE then the index of the term within
 *	the PARMSET will be written into the TORSION.iIndex, otherwise
 *	PARM_NOT_FOUND will be written.
 *
 *	If terms are found and added then return PARM_FOUND_TERMS,
 *	otherwise return PARM_NOT_FOUND.
 *
 *	The caller is responsible for making sure that the TORSION is 
 *	valid.
 */
int
iParmSetFindProperTerms( PARMSET psLib, TORSION tTorsion, BOOL bUseIndex,
		char *sType1, char *sType2, char *sType3, char *sType4 )
{
STRING		s1, s2, s3, s4;

                /* First look for specific parameters */

    strcpy( s1, sType1 );
    strcpy( s2, sType2 );
    strcpy( s3, sType3 );
    strcpy( s4, sType4 );

    zParmSetOrderTorsionAtoms( s1, s2, s3, s4 );

    return( zbParmSetBuildTorsion( psLib, 
    				s1, s2, s3, s4, 
    				tTorsion, bUseIndex ) );
}




/*
 *	iParmSetFindImproperTerms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Search for all improper torsion terms in the PARMSET
 *	that match the ATOM types passed.
 *	The terms are added to the caller's tTorsion.
 *	The terms are only added to the tTorsion if there
 *	is not already a term for that multiplicity or, the term in the
 *	tTorsion is less specific than the term found within the PARMSET.
 *
 *	If terms are found and added then return PARM_FOUND_TERMS,
 *	otherwise return PARM_NOT_FOUND.
 *
 *	If bUseIndex is TRUE then write the index of the term within
 *	the PARMSET into the tTorsion.iIndex field, otherwise write
 *	PARM_NOT_FOUND.
 *
 *	The caller is responsible for making sure that the tTorsion is 
 *	valid.
 */
int
iParmSetFindImproperTerms( PARMSET psLib, TORSION tTorsion, BOOL bUseIndex,
			char *sType1, char *sType2, char *sType3, char *sType4 )
{
STRING		s1, s2, s3, s4;
orderStr	sOrder;

                /* First look for specific parameters */

    strcpy( s1, sType1 );
    strcpy( s2, sType2 );
    strcpy( s3, sType3 );
    strcpy( s4, sType4 );
    strcpy( sOrder, "0123" );

    zParmSetOrderImproperAtoms( s1, s2, s3, s4, sOrder );
    
    if ( zbParmSetBuildImproperTorsion( psLib, 
    					   s1, s2, s3, s4, tTorsion,
					   bUseIndex ) )
	return(PARM_FOUND_TERMS);
    return(PARM_NOT_FOUND);
}





/*
 *      iParmSetFindHBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Search for an HBond parameter in the ParmSet
 *      and return the index if it is found
 *      Otherwise return PARM_NOT_FOUND.
 */
int
iParmSetFindHBond( PARMSET psLib, char *sType1, char *sType2 )
{
HBONDPARMt	*hpPHBond;
int             i, iMax;
BOOL            bFoundOne;
STRING		s1, s2;

    iMax = iVarArrayElementCount( psLib->vaHBonds );
    if ( !iMax )
	return(PARM_NOT_FOUND);

    strcpy( s1, sType1 );
    strcpy( s2, sType2 );
    zParmSetOrderBondAtoms( s1, s2 );

    bFoundOne = FALSE;
    hpPHBond = PVAI( psLib->vaHBonds, HBONDPARMt, 0 );
    for ( i=0; i<iMax; hpPHBond++, i++ ) {
        if ( strcmp( hpPHBond->sType1, s1 ) == 0 ) {
            if ( strcmp( hpPHBond->sType2, s2 ) == 0 ) {
                    bFoundOne = TRUE;
                    break;
            }
        }
    }
    if ( bFoundOne ) {
	MESSAGE(( "-HBond Parameter %s - %s\n", 
		sType1, sType2 ));
	return(i);
    } 
    return(PARM_NOT_FOUND);
}





/*
 *      iParmSetFindNBEdit
 *
 *	Author:	David S. Cerutti (2013)
 *
 *      Search for a non-bonded adjustment in the ParmSet
 *      and return the index if it is found
 *      Otherwise return PARM_NOT_FOUND.
 */
int
iParmSetFindNBEdit( PARMSET psLib, char *sType1, char *sType2 )
{
NBEDITt         *hpPNBEdit;
int             i, iMax;
BOOL            bFoundOne;
STRING		s1, s2;

    iMax = iVarArrayElementCount( psLib->vaNBEdits );
    if ( !iMax )
	return(PARM_NOT_FOUND);

    strcpy( s1, sType1 );
    strcpy( s2, sType2 );
    zParmSetOrderBondAtoms( s1, s2 );

    bFoundOne = FALSE;
    hpPNBEdit = PVAI( psLib->vaNBEdits, NBEDITt, 0 );
    for ( i=0; i<iMax; hpPNBEdit++, i++ ) {
        if ( strcmp( hpPNBEdit->sType1, s1 ) == 0 ) {
            if ( strcmp( hpPNBEdit->sType2, s2 ) == 0 ) {
                    bFoundOne = TRUE;
                    break;
            }
        }
    }
    if ( bFoundOne ) {
	MESSAGE(( "-HBond Parameter %s - %s\n", sType1, sType2 ));
	return(i);
    }
    return(PARM_NOT_FOUND);
}





/*
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *
 *	TORSION routines.
 *
 */

/*
 *    tParmSetTORSIONCreate
 *
 *    Author: Christian Schafmeister (1991)
 *
 *    Create an empty TORSION.
 *	This is a routine (& not a #define) to hide TORSION_MATCHt.
 */
TORSION
tParmSetTORSIONCreate()
{
    return(vaVarArrayCreate(sizeof(TORSION_MATCHt)));
}


/*
 *	ParmSetTORSIONTerm
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the parameters associated with the TORSION element.
 */
void
ParmSetTORSIONTerm( TORSION tTorsion, int iTorsionIndex, int *iPParmSetIndex,
		   char *cPType1, char *cPType2, char *cPType3, char *cPType4,
		   int *iPN, double *dPKp, double *dPP0, double *dPScEE,
		   double *dPScNB, char *sDesc )
{
TORSION_MATCHt	*tmPCur;

    tmPCur = PVAI( tTorsion, TORSION_MATCHt, iTorsionIndex );
    *iPParmSetIndex = tmPCur->iIndex;
    strcpy( cPType1, tmPCur->tpTorsion.sType1 );
    strcpy( cPType2, tmPCur->tpTorsion.sType2 );
    strcpy( cPType3, tmPCur->tpTorsion.sType3 );
    strcpy( cPType4, tmPCur->tpTorsion.sType4 );
    *iPN = tmPCur->tpTorsion.iN;
    *dPKp = tmPCur->tpTorsion.dKp;
    *dPP0 = tmPCur->tpTorsion.dP0;
    *dPScEE = tmPCur->tpTorsion.dScEE;
    *dPScNB = tmPCur->tpTorsion.dScNB;
    strcpy(sDesc, tmPCur->tpTorsion.sDesc);
}



/*
 *	ParmSetTORSIONAddProperTerm
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add the term to the TORSION.
 */
BOOL
bParmSetTORSIONAddProperTerm( TORSION tTorsion,
	char *cPType1, char *cPType2, char *cPType3, char *cPType4,
	int iN, double dKp, double dP0, double dScEE, double dScNB,
	char *sDesc )
{
TORSIONPARMt	tpTorsion;

    strcpy( tpTorsion.sType1, cPType1 );
    strcpy( tpTorsion.sType2, cPType2 );
    strcpy( tpTorsion.sType3, cPType3 );
    strcpy( tpTorsion.sType4, cPType4 );
    zParmSetOrderTorsionAtoms( tpTorsion.sType1,
    			    tpTorsion.sType2,
			    tpTorsion.sType3,
			    tpTorsion.sType4 );

    tpTorsion.iN = iN;
    tpTorsion.dKp = dKp;
    tpTorsion.dP0 = dP0;
    tpTorsion.dScEE = dScEE;
    tpTorsion.dScNB = dScNB;
    strcpy(tpTorsion.sDesc, sDesc);
    strcpy( tpTorsion.sOrder, "0123" );
    
    if (zbParmSetAddToTorsion( tTorsion, 0, &tpTorsion, FALSE )
		!= PARM_NOT_FOUND )
        return(TRUE);
    return(FALSE);
}



/*
 *	ParmSetTORSIONAddImproperTerm
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add the term to the improper TORSION.
 */
BOOL
bParmSetTORSIONAddImproperTerm( TORSION tTorsion,
	char *cPType1, char *cPType2, char *cPType3, char *cPType4,
	int iN, double dKp, double dP0, double dScEE, double dScNB, char *sDesc )
{
TORSIONPARMt	tpTorsion;
orderStr	sOrder;

    strcpy( tpTorsion.sType1, cPType1 );
    strcpy( tpTorsion.sType2, cPType2 );
    strcpy( tpTorsion.sType3, cPType3 );
    strcpy( tpTorsion.sType4, cPType4 );
    strcpy( sOrder, "0123" );
    zParmSetOrderImproperAtoms( tpTorsion.sType1,
    			    tpTorsion.sType2,
			    tpTorsion.sType3,
			    tpTorsion.sType4, sOrder );
    tpTorsion.iN = iN;
    tpTorsion.dKp = dKp;
    tpTorsion.dP0 = dP0;
    tpTorsion.dScEE = 0.;
    tpTorsion.dScNB = 0.;
    strcpy(tpTorsion.sDesc, sDesc);
    strcpy( tpTorsion.sOrder, sOrder );
    
    if (zbParmSetAddToTorsion( tTorsion, 0, &tpTorsion, FALSE )
		!= PARM_NOT_FOUND )
        return(TRUE);
    return(FALSE);
}





/*
 *	ParmSetImproperOrderAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Order the ATOM _types_ and their indices 
 *
 *	Since AMBER impropers require an explicit ordering
 *	of atoms around the central atom, this routine must
 *	reorder the improper found by 'loop.c' according to
 *	the order in the torsion parameter field 'sOrder'.
 *	Also, the list is secondarily sorted by atom number.
 *	The central atom of an improper is always in the
 *	3rd place (slot 2 of 0..3). The other three atoms are 
 *	referred to as 'peripheral atoms'.
 *
 */
void
ParmSetImproperOrderAtoms( TORSION tTorsion, int iTorsionIndex,
	char *cPaTypes[4], int iaIndexes[4] )
{
TORSION_MATCHt	*tmPCur;
char		*cPaTempTypes[4];
int		iaTempIndexes[4];
int		i, j, iPut;
char		*cPTemp;
char		*cPParm[4];
int		iTemp, iBetter, iBetterIndex;

    tmPCur = PVAI( tTorsion, TORSION_MATCHt, iTorsionIndex );

    cPParm[0] = tmPCur->tpTorsion.sType1;
    cPParm[1] = tmPCur->tpTorsion.sType2;
    cPParm[2] = tmPCur->tpTorsion.sType3;
    cPParm[3] = tmPCur->tpTorsion.sType4;

    /*
     *  get exact match for slot 0 if necc
     */
    if ( strcmp( cPParm[0], WILD_CARD_TYPE) != 0  ) {

	/* 
	 *  There's no wild card, so must have exact match of type for 0th slot
	 */
	if (strcmp( cPParm[0], cPaTypes[0] ) != 0) {

	    /*
	     *  must dig it up
	     */
	    if ( strcmp( cPParm[0], cPaTypes[1] ) == 0  &&
	         strcmp( cPParm[1], cPaTypes[1] ) != 0 ) {
		/*
		 *  atom 1 matches 0th slot and mismatches or wild in own
		 */
		SWAP( cPaTypes[1], cPaTypes[0], cPTemp );
		SWAP( iaIndexes[1], iaIndexes[0], iTemp );
	    } else if ( strcmp( cPParm[0], cPaTypes[3] ) == 0  &&
		        strcmp( cPParm[3], cPaTypes[3] ) != 0 ) {
		/*
		 *  atom 3 matches 0th slot and mismatches or wild in own
		 */
		SWAP( cPaTypes[3], cPaTypes[0], cPTemp );
		SWAP( iaIndexes[3], iaIndexes[0], iTemp );
	    } else
		DFATAL(( "Could not order atoms, could not find 0th type: %s\n",
				cPParm[0] ));
	}
    }

    /*
     *  get exact match for slot 1 if necc
     */
    if ( strcmp( cPParm[1], WILD_CARD_TYPE) != 0  ) {

	/* 
	 *  no wild card, so must find exact match for 1st slot
	 */

	if (strcmp( cPParm[1], cPaTypes[1] ) != 0) {

	    /*
	     *  must dig it up
	     */
	    if ( strcmp( cPParm[1], cPaTypes[0] ) == 0  &&
	         strcmp( cPParm[0], WILD_CARD_TYPE) == 0 ) {
		/*
		 *  atom 0 matches 1st slot and wild in own
		 */
		SWAP( cPaTypes[1], cPaTypes[0], cPTemp );
		SWAP( iaIndexes[1], iaIndexes[0], iTemp );
	    } else if ( strcmp( cPParm[1], cPaTypes[3] ) == 0  &&
		        strcmp( cPParm[3], cPaTypes[3] ) != 0 ) {
		/*
		 *  atom 3 matches 1st slot and mismatches or wild in own
		 */
		SWAP( cPaTypes[3], cPaTypes[1], cPTemp );
		SWAP( iaIndexes[3], iaIndexes[1], iTemp );
	    } else
		DFATAL(( "Could not order atoms, could not find 1st type: %s\n",
				cPParm[1] ));
	}
    }

    /*
     *  get exact match for (last) slot 3 if necc
     */
    if ( strcmp( cPParm[3], WILD_CARD_TYPE ) != 0  ) {

	/* 
	 *  no wild card, so must find exact match for 3rd slot
	 */
	if (strcmp( cPParm[3], cPaTypes[3] ) != 0) {

	    /*
	     *  must dig it up
	     */
	    if ( strcmp( cPParm[3], cPaTypes[0] ) == 0  &&
	         strcmp( cPParm[0], WILD_CARD_TYPE ) == 0 ) {
		/*
		 *  atom 0 matches 3rd slot and wild in own
		 */
		SWAP( cPaTypes[3], cPaTypes[0], cPTemp );
		SWAP( iaIndexes[3], iaIndexes[0], iTemp );
	    } else if ( strcmp( cPParm[3], cPaTypes[1] ) == 0  &&
		        strcmp( cPParm[1], WILD_CARD_TYPE ) == 0 ) {
		/*
		 *  atom 1 matches 3rd slot and wild in own
		 */
		SWAP( cPaTypes[3], cPaTypes[1], cPTemp );
		SWAP( iaIndexes[3], iaIndexes[1], iTemp );
	    } else
		DFATAL(( "Could not order atoms, could not find 3rd type: %s\n",
				cPParm[3] ));
	}
    }

    /*
     *  update the order to match the result of the sortings
     */
		/* Now fix the order of the ATOMs to what it */
		/* was in the original parameter set */
    iPut = tmPCur->tpTorsion.sOrder[0]-'0';
    cPaTempTypes[iPut] = cPaTypes[0];
    iaTempIndexes[iPut] = iaIndexes[0];

    iPut = tmPCur->tpTorsion.sOrder[1]-'0';
    cPaTempTypes[iPut] = cPaTypes[1];
    iaTempIndexes[iPut] = iaIndexes[1];

    iPut = tmPCur->tpTorsion.sOrder[2]-'0';
    cPaTempTypes[iPut] = cPaTypes[2];
    iaTempIndexes[iPut] = iaIndexes[2];

    iPut = tmPCur->tpTorsion.sOrder[3]-'0';
    cPaTempTypes[iPut] = cPaTypes[3];
    iaTempIndexes[iPut] = iaIndexes[3];

    /*
     *  having finished with parameter set-based ordering,
     *	sort the peripheral wild card atoms by topological order 
     */
    for (i=0; i<4; i++) {
	if ( i == 2 )	/* fixed atom */
	    continue;
	if (strcmp( cPParm[i], WILD_CARD_TYPE ) != 0 )   /* not wild card */
	    continue;

	/*
	 *  find the least topological order
	 *	wild card atom ahead
	 */
	iBetter = -1;
	iBetterIndex = iaTempIndexes[i];
	for (j=i+1; j<4; j++) {
		if ( j == 2 )	/* fixed atom */
	    	    continue;
		if ( strcmp(cPParm[j], WILD_CARD_TYPE ) != 0 )
		    continue;
		if ( iaTempIndexes[j] < iBetterIndex ) {
			iBetter = j;
			iBetterIndex = iaTempIndexes[j];
		} 
	}
	if ( iBetter != -1 ) {
		SWAP( cPaTempTypes[i], cPaTempTypes[iBetter], cPTemp );
		SWAP( iaTempIndexes[i], iaTempIndexes[iBetter], iTemp );
	}
    }
    
    /*
     *  sort same-type peripheral atoms by topological order
     */
    for (i=0; i<4; i++) {
	if ( i == 2 )	/* center atom */
	    continue;
	/*
	 *  find least-index of same type ahead
	 */
	iBetter = -1;
	iBetterIndex = iaTempIndexes[i];
	for (j=i+1; j<4; j++) {
	    if ( j == 2 )	/* center atom */
	        continue;
	    if ( strcmp( cPaTempTypes[i], cPaTempTypes[j] ) != 0 )
		continue;
	    /*
	     *  atoms are the same type, so check order
	     */
	    if ( iaTempIndexes[j] < iBetterIndex ) {
		iBetter = j;
		iBetterIndex = iaTempIndexes[j];
	    }
	}
	if ( iBetter != -1 ) {
		SWAP( cPaTempTypes[i], cPaTempTypes[iBetter], cPTemp );
		SWAP( iaTempIndexes[i], iaTempIndexes[iBetter], iTemp );
	}
    }
    /*
     *  finally, set the indexes used by the caller
     */
    for ( i=0; i<4; i++ ) {
	cPaTypes[i] = cPaTempTypes[i];
	iaIndexes[i] = iaTempIndexes[i];
    }
}




/*
 *      bParmSetCapableOfHBonding
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the atom type is capable of being involved
 *      in a HBond.
 */
BOOL
bParmSetCapableOfHBonding( PARMSET psParms, char *sType )
{
HBONDPARMt     *hbPCur;
int		iCount, iTotal;

	/* If there are no HBONDS then nothing can HBOND */

    if ( iVarArrayElementCount(psParms->vaHBonds) == 0 ) {
	return(FALSE);
    }
    iTotal = iVarArrayElementCount(psParms->vaHBonds);

    hbPCur = PVAI( psParms->vaHBonds, HBONDPARMt, 0 );
    for ( iCount = 0; iCount < iTotal; iCount++) {
        if ( strcmp( hbPCur->sType1, sType ) == 0 ) return(TRUE);
        if ( strcmp( hbPCur->sType2, sType ) == 0 ) return(TRUE);
        hbPCur++;
    }
    return(FALSE);
}






/*
 *      ParmSetAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return an atom parameter in the ParmSet
 */
void
ParmSetAtom( PARMSET psLib, int i, char *sType, double *dPMass, 
	double *dPPolar, double *dPEpsilon, double *dPR, double *dPEpsilon14,
    double *dPR14, double *dPScreenF, int *iPElement, int *iPHybridization,
    char *sDesc )
{
ATOMPARMt      *apPAtom;

    if ( !iVarArrayElementCount( psLib->vaAtoms ) ) {

	/*
	 *  default values
	 */
	VP0(( "WARNING - using default atom values (NOELEMENT)\n" ));
	strcpy( sType, WILD_CARD_TYPE );
	*dPMass = 0.0;
	*dPPolar = -1.0;
	*dPEpsilon = 0.0;
	*dPR = 0.0;
	*dPEpsilon14 = 0.0;
	*dPR14 = 0.0;
	*dPScreenF = 0.0;
	*iPElement = NOELEMENT;
	*iPHybridization = 0;
	strcpy( sDesc, "??" );
	return;
    }
    apPAtom = PVAI( psLib->vaAtoms, ATOMPARMt, i );
    strcpy( sType, apPAtom->sType);
    *dPMass = apPAtom->dMass;
    *dPPolar = apPAtom->dPolar;
    *dPEpsilon = apPAtom->dEpsilon;
    *dPR = apPAtom->dR;
    *dPEpsilon14 = apPAtom->dEpsilon14;
    *dPR14 = apPAtom->dR14;
    *dPScreenF = apPAtom->dScreenF;
    *iPElement = apPAtom->iElement;
    *iPHybridization = apPAtom->iHybridization;
    strcpy( sDesc, apPAtom->sDesc );
}




/*
 *      ParmSetBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a bond parameter in the ParmSet
 */
void
ParmSetBond( PARMSET psLib, int i, char *sType1, char *sType2, 
	double *dPKb, double *dPR0, char *sDesc )
{
BONDPARMt      *bpPBond;

    if ( !iVarArrayElementCount( psLib->vaBonds ) ) {

	/*
	 *  default values
	 */
	VP0(( "WARNING - using default bond values (0)\n" ));
	strcpy( sType1, WILD_CARD_TYPE );
	strcpy( sType2, WILD_CARD_TYPE );
	*dPKb = 0.0;
	*dPR0 = 0.0;
	strcpy( sDesc, "??" );
	return;
    }
    bpPBond = PVAI( psLib->vaBonds, BONDPARMt, i );
    strcpy( sType1, bpPBond->sType1 );
    strcpy( sType2, bpPBond->sType2 );
    *dPKb = bpPBond->dKb;
    *dPR0 = bpPBond->dR0;
    strcpy( sDesc, bpPBond->sDesc );
}





/*
 *      ParmSetAngle
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a angle parameter in the ParmSet
 */
void
ParmSetAngle( PARMSET psLib, int i, char *sType1, char *sType2, char *sType3,
	double *dPKt, double *dPT0, double *dPTkub, double *dPRkub, char *sDesc )
{
ANGLEPARMt     *apPAngle;

    if ( !iVarArrayElementCount( psLib->vaAngles ) ) {
	/*
	 *  default values
	 */
	VP0(( "WARNING - using default angle values (0)\n" ));
	strcpy( sType1, WILD_CARD_TYPE );
	strcpy( sType2, WILD_CARD_TYPE );
	strcpy( sType3, WILD_CARD_TYPE );
	*dPKt = 0.0;
	*dPT0 = 0.0;
	strcpy( sDesc, "??" );
	return;
    }
    apPAngle = PVAI( psLib->vaAngles, ANGLEPARMt, i );
    strcpy( sType1, apPAngle->sType1 );
    strcpy( sType2, apPAngle->sType2 );
    strcpy( sType3, apPAngle->sType3 );
    *dPKt = apPAngle->dKt;
    *dPT0 = apPAngle->dT0;
    *dPTkub = apPAngle->dTkub;
    *dPRkub = apPAngle->dRkub;
    strcpy( sDesc, apPAngle->sDesc );
}


/*
 *      ParmSetTorsion
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the i'th torsion parameter in the ParmSet
 */
void
ParmSetTorsion( PARMSET psLib, int i, 
	char *sType1, char *sType2, char *sType3, char *sType4,
	int *iPN, double *dPKp, double *dPP0, double *dPScEE,
	double *dPScNB, char *sDesc)
{
TORSIONPARMt   *tpPTorsion;

    if ( !iVarArrayElementCount( psLib->vaTorsions ) ) {
	/*
	 *  default values
	 */
	VP0(( "WARNING - using default torsion values (0)\n" ));
	strcpy( sType1, WILD_CARD_TYPE );
	strcpy( sType2, WILD_CARD_TYPE );
	strcpy( sType3, WILD_CARD_TYPE );
	strcpy( sType4, WILD_CARD_TYPE );
	*iPN  = 0;
	*dPKp = 0;
	*dPP0 = 0;
	strcpy( sDesc, "??" );
	return;
    }
    tpPTorsion = PVAI( psLib->vaTorsions, TORSIONPARMt, i );
    strcpy( sType1, tpPTorsion->sType1 );
    strcpy( sType2, tpPTorsion->sType2 );
    strcpy( sType3, tpPTorsion->sType3 );
    strcpy( sType4, tpPTorsion->sType4 );
    *iPN  = tpPTorsion->iN;
    *dPKp = tpPTorsion->dKp;
    *dPP0 = tpPTorsion->dP0;
    *dPScEE = tpPTorsion->dScEE;
    *dPScNB = tpPTorsion->dScNB;
    strcpy( sDesc, tpPTorsion->sDesc );
}


/*
 *      ParmSetImproper
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return an improper parameter in the ParmSet
 */
void    
ParmSetImproper( PARMSET psLib, int i, 
	char *sType1, char *sType2, char *sType3, char *sType4,
	int *iPN, double *dPKp, double *dPP0, char *sDesc)
{
TORSIONPARMt   *tpPImproper;

    if ( !iVarArrayElementCount( psLib->vaImpropers ) ) {
	/*
	 *  default values
	 */
	VP0(( "WARNING - using default improper torsion values (0)\n" ));
	strcpy( sType1, WILD_CARD_TYPE );
	strcpy( sType2, WILD_CARD_TYPE );
	strcpy( sType3, WILD_CARD_TYPE );
	strcpy( sType4, WILD_CARD_TYPE );
	*iPN  = 0;
	*dPKp = 0;
	*dPP0 = 0;
	strcpy( sDesc, "??" );
	return;
    }
    tpPImproper = PVAI( psLib->vaImpropers, TORSIONPARMt, i );
    strcpy( sType1, tpPImproper->sType1 );
    strcpy( sType2, tpPImproper->sType2 );
    strcpy( sType3, tpPImproper->sType3 );
    strcpy( sType4, tpPImproper->sType4 );
    *iPN  = tpPImproper->iN;
    *dPKp = tpPImproper->dKp;
    *dPP0 = tpPImproper->dP0;
    strcpy( sDesc, tpPImproper->sDesc );
}



/*
 *      ParmSetHBond
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return an HBond parameter in the ParmSet
 */
void
ParmSetHBond( PARMSET psLib, int i, char *sType1, char *sType2, 
	double *dPA, double *dPB, char *sDesc )
{
HBONDPARMt     *hpPHBond;

    if ( !iVarArrayElementCount( psLib->vaHBonds ) ) {
	/*
	 *  default values
	 */
	VP0(( "WARNING - using default hbond values (0)\n" ));
	strcpy( sType1, WILD_CARD_TYPE );
	strcpy( sType2, WILD_CARD_TYPE );
	*dPA = 0;
	*dPB = 0;
	strcpy( sDesc, "??" );
	return;
    }
    hpPHBond = PVAI( psLib->vaHBonds, HBONDPARMt, i );
    strcpy( sType1, hpPHBond->sType1 );
    strcpy( sType2, hpPHBond->sType2 );
    *dPA = hpPHBond->dA;
    *dPB = hpPHBond->dB;
    strcpy( sDesc, hpPHBond->sDesc );
}



/*
 *      ParmSetNBEdit
 *
 *	Author:	David S. Cerutti (2013)
 *
 *      Return a non-bonded pair edit in the ParmSet
 */
void
ParmSetNBEdit( PARMSET psLib, int i, char *sType1, char *sType2, 
	       double *dEI, double *dEJ, double *dRI, double *dRJ,
	       char *sDesc )
{
NBEDITt     *hpPNBEdit;

    if ( !iVarArrayElementCount( psLib->vaNBEdits ) ) {
	/*
	 *  default values
	 */
	VP0(( "WARNING - using zeros for edited Lennard-Jones values (0)\n" ));
	strcpy( sType1, WILD_CARD_TYPE );
	strcpy( sType2, WILD_CARD_TYPE );
	*dEI = 0.0;
	*dEJ = 0.0;
	*dRI = 0.0;
	*dRJ = 0.0;
	strcpy( sDesc, "??" );
	return;
    }
    hpPNBEdit = PVAI( psLib->vaNBEdits, NBEDITt, i );
    strcpy( sType1, hpPNBEdit->sType1 );
    strcpy( sType2, hpPNBEdit->sType2 );
    *dEI = hpPNBEdit->dEI;
    *dEJ = hpPNBEdit->dEJ;
    *dRI = hpPNBEdit->dRI;
    *dRJ = hpPNBEdit->dRJ;
    strcpy( sDesc, hpPNBEdit->sDesc );
}


/*************************************

	ParmSetUpdatexxxx - These functions change the values stored in
	the parmsets parameter tables.  Only the values passed to the
	function are updated so that the parameters that the calling
	function does not want to update are not, if no value is passed.
	
	
	David Rivkin's Additions
	14 August 1992
	
**************************************/


/*
 *      ParmSetUpdateAtom
 *
 *	Author:	David Rivkin
 *
 *      Return an atom parameter in the ParmSet
 */
void
ParmSetUpdateAtom( PARMSET psLib, int i, char *sType, 
	double *dPMass, double *dPPolar, double *dPEpsilon, double *dPR, 
	double *dPScreenF,
	int *iPElement, int *iPHybrid, char *sDescription )
{
ATOMPARMt	*apPAtom;

    apPAtom = PVAI( psLib->vaAtoms, ATOMPARMt, i );

    if (       sType != (char*)NULL  )   strcpy( apPAtom->sType, sType);
    if (      dPMass != (double*)NULL)   apPAtom->dMass = *dPMass;
    if (     dPPolar != (double*)NULL)	apPAtom->dPolar = *dPPolar;
    if (   dPEpsilon != (double*)NULL)	apPAtom->dEpsilon = *dPEpsilon;
    if (         dPR != (double*)NULL)   apPAtom->dR = *dPR;
    if (   dPScreenF != (double*)NULL)  apPAtom->dScreenF = *dPScreenF;
    if (   iPElement != (int*)NULL   )	apPAtom->iElement = *iPElement;
    if (    iPHybrid != (int*)NULL   )	apPAtom->iHybridization = *iPHybrid;
    if (sDescription != (char*)NULL  )	strcpy(apPAtom->sDesc, sDescription);
}




/*
 *      ParmSetUpdateBond
 *
 *	Author:	David Rivkin (1992)
 *	Modified: Christian Schafmeister (Nov 1992)
 *			Atom types have to be pre-ordered.
 *
 *      Return a bond parameter in the ParmSet
 */
void
ParmSetUpdateBond( PARMSET psLib, int i, char *sType1, char *sType2, 
	double *dPKb, double *dPR0, char *sDescription)
{
BONDPARMt      *bpPBond;

    bpPBond = PVAI( psLib->vaBonds, BONDPARMt, i );
    if(      sType1 != (char*)NULL  ) strcpy(bpPBond->sType1, sType1 );
    if(      sType2 != (char*)NULL  ) strcpy(bpPBond->sType2, sType2 );
    if(        dPKb != (double*)NULL) bpPBond->dKb = *dPKb;
    if(        dPR0 != (double*)NULL) bpPBond->dR0 = *dPR0;
    if(sDescription != (char*)NULL  ) strcpy(bpPBond->sDesc, sDescription);
    if ( sType1 || sType2 )
    	zParmSetOrderBondAtoms( bpPBond->sType1, bpPBond->sType2 );
}





/*
 *      ParmSetUpdateAngle
 *
 *	Author:	David Rivkin (1992)
 *	Modified: Christian Schafmeister (Nov 1992)
 *			Atom types have to be pre-ordered.
 *
 *      Return a angle parameter in the ParmSet
 */
void
ParmSetUpdateAngle( PARMSET psLib, int i, 
	char *sType1, char *sType2, char *sType3, 
	double *dPKt, double *dPT0, char *sDescription )
{
ANGLEPARMt     *apPAngle;

    apPAngle = PVAI( psLib->vaAngles, ANGLEPARMt, i );
    if(      sType1 != (char*)NULL  ) strcpy( apPAngle->sType1, sType1 );
    if(      sType2 != (char*)NULL  ) strcpy( apPAngle->sType2, sType2 );
    if(      sType3 != (char*)NULL  ) strcpy( apPAngle->sType3, sType3 );
    if(        dPKt != (double*)NULL) apPAngle->dKt = *dPKt;
    if(        dPT0 != (double*)NULL) apPAngle->dT0 = *dPT0;
    if(sDescription != (char*)NULL  ) strcpy(apPAngle->sDesc, sDescription);

    if ( sType1  ||  sType2 )
    	zParmSetOrderAngleAtoms( apPAngle->sType1,
				 apPAngle->sType2,
				 apPAngle->sType3 );
}





/*
 *      ParmSetUpdateTorsion
 *
 *	Author:	David Rivkin (1992)
 *	Modified: Christian Schafmeister (Nov 1992)
 *			Atom types have to be pre-ordered.
 *
 *      Return a torsion parameter in the ParmSet
 */
void
ParmSetUpdateTorsion( PARMSET psLib, int i, 
	char *sType1, char *sType2, char *sType3, char *sType4,
	int *iPN, double *dPKp, double *dPP0, double *dPScEE,
	double *dPScNB, char *sDescription)
{
TORSIONPARMt	*tpPTorsion;
orderStr	sOrder;

    tpPTorsion = PVAI( psLib->vaTorsions, TORSIONPARMt, i );

    if (      sType1 != (char*)NULL  ) strcpy( tpPTorsion->sType1, sType1 );
    if (      sType2 != (char*)NULL  ) strcpy( tpPTorsion->sType2, sType2 );
    if (      sType3 != (char*)NULL  ) strcpy( tpPTorsion->sType3, sType3 );
    if (      sType4 != (char*)NULL  ) strcpy( tpPTorsion->sType4, sType4 );
    if (         iPN != (int*)NULL   ) tpPTorsion->iN = *iPN;
    if (        dPKp != (double*)NULL) tpPTorsion->dKp = *dPKp;
    if (        dPP0 != (double*)NULL) tpPTorsion->dP0 = *dPP0;
    if (      dPScEE != (double*)NULL) tpPTorsion->dScEE = *dPScEE;
    if (      dPScNB != (double*)NULL) tpPTorsion->dScNB = *dPScNB;
    if (sDescription != (char*)NULL  ) strcpy(tpPTorsion->sDesc, sDescription);

    strcpy( sOrder, "0123" );
    zParmSetOrderTorsionAtoms( tpPTorsion->sType1,
			  tpPTorsion->sType2,
			  tpPTorsion->sType3,
			  tpPTorsion->sType4 );
    strcpy( tpPTorsion->sOrder, sOrder );
}


/*
 *      ParmSetUpdateImproper
 *
 *	Author:	David Rivkin (1992)
 *	Modified: Christian Schafmeister (Nov 1992)
 *			Atom types have to be pre-ordered.
 *
 *      Return a torsion parameter in the ParmSet
 */
void
ParmSetUpdateImproper( PARMSET psLib, int i, 
	char *sType1, char *sType2, char *sType3, char *sType4,
	int *iPN, double *dPKp, double *dPP0, double *dScEE, double *dScNB, char *sDescription)
{
TORSIONPARMt   *tpPTorsion;
orderStr	sOrder;

    tpPTorsion = PVAI( psLib->vaImpropers, TORSIONPARMt, i );
    if (      sType1 != (char*)NULL  ) strcpy( tpPTorsion->sType1, sType1 );
    if (      sType2 != (char*)NULL  ) strcpy( tpPTorsion->sType2, sType2 );
    if (      sType3 != (char*)NULL  ) strcpy( tpPTorsion->sType3, sType3 );
    if (      sType4 != (char*)NULL  ) strcpy( tpPTorsion->sType4, sType4 );
    if (         iPN != (int*)NULL   ) tpPTorsion->iN = *iPN;
    if (        dPKp != (double*)NULL) tpPTorsion->dKp = *dPKp;
    if (        dPP0 != (double*)NULL) tpPTorsion->dP0 = *dPP0;
    if (sDescription != (char*)NULL  ) strcpy(tpPTorsion->sDesc, sDescription);

    strcpy( sOrder, "0123" );
    zParmSetOrderImproperAtoms( tpPTorsion->sType1,
			  tpPTorsion->sType2,
			  tpPTorsion->sType3,
			  tpPTorsion->sType4, sOrder );
    strcpy( tpPTorsion->sOrder, sOrder );
}



/*
 *      ParmSetUpdateHBond
 *
 *	Author:	David Rivkin (1992)
 *	Modified: Christian Schafmeister (Nov 1992)
 *			Atom types have to be pre-ordered.
 *
 *      Return an HBond parameter in the ParmSet
 */
void
ParmSetUpdateHBond( PARMSET psLib, int i, char *sType1, char *sType2, 
	double *dPA, double *dPB, char *sDescription )
{
HBONDPARMt     *hpPHBond;

    hpPHBond = PVAI( psLib->vaHBonds, HBONDPARMt, i );
    if(      sType1 != (char*)NULL  )	strcpy( hpPHBond->sType1, sType1 );
    if(      sType2 != (char*)NULL  )	strcpy( hpPHBond->sType2, sType2 );
    if(         dPA != (double*)NULL) 	hpPHBond->dA = *dPA;
    if(         dPB != (double*)NULL) 	hpPHBond->dB = *dPB;
    if(sDescription != (char*)NULL  )	strcpy(hpPHBond->sDesc, sDescription);

    if ( sType1 || sType2 )
    	zParmSetOrderBondAtoms( hpPHBond->sType1, hpPHBond->sType2 );
}


/*
 *	ParmSetNewAtoms
 *
 *	Author:	David Rivkin (1992)
 *
 *	Destroy the old parameters and create a new set that is empty
 *	of that can hold the iCount number of parameters
 *
 */
 
void
ParmSetNewAtoms( PARMSET psParmSet, int iCount )
{

    VarArrayDestroy( &psParmSet->vaAtoms );
    psParmSet->vaAtoms = vaVarArrayCreate( sizeof( ATOMPARMt ));
    VarArraySetSize( psParmSet->vaAtoms, iCount );
    MESSAGE(( "Atom parameters size changed to %i\n", iCount ));
}


/*
 *	ParmSetNewBonds
 *
 *	Author:	David Rivkin (1992)
 *
 *	Destroy the old parameters and create a new set that is empty
 *	of that can hold the iCount number of parameters
 *
 */
 
void 
ParmSetNewBonds( PARMSET psParmSet, int iCount )
{
    VarArrayDestroy( &psParmSet->vaBonds );
    psParmSet->vaBonds = vaVarArrayCreate( sizeof( BONDPARMt ));
    VarArraySetSize( psParmSet->vaBonds, iCount );
    MESSAGE(( "Bond parameters size changed to %i\n", iCount ));
}

/*
 *	ParmSetNewAngles
 *
 *	Author:	David Rivkin (1992)
 *
 *	Destroy the old parameters and create a new set that is empty
 *	of that can hold the iCount number of parameters
 *
 */
 
void
ParmSetNewAngles( PARMSET psParmSet, int iCount )
{
    VarArrayDestroy( &psParmSet->vaAngles );
    psParmSet->vaAngles = vaVarArrayCreate( sizeof( ANGLEPARMt ));
    VarArraySetSize( psParmSet->vaAngles, iCount );
    MESSAGE(( "Angle parameters size changed to %i\n", iCount ));
}
   	
    	
/*
 *	ParmSetNewTorsions
 *
 *	Author:	David Rivkin (1992)
 *
 *	Destroy the old parameters and create a new set that is empty
 *	of that can hold the iCount number of parameters
 *
 */
 
void
ParmSetNewTorsions( PARMSET psParmSet, int iCount )
{
    VarArrayDestroy( &psParmSet->vaTorsions );
    psParmSet->vaTorsions = vaVarArrayCreate( sizeof( TORSIONPARMt ));
    VarArraySetSize( psParmSet->vaTorsions, iCount );
    MESSAGE(( "Torsion parameters size changed to %i\n", iCount ));
}  



/*
 *	ParmSetNewHBonds
 *
 *	Author:	David Rivkin (1992)
 *
 *	Destroy the old parameters and create a new set that is empty
 *	of that can hold the iCount number of parameters
 *
 */
 
void
ParmSetNewHBonds( PARMSET psParmSet, int iCount )
{
    VarArrayDestroy( &psParmSet->vaHBonds );
    psParmSet->vaHBonds = vaVarArrayCreate( sizeof( HBONDPARMt ));
    VarArraySetSize( psParmSet->vaHBonds, iCount );
    MESSAGE(( "Hydgrogen Bond parameters size changed to %i\n", iCount ));
}



/*
 *	ParmSetNewImpropers
 *
 *	Author:	David Rivkin (1992)
 *
 *	Create a new Torsion parameter array that will hold the expansion/contraction
 *
 */
 
void 
ParmSetNewImpropers( PARMSET psParmSet, int iCount )
{
VARARRAY	vaTemp;

	/* Create a new VARARRAY to hold the torsions */
    vaTemp = vaVarArrayCreate( sizeof( TORSIONPARMt ));
    VarArraySetSize( vaTemp, iCount );
    
    VarArrayDestroy( &psParmSet->vaImpropers );
    psParmSet->vaImpropers = vaTemp;
    MESSAGE(("improper parameters size changed to %i\n", iCount ));
}

int
iParmSetProperCount( PARMSET psParmSet )
{
    return( iVarArrayElementCount( psParmSet->vaTorsions ) );
}

int
iParmSetImproperCount( PARMSET psParmSet )
{
    return( iVarArrayElementCount( psParmSet->vaImpropers ) );
}

#if 0
static void
printimp( PARMSET psParmSet )
{
int		i, count;
TORSIONPARMt	*tpPCur;

	count = iVarArrayElementCount( psParmSet->vaImpropers );
        tpPCur = PVAI( psParmSet->vaImpropers, TORSIONPARMt, 0 );
	for ( i=0; i<count; i++, tpPCur++ ) {
		fprintf(stderr, " %s %s %s %s   %d   %f %f   %s   %s\n",
			tpPCur->sType1, tpPCur->sType2, 
			tpPCur->sType3, tpPCur->sType4,
			tpPCur->iN, tpPCur->dKp, tpPCur->dP0/DEGTORAD,
			tpPCur->sOrder, tpPCur->sDesc );
	}
}

static void
printtors( PARMSET psParmSet )
{
int		i, count;
TORSIONPARMt	*tpPCur;

	count = iVarArrayElementCount( psParmSet->vaTorsions );
        tpPCur = PVAI( psParmSet->vaTorsions, TORSIONPARMt, 0 );
	for ( i=0; i<count; i++, tpPCur++ ) {
		fprintf(stderr, " %d  %s %s %s %s   %d   %f %f   %s\n", i+1,
			tpPCur->sType1, tpPCur->sType2, 
			tpPCur->sType3, tpPCur->sType4,
			tpPCur->iN, tpPCur->dKp, tpPCur->dP0/DEGTORAD, 
			tpPCur->sDesc );
	}
}
#endif
