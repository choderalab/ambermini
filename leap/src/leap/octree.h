/*
 *	File:	octree.h
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
 *     Designed by:    Bill Ross					*
 *     Author:         Bill Ross					*
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Bill Ross						*
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *	Description:
 *		Handle octrees
 *
 */

# ifndef OCTREE_H
# define	OCTREE_H

typedef struct oc {
	int		iNodeNum;
	int		iStatus;
	int		iDepth;
	VECTOR		vCorner;
	int             iAtoms;
	ATOM		*PaAtomList;
	float		*PfCharges;	/* float to save space */
	struct oc	*PonChildren;
} OCTNODEt;

typedef OCTNODEt	*OCTNODE;

typedef	struct	{
	int		iType;	
	int		iMaxDepth;
	int		iTreePoints;
	int		iDielectric;
	double		dGridSize;
	int		*PiDensities;
	float		*PfCharges;	/* float to save space */
	double		*PdHalfEdges;
	double		*PdHalfDiagonals;
	VARARRAY        vaAtoms;
	OCTNODEt	onHead;
} OCTREEt;

typedef OCTREEt		*OCTREE;

/* OCTREEt.iOctType types */

#define OCT_SHELL 		1
#define OCT_INTERIOR_SOLUTE 	2
#define OCT_INTERIOR_SOLVENT	3

/* atom selection options */

#define AT_OCTREE	1

/* dielectric options */

#define DIEL_R2		1
#define DIEL_R		2

/* types of printing */

#define COLOR_CUT	1
#define COLOR_RANGE	2
#define COLOR_DEPTH	3
#define COLOR_NONE	4

extern OCTREE	octOctTreeCreate(UNIT uUnit, int iType, double dGridSpace, 
			double dAddExtent, double dShellExtent, int iIncludeSolvent);
extern void	OctTreeDestroy( OCTREE *PoctTree );

extern void	OctTreeInitCharges(OCTREE octTree, int iAtomOption, int iDielectric, 
			double dCutDist, VECTOR *PvMin, VECTOR *PvMax);	

extern void	OctTreeDescribe();	/* ( OCTREE ) */
extern void	OctTreeDeleteSphere(OCTREE octTree, VECTOR *PvPoint, double dRadius);
extern void	OctTreeUpdateCharge(OCTREE octTree, VECTOR *PvNewPoint, 
			float fCharge, double dCutDist, VECTOR *PvMax, VECTOR *PvMin);
extern RESIDUE	rOctTreeCheckSolvent(OCTREE octTree, VECTOR *PvPoint);

extern void	OctTreePrintGrid(OCTREE octTree, char *sFileName, int iColor);

# endif
