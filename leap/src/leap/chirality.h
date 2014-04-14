/*
 *	File:	chirality.h
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
 *	Description:
 *		Routines to calculate chirality around an ATOM, and
 *		to calculate the orientation of an ATOM with respect to
 *		a plane formed by three other ATOMs.
 *
 *		Chirality is measured using at least three ATOMs
 *
 *			B     A
 *			 \   /
 *			  \ /
 *			   C
 *			   |
 *			   |
 *			   D	Out or into the plane?
 *
 *		The ATOM (D) can be placed either above the plane defined
 *		by (A)-(C)-(D) or below.  The chirality of (C) is calculated
 *		by ordering the ATOMs around it and then crossing (A)-(C) into
 *		(B)-(C) and doing a dot product with (D)-(C).  If the result
 *		is positive then the chirality is +1.0 otherwise -1.0.
 *		The ordering of the ATOMs is determined by sorting them
 *		using the ATOM ID as the key.
 *
 *		Often it is important to calculate the orientation of
 *		an ATOM with respect to three others.  This can be done
 *		by transforming the chirality calculated with respect
 *		to the above ordering into an -orientation- with respect
 *		to a different ordering.
 *
 *		Chirality is defined as the value obtained when the ordering
 *			of the ATOMs is with respect to increasing ATOM ID.
 *		Orientation is defined as the value obtained when the
 *			ordering is arbitrary.
 */

/*  chirality.c  */

extern void		ChiralityOrderNeighbors( ATOM aAtom,
				ATOM *aPAtomA, ATOM *aPAtomB,
				ATOM *aPAtomC, ATOM *aPAtomD );
extern double		dChiralityForAtom( ATOM aAtom );
extern double		dChiralityToOrientation( double dChirality, 
				ATOM aCenter,
				ATOM aAtomA, ATOM aAtomB, 
				ATOM aAtomC, ATOM aAtomD );
extern double		dChiralityFromOrientation( double dOrient, 
				ATOM aCenter,
				ATOM aAtomA, ATOM aAtomB, 
				ATOM aAtomC, ATOM aAtomD );

