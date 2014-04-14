/*
 *      File: classes.h
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
 *      Includes all of the header files required for OOP programming.
 */


#ifndef	CLASSES_H
#define CLASSES_H

#include        "basics.h"
#include        "vector.h"
#include        "varArray.h"

#include        "objekt.h"
#include        "assoc.h"
#include        "oString.h"
#include        "oInteger.h"
#include        "oDouble.h"
#include        "collection.h"
#include        "list.h"
#include        "byteArray.h"

#include        "parmSet.h"
#include        "container.h"
#include        "atom.h"
#include        "residue.h"
#include        "molecule.h"
#include        "unit.h"
#include        "internal.h"

#include        "loop.h"
#include	"parmLib.h"  /* Some class stuff requires PARMLIBRARYs */


#endif  /* CLASSES_H */
