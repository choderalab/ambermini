/*
 *      File: oDouble.h
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
 *              ODOUBLE
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              ODOUBLE is an object that has the same properties
 *              as normal doubles.  It is only declared to
 *              make passing arguments to functions as homogenous as
 *              possible.  It should not be used for any other purpose. 
 *
 *
 */
 
#ifndef ODOUBLE_H
#define ODOUBLE_H

/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/


typedef struct  {
	OBJEKTt         oSuper;
	double          dValue;
} ODOUBLEt;

typedef ODOUBLEt	*ODOUBLE;



/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */


extern ODOUBLE		odODoubleCreate();
extern void		ODoubleDestroy(ODOUBLE *odPDbl);
extern void		ODoubleDescribe(ODOUBLE odDbl);
extern ODOUBLE		odODoubleDuplicate( ODOUBLE od );

#define dODouble( OD )  ( ((ODOUBLE)OD)->dValue )
#define ODoubleSet( OD, V ) ( ((ODOUBLE)OD)->dValue = (V) )


#endif /* ODOUBLE_H */
