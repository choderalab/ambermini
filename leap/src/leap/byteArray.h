/*
 *      File: byteArray.h
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
 *              BYTEARRAY
 *      Superclass: 
 *              OBJEKT
 *
 *
 */

#ifndef BYTEARRAY_H
#define BYTEARRAY_H


/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/

typedef struct  {
	OBJEKTt         oObject;
	int             iDataSize;
} BYTEARRAYHEADERt;
                
typedef struct  {
	BYTEARRAYHEADERt        bahHeader;
	char                    baData[1];
} BYTEARRAYt;


typedef BYTEARRAYt      *BYTEARRAY;


/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/


#define iByteArraySize(ba)      (((BYTEARRAY)ba)->bahHeader.iDataSize)
#define PByteArray(ba)          (((BYTEARRAY)ba)->baData)


/*  byteArray.c  */

extern BYTEARRAY	baByteArrayCreate( int iSize );
extern void		ByteArrayDestroy( BYTEARRAY *baPByteArray );
extern void		ByteArrayDescribe( BYTEARRAY baByteArray );
extern BYTEARRAY	baByteArrayDuplicate( BYTEARRAY baOld );


#endif /* BYTEARRAY_H */
