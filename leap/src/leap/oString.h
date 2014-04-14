/*
 *      File: oString.h
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
 *              OSTRING
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              An OSTRING is a object which contains strings.
 *              It has all of the same properties as normal C strings
 *              but must be accessed as an OBJEKT.  OSTRINGs are 
 *              not supposed to replace C strings, they are only provided
 *              to maintain consistancy in the command line interface.
 *              Since the command line interpreter will be dealing with
 *              MOLECULEs, UNITs, RESIDUEs, ATOMs, LISTs, and also
 *              integers, doubles, and strings there should be
 *              a consistant way of handling these things.
 *
 *
 */
 
#ifndef OSTRING_H
#define OSTRING_H

/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/

typedef struct  {
	OBJEKTt         oSuper;
	int             iMaxLen;
	char		*sString;
} OSTRINGt;

typedef OSTRINGt	*OSTRING;



/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */


extern OSTRING	osOStringCreate();
extern void	OStringDestroy(OSTRING *osPStr);
extern void	OStringDescribe(OSTRING osStr);

extern void	OStringCopy(OSTRING osDest, OSTRING osSrc);
extern void	OStringDefine(OSTRING osStr, char *sStr);
extern void	OStringConcat(OSTRING osA, OSTRING osB);
extern OSTRING	osOStringDuplicate( OSTRING os );

#define iOStringCompare( OS1, OS2 )     \
( strcmp( ((OSTRING)OS1)->sString, ((OSTRING)OS2)->sString )
#define sOString( OS )  ( ((OSTRING)OS)->sString ? ((OSTRING)OS)->sString : "" )


#endif /* OSTRING_H */
