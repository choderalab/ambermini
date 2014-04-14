/*
 *      File:   assoc.h
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
 *              ASSOC
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              An ASSOC (short for association) is an object that
 *              stores a short name and an OBJEKT.
 *              It associates the name with the OBJEKT.
 *              This is used by the command line interface to
 *              store objects with the names that are bound to them.
 *
 */
 

#ifndef ASSOC_H
#define ASSOC_H


/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/

typedef struct  {
        OBJEKTt         oSuper;
        STRING          sName;
        OBJEKT          oObj;
} ASSOCt;

typedef ASSOCt  *ASSOC;




/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primitive superclass
        of this class that it will eventually make it into these routines.
*/


#define AssocSetName(a,n)       (strcpy(((ASSOC)(a))->sName,n))
#define sAssocName(a)           (((ASSOC)(a))->sName)
#define AssocSetObject(a,o)     { OBJEKT zzo; zzo = (OBJEKT)(o);\
                                ((ASSOC)(a))->oObj = zzo; \
                                REF(zzo); }
#define oAssocObject(a)         ( ((ASSOC)(a))->oObj )


/*  assoc.c  */

extern ASSOC            aAssocCreate();
extern void             AssocDestroy( ASSOC *aPAssoc );
extern void             AssocDescribe( ASSOC aAssoc );


#endif /* ASSOC_H */
