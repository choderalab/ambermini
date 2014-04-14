/*
 *      File:   nVector.h
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
 *      Description:
 *              Maintain N dimensional vectors.
 */



#include	"basics.h"

#include        "nVector.h"



/*
 *      nvNVectorCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Allocate memory for the NVECTOR and initialize it
 *      to zero.
 */
NVECTOR
nvNVectorCreate( int iElements )
{
NVECTOR nvNew;

                /* allocate memory for iElements doubles and */
                /* one integer to store the number of elements */


    MESSAGE(( "|||About to MALLOC nvNew\n" ));
    TESTMEMORY();
    MESSAGE(( "|||Passed\n" ));
    MALLOC( nvNew, NVECTOR, sizeof(NVECTORt) + sizeof(double)*iElements-1 );
    MESSAGE(( "|||Just MALLOC nvNew\n" ));
    TESTMEMORY();
    MESSAGE(( "|||Passed\n" ));

    nvNew->iSize = iElements;
    MESSAGE(( "|||Just set number of elements\n" ));
    TESTMEMORY();
    MESSAGE(( "|||Passed\n" ));
    NVectorZero(nvNew);
    MESSAGE(( "|||Just NVectorZero\n" ));
    TESTMEMORY();
    MESSAGE(( "|||Passed\n" ));
    return(nvNew);
}


/*
 *      NVectorDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy an NVECTOR
 */
void    
NVectorDestroy( NVECTOR *nvPX )
{
    FREE( *nvPX );
    *nvPX = NULL;
}



/*
 *      nvNVectorCopy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a new vector and copy the contents of nvVector into it.
 */
NVECTOR 
nvNVectorCopy( NVECTOR nvVector )
{
int             i;
NVECTOR         nvNew;

    nvNew = nvNVectorCreate(iNVectorSize(nvVector));
    for ( i=0; i<iNVectorSize(nvVector); i++ ) 
        NVectorSetElement( nvNew, i, dNVectorElement( nvVector, i ) );
    return(nvNew);
}





/*
 *      NVectorCopy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Copy the contents of nvVector into another.
 */
void    
NVectorCopy( NVECTOR nvDest, NVECTOR nvSource )
{
int             i;

    nvDest->iSize = nvSource->iSize;
    for ( i=0; i<iNVectorSize(nvSource); i++ ) 
        NVectorSetElement( nvDest, i, dNVectorElement( nvSource, i ) );
}



/*
 *      NVectorZero
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the vector with zeros
 */
void
NVectorZero( NVECTOR nvVector )
{
    NVectorFill( nvVector, 0.0 );
}



/*
 *      NVectorFill
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the vector with a number.
 */
void
NVectorFill( NVECTOR nvVector, double dValue )
{
int             i;

    for ( i=0; i<iNVectorSize(nvVector); i++ ) 
        NVectorSetElement( nvVector, i, dValue );
}


/*
 *      NVectorDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe the NVECTOR to stdout.
 */
void
NVectorDescribe( NVECTOR nvVector )
{
int             i;

    PRINTF(( "NVECTOR length=%d\n", iNVectorSize(nvVector) ));
    for ( i=0; i<iNVectorSize(nvVector); i++ ) 
        PRINTF_NO_PREFIX(( "---[%d] = %lf\n", i, dNVectorElement(nvVector,i) ));
    PRINTF_NO_PREFIX(( "\n" ));
}


/*
 *      NVectorAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add two vectors and return the result in a third.
 *      A = X + Y
 */
void
NVectorAdd( NVECTOR nvA, NVECTOR nvX, NVECTOR nvY )
{
int             i;

    for ( i=0; i<iNVectorSize(nvX); i++ )
        NVectorSetElement( nvA, i, 
                dNVectorElement(nvX,i) + dNVectorElement(nvY,i) );
}



/*
 *      NVectorSub
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Subtract one vector from another and return the result in a third.
 *      A = X - Y
 */
void
NVectorSub( NVECTOR nvA, NVECTOR nvX, NVECTOR nvY )
{
int             i;

    for ( i=0; i<iNVectorSize(nvX); i++ )
        NVectorSetElement( nvA, i, 
                dNVectorElement(nvX,i) - dNVectorElement(nvY,i) );
}






/*
 *      NVectorTimesScalar
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Multiply a vector by a scalar and put the result in another vector.
 *      A = X * s 
 */
void
NVectorTimesScalar( NVECTOR nvA, NVECTOR nvX, double dS )
{
int             i;

    for ( i=0; i<iNVectorSize(nvX); i++ )
        NVectorSetElement( nvA, i, 
                dNVectorElement(nvX,i) * dS );
}






/*
 *      dNVectorDot
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the DOT product of two vectors.
 *      return:  X . Y 
 */
double
dNVectorDot( NVECTOR nvX, NVECTOR nvY )
{
int             i;
double          dDot;

    dDot = 0.0;
    for ( i=0; i<iNVectorSize(nvX); i++ )
        dDot += dNVectorElement( nvX, i ) * dNVectorElement( nvY, i );
    return(dDot);
}





/*
 *      dNVectorLen
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the length of the vector.
 */
double
dNVectorLen( NVECTOR nvX )
{
double  dDot;

    dDot = dNVectorDot( nvX, nvX );
    return(sqrt(dDot));
}

