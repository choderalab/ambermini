/*
 *      File: object.c
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
 *
 *
 *      Description:
 *              Contains routines that are common to all OBJEKTS and
 *              subclasses.
 *
 *              ALL classes must insert code into these
 *              routines to call their own creator/destructor routines.
 *
 *      OBJEKT HIERARCHY
 *      ----------------
 *
 *      OBJEKT
 *              ASSOC
 *              COLLECTION
 *                      LIST
 *              BYTEARRAY
 *              CONTAINER
 *                      UNIT
 *                      MOLECULE
 *                      RESIDUE
 *                      ATOM
 *              INTERNAL
 *              PARMSET
 *
 *
 */
 



#include	"basics.h"

#include        "classes.h"




#define DEFAULTSIZE     -1




/*
 *      oCreateSize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to create objects.
 *
 *      Return:
 *              Return a pointer to the object created.
 */
OBJEKT
oCreateSize( int iType, int iSize )
{
OBJEKT  o;

    switch ( iType ) {
    
        case OBJEKTid:
            PRINTF(( "A raw OBJEKT is being created!!\n" ));
            MALLOC( o, OBJEKT, sizeof(OBJEKTt) );
            break;

        case ASSOCid:
            o = (OBJEKT)aAssocCreate();
            break;
        
        case OINTEGERid:
            o = (OBJEKT)oiOIntegerCreate();
            break;

        case ODOUBLEid:
            o = (OBJEKT)odODoubleCreate();
            break;

        case OSTRINGid:
            o = (OBJEKT)osOStringCreate();
            break;

        case COLLECTIONid:   
            case LISTid:
            o = (OBJEKT)cCollectionCreate(iType);
            break;
            
        case BYTEARRAYid:
            o = (OBJEKT)baByteArrayCreate( iSize );
            break;
            
        case ATOMid:
        case RESIDUEid:
        case UNITid:
        case MOLECULEid:
        case CONTAINERid:
            o = (OBJEKT)cContainerCreate(iType);
            break;

        case INTERNALid:
            o = (OBJEKT)iInternalCreate();
            break;

        case PARMSETid:
            o = (OBJEKT)psParmSetCreate();
            break;
                    
        default:
            DFATAL( ("Unknown object type=%d being created with size %d!", iType,
                        iSize ) );
            break;
    }
    
        /* Define the type. */
        
    o->cObjType = iType;
    o->iReferences = 1;    
    return(o);
}




/*
 *      oCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to create objects.
 *      This routine simply calls oCreateSize with size of DEFAULTSIZE.
 *
 *      Return:
 *              Return a pointer to the object created.
 */
OBJEKT
oCreate( int iType )
{
OBJEKT  o;

    o = oCreateSize( iType, DEFAULTSIZE );
                
    return(o);
}







/*
 *      Destroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to destroy objects.
 *      After the call to Destroy, the object will be destroyed and
 *      the value pointed to by oObject will be undefined.
 *
 *      Arguments:
 *      oObject -       A pointer to the object.
 */
void
Destroy( OBJEKT  *oPObject )
{
    if ( *oPObject == NULL ) 
	return;

    MESSAGE(( "Destroying object: %c\n", iObjectType(*oPObject) ));

    switch ( iObjectType(*oPObject) ) {
        case OBJEKTid:
            PRINTF(( "A raw OBJEKT is being destroyed!" ));
            FREE( *oPObject );
            break;

        case ASSOCid:
            AssocDestroy( (ASSOC *)oPObject );
            break;

        case OINTEGERid:
            OIntegerDestroy( (OINTEGER *)oPObject );
            break;

        case ODOUBLEid:
            ODoubleDestroy( (ODOUBLE *)oPObject );
            break;

        case OSTRINGid:
            OStringDestroy( (OSTRING *)oPObject );
            break;

        case COLLECTIONid:
        case LISTid:
            CollectionDestroy( (COLLECTION *)oPObject );
            break;

        case BYTEARRAYid:
            ByteArrayDestroy( (BYTEARRAY *)oPObject );
            break;
            
        case ATOMid:
        case RESIDUEid:
        case UNITid:
        case MOLECULEid:
        case CONTAINERid:
                ContainerDestroy( (CONTAINER *)oPObject );
                break;

        case INTERNALid:
            InternalDestroy((INTERNAL *)oPObject);
            break;
            
        case PARMSETid:
            ParmSetDestroy((PARMSET *)oPObject);
            break;
            
        default:
            DFATAL(("Unknown object type=%c/%d addr %p being destroyed!", 
		iObjectType(*oPObject), iObjectType(*oPObject),
		*oPObject ));
            break;
    }
    
}




/*
 *	oCopy
 *
 *	Create a copy of the OBJEKT.
 *	The OBJEKT that is returned is an identical copy.
 *
 *	THIS IS THE ROUTINE THAT APPLICATIONS MUST USE
 *	TO COPY OBJEKTS, 
 *
 *		NOT oObjectDuplicate.
 */
OBJEKT
oCopy( OBJEKT oCur )
{
OBJEKT	oNew;

		/* First duplicate the OBJEKT */

    MESSAGE( ("copying type %s\n", sObjectType(oCur) ));

    oNew = oObjectDuplicate(oCur);

		/* Then reset the pointers if it is a CONTAINER */

    switch ( iObjectType(oNew) ) {
        case ATOMid:
        case RESIDUEid:
        case UNITid:
        case MOLECULEid:
        case CONTAINERid:
            ContainerResetAllCopyPointers( (CONTAINER)oNew );
            break;
        case OBJEKTid:
        case OINTEGERid:
        case ODOUBLEid:
        case OSTRINGid:
        case COLLECTIONid:
        case LISTid:
        case BYTEARRAYid:
        case PARMSETid:
            break;
        default:
            VP0( ("Unknown object type=%s being copied!",
                 sObjectType(oCur)) );
	    Destroy( &oNew );
	    return(NULL);
            break;
    }

    return(oNew);
}
	


/*
 *      Describe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to return a descriptor
 *      string for the object.
 *      It is the callers responsibility to ensure that there is room
 *      enough in the buffer to describe the object.
 *
 *      Arguments:
 *      oObject -       The object to describe.
 *
 */
void
Describe( OBJEKT oObject )
{

    if ( oObject == NULL ) {
        VP0(( "--NULL--\n" ));
        return;
    }

    switch ( iObjectType(oObject) ) {
        case OBJEKTid:
            VP0(( "OBJEKT\n" ));
            break;
            
        case ASSOCid:
            AssocDescribe( (ASSOC)oObject );
            break;

        case OINTEGERid:
            OIntegerDescribe( (OINTEGER)oObject );
            break;

        case ODOUBLEid:
            ODoubleDescribe( (ODOUBLE)oObject );
            break;

        case OSTRINGid:
            OStringDescribe( (OSTRING)oObject );
            break;

        case COLLECTIONid:
        case LISTid:
            CollectionDescribe( (COLLECTION)oObject );
            break;

        case BYTEARRAYid:
            ByteArrayDescribe( (BYTEARRAY)oObject );
            break;

        case ATOMid:
        case RESIDUEid:
        case UNITid:
        case MOLECULEid:
        case CONTAINERid:
                ContainerDescribe( (CONTAINER)oObject );
                break;
                
        case INTERNALid:
            InternalDescribe( (INTERNAL)oObject );
            break;
            
        case PARMSETid:
            ParmSetDescribe( (PARMSET)oObject );
            break;
            
        default:
            DFATAL( ("Unknown object id=%c being described!",
                 iObjectType(oObject)) );
            break;
    }
    
}



/*
 *      bObjectInClass
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the object has iClass as one of
 *      its superclasses.
 *      The ENTIRE class structure has to be represented here.
 */
BOOL
bObjectInClass( OBJEKT oObject, int iClass )
{
int             iObjClass;

    iObjClass = iObjectType(oObject);
    
    switch ( iClass ) {
        case OBJEKTid:
            return(TRUE);
        case ASSOCid:
            if ( iObjClass==ASSOCid ) return(TRUE);
            break;
        case OINTEGERid:
            if ( iObjClass==OINTEGERid ) return(TRUE);
            break;
        case ODOUBLEid:
            if ( iObjClass==ODOUBLEid ) return(TRUE);
            break;
        case OSTRINGid:
            if ( iObjClass==OSTRINGid ) return(TRUE);
            break;
        case BYTEARRAYid:
            if ( iObjClass==BYTEARRAYid ) return(TRUE);
            break;
        case COLLECTIONid:
            switch ( iObjClass ) {
                case COLLECTIONid:
                case LISTid:
                    return(TRUE);
            }
            break;
        case LISTid:
            switch ( iObjClass ) {
                case LISTid:
                    return(TRUE);
            }
            break;
        case CONTAINERid:
            switch ( iObjClass ) {
                case ATOMid:
                case RESIDUEid:
                case UNITid:
                case MOLECULEid:
                case CONTAINERid:
                    return(TRUE);
            }
            break;
        case UNITid:
            switch ( iObjClass ) {
                case UNITid:
                    return(TRUE);
            }
            break;
        case MOLECULEid:
            switch ( iObjClass ) {
                case MOLECULEid:
                    return(TRUE);
            }
            break;
        case RESIDUEid:
            switch ( iObjClass ) {
                case RESIDUEid:
                    return(TRUE);
            }
            break;
        case ATOMid:
            switch ( iObjClass ) {
                case ATOMid:
                    return(TRUE);
            }
            break;
        case INTERNALid:
            switch ( iObjClass ) {
                case INTERNALid:
                    return(TRUE);
            }
            break;
        case PARMSETid:
            switch ( iObjClass ) {
                case PARMSETid:
                    return(TRUE);
            }
    }
    return(FALSE);
}
          
            

/*
 *      oObjectDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Call the specific routine for the class to duplicate objects
 *	and update objekt-level attributes (iReferences).  The
 *	<class>Duplicate() routines must ONLY be called from this one
 *	to guarantee the objekt-level stuff is right.
 *      The entire OBJEKT MUST be duplicated, including OBJEKTS contained
 *      within.  Internal pointers must also be pointed to within
 *      the duplicate object.
 *
 *      Return:
 *              Return a pointer to the duplicated object.
 */
OBJEKT
oObjectDuplicate( OBJEKT oOld )
{
OBJEKT  o;

    switch ( iObjectType(oOld) ) {
    
        case ATOMid:
        case RESIDUEid:
        case UNITid:
        case MOLECULEid:
        case CONTAINERid:
            o = (OBJEKT)cContainerDuplicate((CONTAINER)oOld);
            break;

        case COLLECTIONid:   
	case LISTid:
            o = (OBJEKT)cCollectionDuplicate((COLLECTION)oOld);
            break;
            
        case OBJEKTid:
            DFATAL( ( "A raw OBJEKT is being duplicated!!" ) );
            break;
 
        case OINTEGERid:
            o = (OBJEKT)oiOIntegerDuplicate( (OINTEGER)oOld );
            break;

        case ODOUBLEid:
            o = (OBJEKT)odODoubleDuplicate( (ODOUBLE)oOld );
            break;

        case OSTRINGid:
            o = (OBJEKT)osOStringDuplicate( (OSTRING)oOld );
            break;

        case BYTEARRAYid:
            o = (OBJEKT)baByteArrayDuplicate((BYTEARRAY)oOld);
            break;
            
        case PARMSETid:
            o = (OBJEKT)psParmSetDuplicate((PARMSET)oOld);
            break;
            
            
        default:
            DFATAL( ("Unknown object id=%c being duplicated", 
                        iObjectType(oOld) ) );
            break;
    }
    o->iReferences = 1;
    
    return(o);
}



/*
 *      sObjectIndexType
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a pointer to a string that represents the type.
 */
char *
sObjectIndexType( int iType )
{

    switch ( iType ) {
        case OBJEKTid:          return("Object");
        case UNITid:            return("Unit");
        case MOLECULEid:        return("Molecule");
        case RESIDUEid:         return("Residue");
        case ATOMid:            return("Atom");
        case ASSOCid:           return("Association");
        case COLLECTIONid:      return("Collection");
        case LISTid:            return("List");
        case BYTEARRAYid:       return("ByteArray");
        case INTERNALid:        return("Internal");
        case PARMSETid:         return("ParameterSet");
        case OSTRINGid:         return("String");
        case ODOUBLEid:         return("Double");
        case OINTEGERid:        return("Integer");
        default: return("?? Unknown type ??");
    }
    return("?? Unknown type ??");	/* for lint */
}



/*
 *	sObjectType
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a string indicating the type of the OBJEKT,.
 */
char *
sObjectType( OBJEKT oObj )
{
    return(sObjectIndexType(iObjectType(oObj)));
}



/*
 *      bObjektWarnType
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Check the type of the argument, print an error and return FALSE
 *      if it is not the correct type.
 */
BOOL
bObjektWarnType( OBJEKT oObj, int iType )
{

    if ( iObjectType(oObj) == iType ) return(TRUE);
    VP0(( "The value must be of the type: %s\n", sObjectIndexType(iType) ));
    return(FALSE);
}





