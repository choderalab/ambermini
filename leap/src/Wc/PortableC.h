#include "COPY.h"
#ifndef _PortableC_h
#define _PortableC_h

/*
* SCCS_data: %Z% %M% %I% %E% %U%
*
* This include file can be used to support portability between K&R C, ANSI C,
* and C++ environments.  This include file should be included at the beginning
* of all include files which use these portability macros.
*
* Four macros are provided:
*	portable prototypes		_()
*	portable varargs		Va_start()
*	portable external linkage	BEGIN_NOT_Cxx ... END_NOT_Cxx
*	portable const			ANSI_const
*
* In order to acheive maximum portabilty, ALL of the following techniques
* should be used.  Your portability decreases whenever any of these techniques
* are not used:
*
* 1) ALWAYS USE PROTOTYPE DECLARATIONS -- This allows your compiler to perform
*    argument type checking, if it can.
*
* 2) USE K&R STYLE DEFINITIONS EXCEPT WITH `...' -- This ensures you do not use
*    non-portable arguments to functions.
*
* 3) ALWAYS USE Va_start() IN VARARGS FUNCTIONS -- This macro works for both
*    K&R and ANSI/C++ environments
*
* 4) ALWAYS EXPORT C INTERFACES -- No C++ interface should be provided by any
*    library, even if the library is implemented in C++ for use by C++ clients:
*    otherwise, ANY change to ANY member (data or function, public or private)
*    of ANY class (base or derived) will force ALL clients to RE-COMPILE not
*    simply re-link!
*
* 5) Use ANSI_const for constant argument declarations - e.g., functions
*    invoked by qsort() should have a prototype like this:
*	int QsortCompare _(( ANSI_const void* left, ANSI_const void* right ));

*--------------------------------
* Portabilty Techniques In Detail
*--------------------------------
*
* 1) ALWAYS USE PROTOTYPE DECLARATIONS
*	You should ALWAYS provide prototype declarations for EVERY function,
*	even static functions.  These prototype declarations should ALWAYS be
*	seen by compilers BEFORE any invocation of the function.  This allows
*	your compiler to perform argument type checking, if it can.
*
*	For example:
*
*	    extern void WcAttachThisToWidget    _(( XtPointer, char*, Widget ));
*	    static String XmpxMakeEntireMessage _(( String, va_list ));
*	    static void XmpxPrint               _(( String, ... ));
*
*	NOTE: There MUST be whitespace between the function name and the `_'
*	and there MUST NOT be whitespace between the `_' and the `('.  These
*	DO NOT WORK as expected:
*
*	    static void WillNotWork_(( String, ... ));
*	    static void WillNotWork_ (( String, ... ));
*	    static void WillNotWork _ (( String, ... ));
*	    

* 2) USE K&R STYLE DEFINITIONS EXCEPT WITH `...'
*	This ensures you do not use non-portable arguments to functions.  If
*	you try to declare a function which takes non-portable arguments, such
*	as a short rather than an int, or a float rather than a double, then
*	you will get something like the following from ANSI C or C++ compilers:
*		ERROR: Function definition does not match prototype
*
*	For example:
*
*	    extern void WcAttachThisToWidget( this, className, widget )
*		XtPointer this;
*		char*     className;
*		Widget    widget;
*	    {
*		...etc...
*	    }
*
*	    static String XmpxMakeEntireMessage( string, argPtr )
*		String string;
*		va_list argPtr;
*	    {
*		...etc...
*	    }
*
*	    ****************************************************************
*	    ** THE ONLY TIME YOU USE BOTH K&R AND ANSI-STYLE DEFINITIONS ***
*	    ****************************************************************
*	    #if NeedFunctionPrototypes
*	    static void XmpxPrint( String msg, ... )
*	    #else
*	    static void XmpxPrint( msg, va_alist )
*	        String      msg;
*	        va_dcl
*	    #endif
*	    {
*	        va_list     argPtr;
*
*	        Va_start( argPtr, msg );      // point at 1st unnamed arg
*
*		...etc...
*	    }
*

* 3) ALWAYS USE Va_start() IN VARARGS FUNCTIONS
*	K&R uses a different macro for beginning processing of variable
*	argument lists.  The Va_start() macro works for K&R as well as newer
*	environments.  The example above shows the use of Va_start().
*
*	Also, note that this include file imports the appropriate varags
*	include files for K&R and newer environments.
*

* 4) ALWAYS EXPORT C INTERFACES
*	The most obvious reason for this is so you can export capability from
*	your C++ code to clients which are not written in C++.
*
*	However, there are severe problems with even exporting C++ capabilities
*	to C++ clients:
*	
*	First, C++ interfaces are very inconsistent across compiler vendors,
*	and even between different versions from the same vendor.  Any given
*	C++ interface can only be reliably used by clients using the same
*	version of the same compiler as that used to build the library.
*
*	Worse, the C++ interfaces are only transparent at the source code
*	level: the client compiler must know the binary layout of each class.
*	The binary layout of classes is a function of the public AND private
*	parts of the class, and the existence of virtual functions.  
*
*	If even a private virtual function of a base class of an exported
*	class is changed, added, or deleted, ALL method invocations of ALL
*	derived classes will need to be RE-COMPILED, not simply re-linked.
*
*	If even a private member of a base class of an exported class is
*	changed, added, or deleted, ALL member accesses of ALL derived classes
*	will need to be RE-COMPILED, not simply re-linked.
*
*	In your public interface include files:
*
*	    BEGIN_NOT_Cxx
*	    ... all external interface declarations ...
*	    END_NOT_Cxx
*
*	In your source files which implement the public interfaces:
*
*	    BEGIN_NOT_Cxx
*	    ... all external interface definitions ...
*	    END_NOT_Cxx

* 5) USE ANSI_const FOR CONSTANT ARGUMENT DECLARATIONS
*	ANSI and C++ compilers accept the `const' keyword to declare read-only
*	arguments to functions.  qsort() is a commonly used capability which
*	requires declaring a function with `const' arguments.  In order to
*	allow the same source code to compile cleanly on all C environments,
*	you need to declare and define such functions using the ANSI_const
*	macro.  For example:
*
*	    int Compare _(( ANSI_const void* left, ANSI_const void* right ));
*
*	    int Compare( left, right )
*		ANSI_const void* left;
*		ANSI_const void* right;
*	    {
*		char* leftStr  = *((char**)left);
*		char* rightStr = *((char**)right);
*		return strcmp( leftStr, rightStr );
*	    }
*/	

#ifndef NeedFunctionPrototypes
#if defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus) || defined(FUNCPROTO)
#define NeedFunctionPrototypes 1
#else
#define NeedFunctionPrototypes 0
#endif /* __STDC__ */
#endif /* notdef NeedFunctionPrototypes */

/* Macro for ANSI or K&R external declarations.  Declare them like this:
**
**      int foo _(( int, char* ));
**
***************** DO NOT forget whitespace before the '_' !! *****************
*/
#ifndef _
#if NeedFunctionPrototypes
#define	_(argsInParens)	argsInParens	/* ANSI --> int foo ( int, char* ); */
#else
#define	_(argsInParens)	()		/* K&R  --> int foo ();             */
#endif
#endif

#if NeedFunctionPrototypes
#include <stdarg.h>
#define Va_start(a,b) va_start(a,b)
#else
#include <varargs.h>
#define Va_start(a,b) va_start(a)
#endif /* NeedFunctionPrototypes */

#ifdef __cplusplus  /* for C++ V2.0 */
#define BEGIN_NOT_Cxx	extern "C" {
#define END_NOT_Cxx	}
#else
#define BEGIN_NOT_Cxx
#define END_NOT_Cxx
#endif

#if NeedFunctionPrototypes
#define ANSI_const const
#else
#define ANSI_const
#endif

#endif /* _PortableC_h */
