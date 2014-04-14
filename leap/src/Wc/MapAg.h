#ifndef _MapAg_h_
#define _MapAg_h_
#include "COPY.h"

/*
* SCCS_data: %Z% %M% %I% %E% %U%
*
* Mapping Agent - MapAg.h - Public Interface
*
* This include file describes an agent for associating arbitrary data.  In all
* cases, the data consists of char* (pointers to anything).
*
* An adaptative hashing algorithm is used: when multiple collisions occur,
* an attempt is made to modify the hashing algorithm to reduce collisions.
*
* NOTE: Most methods are invoked using macros which provide the notational
* convenience of casting the arguments to char* so the client does not need to.
* The functions invoked by these macros are public, of course.
*
* Public Interface:
*
*	MapAg			Opaque Pointer to Mapping Agent Object.
*
*	MapAg_New()		Create a new mapping agent.
*
*	MapAg_Free()		Destroy an existing map database.
*
*	MapAg_Define()		Set mapping of (a,b,c) to data - replaces any 
*				existing mapping of (a,b,c) to anything.
*
*	MapAg_Find()		Find data associated with (a,b,c).
*
*	MapAg_Forget()		Forget association from (a,b,c) to anything.
*
*	MapAg_Purge()		Forget all associations to a `data' pointer.
*
*
* INCOMPATIBILITIES WITH EARLIER VERSIONS:
* To support encapsulation, MapAg_STATIC is no longer publically supported.
*
* Clients must re-compile to use this version, but this version supports
* binary compatibility with future enhancements.
*
*******************************************************************************
*/

#include "PortableC.h"		/* _() macro for portable prototypes */

BEGIN_NOT_Cxx

/*  -- Opaque Pointer to MapAg Object
*******************************************************************************
*/
typedef struct _MapAgRec* MapAg;


/*  -- MapAg Public Methods
*******************************************************************************
*/
extern MapAg MapAg_New      _(());		/* Create a new agent        */
extern void  MapAg_DoFree   _(( MapAg ));	/* Destroy an existing agent */

extern void  MapAg_DoDefine _(( MapAg, char* a, char* b, char* c, char* data ));
extern char* MapAg_DoFind   _(( MapAg, char* a, char* b, char* c ));
extern void  MapAg_DoForget _(( MapAg, char* a, char* b, char* c ));
extern void  MapAg_DoPurge  _(( MapAg, char* data ));

extern void	MapAg_AssertLooksOk _(( MapAg, char*, int ));
/****** usage:	MapAg_AssertLooksOk( agent, __FILE__, __LINE__ ); ******/


/*  -- Macros for Invoking Mapping Agent Functions
*******************************************************************************
    These macros perform casts of arguments so the client does not need to.
    This is only for notational convenience.  If these macros were C functions,
    their prototypes would be:

	void  MapAg_Define _(( MapAg, char* a, char* b, char* c, char* data ));
	char* MapAg_Find   _(( MapAg, char* a, char* b, char* c ));
	void  MapAg_Forget _(( MapAg, char* a, char* b, char* c ));
	void  MapAg_Purge  _(( MapAg, char* data ));
	void  MapAg_Free   _(( MapAg ));

    The assertions are enabled by defining ASSERTIONS (see assert(3)).  In
    this case, the mapping agent argument is referenced twice, so you must NOT
    use side effects!
*/

#ifndef ASSERTIONS 
#define MapAg_Define(ag,a,b,c,d) \
	(MapAg_DoDefine((ag),(char*)(a),(char*)(b),(char*)(c),(char*)(d)))

#define MapAg_Find(ag,a,b,c) \
	(MapAg_DoFind((ag),(char*)(a),(char*)(b),(char*)(c)))

#define MapAg_Forget(ag,a,b,c) \
	(MapAg_DoForget((ag),(char*)(a),(char*)(b),(char*)(c)))

#define MapAg_Purge(ag,data) \
	(MapAg_DoPurge((ag),(char*)(data)))

#define MapAg_Free(ag) \
	(MapAg_DoFree(ag))

#else
#define MapAg_Define(ag,a,b,c,d) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	(MapAg_DoDefine((ag),(char*)(a),(char*)(b),(char*)(c),(char*)(d))))

#define MapAg_Find(ag,a,b,c) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	(MapAg_DoFind((ag),(char*)(a),(char*)(b),(char*)(c))))

#define MapAg_Forget(ag,a,b,c) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	(MapAg_DoForget((ag),(char*)(a),(char*)(b),(char*)(c)))

#define MapAg_Purge(ag,data) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	(MapAg_DoPurge((ag),(char*)(data))))

#define MapAg_Free(ag) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	(MapAg_DoFree(ag)))
#endif

END_NOT_Cxx

#endif /* _MapAg_h_ */
