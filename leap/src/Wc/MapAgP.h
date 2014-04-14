#ifndef _MapAgP_h_
#define _MapAgP_h_
#include "COPY.h"

/*
* SCCS_data: %Z% %M% %I% %E% %U%
*
* Mapping Agent - MapAgP.h - Private Interface
*
* This include file provides declarations private to the implementation
* of Mapping Agents.
*
* If a client needs MapAg_STATIC (a macro to statically allocate and initialize
* Mapping Agent) then this file must be included.  Note that MapAg_STATIC
* breaks the concept of binary encapsulation: in other words, a client using
* MapAg_STATIC must RE-COMPILE with every new version of MapAg.
*
* An agent has 4 states: Initial and Normal, InitialStatic and NormalStatic.
* If an agent is created using MapAg_New() - which is the normal and
* preferred approach - then the agent will transition between the Initial
* and Normal stated.  If an agent is statically created using MapAg_STATIC
* then the agent will transition between the InitialStatic and NormalStatic
* states.  The states provide the following behaviors:
*
*       Initial:
*		MapAg_Transition_InitialToNormal,
*               MapAg_Define_Initial, MapAg_Find_Initial, 
*               MapAg_Forget_Initial, MapAg_Purge_Initial,
*		MapAg_Resize_Initial, MapAg_Free_Dynamic
*       Normal:
*		MapAg_Transition_NormalToInitial,
*               MapAg_Define_Normal, MapAg_Find_Normal, 
*               MapAg_Forget_Normal, MapAg_Purge_Normal,
*		MapAg_Resize_Normal, MapAg_Free_Dynamic
*       InitialStatic:
*		MapAg_Transition_InitialToNormalStatic,
*               MapAg_Define_Initial, MapAg_Find_Initial, 
*               MapAg_Forget_Initial, MapAg_Purge_Initial,
*		MapAg_Resize_Initial, MapAg_Free_Static
*       NormalStatic:
*		MapAg_Transition_NormalToInitialStatic,
*               MapAg_Define_Normal, MapAg_Find_Normal, 
*               MapAg_Forget_Normal, MapAg_Purge_Normal,
*		MapAg_Resize_Normal, MapAg_Free_Static
*
* State Transitions:
*-------------------
*
* Creating a new agent results in an agent in the Initial state.
*
* Defining a new mapping to an Initial agent causes the transition to Normal
* if sucessful (i.e., all the allocations succeed).  If any malloc fails, then
* XtError is called.
*
* From the Normal state, Forgeting the last mapping causes the agent to
* revert to the initial state.
*
* State Implementation:
*----------------------
*
* The technique used is to declare a state object type which will contain the
* pointers to the functions which implement the behaviors.  Four such objects
* are statically allocated and initialized, one for each state.  A state
* transition then simply involves setting the value of the MapAg instance
* member.
*
*******************************************************************************
*/

#include "MapAg.h"

BEGIN_NOT_Cxx

/*  -- Map Object Declaration - Private Object
*******************************************************************************
    Mapping agents implement a hashing scheme.  The a,b,c values are
    hashed, the result is the index into an array of pointers to Map
    objects.  Map objects form a linked list to take care of collisions.
    Note that we store the a,b,c values in the Map object to be certain
    we find the appropriate data.
*/
typedef struct _MapRec	{	/* Stores one entry. */
    char* 		a;
    char*		b;
    char*		c;
    char*		data;
    struct _MapRec*	next;
} MapRec, *Map;


/*  -- State Method Type Declarations
*******************************************************************************
*/

typedef void  (*TransitionMethod) _(( MapAg ));
typedef void  (*DefineMethod)	  _(( MapAg, char*, char*, char*, char* ));
typedef char* (*FindMethod)	  _(( MapAg, char*, char*, char* ));
typedef void  (*ForgetMethod)	  _(( MapAg, char*, char*, char* ));
typedef void  (*PurgeMethod)	  _(( MapAg, char* ));
typedef void  (*ResizeMethod)	  _(( MapAg ));
typedef void  (*FreeMethod)	  _(( MapAg ));

/*  -- State Object Declaration - Private Object
*******************************************************************************
   Every MapAg instance will point to a State object in order to get its
   state dependent behavior.
*/
typedef struct _StateRec {
    TransitionMethod	Transition;
    DefineMethod	Define;	/* Add/change the data mapped to a,b,c	*/
    FindMethod		Find;	/* find data mapped to a,b,c		*/
    ForgetMethod	Forget;	/* forget any mapping of a,b,c		*/
    PurgeMethod		Purge;	/* forget any mapping to data		*/
    ResizeMethod	Resize;	/* PRIVATE resize table method 		*/
    FreeMethod		Free;	/* Free a MapAg and all its mappings	*/
} StateRec, *State;


/*  -- Mapping Agent Object Declaration
*******************************************************************************
*/
typedef struct _MapAgRec {	/* Stores hash table for mapping.	*/
    State	state;		/* Current state			*/
    Map*	table;		/* Pointer to hash table of Maps.	*/
    int		mask;		/* Current size of hash table minus 1.	*/
    int		collisionLimit;	/* keep collisions under this limit	*/
    int		numMaps;	/* Maps currently in hash table.	*/
    int		shiftA;		/* used by Hash: num insig LSBs in `a'	*/
    int		shiftB;		/* used by Hash: num insig LSBs in `b'	*/
    int		shiftC;		/* used by Hash: num insig LSBs in `c'	*/
} MapAgRec;

/*  -- Macro for Static Initialization
*******************************************************************************
    This MUST be kept compatible with MapAg_New()
*/
extern StateRec State_initialStaticRec;

#define MapAg_STATIC {\
/* State state		*/ &State_initialStaticRec, \
/* Map*  table		*/ (Map*)0, \
/* int   mask		*/ (int)0,\
/* int   collisionLimit	*/ (int)0,\
/* int   numMaps	*/ (int)0,\
/* int   shiftA		*/ (int)0,\
/* int   shiftB		*/ (int)0,\
/* int   shiftC		*/ (int)0 \
}

/* -- In-line Method Invocations
*******************************************************************************
   This breaks binary encapsulation!  However, any code using this private
   include file already does not use binary encapsulation.

   Note that we also provide invocation macros for the private methods
   Resize and Transition.
*/
#undef MapAg_Define
#undef MapAg_Find
#undef MapAg_Forget
#undef MapAg_Purge
#undef MapAg_Free

#ifndef ASSERTIONS 
#define MapAg_Define(ag,a,b,c,d) \
	((ag)->state->Define((ag),(char*)(a),(char*)(b),(char*)(c),(char*)(d)))

#define MapAg_Find(ag,a,b,c) \
	((ag)->state->Find((ag),(char*)(a),(char*)(b),(char*)(c)))

#define MapAg_Forget(ag,a,b,c) \
	((ag)->state->Forget((ag),(char*)(a),(char*)(b),(char*)(c)))

#define MapAg_Purge(ag,data) \
	((ag)->state->Purge((ag),(char*)(data)))

#define MapAg_Free(ag) \
	((ag)->state->Free(ag))

#define MapAg_Resize(ag) \
	((ag)->state->Resize(ag))

#define MapAg_Transition(ag) \
	((ag)->state->Transition(ag))

#else
#define MapAg_Define(ag,a,b,c,d) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	((ag)->state->Define((ag),(char*)(a),(char*)(b),(char*)(c),(char*)(d))))

#define MapAg_Find(ag,a,b,c) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	((ag)->state->Find((ag),(char*)(a),(char*)(b),(char*)(c))))

#define MapAg_Forget(ag,a,b,c) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	((ag)->state->Forget((ag),(char*)(a),(char*)(b),(char*)(c)))

#define MapAg_Purge(ag,data) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	((ag)->state->Purge((ag),(char*)(data))))

#define MapAg_Free(ag) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	((ag)->state->Free(ag)))

#define MapAg_Resize(ag) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	((ag)->state->Resize(ag)))

#define MapAg_Transition(ag) \
	(MapAg_AssertLooksOk((ag),__FILE__,__LINE__), \
	((ag)->state->Transition(ag)))

#endif

END_NOT_Cxx

#endif /* _MapAgP_h_ */
