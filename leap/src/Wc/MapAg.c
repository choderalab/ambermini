#include "COPY.h"

/*
* SCCS_data: %Z% %M% %I% %E% %U%
*
* Mapping Agent - MapAg.c
*
* This module implements an agent for associating arbitrary data.  In all
* cases, the data consists of char* (pointers to anything).
*
* The hashing algorithm is tailored  based on the initial a,b,c passed in to 
* MapAg_Define_Initial.  The hashing algorithm may change at resize time
* if too many collisions are occurring.
*
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "MapAgP.h"

BEGIN_NOT_Cxx

/*  -- Private Method Declarations
*******************************************************************************
*/

static void  MapAg_Transition_InitialToNormal       _(( MapAg ));
static void  MapAg_Transition_NormalToInitial       _(( MapAg ));
static void  MapAg_Transition_InitialToNormalStatic _(( MapAg ));
static void  MapAg_Transition_NormalToInitialStatic _(( MapAg ));

static void  MapAg_Define_Initial  _(( MapAg, char*, char*, char*, char* ));
static void  MapAg_Define_Normal   _(( MapAg, char*, char*, char*, char* ));

static char* MapAg_Find_Initial    _(( MapAg, char*, char*, char* ));
static char* MapAg_Find_Normal     _(( MapAg, char*, char*, char* ));

static void  MapAg_Forget_Initial  _(( MapAg, char*, char*, char* ));
static void  MapAg_Forget_Normal   _(( MapAg, char*, char*, char* ));

static void  MapAg_Purge_Initial   _(( MapAg, char* ));
static void  MapAg_Purge_Normal    _(( MapAg, char* ));

static void  MapAg_Resize_Initial  _(( MapAg ));
static void  MapAg_Resize_Normal   _(( MapAg ));

static void  MapAg_Free_Initial       _((MapAg));
static void  MapAg_Free_Normal        _((MapAg));
static void  MapAg_Free_InitialStatic _((MapAg));
static void  MapAg_Free_NormalStatic  _((MapAg));


/* -- Statically Allocate and Initialize State Objects
*******************************************************************************
*/

static StateRec State_initialRec = {
    MapAg_Transition_InitialToNormal,
    MapAg_Define_Initial,
    MapAg_Find_Initial,
    MapAg_Forget_Initial,
    MapAg_Purge_Initial,
    MapAg_Resize_Initial,
    MapAg_Free_Initial
};

static State State_initial = &State_initialRec;

static StateRec State_normalRec = {
    MapAg_Transition_NormalToInitial,
    MapAg_Define_Normal,
    MapAg_Find_Normal,
    MapAg_Forget_Normal,
    MapAg_Purge_Normal,
    MapAg_Resize_Normal,
    MapAg_Free_Normal
};

static State State_normal = &State_normalRec;

/* extern for MapAg_STATIC macro
*/
StateRec State_initialStaticRec = {
    MapAg_Transition_InitialToNormalStatic,
    MapAg_Define_Initial,
    MapAg_Find_Initial,
    MapAg_Forget_Initial,
    MapAg_Purge_Initial,
    MapAg_Resize_Initial,
    MapAg_Free_InitialStatic
};

static State State_initialStatic = &State_initialStaticRec;

static StateRec State_normalStaticRec = {
    MapAg_Transition_NormalToInitialStatic,
    MapAg_Define_Normal,
    MapAg_Find_Normal,
    MapAg_Forget_Normal,
    MapAg_Purge_Normal,
    MapAg_Resize_Normal,
    MapAg_Free_NormalStatic
};

static State State_normalStatic = &State_normalStaticRec;


/* -- Forward Declarations of All Static Functions
==================================================
*/
static char* Malloc _((int bytes));
static char* Calloc _((int bytes, int count));
static void  Free   _((char*));

static void MapAg_FirstGuessAtHashAlgorithm _((MapAg,char*,char*,char*));
static void MapAg_RefineHashAlgorithm _((MapAg,int));
static int MapAg_ReHashUnderCollisionLimit _((MapAg,int));
static void MapAg_ReHash _((MapAg,int,int));
static int MapAg_InInitialState _(( MapAg ));
static int MapAg_InNormalState _(( MapAg ));
static void MapAg_FreeTable _((MapAg));


/* -- Heap Allocation and Free
==============================
   Error if heap allocation fails
*/
static char* Malloc( bytes )
    int bytes;
{
    char* allocated = (char*)malloc( (unsigned)bytes );
    if ( allocated == NULL )
    {
        fprintf( stderr, "MapAg could not allocate %d bytes.\n", bytes );
        abort();
        exit(1);
    }
    return allocated;
}

static char* Calloc( bytes, count )
    int bytes, count;
{
    char* allocated = (char*)calloc( (unsigned)bytes, (unsigned)count );

    if ( allocated == NULL )
    {
        fprintf( stderr, "MapAg could not allocate %d x %d bytes.\n",
                 bytes, count );
        abort();
        exit(1);
    }
    return allocated;
}

static void Free( allocated )
    char* allocated;
{
    if ( allocated != NULL )
        free( allocated );
}



/* -- Hash Method - Find Bucket
================================================
    The algorithm REQUIRES that hashSize be a power of 2!
    Mask is all ones: (some power of 2) - 1
    Shifts are determined based on initial a,b,c passed to MapAg_Define_Initial
*/
#define MapAg_HashBucket(ag,a,b,c) \
  ( ( ( (intptr_t)(a) >> (ag)->shiftA ) +  \
      ( (intptr_t)(b) >> (ag)->shiftB ) +  \
      ( (intptr_t)(c) >> (ag)->shiftC ) ) & (ag)->mask )

/* Determines shifts we use within the hashing algorithm above
*/
static void MapAg_FirstGuessAtHashAlgorithm( ag, a, b, c )
    MapAg       ag;
    char        *a, *b, *c;
{
    if ( a == (char*)0 )
    {
        ag->shiftA = 0;
    }
    else
    {
        /* Find the number of least significant 0 bits in a
        */
        int i = (intptr_t)a;
        for ( ag->shiftA = 1 ; ag->shiftA < 8 * (sizeof(char*)) ; ag->shiftA++ )
        {
            /* shift right some bits, then left again, see if we lose any
            */
            if ( i != (( i >> ag->shiftA ) << ag->shiftA ) )
            {
                ag->shiftA--;   /* cannot mask so many bits */
                break;          /* out of for loop */
            }
        }
    }

    if ( b == (char*)0 )
    {
        ag->shiftB = 0;
    }
    else
    {
        int i = (intptr_t)b;
        for ( ag->shiftB = 1 ; ag->shiftB < 8 * (sizeof(char*)) ; ag->shiftB++ )
        {
            if ( i != (( i >> ag->shiftB ) << ag->shiftB ) )
            {
                ag->shiftB--;
                break;
            }
        }
    }

    if ( c == (char*)0 )
    {
        ag->shiftC = 0;
    }
    else
    {
        int i = (intptr_t)c;
        for ( ag->shiftC = 1 ; ag->shiftC < 8 * (sizeof(char*)) ; ag->shiftC++ )
        {
            if ( i != (( i >> ag->shiftC ) << ag->shiftC ) )
            {
                ag->shiftC-- ;
                break;
            }
        }
    }
}

/* -- Refine the Hashing Algorithm
==================================
    Obviously, this is heuristic.  The basic assumption (or hope) is that
    the collisions are what we need to worry about.  The given hash bucket
    is the one which has too many collisions.  The first cut is to try and
    find a new algorithm which prevents these collisions.  We then do a
    re-layout using an algorithm suggested by these collisions.  If we
    never hit the collision limit during the re-layout, then we are done.
*/
 
static void MapAg_RefineHashAlgorithm( agent, bkt )
    MapAg agent;
    int   bkt;
{
    Map map;

    int oldShiftA = agent->shiftA;
    int oldShiftB = agent->shiftB;
    int oldShiftC = agent->shiftC;

    for ( map = agent->table[bkt] ; map ; map = map->next )
    {
        MapAg_FirstGuessAtHashAlgorithm( agent, map->a, map->b, map->c );

        if ( MapAg_ReHashUnderCollisionLimit( agent, bkt ) )
            return;
    }

    /* OK, a hashing algorithm based on any of the (a,b,c) in the
     * original collision set is WORSE than the original algorithm.
     * Not bloody likely, but OK, stranger things have happened.
     */

    /* Nothing improved the situation.  Punt, revert to the old behavior.
     */
    agent->shiftA = oldShiftA;
    agent->shiftB = oldShiftB;
    agent->shiftC = oldShiftC;
}

static int MapAg_ReHashUnderCollisionLimit( agent, bkt )
    MapAg agent;
    int   bkt;
{
    Map  map;
    int  numOldCollisions, numNewCollisions, size, *collisions;

    numOldCollisions = numNewCollisions = 0;
    for ( map = agent->table[bkt] ; map ; map = map->next )
    {
        ++numOldCollisions;
        if ( MapAg_HashBucket( agent, map->a, map->b, map->c ) == bkt )
            ++numNewCollisions;
    }

    if ( numNewCollisions == numOldCollisions )
        return 0;       /* this algorithm does not help */

    /* This new algorithm does help in this specific set of
    ** collisions.  See if it works for all current mappings:
    ** Look at every map, and compute its new hash location.  Count
    ** the number of new collisions at each hash location.  If any
    ** exceeds the limit, this algorithm is no good.
    */
    size = 1 + agent->mask;
    collisions = (int*)Calloc( size, sizeof(int) );
    for ( bkt = 0 ; bkt < size ; ++bkt )
    {
        for ( map = agent->table[bkt] ; map ; map = map->next )
        {
            int newBkt = MapAg_HashBucket( agent, map->a, map->b, map->c );
            collisions[newBkt] += 1;

            if ( collisions[newBkt] == agent->collisionLimit )
            {
                Free( (char*)collisions );
                return 0;       /* this algorithm does not help */
            }
        }
    }
    Free( (char*)collisions );

    /* This is an improved algorithm.  Go ahead and actually re-layout
    */
    MapAg_ReHash( agent, size, size );
    return 1;
}

static void MapAg_ReHash( agent, newSize, oldSize )
    MapAg agent;
    int   newSize, oldSize;
{
    int  oldBkt, newBkt;
    Map* newTable = (Map*)Calloc( newSize, sizeof(Map) );
    Map* oldTable = agent->table;

    agent->table = newTable;

    for ( oldBkt = 0 ; oldBkt < oldSize ; oldBkt++ )
    {
        /* Each non-empty bucket points to a Map, each Map is a link
        ** in a chain.  Need to go down chain and re-hash each Map.
        */
        Map map, next;

        for ( map = oldTable[oldBkt] ; (Map)0 != map ; map = next )
        {
            /* 'map' gets pointed to by bucket in new hash table, what
            ** was in that bucket becomes the next ptr of 'map'.
            */
            newBkt = MapAg_HashBucket( agent, map->a, map->b, map->c );
            next = map->next;
            map->next = agent->table[newBkt];
            agent->table[newBkt] = map;
        }
    }
    Free( (char*)oldTable );
}

/* -- MapAg Sanity Checks
==========================
    Each check is an individual if statement so its easier to see the problem
*/

static int MapAg_InInitialState( agent )
    MapAg       agent;
{
    if ( agent->state   != State_initial
      && agent->state   != State_initialStatic  )       return 0;
    if ( agent->table   != (Map*)0              )       return 0;
    if ( agent->mask    != 0                    )       return 0;
    if ( agent->numMaps != 0                    )       return 0;
    if ( agent->shiftA  != 0                    )       return 0;
    if ( agent->shiftB  != 0                    )       return 0;
    if ( agent->shiftC  != 0                    )       return 0;
    return 1;
}

static int MapAg_InNormalState( agent )
    MapAg       agent;
{
    if ( agent->state   != State_normal
      && agent->state   != State_normalStatic   )       return 0;
    if ( agent->table   == (Map*)0              )       return 0;
    if ( 0 != (agent->mask & (agent->mask + 1)) )       return 0;
    if ( agent->numMaps <=  0                   )       return 0;
    return 1;
}

extern void MapAg_AssertLooksOk( agent, file, line )
    MapAg       agent;
    char*       file;
    int         line;
{
    if (   !MapAg_InNormalState(  agent )
        && !MapAg_InInitialState( agent ) )
    {
        fprintf(stderr, "Assertion failed: file \"%s\", line %d\n", file, line);
        abort();
        exit(1);
    }
}

/* -- Transition Between States
===============================
   Initial state is always zero filled but with the state pointer set.
*/

static void MapAg_Transition_NormalToInitial( agent )
    MapAg agent;
{
    Free( (char*)agent->table );
    agent->table = (Map*)0;
    agent->mask  = 0;
    agent->shiftA = agent->shiftB = agent->shiftC = 0;
    agent->state = State_initial;
}

static void MapAg_Transition_NormalToInitialStatic( agent )
    MapAg agent;
{
    Free( (char*)agent->table );
    agent->table = (Map*)0;
    agent->mask  = 0;
    agent->shiftA = agent->shiftB = agent->shiftC = 0;
    agent->state = State_initialStatic;
}

static void MapAg_Transition_InitialToNormal( agent )
    MapAg agent;
{
    agent->state = State_normal;
}

static void MapAg_Transition_InitialToNormalStatic( agent )
    MapAg agent;
{
    agent->state = State_normalStatic;
}

/* -- MapAg "Find" Methods
==========================
*/

/*ARGSUSED*/
static char* MapAg_Find_Initial( agent, a, b, c )
    MapAg       agent;
    char        *a, *b, *c;
{
    return (char*)0;
}

static char* MapAg_Find_Normal( agent, a, b, c )
    MapAg       agent;
    char        *a, *b, *c;
{
    Map map = agent->table[ MapAg_HashBucket( agent, a, b, c ) ];

    for ( ; (Map)0 != map ; map = map->next )
        if ( (map->a == a) && (map->b == b) && (map->c == c) )
            return map->data;

    return (char*)0;
}


/*  -- Make new definition of a,b,c to d
========================================================
    Might need to grow hash table.  
    New map is head of chain if hash collision.
*/

static void MapAg_Define_Initial( agent, a, b, c, data )
    MapAg agent;
    char  *a, *b, *c, *data;
{
    Map newMap = (Map)Malloc( sizeof(MapRec) );  /* careful! must initialize! */

    /* Do initial allocation of table
    */
    MapAg_Resize(agent);

    newMap->a = a; newMap->b = b; newMap->c = c; newMap->data = data;
    newMap->next = (Map)0;

    /* Change behaviors to normal behaviors
    */
    MapAg_Transition( agent );

    /* Determine a hashing function based on this first mapping
    */
    MapAg_FirstGuessAtHashAlgorithm( agent, a, b, c );

    agent->table[ MapAg_HashBucket( agent, a, b, c ) ] = newMap;

    agent->numMaps = 1;
}

static void MapAg_Define_Normal( agent, a, b, c, data )
    MapAg agent;
    char  *a, *b, *c, *data;
{
    Map map;
    int bkt;
    int numCollisions = 0;

    bkt = MapAg_HashBucket( agent, a, b, c );

    for ( map = agent->table[ bkt ] ; (Map)0 != map ; map = map->next )
    {
        if ( (map->a == a) && (map->b == b) && (map->c == c) )
        {
            /* Replace existing mapping
            */
            map->data = data;
            return;
        }
        ++numCollisions;
    }

    /* No map at this bucket, or no map at this bucket with same a,b,c
    ** If hash table is too full, resize it.   If too many collisions, try a
    ** new hashing function.  In either of these cases, we need to
    ** re-compute the hash bucket.
    */
    if ( agent->mask <= agent->numMaps )
    {
        MapAg_Resize(agent);
        bkt = MapAg_HashBucket( agent, a, b, c );
    }
    else if ( agent->collisionLimit <= numCollisions )
    {
        MapAg_RefineHashAlgorithm( agent, bkt );
        bkt = MapAg_HashBucket( agent, a, b, c );
    }

    /* Make a new map, initialize it, put it at head of chain from its bucket.
    */
    map = (Map)Malloc( sizeof(MapRec) );
    map->a    = a;
    map->b    = b;
    map->c    = c;
    map->data = data;
    map->next = agent->table[bkt];
    agent->table[bkt] = map;
    agent->numMaps++;
}


/*  -- MapAg Resize Methods
===============================================================================
    When nothing is in the MapAg, MapAg_Resize_Initial is the Resize method
    of the MapAg.  This initial method allocates the hash table.  The table
    size must ALWAYS be a power of two, as we use a mask in the hashing
    function which must be all ones.

    Its a good idea if the initial table size is small, so there are many
    chances early on to refine the hashing algorithm.

    MapAg_Resize_Normal doubles the size of the hash table and re-distributes
    the entries.  It also doubles the size of the collision limit.
*/

#define MAPAG_InitTableSize 8
#define MAPAG_InitCollisionLimit 1

static void MapAg_Resize_Initial(agent)
    MapAg agent;
{
    agent->table = (Map*)Calloc( MAPAG_InitTableSize, sizeof(Map) );

    agent->mask  = MAPAG_InitTableSize-1;
    agent->collisionLimit = MAPAG_InitCollisionLimit;
}

/* Note we DO NOT need to consider collisions during re-sizing: the resized
 * (bigger) table will certainly have the same or less collisions as the old
 * table.  If we did not hit the collision limit before, we will not during the
 * re-hashing.  This is certain, because the bucket calculation will distribute
 * the same f(a,b,c) across more hash table slots.
 */
static void MapAg_Resize_Normal(agent)
    MapAg       agent;
{
    int         oldSize  = agent->mask + 1;
    int         newSize  = 2*oldSize;

    agent->mask = newSize-1;            /* size is now doubled ...      */
    agent->collisionLimit *= 2;         /* so collision limit doubles   */

    MapAg_ReHash( agent, newSize, oldSize );
}

/*  -- Drop Map from MapAg
==========================
*/

/*ARGSUSED*/
static void MapAg_Forget_Initial( agent, a, b, c )
    MapAg agent;
    char  *a, *b, *c;
{
    return;
}

static void MapAg_Forget_Normal( agent, a, b, c )
    MapAg agent;
    char  *a, *b, *c;
{
    int bkt  = MapAg_HashBucket( agent, a, b, c );
    Map map  = agent->table[bkt];
    Map prev = (Map)0;

    for ( ; map ; prev = map, map = map->next )
    {
        if ( (map->a == a) && (map->b == b) && (map->c == c) )
        {
            if ((Map)0 == prev)
                agent->table[bkt] = map->next;
            else
                prev->next = map->next;
            Free( (char*)map );
            if ( --(agent->numMaps) == 0 )
            {
                MapAg_Transition( agent );
            }
            return;
        }
    }
}


/*  -- MapAg Purge Data Methods
===============================================================================
    Purge all mappings to this data from the hash table.  We simply do a
    complete search of the hash table, and remove all Map objects which point
    to this data item.
*/

/*ARGSUSED*/
static void MapAg_Purge_Initial( agent, data )
    MapAg agent;
    char* data;
{
    return;
}

static void MapAg_Purge_Normal( agent, data )
    MapAg agent;
    char* data;
{
    int bkt;

    for ( bkt = 0 ; bkt <= agent->mask ; bkt++ )
    {
        Map prev, map, next;

        prev = (Map)0;
        for ( map = agent->table[bkt] ; map ; map = next )
        {
            next = map->next;

            if ( map->data == data )
            {
                /* Get rid of this Map
                */
                if ( (Map)0 == prev )
                {
                    agent->table[bkt] = next;
                }
                else
                {
                    prev->next = next;
                }

                Free( (char*)map );

                if ( --(agent->numMaps) == 0 )
                {
                    /* Just deleted the last Map for this Agent.
                    */
                    MapAg_Transition( agent );
                    return;
                }
            }
            else
            {
                prev = map;
            }
        }
    }
}


/*  -- Private: Free all maps, free the hash table
--------------------------------------------------
*/

static void MapAg_FreeTable( agent )
    MapAg agent;
{
    int bkt;

    for ( bkt = 0 ; bkt <= agent->mask ; bkt++ )
    {
        Map map, next;
        for ( map = agent->table[bkt] ; map ; map = next )
        {
            next = map->next;
            Free( (char*)map );
        }
    }
    Free( (char*)agent->table );
}

static void MapAg_Free_Initial( agent )
    MapAg agent;
{
    Free( (char*)agent );
}

static void MapAg_Free_Normal( agent )
    MapAg agent;
{
    MapAg_FreeTable( agent );
    Free( (char*)agent );
}

/*ARGSUSED*/
static void MapAg_Free_InitialStatic( agent )
    MapAg agent;
{
    return;
}

static void MapAg_Free_NormalStatic( agent )
    MapAg agent;
{
    /* We cannot free the MapAg itself, since it is statically allocated
     * Therefore, just re-initialize agent.
     */
    MapAg_Transition( agent );
}


/*  -- Mapping Agent Public Functions
*******************************************************************************
*/

/*  -- Create a new Mapping Agent
=======================================================================
*/
extern MapAg MapAg_New()
{
    MapAg agent = (MapAg)Calloc( sizeof(MapAgRec), 1 );

    agent->state = State_initial;

    MapAg_AssertLooksOk( agent, __FILE__, __LINE__ );
    return agent;
}

/*  -- Invoke state specific behaviors
======================================
*/
extern void MapAg_DoDefine( agent, a, b, c, data )
    MapAg agent;
    char  *a, *b, *c, *data;
{
    agent->state->Define( agent, a, b, c, data );
}

extern char* MapAg_DoFind( agent, a, b, c )
    MapAg agent;
    char  *a, *b, *c;
{
    return agent->state->Find( agent, a, b, c );
}

extern void  MapAg_DoForget( agent, a, b, c )
    MapAg agent;
    char  *a, *b, *c;
{
    agent->state->Forget( agent, a, b, c );
}

extern void MapAg_DoPurge( agent, data )
    MapAg agent;
    char* data;
{
    agent->state->Purge( agent, data );
}

extern void MapAg_DoFree( agent )
    MapAg agent;
{
    agent->state->Free( agent );
}

END_NOT_Cxx
