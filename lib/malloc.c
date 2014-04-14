#include <stdio.h>
/*
 *  malloc.c - cmalloc() and crealloc() wrappers for Fortran
 *
 * 	Copyright: (C) OXFORD MOLECULAR LTD, 1995
 */

/* on SGI, MACHINE's LOADLIB should use -lmalloc for reliability */

/* FORTRAN/C interface-specific naming */
#ifdef CLINK_PLAIN
# define cmalloc_ cmalloc
# define crealloc_ crealloc
#endif
#ifdef CLINK_CAPS
# define cmalloc_ CMALLOC
# define crealloc_ CREALLOC
#endif
#ifdef SPARC
/* sparc is weird */
#define _cmalloc_ cmalloc
#define _crealloc_ crealloc
#endif

#include <sys/types.h>
#include <stdlib.h>
#include <malloc.h>

/*****************************************************************************
 *      Name: cmalloc_()
 *  Function: FORTRAN wrapper utility to the C malloc function. A NULL (0)   *
 *            value will result if the request could not be satisfied.       *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1995                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 12/01/95                                                       *
 *---------------------------------------------------------------------------*
 *    Inputs: int *NBytes  - The number of bytes required to be allocated    *
 *   Outputs: NONE                                                           *
 *   Returns: void *ArrPtr - Pointer to the block of memory just allocated   *
 * Externals: NONE                                                           *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/


void *
#ifdef __STDC__
cmalloc_(int *NBytes)
#else
cmalloc_(NBytes)
int	*NBytes;
#endif
{   
	void	*ArrPtr; /* Pointer to the block of allocated memory */

	/* 
	 *  Check that the number of bytes requested is sensible 
	 */
	if (*NBytes > 0) {
		ArrPtr = malloc((size_t)(*NBytes));
	} else {
		ArrPtr = NULL;
	}
	/*
	if (ArrPtr) fprintf(stderr, "mallocd %d\n", *NBytes);
	else fprintf(stderr, "didn't malloc %d\n", *NBytes);
	*/
      
	return(ArrPtr);
}
/*****************************************************************************
 *      Name: crealloc_()
 *  Function: FORTRAN wrapper utility to the C crealloc function. A NULL (0)  *
 *            value will result if the request could not be satisfied.       *
 * Copyright: (C) OXFORD MOLECULAR LTD, 1995                                 *
 *---------------------------------------------------------------------------*
 *    Author: Robert Williams                                                *
 *      Date: 12/01/95                                                       *
 *---------------------------------------------------------------------------*
 *    Inputs: void *ArrPtr - Pointer to the block of memory to creallocate    *
 *            int *NBytes  - The number of bytes required to be allocated    *
 *   Outputs: NONE                                                           *
 *   Returns: void *ArrPtr - Pointer to the block of memory just allocated   *
 * Externals: NONE                                                           *
 *---------------------------------------------------------------------------*
 * MODIFICATION RECORD:                                                      *
 * DD/MM/YY   Initials   Comments                                            *
 *****************************************************************************/


void *
#ifdef __STDC__
crealloc_(void **ArrPtr, int *NBytes)
#else
crealloc_(ArrPtr, NBytes)
void	**ArrPtr;
int	*NBytes;
#endif
{   
	void *NewPtr; /* Pointer to new memory */

	/* 
	 *  Check that the number of bytes requested is sensible 
	 */
	if (*NBytes > 0) {
		NewPtr = realloc(*ArrPtr, (size_t)(*NBytes));
	} else {
		NewPtr = NULL;
	}
      
	return(NewPtr);
}

