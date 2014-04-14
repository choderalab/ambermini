/*
 *      File: varArray.c
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
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              A VARARRAY is an array that can grow in size.
 *              It is represented internally as a single block of
 *              memory that is REALLOCed periodically to allow growth.
 *              
 *              The LINKEDLIST allows the caller to:
 *                      create new VARARRAY 
 *                      destroy VARARRAYs 
 *                      define the size of elements 
 *                      obtain pointers to the start of indexed elements
 *                      add objects to the . 
 *
 *
 *       Body of varArray.[ch] modified by Vladimir Romanovski (1994).  
 *
 *
 ************************************************************************/



#include	"basics.h"
#include        "varArray.h"




#define element(a,i) ((char*)((a)->data + (i)*(a)->size))
#define PORTION  10
#define SLOT_FOR_COUNT(count) \
  ( count > 0 ? (((count + PORTION -1) / PORTION ) * PORTION) : PORTION)


/*-----------------------------------------------------
 *               iVarArrayPointerToIndex
 *
 *
 */
int
iVarArrayPointerToIndex( VARARRAY header, char *data)
{
  if ( header == NULL || data == NULL){
    DFATAL(( " iVarArrayPointerToIndex: VARARRAY or Data is NULL" ));
  }
  return( (data - header->data) / header->size );
}

/*-----------------------------------------------------
 *                iVarArrayElementSize
 *
 *
 */
int
iVarArrayElementSize( VARARRAY header)
{
  if ( header == NULL){
    DFATAL(( " iVarArrayElementSize: VARARRAY is NULL" ));
  }
  return( header->size );
}

/*------------------------------------------------------
 *                iVarArrayElementCount
 *
 *
 */
int
iVarArrayElementCount( VARARRAY header)
{
  if ( header == NULL)
	return( 0 );
  return( header->count );
}

/*-----------------------------------------------------
 *               PVarArrayIndex
 */
char *
PVarArrayIndex( VARARRAY header, int pos)
{
  if ( header == NULL){
    DFATAL(( " PVarArrayIndex: VARARRAY is NULL" ));
  }
  return( element( header,pos));

}

/*-----------------------------------------------------
 *      vaVarArrayCreate
 *
 *
 *      Create a new VARARRAY and initialize it.
 *      The caller must initialize the size of the elements
 *      of the VARARRAY when they create it.
 */

VARARRAY 
vaVarArrayCreate( int size )
{
  VARARRAY new;

  MALLOC(new, VARARRAY, sizeof(HeaderStruct));
  
  if (new == NULL) {
    DFATAL(( "vaVarArrayCreate: not enough memory " ));
  }
  MALLOC(new->data, char *, size *PORTION );

  if (new->data == NULL) {
    DFATAL(( "vaVarArrayCreate: not enough memory " ));
  }
  new->count = 0;
  new->size  = size;
  new->slot  = PORTION;
  
  return( new );
}

/*-----------------------------------------------------
 *      VarArrayDestroy
 *
 *      Destroy the VarArray, after this call the pointer will
 *      be undefined.
 */


void
VarArrayDestroy( VARARRAY *header )
{

  if ( *header == NULL){
    DFATAL(( " VarArrayDestroy: VARARRAY is NULL" ));
  }
  if ( (*header)->data != NULL ) FREE( (*header)->data );
  
  FREE( *header );
  *header = NULL;
}



/*-----------------------------------------------------
 *      VarArrayAdd
 *
 *      Add one element to the VARARRAY and copy the data into it.
 *      This will require REALLOCing the array and probably changing
 *      its address.
 *
 */
void
VarArrayAdd( VARARRAY header, GENP data )
{
  if ( header == NULL) {
    DFATAL(( " VarArrayAdd: VARARRAY is NULL" ));
  }
  if ( header->count == header->slot ) {
    header->slot += PORTION;
    REALLOC(header->data , char*, header->data, header->size * header->slot );
  }

  memcpy(element(header,header->count), (char*)data, header->size);

  header->count++;
}


/*-----------------------------------------------------
 *      vaVarArrayCopy
 *	Copy the VARARRAY.
 *
 *
 */
VARARRAY 
vaVarArrayCopy( VARARRAY header )
{
  VARARRAY	new;

  if ( header == NULL) {
    DFATAL(( " VarArrayAdd: VARARRAY is NULL" ));
  }

  MALLOC(new, VARARRAY, sizeof(HeaderStruct) );

  new->size  = header->size;
  new->count = header->count;
  new->slot  = header->slot;  

  MALLOC(new->data, char*, new->size * header->slot);

  memcpy(new->data, header->data, new->size * header->slot);
  
  return( new );
}

/*-----------------------------------------------------
 *      vaVarArrayCopy2
 *	Copy 2 VARARRAYs into one.
 *
 *
 */
VARARRAY
vaVarArrayCopy2( VARARRAY header1, VARARRAY header2 )
{
  VARARRAY	new;
  int	  copysize;

  if ( header1 == NULL  ||  header2 == NULL ) {
    DFATAL(( " VarArrayCopy2: VARARRAY is NULL" ));
  }
  if ( header1->size != header2->size )
    DFATAL(( " VarArrayCopy2: header sizes different\n" ));

  MALLOC(new, VARARRAY, sizeof(HeaderStruct) );

  new->size  = header1->size;
  new->count = header1->count + header2->count;
  new->slot  = SLOT_FOR_COUNT(new->count);

  MALLOC(new->data, char*, new->size * new->slot);

  copysize = new->size * header1->count;
  memcpy(new->data, header1->data, copysize);
  memcpy(new->data+copysize, header2->data, new->size * header2->count);
  
  return( new );

}

  
void
VarArrayInsertBeforeMore( VARARRAY header, int pos, int num )
{
	int	shift, nslot;
	char	*h;

  	if ( header == NULL)
		DFATAL(( " VarArrayInsertBeforeMore: VARARRAY is NULL" ));

	if ( (pos >= header->count) && (pos < 0 ) )
		DFATAL(( " VarArrayInsertBeforeMore: position=%d",pos ));

	/*
	 *  grow array if necc
	 */
	nslot = SLOT_FOR_COUNT(header->count + num);
  
	if ( header->slot != nslot ) {
	    REALLOC(header->data , char*, header->data, header->size * nslot );
	    header->slot = nslot;
	}

	/*
	 *  update item count
	 */
	header->count += num;

	/*
	 *  open up insert space by shuffling remainder down
	 */
	shift = header->size * num;
	h = element(header, pos);
 	memmove( h+shift, h, (header->count - num - pos) * header->size );
}


/*-----------------------------------------------------
 *      VarArrayInsertBefore
 *
 *	Add one element to the VARARRAY and move all of the data
 *	at index iPos and beyond up one element.
 *	Copy the data at data into the new element that
 *	has been opened up.
 *
 */
void 
VarArrayInsertBefore( VARARRAY header, int pos, GENP data )
{
  VarArrayInsertBeforeMore(header, pos, 1);
  memcpy(element(header,pos), (char*)data, header->size);
}



/*-----------------------------------------------------
 *	VarArrayDelete
 *
 *	Remove an element in a VARARRAY.
 *	Move the data below the one to be removed, up one
 *
 */
void
VarArrayDeleteMore( VARARRAY header, int pos, int num)
{
  int      shift, nslot;
  register char* i, *h;
  
  if ( header == NULL) {
    DFATAL(( " VarArrayDelete: VARARRAY is NULL" ));
  }
  if ( ((pos+num) > header->count) || (pos < 0 ) || (num < 1)){
    DFATAL(( " VarArrayDelete: position=%-5d num=%-5d count=%-5d",
	    pos,num,header->count ));
  }
  header->count -= num;
  
  shift = num * header->size;

  for (i = element( header, pos), h = element( header, header->count);
       i < h  ;
       i++)
    *i = *(i + shift);

  nslot = SLOT_FOR_COUNT( header->count );

  if ( header->slot != nslot ){
    REALLOC(header->data , char*, header->data, header->size * nslot );
    header->slot = nslot;
  }
}


/*-----------------------------------------------------
 *      VarArraySetSize
 *
 *      Change the size of the array in terms of elements.
 *      The size of the array will be adjusted so that it
 *      can contain iElements elements.
 *      All previous contents of the VARARRAY are still there,
 *      unless the VARARRAY was made smaller, then the tail   is lost.
 */
void
VarArraySetSize( VARARRAY header, int ncount )
{
  int nslot;

  if ( header == NULL){
    DFATAL(( " VarArraySetSize: VARARRAY is NULL" ));
  }
  
  if ( ncount < 0  ) {
    DFATAL(( " VarArraySetSize: elements=%5d",ncount ));
  }
  nslot = SLOT_FOR_COUNT(ncount);
  header->count = ncount;

  if ( header->slot != nslot ) {
    REALLOC(header->data , char*, header->data, header->size * nslot );
    header->slot  = nslot;
  }
}


/*-----------------------------------------------------
 *      PVarArrayDebugIndex
 *
 *      Return a pointer to the element within the VARARRAY, but
 *      first check the bounds.  Report an  if there
 *      is an out of bound access.
 */

GENP
PVarArrayDebugIndex( VARARRAY header, int pos, char *file, int line )
{
  if ( header == NULL ) {
    DFATAL(( "Attempting to access an invalid VARARRAY (%s line %d).",
	file, line ));
  }
  if ( header->count == 0 ) {
    DFATAL(( "Attempting to access a no-data VARARRAY (%s line %d).",
	file, line ));
  }
  if ( pos < 0 || pos >= header->count ) {
    DFATAL(( "Attempted to access element: %d in a VARARRAY of size: %d",
	    pos, header->count ));
  }
  return( element( header, pos) );
}

