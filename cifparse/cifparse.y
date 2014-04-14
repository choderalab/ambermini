/*

COPYRIGHT 1996 Rutgers - The State University of New Jersey

This software is provided WITHOUT WARRANTY OF MERCHANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED.  RUTGERS MAKE NO REPRESENTATION OR WARRANTY THAT THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT OR OTHER
PROPRIETARY RIGHT.

The user of this software shall indemnify, hold harmless and defend
Rutgers, its governors, trustees, officers, employees, students,
agents and the authors against any and all claims, suits,
losses, liabilities, damages, costs, fees, and expenses including
reasonable attorneys' fees resulting from or arising out of the
use of this software.  This indemnification shall include, but is
not limited to, any and all claims alleging products liability.

This software may be used only for not-for-profit educational
and research purposes.

*/

/* 
  FILE:         ndb_cifparse_parser.y

  This file is part of the NDBQUERY application,
  a program component of the Nucleic Acid Database Project.

  H. M. Berman, W. K. Olson, D. L. Beveridge, J. K. Westbrook, A. Gelbin,
  T. Demeny, S. H. Shieh, A. R. Srinivasan, and B. Schneider. (1992).
  The Nucleic Acid Database: A Comprehensive Relational Database of 
  Three-Dimensional Structures of Nucleic Acids.  Biophys J., 63, 751-
  759.

  Questions about this software should be directed to:

              ndbadmin@ndbserver.rutgers.edu

  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

/*$$LICENSE-NDB$$*/

%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cifparse.h"
int yylex();
void yyerror();

static int curItemNo, curValueNo, itemIndex;
static int  *fieldList = NULL;

%}
%union {
	char TempBuffer[MAXVALUELENGTH+1];
}
%token <TempBuffer> ITEMNAME
%token <TempBuffer> VALUE
%token <TempBuffer> LOOP
%token <TempBuffer> DATABLOCK
%token <TempBuffer> UNKNOWN
%token <TempBuffer> MISSING
%type  <TempBuffer> ItemName Value 
%%

Datablocks :     Lines
             |   Datablocks  Datablock
             ;
Datablock :    DatablockName Lines 
	  ;

Lines:  /* empty */
  | Lines Line
  ;

Line:       ItemName Value  /* instance of a non-looped item name value pair */
{
  ndb_cif_process_item_name_value_pair();
}
  |   ItemNameList ValueList
  ;


ItemNameList:  

  Loop  ItemName           /*  the beginning of a loop_ */
{
  ndb_cif_process_loop_declaration();
}
  
  | ItemNameList ItemName 
{
  ndb_cif_process_item_name_list();
}
  ;


ValueList:  Value 
{
  ndb_cif_rewind_column();
  ndb_cif_process_value_list();
}
  | ValueList Value
{
  ndb_cif_process_value_list();
}
;

ItemName : ITEMNAME  
{
  /*  sprintf(TempKeyword,"%s",$1); */
  strncpy(TempKeyword,$1,MxNameLen);
}
;

Loop :    LOOP
{
  curItemNo = 0;  curValueNo = 0;  itemIndex = 0;
}
;

Value:   VALUE 
{
  /*  sprintf(TempValue,"%s",$1); */
  strncpy(TempValue,$1,MAXVALUELENGTH);
}
        | UNKNOWN 
{
  strcpy(TempValue,"");
}
        | MISSING
{
  strcpy(TempValue,"");
}
;

DatablockName: DATABLOCK
{
  int idatablock;
  idatablock = ndb_cif_get_datablock_id(&($1)[5]);
  if (idatablock == 0)
    ndb_cif_new_datablock(&($1)[5]);
  else {
    ndb_cif_move_datablock(&($1)[5]);
    ndb_cif_reset_datablock();
  }
}
;

%%

void yyerror(s)
char *s;
/* ---------------------------------------------------------------------------
 * Purpose:  yyerror()
 * 
             Report the location of a parsing error...
 -----------------------------------------------------------------------------*/
 
{
  fprintf( stderr,"%s near line %d\n", s, lineNo);
}

void ndb_cif_process_item_name_list()
/* ----------------------------------------------------------------------
   Purpose: ndb_cif_process_item_name_list()

            Registers the item keyword for the the current item in the 
	    current category.  Maintains an index array of "valid" keyword 
	    names in fieldList[].  This array is used as an indirectly 
	    reference between keywords and values ...  
 * ----------------------------------------------------------------------*/
{
  int prevCategoryId, categoryId;
  char itemKeyword[MxNameLen], categoryName[MxNameLen];

/*
  fprintf( stderr,"Processing item name list at line %d keyword %s\n", lineNo, TempKeyword);
*/
  prevCategoryId = ndb_cif_current_category();
  categoryId = ndb_cif_get_category_name_from_item_name(categoryName, TempKeyword);


  fieldList = (int *) realloc(fieldList, (curItemNo+1) * sizeof(int));
  if (categoryId != 0 && categoryId == prevCategoryId) {
    ndb_cif_get_item_keyword_from_item_name(itemKeyword, TempKeyword);
    ndb_cif_put_item_keyword(itemKeyword);
    itemIndex++;
    fieldList[curItemNo] = itemIndex;
  }
  else {
    fprintf( stderr,"Syntax error line %d with %s\n", lineNo, TempKeyword);
    fieldList[curItemNo] = -1;
  }
  curItemNo++;

}

void ndb_cif_process_value_list()
/* ----------------------------------------------------------------------
     Purpose:  ndb_cif_process_value_list()

               Add the current value to the appropriate column in the 
               the current row.  Start a new row if necessary.
 * ----------------------------------------------------------------------*/
{
  int fieldNo;

/*
  fprintf( stderr,"Processing value list at line %d value %s\n", lineNo, TempValue);
*/
  if (fieldList[curValueNo] != -1) {
    fieldNo = fieldList[curValueNo];
    if (fieldNo == 1) ndb_cif_new_row();
    if (fieldNo != 0) ndb_cif_put_item_value(fieldList[curValueNo], TempValue);
  }
  curValueNo++;

  if (curValueNo == curItemNo) curValueNo = 0;
}

void ndb_cif_process_item_name_value_pair()
/* ----------------------------------------------------------------------
      Purpose: ndb_cif_process_item_name_value_pair()

               Assign the current value to its associated item name.
 * ----------------------------------------------------------------------*/
{
  int categoryId, tmpCategoryId;
  char categoryName[MxNameLen], itemKeyword[MxNameLen];

/*
  fprintf( stderr,"Processing item name value pair at line %d value %s\n", lineNo, TempKeyword);
*/

  tmpCategoryId = ndb_cif_current_category();
  categoryId = ndb_cif_get_category_name_from_item_name(categoryName,TempKeyword);
  curItemNo  = 1;
  curValueNo = 0;
  itemIndex  = 1;
  if (categoryId == 0) {
    if (strcmp(categoryName,"")) {
      ndb_cif_new_category(categoryName);
    }
    else {
      fprintf( stderr, "Missing category name component in %s at line %d\n", 
	     TempKeyword, lineNo);
      return ;
    }
  }
  else if (tmpCategoryId != categoryId) {
    fprintf( stderr,"Category conflict between %s and %s at line %d\n",  
		    cifFiles.datablocks[0].categories[tmpCategoryId].categoryName,
		    categoryName,
		    lineNo);
    return;
  }

  ndb_cif_get_item_keyword_from_item_name(itemKeyword, TempKeyword);
  ndb_cif_put_item_keyword(itemKeyword);
  fieldList = (int *) realloc(fieldList, sizeof(int));
  fieldList[0] = ndb_cif_current_col();
  ndb_cif_process_value_list();
}


void ndb_cif_process_loop_declaration()
/* ----------------------------------------------------------------------
     Purpose: ndb_cif_process_loop_declaration()

              Handles initialization for a new loop, by creating a new 
              category and adding the current item name to this category.
 * ---------------------------------------------------------------------- */
{
  char  categoryName[MxNameLen];
  int   categoryId;

/*
  fprintf( stderr,"Processing loop declaration at line %d value %s\n", lineNo, TempKeyword);
*/

  categoryId = ndb_cif_get_category_name_from_item_name(categoryName, TempKeyword);
  if (categoryId == 0) {
    ndb_cif_new_category(categoryName);
    ndb_cif_process_item_name_list();
  }
  else {
    fprintf( stderr,"Duplicate category name %s at line %d\n",categoryName, lineNo);
  }
}

#if 0
void yyprint(FILE *file, int type, YYSTYPE value)
{
  fprintf(file,"\nType  = %d\n",type);
  fprintf(file,"Value = %s\n",value.TempBuffer);
}
#endif

