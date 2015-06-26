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

Modified by D.A. Case, 12/97: remove unused code, rename the files,
   allow the lexer and parser to co-exist in a program with other lexers
   and parsers.

*/

/* 

  This file is part of the NDBQUERY application,
  a program component of the Nucleic Acid Database Project.

  H. M. Berman, W. K. Olson, D. L. Beveridge, J. K. Westbrook, A. Gelbin,
  T. Demeny, S. H. Shieh, A. R. Srinivasan, and B. Schneider. (1992).
  The Nucleic Acid Database: A Comprehensive Relational Database of 
  Three-Dimensional Structures of Nucleic Acids.  Biophys J., 63, 751-
  759.

  Questions about this software should be directed to:

              ndbadmin@ndbserver.rutgers.edu

  PURPOSE:       This collection of routines provides a simple interface 
                 to the DDL 2.1 compliant CIF files.

  Modified for incorporation into NAB by D.A. Case, 12/97.

*/

/*$$LICENSE-NDB$$*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define CIF_GLOBAL
#include "cifparse.h"

static char null_char[2] = "?"; /* The character output for null data values */

void ndb_cif_print_item_name(FILE *fp, char *itemName, int *linePos);
void ndb_cif_print_item_value(FILE *fp, char *itemValue, int *linePos);
void ndb_cif_print_formatted_item_value(FILE *outfile, int *offset, char *value, int newline_end, int fieldLength);

char ndb_cif_set_null_char(new_null)
char new_null;
     /* Change the CIF character identifying a NULL data value .
      */
{
  char old_null;
  old_null = null_char[0];
  null_char[0] = new_null;
  return old_null;
}

int ndb_cif_read_file(fp)
FILE *fp;
{


  if (fp == NULL ) {
    return 0;
  }
  cifpin = fp;
  lineNo = 0;
  rewind(fp);
  if (cifpparse() != 0) { /* YYACCEPT */
    ndb_cif_rewind_datablock();
    return FALSE;
  }
  return cifFiles.curDatablock+1;

}

int ndb_cif_write_file(fp)
FILE *fp;
{
  int i, datablockId, categoryId, rowId, linePos, numColumn, numRow;
  char itemValue[YYLMAX], itemName[MxNameLen], datablockName[MxNameLen];


  /* Open input and output files */
  if (fp == NULL ) {
    return 0;
  }

  datablockId = ndb_cif_rewind_datablock();
  do {
    ndb_cif_current_datablock_name(datablockName);
    fprintf( stderr, "Writing data block %s[%d of %d]\n",
	   datablockName,
	   ndb_cif_current_datablock(),
	   ndb_cif_count_datablock());
    fprintf(fp, "data_%s\n", datablockName);
    do {
      categoryId = ndb_cif_current_category();
      numColumn = ndb_cif_count_column();
      numRow = ndb_cif_count_row();
      fprintf(fp, "#\n");
      if (numRow <=1) {
	for (i=0; i< numColumn; i++) {
	  ndb_cif_get_item_name(i+1, itemName);
	  ndb_cif_print_item_name(fp, itemName, &linePos);
	  ndb_cif_get_item_value(i+1, itemValue, YYLMAX);
	  ndb_cif_print_item_value(fp, itemValue, &linePos);
	  if (linePos != 0) fprintf(fp, "\n");
	}
      }
      else {
	fprintf(fp, "loop_\n");
	for (i=0; i< numColumn; i++) {
	  ndb_cif_get_item_name(i+1, itemName);
	  ndb_cif_print_item_name(fp, itemName, &linePos);
	  fprintf(fp, "\n");
	}
	do {
	  linePos = 0;
	  for (i=0; i< numColumn; i++) {
	    ndb_cif_get_item_value(i+1, itemValue, YYLMAX);
	    ndb_cif_print_item_value(fp, itemValue, &linePos);
	  }
	  if (linePos != 0) fprintf(fp, "\n");
	  rowId = ndb_cif_next_row();
	}  while (rowId);
      }
      fflush(fp);
      categoryId = ndb_cif_next_category();
    } while (categoryId);
    fprintf(fp, "#\n");
    fflush(fp);
    datablockId = ndb_cif_next_datablock();
  } while (datablockId);
  ndb_cif_rewind_datablock();
  return(1);
}


void ndb_cif_print_item_name(fp, itemName, linePos)
FILE *fp;
char *itemName;
int *linePos;
{
  fprintf(fp, "%s  ", itemName);
  *linePos = strlen(itemName) +2;
}


void ndb_cif_print_item_value(fp, itemValue, linePos) 
FILE *fp;
char *itemValue;
int *linePos;
{
  int i, str_len, multipleLine, multipleWord, embeddedQuotes, len;

  /* NULL VALUE */
  if (itemValue == NULL || !strcmp(itemValue,"")) { 
    if (*linePos < MxNameLen - 1) {
      (*linePos) +=1;
      fputs(null_char, fp);
      if (*linePos < MxNameLen -2) {
	fprintf(fp, "  ");
	(*linePos) +=2;
      }
      else {
	fprintf(fp, "\n");
	(*linePos) = 0;
      }
    }
    else {
      fprintf(fp, "\n?  ");
      (*linePos) = 3;
    }
    return;
  }

  str_len = strlen(itemValue);
  multipleLine = FALSE;
  multipleWord = FALSE;
  embeddedQuotes = FALSE;
  len = 0;
  for (i=0; i<str_len; i++) {
    if (itemValue[i] == ' ') multipleWord = TRUE;
    else if (itemValue[i] == '\n') multipleLine =TRUE ;
    else if (itemValue[i] == '\'' || itemValue[i] == '\"') embeddedQuotes =TRUE ;
    if (itemValue[i] == '\t') len +=8;
    else                      len ++;
  }
  if (embeddedQuotes && multipleWord) multipleLine = TRUE;
  if (embeddedQuotes) multipleWord = TRUE;

  if (len > MxNameLen || multipleLine) {
    if (*linePos != 0)
      fprintf(fp, "\n");
    fprintf(fp, ";%s", itemValue);
    if (itemValue[str_len-1] != '\n')
      fprintf(fp, "\n");
    fprintf(fp, ";\n");
    *linePos = 0;
  }
  else {
    if ((!multipleWord && len + (*linePos) > MxNameLen) ||
	(multipleWord && len + 2 + (*linePos) > MxNameLen)) 
      fprintf(fp, "\n");
    if (multipleWord) {
      fprintf(fp, "\'%s\'", itemValue);
      (*linePos) += (len + 2);
    }
    else {
      fprintf(fp, "%s", itemValue);
      (*linePos) += len ;
    }
    if (*linePos +2 < MxNameLen) {
      fprintf(fp, "  ");
      (*linePos) +=2;
    }
    else {
      fprintf(fp, "\n");
      (*linePos) = 0;
    }
  }
	   
}


void ndb_cif_write_category(fp)
FILE *fp;
{
  int icol, irow, datablockId, categoryId, numColumn, numRow, linePos;
  int str_len1, str_len, i;
  char itemName[MxNameLen];
  NdbCifRowFormat *pRow;
  NdbCifCategoryFormat *pCat;
  static int *colwidth = NULL;

  datablockId = ndb_cif_current_datablock();
  categoryId = ndb_cif_current_category();
  numColumn = ndb_cif_count_column();
  numRow = ndb_cif_count_row();
  fprintf(fp, "#\n");
  pCat = &cifFiles.datablocks[datablockId-1].categories[categoryId-1];
  if (numRow <=1) {
    if (numRow == 1)
      pRow = &pCat->rows[0];

    for (icol=0; icol< numColumn; icol++) {
      ndb_cif_get_item_name(icol+1, itemName);
      ndb_cif_print_item_name(fp, itemName, &linePos);
	
      if (numRow == 0 || pRow->columns[icol] == NULL) 
	ndb_cif_print_formatted_item_value(fp, &linePos, "" , TRUE, 1);
      else 
	ndb_cif_print_formatted_item_value(fp, &linePos, pRow->columns[icol], TRUE, 1);
    }
  }
  else {
    fprintf(fp, "loop_\n");

    /* counting colwidth */
    colwidth = (int *)   realloc(colwidth, sizeof(int) *numColumn);
    for (icol=0; icol< numColumn; icol++) 
      colwidth[icol] = 1;
    ndb_cif_rewind_row();
    for (icol=0; icol< numColumn; icol++) {
      for (irow=0; irow< numRow; irow++, ndb_cif_next_row()) {
	pRow = &pCat->rows[irow];
	if (pRow != NULL && pRow->columns[icol] != NULL) {
	  str_len = strlen(pRow->columns[icol]);
	  str_len1 = 0;
	  for (i=0; i<str_len+1; i++)
	    if (pRow->columns[icol][i] == '\n' || pRow->columns[icol][i] == '\0') {
	      if (colwidth[icol] < i - str_len1 +1)
		colwidth[icol] = i-str_len1;
	      str_len1 = i;
	    }
	}
      }
      ndb_cif_rewind_row();
    }

    for (icol=0; icol< numColumn; icol++) {
      ndb_cif_get_item_name(icol+1, itemName);
      ndb_cif_print_item_name(fp, itemName, &linePos);
      fprintf(fp, "\n");
    }
    irow = 1;
    ndb_cif_rewind_row();
    do {
      pRow = &pCat->rows[irow-1];
      linePos = 0;
      for (icol=0; icol< numColumn; icol++) {
	if (icol == numColumn-1)
	  ndb_cif_print_formatted_item_value(fp, &linePos, pRow->columns[icol], TRUE,
					     colwidth[icol]);
	else
	  ndb_cif_print_formatted_item_value(fp, &linePos, pRow->columns[icol], FALSE,
					     colwidth[icol]);
      }
      irow = ndb_cif_next_row();
    }  while (irow);
  }
  fflush(fp);
}

void ndb_cif_print_formatted_item_value(outfile, offset, 
					value, newline_end, fieldLength)
FILE *outfile;
int *offset;
char *value;
int newline_end;
int fieldLength;
/* ---------------------------------------------------------------------
                       print_cif_value
   Rule: 1. If there're single quotes or newlines in value or it exceed 80 
            characters, then it is assumed to be mutiple line value
	    put semi-colon at the begining and end of the value
	    definitely, the *offset will be set to 0
            
         2. If it is less than 80 character there're spaces in the value,
	    then single quotes are put around the value
         3. Single word

	 For 2 and 3, if it will exceed 80 characters after printing out the
         value has to be checked. and *offset has to be updated. 
         If newline_end is true, only multiple-word and single-word values
         need to print a new line
	 
   Parameters: FILE *outfile: file pointer 
               int  *offset:  to keep track the current position of column
	       char *value:   the value to be printed
	       int newline_end: Should a new line to be printed at the end
                                or not. 
 ---------------------------------------------------------------------*/
{
  int str_len, i, pos, len, newline;
  char tmpstring[MxNameLen], fieldFormat[MxNameLen];

  if (value == NULL || !strcmp("",value)) {
    value = (char *) calloc(2, sizeof(char));
    strcpy(value, null_char);
  }

  str_len = strlen(value);
  for (newline=0; newline<str_len; newline++)
    if (value[newline] == '\'' || value[newline] == '\n') break;

  /* Multiple line value, or there're '\n' or '\' in the value */
  if (str_len > 76 || newline != str_len) {
    if (*offset != 0)
      fprintf(outfile, "\n");
    fprintf(outfile, ";\n");
    pos = 0;
    while (pos < str_len) {
      memset(tmpstring, 0, MxNameLen);
      strncpy(tmpstring, &value[pos], MxNameLen-1);
      if ( str_len -pos > MxNameLen)
	len = 78;
      else
	len = str_len - pos;

      i = 0;
      if (newline != str_len && value[newline] == '\n') 
	for (i= len; i> 0; i--) 
	  if (tmpstring[i] == '\n') break;
      /* Did not find newline, in these 80 character, try to find break point */
      if (i == 0) {
	for (i= len; i> 0; i--) 
	  if (!isalnum(tmpstring[i]) && tmpstring[i] != '.') break;
	if (i == 0)
	  i = len;
      }
      memset(tmpstring, 0, MxNameLen);
      if (value[i+pos] == '\n') 
	strncpy(tmpstring, &value[pos], i);
      else 
	strncpy(tmpstring, &value[pos], i+1);
      pos +=(i+1);
      fprintf(outfile, "%s\n",tmpstring);
    }
    fprintf(outfile, ";\n");
    *offset = 0;
  }
  else {

    /* Check if after adding the new value, it will exceed 80 character */
    if (str_len + (*offset) > 75 || fieldLength + (*offset) > 75)  {
      fprintf(outfile, "\n");
      *offset = 0;
    }

    /* Check if it is mutiple word value */
    for (i=0; i<str_len; i++)
      if (value[i] == ' ') break;

    /* 
     * if value is mutiple words put quotes in the front and at the end, 
     * and if it is not at the beginning of the line, leave some spaces 
     * from previous value
     */
    if (i == str_len) {
      if (*offset != 0) 
	sprintf(tmpstring, " %s\0", value);
      else 
	strcpy(tmpstring, value);
    }
    else {
      if (*offset != 0) 
	sprintf(tmpstring, " \'%s\'\0", value);
      else
	sprintf(tmpstring, "\'%s\'\0", value);
    }
    if (str_len < fieldLength) { 
      if ((fieldLength+3) > 80) str_len +=3;
      else str_len = fieldLength+3;
    } else str_len +=3;

    sprintf(fieldFormat,"%%-%ds",str_len);
    fprintf(outfile, fieldFormat, tmpstring);

    (*offset) += (str_len);

    /* If it need a newline at the end of value, put it */
    if (newline_end) fprintf(outfile, "\n");
  }
}


/* General rules for return values:
 * Return values to user 
 *
 * 0 :  Failure
 * 1 :  Success  or
 * >=1: Index of that has been assigned to the request or current index
 *      This index starts from 1 (which is one more than C-based index)
 */


/* ----------------------------------------------------------------
   initialization
   ---------------------------------------------------------------*/

int ndb_cif_init()
{
    int i;
  cifFiles.numDatablock = 0;
  cifFiles.curDatablock = -1;
  cifFiles.allDatablock = CifDefaultSpace ;
  cifFiles.datablocks = (NdbCifDatablockFormat *)
    calloc(cifFiles.allDatablock, sizeof(NdbCifDatablockFormat));
  for (i = 0; i <   cifFiles.allDatablock; i++) {
      cifFiles.datablocks[i].allCategory = 0;
      cifFiles.datablocks[i].numCategory = 0;      
      cifFiles.datablocks[i].curCategory = -1;      
      memset(cifFiles.datablocks[i].datablockName,0,MxNameLen);
  }
  return 1;
}

/* ----------------------------------------------------------------
   close
   ---------------------------------------------------------------*/
int ndb_cif_close()
{
  int i, j, k, l;
  if (cifFiles.datablocks !=NULL) {
    for (i=cifFiles.allDatablock-1; i>=0; i--) {
      for (j=cifFiles.datablocks[i].allCategory-1; j>=0; j--) {
	if (cifFiles.datablocks[i].categories[j].rows != NULL) {
	  for (k= cifFiles.datablocks[i].categories[j].allRow-1; k>=0; k--) {
	    if (cifFiles.datablocks[i].categories[j].rows[k].columns != NULL) {
	      for (l= cifFiles.datablocks[i].categories[j].allCol-1; l>=0; l--) {
		if (cifFiles.datablocks[i].categories[j].rows[k].columns[l] != NULL)
		  free(cifFiles.datablocks[i].categories[j].rows[k].columns[l]);
		cifFiles.datablocks[i].categories[j].rows[k].columns[l] = NULL;
	      }
	      free(cifFiles.datablocks[i].categories[j].rows[k].columns);
	      cifFiles.datablocks[i].categories[j].rows[k].columns = NULL;
	    }
	  }
	  free(cifFiles.datablocks[i].categories[j].rows);
	  cifFiles.datablocks[i].categories[j].rows = NULL;

	  for (l= cifFiles.datablocks[i].categories[j].allCol-1; l>=0; l--) {
	    if (cifFiles.datablocks[i].categories[j].colNames[l] != NULL)
	      free(cifFiles.datablocks[i].categories[j].colNames[l]);
	    cifFiles.datablocks[i].categories[j].colNames[l] = NULL;
	  }
	  free(cifFiles.datablocks[i].categories[j].colNames);
	  cifFiles.datablocks[i].categories[j].colNames = NULL;
	}
      }
      free(cifFiles.datablocks[i].categories);
      cifFiles.datablocks[i].categories = NULL;
    }
    free(cifFiles.datablocks);
    cifFiles.datablocks = NULL;
  }
    
  cifFiles.numDatablock = 0;
  cifFiles.curDatablock = -1;
  cifFiles.allDatablock = 0;

  return 1;
}
/* ----------------------------------------------------------------
   allocation 
   ---------------------------------------------------------------*/
int ndb_cif_new_datablock(datablockName)
char *datablockName;
{
  int i;
  if (cifFiles.numDatablock ==  cifFiles.allDatablock) {
    cifFiles.allDatablock *= 2;
    cifFiles.datablocks = (NdbCifDatablockFormat *)
      realloc(cifFiles.datablocks, 
	      cifFiles.allDatablock * sizeof(NdbCifDatablockFormat));
    for (i=cifFiles.numDatablock; i < cifFiles.allDatablock; i++) {
      cifFiles.datablocks[i].categories = NULL;
      cifFiles.datablocks[i].allCategory = 0;
      cifFiles.datablocks[i].numCategory = 0;
      cifFiles.datablocks[i].curCategory = -1;
      memset(cifFiles.datablocks[i].datablockName,0,MxNameLen);
    }
  }
  cifFiles.curDatablock = cifFiles.numDatablock;
  cifFiles.numDatablock++;
  strcpy(cifFiles.datablocks[cifFiles.curDatablock].datablockName, datablockName);

  /* Allocate spaces for category */
  cifFiles.datablocks[cifFiles.curDatablock].allCategory = CifDefaultSpace ;
  cifFiles.datablocks[cifFiles.curDatablock].categories = (NdbCifCategoryFormat *)
    calloc(cifFiles.datablocks[cifFiles.curDatablock].allCategory,
	   sizeof(NdbCifCategoryFormat));

  for (i=0; i< cifFiles.datablocks[cifFiles.curDatablock].allCategory; i++) {
    cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows = NULL;
    cifFiles.datablocks[cifFiles.curDatablock].categories[i].colNames = NULL;
  }

  cifFiles.datablocks[cifFiles.curDatablock].numCategory = 0;
  cifFiles.datablocks[cifFiles.curDatablock].curCategory = -1;

  return cifFiles.curDatablock +1;
}

int  ndb_cif_put_datablock_name(datablockName)
char *datablockName;
{
  if (cifFiles.numDatablock == 0)
    ndb_cif_new_datablock(datablockName);
  else
    strcpy(cifFiles.datablocks[cifFiles.curDatablock].datablockName,datablockName);
  return(1);
}


int ndb_cif_new_category(categoryName)
 char *categoryName;
{
  int categoryNo,i;

  if (cifFiles.numDatablock == 0)
    ndb_cif_new_datablock("");

  if (cifFiles.datablocks[cifFiles.curDatablock].numCategory ==  
      cifFiles.datablocks[cifFiles.curDatablock].allCategory) {

    cifFiles.datablocks[cifFiles.curDatablock].allCategory *=2;
    cifFiles.datablocks[cifFiles.curDatablock].categories = 
      (NdbCifCategoryFormat *)
      realloc(cifFiles.datablocks[cifFiles.curDatablock].categories, 
	      cifFiles.datablocks[cifFiles.curDatablock].allCategory * 
	      sizeof(NdbCifCategoryFormat));
    for (i=cifFiles.datablocks[cifFiles.curDatablock].numCategory; 
	 i < cifFiles.datablocks[cifFiles.curDatablock].allCategory; i++) {
      cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows = NULL;
      cifFiles.datablocks[cifFiles.curDatablock].categories[i].colNames = NULL;
    }
      
  }
  categoryNo = cifFiles.datablocks[cifFiles.curDatablock].numCategory;
  cifFiles.datablocks[cifFiles.curDatablock].curCategory = categoryNo;
  cifFiles.datablocks[cifFiles.curDatablock].numCategory++;
  cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].numRow = 0;
  cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].numCol = 0;

  /* Allocate space for Rows */
  cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].allRow = 
    CifDefaultSpace;
  cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].rows= 
    (NdbCifRowFormat *) calloc(CifDefaultSpace, sizeof(NdbCifRowFormat));

  cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].curRow = -1;

  /* Allocate space for Column Names */
  cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].allCol = 
    CifDefaultSpace;
  for (i=0; i<CifDefaultSpace; i++) {
    cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].rows[i].columns
      = (char **)calloc(CifDefaultSpace, sizeof(char *));
  }
  cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].colNames = 
    (char **) calloc(CifDefaultSpace, sizeof(char *));
  for (i=0; i<CifDefaultSpace; i++)
    cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].colNames[i] = NULL;
  cifFiles.datablocks[cifFiles.curDatablock].categories[categoryNo].curCol = -1;
  
  strcpy(cifFiles.datablocks[cifFiles.curDatablock].
	 categories[categoryNo].categoryName, categoryName);

  return cifFiles.datablocks[cifFiles.curDatablock].curCategory +1;
}

int ndb_cif_new_row()
{
  int curCategory, i, curRow, curDatablock;
  NdbCifCategoryFormat *pCategory;
  NdbCifRowFormat *pRow;
  
  if (cifFiles.numDatablock == 0 || cifFiles.curDatablock < 0)
    return 0;

  if (cifFiles.datablocks[cifFiles.curDatablock].numCategory == 0)
    return 0;

  curDatablock = cifFiles.curDatablock;
  curCategory = cifFiles.datablocks[curDatablock].curCategory;
  if (curCategory < 0) return 0;
  pCategory = &cifFiles.datablocks[curDatablock].categories[curCategory];
  if (pCategory == NULL) return 0;

  if (pCategory->numRow == pCategory->allRow) {
    pCategory->allRow *=2;
    pCategory->rows = (NdbCifRowFormat *) 
      realloc(pCategory->rows, pCategory->allRow* sizeof(NdbCifRowFormat));
    for (i = pCategory->numRow ; i< pCategory->allRow ; i++)
      pCategory->rows[i].columns = NULL;
  }
  pCategory->curRow = pCategory->numRow;
  curRow = pCategory->curRow;
  pRow = &pCategory->rows[curRow];
  pCategory->numRow ++;
  pRow->columns = (char **) calloc(pCategory->allCol, sizeof(char *));

  for (i=0; i< pCategory->allCol; i++)
    pRow->columns[i] = (char *) calloc(1, sizeof(char));
  
  return pCategory->curRow +1;
}

int ndb_cif_insert_new_row(rowNo)
int rowNo;
{
  int curCategory, i, j, curDatablock, curRow;
  NdbCifCategoryFormat *pCategory;
  NdbCifRowFormat pRow;

  curDatablock = cifFiles.curDatablock;
  if (cifFiles.numDatablock == 0 || cifFiles.curDatablock < 0)
    return 0;

  if (cifFiles.datablocks[curDatablock].numCategory == 0 )
    return 0;

  curCategory = cifFiles.datablocks[curDatablock].curCategory;
  if (curCategory < 0) return 0;
  pCategory = &cifFiles.datablocks[curDatablock].categories[curCategory];
  if (pCategory == NULL) return 0;
  else if (rowNo > pCategory->numRow)
    return 0;

  curRow = ndb_cif_new_row();

  /* 
   * Move the new allocated to correct space
   */
  memcpy(&pRow, &pCategory->rows[curRow-1], sizeof(NdbCifRowFormat));
  for (i=pCategory->numRow-1; i >= rowNo; i--)
    for (j=0 ; j< pCategory->numCol; j++) 
      pCategory->rows[i].columns[j] = pCategory->rows[i-1].columns[j];

  memcpy(&pCategory->rows[rowNo-1], &pRow, sizeof(NdbCifRowFormat));
  pCategory->curRow = rowNo;
  return pCategory->curRow +1;
}

/*
 * Rewind datablock, category, row or field to first one.
 */
int ndb_cif_rewind_datablock()
{
  int curCategory;

  if (cifFiles.numDatablock <= 0) return 0;

  cifFiles.datablocks[cifFiles.curDatablock].curCategory = 0;
  do {
    ndb_cif_rewind_category();
    curCategory = ndb_cif_next_category();
  } while (curCategory);
  cifFiles.datablocks[cifFiles.curDatablock].curCategory = 0;

  return cifFiles.curDatablock+1;

}

int ndb_cif_rewind_category()
{
  if (cifFiles.numDatablock == 0 ||
      cifFiles.datablocks[cifFiles.curDatablock].numCategory == 0)
    return 0;
  ndb_cif_rewind_row();
  ndb_cif_rewind_column();
  return cifFiles.datablocks[cifFiles.curDatablock].curCategory +1;
}

int ndb_cif_rewind_row()
{
  if (cifFiles.numDatablock == 0 ||
      cifFiles.datablocks[cifFiles.curDatablock].numCategory == 0 ||
      cifFiles.datablocks[cifFiles.curDatablock].
      categories[cifFiles.datablocks[cifFiles.curDatablock].curCategory].numRow == 0)
    return 0;
  cifFiles.datablocks[cifFiles.curDatablock].
    categories[cifFiles.datablocks[cifFiles.curDatablock].curCategory].curRow = 0;
  ndb_cif_rewind_column();
  return 1;
}

int ndb_cif_rewind_column()
{
  int curCategory;
  if (cifFiles.numDatablock == 0 ||
      cifFiles.datablocks[cifFiles.curDatablock].numCategory == 0 ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory < 0 ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory)
    return 0;
  curCategory = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  if (cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].numRow == 0)
    return 0;
  cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curCol = 0;
  return cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curCol+ 1;
}

/*
 * Reset/clear up datablock, category
 */
int ndb_cif_reset_datablocks()
{
  int ret;
  ret = ndb_cif_rewind_datablock();
  while (ret) {
    ret = ndb_cif_next_datablock();
    if (ret) ndb_cif_reset_datablock();
  } 
  return ret;
}

int ndb_cif_reset_datablock()
{
  int ret;
  
  if (cifFiles.numDatablock == 0) return 0;

  ret = ndb_cif_rewind_datablock();
  while (ret) {
    ret = ndb_cif_remove_category();
  } 
  return ret;
}

int ndb_cif_reset_datablock_by_id(datablockId)
int datablockId;
{
  int i, ret;

  if (datablockId <=0 ||
      cifFiles.numDatablock == 0 ||
      cifFiles.numDatablock < datablockId) return 0;

  for (i=0; i<cifFiles.datablocks[datablockId-1].numCategory; i++)
    ret  = ndb_cif_reset_category_by_id(datablockId, i+1);
  return ret;

}

int ndb_cif_reset_category()
{
  int i, nrows;
  
  if (cifFiles.numDatablock == 0 ||
      cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].numCategory == 0  ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory )
    return 0;
  

  nrows = cifFiles.datablocks[cifFiles.curDatablock].
    categories[cifFiles.datablocks[cifFiles.curDatablock].curCategory].numRow;
  for (i=nrows-1; i >=0 ; i--) 
    ndb_cif_remove_row_by_id(cifFiles.curDatablock+1,
				   cifFiles.datablocks[cifFiles.curDatablock].
				   curCategory+1, i+1);

  return 1;
}

int ndb_cif_reset_category_by_id(datablockId, categoryId)
int datablockId;
int categoryId;
{
  int ret;
  
  if (datablockId <=0 || categoryId <=0 || 
      cifFiles.numDatablock < datablockId ||
      cifFiles.datablocks[datablockId-1].numCategory < categoryId)
    return 0;
  
  do {
    ret = ndb_cif_remove_row_by_id(datablockId, categoryId, 1);
  } while (ret);
  return 1;
}

/*
 * Remove datablock, category, row or field to first one.
 */

int ndb_cif_remove_datablock()
{

  if (cifFiles.numDatablock == 0) return 0;

  /* Reset free space allocated to this datablock, and memset only copy stuff, not the
   * stuff come with the pointers
   */
  ndb_cif_reset_datablock();

  if (cifFiles.numDatablock > cifFiles.curDatablock) {
    memcpy(&cifFiles.datablocks[cifFiles.curDatablock], 
	   &cifFiles.datablocks[cifFiles.curDatablock+1], 
	   sizeof(NdbCifDatablockFormat) * 
	   (cifFiles.numDatablock - cifFiles.curDatablock - 1));
  }
  cifFiles.numDatablock --;

  if (cifFiles.numDatablock == cifFiles.curDatablock)   cifFiles.curDatablock --;

  return cifFiles.curDatablock+1;
}

int ndb_cif_remove_datablock_by_id(datablockId)
int datablockId;
{
  if (datablockId <= 0 ||
      datablockId > cifFiles.numDatablock)
    return 0;
  if (datablockId-1 == cifFiles.curDatablock)
    return ndb_cif_remove_datablock();

  /* Reset free spaces with this datablock, and memset only copy stuff, not the
   * stuff come with the pointers
   */
  ndb_cif_reset_datablock_by_id(datablockId);
  if (cifFiles.numDatablock > cifFiles.curDatablock) {
    memcpy(&cifFiles.datablocks[datablockId-1], 
	   &cifFiles.datablocks[datablockId], 
	   sizeof(NdbCifDatablockFormat) * 
	   (cifFiles.numDatablock - datablockId));
  }
  cifFiles.numDatablock --;
  return cifFiles.curDatablock +1;
}

int ndb_cif_remove_datablock_by_name(datablockName)
char *datablockName;
{
  int datablockId;

  datablockId = ndb_cif_get_datablock_id(datablockName);
  if (datablockId == 0) return 0;
  else 
    return (ndb_cif_remove_datablock_by_id(datablockId));
}


int ndb_cif_remove_category()
{
  int categoryId;

  if (cifFiles.numDatablock == 0 ||
      cifFiles.datablocks[cifFiles.curDatablock].numCategory == 0) 
    return 0;

  categoryId = cifFiles.datablocks[cifFiles.curDatablock].curCategory+1;
  ndb_cif_reset_category();

  memcpy(&cifFiles.datablocks[cifFiles.curDatablock].categories[categoryId-1],
	 &cifFiles.datablocks[cifFiles.curDatablock].categories[categoryId],
	 sizeof(NdbCifCategoryFormat) * 
	 (cifFiles.datablocks[cifFiles.curDatablock].numCategory - categoryId+1));

  cifFiles.datablocks[cifFiles.curDatablock].numCategory --;
  
  if (cifFiles.datablocks[cifFiles.curDatablock].numCategory == 
      cifFiles.datablocks[cifFiles.curDatablock].curCategory)
    cifFiles.datablocks[cifFiles.curDatablock].curCategory--;
  
  return cifFiles.datablocks[cifFiles.curDatablock].curCategory +1;
}

int ndb_cif_remove_category_by_id(datablockId, categoryId)
int datablockId;
int categoryId;
{

  if (datablockId <=0 || categoryId <=0 ||
      datablockId > cifFiles.numDatablock ||
      categoryId > cifFiles.datablocks[cifFiles.curDatablock].numCategory) 
    return 0;

  ndb_cif_reset_category_by_id(datablockId, categoryId);
  if (datablockId == cifFiles.curDatablock &&
      categoryId  == cifFiles.datablocks[cifFiles.curDatablock].curCategory)
    return ndb_cif_remove_category();
  memcpy(&cifFiles.datablocks[cifFiles.curDatablock].
	 categories[categoryId],
	 &cifFiles.datablocks[cifFiles.curDatablock].
	 categories[categoryId+1],
	 sizeof(NdbCifCategoryFormat) * 
	 (cifFiles.datablocks[cifFiles.curDatablock].numCategory - categoryId-1));
  cifFiles.datablocks[cifFiles.curDatablock].numCategory --;
  return cifFiles.datablocks[cifFiles.curDatablock].curCategory +1;
}

int ndb_cif_remove_category_by_name(datablockName, categoryName)
char *datablockName;
char *categoryName;
{
  int datablockId, categoryId;
  
  datablockId = ndb_cif_get_datablock_id(datablockName);
  if (datablockId == 0) 
    return 0;
  categoryId = ndb_cif_get_category_id(datablockName, categoryName);
  if (categoryId == 0) return 0;
  else return ndb_cif_remove_category_by_id(datablockId, categoryId);
}

int ndb_cif_remove_row()
{
  int datablockId, categoryId, rowId;

  if (cifFiles.numDatablock <= 0 ||
      cifFiles.curDatablock < 0 ||
      cifFiles.curDatablock >= cifFiles.numDatablock) 
    return 0;
  datablockId = cifFiles.curDatablock;

  if (cifFiles.datablocks[datablockId].numCategory <=0 ||
      cifFiles.datablocks[datablockId].curCategory < 0 ||
      cifFiles.datablocks[datablockId].curCategory >=
      cifFiles.datablocks[datablockId].numCategory )
    return 0;

  categoryId  = cifFiles.datablocks[datablockId].curCategory;
  rowId = cifFiles.datablocks[datablockId].categories[categoryId].curRow; 
    
  return ndb_cif_remove_row_by_id(datablockId+1, categoryId+1, rowId+1);
				 
}

int ndb_cif_remove_row_by_id(datablockId, categoryId, rowId)
int datablockId;
int categoryId;
int rowId;
{

  NdbCifCategoryFormat  *pCategory;
  int i, j;

  if (datablockId <=0 || categoryId <=0 || rowId <=0 ||
      cifFiles.numDatablock < datablockId ||
      cifFiles.datablocks[datablockId-1].numCategory < categoryId ||
      cifFiles.datablocks[datablockId-1].categories[categoryId-1].numRow < rowId)
    return 0;
  
  pCategory = &cifFiles.datablocks[datablockId-1].categories[categoryId-1];

  
  /* 
   * Move the new allocated to correct space
   */
  for (i= 0; i < pCategory->numCol; i++)
    if (pCategory->rows[rowId-1].columns[i] != NULL) {
      free(pCategory->rows[rowId-1].columns[i]);
      pCategory->rows[rowId-1].columns[i] = NULL;
    }
  
  for (i= rowId-1; i < pCategory->numRow-1; i++) {
    for (j=0 ; j< pCategory->numCol; j++) 
      pCategory->rows[i].columns[j] = pCategory->rows[i+1].columns[j];
  }
  for (j=0 ; j< pCategory->numCol; j++) 
    pCategory->rows[pCategory->numRow-1].columns[j] = NULL;
  pCategory->numRow--;
  if (pCategory->curRow == pCategory->numRow)
    pCategory->curRow --;
  return pCategory->curRow+1;
}



int ndb_cif_next_datablock()
{
  if (cifFiles.curDatablock < cifFiles.numDatablock-1) {
    cifFiles.curDatablock++;
    cifFiles.datablocks[cifFiles.curDatablock].curCategory = 0;
    ndb_cif_rewind_category();
    return cifFiles.curDatablock+1;  
  }  else 
    return 0;
}

int ndb_cif_next_category()
{
  if (cifFiles.datablocks == NULL) return 0;
  
  if (cifFiles.datablocks[cifFiles.curDatablock].curCategory < 
      cifFiles.datablocks[cifFiles.curDatablock].numCategory-1) {

    cifFiles.datablocks[cifFiles.curDatablock].curCategory++;
    ndb_cif_rewind_row();
    return cifFiles.datablocks[cifFiles.curDatablock].curCategory+1;
  }  else 
    return 0;
}

int ndb_cif_next_row()
{
  NdbCifCategoryFormat *pCategory;
  int curCategory;

  if (cifFiles.numDatablock == 0)  return 0;
  
  if (cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory) 
    return 0;

  curCategory = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  pCategory = &cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory];

  if (pCategory->curRow >= 0 && pCategory->curRow < pCategory->numRow - 1) {
    ++pCategory->curRow;
    return pCategory->curRow+1;
  }
  else return 0;
}

int ndb_cif_next_column()
{
  int curCategory;
  if (cifFiles.datablocks == NULL)
    return 0;
  
  if (cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory) 
    return 0;
  curCategory = cifFiles.datablocks[cifFiles.curDatablock].curCategory;

  if (cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curCol <
      cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].numCol-1) {
    cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curCol++;
    return cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curCol+1;
  }
  else return 0;
}

int ndb_cif_move_datablock(datablockName)
char *datablockName;
{
  int i;
  if (cifFiles.datablocks == NULL)
    return 0;
		 
  for (i=0; i< cifFiles.numDatablock; i++)
    if (!strcmp(cifFiles.datablocks[i].datablockName,datablockName)) break;
  
  if (i != cifFiles.numDatablock) {
    cifFiles.curDatablock = i;
    return i+1;
  }
  else 
    return 0;
}

int ndb_cif_move_datablock_by_id(datablockId)
int datablockId;
{

  if (datablockId < 0 ||
      datablockId > cifFiles.numDatablock)
    return 0;
		 
  cifFiles.curDatablock = datablockId-1;
  return datablockId;
}

int ndb_cif_move_category_by_name(datablockName, categoryName)
char *datablockName;
char *categoryName;
{
  int i, curDatablock;

  if ((curDatablock = ndb_cif_move_datablock(datablockName)) == 0) 
    return 0;

  for (i=0; i< cifFiles.datablocks[curDatablock-1].numCategory; i++)
    if (!strcmp(cifFiles.datablocks[curDatablock-1].categories[i].categoryName,
		categoryName)) break;
  
  if (i != cifFiles.datablocks[curDatablock-1].numCategory) {
    cifFiles.datablocks[cifFiles.curDatablock].curCategory = i;
/*
    cifFiles.datablocks[cifFiles.curDatablock].categories[i].curRow = 0;
    cifFiles.datablocks[cifFiles.curDatablock].categories[i].curCol= 0;
*/
    return i+1;
  }
  else 
    return 0;
}

int ndb_cif_move_category_by_id(datablockId, categoryId)
int datablockId;
int categoryId;
{
  int  curDatablock;

  if ((curDatablock = ndb_cif_move_datablock_by_id(datablockId)) == 0) 
    return 0;

  if (categoryId <0 || categoryId > cifFiles.datablocks[curDatablock-1].numCategory)
    return 0;


  cifFiles.datablocks[curDatablock-1].curCategory = categoryId-1;
  return categoryId;

}

int ndb_cif_move_category(categoryId)
int categoryId;
{

  if (categoryId <0 || categoryId > cifFiles.datablocks[cifFiles.curDatablock].numCategory)
    return 0;


  cifFiles.datablocks[cifFiles.curDatablock].curCategory = categoryId-1;
  return categoryId;

}

int ndb_cif_move_row_by_name(datablockName, categoryName, rowId)
char *datablockName;
char *categoryName;
int rowId;
{
  int curCategory, curDatablock;

  if ((curCategory = ndb_cif_move_category_by_name(datablockName, categoryName)) == 0) 
    return 0;
  curDatablock = cifFiles.curDatablock;
  
  if (rowId < 0 ||
      rowId > cifFiles.datablocks[curDatablock].categories[curCategory].numRow)
    return 0;

  cifFiles.datablocks[curDatablock].categories[curCategory].curRow = rowId-1;
  return rowId+1;
}

int ndb_cif_move_row_by_id(datablockId, categoryId, rowId)
int datablockId;
int categoryId;
int rowId;
{
  if (datablockId <=0 || categoryId <=0 || rowId <=0 ||
      ndb_cif_move_category_by_id(datablockId, categoryId) == 0) 
    return 0;

  if (rowId <0 || 
      rowId > cifFiles.datablocks[datablockId-1].categories[categoryId-1].numRow)
    return 0;

  cifFiles.datablocks[datablockId-1].categories[categoryId-1].curRow = rowId -1;
  return rowId;

}

int ndb_cif_move_row(rowId)
int rowId;
{
  int  curDatablock, curCategory;

  curDatablock = cifFiles.curDatablock;
  curCategory= cifFiles.datablocks[curDatablock].curCategory;

  return  ndb_cif_move_row_by_id(curDatablock+1, curCategory+1, rowId);

}

int ndb_cif_count_datablock()
{
  return cifFiles.numDatablock;
}

int ndb_cif_count_category()
{
  if (cifFiles.curDatablock >= cifFiles.numDatablock)
    return 0;

  return cifFiles.datablocks[cifFiles.curDatablock].numCategory;
}


int ndb_cif_count_row()
{
  int curDatablock, curCategory;

  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory )
    return 0;

  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  return cifFiles.datablocks[curDatablock].categories[curCategory].numRow;
}


int ndb_cif_count_column()
{
  int curDatablock, curCategory;

  if (cifFiles.curDatablock < 0 ||
      cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory < 0 ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory )
    return 0;

  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  return cifFiles.datablocks[curDatablock].categories[curCategory].numCol;
}


int ndb_cif_current_datablock()
{
  return cifFiles.curDatablock+1;
}

int ndb_cif_current_datablock_name(datablockName)
char *datablockName;
{
  if (cifFiles.curDatablock < 0 ||
      cifFiles.curDatablock >= cifFiles.numDatablock)
    return 0;
  strcpy(datablockName, cifFiles.datablocks[cifFiles.curDatablock].datablockName);
  return cifFiles.curDatablock+1;
}

int ndb_cif_current_category()
{
  if (cifFiles.curDatablock >= cifFiles.numDatablock)
    return 0;

  return cifFiles.datablocks[cifFiles.curDatablock].curCategory+1;
}

int ndb_cif_current_category_name(categoryName)
char *categoryName;
{
  if (cifFiles.curDatablock >= cifFiles.numDatablock)
    return 0;

  strcpy(categoryName, cifFiles.datablocks[cifFiles.curDatablock].categories
	 [cifFiles.datablocks[cifFiles.curDatablock].curCategory].categoryName);
  return cifFiles.datablocks[cifFiles.curDatablock].curCategory+1;
}

int ndb_cif_current_row()
{

  int curDatablock, curCategory;

  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory )
    return 0;

  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  return cifFiles.datablocks[curDatablock].categories[curCategory].curRow+1;
}

int ndb_cif_current_col()
{
  int curDatablock, curCategory;

  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory )
    return 0;

  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  return cifFiles.datablocks[curDatablock].categories[curCategory].curCol+1;
}

int ndb_cif_put_item_keyword(itemKeyword)
char *itemKeyword;
{ 
  NdbCifCategoryFormat *pCategory;
  int curDatablock, curCategory, i, j;

  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory )
    return 0;

  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[curDatablock].curCategory;
  pCategory  = &cifFiles.datablocks[curDatablock].categories[curCategory];

  if (pCategory->numCol == pCategory->allCol) {
    pCategory->allCol *=2;
    pCategory->colNames = (char **) 
      realloc(pCategory->colNames, pCategory->allCol* sizeof(char *));
    for (i=pCategory->numCol; i<pCategory->allCol; i++)
      pCategory->colNames[i] = NULL;
    for (i=0; i<pCategory->allRow; i++) {
      pCategory->rows[i].columns = (char **)
	realloc(pCategory->rows[i].columns, pCategory->allCol * sizeof(char**));
      for (j=pCategory->numCol; j<pCategory->allCol; j++) {
	pCategory->rows[i].columns[j] = NULL;
      }
    }
  }
  pCategory->curCol = pCategory->numCol;
  pCategory->numCol++;
      
  pCategory->colNames[pCategory->curCol] = (char *)
    realloc(pCategory->colNames[pCategory->curCol], 
	    sizeof(char) * (strlen(itemKeyword)+1));
  strcpy(pCategory->colNames[pCategory->curCol], itemKeyword);
  return pCategory->curCol+1;
}

int ndb_cif_remove_item_keyword(itemKeyword)     
char *itemKeyword;
/* Given the name of a column, remove the keyword and all of
   the values from the current category.
   */     
{ 
  NdbCifCategoryFormat *pCategory;
  int curDatablock, curCategory, i, j;
  int ind;
  
  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory )
    return 0;

  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[curDatablock].curCategory;
  pCategory  = &cifFiles.datablocks[curDatablock].categories[curCategory];

  for (i = 0; i < 	pCategory->numCol; i++)
    if (!strcmp(pCategory->colNames[i], itemKeyword)) break;
  if (i == pCategory->numCol)
    return 0;

  ind = i;
  free(pCategory->colNames[ind]);
  for (j = ind; j < pCategory->numCol - 1; j++)
    pCategory->colNames[j] = pCategory->colNames[j + 1]; 
  pCategory->colNames[pCategory->numCol - 1] = NULL;
  
  for (i=0; i<pCategory->allRow; i++) 
    {
      for (j = ind; j < pCategory->numCol - 1; j++) 
	pCategory->rows[i].columns[j] = pCategory->rows[i].columns[j + 1];
      pCategory->rows[i].columns[pCategory->numCol - 1] = NULL;
    }
  return pCategory->numCol--; 
}


int ndb_cif_get_item_name(colId, itemName)
int colId;
char *itemName;
{
  int curDatablock, curCategory;
  NdbCifCategoryFormat *pCategory;
  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory)
    return 0;
  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;

  if (colId > cifFiles.datablocks[curDatablock].categories[curCategory].numCol ||
      cifFiles.datablocks[curDatablock].categories[curCategory].curRow >=
      cifFiles.datablocks[curDatablock].categories[curCategory].numRow)
    return 0;

  pCategory = & cifFiles.datablocks[curDatablock].categories[curCategory];
  sprintf(itemName, "_%s.%s", pCategory->categoryName, 
	  pCategory->colNames[colId-1]);
  return TRUE;
}

int ndb_cif_get_item_shname(colId, itemName)
int colId;
char *itemName;
     /* Retur the name minus the category name (the SHort name. 
      */     
{
  int curDatablock, curCategory;
  NdbCifCategoryFormat *pCategory;
  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory)
    return 0;
  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;

  if (colId > cifFiles.datablocks[curDatablock].categories[curCategory].numCol ||
      cifFiles.datablocks[curDatablock].categories[curCategory].curRow >=
      cifFiles.datablocks[curDatablock].categories[curCategory].numRow)
    return 0;

  pCategory = & cifFiles.datablocks[curDatablock].categories[curCategory];
  sprintf(itemName, "%s", pCategory->colNames[colId-1]);
  return TRUE;
}

int ndb_cif_get_item_value(colId, fieldValue, maxFieldLen)
int colId;
char *fieldValue;
int maxFieldLen;
{
  NdbCifRowFormat *pRow;
  int curDatablock, curCategory, curRow;

/*  memset(fieldValue, NULL, maxFieldLen); */
  fieldValue[0] = '\0';
  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory)
    return 0;
  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  if (colId > cifFiles.datablocks[curDatablock].categories[curCategory].numCol ||
      cifFiles.datablocks[curDatablock].categories[curCategory].numRow == 0 ||
      cifFiles.datablocks[curDatablock].categories[curCategory].curRow >=
      cifFiles.datablocks[curDatablock].categories[curCategory].numRow)
    return 0;

  curRow  = cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curRow;
  if (curRow == -1) return 0;
  pRow = &cifFiles.datablocks[curDatablock].categories[curCategory].rows[curRow];

  if (pRow->columns[colId-1] != NULL) 
    strncpy(fieldValue, pRow->columns[colId-1], maxFieldLen-1);
  fieldValue[maxFieldLen-1] = '\0';
  return TRUE;
}

char *ndb_cif_copy_item_value(colId)
int colId;
  /* In the current row,return a newly allocated string that
     contains a copy of the valule in  colID.
   */
{
 NdbCifRowFormat *pRow;
 int len, curDatablock, curCategory, curRow;
 char *ret = NULL;
 
 if (cifFiles.curDatablock >= cifFiles.numDatablock ||
     cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
     cifFiles.datablocks[cifFiles.curDatablock].numCategory)
   return NULL;
 curDatablock = cifFiles.curDatablock;
 curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
 if (colId > cifFiles.datablocks[curDatablock].categories[curCategory].numCol ||
     cifFiles.datablocks[curDatablock].categories[curCategory].numRow == 0 ||
     cifFiles.datablocks[curDatablock].categories[curCategory].curRow >=
     cifFiles.datablocks[curDatablock].categories[curCategory].numRow)
   return NULL;
 
 curRow  = cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curRow;
 if (curRow == -1) return NULL;
 pRow = &cifFiles.datablocks[curDatablock].categories[curCategory].rows[curRow];
 
 if (pRow->columns[colId-1] != NULL) 
   {
    len = strlen(pRow->columns[colId-1]);
    ret = (char *) calloc(len + 1, sizeof(char));
    strcpy(ret, pRow->columns[colId-1]);
    return ret;
   }
 return NULL;
}


int ndb_cif_output_item(fp, colId)
FILE *fp;
int colId;
  /* In the current row, output the value in colID to port fp.
   */
{
 int i;
 NdbCifRowFormat *pRow;
 int curDatablock, curCategory, curRow;
 
 /*  memset(fieldValue, NULL, maxFieldLen); */
 if (cifFiles.curDatablock >= cifFiles.numDatablock ||
     cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
     cifFiles.datablocks[cifFiles.curDatablock].numCategory)
   return 0;
 curDatablock = cifFiles.curDatablock;
 curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
 if (colId > cifFiles.datablocks[curDatablock].categories[curCategory].numCol ||
     cifFiles.datablocks[curDatablock].categories[curCategory].numRow == 0 ||
     cifFiles.datablocks[curDatablock].categories[curCategory].curRow >=
     cifFiles.datablocks[curDatablock].categories[curCategory].numRow)
   return 0;
 
 curRow  = cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curRow;
 if (curRow == -1) return 0;
 pRow = &cifFiles.datablocks[curDatablock].categories[curCategory].rows[curRow];
 
 if (pRow->columns[colId-1] != NULL) 
   {
    i = 0;
    while (pRow->columns[colId-1][i] != '\0')
      fputc(pRow->columns[colId-1][i++], fp);
    return TRUE;
   }
 else 
   return 0;
}

int ndb_cif_item_value_strncmp(catId, colId, rowId, startPos,
			       fieldValue, maxFieldLen)
int catId;
int colId;
int rowId;
int startPos;
char *fieldValue;
int maxFieldLen;
{
 NdbCifRowFormat *pRow;
 int curDatablock;
 
 if (cifFiles.curDatablock >= cifFiles.numDatablock ||
     cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
     cifFiles.datablocks[cifFiles.curDatablock].numCategory)
   return -1;
 curDatablock = cifFiles.curDatablock;
 
 if (catId <= 0 || catId > cifFiles.datablocks[curDatablock].numCategory)
   return -1;
 if (colId <= 0 ||
     colId > cifFiles.datablocks[curDatablock].categories[catId-1].numCol ||
     rowId <= 0 ||
     rowId > cifFiles.datablocks[curDatablock].categories[catId-1].numRow )
   return -1;
 
 pRow = &cifFiles.datablocks[curDatablock].categories[catId-1].rows[rowId-1];
 
 /*  if (startPos <= 0 || startPos > strlen(pRow->columns[colId-1])) return -1; */
 
 if (pRow->columns[colId-1] != NULL) 
   return (strncmp(fieldValue, &pRow->columns[colId-1][startPos-1], maxFieldLen));
 return(-1);
}

int ndb_cif_item_values_strncmp(catId1, colId1, rowId1, startPos1, 
				catId2, colId2,  rowId2,  startPos2, maxFieldLen)
int catId1;
int colId1;
int rowId1;
int startPos1;
int catId2;
int colId2;
int rowId2;
int startPos2;
int maxFieldLen;
{
  NdbCifRowFormat *pRow1, *pRow2;
  int curDatablock;

  if (cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory)
    return -1;
  curDatablock = cifFiles.curDatablock;

  if (catId1 <= 0 || catId1 > cifFiles.datablocks[curDatablock].numCategory ||
      catId2 <= 0 || catId2 > cifFiles.datablocks[curDatablock].numCategory )
    return -1;
  if (colId1 <= 0 ||
      colId1 > cifFiles.datablocks[curDatablock].categories[catId1-1].numCol ||
      rowId1 <= 0 ||
      rowId1 > cifFiles.datablocks[curDatablock].categories[catId1-1].numRow ||
      colId2 <= 0 ||
      colId2 > cifFiles.datablocks[curDatablock].categories[catId2-1].numCol ||
      rowId2 <= 0 ||
      rowId2 > cifFiles.datablocks[curDatablock].categories[catId2-1].numRow)
    return -1;

  pRow1 = &cifFiles.datablocks[curDatablock].categories[catId1-1].rows[rowId1-1];
  pRow2 = &cifFiles.datablocks[curDatablock].categories[catId2-1].rows[rowId2-1];
/*
  if (startPos1 <= 0 || startPos1 > strlen(pRow1->columns[colId1-1]) ||
      startPos2 <= 0 || startPos2 > strlen(pRow2->columns[colId2-1]))
    return -1;
  */    
  return strncmp(&pRow1->columns[colId1-1][startPos1-1],
		 &pRow2->columns[colId2-1][startPos2-1], maxFieldLen);
}

/*
   The following two functions are like ndb_cif_item_values_strncmp
   but return 0 or i if the string(s) are found in the given 
   columns of the given row,using strcmp for comparison, and not specifying
   the starting column of the string searched against.
   */

int ndb_cif_item_row_1_key(colId, fieldValue)
int colId;
char *fieldValue;
  /* In the current category, examine all the rows and return the row number if fieldValue
     is found, and also also set the current row to i. If the row is not found return 0.
     */
{
 int i;
 int curDatablock, curCategory;
 NdbCifRowFormat *pRow;
 NdbCifCategoryFormat *pCategory;

 if (cifFiles.numDatablock == 0)  return 0;
 if (cifFiles.curDatablock >= cifFiles.numDatablock ||
     (cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory))
   return 0;

  curDatablock = cifFiles.curDatablock;
 curCategory = cifFiles.datablocks[curDatablock].curCategory;
 pCategory = &cifFiles.datablocks[curDatablock].categories[curCategory];
                  
 if (colId <= 0 ||
     colId >  pCategory->numCol)
   return 0;

 for (i = 0; i < pCategory->numRow; i++)
   {
    pRow = &pCategory->rows[i];
    if ((pRow->columns[colId-1] != NULL) 
      && !strcmp(fieldValue, pRow->columns[colId-1]))
      {
       pCategory->curRow = i;    
       return i + 1;
      }
   }
 return 0;
}

int ndb_cif_item_row_2_keys(colId1, fieldValue1, colId2, fieldValue2)
int colId1;
char *fieldValue1;
int colId2;
char *fieldValue2;
  /* In the current category, examine all the rows and return the row number if fieldValue1 
     and 2 are found in the respecive columns of row i, and also also set the current
     row to i. If the row is not found, return 0.
     */
{
 int i;
 int curDatablock, curCategory;
 NdbCifRowFormat *pRow;
 NdbCifCategoryFormat *pCategory;

 if (cifFiles.numDatablock == 0)  return 0;
 if (cifFiles.curDatablock >= cifFiles.numDatablock ||
     (cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory))
   return 0;

  curDatablock = cifFiles.curDatablock;
 curCategory = cifFiles.datablocks[curDatablock].curCategory;
 pCategory = &cifFiles.datablocks[curDatablock].categories[curCategory];
                  
 if ((colId1 <= 0) || (colId1 >  pCategory->numCol) || 
     (colId2 <= 0) || (colId2 >  pCategory->numCol))
   return 0;

 for (i = 0; i < pCategory->numRow; i++)
   {
    pRow = &pCategory->rows[i];
    if ((pRow->columns[colId1-1] != NULL) && (pRow->columns[colId2-1] != NULL))
      if (!strcmp(fieldValue1, pRow->columns[colId1-1]))
	if (!strcmp(fieldValue2, pRow->columns[colId2-1]))
	  {
	   pCategory->curRow = i;    
	   return i + 1;
	  }
   }
 return 0;
}

int ndb_cif_item_row_3_keys(colId1, fieldValue1, colId2, fieldValue2,
			    colId3, fieldValue3)
int colId1;
char *fieldValue1;
int colId2;
char *fieldValue2;
int colId3;
char *fieldValue3;

  /* In the current category, examine all the rows and return the row number if fieldValue1 
     2 and 3 are found in the respective columns of row i, and also also set the current
     row to i. If the row is not found, return 0.
     */
{
 int i;
 int curDatablock, curCategory;
 NdbCifRowFormat *pRow;
 NdbCifCategoryFormat *pCategory;

 if (cifFiles.numDatablock == 0)  return 0;
 if (cifFiles.curDatablock >= cifFiles.numDatablock ||
     (cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory))
   return 0;

  curDatablock = cifFiles.curDatablock;
 curCategory = cifFiles.datablocks[curDatablock].curCategory;
 pCategory = &cifFiles.datablocks[curDatablock].categories[curCategory];
                  
 if ((colId1 <= 0) || (colId1 >  pCategory->numCol) || 
     (colId2 <= 0) || (colId2 >  pCategory->numCol) ||
     (colId3 <= 0) || (colId3 >  pCategory->numCol))
   return 0;

 for (i = 0; i < pCategory->numRow; i++)
   {
    pRow = &pCategory->rows[i];
    if ((pRow->columns[colId1-1] != NULL) && (pRow->columns[colId2-1] != NULL))
      if (!strcmp(fieldValue1, pRow->columns[colId1-1]))
	if (!strcmp(fieldValue2, pRow->columns[colId2-1]))
	  if (!strcmp(fieldValue3, pRow->columns[colId3-1]))
	  {
	   pCategory->curRow = i;    
	   return i + 1;
	  }
   }
 return 0;
}

int ndb_cif_item_row_4_keys(colId1, fieldValue1, colId2, fieldValue2,
			    colId3, fieldValue3, colId4, fieldValue4)
int colId1, colId2, colId3, colId4;
char *fieldValue1, *fieldValue2, *fieldValue3, *fieldValue4;
  /* In the current category, examine all the rows and return the row number if fieldValue1 
     2, 3 and 4 are found in the respective columns of row i, and also also set the current
     row to i. If the row is not found, return 0.
     */
{
 int i;
 int curDatablock, curCategory;
 NdbCifRowFormat *pRow;
 NdbCifCategoryFormat *pCategory;

 if (cifFiles.numDatablock == 0)  return 0;
 if (cifFiles.curDatablock >= cifFiles.numDatablock ||
     (cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory))
   return 0;

  curDatablock = cifFiles.curDatablock;
 curCategory = cifFiles.datablocks[curDatablock].curCategory;
 pCategory = &cifFiles.datablocks[curDatablock].categories[curCategory];
                  
 if ((colId1 <= 0) || (colId1 >  pCategory->numCol) || 
     (colId2 <= 0) || (colId2 >  pCategory->numCol) ||
     (colId3 <= 0) || (colId3 >  pCategory->numCol) ||
     (colId4 <= 0) || (colId4 >  pCategory->numCol))
   return 0;

 for (i = 0; i < pCategory->numRow; i++)
   {
    pRow = &pCategory->rows[i];
    if ((pRow->columns[colId1-1] != NULL) && (pRow->columns[colId2-1] != NULL))
      if (!strcmp(fieldValue1, pRow->columns[colId1-1]))
	if (!strcmp(fieldValue2, pRow->columns[colId2-1]))
	  if (!strcmp(fieldValue3, pRow->columns[colId3-1]))
	    if (!strcmp(fieldValue4, pRow->columns[colId4-1]))
	      {
	       pCategory->curRow = i;    
	       return i + 1;
	      }
   }
 return 0;
}

int ndb_cif_put_item_value(colId, fieldValue)
int colId;
char *fieldValue;
{
  NdbCifRowFormat *pRow;
  int curDatablock, curCategory, curRow, str_len;
  if (colId <= 0 ){
    fprintf(stderr, "Error: index of %d passed to ndb_cif_put_item_value.\n",
	    colId);
    exit(1);
  }
  if (cifFiles.curDatablock< 0 ||
      cifFiles.curDatablock >= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory <0 ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory >=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory)
    return 0;
  curDatablock = cifFiles.curDatablock;
  curCategory  = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  if (colId > cifFiles.datablocks[curDatablock].categories[curCategory].numCol ||
      cifFiles.datablocks[curDatablock].categories[curCategory].curRow >=
      cifFiles.datablocks[curDatablock].categories[curCategory].numRow)
    return 0;

  curRow  = cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].curRow;
  if (curRow == -1) {
    curRow = ndb_cif_new_row();
    curRow --;
  }
  pRow = &cifFiles.datablocks[curDatablock].categories[curCategory].rows[curRow];

  str_len = strlen(fieldValue);
  if (pRow->columns[colId-1] == NULL || str_len > strlen(pRow->columns[colId-1]))
    pRow->columns[colId-1] = (char *)
      realloc( pRow->columns[colId-1], (str_len+1)*sizeof(char));
  strcpy(pRow->columns[colId-1], fieldValue);
  return TRUE;
}

int ndb_cif_get_category_name_from_item_name(categoryName, itemName)
char *categoryName;
char *itemName;
{
  int i, str_len;

  str_len = strlen(itemName);
    for (i=0; i<str_len; i++)
      if (itemName[i] == '.') break;
  if (i == str_len)  {
    categoryName[0] = '\0';
    return 0;
  }
  strncpy(categoryName, &itemName[1],i-1);    
  categoryName[i-1] = '\0';
  if (cifFiles.curDatablock >= cifFiles.numDatablock) return 0;

  for (i=0; i< cifFiles.datablocks[cifFiles.curDatablock].numCategory; i++)
    if (!strcmp(categoryName, 
		cifFiles.datablocks[cifFiles.curDatablock].categories[i].categoryName))
      break;
  if (i == cifFiles.datablocks[cifFiles.curDatablock].numCategory) return 0;
  else return i+1;

}

int ndb_cif_get_item_keyword_from_item_name(itemKeyword, itemName)
char *itemKeyword;
char *itemName;
{
  int i, str_len, curCategory;

  str_len = strlen(itemName);
    for (i=0; i<str_len; i++)
      if (itemName[i] == '.') break;
  if (i == str_len)  {
    itemKeyword[0] = '\0';
    return 0;
  }
  strcpy(itemKeyword, &itemName[i+1]);    
  if (cifFiles.curDatablock <= cifFiles.numDatablock ||
      cifFiles.datablocks[cifFiles.curDatablock].curCategory <=
      cifFiles.datablocks[cifFiles.curDatablock].numCategory)
    return 0;

  curCategory = cifFiles.datablocks[cifFiles.curDatablock].curCategory;
  for (i=0; i< cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].numCol; i++)
    if (!strcmp(itemKeyword, 
		cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].
		colNames[i]))
      break;
  if (i == cifFiles.datablocks[cifFiles.curDatablock].categories[curCategory].numCol) 
    return 0;
  else return i+1;

}

int ndb_cif_get_datablock_id(datablockName)
char *datablockName;
{
  int i;

  if (cifFiles.datablocks == NULL ||
      cifFiles.curDatablock >= cifFiles.numDatablock)
    return 0;
		 
 
  for (i=0; i< cifFiles.numDatablock; i++)
    if (!strcmp(cifFiles.datablocks[i].datablockName, datablockName)) break;
  if (i != cifFiles.numDatablock)
    return i +1;
  else 
    return 0;
}

int ndb_cif_get_category_id(datablockName, categoryName)
char *datablockName;
char *categoryName;
{
  int i, datablockId;

  datablockId = ndb_cif_get_datablock_id(datablockName);
  if (datablockId != 0) {
    for (i=0; i< cifFiles.datablocks[datablockId-1].numCategory; i++)
      if (!strcmp(cifFiles.datablocks[datablockId-1].categories[i].categoryName,
		  categoryName)) break;
    if (i != cifFiles.datablocks[datablockId-1].numCategory) 
      return i+1;
    else 
      return 0;
  }
  else return 0;
}

int ndb_cif_get_column_id(datablockName, categoryName, itemKeyword)
char *datablockName;
char *categoryName;
char *itemKeyword;
{
  int i, datablockId, catId;

  datablockId = ndb_cif_get_datablock_id(datablockName);
  if (datablockId != 0)
  catId = ndb_cif_get_category_id(datablockName, categoryName);
  if (datablockId != 0 && catId != 0) {
    for (i=0; i< cifFiles.datablocks[datablockId-1].categories[catId-1].numCol; i++)
      if (!strcmp(cifFiles.datablocks[datablockId-1].categories[catId-1].colNames[i],
		  itemKeyword)) break;
    if (i != cifFiles.datablocks[datablockId-1].categories[catId-1].numCol)
      return i+1;
    else 
      return 0;
  }
  else return 0;
}


void ndb_cif_print_datablock(fp)
FILE *fp;
{
    int i, j, k;
    if (cifFiles.numDatablock == 0 ||
	cifFiles.curDatablock < 0 || 
	cifFiles.curDatablock >= cifFiles.numDatablock)
	return;

    fprintf(fp, "Data Block: %s\n", cifFiles.datablocks[cifFiles.curDatablock].datablockName);

    for (i=0; i< cifFiles.datablocks[cifFiles.curDatablock].numCategory; i++) {
	fprintf(fp,"Category: %s\n", cifFiles.datablocks[cifFiles.curDatablock].categories[i].categoryName);

	for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) 
	    fprintf(fp, "%s  ", cifFiles.datablocks[cifFiles.curDatablock].categories[i].colNames[j]);
	fprintf(fp, "\n");
	for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) 
	    fprintf(fp, "-------");
	fprintf(fp, "\n");
	for (k=0; k< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numRow; k++)  {
	    for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) { 
		if (cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j] == NULL ||
		    !strcmp(cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j],""))
	  
		    fprintf(fp, "%s ", "(null)");
		else
		    fprintf(fp, "%s ", cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j]);
	    }
	    fprintf(fp, "\n");
	}
    }
}


void ndb_cif_pretty_print_datablock(fp)
FILE *fp;
{
    int i, j, k, l, len;
    int *cwidth;

    if (cifFiles.numDatablock == 0 ||
	cifFiles.curDatablock < 0 || 
	cifFiles.curDatablock >= cifFiles.numDatablock)
	return;


    fprintf(fp, "\n\n--------------------------------------------------------\n");
    fprintf(fp, "Data Block: %s\n", cifFiles.datablocks[cifFiles.curDatablock].datablockName);


    for (i=0; i< cifFiles.datablocks[cifFiles.curDatablock].numCategory; i++) {
	cwidth = (int *) calloc(cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol, sizeof(int));
	for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) {
	    cwidth[j] = strlen(cifFiles.datablocks[cifFiles.curDatablock].categories[i].colNames[j]);
	    if (cwidth[j] < 10) cwidth[j] = 10;
	}
	for (k=0; k< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numRow; k++)  {
	    for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) { 
		if ((cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j] != NULL) &&
		     (strlen(cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j]) > 0)) {
		    if ((strlen(cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j]) > cwidth[j]) &&
                         (strlen(cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j])  > 10)) { 
			cwidth[j] = strlen(cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j]);
		    }
		}
	    }
	}
	for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) { 
	    cwidth[j] += 4;
	}
	fprintf(fp, "\n\n--------------------------------------------------------\n");	
	fprintf(fp,"Category: %s\n", cifFiles.datablocks[cifFiles.curDatablock].categories[i].categoryName);

	for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) {
	    fprintf(fp, "%s", cifFiles.datablocks[cifFiles.curDatablock].categories[i].colNames[j]);
	    len = cwidth[j] - strlen(cifFiles.datablocks[cifFiles.curDatablock].categories[i].colNames[j]);
	    for (l=0; l < len; l++)
		fprintf(fp," ");
	}
	fprintf(fp, "\n");
	for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++)  {
	    for (l=0; l < cwidth[j]-4; l++)
		fprintf(fp,"-");
	    for (l=0; l < 4; l++)
		fprintf(fp," ");
	}
	fprintf(fp, "\n");
	
	for (k=0; k< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numRow; k++)  {
	    for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) { 
		if (cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j] == NULL ||
		    !strcmp(cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j],"")) {
		    
		    fprintf(fp, "%s", "(null)");
		    len = cwidth[j] - 6;
		    for (l=0; l < len; l++)
			fprintf(fp," ");
		} else {
		    fprintf(fp, "%s", cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j]);
		    len = cwidth[j] - strlen(cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j]);
		    for (l=0; l < len; l++)
			fprintf(fp," ");
		}
	    }
	    fprintf(fp, "\n");
	}
	free(cwidth);
    }

}


void ndb_cif_print_category(fp, category)
FILE *fp;
char *category;
{
  int i, j, k;
  if (cifFiles.numDatablock == 0 ||
      cifFiles.curDatablock < 0 || 
      cifFiles.curDatablock >= cifFiles.numDatablock)
    return;
  for (i=0; i< cifFiles.datablocks[cifFiles.curDatablock].numCategory; i++)
    if (!strcmp(cifFiles.datablocks[cifFiles.curDatablock].categories[i].categoryName, category))
      break;
  if (i == cifFiles.datablocks[cifFiles.curDatablock].numCategory) return;
  fprintf(fp,"Category: %s\n", cifFiles.datablocks[cifFiles.curDatablock].categories[i].categoryName);

  for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) 
    fprintf(fp, "%s  ", cifFiles.datablocks[cifFiles.curDatablock].categories[i].colNames[j]);
  fprintf(fp, "\n");
  for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) 
    fprintf(fp, "-------");
  fprintf(fp, "\n");
  for (k=0; k< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numRow; k++)  {
    for (j=0; j< cifFiles.datablocks[cifFiles.curDatablock].categories[i].numCol; j++) { 
      fprintf(fp, "*%s* ", cifFiles.datablocks[cifFiles.curDatablock].categories[i].rows[k].columns[j]);
    }
    fprintf(fp, "\n");
  }
}


void ndb_cif_print_datablocks(fp)
FILE *fp;
{
  int i, j, k, n;
  if (cifFiles.numDatablock == 0 ||
      cifFiles.curDatablock < 0 || 
      cifFiles.curDatablock >= cifFiles.numDatablock)
    return;

  for (n=0; n < cifFiles.numDatablock; n++) {
    fprintf(fp, "Data Block: %d %s\n", n, cifFiles.datablocks[n].datablockName);
    
    for (i=0; i< cifFiles.datablocks[n].numCategory; i++) {
      fprintf(fp,"Category: %s\n", cifFiles.datablocks[n].categories[i].categoryName);
      
      for (j=0; j< cifFiles.datablocks[n].categories[i].numCol; j++) 
	fprintf(fp, "%s  ", cifFiles.datablocks[n].categories[i].colNames[j]);
      fprintf(fp, "\n");
      for (j=0; j< cifFiles.datablocks[n].categories[i].numCol; j++) 
	fprintf(fp, "-------");
      fprintf(fp, "\n");
      for (k=0; k< cifFiles.datablocks[n].categories[i].numRow; k++)  {
	for (j=0; j< cifFiles.datablocks[n].categories[i].numCol; j++) { 
	  fprintf(fp, "%s ", cifFiles.datablocks[n].categories[i].rows[k].columns[j]);
	}
	fprintf(fp, "\n");
      }
    }
  }
}  

int get_column_index(blockIndex, categoryIndex, columnName)
int blockIndex;
int categoryIndex;
char *columnName;
  /*
   *  Return the index of for 'columnName' in the data block and category 
   *  specified by 'blockIndex' and 'categoryIndex'.
   */
{
    int icol, k;
    icol = -1;
    for (k=0; k < cifFiles.datablocks[blockIndex].categories[categoryIndex].numCol; k++) {
        if (!strcmp(columnName,cifFiles.datablocks[blockIndex].categories[categoryIndex].colNames[k])) {
            icol = k;
            break;
        }
    }
    return icol;
}

int get_category_index(blockIndex, categoryName)
int blockIndex;
char *categoryName;
  /*
   *  Return the index of for 'categoryName' in the data block 
   *  specified by 'blockIndex'.
   */
{
    int icat, i;
    icat = -1;
    for (i=0; i< cifFiles.datablocks[blockIndex].numCategory; i++) {
        if (!strcmp(categoryName,cifFiles.datablocks[blockIndex].categories[i].categoryName)) {
            icat = i;
            break;
        }
    }
    return(icat);
}
