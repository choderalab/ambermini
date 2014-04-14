%{ 
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
  FILE:         ndb_cifparse_lex.l

  This file is part of the NDBQUERY application,
  a program component of the Nucleic Acid Database Project.

  H. M. Berman, W. K. Olson, D. L. Beveridge, J. K. Westbrook, A. Gelbin,
  T. Demeny, S. H. Shieh, A. R. Srinivasan, and B. Schneider. (1992).
  The Nucleic Acid Database: A Comprehensive Relational Database of 
  Three-Dimensional Structures of Nucleic Acids.  Biophys J., 63, 751-
  759.

  Questions about this software should be directed to:

              ndbadmin@ndbserver.rutgers.edu

  PURPOSE:    DDL 2.1 compliant CIF file lexer ...

  Modified to allow multiple lexers in single program, support flex, etc.
    by D.A. Case,  12/97

*/

/*$$LICENSE-NDB$$*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cifparse.h"
#include "cifp.tab.h"

/* 
 *  Local prototypes
 */
static void alt_yymore();
static void delete_punctuation();
static char tbuf[MAXVALUELENGTH+1];
/*
 *  Local static variables ...
 */

static int   isText = FALSE;
static int   i, j;
static int   buflen;

%}
%%
[\n];|[\n]    {                                       /* New Lines */
          lineNo++;
          if (yytext[yyleng-1] == ';') {
             if (isText == TRUE) {          /* end of text value */
   	        strncat(tbuf,yytext,MAXVALUELENGTH);
		buflen = strlen(tbuf);
	        delete_punctuation(yylval.TempBuffer,tbuf,&buflen);
                isText = FALSE;
                return(VALUE);
	     } else {  /* text value begins */
	       isText = TRUE;
	       strcpy(tbuf,"");
	     }
	  }
          if (isText == TRUE) alt_yymore();
       }

^;|^;([^;])+      {                      /* multi-line text  value */    
      if (yytext[0] == ';' && yytext[yyleng-1] != '\n') {  /* text value begins*/
	 isText = TRUE;
	 strncpy(tbuf,yytext,MAXVALUELENGTH);
      }
      else {  
          isText = FALSE;
	  strncat(tbuf,yytext,MAXVALUELENGTH);
	  buflen = strlen(tbuf);
          for (i=0; i<buflen; i++)
	       if (tbuf[i] == '\n') lineNo++;

	  delete_punctuation(yylval.TempBuffer,tbuf,&buflen); 
          return(VALUE);
        } 
    }

[\t ]+      {                                    /* white space */ 
      if (isText) alt_yymore();
       else {
	 for (i=0; i<yyleng; i++)
	   if (yytext[i] == '\n') lineNo++;
       }
    }

[Dd][Aa][Tt][Aa][_][^\n\t ]*   {                 /* data_ */ 
      if (isText) alt_yymore();
      else {
        strncpy(yylval.TempBuffer,yytext,MAXVALUELENGTH);
        return (DATABLOCK);
      }
    }

[Ll][Oo][Oo][Pp][_]      {                        /* loop_ */
      if (isText) alt_yymore();
      else return (LOOP);
    }

[?]|[.]  {                                        /* Unknown */
      if (isText) alt_yymore();
      else return (UNKNOWN);
    }

[#][^\n]*[\n]       {                             /* COMMENT */
      if (isText) alt_yymore();
      else lineNo++;
    }

[_][^\n \t]+  {                /* Item names start with _ end with spaces*/
      if (isText) {
         alt_yymore();
      } else {
        /* If the beginning of text is in buffer yytext */
         for (i=0; i<yyleng; i++) 
           if (yytext[i] == '_') break;
	 strncpy(yylval.TempBuffer,yytext,MAXVALUELENGTH);
	 return(ITEMNAME);
      }
   }


[^\n\t _][^\n \t]*    {  /* unquoted or mixed quote strings strings */
     if (!isText) {
        j=0;
        for (i=yyleng-1; i >= 0; i--) {
            if ( yytext[i] == '\'' || yytext[i] == '\"') {
               yytext[i]='\0';
               break;
            } else 
             break;
        }
        for (i=0; i < yyleng; i++) {
             if (yytext[i] == '\'' || yytext[i] == '\"') {
                j++;
                break;
            } else
               break;
        }
        strncpy(yylval.TempBuffer, &yytext[j],MAXVALUELENGTH);
	return(VALUE);
     }
     else {
       alt_yymore();
     }
   }


[\'][^'\n]*[\'\n][ \t\n]   {  /* Single Quoted Strings */

     if (!isText) {

        j=0;
        for (i=yyleng-1; i >= 0; i--) {
            if ( yytext[i] == ' ' || yytext[i] == '\t' ||  yytext[i] == '\n')
                 yytext[i]='\0';
            else if ( yytext[i] == '\'') {
               yytext[i]='\0';
               break;
            } else 
             break;
        }
        for (i=0; i < yyleng; i++) {
            if ( yytext[i] == ' ' || yytext[i] == '\t') 
                j++;
             else if (yytext[i] == '\'') {
                j++;
                break;
            } else
               break;
        }
        strncpy(yylval.TempBuffer, &yytext[j],MAXVALUELENGTH);
	return(VALUE);
     }
     else {

        alt_yymore();
        if (yytext[yyleng-1] == '\n') unput('\n');
     }
   }


[\"][^"\n]*[\"\n][ \t\n]    {  /* Double Quoted Strings */

     if (!isText) {
        j=0;
        for (i=yyleng-1; i >= 0; i--) {
            if ( yytext[i] == ' ' || yytext[i] == '\t' ||  yytext[i] == '\n')
                 yytext[i]='\0';
            else if ( yytext[i] == '\"') {
               yytext[i]='\0';
               break;
            } else 
             break;
        }
        for (i=0; i < yyleng; i++) {
            if ( yytext[i] == ' ' || yytext[i] == '\t') 
                j++;
             else if (yytext[i] == '\"') {
                j++;
                break;
            } else
               break;
        }
        strncpy(yylval.TempBuffer, &yytext[j], MAXVALUELENGTH);
	return(VALUE);
     }
     else {
       alt_yymore();
        if (yytext[yyleng-1] == '\n') unput('\n');
     }
   }

%%
/*---------------------------------------------------------------------------
 *                     delete_punctuation
 * Purpose:                                                                  
 *         Deletes leading and trailing semicolons and quotes on text strings.
 * Input :                                                                   
 *         char *ibuf:   The string which yylex finds                  
 *         char *obuf:   unquoted string
 *         int  *bufLen: The sting length of ibuf                   
 *---------------------------------------------------------------------------*/
static void delete_punctuation(obuf, ibuf, bufLen)
char *obuf;
char *ibuf;
int *bufLen;
{
  int i;
  if (*bufLen == 0) return;
  if (*bufLen > MAXVALUELENGTH) {
    fprintf( stderr,"String value exceeds maximum length %d > %d",*bufLen,MAXVALUELENGTH);
    fprintf( stderr,"You must reset the CIFPARSE definition of MAXVALUELENGTH greater than %d to process this file",*bufLen);
    exit(1);
  }

  for (i=0; i<*bufLen; i++)
    if (ibuf[i] == ';') break;
  strcpy(obuf, &ibuf[i+1]);
  (*bufLen)-=i;
  
  for (i= *bufLen-1; i >=0; i--)
    if (obuf[i] == ';') break;
  if (obuf[i-1] == '\n') i--;			   
  (*bufLen)-= (*bufLen-i);
  obuf[*bufLen] = '\0';
}


int yywrap()
{
   return(1);
}

static void alt_yymore()
{
   buflen = strlen(tbuf);
   if (buflen + strlen(yytext) > MAXVALUELENGTH) {
     fprintf( stderr,"Exceeded maximum buffer size %d starting at line %d",MAXVALUELENGTH,lineNo);
     fprintf( stderr,"Increase the MAXVALUELENGTH to read this file");
   } else 
     strcat(tbuf, yytext);
}

