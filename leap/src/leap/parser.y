/*
 *      File:   parser.y
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
 *             David A. Rivkin                                          *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      A YACC program for parsing commands to LEaP.
 *
 *      The syntax is described in this file.
 *      Comments are handled in the routine cGetChar, which skips over
 *      all text between '#' and '\n' inclusive.
 *
 *      The syntax has the form:        [] delimit necessary syntax elements
 *
 *      [command] [arg1], [arg2], [arg3], ... ;
 *      [variable] = [expression];
 *      [variable] = ;                           Clears a variable 
 *      
 *
 *      The parser requires all commands to be ended with a semicolon ';'
 *      Comments are started with a '#' and terminated with a '\n'.
 *      Strings can be delimited with double quotes '"', or
 *      started with a single '$' character and ended with a comma, space,
 *      curly bracket '}', semicolon, equals sign, etc.
 *
        
 */


%token  LVARIABLE
%token  LSTRING
%token  LNUMBER
%token  LASSIGN 
%token  LENDOFCOMMAND
%token  LOPENLIST
%token  LCLOSELIST
%token  LOPENPAREN
%token  LCLOSEPAREN
%token  LQUIT
%token  LCOMMA
%token  LDOT
%token  LCOMMAND
%token  LDUMMY
%token  LNULL
%token  LNOTSINGLECHAR

%{
#include	<unistd.h>
#include        "basics.h"

#include        "classes.h"

#include        "dictionary.h"
#include        "parmLib.h"

#include        "commands.h"
#include	"block.h"
#include	"parser.h"

#include        "leap.h"
#include        "block.h"

#include        "help.h"

#define         MESSAGEFILTER   MESSPARSER




#define         NULLSTR         "null"

/*
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        GLOBAL VARIABLES

*/

                /* The global variable GfCurrentInput defines the file */
                /* from which input is currently read.  When the file */
                /* is empty, then switch back to stdin */

#define	MAXINPUT	1000
#define	MAXINPUTFILES	10		/* Maximum 10 input files */
					/* can be open at once */


char            GsInputLine[MAXINPUT] = "";
BOOL		GbLastLine = FALSE;
BOOL		bCmdDeleteObj;
int             GiInputPos = 0;
PARMLIB		GplAllParameters;
RESULTt		GrMainResult;
BLOCK		GbCommand = NULL;
BLOCK		GbExecute = NULL;
int		GiClipPrompts = 0;
BOOL		GbGraphicalEnvironment;
STRING		GsProgramName;

extern int	iMemDebug;

static	STRING	*SbFirstSourceFiles = NULL;
static	int	iFirstSource = 0;
static	BOOL	SbUseStartup = TRUE;



/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *
 *	The following is used by the parser to maintain a stack of
 *	files where input is received from.  If the file is NULL
 *	then input is received from the main program in the form
 *	of BLOCKS.  The main program is then responsible for reading
 *	the stdin ( in the command line interface ) or for gathering
 *	keypress events from X-Windows ( in the graphical interface ).
 *
 */


int             GiInputFileStackPos = 0;
FILE*           GfaInputFileStack[MAXINPUTFILES];



/*
 *----------------------------------------------------------------
 *
 *	Not quite GLOBAL variables used by the parser.
 */

                                /* Arguments to routines are passed through */
                                /* an array */
#define MAXARGS         10
#define MAXLISTNEXT     10


ATOM            aDummy;
ASSOC           aaArgs[MAXARGS];
int             iArgCount, i;

                                /* List stuff is used for input of nested */
                                /* lists */
#define MAXLISTNEST     10
ASSOC           aaLists[MAXLISTNEST];
int             iCurrentList = -1;
#define PUSHLIST()      iCurrentList++
#define POPLIST()       iCurrentList--
#define CURRENTLIST     aaLists[iCurrentList]


OBJEKT          o0;
double          dTemp;
ASSOC           aAssoc;
STRING          sTemp;
BOOL		bQuit = FALSE;
BOOL		bCommandFound = FALSE;

                /* There seems to be a problem with YACC not properly */
                /* declaring yylval and yyval */
typedef struct  {
	ASSOC		aVal;
	double		dVal;
	STRING		sVal;
	FUNCTION	fCallback;
} YYSTYPEt;

#define YYSTYPE YYSTYPEt


extern  OBJEKT  oGetObject();           /* ( STRING ) */
extern  int     yyparse();
%}


%start  input





%%
/*------------------------------------------------------------

        RULES
*/

/*
 *	Bogus function for xaUtilMessageFilter
 *
void
yyparse()
{
*/

input   :       line
		{
			return 0;
		}
		;

line    :       LENDOFCOMMAND
	|	instruct
                        {
                        bCommandFound = FALSE;
                        }
        |       error LENDOFCOMMAND
                        {
                            VP0(( "\n" ));
                            yyerrok;
                        }
        ;

instruct:       assign LENDOFCOMMAND
        |       function LENDOFCOMMAND
        ;

assign  :       LVARIABLE LASSIGN express 
                        {
                                /* Set the value of the variable */
                            aAssoc = $<aVal>3;
			    if ( aAssoc != NULL ) {
                                MESSAGE(( "Assigning a value to %s\n", 
                                                        $<sVal>1 ));
                                VariableSet( $<sVal>1, oAssocObject(aAssoc) );
				MESSAGE(( "DEREF (assign) - %s\n",
							sAssocName(aAssoc) ));
                                DEREF( aAssoc );
			    } else {
				MESSAGE(("Not assigning value to %s - rmving\n",
							$<sVal>1 ));
				VariableRemove( $<sVal>1 );
			    }

                        }
        |       LVARIABLE LASSIGN 
                        {
                            MESSAGE(( "Removing variable %s\n", $<sVal>1 ));
                            VariableRemove( $<sVal>1 );
                        }
        ;

express :       rawexp
        |       function
        ;

rawexp  :       LOPENLIST 
                        {
                            PUSHLIST();
                                /* Create an ASSOC for the list */
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, oCreate(LISTid) );
                            CURRENTLIST = aAssoc;
                        }
                    elements LCLOSELIST 
                        {
                            $<aVal>$ = $<aVal>3;
                            POPLIST();
                        }
        |       LNUMBER
                        {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            o0 = oCreate(ODOUBLEid);
                            ODoubleSet( o0, $<dVal>1 );
                            AssocSetObject( aAssoc, o0 );
                            $<aVal>$ = aAssoc;
                        }
        |       LSTRING
                        {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            o0 = oCreate(OSTRINGid);
                            OStringDefine( (OSTRING) o0, $<sVal>1 );
                            AssocSetObject( aAssoc, o0 );
			    DEREF( o0 );	/* keeps count = 1 */
                            $<aVal>$ = aAssoc;
                        }
        |       LVARIABLE
                        {
			    OBJEKT	o = oGetObject($<sVal>1);
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, $<sVal>1 );
                            AssocSetObject( aAssoc, o ); /* REF's o */
                            $<aVal>$ = aAssoc;
                        }
        |       LDUMMY
                        {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, aDummy );
                            $<aVal>$ = aAssoc;
                        }
        |       LNULL
                        {
                            MESSAGE(( "Parsed a null\n" ));
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, NULL );
                            $<aVal>$ = aAssoc;
                        }
        ;


function:       cmdname 
                        { iArgCount = 0;
                        } 
                    args  
                        {
                                /* Execute the command */
			    MESSAGE(( "executing function\n"));
			    bCmdDeleteObj = FALSE;
                            o0 = $<fCallback>1( iArgCount, aaArgs );
			    if ( o0 != NULL ) {
                            	aAssoc = (ASSOC)oCreate(ASSOCid);
                            	AssocSetObject( aAssoc, o0 );
                            	$<aVal>$ = aAssoc;
			    } else {
				MESSAGE(( "func == NULL---\n"));
				$<aVal>$ = NULL;
			    }

                                /* DEREF each of the arguments */

                            for ( i=0; i<iArgCount; i++ ) {
				if ( bCmdDeleteObj ) {
					MESSAGE(( "bCmdDeleteObj---\n" ));
					DEREF( oAssocObject( aaArgs[i] ) );
				}
				MESSAGE(( "DEREF (function) - %s\n",
						sAssocName(aaArgs[i])));
                                DEREF( aaArgs[i] );
                            }
                        }
        ;


elements:       /* nothing */
        |       elements rawexp 
                        {
                                /* Get the element and add it to the list */
                            MESSAGE(( "Adding to list:\n" ));
                            ListAddToEnd( (LIST)oAssocObject(CURRENTLIST), 
							(OBJEKT) $<aVal>2 );
                            $<aVal>$ = CURRENTLIST;
                        } 
        ;

cmdname :	LCOMMAND
        ;

args    :       /* Nothing */
        |       arg
        |       args arg
        ;

arg     :       rawexp  
                        {
                            aaArgs[iArgCount++] = $<aVal>1;
                        }
        ;



/*
 *
 *	Bogus stuff for xaUtilMessageFilter
 *
}
 */



%%
/*------------------------------------------------------------

        ROUTINES

*/

static  BOOL    SbGotUngetc = FALSE;
static  char    ScUngetc;


/*
 *      yyerror
 *
 *      Respond to errors.
 */
int
yyerror( char *sStr )
{
    VP0(( "ERROR: %s\n", sStr ));
    return 1;
}

FILE *
fINPUTFILE()            
{
	if ( GiInputFileStackPos < 0 )
		return(NULL);
	return( GfaInputFileStack[GiInputFileStackPos] );
}



/*
 *	zbGetLine
 *
 *	Get the next line of input from the current input source.
 *	This may be GbCommand or GbExecute depending if
 *	the input is coming from the user or from an execute file.
 *	Return FALSE if there is no more lines to be received from 
 *	the BLOCK.
 */
BOOL
zbGetLine( char *sLine, BOOL *bPFromExecute )
{
BOOL		bGotBlock;
char		c;

    if ( fINPUTFILE() == NULL ) {
	*bPFromExecute = FALSE;
	return(bBlockReadLine( GbCommand, sLine ));
    }

    *bPFromExecute = TRUE;

		/* If there is no line in the execute BLOCK */
		/* then read in another block from the current file */

    if ( bBlockEndOfRead(GbExecute) ) {
	    BlockEmpty( GbExecute );
	    bGotBlock = FALSE;
	    while ( !feof(fINPUTFILE()) && !bGotBlock ) {
		c = fgetc(fINPUTFILE());
		if ( feof(fINPUTFILE()) ) break;
		bGotBlock = bBlockAddChar( GbExecute, c );
	    }

		/* If a complete BLOCK was not read then */
		/* append a final '\n' character, if that doesn't */
		/* make this a complete BLOCK then there is an error */
		/* in the input, also we are at the end of the file */
		/* so pop the file off the stack and say that we have */
		/* a complete BLOCK */

	    if ( !bGotBlock ) {
	    	bGotBlock = bBlockAddChar( GbExecute, '\n' );
	    	if ( fINPUTFILE() != NULL )
	    		fclose( fINPUTFILE() );
	    	INPUTPOPFILE();
	    }
    }
    return(bBlockReadLine( GbExecute, sLine ));
}


		
	


/*
 *      zcGetChar
 *
 *      Get the next character from the line buffer.
 *	If there are no more characters in the line buffer and
 *	GbLastLine is TRUE then return '\0', otherwise if the
 *	line buffer is empty, fill it and return the next character
 *	in it.
 */
char
zcGetChar()
{
char            c;
BOOL		bFromExecute;

                /* If there is a pushed character then return it */

    if ( SbGotUngetc ) {
        SbGotUngetc = FALSE;
        c = ScUngetc;
        goto DONE;
    }

		/* Now if the input line is empty then fill it */
    if ( GsInputLine[GiInputPos] == '\0' ) {
	if ( GbLastLine ) {
	    c = '\0';
	    GbLastLine = FALSE;
	    goto DONE;
	}
	GbLastLine = zbGetLine( GsInputLine, &bFromExecute );
	GiInputPos = 0;
	if ( bFromExecute ) {
	    for ( i=GiClipPrompts; i<=iINPUTSTACKDEPTH(); i++ ) VP2(( ">" ));
	    VP2(( " " ));
	    VP2(( "%s", GsInputLine ));
	} else {
	    VPLOG(( "> %s", GsInputLine ));
	}
    }

    c = GsInputLine[GiInputPos++];

DONE:
    return(c);

}



/*
 *      zUngetc
 *
 *      Push one character back to the file.
 */
void
zUngetc( char c )
{
    SbGotUngetc = TRUE;
    ScUngetc = c;
}






/*
 *      cGetChar
 *
 *      Return the next character, skipping over comments
 *      which start with '#' and end with '\n'.
 */
char
cGetChar()
{
char    c;

	while (1) {
        	c = zcGetChar();
        	if ( c == '#' )
            		while ( ( c=zcGetChar() ) != '\n' )  /* Nothing */ ;
		else
			break;
        } 
	return(c);
}

 
int
iOneCharToken( int c )
{
	switch(c) {
		case ';':
        		return(LENDOFCOMMAND);
		case '=':
        		return(LASSIGN);
		case '*':
        		return(LDUMMY);
		case '(':
        		return(LOPENPAREN);
		case ')':
        		return(LCLOSEPAREN);
		case '{':
        		return(LOPENLIST);
		case '}':
        		return(LCLOSELIST);
		default:
        		return(LNOTSINGLECHAR);
	}
}

/*
 *      yylex
 *
 *      Lexical analyzer for the parser.
 *      Read characters from stdin and return the token types
 *      read and place the value read into the global UNION
 *      yylval.
 *
 *      The things that it recognizes are:
 *              LDOUBLE         [-]###.###E## or ###
 *              LSTRING         "xx xx xxx" or '$'everything up to ' ' ',' ';'
 *              commands        xxxxxxx
 *              LVARIABLE       xxxxxxx which are not commands
 *              LTERMINATOR     ;
 *              LASSIGN         =
 *
 *	Modified 17 November 1992 - David A. Rivkin
 *		Added checking the alias table for command matches.
 *	Total rewrite October 1993 - Bill Ross
 *
 */
int
yylex()
{
STRING          sStr;
int             j, iMax, tok;
BOOL            bGotExp, bGotDot;
char            c;
STRING		sCmd;
STRING		sPossibleCmd;

                /* Skip over blanks, tabs, end of lines etc */

    while ( (c=cGetChar())==' '  ||  c=='\t'  ||  c=='\n'  ||  c == ',' );

    if ( c == '\0' ) 
	return(LENDOFCOMMAND);

    /*
     *  Check the 1-character possibilities: , ; = * ( ) { }
     */
    tok = iOneCharToken( c );
    if ( tok != LNOTSINGLECHAR ) {
        MESSAGE(( "Parsed /%c/\n", c ));
	return(tok);
    }

    /*
     *  it isn't a 1-char thing; read in the rest 
     *	and push back the terminating char
     */
    sStr[0] = c;
    for (j=1;;j++) {

	if ( j >= sizeof(STRING) )
	    DFATAL(( "string too long" ));

	c = cGetChar();
	/*
	 *  NULL terminates anything (?)
	 */
	if ( c == '\0' ) {
	    	sStr[j] = '\0';
	    	break;
	}
	/*
	 *  allow anything inside quotes; chop closing quote
	 */
	if ( sStr[0] == '"' ) {
		if ( c == '"' ) {
			sStr[j] = '\0';
			break;
		}
		sStr[j] = c;
		continue;
	}
	/*
	 *  whitespace is a delimiter outside of quotes
	 */
	if ( c == ' '  ||  c == '\t'  ||  c == '\n'  ||  c == ',' ) {
	    	sStr[j] = '\0';
	    	break;
	}
	/*
	 *  special case for $-type strings: allow embedded single-char
	 *	tokens, except ';'
	 */
	if ( sStr[0] == '$' ) {
		if ( c == ';' ) {
			zUngetc( c );
	    		sStr[j] = '\0';
	    		break;
		}
		sStr[j] = c;
		continue;
	}
	tok = iOneCharToken( c );
	if ( tok != LNOTSINGLECHAR ) {
		zUngetc( c );
	    	sStr[j] = '\0';
	    	break;
	}
	sStr[j] = c;
    }

    /*
     *  see if it's a number
     */
    bGotExp = FALSE;
    bGotDot = FALSE;
    if ( isdigit(sStr[0]) || sStr[0] == '-' || sStr[0] == '+' || 
				( sStr[0] == '.' && isdigit(sStr[1]) ) ) {

        for ( j=0; j<sizeof(STRING); j++ ) {
            MESSAGE(( "Thinking NUMBER got: %c\n", sStr[j] ));
	    switch ( sStr[j] ) {
		case '\0':
        	    if ( sscanf( sStr, "%lf", &yylval.dVal ) != 1 ) {
			VP0(( " Couldn't scan NUMBER from (%s)\n", sStr ));
			return(LNULL);
		    }
        	    MESSAGE(( "Parsed a number: %lf\n", yylval.dVal ));
        	    return(LNUMBER);
		case '.':
		    if ( bGotDot ) {
			VP0(( "(Multiple '.' in NUMBER-like thing (%s))\n", 
								sStr ));
			goto notnum;
		    }
		    if ( bGotExp ) {
			VP0(( 
			 "('.' follows exponent in NUMBER-like thing (%s))\n",
								sStr ));
			goto notnum;
		    }
        	    bGotDot = TRUE;
		    break;
		case 'e':
		case 'E':
            	    if ( bGotExp ) {
			VP0(( "(Multiple 'e' in NUMBER-like thing (%s))\n", 
								sStr ));
			goto notnum;
		    }
                    bGotExp = TRUE;
                    break;
		case '+':
		case '-':
                    break;
		default:
            	    if ( !isdigit(sStr[j]) ) {
			goto notnum;
		    }
            	    break;
	    }
	}
    }
notnum:
    /* 
     *  see if it's a string in quotes
     */
    if ( sStr[0] == '"' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE(( "Parsed a STRING: %s\n", sStr ));
        return(LSTRING);
    }

    /* 
     *  see if it's a string prefixed w/ '$'
     */
    if ( sStr[0] == '$' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE(( "Parsed a STRING: %s\n", sStr ));
        return(LSTRING);
    }

                /* LASTLY!!!!!!!! */
    /* 
     *  see if it's a variable/command 
     */
    strcpy( yylval.sVal, sStr );
    strcpy( sPossibleCmd, sStr );
    StringLower( sPossibleCmd );

    		/* Check if there is an alias that is an exact match */
    if ( (iMax = iVarArrayElementCount( GvaAlias )) ) {
	ALIAS		aAlias;
	aAlias = PVAI( GvaAlias, ALIASt, 0 );
	for ( i=0; i<iMax; i++, aAlias++ ) {
	    if ( strcmp( aAlias->sName, sPossibleCmd ) == 0 ) {
	    	strcpy( sPossibleCmd, aAlias->sCommand );
	    }
        }
    }
                /* Check if there is an exact match of the command */
                /* If a command has already been found for this input
                	line, then do not consider the string a command
                	but rather as a STRING variable */
                	
    if ( !bCommandFound ) {
	for ( j=0; strlen(cCommands[j].sName) != 0; j++ ) {
	    strcpy( sCmd, cCommands[j].sName );
	    StringLower( sCmd );
            if ( strcmp( sCmd, sPossibleCmd ) == 0 ) {
		yylval.fCallback = cCommands[j].fCallback;
		MESSAGE(( "Parsed a command: %s\n", sStr ));
		bCommandFound = TRUE;
		return(LCOMMAND);
	    }
        }
    }


                /* If the variable name is null then return LNULL */

    if ( strcmp( sStr, NULLSTR ) == 0 ) 
	return(LNULL);
    
                /* Return the variable name */

    strcpy( yylval.sVal, sStr );
    MESSAGE(( "Parsed a variable: %s\n", sStr ));
    return(LVARIABLE);

}




/*
 *      oGetObject
 *
 *      If the string is a variable then return the OBJEKT that
 *      it is attached to, otherwise if the string is
 *      a string like: 'unit.mol.res.atom' then parse the
 *      individual names and search the CONTAINERS for
 *      the subcontainers.
 */
OBJEKT
oGetObject( char *sName )
{
CONTAINER       cCont[5];
int             j, k, iSeq;
OBJEKT          oObj;
STRING          sLine, sHead;
BOOL		bDot, bAt, bPdbSeq;
STRING		sGroup;
LIST		lGroup;
OSTRING		osString;
LOOP		lRes;
RESIDUE		rRes;

    oObj = oVariable( sName );
    if ( oObj != NULL ) return(oObj);

        /* Now try to parse the name */
    strcpy( sLine, sName );
    bDot = FALSE;
    bAt = FALSE;
    bPdbSeq = FALSE;
    for ( k=0; k<strlen(sName); k++ ) {
	if ( sLine[k] == '.' ) {
	    sLine[k]=' ';
	    bDot = TRUE;
/*fprintf(stderr, "GOTDOT\n"); */
	} else if ( sLine[k] == '@' ) {
	    sLine[k] = ' ';
	    bAt = TRUE;
	} else if ( sLine[k] == '%' ) {
	    if ( !bDot ) {
	        sLine[k] = ' ';
	        bPdbSeq = TRUE;
	    }
	}
    }


    sRemoveFirstString( sLine, sHead );
    cCont[0] = (CONTAINER)oVariable(sHead);
    if ( cCont[0] == NULL ) {
/* fprintf(stderr, "STRING %s\n", sName); */
    	/* It is not an object variable so...
    	   return the whole string as a OSTRING */
	goto String;
    }
    
/*fprintf(stderr, "VARIABLE %c %c\n", bAt, bPdbSeq);*/
    if ( bPdbSeq ) {
	sRemoveLeadingSpaces( sLine );
	sRemoveFirstString( sLine, sHead );
	if ( strlen(sHead) == 0 ) {
	    goto String;
	}
	if ( isdigit(sHead[0]) ) {

			/* Make sure the rest are digits */
			/* If not return NULL */
	    for ( j=1; j<strlen(sHead); j++ ) {
		if ( !isdigit(sHead[j]) ) {
			    /* It is not an object variable so...
			       return the whole string as a OSTRING */
		    goto String;
		}
	    }
	
			/* Find the PDB sequence number */

	    cCont[1] = NULL;	
	    iSeq = atoi(sHead);
	    lRes = lLoop( (OBJEKT)cCont[0], RESIDUES );
	    while ( (rRes = (RESIDUE)oNext(&lRes)) ) {
		if ( iResiduePdbSequence(rRes)==iSeq ) {
		    cCont[1] = (CONTAINER) rRes;
		}
	    }

	    if ( !cCont[1] ) {
		goto String;
	    }

	    if ( bDot ) {
		sRemoveLeadingSpaces( sLine );
		sRemoveFirstString( sLine, sHead );
		if ( strlen(sHead) == 0 ) {
		    goto String;
		}
		if ( isdigit(sHead[0]) ) {

				/* Make sure the rest are digits */
				/* If not return NULL */
		    for ( j=1; j<strlen(sHead); j++ ) {
			if ( !isdigit(sHead[j]) ) {
				    /* It is not an object variable so...
				       return the whole string as a OSTRING */
			    goto String;
			}
		    }
			
		    iSeq = atoi(sHead);
		    cCont[2] = cContainerFindSequence( cCont[1], 
						DIRECTCONTENTS, iSeq );
		} else {
		    cCont[2] = 
			cContainerFindName( cCont[1], DIRECTCONTENTS, sHead );
		}
		return((OBJEKT)cCont[2]);
	    } else {
		return((OBJEKT)cCont[1]);
	    }
	} else {
	    goto String;
	}
    }


    if ( bDot ) {
	k = 0;
	do {
	    sRemoveLeadingSpaces( sLine );
	    sRemoveFirstString( sLine, sHead );
	    if ( strlen(sHead) == 0 ) break;
	    k++;
	    if ( cCont[k-1]->oHeader.cObjType == PARMSETid ) {
		/*  semi-HACK - parmsets don't have contents in this sense */
		cCont[k] = NULL;
	    } else if ( isdigit(sHead[0]) ) {

			    /* Make sure the rest are digits */
			    /* If not return NULL */
		for ( j=1; j<strlen(sHead); j++ ) {
		    if ( !isdigit(sHead[j]) ) {
   				/* It is not an object variable so...
			    	   return the whole string as a OSTRING */
			goto String;
    		    }
		}
		iSeq = atoi(sHead);
		cCont[k] = cContainerFindSequence( cCont[k-1], 
					    DIRECTCONTENTS, iSeq );
	    } else {
		cCont[k] = 
			cContainerFindName( cCont[k-1], DIRECTCONTENTS, sHead );
	    }
	    if ( cCont[k] == NULL ) break;
	} while ( strlen(sHead) != 0 ) ;
        if ( cCont[k] == NULL ) {
    		/* It is not an object variable so...
    		   return the whole string as a OSTRING */
		goto String;
    	}
	return((OBJEKT)cCont[k]);
    }

			/* If group notation then return the group */

    if ( bAt ) {
	sRemoveLeadingSpaces( sLine );
	sRemoveFirstString( sLine, sGroup );
	if ( iObjectType(cCont[0]) != UNITid ) {
    		/* It is not an object variable so...
    		   return the whole string as a OSTRING */
		goto String;
	}

	lGroup = lUnitGroup( (UNIT)cCont[0], sGroup );
	return((OBJEKT)lGroup);
    }

String:

   /* 
    *  It is not an object variable so...
    *	   return the whole string as a OSTRING
    *	   -- need to set refs to 0 so that
    *	      it will be freed later.. HACK
    */

    osString = (OSTRING)oCreate(OSTRINGid);
    OStringDefine( osString, sName );
    ((OBJEKT)(osString))->iReferences = 0;
    return((OBJEKT)osString );
}





    
/*
 *================================================================
 *
 *	Public routines
 */




/*
 *	ParseArguments
 *
 *	Parse the arguments that the user passes to LEaP from
 *	the command line arguments.
 */
void
ParseArguments( int argc, char *argv[] )
{
char		c;
extern	char	*optarg;

    while ( (c = getopt( argc, argv, "hsI:f:" )) != (char)(EOF) ) {
	switch (c) {
	    case 'h':
		printf( "Usage: %s [options]\n", argv[0] );
		printf( "Options:\n" );
		printf( " -h         Generate this message.\n" );
		printf( " -s         Ignore %s startup file.\n", LEAPRC );
		printf( " -I {dir}   Add {dir} to search path.\n" );
		printf( " -f {file}  Source {file}.\n" );
		exit(1);
	    case 's':
		printf( "-s: Ignoring startup file: %s\n", LEAPRC );
		SbUseStartup = FALSE;
		break;
	    case 'I':
		printf( "-I: Adding %s to search path.\n", optarg );
		BasicsAddDirectory( optarg, 1 );
		break;
	    case 'f':
		printf( "-f: Source %s.\n", optarg );
		if ( iFirstSource == 0 ) {
			MALLOC( SbFirstSourceFiles, STRING *, sizeof(STRING) );
			iFirstSource = 1;
		} else {
			iFirstSource++;
			REALLOC( SbFirstSourceFiles, STRING *, SbFirstSourceFiles,
					iFirstSource * sizeof(STRING));
		}
		strcpy( SbFirstSourceFiles[iFirstSource-1], optarg );
		break;
	}
    }
}






/*
 *	ParseInit
 *
 *	Initialize the parser.
 *	If SbStartup is TRUE the execute the LEAPRC script.
 */
void
ParseInit( RESULTt *rPResult )
{
FILE	*fStartup;
int	iFile;

    VP0(( "\nWelcome to LEaP!\n" ));

#ifdef  DEBUG
    VP0(( "LEaP is running in DEBUG mode!\n" ));
#endif

		/* Initialize memory manager debugging */
    INITMEMORYDEBUG();

    HelpInitialize();

                /* Initialize the first file in the stack to be stdin */
                
    GfaInputFileStack[0] = NULL;
    GbExecute = bBlockCreate();

                /* Create a few OBJEKTs that will be used by the parser */

    aDummy = (ATOM)oCreate(ATOMid);
    ContainerSetName( aDummy, "DUMMY" );
    GplAllParameters = plParmLibCreate();

    VariablesInit();
    GrMainResult.iCommand = CNONE;
    rPResult->iCommand = CNONE;
    
                /* Parse the LEAPRC file if bUseStartup is TRUE */

    if ( SbUseStartup ) {
        fStartup = FOPENNOCOMPLAIN( LEAPRC, "r" );
        if ( fStartup == NULL ) {
	    VP0(( "(no %s in search path)\n", LEAPRC ));
	} else {
	    /*
	     *  source the leaprc
	     */
	    VP0(( "Sourcing %s: %s\n", LEAPRC, GsBasicsFullName ));
	    INPUTPUSHFILE( fStartup );
	    GiClipPrompts = 1;
	    while ( fINPUTFILE() != NULL ) {
	        yyparse();
		if ( GrMainResult.iCommand == CQUIT ) {
	    	    if ( fINPUTFILE() != NULL )
	    	    	fclose( fINPUTFILE() );
	    	    INPUTPOPFILE();
		}
	    }
	    *rPResult = GrMainResult;
	    GiClipPrompts = 0;
	}
    }

		/* Parse the first source file specified */
		/* on the command line using the -f option */

    for (iFile=0; iFile<iFirstSource; iFile++) {
	fStartup = FOPENCOMPLAIN( SbFirstSourceFiles[iFile], "r" );
	if ( fStartup != NULL ) {
	    VP0(( "Sourcing: %s\n", GsBasicsFullName ));
	    INPUTPUSHFILE( fStartup );
	    GiClipPrompts = 1;
	    while ( fINPUTFILE() != NULL ) {
		yyparse();
		if ( GrMainResult.iCommand == CQUIT ) {
	    	    if ( fINPUTFILE() != NULL )
	    	    	fclose( fINPUTFILE() );
	    	    INPUTPOPFILE();
		}
	    }
	    *rPResult = GrMainResult;
	    GiClipPrompts = 0;
	}
    }
    if ( SbFirstSourceFiles )
	FREE( SbFirstSourceFiles );
}




/*
 *	ParseBlock
 *
 *	Parse a BLOCK containing one complete command.
 *	Return in rResult the result of the command.
 */
void
ParseBlock( BLOCK bBlock, RESULTt *rPResult )
{
		/* Set up the BLOCK from which to read the command */

    MESSAGE(( "Parsing block: %s\n", sBlockText(bBlock) ));
    GbCommand = bBlock;
    BlockResetRead( GbCommand );

    GrMainResult.iCommand = CNONE;

		/* Parse the BLOCK */
		/* Keep parsing as long as the 'execute' command */
		/* keeps setting the GLOBAL variable GbContinueParsing */

    do {
	yyparse();
	if ( GrMainResult.iCommand == CQUIT ) {
	    if ( fINPUTFILE() != NULL )
	    	fclose( fINPUTFILE() );
	    INPUTPOPFILE();
	}
    } while ( fINPUTFILE() != NULL );

	/* Reset the bCommandFound variable as a new command may be
	   available in a new input */
    bCommandFound = FALSE;
    
	/* Copy the RESULT from the global result variable */

    *rPResult = GrMainResult;

}





/*
 *	ParseShutdown
 *
 *	Shutdown the parser, release all variables setup in ParseInit
 *	Only needed if debugging memory mgt, since prog mem is all
 *	freed when process exits anyway.
 */
void
ParseShutdown()
{
	if ( !iMemDebug )
		return;

	Destroy( (OBJEKT *)&aDummy );
	VariablesDestroy();
	ParmLibDestroy( &GplAllParameters );

	BlockDestroy( &GbExecute );

	HelpShutdown();
}
