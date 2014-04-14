/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with cifp or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum cifptokentype {
     ITEMNAME = 258,
     VALUE = 259,
     LOOP = 260,
     DATABLOCK = 261,
     UNKNOWN = 262,
     MISSING = 263
   };
#endif
/* Tokens.  */
#define ITEMNAME 258
#define VALUE 259
#define LOOP 260
#define DATABLOCK 261
#define UNKNOWN 262
#define MISSING 263




/* Copy the first part of user declarations.  */
#line 45 "cifparse.y"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cifparse.h"
int cifplex();
void cifperror();

static int curItemNo, curValueNo, itemIndex;
static int  *fieldList = NULL;



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 58 "cifparse.y"
{
	char TempBuffer[MAXVALUELENGTH+1];
}
/* Line 193 of yacc.c.  */
#line 130 "y.tab.c"
	YYSTYPE;
# define cifpstype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 143 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 cifptype_uint8;
#else
typedef unsigned char cifptype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 cifptype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char cifptype_int8;
#else
typedef short int cifptype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 cifptype_uint16;
#else
typedef unsigned short int cifptype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 cifptype_int16;
#else
typedef short int cifptype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined cifpoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined cifpoverflow || YYERROR_VERBOSE */


#if (! defined cifpoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union cifpalloc
{
  cifptype_int16 cifpss;
  YYSTYPE cifpvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union cifpalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (cifptype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T cifpi;				\
	  for (cifpi = 0; cifpi < (Count); cifpi++)	\
	    (To)[cifpi] = (From)[cifpi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T cifpnewbytes;						\
	YYCOPY (&cifpptr->Stack, Stack, cifpsize);				\
	Stack = &cifpptr->Stack;						\
	cifpnewbytes = cifpstacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	cifpptr += cifpnewbytes / sizeof (*cifpptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   20

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  9
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  11
/* YYNRULES -- Number of rules.  */
#define YYNRULES  18
/* YYNRULES -- Number of states.  */
#define YYNSTATES  23

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   263

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? cifptranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const cifptype_uint8 cifptranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const cifptype_uint8 cifpprhs[] =
{
       0,     0,     3,     5,     8,    11,    12,    15,    18,    21,
      24,    27,    29,    32,    34,    36,    38,    40,    42
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const cifptype_int8 cifprhs[] =
{
      10,     0,    -1,    12,    -1,    10,    11,    -1,    19,    12,
      -1,    -1,    12,    13,    -1,    16,    18,    -1,    14,    15,
      -1,    17,    16,    -1,    14,    16,    -1,    18,    -1,    15,
      18,    -1,     3,    -1,     5,    -1,     4,    -1,     7,    -1,
       8,    -1,     6,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const cifptype_uint8 cifprline[] =
{
       0,    70,    70,    71,    73,    76,    77,    80,    84,    90,
      95,   102,   107,   113,   120,   126,   131,   135,   141
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const cifptname[] =
{
  "$end", "error", "$undefined", "ITEMNAME", "VALUE", "LOOP", "DATABLOCK",
  "UNKNOWN", "MISSING", "$accept", "Datablocks", "Datablock", "Lines",
  "Line", "ItemNameList", "ValueList", "ItemName", "Loop", "Value",
  "DatablockName", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const cifptype_uint16 cifptoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const cifptype_uint8 cifpr1[] =
{
       0,     9,    10,    10,    11,    12,    12,    13,    13,    14,
      14,    15,    15,    16,    17,    18,    18,    18,    19
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const cifptype_uint8 cifpr2[] =
{
       0,     2,     1,     2,     2,     0,     2,     2,     2,     2,
       2,     1,     2,     1,     1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const cifptype_uint8 cifpdefact[] =
{
       5,     0,     2,     1,    18,     3,     5,    13,    14,     6,
       0,     0,     0,     4,    15,    16,    17,     8,    10,    11,
       7,     9,    12
};

/* YYDEFGOTO[NTERM-NUM].  */
static const cifptype_int8 cifpdefgoto[] =
{
      -1,     1,     5,     2,     9,    10,    17,    11,    12,    19,
       6
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -11
static const cifptype_int8 cifppact[] =
{
     -11,     0,    12,   -11,   -11,   -11,   -11,   -11,   -11,   -11,
       1,     6,    -1,    12,   -11,   -11,   -11,     6,   -11,   -11,
     -11,   -11,   -11
};

/* YYPGOTO[NTERM-NUM].  */
static const cifptype_int8 cifppgoto[] =
{
     -11,   -11,   -11,    -3,   -11,   -11,   -11,     8,   -11,   -10,
     -11
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const cifptype_uint8 cifptable[] =
{
       3,    20,     7,    13,     7,    14,     4,    22,    15,    16,
      14,     0,     0,    15,    16,     7,     0,     8,    18,     0,
      21
};

static const cifptype_int8 cifpcheck[] =
{
       0,    11,     3,     6,     3,     4,     6,    17,     7,     8,
       4,    -1,    -1,     7,     8,     3,    -1,     5,    10,    -1,
      12
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const cifptype_uint8 cifpstos[] =
{
       0,    10,    12,     0,     6,    11,    19,     3,     5,    13,
      14,    16,    17,    12,     4,     7,     8,    15,    16,    18,
      18,    16,    18
};

#define cifperrok		(cifperrstatus = 0)
#define cifpclearin	(cifpchar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto cifpacceptlab
#define YYABORT		goto cifpabortlab
#define YYERROR		goto cifperrorlab


/* Like YYERROR except do call cifperror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto cifperrlab

#define YYRECOVERING()  (!!cifperrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (cifpchar == YYEMPTY && cifplen == 1)				\
    {								\
      cifpchar = (Token);						\
      cifplval = (Value);						\
      cifptoken = YYTRANSLATE (cifpchar);				\
      YYPOPSTACK (1);						\
      goto cifpbackup;						\
    }								\
  else								\
    {								\
      cifperror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `cifplex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX cifplex (YYLEX_PARAM)
#else
# define YYLEX cifplex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (cifpdebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (cifpdebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      cifp_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
cifp_symbol_value_print (FILE *cifpoutput, int cifptype, YYSTYPE const * const cifpvaluep)
#else
static void
cifp_symbol_value_print (cifpoutput, cifptype, cifpvaluep)
    FILE *cifpoutput;
    int cifptype;
    YYSTYPE const * const cifpvaluep;
#endif
{
  if (!cifpvaluep)
    return;
# ifdef YYPRINT
  if (cifptype < YYNTOKENS)
    YYPRINT (cifpoutput, cifptoknum[cifptype], *cifpvaluep);
# else
  YYUSE (cifpoutput);
# endif
  switch (cifptype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
cifp_symbol_print (FILE *cifpoutput, int cifptype, YYSTYPE const * const cifpvaluep)
#else
static void
cifp_symbol_print (cifpoutput, cifptype, cifpvaluep)
    FILE *cifpoutput;
    int cifptype;
    YYSTYPE const * const cifpvaluep;
#endif
{
  if (cifptype < YYNTOKENS)
    YYFPRINTF (cifpoutput, "token %s (", cifptname[cifptype]);
  else
    YYFPRINTF (cifpoutput, "nterm %s (", cifptname[cifptype]);

  cifp_symbol_value_print (cifpoutput, cifptype, cifpvaluep);
  YYFPRINTF (cifpoutput, ")");
}

/*------------------------------------------------------------------.
| cifp_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
cifp_stack_print (cifptype_int16 *bottom, cifptype_int16 *top)
#else
static void
cifp_stack_print (bottom, top)
    cifptype_int16 *bottom;
    cifptype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (cifpdebug)							\
    cifp_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
cifp_reduce_print (YYSTYPE *cifpvsp, int cifprule)
#else
static void
cifp_reduce_print (cifpvsp, cifprule)
    YYSTYPE *cifpvsp;
    int cifprule;
#endif
{
  int cifpnrhs = cifpr2[cifprule];
  int cifpi;
  unsigned long int cifplno = cifprline[cifprule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     cifprule - 1, cifplno);
  /* The symbols being reduced.  */
  for (cifpi = 0; cifpi < cifpnrhs; cifpi++)
    {
      fprintf (stderr, "   $%d = ", cifpi + 1);
      cifp_symbol_print (stderr, cifprhs[cifpprhs[cifprule] + cifpi],
		       &(cifpvsp[(cifpi + 1) - (cifpnrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (cifpdebug)				\
    cifp_reduce_print (cifpvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int cifpdebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef cifpstrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define cifpstrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
cifpstrlen (const char *cifpstr)
#else
static YYSIZE_T
cifpstrlen (cifpstr)
    const char *cifpstr;
#endif
{
  YYSIZE_T cifplen;
  for (cifplen = 0; cifpstr[cifplen]; cifplen++)
    continue;
  return cifplen;
}
#  endif
# endif

# ifndef cifpstpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define cifpstpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
cifpstpcpy (char *cifpdest, const char *cifpsrc)
#else
static char *
cifpstpcpy (cifpdest, cifpsrc)
    char *cifpdest;
    const char *cifpsrc;
#endif
{
  char *cifpd = cifpdest;
  const char *cifps = cifpsrc;

  while ((*cifpd++ = *cifps++) != '\0')
    continue;

  return cifpd - 1;
}
#  endif
# endif

# ifndef cifptnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for cifperror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from cifptname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
cifptnamerr (char *cifpres, const char *cifpstr)
{
  if (*cifpstr == '"')
    {
      YYSIZE_T cifpn = 0;
      char const *cifpp = cifpstr;

      for (;;)
	switch (*++cifpp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++cifpp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (cifpres)
	      cifpres[cifpn] = *cifpp;
	    cifpn++;
	    break;

	  case '"':
	    if (cifpres)
	      cifpres[cifpn] = '\0';
	    return cifpn;
	  }
    do_not_strip_quotes: ;
    }

  if (! cifpres)
    return cifpstrlen (cifpstr);

  return cifpstpcpy (cifpres, cifpstr) - cifpres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
cifpsyntax_error (char *cifpresult, int cifpstate, int cifpchar)
{
  int cifpn = cifppact[cifpstate];

  if (! (YYPACT_NINF < cifpn && cifpn <= YYLAST))
    return 0;
  else
    {
      int cifptype = YYTRANSLATE (cifpchar);
      YYSIZE_T cifpsize0 = cifptnamerr (0, cifptname[cifptype]);
      YYSIZE_T cifpsize = cifpsize0;
      YYSIZE_T cifpsize1;
      int cifpsize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *cifparg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int cifpx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *cifpfmt;
      char const *cifpf;
      static char const cifpunexpected[] = "syntax error, unexpected %s";
      static char const cifpexpecting[] = ", expecting %s";
      static char const cifpor[] = " or %s";
      char cifpformat[sizeof cifpunexpected
		    + sizeof cifpexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof cifpor - 1))];
      char const *cifpprefix = cifpexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int cifpxbegin = cifpn < 0 ? -cifpn : 0;

      /* Stay within bounds of both cifpcheck and cifptname.  */
      int cifpchecklim = YYLAST - cifpn + 1;
      int cifpxend = cifpchecklim < YYNTOKENS ? cifpchecklim : YYNTOKENS;
      int cifpcount = 1;

      cifparg[0] = cifptname[cifptype];
      cifpfmt = cifpstpcpy (cifpformat, cifpunexpected);

      for (cifpx = cifpxbegin; cifpx < cifpxend; ++cifpx)
	if (cifpcheck[cifpx + cifpn] == cifpx && cifpx != YYTERROR)
	  {
	    if (cifpcount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		cifpcount = 1;
		cifpsize = cifpsize0;
		cifpformat[sizeof cifpunexpected - 1] = '\0';
		break;
	      }
	    cifparg[cifpcount++] = cifptname[cifpx];
	    cifpsize1 = cifpsize + cifptnamerr (0, cifptname[cifpx]);
	    cifpsize_overflow |= (cifpsize1 < cifpsize);
	    cifpsize = cifpsize1;
	    cifpfmt = cifpstpcpy (cifpfmt, cifpprefix);
	    cifpprefix = cifpor;
	  }

      cifpf = YY_(cifpformat);
      cifpsize1 = cifpsize + cifpstrlen (cifpf);
      cifpsize_overflow |= (cifpsize1 < cifpsize);
      cifpsize = cifpsize1;

      if (cifpsize_overflow)
	return YYSIZE_MAXIMUM;

      if (cifpresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *cifpp = cifpresult;
	  int cifpi = 0;
	  while ((*cifpp = *cifpf) != '\0')
	    {
	      if (*cifpp == '%' && cifpf[1] == 's' && cifpi < cifpcount)
		{
		  cifpp += cifptnamerr (cifpp, cifparg[cifpi++]);
		  cifpf += 2;
		}
	      else
		{
		  cifpp++;
		  cifpf++;
		}
	    }
	}
      return cifpsize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
cifpdestruct (const char *cifpmsg, int cifptype, YYSTYPE *cifpvaluep)
#else
static void
cifpdestruct (cifpmsg, cifptype, cifpvaluep)
    const char *cifpmsg;
    int cifptype;
    YYSTYPE *cifpvaluep;
#endif
{
  YYUSE (cifpvaluep);

  if (!cifpmsg)
    cifpmsg = "Deleting";
  YY_SYMBOL_PRINT (cifpmsg, cifptype, cifpvaluep, cifplocationp);

  switch (cifptype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int cifpparse (void *YYPARSE_PARAM);
#else
int cifpparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int cifpparse (void);
#else
int cifpparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int cifpchar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE cifplval;

/* Number of syntax errors so far.  */
int cifpnerrs;



/*----------.
| cifpparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
cifpparse (void *YYPARSE_PARAM)
#else
int
cifpparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
cifpparse (void)
#else
int
cifpparse ()

#endif
#endif
{
  
  int cifpstate;
  int cifpn;
  int cifpresult;
  /* Number of tokens to shift before error messages enabled.  */
  int cifperrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int cifptoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char cifpmsgbuf[128];
  char *cifpmsg = cifpmsgbuf;
  YYSIZE_T cifpmsg_alloc = sizeof cifpmsgbuf;
#endif

  /* Three stacks and their tools:
     `cifpss': related to states,
     `cifpvs': related to semantic values,
     `cifpls': related to locations.

     Refer to the stacks thru separate pointers, to allow cifpoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  cifptype_int16 cifpssa[YYINITDEPTH];
  cifptype_int16 *cifpss = cifpssa;
  cifptype_int16 *cifpssp;

  /* The semantic value stack.  */
  YYSTYPE cifpvsa[YYINITDEPTH];
  YYSTYPE *cifpvs = cifpvsa;
  YYSTYPE *cifpvsp;



#define YYPOPSTACK(N)   (cifpvsp -= (N), cifpssp -= (N))

  YYSIZE_T cifpstacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE cifpval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int cifplen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  cifpstate = 0;
  cifperrstatus = 0;
  cifpnerrs = 0;
  cifpchar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  cifpssp = cifpss;
  cifpvsp = cifpvs;

  goto cifpsetstate;

/*------------------------------------------------------------.
| cifpnewstate -- Push a new state, which is found in cifpstate.  |
`------------------------------------------------------------*/
 cifpnewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  cifpssp++;

 cifpsetstate:
  *cifpssp = cifpstate;

  if (cifpss + cifpstacksize - 1 <= cifpssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T cifpsize = cifpssp - cifpss + 1;

#ifdef cifpoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *cifpvs1 = cifpvs;
	cifptype_int16 *cifpss1 = cifpss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if cifpoverflow is a macro.  */
	cifpoverflow (YY_("memory exhausted"),
		    &cifpss1, cifpsize * sizeof (*cifpssp),
		    &cifpvs1, cifpsize * sizeof (*cifpvsp),

		    &cifpstacksize);

	cifpss = cifpss1;
	cifpvs = cifpvs1;
      }
#else /* no cifpoverflow */
# ifndef YYSTACK_RELOCATE
      goto cifpexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= cifpstacksize)
	goto cifpexhaustedlab;
      cifpstacksize *= 2;
      if (YYMAXDEPTH < cifpstacksize)
	cifpstacksize = YYMAXDEPTH;

      {
	cifptype_int16 *cifpss1 = cifpss;
	union cifpalloc *cifpptr =
	  (union cifpalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (cifpstacksize));
	if (! cifpptr)
	  goto cifpexhaustedlab;
	YYSTACK_RELOCATE (cifpss);
	YYSTACK_RELOCATE (cifpvs);

#  undef YYSTACK_RELOCATE
	if (cifpss1 != cifpssa)
	  YYSTACK_FREE (cifpss1);
      }
# endif
#endif /* no cifpoverflow */

      cifpssp = cifpss + cifpsize - 1;
      cifpvsp = cifpvs + cifpsize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) cifpstacksize));

      if (cifpss + cifpstacksize - 1 <= cifpssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", cifpstate));

  goto cifpbackup;

/*-----------.
| cifpbackup.  |
`-----------*/
cifpbackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  cifpn = cifppact[cifpstate];
  if (cifpn == YYPACT_NINF)
    goto cifpdefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (cifpchar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      cifpchar = YYLEX;
    }

  if (cifpchar <= YYEOF)
    {
      cifpchar = cifptoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      cifptoken = YYTRANSLATE (cifpchar);
      YY_SYMBOL_PRINT ("Next token is", cifptoken, &cifplval, &cifplloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  cifpn += cifptoken;
  if (cifpn < 0 || YYLAST < cifpn || cifpcheck[cifpn] != cifptoken)
    goto cifpdefault;
  cifpn = cifptable[cifpn];
  if (cifpn <= 0)
    {
      if (cifpn == 0 || cifpn == YYTABLE_NINF)
	goto cifperrlab;
      cifpn = -cifpn;
      goto cifpreduce;
    }

  if (cifpn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (cifperrstatus)
    cifperrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", cifptoken, &cifplval, &cifplloc);

  /* Discard the shifted token unless it is eof.  */
  if (cifpchar != YYEOF)
    cifpchar = YYEMPTY;

  cifpstate = cifpn;
  *++cifpvsp = cifplval;

  goto cifpnewstate;


/*-----------------------------------------------------------.
| cifpdefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
cifpdefault:
  cifpn = cifpdefact[cifpstate];
  if (cifpn == 0)
    goto cifperrlab;
  goto cifpreduce;


/*-----------------------------.
| cifpreduce -- Do a reduction.  |
`-----------------------------*/
cifpreduce:
  /* cifpn is the number of a rule to reduce with.  */
  cifplen = cifpr2[cifpn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  cifpval = cifpvsp[1-cifplen];


  YY_REDUCE_PRINT (cifpn);
  switch (cifpn)
    {
        case 7:
#line 81 "cifparse.y"
    {
  ndb_cif_process_item_name_value_pair();
}
    break;

  case 9:
#line 91 "cifparse.y"
    {
  ndb_cif_process_loop_declaration();
}
    break;

  case 10:
#line 96 "cifparse.y"
    {
  ndb_cif_process_item_name_list();
}
    break;

  case 11:
#line 103 "cifparse.y"
    {
  ndb_cif_rewind_column();
  ndb_cif_process_value_list();
}
    break;

  case 12:
#line 108 "cifparse.y"
    {
  ndb_cif_process_value_list();
}
    break;

  case 13:
#line 114 "cifparse.y"
    {
  /*  sprintf(TempKeyword,"%s",$1); */
  strncpy(TempKeyword,(cifpvsp[(1) - (1)].TempBuffer),MxNameLen);
}
    break;

  case 14:
#line 121 "cifparse.y"
    {
  curItemNo = 0;  curValueNo = 0;  itemIndex = 0;
}
    break;

  case 15:
#line 127 "cifparse.y"
    {
  /*  sprintf(TempValue,"%s",$1); */
  strncpy(TempValue,(cifpvsp[(1) - (1)].TempBuffer),MAXVALUELENGTH);
}
    break;

  case 16:
#line 132 "cifparse.y"
    {
  strcpy(TempValue,"");
}
    break;

  case 17:
#line 136 "cifparse.y"
    {
  strcpy(TempValue,"");
}
    break;

  case 18:
#line 142 "cifparse.y"
    {
  int idatablock;
  idatablock = ndb_cif_get_datablock_id(&((cifpvsp[(1) - (1)].TempBuffer))[5]);
  if (idatablock == 0)
    ndb_cif_new_datablock(&((cifpvsp[(1) - (1)].TempBuffer))[5]);
  else {
    ndb_cif_move_datablock(&((cifpvsp[(1) - (1)].TempBuffer))[5]);
    ndb_cif_reset_datablock();
  }
}
    break;


/* Line 1267 of yacc.c.  */
#line 1432 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", cifpr1[cifpn], &cifpval, &cifploc);

  YYPOPSTACK (cifplen);
  cifplen = 0;
  YY_STACK_PRINT (cifpss, cifpssp);

  *++cifpvsp = cifpval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  cifpn = cifpr1[cifpn];

  cifpstate = cifppgoto[cifpn - YYNTOKENS] + *cifpssp;
  if (0 <= cifpstate && cifpstate <= YYLAST && cifpcheck[cifpstate] == *cifpssp)
    cifpstate = cifptable[cifpstate];
  else
    cifpstate = cifpdefgoto[cifpn - YYNTOKENS];

  goto cifpnewstate;


/*------------------------------------.
| cifperrlab -- here on detecting error |
`------------------------------------*/
cifperrlab:
  /* If not already recovering from an error, report this error.  */
  if (!cifperrstatus)
    {
      ++cifpnerrs;
#if ! YYERROR_VERBOSE
      cifperror (YY_("syntax error"));
#else
      {
	YYSIZE_T cifpsize = cifpsyntax_error (0, cifpstate, cifpchar);
	if (cifpmsg_alloc < cifpsize && cifpmsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T cifpalloc = 2 * cifpsize;
	    if (! (cifpsize <= cifpalloc && cifpalloc <= YYSTACK_ALLOC_MAXIMUM))
	      cifpalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (cifpmsg != cifpmsgbuf)
	      YYSTACK_FREE (cifpmsg);
	    cifpmsg = (char *) YYSTACK_ALLOC (cifpalloc);
	    if (cifpmsg)
	      cifpmsg_alloc = cifpalloc;
	    else
	      {
		cifpmsg = cifpmsgbuf;
		cifpmsg_alloc = sizeof cifpmsgbuf;
	      }
	  }

	if (0 < cifpsize && cifpsize <= cifpmsg_alloc)
	  {
	    (void) cifpsyntax_error (cifpmsg, cifpstate, cifpchar);
	    cifperror (cifpmsg);
	  }
	else
	  {
	    cifperror (YY_("syntax error"));
	    if (cifpsize != 0)
	      goto cifpexhaustedlab;
	  }
      }
#endif
    }



  if (cifperrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (cifpchar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (cifpchar == YYEOF)
	    YYABORT;
	}
      else
	{
	  cifpdestruct ("Error: discarding",
		      cifptoken, &cifplval);
	  cifpchar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto cifperrlab1;


/*---------------------------------------------------.
| cifperrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
cifperrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label cifperrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto cifperrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (cifplen);
  cifplen = 0;
  YY_STACK_PRINT (cifpss, cifpssp);
  cifpstate = *cifpssp;
  goto cifperrlab1;


/*-------------------------------------------------------------.
| cifperrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
cifperrlab1:
  cifperrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      cifpn = cifppact[cifpstate];
      if (cifpn != YYPACT_NINF)
	{
	  cifpn += YYTERROR;
	  if (0 <= cifpn && cifpn <= YYLAST && cifpcheck[cifpn] == YYTERROR)
	    {
	      cifpn = cifptable[cifpn];
	      if (0 < cifpn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (cifpssp == cifpss)
	YYABORT;


      cifpdestruct ("Error: popping",
		  cifpstos[cifpstate], cifpvsp);
      YYPOPSTACK (1);
      cifpstate = *cifpssp;
      YY_STACK_PRINT (cifpss, cifpssp);
    }

  if (cifpn == YYFINAL)
    YYACCEPT;

  *++cifpvsp = cifplval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", cifpstos[cifpn], cifpvsp, cifplsp);

  cifpstate = cifpn;
  goto cifpnewstate;


/*-------------------------------------.
| cifpacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
cifpacceptlab:
  cifpresult = 0;
  goto cifpreturn;

/*-----------------------------------.
| cifpabortlab -- YYABORT comes here.  |
`-----------------------------------*/
cifpabortlab:
  cifpresult = 1;
  goto cifpreturn;

#ifndef cifpoverflow
/*-------------------------------------------------.
| cifpexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
cifpexhaustedlab:
  cifperror (YY_("memory exhausted"));
  cifpresult = 2;
  /* Fall through.  */
#endif

cifpreturn:
  if (cifpchar != YYEOF && cifpchar != YYEMPTY)
     cifpdestruct ("Cleanup: discarding lookahead",
		 cifptoken, &cifplval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (cifplen);
  YY_STACK_PRINT (cifpss, cifpssp);
  while (cifpssp != cifpss)
    {
      cifpdestruct ("Cleanup: popping",
		  cifpstos[*cifpssp], cifpvsp);
      YYPOPSTACK (1);
    }
#ifndef cifpoverflow
  if (cifpss != cifpssa)
    YYSTACK_FREE (cifpss);
#endif
#if YYERROR_VERBOSE
  if (cifpmsg != cifpmsgbuf)
    YYSTACK_FREE (cifpmsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (cifpresult);
}


#line 154 "cifparse.y"


void cifperror(s)
char *s;
/* ---------------------------------------------------------------------------
 * Purpose:  cifperror()
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
void cifpprint(FILE *file, int type, YYSTYPE value)
{
  fprintf(file,"\nType  = %d\n",type);
  fprintf(file,"Value = %s\n",value.TempBuffer);
}
#endif


