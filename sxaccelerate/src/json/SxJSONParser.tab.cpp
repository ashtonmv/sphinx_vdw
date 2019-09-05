/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         SxJSONParser_parse
#define yylex           SxJSONParser_lex
#define yyerror         SxJSONParser_error
#define yydebug         SxJSONParser_debug
#define yynerrs         SxJSONParser_nerrs


/* Copy the first part of user declarations.  */
#line 9 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:339  */

#  define YYDEBUG 1
#  include <SxMacroLib.h>
#  include <SxString.h>
#  include <SxGProps.h>
class SxJSONParser;
#  include <SxJSONParser.tab.hpp>
#  define SXPARSER_TYPE SxJSONParser
#  include <SxJSONParser.h>
#  define ast parserPtr



#line 86 "SxJSONParser.tab.cpp" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* In a future release of Bison, this section will be replaced
   by #include "SxJSONParser.tab.hpp".  */
#ifndef YY_SXJSONPARSER_SXJSONPARSER_TAB_HPP_INCLUDED
# define YY_SXJSONPARSER_SXJSONPARSER_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int SxJSONParser_debug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    TK_INTEGER = 258,
    TK_FLOAT = 259,
    TK_BOOL = 260,
    TK_NULL = 261,
    TK_STRING = 262
  };
#endif
/* Tokens.  */
#define TK_INTEGER 258
#define TK_FLOAT 259
#define TK_BOOL 260
#define TK_NULL 261
#define TK_STRING 262

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 24 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:355  */

   ssize_t   intVal;
   bool      boolVal;
   double    floatVal;
   SxString *strPtr;
   SxGProps *nPtr;

#line 148 "SxJSONParser.tab.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif



int SxJSONParser_parse (SxJSONParser* parserPtr);

#endif /* !YY_SXJSONPARSER_SXJSONPARSER_TAB_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 178 "SxJSONParser.tab.cpp" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
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
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
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
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
             && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  5
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   29

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  14
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  11
/* YYNRULES -- Number of rules.  */
#define YYNRULES  22
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  35

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   262

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    10,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    11,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    12,     2,    13,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     8,     2,     9,     2,     2,     2,     2,
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
       5,     6,     7
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    62,    62,    67,    74,    80,    88,   103,   119,   130,
     139,   147,   167,   190,   196,   200,   204,   211,   217,   223,
     229,   237,   244
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "TK_INTEGER", "TK_FLOAT", "TK_BOOL",
  "TK_NULL", "TK_STRING", "'{'", "'}'", "','", "':'", "'['", "']'",
  "$accept", "main", "object", "members", "member", "array", "elements",
  "element", "value", "push", "pop", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   123,   125,
      44,    58,    91,    93
};
# endif

#define YYPACT_NINF -18

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-18)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
     -18,     6,   -18,   -18,    13,   -18,    15,    -2,     2,   -18,
      17,   -18,   -18,   -18,   -18,   -18,   -18,   -18,   -18,   -18,
      -3,   -18,   -18,    12,   -18,   -18,    16,   -18,    12,   -18,
     -18,   -18,   -18,   -18,   -18
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
      21,     0,     2,     3,     0,     1,     0,    21,     0,    22,
       0,     6,    17,    18,    19,    20,    16,    22,    14,    15,
       0,    11,    13,    21,     4,    22,     0,     9,    21,    22,
       8,     5,     7,    12,    10
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -18,   -18,    20,   -18,     3,    28,   -18,   -14,   -18,   -18,
     -17
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    18,    10,    11,    19,    20,    21,    22,     4,
      24
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      27,    12,    13,    14,    15,    16,     5,    28,    31,    30,
      29,    17,    34,    23,    33,    12,    13,    14,    15,    16,
       2,     6,     8,     8,     9,     7,    25,    26,     3,    32
};

static const yytype_uint8 yycheck[] =
{
      17,     3,     4,     5,     6,     7,     0,    10,    25,    23,
      13,    13,    29,    11,    28,     3,     4,     5,     6,     7,
       0,     8,     7,     7,     9,    12,     9,    10,     0,    26
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    15,    16,    19,    23,     0,     8,    12,     7,     9,
      17,    18,     3,     4,     5,     6,     7,    13,    16,    19,
      20,    21,    22,    11,    24,     9,    10,    24,    10,    13,
      21,    24,    18,    21,    24
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    14,    15,    15,    16,    16,    17,    17,    18,    19,
      19,    20,    20,    21,    22,    22,    22,    22,    22,    22,
      22,    23,    24
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     4,     5,     1,     3,     3,     4,
       5,     1,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     0,     0
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (&yylloc, parserPtr, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (0)
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL

/* Print *YYLOCP on YYO.  Private, do not rely on its existence. */

YY_ATTRIBUTE_UNUSED
static unsigned
yy_location_print_ (FILE *yyo, YYLTYPE const * const yylocp)
{
  unsigned res = 0;
  int end_col = 0 != yylocp->last_column ? yylocp->last_column - 1 : 0;
  if (0 <= yylocp->first_line)
    {
      res += YYFPRINTF (yyo, "%d", yylocp->first_line);
      if (0 <= yylocp->first_column)
        res += YYFPRINTF (yyo, ".%d", yylocp->first_column);
    }
  if (0 <= yylocp->last_line)
    {
      if (yylocp->first_line < yylocp->last_line)
        {
          res += YYFPRINTF (yyo, "-%d", yylocp->last_line);
          if (0 <= end_col)
            res += YYFPRINTF (yyo, ".%d", end_col);
        }
      else if (0 <= end_col && yylocp->first_column < end_col)
        res += YYFPRINTF (yyo, "-%d", end_col);
    }
  return res;
 }

#  define YY_LOCATION_PRINT(File, Loc)          \
  yy_location_print_ (File, &(Loc))

# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, Location, parserPtr); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, SxJSONParser* parserPtr)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (yylocationp);
  YYUSE (parserPtr);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  switch (yytype)
    {
          case 3: /* TK_INTEGER  */
#line 48 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "TK_INTEGER, val='%zd'", ((*yyvaluep).intVal)); }
#line 749 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 4: /* TK_FLOAT  */
#line 49 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "TK_FLOAT, val='%f'", ((*yyvaluep).floatVal)); }
#line 755 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 5: /* TK_BOOL  */
#line 50 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "TK_BOOL, val='%d'", ((*yyvaluep).boolVal)); }
#line 761 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 7: /* TK_STRING  */
#line 47 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "TK_STRING, val='%s'", ((*yyvaluep).strPtr)->ascii()); }
#line 767 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 15: /* main  */
#line 52 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "main"); }
#line 773 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 16: /* object  */
#line 53 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "object"); }
#line 779 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 17: /* members  */
#line 54 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "members"); }
#line 785 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 18: /* member  */
#line 55 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "member"); }
#line 791 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 19: /* array  */
#line 56 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "array"); }
#line 797 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 20: /* elements  */
#line 57 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "elements"); }
#line 803 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 21: /* element  */
#line 58 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "element"); }
#line 809 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;

    case 22: /* value  */
#line 59 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:684  */
      { fprintf (yyoutput, "value"); }
#line 815 "SxJSONParser.tab.cpp" /* yacc.c:684  */
        break;


      default:
        break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, SxJSONParser* parserPtr)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  YY_LOCATION_PRINT (yyoutput, *yylocationp);
  YYFPRINTF (yyoutput, ": ");
  yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp, parserPtr);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, SxJSONParser* parserPtr)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                       , &(yylsp[(yyi + 1) - (yynrhs)])                       , parserPtr);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, yylsp, Rule, parserPtr); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
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

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, SxJSONParser* parserPtr)
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);
  YYUSE (parserPtr);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yytype)
    {
          case 7: /* TK_STRING  */
#line 37 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1257  */
      { delete ((*yyvaluep).strPtr); }
#line 1165 "SxJSONParser.tab.cpp" /* yacc.c:1257  */
        break;


      default:
        break;
    }
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (SxJSONParser* parserPtr)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

/* Location data for the lookahead symbol.  */
static YYLTYPE yyloc_default
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  = { 1, 1, 1, 1 }
# endif
;
YYLTYPE yylloc = yyloc_default;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.
       'yyls': related to locations.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls;
    YYLTYPE *yylsp;

    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[3];

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yylsp = yyls = yylsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  yylsp[0] = yylloc;
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;
        YYLTYPE *yyls1 = yyls;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yyls1, yysize * sizeof (*yylsp),
                    &yystacksize);

        yyls = yyls1;
        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
        YYSTACK_RELOCATE (yyls_alloc, yyls);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex (&yylval, &yylloc, scannerPtr);
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  *++yylsp = yylloc;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location.  */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 62 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("main#1");
                  SxGProps *root = parserPtr->getRoot ();
                  ast->addEdge (root, (yyval.nPtr));
               }
#line 1463 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 3:
#line 67 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("main#2");
                  SxGProps *root = parserPtr->getRoot ();
                  ast->addEdge (root, (yyval.nPtr));
               }
#line 1473 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 4:
#line 74 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("object#1");
                  (yyval.nPtr) = (yyvsp[-3].nPtr);
                  (yyval.nPtr)->setProperty ("__sx_Key", "undefined");
                  (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[-2]),(yylsp[-1])));
               }
#line 1484 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 5:
#line 80 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("object#2");
                  (yyval.nPtr) = (yyvsp[-2].nPtr);
                  (yyval.nPtr)->setProperty ("__sx_Key", "undefined");
                  (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[-3]),(yylsp[-1])));
               }
#line 1495 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 6:
#line 88 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("members#1: " << (yyvsp[0].nPtr));
                  (yyval.nPtr) = parserPtr->peek ();
                  SxVariantType::DataType chT =
                  (SxVariantType::DataType)(yyvsp[0].nPtr)->getProperty ("__sx_Value").getType ();

                  if (!(chT == SxVariantType::List || chT == SxVariantType::Group))
                  {
                     (yyval.nPtr)->setProperty ((yyvsp[0].nPtr)->getProperty ("__sx_Key").getString (),
                                      (yyvsp[0].nPtr)->getProperty ("__sx_Value"));
                     ast->getAst ()->removeElement (*(yyvsp[0].nPtr));
                  } else {
                     ast->addEdge ((yyval.nPtr), (yyvsp[0].nPtr));
                  }
               }
#line 1515 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 7:
#line 103 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("members#2: " << (yyvsp[-2].nPtr) << "-" << (yyvsp[0].nPtr));
                  (yyval.nPtr) = parserPtr->peek ();
                  SxVariantType::DataType chT =
                  (SxVariantType::DataType)(yyvsp[0].nPtr)->getProperty ("__sx_Value").getType ();

                  if (!(chT == SxVariantType::List || chT == SxVariantType::Group)) {
                     (yyval.nPtr)->setProperty ((yyvsp[0].nPtr)->getProperty ("__sx_Key").getString (),
                                      (yyvsp[0].nPtr)->getProperty ("__sx_Value"));
                     ast->getAst ()->removeElement (*(yyvsp[0].nPtr));
                  } else {
                     ast->addEdge ((yyval.nPtr), (yyvsp[0].nPtr));
                  }
               }
#line 1534 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 8:
#line 119 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("member#1");
                  (yyval.nPtr) = (yyvsp[0].nPtr);
                  if (!ast->validateKey (*(yyvsp[-2].strPtr)))
                     SxJSONParser_error (&(yylsp[-2]), ast, SxString("Invalid key: ")+ *(yyvsp[-2].strPtr));
                  (yyval.nPtr)->setProperty ("__sx_Key", *(yyvsp[-2].strPtr));
                  (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[-2]),(yylsp[0])));
                  delete (yyvsp[-2].strPtr);
               }
#line 1548 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 9:
#line 130 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("array#1");
                  (yyval.nPtr) = (yyvsp[-3].nPtr);
                  (yyval.nPtr)->setProperty ("__sx_Key", "undefined");
                  SxVariant v;
                  v.setType ((int)SxVariantType::List);
                  (yyval.nPtr)->setProperty ("__sx_Value", v);
                  (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[-2]),(yylsp[-1])));
               }
#line 1562 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 10:
#line 139 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("array#2");
                  (yyval.nPtr) = (yyvsp[-4].nPtr);
                  (yyval.nPtr)->setProperty ("__sx_Key", "undefined");
                  (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[-3]),(yylsp[-1])));
               }
#line 1573 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 11:
#line 147 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("elements#1: " << (yyvsp[0].nPtr));
                  (yyval.nPtr) = parserPtr->peek ();
                  if ((yyval.nPtr)->getProperty ("__sx_Value").getType () != SxVariantType::List) {
                     SxVariant v;
                     v.setType ((int)SxVariantType::List);
                     v.set (SxList<SxVariant>());
                     (yyval.nPtr)->setProperty ("__sx_Value", v);
                  }
                  SxVariantType::DataType chT =
                  (SxVariantType::DataType)(yyvsp[0].nPtr)->getProperty ("__sx_Value").getType ();

                  if (!(chT == SxVariantType::List || chT == SxVariantType::Group))
                  {
                     (yyval.nPtr)->getProperty ("__sx_Value").append ((yyvsp[0].nPtr)->getProperty ("__sx_Value"));
                     ast->getAst ()->removeElement (*(yyvsp[0].nPtr));
                  } else {
                     ast->addEdge ((yyval.nPtr), (yyvsp[0].nPtr));
                  }
               }
#line 1598 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 12:
#line 167 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                  SX_DBG_YACC ("elements#2: " << (yyvsp[-2].nPtr) << "-" << (yyvsp[0].nPtr));
                  (yyval.nPtr) = parserPtr->peek ();

                  if ((yyval.nPtr)->getProperty ("__sx_Value").getType () != SxVariantType::List) {
                     SxVariant v;
                     v.setType ((int)SxVariantType::List);
                     (yyval.nPtr)->setProperty ("__sx_Value", v);
                  }

                  SxVariantType::DataType chT =
                  (SxVariantType::DataType)(yyvsp[0].nPtr)->getProperty ("__sx_Value").getType ();

                  if (!(chT == SxVariantType::List || chT == SxVariantType::Group))
                  {
                     (yyval.nPtr)->getProperty ("__sx_Value").append ((yyvsp[0].nPtr)->getProperty ("__sx_Value"));
                     ast->getAst ()->removeElement (*(yyvsp[0].nPtr));
                  } else {
                     ast->addEdge ((yyval.nPtr), (yyvsp[0].nPtr));
                  }
               }
#line 1624 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 13:
#line 190 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_YACC ("element#1");
                    (yyval.nPtr) = (yyvsp[0].nPtr);
                 }
#line 1633 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 14:
#line 196 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_YACC ("value#1");
                    (yyval.nPtr) = (yyvsp[0].nPtr);
                 }
#line 1642 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 15:
#line 200 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_YACC ("value#2");
                    (yyval.nPtr) = (yyvsp[0].nPtr);
                 }
#line 1651 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 16:
#line 204 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_YACC ("value#3");
                    (yyval.nPtr) = ast->addNode (SxVariantType::String);
                    (yyval.nPtr)->setProperty ("__sx_Value", *(yyvsp[0].strPtr));
                    (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[0])));
                    delete (yyvsp[0].strPtr);
                 }
#line 1663 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 17:
#line 211 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_YACC ("value#4");
                    (yyval.nPtr) = ast->addNode (SxVariantType::Int);
                    (yyval.nPtr)->setProperty ("__sx_Value", (yyvsp[0].intVal));
                    (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[0])));
                 }
#line 1674 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 18:
#line 217 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_YACC ("value#5");
                    (yyval.nPtr) = ast->addNode (SxVariantType::Double);
                    (yyval.nPtr)->setProperty ("__sx_Value", (yyvsp[0].floatVal));
                    (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[0])));
                 }
#line 1685 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 19:
#line 223 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_YACC ("value#6");
                    (yyval.nPtr) = ast->addNode (SxVariantType::Bool);
                    (yyval.nPtr)->setProperty ("__sx_Value", (yyvsp[0].boolVal));
                    (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[0])));
                 }
#line 1696 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 20:
#line 229 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_YACC ("value#6");
                    (yyval.nPtr) = ast->addNode (SxVariantType::String);
                    (yyval.nPtr)->setProperty ("__sx_Value", SxString (""));
                    (yyval.nPtr)->getProperty ("__sx_Value").setTag (SXTAG_YACC((yylsp[0])));
                 }
#line 1707 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 21:
#line 237 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    (yyval.nPtr) = ast->push (SxVariantType::Group);
                    (yyval.nPtr)->setProperty ("__sx_Key", "undefined");
                    SX_DBG_MSG ("push: " << (yyval.nPtr));
                 }
#line 1717 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;

  case 22:
#line 244 "/home/freysoldt/devel/sphinx/sxaccelerate/src/json/SxJSONParser.ypp" /* yacc.c:1646  */
    {
                    SX_DBG_MSG ("pop");
                    (yyval.nPtr) = parserPtr->pop ();
                 }
#line 1726 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
    break;


#line 1730 "SxJSONParser.tab.cpp" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (&yylloc, parserPtr, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (&yylloc, parserPtr, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }

  yyerror_range[1] = yylloc;

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, &yylloc, parserPtr);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  yyerror_range[1] = yylsp[1-yylen];
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, yylsp, parserPtr);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  yyerror_range[2] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the lookahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
  *++yylsp = yyloc;

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (&yylloc, parserPtr, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, &yylloc, parserPtr);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, yylsp, parserPtr);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
