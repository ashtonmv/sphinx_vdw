/* A Bison parser, made by GNU Bison 2.7.  */

/* Bison interface for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2012 Free Software Foundation, Inc.
   
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

#ifndef YY_SXPARSER_SXPARSER_TAB_HPP_INCLUDED
# define YY_SXPARSER_SXPARSER_TAB_HPP_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int SxParser_debug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     STR = 259,
     VEC = 260,
     VAR = 261,
     FUNC = 262,
     STRFUNC = 263,
     LIST = 264,
     NEG = 265
   };
#endif
/* Tokens.  */
#define NUM 258
#define STR 259
#define VEC 260
#define VAR 261
#define FUNC 262
#define STRFUNC 263
#define LIST 264
#define NEG 265



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 2058 of yacc.c  */
#line 38 "SxParser.ypp"

   char     *str;
   double    val;
   SxSymbol *sym;


/* Line 2058 of yacc.c  */
#line 84 "SxParser.tab.hpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE SxParser_lval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int SxParser_parse (void *YYPARSE_PARAM);
#else
int SxParser_parse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int SxParser_parse (void);
#else
int SxParser_parse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_SXPARSER_SXPARSER_TAB_HPP_INCLUDED  */
