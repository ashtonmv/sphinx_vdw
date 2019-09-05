#ifndef _SX_PARSER_BASE_H_
#define _SX_PARSER_BASE_H_


#include <SxException.h>
#include <SxFunction.h>
#include <SxString.h>
#include <SxSed.h>
#include <SxSortedList.h>
#include <SxParserKit.h>
#include <iostream>
#include <sstream>



/** \brief This class ...

    \author Sixten Boeck */
class SX_EXPORT_PARSER_KIT SxParserBase
{
   public:

      typedef SxFunction<char,void>  ReaderFunction;
      void *scannerPtr;          // .ypp: lex-param
      std::istream *inStream;
      ReaderFunction readerCB;

      SxParserBase ();
      virtual ~SxParserBase ();

      void setLineOffset (ssize_t);
      void setSearchPath (const SxString &);
      void setMaxIncludeDepth (ssize_t);

      void enableTraces (bool state=true);

      int readString (const SxString &,
                      const SxString &label="(stdin)", ssize_t maxErrors=1);
      int readStream (std::istream *,
                      const SxString &label="(cin)", ssize_t maxErrors=1);
      int readFile   (const SxString &, ssize_t maxErrors=1);
      int readFunction (const ReaderFunction &,
                        const SxString &label="(in)", ssize_t maxErrors=1);

      virtual bool resume () const;
      virtual void errorHandler (const SxString &msg,
                                 ssize_t line0, ssize_t col0,
                                 ssize_t row1, ssize_t line1);

      SxString getInFile () const { return inFile; }
      void     setInFile (const SxString &);

      // --- Lex input file buffer, only for internal usage
      ssize_t getLineOffset () const { return lineOffset; }
      void validateIncludes (const SxString &) const;
      bool pushIncludes (const SxString &, int, int, int, int);
      bool popInclude (int *, int *, int *, int *);
      void setIncLoc (int, int, int, int);
      bool processNextInclude ();
      SxSortedList<SxString> resolveFiles (const SxString &) const;

      // --- Lex stack, only for internal usage (SXPARSER_PUSH_STATE)
      void pushLex (int, int, int);
      void repushLex (int, int);
      void popLex (int *, int *);
      int getLexState ();
      void appendLex (const SxString &);
      SxList<SxString> collectLex ();

      SxString getLexTag ();

      // --- Lex indent/unindent management
      void setIndent (int); // modify indent counter of current line
      int getIndent ();     // return curLineIndent
      void updateIndent ();
      int indentDiff ();   // return curLineIndent - prevLineIndent
      bool unindent (); // return true if --prevLineIndent > curLineIndent

      // --- Yacc search path
      SxList<SxString> getSearchPath () const { return searchPath; }

   protected:

      class Error {
         public:
            ssize_t line0, col0, line1, col1;
            SxString msg;
            Error (ssize_t line0_=-1, ssize_t col0_=-1,
                   ssize_t line1_=-1, ssize_t col1_=-1,
                   const SxString &msg_="")
               : line0(line0_), col0(col0_), line1(line1_), col1(col1_),
                 msg(msg_) { /* empty */ }
           ~Error () { /* empty */ }
      };
      SxList<Error> errors;

      ssize_t lineOffset;
      SxString inFile;
      ssize_t maxErrors;
      bool traces;

      // --- Lex input buffers
      SxList<SxString>          lexIncludes;    // filename stack
      SxList<SxString>          lexIncStack;    // not yet opened files
      SxList<std::istream *>    lexStreams;
      SxList<SxArray<int> >     lexLocs;
      SxArray<int>              lexIncLoc;      // location of wildcard
      // --- Lex state stack
      SxList<int>               lexStates;
      SxList<SxList<SxString> > lexStack;
      SxList<int>      lexCols, lexLines;
      int curLineIndent, prevLineIndent;

      SxList<SxString> searchPath;
      ssize_t maxIncludeDepth;


      virtual void initScanner (bool)=0;
      virtual void destroyScanner ()=0;
      virtual int parse ()=0;

      void handleExceptions ();
};

#endif /* _SX_PARSER_H_ */


// --- convert SXPARSER_TYPE##_error to its actual value, SxFooParser_error
#define SXPARSER_FUNC_WRAPPER(x,func) x ## func
#define SXPARSER_FUNC(x,func) SXPARSER_FUNC_WRAPPER(x,func)

#include <SxString.h>

// --- included from lexer *.lpp
#ifdef FLEX_SCANNER
#  ifndef SX_PARSER_LEX_H
#     define SX_PARSER_LEX_H
#     ifndef SXPARSER_TYPE
#        error "Including SxParserBase.h from lexer requires SXPARSER_TYPE to be defined."
#     endif /* SX_PARSER_TYPE */

#     define SXTAG_LEX()  yyextra->getLexTag () + "-"                         \
                          + yylloc->last_line + "."                           \
                          + yylloc->last_column;

#     define YY_USER_INIT  yylloc->last_line += yyextra->getLineOffset();
#     define YY_USER_ACTION                                                   \
         yylloc->first_line   = yylloc->last_line;                            \
         yylloc->first_column = yylloc->last_column;                          \
         for(int i = 0; yytext[i] != '\0'; i++)  {                            \
            unsigned char c = yytext[i];                                      \
             if(c == '\n') {                                                  \
                 ++yylloc->last_line;                                         \
                 yylloc->last_column = 1;                                     \
             }  else  {                                                       \
                if (c <= 0x7F)  /* 7bit ASCII */                              \
                   ++yylloc->last_column;                                     \
                else if (c <= 0xBF) /* continuing UTF8 */                     \
                   /* do nothing */;                                          \
                else if (c <= 0xDF) /* begin 2-byte seq. */                   \
                   yylloc->last_column += 1;                                  \
                else if (c <= 0xEF) /* begin 3-byte seq. */                   \
                   yylloc->last_column += 1;                                  \
                else                /* begin 4-byte seq. */                   \
                   yylloc->last_column += 1;                                  \
             }                                                                \
         }

#     define YY_INPUT(buf,result,max_size) {                                  \
         if (yyextra->readerCB) {                                             \
            char c = yyextra->readerCB ();                                    \
            if (c == '\r')  c = yyextra->readerCB ();                         \
            if (c == '\0')  result = YY_NULL;                                 \
            else { buf[0] = c; result = 1; }                                  \
         } else {                                                             \
            char c = yyextra->inStream->get();                                \
            if (c == '\r')  c = yyextra->inStream->get ();                    \
            if (yyextra->inStream->eof()) result = YY_NULL;                   \
            else { buf[0] = c; result = 1; }                                  \
         }                                                                    \
      }

#     define YY_EXTRA_TYPE SXPARSER_TYPE *
#     define SXPARSER_PREFIX()     SXPARSER_TYPE ## _

      class YYLTYPE;
      class SXPARSER_TYPE;

      extern "C" bool SXPARSER_FUNC(SXPARSER_TYPE,_error) (YYLTYPE *, SXPARSER_TYPE *, const SxString&);
#     define SXPARSER_SEND(token) {                                           \
         SX_DBG_LEX("return " #token);                                        \
         return token; }
#     define SXPARSER_UNPUT(c)  {                                             \
         SX_DBG_LEX("unput ('" <<SxString(c).substitute('\n', "\\n")<<"')");  \
         if (c == '\n')  { yylloc->last_column = 999; --yylloc->first_line;   \
                           --yylloc->last_line; }                             \
         else            { --yylloc->last_column; }                           \
         yyunput (c, yytext, yyscanner); }
#     define SXPARSER_PUSH_STATE(state) {                                     \
         yyextra->pushLex(state,yylloc->first_line, yylloc->first_column);    \
         SX_DBG_LEX("push state: #" << state);                                \
         yy_push_state (state, yyscanner);  }
#     define SXPARSER_REPUSH_STATE() {                                        \
         yyextra->repushLex(yylloc->first_line, yylloc->first_column);        \
         yy_push_state (yyextra->getLexState(), yyscanner);  }
#     define SXPARSER_POP_STATE()       {                                     \
         yyextra->popLex(&yylloc->first_line, &yylloc->first_column);         \
         struct yyguts_t *yyg_ = (struct yyguts_t*)yyscanner;                 \
         if ( yyg_->yy_start_stack_ptr <= 0 )  {                              \
           SX_DBG_LEX("POP STATE: <INITIAL>");                                \
           BEGIN (INITIAL);                                                   \
         }  else  {                                                           \
            SX_DBG_MSG("LEX: pop state: #" << yy_top_state (yyscanner));      \
            yy_pop_state (yyscanner);                                         \
         }                                                                    \
      }
#     define SXPARSER_REPLACE_STATE(state)                                    \
         /* --- POP_STATE */                                                  \
         yyextra->popLex (&yylloc->first_line, &yylloc->first_column);        \
         struct yyguts_t *yyg_ = (struct yyguts_t*)yyscanner;                 \
         if ( yyg_->yy_start_stack_ptr <= 0 )  BEGIN (INITIAL);               \
         else yy_pop_state (yyscanner);                                       \
         /* --- PUSH STATE */                                                 \
         yyextra->pushLex (state,yylloc->first_line, yylloc->first_column);   \
         yy_push_state (state, yyscanner)
#     define SXPARSER_COLLECT_STATE()                                         \
         yyextra->collectLex()
#     define SXPARSER_APPEND_STATE(text)                                      \
         yyextra->appendLex(text)

#     define SXPARSER_SKIP_WHITESPACE() {                                     \
         while (*yytext == ' ' || *yytext == '\t')  {                         \
            ++yylloc->first_column; ++yytext;                                 \
         }                                                                    \
      }
#     define SXPARSER_FATAL(msg) { \
           SXPARSER_FUNC(SXPARSER_TYPE,_error) (yylloc, yyextra, msg);        \
           SxParserBase *scannerPtr = reinterpret_cast<SxParserBase *>(yyextra); \
           yyterminate ();                                                    \
      }
#     define SXPARSER_ERROR(msg) { \
           SXPARSER_FUNC(SXPARSER_TYPE,_error) (yylloc, yyextra, msg);        \
           SxParserBase *scannerPtr = reinterpret_cast<SxParserBase *>(yyextra); \
           if (!scannerPtr->resume())  yyterminate ();  \
      }
#     define SXPARSER_VALIDATE_INCLUDE(pattern)  {                            \
         yyextra->validateIncludes (pattern);                                 \
      }
#     define SXPARSER_PUSH_INCLUDE(pattern)  {                                \
         yyextra->setIncLoc (yylloc->first_line, yylloc->first_column,        \
                             yylloc->last_line,  yylloc->last_column);        \
         yyextra->pushIncludes (pattern,                                      \
                                yylloc->first_line,                           \
                                yylloc->first_column,                         \
                                yylloc->last_line,                            \
                                yylloc->last_column);                         \
         SXPARSER_REPLACE_STATE(INITIAL);                                     \
         yylloc->last_line = 1;                                               \
      }
#     define SXPARSER_POP_INCLUDE_OR_EXIT()  {                                \
         if (!yyextra->processNextInclude ())  {                              \
            if (yyextra->popInclude (&yylloc->first_line, &yylloc->first_column, \
                                     &yylloc->last_line,  &yylloc->last_column)) \
            {                                                                    \
               SXPARSER_POP_STATE ();                                            \
            } else { yyterminate (); }                                           \
         } else { cout << "XXX\n"; }\
      }

#     define SXPARSER_FOOTER                      \
         extern int SXPARSER_FUNC(SXPARSER_TYPE,_lex_init) (void *);          \
         extern int SXPARSER_FUNC(SXPARSER_TYPE,_set_extra) (void *, void *); \
         extern int SXPARSER_FUNC(SXPARSER_TYPE,_lex_destroy) (void *);       \
         extern int SXPARSER_FUNC(SXPARSER_TYPE,_debug);                      \
         void SXPARSER_TYPE::initScanner(bool trace=false) {                  \
            SXPARSER_FUNC(SXPARSER_TYPE,_lex_init)(&scannerPtr);              \
            SXPARSER_FUNC(SXPARSER_TYPE,_set_extra)(this, scannerPtr);        \
            SXPARSER_FUNC(SXPARSER_TYPE,_debug) = (trace == true);            \
         }                                                                    \
         void SXPARSER_TYPE::destroyScanner() {                               \
            SXPARSER_FUNC(SXPARSER_TYPE,_lex_destroy)(scannerPtr);            \
         }

#     define SX_DBG_LEX(msg)                                                  \
         SX_DBG_MSG ("LEX[" << lexerPtr->getLexState() << "]: "               \
                     << yylloc->first_line << "." << yylloc->first_column     \
                     << "-"                                                   \
                     << yylloc->last_line << "." << yylloc->last_column       \
                     << " " << msg                                            \
                     << "[" << SxString(yytext).substitute("\n", "\\n")       \
                                               .substitute("\t", "\\t")       \
                     << "]");

//#     define parserPtr ((SxQParser *)yyextra)

#  endif /* SX_PARSER_LEX_H */
#endif /* FLEX_SCANNER */



// --- included from yacc *.ypp
#ifdef yyparse
#  ifndef SX_PARSER_YACC_H
#     define SX_PARSER_YACC_H

#     ifndef SXPARSER_TYPE
#        error "Including SxParserBase.h from YACC requires SXPARSER_TYPE to be defined."
#     endif /* SX_PARSER_TYPE */

   int SXPARSER_FUNC(SXPARSER_TYPE,_lex)(YYSTYPE *lvalp, YYLTYPE *location, void *scannerPtr);

   class YYLTYPE;
   class SXPARSER_TYPE;
   extern "C" void SXPARSER_FUNC(SXPARSER_TYPE,_error)(YYLTYPE *loc, SXPARSER_TYPE *pPtr,
                      const SxString &msg)
   {
      SxParserBase *parser = reinterpret_cast<SxParserBase *>(pPtr);
      parser->errorHandler ( msg,
                 loc->first_line, loc->first_column,
                 loc->last_line,  loc->last_column);
   }

#  define SX_DBG_YACC(msg)                                                    \
      SX_DBG_MSG ("YACC: "                                                    \
                  << " " << msg);
#  define SXPARSER_ERROR_HANDLER() {                                          \
      if (parserPtr->resume ())  { yyclearin; }                               \
      else                       { YYABORT;   }                               \
   }

#  define _SxTag1(iLoc1)                                                      \
      parserPtr->getInFile() + ":"                                            \
      + iLoc1.first_line + "." + iLoc1.first_column + "-"                     \
      + iLoc1.last_line  + "." + (iLoc1.last_column-1)
#  define _SxTag2(iLoc1,iLoc2)                                                \
      parserPtr->getInFile() + ":"                                            \
      + iLoc1.first_line + "." + iLoc1.first_column + "-"                     \
      + iLoc2.last_line  + "." + (iLoc2.last_column-1)

# define SXTAG_YACC(...)  SX_VMACRO_FUNC (_SxTag, __VA_ARGS__)

#  define scannerPtr parserPtr->scannerPtr

#  endif /* SX_PARSER_YACC_H */
#endif /* yyparse */
