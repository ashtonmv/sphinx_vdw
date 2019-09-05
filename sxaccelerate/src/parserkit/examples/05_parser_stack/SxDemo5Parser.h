// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_Demo5_PARSER_H_
#define _SX_Demo5_PARSER_H_

#include <SxDemo5.h>
#include <SxParserBase.h>
#include <SxDemo5Ast.h>

/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....
    */
class SX_EXPORT_DEMO_5 SxDemo5Parser : public SxParserBase,
                                       public SxDemo5Ast
{
   public:

      SxDemo5Parser ();
     ~SxDemo5Parser ();

      // --- yacc stack interface
      void push (ssize_t);
      SxDemo5AstNode &pop ();
      SxDemo5AstNode &getCurrent ();
      SxDemo5AstNode &getParent ();


      // --- defined in SXPARSER_FOOTER
      virtual void initScanner (bool);
      virtual void destroyScanner ();


   protected:

      // --- yacc stack
      SxList<ssize_t> stack;

      virtual int parse ();

};

#endif /* _SX_Demo5_PARSER_H_ */
