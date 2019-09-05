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

#ifndef _SX_Demo4_PARSER_H_
#define _SX_Demo4_PARSER_H_

#include <SxDemo4.h>
#include <SxParserBase.h>
#include <SxDemo4Ast.h>

/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....
    */
class SX_EXPORT_DEMO_4 SxDemo4Parser : public SxParserBase,
                                       public SxDemo4Ast
{
   public:

      SxDemo4Parser ();
     ~SxDemo4Parser ();

      // --- yacc stack interface
      void push (ssize_t);
      SxDemo4AstNode &pop ();
      SxDemo4AstNode &getCurrent ();
      SxDemo4AstNode &getParent ();


      // --- defined in SXPARSER_FOOTER
      virtual void initScanner (bool);
      virtual void destroyScanner ();


   protected:

      // --- yacc stack
      SxList<ssize_t> stack;

      virtual int parse ();

};

#endif /* _SX_Demo4_PARSER_H_ */
