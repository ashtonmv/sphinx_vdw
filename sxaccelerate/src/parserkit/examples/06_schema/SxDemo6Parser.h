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

#ifndef _SX_Demo6_PARSER_H_
#define _SX_Demo6_PARSER_H_

#include <SxDemo6.h>
#include <SxParserBase.h>
#include <SxDemo6Ast.h>

/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....
    */
class SX_EXPORT_DEMO_6 SxDemo6Parser : public SxParserBase,
                                       public SxDemo6Ast
{
   public:

      SxDemo6Parser ();
     ~SxDemo6Parser ();

      // --- yacc stack interface
      void push (ssize_t);
      SxDemo6AstNode &pop ();
      SxDemo6AstNode &getCurrent ();
      SxDemo6AstNode &getParent ();


      // --- defined in SXPARSER_FOOTER
      virtual void initScanner (bool);
      virtual void destroyScanner ();


   protected:

      // --- yacc stack
      SxList<ssize_t> stack;

      virtual int parse ();

};

#endif /* _SX_Demo6_PARSER_H_ */
