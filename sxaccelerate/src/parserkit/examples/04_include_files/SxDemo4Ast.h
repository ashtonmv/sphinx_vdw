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

#ifndef _SX_DEMO_4_AST_H_
#define _SX_DEMO_4_AST_H_

#include <SxDemo4.h>
#include <SxString.h>
#include <SxGraph.h>

class SX_EXPORT_DEMO_4 SxDemo4AstNode
{
   public:

      ssize_t id;

      SxDemo4AstNode () { }
      SxDemo4AstNode (ssize_t) { }

      bool operator== (const SxDemo4AstNode &) const { return false; }
};

inline size_t sxHash (const SxDemo4AstNode &in)
{
   return SxHashFunction::hash (in.id);
}



/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author John Doe, johndoe@example.org */
class SX_EXPORT_DEMO_4 SxDemo4Ast

{
   public:

      SxDemo4Ast ();
     ~SxDemo4Ast ();

   protected:

     ssize_t idx;
     SxGraph<SxDemo4AstNode> ast;

};

#endif /* _SX_DEMO_4_AST_H_ */
