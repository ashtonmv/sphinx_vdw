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

#ifndef _SX_DEMO_6_AST_H_
#define _SX_DEMO_6_AST_H_

#include <SxDemo6.h>
#include <SxString.h>
#include <SxGraph.h>
#include <SxVariant.h>

class SX_EXPORT_DEMO_6 SxDemo6AstNode
{
   public:

      ssize_t id;

      enum Type { Group,
                  KeyVal,
                  StringVal,
                  NumberVal,
                  ArrayVal,
                  TrueVal,
                  FalseVal,
                  NullVal
      } type;

      SxVariant data;
      int line0, col0, line1, col1;

      SxDemo6AstNode ();
      SxDemo6AstNode (ssize_t id_, const Type &type_=Group);
      SxDemo6AstNode (ssize_t id_, const Type &type_, const SxVariant &data_);
     ~SxDemo6AstNode ();

      bool operator== (const SxDemo6AstNode &) const;
};

std::ostream &operator<< (std::ostream &s, const SxDemo6AstNode &in);

inline size_t sxHash (const SxDemo6AstNode &in)
{
   return SxHashFunction::hash (in.id);
}



/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author John Doe, johndoe@example.org */
class SX_EXPORT_DEMO_6 SxDemo6Ast

{
   public:

      SxDemo6Ast ();
     ~SxDemo6Ast ();

      ssize_t addNode (const SxDemo6AstNode::Type &);
      ssize_t addNode (const SxDemo6AstNode::Type &,
                       const SxVariant &);
      void addEdge (ssize_t, ssize_t);
      void prependEdge (ssize_t, ssize_t);

      SxDemo6AstNode &getNode (ssize_t);

      const SxGraph<SxDemo6AstNode> &getAst () const;

      void printAst () const;

   protected:

     ssize_t idx;
     SxGraph<SxDemo6AstNode> ast;

};

#endif /* _SX_DEMO_6_AST_H_ */
