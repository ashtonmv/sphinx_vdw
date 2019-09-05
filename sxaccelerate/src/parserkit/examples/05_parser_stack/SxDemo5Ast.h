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

#ifndef _SX_DEMO_5_AST_H_
#define _SX_DEMO_5_AST_H_

#include <SxDemo5.h>
#include <SxString.h>
#include <SxGraph.h>
#include <SxVariant.h>

class SX_EXPORT_DEMO_5 SxDemo5AstNode
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

      SxDemo5AstNode ();
      SxDemo5AstNode (ssize_t id_, const Type &type_=Group);
      SxDemo5AstNode (ssize_t id_, const Type &type_, const SxVariant &data_);
     ~SxDemo5AstNode ();

      bool operator== (const SxDemo5AstNode &) const;
};

std::ostream &operator<< (std::ostream &s, const SxDemo5AstNode &in);

inline size_t sxHash (const SxDemo5AstNode &in)
{
   return SxHashFunction::hash (in.id);
}



/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author John Doe, johndoe@example.org */
class SX_EXPORT_DEMO_5 SxDemo5Ast

{
   public:

      SxDemo5Ast ();
     ~SxDemo5Ast ();

      ssize_t addNode (const SxDemo5AstNode::Type &);
      ssize_t addNode (const SxDemo5AstNode::Type &,
                       const SxVariant &);
      void addEdge (ssize_t, ssize_t);
      void prependEdge (ssize_t, ssize_t);

      SxDemo5AstNode &getNode (ssize_t);

      const SxGraph<SxDemo5AstNode> &getAst () const;

      void printAst () const;

   protected:

     ssize_t idx;
     SxGraph<SxDemo5AstNode> ast;

};

#endif /* _SX_DEMO_5_AST_H_ */
