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

#ifndef _SX_GQ_PATTERN_H_
#define _SX_GQ_PATTERN_H_

#include <SxPtr.h>
#include <SxGraph.h>
#include <SxGProps.h>
#include <stddef.h>

class SxGQExprBase;

/** \brief Node type for Search Pattern Graph

    \b SxGQPattern = SPHInX Search Pattern Class

     SxGQPattern class represent a node in the
     query pattern graph. A pattern graph is
     used in graph pattern matching where a
     sub-graph matching the pattern graph
     is to be found inside a given larger
     graph.

     Each SxGQPattern object contains a pointer
     to the expression that must evaluate to
     true, in order to successfully continue
     pattern matching.
 */

class SxGQPattern
{
   typedef SxPtr<SxUniqueList<ssize_t> > Selection;
   typedef SxPtr<SxList<Selection> >     SelSet;

   public:
      enum RelationType {OrderedDirect, OrderedIndirect,
                         UnorderedDirect, UnorderedIndirect};

      SxGQPattern ();
      SxGQPattern (ssize_t id_);
      SxGQPattern (ssize_t id_, const SxPtr<SxGQExprBase> &expr_);
      SxGQPattern (const SxGQPattern &in);
     ~SxGQPattern ();
      SxGQPattern &operator= (const SxGQPattern &in);
      bool operator== (const SxGQPattern &in) const;

      void setRelType (int rType);
      RelationType getRelType () const;
      ssize_t getId () const;
      void setExpr (const SxPtr<SxGQExprBase> &expr_);

      // eval function finds a single match that satisfy the
      // given expression according to specified operator
      bool eval (const SxGraph<SxGProps>::ConstIterator &it,
                 const Selection &sel) const;

      // matchAll function finds all possible matches that satisfy the
      // given expression according to specified operator
      bool matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                     const SelSet &sels) const;

   protected:
      ssize_t id;
      SxPtr<SxGQExprBase> expr;
      RelationType relType;
};

size_t sxHash (const SxGQPattern &in);
#endif /* _SX_ASG_NODE_H_ */
