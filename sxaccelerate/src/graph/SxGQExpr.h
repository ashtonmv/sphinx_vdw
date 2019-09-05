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

#ifndef _SX_GQ_EXPR_H_
#define _SX_GQ_EXPR_H_

#include <SxGQExprBase.h>

/** \brief Graph Query Basic Expression Class

    \b SxGQExpr = SPHInX Graph Query Basic Expression Class

     SxGQExpr class represents a basic
     expression in graph query. In the
     AST this class represents the leaf
     nodes. It matches left and right
     expression for equality or
     non-equality depending on OpType.

\code
#include <SxGQueryNode.h>
#include <SxGQExpr.h>

SxGQExpr expr = (N("propName") == "propVal");

or

SxGQExpr expr = (N("propName") != "propVal");

or

// expr will evaluate to true, if the
// node contains a property with name
// 'propName' with any value.

SxGQExpr expr = (N("propName").any());

\endcode
 */
class SxGQExpr : public SxGQExprBase
{
   public:
      SxGQExpr ();
      SxGQExpr (const SxString &propName_,
                const SxVariant &right_,
                const OpType &type_);
      SxGQExpr (const SxString &propName_,
                const SxVariant &right_,
                const SxString &capName_,
                const OpType &type_, bool isCap_);

      // eval function finds the match that satisfy the
      // given expression according to specified operator
      bool eval (const SxGraph<SxGProps>::ConstIterator &it,
                 const Selection &sel) const override;
      bool matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                     const SelSet &sels) const override;

      SxPtr<SxList<SxPtr<SxGQExprBase> > > firsts () const override;
      SxPtr<SxList<SxPtr<SxGQExprBase> > > lasts () const override;
      SxPtr<SxGQExprBase> first() const override;
      SxPtr<SxGQExprBase> last() const override;

      void setLast (const SxPtr<SxGQExprBase> &p) override;
      void setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst) override;

      bool isOp () const override;
      OpType getOp () const override;
      OpType getRightOp () const override;

      bool isCaptured () const;
      size_t getHash () const override;

      void makeGraph (SxGraph<SxGQPattern> *g) const override;
      SxGQPattern getGraphNode () const override;

   protected:

      std::ostream &print (std::ostream &s) const override;

      SxString propName;
      SxString captureName;
      bool captured;
      SxVariant right;
      ssize_t id;
      OpType opType;

};


#endif /*_SX_GQ_EXPR_H_*/
