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

#ifndef _SX_GQ_OR_H_
#define _SX_GQ_OR_H_

#include <SxGQExprBase.h>

/** \brief Graph Query Binary OR operator

    \b SxGQOr = SPHInX Graph Query Binary OR Operator Class

     SxGQOr class represents a binary OR
     operator between two nodes in a query.
     A binary OR evaluatese to true when
     either left or right operands evaluate
     to true.

\code
N1 || N2 || N3 ...

where N represents the class sx::N

\endcode
 */
class SxGQOr : public SxGQExprBase
{

   public:
      SxGQOr ();
      SxGQOr (const SxPtr<SxGQExprBase> &l, const SxPtr<SxGQExprBase> &r);
     ~SxGQOr ();

      bool eval (const SxGraph<SxGProps>::ConstIterator &it,
                 const Selection &sel) const override;
      bool matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                     const SelSet &sels) const override;

      SxPtr<SxList<SxPtr<SxGQExprBase> > > firsts () const override;
      SxPtr<SxList<SxPtr<SxGQExprBase> > > lasts () const override;
      SxPtr<SxGQExprBase> first () const override;
      SxPtr<SxGQExprBase> last () const override;
      void setLast (const SxPtr<SxGQExprBase> &p) override;
      void setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst) override;

      bool isOp () const override;
      OpType getOp () const override;
      OpType getRightOp () const override;
      size_t getHash () const override;

      void makeGraph (SxGraph<SxGQPattern> *g) const override;
      SxGQPattern getGraphNode () const override;

   protected:
      std::ostream &print (std::ostream &s) const override;

      SxPtr<SxGQExprBase> left;
      SxPtr<SxGQExprBase> right;
};

#endif /* _SX_GQ_OR_H_ */
