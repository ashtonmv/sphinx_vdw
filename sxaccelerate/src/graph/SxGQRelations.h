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

#ifndef _SX_GQ_RELATIONS_H_
#define _SX_GQ_RELATIONS_H_

#include <SxGQExprBase.h>

/** \brief Graph Query Parent-Child relationship expression

    \b SxGQRelations = SPHInX Graph Query Parent-Child Relationship Class

     SxGQRelations class stores the
     representation for relationship
     between parent and child nodes.

\code

N1 ^ ( N2 < N3 )

where N1 is the parent of N2 and N3

\endcode
 */
template<class FromType,class ToType>
class SxGQRelations : public SxGQExprBase
{
   public:
      SxGQRelations ();
      SxGQRelations (const SxPtr<FromType> &frm_, size_t min_=1, size_t max_=1);
      SxGQRelations (const SxPtr<FromType> &frm_, const SxPtr<ToType> &to_,
                     size_t min_=1, size_t max_=1);
     ~SxGQRelations ();

      bool matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                     const SelSet &sels) const override;

      // eval function finds match for given parent/child
      // pattern
      bool eval (const SxGraph<SxGProps>::ConstIterator &it,
                 const Selection &sel) const override;

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

      SxPtr<FromType> from;
      SxPtr<ToType>   to;
      size_t minDist;
      size_t maxDist;

};

#include <SxGQRelations.hpp>

#endif /*_SX_GQ_RELATIONS_H_*/
