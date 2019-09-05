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

#include <SxGQOr.h>

SxGQOr::SxGQOr () { }
SxGQOr::SxGQOr (const SxPtr<SxGQExprBase> &l, const SxPtr<SxGQExprBase> &r)
                : left(l), right(r) { }
SxGQOr::~SxGQOr () { }

bool SxGQOr::eval (const SxGraph<SxGProps>::ConstIterator &it,
                   const Selection &sel) const
{
   SX_CHECK (it.isValid ());
   SX_CHECK (left.getPtr ());
   SX_CHECK (right.getPtr ());

   bool res = left->eval (it, sel);
   if (!res) {
      res = right->eval (it, sel);
   }
   return res;
}

bool SxGQOr::matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                       const SelSet &sels) const
{
   SX_CHECK (it.isValid ());
   SX_CHECK (left.getPtr ());
   SX_CHECK (right.getPtr ());

   bool res = left->matchAll (it, sels);
   if (!res) {
      res = right->matchAll (it, sels);
   }
   return res;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQOr::firsts () const
{
   SX_EXIT;
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   eList->append (first ());
   return eList;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQOr::lasts () const
{
   if (right->isOp ())
      return right->lasts ();
   else {
      SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
         SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();
      eList->append (right);
      return eList;
   }
}

SxPtr<SxGQExprBase> SxGQOr::first () const
{
   SX_CHECK (left.getPtr ());
   if (left->isOp())
      return left->first ();
   else
      return left;
}

SxPtr<SxGQExprBase> SxGQOr::last () const
{
   SX_CHECK (right.getPtr ());
   if (right->isOp ())
      return right->last ();
   else
      return right;
}

void SxGQOr::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_CHECK (right.getPtr ());
   if (right->isOp ())
      right->setLast (p);
   else
      right = p;
}

void SxGQOr::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_CHECK (right.getPtr ());
   SX_CHECK (right->isOp ());
   right->setLasts (lst);
}

bool SxGQOr::isOp () const
{
   return false;
}

SxGQExprBase::OpType SxGQOr::getOp () const
{
   return SxGQExprBase::OpType::OR;
}

SxGQExprBase::OpType SxGQOr::getRightOp () const
{
   return right->getRightOp ();
}


size_t SxGQOr::getHash () const
{
   return ( left->getHash ()
         + SxHashFunction::hash (SxString("||"))
         + right->getHash () );
}

void SxGQOr::makeGraph (SxGraph<SxGQPattern> *g) const
{
   SxPtr<SxGQOr> expr = SxPtr<SxGQOr>::create (*this);
   SxGQPattern n = getGraphNode ();
   g->createNode (n);
   g->begin (n)->setExpr (expr);
}

SxGQPattern SxGQOr::getGraphNode () const
{
   SX_CHECK (!left->isOp ());
   SxPtr<SxGQOr> expr = SxPtr<SxGQOr>::create (*this);
   SxGQPattern n((ssize_t)getHash ());
   n.setExpr (expr);
   return n;
}

std::ostream &SxGQOr::print (std::ostream &s) const
{
   s << "\"";
   s << *left;
   s << "||";
   s << *right;
   s << "\"";
   return s;
}
