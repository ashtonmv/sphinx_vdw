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

#include <SxGQExprList.h>

SxGQExprList::SxGQExprList () : sibType(OrderedDirect) {}

SxGQExprList::SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                            SxGQExprList::SiblingType sType)
                            : sibType(sType)
{
   exprList << e1;
}

SxGQExprList::SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                            const SxPtr<SxGQExprBase> &e2,
                            SxGQExprList::SiblingType sType)
                            : sibType(sType)
{
   exprList << e1 << e2;
}

SxGQExprList::SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                            const SxPtr<SxGQExprList> &e2,
                            SxGQExprList::SiblingType sType)
                            : sibType(sType)
{
   SX_CHECK (sibType == e2->sibType);
   exprList << e1 << e2->exprList;
}

SxGQExprList::SxGQExprList (const SxPtr<SxGQExprList> &e1,
                            const SxPtr<SxGQExprBase> &e2,
                            SxGQExprList::SiblingType sType)
                            : sibType(sType)
{
   SX_CHECK (sibType == e1->sibType);
   exprList << e1->exprList << e2;
}

SxGQExprList::~SxGQExprList () { }

void SxGQExprList::append (const SxGQExprList &lst)
{
   exprList << lst.exprList;
}

void SxGQExprList::append (const SxPtr<SxGQExprBase> &e)
{
   SX_CHECK (e.getPtr ());
   exprList << e;
}

void SxGQExprList::append (const SxPtr<SxGQExprList> &lst)
{
   SX_CHECK (lst.getPtr ());
   append (*lst);
}

SxGQExprList::SiblingType SxGQExprList::getSibType () const
{
   return sibType;
}

const SxList<SxPtr<SxGQExprBase> > &SxGQExprList::getExprList () const
{
   return exprList;
}

bool SxGQExprList::matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                             const SelSet &sels) const
{
   SX_UNUSED (it, sels);
   SX_EXIT;
   return false;
}

bool SxGQExprList::eval (const SxGraph<SxGProps>::ConstIterator &it,
                         const Selection &sel) const
{
   SX_UNUSED (it,sel);
   SX_EXIT;
   return false;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQExprList::firsts () const
{
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   for (auto listIt = exprList.begin();listIt != exprList.end(); ++listIt)
   {
      eList->append (*listIt);
   }
   return eList;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQExprList::lasts () const
{
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   for (auto listIt = exprList.begin();listIt != exprList.end(); ++listIt)
   {
      if ((*listIt)->isOp ())
         eList->append ((*listIt)->last ());
      else
         eList->append (*listIt);
   }
   return eList;
}

void SxGQExprList::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   SX_CHECK (exprList.getSize () == lst->getSize());
   auto lastsIt = lst->begin ();
   for (auto listIt = exprList.begin();listIt != exprList.end(); ++listIt)
   {
      if ((*listIt)->isOp ())
         (*listIt)->setLast (*lastsIt);
      else
         *listIt = *lastsIt;
      ++lastsIt;
   }
}

SxPtr<SxGQExprBase> SxGQExprList::first () const
{
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   if (exprList.first ()->isOp())
      return exprList.first ()->first ();
   else
      return exprList.first ();
}

SxPtr<SxGQExprBase> SxGQExprList::last () const
{
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   if (exprList.last ()->isOp ())
      return exprList.last ()->last ();
   else
      return exprList.last ();
}
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

void SxGQExprList::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   if (exprList.last ()->isOp ())
      exprList.last ()->setLast (p);
   else
      exprList.last () = p;
}

bool SxGQExprList::isOp () const
{
   return true;
}

SxGQExprBase::OpType SxGQExprList::getOp () const
{
   return SxGQExprBase::OpType::Sibling;
}

SxGQExprBase::OpType SxGQExprList::getRightOp () const
{
   if (exprList.getSize () > 0)
      return SxGQExprBase::OpType::Sibling;
   else
      if (exprList.last ()->isOp ())
         return exprList.last ()->getOp ();
      else
         return SxGQExprBase::OpType::None;
}

size_t SxGQExprList::getHash () const
{
   SX_EXIT;
   return 0;
}

void SxGQExprList::makeGraph (SxGraph<SxGQPattern> *g) const
{
   for (auto listIt = exprList.begin();listIt != exprList.end(); ++listIt)
   {
      if ((*listIt)->isOp ()) {
         (*listIt)->makeGraph (g);
      } else {
         g->createNode ((*listIt)->getGraphNode ());
      }
   }
}

SxGQPattern SxGQExprList::getGraphNode () const
{
   SX_EXIT;
}

std::ostream &SxGQExprList::print (std::ostream &s) const
{
   for (auto listIt = exprList.begin();listIt != exprList.end(); ++listIt)
   {
      s << "\n";
      s << *(*listIt);
      if (!(*listIt)->isOp())  s << ";\n";
   }
   return s;
}

