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

template<class FT,class TT>
SxGQRelations<FT,TT>::SxGQRelations () {}

template<class FT,class TT>
SxGQRelations<FT,TT>::SxGQRelations (const SxPtr<FT> &frm_,
                                     size_t min_, size_t max_)
                                   : from(frm_), minDist(min_), maxDist(max_)
                                     { }

template<class FT,class TT>
SxGQRelations<FT,TT>::SxGQRelations (const SxPtr<FT> &frm_,
                                     const SxPtr<TT> &to_,
                                     size_t min_, size_t max_)
                                     : from(frm_), to(to_), minDist(min_),
                                       maxDist(max_) { }

template<class FT,class TT>
SxGQRelations<FT,TT>::~SxGQRelations () { }

template<class FT,class TT>
bool SxGQRelations<FT,TT>::matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                                     const SelSet &sels) const
{
   SX_UNUSED (it, sels);
   SX_EXIT;
   return false;
}

// eval function finds match for given parent/child
// pattern
template<class FT,class TT>
bool SxGQRelations<FT,TT>::eval (const SxGraph<SxGProps>::ConstIterator &it,
                                 const Selection &sel) const
{
   SX_UNUSED (it, sel);
   SX_EXIT;
   return false;
}

template<class FT,class TT>
SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQRelations<FT,TT>::firsts () const
{
   SX_CHECK (from.getPtr ());
   return from->firsts ();
}

template<class FT,class TT>
SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQRelations<FT,TT>::lasts () const
{
   SX_CHECK (to.getPtr ());
   if (to->isOp ())
      return to->lasts ();
   else {
      SxPtr<SxList<SxPtr<SxGQExprBase> > > res =
         SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();
      res->append (to);
      return res;
   }
}

template<class FT,class TT>
SxPtr<SxGQExprBase> SxGQRelations<FT,TT>::first () const
{
   SX_CHECK (from.getPtr ());
   if (from->isOp ())
      return from->first ();
   else
      return from;
}

template<class FT,class TT>
SxPtr<SxGQExprBase> SxGQRelations<FT,TT>::last () const
{
   SX_CHECK (to.getPtr ());
   if (to->isOp ())
      return to->last ();
   else
      return to;
}

template<class FT,class TT>
void SxGQRelations<FT,TT>::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_CHECK (to.getPtr ());
   if (to->isOp ())
      to->setLast (p);
   else
      to = p;
}

template<class FT,class TT>
void SxGQRelations<FT,TT>::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_CHECK (to.getPtr ());
   SX_CHECK (to->isOp ());
   to->setLasts (lst);
}

template<class FT,class TT>
bool SxGQRelations<FT,TT>::isOp () const
{
   return true;
}

template<class FT,class TT>
SxGQExprBase::OpType SxGQRelations<FT,TT>::getOp () const
{
   return SxGQExprBase::OpType::Relation;
}

template<class FT,class TT>
SxGQExprBase::OpType SxGQRelations<FT,TT>::getRightOp () const
{
   return to->getRightOp ();
}


template<class FT,class TT>
size_t SxGQRelations<FT,TT>::getHash () const
{
   SX_EXIT;
   return 0;
}

template<class FT,class TT>
void SxGQRelations<FT,TT>::makeGraph (SxGraph<SxGQPattern> *g) const
{
   SxPtr<SxGQExprBase> last;
   if (from->isOp ()) {
      from->makeGraph (g);
      last = from->last ();
   } else {
      last = from;
   }

   SxGQPattern parentNode  = last->getGraphNode ();

   g->createNode (parentNode);

   // if node is already in graph, the new one will not be
   // inserted, hence to make the changes to the node present
   // in the graph use the returned iterator.

   g->begin(parentNode)->setRelType (to->getSibType ());
   g->begin(parentNode)->setExpr (last);

   for (auto listIt = to->getExprList().begin();
         listIt != to->getExprList().end(); ++listIt)
   {
      SxPtr<SxGQExprBase> child;
      if((*listIt)->isOp ())
         child = (*listIt)->first ();
      else
         child = (*listIt);

      SxGQPattern childNode = child->getGraphNode ();

      g->createNode (childNode);

      g->begin (childNode)->setExpr (child);

      g->createEdge (parentNode, childNode);

      if((*listIt)->isOp ())
         (*listIt)->makeGraph (g);
   }

}

template<class FT,class TT>
SxGQPattern SxGQRelations<FT,TT>::getGraphNode () const
{
   SX_EXIT;
}

template<class FT,class TT>
std::ostream &SxGQRelations<FT,TT>::print (std::ostream &s) const
{
   if (from->isOp ())
      s << *from;

   for (auto listIt = to->getExprList().begin();
         listIt != to->getExprList().end(); ++listIt)
   {
      s << "\n";
      if(from->isOp ())
         s << *(from->last ());
      else
         s << *from;
      s << " -> ";
      if((*listIt)->isOp ())
         s << *((*listIt)->first ());
      else
         s << *(*listIt);

      s << ";\n";
      if((*listIt)->isOp ())
         s << *(*listIt);
   }
   return s;
}

