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

#include <SxGQExpr.h>

SxGQExpr::SxGQExpr ()
   : captured(false),
     id(-1),
     opType(OpType::None)
{
   // empty
}

SxGQExpr::SxGQExpr (const SxString &propName_,
                    const SxVariant &right_,
                    const OpType &type_)
{
   SX_CHECK (propName_!= "");
   propName = propName_;
   right    = right_;
   opType   = type_;
   id       = (ssize_t)getHash ();
   captured = false;
}

SxGQExpr::SxGQExpr (const SxString &propName_,
                    const SxVariant &right_,
                    const SxString &capName_,
                    const OpType &type_,
                    bool isCap_)
{
   SX_CHECK (propName_!= "");
   propName     = propName_;
   captureName  = capName_;
   captured     = isCap_;
   right        = right_;
   opType       = type_;
   id           = (ssize_t)getHash ();
}

// eval function finds the match that satisfy the
// given expression according to specified operator
bool SxGQExpr::eval (const SxGraph<SxGProps>::ConstIterator &it,
                     const Selection &sel) const
{
   if (!it.isValid ()) return false;
   bool res = false;
   switch (opType) {
      case SxGQExpr::OpType::Any:
         if (it->hasProperty(propName)) {
            res = true;
            sel->append (it.getIdx ());
         }
         return res;
      case SxGQExpr::OpType::Equal:
         res = (it->hasProperty (propName) && it->getProperty(propName) == right);

         if (res == true) sel->append (it.getIdx ());
         if (res == false)
            SX_DBG_MSG ("SxGQExpr::Equal 'False' for Idx: "+
                  SxString(it.getIdx ())+
                  " == "+right.toString ());
         return res;
      case SxGQExpr::OpType::NotEqual:
         res = !(it->hasProperty (propName) && it->getProperty(propName) == right);
         if (res == true) sel->append (it.getIdx ());
         if (res == false)
            SX_DBG_MSG ("SxGQExpr::NotEqual 'False' for Idx: "+
                  SxString(it.getIdx ())+
                  " == "+right.toString ());
         return res;
      default:
         return false;
   }
   return false;
}

bool SxGQExpr::matchAll (const SxGraph<SxGProps>::ConstIterator &it,
                         const SelSet &sels) const
{
   Selection sel = Selection::create ();
   if ( eval (it, sel) ) {

      SX_CHECK (sel->getSize ()>0, sel->getSize());

      // if empty, then add as first
      if (sels->getSize () == 0) {
         sels->append (sel);
         return true;
      }

      // cross of n x 1 = n
      for (auto it1 = sels->begin ();it1 != sels->end (); ++it1) {
         (*it1)->append (*sel);
      }
      return true;
   }
   return false;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQExpr::firsts () const
{
   SX_EXIT;
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   eList->append (first ());
   return eList;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQExpr::lasts () const
{
   SX_EXIT;
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   eList->append (last ());
   return eList;
}

SxPtr<SxGQExprBase> SxGQExpr::first() const
{
   SX_EXIT;
   SxPtr<SxGQExpr> e = SxPtr<SxGQExpr>::create (*this);
   return e;
}

SxPtr<SxGQExprBase> SxGQExpr::last() const
{
   SX_EXIT;
   SxPtr<SxGQExpr> e = SxPtr<SxGQExpr>::create (*this);
   return e;
}

void SxGQExpr::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_UNUSED (p);
   SX_EXIT;
}

void SxGQExpr::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_UNUSED (lst);
   SX_EXIT;
}

bool SxGQExpr::isOp () const
{
   return false;
}

SxGQExprBase::OpType SxGQExpr::getOp () const
{
   return SxGQExprBase::OpType::None;
}

SxGQExprBase::OpType SxGQExpr::getRightOp () const
{
   return SxGQExprBase::OpType::None;
}

bool SxGQExpr::isCaptured () const
{
   return captured;
}

size_t SxGQExpr::getHash () const
{
   size_t res = ( SxHashFunction::hash (propName)
                + SxHashFunction::hash (opType) );
   if (right.isInitialized ())
      res += SxHashFunction::hash (right.toString());
   return res;
}

void SxGQExpr::makeGraph (SxGraph<SxGQPattern> *g) const
{
   SxPtr<SxGQExpr> expr = SxPtr<SxGQExpr>::create (*this);
   SxGQPattern n = getGraphNode ();
   g->createNode (n);
   g->begin (n)->setExpr (expr);
}

SxGQPattern SxGQExpr::getGraphNode () const
{
   SxPtr<SxGQExpr> expr = SxPtr<SxGQExpr>::create (*this);
   SxGQPattern n((ssize_t)getHash ());
   n.setExpr (expr);
   return n;
}

std::ostream &SxGQExpr::print (std::ostream &s) const
{
   s << "\'N(" << propName << ")";
   if (opType == SxGQExprBase::OpType::NotEqual)  s << "!=";
   else                                           s << "==";
   if (opType == SxGQExprBase::OpType::Any)  s << "?";
   else                                      s << right.toString ();
   s << "\'";
   return s;
}
