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

#ifndef _SX_G_QUERY_H_
#define _SX_G_QUERY_H_

#include <SxGQExprBase.h>
#include <SxN.h>
#include <SxGQExpr.h>
#include <SxGQExprList.h>
#include <SxGQAnd.h>
#include <SxGQOr.h>
#include <SxGQRelations.h>


/** \brief Graph Query Class

    \b SxGQuery = SPHInX Graph Query Class

    SxGQuery is the top level class used to parse
    then store the given query. It is also used
    to find matching graph patterns.

 */
class SxGQuery
{

   public:
      typedef SxPtr<SxUniqueList<ssize_t> > Selection;
      typedef SxPtr<SxList<Selection> >  SelSet;

      SxGQuery ();
      SxGQuery (const SxGQuery &in);
      SxGQuery (const sx::N &n);
      SxGQuery (const SxPtr<SxGQExprBase> &r);
      SxGQuery (const SxPtr<SxGQExprList> &r);
     ~SxGQuery ();

      SxGQuery &operator= (const SxGQuery &in);

      // operator= allows to assign a SxGQExprBase object to SxGQuery
      SxGQuery &operator= (const SxPtr<SxGQExprBase> &r);

      SxGQuery &operator= (const SxPtr<SxGQExprList> &r);

      void makeGraph ();

      // returns the matched selections of nodes' Idx.
      SelSet matchAll (const SxPtr<SxGraph<SxGProps> > &g,
                       const SxGraph<SxGProps>::Iterator &gIt);
      SelSet matchAll (const SxPtr<SxGraph<SxGProps> > &g);

      Selection match (const SxPtr<SxGraph<SxGProps> > &g,
                       const SxGraph<SxGProps>::Iterator &gIt);
      Selection  match (const SxPtr<SxGraph<SxGProps> > &g);

      bool matchOneIncomingOnce (const SxGraph<SxGQPattern>::ConstIterator &it,
                                 const SxGraph<SxGProps>::ConstIterator &gIt,
                                 Selection &sel,
                                 Selection &visited) const;

      bool matchAllIncomingOnce (const SxGraph<SxGQPattern>::ConstIterator &it,
                                 const SxGraph<SxGProps>::ConstIterator &gIt,
                                 Selection &sel,
                                 Selection &visited) const;


      // Single match function that evaluates current
      // pattern node and recursively checks all incoming
      // and outgoing edges too.
      bool evalCurrent (const SxGraph<SxGQPattern>::ConstIterator &it,
                        const SxGraph<SxGProps>::ConstIterator &gIt,
                        Selection &sel,
                        Selection &visited) const;


      // match all in data graph for the given incoming node in pattern graph
      bool matchAllOneIncoming (const SxGraph<SxGQPattern>::ConstIterator &it,
                                const SxGraph<SxGProps>::ConstIterator &gIt,
                                SelSet &sels,
                                Selection &visited) const;

      bool matchAllIncoming (const SxGraph<SxGQPattern>::ConstIterator &it,
                             const SxGraph<SxGProps>::ConstIterator &gIt,
                             SelSet &sels,
                             Selection &visited) const;

      // match all function that evaluates current
      // pattern node and recursively checks all incoming
      // and outgoing edges too.
      bool matchAllCurrent (const SxGraph<SxGQPattern>::ConstIterator &it,
                            const SxGraph<SxGProps>::ConstIterator &gIt,
                            SelSet &sels,
                            Selection &visited) const;

      inline void resizeSelections (SelSet &sels,
                                    ssize_t setSize, ssize_t selSize) const;

      inline SelSet crossSelections (SelSet &sels1, SelSet &sels2) const;

      // match all ordered-direct starting from given child node in data graph
      inline bool matchAllOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                              const SxGraph<SxGProps>::ConstIterator &gIt,
                              SelSet &sels, ssize_t childIdx,
                              Selection &visited) const;

      // match all child nodes of pattern graph ordered-direct (default)
      inline bool matchAllChildrenOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                      const SxGraph<SxGProps>::ConstIterator &gIt,
                                      SelSet &sels,
                                      Selection &visited) const;

      // match all child nodes of pattern graph unordered-indirect
      bool matchAllChildrenUI (const SxGraph<SxGQPattern>::ConstIterator &it,
                               const SxGraph<SxGProps>::ConstIterator &pIt,
                               SelSet &sels,
                               Selection &visited) const;


      // match all child nodes of pattern graph ordered-indirect
      inline bool matchAllChildrenOI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                      const SxGraph<SxGProps>::ConstIterator &gIt,
                                      SelSet &sels,
                                      Selection &visited) const;

      // match all child nodes of pattern graph unordered-direct
      inline bool matchAllChildrenUD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                      const SxGraph<SxGProps>::ConstIterator &gIt,
                                      SelSet &sels,
                                      Selection &visited) const;

      // recursively find one match starting at given child node unordered-direct
      bool matchOnceUDR (const SxGraph<SxGQPattern>::ConstIterator &it,
                         const SxGraph<SxGProps>::ConstIterator &gIt,
                         ssize_t childIdx, ssize_t exprChildIdx,
                         Selection &sel,
                         Selection &visited) const;

      // go through all child nodes and find one match unordered-direct
      inline bool matchOnceChildrenUD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                       const SxGraph<SxGProps>::ConstIterator &gIt,
                                       Selection &sel,
                                       Selection &visited) const;

      // go through all child nodes and find one match unordered-indirect
      inline bool matchOnceChildrenUI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                       const SxGraph<SxGProps>::ConstIterator &gIt,
                                       Selection &sel,
                                       Selection &visited) const;

      // find one match for ordered and direct
      bool matchOnceOD (const SxGraph<SxGQPattern>::ConstIterator &pIt,
                        const SxGraph<SxGProps>::ConstIterator &gIt,
                        Selection &sel,
                        ssize_t childIdx,
                        Selection &visited) const;

      // match once ordered and direct
      inline bool matchOnceChildrenOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                       const SxGraph<SxGProps>::ConstIterator &gIt,
                                       Selection &sel,
                                       Selection &visited) const;

      // find one match ordered-indirect
      inline bool matchOnceChildrenOI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                       const SxGraph<SxGProps>::ConstIterator &gIt,
                                       Selection &sel,
                                       Selection &visited) const;

      friend std::ostream &operator<< (std::ostream &s, const SxGQuery &q);

   protected:
      SxPtr<SxGQExprBase> rootExpr;
      SxGraph<SxGQPattern> patternGraph;
      SelSet selections;
};

      std::ostream &operator<< (std::ostream &s, const SxGQuery &q)
      {
         s << "\n\ndigraph G {\n";
         s << *(q.rootExpr);
         s << "\n}\n\n";
         return s;
      }


// -----------------------------------------------------------------------------

// returns SxGQExprBase representing an Equal condition
template<class T>
SxPtr<SxGQExprBase>
operator== (const sx::N &n_, const T &val_)
{
   SxPtr<SxGQExpr> e1 = SxPtr<SxGQExpr>::create (n_.getPropName (), SxVariant(val_),
                                                 SxGQExpr::OpType::Equal);
   return e1;
}

// returns SxGQExprBase representing a not-Equal condition
template<class T>
SxPtr<SxGQExprBase>
operator!= (const sx::N &n_, const T &val_)
{
   SxPtr<SxGQExpr> e1 = SxPtr<SxGQExpr>::create (n_.getPropName (), SxVariant(val_),
                                                 SxGQExpr::OpType::NotEqual);
   return e1;
}

// -----------------------------------------------------------------------------

// returns SxGQExprBase representing a binary AND condition
SxPtr<SxGQExprBase>
operator&& (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQAnd> expr1 = SxPtr<SxGQAnd>::create (e1, e2);
   return expr1;
}

// returns SxGQExprBase representing a binary OR condition
SxPtr<SxGQExprBase>
operator|| (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQOr> expr1 = SxPtr<SxGQOr>::create (e1, e2);
   return expr1;
}

// -----------------------------------------------------------------------------
// Ordered Direct

// joins two child nodes represented by two SxGQExprBase objects
SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::OrderedDirect);
   return ptr;
}

// joins two child nodes represented by SxGQExprList and SxGQExprBase
SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprList> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   e1->append (e2);
   return e1;
}

// joins two child nodes represented by SxGQExprBase and SxGQExprList
SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprList> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::OrderedDirect);
   return ptr;
}

// joins two child nodes represented by two SxGQExprList objects
SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprList> &e1,
            const SxPtr<SxGQExprList> &e2)
{
   e1->append (e2);
   return e1;
}

// -----------------------------------------------------------------------------
// Unordered Indirect

// joins two child nodes represented by two SxGQExpr objects
SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::UnorderedIndirect);
   return ptr;
}

// joins two child nodes represented by SxGQExprList and SxGQExprBase
SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprList> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   e1->append (e2);
   return e1;
}

// joins two child nodes represented by SxGQExprBase and SxGQExprList
SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprList> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::UnorderedIndirect);
   return ptr;
}

// joins two child nodes represented by two SxGQExprList objects
SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprList> &e1,
           const SxPtr<SxGQExprList> &e2)
{
   e1->append (e2);
   return e1;
}

// -----------------------------------------------------------------------------
// Unordered Direct

// joins two child nodes represented by two SxGQExpr objects
SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::UnorderedDirect);
   return ptr;
}

// joins two child nodes represented by SxGQExprList and SxGQExprBase
SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprList> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   e1->append (e2);
   return e1;
}

// joins two child nodes represented by SxGQExprBase and SxGQExprList
SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprList> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::UnorderedDirect);
   return ptr;
}

// joins two child nodes represented by two SxGQExprList objects
SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprList> &e1,
            const SxPtr<SxGQExprList> &e2)
{
   e1->append (e2);
   return e1;
}

// -----------------------------------------------------------------------------
// Ordered Indirect

// joins two child nodes represented by two SxGQExpr objects
SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::OrderedIndirect);
   return ptr;
}

// joins two child nodes represented by SxGQExprList and SxGQExprBase
SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprList> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   e1->append (e2);
   return e1;
}

// joins two child nodes represented by SxGQExprBase and SxGQExprList
SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprList> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::OrderedIndirect);
   return ptr;
}

// joins two child nodes represented by two SxGQExprList objects
SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprList> &e1,
           const SxPtr<SxGQExprList> &e2)
{
   e1->append (e2);
   return e1;
}

// -----------------------------------------------------------------------------

// joins SxGQExprBase(parent) with SxGQExprList(children) as SxGQRelations
SxPtr<SxGQExprBase>
operator^ (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprList> &elist)
{
   if (e1->isOp ()) {
      if (e1->getOp () == SxGQExprBase::OpType::Relation) {
         if (e1->getRightOp () == SxGQExprBase::OpType::Relation) {
            SxPtr<SxGQExprBase> p = e1->last ();
            SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (*elist);
            SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> > ptr =
            SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> >::create (p, child);
            e1->setLast (ptr);
            return e1;
         } else if (e1->getRightOp () == SxGQExprBase::OpType::Sibling) {
            SxPtr<SxList<SxPtr<SxGQExprBase> > > lst = e1->lasts ();
            auto lastsIt = lst->begin ();

            SxPtr<SxList<SxPtr<SxGQExprBase> > > newLasts =
            SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();

            for (; lastsIt != lst->end (); ++lastsIt) {
               // deep copy the sibling list, to avoid any problems due to shared obj
               SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (*elist);
               SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> > ptr =
               SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> >::create (*lastsIt, child);
               newLasts->append (ptr);
            }

            e1->setLasts (newLasts);
            return e1;
         }

      } else if (e1->getOp () == SxGQExprBase::OpType::Sibling) {
         SxPtr<SxList<SxPtr<SxGQExprBase> > > lst = e1->lasts ();
         auto lastsIt = lst->begin ();

         SxPtr<SxList<SxPtr<SxGQExprBase> > > newLasts =
         SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();

         for (; lastsIt != lst->end (); ++lastsIt) {
            // deep copy the sibling list, to avoid any problems due to shared obj
            SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (*elist);
            SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> > ptr =
            SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> >::create (*lastsIt, child);
            newLasts->append (ptr);
         }

         e1->setLasts (newLasts);
         return e1;
      } else {
         // currently no other operator has high enough
         // precedence to appear as parent
         SX_EXIT;
      }
   } else {
      SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> > ptr =
      SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> >::create (e1, elist);
      return ptr;
   }
   SX_EXIT;
   return e1;
}

// joins SxGQExprBase(parent) with SxGQExprBase(children) as SxGQRelations
SxPtr<SxGQExprBase>
operator^ (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   if (e1->isOp ()) {
      if (e1->getOp () == SxGQExprBase::OpType::Relation) {
         if (e1->getRightOp () == SxGQExprBase::OpType::Relation) {
            SxPtr<SxGQExprBase> p = e1->last ();
            SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (e2);
            SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> > ptr =
            SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> >::create (p, child);
            e1->setLast (ptr);
            return e1;
         } else if (e1->getRightOp () == SxGQExprBase::OpType::Sibling) {
            SxPtr<SxList<SxPtr<SxGQExprBase> > > lst = e1->lasts ();
            auto lastsIt = lst->begin ();

            SxPtr<SxList<SxPtr<SxGQExprBase> > > newLasts =
            SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();
            for (; lastsIt != lst->end (); ++lastsIt) {
               SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (e2);
               SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> > ptr =
               SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> >::create (*lastsIt, child);
               newLasts->append (ptr);
            }
            e1->setLasts (newLasts);
            return e1;
         }

      } else if (e1->getOp () == SxGQExprBase::OpType::Sibling) {
         SxPtr<SxList<SxPtr<SxGQExprBase> > > lst = e1->lasts ();
         auto lastsIt = lst->begin ();

         SxPtr<SxList<SxPtr<SxGQExprBase> > > newLasts =
         SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();
         for (; lastsIt != lst->end (); ++lastsIt) {
            SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (e2);
            SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> > ptr =
            SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> >::create (*lastsIt, child);
            newLasts->append (ptr);
         }
         e1->setLasts (newLasts);
         return e1;
      } else {
         // currently no other operator has high enough
         // precedence to appear as parent
         SX_EXIT;
      }
   } else {
      SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (e2);
      SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> > ptr =
      SxPtr<SxGQRelations<SxGQExprBase,SxGQExprList> >::create (e1, child);
      return ptr;
   }
   SX_EXIT;
   return e1;
}

// -----------------------------------------------------------------------------


#endif /*_SX_G_QUERY_H_*/
