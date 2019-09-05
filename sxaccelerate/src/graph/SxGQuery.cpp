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

#include <SxGQuery.h>

SxGQuery::SxGQuery () {}

SxGQuery::SxGQuery (const SxGQuery &in)
{
   this->rootExpr = in.rootExpr;
   makeGraph ();
}

SxGQuery::SxGQuery (const sx::N &n)
{
   SxPtr<SxGQExprBase> e = n;
   this->rootExpr = e;
   makeGraph ();
}

SxGQuery::SxGQuery (const SxPtr<SxGQExprBase> &r)
{
   SX_CHECK (r.getPtr ());
   rootExpr = r;
   makeGraph ();
}

SxGQuery::SxGQuery (const SxPtr<SxGQExprList> &r)
{
   SX_CHECK (r.getPtr ());
   rootExpr = r;
   makeGraph ();
}

SxGQuery::~SxGQuery () { }

SxGQuery &SxGQuery::operator= (const SxGQuery &in)
{
   if (this == &in) return *this;
   this->rootExpr = in.rootExpr;
   patternGraph = SxGraph<SxGQPattern>();
   makeGraph ();
   return *this;
}

// operator= allows to assign a SxGQExprBase object to SxGQuery
SxGQuery &SxGQuery::operator= (const SxPtr<SxGQExprBase> &r)
{
   SX_CHECK (r.getPtr ());
   this->rootExpr = r;
   patternGraph = SxGraph<SxGQPattern>();
   makeGraph ();
   return *this;
}

SxGQuery &SxGQuery::operator= (const SxPtr<SxGQExprList> &r)
{
   SX_CHECK (r.getPtr ());
   this->rootExpr = r;
   patternGraph = SxGraph<SxGQPattern>();
   makeGraph ();
   return *this;
}

void SxGQuery::makeGraph ()
{
   if (rootExpr->isOp ()) {
      rootExpr->makeGraph (&patternGraph);
   } else {
      patternGraph.createNode (rootExpr->getGraphNode ());
   }
}

// returns the matched selections of nodes' Idx.
SxGQuery::SelSet SxGQuery::matchAll (const SxPtr<SxGraph<SxGProps> > &g,
                           const SxGraph<SxGProps>::Iterator &gIt)
{
   selections = SxGQuery::SelSet::create ();
   SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();
   bool res = false, finalRes = false;
   SxGraph<SxGProps>::ConstIterator tIt;
   // list for visited pattern nodes overall
   SxGQuery::Selection gVisited = SxGQuery::Selection::create ();
   SxGQuery::Selection visited = SxGQuery::Selection::create ();
   SxGQuery::Selection fullVisited = SxGQuery::Selection::create ();

   auto it = gIt;
   // this iterator goes over all nodes
   // irrespective of edges. It allows to
   // cover all connected components of the
   // graph.
   auto patternIt = patternGraph.begin (sx::Undefined, (*patternGraph.begin ()));
   while (patternIt.isValid ()) {
      while (it.isValid ()) {
         SX_CHECK (selections.getPtr ());
         SX_CHECK (tmpSels.getPtr ());
         visited = SxGQuery::Selection::create ();
         // it is important to start with fresh iterator
         // in order to cover all possible paths to a node
         tIt = g->begin (*it);
         res = matchAllCurrent (patternGraph.begin (*patternIt), tIt, tmpSels, visited);
         if (res == true) {
            // at least one match found during this loop
            finalRes = true;
            // store the full matched set of pattern nodes
            fullVisited = visited;
            selections->append (std::move(*tmpSels));
         }
         ++it;
      }
      if (!finalRes) {
         return SxGQuery::SelSet::create ();
      }
      gVisited->append (*fullVisited);
      finalRes = false;
      it = gIt;
      // skip the visited pattern nodes
      while (patternIt.isValid () && gVisited->contains (patternIt->getId ()))
         ++patternIt;
   }
   return selections;
}

SxGQuery::SelSet SxGQuery::matchAll (const SxPtr<SxGraph<SxGProps> > &g)
{
   return matchAll (g, g->begin (sx::Undefined, (*g->begin ())));
}

SxGQuery::Selection SxGQuery::match (const SxPtr<SxGraph<SxGProps> > &g,
                           const SxGraph<SxGProps>::Iterator &gIt)
{

   bool res = false;
   SxGQuery::Selection sel = SxGQuery::Selection::create ();
   SxGQuery::Selection visited = SxGQuery::Selection::create ();
   auto it = gIt;
   SxGraph<SxGProps>::ConstIterator tIt;

   auto patternIt = patternGraph.begin (sx::Undefined, (*patternGraph.begin ()));
   while (patternIt.isValid ()) {
      while (it.isValid ()) {
         tIt = g->begin (*it);
         res = evalCurrent (patternGraph.begin (*patternIt), tIt, sel, visited);
         if (res == true)  break;
         ++it;
      }
      if (!res) {
         sel->removeAll ();
         break;
      }
      res = false;
      it = gIt;
      // skip the visited pattern nodes
      while (patternIt.isValid () && visited->contains (patternIt->getId ()))
         ++patternIt;
   }

   return sel;
}

SxGQuery::Selection SxGQuery::match (const SxPtr<SxGraph<SxGProps> > &g)
{
   return match (g, g->begin (sx::Undefined, (*g->begin ())));
}


bool SxGQuery::matchOneIncomingOnce (const SxGraph<SxGQPattern>::ConstIterator &it,
                                     const SxGraph<SxGProps>::ConstIterator &gIt,
                                     SxGQuery::Selection &sel,
                                     SxGQuery::Selection &visited) const
{
   bool finalRes = false;
   for (ssize_t i = 0; i < gIt.getSizeIn (); ++i) {
      if(evalCurrent (it, gIt.in (i), sel, visited)) {
         finalRes = true; // atleast one match found
         break;
      }
   }
   return finalRes;
}

bool SxGQuery::matchAllIncomingOnce (const SxGraph<SxGQPattern>::ConstIterator &it,
                                     const SxGraph<SxGProps>::ConstIterator &gIt,
                                     SxGQuery::Selection &sel,
                                     SxGQuery::Selection &visited) const
{
   SX_CHECK (it.isValid ());
   SX_CHECK (gIt.isValid ());

   SxGQuery::Selection resSel  = SxGQuery::Selection::create ();
   bool res = true, visitable = false;

   for (ssize_t i = 0; i < it.getSizeIn (); ++i) {
      auto inIt = it.in (i);
      if (!visited->contains (inIt->getId ())) {
         visitable = true;
         res = matchOneIncomingOnce (inIt, gIt, resSel, visited);
         if (!res)  break;
      }
   }
   if (res && visitable) {
      resSel->append (*sel);
      sel = resSel;
   }
   return res;
}

// Single match function that evaluates current
// pattern node and recursively checks all incoming
// and outgoing edges too.
bool SxGQuery::evalCurrent (const SxGraph<SxGQPattern>::ConstIterator &it,
                            const SxGraph<SxGProps>::ConstIterator &gIt,
                            SxGQuery::Selection &sel,
                            SxGQuery::Selection &visited) const
{
   SX_CHECK (it.isValid ());
   SX_CHECK (gIt.isValid ());

   // if this node has already been visited
   // still evaluate it so that calling function
   // can know success/failure. But don't append
   // it to the result.
   if (visited->contains (it->getId ())) {
      SxGQuery::Selection tmpSel = SxGQuery::Selection::create ();
      return it->eval (gIt, tmpSel);
   }
   bool res = false;
   ssize_t siz = sel->getSize ();
   res = it->eval (gIt, sel);

   if (res) {
      visited->append (it->getId ());
      res = matchAllIncomingOnce (it, gIt, sel, visited);
   }

   if (res && it.getSizeOut () > 0) {
      switch(it->getRelType ()) {
         case SxGQPattern::RelationType::OrderedDirect:
            res = matchOnceChildrenOD (it, gIt, sel, visited);
            break;
         case SxGQPattern::RelationType::OrderedIndirect:
            res = matchOnceChildrenOI (it, gIt, sel, visited);
            break;
         case SxGQPattern::RelationType::UnorderedDirect:
            res = matchOnceChildrenUD (it, gIt, sel, visited);
            break;
         case SxGQPattern::RelationType::UnorderedIndirect:
            res = matchOnceChildrenUI (it, gIt, sel, visited);
            break;
      }
   }
   if (!res)  sel->resize (siz);
   return res;
}

// match all in data graph for the given incoming node in pattern graph
bool SxGQuery::matchAllOneIncoming (const SxGraph<SxGQPattern>::ConstIterator &it,
                                    const SxGraph<SxGProps>::ConstIterator &gIt,
                                    SxGQuery::SelSet &sels,
                                    SxGQuery::Selection &visited) const
{
   SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
   bool finalRes = false;
   SxGQuery::Selection tmpVisited;
   for (ssize_t i = 0; i < gIt.getSizeIn (); ++i) {
      SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();
      tmpVisited = SxGQuery::Selection::create (*visited);
      if(matchAllCurrent (it, gIt.in (i), tmpSels, tmpVisited)) {
         finalRes = true; // atleast one match found
         resSet->append (std::move(*crossSelections (tmpSels, sels)) );
      }
   }
   visited = tmpVisited; // update visited list
   if (finalRes)  sels = resSet;
   return finalRes;
}


bool SxGQuery::matchAllIncoming (const SxGraph<SxGQPattern>::ConstIterator &it,
                                 const SxGraph<SxGProps>::ConstIterator &gIt,
                                 SxGQuery::SelSet &sels,
                                 SxGQuery::Selection &visited) const
{
   SX_CHECK (it.isValid ());
   SX_CHECK (gIt.isValid ());

   SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
   resSet= crossSelections (resSet, sels);

   bool res = true, visitable = false;

   for (ssize_t i = 0; i < it.getSizeIn (); ++i) {
      auto inIt = it.in (i);
      if (!visited->contains (inIt->getId ())) {
         visitable = true;
         res = matchAllOneIncoming (inIt, gIt, resSet, visited);
         if (!res)  break;
      }
   }
   if (visitable && res)  sels = resSet;
   return res;
}

// match all function that evaluates current
// pattern node and recursively checks all incoming
// and outgoing edges too.
bool SxGQuery::matchAllCurrent (const SxGraph<SxGQPattern>::ConstIterator &it,
                                const SxGraph<SxGProps>::ConstIterator &gIt,
                                SxGQuery::SelSet &sels,
                                SxGQuery::Selection &visited) const
{
   SX_CHECK (it.isValid ());
   SX_CHECK (gIt.isValid ());

   // if this node has already been visited
   // still evaluate it so that calling function
   // can know success/failure. But don't append
   // it to the result.
   if (visited->contains (it->getId ())) {
      SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();
      return it->matchAll(gIt, tmpSels);
   }

   bool res = false;
   ssize_t setSize = sels->getSize ();
   ssize_t siz = (setSize == 0)? 0 : sels->first()->getSize ();

   res = it->matchAll (gIt, sels);

   if (res) {
      visited->append (it->getId ());
      res = matchAllIncoming (it, gIt, sels, visited);
   }

   if (res && it.getSizeOut () > 0) {
      switch(it->getRelType ()) {
         case SxGQPattern::RelationType::OrderedDirect:
            res = matchAllChildrenOD (it, gIt, sels, visited);
            break;
         case SxGQPattern::RelationType::OrderedIndirect:
            res = matchAllChildrenOI (it, gIt, sels, visited);
            break;
         case SxGQPattern::RelationType::UnorderedDirect:
            res = matchAllChildrenUD (it, gIt, sels, visited);
            break;
         case SxGQPattern::RelationType::UnorderedIndirect:
            res = matchAllChildrenUI (it, gIt, sels, visited);
            break;
      }
   }

   if (!res)  resizeSelections (sels, setSize, siz);

   return res;
}

void SxGQuery::resizeSelections (SxGQuery::SelSet &sels,
                                 ssize_t setSize,
                                 ssize_t selSize) const
{
   sels->resize (setSize);

   for (auto it = sels->begin ();it != sels->end (); ++it) {
      (*it)->resize (selSize);
   }
}

SxGQuery::SelSet SxGQuery::crossSelections (SxGQuery::SelSet &sels1,
                                            SxGQuery::SelSet &sels2) const
{
   SX_CHECK (sels1.getPtr ());
   SX_CHECK (sels2.getPtr ());
   SxGQuery::SelSet res = SxGQuery::SelSet::create ();

   if (sels1->getSize () == 0) {
      for (auto it2 = sels2->begin ();it2 != sels2->end (); ++it2) {
         SxGQuery::Selection ptr = SxGQuery::Selection::create ();
         ptr->append ((*(*it2)));
         res->append (ptr);
      }
   } else if (sels2->getSize () == 0) {
      for (auto it1 = sels1->begin ();it1 != sels1->end (); ++it1) {
         SxGQuery::Selection ptr = SxGQuery::Selection::create ();
         ptr->append ((*(*it1)));
         res->append (ptr);
      }
   } else {
      for (auto it1 = sels1->begin ();it1 != sels1->end (); ++it1) {
         for (auto it2 = sels2->begin ();it2 != sels2->end (); ++it2) {
            SxGQuery::Selection ptr = SxGQuery::Selection::create ();
            ptr->append ((*(*it1)));
            ptr->append ((*(*it2)));
            res->append (ptr);
         }
      }
   }

   SX_CHECK (res.getPtr ());
   return res;
}

// match all ordered-direct starting from given child node in data graph
bool SxGQuery::matchAllOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                           const SxGraph<SxGProps>::ConstIterator &gIt,
                           SxGQuery::SelSet &sels, ssize_t childIdx,
                           SxGQuery::Selection &visited) const
{
   SX_CHECK (it.getSizeOut () > 0, it.getSizeOut ());

   bool res = true;
   ssize_t setSize = sels->getSize ();
   ssize_t siz = (setSize == 0)? 0 : sels->first()->getSize ();
   ssize_t outSize = gIt.getSizeOut ();

   auto currIt = gIt.out (childIdx);

   if ((outSize-childIdx) < it.getSizeOut ()) {
      return false;
   }

   for (ssize_t j = 0; j < it.getSizeOut (); ++j) {

      res = matchAllCurrent (it.out(j), gIt.out (childIdx), sels, visited);
      if (!res)  break;
      SX_CHECK (childIdx < outSize, childIdx, outSize);
      ++childIdx;
   }

   if (!res)  resizeSelections (sels, setSize, siz);
   return res;
}

// match all child nodes of pattern graph ordered-direct (default)
bool SxGQuery::matchAllChildrenOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                   const SxGraph<SxGProps>::ConstIterator &gIt,
                                   SxGQuery::SelSet &sels,
                                   SxGQuery::Selection &visited) const
{
   SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
   bool finalRes = false;
   ssize_t outSize = gIt.getSizeOut ();
   SxGQuery::Selection tmpVisited = SxGQuery::Selection::create (*visited);
   for (ssize_t c = 0; c < outSize; ++c) {
      SxGQuery::SelSet childSels = SxGQuery::SelSet::create ();
      tmpVisited = SxGQuery::Selection::create (*visited);

      if (matchAllOD (it, gIt, childSels, c, tmpVisited)) {
         finalRes = true;
         resSet->append (std::move(*crossSelections(sels, childSels)));
      }
   }
   // update the visited list once all matches
   // for a particular node have been found.
   visited = tmpVisited;
   if (finalRes)  sels = resSet;
   return finalRes;
}

// match all child nodes of pattern graph unordered-indirect
bool SxGQuery::matchAllChildrenUI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                   const SxGraph<SxGProps>::ConstIterator &pIt,
                                   SxGQuery::SelSet &sels,
                                   SxGQuery::Selection &visited) const
{
   SxGQuery::Selection prevMatches = SxGQuery::Selection::create ();
   SxPtr<SxList<SxGQuery::SelSet> > prevSels = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   ssize_t outSize = pIt.getSizeOut ();
   SxGQuery::Selection tmpVisited = SxGQuery::Selection::create (*visited);
   bool finalRes = false;
   bool res = false;

   SxGQuery::Selection currMatches   = SxGQuery::Selection::create ();
   SxPtr<SxList<SxGQuery::SelSet> > currSels = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   for (ssize_t exprIdx = 0; exprIdx < it.getSizeOut (); ++exprIdx) {

      auto chIt = it.out (exprIdx);
      finalRes = false;
      res      = false;

      for (ssize_t cIdx = 0; cIdx < outSize; ++cIdx) {

         SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
         tmpVisited = SxGQuery::Selection::create (*visited);
         res = matchAllCurrent (chIt, pIt.out (cIdx), resSet, tmpVisited);
         if (res == true) {
            // atleast one match found for curr expr
            finalRes = true;
            currMatches->append (cIdx);
            SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();

            if (prevMatches->getSize () > 0) {

               auto pMIt = prevMatches->begin ();
               auto pSIt = prevSels->begin ();
               for (;pMIt != prevMatches->end (); (++pMIt,++pSIt)) {
                  tmpSels->append (std::move( *(crossSelections (*pSIt, resSet)) ) );
               }

            } else {
               // prevMatches is empty i-e. we are in first expression
               *tmpSels = *resSet;
            }

            SX_CHECK (tmpSels->getSize () > 0, tmpSels->getSize ());
            currSels->append (tmpSels);

         }

      }

      visited = tmpVisited;

      // if no match found then return
      if (finalRes == false) {
         return false;
      }

      *prevMatches = std::move(*currMatches);
      *prevSels = std::move(*currSels);

   }

   SxGQuery::SelSet rSet = SxGQuery::SelSet::create ();
   for (auto sIt = prevSels->begin (); sIt != prevSels->end (); ++sIt) {
      rSet->append ( std::move( *(crossSelections (sels, *sIt)) ) );
   }
   sels = rSet;
   return true;
}

// match all child nodes of pattern graph ordered-indirect
bool SxGQuery::matchAllChildrenOI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                   const SxGraph<SxGProps>::ConstIterator &gIt,
                                   SxGQuery::SelSet &sels,
                                   SxGQuery::Selection &visited) const
{
   SxGQuery::Selection prevMatches = SxGQuery::Selection::create ();
   SxPtr<SxList<SxGQuery::SelSet> > prevSels = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   ssize_t cIdx = 0;
   ssize_t outSize = gIt.getSizeOut ();
   ssize_t exprsLeft = it.getSizeOut ();
   SxGQuery::Selection tmpVisited;
   bool finalRes = false;
   bool res = false;

   SxGQuery::Selection currMatches   = SxGQuery::Selection::create ();
   SxPtr<SxList<SxGQuery::SelSet> > currSels = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   for (ssize_t exprIdx = 0; exprIdx < it.getSizeOut (); ++exprIdx) {
      auto chIt = it.out (exprIdx);

      finalRes = false;
      res      = false;

      for (; cIdx < outSize; ++cIdx) {
         // if number of siblings left to process
         // are less than the number of expresssions
         // left then stop
         if ((outSize-cIdx) < exprsLeft)
            break;

         SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
         tmpVisited = SxGQuery::Selection::create (*visited);
         res = matchAllCurrent (chIt, gIt.out (cIdx), resSet, tmpVisited);
         if (res == true) {
            // atleast one match found for curr expr
            finalRes = true;
            currMatches->append (cIdx);
            SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();

            if (prevMatches->getSize () > 0) {

               auto pMIt = prevMatches->begin ();
               auto pSIt = prevSels->begin ();
               for (;pMIt != prevMatches->end (); (++pMIt,++pSIt)) {

                  if (*pMIt < cIdx) {
                     tmpSels->append (std::move( *(crossSelections (*pSIt, resSet)) ) );
                  } else {
                     // childIdx's in prevMatches are in sorted order
                     // so once (*pMIt >= cIdx), the rest of them will
                     // also be greater than cIdx. Due to "Ordered"
                     // requirement those won't match. so break loop.
                     break;
                  }

               }
            } else {
               // prevMatches is empty i-e. we are in first expression
               *tmpSels = *resSet;
            }

            SX_CHECK (tmpSels->getSize () > 0, tmpSels->getSize ());
            currSels->append (tmpSels);
         }

      }
      visited = tmpVisited;
      // if no match found then return
      if (finalRes == false)  return false;

      --exprsLeft;
      *prevMatches = std::move (*currMatches);
      *prevSels = std::move (*currSels);

      // start next loop from firstChildIdx+1
      // that matched in prev iteration.
      SX_CHECK (prevMatches->getSize() > 0, prevMatches->getSize ());
      cIdx = prevMatches->first ()+1;

   }

   SxGQuery::SelSet rSet = SxGQuery::SelSet::create ();
   for (auto sIt = prevSels->begin (); sIt != prevSels->end (); ++sIt) {
      rSet->append ( std::move( *(crossSelections (sels, *sIt)) ) );
   }
   sels = rSet;
   return true;

}

// match all child nodes of pattern graph unordered-direct
bool SxGQuery::matchAllChildrenUD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                   const SxGraph<SxGProps>::ConstIterator &gIt,
                                   SxGQuery::SelSet &sels,
                                   SxGQuery::Selection &visited) const
{

   SX_CHECK (it.getSizeOut () > 0, it.getSizeOut ());

   SxGQuery::Selection prevMatches = SxGQuery::Selection::create ();
   SxGQuery::Selection tmpVisited = SxGQuery::Selection::create (*visited);
   SxPtr<SxList<SxGQuery::SelSet> > prevSels =
                                    SxPtr<SxList<SxGQuery::SelSet> >::create ();
   ssize_t outSize = gIt.getSizeOut ();
   ssize_t cIdx    = 0;
   ssize_t exprIdx = 0;

   bool finalRes = false;
   bool res = false;
   bool res1 = false, res2 = false;

   auto exprIt = it.out (exprIdx);

   // go over the childList to find all matches for first expr
   for (; cIdx < outSize; ++cIdx) {

      SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
      tmpVisited = SxGQuery::Selection::create (*visited);
      res = matchAllCurrent (exprIt, gIt.out (cIdx), resSet, tmpVisited);
      if (res == true) {
         // at least one match found for curr expr
         finalRes = true;
         prevMatches->append (cIdx);
         prevSels->append (resSet);
      }

   }
   visited = tmpVisited;

   if (!finalRes)  return false;
   exprIdx++;

   for (; exprIdx < it.getSizeOut (); ++exprIdx) {

      exprIt = it.out (exprIdx);

      SxGQuery::Selection currMatches = SxGQuery::Selection::create ();
      SxPtr<SxList<SxGQuery::SelSet> > currSels =
                                    SxPtr<SxList<SxGQuery::SelSet> >::create ();
      finalRes = false;
      res      = false;
      auto pMIt = prevMatches->begin ();
      auto pSIt = prevSels->begin ();

      for (; pMIt != prevMatches->end (); ++pMIt, ++pSIt) {

         cIdx = *pMIt;
         res1 = false; res2 = false;
         if ((cIdx-1) >= 0) {
            auto prevIt = gIt.out ((cIdx-1));
            SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
            tmpVisited = SxGQuery::Selection::create (*visited);
            res1 = matchAllCurrent (exprIt, prevIt, resSet, tmpVisited);
            if (res1) {
               currMatches->append ((cIdx-1));
               currSels->append ( crossSelections (*pSIt, resSet) );
            }
         }

         if ((cIdx+1) < outSize) {
            auto nextIt = gIt.out ((cIdx+1));
            SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
            tmpVisited = SxGQuery::Selection::create (*visited);
            res2 = matchAllCurrent (exprIt, nextIt, resSet, tmpVisited);
            if (res2) {
               currMatches->append ((cIdx+1));
               currSels->append ( crossSelections (*pSIt, resSet) );
            }
         }

         res = res1 || res2;

         // if any one of the prevMatched child results in
         // successful match for next sibling/expr, it is enough
         // to continue with matching next exprs.
         if (res)  finalRes = true;

      }
      visited = tmpVisited;

      // if no match found then return
      if (finalRes == false)  return false;

      prevMatches = currMatches;
      prevSels = currSels;

   }
   SxGQuery::SelSet rSet = SxGQuery::SelSet::create ();
   for (auto sIt = prevSels->begin (); sIt != prevSels->end (); ++sIt) {
      rSet->append ( std::move( *(crossSelections (sels, *sIt)) ) );
   }
   sels = rSet;
   return true;

}


// recursively find one match starting at given child node unordered-direct
bool SxGQuery::matchOnceUDR (const SxGraph<SxGQPattern>::ConstIterator &it,
                             const SxGraph<SxGProps>::ConstIterator &gIt,
                             ssize_t childIdx, ssize_t exprChildIdx,
                             SxGQuery::Selection &sel,
                             SxGQuery::Selection &visited) const
{
   bool res;
   ssize_t siz = sel->getSize ();
   res = evalCurrent (it.out (exprChildIdx), gIt.out (childIdx), sel, visited);

   if (res == true) {
      ++exprChildIdx;
      if (exprChildIdx < it.getSizeOut ()) {
         bool res1 = false, res2 = false;

         if ((childIdx+1) < gIt.getSizeOut ()) {
            res1 = matchOnceUDR (it, gIt, (childIdx+1), exprChildIdx, sel, visited);
            if (res1)  return true;
            else       sel->resize (siz);
         }

         if ((childIdx-1) >= 0) {
            res2 = matchOnceUDR (it, gIt, (childIdx-1), exprChildIdx, sel, visited);
            if (res2)  return true;
            else       sel->resize (siz);
         }

         res = res1 || res2;
         // incase could not match either of next or prev neighbor child
         if (!res)  sel->resize (siz);

      } else {
         // path completed
         return true;
      }
   }

   return res;
}

// go through all child nodes and find one match unordered-direct
bool SxGQuery::matchOnceChildrenUD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                    const SxGraph<SxGProps>::ConstIterator &gIt,
                                    SxGQuery::Selection &sel,
                                    SxGQuery::Selection &visited) const
{
   for (ssize_t i=0; i < gIt.getSizeOut (); ++i) {
      if (matchOnceUDR (it, gIt, i, 0, sel, visited))
         return true;
   }
   return false;
}

// go through all child nodes and find one match unordered-indirect
bool SxGQuery::matchOnceChildrenUI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                    const SxGraph<SxGProps>::ConstIterator &gIt,
                                    SxGQuery::Selection &sel,
                                    SxGQuery::Selection &visited) const
{
   ssize_t siz = sel->getSize ();
   ssize_t outSize = gIt.getSizeOut ();
   ssize_t childIdx = 0;
   ssize_t exprChildIdx = 0;
   bool finalRes = false;

   auto chIt = it.out (exprChildIdx);
   for (; childIdx < outSize; ++childIdx) {
      if (evalCurrent (chIt, gIt.out (childIdx), sel, visited)) {
         // go to next expr
         exprChildIdx++;
         if (exprChildIdx >= it.getSizeOut ()) {
            finalRes = true;
            break;
         }

         chIt = it.out (exprChildIdx);

         // loop will increment it to 0
         childIdx = -1;
      }
   }
   if (!finalRes)  sel->resize (siz);

   return finalRes;
}

// find one match for ordered and direct
bool SxGQuery::matchOnceOD (const SxGraph<SxGQPattern>::ConstIterator &pIt,
                            const SxGraph<SxGProps>::ConstIterator &gIt,
                            SxGQuery::Selection &sel,
                            ssize_t childIdx,
                            SxGQuery::Selection &visited) const
{
   ssize_t siz = sel->getSize ();
   if ((gIt.getSizeOut ()-childIdx) < pIt.getSizeOut ()) {
      return false;
   }

   for (ssize_t j=0; j < pIt.getSizeOut (); ++j) {
      auto chIt = pIt.out (j);
      if (!evalCurrent (chIt, gIt.out (childIdx), sel, visited)) {
         sel->resize (siz);
         return false;
      }

      ++childIdx;
   }

   return true;
}

// match once ordered and direct
bool SxGQuery::matchOnceChildrenOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                    const SxGraph<SxGProps>::ConstIterator &gIt,
                                    SxGQuery::Selection &sel,
                                    SxGQuery::Selection &visited) const
{
   for (ssize_t i=0; i < gIt.getSizeOut (); ++i) {
      if (matchOnceOD (it, gIt, sel, i, visited)) return true;
   }

   return false;
}

// find one match ordered-indirect
bool SxGQuery::matchOnceChildrenOI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                    const SxGraph<SxGProps>::ConstIterator &gIt,
                                    SxGQuery::Selection &sel,
                                    SxGQuery::Selection &visited) const
{
   ssize_t siz = sel->getSize ();
   ssize_t exprsLeft = it.getSizeOut ();
   ssize_t outSize = gIt.getSizeOut ();

   bool finalRes = false;
   ssize_t i = 0;
   auto chIt = it.out (i);
   for (ssize_t j=0; j < gIt.getSizeOut (); ++j) {

      if ((outSize-j) < exprsLeft)  break;

      if (evalCurrent (chIt, gIt.out (j), sel, visited)) {
         ++i;
         if (i >= it.getSizeOut ()) {
            finalRes = true;
            break;
         }

         chIt = it.out (i);
         exprsLeft--;
      }
   }

   if (!finalRes)  sel->resize (siz);

   return finalRes;
}

