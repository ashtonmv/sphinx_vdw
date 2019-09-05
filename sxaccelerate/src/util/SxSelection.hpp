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

template<class C>
SxSelection<C>::SxSelection ()
{ }

template<class C>
SxSelection<C>::SxSelection (const SxCList<typename SxSelection<C>::SelIdx> &in)
{
   SX_CHECK (in.getSize () >= 0, in.getSize ());
   selected = SxPtr<typename SxSelection<C>::SelStorage>::create (in);
}

template<class C>
SxSelection<C>::SxSelection (const SxPtr<typename
                             SxSelection<C>::SelStorage> &selected_)
   : selected(selected_)
{ }

template<class C>
SxSelection<C>::~SxSelection ()
{ }

template<class C>
void SxSelection<C>::setContainer (SxPtr<C> ptr)
{
   SX_CHECK (ptr.getPtr ());
   container = ptr;
}

template<class C>
ssize_t SxSelection<C>::getSize () const
{
   SX_CHECK (selected.getPtr ());
   return selected->getSize ();
}

template<class C>
SxPtr<C> SxSelection<C>::getContainer () const
{
   return container;
}


// sort using bubble sort [it, lastIt]
template<class C>
template<class Iterator1,class Iterator2>
void SxSelection<C>::selSort (const Iterator1 &it1_, const Iterator1 &lastIt1_,
                              const Iterator2 &it2_, const Iterator2 &lastIt2_,
                              SxCBoundPtr<int, const decltype(Iterator1())&,
                              const decltype(Iterator1())&> comp)
{
   SX_CHECK(it1_.isValid ());
   SX_CHECK(lastIt1_.isValid ());// should point to a valid element

   auto it1 = it1_; auto lastIt1 = lastIt1_;
   auto it2 = it2_; auto lastIt2 = lastIt2_;

   if (it1 == lastIt1) return;

   Iterator1 tmpIt = it1, tmpIt2 = it1, itPrev = it1;
   Iterator2 tmpIt21 = it2, tmpIt22 = it2, itPrev2 = it2;

   bool forward = it1.isForward ();

   if (forward)  { ++lastIt1; ++lastIt2; }
   else          { --lastIt1; --lastIt2; }

   while (tmpIt2 != lastIt1) {
      SX_CHECK(tmpIt2.isValid ());
      tmpIt  = it1;  tmpIt21 = it2;
      itPrev = it1;  itPrev2 = it2;

      if (forward)  { ++tmpIt; ++tmpIt21; }
      else          { --tmpIt; --tmpIt21; }

      SX_CHECK(tmpIt.isValid ());

      while (tmpIt != lastIt1) {
         SX_CHECK(tmpIt.isValid ());
         SX_CHECK(itPrev.isValid ());

         if (comp(tmpIt,itPrev) == -1) {
            auto val = *tmpIt;
            *tmpIt   = *itPrev;
            *itPrev  = val;

            auto val1 = *tmpIt21;
            *tmpIt21  = *itPrev2;
            *itPrev2  = val1;
         }
         if (forward) {
            ++tmpIt;   ++itPrev;
            ++tmpIt21; ++itPrev2;
         } else {
            --tmpIt;   --itPrev;
            --tmpIt21; --itPrev2;
         }
      }

      if (forward)  { ++tmpIt2; ++tmpIt22; }
      else          { --tmpIt2; --tmpIt22; }
   }
}

template<class C>
template<class Iterator1,class Iterator2>
void SxSelection<C>::selSort (const Iterator1 &it1, const Iterator1 &lastIt1,
                              const Iterator2 &it2, const Iterator2 &lastIt2)
{
   selSort (it1, lastIt1, it2, lastIt2,
            [](const Iterator1 &itA, const Iterator1 &itB)->int
            { if (*itA < *itB)        return -1;
              else if (*itA == *itB)  return  0;
              else                    return  1;
            });
}

template<class C>
void SxSelection<C>::removeDuplicates(SxCList<typename SxSelection<C>::TElem> *lst,
                                      SxCList<typename SxSelection<C>::SelIdx> *lstIdx)
{
   if (lst->getSize () == 0) return;

   typedef typename SxCList<TElem>::Node   ListNode;
   typedef typename SxCList<SelIdx>::Node  IdxNode;

   ListNode *ptr   = lst->firstElement;
   IdxNode *idxPtr = lstIdx->firstElement;

   while (ptr->next != NULL) {
      if (ptr->elem == ptr->next->elem) {
         lst->removeItem (ptr->next);
         lstIdx->removeItem (idxPtr->next);
         continue;
      }
      ptr    = ptr->next;
      idxPtr = idxPtr->next;
   }

}


// Compute Union
template<class C>
SxSelection<C> SxSelection<C>::operator| (const SxSelection<C> &in)
{

   // they must be representing the same object
   SX_CHECK (container == in.container);

   // generate a separate list of selected items
   // and their corresponding indices
   SxCList<TElem>  elemList;
   SxCList<SelIdx> listIdx;

   for (auto it = begin (); it != end (); ++it) {
      elemList.append (*it);
      listIdx.append (it.getSelIdx ());
   }

   SxCList<TElem>  inElemList;
   SxCList<SelIdx> inListIdx;

   for (auto it1 = in.begin (); it1 != in.end (); ++it1) {
      inElemList.append (*it1);
      inListIdx.append (it1.getSelIdx ());
   }

   // sort both selections along with their idx list
   if (elemList.getSize () > 0)
      selSort (elemList.begin (), elemList.fromLast (),
               listIdx.begin (), listIdx.fromLast ());

   if (inElemList.getSize () > 0)
      selSort (inElemList.begin (), inElemList.fromLast (),
               inListIdx.begin (), inListIdx.fromLast ());

   removeDuplicates (&elemList, &listIdx);
   removeDuplicates (&inElemList, &inListIdx);

   auto inIt = inElemList.begin (); auto inIdxIt = inListIdx.begin ();
   auto it   = elemList.begin ();   auto idxIt   = listIdx.begin ();


   SxCList<SelIdx> resList;
   while (inIt.isValid () && it.isValid ()) {
      if (*it < *inIt) {
         resList.append (*idxIt);
         ++it; ++idxIt;
      } else if (*inIt < *it) {
         resList.append (*inIdxIt);
         ++inIt; ++inIdxIt;
      } else {
         resList.append (*inIdxIt);
         ++inIt;    ++it;
         ++inIdxIt; ++idxIt;
      }
   }

   while (inIt.isValid ()) {
      resList.append (*inIdxIt);
      ++inIt; ++inIdxIt;
   }

   while (it.isValid ()) {
      resList.append (*idxIt);
      ++it; ++idxIt;
   }

   SxSelection<typename Container::SelContainer> sel(resList);
   sel.container = this->getContainer ();
   return sel;
}

// Compute Intersection
template<class C>
SxSelection<C> SxSelection<C>::operator& (const SxSelection<C> &in)
{
   // they must be representing the same object
   SX_CHECK (container == in.container);

   // generate a separate list of selected items
   // and their corresponding indices
   SxCList<TElem>  elemList;
   SxCList<SelIdx> listIdx;

   for (auto it = begin (); it != end (); ++it) {
      elemList.append (*it);
      listIdx.append (it.getSelIdx ());
   }

   SxCList<TElem>  inElemList;
   SxCList<SelIdx> inListIdx;

   for (auto it1 = in.begin (); it1 != in.end (); ++it1) {
      inElemList.append (*it1);
      inListIdx.append (it1.getSelIdx ());
   }

   // sort both selections along with their idx list
   if (elemList.getSize () > 0)
      selSort (elemList.begin (), elemList.fromLast (),
               listIdx.begin (), listIdx.fromLast ());

   if (inElemList.getSize () > 0)
      selSort (inElemList.begin (), inElemList.fromLast (),
               inListIdx.begin (), inListIdx.fromLast ());


   removeDuplicates (&elemList, &listIdx);
   removeDuplicates (&inElemList, &inListIdx);

   auto inIt = inElemList.begin (); auto inIdxIt = inListIdx.begin ();
   auto it   = elemList.begin ();   auto idxIt   = listIdx.begin ();


   SxCList<SelIdx> resList;
   while (inIt.isValid () && it.isValid ()) {
      if (*it < *inIt) {
         ++it; ++idxIt;
      } else if (*inIt < *it) {
         ++inIt; ++inIdxIt;
      } else {
         resList.append (*inIdxIt);
         ++inIt;    ++it;
         ++inIdxIt; ++idxIt;
      }
   }

   SxSelection<typename Container::SelContainer> sel(resList);
   sel.container = this->getContainer ();
   return sel;
}

// Compute Complement
template<class C>
SxSelection<C> SxSelection<C>::operator- (const SxSelection<C> &in)
{
   // they must be representing the same object
   SX_CHECK (container == in.container);

   // generate a separate list of selected items
   // and their corresponding indices
   SxCList<TElem>  elemList;
   SxCList<SelIdx> listIdx;

   for (auto it = begin (); it != end (); ++it) {
      elemList.append (*it);
      listIdx.append (it.getSelIdx ());
   }

   SxCList<TElem>  inElemList;
   SxCList<SelIdx> inListIdx;

   for (auto it1 = in.begin (); it1 != in.end (); ++it1) {
      inElemList.append (*it1);
      inListIdx.append (it1.getSelIdx ());
   }

   // sort both selections along with their idx list
   if (elemList.getSize () > 0)
      selSort (elemList.begin (), elemList.fromLast (),
               listIdx.begin (), listIdx.fromLast ());

   if (inElemList.getSize () > 0)
      selSort (inElemList.begin (), inElemList.fromLast (),
               inListIdx.begin (), inListIdx.fromLast ());

   removeDuplicates (&elemList, &listIdx);
   removeDuplicates (&inElemList, &inListIdx);

   auto inIt = inElemList.begin (); auto inIdxIt = inListIdx.begin ();
   auto it   = elemList.begin ();   auto idxIt   = listIdx.begin ();


   SxCList<SelIdx> resList;
   while (inIt.isValid () && it.isValid ()) {
      if (*it < *inIt) {
         resList.append (*idxIt);
         ++it; ++idxIt;
      } else if (*inIt < *it) {
         ++inIt; ++inIdxIt;
      } else {
         ++inIt;    ++it;
         ++inIdxIt; ++idxIt;
      }
   }

   while (it.isValid ()) {
      resList.append (*idxIt);
      ++it; ++idxIt;
   }

   SxSelection<typename Container::SelContainer> sel(resList);
   sel.container = this->getContainer ();
   return sel;
}

// Compute Symmetric Complement
template<class C>
SxSelection<C> SxSelection<C>::operator!= (const SxSelection<C> &in)
{
   // they must be representing the same object
   SX_CHECK (container == in.container);

   // generate a separate list of selected items
   // and their corresponding indices
   SxCList<TElem>  elemList;
   SxCList<SelIdx> listIdx;

   for (auto it = begin (); it != end (); ++it) {
      elemList.append (*it);
      listIdx.append (it.getSelIdx ());
   }

   SxCList<TElem>  inElemList;
   SxCList<SelIdx> inListIdx;

   for (auto it1 = in.begin (); it1 != in.end (); ++it1) {
      inElemList.append (*it1);
      inListIdx.append (it1.getSelIdx ());
   }

   // sort both selections along with their idx list
   if (elemList.getSize () > 0)
      selSort (elemList.begin (), elemList.fromLast (),
               listIdx.begin (), listIdx.fromLast ());

   if (inElemList.getSize () > 0)
      selSort (inElemList.begin (), inElemList.fromLast (),
               inListIdx.begin (), inListIdx.fromLast ());


   removeDuplicates (&elemList, &listIdx);
   removeDuplicates (&inElemList, &inListIdx);

   auto inIt = inElemList.begin (); auto inIdxIt = inListIdx.begin ();
   auto it   = elemList.begin ();   auto idxIt   = listIdx.begin ();


   SxCList<SelIdx> resList;
   while (inIt.isValid () && it.isValid ()) {
      if (*it < *inIt) {
         resList.append (*idxIt);
         ++it; ++idxIt;
      } else if (*inIt < *it) {
         resList.append (*inIdxIt);
         ++inIt; ++inIdxIt;
      } else {
         ++inIt;    ++it;
         ++inIdxIt; ++idxIt;
      }
   }

   while (inIt.isValid ()) {
      resList.append (*inIdxIt);
      ++inIt; ++inIdxIt;
   }

   while (it.isValid ()) {
      resList.append (*idxIt);
      ++it; ++idxIt;
   }

   SxSelection<typename Container::SelContainer> sel(resList);
   sel.container = this->getContainer ();
   return sel;
}

template<class C>
typename C::Iterator
SxSelection<C>::getIterator (const typename SelStorage::Iterator &it)
{
   SX_CHECK (it.isValid ());
   SX_CHECK (container.getPtr ());
   return container->getIterator (*it);
}

template<class C>
typename C::ConstIterator
SxSelection<C>::getIterator (const typename SelStorage::ConstIterator &it) const
{
   SX_CHECK (it.isValid ());
   SX_CHECK (container.getPtr ());
   return container->getIterator (*it);
}

template<class C>
typename C::ConstIterator
SxSelection<C>::getConstIterator (const Iterator &it) const
{
   SX_CHECK (it.isValid ());
   SX_CHECK (container.getPtr ());
   return container->getIterator (it.getSelIdx ());
}

// Iterators
template<class C>
typename SxSelection<C>::ConstIterator
SxSelection<C>::begin () const
{
   SX_CHECK (selected.getPtr ());
   return ConstIterator (sx::Forward, this,
                         selected->begin ());
}

template<class C>
typename SxSelection<C>::Iterator
SxSelection<C>::begin ()
{
   SX_CHECK (selected.getPtr ());
   return Iterator (sx::Forward, this,
                    selected->begin ());
}

template<class C>
typename SxSelection<C>::ConstIterator
SxSelection<C>::begin (ssize_t idx_, sx::Direction dir_) const
{
   SX_CHECK (selected.getPtr ());
   if (idx_ >=0 && idx_ < getSize())
      return ConstIterator (dir_, this,
                            selected->begin(idx_) );
   else
      return ConstIterator (dir_, this);
}

template<class C>
typename SxSelection<C>::Iterator
SxSelection<C>::begin (ssize_t idx_, sx::Direction dir_)
{
   SX_CHECK (selected.getPtr ());
   if (idx_ >= 0 && idx_ < getSize())
      return Iterator (dir_, this,
                       selected->begin(idx_) );
   else
      return Iterator (dir_, this);
}

template<class C>
typename SxSelection<C>::ConstIterator
SxSelection<C>::end () const
{
   SX_CHECK (selected.getPtr ());
   return typename SxSelection<C>::ConstIterator (sx::Forward, this);
}

template<class C>
typename SxSelection<C>::Iterator
SxSelection<C>::end ()
{
   SX_CHECK (selected.getPtr ());
   return typename SxSelection<C>::Iterator (sx::Forward, this);
}

template<class C>
typename SxSelection<C>::ConstIterator
SxSelection<C>::fromLast () const
{
   SX_CHECK (selected.getPtr ());
   return typename SxSelection<C>::ConstIterator (sx::Backward, this,
                                                  selected->fromLast ());
}

template<class C>
typename SxSelection<C>::Iterator
SxSelection<C>::fromLast ()
{
   SX_CHECK (selected.getPtr ());
   return typename SxSelection<C>::Iterator (sx::Backward, this,
                                             selected->fromLast ());
}

template<class C>
typename SxSelection<C>::ConstIterator
SxSelection<C>::toFirst () const
{
   SX_CHECK (selected.getPtr ());
   return typename SxSelection<C>::ConstIterator (sx::Backward, this);
}

template<class C>
typename SxSelection<C>::Iterator
SxSelection<C>::toFirst ()
{
   SX_CHECK (selected.getPtr ());
   return typename SxSelection<C>::Iterator (sx::Backward, this);
}

