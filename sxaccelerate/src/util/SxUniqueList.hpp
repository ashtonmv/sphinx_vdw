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

template<class T, class H>
SxUniqueList<T,H>::SxUniqueList ()
   : SxList<T> ()
{
   // empty
}

template<class T, class H>
SxUniqueList<T,H>::SxUniqueList (const SxList<T> &list)
   : SxList<T> ()
{
   append (list);
}

template<class T, class H>
SxUniqueList<T,H>::SxUniqueList (const SxUniqueList<T> &list)
   : SxList<T> ()
{
   append (list);
}

template<class T, class H>
SxUniqueList<T,H>::SxUniqueList (SxUniqueList<T> &&list)
   : SxList<T> (std::move (list))
{
   ssize_t idx;
   typename SxList<T>::Node *node = this->firstElement;
   while (node) {
      idx = -1;
      if (!table.contains (node->elem, &idx)) {
         table.elements[idx] = node;
      }
      node = node->next;
   }
}

template<class T, class H>
SxUniqueList<T,H>::SxUniqueList (const std::initializer_list<T> &list)
   : SxList<T> ()
{
   append (list);
}

template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::operator= (const SxList<T> &list)
{
   removeAll ();
   return append (list);
}

template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::operator= (const SxUniqueList<T> &list)
{
   if (this == &list)  return *this;
   removeAll ();
   return append (list);
}

template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::operator= (SxUniqueList<T> &&list)
{
   if (this == &list)  return *this;
   table.removeAll ();
   SxList<T>::operator= (std::move(list));
   ssize_t idx;
   typename SxList<T>::Node *node = this->firstElement;
   while (node) {
      idx = -1;
      if (!table.contains (node->elem, &idx)) {
         table.elements[idx] = node;
      }
      node = node->next;
   }
   return *this;
}

template<class T, class H>
SxUniqueList<T,H> &
SxUniqueList<T,H>::operator= (const std::initializer_list<T> &list)
{
   removeAll ();
   return append (list);
}


template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::operator<< (const T &elem)
{
   return append (elem);
}

template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::operator<< (const SxList<T> &list)
{
   return append (list);
}

template<class T, class H>
SxUniqueList<T,H> &
SxUniqueList<T,H>::operator<< (const std::initializer_list<T> &list)
{
   return append (list);
}

//template<class T, class H>
//ssize_t SxUniqueList<T,H>::operator|= (const T &elem)
//{
//   ssize_t idx;
//   if (!table.contains (elem, &idx))  {
//      SxList<T>::append (elem);
//      table.elements[idx] = this->lastElement;
//      idx = this->getSize() - 1;
//   }  else  {
//      idx = 0;
//      typename SxList<T>::Node *ptr = this->firstElement;
//      while (ptr)  {
//         if (ptr->element == elem)  {
//            return idx;
//         }
//         ptr = ptr->next;
//         idx++;
//      }
//   }
//   return idx;
//}

template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::append (const T &elem)
{
   ssize_t idx = -1;
   if (!table.contains (elem, &idx))  {
      SxList<T>::append (elem);
      table.elements[idx] = this->lastElement;
   }
   return *this;
}

template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::append (const SxList<T> &list)
{
   ssize_t idx = -1;
   typename SxList<T>::ConstIterator it;
   for (it = list.begin (); it != list.end (); ++it) {
      if (!table.contains (*it, &idx))  {
         SxList<T>::append (*it);
         table.elements[idx] = this->lastElement;
      }
   }
   return *this;
}

template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::insert (ssize_t newPos, const T &elem)
{
   ssize_t idx;
   if (!table.contains (elem, &idx))  {
      SxList<T>::insert (newPos, elem);
      typename SxList<T>::Node *ptr = this->firstElement;
      for (ssize_t i = 0; i < newPos; i++)  {
         ptr = ptr->next;
      }
      table.elements[idx] = ptr;
   }
   return *this;
}

template<class T, class H>
SxUniqueList<T,H> &
SxUniqueList<T,H>::append (const std::initializer_list<T> &list)
{
   ssize_t idx = -1;
   typename SxList<T>::ConstIterator it;
   for (const T &elem : list) {
      if (!table.contains (elem, &idx))  {
         SxList<T>::append (elem);
         table.elements[idx] = this->lastElement;
      }
   }
   return *this;
}


template<class T, class H>
SxUniqueList<T,H> &SxUniqueList<T,H>::prepend (const T &elem)
{
   ssize_t idx;
   if (!table.contains (elem, &idx))  {
      SxList<T>::prepend (elem);
      table.elements[idx] = this->firstElement;
   }
   return *this;
}

template<class T, class H>
ssize_t SxUniqueList<T,H>::findPos (const T &elem) const
{
   ssize_t idx = table.findPos (elem);
   if (idx >= 0)  {
      typename SxList<T>::Node *ptr = table.elements[idx];
      idx = 0;
      while ((ptr = ptr->prev))  {
         idx++;
      }
   }

   return idx;
}

template<class T, class H>
bool SxUniqueList<T,H>::contains (const T &elem) const
{
   return table.contains (elem);
}

template<class T, class H>
void SxUniqueList<T,H>::remove (ssize_t idx)
{
   SX_CHECK (idx >= 0 && idx < this->size, idx, this->size);

   typename SxList<T>::Node *ptr = this->firstElement;
   for (ssize_t i = 0; i < idx; i++)  {
      ptr = ptr->next;
   }

   remove (ptr->element);
}

template<class T, class H>
void SxUniqueList<T,H>::removeElement (const T &elem)
{
   ssize_t idx = table.findPos (elem);
   if (idx >= 0)  {
      SxList<T>::removeItem (table.elements[idx]);
      table.remove (idx);
   }
}

template<class T, class H>
void SxUniqueList<T,H>::removeFirst ()
{
   if (this->firstElement)  {
      removeElement (this->firstElement->elem);
   }
}

template<class T, class H>
void SxUniqueList<T,H>::removeLast ()
{
   if (this->lastElement)  {
      removeElement (this->lastElement->elem);
   }
}

template<class T, class H>
void SxUniqueList<T,H>::removeAll ()
{
   SxList<T>::removeAll ();
   table.removeAll ();
}

template<class T, class H>
void SxUniqueList<T,H>::print () const
{
   table.print ();

   std::cout << "size=" << this->size << std::endl;

   typename SxList<T>::Node *ptr = this->firstElement;
   ssize_t idx = 0;
   while (ptr)  {
      std::cout << idx << ":" << (void*)ptr << "='" << ptr->element << "' ";
      ptr = ptr->next;
      idx++;
   }
   if (this->size > 0)  {
      std::cout << std::endl;
   }
}


/** \brief List of elements, where each element is unique

   \b SxUniqueList = SFHIngX List containing unique elements only

    \author C. Freysoldt, freyso@fhi-berlin.mpg.de
    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template <class T>
class SxUniqueList<T,SxNull> : public SxList<T> 
{
   public:
      // -- Constructor 
      SxUniqueList ();
      SxUniqueList (const SxList<T> &);
      SxUniqueList (const SxUniqueList<T> &);
      SxUniqueList (SxUniqueList<T> &&);

      // -- Assignment
      SxUniqueList<T,SxNull> &operator= (const SxList<T> &);
      SxUniqueList<T,SxNull> &operator= (const SxUniqueList<T> &);
      SxUniqueList<T,SxNull> &operator= (SxUniqueList<T> &&);

      // -- Insertion
      inline SxUniqueList<T,SxNull> &operator<< (const T &);
      inline SxUniqueList<T,SxNull> &operator<< (const SxList<T> &);

      SxUniqueList<T,SxNull> &append (const T &);
      SxUniqueList<T,SxNull> &append (const SxList<T> &);
      SxUniqueList<T,SxNull> &insert (ssize_t newPos, const T &);
      SxUniqueList<T,SxNull> &prepend (const T &);

      void print () const;
};

template<class T>
SxUniqueList<T,SxNull>::SxUniqueList ()
   : SxList<T> ()
{
   // empty
}

template<class T>
SxUniqueList<T,SxNull>::SxUniqueList (const SxList<T> &list)
   : SxList<T> ()
{
   append (list);
}

template<class T>
SxUniqueList<T,SxNull>::SxUniqueList (const SxUniqueList<T> &list)
   : SxList<T> (list)
{ }

template<class T>
SxUniqueList<T,SxNull>::SxUniqueList (SxUniqueList<T> &&list)
   : SxList<T> (std::move(list))
{ }

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::operator= (const SxList<T> &l)
{
   SxList<T>::removeAll ();
   return append (l);
}

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::operator= (
   const SxUniqueList<T> &list)
{
   SxList<T>::removeAll ();
   return append (list);
}

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::operator= (
   SxUniqueList<T> &&list)
{
   SxList<T>::operator=(std::move(list));
   return *this;
}

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::operator<< (const T &elem)
{
   return append (elem);
}

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::operator<< (const SxList<T> &l)
{
   return append (l);
}

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::append (const T &elem)
{
   if (!SxList<T>::contains (elem))  {
      SxList<T>::append (elem);
   }
   return *this;
}

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::append (const SxList<T> &list)
{
   typename SxList<T>::ConstIterator it = list.begin ();
   typename SxList<T>::ConstIterator itEnd = list.end ();
   for (; it != itEnd; ++it) {
      if (!SxList<T>::contains (*it))  {
         SxList<T>::append (*it);
      }
   }
   return *this;
}

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::insert (ssize_t newPos,
                                                        const T &elem)
{
   if (!SxList<T>::contains (elem))  {
      SxList<T>::insert (newPos, elem);
   }
   return *this;
}

template<class T>
SxUniqueList<T,SxNull> &SxUniqueList<T,SxNull>::prepend (const T &elem)
{
   if (!SxList<T>::contains (elem))  {
      SxList<T>::prepend (elem);
   }
   return *this;
}

template<class T>
void SxUniqueList<T,SxNull>::print () const
{
   typename SxList<T>::Node *ptr = this->firstElement;
   ssize_t idx = 0;
   while (ptr)  {
      std::cout << idx << ":" << (void*)ptr << "='" << ptr->element << "' ";
      ptr = ptr->next;
      idx++;
   }
   if (this->size > 0)  {
      std::cout << std::endl;
   }
}
