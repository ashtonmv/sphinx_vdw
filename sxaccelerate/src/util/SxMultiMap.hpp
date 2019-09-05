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

template<class K, class V, class H>
SxMultiMap<K,V,H>::SxMultiMap () //: SxMemConsumer(), nBytes(0)
{
   // empty
}

template<class K, class V, class H>
SxMultiMap<K,V,H>::SxMultiMap (const SxMultiMap<K,V,H> &map_)
{
   append (map_, true);
}

template<class K, class V, class H>
SxMultiMap<K,V,H> &SxMultiMap<K,V,H>::operator= (const SxMultiMap<K,V,H> &map_)
{
   removeAll ();
   return append (map_, true);
}

template<class K, class V, class H>
SxMultiMap<K,V,H>::~SxMultiMap ()
{
   removeAll ();
}

template<class K, class V, class H>
V &SxMultiMap<K,V,H>::operator() (const K &key_)
{
   return insert (key_, V (), true);
}

template<class K, class V, class H>
const V &SxMultiMap<K,V,H>::operator() (const K &key_) const
{
   return getValues(key_).first ();
}

template<class K, class V, class H>
SxList<V> &SxMultiMap<K,V,H>::getValues (const K &key_)
{
   ssize_t idx = table.findPos (key_);
   if (idx >= 0)  {
      return table.value(idx)->element.values;
   }

   throw SxException("The map has no elements for this key.",__FILE__,__LINE__);
}

template<class K, class V, class H>
const SxList<V> &SxMultiMap<K,V,H>::getValues (const K &key_) const
{
   ssize_t idx = table.findPos (key_);
   if (idx >= 0)  {
      return table.value(idx)->element.values;
   }

   throw SxException("The map has no elements for this key.",__FILE__,__LINE__);
}

template<class K, class V, class H>
ssize_t SxMultiMap<K,V,H>::findKey (const K &key_) const
{
   return table.findPos (key_);
}

template<class K, class V, class H>
SxList<V> &SxMultiMap<K,V,H>::getValuesIdx (ssize_t idx_)
{
   return table.value(idx_)->element.values;
}

template<class K, class V, class H>
const SxList<V> &SxMultiMap<K,V,H>::getValuesIdx (ssize_t idx_) const
{
   return table.value(idx_)->element.values;
}

template<class K, class V, class H>
V &SxMultiMap<K,V,H>::insert (const K &key_, const V &value_,bool overwrite_)
{
   ssize_t idx;
   if (!table.contains (key_, &idx))  {
      nodes.append (Node ());
      Node &node = nodes.lastElement->element;
      node.key = key_;
      node.values.append (value_);
      keys.append (key_);
      table.elements[idx] = keys.lastElement;
      table.value(idx) = nodes.lastElement;

      sequence.append (Pair ());
      Pair &pair = sequence.lastElement->element;
      pair.node = nodes.lastElement;
      pair.value = node.values.lastElement;
      node.sequence.append (sequence.lastElement);

      return node.values.lastElement->element;
   }  else  {
      Node &node = table.value(idx)->element;
      if (overwrite_)  {
         // --- overwrite value of existing key
         if (node.values.getSize () > 1)  {
            node.values.removeAll ();
            node.values << value_;

            typename SxList<typename SxList<Pair>::List*>::ConstIterator it;
            typename SxList<typename SxList<Pair>::List*>::ConstIterator itEnd;

            it = node.sequence.begin ();
            itEnd = node.sequence.end ();

            for (++it; it != itEnd; ++it)  {
               sequence.removeItem (*it);
            }

            typename SxList<Pair>::List* ptr = node.sequence.first ();
            node.sequence.removeAll ();
            node.sequence << ptr;
            ptr->element.value = node.values.lastElement;
         }
      }  else  {
         // --- append to list of values with the same key
         node.values.append (value_);

         sequence.append (Pair ());
         Pair &pair = sequence.lastElement->element;
         pair.node = table.value(idx);
         pair.value = node.values.lastElement;

         node.sequence.append (sequence.lastElement);
      }
      return node.values.lastElement->element;
   }
}

template<class K, class V, class H>
SxMultiMap<K,V,H> &SxMultiMap<K,V,H>::append (const K &key_, const V &value_)
{
   insert (key_, value_, false);
   return *this;
}

template<class K, class V, class H>
ssize_t SxMultiMap<K,V,H>::getSize () const
{
   return sequence.getSize ();
}

template<class K, class V, class H>
bool SxMultiMap<K,V,H>::containsKey (const K &key_) const
{
   return table.findPos (key_) >= 0;
}

template<class K, class V, class H>
bool SxMultiMap<K,V,H>::containsValue (const V &value_) const
{
   typename SxList<Pair>::ConstIterator it = sequence.begin ();
   typename SxList<Pair>::ConstIterator itEnd = sequence.end ();

   for (; it != itEnd; ++it)  {
      if ((*it).value->element == value_)  {
         return true;
      }
   }
   return false;
}

template<class K, class V, class H>
bool SxMultiMap<K,V,H>::removeKey (const K &key_)
{
   ssize_t idx = table.findPos (key_);
   if (idx >= 0)  {
      typename SxList<typename SxList<Pair>::List*>::ConstIterator it;
      typename SxList<typename SxList<Pair>::List*>::ConstIterator itEnd;

      it = table.value(idx)->element.sequence.begin ();
      itEnd = table.value(idx)->element.sequence.end ();

      for (; it != itEnd; ++it)  {
         sequence.removeItem (*it);
      }

      nodes.removeItem (table.value(idx));
      keys.removeItem (table.elements[idx]);
      table.remove (idx);
      return true;
   }

   return false; // errror
}

template<class K, class V, class H>
bool SxMultiMap<K,V,H>::removeValue (const V &value_)
{
   typename SxList<Pair>::Iterator it = sequence.begin ();
   typename SxList<Pair>::Iterator itEnd = sequence.end ();

   for (; it != itEnd; ++it)  {
      if ((*it).value->element == value_)  {
         Node &node = (*it).node->element;
         if (node.values.getSize () < 2)  {
            return removeKey (node.key);
         }  else  {
            // --- remove one value
            node.values.removeItem ((*it).value);
            node.sequence.removeElement (it.ptr);
            sequence.removeItem (it.ptr);
         }
         return true;
      }
   }
   return false;
}

template<class K, class V, class H>
void SxMultiMap<K,V,H>::removeAll ()
{
   table.removeAll ();
   keys.removeAll ();
   nodes.removeAll ();
   sequence.removeAll ();
}

template<class K, class V, class H>
V &SxMultiMap<K,V,H>::first ()
{
   SX_CHECK (getSize() > 0, getSize());
   return sequence.first().value->element;
}

template<class K, class V, class H>
const V &SxMultiMap<K,V,H>::first () const
{
   SX_CHECK (getSize() > 0, getSize());
   return sequence.first().value->element;
}

template<class K, class V, class H>
V &SxMultiMap<K,V,H>::last ()
{
   SX_CHECK (getSize() > 0, getSize());
   return sequence.last().value->element;
}

template<class K, class V, class H>
const V &SxMultiMap<K,V,H>::last () const
{
   SX_CHECK (getSize() > 0, getSize());
   return sequence.last().value->element;
}

template<class K, class V, class H>
const SxList<K> &SxMultiMap<K,V,H>::getKeys () const
{
   return keys;
}

//template<class K, class V, class H>
//SxList<K> SxMultiMap<K,V,H>::getKeys () const
//{
//   SxList<K> keyList;
//
//   typename SxList<Pair>::ConstIterator it = sequence.begin ();
//   typename SxList<Pair>::ConstIterator itEnd = sequence.end ();
//
//   for (; it != itEnd; ++it)  {
//      keyList << (*it).node->element.key;
//   }
//
//   return keyList;
//}

template<class K, class V, class H>
SxList<V> SxMultiMap<K,V,H>::getValues () const
{
   SxList<V> valueList;

   typename SxList<Pair>::ConstIterator it = sequence.begin ();
   typename SxList<Pair>::ConstIterator itEnd = sequence.end ();

   for (; it != itEnd; ++it)  {
      valueList << (*it).value->element;
   }

   return valueList;
}

template<class K, class V, class H>
typename SxMultiMap<K,V,H>::ConstIterator SxMultiMap<K,V,H>::begin () const
{
   return ConstIterator (sequence.begin ());
}

template<class K, class V, class H>
typename SxMultiMap<K,V,H>::Iterator SxMultiMap<K,V,H>::begin ()
{
   return Iterator (sequence.begin ());
}

template<class K, class V, class H>
typename SxMultiMap<K,V,H>::ConstIterator SxMultiMap<K,V,H>::end () const
{
   return ConstIterator (sequence.end ());
}

template<class K, class V, class H>
typename SxMultiMap<K,V,H>::Iterator SxMultiMap<K,V,H>::end ()
{
   return Iterator (sequence.end ());
}

template<class K, class V, class H>
SxMultiMap<K,V,H> &SxMultiMap<K,V,H>::append (const SxMultiMap<K,V,H> &map,
                                              bool                    overwrite)
{
   ConstIterator it;
   for (it = map.begin(); it != map.end(); ++it)  {
      if (overwrite || !containsKey (it.getKey ()))  {
         append (it.getKey (), it.getValue ());
      }
   }
   return *this;
}

template<class K, class V, class H>
void SxMultiMap<K,V,H>::print () const
{
   typename SxList<Pair>::ConstIterator it = sequence.begin ();
   typename SxList<Pair>::ConstIterator itEnd = sequence.end ();

   for (; it != itEnd; ++it)  {
      std::cout << (*it).node->element.key << " = "
                << (*it).value->element << std::endl;
   }
}

// --------------------------------------------------------------------------

template<class K, class V, class H>
SxMultiMap<K,V,H>::ConstIterator::ConstIterator ()
{
   // empty
}


template<class K, class V, class H>
SxMultiMap<K,V,H>::ConstIterator::ConstIterator (
      const typename SxList<Pair>::ConstIterator &it_)
   : it(it_)
{
   // empty
}


template<class K, class V, class H>
SxMultiMap<K,V,H>::ConstIterator::~ConstIterator ()
{
   // empty
}



template<class K, class V, class H>
bool 
SxMultiMap<K,V,H>::ConstIterator::operator== (
      const typename SxMultiMap<K,V,H>::ConstIterator &in) const
{
   return it == in.it;
}


template<class K, class V, class H>
bool SxMultiMap<K,V,H>::ConstIterator::operator!= (
      const typename SxMultiMap<K,V,H>::ConstIterator &in) const
{
   return !operator== (in);
}


template<class K, class V, class H>
const typename SxMultiMap<K,V,H>::ConstIterator
SxMultiMap<K,V,H>::ConstIterator::operator++ (int)
{
   typename SxMultiMap<K,V,H>::ConstIterator copy (it);
   ++it;
   return copy;
}

template<class K, class V, class H>
typename SxMultiMap<K,V,H>::ConstIterator
&SxMultiMap<K,V,H>::ConstIterator::operator++ ()
{
   ++it;
   return *this;
}


template<class K, class V, class H>
const K &SxMultiMap<K,V,H>::ConstIterator::getKey () const
{
   return it->node->element.key;
}

template<class K, class V, class H>
const V &SxMultiMap<K,V,H>::ConstIterator::getValue () const
{
   return it->value->element;
}

//template<class K, class V, class H>
//const SxList<V> &SxMultiMap<K,V,H>::ConstIterator::getValues () const
//{
//   return it->node->element.values;
//}

template<class K, class V, class H>
const V &SxMultiMap<K,V,H>::ConstIterator::operator-> () const
{
   return it->value->element;
}


// --------------------------------------------------------------------------

template<class K, class V, class H>
SxMultiMap<K,V,H>::Iterator::Iterator ()
{
   // empty
}


template<class K, class V, class H>
SxMultiMap<K,V,H>::Iterator::Iterator (
   const typename SxList<Pair>::Iterator &it_)
   : it(it_)
{
   // empty
}


template<class K, class V, class H>
SxMultiMap<K,V,H>::Iterator::~Iterator ()
{
   // empty
}


template<class K, class V, class H>
bool SxMultiMap<K,V,H>::Iterator::operator== (
   const typename SxMultiMap<K,V,H>::Iterator &in)
{
   return it == in.it;
}


template<class K, class V, class H>
bool SxMultiMap<K,V,H>::Iterator::operator!= (
   const typename SxMultiMap<K,V,H>::Iterator &in)
{
   return !operator== (in);
}


template<class K, class V, class H>
const typename SxMultiMap<K,V,H>::Iterator
SxMultiMap<K,V,H>::Iterator::operator++ (int)
{
   typename SxMultiMap<K,V,H>::Iterator copy (it);
   ++it;
   return copy;
}

template<class K, class V, class H>
typename SxMultiMap<K,V,H>::Iterator &SxMultiMap<K,V,H>::Iterator::operator++ ()
{
   ++it;
   return *this;
}


template<class K, class V, class H>
const K &SxMultiMap<K,V,H>::Iterator::getKey ()
{
   return it->node->element.key;
}


template<class K, class V, class H>
V &SxMultiMap<K,V,H>::Iterator::getValue ()
{
   return it->value->element;
}

//template<class K, class V, class H>
//SxList<V> &SxMultiMap<K,V,H>::Iterator::getValues ()
//{
//   return it->node->element.values;
//}

template<class K, class V, class H>
V &SxMultiMap<K,V,H>::Iterator::operator-> ()
{
   return it->value->element;
}



// --------------------------------------------------------------------------
//
//template<class K, class V, class H>
//size_t getNBytes (const SxMultiMap<K,V,H> &in)
//{
//   return in.getNBytes ();
//}
