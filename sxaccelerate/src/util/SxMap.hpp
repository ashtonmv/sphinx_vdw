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
SxMap<K,V,H>::SxMap () : SxThis<SxMap<K,V,H> > ()
{
   // empty
}

template<class K, class V, class H>
SxMap<K,V,H>::SxMap (const SxMap<K,V,H> &map_)
   : SxThis<SxMap<K,V,H> > ()
{
   append (map_, true);
}

template<class K, class V, class H>
SxMap<K,V,H>::SxMap (const std::initializer_list<SxPair<K,V> > &list_)
   : SxThis<SxMap<K,V,H> > ()
{
   append (list_, true);
}

template<class K, class V, class H>
SxMap<K,V,H>::SxMap (SxMap<K,V,H> &&map_)
   : SxThis<SxMap<K,V,H> > ()
{
   keys = std::move (map_.keys);
   values = std::move (map_.values);
   ssize_t idx;
   typename SxList<K>::Node *kNode = keys.firstElement;
   typename SxList<V>::Node *vNode = values.firstElement;
   while (kNode) {
      idx = -1;
      K key = kNode->elem;
      if (!table.contains (key, &idx)) {
         table.elements[idx] = kNode;
         table.value(idx) = vNode;
      }
      kNode = kNode->next;
      vNode = vNode->next;
   }
}

template<class K, class V, class H>
SxMap<K,V,H> &SxMap<K,V,H>::operator= (const SxMap<K,V,H> &map_)
{
   if (this == &map_)  return *this;
   removeAll ();
   append (map_, true);
   return *this;
}

template<class K, class V, class H>
SxMap<K,V,H> &SxMap<K,V,H>::operator= (SxMap<K,V,H> &&map_)
{
   if (this == &map_)  return *this;
   removeAll ();
   keys = std::move (keys);
   values = std::move (values);
   ssize_t idx;
   typename SxList<K>::Node *kNode = keys.firstElement;
   typename SxList<V>::Node *vNode = values.firstElement;
   while (kNode) {
      idx = -1;
      K key = kNode->elem;
      if (!table.contains (key, &idx)) {
         table.elements[idx] = kNode;
         table.value(idx) = vNode;
      }
      kNode = kNode->next;
      vNode = vNode->next;
   }
   return *this;
}

template<class K, class V, class H>
SxMap<K,V,H> &
SxMap<K,V,H>::operator= (const std::initializer_list<SxPair<K,V> > &list_)
{
   removeAll ();
   append (list_, true);
   return *this;
}

template<class K, class V, class H>
SxMap<K,V,H>::~SxMap ()
{
   removeAll ();
}

template<class K, class V, class H>
inline V &SxMap<K,V,H>::operator() (const K &key_)
{
   ssize_t idx = -1;
   if (!table.contains (key_, &idx))  {
      keys.append (key_);
      table.elements[idx] = keys.lastElement;
      values.append (V ());
      table.value(idx) = values.lastElement;
      return values.lastElement->elem;
   }  else  {
      return table.value(idx)->elem;
   }
}

template<class K, class V, class H>
inline const V &SxMap<K,V,H>::operator() (const K &key_) const
{
   ssize_t idx = table.findPos (key_);
   if (idx >= 0)  {
      return table.value(idx)->elem;
   }
   throw SxException("The map has no element for this key.",__FILE__,__LINE__);
}

template<class K, class V, class H>
ssize_t SxMap<K,V,H>::findKey (const K &key_) const
{
   return table.findPos (key_);
}

template<class K, class V, class H>
V &SxMap<K,V,H>::getValueIdx (ssize_t idx_)
{
   return table.value(idx_)->elem;
}

template<class K, class V, class H>
const V &SxMap<K,V,H>::getValueIdx (ssize_t idx_) const
{
   return table.value(idx_)->elem;
}

template<class K, class V, class H>
ssize_t SxMap<K,V,H>::getSize () const
{
   return keys.getSize ();
}

template<class K, class V, class H>
inline bool SxMap<K,V,H>::containsKey (const K &key_) const
{
   return table.findPos (key_) >= 0;
}

template<class K, class V, class H>
inline bool SxMap<K,V,H>::hasKey (const K &key_) const
{
   return table.findPos (key_) >= 0;
}

template<class K, class V, class H>
bool SxMap<K,V,H>::containsValue (const V &value_) const
{
   typename SxList<V>::ConstIterator it = values.begin ();
   typename SxList<V>::ConstIterator itEnd = values.end ();
   for (; it != itEnd; ++it)  {
      if (*it == value_) return true;
   }
   return false;
}

template<class K, class V, class H>
bool SxMap<K,V,H>::hasValue (const V &value_) const
{
   typename SxList<V>::ConstIterator it = values.begin ();
   typename SxList<V>::ConstIterator itEnd = values.end ();
   for (; it != itEnd; ++it) {
      if (*it == value_) return true;
   }
   return false;
}

template<class K, class V, class H>
bool SxMap<K,V,H>::removeKey (const K &key_)
{
   ssize_t idx = table.findPos (key_);
   if (idx >= 0)  {
      values.removeItem (table.value(idx));
      keys.removeItem (table.elements[idx]);
      table.remove (idx);
      return true;
   }
   return false;
}

template<class K, class V, class H>
bool SxMap<K,V,H>::removeValue (const V &value_)
{
   typename SxList<V>::ConstIterator it = values.begin ();
   typename SxList<K>::ConstIterator itK = keys.begin ();
   typename SxList<V>::ConstIterator itEnd = values.end ();

   for (; it != itEnd; ++it, ++itK)  {
      if (*it == value_)  {
         return removeKey (*itK);
      }
   }
   return false;
}

template<class K, class V, class H>
void SxMap<K,V,H>::removeItem (typename SxMap<K,V,H>::Iterator *it)
{
   if((*it).isForward ()) {
      SxMap<K,V,H>::Iterator tmpIt = (*it)++;
      removeKey (tmpIt.getKey ());
   }
   else {
      SxMap<K,V,H>::Iterator tmpIt = (*it)--;
      removeKey (tmpIt.getKey ());
   }
}

template<class K, class V, class H>
void SxMap<K,V,H>::removeAll ()
{
   table.removeAll ();
   keys.removeAll ();
   values.removeAll ();
}

template<class K, class V, class H>
const SxList<K> &SxMap<K,V,H>::getKeys () const
{
   return keys;
}

template<class K, class V, class H>
const SxList<V> &SxMap<K,V,H>::getValues () const
{
   return values;
}

template<class K, class V, class H>
typename SxMap<K,V,H>::ConstIterator SxMap<K,V,H>::begin () const
{
   return ConstIterator (sx::Forward, this, keys.begin (), values.begin ());
}

template<class K, class V, class H>
typename SxMap<K,V,H>::Iterator SxMap<K,V,H>::begin ()
{
   return Iterator (sx::Forward, this, keys.begin (), values.begin ());
}

template<class K, class V, class H>
typename SxMap<K,V,H>::ConstIterator
SxMap<K,V,H>::begin (const K &key, sx::Direction dir_) const
{
   ConstIterator it = begin ();
   for (; it != end (); ++it) {
      if (it.getKey () == key) {
         if (dir_ == sx::Backward)
            it = it.backward ();
         break;
      }
   }
   return it;
}

template<class K, class V, class H>
typename SxMap<K,V,H>::Iterator 
SxMap<K,V,H>::begin (const K &key, sx::Direction dir_)
{
   Iterator it = begin ();
   for (; it != end (); ++it) {
      if (it.getKey () == key) {
         if (dir_ == sx::Backward)
            it = it.backward ();
         break;
      }
   }
   return it;
}

template<class K, class V, class H>
typename SxMap<K,V,H>::ConstIterator SxMap<K,V,H>::end () const
{
   return ConstIterator (sx::Forward, this, keys.end (), values.end ());
}

template<class K, class V, class H>
typename SxMap<K,V,H>::Iterator SxMap<K,V,H>::end ()
{
   return Iterator (sx::Forward, this, keys.end (), values.end ());
}

template<class K, class V, class H>
typename SxMap<K,V,H>::ConstIterator SxMap<K,V,H>::fromLast () const
{
   return ConstIterator (sx::Backward, this, keys.fromLast (), values.fromLast ());
}

template<class K, class V, class H>
typename SxMap<K,V,H>::Iterator SxMap<K,V,H>::fromLast ()
{
   return Iterator (sx::Backward, this, keys.fromLast (), values.fromLast ());
}

template<class K, class V, class H>
typename SxMap<K,V,H>::ConstIterator SxMap<K,V,H>::toFirst () const
{
   return ConstIterator (sx::Backward, this, keys.toFirst (), values.toFirst ());
}

template<class K, class V, class H>
typename SxMap<K,V,H>::Iterator SxMap<K,V,H>::toFirst ()
{
   return Iterator (sx::Backward, this, keys.toFirst (), values.toFirst ());
}

template<class K, class V, class H>
SxMap<K,V,H> &SxMap<K,V,H>::append (const SxMap<K,V,H> &map, bool overwrite)
{
   ConstIterator it;
   for (it = map.begin(); it != map.end(); ++it)  {
      if (overwrite || !containsKey (it.getKey ()))  {
         operator()(it.getKey ()) = it.getValue ();
      }
   }
   return *this;
}

template<class K, class V, class H>
SxMap<K,V,H> &
SxMap<K,V,H>::append (const std::initializer_list<SxPair<K,V> > &list_,
                      bool overwrite)
{
   for (const SxPair<K,V> &pair : list_)  {
      if (overwrite || !containsKey (pair.key))  {
         operator()(pair.key) = pair.value;
      }
   }
   return *this;
}

template<class K, class V, class H>
typename SxMap<K,V,H>::Iterator
SxMap<K,V,H>::getIterator (const typename SxMap<K,V,H>::SelIdx &idx)
{
   return begin (*idx);
}

template<class K, class V, class H>
typename SxMap<K,V,H>::ConstIterator
SxMap<K,V,H>::getIterator (const typename SxMap<K,V,H>::SelIdx &idx) const
{
   return begin (*idx);
}

template<class K, class V, class H>
void SxMap<K,V,H>::print () const
{
   ConstIterator it;
   for (it = begin(); it != end(); ++it)  {
      std::cout << it.getKey () << " = " << it.getValue () << std::endl;
   }
}

template<class K, class V, class H>
SxSelection<SxMap<K,V,H> >
SxMap<K,V,H>::where (SxCBoundPtr<bool,typename SxMap<K,V,H>::ConstIterator> f)
{
   SxCList<SelIdx> lst;
   ConstIterator it = begin ();
   for (;it != end (); ++it) {
      if (f (it)) {
         lst.append (&(it.getKey ()));
      }
   }
   SxSelection<SxMap<K,V,H> > sel(lst);
   sel.container = this->getThis ();
   return sel;
}

// --------------------------------------------------------------------------
/** \brief Maps keys to values

    \b SxMap = SPHInX Key-Value Map

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template<class K, class V>
class SxMap<K,V,SxNull> : public SxThis<SxMap<K,V,SxNull> >
{
   public:
      typedef K  Key;
      typedef V  Value;
      typedef SxMap<K,V,SxNull> Container;

      // --- SxSelection:
      typedef V               TElem;
      typedef const K        *SelIdx;
      typedef SxList<SelIdx>  SelStorage;
      typedef Container       SelContainer;


      template<class SKey,class SValue, class SContainer, class IT,
               class ITK, class ITV>
      class State
      {
         SX_ITERATOR_STATE
         // friend to State<const ...>
         template<class CK,class CV,class CC,class CIT,class CITK, class CITV>
         friend class State;

         public:

            State () {}

            State (sx::Direction dir_,
                   SContainer *container_, ITK itK_,
                   ITV itV_)
               :   dir(dir_), container(container_), itK(itK_), itV(itV_) { }

            State (const State &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), itK(in.itK), itV(in.itV)
            {
               SX_UNUSED (cmode_);
            }

            State (State &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container),
                 itK(std::move(in.itK)), itV(std::move(in.itV))
            {
               SX_UNUSED (cmode_);
               in.dir = sx::Undefined; in.container = NULL;
               in.itK = ITK(); in.itV = ITV();
            }

            // non-const to const cast
            template<class CK,class CV, class CC,class CI,class CITK,class CITV>
            State (const State<CK,CV,CC,CI,CITK,CITV> &in)
                  : dir(in.dir), container(in.container), itK(in.itK),
                    itV(in.itV) { }
            SKey &getKey () { return *itK; }
            SValue &getValue () { return *itV; }
            typename Container::SelIdx getSelIdx () const {
               return &(*itK);
            }

         protected:
            SContainer *container;
            ITK itK;
            ITV itV;

            // --- SX_ITERATOR / SX_CONST_ITERATOR callbacks
            void copy (const IT &in_) {
               container = in_.container; dir = in_.dir;
               itK = in_.itK; itV = in_.itV;
            }
            void move (IT &&in_) {
               dir = in_.dir; container = in_.container;
               itK = std::move(in_.itK); itV = std::move(in_.itV);
               in_.dir = sx::Undefined; in_.container = NULL;
               in_.itK = ITK(); in_.itV = ITV();
            }
            void next () {
               if (itK.isValid ()) {
                  ++itK;
                  ++itV;
               }
            }
            void prev () {
               if (itK.isValid ()) {
                  --itK;
                  --itV;
               }
            }
            bool valid () const { return (itK.isValid () && itV.isValid ()); }
            SValue &getRef () { SX_CHECK (itV.isValid ()); return  *itV; }
            SValue *getPtr () { SX_CHECK (itV.isValid ()); return &(*itV); }
            bool equal (const IT &in) const { return itK == in.itK; }

      };

      class Iterator;

      class ConstIterator
            : public State<const K, const V, const Container, ConstIterator,
                           typename SxList<K>::ConstIterator,
                           typename SxList<V>::ConstIterator>
      {
         SX_CONST_ITERATOR_NO_LAMBDAS(V,Container,ConstIterator)
         public:
            /** \brief Default constructor */
            ConstIterator () : State<const K, const V, const Container,
                                     ConstIterator,
                                     typename SxList<K>::ConstIterator,
                                     typename SxList<V>::ConstIterator
                                    > () { };
            /** \brief Typecast constructor */
            ConstIterator (const Iterator &it)
                          : State<const K, const V, const Container,
                                  ConstIterator,
                                  typename SxList<K>::ConstIterator,
                                  typename SxList<V>::ConstIterator
                                 > (it) { };
            /** \brief Creates a new iterator */
            ConstIterator (sx::Direction dir_,
                           const Container *c,
                           typename SxList<K>::ConstIterator itK_,
                           typename SxList<V>::ConstIterator itV_)
                          : State<const K, const V, const Container,
                                  ConstIterator,
                                  typename SxList<K>::ConstIterator,
                                  typename SxList<V>::ConstIterator
                                 > (dir_, c, itK_, itV_) { };

            ConstIterator (const ConstIterator &in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const K, const V, const Container,ConstIterator,
                       typename SxList<K>::ConstIterator,
                       typename SxList<V>::ConstIterator> (in, cmode_) { }

            ConstIterator (ConstIterator &&in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const K, const V, const Container, ConstIterator,
                       typename SxList<K>::ConstIterator,
                       typename SxList<V>::ConstIterator>
                       (std::move(in), cmode_) { }

            /** \brief Destroy iterator */
           ~ConstIterator () { };
      };

      class Iterator
            : public State< const K, V, Container, Iterator,
                            typename SxList<K>::Iterator,
                            typename SxList<V>::Iterator>

      {
         SX_ITERATOR_NO_LAMBDAS(V,Container,Iterator)
         public:
            /** \brief Default constructor */
            Iterator () : State<const K, V, Container, Iterator,
                                typename SxList<K>::Iterator,
                                typename SxList<V>::Iterator
                               > () { };
            /** \brief Creates a new iterator */
            Iterator (sx::Direction dir_,
                      Container *c,
                      typename SxList<K>::Iterator itK_,
                      typename SxList<V>::Iterator itV_)
                      : State<const K, V, Container, Iterator,
                              typename SxList<K>::Iterator,
                              typename SxList<V>::Iterator
                             > (dir_, c, itK_, itV_) { };

            Iterator (const Iterator &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const K, V, Container,Iterator,
                       typename SxList<K>::Iterator,
                       typename SxList<V>::Iterator> (in, cmode_) { }

            Iterator (Iterator &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const K, V, Container, Iterator,
                       typename SxList<K>::Iterator,
                       typename SxList<V>::Iterator> (std::move(in), cmode_) { }

            /** \brief Destroy iterator */
            ~Iterator () { };

            friend class ConstIterator;

      };


      SxMap ();
      SxMap (const SxMap<K,V,SxNull> &);
      SxMap (SxMap<K,V,SxNull> &&);
      SxMap (const std::initializer_list<SxPair<K,V> > &);
      ~SxMap ();

      inline V &operator() (const K &key);
      inline const V &operator() (const K &key) const;

      SxMap<K,V,SxNull> &operator= (const SxMap<K,V,SxNull> &);
      SxMap<K,V,SxNull> &operator= (SxMap<K,V,SxNull> &&);
      SxMap<K,V,SxNull> &operator= (const std::initializer_list<SxPair<K,V>> &);

      ssize_t findKey (const K &) const;
            V &getValueIdx (ssize_t idx);
      const V &getValueIdx (ssize_t idx) const;

      inline Iterator      getIterator (const SelIdx &idx);
      inline ConstIterator getIterator (const SelIdx &idx) const;

      bool containsKey (const K &) const;
      bool hasKey (const K &) const;
      bool containsValue (const V &) const;
      bool hasValue (const V &) const;

      ssize_t getSize () const { return keys.getSize (); }

      bool removeKey (const K &);
      bool removeValue (const V &);
      void removeItem (typename SxMap<K,V,SxNull>::Iterator *it);
      void removeAll ();

      const SxList<K> &getKeys () const { return keys; }
      const SxList<V> &getValues () const { return values; }

      typename SxMap<K,V,SxNull>::ConstIterator begin () const;
      typename SxMap<K,V,SxNull>::Iterator begin ();

      typename SxMap<K,V,SxNull>::ConstIterator begin
               (const K &key, sx::Direction dir_ = sx::Forward) const;
      typename SxMap<K,V,SxNull>::Iterator begin
               (const K &key, sx::Direction dir_ = sx::Forward);

      typename SxMap<K,V,SxNull>::ConstIterator end () const;
      typename SxMap<K,V,SxNull>::Iterator end ();

      typename SxMap<K,V,SxNull>::ConstIterator fromLast () const;
      typename SxMap<K,V,SxNull>::Iterator fromLast ();
      typename SxMap<K,V,SxNull>::ConstIterator toFirst () const;
      typename SxMap<K,V,SxNull>::Iterator toFirst ();

      SxMap<K,V,SxNull> &append (const SxMap<K,V,SxNull> &map,
                                 bool                    overwrite=true);

      SxMap<K,V,SxNull> &append (const std::initializer_list<SxPair<K,V>> &list,
                                 bool overwrite=true);

      void print () const;

      SxSelection<SxMap<K,V,SxNull> > where (SxCBoundPtr<bool,
                                             typename SxMap<K,V,SxNull>::ConstIterator> f);


   protected:
      SxList<K> keys;
      SxList<V> values;
};

template<class K, class V>
SxMap<K,V,SxNull>::SxMap () : SxThis<SxMap<K,V,SxNull> > ()
{
   // empty
}

template<class K, class V>
SxMap<K,V,SxNull>::SxMap (const SxMap<K,V,SxNull> &map_)
   : SxThis<SxMap<K,V,SxNull> > ()
{
   append (map_, true);
}

template<class K, class V>
SxMap<K,V,SxNull>::SxMap (SxMap<K,V,SxNull> &&map_)
   : SxThis<SxMap<K,V,SxNull> > ()
{
   keys   = std::move (map_.keys);
   values = std::move (map_.values);
}

template<class K, class V>
SxMap<K,V,SxNull>::SxMap (const std::initializer_list<SxPair<K,V> > &list_)
   : SxThis<SxMap<K,V,SxNull> > ()
{
   append (list_, true);
}

template<class K, class V>
SxMap<K,V,SxNull> &SxMap<K,V,SxNull>::operator= (const SxMap<K,V,SxNull> &map_)
{
   if (this == &map_)  return *this;
   removeAll ();
   return append (map_, true);
}

template<class K, class V>
SxMap<K,V,SxNull> &SxMap<K,V,SxNull>::operator= (SxMap<K,V,SxNull> &&map_)
{
   if (this == &map_)  return *this;
   keys   = std::move (map_.keys);
   values = std::move (map_.values);
}

template<class K, class V>
SxMap<K,V,SxNull> &
SxMap<K,V,SxNull>::operator= (const std::initializer_list<SxPair<K,V> > &list_)
{
   removeAll ();
   return append (list_, true);
}

template<class K, class V>
SxMap<K,V,SxNull>::~SxMap ()
{
   removeAll ();
}

template<class K, class V>
inline V &SxMap<K,V,SxNull>::operator() (const K &key_)
{
   typename SxList<K>::Iterator itKey = keys.begin ();
   typename SxList<K>::Iterator itEnd = keys.end ();
   typename SxList<V>::Iterator itValue = values.begin ();

   for (; itKey != itEnd; ++itKey, ++itValue)  {
      if (*itKey == key_)  {
         return *itValue;
      }
   }

   keys << key_;
   values << V();

   return values.lastElement->elem;
}

template<class K, class V>
inline const V &SxMap<K,V,SxNull>::operator() (const K &key_) const
{
   typename SxList<K>::ConstIterator itKey = keys.begin ();
   typename SxList<K>::ConstIterator itEnd = keys.end ();
   typename SxList<V>::ConstIterator itValue = values.begin ();

   for (; itKey != itEnd; ++itKey, ++itValue)  {
      if (*itKey == key_)  {
         return *itValue;
      }
   }

   throw SxException("The map has no element for this key.",__FILE__,__LINE__);
}

template<class K, class V>
ssize_t SxMap<K,V,SxNull>::findKey (const K &key_) const
{
   return keys.findPos (key_);
}

template<class K, class V>
V &SxMap<K,V,SxNull>::getValueIdx (ssize_t idx_)
{
   return values(idx_);
}

template<class K, class V>
const V &SxMap<K,V,SxNull>::getValueIdx (ssize_t idx_) const
{
   return values(idx_);
}

template<class K, class V>
typename SxMap<K,V,SxNull>::Iterator
SxMap<K,V,SxNull>::getIterator (const typename SxMap<K,V,SxNull>::SelIdx &idx)
{
   return begin (*idx);
}

template<class K, class V>
typename SxMap<K,V,SxNull>::ConstIterator
SxMap<K,V,SxNull>::getIterator (const typename
                                SxMap<K,V,SxNull>::SelIdx &idx) const
{
   return begin (*idx);
}

template<class K, class V>
inline bool SxMap<K,V,SxNull>::containsKey (const K &key_) const
{
   return keys.contains (key_);
}

template<class K, class V>
inline bool SxMap<K,V,SxNull>::hasKey (const K &key_) const
{
   return keys.contains (key_);
}

template<class K, class V>
bool SxMap<K,V,SxNull>::containsValue (const V &value_) const
{
   return values.contains (value_);
}

template<class K, class V>
bool SxMap<K,V,SxNull>::hasValue (const V &value_) const
{
   return values.contains (value_);
}

template<class K, class V>
bool SxMap<K,V,SxNull>::removeKey (const K &key_)
{
   ssize_t idx = keys.findPos (key_);
   if (idx >= 0)  {
      keys.remove (idx);
      values.remove (idx);
      return true;
   }
   return false;
}

template<class K, class V>
bool SxMap<K,V,SxNull>::removeValue (const V &value_)
{
   ssize_t idx = values.findPos (value_);
   if (idx >= 0)  {
      keys.remove (idx);
      values.remove (idx);
      return true;
   }
   return false;
}

template<class K,class V>
void SxMap<K,V,SxNull>::removeItem (typename SxMap<K,V,SxNull>::Iterator *it)
{
   if((*it).isForward ()) {
      SxMap<K,V,SxNull>::Iterator tmpIt = (*it)++;
      removeKey (tmpIt.getKey ());
   } else {
      SxMap<K,V,SxNull>::Iterator tmpIt = (*it)--;
      removeKey (tmpIt.getKey ());
   }
}

template<class K, class V>
void SxMap<K,V,SxNull>::removeAll ()
{
   keys.removeAll ();
   values.removeAll ();
}

template<class K, class V>
typename SxMap<K,V,SxNull>::ConstIterator SxMap<K,V,SxNull>::begin () const
{
   return ConstIterator (sx::Forward, this, keys.begin (), values.begin ());
}

template<class K, class V>
typename SxMap<K,V,SxNull>::Iterator SxMap<K,V,SxNull>::begin ()
{
   return Iterator (sx::Forward, this, keys.begin (), values.begin ());
}

template<class K, class V>
typename SxMap<K,V,SxNull>::ConstIterator 
SxMap<K,V,SxNull>::begin (const K &key, sx::Direction dir_) const
{
   ConstIterator it = begin ();
   for (; it != end (); ++it) {
      if (it.getKey () == key) {
         if (dir_ == sx::Backward)
            it = it.backward ();
         break;
      }
   }
   return it;
}

template<class K, class V>
typename SxMap<K,V,SxNull>::Iterator 
SxMap<K,V,SxNull>::begin (const K &key, sx::Direction dir_)
{
   Iterator it = begin ();
   for (; it != end (); ++it) {
      if (it.getKey () == key) {
         if (dir_ == sx::Backward)
            it = it.backward ();
         break;
      }
   }
   return it;
}

template<class K, class V>
typename SxMap<K,V,SxNull>::ConstIterator SxMap<K,V,SxNull>::end () const
{
   return ConstIterator (sx::Forward, this, keys.end (), values.end ());
}

template<class K, class V>
typename SxMap<K,V,SxNull>::Iterator SxMap<K,V,SxNull>::end ()
{
   return Iterator (sx::Forward, this, keys.end (), values.end ());
}

template<class K, class V>
typename SxMap<K,V,SxNull>::ConstIterator SxMap<K,V,SxNull>::fromLast () const
{
   return ConstIterator (sx::Backward, this, keys.fromLast (), values.fromLast ());
}

template<class K, class V>
typename SxMap<K,V,SxNull>::Iterator SxMap<K,V,SxNull>::fromLast ()
{
   return Iterator (sx::Backward, this, keys.fromLast (), values.fromLast ());
}

template<class K, class V>
typename SxMap<K,V,SxNull>::ConstIterator SxMap<K,V,SxNull>::toFirst () const
{
   return ConstIterator (sx::Backward, this, keys.toFirst (), values.toFirst ());
}

template<class K, class V>
typename SxMap<K,V,SxNull>::Iterator SxMap<K,V,SxNull>::toFirst ()
{
   return Iterator (sx::Backward, this, keys.toFirst (), values.toFirst ());
}

template<class K, class V>
SxMap<K,V,SxNull> &SxMap<K,V,SxNull>::append (const SxMap<K,V,SxNull> &map_,
                                              bool                   overwrite_)
{
   ConstIterator it = map_.begin ();
   ConstIterator itEnd = map_.end ();
   for (; it != itEnd; ++it)  {
      if (overwrite_ || !containsKey (it.getKey ()))  {
         operator()(it.getKey ()) = it.getValue ();
      }
   }
   return *this;
}

template<class K, class V>
SxMap<K,V,SxNull> &
SxMap<K,V,SxNull>::append (const std::initializer_list<SxPair<K,V> > &list_,
                           bool                                      overwrite_)
{
   for (const SxPair<K,V> &pair : list_)  {
      if (overwrite_ || !containsKey (pair.key))  {
         operator()(pair.key) = pair.value;
      }
   }
   return *this;
}

template<class K, class V>
void SxMap<K,V,SxNull>::print () const
{
   ConstIterator it;
   for (it = begin(); it != end(); ++it)  {
      std::cout << it.getKey () << " = " << it.getValue () << std::endl;
   }
}

template<class K, class V>
SxSelection<SxMap<K,V,SxNull> >
SxMap<K,V,SxNull>::where (SxCBoundPtr<bool,
                          typename SxMap<K,V,SxNull>::ConstIterator> f)
{
   SxCList<SelIdx> lst;
   ConstIterator it = begin ();
   for (;it != end (); ++it) {
      if (f (it)) {
         lst.append (&(it.getKey()));
      }
   }
   SxSelection<SxMap<K,V,SxNull> > sel(lst);
   sel.container = this->getThis ();
   return sel;
}
