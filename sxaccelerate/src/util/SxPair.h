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

#ifndef _SX_PAIR_H_
#define _SX_PAIR_H_


/** \brief Key Value pair

    Compare operations of SxPairs are based only on the keys.

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template <class K, class V>
class SxPair
{
   public:
      K key;
      V value;

      SxPair ();
      SxPair (const K&, const V&);
      ~SxPair ();

      inline SxPair<K,V> &operator=  (const SxPair<K,V> &);
      inline bool operator== (const SxPair<K,V> &) const;
      inline bool operator!= (const SxPair<K,V> &) const;
      inline bool operator<  (const SxPair<K,V> &) const;
      inline bool operator>  (const SxPair<K,V> &) const;
};

template<class K, class V>
SxPair<K,V>::SxPair ()
   : key(K()), value(V())
{
   // empty
}

template<class K, class V>
SxPair<K,V>::SxPair (const K &key_, const V &value_)
   : key(key_),
     value(value_)
{
   // empty
}

template<class K, class V>
SxPair<K,V>::~SxPair ()
{
   // empty
}

template<class K, class V>
SxPair<K,V> &SxPair<K,V>::operator= (const SxPair<K,V> &in)
{
   if (this == &in) return *this;
   key   = in.key;
   value = in.value;
   return *this;
}

template<class K, class V>
bool SxPair<K,V>::operator== (const SxPair<K,V> &in) const
{
   return (key == in.key);
}

template<class K, class V>
bool SxPair<K,V>::operator!= (const SxPair<K,V> &in) const
{
   return (key != in.key);
}

template<class K, class V>
bool SxPair<K,V>::operator< (const SxPair<K,V> &in) const
{
   return (key < in.key);
}

template<class K, class V>
bool SxPair<K,V>::operator> (const SxPair<K,V> &in) const
{
   return (key > in.key);
}

template<class K, class V>
std::ostream &operator<< (std::ostream &s, const SxPair<K,V> &in)
{
   s << "[" << in.key << ", " << in.value << "]";
   return s;
}

template<class K, class V>
size_t getNBytes (const SxPair<K,V> &in)
{
   size_t nBytes = sizeof (in)
                 + getNBytes (in.key)   - sizeof (in.key)
                 + getNBytes (in.value) - sizeof (in.value);

   return nBytes;
}

#endif /* _SX_PAIR_H_ */
