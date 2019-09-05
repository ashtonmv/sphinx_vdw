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

#ifndef _SX_MAP_H_
#define _SX_MAP_H_

#include <SxException.h>
#include <SxHashTable.h>
#include <SxSelection.h>
#include <SxPtr.h>
#include <SxPair.h>

/** \brief Maps keys to values

    \b SxMap = SPHInX Key-Value Map

    This template class creates key-value pairs and manages them.
    The data types of both key and value are template arguments. It can
    be instantiated for all data types which define
    - a default constructor
    - a copy constructor
    - an operator==

    \par Example:

\code
   SxMap<SxString,double> hash;
   hash("PI") = 3.14;
   hash("2PI") = 6.28;
   SxMap<SxString,double>::Iterator it;
   for (it = hash.begin(); it != hash.end(); ++it)  {
      cout << it->getKey() << " = " << it->getValue() << endl;
   }
   cout << "2pi = " << hash.getValue("2pi") << endl;
\endcode

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template<class K, class V, class H=SxHashFunction>
class SxMap : public SxThis<SxMap<K,V,H> >
{
   public:

      typedef K  Key;
      typedef V  Value;
      typedef SxMap<K,V,H> Container;

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
                 : dir(dir_), container(container_), itK(itK_), itV(itV_) { }

            State (const State &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), itK(in.itK), itV(in.itV)
            {
               SX_UNUSED (cmode_);
            }

            State (State &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), itK(std::move(in.itK)),
                 itV(std::move(in.itV))
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

      /** \brief Iterator for SxMap objecs

          Using this class it is possible to iterate over all key-value-pairs
          of an SxMap object.

          \par Example
\code
   SxMap<SxString,double> hash;
   hash("PI") = 3.14;
   hash("2PI") = 6.28;
   SxMap<SxString,double>::Iterator it;
   for (it = hash.begin(); it != hash.end(); ++it)  {
      cout << it.getKey << "=" << it.getValue();
   }
\endcode
           If all keys/values have to be accessed the iteration is significantly
           faster then a direct access.  */


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
            : public State<const K, V, Container, Iterator,
                           typename SxList<K>::Iterator,
                           typename SxList<V>::Iterator>

      {
         SX_ITERATOR_NO_LAMBDAS(V,Container,Iterator)
         public:
            /** \brief Default constructor */
            Iterator () : State<const K, V, Container, Iterator,
                                typename SxList<K>::Iterator,
                                typename SxList<V>::Iterator
                               > () { }

            /** \brief Creates a new iterator */
            Iterator (sx::Direction dir_,
                      Container *c,
                      typename SxList<K>::Iterator itK_,
                      typename SxList<V>::Iterator itV_)
                      : State<const K, V, Container, Iterator,
                              typename SxList<K>::Iterator,
                              typename SxList<V>::Iterator
                             > (dir_, c, itK_, itV_) { }

            Iterator (const Iterator &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const K, V, Container,Iterator,
                       typename SxList<K>::Iterator,
                       typename SxList<V>::Iterator> (in, cmode_) { }

            Iterator (Iterator &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const K, V, Container, Iterator,
                       typename SxList<K>::Iterator,
                       typename SxList<V>::Iterator> (std::move(in), cmode_)
                       { }

            /** \brief Destroy iterator */
            ~Iterator () { };

            friend class ConstIterator;

      };

      /** \brief Default constructor
       */
      SxMap ();
      SxMap (const SxMap<K,V,H> &);
      SxMap (SxMap<K,V,H> &&);
      SxMap (const std::initializer_list<SxPair<K,V> > &);

      /** \brief Destructor
        */
      ~SxMap ();

      /** \brief Assign a key.

          This is the default way of assigning a key.

          \par Example:
\code
   SxMap<SxString,double> hash;
   hash("PI") = 3.14;
\endcode */
      inline V &operator() (const K &key);
      inline const V &operator() (const K &key) const;

      // -- Assignment
      SxMap<K,V,H> &operator= (const SxMap<K,V,H> &);
      SxMap<K,V,H> &operator= (SxMap<K,V,H> &&);
      SxMap<K,V,H> &operator= (const std::initializer_list<SxPair<K,V> > &);

      /** \brief Find hash table index of a key.

          \par Example:
\code
   ssize_t idx = table.findKey ("PI");
   if (idx >= 0)  {
      return table.getValueIdx (idx);
   }
\endcode

          Which is faster than:
\code
   if (table.containsKey ("PI"))  {
      return table("PI");
   }
\endcode */
      ssize_t findKey (const K &) const;
            V &getValueIdx (ssize_t idx);
      const V &getValueIdx (ssize_t idx) const;

      /** return container iterator/ConstIterator */
      inline Iterator      getIterator (const SelIdx &);
      inline ConstIterator getIterator (const SelIdx &) const;

      /** \return number of key-value pairs */
      ssize_t getSize () const;

      /** \brief Checks whether a certain key is present in the hash */
      bool containsKey (const K &) const;
      bool hasKey (const K &) const;

      /** \brief Checks whether a certain value is present in the hash */
      bool containsValue (const V &) const;
      bool hasValue (const V &) const;

      /** \brief Removes a key-value pair
        */
      bool removeKey (const K &);

      /** \brief Removes a key-value pair by value
        */
      bool removeValue (const V &);

      void removeItem (typename SxMap<K,V,H>::Iterator *it);

      /** \brief Removes all key-value pairs
        */
      void removeAll ();

      /** \brief Returns the list of all defined keys in the map.

          \par Example:
\code
   SxMap<SxString,SxString> env;
   env("SX_INCLUDE_PATH") = "/home/johndoe";
   env("PATH") = "/bin:/usr/bin";
   cout << env.getKeys() << endl;
\endcode
         yields the output "SX_INCLUDE_PATH" and "PATH".
         \sa getValues */
      const SxList<K> &getKeys () const;

      /** \brief Returns the list of all defined values in the map.

          \par Example:
\code
   SxMap<SxString,SxString> env;
   env("SX_INCLUDE_PATH") = "/home/johndoe";
   env("PATH") = "/bin:/usr/bin";
   cout << env.getValues() << endl;
\endcode
         yields the output "/home/johndoe" and "/bin:/usr/bin".
         \sa getKeys
       */
      const SxList<V> &getValues () const;

      /** \brief Create iterator pointing to the begin of the map

          Use this function find the first key-value pair.
          \par Example:
\code
   SxMap<SxString,SxString> hash = ...;
   SxMap<SxString,SxString>::Iterator it;
   for (it = hash.begin(); it != hash.end(); ++it)  {
      cout << (*it).getKey() << " = " << (*it).getValue () << endl;
   }
\endcode */
      typename SxMap<K,V,H>::ConstIterator begin () const;
      typename SxMap<K,V,H>::Iterator begin ();

      typename SxMap<K,V,H>::ConstIterator begin
               (const K &key, sx::Direction dir_ = sx::Forward) const;
      typename SxMap<K,V,H>::Iterator begin
               (const K &key, sx::Direction dir_ = sx::Forward);

      /** \brief Create iterator pointing to the end of the map

          Use this function find the last key-value pair.
          \par Example:
\code
   SxMap<SxString,SxString> hash = ...;
   SxMap<SxString,SxString>::Iterator it;
   for (it = hash.begin(); it != hash.end(); ++it)  {
      cout << (*it).getKey() << " = " << (*it).getValue () << endl;
   }
\endcode */
      typename SxMap<K,V,H>::ConstIterator end () const;
      typename SxMap<K,V,H>::Iterator end ();


      typename SxMap<K,V,H>::ConstIterator fromLast () const;
      typename SxMap<K,V,H>::Iterator fromLast ();

      typename SxMap<K,V,H>::ConstIterator toFirst () const;
      typename SxMap<K,V,H>::Iterator toFirst ();


      /** \brief append map to existing map

          This function extends an existing map. It appends key/value pairs of
          the existing map. If a key exists already in the current map it can 
          either be overwritten or ignored. */
      SxMap<K,V,H> &append (const SxMap<K,V,H> &map, bool overwrite=true);

      SxMap<K,V,H> &append (const std::initializer_list<SxPair<K,V> > &list_,
                            bool overwrite=true);

      /** \brief print the content of a map

          This function can be used to print all key-value pairs to
          stdout. Useful for debugging purposes. */
      void print () const;

      SxSelection<SxMap<K,V,H> > where (SxCBoundPtr<bool,
                                        typename SxMap<K,V,H>::ConstIterator> f);


   protected:
      SxHashTable<K,V,H> table;
      SxList<K> keys;
      SxList<V> values;
};

#include <SxMap.hpp>

#endif /* _SX_MAP_H_ */
