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

#ifndef _SX_MULTI_MAP_H_
#define _SX_MULTI_MAP_H_

#include <SxException.h>
#include <SxHashTable.h>

/** \brief Maps keys to values

    \b SxMultiMap = SPHInX Key-Value Map

    This template class creates key-value pairs and manages them.
    The data types of both key and value are template arguments. It can
    be instantiated for all data types which define
    - a default constructor
    - a copy constructor
    - an operator==
    
    \par Example:

\code
   SxMultiMap<SxString,double> hash;
   hash("PI") = 3.14;
   hash("2PI") = 6.28;
   SxMultiMap<SxString,double>::Iterator it;
   for (it = hash.begin(); it != hash.end(); ++it)  {
      cout << it->getKey() << " = " << it->getValue() << endl;
   }
   cout << "2pi = " << hash("2pi") << endl;
\endcode

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template<class K, class V, class H=SxHashFunction>
class SxMultiMap
{
   protected:
      class Pair;
      
   public:
   
      typedef K  Key;
      typedef V  Value;
      
      /** \brief Iterator for SxMultiMap objecs

          Using this class it is possible to iterate over all key-value-pairs 
          of an SxMultiMap object.

          \par Example
\code
   SxMultiMap<SxString,double> hash;
   hash("PI") = 3.14;
   hash("2PI") = 6.28;
   SxMultiMap<SxString,double>::Iterator it;
   for (it = hash.begin(); it != hash.end(); ++it)  {
      cout << it.getKey << "=" << it.getValue();
   }
\endcode
           If all keys/values have to be accessed the iteration is significantly
           faster then a direct access.  */
      class ConstIterator
      {
         public:
            /** \brief Default constructor */
            ConstIterator ();
            /** \brief Creates a new iterator */
            ConstIterator (const typename SxList<Pair>::ConstIterator &);
            /** \brief Destroy iterator */
           ~ConstIterator ();

            /** \brief Check if two iterators are equal */
            bool operator== (const typename SxMultiMap<K,V,H>::ConstIterator &)
                 const;
            /** \brief Check if two iterators are not equal */
            bool operator!= (const typename SxMultiMap<K,V,H>::ConstIterator &)
                 const;
            /** \brief Postfix incremental operator of an iterator */
            const typename SxMultiMap<K,V,H>::ConstIterator  operator++ (int);
            /** \brief Prefix incremental operator of an iterator */
                  typename SxMultiMap<K,V,H>::ConstIterator &operator++ ();

            /** \brief Extract the key from the key-value pair */
            const K &getKey () const;
            /** \brief Extract the value from the key-value pair */
            const V &getValue () const;
            /** \brief Dereferencing operator */
            const V &operator->() const;

         protected:
            typename SxList<Pair>::ConstIterator it;
      };
      
      class Iterator
      {
         public:
            /** \brief Default constructor */
            Iterator ();
            /** \brief Creates a new iterator */
            Iterator (const typename SxList<Pair>::Iterator &);
            /** \brief Destroy iterator */
           ~Iterator ();

            /** \brief Check if two iterators are equal */
            bool operator== (const typename SxMultiMap<K,V,H>::Iterator &);
            /** \brief Check if two iterators are not equal */
            bool operator!= (const typename SxMultiMap<K,V,H>::Iterator &);
            /** \brief Postfix incremental operator of an iterator */
            const typename SxMultiMap<K,V,H>::Iterator  operator++ (int);
            /** \brief Prefix incremental operator of an iterator */
                  typename SxMultiMap<K,V,H>::Iterator &operator++ ();

            /** \brief Extract the key from the key-value pair */
            const K &getKey ();
            /** \brief Extract the value from the key-value pair */
            V &getValue ();
            /** \brief Dereferencing operator */
            V &operator->();

         protected:
            typename SxList<Pair>::Iterator it;
      };
      

      /** \brief Default constructor 
       */
      SxMultiMap ();
      SxMultiMap (const SxMultiMap<K,V,H> &);

      /** \brief Destructor 
        */
      ~SxMultiMap ();

      /** \brief Assign a key.

          This is the default way of assigning a key.

          \par Example: 
\code          
   SxMultiMap<SxString,double> hash;
   hash("PI") = 3.14;
\endcode */
      inline V &operator() (const K &key);
      inline const V &operator() (const K &key) const;
      
      // -- Assignment
      SxMultiMap<K,V,H> &operator= (const SxMultiMap<K,V,H> &);
      
      /** Creates 1:N key:value multimap as key:value1,value2,.. */
      SxMultiMap<K,V,H> &append (const K &key, const V &value);
      SxList<V> &getValues (const K &key);
      const SxList<V> &getValues (const K &key) const;
      
      ssize_t findKey (const K &) const;
      SxList<V> &getValuesIdx (ssize_t idx);
      const SxList<V> &getValuesIdx (ssize_t idx) const;
            
      /** \return number of key-value pairs */
      ssize_t getSize () const;
      
      /** \brief Checks whether a certain key is present in the hash */
      bool containsKey (const K &) const;

      /** \brief Checks whether a certain value is present in the hash */
      bool containsValue (const V &) const;

      /** \brief Removes a key-value pair 
        */
      bool removeKey (const K &);

      /** \brief Removes a key-value pair by value
        */
      bool removeValue (const V &);
            
      /** \brief Removes all key-value pairs
        */
      void removeAll ();
      
      // Accessfunctions
      inline V &first ();
      inline const V &first () const;
      inline V &last ();
      inline const V &last () const;

      /** \brief Returns the list of all defined keys in the map.

          \par Example:
\code          
   SxMultiMap<SxString,SxString> env;
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
   SxMultiMap<SxString,SxString> env;
   env("SX_INCLUDE_PATH") = "/home/johndoe";
   env("PATH") = "/bin:/usr/bin";
   cout << env.getValues() << endl;
\endcode
         yields the output "/home/johndoe" and "/bin:/usr/bin".  
         \sa getKeys
       */
      SxList<V> getValues () const;

      /** \brief Create iterator pointing to the begin of the map

          Use this function find the first key-value pair.
          \par Example:
\code
   SxMultiMap<SxString,SxString> hash = ...;
   SxMultiMap<SxString,SxString>::Iterator it;
   for (it = hash.begin(); it != hash.end(); ++it)  {
      cout << (*it).getKey() << " = " << (*it).getValue () << endl;
   }
\endcode */
      typename SxMultiMap<K,V,H>::ConstIterator begin () const;
      typename SxMultiMap<K,V,H>::Iterator begin ();


      /** \brief Create iterator pointing to the end of the map

          Use this function find the last key-value pair.
          \par Example:
\code
   SxMultiMap<SxString,SxString> hash = ...;
   SxMultiMap<SxString,SxString>::Iterator it;
   for (it = hash.begin(); it != hash.end(); ++it)  {
      cout << (*it).getKey() << " = " << (*it).getValue () << endl;
   }
\endcode */
      typename SxMultiMap<K,V,H>::ConstIterator end () const;
      typename SxMultiMap<K,V,H>::Iterator end ();

      /** \brief append map to existing map

          This function extends an existing map. It appends key/value pairs of
          the existing map. If a key exists already in the current map it can 
          either be overwritten or ignored. */
      SxMultiMap<K,V,H> &append (const SxMultiMap<K,V,H> &map,
                                 bool                    overwrite=true);

      /** \brief print the content of a map 

          This function can be used to print all key-value pairs to
          stdout. Useful for debugging purposes. */
      void print () const;

   protected:
      class Node {
         public:
            K         key;
            SxList<V> values;
            SxList<typename SxList<Pair>::List*> sequence;
      };
      class Pair {
         public:
            typename SxList<Node>::List* node;
            typename SxList<V>::List*    value;
      };
      
      SxHashTable<K,Node,H> table;

      SxList<K>    keys;
      SxList<Node> nodes;
      SxList<Pair> sequence;
      
      V &insert (const K &key, const V &value, bool overwrite);
};

#include <SxMultiMap.hpp>

#endif /* _SX_MULTI_MAP_H_ */
