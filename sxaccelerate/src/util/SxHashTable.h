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

#ifndef _SX_HASH_TABLE_H_
#define _SX_HASH_TABLE_H_

#include <SxList.h>
#include <SxHashFunction.h>
#include <SxFunctor.h> /* SxNull */

/** \brief Hash table

    \b SxHashTable = SFHIngX Hash table

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template <class K, class V, class H>
class SxHashTable : public SxArray<typename SxList<K>::Node*>
{
   public:
      // --- Constructor 
      SxHashTable ();
      
      // -- Insertion
      bool contains (const K &, ssize_t *);
      inline typename SxList<V>::Node* &value (ssize_t idx)
         { return values.elements[idx]; }
      inline typename SxList<V>::Node* &value (ssize_t idx) const
         { return values.elements[idx]; }
      
      // -- Search
      ssize_t findPos (const K &) const;
      bool contains (const K &) const;
            
      // -- Deletion
      inline void remove (ssize_t);
      inline void removeAll ();
      
      void resize (ssize_t newSize);
      
      void print () const;
      
   protected:
      ssize_t tableSize;
      ssize_t tableUsed;
      ssize_t tableBound;
      
      SxArray<typename SxList<V>::Node*> values;
      
      typename SxList<K>::Node deleted;
};

/** \brief Hash table

    \b SxHashTable = SFHIngX Hash table

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template <class K, class H>
class SxHashTable<K,SxNull,H> : public SxArray<typename SxList<K>::Node*>
{
   public:
      // --- Constructor 
      SxHashTable ();
      
      // -- Insertion
      bool contains (const K &, ssize_t *);
      
      // -- Search
      ssize_t findPos (const K &) const;
      bool contains (const K &) const;
            
      // -- Deletion
      inline void remove (ssize_t);
      inline void removeAll ();
      
      void resize (ssize_t newSize);
      
      void print () const;
      
   protected:
      ssize_t tableSize;
      ssize_t tableUsed;
      ssize_t tableBound;
      
      typename SxList<K>::Node deleted;
};

#include <SxHashTable.hpp>

#endif /* _SX_HASH_TABLE_H_ */
