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

#ifndef _SX_SET_2_H_
#define _SX_SET_2_H_

#include <SxArray.h>
#include <SxException.h>
#include <SxHashFunction.h>

/** \brief Set of elements

    This Hash table does not keep the insert order. The order of elements
    might change after any append operation which can cause re-hashing of
    the whole table.

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template <class K, class H=SxHashFunction>
class SxSet
{
   public:
      // --- Constructor
      SxSet ();
      ~SxSet ();

      // -- Insertion
      inline SxSet<K,H> &operator<< (const K &);
      inline SxSet<K,H> &operator<< (const SxSet<K,H> &);

      /** Returns true only if there was no such element before. */
      bool append (const K &);
      bool append (const SxSet<K,H> &);

      // -- Search
      bool contains (const K &) const;

      // -- Deletion
      inline void removeElement (const K &);
      inline void removeAll ();

      // -- Size
      inline ssize_t getSize () const { return size; }
      size_t getNBytes () const;

      /** Converts this set to a SxList object.
          Does not keep the insert order of elements. */
      inline SxList<K> toList () const;

   protected:
      ssize_t size;  //! Number of valid elements
      ssize_t nUsed; //! Number of valid and deleted elements
      ssize_t bound; //! Limit Number of elements for the next resize

      SxArray<K>    table; //! The elements spread at hashed indices
      SxArray<char> flag;  //! 'u':Used, 'e':Empty, 'd':Deleted

      ssize_t findPos (const K &) const;

      void resize (ssize_t newSize);
};

template<class K, class H>
SxSet<K,H>::SxSet ()
   : size(0),
     nUsed(0),
     bound(0)
{
   // empty
}

template<class K, class H>
SxSet<K,H>::~SxSet ()
{
   removeAll ();
}

template<class K, class H>
SxSet<K,H> &SxSet<K,H>::operator<< (const K &in)
{
   append (in);
   return *this;
}

template<class K, class H>
SxSet<K,H> &SxSet<K,H>::operator<< (const SxSet<K,H> &in)
{
   append (in);
   return *this;
}

template<class K, class H>
bool SxSet<K,H>::append (const K &elem_)
{
   if (nUsed >= bound)  {
      if (size * 2 < table.size)  {
         resize (table.size - 1);
      }  else  {
         resize (table.size + 1);
      }
   }

   size_t h = H::hash (elem_);
   ssize_t idx = (ssize_t)h % table.size;
   ssize_t step = (ssize_t)h % (table.size - 1) + 1;
   ssize_t lastIdx = idx;
   ssize_t replaceIdx = table.size;

   do {
      if (flag.elements[idx] == 'e')  { // Empty
         // --- insert new element
         if (replaceIdx != table.size)  {
            idx = replaceIdx;
         }  else  {
            nUsed++;
         }
         size++;
         table.elements[idx] = elem_;
         flag.elements[idx] = 'u'; // Used
         return true;
      }  else if (flag.elements[idx] == 'd')  { // Deleted
         replaceIdx = idx;
      }  else if (table.elements[idx] == elem_)  {
         return false;
      }
      idx = (idx + step) % table.size;
   } while (idx != lastIdx);

   SX_DBG_MSG ("SxSet::append: Can't insert");
   return false;
}

template<class K, class H>
bool SxSet<K,H>::append (const SxSet<K,H> &in)
{
   bool res = false;
   for (ssize_t i=0; i < in.table.size; ++i)  {
      if (in.flag.elements[i] == 'u' && append (in.table.elements[i]))  {
         res = true;
      }
   }

   return res;
}

template<class K, class H>
ssize_t SxSet<K,H>::findPos (const K &elem_) const
{
   if (size > 0)  {
      size_t h = H::hash (elem_);
      ssize_t idx = (ssize_t)h % table.size;
      ssize_t step = (ssize_t)h % (table.size - 1) + 1;
      ssize_t lastIdx = idx;

      do {
         if (flag.elements[idx] == 'e')  { // Empty
            return -1;
         }  else if (flag.elements[idx]  == 'u' && // Used
                     table.elements[idx] == elem_)
         {
            return idx;
         }
         idx = (idx + step) % table.size;
      } while (idx != lastIdx);
   }
   return -1;
}

template<class K, class H>
bool SxSet<K,H>::contains (const K &elem_) const
{
   return findPos (elem_) >= 0;
}

template<class K, class H>
void SxSet<K,H>::removeElement (const K &elem_)
{
   ssize_t idx = findPos (elem_);
   if (idx >= 0)  {
      flag.elements[idx] = 'd'; // Deleted
      table.elements[idx] = K ();
      size--;
   }
}

template<class K, class H>
void SxSet<K,H>::removeAll ()
{
   table.resize (0);
   flag.resize (0);
   size = 0;
   nUsed = 0;
   bound = 0;
}

template<class K, class H>
SxList<K> SxSet<K,H>::toList () const
{
   SxList<K> res;
   for (ssize_t i=0; i < table.size; ++i)  {
      if (flag.elements[i] == 'u')  {
         res << table.elements[i];
      }
   }
   return res;
}

template<class K, class H>
void SxSet<K,H>::resize (ssize_t newSize_)
{
   SX_CHECK (newSize_ < 1610612741, newSize_);

   // --- http://planetmath.org/encyclopedia/GoodHashTablePrimes.html
   ssize_t prime[30] = {0, 3, 11, 23, 53, 97, 193, 389, 769, 1543, 3079, 6151,
                        12289, 24593, 49157, 98317, 196613, 393241, 786433,
                        1572869, 3145739, 6291469, 12582917, 25165843, 50331653,
                        100663319, 201326611, 402653189, 805306457, 1610612741};
   ssize_t idx = 0;
   while (prime[idx] <= newSize_)  {
      idx++;
   }
   newSize_ = prime[idx];
   ssize_t newBound = (ssize_t)((double)newSize_ * 0.77 + 0.5);
   ssize_t n = table.size;
   ssize_t i;

   if (size < newBound)  {
      SxArray<K> newTable(newSize_);
      SxArray<char> newFlag(newSize_);
      if (newSize_ > 0)  {
         newFlag.set ('e'); // Empty
      }

      for (i = 0; i < n; i++)  {
         if (flag.elements[i] == 'u')  {
            size_t h = H::hash (table.elements[i]);
            idx = (ssize_t)h % newSize_;
            ssize_t step = (ssize_t)h % (newSize_ - 1) + 1;
            ssize_t lastIdx = idx;

            while (newFlag.elements[idx] != 'e')  {
               idx = (idx + step) % newSize_;
               if (idx == lastIdx)  {
                  SX_EXIT;
               }
            }
            newTable.elements[idx] = table.elements[i];
            newFlag.elements[idx] = 'u'; //Used;
         }
      }

      nUsed = size;
      bound = newBound;

      table.resize (newSize_, false);
      size_t len = sizeof(K) * static_cast<size_t>(table.size);
      ::memcpy (table.elements, newTable.elements, len);

      flag.resize (newSize_, false);
      len = static_cast<size_t>(flag.size);
      ::memcpy (flag.elements, newFlag.elements, len);
   }
}

template<class K, class H>
size_t SxSet<K,H>::getNBytes () const
{
   size_t nBytes = sizeof(*this)
                 + ::getNBytes (table) - sizeof (table)
                 + ::getNBytes (flag)  - sizeof (flag);

   return nBytes;
}

template<class K, class H>
size_t getNBytes (const SxSet<K,H> &in)
{
   return in.getNBytes ();
}

template<class K, class H>
std::ostream &operator<< (std::ostream &s, const SxSet<K,H> &in)
{
   s << in.toList ();
   return s;
}

#endif /* _SX_SET_2_H_ */
