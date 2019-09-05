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

#ifndef _SX_SET_H_
#define _SX_SET_H_

#include <SxList.h>
#include <SxHashTable.h>


/** \brief List of elements, where each element is unique

   \b SxUniqueList = SFHIngX List containing unique elements only

    If SxUniqueList<T> fails to compile use SxUniqueList<T,SxNull>
    or add a hash function to SxHashFunction.h.

    The order of keys in SxUniqueList is according to the insert order.

    \author C. Freysoldt, freyso@fhi-berlin.mpg.de
    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template <class T, class H=SxHashFunction>
class SxUniqueList : public SxList<T> 
{
   public:
      /// Constructor 
      SxUniqueList ();
      SxUniqueList (const SxList<T> &);
      SxUniqueList (const SxUniqueList<T> &);
      SxUniqueList (SxUniqueList<T> &&);
      SxUniqueList (const std::initializer_list<T> &);

      // -- Assignment
      SxUniqueList<T,H> &operator= (const SxList<T> &);
      SxUniqueList<T,H> &operator= (const SxUniqueList<T> &);
      SxUniqueList<T,H> &operator= (SxUniqueList<T> &&);
      SxUniqueList<T,H> &operator= (const std::initializer_list<T> &);

      // -- Insertion
      inline SxUniqueList<T,H> &operator<< (const T &);
      inline SxUniqueList<T,H> &operator<< (const SxList<T> &);
      inline SxUniqueList<T,H> &operator<< (const std::initializer_list<T> &);

      SxUniqueList<T,H> &append (const T &);
      SxUniqueList<T,H> &append (const SxList<T> &);
      SxUniqueList<T,H> &append (const std::initializer_list<T> &);
      SxUniqueList<T,H> &insert (ssize_t newPos, const T &);
      SxUniqueList<T,H> &prepend (const T &);

      // -- Search
      ssize_t findPos (const T &) const; // slow
      bool contains (const T &) const;

      // -- Deletion
      inline void remove (ssize_t);
      inline void removeElement (const T &);
      inline void removeFirst ();
      inline void removeLast ();
      inline void removeAll ();

      void print () const;

   protected:
      SxHashTable<T,SxNull,H> table;
};

#include <SxUniqueList.hpp>

#endif /* _SX_SET_H_ */
