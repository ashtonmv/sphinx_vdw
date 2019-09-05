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

#ifndef _SX_ARRAY_LIST_H_
#define _SX_ARRAY_LIST_H_

#include <SxString.h>

/** \brief Dynamic Array

    Growable array with sequence 0,1,2,4,8,16...

    Idx=int can allow to create an array in 16 bytes instead of 24 bytes.

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template <class T, class Idx=ssize_t>
class SxArrayList
{
   public:
      T *elements;
      Idx size;
      Idx sizeAlloc;

      SxArrayList ();
      explicit SxArrayList (ssize_t);
      SxArrayList (const SxArrayList<T,Idx> &);
      SxArrayList (SxArrayList<T,Idx> &&);
      ~SxArrayList ();

      // -- Assignment
      SxArrayList<T,Idx> &operator=  (const SxArrayList<T,Idx> &);
      SxArrayList<T,Idx> &operator=  (SxArrayList<T,Idx> &&);

      /** Standard dereferencing operator for non-const objectes. */
      T &operator() (ssize_t);
      /** Standard dereferencing operator for const objectes. */
      const T &operator() (ssize_t) const;

      // -- Insertion
      ssize_t append (const T &);
      ssize_t prepend (const T &);
      ssize_t binaryInsert (const T &);
      void insert (ssize_t newPos, const T &);

      // -- Search
      ssize_t findPos (const T &, ssize_t from=0) const;
      ssize_t binarySearch (const T &) const;
      bool contains (const T &) const;

      // -- Deletion
      void remove (ssize_t idx);
      void remove (ssize_t idx, ssize_t count);
      void removeAll ();

      // -- Size 
      inline ssize_t getSize () const { return size; }
      void resize (ssize_t newSize, bool keep=false);
      void extend ();
};

template<class T, class Idx>
SxArrayList<T,Idx>::SxArrayList ()
   : elements(NULL),
     size(0),
     sizeAlloc(0)
{
   // empty
}

template<class T, class Idx>
SxArrayList<T,Idx>::SxArrayList (ssize_t n)
{
   SX_CHECK (n >= 0, n);

   if (n == 0)  {
      sizeAlloc = 0;
      elements  = NULL;
   }  else  {
      sizeAlloc = static_cast<Idx>(n);
      elements  = new T [sizeAlloc];
   }
   size = static_cast<Idx>(n);
}

template<class T, class Idx>
SxArrayList<T,Idx>::SxArrayList (const SxArrayList<T,Idx> &in)
{
   SX_CHECK (in.size >= 0, in.size);

   sizeAlloc = in.sizeAlloc;
   if (sizeAlloc == 0)  {
      elements = NULL;
   }  else  {
      elements   = new T [sizeAlloc];
      T *destPtr = elements;
      T *srcPtr  = in.elements;
      for (Idx i=0; i < in.size; i++)  {
         *destPtr++ = *srcPtr++;
      }
   }
   size = in.size;
}

template<class T, class Idx>
SxArrayList<T,Idx>::SxArrayList (SxArrayList<T,Idx> &&in)
{
   SX_CHECK (in.size >= 0, in.size);

   sizeAlloc = in.sizeAlloc;
   if (sizeAlloc == 0)  {
      elements = NULL;
   }  else  {
      elements   = in.elements;
   }
   size = in.size;
   in.size = 0;
   in.sizeAlloc = 0;
   in.elements = NULL;
}

template<class T, class Idx>
SxArrayList<T,Idx>::~SxArrayList ()
{
   if (sizeAlloc)  delete [] elements;
   elements = NULL;
   size = 0;
   sizeAlloc = 0;
}

template<class T, class Idx>
SxArrayList<T,Idx> &SxArrayList<T,Idx>::operator= (
   const SxArrayList<T,Idx> &in)
{
   if (&in == this)  return *this;

   if (sizeAlloc)  delete [] elements;

   sizeAlloc = in.sizeAlloc;
   if (sizeAlloc == 0)  {
      elements = NULL;
   }  else  {
      elements   = new T [sizeAlloc];
      T *destPtr = elements;
      T *srcPtr  = in.elements;
      for (Idx i=0; i < in.size; i++)  {
         *destPtr++ = *srcPtr++;
      }
   }
   size = in.size;

   return *this;
}

template<class T, class Idx>
SxArrayList<T,Idx> &SxArrayList<T,Idx>::operator= (
   SxArrayList<T,Idx> &&in)
{
   if (&in == this)  return *this;

   if (sizeAlloc)  delete [] elements;

   sizeAlloc = in.sizeAlloc;
   if (sizeAlloc == 0)  {
      elements = NULL;
   }  else  {
      elements   = in.elements;
   }
   size = in.size;
   in.size = 0;
   in.sizeAlloc = 0;
   in.elements = NULL;

   return *this;
}

template<class T, class Idx>
T &SxArrayList<T,Idx>::operator() (ssize_t i)
{
   SX_CHECK (i >= 0 && i < static_cast<ssize_t>(size), i, size);
   return elements[i];
}


template<class T, class Idx>
const T &SxArrayList<T,Idx>::operator() (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < static_cast<ssize_t>(size), i, size);
   return elements[i];
}

template<class T, class Idx>
ssize_t SxArrayList<T,Idx>::append (const T &in)
{
   ssize_t idx = size;
   extend ();
   elements[idx] = in;
   return idx;
}

template<class T, class Idx>
ssize_t SxArrayList<T,Idx>::prepend (const T &in)
{
   insert (0, in);
   return 0;
}

template<class T, class Idx>
ssize_t SxArrayList<T,Idx>::binaryInsert (const T &in)
{
   ssize_t idx = size;
   extend ();
   // --- find sorted position
   while (idx > 0 && elements[idx - 1] > in)  {
      elements[idx] = elements[idx - 1];
      --idx;
   }
   elements[idx] = in;
   return idx;
}


template<class T, class Idx>
void SxArrayList<T,Idx>::insert (ssize_t newPos, const T &in)
{
   ssize_t idx = size;
   SX_CHECK (newPos <= idx, newPos, idx);
   extend ();
   for (ssize_t i=idx; i > newPos; --i)  elements[i] = elements[i - 1];
   elements[newPos] = in;
}

template<class T, class Idx>
ssize_t SxArrayList<T,Idx>::findPos (const T &in, ssize_t from) const
{
   SX_CHECK (from >= 0 && from <= size, from, size);

   ssize_t n = getSize ();
   for (ssize_t i=from; i < n; i++)  {
      if (elements[i] == in) return i;
   }
   return -1;
}

template<class T, class Idx>
bool SxArrayList<T,Idx>::contains (const T &in) const
{
   return (findPos (in) >= 0);
}

template<class T, class Idx>
ssize_t SxArrayList<T,Idx>::binarySearch (const T &in) const
{
   ssize_t l = 0;
   ssize_t r = size - 1;
   while (r >= l)  {
      const ssize_t m = (l + r) / 2;
      if (in == elements[m]) return m;
      if (in <  elements[m]) r = m - 1;
      else                   l = m + 1;
   }
   return -1;
}

template<class T, class Idx>
void SxArrayList<T,Idx>::remove (ssize_t idx)
{
   SX_CHECK (idx >= 0 && idx < static_cast<ssize_t>(size), idx, size);

   ssize_t n = size - 1;
   for (ssize_t i=idx; i < n; ++i)  elements[i] = elements[i + 1];
   elements[n] = T(); // destroy
   --size;
}

template<class T, class Idx>
void SxArrayList<T,Idx>::remove (ssize_t idx, ssize_t count)
{
   SX_CHECK (count > 0, count);
   SX_CHECK (idx >= 0 && idx + count <= static_cast<ssize_t>(size),
                  idx, size);

   ssize_t n = size;
   ssize_t nSwap = n - count - idx;
   for (ssize_t i=idx; i < idx + nSwap; ++i)  elements[i] = elements[i + count];
   for (ssize_t i=idx + nSwap; i < n; ++i)    elements[i] = T(); // destroy
   size -= static_cast<Idx>(count);
}

template<class T, class Idx>
void SxArrayList<T,Idx>::removeAll ()
{
   resize (0);
}

template<class T, class Idx>
void SxArrayList<T,Idx>::resize (ssize_t newSize, bool keep)
{
   SX_CHECK (newSize >= 0, newSize);
   if (newSize == static_cast<ssize_t>(sizeAlloc))  {
      size = sizeAlloc;
      return;
   }

   T *buffer = NULL;
   if (newSize > 0)  {
      buffer = new T [newSize];
      if (keep)  {
         ssize_t n = size;
         ssize_t min = (n < newSize) ? n : newSize;
         for (ssize_t i=0; i < min; ++i)   buffer[i] = elements[i];
      }
   }
   if (sizeAlloc)  delete [] elements;
   size      = static_cast<Idx>(newSize);
   sizeAlloc = static_cast<Idx>(newSize);
   elements  = buffer;
   buffer    = NULL;
}

template<class T, class Idx>
void SxArrayList<T,Idx>::extend ()
{
   if (size == sizeAlloc)  {
      if (sizeAlloc < 1) sizeAlloc = 1;
      else               sizeAlloc *= 2;
      T *buffer = new T [sizeAlloc];
      ssize_t n = size;
      for (ssize_t i=0; i < n; ++i)   buffer[i] = elements[i];
      if (elements)  delete [] elements;
      elements  = buffer;
      buffer    = NULL;
   }
   ++size;
}

template<class T, class Idx>
std::ostream &operator<< (std::ostream &s, const SxArrayList<T,Idx> &in)
{
   s << "[";
   ssize_t n = in.getSize ();
   for (ssize_t i=0; i < n; ++i)  {
      if (i != 0)  s << ", ";
      s << in(i);
   }
   s << "];";
   return s;
}

template<class T, class Idx>
size_t getNBytes (const SxArrayList<T,Idx> &in)
{
   size_t nBytes = sizeof (in);
   ssize_t n = in.sizeAlloc;
   for (ssize_t i=0; i < n; ++i)
      nBytes += ::getNBytes (in.elements[i]);

   return nBytes;
}

#endif /* _SX_ARRAY_LIST_H_ */
