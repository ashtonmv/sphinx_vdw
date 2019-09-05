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

#ifndef _SX_HASHED_ARRAY_H_
#define _SX_HASHED_ARRAY_H_

#include <SxString.h>
//#include <SxTypeDefs.h>

/** \brief Hashed Array Tree

   http://en.wikipedia.org/wiki/Hashed_array_tree

   Two-level Array
   Larson, Per-Åke, "Dynamic Hash Tables",
   Communications of the ACM, 31(4), 1988, p 446–457

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */
template <class T>
class SxHashedArray
{
   public:
      T **map;               //! Directory of pointers to blocks.
      ssize_t mapSize;       //! Number of blocks (second level arrays).
      ssize_t lastBlockIdx;  //! The last index in a block.
      ssize_t shift;         //! Ammount of shift to index the first level.
      ssize_t size;          //! Number of used elements.

      /** This constructor allocates a specified number of elements. The array
          is <em> not </em> being preset wit any value (such as 0.)! */
      SxHashedArray (ssize_t size=0, ssize_t blockSize=512);
      /** Copy constructor. */
      SxHashedArray (const SxHashedArray<T> &);
      SxHashedArray (SxHashedArray<T> &&);
      /** Destructor. */
      ~SxHashedArray ();

      /** Standard assignment operator. */
      SxHashedArray<T> &operator= (const SxHashedArray<T> &);
      SxHashedArray<T> &operator= (SxHashedArray<T> &&);

      /** Standard dereferencing operator for non-const objectes. */
      T &operator() (ssize_t);
      /** Standard dereferencing operator for const objectes. */
      const T &operator() (ssize_t) const;

      /** Initialize the (previously allocated) array with a certain value. */
      void set (const T &);

      // -- Insertion
      ssize_t append (const T &);
      ssize_t append (const T &, ssize_t count);

      // --Size 
      /** Returns the number of used elements in the array. */
      inline ssize_t getSize () const { return size; }

      /** Resize the array. New array keeps the old elements by default. */
      void resize (ssize_t newSize);

      // -- Deletion
      void removeLast ();
      void removeAll ();

   protected:
      void resizeMap (ssize_t mapSize);
      void copy (const SxHashedArray<T> &);
};

template<class T>
SxHashedArray<T>::SxHashedArray (ssize_t size_, ssize_t blockSize_)
{
   SX_CHECK (blockSize_ > 0, blockSize_);

   map = NULL;
   mapSize = 0;
   lastBlockIdx = 0;
   shift = 0;
   size = 0;

   // --- blockSize is the nearest power of 2
   while (lastBlockIdx < blockSize_)  {
      shift++;
      lastBlockIdx = 1 << shift;
   }
   --lastBlockIdx;
   resize (size_);
}

template<class T>
SxHashedArray<T>::SxHashedArray (const SxHashedArray<T> &in)
{
   map = NULL;
   mapSize = 0;
   lastBlockIdx = in.lastBlockIdx;
   shift = in.shift;
   size = 0;

   copy (in);
}

template<class T>
SxHashedArray<T>::SxHashedArray (SxHashedArray<T> &&in)
{
   map = in.map;
   mapSize = in.mapSize;
   lastBlockIdx = in.lastBlockIdx;
   shift = in.shift;
   size = in.size;

   in.map = NULL;
   in.mapSize = 0;
   in.lastBlockIdx = 0;
   in.shift = 0;
   in.size = 0;
}

template<class T>
SxHashedArray<T>::~SxHashedArray ()
{
   removeAll ();
}

template<class T>
SxHashedArray<T> &SxHashedArray<T>::operator= (const SxHashedArray<T> &in)
{
   if (&in == this)  return *this;
   copy (in);
   return *this;
}

template<class T>
SxHashedArray<T> &SxHashedArray<T>::operator= (SxHashedArray<T> &&in)
{
   if (&in == this)  return *this;

   resize (0);

   map = in.map;
   mapSize = in.mapSize;
   lastBlockIdx = in.lastBlockIdx;
   shift = in.shift;
   size = in.size;

   in.map = NULL;
   in.mapSize = 0;
   in.lastBlockIdx = 0;
   in.shift = 0;
   in.size = 0;

   return *this;
}

template<class T>
void SxHashedArray<T>::copy (const SxHashedArray<T> &in)
{
   resize (in.getSize ());
   SX_CHECK (size == in.size, size, in.size);
   if (size == 0) return;

   if (lastBlockIdx == in.lastBlockIdx)  {
      SX_CHECK (mapSize == in.mapSize, mapSize, in.mapSize);
      for (ssize_t i = 0; i < mapSize - 1; i++)  {
         T *dst = map[i];
         const T *src = in.map[i];
         for (ssize_t j=0; j <= lastBlockIdx; j++)
            dst[j] = src[j];
      }
      if (mapSize > 0)  {
         T *dst = map[mapSize - 1];
         const T *src = in.map[mapSize - 1];
         ssize_t lastIdx = (size - 1) & lastBlockIdx;
         for (ssize_t j=0; j <= lastIdx; j++)
            dst[j] = src[j];
      }
   }  else  {
      for (ssize_t i = 0; i < size; i++)  {
         operator()(i) = in(i);
      }
   }
}

template<class T>
T &SxHashedArray<T>::operator() (ssize_t i)
{
   SX_CHECK (i >= 0 && i < size, i, size);
   return *(map[i >> shift] + (i & lastBlockIdx));
}

template<class T>
const T &SxHashedArray<T>::operator() (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < size, i, size);
   return *(map[i >> shift] + (i & lastBlockIdx));
}

template<class T>
void SxHashedArray<T>::set (const T &in)
{
   SX_CHECK (size > 0);
   for (ssize_t i = 0; i < mapSize - 1; i++)  {
      T *block = map[i];
      SX_CHECK (block);
      for (ssize_t j=0; j <= lastBlockIdx; j++)
         block[j] = in;
   }
   if (mapSize > 0)  {
      T *block = map[mapSize - 1];
      SX_CHECK (block);
      ssize_t lastIdx = (size - 1) & lastBlockIdx;
      for (ssize_t j=0; j <= lastIdx; j++)
         block[j] = in;
   }
}

template<class T>
ssize_t SxHashedArray<T>::append (const T &in)
{
   resize (size + 1);
   operator()(size - 1) = in;
   return size;
}

template<class T>
ssize_t SxHashedArray<T>::append (const T &in, ssize_t count)
{
   SX_CHECK (count > 0, count);
   ssize_t i;
   ssize_t n = getSize ();

   resize (n + count);
   for (i=0; i < count; i++) {
      operator()(n + i) = in;
   }
   return size;
}

template<class T>
void SxHashedArray<T>::removeLast ()
{
   resize (size - 1);
}

template<class T>
void SxHashedArray<T>::removeAll ()
{
   resize (0);
}

template<class T>
void SxHashedArray<T>::resizeMap (ssize_t newSize)
{
   SX_CHECK (newSize >= 0, newSize);
   if (newSize == mapSize)  return;

   T **tmp =  NULL;
   if (newSize > 0)  {
      tmp =  new T*[(size_t)newSize];
      ssize_t min = (mapSize < newSize) ? mapSize : newSize;
      for (ssize_t i = 0; i < min; i++)  {
         tmp[i] = map[i];
      }
      for (ssize_t i = min; i < newSize; i++)  {
         tmp[i] = new T [(size_t)(lastBlockIdx + 1)];
      }
   }
   for (ssize_t i = newSize; i < mapSize; i++)  {
      SX_CHECK (map[i]);
      delete [] map[i];
   }
   if (mapSize > 0) delete [] map;
   map = tmp;
   mapSize = newSize;
}

template<class T>
void SxHashedArray<T>::resize (ssize_t newSize)
{
   SX_CHECK (newSize >= 0, newSize);
   if (newSize == size) return;

   if (newSize > 0)  {
      ssize_t oldMapSize = mapSize;
      ssize_t lastMapIdx = (newSize - 1) >> shift;
      resizeMap (lastMapIdx + 1);
      if (newSize < size)  {
         // --- call the destructors of removed elements in the last block
         ssize_t blockIdx = (newSize - 1) & lastBlockIdx;
         ssize_t oldBlockIdx = (size - 1) & lastBlockIdx;
         if (mapSize < oldMapSize)  {
            oldBlockIdx = lastBlockIdx;
         }
         if (blockIdx != oldBlockIdx)  {
            T *block = map[lastMapIdx];
            SX_CHECK (block);
            T empty;
            for (ssize_t i = blockIdx; i < oldBlockIdx; i++)  {
               // --- destroy map[lastMapIdx] at idx[i + 1];
               block[i +1] = empty;
            }
         }
      }
   }  else  {
      resizeMap (0);
   }
   size = newSize;
}

template<class T>
std::ostream& operator<< (std::ostream &s, const SxHashedArray<T> &in)
{
   s << "[";
   for (ssize_t i=0; i < in.getSize(); ++i)  {
      if (i != 0)  s << ", ";
      s << in(i);
   }
   s << "];";
   return s;
}

template<class T>
size_t getNBytes (const SxHashedArray<T> &in)
{
   size_t nBytes = sizeof (in);

   // --- elementwisely evaluation in case of something like
   //     SxHashedArray<SxList> or SxHashedArray<SxArray>
   for (ssize_t i = 0; i < in.mapSize; i++)  {
      nBytes += sizeof (T *);
      const T *block = in.map[i];
      SX_CHECK (block);
      for (ssize_t j=0; j <= in.lastBlockIdx; j++)
         nBytes += getNBytes (block[j]);
   }

   return nBytes;
}

#endif /* _SX_HASHED_ARRAY_H_ */
