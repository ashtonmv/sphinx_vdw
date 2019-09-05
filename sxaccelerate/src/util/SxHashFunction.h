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

#ifndef _SX_HASH_FUNCTION_H_
#define _SX_HASH_FUNCTION_H_

//#include <SxString.h>
#include <SxArray.h>
#include <SxPtr.h>

/** \brief Hash Function

    \b SxHashFunction = SPHInX Hash Functions

\code
   cout << SxHashFunction::hash ("key") << endl;
\endcode

    \par Custom hash key: 
\code
#include <SxUniqueList.h>

class A
{
   public:
      int data;
   
      A () : data(0) {  }
      A (int data_) : data(data_) {  }
      
      bool operator== (const A &in) const  {
         return (in.data == data);
      }
};

class Hash
{
   public:
      static size_t hash (const A &in)
      {
         return (size_t)in.data;
      }
};

std::ostream &operator<< (std::ostream &s, const A &in)
{
   s << in.data;
   return s;
}

int main (int, char **)
{
   SxUniqueList<A, Hash> list;
   
   list << A(1) << A(1) << A(3) << A(2);
   
   std::cout << list << std::endl;
   
   // prints
   //  0: 1
   //  1: 3
   //  2: 2
   
   return 0;
}
\endcode

    \author Sixten Boeck, boeck@gemmantics.com
    \author Vaclav Bubnik, bubnik@gemmantics.com */

static uint32_t constexpr operator"" _SX (const char *s, size_t len);

class SxHashFunction
{
   public:
      
      static size_t hash (int key) { return (size_t)key; }
      static size_t hash (unsigned int key) { return (size_t)key; }
      static size_t hash (long key) { return (size_t)key; }
      static size_t hash (unsigned long key) { return (size_t)key; }
      static size_t hash (long long key) { return (size_t)key; }
      static size_t hash (unsigned long long key) { return (size_t)key; }
      static size_t hash (float);
      static size_t hash (double);
      static size_t hash (const SxArray<char> &);
      template<class T> static size_t hash (const SxPtr<T> &);
      template<class T> static size_t hash (const T *);
      
      static size_t hash (double *key, ssize_t n);
      static size_t hash (const SxArray<double> &);
      static size_t hash (const SxList<double> &);
      
      // --- combine
      static size_t combine (size_t prev, size_t h);

      // --- hash functions
      static unsigned int murmur (const void *key, int len, unsigned int seed);
      static uint32_t constexpr jenkinsHash (const char *key, size_t len);
};

inline size_t SxHashFunction::hash (float key)
{
   return murmur (reinterpret_cast<const void*>(&key), (int)sizeof(key), 12345);
}

inline size_t SxHashFunction::hash (double key)
{
   return murmur (reinterpret_cast<const void*>(&key), (int)sizeof(key), 12345);
}
      
inline size_t SxHashFunction::hash (const SxArray<char> &key)
{
   return murmur ((const void *)key.elements, (int)key.getSize (), 12345);
}

template<class T>
inline size_t SxHashFunction::hash (const SxPtr<T> &key)
{
   return murmur ((const void *)key.getPtr (), (int)sizeof(void*), 12345);
}

template<class T>
inline size_t SxHashFunction::hash (const T *key)
{
   return murmur ((const void *)key, (int)sizeof(void*), 12345);
}

inline size_t SxHashFunction::hash (double *key, ssize_t n)
{
   size_t h = 12345;
   for (ssize_t i = 0; i < n; ++i)  {
      h ^= hash (key[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
   }
   return h;
}

inline size_t SxHashFunction::hash (const SxArray<double> &key)
{
   return hash (key.elements, key.getSize ());
}

inline size_t SxHashFunction::hash (const SxList<double> &key)
{
   size_t h = 12345;
   SxList<double>::ConstIterator it;
   for (it = key.begin (); it != key.end (); ++it) {
      h ^= hash (*it) + 0x9e3779b9 + (h << 6) + (h >> 2);
   }
   return h;
}

inline size_t SxHashFunction::combine (size_t prev, size_t h)
{
   // --- golden ratio
   //     http://burtleburtle.net/bob/hash/doobs.html
   return (prev ^= h + 0x9e3779b9 + (prev << 6) + (prev >> 2));
}

// --- MurmurHashNeutral2, by Austin Appleby
//     LICENSE: public domain / MIT
//     Same as MurmurHash2, but endian- and alignment-neutral.
//     Half the speed though, alas.
inline unsigned int SxHashFunction::murmur (const void   *key,
                                            int          len,
                                            unsigned int seed)
{
   SX_CHECK (len >= 0);
   const unsigned int m = 0x5bd1e995;
   const int r = 24;

   unsigned int h = seed ^ static_cast<unsigned int>(len);

   const unsigned char * data = (const unsigned char *)key;

   while(len >= 4)
   {
      unsigned int k;

      k  = data[0];
      k |= static_cast<unsigned int>(data[1]) << 8;
      k |= static_cast<unsigned int>(data[2]) << 16;
      k |= static_cast<unsigned int>(data[3]) << 24;

      k *= m; 
      k ^= k >> r; 
      k *= m;

      h *= m;
      h ^= k;

      data += 4;
      len -= 4;
   }
   
   if (len == 3)  { h ^= static_cast<unsigned int>(data[2]) << 16; }
   if (len >= 2)  { h ^= static_cast<unsigned int>(data[1]) << 8;  }
   if (len >= 1)  { h ^= data[0]; h *= m; }

   h ^= h >> 13;
   h *= m;
   h ^= h >> 15;

   return h;
}

// ref.: https://en.wikipedia.org/wiki/Jenkins_hash_function
uint32_t constexpr SxHashFunction::jenkinsHash (const char *key,
                                                       size_t len)
{
   uint32_t i = 0;
   uint32_t hash = 0;
   while (i != len) {
      hash += (uint32_t)key[i++];
      hash += hash << 10;
      hash ^= hash >> 6;
   }
   hash += hash << 3;
   hash ^= hash >> 11;
   hash += hash << 15;
   return hash;
}

// --- user defined literal
static uint32_t constexpr operator"" _SX (const char *s, size_t len)
{
   return SxHashFunction::jenkinsHash (s, len);
}

#endif /* _SX_HASH_FUNCTION_H_ */
