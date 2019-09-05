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

#ifndef _SX_ALLOC_CACHE_H_
#define _SX_ALLOC_CACHE_H_

#include <SxUtil.h>
#include <SxError.h>

// small is already defined for Windows 8, 8.1 as char
#ifdef small
#   undef small
#endif

/** \brief This is a cache for memory allocations

    The purpose of this class is to temporarily retain memory allocations
    within a thread when memory of the same size must be sequentially allocated
    and deallocated, and the allocation/deallocation events are beyond the
    realm of temporaries and C++ move semantics. One important example is the
    use of temporary vectors in a openmp-parallelized loop. Obviously, each
    thread needs it's separate memory, which should be allocated once before
    the loop. On the other hand, the temporaries are created within the loop,
    so they cannot survive from one loop iteration to the next.

    This class solves the problem by providing a temporary, thread-local
    allocation cache that can be accessed globally from anywhere in the code.
    The cache is not thread-safe, in order to make it fast. To make it
    thread-safe, use thread-local caches via the macro
    \code
    SX_ALLOC_CACHE;
    \endcode

    You can access the thread's local cache (if set) via
    SxAllocation::threadCache, or via the following global helper functions
    which use the cache if available, and fall back to
    malloc/posix_memalign/free otherwise.

    To get an allocation from the cache, use
    \code
    void *mem  = SxAllocation::get (nBytes);
    void *mem2 = SxAllocation::getAligned (nBytes, 16);
    \endcode
    If no suitable allocation is available, new memory will be allocated

    To return allocated memory to the cache, use
    \code
    SxAllocation::retain (ptr, nBytes);
    \endcode
    If the maximum number of allocations in the cache is reached, the oldest
    retained allocation is freed (in the expectation that more recently
    returned allocations might be hot in the CPU cache for subsequent
    allocation requests).

    There is no restrictions where ptr has been allocated, but be aware of
    false sharing if small memory allocations are transfered between threads
    (i.e. separate memory allocations sitting on the same cache line can cause
    delays if accessed concurrently). To circumvent that in critical routines
    (frequent access), allocate memory aligned to a cache line boundary
    (=64 on Intel).

    \author C. Freysoldt freysoldt@mpie.de */
class SX_EXPORT_UTIL alignas(64) SxAllocCache
{
   public:
      /** Constructor
        @param maxSize_  max. allocation to cache
        @param nSlots    number of slots to use, <= 16
        */
      SxAllocCache (size_t maxSize_ = 2097152, unsigned char nSlots = 16);
      ~SxAllocCache ();

      // --- these functions are inline, to help branch prediction and
      //     compiler optimization
      /// Get a memory allocation
      inline void *get (size_t nBytes);
      /// Get an aligned memory allocation
      inline void *getAligned (size_t nBytes, size_t align);
      /// Return a memory allocation
      inline void retain (void *ptr, size_t nBytes);

   protected:
      /// Large allocations
      class {
         public:
            /// Large allocations: pointer
            void *ptr;
            /// Large allocations: sizes
            size_t nBytes;
      } large[16];
      /// Small allocations
      void* small[6][16];
      /// which slot to use for small sizes (1/2 byte for 1..64)
      char slotSmall[32];

   public:
      /// Maximum size to cache
      size_t maxSize;

   protected:
      /// Number of large slots to use
      const unsigned char nSlotsLarge;

      /// Large allocations: most recently returned allocation
      unsigned char largeHead;
      /// Large allocations: number of allocs
      unsigned char nLarge;

      /** \brief Number of available small allocations
        \note The last element used for uncached sizes (always 0).
      */
      unsigned char nsmall[6+1];
      /// Number of currently used small slots
      unsigned char smallUsed;

      /// Get a memory allocation
      void *getLarge (size_t nBytes, size_t align = 16);
      /// Return a memory allocation
      void retainLarge (void *ptr, size_t nBytes);
      /// Return a memory allocation
      void retainSmall (void *ptr, size_t nBytes);
      /// Auxiliary function: start a new small slot
      inline void newSmallSlot(int iSlot, size_t nBytes, void *ptr);
      /// Auxiliary function: slot id for a given small size
      inline int getSmallId (size_t nBytes)  {
         int id = nBytes & 63;
         return ((slotSmall[id>>1]) >> ((id & 1)*4)) & 0xf;
      }
};

/// This is the namespace class for thread-specific global caches
/// with limited lifetime as created by SxAllocation::TemporaryCache
class SxAllocation
{
   public:
      static thread_local SxAllocCache *threadCache;
      class TemporaryCache {
         private:
            /// The new temporary cache
            SxAllocCache newCache;
            /// Address of the old cache
            SxAllocCache *oldCache;
         public:
            TemporaryCache (size_t maxSize = 2097152, unsigned char nSlots = 16)
               : newCache (maxSize, nSlots)
            {
               oldCache = threadCache;
               threadCache = &newCache;
            }
            ~TemporaryCache ()  {
               // reactivate old cache
               threadCache = oldCache;
            }
      };

      /// Get a memory allocation
      static inline void *get (size_t nBytes)  {
         if (threadCache) return threadCache->get (nBytes);
         void *res = malloc (nBytes);
         if (!res) sxOutOfMemoryHandler ();
         return res;
      }
      /// Get an aligned memory allocation
      static inline void *getAligned (size_t nBytes, size_t align)  {
         if (threadCache) return threadCache->getAligned (nBytes, align);
         void *res = NULL;
#        ifdef WIN32
            res = _aligned_malloc (nBytes, align);
            if (res == NULL) sxOutOfMemoryHandler ();
#        else
            if (posix_memalign(&res, align, nBytes) != 0)  {
               sxOutOfMemoryHandler ();
            }
#        endif
         return res;
      }
      /// Return a memory allocation for caching
      static inline void retain (void *ptr, size_t nBytes)  {
         if (threadCache)
            threadCache->retain (ptr, nBytes);
         else
            free (ptr);
      }
};

// we purposefully risk name clashes: temporary allocation caches should be
// introduced at infrequent intervals, and definitely not multiple times
// in the same scope
#define SX_ALLOC_CACHE SxAllocation::TemporaryCache _myTempCache

/// Get a memory allocation
inline void *SxAllocCache::get (size_t nBytes)
{
   if (nBytes > 64) return getLarge(nBytes);
   int id = getSmallId (nBytes);
   if (nsmall[id])  {
      return small[id][int(--nsmall[id])];
   }
   // fallback:
   void *res = malloc(nBytes);
   if (!res) sxOutOfMemoryHandler ();
   return res;
}

/// Get an aligned memory allocation
inline void *SxAllocCache::getAligned (size_t nBytes, size_t align)
{
   if (nBytes > 64) return getLarge(nBytes, align);
   int id = getSmallId (nBytes);
   SX_CHECK(((align-1) & align) == 0, align); //must be power of 2
   size_t alignmask = align-1;
   int j = nsmall[id];
   while (j)  {
      j--;
      void *res = small[id][j];
      if ((reinterpret_cast<size_t>(res) & alignmask) == 0) {
         nsmall[id]--;
         small[id][j] = small[id][int(nsmall[id])];
         return res;
      }
   }
   // fallback:
   void *res = NULL;
#  ifdef WIN32
      res = _aligned_malloc (nBytes, align);
      if (res == NULL) sxOutOfMemoryHandler ();
#  else
      if (posix_memalign(&res, align, nBytes) != 0) sxOutOfMemoryHandler ();
#  endif

   return res;
}

/// Return a memory allocation
inline void SxAllocCache::retain (void *ptr, size_t nBytes)
{
   SX_CHECK (ptr);
   SX_CHECK (nBytes, nBytes);
   if (nBytes > 64) {
      if (nBytes > maxSize)  {
         free(ptr);
         return;
      }
      retainLarge(ptr, nBytes);
      return;
   }
   int id = getSmallId (nBytes);
   if (id < 6)  {
      if (nsmall[id] == 16)  {
         free (small[id][0]);
         for (int i = 0; i < 15; i++)
            small[id][i] = small[id][i+1];
         small[id][15] = ptr;
         return;
      }
      small[id][int(nsmall[id]++)] = ptr;
   } else {
      retainSmall(ptr, nBytes);
   }
}


#endif /* _SX_ALLOC_CACHE_H_ */
