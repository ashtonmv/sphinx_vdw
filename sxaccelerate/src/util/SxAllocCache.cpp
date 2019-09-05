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

#include <SxAllocCache.h>

#include <atomic>
#include <SxString.h>

// --- quick hack: class to create one output file per thread for logging
//class SxMyLog {
//   public:
//      FILE *fp;
//      SxMyLog () {
//         static std::atomic<int> g_id(0);
//         int id = g_id++;
//         fp = fopen (("t" + SxString(id) + ".log").ascii (), "w");
//      }
//      ~SxMyLog () { fclose (fp); }
//};
//namespace { SxMyLog thread_local tLog; }

thread_local SxAllocCache *SxAllocation::threadCache = NULL;

SxAllocCache::SxAllocCache (size_t maxSize_, unsigned char nSlots)
: maxSize (maxSize_), nSlotsLarge(nSlots),
  largeHead(0), nLarge(0), smallUsed(0)
{
   for (int i = 0; i < 16; i++) {
      large[i].nBytes = 0;
      large[i].ptr = NULL;
   }
   for (int i = 0; i < 6+1; i++) nsmall[i] = 0;
   for (int i = 0; i < 32; i++) slotSmall[i] = 0x66;
//   fprintf(tLog.fp, "Temp cache @%p\n", (void*)this); fflush(tLog.fp);
}

void *SxAllocCache::getLarge(size_t nBytes, size_t align)
{
   size_t alignMask = align - 1;
   SX_CHECK((alignMask & align) == 0, align); // must be power of 2
   // search large slots
   if (nBytes <= maxSize)  {
      for (int i = 0, ix = largeHead; i < nLarge; i++, ix++)  {
         ix &= 0xf; // cyclic
         if (   large[ix].nBytes == nBytes
             && (reinterpret_cast<size_t>(large[ix].ptr) & alignMask) == 0)
         {
            void *res = large[ix].ptr;
            if (i * 2 < nLarge)  {
               // move lower slots to close gap at ix
               int jx = (ix - 1)&0xf;
               while (i--)  {
                  large[ix] = large[jx];
                  ix=jx;
                  jx--; jx &=0xf;
               }
               large[ix].nBytes=0;
#              ifndef NDEBUG
                  large[ix].ptr = NULL;
#              endif
               largeHead++;
               largeHead &=0xf;
            } else {
               // move upper slots to close gap at ix
               int jx = (ix + 1)&0xf;
               while (++i<nLarge)  {
                  large[ix] = large[jx];
                  large[ix] = large[jx];
                  ix=jx;
                  jx++; jx &=0xf;
               }
               large[ix].nBytes = 0;
#              ifndef NDEBUG
                  large[ix].ptr = NULL;
#              endif
            }
            nLarge--;
//            fprintf(tLog.fp, "%p get %ld -> %p\n", (void*)this, nBytes, res); fflush(tLog.fp);
            return res;
         }
      }
   }
   // do a real allocation
   void *res = NULL;
#  ifdef WIN32
      res = _aligned_malloc (nBytes, align);
      if (res == NULL) sxOutOfMemoryHandler ();
#  else
      if (posix_memalign(&res, align, nBytes) != 0) sxOutOfMemoryHandler ();
#  endif

   return res;
}

void SxAllocCache::retainLarge(void *ptr, size_t nBytes)
{
   // search large slots
   if (nLarge == nSlotsLarge)  {
      int pos = (largeHead + nLarge-1) & 0xf;
      free (large[pos].ptr);
      large[pos].nBytes = 0;
      nLarge--;
   }
   largeHead--;
   largeHead &= 0xf;
   large[int(largeHead)].ptr    = ptr;
   large[int(largeHead)].nBytes = nBytes;
   nLarge++;
   //fprintf(tLog.fp, "%p ret %ld -> %p\n", (void*)this, nBytes, ptr); fflush(tLog.fp);
}

inline void SxAllocCache::newSmallSlot(int iSlot, size_t nBytes, void *ptr)
{
   int id = nBytes & 63;
   char &byte = slotSmall[id>>1];
   byte = char(byte + ((iSlot-6)<< (4*(id&1))));
   nsmall[iSlot] = 1;
   small[iSlot][0] = ptr;
}

void SxAllocCache::retainSmall (void *ptr, size_t nBytes)
{
   if (smallUsed == 6)  {
      for (int i = 0; i < 6; i++)  {
         if (nsmall[i] == 0)  {
            // find which size this slot was used for
            // and reset the slot for that size to magic id=7
            for (int k = 0; k < 32; k++)  {
               if ((slotSmall[k] & 0x0f) == i)  {
                  slotSmall[k] = char((slotSmall[k] & 0xf0) | 0x06);
                  break;
               }
               if ((slotSmall[k] & 0xf0) == i << 4)  {
                  slotSmall[k] = char((slotSmall[k] & 0x0f) | 0x60);
                  break;
               }
            }
            // use this free slot for new size nBytes
            newSmallSlot(i, nBytes, ptr);
            return;
         }
      }
#     ifndef NDEBUG
      //std::cout << "allocation cache small size slot overflow" << std::endl;
#     endif
      //free (ptr); return; // give up
      // --- empty slot 0 and move down all others
      for (int iPtr = 0; iPtr < nsmall[0]; iPtr++)
         free(small[0][iPtr]);
      for (int iSlot = 1; iSlot < 6; iSlot++)  {
         nsmall[iSlot-1] = nsmall[iSlot];
         for (int iPtr = 0; iPtr < nsmall[iSlot]; iPtr++)
            small[iSlot-1][iPtr] = small[iSlot][iPtr];
      }
      // update slot ids
      for (int k = 0; k < 32; k++)  {
         if ((slotSmall[k] & 0x0f) == 0)  {
            slotSmall[k] = char((slotSmall[k] & 0xf0) | 0x06);
         } else if ((slotSmall[k] & 0x0f) < 0x06) {
            slotSmall[k] = char(slotSmall[k] - 1);
         }
         if ((slotSmall[k] & 0xf0) == 0)  {
            slotSmall[k] = char((slotSmall[k] & 0x0f) | 0x60);
         } else if ((slotSmall[k] & 0xf0) < 0x60)  {
            slotSmall[k] = char(slotSmall[k] - 0x10);
         }
      }
      newSmallSlot(5, nBytes, ptr);
      return;
   }
   // use next free slot for new size=nBytes
   newSmallSlot(smallUsed++, nBytes, ptr);
}

SxAllocCache::~SxAllocCache ()
{
   for (int i = 0; i < 16; i++)  {
      if (large[i].nBytes) free(large[i].ptr);
#     ifndef NDEBUG
         large[i].ptr = NULL;
#     endif
   }
   for (int i = 0; i < smallUsed; i++) {
      for (int j = nsmall[i]; j; )  {
         free (small[i][--j]);
#        ifndef NDEBUG
            small[i][j] = NULL;
#        endif
      }
   }
}
