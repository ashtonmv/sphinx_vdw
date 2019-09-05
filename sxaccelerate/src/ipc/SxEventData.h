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

#ifndef _SX_EVENT_DATA_H_
#define _SX_EVENT_DATA_H_

#include <SxIPC.h>
#include <SxHashFunction.h>
#include <SxArray.h>
#include <SxString.h>
#include <type_traits>

namespace sx {
   enum EventData : uint8_t
   {
      Int    = 1,
      String = 2
   };
}

class SxEvent;

/** \brief

   \b SxEventData = Event Data Storage

   A storage class for event data send over the SxEventBus.
   Potential users of the event bus will use SxEvent for
   retreiving event data from the listeners.

*/

class SX_EXPORT_IPC SxEventData
{
   public:

      friend class SxEvent;

      SxEventData ();
      SxEventData (const SxEventData &in_);
     ~SxEventData ();

      SxEventData &operator= (const SxEventData &in_);

      template<typename T1>
      SxEventData (const T1 &t1);

      template<typename T1, typename T2>
      SxEventData (const T1 &t1, const T2 &t2);

      template<typename T1, typename T2, typename T3>
      SxEventData (const T1 &t1, const T2 &t2, const T3 &t3);

      template<typename T1, typename T2, typename T3, typename T4>
      SxEventData (const T1 &t1, const T2 &t2,
                   const T3 &t3, const T4 &t4);

      template<typename T1, typename T2, typename T3,
               typename T4, typename T5>
      SxEventData (const T1 &t1, const T2 &t2, const T3 &t3,
                   const T4 &t4, const T5 &t5);

      template<typename T1, typename T2, typename T3,
               typename T4, typename T5, typename T6>
      SxEventData (const T1 &t1, const T2 &t2, const T3 &t3,
                   const T4 &t4, const T5 &t5, const T6 &t6);

      void setEventHash (uint32_t evHash) { eventHash = evHash; };

      uint64_t getLayout ()
      {
         uint64_t l = 0ULL;
         for (uint64_t i = 0; i < nElem; ++i) {
            uint8_t value = layout[i];
            l |= ((uint64_t)value) << (i * 8);
         }
         return l;
      }

      static uint64_t getLayout (const SxList<uint8_t> &ll)
      {
         SX_TRACE ();
         SX_CHECK (ll.getSize () < 7, ll.getSize ());
         uint64_t l = 0ULL;
         uint64_t i = 0ULL;
         for (uint8_t value : ll) {
            l |= (uint64_t)value << (i * 8);
            ++i;
         }
         return l;
      }

      static uint64_t getLayout (const SxArray<uint8_t> &a)
      {
         SX_TRACE ();
         SxList<uint8_t> l;
         for (uint8_t v : a) l << v;
         return getLayout(l);
      }

   protected:

      uint8_t iElem;      // current layout position
      uint8_t nElem;      // how many types
      uint8_t layout[6];  // event data layout (types)
      size_t offset;      // cursor position in bytes
      size_t capacity;    // maximum storage size in bytes
      uint8_t *data;      // byte array
      uint32_t eventHash; // hashed event name
      uint32_t unused;    // currently unused, manual padding

      inline void resize (size_t missing);

      inline void read (int64_t *i);
      inline void read (SxString *str);

      inline void write (int64_t i, size_t idx);
      inline void write (const SxString &str, size_t idx);
      inline void write (const char *str, size_t idx);

      inline ssize_t getSize (const int64_t &t);
      inline ssize_t getSize (const SxString &str);
      inline ssize_t getSize (const char *str);
};

template<typename T1>
SxEventData::SxEventData (const T1 &t1)
   :  iElem (0),
      nElem (1),
      offset (0),
      capacity (0),
      data (NULL),
      eventHash(0)
{
   SX_TRACE ();
   capacity = getSize (t1);
   //data = (uint8_t*)aligned_alloc (8, capacity);
   data = (uint8_t*) new int64_t[capacity/sizeof (int64_t) + 1];
   write (t1, 0);
   offset = 0;
}

template<typename T1, typename T2>
SxEventData::SxEventData (const T1 &t1, const T2 &t2)
   :  iElem (0),
      nElem (2),
      offset (0),
      capacity (0),
      data (NULL),
      eventHash (0)
{
   SX_TRACE ();
   capacity = getSize (t1) + getSize (t2);
   //data = (uint8_t*)aligned_alloc (8, capacity);
   data = (uint8_t*) new int64_t[capacity/sizeof (int64_t) + 1];
   write (t1, 0);
   write (t2, 1);
   offset = 0;
}

template<typename T1, typename T2, typename T3>
SxEventData::SxEventData (const T1 &t1, const T2 &t2, const T3 &t3)
   :  iElem (0),
      nElem (3),
      offset (0),
      capacity (0),
      data (NULL),
      eventHash (0)
{
   SX_TRACE ();
   capacity = getSize (t1) + getSize (t2) + getSize (t3);
   //data = (uint8_t*)aligned_alloc (8, capacity);
   data = (uint8_t*) new int64_t[capacity/sizeof (int64_t) + 1];
   write (t1, 0);
   write (t2, 1);
   write (t3, 2);
   offset = 0;
}

template<typename T1, typename T2, typename T3, typename T4>
SxEventData::SxEventData (const T1 &t1,
                          const T2 &t2,
                          const T3 &t3,
                          const T4 &t4)
   :  iElem (0),
      nElem (4),
      offset (0),
      capacity (0),
      data (NULL),
      eventHash (0)
{
   SX_TRACE ();
   capacity = getSize (t1) + getSize (t2) + getSize (t3) + getSize (t4);

   //data = (uint8_t*)aligned_alloc (8, capacity);
   data = (uint8_t*) new int64_t[capacity/sizeof (int64_t) + 1];
   write (t1, 0);
   write (t2, 1);
   write (t3, 2);
   write (t4, 3);
   offset = 0;
}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
SxEventData::SxEventData (const T1 &t1,
                          const T2 &t2,
                          const T3 &t3,
                          const T4 &t4,
                          const T5 &t5)
   :  iElem (0),
      nElem (5),
      offset (0),
      capacity (0),
      data (NULL),
      eventHash (0)
{
   SX_TRACE ();
   capacity = getSize (t1) + getSize (t2) + getSize (t3)
            + getSize (t4) + getSize (t5);
   //data = (uint8_t*)aligned_alloc (8, capacity);
   data = (uint8_t*) new int64_t[capacity/sizeof (int64_t) + 1];
   write (t1, 0);
   write (t2, 1);
   write (t3, 2);
   write (t4, 3);
   write (t5, 4);
   offset = 0;
}

template<typename T1, typename T2, typename T3,
         typename T4, typename T5, typename T6>
SxEventData::SxEventData (const T1 &t1,
                          const T2 &t2,
                          const T3 &t3,
                          const T4 &t4,
                          const T5 &t5,
                          const T6 &t6)
   :  iElem (0),
      nElem (6),
      offset (0),
      capacity (0),
      data (NULL),
      eventHash (0)
{
   SX_TRACE ();
   capacity = getSize (t1) + getSize (t2) + getSize (t3)
            + getSize (t4) + getSize (t5) + getSize (t6);

   //data = (uint8_t*)aligned_alloc (8, capacity);
   data = (uint8_t*) new int64_t[capacity/sizeof (int64_t) + 1];
   write (t1, 0);
   write (t2, 1);
   write (t3, 2);
   write (t4, 3);
   write (t5, 4);
   write (t6, 5);
   offset = 0;
}

inline ssize_t SxEventData::getSize (const int64_t &t)
{
   SX_UNUSED (t);
   return sizeof (int64_t);
}

inline ssize_t SxEventData::getSize (const SxString &str)
{
   SX_TRACE ();
   return (ssize_t)(str.getSize () + 1);
}

inline ssize_t SxEventData::getSize (const char *str)
{
   SX_TRACE ();
   return (ssize_t)(strlen (str) + 1);
}

void SxEventData::read (int64_t *i)
{
   SX_TRACE ();
   SX_CHECK (offset + sizeof (int64_t) <= capacity, offset, capacity);
   *i = *((int64_t*)(data + offset));
   offset += sizeof (int64_t);
}

void SxEventData::read (SxString *str)
{
   SX_TRACE ();
   *str = SxString ((const char *)&data[offset]);
   offset += str->getSize () + 1;
}

void SxEventData::write (int64_t i, size_t layoutIdx)
{
   SX_TRACE ();
   SX_CHECK (layoutIdx < nElem, layoutIdx, nElem);
   layout[layoutIdx] = sx::Int;
   memcpy (data + offset, &i, sizeof (int64_t));
   offset += sizeof (int64_t);
}

void SxEventData::write (const SxString &str, size_t layoutIdx)
{
   SX_TRACE ();
   SX_CHECK (layoutIdx < nElem, layoutIdx, nElem);
   layout[layoutIdx] = sx::String;
   for (const char &c : str) {
      data[offset++] = c;
   }
}

void SxEventData::write (const char *str, size_t layoutIdx)
{
   SX_TRACE ();
   SX_CHECK (layoutIdx < nElem, layoutIdx, nElem);
   layout[layoutIdx] = sx::String;
   // strcpy, and increasing offset
   while ((data[offset++] = *str++)) { }
}

#endif /* _SX_EVENT_DATA_H_ */
