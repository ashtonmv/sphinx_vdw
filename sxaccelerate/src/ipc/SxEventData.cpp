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

#include <SxEventData.h>

SxEventData::SxEventData ()
   :  iElem (0),
      nElem (0),
      offset (0),
      capacity (0),
      data (NULL),
      eventHash (0)
{
   SX_TRACE ();
   for (size_t idx = 0; idx < 6; ++idx)
      layout[idx] = 0;
}

SxEventData::~SxEventData ()
{
   if (data) free (data);
}

SxEventData::SxEventData (const SxEventData &in)
   :  iElem (in.iElem),
      nElem (in.nElem),
      offset (0),
      capacity (in.capacity),
      data (NULL),
      eventHash (in.eventHash)
{
   for (size_t idx = 0; idx < 6; ++idx)
      layout[idx] = in.layout[idx];

   //data = (uint8_t*)aligned_alloc (8, capacity);
   data = (uint8_t*) new int64_t[capacity/sizeof (int64_t) + 1];
   memcpy (data, in.data, capacity);
}

SxEventData &SxEventData::operator= (const SxEventData &in)
{
   if (this == &in) return *this;
   if (data) free (data);
   iElem = in.iElem;
   nElem = in.nElem;
   memcpy (layout, in.layout, 6 * sizeof (uint8_t));
   offset = 0;
   capacity = in.capacity;
   //data = (uint8_t*)aligned_alloc (8, capacity);
   data = (uint8_t*) new int64_t[capacity/sizeof (int64_t) + 1];
   memcpy (data, in.data, capacity);
   eventHash = in.eventHash;
   return *this;
}
