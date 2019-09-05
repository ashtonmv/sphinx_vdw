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

#ifndef _SX_EVENT_H_
#define _SX_EVENT_H_

#include <SxIPC.h>
#include <SxEventData.h>

class SX_EXPORT_IPC SxEvent
{
   public:

      SxEvent ();
      SxEvent (const SxPtr<SxEventData> &data_);
      SxEvent (const SxEvent &in_);
      SxEvent &operator= (const SxEvent &in_);
     ~SxEvent ();

      bool operator== (uint32_t evHash)
      {
         SX_TRACE ();
         SX_CHECK (data.getPtr ());
         return (data.getPtr ()->eventHash == evHash);
      }

      uint32_t getHash ()
      {
         SX_TRACE ();
         SX_CHECK (data.getPtr ());
         return data.getPtr ()->eventHash;
      }

      uint64_t getLayout ()
      {
         SX_TRACE ();
         SX_CHECK (data.getPtr ());
         return data.getPtr ()->getLayout ();
      }

      SxString shiftString ();
      int64_t  shiftInt ();
      SxEvent &operator>> (SxString &str);
      SxEvent &operator>> (int64_t &i);

   protected:
      SxPtr<SxEventData> data;
};


#endif /* _SX_EVENT_H_ */
