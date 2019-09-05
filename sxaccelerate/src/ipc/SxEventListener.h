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

#ifndef _SX_EVENT_LISTENER_H
#define _SX_EVENT_LISTENER_H

#include <SxIPC.h>
#include <SxTime.h>
#include <SxMutex.h>
#include <SxMap.h>
#include <SxArray.h>
#include <SxString.h>
#include <SxEvent.h>
#include <atomic>

/** \brief

   \b SxEventData = Event Listener/Receiver

   A event listener can subscribe to any event of interest.
   The thread in which 'waitForEvents' (and similar) is used
   is put to sleep (with a set sleep time or) until an event is received
   or the event bus is terminated. The events are stored in the order they were
   received.
*/

class SX_EXPORT_IPC SxEventListener
{
   public:

      SxEventListener ();
      SxEventListener (const SxArray<SxString> &evNames,
                       const SxArray<uint32_t> &evHashes);
      SxEventListener (
         const std::initializer_list<SxPair<const char*, uint32_t> > &initList
      );

     ~SxEventListener ();

      void setName (const SxString &name_) { SX_MUTEX (mtx) { name = name_; } }

      SxString getName () const;

      SxString getEventName (uint32_t evHash) const;

      void setEvents (const SxArray<SxString> &evNames,
                      const SxArray<uint32_t> &evHashes);

      void setSubscribeTag (uint32_t event, const SxString &sxTag);

      SxString getSubscribeTag (uint32_t evHash) const;

      template<uint32_t event>
      void set (const SxEvent &args);

      size_t  getSize () const { return totalSize.load (std::memory_order_relaxed); }
      SxEvent get ();
      SxEvent getLast ();
      SxList<SxEvent> getAll ();

      bool hasData (uint32_t event) const;

      template<uint32_t event>
      void waitForEvents (
         bool *terminateCalled, size_t nEvents = 1, unsigned int sleepMS = 10
      );

      // --- sleeps until atleast one event has data
      void waitForEvents_Or  (
         SxArray<uint32_t> events, bool *terminateCalled, unsigned int sleepMS = 10
      );

      // --- sleeps until each event has data
      void waitForEvents_And (
         SxArray<uint32_t> events, bool *terminateCalled, unsigned int sleepMS = 10
      );

      bool hasEventAssociation (uint32_t event) const;

      void setUnsubscribed () { SX_MUTEX (mtx) { unsubscribed = true; } }

      // --- can be used to set the terminateSignal value, kill by default
      static void terminate (bool kill = true);

      // --- copy constructors shall not be usable here and in inherited classes
      //     Since the event bus stores pointers to objects.
      SxEventListener (const SxEventListener &in_) = delete;
      SxEventListener &operator= (const SxEventListener &in_) = delete;

   protected:

      mutable SxMutex mtx;
      SxString name;
      bool unsubscribed;
      SxMap<uint32_t, SxString> eventNames;
      SxMap<uint32_t, SxString> subscrTags;
      std::atomic<size_t> totalSize;
      SxList<SxEvent> storage;

      mutable SxMutex waitMtx;
      SxMap<uint32_t, std::atomic<size_t>* > waits;

      static std::atomic<bool> terminateSignal;
      std::atomic<size_t> *getWaitSize (uint32_t event) const;
};

template<uint32_t event>
void SxEventListener::set (const SxEvent &args)
{
   SX_TRACE ();
   std::atomic<size_t> *size = getWaitSize (event);
   SX_MUTEX (mtx)  {
      SX_CHECK (eventNames.hasKey (event), event);
      storage.append (args);
      size->fetch_add (1, std::memory_order_relaxed);
      totalSize.fetch_add (1, std::memory_order_relaxed);
   }
}

template<uint32_t event>
void SxEventListener::waitForEvents (
   bool *terminateCalled, size_t nEvents, unsigned int sleepMS)
{
   SX_TRACE ();
   SX_CHECK (terminateCalled, event);
   SX_CHECK (eventNames.hasKey (event), event);

   static std::atomic<bool> &killSig = SxEventListener::terminateSignal;
   std::atomic<size_t> *size = getWaitSize (event);
   size_t n = size->load (std::memory_order_relaxed);
   while (n < nEvents && !killSig.load (std::memory_order_relaxed)) {
      SxTime::msleep (sleepMS);
      n = size->load (std::memory_order_relaxed);
   }
   *terminateCalled = killSig.load (std::memory_order_relaxed);
}

#endif // _SX_EVENT_LISTENER_H
