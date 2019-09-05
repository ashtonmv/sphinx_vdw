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

#include <SxEventBus.h>

std::atomic<bool> SxEventListener::terminateSignal(false);

void SxEventListener::terminate (bool kill)
{
   SX_TRACE ();
   terminateSignal.store (kill, std::memory_order_relaxed);
}

SxEventListener::SxEventListener (
   const std::initializer_list<SxPair<const char*, uint32_t> > &initList)
   :  unsubscribed (true),
      totalSize(0)
{
   SX_TRACE ();
   for (const SxPair<const char*, uint32_t> &p : initList) {
      // --- spelling mistake, check the event names
      SX_CHECK (SxHashFunction::jenkinsHash (p.key, strlen (p.key)) == p.value,
                p.key);
      eventNames(p.value) = SxString (p.key);
      storage = SxList<SxEvent> ();
      waits(p.value) = new std::atomic<size_t>(0);
      SX_CHECK (waits(p.value));
   }
}

SxEventListener::SxEventListener ()
   :  unsubscribed(true),
      totalSize(0)
{
   // empty
}

SxEventListener::SxEventListener(const SxArray<SxString> &evNames,
                                 const SxArray<uint32_t> &evHashes)
   :  unsubscribed(true),
      totalSize(0)
{
   SX_TRACE ();
   SX_CHECK (evNames.getSize () == evHashes.getSize (),
             evNames.getSize (), evHashes.getSize ());

   for (ssize_t idx = 0; idx < evNames.getSize (); ++idx) {
      // --- spelling mistake, check the event names
      SX_CHECK (SxHashFunction::jenkinsHash (evNames(idx).elements,
                                             evNames(idx).getSize ())
                 == evHashes(idx), evNames(idx));

      uint32_t evHash = evHashes(idx);
      eventNames(evHash) = evNames(idx);
      storage = SxList<SxEvent> ();
      waits(evHash) = new std::atomic<size_t>(0);
   }
}

SxEventListener::~SxEventListener ()
{
   SX_TRACE ();
   SX_MUTEX (mtx) {
      for (auto it = waits.begin (); it != waits.end (); ++it) {
         delete it.getValue ();
      }
      SX_CHECK (unsubscribed); //-> SxEventBus::unsubscribeAll (&listener);
   }
}

void SxEventListener::setEvents (const SxArray<SxString> &evNames,
                                 const SxArray<uint32_t> &evHashes)
{
   SX_TRACE ();
   SX_CHECK (evNames.getSize () == evHashes.getSize (),
             evNames.getSize (), evHashes.getSize ());

   SX_MUTEX (mtx) {
      for (ssize_t idx = 0; idx < evNames.getSize (); ++idx) {

         uint32_t evHash = evHashes(idx);

         // --- spelling mistake, check the event names
         SX_CHECK (SxHashFunction::jenkinsHash (evNames(idx).elements,
                                                evNames(idx).getSize ())
                   == evHashes(idx), evNames(idx));

         // --- is the event name already there?
         SX_CHECK (!eventNames(evHash), evNames(idx), evHash);
         eventNames(evHash) = evNames(idx);
         waits(evHash) = new std::atomic<size_t>(0);
      }
   }
}

std::atomic<size_t> *SxEventListener::getWaitSize (uint32_t event) const
{
   SX_TRACE ();
   SxMutexBlock lock (&waitMtx);
   SX_CHECK (waits.hasKey (event), event);
   std::atomic<size_t>* waiter = waits(event);
   return waiter;
}

bool SxEventListener::hasData (uint32_t event) const
{
   return getWaitSize(event)->load (std::memory_order_relaxed) > 0 ? true : false;
}

SxEvent SxEventListener::get ()
{
   SX_TRACE ();
   SxEvent evData;
   SX_MUTEX (mtx) {
      evData = storage.first ();
      storage.removeFirst ();
   }
   std::atomic<size_t> *size = getWaitSize (evData.getHash ());
   size->fetch_sub (1, std::memory_order_relaxed);
   totalSize.fetch_sub (1, std::memory_order_relaxed);
   return evData;
}

SxList<SxEvent> SxEventListener::getAll ()
{
   SX_TRACE ();
   SxList<SxEvent> evData;
   SX_MUTEX (mtx) {
      evData = std::move (storage);
      storage = SxList<SxEvent>();
      SX_MUTEX (waitMtx) {
         for (auto it = waits.begin (); it != waits.end (); ++it) {
            it.getValue ()->store (0, std::memory_order_relaxed);
         }
      }
      totalSize.store (0, std::memory_order_relaxed);
   }
   return evData;
}

SxEvent SxEventListener::getLast ()
{
   SX_TRACE ();
   SxEvent evData;
   SX_MUTEX (mtx) {
      evData = storage.last ();
      storage.removeAll ();
      SX_MUTEX (waitMtx) {
         for (auto it = waits.begin (); it != waits.end (); ++it) {
            it.getValue ()->store (0, std::memory_order_relaxed);
         }
      }
      totalSize.store (0, std::memory_order_release);
   }
   return evData;
}

void SxEventListener::waitForEvents_Or (
   SxArray<uint32_t> events, bool *terminateCalled, unsigned int sleepMS)
{
   SX_TRACE ();

   static std::atomic<bool> &killSig = SxEventListener::terminateSignal;

   SX_CHECK (terminateCalled, events.getSize () > 0, events.getSize ());

   SxArray<std::atomic<size_t> *> waitConds;
   SX_MUTEX (mtx) {

#     ifndef NDEBUG
      for (uint32_t ev : events)
         SX_CHECK (eventNames.hasKey (ev), ev);
#     endif

      waitConds = SxArray<std::atomic<size_t> *>(events.getSize ());
      waitConds.set (NULL);
      for (ssize_t idx = 0; idx < events.getSize (); ++idx) {
         waitConds(idx) = getWaitSize (events(idx));
      }

   }

   size_t n = 0;
   for (const std::atomic<size_t> *w : waitConds)
      n += w->load (std::memory_order_relaxed);

   while (n <= 0 && !killSig.load (std::memory_order_relaxed)) {
      SxTime::msleep (sleepMS);
      for (const std::atomic<size_t> *w : waitConds)
         n += w->load (std::memory_order_relaxed);
   }
   *terminateCalled = killSig.load (std::memory_order_relaxed);
}

void SxEventListener::waitForEvents_And (
   SxArray<uint32_t> events, bool *terminateCalled, unsigned int sleepMS)
{
   SX_TRACE ();

   static std::atomic<bool> &killSig = SxEventListener::terminateSignal;

   SX_CHECK (terminateCalled, events.getSize () > 0, events.getSize ());

   SxArray<std::atomic<size_t> *> waitConds;
   SX_MUTEX (mtx) {

#     ifndef NDEBUG
         for (uint32_t ev : events)
            SX_CHECK (eventNames.hasKey (ev), ev);
#     endif

      waitConds = SxArray<std::atomic<size_t> *>(events.getSize ());
      waitConds.set (NULL);
      for (ssize_t idx = 0; idx < events.getSize (); ++idx) {
         waitConds(idx) = getWaitSize (events(idx));
      }

   }

   size_t n = 0;
   for (const std::atomic<size_t> *w : waitConds)
      n += w->load (std::memory_order_relaxed) > 0 ? 1 : 0;
   size_t waitN = events.getSize ();
   while (!killSig.load (std::memory_order_relaxed)) {
      SxTime::msleep (sleepMS);
      for (const std::atomic<size_t> *w : waitConds)
         n += w->load (std::memory_order_relaxed) > 0 ? 1 : 0;
      if (n < waitN) n = 0;
      else           break;
   }
   *terminateCalled = killSig.load (std::memory_order_relaxed);
}

bool SxEventListener::hasEventAssociation (uint32_t event) const
{
   SX_TRACE ();
   bool hasEvent = false;
   SX_MUTEX (mtx) {
      hasEvent = eventNames.hasKey (event);
   }
   return hasEvent;
}

SxString SxEventListener::getName () const
{
   SX_TRACE ();
   SxString res;
   SX_MUTEX (mtx) {
      res = name;
   }
   return res;
}

SxString SxEventListener::getEventName (uint32_t evHash) const
{
   SX_TRACE ();
   SxString res;
   SX_MUTEX (mtx) {
      res = eventNames(evHash);
   }
   return res;
}

void SxEventListener::setSubscribeTag (uint32_t event, const SxString &sxTag)
{
   SX_TRACE ();
   SX_MUTEX (mtx) {
      unsubscribed = false;
      subscrTags(event) = sxTag;
   }
}

SxString SxEventListener::getSubscribeTag (uint32_t evHash) const
{
   SX_TRACE ();
   SxString res;
   SX_MUTEX (mtx) {
      res = subscrTags(evHash);
   }
   return res;
}
