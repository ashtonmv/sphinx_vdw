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

SxEventBus::SxEventBus ()
   : listenerMtx(SxMutex ())
{
   SX_TRACE ();
}

SxEventBus::~SxEventBus ()
{
   SX_TRACE ();
}

SxEventBus& SxEventBus::getGlobalObj ()
{
   SX_TRACE ();
   static SxEventBus globalBus;
   return globalBus;
}

void SxEventBus::terminate ()
{
   SX_TRACE ();
   static SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx)  {
      SxEventListener::terminate ();
      eventBus.listener.removeAll ();
   }
}

void SxEventBus::unsubscribeAll (SxEventListener *listenerPtr)
{
   SX_TRACE ();
   // --- this is costly, it will stalles the bus
   static SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx) {
      ListenerLookup &listener = getGlobalObj ().listener;
      for (auto it = listener.begin (); it != listener.end (); ++it)  {
         SxArray<SxEventListener*> &store = it.getValue ();
         ssize_t pos = store.findPos (listenerPtr);
         if (pos == -1) continue;
         store(pos) = NULL;
      }
      listenerPtr->setUnsubscribed ();
   }
}

std::ostream &operator<< (std::ostream &stream, SxEventBus &bus)
{
   SX_MUTEX (bus.listenerMtx) {

      const SxMap<uint32_t, SxArray<SxEventListener*> > &listener
         = bus.listener;
      stream << "Subscribed listener:" << endl;
      for (auto lstIt = listener.begin (); lstIt != listener.end (); ++lstIt)  {
         uint32_t evHash = lstIt.getKey ();
         for (SxEventListener *lst : lstIt.getValue ()) {
            if (lst == NULL) continue;
            const SxString &evName = lst->getEventName (evHash);
            const SxString &sxTag  = lst->getSubscribeTag (evHash);
            const SxString &lstName = lst->getName ();

            stream << "0x" << std::hex << evHash  << " | "
                   << evName  << "->"
                   << lstName << "@"
                   << sxTag   << endl;
         }
      }

      stream << "Orphaned events:\n";
      for (auto it = bus.unreceived.begin (); it != bus.unreceived.end (); ++it) {
         stream << *it << endl;
      }
      stream << std::dec;
   }
   return stream;
}

void SxEventBus::subscribe (uint32_t event,
                            SxEventListener *listenerPtr,
                            const SxArray<uint8_t> &layout)
{
   SX_TRACE ();
   static SxEventBus &eventBus = getGlobalObj ();

   SxMutexBlock (&eventBus.listenerMtx);
   const ssize_t RESIZE_COUNT = 100;
   static ListenerLookup &listener = eventBus.listener;

   // --- the listener is not associated with this event
   SX_CHECK (listenerPtr->hasEventAssociation(event), event,
             listenerPtr->getName (), listenerPtr->getSubscribeTag (event));

   // --- dublicate subscribe
   SX_CHECK (!listener(event).contains (listenerPtr), event,
             listenerPtr->getName (), listenerPtr->getSubscribeTag (event));

   if (listener(event).getSize () == 0) {
      listener(event).resize (RESIZE_COUNT);
      listener(event).set (NULL);
   }

   ssize_t nextNull = listener(event).findPos (NULL);
   if (nextNull == -1)  {
      ssize_t oldSize = listener(event).getSize ();
      listener(event).resize (RESIZE_COUNT + oldSize, true);

      for (ssize_t idx = 0; idx < RESIZE_COUNT - 1; ++idx)
         listener(event)(idx + oldSize) = NULL;

      listener(event)(oldSize) = listenerPtr;
   } else {
      listener(event)(nextNull) = listenerPtr;
   }

   // --- just in debug?!
   uint64_t layoutHash = SxEventData::getLayout (layout);
   if (eventBus.layoutLookup.hasKey (event)) {
      // --- event subscribed with different layout, BUG!
      SX_CHECK (eventBus.layoutLookup(event) == layoutHash);
   } else {
      eventBus.layoutLookup(event) = layoutHash;
   }
}

void SxEventBus::unsubscribe (uint32_t event, SxEventListener *listenerPtr)
{
   SX_TRACE ();
   SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx) {
      static ListenerLookup &listener = eventBus.listener;

      // --- the listener is not associated with this event
      SX_CHECK (listenerPtr->hasEventAssociation(event), event,
                listenerPtr->getName (), listenerPtr->getSubscribeTag (event));

      ssize_t pos = listener(event).findPos (listenerPtr);
      SX_CHECK (pos >= 0, event, listenerPtr->getName ()); // listener not subscribed
      listener(event)(pos) = NULL;
   }
}
