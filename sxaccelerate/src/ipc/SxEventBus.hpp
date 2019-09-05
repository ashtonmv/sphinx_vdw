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

template<uint32_t event, typename T1>
void SxEventBus::pub (const T1 &t1,
                      const SxString &pubTag,
                      const char *evName)
{
   SX_TRACE ();
   SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx) {
      uint32_t count = 0;
      const SxArray<SxEventListener*> &listener = eventBus.listener(event);
      for (SxEventListener *l : listener) {
         if (l == NULL) continue;
         SxPtr<SxEventData> data = SxPtr<SxEventData>::create (t1);
         data->setEventHash (event);
         SX_CHECK (data->getLayout () == eventBus.layoutLookup(event),
                   data->getLayout (),   eventBus.layoutLookup(event), event, pubTag);
         l->set<event> (SxEvent (data));
         ++count;
      }
      if (count == 0) {
         SxString msg = SxString::sprintf("0x%x", event) + " | " + SxString (evName) +
                        "@" + pubTag;
         eventBus.unreceived.append (msg);
      }
   }
}

template<uint32_t event, typename T1, typename T2>
void SxEventBus::pub (const T1 &t1,
                      const T2 &t2,
                      const SxString &pubTag,
                      const char *evName)
{
   SX_TRACE ();
   SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx) {
      uint32_t count = 0;
      const SxArray<SxEventListener*> &listener = eventBus.listener(event);
      for (SxEventListener *l : listener) {
         if (l == NULL) continue;
         SxPtr<SxEventData> data = SxPtr<SxEventData>::create (t1, t2);
         data->setEventHash (event);
         SX_CHECK (data->getLayout () == eventBus.layoutLookup(event),
                   data->getLayout (),   eventBus.layoutLookup(event), event, pubTag);
         l->set<event> (SxEvent (data));
         ++count;
      }

      if (count == 0) {
         SxString msg = SxString::sprintf("0x%x", event) + " | " + SxString (evName) +
                        "@" + pubTag;
         eventBus.unreceived.append (msg);
      }
   }
}

template<uint32_t event, typename T1, typename T2, typename T3>
void SxEventBus::pub (const T1 &t1,
                      const T2 &t2,
                      const T3 &t3,
                      const SxString &pubTag,
                      const char *evName)
{
   SX_TRACE ();
   SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx) {
      uint32_t count = 0;
      const SxArray<SxEventListener*> &listener = eventBus.listener(event);
      for (SxEventListener *l : listener) {
         if (l == NULL) continue;
         SxPtr<SxEventData> data = SxPtr<SxEventData>::create (t1, t2, t3);
         data->setEventHash (event);
         SX_CHECK (data->getLayout () == eventBus.layoutLookup(event),
                   data->getLayout (),   eventBus.layoutLookup(event), event, pubTag);
         l->set<event> (SxEvent (data));
         ++count;
      }

      if (count == 0) {
         SxString msg = SxString::sprintf("0x%x", event) + " | " + SxString (evName) +
                        "@" + pubTag;
         eventBus.unreceived.append (msg);
      }
   }
}

template<uint32_t event, typename T1, typename T2, typename T3,
                         typename T4>
void SxEventBus::pub (const T1 &t1,
                      const T2 &t2,
                      const T3 &t3,
                      const T4 &t4,
                      const SxString &pubTag,
                      const char *evName)
{
   SX_TRACE ();
   SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx) {
      uint32_t count = 0;
      const SxArray<SxEventListener*> &listener = eventBus.listener(event);
      for (SxEventListener *l : listener) {
         if (l == NULL) continue;
         SxPtr<SxEventData> data = SxPtr<SxEventData>::create (t1, t2, t3, t4);
         data->setEventHash (event);
         SX_CHECK (data->getLayout () == eventBus.layoutLookup(event),
                   data->getLayout (),   eventBus.layoutLookup(event), event, pubTag);
         l->set<event> (SxEvent (data));
         ++count;
      }

      if (count == 0) {
         SxString msg = SxString::sprintf("0x%x", event) + " | " + SxString (evName) +
                        "@" + pubTag;
         eventBus.unreceived.append (msg);
      }
   }
}

template<uint32_t event, typename T1, typename T2, typename T3,
                         typename T4, typename T5>
void SxEventBus::pub (const T1 &t1,
                      const T2 &t2,
                      const T3 &t3,
                      const T4 &t4,
                      const T5 &t5,
                      const SxString &pubTag,
                      const char *evName)
{
   SX_TRACE ();
   SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx) {
      uint32_t count = 0;
      const SxArray<SxEventListener*> &listener = eventBus.listener(event);
      for (SxEventListener *l : listener) {
         if (l == NULL) continue;
         SxPtr<SxEventData> data = SxPtr<SxEventData>::create (t1, t2, t3, t4, t5);
         data->setEventHash (event);
         SX_CHECK (data->getLayout () == eventBus.layoutLookup(event),
                   data->getLayout (),   eventBus.layoutLookup(event), event, pubTag);
         l->set<event> (SxEvent (data));
         ++count;
      }

      if (count == 0) {
         SxString msg = SxString::sprintf("0x%x", event) + " | " + SxString (evName) +
                        "@" + pubTag;
         eventBus.unreceived.append (msg);
      }
   }
}

template<uint32_t event, typename T1, typename T2, typename T3,
                         typename T4, typename T5, typename T6>
void SxEventBus::pub (const T1 &t1,
                      const T2 &t2,
                      const T3 &t3,
                      const T4 &t4,
                      const T5 &t5,
                      const T6 &t6,
                      const SxString &pubTag,
                      const char *evName)
{
   SX_TRACE ();
   SxEventBus &eventBus = getGlobalObj ();
   SX_MUTEX (eventBus.listenerMtx) {
      uint32_t count = 0;
      const SxArray<SxEventListener*> &listener = eventBus.listener(event);
      for (SxEventListener *l : listener) {
         if (l == NULL) continue;
         SxPtr<SxEventData> data = SxPtr<SxEventData>::create (t1, t2, t3, t4, t5, t6);
         data->setEventHash (event);
         SX_CHECK (data->getLayout () == eventBus.layoutLookup(event),
                   data->getLayout (),   eventBus.layoutLookup(event), event, pubTag);
         l->set<event> (SxEvent (data));
         ++count;
      }

      if (count == 0) {
         SxString msg = SxString::sprintf("0x%x", event) + " | " + SxString (evName) +
                        "@" + pubTag;
         eventBus.unreceived.append (msg);
      }
   }
}

template<uint32_t event>
inline void SxEventBus::subscribe (SxEventListener *listenerPtr,
                            const SxList<uint8_t> &layout)
{
   SX_TRACE ();
   subscribe (event, listenerPtr, layout);
}

template<uint32_t event>
inline void SxEventBus::unsubscribe (SxEventListener *listenerPtr)
{
   SX_TRACE ();
   unsubscribe (event, listenerPtr);
}

