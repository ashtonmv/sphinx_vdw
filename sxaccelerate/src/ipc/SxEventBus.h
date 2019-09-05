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

#ifndef _SX_EVENT_BUS_H_
#define _SX_EVENT_BUS_H_

#include <SxIPC.h>
#include <SxMacroLib.h>
#include <SxEventListener.h>
#include <SxUniqueList.h>

/**
    \b SxEventBus = Global Event Bus

    A global event bus for forwarding
    events to the respective listeners.
    Anyone can publish; subscribe a
    SxEventListener in order to receive events.
*/

class SX_EXPORT_IPC SxEventBus
{
   public:

     ~SxEventBus ();

      friend std::ostream &operator<< (std::ostream &s, SxEventBus &bus);

      static SxEventBus &getGlobalObj ();

      template<uint32_t event,typename T1>
      static void pub (const T1 &t1,
                       const SxString &pubTag,
                       const char *evName);

      template<uint32_t event, typename T1, typename T2>
      static void pub (const T1 &t1,
                       const T2 &t2,
                       const SxString &pubTag,
                       const char *evName);

      template<uint32_t event, typename T1, typename T2, typename T3>
      static void pub (const T1 &t1,
                       const T2 &t2,
                       const T3 &t3,
                       const SxString &pubTag,
                       const char *evName);

      template<uint32_t event, typename T1, typename T2, typename T3,
                               typename T4>
      static void pub (const T1 &t1,
                       const T2 &t2,
                       const T3 &t3,
                       const T4 &t4,
                       const SxString &pubTag,
                       const char *evName);

      template<uint32_t event, typename T1, typename T2, typename T3,
                               typename T4, typename T5>
      static void pub (const T1 &t1,
                       const T2 &t2,
                       const T3 &t3,
                       const T4 &t4,
                       const T5 &t5,
                       const SxString &pubTag,
                       const char *evName);

      template<uint32_t event, typename T1, typename T2, typename T3,
                               typename T4, typename T5, typename T6>
      static void pub (const T1 &t1,
                       const T2 &t2,
                       const T3 &t3,
                       const T4 &t4,
                       const T5 &t5,
                       const T6 &t6,
                       const SxString &pubTag,
                       const char *evName);

      template<uint32_t event>
      static inline void subscribe (SxEventListener * listenerPtr,
                                    const SxList<uint8_t> &layout);

      template<uint32_t event>
      static inline void unsubscribe (SxEventListener *listenerPtr);

      static void unsubscribe (uint32_t evHash,
                               SxEventListener *listenerPtr);

      static void subscribe (uint32_t evHash,
                             SxEventListener *listenerPtr,
                             const SxArray<uint8_t> &layout);

      static void unsubscribeAll (SxEventListener *listenerPtr);

      static void terminate ();

   protected:
      using ListenerLookup = SxMap<uint32_t, SxArray<SxEventListener*> >;
      SxMutex listenerMtx;
      ListenerLookup listener;

      // --- possibly remove it in release
      SxMap<uint32_t, uint64_t> layoutLookup;
      SxUniqueList<SxString> unreceived;

      SxEventBus ();
};

std::ostream &operator<< (std::ostream &s, SxEventBus &bus);

// --- Subscribing and Publish Macro
#define SX_EVENT_BUS_SUB(EvName, layout, listener)                            \
      listener.setSubscribeTag (EvName ""_SX, SX_TAG);                        \
      SxEventBus::subscribe<EvName ""_SX>(&listener, SxList<uint8_t>() << layout)

#define _SX_EVENT_BUS_PUB2(EvName, arg1)                                      \
   SxEventBus::pub<EvName ""_SX> ((arg1), SX_TAG, EvName)

#define _SX_EVENT_BUS_PUB3(EvName, arg1, arg2)                                \
   SxEventBus::pub<EvName ""_SX> ((arg1), (arg2), SX_TAG, EvName)

#define _SX_EVENT_BUS_PUB4(EvName, arg1, arg2, arg3)                          \
   SxEventBus::pub<EvName ""_SX> ((arg1), (arg2), (arg3), SX_TAG, EvName)

#define _SX_EVENT_BUS_PUB5(EvName, arg1, arg2, arg3, arg4)                    \
   SxEventBus::pub<EvName ""_SX> ((arg1), (arg2), (arg3), (arg4), SX_TAG, EvName)

#define SX_EVENT_BUS_PUB6(EvName, arg1, arg2, arg3, arg4, arg5)               \
   SxEventBus::pub<EvName ""_SX> ((arg1), (arg2), (arg3), (arg4),             \
                                  (arg5), SX_TAG, EvName)

#define _SX_EVENT_BUS_PUB6(EvName, arg1, arg2, arg3, arg4, arg5, arg6)        \
   SxEventBus::pub<EvName ""_SX> ((arg1), (arg2), (arg3), (arg4),             \
                                  (arg5), (arg6), SX_TAG, EvName)

#define SX_EVENT_BUS_PUB(...) SX_VMACRO(_SX_EVENT_BUS_PUB, __VA_ARGS__)

#include <SxEventBus.hpp>

#endif /* _SX_EVENT_BUS_H */
