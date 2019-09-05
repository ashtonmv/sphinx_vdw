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

#ifndef _SX_SIGSLOTS_H_
#define _SX_SIGSLOTS_H_

#include <SxSignals.h>
#include <SxPtr.h>
#include <SxSlot.h>
#ifndef NDEBUG
#   include <SxError.h>
#endif
/** \brief Type safe signal slot mechanism

    connect signals and slots:

\code
#include <SxSigSlots.h>

   class MySlider
{
   public:
      MySlider () { }

   public slots:

      void slotABC (double val)
      {
         printf ("MySlider::slotABC received: %g\n", val);
      }

   signals:

      SxSignal<void,double> SX_SIGNAL(valueChanged);
};




class MyIsosurface
{
   public:
      MyIsosurface () { }

      void generateIso (double threshold)
      {
         printf ("MyIsosurface::slotABC received: %g\n", threshold);
      }

};


int main ()
{
   MySlider slider;
   SxPtr<MyIsosurface> isoPtr;

   {
      SxPtr<MyIsosurface> iso = SxPtr<MyIsosurface>::create ();
      sxconnect (&slider, valueChanged, iso, MyIsosurface::generateIso);

      isoPtr = iso;
      slider.valueChanged.send (5.0);
   }
\endcode
    \author Sixten Boeck, boeck@mpie.de */

template<class T>
struct SxRemoveRef { typedef T type; };

template<class T>
struct SxRemoveRef<T&> { typedef T type; };

template<class T>
struct SxRemoveRef<T&&> { typedef T type; };

#   define SX_GET_SLOT(PTR,SLOT) \
           PTR->SLOT##Link

#   define SX_SLOT(SLOT, ...) \
           SxSlot<__VA_ARGS__> SLOT##Link; \
           void SLOT (__VA_ARGS__);

#   define SX_LINK_SLOT(OBJPTR,SLOT)                \
           SLOT##Link(decltype(SLOT##Link)::create( \
           OBJPTR,&SxRemoveRef<decltype(*OBJPTR)>::type::SLOT))

#ifdef NDEBUG
#   define sxconnect(SIGNAL,SLOT)               \
       SxSigConnector::registerPtr(SIGNAL,SLOT);\
       SIGNAL.connect(&SLOT,                    \
                      SLOT.cbCheck.getCB());
#else
#   define sxconnect(SIGNAL,SLOT)                          \
       SxSigConnector::registerPtr(SIGNAL,SLOT);           \
       SIGNAL.setSigName(#SIGNAL);                         \
       SIGNAL.connect(&SLOT,                               \
                      SLOT.cbCheck.getCB(), #SLOT, SX_TAG);
#endif /* NDEBUG */

#define sxdisconnect(SIGNAL,SLOT)              \
   SxSigConnector::deregisterPtr(SIGNAL,SLOT); \
   SIGNAL.disconnect (&SLOT,                   \
                      SLOT.cbCheck.getCB());


/* \brief Signal connector

   \b SxSigConnector = S/PHI/nX Signal Slot Connector

   This wrapper class is used by the sxconnect macro.

   \author Sixten Boeck, boeck@mpie.de */
struct SxSigConnector
{

   template<class T>
   static void *getPtr(T in)
   {
      SX_CHECK (in == 0, in);
      return NULL;
   }

   template<class T>
   static T *getPtr(T *p)
   {
      SX_CHECK (p);
      return p;
   }

   template<class T>
   static T *getPtr(const SxPtr<T> &p)
   {
      SX_CHECK (p.getPtr());
      return p.getPtr();
   }

   template<class T>
   static void registerPtr (T in) { SX_CHECCK (in == 0, in); }

   template<class T>
   static void registerPtr (T *) { /* empty */ }

   template<class SIGNAL,class SLOT>
   static void registerPtr (SIGNAL &signal, SLOT &slot)
   {
      slot.registerSignal (signal);
      signal.registerLink (&slot, &SLOT::deregisterSignal);
   }

   template<class T>
   static void deregisterPtr (int) { /* empty */ }

   template<class T>
   static void deregisterPtr (T *) { /* empty */ }

   template<class SIGNAL,class SLOT>
   static void deregisterPtr (SIGNAL &signal, SLOT &slot)
   {
      slot.deregisterSignal (signal);
      signal.deregisterLink (&slot, &SLOT::deregisterSignal);
   }
};



#endif /* _SX_SIGSLOTS_H_ */
