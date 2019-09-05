
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

#ifndef _SX_SIGNALS_H_
#define _SX_SIGNALS_H_

#include <SxCBoundPtr.h>
#include <SxCList.h>
#include <iostream>
#include <string.h>
// --- support QT-like doxygen documentation 
#ifdef BUILD_DOXYGEN
#   define SX_SIGNAL(name)  name()
#else
#   define slots    /* empty */
#   undef signals
#   define signals  public /* empty */
#   define SX_SIGNAL(name)  name
#endif

#   ifndef NDEBUG
       class SxSigDbgInfo
       {
          public:
             SxSigDbgInfo () {
                tag[0] = slot[0] = slotTag[0] = '\0';
             }
             SxSigDbgInfo (const char *slot_, const char *tag_,
                           const char *signalName)
             {
                if (tag_ == NULL) {
                   tag[0] = '\0';
                } else {
                   strncpy (tag, tag_, 239);
                   tag[239] = '\0';
                }
                if (slot_ == NULL) {
                   slot[0] = '\0';
                } else {
                   strncpy (slot, slot_, 79);
                   slot[79] = '\0';
                }
                setSlotTag (signalName);
             }
             ~SxSigDbgInfo () {}
             bool operator == (const SxSigDbgInfo &inf) {
                if ((strncmp(tag, inf.tag, 240) == 0) &&
                    (strncmp(slot,inf.slot, 80) == 0))
                   return true;
                else
                   return false;
             }
             void setSlotTag (const char *signalName) {
                if (slot[0] == '\0') {
                   slotTag[0] = '\0';
                } else {
                   strcpy (slotTag, "'");
                   strcat (slotTag, slot);
                   strcat (slotTag, "' is connected to:\n- '");
                   strcat (slotTag, signalName);
                   strcat (slotTag, "' at ");
                   strcat (slotTag, tag);
                   slotTag[247] = '\0';
                }
             }
             const char *getTag () const {
                return tag;
             }
             const char *getSlot () const {
                return slot;
             }
             const char *getSlotTag () const {
                if (slotTag[0] == '\0')  return NULL;
                else                     return slotTag;
             }
          protected:
             char tag[240];
             char slot[80];
             char slotTag[248];
       };
#   endif /* NDEBUG */

/** \brief Signal/Slot handling

    \b SxSignals = S/PHI/nX Signals

    This class provides the base for a Signal/Slot mechanism in S/PHI/nX. 
\code
#include <SxSigSlots.h>
class MySlider
{
   public:
      MySlider ()
      {
         sxconnect (this, valueChanged, this, MySlider::slotABC);
      }

   public slots:
 
      void slotABC (double val)
      {
         printf ("AAA::slotABC received: %g\n", val);
      }

   public signals:
      SxSignal<void,double> SIGNAL(valueChanged);
};

class MyIsosurface
{
   public:
      MyIsosurface () { }
 
      void generateIso (double threshold)
      {
         printf ("BBB::slotABC received: %g\n", threshold);
      }

};


int main ()
{
   MySlider slider;
   MyIsosurface iso;

   sxconnect (&slider, valueChanged, &iso, MyIsosurface::generateIso);

// slider.valueChanged.connect (&fooIso);


   sxconnect (&slider, valueChanged, NULL, fooIso);
   slider.valueChanged.send (0.5);
}
\endcode
    \author Sixten Boeck, boeck@mpie.de */
// ---------------------------------------------------------------------------
template<class P1=SxNull, class P2=SxNull, class P3=SxNull, class P4=SxNull, class P5=SxNull, class P6=SxNull>
class SxSignal
{
   public:
      typedef SxCBoundPtr<void,P1,P2,P3,P4,P5,P6>  TBoundPtr;
      typedef SxCBoundPtr<void, SxSignal<P1,P2,P3,P4,P5,P6>&> LBoundPtr;

#   ifndef NDEBUG
       void setSigName (const char *sigName) {
          strncpy (signalName, sigName, 79);
          signalName[79] = '\0';
       }
#   endif /* NDEBUG */

      friend std::ostream &operator<< (std::ostream &os,
                                       const SxSignal<P1,P2,P3,P4,P5,P6> &sig_)
      {
#     ifndef NDEBUG
         os << "'" << sig_.signalName
            << "' is connected to: "
            << std::endl;

         SxCList<SxSigDbgInfo>::Node *iPtr = sig_.infoList.firstElement;

         while (iPtr) {
            if (((iPtr->elem).getSlot ())[0] !=  '\0') {
               os << "- '" << (iPtr->elem).getSlot () << "' at "
                  << (iPtr->elem).getTag () << std::endl;
            }
            iPtr = iPtr->next;
         }

#     else
         SX_UNUSED (sig_);
         SX_EXIT;
#     endif
         return os;
      }

      ~SxSignal () {
         typename SxCList<LBoundPtr>::Node *linkPtr = links.firstElement;
         while (linkPtr) {
            (linkPtr->elem)(*this);
            linkPtr = linkPtr->next;
         }
      }
      template<class DST, class Class,class TP1>
      void registerLink (DST *dst, void (Class::* const &slot)(TP1)) {
         LBoundPtr bPtr = LBoundPtr::create (*dst, slot);
         SX_CHECK (!links.contains (bPtr));
         links << bPtr;
      }
      template<class DST, class Class,class TP1>
      void deregisterLink (DST *dst, void (Class::* const &slot)(TP1)) {
         LBoundPtr bPtr = LBoundPtr::create (*dst, slot);
         SX_CHECK (links.contains (bPtr));
         links.removeElement (bPtr);
      }
      bool isConnected () const { return receivers.getSize() > 0; }

      // --- global function slots
      template<class TP1>
      void connect (void *dst, void (*slot)(TP1),
                    const char *slot_ = NULL,
                    const char *tag_ = NULL)
      {
         SX_CHECK (!dst);
         TBoundPtr bPtr = TBoundPtr::create (slot);
         SX_CHECK (!receivers.contains (bPtr));
         receivers << bPtr;
#   ifndef NDEBUG
       infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
       SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
      }
      template<class TP1,class TP2>
      void connect (void *dst, void (*slot)(TP1,TP2),
                    const char *slot_ = NULL,
                    const char *tag_ = NULL)
      {
         SX_CHECK (!dst);
         TBoundPtr bPtr = TBoundPtr::create (slot);
         SX_CHECK (!receivers.contains (bPtr));
         receivers << bPtr;
#   ifndef NDEBUG
         infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
         SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
      }
      template<class TP1,class TP2,class TP3>
      void connect (void *dst, void (*slot)(TP1,TP2,TP3),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         SX_CHECK (!dst);
         TBoundPtr bPtr = TBoundPtr::create (slot);
         SX_CHECK (!receivers.contains (bPtr));
         receivers << bPtr;
#   ifndef NDEBUG
         infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
         SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
      }
      template<class TP1,class TP2,class TP3,class TP4>
      void connect (void *dst, void (*slot)(TP1,TP2,TP3,TP4),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         SX_CHECK (!dst);
         TBoundPtr bPtr = TBoundPtr::create (slot);
         SX_CHECK (!receivers.contains (bPtr));
         receivers << bPtr;
#   ifndef NDEBUG
         infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
         SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
      }
      template<class TP1,class TP2,class TP3,class TP4,class TP5>
      void connect (void *dst, void (*slot)(TP1,TP2,TP3,TP4,TP5),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         SX_CHECK (!dst);
         TBoundPtr bPtr = TBoundPtr::create (slot);
         SX_CHECK (!receivers.contains (bPtr));
         receivers << bPtr;
#   ifndef NDEBUG
         infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
         SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
      }
      template<class TP1,class TP2,class TP3,class TP4,class TP5,class TP6>
      void connect (void *dst, void (*slot)(TP1,TP2,TP3,TP4,TP5,TP6),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         SX_CHECK (!dst);
         TBoundPtr bPtr = TBoundPtr::create (slot);
         SX_CHECK (!receivers.contains (bPtr));
         receivers << bPtr;
#   ifndef NDEBUG
         infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
         SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
      }

      // --- member function slots
      template<class DST,class Class,class TP1>
      void connect (DST *dst, void (Class::* const &slot)(TP1),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2>
      void connect (DST *dst, void (Class::* const &slot)(TP1,TP2),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2,class TP3>
      void connect (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2,class TP3,class TP4>
      void connect (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3,TP4),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2,class TP3,class TP4,
               class TP5>
      void connect (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3,TP4,TP5),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED (slot_, tag_);
#   endif
         }
      }
      template<class DST,class Class,class TP1,class TP2,class TP3,class TP4,
               class TP5,class TP6>
      void connect (DST *dst,
                    void (Class::* const &slot)(TP1,TP2,TP3,TP4,TP5,TP6),
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }

      // --- const member function slots
      template<class DST,class Class,class TP1>
      void connect (DST *dst, void (Class::* const &slot)(TP1)const,
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2>
      void connect (DST *dst, void (Class::* const &slot)(TP1,TP2)const,
                    const char *slot_ = NULL,const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2,class TP3>
      void connect (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3)const,
                    const char *slot_ = NULL,const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2,class TP3,class TP4>
      void connect (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3,TP4)const,
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2,class TP3,class TP4,
               class TP5>
      void connect (DST *dst,
                    void (Class::* const &slot)(TP1,TP2,TP3,TP4,TP5)const,
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }
      template<class DST,class Class,class TP1,class TP2,class TP3,
               class TP4,class TP5,class TP6>
      void connect (DST *dst,
                    void (Class::* const &slot)(TP1,TP2,TP3,TP4,TP5,TP6)const,
                    const char *slot_ = NULL, const char *tag_ = NULL)
      {
         TBoundPtr bPtr = TBoundPtr::create (*dst, slot);
         SX_CHECK (!receivers.contains (bPtr));
         if (!receivers.contains (bPtr)) {
            receivers << bPtr;
#   ifndef NDEBUG
            infoList << SxSigDbgInfo (slot_, tag_, signalName);
#   else
            SX_UNUSED(slot_, tag_);
#   endif /* NDEBUG */
         }
      }

      void slotDeregister (void *p, const char *)
      {
#   ifdef NDEBUG
          typename SxCList<TBoundPtr>::Node *rPtr = receivers.firstElement;
          while (rPtr) {
             if (rPtr->elem.getCallee () == p)  {
                typename SxCList<LBoundPtr>::Node *lPtr = links.firstElement;
                while (lPtr) {
                   if (lPtr->elem.getCallee () == rPtr->elem.getCallee ()) {
                      links.removeItem (lPtr);
                      break;
                   }
                }
                receivers.removeItem (rPtr);
                break;
             }
             rPtr = rPtr->next;
          }
#   else
          typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
          typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
          while (rPtr && infPtr) {
             if (rPtr->elem.getCallee () == p)  {
                typename SxCList<LBoundPtr>::Node *lPtr = links.firstElement;
                while (lPtr) {
                   if (lPtr->elem.getCallee () == rPtr->elem.getCallee ()) {
                      links.removeItem (lPtr);
                      break;
                   }
                }
                receivers.removeItem (rPtr);
                infoList.removeItem (infPtr);
                break;
             }
             rPtr   = rPtr->next;
             infPtr = infPtr->next;
          }
#   endif /* NDEBUG */
      }

#   ifndef NDEBUG
       void removeDebugInfo (const TBoundPtr &bPtr)
       {
          typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
          typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
          while (rPtr && infPtr) {
             if (rPtr->elem == bPtr)  {
                infoList.removeItem (infPtr);
                break;
             }
             rPtr   = rPtr->next;
             infPtr = infPtr->next;
          }
       }
#   endif /* NDEBUG */
      void disconnect (const TBoundPtr &cb)
      {
         SX_CHECK (receivers.contains (cb));
#   ifndef NDEBUG
       removeDebugInfo (cb);
#   endif /* NDEBUG */
         receivers.removeElement (cb);
      }

      template<class DST,class Class,class TP1>
      void disconnect (DST *dst, void (Class::* const &slot)(TP1))
      {
         disconnect (TBoundPtr::create (*dst, slot));
      }

      template<class DST,class Class,class TP1,class TP2>
      void disconnect (DST *dst, void (Class::* const &slot)(TP1,TP2))
      {
         disconnect (TBoundPtr::create (*dst, slot));
      }

      template<class DST,class Class,class TP1,class TP2,class TP3>
      void disconnect (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3))
      {
         disconnect (TBoundPtr::create (*dst, slot));
      }

      template<class DST,class Class,class TP1,class TP2,class TP3,class TP4>
      void disconnect (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3,TP4))
      {
         disconnect (TBoundPtr::create (*dst, slot));
      }

      template<class DST,class Class,
               class TP1,class TP2,class TP3,
               class TP4,class TP5>
      void disconnect (DST *dst,
                       void (Class::* const &slot)(TP1,TP2,TP3,TP4,TP5))
      {
         disconnect (TBoundPtr::create (*dst, slot));
      }

      template<class DST,class Class,
               class TP1,class TP2,class TP3,
               class TP4,class TP5,class TP6>
      void disconnect (DST *dst,
                       void (Class::* const &slot)(TP1,TP2,TP3,TP4,TP5,TP6))
      {
         disconnect (TBoundPtr::create (*dst, slot));
      }

      void send ()
      {
#   ifdef NDEBUG
         typename SxCList<TBoundPtr>::Node *rPtr = receivers.firstElement;
         while (rPtr) {
            (rPtr->elem)(NULL);
            rPtr = rPtr->next;
         }
#   else
         typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
         typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
         while (rPtr && infPtr) {
            (rPtr->elem)((infPtr->elem).getSlotTag ());
            rPtr   = rPtr->next;
            infPtr = infPtr->next;
         }
#   endif
      }

      template<class TP1>
      void send (const TP1 &p1)
      {

#   ifdef NDEBUG
         typename SxCList<TBoundPtr>::Node *rPtr = receivers.firstElement;
         while (rPtr) {
            (rPtr->elem)(p1, NULL);
            rPtr = rPtr->next;
         }
#   else
         typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
         typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
         while (rPtr && infPtr) {
            (rPtr->elem)(p1, (infPtr->elem).getSlotTag ());
            rPtr   = rPtr->next;
            infPtr = infPtr->next;
         }
#   endif
      }

      template<class TP1,class TP2>
      void send (const TP1 &p1, const TP2 &p2)
      {
#   ifdef NDEBUG
         typename SxCList<TBoundPtr>::Node *rPtr = receivers.firstElement;
         while (rPtr) {
            (rPtr->elem)(p1, p2, NULL);
            rPtr = rPtr->next;
         }
#   else
         typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
         typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
         while (rPtr && infPtr) {
            (rPtr->elem)(p1, p2, (infPtr->elem).getSlotTag ());
            rPtr   = rPtr->next;
            infPtr = infPtr->next;
         }
#   endif
      }

      template<class TP1,class TP2,class TP3>
      void send (const TP1 &p1, const TP2 &p2, const TP3 &p3)
      {
#   ifdef NDEBUG
         typename SxCList<TBoundPtr>::Node *rPtr = receivers.firstElement;
         while (rPtr) {
            (rPtr->elem)(p1, p2, p3, NULL);
            rPtr = rPtr->next;
         }
#   else
         typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
         typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
         while (rPtr && infPtr) {
            (rPtr->elem)(p1, p2, p3, (infPtr->elem).getSlotTag ());
            rPtr   = rPtr->next;
            infPtr = infPtr->next;
         }
#   endif
      }

      template<class TP1,class TP2,class TP3, class TP4>
      void send (const TP1 &p1, const TP2 &p2, const TP3 &p3, const TP4 &p4)
      {
#   ifdef NDEBUG
         typename SxCList<TBoundPtr>::Node *rPtr = receivers.firstElement;
         while (rPtr) {
            (rPtr->elem)(p1, p2, p3, p4, NULL);
            rPtr = rPtr->next;
         }
#   else
         typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
         typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
         while (rPtr && infPtr) {
            (rPtr->elem)(p1, p2, p3, p4, (infPtr->elem).getSlotTag ());
            rPtr   = rPtr->next;
            infPtr = infPtr->next;
         }
#   endif
      }

      template<class TP1,class TP2,class TP3, class TP4, class TP5>
      void send (const TP1 &p1, const TP2 &p2, const TP3 &p3, 
                 const TP4 &p4, const TP5 &p5)
      {
#   ifdef NDEBUG
         typename SxCList<TBoundPtr>::Node *rPtr = receivers.firstElement;
         while (rPtr) {
            (rPtr->elem)(p1, p2, p3, p4, p5, NULL);
            rPtr = rPtr->next;
         }
#   else
         typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
         typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
         while (rPtr && infPtr) {
            (rPtr->elem)(p1, p2, p3, p4, p5, (infPtr->elem).getSlotTag ());
            rPtr   = rPtr->next;
            infPtr = infPtr->next;
         }
#   endif
      }

      template<class TP1,class TP2,class TP3, class TP4, class TP5, class TP6>
      void send (const TP1 &p1, const TP2 &p2, const TP3 &p3, 
                 const TP4 &p4, const TP5 &p5, const TP6 &p6)
      {
#   ifdef NDEBUG
         typename SxCList<TBoundPtr>::Node *rPtr = receivers.firstElement;
         while (rPtr) {
            (rPtr->elem)(p1, p2, p3, p4, p5, p6, NULL);
            rPtr = rPtr->next;
         }
#   else
         typename SxCList<SxSigDbgInfo>::Node *infPtr = infoList.firstElement;
         typename SxCList<TBoundPtr>::Node *rPtr   = receivers.firstElement;
         while (rPtr && infPtr) {
            (rPtr->elem)(p1, p2, p3, p4, p5, p6, (infPtr->elem).getSlotTag ());
            rPtr   = rPtr->next;
            infPtr = infPtr->next;
         }
#   endif
      }

   protected:
      SxCList<TBoundPtr> receivers;
      SxCList<LBoundPtr> links;
#   ifndef NDEBUG
       char signalName[80] = {'\0'};
       SxCList<SxSigDbgInfo> infoList;
#   endif /* NDEBUG */
};


#endif /* _SX_SIGNALS_H_ */
