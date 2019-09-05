#ifndef _SX_SLOT_H_
#define _SX_SLOT_H_

#include <SxCBoundPtr.h>
#include <SxCList.h>

template<class P1=SxNull, class P2=SxNull, class P3=SxNull, class P4=SxNull,
         class P5=SxNull, class P6=SxNull>
class SxSlot
{

   private:
      SxSlot () {}
   public:
      typedef SxCBoundPtr<void,P1,P2,P3,P4,P5,P6> TBoundPtr;

      SxSlot (const SxSlot<P1,P2,P3,P4,P5,P6> &slot_) {
         *this = slot_;
      }
     ~SxSlot () {
         SxCList<SigBoundPtr>::Node *ptr = sigBPtrs.firstElement;
         while (ptr) {
            (ptr->elem)(this, "");
            ptr = ptr->next;
         }
      }
      template<class TP1,class TP2,class TP3,
               class TP4,class TP5,class TP6>
      struct SxCBCheck {
         typedef void (SxSlot::*memFunc) (P1,P2,P3,P4,P5,P6);
         memFunc getCB () {
            return &SxSlot<P1,P2,P3,P4,P5,P6>::callback<P1,P2,P3,P4,P5,P6>;
         }
      };
      template<class TP1,class TP2,class TP3,
               class TP4,class TP5>
      struct SxCBCheck<TP1,TP2,TP3,TP4,TP5,SxNull> {
         typedef void (SxSlot::*memFunc) (P1,P2,P3,P4,P5);
         memFunc getCB () {
            return &SxSlot<P1,P2,P3,P4,P5,P6>::callback<P1,P2,P3,P4,P5>;
         }
      };
      template<class TP1,class TP2,class TP3,
               class TP4>
      struct SxCBCheck<TP1,TP2,TP3,TP4,SxNull,SxNull> {
         typedef void (SxSlot::*memFunc) (P1,P2,P3,P4);
         memFunc getCB () {
            return &SxSlot<P1,P2,P3,P4,P5,P6>::callback<P1,P2,P3,P4>;
         }
      };
      template<class TP1,class TP2,class TP3>
      struct SxCBCheck<TP1,TP2,TP3,SxNull,SxNull,SxNull> {
         typedef void (SxSlot::*memFunc) (P1,P2,P3);
         memFunc getCB () {
            return &SxSlot<P1,P2,P3,P4,P5,P6>::callback<P1,P2,P3>;
         }
      };
      template<class TP1,class TP2>
      struct SxCBCheck<TP1,TP2,SxNull,SxNull,SxNull,SxNull> {
         typedef void (SxSlot::*memFunc) (P1,P2);
         memFunc getCB () {
            return &SxSlot<P1,P2,P3,P4,P5,P6>::callback<P1,P2>;
         }
      };
      template<class TP1>
      struct SxCBCheck<TP1,SxNull,SxNull,SxNull,SxNull,SxNull> {
         typedef void (SxSlot::*memFunc) (P1);
         memFunc getCB () {
            return &SxSlot<P1,P2,P3,P4,P5,P6>::callback<P1>;
         }
      };

      template<class TP1>
      void callback (TP1 p1) {
         slotBPtr (p1);
      }
      template<class TP1,class TP2>
      void callback (TP1 p1,TP2 p2) {
         slotBPtr (p1,p2);
      }
      template<class TP1,class TP2,class TP3>
      void callback (TP1 p1,TP2 p2,TP3 p3) {
         slotBPtr (p1,p2,p3);
      }
      template<class TP1,class TP2,class TP3,class TP4>
      void callback (TP1 p1,TP2 p2,TP3 p3,TP4 p4) {
         slotBPtr (p1,p2,p3,p4);
      }
      template<class TP1,class TP2,class TP3,class TP4,class TP5>
      void callback (TP1 p1,TP2 p2,TP3 p3,TP4 p4,TP5 p5) {
         slotBPtr (p1,p2,p3,p4,p5);
      }
      template<class TP1,class TP2,class TP3,class TP4,class TP5,class TP6>
      void callback (TP1 p1,TP2 p2,TP3 p3,TP4 p4,TP5 p5,TP6 p6) {
         slotBPtr (p1,p2,p3,p4,p5,p6);
      }

      template<class DST,class Class,class TP1>
      static SxSlot<P1,P2,P3,P4,P5,P6>
      create (DST *dst, void (Class::* const &slot)(TP1)) {
         SxSlot<P1,P2,P3,P4,P5,P6> slotObj;
         slotObj.slotBPtr = TBoundPtr::create (*dst,slot);
         return slotObj;
      }
      template<class DST,class Class,class TP1,class TP2>
      static SxSlot<P1,P2,P3,P4,P5,P6>
      create (DST *dst, void (Class::* const &slot)(TP1,TP2)) {
         SxSlot<P1,P2,P3,P4,P5,P6> slotObj;
         slotObj.slotBPtr = TBoundPtr::create (*dst,slot);
         return slotObj;
      }
      template<class DST,class Class,class TP1,class TP2,class TP3>
      static SxSlot<P1,P2,P3,P4,P5,P6>
      create (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3)) {
         SxSlot<P1,P2,P3,P4,P5,P6> slotObj;
         slotObj.slotBPtr = TBoundPtr::create (*dst,slot);
         return slotObj;
      }
      template<class DST,class Class,class TP1,
               class TP2,class TP3,class TP4>
      static SxSlot<P1,P2,P3,P4,P5,P6>
      create (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3,TP4)) {
         SxSlot<P1,P2,P3,P4,P5,P6> slotObj;
         slotObj.slotBPtr = TBoundPtr::create (*dst,slot);
         return slotObj;
      }
      template<class DST,class Class,class TP1,
               class TP2,class TP3,class TP4,class TP5>
      static SxSlot<P1,P2,P3,P4,P5,P6>
      create (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3,TP4,TP5)) {
         SxSlot<P1,P2,P3,P4,P5,P6> slotObj;
         slotObj.slotBPtr = TBoundPtr::create (*dst,slot);
         return slotObj;
      }
      template<class DST,class Class,class TP1,
               class TP2,class TP3,class TP4,class TP5,class TP6>
      static SxSlot<P1,P2,P3,P4,P5,P6>
      create (DST *dst, void (Class::* const &slot)(TP1,TP2,TP3,TP4,TP5,TP6)) {
         SxSlot<P1,P2,P3,P4,P5,P6> slotObj;
         slotObj.slotBPtr = TBoundPtr::create (*dst,slot);
         return slotObj;
      }

      void registerSignal (SxSignal<P1,P2,P3,P4,P5,P6> &signal) {
         SigBoundPtr bPtr = SigBoundPtr::create
                            (signal, &SxSignal<P1,P2,P3,P4,P5,P6>::slotDeregister);
         SX_CHECK (!sigBPtrs.contains (bPtr));
         sigBPtrs << bPtr;
      }
      void deregisterSignal (SxSignal<P1,P2,P3,P4,P5,P6> &signal) {
         SigBoundPtr bPtr = SigBoundPtr::create
                            (signal, &SxSignal<P1,P2,P3,P4,P5,P6>::slotDeregister);
         SX_CHECK (sigBPtrs.contains (bPtr));
         sigBPtrs.removeElement(bPtr);
      }
      struct SxCBCheck<P1,P2,P3,P4,P5,P6> cbCheck;

   protected:

      typedef SxCBoundPtr<void,void*,const char*,
                          SxNull,SxNull,SxNull,SxNull> SigBoundPtr;

      TBoundPtr slotBPtr;
      SxCList<SigBoundPtr> sigBPtrs;

};

#endif //_SX_SLOT_H_
