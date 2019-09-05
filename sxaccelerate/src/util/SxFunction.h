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

#ifndef _SX_FUNCTION_H_
#define _SX_FUNCTION_H_

#include <SxPtr.h>
#include <SxCBoundPtr.h>
#include <SxSigSlots.h>


/** \brief Safebound pointer class

    \b SxFunction = S/PHI/nX Bound Pointer/Functors

    Bound pointers is a type- and thread-safe way of handling function
    pointers in S/PHI/nX. In contrast to C function pointers this class
    is also able to point to static and even non-static member functions
    of classes.

    This class is a safe wrapper for the (dangerous) low-level class
    SxCBoundPtr.  In contrast to SxCBoundPtr this class can keep track of
    the life time of the objects in case of pointer-to-members kind of
    bound pointers.

\code
#include <SxFunction.h>

// --- function pointer to a global function
// --- function pointer to a static member function
// --- function pointer to a non-static member function
\endcode


    \author Sixten Boeck, boeck@mpie.de */
template<class RES,
         class P1 = SxNull, class P2 = SxNull, class P3 = SxNull,
         class P4 = SxNull, class P5 = SxNull, class P6 = SxNull
        >
class SxFunction : public SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>
{


   public:

      typedef SxFunction<RES,P1,P2,P3,P4,P5,P6>  TBoundPtr;
      typedef SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6> TCBoundPtr;


   public:

      class Receiver
      {
         public:
            Receiver (bool valid_=false) : valid (valid_) { }
           ~Receiver () {
               if(isValid())  sigDeleting.send (this);
            }
            bool isValid () const { return valid; }

         signals:

            SxSignal<void *, const char *> SX_SIGNAL (sigDeleting);

         public slots:
            void slotDeleteObj (void *, const char *) { valid = false; }

         protected:
            bool valid;
      };
      SxPtr<Receiver> receiver;

      SxFunction ()
         : SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6> (),
           receiver (SxPtr<Receiver>::create (false))
      { }

     ~SxFunction () {  }

      SxFunction (const SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6> &in)
         : SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6> (in),
           receiver (SxPtr<Receiver>::create (true))
      {  }

      friend std::ostream
      &operator<< (std::ostream &os,
                   const SxFunction<RES,P1,P2,P3,P4,P5,P6> &bPtr_)
      {
#     ifndef NDEBUG
         if ((bPtr_.debugInf.getTag ())[0] != '\0') {
            os << "SxFunction : "
               << bPtr_.debugInf.getTag () << "\n";
         }
#     else
         SX_UNUSED (bPtr_);
         SX_EXIT;
#     endif
         return os;
      }
      // --------------------------------------------------------------------

      template<class TRES,class TP1>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (TRES (*f)(TP1),
              const char *tag_ = NULL)
      {
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (f, tag_)
         );
         return bPtr;
      }

      template<class TRES,class TP1,class TP2>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (TRES (*f)(TP1,TP2),
              const char *tag_ = NULL)
      {
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (f, tag_)
         );
         return bPtr;
      }

      template<class TRES,class TP1,class TP2,class TP3>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (TRES (*f)(TP1,TP2,TP3),
              const char *tag_ = NULL)
      {
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (f, tag_)
         );
         return bPtr;
      }

      template<class TRES,class TP1,class TP2,class TP3,class TP4>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (TRES (*f)(TP1,TP2,TP3,TP4),
              const char *tag_ = NULL)
      {
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (f, tag_)
         );
         return bPtr;
      }

      template<class TRES,class TP1,class TP2,class TP3,class TP4,class TP5>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (TRES (*f)(TP1,TP2,TP3,TP4,TP5),
              const char *tag_ = NULL)
      {
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (f, tag_)
         );
         return bPtr;
      }

      template<class TRES,
               class TP1,class TP2,class TP3,
               class TP4,class TP5,class TP6>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (TRES (*f)(TP1,TP2,TP3,TP4,TP5,TP6),
              const char *tag_ = NULL)
      {
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (f, tag_)
         );
         return bPtr;
      }

      // --------------------------------------------------------------------

      template <class Callee>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (Callee c, const char *tag_ = NULL)
      {
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (c, tag_)
         );
         return bPtr;
      }

      // --------------------------------------------------------------------

      template <class Callee,class TRES,class Class,class TP1>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,TRES (Class::* const &f)(TP1),
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,TRES (Class::* const &f)(TP1,TP2),
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,
                class TP1,class TP2,class TP3>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,TRES (Class::* const &f)(TP1,TP2,TP3),
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,
                class TP1,class TP2,class TP3,class TP4>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,TRES (Class::* const &f)(TP1,TP2,TP3,TP4),
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,
                class TP1,class TP2,class TP3,class TP4,class TP5>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,
              TRES (Class::* const &f)(TP1,TP2,TP3,TP4,TP5),
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,
                class TP1,class TP2,class TP3,class TP4,class TP5,class TP6>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,
              TRES (Class::* const &f)(TP1,TP2,TP3,TP4,TP5,TP6),
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      // --------------------------------------------------------------------

      template <class Callee,class TRES,class Class,class TP1>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,TRES (Class::* const &f)(TP1)const,
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,TRES (Class::* const &f)(TP1,TP2)const,
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,
                class TP1,class TP2,class TP3>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,
              TRES (Class::* const &f)(TP1,TP2,TP3)const,
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,
                class TP1,class TP2,class TP3,class TP4>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,
              TRES (Class::* const &f)(TP1,TP2,TP3,TP4)const,
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,
                class TP1,class TP2,class TP3,class TP4,class TP5>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,
              TRES (Class::* const &f)(TP1,TP2,TP3,TP4,TP5)const,
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      template <class Callee,class TRES,class Class,
                class TP1,class TP2,class TP3,class TP4,class TP5,class TP6>
      static SxFunction<RES,P1,P2,P3,P4,P5,P6>
      create (const SxPtr<Callee> &c,
              TRES (Class::* const &f)(TP1,TP2,TP3,TP4,TP5,TP6)const,
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,P1,P2,P3,P4,P5,P6> bPtr;
         bPtr = SxFunction<RES,P1,P2,P3,P4,P5,P6> (
            SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6>::create (*c.getPtr(), f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }


      // --------------------------------------------------------------------

      operator int () const
      {
         return receiver->isValid() && TCBoundPtr::operator int();
      }

      // --------------------------------------------------------------------


      RES operator() () const 
      {
         SX_CHECK (receiver->isValid()); // object has been deleted meanwhile!
         return TCBoundPtr::operator() ();
      }

      RES operator() (P1 p1) const 
      {
         SX_CHECK (receiver->isValid()); // object has been deleted meanwhile!
         return TCBoundPtr::operator() (p1);
      }


      RES operator() (P1 p1, P2 p2) const 
      {
         SX_CHECK (receiver->isValid()); // object has been deleted meanwhile!
         return TCBoundPtr::operator() (p1, p2);
      }

      RES operator() (P1 p1, P2 p2, P3 p3) const 
      {
         SX_CHECK (receiver->isValid()); // object has been deleted meanwhile!
         return TCBoundPtr::operator() (p1, p2, p3);
      }

      RES operator() (P1 p1, P2 p2, P3 p3, P4 p4) const 
      {
         SX_CHECK (receiver->isValid()); // object has been deleted meanwhile!
         return TCBoundPtr::operator() (p1, p2, p3, p4);
      }

      RES operator() (P1 p1, P2 p2, P3 p3, P4 p4, P5 p5) const 
      {
         SX_CHECK (receiver->isValid()); // object has been deleted meanwhile!
         return TCBoundPtr::operator() (p1, p2, p3, p4, p5);
      }

      RES operator() (P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6) const 
      {
         SX_CHECK (receiver->isValid()); // object has been deleted meanwhile!
         return TCBoundPtr::operator() (p1, p2, p3, p4, p5, p6);
      }
};



// SPECIAL CASE: void

template<class RES>
class SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> 
   : public SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
{

   public:

      typedef SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull>  TBoundPtr;
      typedef SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull> TCBoundPtr;

   public:

      class Receiver
      {
         public:
            Receiver (bool valid_=false) : valid (valid_) { }
           ~Receiver () {
               if (isValid())  sigDeleting.send (this);
            }
            bool isValid () const { return valid; }

         signals:
            SxSignal<void *, const char *> SX_SIGNAL (sigDeleting);

         public slots:
            void slotDeleteObj (void *, const char *) { printf ("DELETING\n"); valid = false; }

         protected:
            bool valid;
      };
      SxPtr<Receiver> receiver;

      SxFunction ()
         : SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull> (),
           receiver (SxPtr<Receiver>::create (false))
      { }

     ~SxFunction () {  }

      SxFunction (const SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull> &in)
         : SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull> (in),
           receiver (SxPtr<Receiver>::create (true))
      {  }

      friend std::ostream &operator<< (std::ostream &os,
                                       const SxFunction<RES,void,SxNull,
                                                        SxNull,SxNull,SxNull,
                                                        SxNull> &bPtr_)
      {
#     ifndef NDEBUG
         if ((bPtr_.debugInf.getTag ())[0] != '\0') {
            os << "SxFunction : "
               << bPtr_.debugInf.getTag () << "\n";
         }
#     else
         SX_UNUSED (bPtr_);
         SX_EXIT;
#     endif
         return os;
      }
      // --------------------------------------------------------------------

      template<class TRES>
      static SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull>
      create (TRES (*f)(void),
              const char *tag_ = NULL)
      {
         SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> bPtr;
         bPtr = SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> (
            SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>::create (f, tag_)
         );
         return bPtr;
      }

      // --------------------------------------------------------------------

      template <class Callee>
      static SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull>
      create (Callee c,
              const char *tag_ = NULL)
      {
         SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> bPtr;
         bPtr = SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> (
            SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>::create (c, tag_)
         );
         return bPtr;
      }

      // --------------------------------------------------------------------

      template <class Callee,class TRES,class Class>
      static SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull>
      create (const SxPtr<Callee> &c,TRES (Class::* const &f)(void),
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> bPtr;
         bPtr = SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> (
            SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>::create (*c.getPtr(),
                                                                                f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }

      // --------------------------------------------------------------------

      template <class Callee,class TRES,class Class>
      static SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull>
      create (const SxPtr<Callee> &c,TRES (Class::* const &f)(void)const,
              const char *tag_ = NULL)
      {
         SX_CHECK (c.getPtr());
         SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> bPtr;
         bPtr = SxFunction<RES,void,SxNull,SxNull,SxNull,SxNull,SxNull> (
            SxCBoundPtr<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>::create (*c.getPtr(),
                                                                                f, tag_)
         );
         c.receiver->sigDeleting.connect (bPtr.receiver.getPtr(),
                                 &TBoundPtr::Receiver::slotDeleteObj);
         bPtr.receiver->sigDeleting.connect (c.receiver,
                                   &SxPtr<Callee>::Receiver::slotDeleting);
         return bPtr;
      }


      // --------------------------------------------------------------------

      operator int () const
      {
         return receiver->isValid() && TCBoundPtr::operator int();
      }

      // --------------------------------------------------------------------


      RES operator() () const 
      {
         SX_CHECK (receiver->isValid()); // object has been deleted meanwhile!
         return TCBoundPtr::operator() ();
      }


};

#endif /* _SX_FUNCTION_H_ */
