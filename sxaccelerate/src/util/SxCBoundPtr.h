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

#ifndef _SX_CBOUND_PTR_H_
#define _SX_CBOUND_PTR_H_

#include <SxFunctor.h>


/** \brief Handle function pointers (UNSAFE!!!)

    \b SxCBoundPtr S/PHI/nX (unsafe) C-like bound-pointers

    \par Warning

    The usage of this class is unsafe. If not sure please use SxFunction
    instead! SxFunction monitors the object's lifetime.

    This class provides a convienient way of handle all function pointers
    in S/PHI/nX. Use this class to create function pointers to static
    global functions and to member functions.

\code
// --- functions with signature 'void foo(const SxString &)'
// (1) global function
void foo (const SxString &) { printf ("this is foo\n"); }

// (2) member function
class A
{
   public:
      void foo (const SxString &) { printf ("this is A::foo\n"); }
};




// --- define callbacks with proper signature
SxCBoundPtr<void,const SxString &> cb1, cb2;

cb1 = SxCBoundPtr<void,const SxString &>::create (&foo);

A a;
cb2 = SxCBoundPtr<void,const SxString &>::create (a, &A::foo);

// --- invoke callbacks
cb1 ("abc");    // => foo ("abc")
cb2 ("abc");    // => a.foo ("abc");
\endcode

\par Warning
If the provided object is deleted the callback becomes an invalid
bound pointer. Currently, this is not traced in the DEBUG mode!!!
Better use directly SxSignal instead. S/PHI/nX signals take care
of invalid bound pointers.

\author Sixten Boeck, boeck@mpie.de
*/
template<class RES,
         class P1 = SxNull, class P2 = SxNull, class P3 = SxNull,
         class P4 = SxNull, class P5 = SxNull, class P6 = SxNull
        >
class SxCBoundPtr : public SxFunctor<RES,P1,P2,P3,P4,P5,P6>
{
   public:

      SxCBoundPtr () { }

      template<class Callee>
      SxCBoundPtr (Callee c)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6>(SxLambdaTranslator<Callee,RES,
                                            P1,P2,P3,P4,P5,P6>(c, NULL)) { }

      template<class TRES>
      SxCBoundPtr (const SxFuncTranslator<TRES(*)(),RES,P1,P2,P3,P4,P5,P6> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }


      template<class TRES, class TP1>
      SxCBoundPtr (const SxFuncTranslator<TRES(*)(TP1),RES,P1,P2,P3,P4,P5,P6> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }


      template<class TRES, class TP1,class TP2>
      SxCBoundPtr (const SxFuncTranslator<TRES(*)(TP1,TP2),RES,P1,P2,P3,P4,P5,P6> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }


      template<class TRES, class TP1,class TP2,class TP3>
      SxCBoundPtr (const SxFuncTranslator<TRES(*)(TP1,TP2,TP3),RES,P1,P2,P3,P4,P5,P6> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }


      template<class TRES, class TP1,class TP2,class TP3, class TP4>
      SxCBoundPtr (const SxFuncTranslator<TRES(*)(TP1,TP2,TP3,TP4),RES,P1,P2,P3,P4,P5,P6> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }


      template<class TRES, class TP1,class TP2,class TP3, class TP4,class TP5>
      SxCBoundPtr (const SxFuncTranslator<TRES(*)(TP1,TP2,TP3,TP4,TP5),RES,P1,P2,P3,P4,P5,P6> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }


      template<class TRES,
               class TP1,class TP2,class TP3, class TP4,class TP5,class TP6>
      SxCBoundPtr (const SxFuncTranslator<TRES(*)(TP1,TP2,TP3,TP4,TP5,TP6),RES,P1,P2,P3,P4,P5,P6> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      friend std::ostream &operator<< (std::ostream &os,
                                       const SxCBoundPtr<RES,P1,P2,P3,P4,P5,P6> &cbPtr_)
      {
#     ifndef NDEBUG
         if ((cbPtr_.debugInf.getTag ())[0] != '\0') {
            os << "SxCBoundPtr : "
               << cbPtr_.debugInf.getTag () << "\n";
         }
#     else
         SX_UNUSED (cbPtr_);
         SX_EXIT;
#     endif
         return os;
      }

      // --------------------------------------------------------------------

      template<class Callee>
      SxCBoundPtr (const SxLambdaTranslator<Callee,RES,P1,P2,P3,P4,P5,P6> &in)
        : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      // --------------------------------------------------------------------

      template<class Callee, class TRES, class Class>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(),RES> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1),RES,P1> &in)
        : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2),RES,P1,P2> &in)
        : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2,class TP3>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3),RES,P1,P2,P3> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2,class TP3, class TP4>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4),RES,P1,P2,P3,P4> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2,class TP3, class TP4,class TP5>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4,TP5),RES,P1,P2,P3,P4,P5> &in)
        : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2,class TP3, class TP4,class TP5,class TP6>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4,TP5,TP6),RES,P1,P2,P3,P4,P5,P6> &in)
        : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      // --------------------------------------------------------------------

      template<class Callee, class TRES, class Class>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)()const,RES> &in)
        : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1)const,RES,P1> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2)const,RES,P1,P2> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2,class TP3>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3)const,RES,P1,P2,P3> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2,class TP3, class TP4>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4)const,RES,P1,P2,P3,P4> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2,class TP3, class TP4,class TP5>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4,TP5)const,RES,P1,P2,P3,P4,P5> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }

      template<class Callee, class TRES, class Class,
               class TP1,class TP2,class TP3, class TP4,class TP5,class TP6>
      SxCBoundPtr (const SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4,TP5,TP6)const,RES,P1,P2,P3,P4,P5,P6> &in)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (in) { /* empty */ }



      // --------------------------------------------------------------------

      template<class TRES, class TP1>
      static SxFuncTranslator<TRES (*)(TP1),RES,P1>
      create (TRES (*f)(TP1), const char *tag_ = NULL)
      {
         return SxFuncTranslator<TRES (*)(TP1),RES,P1> (f, tag_);
      }

      // --------------------------------------------------------------------

      template <class Callee>
      static SxLambdaTranslator<Callee,RES,P1,P2,P3,P4,P5,P6>
      create (Callee &c, const char *tag_ = NULL)
      {
         return SxLambdaTranslator<Callee,RES,P1,P2,P3,P4,P5,P6>(c, tag_);
      }

      // --------------------------------------------------------------------

      template <class Callee,class TRES,class Class>
      static SxMemTranslator<Callee,TRES (Class::*)(),RES,P1,P2,P3,P4,P5,P6>
      create (Callee &c,TRES (Class::* const &f)(), const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)();
         return SxMemTranslator<Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1>
      static SxMemTranslator<Callee,TRES (Class::*)(TP1),RES,P1,P2,P3,P4,P5,P6>
      create (Callee &c,TRES (Class::* const &f)(TP1), const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1);
         return SxMemTranslator<Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2>
      static SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2),RES,P1,P2,P3,P4,P5,P6>
      create (Callee &c,TRES (Class::* const &f)(TP1,TP2), const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2);
         return SxMemTranslator<Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2,class TP3>
      static SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3),RES,P1,P2,P3,P4,P5,P6>
      create (Callee &c,TRES (Class::* const &f)(TP1,TP2,TP3),
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2,TP3);
         return SxMemTranslator<Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2,class TP3,class TP4>
      static SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4),RES,P1,P2,P3,P4,P5,P6>
      create (Callee &c,TRES (Class::* const &f)(TP1,TP2,TP3,TP4),
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2,TP3,TP4);
         return SxMemTranslator<Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2,class TP3,class TP4,class TP5>
      static SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4,TP5),RES,P1,P2,P3,P4,P5,P6>
      create (Callee &c,TRES (Class::* const &f)(TP1,TP2,TP3,TP4,TP5),
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2,TP3,TP4,TP5);
         return SxMemTranslator<Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2,class TP3,class TP4,class TP5,class TP6>
      static SxMemTranslator<Callee,TRES (Class::*)(TP1,TP2,TP3,TP4,TP5,TP6),RES,P1,P2,P3,P4,P5,P6>
      create (Callee &c,TRES (Class::* const &f)(TP1,TP2,TP3,TP4,TP5,TP6),
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2,TP3,TP4,TP5,TP6);
         return SxMemTranslator<Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      // --------------------------------------------------------------------

      template <class Callee,class TRES,class Class>
      static SxMemTranslator<const Callee,TRES (Class::*)()const,RES,P1,P2,P3,P4,P5,P6>
      create (const Callee &c, TRES (Class::* const &f)()const,
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)()const;
         return SxMemTranslator<const Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1>
      static SxMemTranslator<const Callee,TRES (Class::*)(TP1)const,RES,P1,P2,P3,P4,P5,P6>
      create (const Callee &c, TRES (Class::* const &f)(TP1)const,
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1)const;
         return SxMemTranslator<const Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2>
      static SxMemTranslator<const Callee,TRES (Class::*)(TP1,TP2)const,RES,P1,P2,P3,P4,P5,P6>
      create (const Callee &c, TRES (Class::* const &f)(TP1,TP2)const,
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2)const;
         return SxMemTranslator<const Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2,class TP3>
      static SxMemTranslator<const Callee,TRES (Class::*)(TP1,TP2,TP3)const,RES,P1,P2,P3,P4,P5,P6>
      create (const Callee &c, TRES (Class::* const &f)(TP1,TP2,TP3)const,
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2,TP3)const;
         return SxMemTranslator<const Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2,class TP3,class TP4>
      static SxMemTranslator<const Callee,TRES (Class::*)(TP1,TP2,TP3,TP4)const,RES,P1,P2,P3,P4,P5,P6>
      create (const Callee &c, TRES (Class::* const &f)(TP1,TP2,TP3,TP4)const,
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2,TP3,TP4)const;
         return SxMemTranslator<const Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2,class TP3,class TP4,class TP5>
      static SxMemTranslator<const Callee,TRES (Class::*)(TP1,TP2,TP3,TP4,TP5)const,RES,P1,P2,P3,P4,P5,P6>
      create (const Callee &c, TRES (Class::* const &f)(TP1,TP2,TP3,TP4,TP5)const,
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2,TP3,TP4,TP5)const;
         return SxMemTranslator<const Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }

      template <class Callee,class TRES,class Class,class TP1,class TP2,class TP3,class TP4,class TP5,class TP6>
      static SxMemTranslator<const Callee,TRES (Class::*)(TP1,TP2,TP3,TP4,TP5,TP6)const,RES,P1,P2,P3,P4,P5,P6>
      create (const Callee &c, TRES (Class::* const &f)(TP1,TP2,TP3,TP4,TP5,TP6)const,
              const char *tag_ = NULL)
      {
         typedef TRES (Class::*MemFunc)(TP1,TP2,TP3,TP4,TP5,TP6)const;
         return SxMemTranslator<const Callee,MemFunc,RES,P1,P2,P3,P4,P5,P6>(c,f,tag_);
      }
};


#endif /* _SX_CBOUND_PTR_H_ */
