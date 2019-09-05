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

#ifndef _SX_FUNCTOR_H_
#define _SX_FUNCTOR_H_

#include <string.h>
#include <SxError.h>
#include <SxConfig.h>
#include <iostream>


#   ifndef NDEBUG
       class SxFunctorDbgInfo
       {
          public:
             SxFunctorDbgInfo () {
                tag[0] = '\0';
             }
             SxFunctorDbgInfo (const char *tag_)
             {
                if (tag_ == NULL) {
                   tag[0] = '\0';
                } else {
                   strncpy (tag, tag_, 239);
                   tag[239] = '\0';
                }
             }
             ~SxFunctorDbgInfo () {}
             bool operator == (const SxFunctorDbgInfo &inf) {
                if ((strncmp(tag, inf.tag, 240) == 0))
                   return true;
                else
                   return false;
             }
             void setTag (const char *tag_) {
                if (tag_ == NULL) {
                   tag[0] = '\0';
                } else {
                   strncpy (tag, tag_, 239);
                   tag[239] = '\0';
                }
             }
             const char * getTag () const {
                return tag;
             }
          protected:
             char tag[240];
       };
#   endif /* NDEBUG */




/** \brief Empty data type to handle variadic template classes

    Use this class to define variatic template classes. For example
\code
template<class A, class B=SxNull,class C=SxNull>
class MyClass { };

template<class A>
class MyClass<A,SxNull,SxNull>
{
   public:
      someCode
};

template<class A, class B>
class MyClass<A,B,SxNull>
{
   public:
      someCode
};
\endcode
\author Sixten Boeck, boeck@mpie.de */
class SxNull { };

// --------------------------------------------------------------------------

/** \brief Basis class for all functors

    \b SxFunctorStorage = S/PHI/nX Functor Storage Base class

       this class serves as a base class for SxLambdaStorage and
       SxBoundStorage classes.
 */

class SxFunctorStorage {
   public:

      /** \brief Type of a pointer to member functions

          According to the C++ standard, all pointer-to-member-functions
          can be byte-converted via reinterpret_cast<> to
          pointer-to-member-functions of a different class and signature.
          Back-conversion yields the original pointer. We use this feature to
          store any pointer-to-member-function in a
          void (SxFunctorStorage::*)() called "MemFunc". It can never be used
          as such, only after type-converting to the original class. This is
          done in SxMemTranslator::metaFunc .
          For more on the intricacies of the pointer-to-member-functions, see

          https://www.codeproject.com/Articles/7150/Member-Function-Pointers-and-the-Fastest-Possible

        */
      typedef void (SxFunctorStorage::*MemFunc)();

      /** \brief Type of a pointer to static global functions
        */
      typedef void (*Func)();


      virtual operator int() const = 0;
      virtual void *getCallee () const = 0;
      virtual bool operator== (const SxFunctorStorage &) const = 0;
      virtual SxFunctorStorage *clone () const = 0;
      virtual ~SxFunctorStorage () {}
};

/** \brief Basis class for lambda functor storage

    \b SxLambdaStorage = S/PHI/nX Lambda Functor Storage class

       this class serves as a storage class for Lambda functor.
 */

template<class Callee>
class SxLambdaStorage : public SxFunctorStorage
{
   public:

      SxLambdaStorage () { }

     /** \param callee_   Lambda object
          \param func_    Pointer to member function, currently unused
      */
      SxLambdaStorage (const Callee &callee_, const MemFunc &func_)
         : callee(callee_)
      {
         memFunc = func_;
      }

      SxLambdaStorage (const SxLambdaStorage<Callee> &in)
         : callee(in.callee)
      {
         memFunc = in.memFunc;
      }

      SxLambdaStorage &operator= (const SxLambdaStorage<Callee> &in)
      {
         if (this == &in) return *this;
         // copy the lambda object
         callee = in.callee;
         memFunc = in.memFunc;
         return *this;
      }

      /** \brief support 'if (callbackPtr)'
        */
      operator int() const override
      {
         return 1; // always true for lambda
      }

      bool operator== (const SxLambdaStorage<Callee> &in)
      {
         if (this == &in) return true;
         return false;
      }

      bool operator== (const SxFunctorStorage &in) const override {
         const SxLambdaStorage<Callee> *fm = dynamic_cast<const SxLambdaStorage<Callee>*>(&in);
         if (!fm)  return false; // could not cast
         if (this == fm) return true;
         return false;
         // lambda comparison causes an error on clang++
         //return (callee == (fm->callee) && memFunc == (MemFunc)(fm->memFunc));
      }

      bool operator!= (const SxLambdaStorage<Callee> &in)
      {
         return !operator== (in);
      }

      SxLambdaStorage<Callee> *clone () const override {
         SxLambdaStorage<Callee> *ptr = new SxLambdaStorage<Callee> (*this);
         return ptr;
      }

      // this function is called from SxSignal
      void *getCallee () const override {
         return NULL;
      }

     ~SxLambdaStorage<Callee> () override {}

      /** \brief copy of lambda object
        */
      Callee callee;

      /** \brief function or member function data
        */
      MemFunc memFunc;
};

/** \brief Storage class for global and Member functor storage

    \b SxBoundStorage = S/PHI/nX Global and Member Functor Storage class

 */

class SxBoundStorage : public SxFunctorStorage
{

   public:

      SxBoundStorage () : callee(0), func(0) {  }

      /** \param func_    Pointer to global function
          */
      SxBoundStorage (Func *func_)
         : callee(0),
           func(func_)
      {  }

     /** \param callee_   Pointer to object or NULL for global functions
          \param func_    Pointer to member function
      */
      SxBoundStorage (const void *callee_, const MemFunc &func_)
         : callee(const_cast<void *>(callee_)),
           func(0)
      {
         SX_CHECK (callee_);
         memFunc = func_;
      }

      SxBoundStorage (const SxBoundStorage &in)
         : callee(in.callee)
      {
         if (callee)  {  // member function pointer
            memFunc = in.memFunc;
         } else {        // global function pointer
            func = in.func;
         }
      }

      SxBoundStorage &operator= (const SxBoundStorage &in)
      {
         if (this == &in) return *this;
         callee = in.callee;
         if (callee)  {  // member function pointer
            memFunc = in.memFunc;
         } else {        // global function pointer
            func = in.func;
         }
         return *this;
      }

      /** \brief support 'if (callbackPtr)'
        */
      operator int() const override
      {
         return (callee == NULL && func != NULL) || (callee != NULL);
      }

      bool operator== (const SxBoundStorage &in)
      {
         if (!callee)  {  // C function
            return (func == in.func && !in.callee);
         }  else  {  // member function
            return (callee == in.callee && memFunc == in.memFunc);
         }
      }

      bool operator== (const SxFunctorStorage &in) const override {
         const SxBoundStorage *fm = dynamic_cast<const SxBoundStorage *>(&in);
         if (!fm)  return false; // could not cast
         if (!callee)  {  // C function
            return (func == fm->func && !fm->callee);
         }  else  {  // member function
            return (callee == fm->callee && memFunc == fm->memFunc);
         }
      }

      bool operator!= (const SxBoundStorage &in)
      {
         return !operator== (in);
      }

      SxBoundStorage *clone () const override {
         SxBoundStorage *ptr = new SxBoundStorage (*this);
         return ptr;
      }

      void *getCallee () const override {
         return callee;
      }

     ~SxBoundStorage () override {}

      /** \brief Pointer to object or NULL for global functions
        */
      void *callee;

      /** \brief function or member function data
        */
      union {
         void *func;
         MemFunc memFunc;
      };
};



// --------------------------------------------------------------------------

/** \brief Basis class for all functors

    \b SxFunctorBase = S/PHI/nX Functor Base class

    This functor class is a container for bound function pointers. It stores
    the object (callee), the member function or the static global function.
    The information is anonymized.

    \author Sixten Boeck, boeck@mpie.de */
class SxFunctorBase
{

#   ifndef NDEBUG
       public:
          SxFunctorDbgInfo debugInf;
#   endif /* NDEBUG */

   public:
      SxFunctorBase () : fBPtr(NULL) { }
      SxFunctorBase (const SxFunctorStorage *fb_, const char *tag_=NULL) : fBPtr(fb_)
      {
#   ifndef NDEBUG
         setDebugInfo (tag_);
#   else
         SX_UNUSED(tag_);
#   endif /* NDEBUG */
      }
      SxFunctorBase (const SxFunctorBase &in)
          : fBPtr(in.fBPtr?in.fBPtr->clone():NULL) {
#   ifndef NDEBUG
         setDebugInfo (in.debugInf.getTag());
#   endif /* NDEBUG */
      }

      SxFunctorBase &operator= (const SxFunctorBase &in) {
         if (!in.fBPtr) return *this;
         if (*this == in) return *this;
         if (fBPtr) delete fBPtr;

         fBPtr = in.fBPtr->clone (); // call corresponding subclass.clone()
#   ifndef NDEBUG
         setDebugInfo (in.debugInf.getTag());
#   endif /* NDEBUG */
         return *this;
      }

#      ifndef NDEBUG
       void setDebugInfo (const char *tag_) {
          debugInf.setTag (tag_);
       }
#      endif /* NDEBUG */

      friend std::ostream &operator<< (std::ostream &os,
                                       const SxFunctorBase &functor_);

      bool operator== (const SxFunctorBase &in) {
         if (!fBPtr || !in.fBPtr)  return false;

         return (*fBPtr == *in.fBPtr);
      }

      bool operator!= (const SxFunctorBase &in) {
         return !(*this == in);
      }

      operator int() const {
         return (fBPtr && int(*fBPtr));
      }

      void *getCallee () const {
         if(!fBPtr)  return NULL;
         return fBPtr->getCallee ();
      }

     ~SxFunctorBase () { if (fBPtr) delete fBPtr; fBPtr = NULL; }

   protected:
      const SxFunctorStorage *fBPtr;
};

inline std::ostream &operator<< (std::ostream &os, const SxFunctorBase &functor_)
{
#  ifndef NDEBUG
      if ((functor_.debugInf.getTag ())[0] != '\0') {
         os << "SxFunctorStorage : "
            << functor_.debugInf.getTag () << "\n";
      }
#  else
      SX_UNUSED (functor_);
      SX_EXIT;
#  endif
   return os;
}

// --------------------------------------------------------------------------

template<class RES, 
         class P1 = SxNull, class P2 = SxNull, class P3 = SxNull, 
         class P4 = SxNull, class P5 = SxNull, class P6 = SxNull,
         class Terminator = SxNull
        >
class SxFunctor : public SxFunctorBase
{
   // empty
};

// --- SxFunctor<RES>
template<class RES>
class SxFunctor<RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctorBase
{
   public:
      SxFunctor () : SxFunctorBase () { /* empty */ }
      SxFunctor (const SxFunctor &in) : SxFunctorBase(in), metaFunc(in.metaFunc)
                                        { }
      SxFunctor &operator= (const SxFunctor &in) {
         SxFunctorBase::operator= (in);
         metaFunc = in.metaFunc;
         return *this;
      }

      bool operator== (const SxFunctor &in) {
         return SxFunctorBase::operator== (in);
      }

      bool operator!= (const SxFunctor &in) {
         return SxFunctorBase::operator!= (in);
      }

      using SxFunctorBase::operator int;

     ~SxFunctor () { }

      RES operator() () const {
         return metaFunc (fBPtr);
      }

   protected:

      // meta function pointer (function pointer + functor basis)
      typedef RES ( *MetaFunc)(const SxFunctorStorage *);

      SxFunctor (MetaFunc mf_, const SxFunctorStorage *fb_, const char *tag_=NULL)
          : SxFunctorBase(fb_, tag_), metaFunc(mf_) {  }

   private:
      MetaFunc metaFunc;  // function pointer + functor basis
};

// --- SxFunctor<RES,P1>
template<class RES, class P1>
class SxFunctor<RES,P1,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctorBase
{
   public:
      SxFunctor () : SxFunctorBase () { /* empty */ }
      SxFunctor (const SxFunctor &in) : SxFunctorBase(in), metaFunc(in.metaFunc)
                                        { }
      SxFunctor &operator= (const SxFunctor &in) {
         SxFunctorBase::operator= (in);
         metaFunc = in.metaFunc;
         return *this;
      }

      bool operator== (const SxFunctor &in) {
         return SxFunctorBase::operator== (in);
      }

      bool operator!= (const SxFunctor &in) {
         return SxFunctorBase::operator!= (in);
      }

      using SxFunctorBase::operator int;

     ~SxFunctor () { }
      RES operator() (P1 p1) const {
            return metaFunc (fBPtr, p1);
      }

   protected:

      // meta function pointer (function pointer + functor basis)
      typedef RES (*MetaFunc)(const SxFunctorStorage *, P1);

      SxFunctor (MetaFunc mf_, const SxFunctorStorage *fb_, const char *tag_=NULL)
          : SxFunctorBase(fb_, tag_), metaFunc(mf_) {  }

   private:
      MetaFunc metaFunc;  // function pointer + functor basis
};

// --- SxFunctor<RES,P1,P2>
template<class RES, class P1, class P2>
class SxFunctor<RES,P1,P2,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctorBase
{
   public:
      SxFunctor () : SxFunctorBase () { /* empty */ }
      SxFunctor (const SxFunctor &in) : SxFunctorBase(in), metaFunc(in.metaFunc)
                                        { }
      SxFunctor &operator= (const SxFunctor &in) {
         SxFunctorBase::operator= (in);
         metaFunc = in.metaFunc;
         return *this;
      }

      bool operator== (const SxFunctor &in) {
         return SxFunctorBase::operator== (in);
      }

      bool operator!= (const SxFunctor &in) {
         return SxFunctorBase::operator!= (in);
      }

      using SxFunctorBase::operator int;

     ~SxFunctor () { }

      RES operator() (P1 p1, P2 p2) const {
         return metaFunc (fBPtr, p1, p2);
      }

   protected:

      // meta function pointer (function pointer + functor basis)
      typedef RES (*MetaFunc)(const SxFunctorStorage *, P1, P2);

      SxFunctor (MetaFunc mf_, const SxFunctorStorage *fb_, const char *tag_=NULL)
          : SxFunctorBase(fb_, tag_), metaFunc(mf_) {  }

   private:
      MetaFunc metaFunc;  // function pointer + functor basis
};

// --- SxFunctor<RES,P1,P2,P3>
template<class RES, class P1, class P2, class P3>
class SxFunctor<RES,P1,P2,P3,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctorBase
{
   public:
      SxFunctor () : SxFunctorBase () { /* empty */ }
      SxFunctor (const SxFunctor &in) : SxFunctorBase(in), metaFunc(in.metaFunc)
                                        { }
      SxFunctor &operator= (const SxFunctor &in) {
         SxFunctorBase::operator= (in);
         metaFunc = in.metaFunc;
         return *this;
      }

      bool operator== (const SxFunctor &in) {
         return SxFunctorBase::operator== (in);
      }

      bool operator!= (const SxFunctor &in) {
         return SxFunctorBase::operator!= (in);
      }

      using SxFunctorBase::operator int;

     ~SxFunctor () { }

      RES operator() (P1 p1, P2 p2, P3 p3) const {
         return metaFunc (fBPtr, p1, p2, p3);
      }

   protected:

      // meta function pointer (function pointer + functor basis)
      typedef RES ( *MetaFunc)(const SxFunctorStorage *, P1, P2, P3);

      SxFunctor (MetaFunc mf_, const SxFunctorStorage *fb_, const char *tag_=NULL)
          : SxFunctorBase(fb_, tag_), metaFunc(mf_) {  }

   private:
      MetaFunc metaFunc;  // function pointer + functor basis
};

// --- SxFunctor<RES,P1,P2,P3,P4>
template<class RES, class P1, class P2, class P3, class P4>
class SxFunctor<RES,P1,P2,P3,P4,SxNull,SxNull,SxNull>
   : public SxFunctorBase
{
   public:
      SxFunctor () : SxFunctorBase () { /* empty */ }
      SxFunctor (const SxFunctor &in) : SxFunctorBase(in), metaFunc(in.metaFunc)
                                        { }
      SxFunctor &operator= (const SxFunctor &in) {
         SxFunctorBase::operator= (in);
         metaFunc = in.metaFunc;
         return *this;
      }

      bool operator== (const SxFunctor &in) {
         return SxFunctorBase::operator== (in);
      }

      bool operator!= (const SxFunctor &in) {
         return SxFunctorBase::operator!= (in);
      }

      using SxFunctorBase::operator int;

     ~SxFunctor () { }
      RES operator() (P1 p1, P2 p2, P3 p3, P4 p4) const {
        return metaFunc (fBPtr, p1, p2, p3, p4);
      }

   protected:

      // meta function pointer (function pointer + functor basis)
      typedef RES ( *MetaFunc)(const SxFunctorStorage *, P1, P2, P3 , P4);

      SxFunctor (MetaFunc mf_, const SxFunctorStorage *fb_, const char *tag_=NULL)
          : SxFunctorBase(fb_, tag_), metaFunc(mf_) {  }

   private:
      MetaFunc metaFunc;  // function pointer + functor basis
};

// --- SxFunctor<RES,P1,P2,P3,P4,P5>
template<class RES, class P1, class P2, class P3, class P4, class P5>
class SxFunctor<RES,P1,P2,P3,P4,P5,SxNull,SxNull>
   : public SxFunctorBase
{
   public:
      SxFunctor () : SxFunctorBase () { /* empty */ }
      SxFunctor (const SxFunctor &in) : SxFunctorBase(in), metaFunc(in.metaFunc)
                                        { }
      SxFunctor &operator= (const SxFunctor &in) {
         SxFunctorBase::operator= (in);
         metaFunc = in.metaFunc;
         return *this;
      }

      bool operator== (const SxFunctor &in) {
         return SxFunctorBase::operator== (in);
      }

      bool operator!= (const SxFunctor &in) {
         return SxFunctorBase::operator!= (in);
      }

      using SxFunctorBase::operator int;

     ~SxFunctor () { }

      RES operator() (P1 p1, P2 p2, P3 p3, P4 p4, P5 p5) const {
            return metaFunc (fBPtr, p1, p2, p3, p4, p5);
      }

   protected:

      // meta function pointer (function pointer + functor basis)
      typedef RES ( *MetaFunc)(const SxFunctorStorage &, P1, P2, P3, P4, P5);

      SxFunctor (MetaFunc mf_, const SxFunctorStorage *fb_, const char *tag_=NULL)
          : SxFunctorBase(fb_, tag_), metaFunc(mf_) {  }

   private:
      MetaFunc metaFunc;  // function pointer + functor basis
};

// --- SxFunctor<RES,P1,P2,P3,P4,P5,P6>
template<class RES, class P1, class P2, class P3, class P4, class P5, class P6>
class SxFunctor<RES,P1,P2,P3,P4,P5,P6,SxNull>
   : public SxFunctorBase
{
   public:
      SxFunctor () : SxFunctorBase () { /* empty */ }
      SxFunctor (const SxFunctor &in) : SxFunctorBase(in), metaFunc(in.metaFunc)
                                        { }
      SxFunctor &operator= (const SxFunctor &in) {
         SxFunctorBase::operator= (in);
         metaFunc = in.metaFunc;
         return *this;
      }

      bool operator== (const SxFunctor &in) {
         return SxFunctorBase::operator== (in);
      }

      bool operator!= (const SxFunctor &in) {
         return SxFunctorBase::operator!= (in);
      }

      using SxFunctorBase::operator int;

     ~SxFunctor () { }


      RES operator() (P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6) const {
            return metaFunc (fBPtr, p1, p2, p3, p4, p5, p6);
      }

   protected:

      // meta function pointer (function pointer + functor basis)
      typedef RES ( *MetaFunc)(const SxFunctorStorage &, P1, P2, P3, P4, P5, P6);

      SxFunctor (MetaFunc mf_, const SxFunctorStorage *fb_, const char *tag_=NULL)
          : SxFunctorBase(fb_, tag_), metaFunc(mf_) {  }

   private:
      MetaFunc metaFunc;  // function pointer + functor basis
};
//----------------------------------------------------------------------------

// translate function pointers to functors

// SxFuncTranslator<Func, RES,...>
template<class Func, class RES, 
         class P1 = SxNull, class P2 = SxNull, class P3 = SxNull, 
         class P4 = SxNull, class P5 = SxNull, class P6 = SxNull,
         class Terminator = SxNull
        >
class SxFuncTranslator : public SxFunctor<RES,P1,P2,P3,P4,P5,P6,Terminator>
{ /* empty */ };

// SxFuncTranslator<Func, RES>
template<class Func, class RES>
class SxFuncTranslator<Func,RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES>
{
   public:
      SxFuncTranslator (Func func_, const char *tag_ = NULL)
         : SxFunctor<RES> (metaFunc,
                          (new SxBoundStorage(
                          (SxFunctorStorage::Func *)const_cast<void *>(
                          (const void *)func_) )), tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*>(functor_);
         SX_CHECK (ftr->func);
#        ifdef WIN32
            return (reinterpret_cast<RES (*)()>(ftr->func))();
#        else
            return (Func (ftr->func))();
#        endif
      }
};

// SxFuncTranslator<Func,RES,P1>
template<class Func, class RES, class P1>
class SxFuncTranslator<Func,RES,P1,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1>
{
   public:

      SxFuncTranslator (RES (*func_)(P1), const char *tag_ = NULL)
         : SxFunctor<RES,P1> (metaFunc,
                             (new SxBoundStorage(
                             (SxFunctorStorage::Func *)const_cast<void *>(
                             (const void *)func_) )), tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, P1 p1)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*>(functor_);
         SX_CHECK (ftr->func);
#        ifdef WIN32
            return (reinterpret_cast<RES (*)(P1)>(ftr->func))(p1);
#        else
            return (Func(ftr->func))(p1);
#        endif
      }
};

// SxFuncTranslator<Func,RES,P1,P2>
template<class Func, class RES, class P1, class P2>
class SxFuncTranslator<Func,RES,P1,P2,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor< RES,P1,P2>
{
   public:
      SxFuncTranslator (Func func_, const char *tag_ = NULL) 
         : SxFunctor<RES,P1,P2> (metaFunc,
                                (new SxBoundStorage(
                                (SxFunctorStorage::Func *)const_cast<void *>(
                                (const void *)func_))), tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, P1 p1, P2 p2)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*>(functor_);
         SX_CHECK (ftr->func);
#        ifdef WIN32
            return (reinterpret_cast<RES (*)(P1,P2)>(ftr->func))(p1,p2);
#        else
            return (Func (ftr->func))(p1, p2);
#        endif
      }
};

// SxFuncTranslator<Func,RES,P1,P2,P3>
template<class Func, class RES, class P1, class P2, class P3>
class SxFuncTranslator<Func,RES,P1,P2,P3,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3>
{
   public:
      SxFuncTranslator (Func func_, const char *tag_ = NULL) 
         : SxFunctor<RES,P1,P2,P3> (metaFunc,
                                   (new SxBoundStorage(
                                   (SxFunctorStorage::Func *)const_cast<void *>(
                                   (const void *)func_))), tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, P1 p1, P2 p2, P3 p3)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*>(functor_);
         SX_CHECK (ftr->func);
#        ifdef WIN32
            return (reinterpret_cast<RES (*)(P1,P2,P3)>(ftr->func))(p1,p2,p3);
#        else
            return (Func (ftr->func))(p1, p2, p3);
#        endif
      }
};

// SxFuncTranslator<Func,RES,P1,P2,P3,P4>
template<class Func, class RES, class P1, class P2, class P3, class P4>
class SxFuncTranslator<Func,RES,P1,P2,P3,P4,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4>
{
   public:
      SxFuncTranslator (Func func_, const char *tag_ = NULL) 
         : SxFunctor<RES,P1,P2,P3,P4> (metaFunc,
                                      (new SxBoundStorage(
                                      (SxFunctorStorage::Func *)const_cast<void *>(
                                      (const void *)func_))), tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*>(functor_);
         SX_CHECK (ftr->func);
#        ifdef WIN32
            return (reinterpret_cast<RES (*)(P1,P2,P3,P4)>(ftr->func))(p1,p2,p3,p4);
#        else
            return (Func (ftr->func))(p1, p2, p3, p4);
#        endif
      }
};

// SxFuncTranslator<Func,RES,P1,P2,P3,P4,P5>
template<class Func, class RES, class P1, class P2, class P3, class P4, class P5>
class SxFuncTranslator<Func,RES,P1,P2,P3,P4,P5,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4,P5>
{
   public:
      SxFuncTranslator (Func func_, const char *tag_ = NULL) 
         : SxFunctor<RES,P1,P2,P3,P4,P5> (metaFunc,
                                         (new SxBoundStorage(
                                         (SxFunctorStorage::Func *)const_cast<void *>(
                                         (const void *)func_))), tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4, P5 p5)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*>(functor_);
         SX_CHECK (ftr->func);
#        ifdef WIN32
            return (reinterpret_cast<RES (*)(P1,P2,P3,P4,P5)>(ftr->func))(p1,p2,p3,p4,p5);
#        else
            return (Func (ftr->func))(p1, p2, p3, p4, p5);
#        endif
      }
};

// SxFuncTranslator<Func,RES,P1,P2,P3,P4,P5,P6>
template<class Func, class RES, class P1, class P2, class P3, class P4, class P5, class P6>
class SxFuncTranslator<Func,RES,P1,P2,P3,P4,P5,P6,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4,P5,P6>
{
   public:
      SxFuncTranslator (Func func_, const char *tag_ = NULL) 
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (metaFunc,
                                            (new SxBoundStorage(
                                            (SxFunctorStorage::Func *)const_cast<void *>(
                                            (const void *)func_))), tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*>(functor_);
         SX_CHECK (ftr->func);
#        ifdef WIN32
            return (reinterpret_cast<RES (*)(P1,P2,P3,P4,P5,P6)>(ftr->func))(p1,p2,p3,p4,p5,p6);
#        else
            return (Func (ftr->func))(p1, p2, p3, p4, p5, p6);
#        endif
      }
};



// ---------------------------------------------------------------------------

/*
 SxIsCallable allows to check at compile time whether the type T has
 an public overloaded operator() function with specific signature.
 */
struct SxIsCallable {

   template <typename T, typename RES>
   class ZeroArgs {
      template <typename U, RES (U::*)()> class Check;
      template <typename U, RES (U::*)()const> class ConstCheck;
      template <typename U> static char func(Check<U, &U::operator()> *);
      template <typename U> static char func(ConstCheck<U, &U::operator()> *);
      template <typename U> static int func(...);
      public:
         static const bool value = sizeof(func<T>(0)) == sizeof(char);
   };

   template <typename T, typename RES, typename P1>
   class OneArgs {
      template <typename U, RES (U::*)(P1)> class Check;
      template <typename U, RES (U::*)(P1)const> class ConstCheck;
      template <typename U> static char func(Check<U, &U::operator()> *);
      template <typename U> static char func(ConstCheck<U, &U::operator()> *);
      template <typename U> static int func(...);
      public:
         static const bool value = sizeof(func<T>(0)) == sizeof(char);
   };

   template <typename T, typename RES, typename P1, typename P2>
   class TwoArgs {
      template <typename U, RES (U::*)(P1, P2)> class Check;
      template <typename U, RES (U::*)(P1, P2)const> class ConstCheck;
      template <typename U> static char func(Check<U, &U::operator()> *);
      template <typename U> static char func(ConstCheck<U, &U::operator()> *);
      template <typename U> static int func(...);
      public:
         static const bool value = sizeof(func<T>(0)) == sizeof(char);
   };

   template <typename T, typename RES, typename P1, typename P2, typename P3>
   class ThreeArgs {
      template <typename U, RES (U::*)(P1, P2, P3)> class Check;
      template <typename U, RES (U::*)(P1, P2, P3)const> class ConstCheck;
      template <typename U> static char func(Check<U, &U::operator()> *);
      template <typename U> static char func(ConstCheck<U, &U::operator()> *);
      template <typename U> static int func(...);
      public:
         static const bool value = sizeof(func<T>(0)) == sizeof(char);
   };

   template <typename T, typename RES, typename P1, typename P2, typename P3,
             typename P4>
   class FourArgs {
      template <typename U, RES (U::*)(P1, P2, P3, P4)> class Check;
      template <typename U, RES (U::*)(P1, P2, P3, P4)const> class ConstCheck;
      template <typename U> static char func(Check<U, &U::operator()> *);
      template <typename U> static char func(ConstCheck<U, &U::operator()> *);
      template <typename U> static int func(...);
      public:
         static const bool value = sizeof(func<T>(0)) == sizeof(char);
   };

   template <typename T, typename RES, typename P1, typename P2, typename P3,
             typename P4, typename P5>
   class FiveArgs {
      template <typename U, RES (U::*)(P1, P2, P3, P4, P5)> class Check;
      template <typename U, RES (U::*)(P1, P2, P3, P4, P5)const> class ConstCheck;
      template <typename U> static char func(Check<U, &U::operator()> *);
      template <typename U> static char func(ConstCheck<U, &U::operator()> *);
      template <typename U> static int func(...);
      public:
         static const bool value = sizeof(func<T>(0)) == sizeof(char);
   };

   template <typename T, typename RES, typename P1, typename P2, typename P3,
             typename P4, typename P5, typename P6>
   class SixArgs {
      template <typename U, RES (U::*)(P1, P2, P3, P4, P5, P6)> class Check;
      template <typename U, RES (U::*)(P1, P2, P3, P4, P5, P6)const> class ConstCheck;
      template <typename U> static char func(Check<U, &U::operator()> *);
      template <typename U> static char func(ConstCheck<U, &U::operator()> *);
      template <typename U> static int func(...);
      public:
         static const bool value = sizeof(func<T>(0)) == sizeof(char);
   };
};


/*
 SxCallWrapper allows compile time conditional call to an object's operator()
 function if the condition is set to true. If condition is false then function
 is not called.
 */
template<bool Condition>
class SxCallWrapper {
   public:
      template<typename T, typename RES = void, typename... Args>
      static void call (T obj, Args... args) {
         SX_UNUSED(obj);
         SX_UNUSED(args...);
         SX_EXIT;
      }
};

template<>
class SxCallWrapper<true> {
   public:
      template<typename T, typename RES = void, typename... Args>
      static RES call (T obj, Args... args) {
         return obj(args...);
      }
};

// ------------------------------------------------------------------------

template<class Callee, class RES, 
         class P1 = SxNull, class P2 = SxNull, class P3 = SxNull,
         class P4 = SxNull, class P5 = SxNull, class P6 = SxNull,
         class Terminator = SxNull
        >
class SxLambdaTranslator : public SxFunctor<RES,P1,P2,P3,P4,P5,P6,Terminator>
{ /* empty */ };


template<class Callee, class RES>
class SxLambdaTranslator<Callee,RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES>
{
   public:
      SxLambdaTranslator (Callee &callee_, const char *tag_ = NULL)
         : SxFunctor<RES> (metaFunc,
                           (new SxLambdaStorage<Callee>(callee_,
                            SxFunctorStorage::MemFunc (nullptr))), tag_)
      { }

      static RES metaFunc (const SxFunctorStorage *functor_)
      {
         SX_CHECK (functor_);
         const SxLambdaStorage<Callee> *ftr = dynamic_cast<const SxLambdaStorage<Callee>*>(functor_);
         SX_CHECK (ftr);
         //Callee *callee = (Callee ftr->callee;
         static const bool ZeroArgs = SxIsCallable::ZeroArgs<Callee, RES>::value;
         if (ZeroArgs) {
            return SxCallWrapper<ZeroArgs>::template
                   call<Callee, RES>(ftr->callee);
         } else { SX_EXIT; }
      }
};

template<class Callee, class RES,
         class P1>
class SxLambdaTranslator<Callee,RES,P1,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1>
{
   public:
      SxLambdaTranslator (Callee &callee_, const char *tag_ = NULL)
         : SxFunctor<RES,P1> (metaFunc,
                              (new SxLambdaStorage<Callee>(callee_,
                               SxFunctorStorage::MemFunc (nullptr))), tag_)
      { }

      static RES metaFunc (const SxFunctorStorage *functor_, P1 p1)
      {
         SX_CHECK (functor_);
         const SxLambdaStorage<Callee> *ftr = dynamic_cast<const SxLambdaStorage<Callee>*>(functor_);
         SX_CHECK (ftr);
         //Callee *callee = (Callee *)functor_.callee;
         static const bool OneArgs = SxIsCallable::OneArgs<Callee, RES,
                                                         P1>::value;
         if (OneArgs) {
            return SxCallWrapper<OneArgs>::template
                   call<Callee, RES, P1>(ftr->callee, p1);
         } else { SX_EXIT; }
      }
};

template<class Callee, class RES, 
         class P1, class P2>
class SxLambdaTranslator<Callee,RES,P1,P2,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2>
{
   public:
      SxLambdaTranslator (Callee &callee_, const char *tag_ = NULL)
         : SxFunctor<RES,P1,P2> (metaFunc,
                                 (new SxLambdaStorage<Callee>(callee_,
                                  SxFunctorStorage::MemFunc (nullptr))), tag_)
      { }

      static RES metaFunc (const SxFunctorStorage *functor_, P1 p1, P2 p2)
      {
         SX_CHECK (functor_);
         const SxLambdaStorage<Callee> *ftr = dynamic_cast<const SxLambdaStorage<Callee>*>(functor_);
         SX_CHECK (ftr);
         //Callee *callee = (Callee *)functor_.callee;
         static const bool TwoArgs = SxIsCallable::TwoArgs<Callee, RES,
                                                         P1, P2>::value;
         if (TwoArgs) {
            return SxCallWrapper<TwoArgs>::template
                   call<Callee, RES, P1, P2>(ftr->callee, p1, p2);
         } else { SX_EXIT; }
      }
};

template<class Callee, class RES, 
         class P1, class P2, class P3>
class SxLambdaTranslator<Callee,RES,P1,P2,P3,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3>
{
   public:
      SxLambdaTranslator (Callee &callee_, const char *tag_ = NULL)
         : SxFunctor<RES,P1,P2,P3> (metaFunc,
                                    (new SxLambdaStorage<Callee>(callee_,
                                     SxFunctorStorage::MemFunc (nullptr))), tag_)
      { }

      static RES metaFunc (const SxFunctorStorage *functor_,
                           P1 p1, P2 p2, P3 p3)
      {
         SX_CHECK (functor_);
         const SxLambdaStorage<Callee> *ftr = dynamic_cast<const SxLambdaStorage<Callee>*>(functor_);
         SX_CHECK (ftr);
         //Callee *callee = (Callee *)functor_.callee;
         static const bool ThreeArgs = SxIsCallable::ThreeArgs<Callee, RES, P1, P2,
                                                    P3>::value;
         if (ThreeArgs) {
            return SxCallWrapper<ThreeArgs>::template
                   call<Callee, RES, P1, P2, P3>(ftr->callee, p1, p2, p3);
         } else { SX_EXIT; }
      }
};

template<class Callee, class RES, 
         class P1, class P2, class P3, class P4>
class SxLambdaTranslator<Callee,RES,P1,P2,P3,P4,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4>
{
   public:
      SxLambdaTranslator (Callee &callee_, const char *tag_ = NULL)
         : SxFunctor<RES,P1,P2,P3,P4> (metaFunc,
                                       (new SxLambdaStorage<Callee>(callee_,
                                        SxFunctorStorage::MemFunc (nullptr))), tag_)
      { }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4)
      {
         SX_CHECK (functor_);
         const SxLambdaStorage<Callee> *ftr = dynamic_cast<const SxLambdaStorage<Callee>*>(functor_);
         SX_CHECK (ftr);
         //Callee *callee = (Callee *)functor_.callee;
         static const bool FourArgs = SxIsCallable::FourArgs<Callee, RES, P1,
                                                           P2, P3, P4>::value;
         if (FourArgs) {
            return SxCallWrapper<FourArgs>::template
                   call<Callee, RES, P1, P2, P3, P4>(ftr->callee, p1, p2,
                                                     p3, p4);
         } else { SX_EXIT; }
      }
};

template<class Callee, class RES, 
         class P1, class P2, class P3, class P4, class P5>
class SxLambdaTranslator<Callee,RES,P1,P2,P3,P4,P5,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4,P5>
{
   public:
      SxLambdaTranslator (Callee &callee_, const char *tag_ = NULL)
         : SxFunctor<RES,P1,P2,P3,P4,P5> (metaFunc,
                                          (new SxLambdaStorage<Callee>(callee_,
                                           SxFunctorStorage::MemFunc (nullptr))), tag_)
      { }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4, P5 p5)
      {
         SX_CHECK (functor_);
         const SxLambdaStorage<Callee> *ftr = dynamic_cast<const SxLambdaStorage<Callee>*>(functor_);
         SX_CHECK (ftr);
         //Callee *callee = (Callee *)functor_.callee;
         static const bool FiveArgs = SxIsCallable::FiveArgs<Callee, RES, P1, P2,
                                                           P3, P4, P5>::value;
         if (FiveArgs) {
            return SxCallWrapper<FiveArgs>::template
                   call<Callee, RES, P1, P2, P3, P4, P5>(ftr->callee,
                                                         p1, p2, p3,
                                                         p4, p5);
         } else { SX_EXIT; }
      }
};

template<class Callee, class RES, 
         class P1, class P2, class P3, class P4, class P5, class P6>
class SxLambdaTranslator<Callee,RES,P1,P2,P3,P4,P5,P6,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4,P5,P6>
{
   public:
      SxLambdaTranslator (Callee &callee_, const char *tag_ = NULL)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (metaFunc,
                                             (new SxLambdaStorage<Callee>(callee_,
                                              SxFunctorStorage::MemFunc (nullptr))), tag_)
           { }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6)
      {
         SX_CHECK (functor_);
         const SxLambdaStorage<Callee> *ftr = dynamic_cast<const SxLambdaStorage<Callee>*>(functor_);
         SX_CHECK (ftr);
         //Callee *callee = (Callee *)functor_.callee;
         static const bool SixArgs = SxIsCallable::SixArgs<Callee, RES, P1, P2,
                                                         P3, P4, P5, P6>::value;
         if (SixArgs) {
            return SxCallWrapper<SixArgs>::template
                   call<Callee, RES, P1, P2, P3, P4, P5, P6>(ftr->callee,
                                                             p1, p2, p3,
                                                             p4, p5, p6);
         } else { SX_EXIT; }
      }
};

// ------------------------------------------------------------------------

// translate member function pointers to functors
// This implementation allows to slightly modify the signature on the fly
// The bound pointer's calling signature is RES func(P1, ... P6)
// The member function's calling signature is TRES A::func (TP1, ..., TP6)
// i.e., P1 must be convertible to TP1, P2 to TP2, ...
//       and TRES must be convertible to RES
// The conversion and possible compile time errors occur in metaFunc.
// SxCBoundPtr takes care of getting the number of parameters right,
// and puts the member function's signature into SxMemTranslator's
// MemberFunc template parameter

// SxMemTranslator<Callee, MemberFunc, RES, ...>
template<class Callee, class MemberFunc, class RES, 
         class P1 = SxNull, class P2 = SxNull, class P3 = SxNull,
         class P4 = SxNull, class P5 = SxNull, class P6 = SxNull,
         class Terminator = SxNull
        >
class SxMemTranslator : public SxFunctor<RES,P1,P2,P3,P4,P5,P6,Terminator>
{ /* empty */ };


// SxMemTranslator<Callee, MemberFunc, RES>
template<class Callee, class MemberFunc, class RES>
class SxMemTranslator<Callee,MemberFunc,RES,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES>
{
   public:
      SxMemTranslator (Callee &callee_, const MemberFunc &memFunc_, const char *tag_=NULL)
         : SxFunctor<RES> (metaFunc, (new SxBoundStorage(&callee_,
                           reinterpret_cast<SxFunctorStorage::MemFunc>(memFunc_))),
                           tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*> (functor_);
         SX_CHECK (ftr);
         Callee *callee = (Callee *)ftr->callee;
         MemberFunc memFunc = reinterpret_cast<MemberFunc>(ftr->memFunc);
         return (callee->*memFunc)();
      }
};

// SxMemTranslator<Callee, MemberFunc, RES, P1>
template<class Callee, class MemberFunc, class RES, 
         class P1>
class SxMemTranslator<Callee,MemberFunc,RES,P1,SxNull,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1>
{
   public:
      SxMemTranslator (Callee &callee_, const MemberFunc &memFunc_, const char *tag_=NULL)
         : SxFunctor<RES,P1> (metaFunc, (new SxBoundStorage(&callee_,
                           reinterpret_cast<SxFunctorStorage::MemFunc>(memFunc_))),
                           tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, P1 p1)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*> (functor_);
         SX_CHECK (ftr);
         Callee *callee = (Callee *)(ftr->callee);
         MemberFunc memFunc = reinterpret_cast<MemberFunc>(ftr->memFunc);
         return (callee->*memFunc)(p1);
      }
};

// SxMemTranslator<Callee, MemberFunc, RES, P1, P2>
template<class Callee, class MemberFunc, class RES, 
         class P1, class P2>
class SxMemTranslator<Callee,MemberFunc,RES,P1,P2,SxNull,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2>
{
   public:
      SxMemTranslator (Callee &callee_, const MemberFunc &memFunc_, const char *tag_=NULL)
         : SxFunctor<RES,P1,P2> (metaFunc, (new SxBoundStorage(&callee_,
                           reinterpret_cast<SxFunctorStorage::MemFunc>(memFunc_))),
                           tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, P1 p1, P2 p2)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*> (functor_);
         SX_CHECK (ftr);
         Callee *callee = (Callee *)(ftr->callee);
         MemberFunc memFunc = reinterpret_cast<MemberFunc>(ftr->memFunc);
         return (callee->*memFunc)(p1, p2);
      }
};

// SxMemTranslator<Callee, MemberFunc, RES, P1, P2, P3>
template<class Callee, class MemberFunc, class RES, 
         class P1, class P2, class P3>
class SxMemTranslator<Callee,MemberFunc,RES,P1,P2,P3,SxNull,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3>
{
   public:
      SxMemTranslator (Callee &callee_, const MemberFunc &memFunc_, const char *tag_=NULL)
         : SxFunctor<RES,P1,P2,P3> (metaFunc, (new SxBoundStorage(&callee_,
                           reinterpret_cast<SxFunctorStorage::MemFunc>(memFunc_))),
                           tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*> (functor_);
         SX_CHECK (ftr);
         Callee *callee = (Callee *)(ftr->callee);
         MemberFunc memFunc = reinterpret_cast<MemberFunc>(ftr->memFunc);
         return (callee->*memFunc)(p1, p2, p3);
      }
};

// SxMemTranslator<Callee, MemberFunc, RES, P1, P2, P3, P4>
template<class Callee, class MemberFunc, class RES, 
         class P1, class P2, class P3, class P4>
class SxMemTranslator<Callee,MemberFunc,RES,P1,P2,P3,P4,SxNull,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4>
{
   public:
      SxMemTranslator (Callee &callee_, const MemberFunc &memFunc_, const char *tag_=NULL)
         : SxFunctor<RES,P1,P2,P3,P4> (metaFunc, (new SxBoundStorage(&callee_,
                           reinterpret_cast<SxFunctorStorage::MemFunc>(memFunc_))),
                           tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*> (functor_);
         SX_CHECK (ftr);
         Callee *callee = (Callee *)(ftr->callee);
         MemberFunc memFunc = reinterpret_cast<MemberFunc>(ftr->memFunc);
         return (callee->*memFunc)(p1, p2, p3, p4);
      }
};

// SxMemTranslator<Callee, MemberFunc, RES, P1, P2, P3, P4, P5>
template<class Callee, class MemberFunc, class RES, 
         class P1, class P2, class P3, class P4, class P5>
class SxMemTranslator<Callee,MemberFunc,RES,P1,P2,P3,P4,P5,SxNull,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4,P5>
{
   public:
      SxMemTranslator (Callee &callee_, const MemberFunc &memFunc_, const char *tag_=NULL)
         : SxFunctor<RES,P1,P2,P3,P4,P5> (metaFunc, (new SxBoundStorage(&callee_,
                           reinterpret_cast<SxFunctorStorage::MemFunc>(memFunc_))),
                           tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4, P5 p5)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*> (functor_);
         SX_CHECK (ftr);
         Callee *callee = (Callee *)(ftr->callee);
         MemberFunc memFunc = reinterpret_cast<MemberFunc>(ftr->memFunc);
         return (callee->*memFunc)(p1, p2, p3, p4, p5);
      }
};

// SxMemTranslator<Callee, MemberFunc, RES, P1, P2, P3, P4, P5, P6>
template<class Callee, class MemberFunc, class RES, 
         class P1, class P2, class P3, class P4, class P5, class P6>
class SxMemTranslator<Callee,MemberFunc,RES,P1,P2,P3,P4,P5,P6,SxNull>
   : public SxFunctor<RES,P1,P2,P3,P4,P5,P6>
{
   public:
      SxMemTranslator (Callee &callee_, const MemberFunc &memFunc_, const char *tag_=NULL)
         : SxFunctor<RES,P1,P2,P3,P4,P5,P6> (metaFunc, (new SxBoundStorage(&callee_,
                           reinterpret_cast<SxFunctorStorage::MemFunc>(memFunc_))),
                           tag_)
      { /* empty */ }

      static RES metaFunc (const SxFunctorStorage *functor_, 
                           P1 p1, P2 p2, P3 p3, P4 p4, P5 p5, P6 p6)
      {
         SX_CHECK (functor_);
         const SxBoundStorage *ftr = dynamic_cast<const SxBoundStorage*> (functor_);
         SX_CHECK (ftr);
         Callee *callee = (Callee *)(ftr->callee);
         MemberFunc memFunc = reinterpret_cast<MemberFunc>(ftr->memFunc);
         return (callee->*memFunc)(p1, p2, p3, p4, p5, p6);
      }
};

#endif /* _SX_FUNCTOR_H_ */
