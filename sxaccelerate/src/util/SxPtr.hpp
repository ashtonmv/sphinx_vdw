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

#include <stdio.h>
#include <SxError.h>

template<class T>
SxPtr<T>::SxPtr ()
   : ptr(NULL),
     selfPtr(false),
     refCounter(NULL),
     receiver(NULL)
{
   // empty
}


template<class T>
SxPtr<T>::SxPtr (T *ptr_, bool selfPtr_)
   : ptr(ptr_),
     selfPtr(selfPtr_),
     refCounter(NULL),
     receiver(NULL)
{
   ref ();
}



template<class T>
template<class B>
SxPtr<T>::SxPtr (const SxPtr<B> &in)
   : ptr (dynamic_cast<T *>(in.ptr)),
     selfPtr(in.selfPtr),
     refCounter(in.refCounter),
     receiver((typename SxPtr<T>::Receiver *)in.receiver)
{
   SX_CHECK (in.ptr == NULL || dynamic_cast<T *>(in.ptr));
   SX_CHECK (in.refCounter);
   ref ();
}



template<class T>
SxPtr<T>::SxPtr (const SxPtr<T> &in)
   : ptr (in.ptr),
     selfPtr(in.selfPtr),
     refCounter(in.refCounter),
     receiver(in.receiver)
{
   if (ptr || refCounter)  {
      SX_CHECK (ptr);
      SX_CHECK (receiver);
      SX_CHECK (refCounter);
      ref ();
   }  else  {
      SX_CHECK (!ptr);
      SX_CHECK (!receiver);
      SX_CHECK (!refCounter);
   }
}



template<class T>
SxPtr<T>::~SxPtr ()
{
   unref ();
}


template<class T>
void SxPtr<T>::ref ()
{
   if (!refCounter)  {
      refCounter  = new std::atomic<int>{0};
      receiver    = new Receiver;
   }  else  {
      SX_CHECK (ptr); // do not reference count NULL pointers
      //(*refCounter)++;
      refCounter->fetch_add (1, std::memory_order_relaxed);
   }
}


template<class T>
void SxPtr<T>::unref ()
{
   if (!refCounter)  return;
   //SX_CHECK (ptr || (!ptr && *refCounter == 0), *refCounter);
   SX_CHECK (ptr || (!ptr && refCounter->load (std::memory_order_relaxed) == 0),
             refCounter->load (std::memory_order_relaxed));

   //(*refCounter)--;
   refCounter->fetch_sub (1, std::memory_order_relaxed);

   if (*refCounter < 0)  {
      receiver->sigDeleting.send (const_cast<void *>((const void *)ptr));
      if (!selfPtr)  delete ptr;
      delete refCounter;
      delete receiver;
      ptr         = NULL;
      refCounter  = NULL;
      receiver    = NULL;
   }
}


template<class T>
void SxPtr<T>::Receiver::slotDeleting (void *r, const char *dbgInfo)
{
   sigDeleting.slotDeregister (r, dbgInfo);
}


template<class T>
SxPtr<T> SxPtr<T>::create ()
{
   T *ptr = new T;
   SxPtr<T> res(ptr);
   return res;
}

template<class T>
template<class A1>
SxPtr<T> SxPtr<T>::create (const A1 &a1)
{
   return SxPtr<T> (new T (a1));
}


template<class T>
template<class A1, class A2>
SxPtr<T> SxPtr<T>::create (const A1 &a1, const A2 &a2)
{
   return SxPtr<T> (new T (a1, a2));
}


template<class T>
template<class A1, class A2, class A3>
SxPtr<T> SxPtr<T>::create (const A1 &a1, const A2 &a2,
                                           const A3 &a3)
{
   return SxPtr<T> (new T (a1, a2, a3));
}


template<class T>
template<class A1, class A2, class A3, class A4>
SxPtr<T> SxPtr<T>::create (const A1 &a1, const A2 &a2,
                                           const A3 &a3, const A4 &a4)
{
   return SxPtr<T> (new T (a1, a2, a3, a4));
}


template<class T>
template<class A1, class A2, class A3, class A4, class A5>
SxPtr<T> SxPtr<T>::create (const A1 &a1, const A2 &a2,
                                           const A3 &a3, const A4 &a4,
                                           const A5 &a5)
{
   return SxPtr<T> (new T (a1, a2, a3, a4, a5));
}


template<class T>
template<class A1, class A2, class A3, class A4, class A5, class A6>
SxPtr<T> SxPtr<T>::create (const A1 &a1, const A2 &a2,
                                           const A3 &a3, const A4 &a4,
                                           const A5 &a5, const A6 &a6)
{
   return SxPtr<T> (new T (a1, a2, a3, a4, a5, a6));
}


template<class T>
SxPtr<T> SxPtr<T>::self (T *ptr_)
{
   return SxPtr<T> (ptr_, true);
}


template<class T>
SxPtr<T> &SxPtr<T>::initFromCPointer (T *inPtr)
{
   SX_CHECK (inPtr);  // use SxPtr<> () to create a NULL pointer
   unref ();
   ptr = inPtr;
   selfPtr = false;
   ref ();
   return *this;
}

template<class T>
T *SxPtr<T>::getPtr () const
{
   return ptr;
}


template<class T>
const T *SxPtr<T>::getConstPtr () const
{
   return ptr;
}




template<class T>
template<class B>
SxPtr<T> &
SxPtr<T>::operator= (const SxPtr<B> &in)
{
   SX_CHECK (in.ptr == NULL || dynamic_cast<T *>(in.ptr));

   if ((const void *)(&in) == (void *)(this))  return *this;

   unref ();

   ptr         = dynamic_cast<T *>(in.ptr);
   selfPtr     = in.selfPtr;
   refCounter  = in.refCounter;
   receiver    = (typename SxPtr<T>::Receiver *)in.receiver;
   if (!ptr) return *this; // allow assigning NULL autopointers
   SX_CHECK (refCounter);

   ref ();


   return *this;
}


template<class T>
SxPtr<T> &
SxPtr<T>::operator= (const SxPtr<T> &in)
{
   if (&in == this)  return *this;

   unref ();

   ptr         = in.ptr;
   selfPtr     = in.selfPtr;
   refCounter  = in.refCounter;
   receiver    = in.receiver;
   if (!ptr) return *this; // allow assigning NULL autopointers
   SX_CHECK (refCounter);

   ref ();


   return *this;
}


template<class T>
bool SxPtr<T>::operator== (const SxPtr<T> &in) const
{
   return (ptr == in.ptr);
}


template<class T>
bool SxPtr<T>::operator!= (const SxPtr<T> &in) const
{
   return (ptr != in.ptr);
}



template<class T>
T &SxPtr<T>::operator* () const
{
   SX_CHECK (ptr);  // has create() been called alread?
   return *ptr;
}

template<class T>
T *SxPtr<T>::operator-> () const
{
   SX_CHECK (ptr);
   return ptr;
}

template<class T>
template <class B>
SxPtr<T>::operator B* ()
{
   SX_CHECK(dynamic_cast<B*>(ptr));
   return static_cast<B*> (ptr);
}

// --------------------------------------------------------------------------
template<class T>
SxConstPtr<T>::SxConstPtr ()
   : SxPtr<const T> ()
{
   // empty
}


// from SxConstPtr<T> t = SxPtr<B>::create ();
// with B being derived from T
template<class T>
template<class B>
SxConstPtr<T>::SxConstPtr (const SxPtr<B> &in)
   : SxPtr<const T> ()
{
   SX_CHECK (in.ptr == NULL || dynamic_cast<const T *>(in.ptr));
   this->ptr         = dynamic_cast<const T*>(in.ptr);
   this->selfPtr     = in.selfPtr;
   this->refCounter  = in.refCounter;
   this->receiver    = (typename SxConstPtr<T>::Receiver *)in.receiver;
   if (this->ptr) this->ref ();
}

// from SxConstPtr<T> t = SxConstPtr<T>::create ();
template<class T>
SxConstPtr<T>::SxConstPtr (const SxPtr<const T> &in)
   : SxPtr<const T> ()
{
   SX_CHECK (in.ptr == NULL || dynamic_cast<const T *>(in.ptr));
   this->ptr         = dynamic_cast<const T *>(in.ptr);
   this->selfPtr     = in.selfPtr;
   this->refCounter  = in.refCounter;
   this->receiver    = (typename SxConstPtr<T>::Receiver *)in.receiver;
   if (this->ptr) this->ref ();
}

// from SxConstPtr<T> t = SxConstPtr<B>::create ();
// with B being derived from T
template<class T>
template<class B>
SxConstPtr<T>::SxConstPtr (const SxPtr<const B> &in)
   : SxPtr<const T> ()
{
   SX_CHECK (in.ptr == NULL || dynamic_cast<const T *>(in.ptr));
   this->ptr         = dynamic_cast<const T *>(in.ptr);
   this->selfPtr     = in.selfPtr;
   this->refCounter  = in.refCounter;
   this->receiver    = (typename SxConstPtr<T>::Receiver *)in.receiver;
   if (this->ptr) this->ref ();
}

// SxConstPtr<T> a;
// SxConstPtr<T> t(a);
template<class T>
SxConstPtr<T>::SxConstPtr (const SxConstPtr<T> &in)
   : SxPtr<const T> ()
{
   this->ptr         = in.ptr,
   this->selfPtr     = in.selfPtr;
   this->refCounter  = in.refCounter;
   this->receiver    = (typename SxConstPtr<T>::Receiver *)in.receiver;
   if (this->ptr) this->ref ();
}

// SxConstPtr<B> b;
// SxConstPtr<T> t(b);
template<class T>
template<class B>
SxConstPtr<T>::SxConstPtr (const SxConstPtr<B> &in)
   : SxPtr<const T> ()
{
   SX_CHECK (in.ptr == NULL || dynamic_cast<const T *>(in.ptr));
   this->ptr         = dynamic_cast<const T *>(in.ptr),
   this->selfPtr     = in.selfPtr;
   this->refCounter  = in.refCounter;
   this->receiver    = (typename SxConstPtr<T>::Receiver *)in.receiver;
   if (this->ptr) this->ref ();
}


// --------------------------------------------------------------------------
template<class T>
SxThis<T>::SxThis () : self (SxPtr<T>::self ((T *)this)) 
{ 
   // empty
}

template<class T>
SxThis<T>::~SxThis ()
{
   // empty
}

template<class T>
SxPtr<T> SxThis<T>::getThis () const 
{
   return self;
}


