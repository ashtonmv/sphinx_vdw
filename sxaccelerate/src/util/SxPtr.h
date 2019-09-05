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

#ifndef _SX_AUTO_POINTER_H_
#define _SX_AUTO_POINTER_H_

#include <SxError.h>
#include <SxSignals.h>
#include <atomic>

template<class T> class SxThis;

template <class T>
class SxConstPtr;

/** \brief Self managing pointer handler


    \b SxPtr = SPHInX Automatic Pointer Handler

    This class allows the usage of the C++ new statement for allocation
    without worrying about deallocation of the memory.
    By using the C++ new statement memory is allocated. Of course, at some
    time it also has to be deleted. Otherwise one would cause a memory leak.
    Unfortunately, developers tend to forget the delete call. Hence, the
    default usage of C++ new is rather dangerous.
    This class takes over the responsibility of deleting the memory.

    \par Example:

    Consider the conventional C++ code
\code
   MyClass *myClass = new MyClass;
   ...
    myClass->memberFunction();
   *myClass = 1;
   ...
   delete myClass;
\endcode

     Using the autopointer it would read
\code
   SxPtr<MyClass> myClass = SxPtr<MyClass>::create ();
   ...
    myClass->memberFunction ();
   *myClass = 1;
   ...
   // no deletion required.
\endcode

   Thus, the usage is exactly the same as for conventional pointers. Only,
   the final deletion is automatically performed.

    \par Misusage 1: Do not use it for arrays or vectors

    Please do not use this class for allocating arrays or vectors.
    This is handled by the SxArray or SxVector classes!

    \par Misusage 2: For constness use SxConstPtr

    Instead of SxPtr<const T> use always SxConstPtr<T>.

    \sa     SxConstPtr
    \author Sixten Boeck, boeck@mpie.de */
template<class T>
class SxPtr
{
   public:

      /** \brief Default constructor

          No memory is allocated by this default constructur. */
      SxPtr ();
      /** \brief Copy constructor

          The copy constructor copies the pointer (not its value) of the
          input argument. As usual in reference counting, the source's
          reference counter is increased afterwards. */
      template<class B>
         SxPtr (const SxPtr<B> &);
         SxPtr (const SxPtr<T> &);

      /** \brief Destructor

          This destructor decreases the reference counter by one. If there
          is no further external object pointing towards it (referenc counter
          is zero) both the pointer and the reference counter are
          deleted. */
     ~SxPtr ();

     /** @{

         \brief Counterpart to the conventional new statement

         When creating a pointer one must nut use the conventional C++
         \b new statement. Instead, one of the create functions can be
         used

         \par Example:
\code
   // int *i = new int;
   SxPtr<int> i = SxPtr<int>::create ();

   // MyClass *obj = new MyClass (1, 2., PI);
   SxPtr<MyClass> obj = SxPtr<MyClass>::create (1, 2., PI);
\endcode
     */
     static SxPtr<T> create ();
     template<class A1>
     static SxPtr<T> create (const A1 &);
     template<class A1, class A2>
     static SxPtr<T> create (const A1 &, const A2 &);
     template<class A1, class A2, class A3>
     static SxPtr<T> create (const A1 &, const A2 &, const A3 &);
     template<class A1, class A2, class A3, class A4>
     static SxPtr<T> create (const A1 &, const A2 &, const A3 &,
                                     const A4 &);
     template<class A1, class A2, class A3, class A4, class A5>
     static SxPtr<T> create (const A1 &, const A2 &, const A3 &,
                                     const A4 &, const A5 &);
     template<class A1, class A2, class A3, class A4, class A5, class A6>
     static SxPtr<T> create (const A1 &, const A2 &, const A3 &,
                                     const A4 &, const A5 &, const A6 &);
     /** @} */


     /** @{
         \brief Take over the management of C-like pointers.

         \b Attention! Use this function with great care! SxPtrs
         should be initialized with the create functions. Under normal
         circumstances there is no reason to initialize SxPtrs
         directly from C-like pointers.
         However, sometimes one has to deal with external libraries which
         return objects alloced using normal memory management. In these
         cases it might be useful to take over the memory control using 
         this function.
         \sa create */
     SxPtr<T> &initFromCPointer (T *);
     /// @}


     /** @{
         \brief Return the actual pointer (USE WITH CARE)

         \b Attention! Use this function with great care! SxPtrs
         should not present their pointer members. Particulary, as soon
         as the returned vector is deleted or otherwisely referenced the
         memory management will break down!
         However, they are situations when dealing with the actual pointer
         is useful. */
     T *getPtr () const;
     const T *getConstPtr () const;
     /** @} */

     /** @{ 
         \brief Assignment operator

         Copy the pointer (not the value) from another auto pointer object. */
     template<class B>
           SxPtr<T> &operator= (const SxPtr<B> &);
           SxPtr<T> &operator= (const SxPtr<T> &);
     /** @} */

     /** @{
         \brief Compare SxPtrs

         The comparison is particluarily useful to check if the pointer
         is set at all.
\code
   if (autoPtr == SxPtr<SomeType> ())  {
      // auto pointer is empty
   }  else  {
      // auto pointer is set
   }
\endcode
   */
     bool operator== (const SxPtr<T> &) const;
     bool operator!= (const SxPtr<T> &) const;
     /** @} */

     /** @{
         \brief Dereferencing operator

         Get the pointer's value. Using as in conventional C
\code
   SxPtr<int> iPtr = new int;
   *i = 1;
   printf ("%d\n", *i);
\endcode
          */
           T &operator* () const;

     /** @{
          \brief Dereferencing operator

          Call a functions member. Using as in conventional C++
\code
   SxPtr<MyClass> myClass = new MyClass;
   myClass->foo ();
\endcode
   */
           T *operator->() const;
     /** @} */

           template<class B>
           explicit operator B* ();


// protected:


     /** \brief Conventional pointer

         This member variable is the actual pointer. */
     T   *ptr;
     bool selfPtr;

   public:
     /** \brief Reference counter

         As the name might imply it counts the number of references, i.e, 
         the number of extrenal objects also pointing to this ptr member.
         Only if the reference counter is 0 the destructor destroys the
         objects. */
     std::atomic<int> *refCounter;

   class Receiver
   {
      signals:

     /** \brief (FOR INTERNAL USAGE ONLY) A reference counted signal emitted when the object is deleted

         \b Warning

         Do no connect this signal using sxconnect. Use SxSigConnector::registerPtr instead.

         This signal is being emitted before the last SxPtr is destroyed
         and the object ptr T * is deleted. This signal is basically used
         for automatic unreferencing of the signal/slot mechanism in S/PHI/nX.
      */
         SxSignal<void *, const char *> SX_SIGNAL(sigDeleting);

      public slots:
         void slotDeleting (void *r, const char *dbgInfo);
   };

   Receiver *receiver;

   protected:

     /** \brief Increase the reference counter by one */
     void ref ();

     /** \brief Decrease the reference counter by one */
     void unref ();

      /** \brief Initialize the auto pointer from a C++ pointer.

          It constructs an auto pointer directly from an conventional
          C++ pointer. This constructor has been protected in order to
          forbid the direct usage of the \b new statement. Instead
          the create member function should be applied

          \par Example
\code
   SxPtr<MyClass> myClass = SxPtr<MyClass>::create ();
\endcode
      */
      SxPtr (T *, bool selfPtr_=false);

      /** \brief create 'this' pointer

          Create a 'this' like auto pointer. In order to access a 'this'
          autopointer simply derive the class from SxThis<T>. */
      static SxPtr<T> self (T *);
      friend class SxThis<T>;

   public:
      /// Type cast to bool (true if pointer is defined)
      operator bool () const { return ptr != NULL; }

      /// Type cast to SxConstPtr reference
      operator const SxConstPtr<T> & () const
      {
         return reinterpret_cast<const SxConstPtr<T> &>(*this);
      }

};


/** \brief const autopointer

    Use this class if the pointer type should be const. Before continuing
    please read first the documentation of class ::SxPtr.

    \par Example:
    Assume that A is the base class of B
\code
   // normal construction
   SxConstPtr<A> a1 = SxConstPtr<A>::create ();

   // down-cast to basis type
   SxConstPtr<A> a2 = SxConstPtr<B>::create ();

   // honor constness of A
   SxPtr<A>   aNonConst = SxPtr<A>::create ();
   SxConstPtr<A> aConst = aNonConst;  // type cast to const

   // honor constness of derived B
   SxPtr<B>   bNonConst = SxPtr<A>::create ();
   SxConstPtr<A> aConst = bNonConst;  // type cast to const
\endcode

    \sa    SxPtr
    \autor Sixten Boeck
 */
template<class T>
class SxConstPtr : public SxPtr<const T>
{
   public:

      /** \brief Default constructor */
      SxConstPtr ();

      /** \brief convert from SxPtr<B>::create ()

          This constructor is basically used in expressions similar to
\code
   SxConstPtr<T> t = SxPtr<B>::create ();
\endcode
      with B being derived from T
      */
      template<class B>
      SxConstPtr (const SxPtr<B> &);

      /** \brief convert from SxPtr<T>::create ()

          This constructor is basically used in expressions similar to
\code
    SxConstPtr<T> t = SxConstPtr<T>::create ();
\endcode
      */
      SxConstPtr (const SxPtr<const T> &);

      /** \brief convert from SxPtr<T>::create ()

          This constructor is basically used in expressions similar to
\code
   SxConstPtr<T> t = SxConstPtr<B>::create ();
\endcode
      with B being derived from T
      */

      template<class B>
      SxConstPtr (const SxPtr<const B> &);

      /** \brief convert from SxPtr<T>::create ()

          This constructor is basically used in expressions similar to
\code
   SxConstPtr<T> a;
   SxConstPtr<T> t(a);
\endcode
      */
      SxConstPtr (const SxConstPtr<T> &);

      /** \brief convert from SxPtr<T>::create ()

          This constructor is basically used in expressions similar to
\code
   SxConstPtr<T> b;
   SxConstPtr<T> t(b);
\endcode
      */
      template<class B>
      SxConstPtr (const SxConstPtr<B> &);

};


/** \brief Create SxPtr<T> from 'this'

    Sometimes a 'this'-like autopointer is handy. In order to create a 'this'
    auto pointer derive the class from SxThis and retrieve the corrresponding
    this autopointer using getThis()

\code
    class A : public SxThis<A>
    {
       public:
          void foo () {
             SxPtr<A> myThisPtr = getThis ();
          }
    };
\endcode

\author Sixten Boeck, boeck@mpie.de
*/
template<class T>
class SxThis
{
   public:

      SxThis ();
     ~SxThis ();

      /** \brief return a 'this' autopointer 
        */
      SxPtr<T> getThis () const;

   private:
      SxPtr<T>      self;
};



#include <SxPtr.hpp>

#endif /* _SX_AUTO_POINTER_H_ */
