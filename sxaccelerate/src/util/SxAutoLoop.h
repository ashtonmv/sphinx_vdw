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

#ifndef _SX_AUTO_LOOP_H_
#define _SX_AUTO_LOOP_H_

#include <SxUtil.h>
#include <SxError.h>

/** \brief Automatic loops

    This is a simple auxiliary class to enable automatic loops over
    containers. The container will provide the actual loop limit.

    The macro SX_LOOP then provides automatic loops from 0 to
    the limit implicitly provided in the loop body. In context other
    than accessing the autoloop-compatible containers, the loop counter
    should behave like a ssize_t.
    \code
    SX_LOOP(i)
    {
       container(i) = i;
    }
    \endcode
    \note SX_LOOP does not allow for empty containers, except if the
    limit is set explicitly:
    \code
    SX_LOOP(iAtom(nAtom)) // OK even if nAtom = 0
    {
       force(iAtom) += ...
    }
    Setting the limit explicitly should otherwise be avoided, unless
    the loop limits cannot be deduced.

    Autoloop support is available for selected containers. To add
    autoloop support for some function "get", returning "MyType",
    add a function
    \code
    MyType get(const SxAutoLoop &it)
    {
       it.setLimit (getSize ());
       return get(it.i);
    }
    \endcode
    to set the autoloop limit to "getSize ()".

    \author Christoph Freysoldt

  */
class SX_EXPORT_UTIL SxAutoLoop
{
   protected:
      mutable ssize_t limit;
   public:
      /// The actual loop iterator
      ssize_t i;

      /// Constructor
      SxAutoLoop ()
         : limit (-1), i(0)
      { }

      /// Initializer
      SxAutoLoop (int i_)
         : limit (-1), i(i_)
      { }

      /** \brief Change i initializer to limit initializer
          \note This construct allows to explicit the upper limit
          in SX_LOOP via
          \code
SX_LOOP(iDir(3))
          \endcode
          In the for initialization statement, the brackets will be seen
          as constructor and will set i=3.
          [We need the i constructor to allow for mixing of SxAutoLoop and
          integer variables in multi-index container classes.]
          In the loop condition (directly after the initialization), the
          brackets will be seen as operator(), where we can set limit=3
          and i=0 (as intended by the notation).
          Note that setting the limits explicitly will also restore the
          proper behavior for limit=0.
        */
      SxAutoLoop& operator() (ssize_t limitIn)
      { 
         if (limit == - 1 && limitIn == i)  {
            limit = limitIn;
            i = 0;
         }
         return *this; 
      }

      /// Initializer
      SxAutoLoop (long i_)
         : limit (-1), i(i_)
      { }

      /// Forbidden copy constructor
      // note: this is forbidden, because setting limits in copies
      // would not change the actual loop limit.
      SxAutoLoop (const SxAutoLoop &)
      {
         SX_EXIT;
      }

      /// Automatic cast to ssize_t for int-like uses of the iterator
      operator ssize_t () const
      {
         return i;
      }

      /// Set/check the limit by an auto-loop compatible container
      void setLimit (ssize_t limit_) const
      {
         SX_CHECK (limit_ > 0, limit_);
         if (limit != -1)  {
            SX_CHECK (limit == limit_, limit, limit_);
         } else {
            limit = limit_;
         }
      }

      /// loop iteration
      void operator++ ()
      {
         // limit must be set somewhere in the loop body
         SX_CHECK (limit > 0, limit);
         ++i;
      }

      /// loop checking
      bool more () const
      {
         return (limit < 0 || i < limit);
      }

      /// Do nothing template for setting the limit
      template <class T>
      static void setLimit (const T&, ssize_t) { }


};

/// Static version of setLimit
template<>
inline void SxAutoLoop::setLimit (const SxAutoLoop &i, ssize_t limit)
{
   i.setLimit (limit);
}

#define SX_LOOP(i) for (SxAutoLoop i; i.more (); ++i)
#define SX_LOOP2(i,j) SX_LOOP(i) SX_LOOP(j) 
#define SX_LOOP3(i,j,k) SX_LOOP(i) SX_LOOP(j) SX_LOOP(k)
#define SX_LOOP4(i,j,k,l) SX_LOOP(i) SX_LOOP(j) SX_LOOP(k) SX_LOOP(l)

#endif /* _SX_AUTO_LOOP_H_ */
