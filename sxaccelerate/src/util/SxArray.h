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
#ifndef _SX_ARRAY_H_
#define _SX_ARRAY_H_

#include <SxUtil.h>
#include <SxError.h>
#include <SxSelection.h>
#include <SxList.h>
#include <SxRange.h>
#include <SxMemConsumer.h>
#include <SxStack.h>
#include <SxIterator.h>
#include <SxContainer.h>
#include <SxAlg.h>
#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include <SxAutoLoop.h>
#include <SxPtr.h>
template<class T> void sort (ssize_t *, const T *, ssize_t);

/** \brief Array container for all data types

    \b SxArray = SPHInX Array Container Class

    This class can be used for not defragmented memory management for all 
    data types other than those defined in the \ref page_vector.
    Note, that in higher SPHInX classes the usage of the C++ statements
    \b new and \b delete are forbidden! Use this class instead, e.g.

    \b wrong:
       \verbatim
       int *intArray = new intArray[5];
       ...
       intArray[2] = ...;
       ...
       delete [] intArray;
    \endverbatim

    \b coorect:
    \verbatim
       #include <SxArray.h>
       ...
       SxArray<int> intArray(5);
       ...
       intArray(2) = ...;
    \endverbatim

    \par Using lists to initialize arrays:
    Sometimes the actual size of an array is not known at construction time,
    e.g., when data have to be taken from an input file. In this case one
    can use SxList objects for the intialization:

\code
#include <SxArray.h>
#include <SxList.h>
...

// --- read number of integers from a file into a list
SxList<int> dataList;
while (filePointer)  {
   dataList << ...;
}

// --- initialize array from list
SxArray<int> data = dataList;
...
\endcode

    \sa SxList
    \author  Sixten Boeck */
template<class T>
class SxArray : public SxMemConsumer,
                public SxThis<SxArray<T> >
{
   public:

      typedef SxArray<T> Container;

      // --- SxSelection:
      typedef T                TElem;
      typedef ssize_t          SelIdx;
      typedef SxList<SelIdx>   SelStorage;
      typedef Container        SelContainer;

      template<class SValue,class Container,class IT>
      class State
      {
         SX_ITERATOR_STATE
         // frient to State<const...>
         template<class CV,class CC,class CIT> friend class State;

         public:
            State (sx::Direction dir_=sx::Forward, Container *c=NULL,
                   ssize_t i=-1)
                   : dir(dir_), container(c), idx(i) { }

            State (const State &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), idx(in.idx) {
               SX_UNUSED (cmode_);
            }

            State (State &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : dir(in.dir), container(in.container), idx(in.idx)
            {
               SX_UNUSED (cmode_);
               in.dir = sx::Undefined; in.container = NULL; in.idx = -1;
            }
            // non-const to const cast
            template<class CV,class CC,class CIT>
            State (const State<CV, CC,CIT> &in)
               : dir(in.dir),container(in.container), idx(in.idx) { }
            ssize_t getIdx () const { return idx; }
            typename Container::SelIdx getSelIdx () const {
               SX_CHECK (container);
               return idx;
            }
         protected:
            Container *container;
            ssize_t idx;

            void copy (const IT &in) {
               dir = in.dir; container = in.container; idx = in.idx;
            }
            void move (IT &&in) {
               dir = in.dir; container = in.container; idx = in.idx;
               in.dir = sx::Undefined; in.container = NULL; in.idx = -1;
            }
            IT insertElem (ssize_t newPos, const SValue &elem) {
               SX_CHECK (container);
               return container->insert (newPos, elem);
            }
            IT prependElem (const SValue &elem) {
               SX_CHECK (container);
               return container->prepend (elem);
            }
            IT appendElem (const SValue &elem) {
               SX_CHECK (container);
               return container->append (elem);
            }
            // --- SX_ITERATOR / SX_CONST_ITERATOR callbacks
            void next () {
               SX_CHECK (container);
               ++idx;
               if (idx > container->getSize()-1)  idx = -1;
            }
            void prev () {
               SX_CHECK (container);
               --idx;
               if (idx < 0)  idx = -1;
            }
            bool valid () const { return idx != -1; }
            SValue &getRef () { SX_CHECK (container); return (*container)(idx); }
            SValue *getPtr () { SX_CHECK (container); return &(*container)(idx); }
            bool equal (const IT &in) const {
               return (idx == in.idx && container == in.container);
            }

      };

      class Iterator : public State<T,Container,Iterator>
      {
         SX_ITERATOR (T,Container,Iterator)
         public:
            Iterator (sx::Direction dir_=sx::Forward, Container *c=NULL,
                      ssize_t i=-1)
                      : State<T, Container,Iterator> (dir_, c, i)
                      { }

            Iterator (const Iterator &in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<T,Container,Iterator> (in, cmode_) { }

            Iterator (Iterator &&in, sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<T,Container,Iterator> (std::move(in), cmode_) { }
      };

      class ConstIterator : public State<const T,const Container,ConstIterator>
      {
         SX_CONST_ITERATOR (T,Container,ConstIterator)
         public:
            ConstIterator (sx::Direction dir_=sx::Forward,
                           const Container *c=NULL, ssize_t i=-1)
               : State<const T,const Container,ConstIterator> (dir_, c, i)
                 { }

            ConstIterator (const ConstIterator &in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const T,const Container,ConstIterator> (in, cmode_)
                 { }

            ConstIterator (ConstIterator &&in,
                           sx::ItCopyMode cmode_ = sx::CopyAll)
               : State<const T,const Container,ConstIterator>
                 (std::move(in), cmode_) { }

            // non-const to const cast
            ConstIterator (const Iterator &in)
               : State<const T,const Container,ConstIterator> (in)
                 { }
      };
      SX_CONTAINER(Container)


      /** array is the pointer to the begin of the data array. Do access data
          from that array use the dereferencing functions. */
      T *elements;
      /** The number of used elements in the memory chunk. */
      ssize_t size;
      /** Real number of preallocated elements in memory. */
      ssize_t sizeAlloc;
      /** The number of elements to allocate in one chunk. Variable. */
      ssize_t chunkSize;

      /** Standard constructor. */
      SxArray ();
      /** This constructor allocates a specified number of elements. The array
          is <em> not </em> being preset wit any value (such as 0.)! */
      explicit SxArray (ssize_t, ssize_t chunkSizeIn=1);
      /** Copy constructor. */
      SxArray (const SxArray<T> &, ssize_t chunkSizeIn=1);
      SxArray (SxArray<T> &&, ssize_t chunkSizeIn=1);


      /** Creates an array (unfragmented memory) from a SxList<T> object, which
          is fragmented.*/
      SxArray (const SxList<T> &, ssize_t chunkSizeIn=1);

      SxArray (const SxCList<T> &);

      /** Creates a new SxArray by copying a specified number of elements 
          from the provided input array. */
      SxArray (const T *, ssize_t, ssize_t chunkSizeIn=1);
      /** Creates a new SxArray by importing a stack.  */
      SxArray (const SxStack<T> &stack, ssize_t chunkSizeIn=1);
      template<long Begin,long End,long Step=1>
      SxArray (const SxRange<Begin,End,Step> &range);
      /** Creates a new SxArray from initializer_list  */
      SxArray (const std::initializer_list<T> &, ssize_t chunkSizeIn=1);

      /** Destructor. */
      ~SxArray ();

      // Iterators
      typename SxArray<T>::ConstIterator begin () const {
         return typename SxArray<T>::ConstIterator (sx::Forward, this,
                                                    getSize() == 0 ? -1 : 0);
      }
      typename SxArray<T>::Iterator begin () {
         return typename SxArray<T>::Iterator (sx::Forward, this, 
                                               getSize() == 0 ? -1 : 0);
      }
      typename SxArray<T>::Iterator
      begin (ssize_t idx_, sx::Direction dir_ = sx::Forward) {
         if (idx_ >=0 && idx_ < getSize())
            return typename SxArray<T>::Iterator (dir_, this,
                                                  idx_ );
         else
            return typename SxArray<T>::Iterator (dir_, this,
                                                  -1 );
      }
      typename SxArray<T>::ConstIterator
      begin (ssize_t idx_, sx::Direction dir_ = sx::Forward) const {
         if (idx_ >=0 && idx_ < getSize())
            return typename SxArray<T>::ConstIterator (dir_, this,
                                                  idx_ );
         else
            return typename SxArray<T>::ConstIterator (dir_, this,
                                                  -1 );
      }
      typename SxArray<T>::ConstIterator end () const {
         return typename SxArray<T>::ConstIterator (sx::Forward, this, -1);
      }
      typename SxArray<T>::Iterator end () {
         return typename SxArray<T>::Iterator (sx::Forward, this, -1);
      }
      typename SxArray<T>::ConstIterator fromLast () const {
         return typename SxArray<T>::ConstIterator (sx::Backward, this,
                                                    getSize() == 0 ? -1 :
                                                    (getSize()-1));
      }
      typename SxArray<T>::Iterator fromLast () {
         return typename SxArray<T>::Iterator (sx::Backward, this,
                                               getSize() == 0 ? -1 :
                                               (getSize()-1));
      }
      typename SxArray<T>::ConstIterator toFirst () const {
         return typename SxArray<T>::ConstIterator (sx::Backward, this, -1);
      }
      typename SxArray<T>::Iterator toFirst () {
         return typename SxArray<T>::Iterator (sx::Backward, this, -1);
      }


      /** Standard assignment operator. */
      SxArray<T> &operator= (const SxArray<T> &);
      /** Move assignment operator. */
      SxArray<T> &operator= (SxArray<T> &&);
      /** Assignment from SxList. */
      SxArray<T> &operator= (const SxList<T> &);
      /** \brief Assignment from SxStack
       */
      SxArray<T> &operator= (const SxStack<T> &);
      /** initializer_list assignment operator. */
      SxArray<T> &operator= (const std::initializer_list<T> &in);

      /** Returns true if arrays are equal. */
      bool       operator== (const SxArray<T> &) const;
      /** Returns true if arrays are not equal. */
      bool       operator!= (const SxArray<T> &) const;
      /** Standard dereferencing operator for non-const objectes. */
      T &operator() (ssize_t);
      /** Standard dereferencing operator for const objectes. */
      const T &operator() (ssize_t) const;
      /** Standard dereferencing operator for non-const objectes. */
      T &operator() (SxAutoLoop &);
      /** Standard dereferencing operator for const objectes. */
      const T &operator() (SxAutoLoop &) const;
      /** Dereferencing function. */
      const T at (ssize_t) const;

      inline Iterator      getIterator (const SelIdx &);
      inline ConstIterator getIterator (const SelIdx &) const;
      inline ConstIterator getConstIterator (const Iterator &) const;

      SxPtr<SxArray<T> > getContainer () const;

      /** Initialize the (previously allocated) array with a certain value. */
      void set (const T &);

      // -- Insertion
      Iterator prepend (const T &, ssize_t count=1);
      Iterator prepend (const SxArray<T> &);
      Iterator prepend (SxArray<T> &&);
      Iterator prepend (const std::initializer_list<T> &);

      Iterator insert (ssize_t newPos, const T &, ssize_t count=1);
      Iterator insert (ssize_t newPos, const SxArray<T> &);
      Iterator insert (ssize_t newPos, SxArray<T> &&);
      Iterator insert (ssize_t newPos, const std::initializer_list<T> &);

      Iterator append (const T &, ssize_t count=1);
      Iterator append (const SxArray<T> &);
      Iterator append (SxArray<T> &&);
      Iterator append (const std::initializer_list<T> &);

      /** Replace is a faster of:
          remove(from, nRemove) followed by insert(from, SxArray()). */
      void replace (ssize_t from, ssize_t nRemove, const SxArray<T> &);



      /** \brief Join List of Arrays into single Array. If the type is
                 trivially copyable 'memcpy' is used, if not all elements
                 are manually copied.

      \code
         SxArray<char> ar1 ({'a', 'b', 'c', 'd'});
         SxArray<char> ar2 ({'1', '2', '3', '4', '5'});
         SxArray<char> ar3 ({'x', 'y', 'z'});

         SxList<SxArray<char> > list;
         list << ar1 << ar2 << ar3;

         SxArray<char> result = SxArray<char>::join (list);
      \endcode

      \param list List of arrays to join

      \return Returns a new SxArray object.
      */
      static SxArray<T> join (const SxList<SxArray<T> > &list);

      // -- Search
      ssize_t findPos (const T &, ssize_t from=0) const;
      bool contains (const T &) const;

      // -- Deletion
      void remove (ssize_t idx, ssize_t count=1);
      void removeElement (const T &);
      void removeItem (typename SxArray<T>::Iterator *it);
      void removeFirst ();
      void removeLast ();
      void removeAll ();


      /** Set new chunk size.
          \param newChunkSize Specifies the amount of elements in one chunk.

          The chunk size is dynamic and can be modified anytime. The amount
          of allocated memory will be updated with first resize().
          \b example:
\code
   SxArray<int> a;

   a.setChunkSize (100); // allocate memory by 100 elements
   a.resize (1); // new [100], 1 element and other 99 preallocated
   a.resize (2); // 2 elements and other 98 preallocated
                 // a.getSize() == 2;

   a.setChunkSize (5); // allocate memory by 5 elements, still 100 allocated
   a.resize (3); // new [5], 3 elements and other 2 preallocated
   a.resize (0); // free preallocated memory, 0 preallocated
\endcode

          The chunk size can be also specified from the constructors. */
      void setChunkSize (ssize_t newChunkSize);

      /** Resize the array by allocating a new array.
         @param keep if set to false the new array after calling has undefined
                     elements.
                     if set to true and the new size is large than the old one
                     the trailing part is being set with 0.
                keep if set to true and the new size is smaller than the old one
                     the rest of the array will be truncated. */
      void resize (ssize_t, bool keep=false);
      /** Sort an array according to an index lookup table */
      void sortByIdx (const SxArray<ssize_t> &);
      /** Create an index array according to the provided values */
      SxArray<ssize_t> getSortIdx () const;

      /** Returns the size of the array */
      inline ssize_t getSize () const;

      /** Returns the number of allocated elements in memory.
          Should not be used except testing purposes. */
      inline ssize_t getSizeAlloc () const;

      /** Converts an array to a SxList object
        */
      inline SxList<T> toList () const;

   protected:
      /** Insert empty space to the array.
          \param idx Specifies the first index.
          \param count Specifies the number of elements to insert. */
      void insertSpace (ssize_t idx, ssize_t count);

      /** Returns the number of elements to allocate in memory.
          \param newSize Specifies the number of elements.
          \return the number of elements to allocate in memory
                  according to the current chunkSize.

          Should not be used except testing purposes. */
      ssize_t getSizeAlloc (ssize_t newSize) const;
};


template<class T>
SxArray<T>::SxArray ()
   : SxMemConsumer (),
     SxThis<SxArray<T> > ()
{
   elements  = NULL;
   size      = 0;
   sizeAlloc = 0;
   chunkSize = 1;
}

template<class T>
SxArray<T>::SxArray (ssize_t n, ssize_t chunkSizeIn)
   : SxMemConsumer (),
     SxThis<SxArray<T> > ()
{
   SX_CHECK (n >= 0, n);
   SX_CHECK (chunkSizeIn > 0, chunkSizeIn);

   chunkSize= chunkSizeIn;
   if (n == 0)  {
      sizeAlloc = 0;
      elements  = NULL;
   }  else  {
      sizeAlloc= getSizeAlloc (n);
      elements = new T [static_cast<size_t>(sizeAlloc)];
   }
   size     = n;

   TRACK_MALLOC (*this, 1);
}


template<class T>
SxArray<T>::SxArray (const SxArray<T> &in, ssize_t chunkSizeIn)
   : SxMemConsumer (in),
     SxThis<SxArray<T> > ()
{
   SX_CHECK (in.size >= 0, in.size);
   SX_CHECK (chunkSizeIn > 0, chunkSizeIn);

   chunkSize = chunkSizeIn;
   sizeAlloc = getSizeAlloc (in.size);

   if (in.size == 0)  {
      elements = NULL;
   }  else  {
      elements   = new T [static_cast<size_t>(sizeAlloc)];
      T *destPtr = elements;
      T *srcPtr  = in.elements;
      for (ssize_t i=0; i < in.size; i++)  {
//        elements[i]  =  in.elements[i];
         *destPtr++    = *srcPtr++;
      }
   }
   size  = in.size;
   TRACK_MALLOC (*this, 1);
}

template<class T>
SxArray<T>::SxArray (SxArray<T> &&in, ssize_t chunkSizeIn)
   : SxMemConsumer (in),
     SxThis<SxArray<T> > ()
{
   SX_CHECK (in.size >= 0, in.size);
   SX_CHECK (chunkSizeIn > 0, chunkSizeIn);

   chunkSize = chunkSizeIn;
   sizeAlloc = getSizeAlloc (in.size);

   if (in.size == 0)  {
      elements = NULL;
   }  else  {
      elements   = in.elements;
   }
   size  = in.size;

   in.elements = NULL;
   in.size = 0;
   in.sizeAlloc = 0;
   in.chunkSize = 1;

   TRACK_MALLOC (*this, 1);
}

template<class T>
SxArray<T>::SxArray (const SxList<T> &in, ssize_t chunkSizeIn)
   : SxMemConsumer (),
     SxThis<SxArray<T> > ()
{
   SX_CHECK (in.getSize() >= 0, in.getSize());
   SX_CHECK (chunkSizeIn > 0, chunkSizeIn);

   ssize_t n = in.getSize();

   chunkSize = chunkSizeIn;
   sizeAlloc = getSizeAlloc (n);

   if (n == 0)  {
      elements = NULL;
   }  else  {
      elements = new T [static_cast<size_t>(sizeAlloc)];
      typename SxList<T>::ConstIterator it = in.begin();
      for (ssize_t i=0; i < n; i++, it++)  elements[i] = *it;
   }
   size  = n;
   TRACK_MALLOC (*this, 1);
}

template<class T>
SxArray<T>::SxArray (const SxCList<T> &in)
   : SxMemConsumer (),
     SxThis<SxArray<T> > ()
{
   SX_CHECK (in.getSize() >= 0, in.getSize());

   ssize_t n = in.getSize();

   chunkSize = 1;
   sizeAlloc = getSizeAlloc (n);

   if (n == 0)  {
      elements = NULL;
   }  else  {
      elements = new T [static_cast<size_t>(sizeAlloc)];
      typename SxCList<T>::Node *srcPtr = in.firstElement;
      size_t i = 0;
      while (srcPtr) {
         elements[i++] = srcPtr->elem;
         srcPtr = srcPtr->next;
      }
   }
   size  = n;
   TRACK_MALLOC (*this, 1);
}


template<class T>
SxArray<T>::SxArray (const T *in, ssize_t n, ssize_t chunkSizeIn)
   : SxMemConsumer (),
     SxThis<SxArray<T> > ()
{
   SX_CHECK (in || (!in && n==0));
   SX_CHECK (chunkSizeIn > 0, chunkSizeIn);

   if (n == 0)  {
      elements  = NULL;
      size      = 0;
      sizeAlloc = 0;
      chunkSize = 1;
   }  else  {

      chunkSize  = chunkSizeIn;
      sizeAlloc  = getSizeAlloc (n);
      T *destPtr = elements = new T [static_cast<size_t>(sizeAlloc)];
      const T *srcPtr = in;
      for (ssize_t i=0; i < n; i++)  {
         //     elements[i]  = in[i]
         *destPtr++    = *srcPtr++;
      }
      size = n;
      TRACK_MALLOC (*this, 1);
   }
}

template <class T>
SxArray<T>::SxArray (const SxStack<T> &stack, ssize_t chunkSizeIn)
   : SxThis<SxArray<T> > (),
     elements(NULL), size(static_cast<ssize_t>(stack.getSize ()))
{
   SX_CHECK (chunkSizeIn > 0, chunkSizeIn);

   chunkSize = chunkSizeIn;
   sizeAlloc = getSizeAlloc (size);

   if (size == 0) return;   
   elements = new T[static_cast<size_t>(sizeAlloc)];
   stack.exportStack (elements, static_cast<size_t>(size));
}

template<class T>
template<long Begin,long End,long Step>
SxArray<T>::SxArray (const SxRange<Begin,End,Step> &range)
   : SxThis<SxArray<T> > ()
{
   ssize_t n = abs ((End - Begin) / Step);
   chunkSize = 1;
   if (n == 0)  {
      sizeAlloc = 0;
      elements  = NULL;
   }  else  {
      sizeAlloc= getSizeAlloc (n);
      elements = new T [static_cast<size_t>(sizeAlloc)];
   }
   size = n;

   typename SxRange<Begin,End,Step>::ConstIterator it = range.begin();
   for (ssize_t i=0; i < n; ++i, ++it)  elements[i] = (T)(*it);

}

template<class T>
SxArray<T>::SxArray (const std::initializer_list<T> &in, ssize_t chunkSizeIn)
   : SxThis<SxArray<T> >()
{
   SX_CHECK (in.size () >= 0, in.size ());
   SX_CHECK (chunkSizeIn > 0, chunkSizeIn);

   ssize_t n = in.size();

   chunkSize = chunkSizeIn;
   sizeAlloc = getSizeAlloc (n);

   if (n == 0)  {
      elements = NULL;
   }  else  {
      elements = new T [static_cast<size_t>(sizeAlloc)];
      ssize_t i = 0;
      for (auto &elem : in)  elements[i++] = elem;
   }
   size  = n;
   TRACK_MALLOC (*this, 1);
}

template<class T>
SxArray<T>::~SxArray ()
{
   if (sizeAlloc)  delete [] elements;
}


template<class T>
SxArray<T> &SxArray<T>::operator= (const SxArray<T> &in)
{
   if (&in == this)  return *this;

   if (sizeAlloc)  delete [] elements;

   chunkSize = in.chunkSize;
   SX_CHECK (chunkSize > 0, chunkSize);
   sizeAlloc = getSizeAlloc (in.size);

   if (in.size == 0)  {
      elements = NULL;
   }  else  {
      T *destPtr = elements  = new T [static_cast<size_t>(sizeAlloc)];
      T *srcPtr  = in.elements;
      for (ssize_t i=0; i < in.size; i++)  {
//        elements[i] =  in.elements[i];
         *destPtr++   = *srcPtr++;
      }
   }
   size = in.size;

   TRACK_MALLOC (*this, 1);
   return *this;
}


template<class T>
SxArray<T> &SxArray<T>::operator= (SxArray<T> &&in)
{
   if (&in == this)  return *this;

   if (sizeAlloc)  delete [] elements;

   chunkSize = in.chunkSize;
   SX_CHECK (chunkSize > 0, chunkSize);
   sizeAlloc = getSizeAlloc (in.size);

   if (in.size == 0)  {
      elements = NULL;
   }  else  {
      elements = in.elements;
   }
   size = in.size;

   in.elements = NULL;
   in.size = 0;
   in.sizeAlloc = 0;
   in.chunkSize = 1;

   TRACK_MALLOC (*this, 1);
   return *this;
}


template<class T>
SxArray<T> &SxArray<T>::operator= (const SxList<T> &in)
{
   ssize_t listSize = in.getSize ();
   resize (listSize);

   if (listSize > 0)  {
      // --- copy list
      typename SxList<T>::ConstIterator it = in.begin ();
      T *destPtr = elements;
      for (ssize_t i = 0; i < listSize; i++, ++it, ++destPtr)  {
         *destPtr = *it;
      }
   }

   TRACK_MALLOC (*this, 1);
   return *this;
}

template<class T>
SxArray<T> &SxArray<T>::operator= (const SxStack<T> &stack)
{
   resize (stack.getSize ());
   stack.exportStack(elements, size);
   return *this;
}

template<class T>
SxArray<T> &SxArray<T>::operator= (const std::initializer_list<T> &in)
{
   ssize_t listSize = in.size ();
   resize (listSize);

   if (listSize > 0)  {
      // --- copy list
      ssize_t i = 0;
      for (auto &elem : in)  elements[i++] = elem;
   }

   TRACK_MALLOC (*this, 1);
   return *this;
}


template<class T>
bool SxArray<T>::operator== (const SxArray<T> &in) const
{
   if (size != in.getSize()) return false;

   for (ssize_t i = 0; i < size; ++i)  {
      if (elements[i] != in(i))  {
         return false;
      }
   }

   return true;
}

template<class T>
bool SxArray<T>::operator!= (const SxArray<T> &in) const
{
   if (size != in.getSize ()) return true;

   for (ssize_t i = 0; i < size; ++i)  {
      if (elements[i] != in(i))  {
         return true;
      }
   }

   return false;
}

template<class T>
T &SxArray<T>::operator() (ssize_t i)
{
   SX_CHECK (i >= 0 && i < size, i, size);
   return elements[i];
}


template<class T>
const T &SxArray<T>::operator() (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < size, i, size);
   return elements[i];
}

template<class T>
T &SxArray<T>::operator() (SxAutoLoop &i)
{
   i.setLimit (size);
   SX_CHECK (i >= 0 && i < size, i.i, size);
   return elements[i];
}


template<class T>
const T &SxArray<T>::operator() (SxAutoLoop &i) const
{
   i.setLimit (size);
   SX_CHECK (i >= 0 && i < size, i.i, size);
   return elements[i];
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::getIterator (const typename SxArray<T>::SelIdx &idx)
{
   SX_CHECK (idx >= 0 && idx < size, idx, size);
   return typename SxArray<T>::Iterator (sx::Forward, this,
                                         idx);
}

template<class T>
typename SxArray<T>::ConstIterator
SxArray<T>::getIterator (const typename SxArray<T>::SelIdx &idx) const
{
   SX_CHECK (idx >= 0 && idx < size, idx, size);
   return typename SxArray<T>::ConstIterator (sx::Forward, this,
                                              idx);
}

template<class T>
typename SxArray<T>::ConstIterator
SxArray<T>::getConstIterator (const typename SxArray<T>::Iterator &it) const
{
   return it;
}

template<class T>
SxPtr<SxArray<T> > SxArray<T>::getContainer () const
{
   return this->getThis ();
}


template<class T>
const T SxArray<T>::at (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < size, i, size);
   return elements[i];
}


template<class T>
void SxArray<T>::set (const T &in)
{
   SX_CHECK (size > 0);
   T *ptr = elements;
   for (ssize_t i=0; i < size; i++)  {
//    elements[i] = in;
      *ptr++      = in;
   }
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::prepend (const T &in, ssize_t count)
{
   SX_CHECK (count > 0, count);

   return insert (0, in, count);
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::prepend (const SxArray<T> &in)
{
   return insert (0, in);
}


template<class T>
typename SxArray<T>::Iterator
SxArray<T>::prepend (SxArray<T> &&in)
{
   insert (0, std::move(in));

   return begin (0);
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::prepend (const std::initializer_list<T> &in)
{
   return insert (0, in);
}

template<class T>
void SxArray<T>::insertSpace (ssize_t idx, ssize_t count)
{
   SX_CHECK (idx >= 0 && idx <= size, idx, size);
   SX_CHECK (count > 0, count);

   ssize_t i;
   ssize_t n = getSize ();
   ssize_t newSize = getSize() + count;
   ssize_t newSizeAlloc = getSizeAlloc (newSize);

   if (newSizeAlloc == sizeAlloc)  {
      // --- shift (in-space compared to memmove)
      for (i=n-1; i >= idx; i--)  {
         elements[i+count] = elements[i];
      }
      size = newSize;
   }  else  {
      // --- create empty space during reallocation
      T *buffer = NULL;
      buffer = new T [newSizeAlloc];
      for (i=0; i < idx; i++)  {
         buffer[i] = elements[i];
      }
      for (i=idx; i < n; i++)  {
         buffer[i + count] = elements[i];
      }

      if (sizeAlloc)  delete [] elements;
      sizeAlloc = newSizeAlloc;
      elements  = buffer;
      buffer    = NULL;
      size      = newSize;
      TRACK_MALLOC (*this, 1);
   }
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::insert (ssize_t newPos, const T &in, ssize_t count)
{
   SX_CHECK (newPos >= 0 && newPos <= size, newPos, size);
   SX_CHECK (count > 0, count);

   ssize_t i;

   insertSpace (newPos, count);

   for (i=0; i < count; i++)  {
      elements[newPos+i] = in;
   }
   return begin (newPos);
}


template<class T>
typename SxArray<T>::Iterator
SxArray<T>::insert (ssize_t newPos, const SxArray<T> &in)
{
   SX_CHECK (newPos >= 0 && newPos <= size, newPos, size);
   SX_CHECK (in.getSize() >= 0, in.getSize());

   // --- possible cast from empty list
   if (in.getSize() == 0) return end ();

   ssize_t i;
   ssize_t count = in.getSize ();

   insertSpace (newPos, count);

   // --- insert the input to the empty space
   for (i=0; i < count; i++)  {
      elements[newPos+i] = in(i);
   }
   return begin (newPos);
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::insert (ssize_t newPos, SxArray<T> &&in)
{
   SX_CHECK (newPos >= 0 && newPos <= size, newPos, size);
   SX_CHECK (in.getSize() >= 0, in.getSize());

   // --- possible cast from empty list
   if (in.getSize() == 0) return end ();

   ssize_t i;
   ssize_t count = in.getSize ();

   insertSpace (newPos, count);

   // --- insert the input to the empty space
   for (i=0; i < count; i++)  {
      elements[newPos+i] = std::move(in(i));
   }
   return begin (newPos);
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::insert (ssize_t newPos, const std::initializer_list<T> &in)
{
   SX_CHECK (newPos >= 0 && newPos <= size, newPos, size);
   SX_CHECK (in.size () >= 0, in.size ());

   // --- possible cast from empty list
   if (in.size () == 0) return end ();

   ssize_t count = in.size ();

   insertSpace (newPos, count);

   ssize_t i = 0;
   // --- insert the input to the empty space
   for (auto &elem : in)  {
      elements[newPos+i] = elem;
      i++;
   }
   return begin (newPos);
}


template<class T>
typename SxArray<T>::Iterator
SxArray<T>::append (const T &in, ssize_t count)
{
   SX_CHECK (count > 0, count);

   ssize_t i;
   ssize_t n = getSize ();

   resize (n+count, true);

   for (i=0; i < count; i++) {
      elements[n+i] = in;
   }

   return begin (size - 1);
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::append (const SxArray<T> &in)
{
   if (in.getSize () > 0)  {
      ssize_t n = getSize ();
      ssize_t m = in.getSize ();
      resize (n + m, true);

      for (ssize_t i=0; i < m; i++)  {
         elements[n+i] = in(i);
      }
   }

   return begin (size - 1);
}


template<class T>
typename SxArray<T>::Iterator
SxArray<T>::append (SxArray<T> &&in)
{
   if (in.getSize () > 0)  {
      ssize_t n = getSize ();
      ssize_t m = in.getSize ();
      resize (n + m, true);

      for (ssize_t i=0; i < m; i++)  {
         elements[n+i] = std::move(in(i));
      }
   }

   return begin (size - 1);
}

template<class T>
typename SxArray<T>::Iterator
SxArray<T>::append (const std::initializer_list<T> &in)
{
   if (in.size () > 0)  {
      ssize_t n = getSize ();
      ssize_t m = in.size ();
      resize (n + m, true);

      ssize_t i = 0;
      for (auto &elem : in)  {
         elements[n+i] = elem;
         i++;
      }
   }

   return begin (size - 1);
}


template<class T>
void SxArray<T>::replace (ssize_t from, ssize_t nRemove, const SxArray<T> &in)
{
   SX_CHECK (from >= 0 && from <= size, from, size);
   SX_CHECK (in.getSize() >= 0, in.getSize());
   SX_CHECK (nRemove >= 0, nRemove);
   SX_CHECK (from + nRemove <= size, from + nRemove, size);

   ssize_t i;
   ssize_t n = getSize();
   ssize_t count = in.getSize();
   ssize_t newSize = getSize() - nRemove + in.getSize();
   ssize_t newSizeAlloc = getSizeAlloc (newSize);
   ssize_t idx = from + nRemove;
   ssize_t shift = from + count - idx;

   if (newSizeAlloc == sizeAlloc)  {
      if (shift < 0)  {
         // --- move left by (shift)
         for (i=idx; i < n; i++)  elements[i+shift] = elements[i];
      } else if (shift > 0)  {
         // --- move right by (shift)
         for (i=n-1; i >= idx; i--)   elements[i+shift] = elements[i];
      }
      size = newSize;
   }  else  {
      // --- create empty space during reallocation
      T *buffer = NULL;
      buffer = new T [newSizeAlloc];
      for (i=0; i < from; i++)  {
         buffer[i] = elements[i];
      }
      for (i=idx; i < n; i++)  {
         buffer[i + shift] = elements[i];
      }

      if (sizeAlloc)  delete [] elements;
      sizeAlloc = newSizeAlloc;
      elements  = buffer;
      buffer    = NULL;
      size      = newSize;
      TRACK_MALLOC (*this, 1);
   }

   // --- insert the input to the empty space
   for (i=0; i < count; i++)  {
      elements[from + i] = in(i);
   }
}

template<class T>
SxArray<T> SxArray<T>::join (const SxList<SxArray<T> > &list)
{
   SxArray<T> result;
   uint64_t size = 0;
   // Get required Array size
   for (const SxArray<T> &elem : list) {
      size += elem.getSize ();
   }
   result.resize (size);
   size = 0;

   // Trivial type. Memcpy.
   if (std::is_trivially_copyable<T>::value)  {

      size = 0;
      uint64_t nBytes = 0;
      // Copy memory directly.
      for (const SxArray<T> &elem : list) {
         nBytes = sizeof (T) * elem.getSize ();
         memcpy ((char *)result.elements + size,
                 (char *)elem.elements,
                 nBytes);
         size += nBytes;
      }

   // Non-trivial type. Element copy.
   } else {
      for (const SxArray<T> &listElem : list)  {
         for (const T &arrayElem : listElem)  {
            result(size) = arrayElem;
            size++;
         }
      }
   }
   return result;
}

template<class T>
ssize_t SxArray<T>::findPos (const T &in, ssize_t from) const
{
   SX_CHECK (from >= 0 && from <= size, from, size);

   ssize_t result = -1; // not found
   ssize_t i;
   ssize_t n = getSize ();

   // --- search for ==
   for (i=from; i < n; i++)  {
      if (elements[i] == in) {
         result = i;
         break;
      }
   }

   return result;
}


template<class T>
bool SxArray<T>::contains (const T &in) const
{
   return (findPos (in) >= 0);
}


template<class T>
void SxArray<T>::remove (ssize_t idx, ssize_t count)
{
   SX_CHECK (idx >= 0 && idx < size, idx, size);
   SX_CHECK (count > 0, count);
   SX_CHECK (idx+count <= size, idx+count, size);

   ssize_t i;
   ssize_t newSize = getSize() - count;
   ssize_t newSizeAlloc = getSizeAlloc (newSize);

   if (newSizeAlloc == sizeAlloc)  {
      // --- shift, (in-space compared to memmove)
      for (i=idx; i < newSize; i++)  {
         elements[i] = elements[i + count];
      }
      size = newSize;
   }  else  {
      // --- remove during reallocation
      T *buffer = NULL;
      if (newSizeAlloc > 0)  { 
         buffer = new T [newSizeAlloc];
         for (i=0; i < idx; i++)  {
            buffer[i] = elements[i];
         }
         for (i=idx; i < newSize; i++)  {
            buffer[i] = elements[i + count];
         }
      }
      if (sizeAlloc)  delete [] elements;
      sizeAlloc = newSizeAlloc;
      elements  = buffer;
      buffer    = NULL;
      size      = newSize;
      TRACK_MALLOC (*this, 1);
   }
}


template<class T>
void SxArray<T>::removeElement (const T &in)
{
   ssize_t idx;

   idx = findPos (in);
   if (idx >= 0)  {
      remove (idx);
   }
}

template<class T>
void SxArray<T>::removeItem (typename SxArray<T>::Iterator *it)
{
   SX_CHECK (it);
   if(it->isForward ()) {
      removeElement (*(*it));
      if (it->getIdx() > (getSize()-1))
         *it = end ();
   }
   else {
      removeElement (*(*it)--);
      if (it->getIdx() < 0)
         *it = toFirst ();
   }
}

template<class T>
void SxArray<T>::removeFirst ()
{
   remove (0);
}


template<class T>
void SxArray<T>::removeLast ()
{
   resize (getSize() - 1, true);
}


template<class T>
void SxArray<T>::removeAll ()
{
   resize (0);
}


template<class T>
void SxArray<T>::setChunkSize (ssize_t newChunkSize)
{
   SX_CHECK (newChunkSize > 0, newChunkSize);
   chunkSize = newChunkSize;
}


template<class T>
inline void SxArray<T>::resize (ssize_t newSize, bool keep)
{
   SX_CHECK (newSize >= 0, newSize);
   if (newSize == size)  return;

   ssize_t newSizeAlloc = getSizeAlloc(newSize);
   if (newSizeAlloc != sizeAlloc)  {
//      std::cout << "resize(" << newSize << ") alloc: " <<
//                   sizeAlloc << " to "<<newSizeAlloc<<std::endl;
      T *buffer = NULL;
      if (newSizeAlloc > 0)  {
         buffer = new T [static_cast<size_t>(newSizeAlloc)];
         if (keep)  {
            ssize_t i, min = (size < newSize) ? size : newSize;
            for (i=0; i < min; i++)   buffer[i] = elements[i];
         }
      }
      if (sizeAlloc)  delete [] elements;
      sizeAlloc = newSizeAlloc;
      elements  = buffer;
      buffer    = NULL;
      TRACK_MALLOC (*this, 1);
   }
   size = newSize;
}

template<>
inline void SxArray<char>::resize (ssize_t newSize, bool keep)
{
   SX_CHECK (newSize >= 0, newSize);
   if (newSize == size)  return;

   ssize_t newSizeAlloc = newSize;
   if (chunkSize > 1)  {
      newSizeAlloc = getSizeAlloc(newSize);
   }
   if (newSizeAlloc != sizeAlloc)  {
      char *buffer = NULL;
      if (newSizeAlloc > 0)  {
         buffer = new char [static_cast<size_t>(newSizeAlloc)];
         if (keep)  {
            ssize_t min = (size < newSize) ? size : newSize;
            // --- TODO: memcpy only for large n
            ::memcpy (buffer, elements, static_cast<size_t>(min));
            //for (i=0; i < min; i++)   buffer[i] = elements[i];
         }
      }
      if (sizeAlloc)  delete [] elements;
      sizeAlloc = newSizeAlloc;
      elements  = buffer;
      buffer    = NULL;
      TRACK_MALLOC (*this, 1);
   }
   size = newSize;
}

template<class T>
ssize_t SxArray<T>::getSize () const
{
   return size;
}


template<class T>
ssize_t SxArray<T>::getSizeAlloc () const
{
   return sizeAlloc;
}


template<class T>
ssize_t SxArray<T>::getSizeAlloc (ssize_t newSize) const
{
   return ((newSize/chunkSize)*chunkSize) + (newSize%chunkSize > 0?chunkSize:0);
}


template<class T>
void SxArray<T>::sortByIdx (const SxArray<ssize_t> &idxArray)
{
   SX_CHECK (idxArray.size > 1, idxArray.size);
   SX_CHECK (size == idxArray.size, size, idxArray.size);

   ssize_t i;

   T *temp = new T [size];
   // --- copy array
   const T *srcPtr = elements;
   T *destPtr      = temp;
   for (i=0; i < size; i++)
//     temp[i]   =  elements[i]
      *destPtr++ = *srcPtr++;

   // --- rearrange elements
   const ssize_t *inPtr = idxArray.elements;
   destPtr = elements;
   for (i=0; i < size; i++)
//     elements[i]  = temp[idxArray.elements[i]];
      *destPtr++ = temp[ *inPtr++ ];

   delete [] temp;
}

template<class T>
SxArray<ssize_t> SxArray<T>::getSortIdx () const
{
   SxArray<ssize_t> sortIdx (getSize());

   if (getSize() == 0)  {
      return sortIdx;
   } else if (getSize() == 1)  {
      sortIdx(0) = 0;
      return sortIdx;
   }

   ::sort (sortIdx.elements, elements, getSize());
   return sortIdx;
}


template<class T>
SxList<T> SxArray<T>::toList () const
{
   SxList<T> res;
   for (ssize_t i=0; i < size; ++i)  res << elements[i];
   return res;
}



// --------------------------------------------------------------------------
template<class T>
SxArray<typename T::TypeMapper::Type> chop (const SxArray<T> &a)
{
   SX_CHECK (a.size > 0, a.size);
   SxArray<typename T::TypeMapper::Type> res; res.resize (a.size);
   for (ssize_t i=0; i < a.size; i++)
      res(i) = a(i).chop();
   return res;
}



template<class T>
SxArray<T> operator- (const SxArray<T> &a, const SxArray<T> &b)
{
   SX_CHECK (a.size == b.size && a.size > 0, a.size, b.size);
   SxArray<T> res; res.resize (a.size);
   for (ssize_t i=0; i < a.size; i++)
      res(i) = a(i) - b(i);
   return res;
}


template<class T>
SxArray<T> operator+ (const SxArray<T> &a, const SxArray<T> &b)
{
   SX_CHECK (a.size == b.size && a.size > 0, a.size, b.size);
   SxArray<T> res; res.resize (a.size);
   for (ssize_t i=0; i < a.size; i++)
      res(i) = a(i) + b(i);
   return res;
}


template<class T>
SxArray<T> operator^ (const SxArray<T> &a, const SxArray<T> &b)
{
   SX_CHECK (a.size == b.size && a.size > 0, a.size, b.size);
   SxArray<T> res; res.resize (a.size);
   for (ssize_t i=0; i < a.size; i++)
      res(i) = a(i) ^ b(i);
   return res;
}


// --------------------------------------------------------------------------
template<class T>
std::ostream& operator<< (std::ostream &s, const SxArray<T> &in)
{
   s << "[";
   for (ssize_t i = 0; i < in.getSize (); ++i)  {
      if (i != 0)  s << ", ";
      s << in(i);
   }
   s << "];";
   return s;
}

// --------------------------------------------------------------------------

template<class T>
size_t getNBytes (const SxArray<T> &in)
{
   size_t nBytes = 0;
   size_t nElem = in.getSizeAlloc ();
   // --- elementwisely evaluation in case of something like
   //     SxArray<SxList> or SxArray<SxArray>
   for (size_t i=0; i < nElem; ++i)
      nBytes += getNBytes (in(i));
   nBytes += sizeof (T *)       // *elements
           + sizeof (ssize_t)   //  size
           + sizeof (ssize_t)   //  sizeAlloc
           + sizeof (ssize_t);  //  chunkSize
   return nBytes;
}



// --- obsolete sort() function
template<class T>
void sort (ssize_t *sortIdx, const T *inValues, ssize_t nElem)
{
   SX_CHECK (sortIdx);
   SX_CHECK (inValues);
   SX_CHECK (nElem >= 1, nElem);

   if (nElem == 1)  return;

   ssize_t i, j, ir, l, ind;
   T q;

   for (i=0; i < nElem; i++)  sortIdx[i] = i;

   l=nElem/2+1;
   ir=nElem;

   for (;;)  {
      if (l > 1)  {
         l--;
         ind = sortIdx[l-1];
         q   = inValues[ind];
      }  else  {
         ind = sortIdx[ir-1];
         q   = inValues[ind];
         sortIdx[ir-1] = sortIdx[0];
         ir--;
         if (ir == 1)  {
            sortIdx[0]=ind;
            return;
         }
      }
      i=l;
      j=l+l;
      while ( j <= ir)  {
         if (j < ir)
            if (inValues[sortIdx[j-1]] < inValues[sortIdx[j+1-1]]) j++;

         if (q < inValues[sortIdx[j-1]])  {
            sortIdx[i-1]=sortIdx[j-1];
            i=j;
            j=j+j;
         }  else  {
            j=ir+1;
         }
      }
      sortIdx[i-1]=ind;
   }

}
#endif /* _SX_ARRAY_H_ */
