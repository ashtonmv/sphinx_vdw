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
#ifndef _SX_SORTED_LIST_H_
#define _SX_SORTED_LIST_H_

#include <SxList.h>
#include <stdio.h>

/** \brief A sorted list.
    \author Christoph Freysoldt */
template<class T>
class SxSortedList : public SxList<T>
{
   public:
      SxSortedList ();
      /** \brief Make a sorted list out of an unsorted one */
      SxSortedList (const SxList<T> &);
      SxSortedList (const SxSortedList<T> &);
      SxSortedList (SxSortedList<T> &&);
      SxSortedList (const std::initializer_list<T> &);
      ~SxSortedList ();

      SxSortedList<T> &operator= (const SxSortedList<T> &);
      SxSortedList<T> &operator= (SxSortedList<T> &&);
      /** \brief Insert item in sorted list, returns new position

          \param in   element to be inserted
          \return position of new element in list

          The return value allows to build sorted lists of complicated objects:
\code
  SxVector3<Double> v;
  double length;
  SxSortedList<double> lengthes;
  SxList<SxVector3<Double> > vectors;
  int pos;
  ...
  // some loop setting v to different values
     length = sqrt(v.absSqr ());
     pos = lengthes.append (length);
     vectors.insert (pos, v);
\endcode
       */
      ssize_t append (const T &in);
      ssize_t append (const std::initializer_list<T> &in);
      /** \brief Insert new element (like #append, but returns the list) */
      inline SxSortedList<T>& operator<< (const T&);
      inline SxSortedList<T>& operator<< (const std::initializer_list<T>&);
};




template<class T>
SxSortedList<T>::SxSortedList () : SxList<T> ()
{
   // empty
}

template <class T>
inline SxSortedList<T> &SxSortedList<T>::operator<< (const T &in)
{
   append (in);
   return *this;
}

template <class T>
inline SxSortedList<T>&
SxSortedList<T>::operator<< (const std::initializer_list<T> &in)
{
   for (const T &elem : in)
      append (elem);
   return *this;
}

template<class T>
SxSortedList<T>::SxSortedList (const SxList<T> &in) : SxList<T> ()
{
   typename SxList<T>::ConstIterator it;
   for (it = in.begin (); it != in.end (); ++it)
      this->append (*it);
}

template<class T>
SxSortedList<T>::SxSortedList (const SxSortedList<T> &in) : SxList<T> (in)
{ }

template<class T>
SxSortedList<T>::SxSortedList (SxSortedList<T> &&in) : SxList<T> (std::move(in))
{ }

template<class T>
SxSortedList<T>::SxSortedList (const std::initializer_list<T> &in)
   : SxList<T> ()
{
   for (const T &elem : in)
      this->append (elem);
}

template<class T>
SxSortedList<T>::~SxSortedList ()
{
   // empty
}

template<class T>
SxSortedList<T> &SxSortedList<T>::operator= (const SxSortedList<T> &in)
{
   if (this == &in)  return *this;
   SxList<T>::operator= (in);
   return *this;
}

template<class T>
SxSortedList<T> &SxSortedList<T>::operator= (SxSortedList<T> &&in)
{
   if (this == &in)  return *this;
   SxList<T>::operator= (std::move(in));
   return *this;
}

template<class T>
ssize_t SxSortedList<T>::append (const T &in)
{
   if (this->getSize() == 0)  {
      ((SxList<T> *)this)->append (in);
      return 0;
   }
   /*
   if (this->getSize() == 1)  {
      if ( this->firstElement->element < in )
         ((SxList<T> *)this)->append (in);
      else
         prepend (in);
      return 0;
   }
   */

   typename SxList<T>::Iterator it = this->begin ();
   ssize_t id = -1;
   ssize_t n = this->getSize();
   for (ssize_t i = 0; i < n; i++, it++)  {
      if ( *it > in )  {
         id = i;  break;
      }
   }
   if ( id == 0 )  {
      ((SxList<T> *)this)->prepend (in);
      return 0;
   }
   if ( id < 0 )  {
      ((SxList<T> *)this)->append (in);
      return n;
   }
   ((SxList<T> *)this)->insert (id, in);
   
   return id;
}

template<class T>
ssize_t SxSortedList<T>::append (const std::initializer_list<T> &in)
{
   ssize_t res = -1;
   for (const T &elem : in) {
      res = this->append (elem);
   }
   return res;
}

#endif /* _SX_SORTED_LIST_H_ */
