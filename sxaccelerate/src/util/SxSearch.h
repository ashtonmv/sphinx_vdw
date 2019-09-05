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

#ifndef _SX_SEARCH_H_
#define _SX_SEARCH_H_

#include <SxError.h>
#include <SxArray.h>

/** \brief Search algorithms

    \b example:
\code
   SxArray<SxString> keywords(6);
   keywords(0) = "break";
   keywords(1) = "case";
   keywords(2) = "class";
   keywords(3) = "delete";
   keywords(4) = "do";
   keywords(5) = "float";
   
   // --- The array has to be sorted in order to get meaningful results
   // --- from binary search algorithm. It is already sorted in this example.
   // keywords.sort (); // simple and slow sorting
   
   std::cout << "break: " << SxSearch<SxString>::binary (keywords, "break") <<
                std::endl;
   std::cout << "do: "    << SxSearch<SxString>::binary (keywords, "do") <<
                std::endl;
   std::cout << "float: " << SxSearch<SxString>::binary (keywords, "float") <<
                std::endl;
\endcode

    
    \ingroup Tools

    \author Vaclav Bubnik, bubnik@mpie.de */
template<class T>
class SxSearch
{
   public:
      /** Binary searching.
          \param v_ Specifies the array. The array has to be sorted.
          \param in_ An element which will be searched.
          \return index of the element or -1 if the element was not found
          
          Binary searching is faster than findPos for n > 6 elements.
          Binary searching is 100x faster than findPos for an array
          with 4000 elements and more. */
      static inline ssize_t binary (const SxArray<T> &v_, const T &in_);
};


/** \brief Search algorithms with custom comparator

    I have found this in C++ book written by Bjarne Stroustrup.

    \b example:
\code
// --- comparator
class SxCmpDouble
{
   public:
      static inline bool equal (double a, double b)
      {
         return (fabs(a - b) < 1.0e-4 ? true : false);
      }
};    

...

SxArray<double> r(2);
r(0) = 1.0;
r(1) = 2.5;
std::cout << "SxSearchCmp: " <<
             SxSearchCmp<double, SxCmpDouble>::linear (r, 2.5) <<
             std::endl;
\endcode    

    \ingroup Tools

    \author Vaclav Bubnik, bubnik@mpie.de */
template<class T, class Comp>
class SxSearchCmp
{
   public:
      /** Linear searching.
          \param v_ Specifies the array. The array can be unsorted.
          \param in_ An element which will be searched.
          \return index of the element or -1 if the element was not found */
      static inline ssize_t linear (const SxArray<T> &v_, const T &in_);
      
      /** Binary searching.
          \param v_ Specifies the array. The array has to be sorted.
          \param in_ An element which will be searched.
          \return index of the element or -1 if the element was not found */
      static inline ssize_t binary (const SxArray<T> &v_, const T &in_);
};


template<class T>
ssize_t SxSearch<T>::binary (const SxArray<T> &v_, const T &in_)
{
   ssize_t result = -1;
   ssize_t i;
   ssize_t first = 0;
   ssize_t last = v_.getSize ();
   
#ifndef NDEBUG
   // --- test if the array is sorted
   for (i=first; i < last-1; i++)  {
      if (v_(i) > v_(i+1))  {
         // --- the array is not sorted
         std::cerr << "SxSearch::binary: The array is not sorted." << std::endl;
         SX_EXIT;
      }
   }
#endif
   
   // --- search for == 
   // --- TODO: can be further optimized
//   last = v_.getSize () - 1;
//   while (first <= last)  {
//      i = (first + last) / 2;
//      if (in_ < v_(i))  {
//         last = i - 1;
//      }
//      else if (in_ > v_(i))  {
//         first = i + 1;
//      }  else  {
//         // --- found at index i
//         result = i;
//         break;                        
//      }
//   }
   
   // --- optimized (worst case time consumption)
   // --- n   32 : 1.94x faster
   // --- n 4096 : 2.0x faster 
   while (first < last)  {
      i = ((size_t)first + (size_t)last) >> 1;
      if (v_(i) < in_)  {
         first = i + 1;
      }  else  {
         last = i;
      }
   }
   if ((first < v_.getSize()) && (v_(first) == in_))  {
      result = first;
   }

   return result;
}


template<>
ssize_t SxSearch<float>::binary (const SxArray<float> &,
                                 const float          &)
{
   // --- can not search for equal floats
   SX_EXIT;
   
   return -1;
}


template<>
ssize_t SxSearch<double>::binary (const SxArray<double> &,
                                  const double          &)
{
   // --- can not search for equal doubles
   SX_EXIT;
   
   return -1;
}


template<class T, class Comp>
ssize_t SxSearchCmp<T, Comp>::linear (const SxArray<T> &v_, const T &in_)
{
   ssize_t result = -1;
   ssize_t i;
   ssize_t n = v_.getSize ();
   
   // --- search for ==   
   for (i=0; i < n; i++)  {
      if (Comp::equal (v_(i), in_)) {
         result = i;
         break;
      }
   }
   
   return result;
}


template<class T, class Comp>
ssize_t SxSearchCmp<T, Comp>::binary (const SxArray<T> &v_, const T &in_)
{
   ssize_t result = -1;
   ssize_t i;
   ssize_t first = 0;
   ssize_t last = v_.getSize ();
   
   while (first < last)  {
      i = ((size_t)first + (size_t)last) >> 1;
      if (Comp::less (v_(i), in_))  {
         first = i + 1;
      }  else  {
         last = i;
      }
   }
   if ((first < v_.getSize()) && Comp::equal (v_(first), in_))  {
      result = first;
   }

   return result;
}

#endif /* _SX_SEARCH_H_ */
