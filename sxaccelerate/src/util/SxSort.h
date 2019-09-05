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

#ifndef _SX_SORT_H_
#define _SX_SORT_H_

#include <SxError.h>
#include <SxArray.h>

/** \brief Sorting template

    TODO: add some comparator template
    TODO: template for any container (vector, array, ...):
          then the current T will be T::typeMapper

    \ingroup Tools

    \author Sixten Boeck, boeck@mpie.de
    \author Vaclav Bubnik, bubnik@mpie.de */
template<class T>
class SxSort
{
   public:
      /** Default constructor. */
      SxSort ();

      /** Quick sort.
          \param v_ The array to sort.
          \param idx_ The array containing the sorting indices.
          \return the indices of the sorted v_.

\code
   SxArray<float> v(3);
   SxArray<size_t> indices; // empty idx
   SxSort<float>  sort;
   int i;

   // --- fill with the values
   v(0) = 3.0f;
   v(1) = 1.0f;
   v(2) = 2.0f;

   // --- sort
   indices = sort.quickSortToIdx (v, &indices);  // presorted
   // or
   indices = sort.quickSortToIdx (v);

   // --- print
   for (i=0; i < indices.getSize(); i++)  {
      cout << v(indices(i)) << endl;
   }

   // --- ouput
   // 1.0f
   // 2.0f
   // 3.0f
\endcode 

      The idx_ is automatically filled as an sorted array with values
      { 0, 1, 2, ... getSize()-1 }
      whenever the array of indices (idx_) has a different size than
      the array with values (v_).

      The more presorted idx_ => faster sorting.
      */
      inline SxArray<ssize_t> quickSortToIdx (const SxArray<T> &v_) const;
      SxArray<ssize_t> &quickSortToIdx (const SxArray<T> &v_,
                                        SxArray<ssize_t> *idxPtr) const;

   protected:
      /** The size of the stack used for iterative quickSort. 
       */
      const ssize_t QSORT_STACK_SIZE;

      /** The minimal number of elements for quickSort,
          insertion sort will be performed otherwise. */
      const ssize_t QSORT_MIN_ELEMENTS;
};

// Default constructor.
template<class T>
SxSort<T>::SxSort ()
   : QSORT_STACK_SIZE (1024),
     QSORT_MIN_ELEMENTS(8)
{
   // empty
}

template<class T>
SxArray<ssize_t> SxSort<T>::quickSortToIdx (const SxArray<T> &v_) const
{
   SxArray<ssize_t> idx (v_.getSize());
   ssize_t i, nElements = idx.getSize();
   for (i=0; i < nElements; i++)  idx(i) = i;

   return quickSortToIdx (v_, &idx);
}

template<class T>
SxArray<ssize_t> &SxSort<T>::quickSortToIdx (const SxArray<T> &v_,
                                             SxArray<ssize_t>  *idxPtr) const
{
   SxArray<ssize_t> &idx = *idxPtr;
   SX_CHECK (v_.getSize() == idx.getSize(),
             v_.getSize(),   idx.getSize());

   ssize_t i, nElements = v_.getSize();

   ssize_t itLeft  = 0;
   ssize_t itRight = nElements-1;

   ssize_t swap;
   T x;
   ssize_t ix;

   ssize_t j, m, n;

   static SxArray<ssize_t> stack(QSORT_STACK_SIZE);
   ssize_t *sp = stack.elements;

   // --- first two bounds to the stack
   *sp++ = itLeft;
   *sp++ = itRight;

   while (1)  { /* A */
      if (sp != stack.elements)  {
         // --- recovery bounds from the stack
         sp -= 2;
         i = sp[0];
         j = sp[1];

         if (j-i >= QSORT_MIN_ELEMENTS) {
            // --- more than minimum for quickSort

            // --- make backup
            itLeft = i;
            itRight = j;

            // --- median
            x = v_(idx(i + ((j - i)/2)));

            do {
               while (x > v_(idx(i))) i++;
               while (v_(idx(j)) > x) j--;

               if (i <= j) {
                  // --- swap the elements
                  swap   = idx(i);
                  idx(i) = idx(j);
                  idx(j) = swap;
                  i++;
                  j--;
               }
            } while (i < j);

            // --- left
            *sp++ = itLeft;
            *sp++ = j;

            // --- right
            *sp++ = i;
            *sp++ = itRight;
         }
         else if (j-i >= 1)  { 
            // --- [2..QSORT_MIN_ELEMENTS] insertion sort
            for (m=i+1; m <= j; m++)  {
               ix = idx(m);
               x  = v_(ix);
               // --- sort x to the beginnig of the subarray [i..j]
               n = m-1;
               while ((n >= i) && (v_(idx(n)) > x)) {
                  idx(n+1) = idx(n);
                  n--;
               }
               idx(n+1) = ix;
            }
         }
      }  else {
         break; /* A */
      }
   }

   // --- sorted indices
   return idx;
}

#endif /* _SX_SORT_H_ */
