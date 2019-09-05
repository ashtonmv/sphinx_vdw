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

#ifndef _SX_SPARSE_UTIL_H_
#define _SX_SPARSE_UTIL_H_

#include <SxError.h>
#include <SxArray.h>
#include <SxSearch.h>
#include <SxBitArray.h>


/** \brief A utility for sparse arrays.

    \b example:
\code
   SxArray<SxString> m;
   SxArray<ssize_t> idx;

   m = SxList<SxString> () << "box" << "line" << "point" << "quad";
   idx.resize (2);
   idx(0) = 1;
   idx(1) = 3;

   // --- prints: [line, quad];
   std::cout << SxSparseUtil<SxString>::gather (m, idx) << std::endl;
\endcode

    \ingroup Tools

    \author Vaclav Bubnik, bubnik@mpie.de */
template<class T>
class SxSparseUtil
{
   public:
      /** Gather values using the indices.
          \param in_ Specifies the input array.
          \param idx_ Specifies the indices of elements to gather.
          \return the array with gathered elements.

          \begin{verbatim}
          in(4) idx(2) -> array(2)
            A      0         A
            B      2         C
            C
            D
          \end{verbatim}

          Gather exists as a special instruction \c lvi in vector processors
          with data parallelism. It is used to load the elements from a sparse
          vector in memory to the vector register file in CPU. */
      static inline SxArray<T> gather (const SxArray<T>       &in_,
                                       const SxArray<ssize_t> &idx_);

      /** Scatter values using the indices.
          \param src_ Specifies the elements to scatter.
          \param idx_ Specifies the indices. Same size as src_.
          \param dst_ Specifies the destination.
          \return modified input with scattered values (the same as dst_) 

          \begin{verbatim}
          src(2) idx(2) dst(4) -> dst(4)
           a       0      A         a
           c       2      B         B
                          C         c
                          D         D
          \end{verbatim}

          Scatter exists as a special instruction \c svi in vector processors
          with data parallelism. It is used to store the elements from the 
          vector register file in CPU to a sparse vector in memory.*/
      static inline SxArray<T> &scatter (const SxArray<T>       &src_,
                                         const SxArray<ssize_t> &idx_,
                                         SxArray<T>             *dst_);

      /** Scatter one value using the indices.
          \param in_ Specifies the element to scatter.
          \param idx_ Specifies the indices.
          \param dst_ Specifies the destination.
          \return modified input with scattered values (the same as dst_) 

          \begin{verbatim}
          in   idx(2) dst(4) -> dst(4)
           x     0      a         x
                 2      B         B
                        c         x
                        D         D
          \end{verbatim}

          This is the scalar version of scattering. It is equivalent to
          \e dst.set(in) limited to selected indices specified in \a idx. */
      static inline SxArray<T> &scatter (const T                &in_,
                                         const SxArray<ssize_t> &idx_,
                                         SxArray<T>             *dst_);

      /** Mask the input array.
          \param in_ Specifies the elements.
          \param mask_ Specifies the mask (a set of 0/1).
          \param pattern_ Specifies the pattern {0, 1} used to indicate
                          presence of elements.
          \return the array with masked elements. Only the elements that have
                  corresponding mask set to the value of \a pattern are copied.

          \begin{verbatim}
          in(4) mask(4) -> array(2)
            A      0         B
            B      1         C
            C      1
            D      0
          \end{verbatim}

          Equal to gather (in, maskToIdx (mask, pattern)), but faster. */
      static inline SxArray<T> mask (const SxArray<T>    &in_,
                                     const SxArray<int>  &mask_,
                                     int                 pattern_);


      /** Convert mask to indices.
          \param mask_ Specifies the mask (a set of 0/1).
          \param pattern_ Specifies the pattern {0, 1} used to indicate
                          presence of indices.
          \return the array with indices.

          \begin{verbatim}
          mask(4) pattern -> array(2)
            0        1         1
            1                  2
            1
            0
          \end{verbatim} */
      static inline SxArray<ssize_t> maskToIdx (const SxArray<int>  &mask_,
                                                int                 pattern_);

      /** Convert mask to indices.
          \param mask_ Specifies the mask.
          \param pattern_ Specifies the pattern {0, 1} used to indicate
                          presence of indices.
          \return the array with indices. */
      static inline SxArray<ssize_t> maskToIdx (const SxBitArray &mask_,
                                                int              pattern_);

      /** Convert indices to mask.
          \param pattern_ Specifies the pattern {0, 1} used to indicate
                          presence of indices.
          \param idx_ Specifies the indices.
          \param dst_ Specifies the destination mask to fill. The mask must
                      contain enough space to hold all indices and will be
                      not resized.
          \return the resulting mask (the same as dst_)

          \begin{verbatim}
          pattern idx(2) dst(4)-> dst(4)
             1       0      x        1
                     1      x        1
                            x        0
                            x        0
          \end{verbatim} */
      static inline SxArray<int> &idxToMask (int                    pattern_,
                                             const SxArray<ssize_t> &idx_,
                                             SxArray<int>           *dst_);

      /** Convert indices to mask.
          \param pattern_ Specifies the pattern {0, 1} used to indicate
                          presence of indices.
          \param idx_ Specifies the indices.
          \param dst_ Specifies the destination mask to fill. The mask must
                      contain enough space to hold all indices and will be
                      not resized.
          \return the resulting mask (the same as dst_) */
      static inline SxBitArray &idxToMask (int                    pattern_,
                                           const SxArray<ssize_t> &idx_,
                                           SxBitArray             *dst_);


      /** Collect all indices of the specified value.
          \param src_ Specifies the source.
          \param in_ Specifies the element to find.
          \return the array with indices to all occurrences of the element.

          \begin{verbatim}
          src(4) in -> array(2)
             u    u      0
             v           2
             u
             w
          \end{verbatim} */
      static inline SxArray<ssize_t> getIdx (const SxArray<T> &src_,
                                             const T          &in_);

      /** Find indices of all specified values.
          \param srcA_ Specifies the source.
          \param srcB_ Specifies the samples to find.
          \param sorted_ true if the source is sorted, faster
          \return the array with indices for each sample, index -1 states
                  that a sample was not found in the source array.

          \begin{verbatim}
          srcA(4) srcB(3) -> array(2)
             a      c           2
             b      d           3
             c      f          -1  (f not found)
             d
          \end{verbatim} */
      static inline SxArray<ssize_t> getIdx (const SxArray<T> &srcA_,
                                             const SxArray<T> &srcB_,
                                             bool             sorted_=false);
};

// --- Gather values using the indices.
template<class T>
SxArray<T> SxSparseUtil<T>::gather (const SxArray<T>       &in_,
                                    const SxArray<ssize_t> &idx_)
{
   SxArray<T> result;
   ssize_t i;
   ssize_t n = idx_.getSize ();

   result.resize (n);

   for (i=0; i < n; i++)  {
      SX_CHECK (idx_(i)>=0 && idx_(i)<in_.getSize(),idx_(i),in_.getSize());
      result(i) = in_(idx_(i));
   }

   return result;
}

// --- Scatter one value using the indices.
template<class T>
SxArray<T> & SxSparseUtil<T>::scatter (const SxArray<T>       &src_,
                                       const SxArray<ssize_t> &idx_,
                                       SxArray<T>             *dst_)
{
   SX_CHECK (dst_);
   SX_CHECK (idx_.getSize() == src_.getSize(),
             idx_.getSize(), src_.getSize());

   ssize_t i;
   ssize_t n = idx_.getSize ();

   for (i=0; i < n; i++)  {
      SX_CHECK (idx_(i) >= 0 && idx_(i) < dst_->getSize(),
                idx_(i), dst_->getSize());
      (*dst_)(idx_(i)) = src_(i);
   }

   return (*dst_);
}

// --- Scatter one value using the indices.
template<class T>
SxArray<T> & SxSparseUtil<T>::scatter (const T                &in_,
                                       const SxArray<ssize_t> &idx_,
                                       SxArray<T>             *dst_)
{
   SX_CHECK (dst_);

   ssize_t i;
   ssize_t n = idx_.getSize ();

   for (i=0; i < n; i++)  {
      SX_CHECK (idx_(i) >= 0 && idx_(i) < dst_->getSize(),
                idx_(i), dst_->getSize());
      (*dst_)(idx_(i)) = in_;
   }

   return (*dst_);
}

// --- Mask the input array.
template<class T>
SxArray<T> SxSparseUtil<T>::mask (const SxArray<T>    &in_,
                                  const SxArray<int>  &mask_,
                                  int                 pattern_)
{
   SX_CHECK (in_.getSize() == mask_.getSize(),
             in_.getSize(), mask_.getSize());

   ssize_t i;
   ssize_t nUsed;
   ssize_t n;
   SxArray<T> result;

   n = mask_.getSize ();
   result.resize (n);

   nUsed = 0;
   if (pattern_ == 1)  {
      for (i=0; i < n; i++)  {
         if (mask_(i) == 1)  {
            result(nUsed) = in_(i);
            nUsed++;
         }
      }
   } else if (pattern_ == 0)  {
      for (i=0; i < n; i++)  {
         if (mask_(i) == 0)  {
            result(nUsed) = in_(i);
            nUsed++;
         }
      }
   }  else  {
      SX_EXIT;
   }

   result.resize (nUsed, true);

   return result;
}

// --- Convert mask to indices.
template<class T>
SxArray<ssize_t> SxSparseUtil<T>::maskToIdx (const SxArray<int> &mask_,
                                             int                pattern_)
{
   ssize_t i;
   ssize_t nUsed;
   ssize_t n;
   SxArray<ssize_t> result;

   n = mask_.getSize ();
   result.resize (n);

   nUsed = 0;
   if (pattern_ == 1)  {
      for (i=0; i < n; i++)  {
         if (mask_(i) == 1)  {
            result(nUsed) = i;
            nUsed++;
         }
      }
   } else if (pattern_ == 0)  {
      for (i=0; i < n; i++)  {
         if (mask_(i) == 0)  {
            result(nUsed) = i;
            nUsed++;
         }
      }
   }  else  {
      SX_EXIT;
   }

   result.resize (nUsed, true);

   return result;
}

// --- Convert mask to indices.
template<class T>
SxArray<ssize_t> SxSparseUtil<T>::maskToIdx (const SxBitArray &mask_,
                                             int              pattern_)
{
   ssize_t i;
   ssize_t nUsed;
   ssize_t n;
   SxArray<ssize_t> result;

   n = (ssize_t)mask_.getSize ();
   result.resize (n);

   nUsed = 0;
   if (pattern_ == 1)  {
      for (i=0; i < n; i++)  {
         if (mask_(i) == true)  {
            result(nUsed) = i;
            nUsed++;
         }
      }
   } else if (pattern_ == 0)  {
      for (i=0; i < n; i++)  {
         if (mask_(i) == false)  {
            result(nUsed) = i;
            nUsed++;
         }
      }
   }  else  {
      SX_EXIT;
   }

   result.resize (nUsed, true);

   return result;
}

// --- Convert indices to mask.
template<class T>
SxArray<int> & SxSparseUtil<T>::idxToMask (int                    pattern_,
                                           const SxArray<ssize_t> &idx_,
                                           SxArray<int>           *dst_)
{
   SX_CHECK (dst_);
   SX_CHECK (pattern_ == 0 || pattern_ == 1, pattern_);

   ssize_t i;
   ssize_t n;

   if ((*dst_).getSize () > 0)  {
      (*dst_).set (1 - pattern_);
   }

   n = idx_.getSize ();

   for (i=0; i < n; ++i)  {
      SX_CHECK (idx_(i) >= 0 && idx_(i) < dst_->getSize(),
                idx_(i), dst_->getSize());
      (*dst_)(idx_(i)) = pattern_;
   }

   return (*dst_);
}

// --- Convert indices to mask.
template<class T>
SxBitArray &SxSparseUtil<T>::idxToMask (int                    pattern_,
                                        const SxArray<ssize_t> &idx_,
                                        SxBitArray             *dst_)
{
   SX_CHECK (dst_);
   SX_CHECK (pattern_ == 0 || pattern_ == 1, pattern_);

   ssize_t i;
   ssize_t n;

   if ((*dst_).getSize () > 0)  {
      (*dst_) = (1 - pattern_);
   }

   n = idx_.getSize ();

   if (pattern_ == 1)  {
      for (i=0; i < n; ++i)  {
         SX_CHECK (idx_(i) >= 0 && idx_(i) < (ssize_t)dst_->getSize(),
                   idx_(i), dst_->getSize());
         (*dst_).setBit (idx_(i));
      }
   }  else  {
      for (i=0; i < n; ++i)  {
         SX_CHECK (idx_(i) >= 0 && idx_(i) < (ssize_t)dst_->getSize(),
                   idx_(i), dst_->getSize());
         (*dst_).unsetBit (idx_(i));
      }
   }

   return (*dst_);
}

// --- Collect all indices of the specified value.
template<class T>
SxArray<ssize_t> SxSparseUtil<T>::getIdx (const SxArray<T>   &src_,
                                          const T            &in_)
{
   ssize_t i;
   ssize_t nUsed;
   ssize_t n;
   SxArray<ssize_t> result;

   n = src_.getSize ();
   result.resize (n);

   nUsed = 0;
   for (i=0; i < n; i++)  {
      if (src_(i) == in_)  {
         result(nUsed) = i;
         nUsed++;
      }
   }

   result.resize (nUsed, true);

   return result;
}

// --- Find indices of all specified values.
template<class T>
SxArray<ssize_t> SxSparseUtil<T>::getIdx (const SxArray<T> &src1_,
                                          const SxArray<T> &src2_,
                                          bool             sorted_)
{
   ssize_t i1;
   ssize_t i2;
   ssize_t n1;
   ssize_t n2;
   SxArray<ssize_t> result;

   n1 = src1_.getSize ();
   n2 = src2_.getSize ();
   result.resize (n2);

   if (sorted_)  {
      // --- slightly faster version (n2 * log n1)
      for (i2=0; i2 < n2; i2++)  {
         result(i2) = SxSearch<T>::binary (src1_, src2_(i2));
      }
   }  else  {
      // --- simple and slow version (n2 * n1)
      for (i2=0; i2 < n2; i2++)  {
         result(i2) = -1;
         for (i1=0; i1 < n1; i1++)  {
            if (src1_(i1) == src2_(i2))  {
               result(i2) = i1;
               break;
            }
         }
      }
   }

   return result;
}

#endif /* _SX_SPARSE_UTIL_H_ */
