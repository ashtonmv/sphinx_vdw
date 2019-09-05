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
// Do not include this file directly. Use wrapper classes, e.g. 
// SxVector<T> or SxDirac<T> instead.
//----------------------------------------------------------------------------

#include <SxPrecision.h>
#include <SxRandom.h>
#include <SxVector3.h>
#include <SxMatrix3.h>
#include <SxList.h>
#include <SxArray.h>
#include <SxError.h>
#include <SxComplex.h>
#include <SxMathLib.h>
#include <SxIdx.h>
//#include <SxFFT.h>       /* data type 'fftw_complex' is defined here */
#include <SxTypeMapper.h>
#include <SxMemConsumer.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <typeinfo>

#include <SxLoopMPI.h>  // LoopMPI
#include <SxAllocCache.h>


#ifndef SXVEC
#  error "Undefined variable SXVEC"
#endif
#ifndef SXAUX
#  error "Undefined variable SXAUX"
#endif


#undef RS    // avoid HPUX <stdlib.h> <--> <veclib.h> problem
#include <stdlib.h>

// --- bitmask manipulators
#define MANAGE_MEMORY      0x01   /* if false destructor doesn't delete array */
#define SUB_INDEX          0x02   /* if true  operator= performs a real copy  */
#define IS_TRIANGULAR      0x04   /* only one triangle is given               */
#define UPPER_RIGHT        0x08   /* relevant only if IS_TRIANGULAR           */ 
#define WILDCARD_32        0x20
#define WILDCARD_64        0x40
#define REFERENCE          0x3c   /* used for reference counting              */

//using std::ios;
//using std::ios_base;
//using std::setw;
//using std::fixed;
//using namespace std;

// --- vectraits
#ifndef SX_VEC_TRAITS
#define SX_VEC_TRAITS
template <class T>
class SxVecTraits {};
#endif

/**
  \ingroup group_num
  \see     SxComplex<T>
  \see     SxError
  \author  Sixten Boeck
  */
template <class T> 
class SXVEC : public SxMemConsumer
{
   public:

      mutable MEMHANDLE      *handle;

      typedef T                            TypeMapper;


      /** The array of elements */
      mutable typename T::Type *elements;
      /** The iterator class has been introduced to access the elements of a 
          vector similar to C like pointers. For iterative methods or 
          algorithms you may speed up the code dramatically. In addition 
          usage of iterators is safer than pointers due to boundary checks.\\ 
          {\bf example:}                                                   \\
          \begin{verbatim}
          SxVector<float> theVector;
          ...
          SxVector<float>::Iterator it;
          for (it = theVector.begin(); it != theVector.end(); it++)  {
             *it = {whatever};
          }
          \end{verbatim}
          will be much faster than                                         \\
          \begin{verbatim}
          SxVector<float> theVector;
          ...
          for (ssize_t i = 0; i != theVector.getSize(); i++)  {
             theVector(i) = {whatever};
          }
          \end{verbatim}
          @see SxVector<T>::begin();
          @see SxVector<T>::end();
          @see SxVector<T>::operator();
       */
      class Iterator
      {
         public:
            inline Iterator (typename T::Type *ptrIn=NULL)  {
               ptr = ptrIn;
            }
            inline typename T::Type &operator* () const  {
               SX_CHECK (ptr);
               return *ptr;
            }
            inline bool operator== (typename SXVEC<T>::Iterator it)  {
               return ( it.ptr == ptr );
            }
            inline bool operator!= (typename SXVEC<T>::Iterator it)  {
               return ( it.ptr != ptr);
            }
            // postfix
            inline const typename SXVEC<T>::Iterator operator++ (int)  {
               typename SXVEC<T>::Iterator copy = Iterator (ptr);
               SX_CHECK (ptr);
               if (ptr)  ptr++;
               return copy;
            }
            // prefix
            inline typename SXVEC<T>::Iterator &operator++ ()  {
               SX_CHECK (ptr);
               if (ptr)  ptr++;
               return *this;
            }
            // postfix
            inline const typename SXVEC<T>::Iterator operator-- (int)  {
               typename SXVEC<T>::Iterator copy = Iterator (ptr);
               if (ptr)  ptr--;
               return copy;
            }
            // prefix
            inline typename SXVEC<T>::Iterator &operator-- ()  {
               if (ptr)  ptr--;
               return *this;
            }

            inline void operator+= (ssize_t i)  {
               SX_CHECK (ptr);
               ptr += i;
            }
            inline typename SXVEC<T>::Iterator operator+  (ssize_t i)  {
               SX_CHECK (ptr);
               typename SXVEC<T>::Iterator res;
               res.ptr = ptr + i;
               return res;
            }
            inline typename SXVEC<T>::Iterator &
               operator=  (const typename SXVEC<T>::Iterator &in)  
            {
               SX_CHECK (in.ptr);
               if ( this == &in )  return *this;
               ptr = in.ptr;
               return *this;
            }

            typename T::Type *ptr;
      };


     
      // --------------------------------------------------------------------
      /**@name Constructors and Destructors */
      //@{
      inline SXVEC ();
      inline SXVEC (const SXVEC<Int> &);
      inline SXVEC (const SXVEC<Float> &);
      inline SXVEC (const SXVEC<Double> &);
      inline SXVEC (const SXVEC<Complex8> &);
      inline SXVEC (const SXVEC<Complex16> &);
      inline SXVEC (const SXVEC<typename T::TReal> &,
                       const SXVEC<typename T::TReal> &);
      inline SXVEC (const SxList<typename T::Type> &);
      inline SXVEC (const SxArray<typename T::Type> &);
      /** \brief Set from imported stack */
             SXVEC (const SxStack<typename T::Type> &); 
      inline SXVEC (void *in, ssize_t nElem);
      inline SXVEC (const SxVector3<T> &);
      inline explicit SXVEC (ssize_t);
      inline          SXVEC (ssize_t, const typename T::Type &);
      inline virtual ~SXVEC ();
      //@}


      // --------------------------------------------------------------------
      /**@name Reference counting */
      //@{
      inline void ref (const SXVEC<T> &);
      inline void unref () const;
      inline void copy (const SXVEC<T> &in);
      inline SXVEC<T> getCopy () const;
      inline void init (ssize_t);
      //@}

      // --------------------------------------------------------------------
      /**@name Dereferencing */
      //@{
      inline       typename T::Type &operator() (ssize_t i);
      inline const typename T::Type &operator() (ssize_t i) const;
      inline       typename T::Type &operator() (SxAutoLoop &i);
      inline const typename T::Type &operator() (SxAutoLoop &i) const;
      inline SXVEC<T>  operator() (const SxIdx &) const;
      inline SXVEC<T>& reshape (ssize_t);
      inline SXVEC<T>& reshape (ssize_t, ssize_t);
      //@}

      // --------------------------------------------------------------------
      /**@name Algebraic operators */
      //@{
      inline SXVEC<T>  &operator= (const typename T::Type &);
      inline SXVEC<T>  &operator= (const SXVEC<T> &);
      inline SXVEC<T>  &operator= (const SxArray<ssize_t> &);
      /** \brief Import stack */
      inline SXVEC<T>  &operator= (const SxStack<typename T::Type> &);
      inline void          operator<<= (const SXVEC<T> &);
      inline void          operator<<  (const SXVEC<T> &);
      inline SXVEC<T>   operator- ()                    const;
      inline void          operator+= (const typename T::Type &);
      inline void          operator-= (const typename T::Type &);
      inline void          operator+= (const SXVEC<T> &);
      inline void          operator-= (const SXVEC<T> &);
      inline void          operator*= (const typename T::Type &);
      inline void          operator*= (const SXVEC<T> &);
      inline void          operator/= (const typename T::Type &);
      inline void          operator/= (const SXVEC<T> &);
      inline void          plus_assign_ax (const typename T::Type &a,
                                           const SXVEC<T> &x);
      //@}

      // --------------------------------------------------------------------
      /**@name Iterator functions */
      //@{
      inline typename SXVEC<T>::Iterator begin () const;
      inline typename SXVEC<T>::Iterator end ()   const;
      //@}

      // --------------------------------------------------------------------
      /**@name Utilities */
      //@{
      inline void set (const typename T::Type &);
      inline void set (const SXVEC<T> &);
      inline void set (const SxStack<typename T::Type> &, size_t stackSize, 
                       ssize_t offset = 0);
      inline typename T::Type minval (int *idx=NULL) const;
      inline typename T::Type maxval (int *idx=NULL) const;
      inline SXVEC<typename T::TReal> real () const;
      inline SXVEC<typename T::TReal> imag () const;
      inline SXVEC<T> sqr () const;
      inline SXVEC<T> cub () const;
      /** Computes \f$ \sum_{i=start}^{\rm end} v_i\f$ 
          \param end The (default) value of -1 is equivalent to getSize()-1;
        */
      inline typename T::Type sum (ssize_t startIdx=0, ssize_t endIdx=-1) const;
      inline typename T::Type sum (ssize_t startIdx, ssize_t endIdx, ssize_t step) const;

      inline typename T::Type chop () const;
      inline typename T::Type product () const;
      inline SXVEC<T> conj () const;
      inline ssize_t getSize () const;
       //\brief reshape + resize, looses data
      inline void reformat (ssize_t); 
      //\brief reshape + resize, looses data
      inline void reformat (ssize_t, ssize_t);
#     ifdef MSVC
      inline void resize (ssize_t, bool keep=false, 
                          const typename T::Type & fillVal = 
                          (T::Type)0.);
#     else
      inline void resize (ssize_t, bool keep=false, 
                          const typename T::Type & fillVal = 
                          (typename T::Type)0.);
#     endif
      inline SxVector3<T> toVector3 () const;
      // conversion to 3x3 matrix (if it is a 3x3 matrix)
      inline operator SxMatrix3<T> () const;
      inline bool isValid () const;
      inline SXVEC<T> abs () const;
      inline SXVEC<typename T::TReal> absSqr () const;
      /** Normalizes a vector. In case of matrices normalization per column will
          be computed */
      inline void normalize ();
      inline void randomize ();
      /// Sum absolute square of elements
      inline typename T::TReal::Type normSqr () const;
      /// Norm of vector
      inline typename T::TReal::Type norm () const;
      SxArray<ssize_t> getSortIdx () const;
      void sortByIdx (const SxArray<ssize_t> &sortIdx);
      SXVEC<T> getSorted (const SxArray<ssize_t> &sortIdx) const;

      typename T::Type integrate (const typename T::Type &step) const;
      /** Converts a matrix in row-major order (C order) to column major order
          (FORTRAN order).
          Note: 
          If SFHIngX is linked against a FORTRAN math library, such as 
          VecLIB or ESSL, you need to use ::toFBLASMajor instead!!! */
      SXVEC<T> toColMajor (const SxVector3<Int> &) const;
      /** Converts a matrix in row-major order (C order) to column major order
          (FORTRAN order).
          Note: 
          If SFHIngX is linked against a FORTRAN math library, such as 
          VecLIB or ESSL, you need to use ::toFBLASMajor instead!!! */
      SXVEC<T> toColMajor (int nx, int ny, int nz) const;
      /** Converts a matrix from FBLAS order to row-major order (C order) */
      SXVEC<T> toRowMajor (const SxVector3<Int> &) const;
      /** Converts a matrix from FBLAS order to row-major order (C order) */
      SXVEC<T> toRowMajor (int nx, int ny, int nz) const;
      /** Converts a matrix from row-major order (C order) to FORTRAN library
          complient transposed column-major order (FBLAS order). If SFHIngX
          is linked agains FORTRAN math libraries, such as VecLIB or ESSL, the
          x and z order is transposed internally. This doesn't affect the 
          programmers interface of the vector class. */
      SXVEC<T> toFBLASMajor (const SxVector3<Int> &) const;
      /** Converts a matrix from row-major order (C order) to FORTRAN library
          complient transposed column-major order (FBLAS order). If SFHIngX
          is linked agains FORTRAN math libraries, such as VecLIB or ESSL, the
          x and z order is transposed internally. This doesn't affect the 
          programmers interface of the vector class. */
      SXVEC<T> toFBLASMajor (int nx, int ny, int nz) const;
      SxArray<typename T::Type> toArray () const;
      //@}

      /**@name Functions working on matrices */
      //@{
      inline ssize_t nRows () const;
      inline ssize_t nCols () const;       
      inline       typename T::Type &operator() (ssize_t r, ssize_t c);
      inline const typename T::Type &operator() (ssize_t r, ssize_t c) const;
      // autoloop wrapper (I1, I2 can be int, ssize_t, or SxAutoLoop)
      template<class I1, class I2>
      inline       typename T::Type &operator() (const I1 &r,
                                                 const I2 &c);
      template<class I1, class I2>
      inline const typename T::Type &operator() (const I1 &r,
                                                 const I2 &c) const;
      inline SXVEC<T> row (ssize_t r) const;
      inline SXVEC<T> col (ssize_t c) const;
      inline SXVEC<T> colRef (ssize_t c);
      inline SXVEC<T> colRef (ssize_t c) const;
      inline SXVEC<T> colRef (SxAutoLoop &c)
      {
         c.setLimit (nCols ());
         return colRef (c.i);
      }
      inline SXVEC<T> colRef (SxAutoLoop &c) const
      {
         c.setLimit (nCols ());
         return colRef (c.i);
      }
      inline SXVEC<T> colRef (ssize_t c, ssize_t startIdx, ssize_t endIdx);
      inline SXVEC<T> getBlock (ssize_t rowOffset, ssize_t colOffset, ssize_t nRowsIn, ssize_t nColsIn) const;
      inline void     setBlock (const SXVEC<T> &block, ssize_t rowOffset, ssize_t colOffset);
      inline SXVEC<T> transpose () const;
      inline SXVEC<T> adjoint () const;
      inline SXVEC<T> inverse () const;
      inline SXVEC<T> solve (const SXVEC<T> &b) const;
      inline SXVEC<T> diag () const;
      inline typename T::Type trace () const;
      inline SXVEC<T> expand () const; 
      SXVEC<T> identity ();
      SXVEC<T> identity (const SXVEC<T> &diagIn) const;
      SXVEC<T> identity (const SXVEC<T> &diagIn, ssize_t n) const;
      SXVEC<T> choleskyDecomposition (enum UPLO = LowerLeft) const;
      /** \brief In-place rotation
          \param rotMat a square rotation matrix, to be applied from the right.
          This is a memory-friendly in-place rotation
          \f[
          M \leftarrow M U
          \f]
          where (this) M is a N x M matrix, and U (rotMat) is a M x M
          square matrix.
          \author C. Freysoldt
        */
      void rotate(const SXVEC<T> &rotMat);
      /** \brief Overlap computation

        This computes the overlap matrix x^T y. It is equivalent to
        \code
        x.adjoint () ^ y
        \endcode
        but avoids the explicit adjoint.
          
      */
      SXVEC<T> overlap (const SXVEC<T> &y) const;

      /** \brief Overlap computation of truncated vectors

        This computes the overlap matrix x^T y restricted to the
        first sumSize elements of each column.

      */
      SXVEC<T> overlap (const SXVEC<T> &y, ssize_t sumSize) const;
      //@}

      /** Print out a vector either in simple output format or in vector
          format (similar to \URL[Mathematica]{http://www.wolfram.com}) */
      inline void print (bool vectorForm=false) const;
      inline void print (std::ostream &, bool vectorForm=false) const;
      /** Recusive print */
      void  println (std::ostream &, bool) const;

#ifdef SXAUXHH
#  include SXAUXHH
#endif      

   protected:
      inline void              setSize    (ssize_t);
      inline typename T::Type *fastNew    (size_t) const;
      inline void              fastDelete (typename T::Type *) const;
      inline void    initMemHandle (const char &paramIn=MANAGE_MEMORY);
      inline void    initMemHandle (const ssize_t &sizeIn, 
                                    const ssize_t &nRowsIn=0, 
                                    const ssize_t &nColsIn=1,
                                    const char &paramIn=MANAGE_MEMORY);
      inline void    copyMemHandle (const MEMHANDLE *);
};

#ifdef NDEBUG
#   define VALIDATE_VECTOR(expr)         ((void)0);
#else
template<class T>
inline void VALIDATE_VECTOR (const SXVEC<T> &in)
{
   const typename T::Type *srcPtr = in.elements;
   
   for (ssize_t i=0; i < in.handle->size; i++, srcPtr++)  {
      SX_CHECK_NUM ((typename T::Type)(*srcPtr));
   }
}
#endif /* NDEBUG */      


#ifndef SKIP_DOXYGEN

#ifdef SXAUXHPP
#  include SXAUXHPP
#endif      

//------------------------------------------------------------------------------
// Standard constructor
//------------------------------------------------------------------------------
template<class T>
SXVEC<T>::SXVEC () : SxMemConsumer(), handle(NULL), elements(NULL)
{
   init (0);
}


//------------------------------------------------------------------------------
// Copy constructor. This differs from C++ standard! Not the elements will be
// copied but the pointers! To copy elements call copy()
//------------------------------------------------------------------------------
template<>
inline SXVEC<Int>::SXVEC (const SXVEC<Int> &in) 
   : SxMemConsumer(in), handle(NULL), elements(NULL)
{
   ref (in);
}


template<>
inline SXVEC<Float>::SXVEC (const SXVEC<Float> &in) 
   : SxMemConsumer(in), handle(NULL), elements(NULL)
{
   ref (in);
}


template<>
inline SXVEC<Double>::SXVEC (const SXVEC<Double> &in) 
   : SxMemConsumer(in), handle(NULL), elements(NULL)
{
   ref (in);
}


template<>
inline SXVEC<Complex8>::SXVEC (const SXVEC<Complex8> &in) 
   : SxMemConsumer(in), handle(NULL), elements(NULL)
{
   ref (in);
}


template<>
inline SXVEC<Complex16>::SXVEC (const SXVEC<Complex16> &in) 
   : SxMemConsumer(in), handle(NULL), elements(NULL)
{
   ref (in);
}


template<class T>
SXVEC<T>::SXVEC (const SxList<typename T::Type> &list) 
   : SxMemConsumer(), handle(NULL), elements(NULL)
{
   SX_CHECK (list.getSize() > 0, list.getSize());
   ssize_t n = list.getSize ();

   init (n);

   typename T::Type *ptr                          = elements;
   typename SxList<typename T::Type>::ConstIterator it = list.begin ();

   for (ssize_t i=0; i < n; i++)  *ptr++ = *it++;
}


template<class T>
SXVEC<T>::SXVEC (const SxArray<typename T::Type> &array) 
   : SxMemConsumer(), handle(NULL), elements(NULL)
{
   //SX_CHECK (array.getSize() > 0, array.getSize()); 0-size well possible
   ssize_t n = array.getSize ();

   init (n);

   typename T::Type *ptr = elements;
   for (ssize_t i=0; i < n; i++)  *ptr++ = array(i);
}

template<class T>
SXVEC<T>::SXVEC (const SxStack<typename T::Type> &stack) 
   : SxMemConsumer(), handle(NULL), elements(NULL)
{
   SX_CHECK (stack.getSize() > 0, stack.getSize());
   init (stack.getSize ());
   set (stack, getSize ());
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

#ifndef SX_GETNAN
#define SX_GETNAN
// --- real nan
template<class T>
inline T getNan ()
{
   return T(sqrt(-1.));
}

// --- complex nan = nan + i nan
template <>
inline SxComplex<double> getNan ()
{
   double rnan = getNan<double> ();
   return SxComplex<double> (rnan, rnan);
}

template <>
inline SxComplex<float> getNan ()
{
   float rnan = getNan<float> ();
   return SxComplex<float> (rnan, rnan);
}
#endif

template<class T>
typename T::Type *SXVEC<T>::fastNew (size_t nElem) const
{
   typename T::Type *res
     = (typename T::Type*)SxAllocation::getAligned (nElem * sizeof(typename T::Type), 16);

#  ifndef NDEBUG   
      // --- In the DEBUG mode we initialize explicitly the entire vector
      //     with NaNs (not-a-number). The macro VALIDATE_VECTOR can so
      //     identify uninitialized vectors.
      typename T::Type  nan = getNan<typename T::Type> ();
      typename T::Type *ptr = res;
      for (size_t i=0; i < nElem; i++)  *ptr++ = nan;
#  endif      
   return res;
}

template<class T>
void SXVEC<T>::fastDelete (typename T::Type *ptr) const
{
   if (handle->size && ptr)
      SxAllocation::retain (ptr, handle->size * sizeof(typename T::Type));
}


template<class T>
void SXVEC<T>::initMemHandle (const char &paramIn)
{
   handle = (MEMHANDLE *)SxAllocation::get (sizeof(MEMHANDLE));

   handle->refCounter = 0;
   handle->size       = 0;
   handle->nRows      = 0;
   handle->nCols      = 1;         // empty column vector
   handle->parameters = paramIn;
   handle->auxData.init ();
}

template<class T>
void SXVEC<T>::initMemHandle (const ssize_t &sizeIn,
                              const ssize_t &nRowsIn,
                              const ssize_t &nColsIn,
                              const char &paramIn)
{
   handle = (MEMHANDLE *)SxAllocation::get (sizeof(MEMHANDLE));

   handle->refCounter = 0;
   handle->size       = sizeIn;
   handle->nRows      = nRowsIn;
   handle->nCols      = nColsIn;
   handle->parameters = paramIn;
   handle->auxData.init ();
}


template<class T>
void SXVEC<T>::copyMemHandle (const MEMHANDLE *handleIn)
{
   if (handle)  
      SxAllocation::retain (handle, sizeof(MEMHANDLE));

   handle = (MEMHANDLE *)SxAllocation::get (sizeof(MEMHANDLE));

   handle->refCounter = 0;
   handle->size       = handleIn->size;
   handle->nRows      = handleIn->nRows;
   handle->nCols      = handleIn->nCols;
   handle->parameters = (char)((handleIn->parameters | MANAGE_MEMORY)
                        & ~(SUB_INDEX));
   handle->auxData.initPtr ();
   handle->auxData    = handleIn->auxData;
}



template<class T>
void SXVEC<T>::ref (const SXVEC<T> &in)
{
#  ifndef WIN32
   SX_CHECK (in.handle->size > 0 && in.elements, in.handle->size);
#  endif /* WIN32 */

   (in.handle->refCounter)++;
   handle   = in.handle;
   elements = in.elements;
}

//------------------------------------------------------------------------------
// Decrease the reference counter. If the counter is zero all dynamically
// allocated memory will be deallocated. A value of '-1' indicates that the
// memory will be managed from outside.
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::unref () const
{
   SX_CHECK (handle);

   if ( handle->refCounter == 0)  {  // destroy members

      if ( handle->parameters & MANAGE_MEMORY )  fastDelete (elements);
      SxAllocation::retain (handle, sizeof(MEMHANDLE));

   }  else  {                 // just decrease reference counter

      (handle->refCounter)--;

   }

   handle = NULL;
}



//------------------------------------------------------------------------------
// Copy constructors with type cast.
//------------------------------------------------------------------------------
template<>
inline SXVEC<Int>::SXVEC (const SXVEC<Float> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Int::Type  *destPtr;
   Float::Type *srcPtr = in.elements;
   //elements = destPtr = new float [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Int::Type)(*srcPtr++);
   }
}

template<>
inline SXVEC<Int>::SXVEC (const SXVEC<Double> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Int::Type  *destPtr;
   Double::Type *srcPtr = in.elements;
   //elements = destPtr = new double [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Int::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Int>::SXVEC (const SXVEC<Complex8> &in) : handle(NULL)
{
   SX_EXIT; // complex->int conversion ??
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Int::Type *destPtr;
   Complex8::Type *srcPtr = in.elements;
   //elements = destPtr   = new double [len];
   elements = destPtr     = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Int::Type)(*srcPtr++).re;
   }
}


template<>
inline SXVEC<Int>::SXVEC (const SXVEC<Complex16> &in) : handle(NULL)
{
   SX_EXIT; // complex->int conversion ??
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Int::Type *destPtr;
   Complex16::Type *srcPtr = in.elements;
   //elements = destPtr    = new double [len];
   elements = destPtr      = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Int::Type)(*srcPtr++).re;
   }
}



template<>
inline SXVEC<Float>::SXVEC (const SXVEC<Int> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Float::Type  *destPtr;
   Int::Type *srcPtr = in.elements;
   //elements = destPtr = new double [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Float::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Float>::SXVEC (const SXVEC<Double> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Float::Type  *destPtr;
   Double::Type *srcPtr = in.elements;
   //elements = destPtr = new double [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Float::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Float>::SXVEC (const SXVEC<Complex8> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Float::Type *destPtr;
   Complex8::Type *srcPtr = in.elements;
   //elements = destPtr   = new double [len];
   elements = destPtr     = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Float::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Float>::SXVEC (const SXVEC<Complex16> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Float::Type *destPtr;
   Complex16::Type *srcPtr = in.elements;
   //elements = destPtr    = new double [len];
   elements = destPtr      = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Float::Type)(*srcPtr++);
   }
}



template<>
inline SXVEC<Double>::SXVEC (const SXVEC<Int> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Double::Type *destPtr;
   Int::Type  *srcPtr = in.elements;
   //elements = destPtr = new double [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Double::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Double>::SXVEC (const SXVEC<Float> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Double::Type *destPtr;
   Float::Type  *srcPtr = in.elements;
   //elements = destPtr = new double [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Double::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Double>::SXVEC (const SXVEC<Complex8> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Double::Type   *destPtr;
   Complex8::Type *srcPtr = in.elements;
   //elements = destPtr   = new double [len];
   elements = destPtr     = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Double::Type)(*srcPtr++);
   }
}

template<>
inline SXVEC<Double>::SXVEC (const SXVEC<Complex16> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Double::Type *destPtr;
   Complex16::Type *srcPtr = in.elements;
   //elements = destPtr    = new double [len];
   elements = destPtr      = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Double::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Complex8>::SXVEC (const SXVEC<Int> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Complex8::Type *destPtr;
   Int::Type *srcPtr        = in.elements;
   //elements = destPtr = new SxComplex8 [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      destPtr->re = static_cast<float>(*srcPtr++);
      destPtr->im = 0.;
      destPtr++; 
   }
}


template<>
inline SXVEC<Complex8>::SXVEC (const SXVEC<Float> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Complex8::Type *destPtr;
   Float::Type *srcPtr        = in.elements;
   //elements = destPtr = new SxComplex8 [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      destPtr->re = *srcPtr++;
      destPtr->im = 0.;
      destPtr++; 
   }
}


template<>
inline SXVEC<Complex8>::SXVEC (const SXVEC<Double> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Complex8::Type *destPtr;
   Double::Type *srcPtr = in.elements;
   //elements = destPtr = new SxComplex8 [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      destPtr->re = (Complex8::Real)*srcPtr++;
      destPtr->im = (Complex8::Real)0.;
      destPtr++; 
   }
}

template<>
inline SXVEC<Complex8>::SXVEC (const SXVEC<Complex16> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Complex8::Type *destPtr;
   Complex16::Type *srcPtr = in.elements;
   //elements = destPtr    = new SxComplex8 [len];
   elements = destPtr      = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      destPtr->re = (Complex8::Real)srcPtr->re;
      destPtr->im = (Complex8::Real)srcPtr->im;
      destPtr++; srcPtr++;
   }
}


template<>
inline SXVEC<Complex16>::SXVEC (const SXVEC<Int> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Complex16::Type *destPtr;
   Int::Type *srcPtr  = in.elements;
   //elements = destPtr = new SxComplex16 [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Float::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Complex16>::SXVEC (const SXVEC<Float> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Complex16::Type *destPtr;
   Float::Type *srcPtr  = in.elements;
   //elements = destPtr = new SxComplex16 [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Double::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Complex16>::SXVEC (const SXVEC<Double> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Complex16::Type *destPtr;
   Double::Type *srcPtr = in.elements;
   //elements = destPtr = new SxComplex16 [len];
   elements = destPtr   = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Double::Type)(*srcPtr++);
   }
}


template<>
inline SXVEC<Complex16>::SXVEC (const SXVEC<Complex8> &in) : handle(NULL)
{
   SX_CHECK     (in.handle);
   SX_CHECK (in.handle->size > 0, in.handle->size);

   const ssize_t &len = in.handle->size;
   Complex16::Type *destPtr;
   Complex8::Type *srcPtr = in.elements;
   //elements = destPtr   = new SxComplex16 [len];
   elements = destPtr     = fastNew (len);
   copyMemHandle (in.handle);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = (Complex16::Type) (*srcPtr++);
   }
}


template<class T>
SXVEC<T>::SXVEC (const SXVEC<typename T::TReal> &reVec,
                 const SXVEC<typename T::TReal> &imVec)
   : handle (NULL)
{
   SX_CHECK      (reVec.handle && imVec.handle);
   SX_CHECK (reVec.handle->size > 0 && imVec.handle->size,
             reVec.handle->size,       imVec.handle->size);
   SX_CHECK (reVec.handle->size    ==  imVec.handle->size,
             reVec.handle->size,       imVec.handle->size);

   const ssize_t &len = reVec.handle->size;
   typename T::Cmplx  *destPtr;
   typename T::Real   *rePtr = reVec.elements;
   typename T::Real   *imPtr = imVec.elements;
   //elements = destPtr   = new SxComplex16 [len];
   elements = destPtr     = fastNew (len);
   initMemHandle (reVec.handle->size, reVec.handle->size, 1, MANAGE_MEMORY);

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = SxComplex<typename T::Real> (*rePtr++, *imPtr++);
   }
}


//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template<class T>
SXVEC<T>::SXVEC (void *in, ssize_t nElem) : handle(NULL)
{
   elements   = (typename T::Type *)in;

   // Exclude memory management in this instance. Array of elements will be
   // administrated from outside
   initMemHandle (nElem, nElem, nElem > 0 ? 1 : 0, 0);
}

template<class T>
SXVEC<T>::SXVEC (const SxVector3<T> &in) : handle(NULL)
{
   elements   = fastNew (3);
   initMemHandle (3, 3, 1, MANAGE_MEMORY);

   for (ssize_t i=0; i < 3; i++)  elements[i] = in.v[i];
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
SXVEC<T>::SXVEC (ssize_t nElements) : handle(NULL), elements (NULL)
{
   init (nElements);
}


template<class T>
SXVEC<T>::SXVEC (ssize_t nElements, const typename T::Type &in) 
   : handle(NULL), elements (NULL)
{
   init (nElements);
   set (in);
}


//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<class T>
SXVEC<T>::~SXVEC ()
{
   unref ();  
}



//------------------------------------------------------------------------------
// The actually copy constructor
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::copy (const SXVEC<T> &in)
{
   SX_CHECK      (handle);
// SX_CHECK (*size == 0 || *size == *in.size, *size, *in.size);

   const ssize_t &len = in.handle->size;

   if (handle->size > 0)  { reformat (handle->size); resize (0); }

   //elements = new T [len];
   elements = fastNew (len);
   copyMemHandle (in.handle);

   /*
   const typename T::Type *srcPtr  = in.elements;
   typename       T::Type *destPtr = elements;
   for (ssize_t i=0; i < len; i++)  {
      *destPtr++  = *srcPtr++;
   }
   */
   SX_CHECK (   in.elements >=    elements + len
             || elements    >= in.elements + len,
             in.elements - elements, len);
   memcpy (elements, in.elements, len * sizeof(typename T::Type));
}

template<class T>
SXVEC<T> SXVEC<T>::getCopy () const
{
   SXVEC<T> result;
   result.copy(*this);

   return result;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------



template<class T>
void SXVEC<T>::setSize (ssize_t newSize)
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size == 0, handle->size);
   SX_CHECK     (!elements);

   init (newSize);
}


//------------------------------------------------------------------------------
// Initialize helper variables and allocate memory.
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::init (ssize_t nElem)
{
   SX_CHECK     (!elements);
   SX_CHECK (nElem >= 0, nElem);

   if ( !nElem )  {
      elements   = NULL;
   }  else  {
      //elements   = new T [nElem];
      elements = fastNew (nElem);
   }

   initMemHandle (nElem, nElem, nElem > 0 ? 1 : 0, MANAGE_MEMORY);
   TRACK_MALLOC (*this, 1);
}


//------------------------------------------------------------------------------
// Reallocates memory for vector.
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::reformat (ssize_t newSize)
{
   SX_CHECK     (handle);
   SX_CHECK (newSize > 0, newSize);
   if (handle->size > 0)  reshape (handle->size); // resize expects 1d vectors
   resize  (newSize);
   reshape (newSize);
}


template<class T>
void SXVEC<T>::reformat (ssize_t newNRows, ssize_t newNCols)
{
   SX_CHECK      (handle);
   SX_CHECK (newNRows > 0 && newNCols > 0, newNRows, newNCols);
   ssize_t newSize = newNRows * newNCols;
   if (handle->size > 0)  reshape (handle->size); // resize expects 1d vectors
   resize  (newSize); 
   reshape (newNRows, newNCols);
}


template<class T>
void SXVEC<T>::resize (ssize_t newSize, bool keep, 
                       const typename T::Type &fillVal
                      )                       
{
   SX_CHECK (newSize >= 0, newSize);
   SX_CHECK (handle);

   SX_CHECK (handle->nCols < 2, handle->nCols);
                              // only allowed for vectors!
   SX_CHECK (handle->parameters & MANAGE_MEMORY); 
                              // handled from oudside!!!

   if (handle->size != newSize )  {
      SX_CHECK (handle->refCounter == 0, handle->refCounter);  
                                 // If not fulfilled: Pointers to this
                                 // instance point to nowhere! 

      if (keep)  {
         ssize_t i, min = newSize < handle->size ? newSize : handle->size;
         //T *temp    = new T [newSize];
         typename T::Type *temp    = fastNew (newSize);
         typename T::Type *srcPtr  = elements;
         typename T::Type *destPtr = temp;

         // --- copy leading part of vector
         for (i=0; i < min; i++) 
            *destPtr++ = *srcPtr++;

         // --- fill rest of vector with fillVal
         if (newSize > min)  {
            destPtr = temp + min;   
            for (i=min; i < newSize; i++)
               *destPtr++ = fillVal;
         }

         //delete [] elements;
         fastDelete (elements);

         elements = temp;

      }  else  {

         //delete [] elements;
         fastDelete (elements);

         if (newSize == 0)
            elements   = NULL;
         else  {
            //elements   = new T [newSize];
            elements   = fastNew (newSize);
         }
      }
      handle->size = newSize;
   }
   handle->nRows = handle->size;
   handle->nCols = 1;
   TRACK_MALLOC (*this, 1);
}


template<class T>
SxVector3<T> SXVEC<T>::toVector3 () const
{
   SX_CHECK (getSize() == 3, getSize());

   return SxVector3<T> (elements[0], elements[1], elements[2]);
}

template<class T>
inline SXVEC<T>::operator SxMatrix3<T> () const
{
   SX_CHECK(nCols () == 3 && nRows () == 3, nCols (), nRows ());
   return SxMatrix3<T> (elements[0], elements[3], elements[6],
                        elements[1], elements[4], elements[7],
                        elements[2], elements[5], elements[8]);
}


template<class T>
bool SXVEC<T>::isValid () const
{
   ssize_t len = getSize();
   if (len < 0)  {
      printf ("ERROR: Vector size is corrupt!\n");
      return false;
   }
   if (len == 0)  {
      printf ("ERROR: Vector is empty!\n");
      return false;
   }

   const typename T::Type *srcPtr = elements;
   for (ssize_t i=0; i < len; i++, srcPtr++)  {
      if (!::isValid(*srcPtr))  {
         printf ("ERROR: Vector is corrupt!\n");
         return false;
      }
   }
   return true;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
SXVEC<T> SXVEC<T>::operator- () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = elements;
#  ifdef USE_OPENMP
      if ( len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
            private(smpIdx) schedule(static)
            for (smpIdx=0; smpIdx < len; smpIdx++)
               destPtr[smpIdx] = -(srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)
         *destPtr++  = -(*srcPtr++);
   }

   return res;
}


template<class T>
inline SXVEC<T> SXVEC<T>::sqr () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<T> res = SXVEC<T> (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   typename T::Type *srcPtr  = elements;
   typename T::Type *destPtr = res.elements;
#  ifdef USE_OPENMP      
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = srcPtr[smpIdx] * srcPtr[smpIdx];
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++, ++srcPtr)
         *destPtr++ = *srcPtr * *srcPtr;
   }

   VALIDATE_VECTOR (res);
   return res;
}


template<>
inline SXVEC<Complex8> SXVEC<Complex8>::sqr () const
{
   SX_EXIT; // should be absSqr () ?!
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Complex8> res = SXVEC<Complex8> (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   Complex8::Type *srcPtr  = elements;
   Complex8::Type *destPtr = res.elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = (srcPtr[smpIdx])^2;
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)  
         *destPtr++ = (*srcPtr++).absSqr ();
   }

   VALIDATE_VECTOR (res);
   return res;
}


template<>
inline SXVEC<Complex16> SXVEC<Complex16>::sqr () const
{
   SX_EXIT; // should be absSqr () ?!
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Complex16> res = SXVEC<Complex16> (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   Complex16::Type *srcPtr  = elements;
   Complex16::Type *destPtr = res.elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = (srcPtr[smpIdx])^2;
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)
         *destPtr++ = (*srcPtr++).absSqr ();
   }

   VALIDATE_VECTOR (res);
   return res;
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
SXVEC<T> SXVEC<T>::abs () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<T> res = SXVEC<T> (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   typename T::Type *srcPtr  = elements;
   typename T::Type *destPtr = res.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = ::fabs((double)srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)  
         *destPtr++ = ::fabs((double)*srcPtr++);
   }

   VALIDATE_VECTOR (res);
   return res;
}



template<>
inline SXVEC<Int::TReal> SXVEC<Int>::absSqr () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Int::TReal> res (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   Int::Type *srcPtr  = elements;
   Int::Real *destPtr = res.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = srcPtr[smpIdx] * srcPtr[smpIdx]; 
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++, ++srcPtr)  
         *destPtr++ = *srcPtr * *srcPtr; 
   }

   VALIDATE_VECTOR (res);
   return res;
}

template<>
inline SXVEC<Float::TReal> SXVEC<Float>::absSqr () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Float::TReal> res (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   Float::Type *srcPtr  = elements;
   Float::Real *destPtr = res.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = srcPtr[smpIdx] * srcPtr[smpIdx]; 
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++, ++srcPtr)
         *destPtr++ = *srcPtr * *srcPtr;
   }

   VALIDATE_VECTOR (res);
   return res;
}


template<>
inline SXVEC<Double::TReal> SXVEC<Double>::absSqr () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Double::TReal> res (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   Double::Type *srcPtr  = elements;
   Double::Real *destPtr = res.elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = srcPtr[smpIdx] * srcPtr[smpIdx]; 
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++, ++srcPtr)  
         *destPtr++ = *srcPtr * *srcPtr; 
   }

   VALIDATE_VECTOR (res);
   return res;
}


template<>
inline SXVEC<Complex8::TReal> SXVEC<Complex8>::absSqr () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Complex8::TReal> res (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   Complex8::Type *srcPtr  = elements;
   Complex8::Real *destPtr = res.elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = (srcPtr[smpIdx])^2; 
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)  
         *destPtr++ = (*srcPtr++).absSqr (); 
   }

   VALIDATE_VECTOR (res);
   return res;
}



template<>
inline SXVEC<Complex16::TReal> SXVEC<Complex16>::absSqr () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Complex16::TReal> res (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   Complex16::Type *srcPtr  = elements;
   Complex16::Real *destPtr = res.elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = (srcPtr[smpIdx])^2; 
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)
         *destPtr++ = (*srcPtr++).absSqr (); 
   }

   VALIDATE_VECTOR (res);
   return res;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
SXVEC<T> SXVEC<T>::cub () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t len = handle->size;
   SXVEC<T> res = SXVEC<T> (len);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;

   typename T::Type *srcPtr  = elements;
   typename T::Type *destPtr = res.elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         typename T::Type *valPtr;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx,valPtr) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)  {
            valPtr = &srcPtr[smpIdx];
            destPtr[smpIdx] = *valPtr * *valPtr * *valPtr;
         }
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++, ++srcPtr)
            *destPtr++ = *srcPtr * *srcPtr * *srcPtr;
      }

   VALIDATE_VECTOR (res);
   return res;
}


// note: see also 'operator() (ssize_t) const'
template<class T>
typename T::Type &SXVEC<T>::operator() (ssize_t i)
{ 
   SX_CHECK (handle);
   SX_CHECK (i >= 0 && i < handle->size, i, handle->size);
   return elements[i]; 
}


// note: see also 'operator() (ssize_t)'
template<class T>
const typename T::Type &SXVEC<T>::operator() (ssize_t i) const
{ 
   SX_CHECK (handle);
   SX_CHECK (i >= 0 && i < handle->size, i, handle->size);
   return elements[i]; 
}

template<class T>
typename T::Type &SXVEC<T>::operator() (SxAutoLoop &i)
{ 
   SX_CHECK (handle);
   SX_CHECK (i >= 0 && i < handle->size, i, handle->size);
   i.setLimit (handle->size);
   return elements[i]; 
}


// note: see also 'operator() (ssize_t)'
template<class T>
const typename T::Type &SXVEC<T>::operator() (SxAutoLoop &i) const
{ 
   SX_CHECK (handle);
   SX_CHECK (i >= 0 && i < handle->size, i, handle->size);
   i.setLimit (handle->size);
   return elements[i]; 
}

template<class T>
SXVEC<T> SXVEC<T>::operator() (const SxIdx &idx) const
{
   SX_CHECK      (handle);
   SX_CHECK (idx.end < handle->size, idx.end, handle->size);
   // temporary objects are not allowed to use with subindeces
   // SX_CHECK      (handle->parameters & MANAGE_MEMORY);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) ); // not yet implemented

   SXVEC<T> res;

   res.elements        = elements + idx.start;
   res.handle->size    = res.handle->nRows = idx.end - idx.start + 1;
   res.handle->nCols   = 1;
   res.handle->auxData = handle->auxData;

  // --- memory will be managed by 'this'
  res.handle->parameters = SUB_INDEX;
//  delete res.refCounter;
//  res.refCounter = refCounter;

   return res;
}



template<class T>
SXVEC<T>& SXVEC<T>::reshape (ssize_t vecLen)
{
   SX_CHECK      (handle);
   SX_CHECK (handle->size == vecLen, handle->size, vecLen);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) ); // not yet implemented
// shape->removeAll()   will be done in SxList<T>::operator=
   handle->size  = handle->nRows = vecLen;
   handle->nCols = 1;
   return *this;
}


template<class T>
SXVEC<T>& SXVEC<T>::reshape (ssize_t newNRows, ssize_t newNCols)
{
   SX_CHECK       (handle);
   SX_CHECK (   (     handle->size == newNRows * newNCols
             && !(handle->parameters & IS_TRIANGULAR))
             || (   newNRows == newNCols
             && handle->size == newNRows*(newNRows+1)/2
             && (handle->parameters & IS_TRIANGULAR)
             ), handle->size, newNRows, newNCols);
// shape->removeAll()   will be done in SxList<T>::operator=
   //*shape = Shape() << newNRows << newNCols;
   handle->nRows = newNRows;
   handle->nCols = newNCols;
   return *this;
}




//------------------------------------------------------------------------------
// Assignment operator (assign vector with scalar)
//------------------------------------------------------------------------------
template<class T>
SXVEC<T> &SXVEC<T>::operator= (const typename T::Type &in)
{
   SX_CHECK (handle);
   SX_CHECK (elements && handle->size > 0, handle->size);
   SX_CHECK     ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   SX_CHECK_NUM (in);

   typename T::Type *destPtr = elements;

   const ssize_t &len = handle->size;
   for (ssize_t i=0; i < len; i++)  {
//    elements[i] = in;
      *destPtr++  = in;
   }

   return *this;
}


//------------------------------------------------------------------------------
// Assignment operator
//------------------------------------------------------------------------------
template<class T>
SXVEC<T> &SXVEC<T>::operator= (const SXVEC<T> &in)
{
   // It's a time critical routine! So avoid statements like 'a = a'!
   // Normally here would be a line like
   //    if ( this == &in )  return *this;
   // would be.
   SX_CHECK (this != &in);

   // --- If LHS is a sub indexed vector, elements will be copied!!! (slower)
   if (   (handle->parameters & SUB_INDEX)
       && (!in.getSize() == 0))   // --- expression "subIdx = SxVec<T> ()" 
                                  //     can release SUB_INDEX again
   {
      SX_CHECK (handle->size == in.handle->size,
                handle->size, in.handle->size);
      SX_CHECK  (handle->refCounter == 0,
                 handle->refCounter);  // otherwise we change more than 'this'

      const ssize_t &len = in.handle->size;
      typename T::Type *srcPtr  = in.elements;
      typename T::Type *destPtr = elements;
      handle->auxData = in.handle->auxData;
      for (ssize_t i=0; i < len; i++)
         *destPtr++ = *srcPtr++;

   // --- copy of pointers (fast)
   }  else  {

      unref ();
//    SX_CHECK (*refCounter == 0); // otherwise the references point to nowhere!

      (in.handle->refCounter)++;
      handle   = in.handle;
      elements = in.elements;
   }
   return *this;
}

// workaround until size_t is used.
template<>
inline SXVEC<Int> &SXVEC<Int>::operator= (const SxArray<ssize_t> &in)
{
   SX_CHECK (getSize() == in.getSize(), getSize(), in.getSize());
   for (ssize_t i=0; i < getSize(); ++i)  {
      SX_CHECK (in(i) < 0x0fffffff);
      elements[i] = (int)(in(i));
   }
   return *this;
}

template<class T>
SXVEC<T>& SXVEC<T>::operator= (const SxStack<typename T::Type> &stack)
{
   ssize_t n = (ssize_t) stack.getSize ();
   if (getSize () != n) resize(n);
   set (stack, n);
   return *this;
}

template<class T>
void SXVEC<T>::set (const SxStack<typename T::Type> &stack, size_t stackSize,
                    ssize_t offset)
{
   SX_CHECK (offset >= 0, offset);
   SX_CHECK ((size_t)getSize () >= size_t(offset) + stackSize,
             getSize (), offset, stackSize);
   // This is a trick to get access to SxStack's protected member
   // function exportStack
   // we can't use friends, because we would need to map
   // double->Double for that
   class StackAccess : public SxStack<typename T::Type> {
      public: 
         void exportStack (typename T::Type *target, size_t n) const  {
            SxStack<typename T::Type>::exportStack (target, n);
         }
   };
   const StackAccess &stackAlias = static_cast<const StackAccess&>(stack);
   stackAlias.exportStack(elements + offset, stackSize);
}

template<class T>
void SXVEC<T>::operator<<= (const SXVEC<T> &in)
{
   SX_CHECK      (handle);
   SX_CHECK      (in.handle->size);
   SX_CHECK (handle->size == in.handle->size, handle->size, in.handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   set (in);
   VALIDATE_VECTOR (*this);
}


template<class T>
void SXVEC<T>::operator<< (const SXVEC<T> &in)
{
   SX_CHECK      (handle);
   SX_CHECK      (in.handle->size);
   SX_CHECK (handle->size >= in.handle->size, handle->size, in.handle->size);
// SX_CHECK  (handle->nCols <  2, shape->nCols);  // not for matrices!
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   handle->auxData = in.handle->auxData;

   ssize_t i, min  = minimum (getSize(), in.getSize());
   ssize_t max     = maximum (getSize(), in.getSize());

   typename T::Type *srcPtr  = in.elements;
   typename T::Type *destPtr = elements;
   // --- copy leading part
   for (i=0; i < min; i++)
      *destPtr++  = *srcPtr++;

   // --- fill rest with zeros
   if (min < max)  {
      destPtr = &elements[min];
      for (i=min; i < max; i++)  
         *destPtr++ = (typename T::Type)0.;
   }
   VALIDATE_VECTOR (*this);
}


//------------------------------------------------------------------------------
// Operator '*='
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::operator*= (const typename T::Type &in)
{
   SX_CHECK      (handle);
   SX_CHECK  (handle->size > 0, handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   SX_CHECK_NUM  (in);

   typename T::Type *destPtr = elements;
   const ssize_t &len = handle->size;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         const typename T::Type *valPtr = &in;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]  *= *valPtr;
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++  *= in;
      }

   VALIDATE_VECTOR (*this);
}


template<class T>
void SXVEC<T>::operator*= (const SXVEC<T> &in)
{
   SX_CHECK      (handle);
   SX_CHECK      (in.handle->size);
   SX_CHECK (in.handle->size == handle->size && handle->size > 0,
             handle->size, in.handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   typename T::Type *destPtr = elements;
   typename T::Type *srcPtr  = in.elements;
   const ssize_t &len = handle->size;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] *= srcPtr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ *= *srcPtr++;
      }
      
   VALIDATE_VECTOR (*this);
}

//------------------------------------------------------------------------------
// Operator '/='
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::operator/= (const typename T::Type &s)
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);
   SX_CHECK (::fabs((Double::Type)s) > 1e-10, (Double::Type)s);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   SX_CHECK_DIV (s);

   typename T::Type d = (typename T::Type)1.0 / s;

   operator*= (d);
}

template<class T>
void SXVEC<T>::operator/= (const SXVEC<T> &in)
{
   SX_CHECK      (handle);
   SX_CHECK      (in.handle);
   SX_CHECK (in.handle->size == handle->size && handle->size > 0,
             handle->size, in.handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   typename T::Type *destPtr = elements;
   typename T::Type *srcPtr  = in.elements;
   const ssize_t &len = handle->size;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] /= srcPtr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ /= *srcPtr++;
      }
      
   VALIDATE_VECTOR (*this);
}

//------------------------------------------------------------------------------
// Operator '+='
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::operator+= (const typename T::Type &in)
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   SX_CHECK_NUM  (in);

   typename T::Type *destPtr = elements;
   const ssize_t &len = handle->size;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const typename T::Type *valPtr = &in;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] += *valPtr;
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++  += in;
      }

   VALIDATE_VECTOR (*this);
}

template<class T>
void SXVEC<T>::operator-= (const typename T::Type &in)
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   SX_CHECK_NUM  (in);

   typename T::Type *destPtr = elements;
   const ssize_t &len = handle->size;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const typename T::Type *valPtr = &in;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] -= *valPtr;
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++  -= in;
      }

   VALIDATE_VECTOR (*this);
}

template<class T>
void SXVEC<T>::operator+= (const SXVEC<T> &in)
{
   SX_CHECK      (handle);
   SX_CHECK (handle->size == in.handle->size, handle->size, in.handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   typename T::Type *srcPtr  = in.elements;
   typename T::Type *destPtr = elements;
   const ssize_t &len = handle->size;

#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] += srcPtr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++  += *srcPtr++;
      }

   VALIDATE_VECTOR (*this);
}


//------------------------------------------------------------------------------
// Y += a*X
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::plus_assign_ax (const typename T::Type &a, const SXVEC<T> &x)
{
   SX_CHECK (getSize() == x.getSize(), getSize(), x.getSize());
   SX_CHECK  (getSize() > 0, getSize());
   
   axpy (elements, a, x.elements, (int)getSize());
}


//------------------------------------------------------------------------------
// Operator '-='
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::operator-= (const SXVEC<T> &in)
{
   SX_CHECK      (handle);
   SX_CHECK (handle->size == in.handle->size, handle->size, in.handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) ); // not yet implemented

   const ssize_t &len = handle->size;
   typename T::Type *srcPtr  = in.elements;
   typename T::Type *destPtr = elements;

#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] -= srcPtr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++  -= *srcPtr++;
      }

   VALIDATE_VECTOR (*this);
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::set (const typename T::Type &value)
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);
   SX_CHECK_NUM (value);

   typename T::Type *destPtr = elements;
   const ssize_t &len = handle->size;
   for (ssize_t i=0; i < len; i++)
      *destPtr++ = value;
}


template<class T>
void SXVEC<T>::set (const SXVEC<T> &in)
{
   SX_CHECK      (handle);
   SX_CHECK      (in.handle->size);
   SX_CHECK  (handle->size > 0, handle->size);
   SX_CHECK (handle->size == in.handle->size, handle->size, in.handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) ); // not yet implemented

   handle->auxData = in.handle->auxData;

   typename T::Type *srcPtr  = in.elements;
   typename T::Type *destPtr = elements;
   const ssize_t &len = handle->size;
   /*
   for (ssize_t i=0; i < len; i++)
      *destPtr++  = *srcPtr++;
   */
   memcpy (destPtr, srcPtr, len * sizeof(typename T::Type));

   VALIDATE_VECTOR (*this);
}


/// \todo: SMP parallelization
template<class T>
typename T::Type SXVEC<T>::minval (int *idx) const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);
   typename T::Type *ptr = elements;
   typename T::Type  res = *ptr++;
   ssize_t pos = 0;
   const ssize_t &len = handle->size;
   for (ssize_t i=1; i < len; i++, ptr++)  {
      if ( *ptr < res )  { res = *ptr; pos = i; }
   }
   SX_CHECK_NUM (res);
   if (idx)  *idx = (int)pos;
   return res;
}

/// \todo: SMP parallelization
template<class T>
typename T::Type SXVEC<T>::maxval (int *idx) const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);
   typename T::Type *ptr = elements;
   typename T::Type  res = *ptr++;
   ssize_t pos = 0;
   const ssize_t &len = handle->size;
   for (ssize_t i=1; i < len; i++, ptr++)  {
      if ( *ptr > res )   {  res = *ptr; pos = i; }
   }
   SX_CHECK_NUM (res);
   if (idx)  *idx = (int)pos;
   return res;
}


template<class T>
SXVEC<typename T::TReal> SXVEC<T>::real () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);
   
   const ssize_t &len = handle->size;
   SXVEC<typename T::TReal> res (len);
   res.handle->nRows = handle->nRows;
   res.handle->nCols = handle->nCols;
   res.handle->auxData = handle->auxData;
   typename T::Real *destPtr = res.elements;
   typename T::Type *srcPtr  = elements;
   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = toReal(*srcPtr++);
   }

   VALIDATE_VECTOR (res);
   return res;
}

template<class T>
SXVEC<typename T::TReal> SXVEC<T>::imag () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);
   
   const ssize_t &len = handle->size;
   SXVEC<typename T::TReal> res (len);
   res.handle->nRows = handle->nRows;
   res.handle->nCols = handle->nCols;
   res.handle->auxData = handle->auxData;
   typename T::Real *destPtr = res.elements;
   typename T::Type *srcPtr  = elements;

   for (ssize_t i=0; i < len; i++)  {
      *destPtr++ = toImag(*srcPtr++);
   }

   VALIDATE_VECTOR (res);
   return res;
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<>
inline SXVEC<Float> SXVEC<Float>::conj () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;

   SXVEC<Float> res (len);
   res.handle->nRows    = handle->nRows;
   res.handle->nCols    = handle->nCols;
   res.handle->auxData  = res.handle->auxData;
   Float::Type *destPtr = res.elements;
   Float::Type *srcPtr  = elements;
   for (ssize_t i=0; i < len; i++)
      *destPtr++ = *srcPtr++;
   
   VALIDATE_VECTOR (res);
   return res;
}


template<>
inline SXVEC<Double> SXVEC<Double>::conj () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Double> res (len);
   res.handle->nRows     = handle->nRows;
   res.handle->nCols     = handle->nCols;
   res.handle->auxData   = handle->auxData;
   Double::Type *destPtr = res.elements;
   Double::Type *srcPtr  = elements;
   for (ssize_t i=0; i < len; i++)
      *destPtr++ = *srcPtr++;
   
   VALIDATE_VECTOR (res);
   return res;
}



template<>
inline SXVEC<Complex8> SXVEC<Complex8>::conj () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Complex8> res (len);
   res.handle->nRows       = handle->nRows;
   res.handle->nCols       = handle->nCols;
   res.handle->auxData     = handle->auxData;
   Complex8::Type *destPtr = res.elements;
   Complex8::Type *srcPtr  = elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = srcPtr[smpIdx].conj();
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)
         *destPtr++ = (*srcPtr++).conj();
   }
   
   VALIDATE_VECTOR (res);
   return res;
}


template<>
inline SXVEC<Complex16> SXVEC<Complex16>::conj () const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size > 0, handle->size);

   const ssize_t &len = handle->size;
   SXVEC<Complex16> res (len);
   res.handle->nRows        = handle->nRows;
   res.handle->nCols        = handle->nCols;
   res.handle->auxData      = handle->auxData;
   Complex16::Type *destPtr = res.elements;
   Complex16::Type *srcPtr  = elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = srcPtr[smpIdx].conj();
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)
         *destPtr++ = (*srcPtr++).conj();
   }

   VALIDATE_VECTOR (res);
   return res;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::normalize ()
{
   SX_CHECK     (handle && elements);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   typename T::Type nrm2, c;
   if (handle->nCols > 1) {  // matrices: per column vector normalization
      ssize_t ic, nr = nRows(), nc = nCols(), offset = 0;
      for (ic=0; ic < nc; ic++, offset+=nr)  {
         nrm2 = norm2 (elements+offset, static_cast<int>(nr));
         SX_CHECK_DIV (nrm2);
         c = (typename T::Type)(1.)/nrm2;
         scale (elements+offset, c, static_cast<int>(nr));
      }
   }  else  {                    // vector
      const ssize_t &len = handle->size;
      nrm2 = norm2 (elements, static_cast<int>(len));
      SX_CHECK_DIV (nrm2);
      c = (typename T::Type)(1.)/nrm2;
      scale (elements, c, static_cast<int>(len));
   }

   VALIDATE_VECTOR (*this);
}

template <class T>
typename T::TReal::Type SXVEC<T>::norm () const
{
   SX_CHECK (handle);
   SX_CHECK (handle->size > 0);
   return norm2(elements, static_cast<int>(handle->size)); // this is the norm
}

template <class T>
typename T::TReal::Type SXVEC<T>::normSqr () const
{
   SX_CHECK (handle);
   SX_CHECK (handle->size > 0);
   typename T::TReal::Type nrm = norm2(elements,static_cast<int>(handle->size));
   return nrm*nrm;
}

//------------------------------------------------------------------------------
// Presets values with uniformly distributed random numbers.
// Suppose random number generator has been initialized before calling.
//------------------------------------------------------------------------------
template<>
inline void SXVEC<Float>::randomize ()
{
   SX_CHECK     (handle && elements);
   SX_CHECK (handle->size > 0, handle->size);

   Float::Type *destPtr = elements;
   const ssize_t &len = handle->size;
   for (ssize_t i=0; i < len; i++)  {
//    elements[i] = (Float::Type)SxRandom::get();
      *destPtr++  = (Float::Type)SxRandom::get();
   }
   normalize ();
}


template<>
inline void SXVEC<Double>::randomize ()
{
   SX_CHECK     (handle && elements);
   SX_CHECK (handle->size > 0, handle->size);

   Double::Type *destPtr = elements;
   const ssize_t &len = handle->size;
   for (ssize_t i=0; i < len; i++)  {
//    elements[i] = (Double::Type)SxRandom::get();
      *destPtr++  = (Double::Type)SxRandom::get();
   }
   normalize ();
}

//------------------------------------------------------------------------------
// Specific SxComplex8 version
// Presets values with uniformly distributed random numbers.
// Suppose random number generator has been initialized before calling.
//------------------------------------------------------------------------------
template<>
inline void SXVEC<Complex8>::randomize ()
{
   SX_CHECK     (handle && elements);
   SX_CHECK (handle->size > 0, handle->size);

   Complex8::Type *destPtr = elements;
   const ssize_t &len = handle->size;
   for (ssize_t i=0; i < len; i++)  {
//    elements[i] = SxComplex8 (SxRandom::get(), SxRandom::get());
      *destPtr++  = Complex8::Type ((float)SxRandom::get(), 
                                    (float)SxRandom::get());
   }
   normalize ();
}

//------------------------------------------------------------------------------
// Specific SxComplex8 version
// Presets values with uniformly distributed random numbers.
// Suppose random number generator has been initialized before calling.
//------------------------------------------------------------------------------
template<>
inline void SXVEC<Complex16>::randomize ()
{
   SX_CHECK     (handle && elements);
   SX_CHECK (handle->size > 0, handle->size);

   Complex16::Type *destPtr = elements;
   const ssize_t &len = handle->size;
   for (ssize_t i=0; i < len; i++)  {
//    elements[i] = SxComplex16 (SxRandom::get(), SxRandom::get());
      *destPtr++  = Complex16::Type (SxRandom::get(), SxRandom::get());
   }
   normalize ();
}



//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
ssize_t SXVEC<T>::getSize () const
{
   SX_CHECK (handle);
   return handle->size;
}






//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
SxArray<ssize_t> SXVEC<T>::getSortIdx () const
{
   SxArray<ssize_t> sortIdx (getSize());

   if (getSize() == 1)  {
      sortIdx(0) = 0; 
      return sortIdx;
   }

   ::sort (sortIdx.elements, elements, getSize());
   return sortIdx;
}


//------------------------------------------------------------------------------
// Sort elements of vector.
//------------------------------------------------------------------------------
template<class T>
void SXVEC<T>::sortByIdx (const SxArray<ssize_t> &sortIdx)
{
   SX_CHECK  (sortIdx.getSize() >= 1, sortIdx.getSize());
   SX_CHECK (getSize() == sortIdx.getSize(),
             getSize(),   sortIdx.getSize());
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   if (sortIdx.getSize() == 1)  return;

   ssize_t i;
   const ssize_t &len = handle->size;

   //T *temp = new T [*size];
   typename T::Type *temp = fastNew (len);
   // --- copy array
   for (i=0; i < len; i++) 
      temp[i] = elements[i];
   // --- rearrange elements
   for (i=0; i < len; i++) 
      elements[i] = temp[sortIdx.elements[i]];

   //delete [] temp;
   fastDelete (temp);
}

template<class T>
SXVEC<T> SXVEC<T>::getSorted (const SxArray<ssize_t> &idx) const
{
   SX_CHECK (getSize () == idx.getSize (), getSize (), idx.getSize ());
   ssize_t n = getSize ();
   SXVEC<T> res(n);
   res.handle->nRows   = handle->nRows;
   res.handle->nCols   = handle->nCols;
   res.handle->auxData = handle->auxData;
#ifdef USE_OPENMP
#  pragma omp parallel for if (n > sxChunkSize)
#endif
   for (ssize_t i = 0; i < n; ++i)
      res.elements[i] = elements[idx(i)];
   VALIDATE_VECTOR (res);
   return res;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<>
inline Int::Type SXVEC<Int>::sum (ssize_t startIdx, ssize_t endIdx) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1 : endIdx;
   Int::Type *ptr = elements;
   Int::Type  res = (Int::Type)0;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(ptr,startIdx,lastIdx) \
         private(smpIdx) schedule(static)           \
         reduction(+:res)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx++)  {
             res += elements[smpIdx];
         }
      } else
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i++)  {
            res += *ptr++;
         }
      }
         
   SX_CHECK_NUM (res);
   return res;
}

template<>
inline Float::Type SXVEC<Float>::sum (ssize_t startIdx, ssize_t endIdx) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1 : endIdx;
   Float::Type *ptr = elements;
   Float::Type  res = (Float::Type)0.;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(ptr,startIdx,lastIdx) \
         private(smpIdx) schedule(static)           \
         reduction(+:res)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx++)  {
             res += elements[smpIdx];
         }
      } else
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i++)  {
            res += *ptr++;
         }
      }
         
   SX_CHECK_NUM (res);
   return res;
}

template<>
inline Double::Type SXVEC<Double>::sum (ssize_t startIdx, ssize_t endIdx) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1: endIdx;
   Double::Type *ptr = elements;
   Double::Type  res = (Double::Type)0.;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(ptr,startIdx,lastIdx) \
         private(smpIdx) schedule(static)           \
         reduction(+:res)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx++)  {
            res += elements[smpIdx];
         }
      } else
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i++)  {
            res += *ptr++;
         }
      }
   SX_CHECK_NUM (res);
   return res;
}

template<>
inline Complex8::Type SXVEC<Complex8>::sum (ssize_t startIdx, ssize_t endIdx) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1 : endIdx;
   Complex8::Type *ptr = elements;
   Complex8::Type  res = (Complex8::Type)0.;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         Complex8::Real resRe=0., resIm=0.;
         ssize_t smpIdx;
#        pragma omp parallel for shared(startIdx,lastIdx) \
         private(smpIdx) schedule(static)       \
         reduction(+:resRe,resIm)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx++)  {
            resRe += elements[smpIdx].re;
            resIm += elements[smpIdx].im;
         }
         res = Complex8::Type (resRe, resIm);
      } else
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i++)  {
            res += *ptr++;
         }
      }
   SX_CHECK_NUM (res);
   return res;
}

template<>
inline Complex16::Type SXVEC<Complex16>::sum (ssize_t startIdx, ssize_t endIdx) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1 : endIdx;
   Complex16::Type *ptr = elements;
   Complex16::Type  res = (Complex16::Type)0.;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         Complex16::Real resRe=0., resIm=0.;
         ssize_t smpIdx;
#        pragma omp parallel for shared(startIdx,lastIdx) \
         private(smpIdx) schedule(static)       \
         reduction(+:resRe,resIm)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx++)  {
            resRe += elements[smpIdx].re;
            resIm += elements[smpIdx].im;
         }
         res = Complex16::Type (resRe, resIm);
      } else
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i++)  {
            res += *ptr++;
         }
      }
   SX_CHECK_NUM (res);
   return res;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<>
inline Float::Type SXVEC<Float>::sum (ssize_t startIdx, ssize_t endIdx, ssize_t step) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      (step != 0);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1 : endIdx;
   Float::Type *ptr = elements;
   Float::Type  res = (Float::Type)0.;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(ptr,startIdx,lastIdx,step) \
         private(smpIdx) schedule(static)                \
         reduction(+:res)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx+=step)  {
            res += elements[smpIdx];
         }
      } else 
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i+=step, ptr+=step)  {
            res += *ptr;
         }
      }
   SX_CHECK_NUM (res);
   return res;
}

template<>
inline Double::Type SXVEC<Double>::sum (ssize_t startIdx, ssize_t endIdx, ssize_t step) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      (step != 0);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1: endIdx;
   Double::Type *ptr = elements;
   Double::Type  res = (Double::Type)0.;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(ptr,startIdx,lastIdx,step) \
         private(smpIdx) schedule(static)                \
         reduction(+:res)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx+=step)  {
            res += elements[smpIdx];
         }
      } else 
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i+=step, ptr+=step)  {
            res += *ptr;
         }
      }
   SX_CHECK_NUM (res);
   return res;
}

template<>
inline Complex8::Type SXVEC<Complex8>::sum (ssize_t startIdx, ssize_t endIdx, ssize_t step) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      (step != 0);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1 : endIdx;
   Complex8::Type *ptr = elements;
   Complex8::Type  res = (Complex8::Type)0.;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         Complex8::Real resRe=0., resIm=0.;
         ssize_t smpIdx;
#        pragma omp parallel for shared(startIdx,lastIdx,step) \
         private(smpIdx) schedule(static)            \
         reduction(+:resRe,resIm)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx+=step)  {
            resRe += elements[smpIdx].re;
            resIm += elements[smpIdx].im;
         }
         res = Complex8::Type (resRe, resIm);
      } else
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i+=step, ptr+=step)  {
            res += *ptr;
         }
      }
   SX_CHECK_NUM (res);
   return res;
}

template<>
inline Complex16::Type 
SXVEC<Complex16>::sum (ssize_t startIdx, ssize_t endIdx, ssize_t step) const
{
   SX_CHECK (startIdx >= 0 && (startIdx <= endIdx || endIdx==-1),
             startIdx, endIdx);
   SX_CHECK (endIdx >= -1 && endIdx < handle->size, endIdx, handle->size);
   SX_CHECK      (step != 0);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t lastIdx = (endIdx == -1) ? handle->size - 1 : endIdx;
   Complex16::Type *ptr = elements;
   Complex16::Type  res = (Complex16::Type)0.;
#  ifdef USE_OPENMP   
      ssize_t len = lastIdx - startIdx;
      if (len > sxChunkSize)  {
         Complex16::Real resRe=0., resIm=0.;
         ssize_t smpIdx;
#        pragma omp parallel for shared(startIdx,lastIdx,step) \
         private(smpIdx) schedule(static)            \
         reduction(+:resRe,resIm)
         for (smpIdx=startIdx; smpIdx<=lastIdx; smpIdx+=step)  {
            resRe += elements[smpIdx].re;
            resIm += elements[smpIdx].im;
         }
         res = Complex16::Type (resRe, resIm);
      } else
#  endif /* USE_OPENMP */      
      {
         ptr += startIdx;
         for (ssize_t i=startIdx; i<=lastIdx; i+=step, ptr+=step)  {
            res += *ptr;
         }
      }
   SX_CHECK_NUM (res);
   return res;
}

#ifdef USE_LOOPMPI
template<>
inline void SxLoopMPI::sum (SXVEC<Int> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum (inout.elements, inout.elements, inout.getSize ());
}

template<>
inline void SxLoopMPI::sum (SXVEC<Double> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum (inout.elements, inout.elements, inout.getSize ());
}

template<>
inline void SxLoopMPI::sum (SXVEC<Complex16> &inout)
{
   // actual implementation is in SxLoopMPI
   SxLoopMPI::sum ((double*)inout.elements, (double*)inout.elements,
                   inout.getSize () * 2);
}

template <>
inline void SxLoopMPI::bcast (SXVEC<Double> &inout, int source)
{
   SxLoopMPI::bcast(inout.elements, inout.getSize(), source);
}

template <>
inline void SxLoopMPI::bcast (SXVEC<Complex16> &inout, int source)
{
   SxLoopMPI::bcast((double*)inout.elements, 2 * inout.getSize(), source);
}
#else
template<> inline void SxLoopMPI::sum (SXVEC<Int> &) { }
template<> inline void SxLoopMPI::sum (SXVEC<Double> &) { }
template<> inline void SxLoopMPI::sum (SXVEC<Complex16> &) { }
template<> inline void SxLoopMPI::bcast (SXVEC<Double> &, int) { }
template<> inline void SxLoopMPI::bcast (SXVEC<Complex16> &, int) { }
#endif

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
typename T::Type SXVEC<T>::chop () const
{
   SX_CHECK (getSize() == 1, getSize());

   return *elements;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template<class T>
typename T::Type SXVEC<T>::product () const
{
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   typename T::Type  res = (typename T::Type)1.;
   typename T::Type *ptr = elements;
   const ssize_t &len = handle->size;
   for (ssize_t i=0; i < len; i++)  {
//    res += elements[i];
      res *= *ptr++;
   }
   SX_CHECK_NUM (res);
   return res;
}


//------------------------------------------------------------------------------
// Simpson integration on logarithmic mesh
//------------------------------------------------------------------------------
template<class T>
typename T::Type SXVEC<T>::integrate (const typename T::Type &step) const
{
   SX_CHECK (handle);
   SX_CHECK (handle->size > 4, handle->size); //check for minimal input values
   SX_CHECK ((Double::Type)step > 0.0, (Double::Type)step);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   typename T::Type evenpartresult=0,result;

   ssize_t n=handle->size,nend=n-1,nt=n/2;
   
	if ( 2 * nt == n )  {	//if number of input grippoints is even, compute the 
									//upper values by the 3/8-Simpson rule
		evenpartresult=(3. * step/8.) * (elements[n-1] + 3. * elements[n-2]
												+ 3. * elements[n-3] + elements[n-4]);
		nend-=3;    
	}
   
	result = (step/3.0) * (4.0 * sum(1, nend-1, 2) + 2.0 * sum(2, nend-2, 2) 
								 + elements[0] + elements[nend]) + evenpartresult;

   SX_CHECK_NUM (result);
   return result;
}

template<class T>
SXVEC<T> SXVEC<T>::toRowMajor (const SxVector3<Int> &dim) const
{
   return toRowMajor (dim(0), dim(1), dim(2));
}

template<class T>
SXVEC<T> SXVEC<T>::toRowMajor (int nx, int ny, int nz) const
{
   SX_CHECK (nx > 0 && ny > 0 && nz > 0, nx, ny, nz);
   SX_CHECK       (handle);
   SX_CHECK  (handle->size == nx*ny*nz, handle->size, nx*ny*nz);
   
   SXVEC<T> res (handle->size);
   res.handle->auxData = handle->auxData;
   typename T::Type *destPtr = res.elements;

   int idx, x, y, z;
   for (x=0; x < nx; x++)  {
      for (y=0; y < ny; y++)  {
         for (z=0; z < nz; z++)  {
//          idx   = (z*ny+y)*nx+x;  // fortran style
            idx   = (x*ny+y)*nz+z;  // c
            *destPtr++ = elements[idx];
         }
      }
   }
   return res;

}

template<class T>
SXVEC<T> SXVEC<T>::toColMajor (const SxVector3<Int> &dim) const
{
   return toColMajor (dim(0), dim(1), dim(2));
}

template<class T>
SXVEC<T> SXVEC<T>::toColMajor (int nx, int ny, int nz) const
{
   SX_CHECK (nx > 0 && ny > 0 && nz > 0, nx, ny, nz);
   SX_CHECK       (handle);
   SX_CHECK  (handle->size == nx*ny*nz, handle->size, nx*ny*nz);
   
   SXVEC<T> res (handle->size);
   res.handle->auxData = handle->auxData;
   typename T::Type *destPtr = res.elements;

   int idx, x, y, z;
   for (z=0; z < nz; z++)  {
      for (y=0; y < ny; y++)  {
         for (x=0; x < nx; x++)  {
            //idx = (z*ny+y)*nx+x;  // fortran
            idx   = (x*ny+y)*nz+z;  // c
            *destPtr++ = elements[idx];
         }
      }
   }

   return res;

}


template<class T>
SXVEC<T> SXVEC<T>::toFBLASMajor (const SxVector3<Int> &dim) const
{
   return toFBLASMajor (dim(0), dim(1), dim(2));
}


template<class T>
SXVEC<T> SXVEC<T>::toFBLASMajor (int nx, int ny, int nz) const
{
   SX_CHECK (nx > 0 && ny > 0 && nz > 0, nx, ny, nz);
   SX_CHECK       (getSize());
   SX_CHECK  (handle->size == nx*ny*nz, handle->size, nx*ny*nz);

   SXVEC<T> res (handle->size);
   typename T::Type *dstPtr = res.elements;

   int idx, x, y, z;
   for (z=0; z < nz; z++)  {
      for (y=0; y < ny; y++)  {
         for (x=0; x < nx; x++)  {
            idx = (z*ny+y)*nx+x;
            *dstPtr++ = elements[idx];
         }
      }
   }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
SxArray<typename T::Type> SXVEC<T>::toArray() const
{
   SxArray<typename T::Type> res(handle->size);

   const ssize_t &len = handle->size;
   for (ssize_t i = 0; i < len; i++)  res(i) = elements[i];

   return res;
}




// ==============================================================================
//
//                          Functions working on matrices
// ==============================================================================
template<class T>
ssize_t SXVEC<T>::nCols () const
{
   SX_CHECK (handle);
   return handle->nCols;
}


template<class T>
ssize_t SXVEC<T>::nRows () const
{
   SX_CHECK (handle);
   return handle->nRows;
}


// Note: see also 'operator() (ssize_t, ssize_t) const'
template<class T>
typename T::Type &SXVEC<T>::operator() (ssize_t r, ssize_t c)
{
   SX_CHECK     (handle);

   SX_CHECK (r >= 0 && r < nRows(), r, nRows());
   SX_CHECK (c >= 0 && c < nCols(), c, nCols());

   // return this->elements[r * nCols() + c];   // row major order
   ssize_t idx = 0;

   if ( !(handle->parameters & IS_TRIANGULAR))  {
      idx = c * nRows() + r;
   }  else  {  // triangular matrix
      if (handle->parameters & UPPER_RIGHT) {
         SX_CHECK (r <= c, r, c);
      //   idx = c*(c+1)/2 + r;
         idx = ( c + r*(2* nRows() - (r+1))/2);
      }  else  {
         SX_EXIT; // lower triangle not allowed 
               // (can't compute conj() of a reference as return variable!!!
         SX_CHECK (r >= c, r, c);
         idx = r*(r+1)/2 + c; // lower triangle
      }
   }
   SX_CHECK (idx < handle->size, idx, handle->size);
   return elements[idx];   // column major order
}



// Note: see also 'operator() (ssize_t, ssize_t)'
template<class T>
const typename T::Type &SXVEC<T>::operator() (ssize_t r, ssize_t c) const
{
   SX_CHECK     (handle);

   SX_CHECK (r >= 0 && r < nRows(), r, nRows());
   SX_CHECK (c >= 0 && c < nCols(), c, nCols());

   // return this->elements[r * nCols() + c];   // row major order
   ssize_t idx = 0;

   if ( !(handle->parameters & IS_TRIANGULAR))  {
      idx = c * nRows() + r;
   }  else  {  // triangular matrix
      if (handle->parameters & UPPER_RIGHT) {
         SX_CHECK (r <= c, r, c);
      //   idx = c*(c+1)/2 + r;
         idx = ( c + r*(2* nRows() - (r+1))/2);
      }  else  {
         SX_EXIT; // lower triangle not allowed 
               // (can't compute conj() of a reference as return variable!!!
         SX_CHECK (r >= c, r, c);
         idx = r*(r+1)/2 + c; // lower triangle
      }
   }
   SX_CHECK (idx < handle->size, idx, handle->size);
   return elements[idx];   // column major order
}

template<class T>
template<class I1, class I2>
inline
typename T::Type &SXVEC<T>::operator() (const I1 &r,
                                        const I2 &c)
{
   SX_CHECK (handle);
   SxAutoLoop::setLimit (r, nRows ());
   SxAutoLoop::setLimit (c, nCols ());
   return operator() ((ssize_t)r, (ssize_t)c);
}

template<class T>
template<class I1, class I2>
inline
const typename T::Type &SXVEC<T>::operator() (const I1 &r,
                                              const I2 &c) const
{
   SX_CHECK (handle);
   SxAutoLoop::setLimit (r, nRows ());
   SxAutoLoop::setLimit (c, nCols ());
   return operator() ((ssize_t)r, (ssize_t)c);
}

template<class T>
SXVEC<T> SXVEC<T>::row (ssize_t r) const
{
   SX_CHECK     (handle);

   SX_CHECK (r >= 0 && r < nRows(), r, nRows());

   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t nCol   = nCols();
   ssize_t stride = nRows();

   SXVEC<T> res (nCol);
   res.handle->auxData = handle->auxData;
   
   typename SXVEC<T>::Iterator srcIt  = this->begin();
   typename SXVEC<T>::Iterator destIt = res.begin();
   srcIt += r;

   for (ssize_t c=0; c < nCol; c++, srcIt += stride, destIt++)  {
      *destIt = *srcIt;
   }

   res.handle->nRows = handle->nCols;
   res.handle->nCols = 1;

   return res;
}


template<class T>
SXVEC<T> SXVEC<T>::col (ssize_t c) const
{
   SX_CHECK     (handle);

   SX_CHECK (c >= 0 && c < nCols(), c, nCols());

   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t nr   = nRows();

   SXVEC<T> res (nr);
   res.handle->auxData = handle->auxData;
   
   typename SXVEC<T>::Iterator srcIt  = this->elements + (c*nr);
   typename SXVEC<T>::Iterator destIt = res.begin();

   for (ssize_t r=0; r < nr; r++, srcIt++, destIt++)  {
      *destIt = *srcIt;
   }

   res.handle->nRows = nr;
   res.handle->nCols = 1;

   return res;
}




template<class T>
SXVEC<T> SXVEC<T>::colRef (ssize_t c)
{
   SX_CHECK     (handle);

   SX_CHECK (c >= 0 && c < nCols(), c, nCols());

   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   SXVEC<T> res = SXVEC<T> (this->elements + (c * nRows()), nRows());
   res.handle->auxData = handle->auxData;
   res.handle->nCols   = 1;

   return res;
}

template<class T>
SXVEC<T> SXVEC<T>::colRef (ssize_t c) const
{
   SX_CHECK     (handle);

   SX_CHECK (c >= 0 && c < nCols(), c, nCols());

   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   SXVEC<T> res = SXVEC<T> (this->elements + (c * nRows()), nRows());

   res.handle->auxData = handle->auxData;
   res.handle->nCols   = 1;

   return res;
}

template<class T>
SXVEC<T> SXVEC<T>::colRef (ssize_t c, ssize_t startIdx, ssize_t endIdx)
{
   SX_CHECK      (handle);
   SX_CHECK (startIdx <= endIdx, startIdx, endIdx);
   SX_CHECK  (startIdx >= 0, startIdx);
   SX_CHECK (endIdx   <  handle->size, endIdx, handle->size);

   SX_CHECK (c >= 0 && c < nCols(), c, nCols());
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   SXVEC<T> res = SXVEC<T> (   this->elements + startIdx 
                                  + (c * nRows()), endIdx+1-startIdx);

   res.handle->nCols   = 1;
   res.handle->auxData = handle->auxData;

   return res;
}

template<class T>
SXVEC<T> SXVEC<T>::getBlock (ssize_t rowOffset, ssize_t colOffset, ssize_t nRowsIn, ssize_t nColsIn) const
{
   SX_CHECK (rowOffset + nRowsIn <= nRows(), rowOffset + nRowsIn, nRows());
   SX_CHECK (colOffset + nColsIn <= nCols(), colOffset + nColsIn, nCols());
   const SXVEC<T> &matrix = *this;
   SXVEC<T> result;
   result.reformat(nRowsIn, nColsIn);
   for (ssize_t iRow = 0; iRow < nRowsIn; iRow++)  {
      for (ssize_t iCol = 0; iCol < nColsIn; iCol++)  {
         result(iRow,iCol) = matrix(rowOffset+iRow,colOffset+iCol);
      }
   }
   result.handle->auxData = handle->auxData;
   return result;
}

template<class T>
void SXVEC<T>::setBlock (const SXVEC<T> &block, ssize_t rowOffset, ssize_t colOffset)
{
   SXVEC<T> &matrix = *this;
   SX_CHECK (rowOffset + block.nRows() <= nRows(), rowOffset + block.nRows(), nRows());
   SX_CHECK (colOffset + block.nCols() <= nCols(), colOffset + block.nCols(), nCols());
   for (ssize_t iRow = 0; iRow < block.nRows(); iRow++)  {
      for (ssize_t iCol = 0; iCol < block.nCols(); iCol++)  {
         matrix(rowOffset + iRow, colOffset + iCol) = block(iRow, iCol);
      }
   } 
}


template<class T>
SXVEC<T> SXVEC<T>::transpose () const
{
   SX_CHECK     (handle);
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   const ssize_t &len = handle->size;
   SXVEC<T> res (len);

   // --- transpose shape information
   res.handle->nRows   = handle->nCols;
   res.handle->nCols   = handle->nRows;
   res.handle->auxData = handle->auxData;
   if (res.handle->nRows == 0)  {  // input vector is a 'unshaped' row vector
      res.handle->nRows = 1;
   }
//   *res.shape =  typename SXVEC<T>::Shape () 
//              << (*this->shape)(1) << (*this->shape)(0);
   const ssize_t &newNCols = handle->nRows;


   // --- transpose elements
   typename SXVEC<T>::Iterator srcIt  = begin(); 
   typename SXVEC<T>::Iterator destIt = res.begin();

   for (ssize_t i=0, x=0; i < len; i++, destIt++)  {
      *destIt = *srcIt;
      srcIt += newNCols; x += newNCols;
      if ( x >= len )  {
         x -= len; x++;
         srcIt = begin(); srcIt += x;
      }
   }

   VALIDATE_VECTOR (res);
   return res;
}


template<>
inline SXVEC<Float> SXVEC<Float>::adjoint () const
{
   SX_CHECK (handle);
   SX_CHECK ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   return transpose ();
}


template<>
inline SXVEC<Double> SXVEC<Double>::adjoint () const
{
   SX_CHECK (handle);
   SX_CHECK ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   return transpose ();
}


template<>
inline SXVEC<Complex8> SXVEC<Complex8>::adjoint () const
{
   SX_CHECK (handle);
   SX_CHECK ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   return conj().transpose ();
}


template<>
inline SXVEC<Complex16> SXVEC<Complex16>::adjoint () const
/*
{
   SX_CHECK (handle);
   SX_CHECK ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   return conj().transpose ();
}
*/
{
   SX_CHECK (handle);
   SX_CHECK ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t len = handle->size;
   SXVEC<Complex16> res (len);

   // --- transpose shape information
   res.handle->nRows   = handle->nCols;
   res.handle->nCols   = handle->nRows;
   res.handle->auxData = handle->auxData;
   if (res.handle->nRows == 0)  {  // input vector is a 'unshaped' row vector
      res.handle->nRows = 1;
   }
   ssize_t newNCols = handle->nRows;

   // --- transpose elements
   SXVEC<Complex16>::Iterator srcIt  = begin(); 
   SXVEC<Complex16>::Iterator destIt = res.begin();

   for (ssize_t i=0, x=0; i < len; i++, destIt++)  {
      (*destIt).re =  (*srcIt).re;
      (*destIt).im = -(*srcIt).im;
      srcIt += newNCols; x += newNCols;
      if ( x >= len )  {
         x -= len; x++;
         srcIt = begin(); srcIt += x;
      }
   }

   VALIDATE_VECTOR (res);
   return res;
}



template<class T>
SXVEC<T> SXVEC<T>::inverse () const
{
   SX_CHECK      (handle);

   
   SXVEC<T> res; 

   if (nRows () == nCols ())  {   // for sym. matrices
      res.copy (*this);  // the input array will be overwritten by LAPACK
      
      if ( !(handle->parameters & IS_TRIANGULAR ) )  { // --- general matrix

         matInverse (res.elements, (int)nRows (), (int)nCols ());

      }  else  {                                // --- triangular matrix
         SX_CHECK (handle->parameters & UPPER_RIGHT);

         // --- change UpperRight -> LowerLeft
         matInverseTri (res.elements, (int)nRows (), LowerLeft);

         res.handle->parameters |= (IS_TRIANGULAR | UPPER_RIGHT);
      }
   }  else  {// non-square, assume full rank matrix
      if (nRows () < nCols ())  {// right inverse defined
         res = transpose ().conj () ^ (*this ^ transpose ().conj ()).inverse ();
      }  else  {// left inverse defined
         res = (transpose ().conj () ^ *this).inverse () ^ transpose ().conj ();
      }
    //res.reshape (nCols (), nRows ()); //no. cols <-> no. rows AUTO?
   }

   VALIDATE_VECTOR (res);
   return res;
}

template<class T>
SXVEC<T> SXVEC<T>::solve (const SXVEC<T> &b) const
{
   // Check input
   SX_CHECK(this->getSize() > 0);
   SX_CHECK(b.getSize() > 0);
   SX_CHECK(this->nRows() == b.nRows(), this->nRows(), b.nRows());

   SXVEC<T> result(this->nCols() * b.nCols());
   result.reshape(this->nCols(),b.nCols());

   SXVEC<T> bCopy = b.getCopy();
   SXVEC<T> thisCopy = this->getCopy();

   solveLinEq(thisCopy.elements, int(thisCopy.nRows()), int(thisCopy.nCols()), 
                 bCopy.elements, int(bCopy.nCols()));

   for(ssize_t c = 0; c < bCopy.nCols(); c++)   {
      for(ssize_t r = 0; r < this->nCols(); r++)   {
         result(r,c) = bCopy(r,c);
      }
   }

   VALIDATE_VECTOR(result);
   return result;
}

template<class T>
SXVEC<T> SXVEC<T>::diag () const
{
   SX_CHECK      (handle);
   SX_CHECK (nRows() == nCols(), nRows(), nCols());

   const ssize_t &len = handle->nRows;
   SXVEC<T> res (len);
   res.handle->auxData = handle->auxData;
   typename SXVEC<T>::Iterator srcIt  = begin();
   typename SXVEC<T>::Iterator destIt = res.begin();
   if ( !(handle->parameters & IS_TRIANGULAR) )  {
      ssize_t offset = len+1;
      for (ssize_t i=0; i < len; i++, srcIt+=offset, destIt++)  
         *destIt = *srcIt;
   }  else  {  // triangular matrix
      ssize_t offset = len+1;
      for (ssize_t i=0; i < len; i++, srcIt+=offset, destIt++)  {
         *destIt = *srcIt;
         offset--;
      }
   }
   VALIDATE_VECTOR (res);

   return res;
}


template<class T>
typename T::Type SXVEC<T>::trace () const
{
   SX_CHECK ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   return diag().sum();
}


template<class T>
SXVEC<T> SXVEC<T>::expand () const
{
   SX_CHECK      (handle);
   SX_CHECK      (handle->parameters & IS_TRIANGULAR ); // for sym.matrices only
   SX_CHECK  (handle->size > 0, handle->size);
   SX_CHECK (nRows() == nCols(), nRows(), nCols());

   const ssize_t &n = handle->nRows;
   ssize_t r, c;

   SXVEC<T> res (n*n); res.reshape (n, n);
   res.handle->auxData = handle->auxData;
   typename SXVEC<T>::Iterator destLeftIt;   // goes column wise
   typename SXVEC<T>::Iterator destRightIt;  // goes row wise
   typename SXVEC<T>::Iterator srcIt = begin();

   for (r=0; r < n; r++)  {
      destLeftIt  = res.begin() + r*n + r;
      destRightIt = destLeftIt;
      for (c=r; c < n; c++)  {
         *destRightIt =  *srcIt;            destRightIt +=n ;
         *destLeftIt  = conjugate(*srcIt);  destLeftIt++;
         srcIt++;
      }
   }

   // --- resulting matrix is not triangular any more
   res.handle->parameters 
      = (res.handle->parameters | IS_TRIANGULAR) ^ IS_TRIANGULAR;
  
   VALIDATE_VECTOR (res); 
   return res;
}
template<>
inline SXVEC<Complex16> SXVEC<Complex16>::expand () const
{
   SX_CHECK      (handle);
   SX_CHECK      (handle->parameters & IS_TRIANGULAR );  // only for symmetric matrices
   SX_CHECK  (handle->size > 0, handle->size);
   SX_CHECK (nRows() == nCols(), nRows(), nCols());

   ssize_t r, c;
   const ssize_t &n = handle->nRows;

   SXVEC<Complex16> res (n*n); res.reshape (n, n);
   res.handle->auxData = handle->auxData;
   SXVEC<Complex16>::Iterator destLeftIt;   // goes column wise
   SXVEC<Complex16>::Iterator destRightIt;  // goes row wise
   SXVEC<Complex16>::Iterator srcIt = begin();

   for (r=0; r < n; r++)  {
      destLeftIt  = res.begin() + r*n + r;
      destRightIt = destLeftIt;
      for (c=r; c < n; c++)  {
         *destRightIt =  *srcIt;          destRightIt += n;
         *destLeftIt  = (*srcIt).conj();  destLeftIt++;
         srcIt++;
      }
   }
   
   // --- resulting matrix is not triangular any more
   res.handle->parameters 
      = (res.handle->parameters | IS_TRIANGULAR) ^ IS_TRIANGULAR;

   VALIDATE_VECTOR (res);
   return res;
}

template<class T>
SXVEC<T> SXVEC<T>::identity ()
{
   SX_CHECK      (handle);
   SX_CHECK  (handle->size > 0, handle->size);

   SXVEC<T> res;

   if ( !(handle->parameters & IS_TRIANGULAR) )  {
      SX_CHECK (nRows() == nCols(), nRows(), nCols());

      const ssize_t &len = handle->size;
      const ssize_t &n   = handle->nRows;
//      SXVEC<T> res (n*n);
      res.resize (n*n);
      res.handle->auxData = handle->auxData;
      res.reshape (n, n);
      res.set (0.);
      typename T::Type *destPtr = res.elements;
      const typename T::Type  c = (typename T::Type)1.;

      const ssize_t stride = nCols() + 1;
      for (ssize_t i=0; i < len; i+=stride, destPtr+=stride)  {
         *destPtr = c;
      }
   }  else  {  // triangular matrix

//      SXVEC<T> res ( n*(n+1)/2 );
      const ssize_t &len = handle->size;
      const ssize_t n    = (ssize_t)(0.5 * sqrt(1. + 8.*(double)len) - 0.5);
                       // from x = n (n+1) / 2
      SX_CHECK (len == n*(n+1)/2, len, n);
      res.resize (len);
      res.handle->auxData = handle->auxData;
      res.set (0.);
      res.handle->parameters |= (char)((handle->parameters & UPPER_RIGHT)
                                + IS_TRIANGULAR);
      res.reshape (n, n);
      const typename T::Type  c = (typename T::Type)1.;

      for (ssize_t i = 0; i < n; i++)  res(i,i) = c;
      /* TODO: fix UPPER_RIGHT-bug in vector class
         the following is for TRUE col-major order

      typename T::Type *destPtr = res.elements;
      if (handle->parameters & UPPER_RIGHT)  {
         ssize_t stride = 1;
         for (ssize_t i = 0; i < n; i++, stride++, destPtr += stride)  {
            *destPtr = c;
         }
      }  else  {  // lower left

         ssize_t stride = n+1;
         for (ssize_t i = 0; i < n; i++, stride--, destPtr += stride)  {
            *destPtr = c;
         }
      }
      */
   }
      
   VALIDATE_VECTOR (res);
   return res;
}

template<class T>
SXVEC<T> SXVEC<T>::identity (const SXVEC<T> &diagIn) const
{
   return identity (diagIn, diagIn.getSize());
}


template<class T>
SXVEC<T> SXVEC<T>::identity (const SXVEC<T> &diagIn, ssize_t n) const
{
   SX_CHECK  (n > 0, n);
   SX_CHECK  (diagIn.getSize() > 0,  diagIn.getSize());
   SX_CHECK (diagIn.getSize() >= n, diagIn.getSize(), n);

   SXVEC<T> res;
   
   if ( !(handle->parameters & IS_TRIANGULAR) )  {
      ssize_t stride = n + 1;
//      SXVEC<T> res (n*n);
      res.resize (n*n);
      res.reshape (n, n);
      res.handle->auxData = handle->auxData;
      res.set (0.); 
      typename T::Type *destPtr = res.elements;
      typename T::Type *srcPtr  = diagIn.elements;

      const ssize_t &len = res.handle->size;
      for (ssize_t i=0; i < len; i+=stride, destPtr+=stride)  {
         *destPtr = *srcPtr++;
      }
   }  else  {  // triangular matrix

//      SXVEC<T> res ( n*(n+1)/2 );
      res.resize ( n*(n+1)/2 );
      res.handle->auxData = handle->auxData;
      res.handle->parameters |= (char)((handle->parameters & UPPER_RIGHT)
                                + IS_TRIANGULAR);
      res.reshape (n,n);
      res.set (0.);
      
      for (ssize_t i = 0; i < n; i++)  res(i,i) = diagIn(i);
      /* TODO: fix UPPER_RIGHT-bug in vector class
         the following is for TRUE col-major order

      typename T::Type *destPtr = res.elements;
      typename T::Type *srcPtr  = diagIn.elements;

      if (handle->parameters & UPPER_RIGHT)  {
         ssize_t stride = 1;
         for (ssize_t i = 0; i < n; i++, stride++, destPtr += stride)  {
            *destPtr = *srcPtr++;
         }
      }  else  {  // lower left

         ssize_t stride = n+1;
         for (ssize_t i = 0; i < n; i++, stride--, destPtr += stride)  {
            *destPtr = *srcPtr++;
         }
      }
      */
   }

   VALIDATE_VECTOR (res);

   return res;
}


template<class T>
SXVEC<T> SXVEC<T>::choleskyDecomposition (enum UPLO uplo) const
{
   SX_CHECK      (handle);
   SX_CHECK      (getSize () > 0);
   SX_CHECK (nRows() == nCols(), nRows(), nCols());
   SX_CHECK      ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented
   SX_CHECK (uplo == UpperRight || uplo == LowerLeft);

   
   SXVEC<T> res; res.copy (*this);

   cholesky (res.elements, uplo, res.elements, (int)nRows());

   ssize_t n = nRows();
   if (uplo == LowerLeft)  {
      // --- clean up upper part of triangle of U
      for (ssize_t c=1; c < n; c++)
         for (ssize_t r=0; r < c; r++) 
            res(r,c) = 0.;
   } else {
      // --- clean up lower part of triangle of L
      for (ssize_t c=0; c < n; c++)
         for (ssize_t r=c+1; r < n; r++) 
            res(r,c) = 0.;
   }

   VALIDATE_VECTOR (res);
   return res;
}

template <class T>
void SXVEC<T>::rotate (const SXVEC<T> &rotMat)
{
   SX_CHECK (getSize () > 0);
   SX_CHECK (rotMat.nCols () == rotMat.nRows (),
             rotMat.nCols (), rotMat.nRows ());
   SX_CHECK (rotMat.getSize () == rotMat.nRows () * rotMat.nRows (),
             rotMat.getSize (), rotMat.nRows ());
   SX_CHECK (nCols () == rotMat.nCols (),
             nCols (), nRows ());
   inPlaceRot (elements, rotMat.elements, (int)nRows (), (int)nCols ());
   VALIDATE_VECTOR (*this);
}

template <class T>
SXVEC<T> SXVEC<T>::overlap (const SXVEC<T> &y) const
{
   SX_CHECK (getSize () > 0);
   SX_CHECK (y.getSize () > 0);
   SX_CHECK (nRows () == y.nRows (), nRows (), y.nRows ());
   if (nCols () == 1 && y.nCols () == 1)  {
      SXVEC<T> res(1);
      res(0) = dot(*this,y);
      VALIDATE_VECTOR (res);
      return res;
   }
   return overlap (y, y.nRows ());
}

template <class T>
SXVEC<T> SXVEC<T>::overlap (const SXVEC<T> &y, ssize_t sumSize) const
{
   SX_CHECK (getSize () > 0, getSize ());
   SX_CHECK (y.getSize () > 0, y.getSize ());
   SX_CHECK (nCols () > 0, nCols ());
   SX_CHECK (y.nCols () > 0, y.nCols ());
   SX_CHECK (sumSize > 0, sumSize);
   SX_CHECK (nRows () >= sumSize, nRows (), sumSize);
   SX_CHECK (y.nRows () >= sumSize, y.nRows (), sumSize);

   SXVEC<T> res(nCols () * y.nCols ());
   res.reshape (nCols (), y.nCols ());
   matovlp (res.elements, elements, y.elements,
            (int)nRows (), (int)nCols (), (int)y.nRows (), (int)y.nCols (),
            (int)sumSize);
   VALIDATE_VECTOR (res);
   return res;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

template<class T>
void SXVEC<T>::print (bool vectorForm) const
{
   print (cout, vectorForm);
   cout << endl;
}

template<class T>
void SXVEC<T>::print (std::ostream &s, bool vectorForm) const
{
   int prec = static_cast<int>(s.precision());
   ios::fmtflags flags = s.flags ();
   if (vectorForm)  {
      s.setf (ios::fixed, ios::floatfield);
      s.precision (3);
      s.fill(' ');
      s << setw(7);
   }

   if (handle->nCols > 1)  {
      // --- before printing, matrices have to be transposed due to
      //     their column major storage order
      SXVEC<T> tMat;
      if ( handle->parameters & IS_TRIANGULAR )  tMat = expand();
      else                                       tMat = *this;
      tMat.println (s, vectorForm);
   }  else  {
      println (s, vectorForm);
   }

   if (vectorForm)  { s << endl; }
   
   // --- restore predefined format
//   s.setf (ios_base::fmtflags(0), ios_base::floatfield); 
   s.precision (prec);
   s.flags (flags);
}


//template<class T>
//void SXVEC<T>::println (std::ostream &s, bool vectorForm) const
//{
//   SX_CHECK     (size);
//   SX_CHECK (*size >= 0, *size);
//   SX_CHECK      ( !(*parameters & IS_TRIANGULAR) );  // not yet implemented
//
//   int width = s.width();
//
//   if (*size == 0)  {
//      s << "(empty)" << endl;
//   }
//
//   // --- recursion loop
//   if (shape->getSize() == 1)  {  // last recursion level reached
//
//      for (int i=0; i < *size-1; i++)  
//         s << setw(width) << elements[i] << ", ";
//      s << setw(width) << elements[*size-1];
//
//   }  else  {                     // call deeper recursion level
//
//      Shape subShape (*shape);
//      subShape.remove (0);
//      Int::Type step = SXVEC<Int> (subShape).product();
//
//      SXVEC<T> v;
//      if ( !vectorForm )  s << "{";
//      for (int d=0; d < (*shape)(0); d++)  {
//         v = SXVEC<T> (elements + d*step, step);
//         v.reshape (subShape);
//         s << setw (width);
//         v.println (s, vectorForm);  // go one level deeper
//         if (d < (*shape)(0)-1)
//            if (vectorForm) {s << endl;} else {s << ", ";}
//      }
//      if (vectorForm) {s << endl;} else {s << "}";}
//   }
//}

template<class T>
void SXVEC<T>::println (std::ostream &s, bool vectorForm) const
{
   SX_CHECK     (handle);
   SX_CHECK (handle->size >= 0, handle->size);
   SX_CHECK     ( !(handle->parameters & IS_TRIANGULAR) );  // not yet implemented

   ssize_t i, r, c;
   int width = static_cast<int>(s.width());

   if ( !vectorForm )  s << "[";
   if (handle->size == 0)  {            // empty

      s << "empty";

   } else if (handle->nCols <= 1)  {     // vector print out

      for (i=0; i < handle->size - 1; i++)  
         s << setw(width) << elements[i] << ", ";
      s << setw(width) << elements[handle->size - 1];

   } else {

      for (r=0; r < handle->nRows; r++)  {  // matrix print out
         if ( !vectorForm )  s << "[";
         for (c=0; c < (handle->nCols)-1; c++)  {
            s << setw(width) << (*this)(r,c) << ", ";
         }
         s << setw(width) << (*this)(r,(handle->nCols)-1);
         if ( vectorForm) {s << endl;} else {s << "]";}
         if (!vectorForm && r < handle->nRows-1)  s << ",";
      }
   }
   if ( !vectorForm )  s << "]";
}


//===========================================================================

template<class T>
typename SXVEC<T>::Iterator SXVEC<T>::begin () const
{
   return Iterator (elements);
// return typename SXVEC<T>::Iterator (elements);
}


template<class T>
typename SXVEC<T>::Iterator SXVEC<T>::end () const
{
   return Iterator (elements + getSize());
// return typename SXVEC<T>::Iterator(elements + getSize ());
}

//===========================================================================
// The following functions are globally defined and will go to a
// SxOperators.h declaration file.
//===========================================================================


// --- scalar product   scp<?> = Vector<A> . Vector<B>
/// \todo: SMP parallelization
template<class A,class B>
inline typename SxTypeCaster<A,B>::Type
dot (const SXVEC<A> &a, const SXVEC<B> &b)
{
   SX_CHECK      (a.handle);
   SX_CHECK      (b.handle);
   SX_CHECK  (a.handle->size > 0, a.handle->size);
   SX_CHECK (a.handle->size == b.handle->size,
             a.handle->size, b.handle->size);
   
//   SXVEC<typename SxTypeCaster<A,B>::TRes> aVec(a);
//   SXVEC<typename SxTypeCaster<A,B>::TRes> aVec(a);
//   SXVEC<typename SxTypeCaster<A,B>::TRes> bVec(b);

   typename SxTypeCaster<A,B>::TRes::Type sum = 0.;
//   typename A::Type *aPtr = a.elements;
//   typename B::Type *bPtr = b.elements;
//   for (ssize_t i=0; i < *a.size; i++, aPtr++, bPtr++)
//      sum += ( conjugate(*aPtr) * *bPtr);

   sum = scalarProduct (a.elements, b.elements, static_cast<int>(a.getSize()));

   SX_CHECK_NUM (sum);
   return sum;
}

// --- Vector<?> = float + Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Float,B>::TRes> 
operator+ (const float &s, const SXVEC<B> &b)
{
   SX_CHECK     (b.handle);
   SX_CHECK (b.handle->size > 0, b.handle->size);

   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<Float,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = b.handle->auxData;

   ssize_t i;
   typename SxTypeCaster<Float,B>::TRes::Type *destPtr = res.elements;
   const typename B::Type *src1Ptr = b.elements;
#  ifdef USE_OPENMP   
      if (len > sxChunkSize)  {
         const float *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] 
               = (typename SxTypeCaster<Float,B>::TRes::Type)(*valPtr)
               + (typename SxTypeCaster<Float,B>::TRes::Type)src1Ptr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {         
         for (i=0; i < len; i++)
            *destPtr++ 
               = (typename SxTypeCaster<Float,B>::TRes::Type)s
               + (typename SxTypeCaster<Float,B>::TRes::Type)*src1Ptr++;
      }

   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = double + Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Double,B>::TRes> 
operator+ (const double &s, const SXVEC<B> &b)
{
   SX_CHECK     (b.handle);
   SX_CHECK (b.handle->size > 0, b.handle->size);

   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<Double,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = b.handle->auxData;

   ssize_t i;
   typename SxTypeCaster<Double,B>::TRes::Type *destPtr = res.elements;
   const typename B::Type *src1Ptr = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const double *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<Double,B>::TRes::Type)(*valPtr)
               + (typename SxTypeCaster<Double,B>::TRes::Type)src1Ptr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {
         for (i=0; i < len; i++)
            *destPtr++ 
               = (typename SxTypeCaster<Double,B>::TRes::Type)s
               + (typename SxTypeCaster<Double,B>::TRes::Type)*src1Ptr++;
      }

   VALIDATE_VECTOR (res);
   return res;
}

// --- Vector<?> = SxComplex8 + Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Complex8,B>::TRes> 
operator+ (const SxComplex8 &s, const SXVEC<B> &b)
{
   SX_CHECK     (b.handle);
   SX_CHECK (b.handle->size > 0, b.handle->size);

   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<Complex8,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = b.handle->auxData;

   ssize_t i;
   typename SxTypeCaster<Complex8,B>::TRes::Type *destPtr = res.elements;
   const typename B::Type *src1Ptr = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const SxComplex8 *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,src1Ptr,valPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<Double,B>::TRes::Type)(*valPtr)
               + (typename SxTypeCaster<Double,B>::TRes::Type)src1Ptr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {
         for (i=0; i < len; i++)
            *destPtr++ = s + *src1Ptr++;
      }

   VALIDATE_VECTOR (res);
   return res;
}

// --- Vector<?> = SxComplex16 + Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Complex16,B>::TRes> 
operator+ (const SxComplex16 &s, const SXVEC<B> &b)
{
   SX_CHECK     (b.handle);
   SX_CHECK (b.handle->size > 0, b.handle->size);

   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<Complex16,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = b.handle->auxData;

   ssize_t i;
   typename SxTypeCaster<Complex16,B>::TRes::Type *destPtr = res.elements;
   const typename B::Type *src1Ptr = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const SxComplex16 *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,src1Ptr,valPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<Double,B>::TRes::Type)(*valPtr)
               + (typename SxTypeCaster<Double,B>::TRes::Type)src1Ptr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {
         for (i=0; i < len; i++)
            *destPtr++ = s + *src1Ptr++;
      }

   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> + float
template<class A>
inline SXVEC<typename SxTypeCaster<A,Float>::TRes> 
operator+ (const SXVEC<A> &a, const float &s)
{
   return (s + a);
}


// --- Vector<?> = Vector<A> + double
template<class A>
inline SXVEC<typename SxTypeCaster<A,Double>::TRes> 
operator+ (const SXVEC<A> &a, const double &s)
{
   return (s + a);
}


// --- Vector<?> = Vector<A> + int
template<class A>
inline SXVEC<typename SxTypeCaster<A,Int>::TRes> 
operator+ (const SXVEC<A> &a, int s)
{
   SX_CHECK     (a.handle->size);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Int>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;
   
   typename SxTypeCaster<A,Int>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const int *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = src1Ptr[smpIdx] + (*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = *src1Ptr++ + s;
      }

   VALIDATE_VECTOR (res);
   return res;
}

// --- Vector<?> = Vector<A> + SxComplex8
template<class A>
inline SXVEC<typename SxTypeCaster<A,Complex8>::TRes> 
operator+ (const SXVEC<A> &a, const SxComplex8 &s)
{
   SX_CHECK     (a.handle->size);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Complex8>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;
   
   typename SxTypeCaster<A,Complex8>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const SxComplex8 *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,src1Ptr,valPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,Complex8>::TRes::Type)src1Ptr[smpIdx]  
               + (typename SxTypeCaster<A,Complex8>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++ 
               = (typename SxTypeCaster<A,Complex8>::TRes::Type)*src1Ptr++  
               + (typename SxTypeCaster<A,Complex8>::TRes::Type)s;
      }

   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> + SxComplex16
template<class A>
inline SXVEC<typename SxTypeCaster<A,Complex16>::TRes> 
operator+ (const SXVEC<A> &a, const SxComplex16 &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Complex16>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Complex16>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const SxComplex16 *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,src1Ptr,valPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,Complex16>::TRes::Type)src1Ptr[smpIdx]  
               + (typename SxTypeCaster<A,Complex16>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++      
               = (typename SxTypeCaster<A,Complex16>::TRes::Type)*src1Ptr++  
               + (typename SxTypeCaster<A,Complex16>::TRes::Type)s;
      }

   VALIDATE_VECTOR (res);
   return res;
}





// --- Vector<?> = Vector<A> + Vector<B>
template<class A,class B>
inline SXVEC<typename SxTypeCaster<A,B>::TRes>
operator+ (const SXVEC<A> &a, const SXVEC<B> &b) 
{
   SX_CHECK      (a.handle && b.handle);
   SX_CHECK (a.handle->size == b.handle->size,
             a.handle->size,   b.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = assign (a.handle->auxData, b.handle->auxData);

   typename SxTypeCaster<A,B>::TRes::Type *destPtr = res.elements;
   typename A::Type *src1Ptr = a.elements;
   typename B::Type *src2Ptr = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,src1Ptr,src2Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,B>::TRes::Type)src1Ptr[smpIdx]  
               + (typename SxTypeCaster<A,B>::TRes::Type)src2Ptr[smpIdx];
         return res;
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++     
               = (typename SxTypeCaster<A,B>::TRes::Type)*src1Ptr++  
               + (typename SxTypeCaster<A,B>::TRes::Type)*src2Ptr++;
      }

   VALIDATE_VECTOR (res);
   return res;
}



// --- Vector<?> = float - Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Float,B>::TRes> 
operator- (const float &s, const SXVEC<B> &b)
{
   SX_CHECK     (b.handle);
   SX_CHECK (b.handle->size > 0, b.handle->size);
   
   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<Float,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = b.handle->auxData;

   typename SxTypeCaster<Float,B>::TRes::Type *destPtr = res.elements;
   const typename B::Type *src1Ptr = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const float *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<Float,B>::TRes::Type)(*valPtr)
               - (typename SxTypeCaster<Float,B>::TRes::Type)src1Ptr[smpIdx];
      } else
#  endif /* USE_OPENMP */         
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ 
               = (typename SxTypeCaster<Float,B>::TRes::Type)s
               - (typename SxTypeCaster<Float,B>::TRes::Type)*src1Ptr++;
      }
      
   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = double - Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Double,B>::TRes> 
operator- (const double &s, const SXVEC<B> &b)
{
   SX_CHECK     (b.handle);
   SX_CHECK (b.handle->size > 0, b.handle->size);

   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<Double,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = b.handle->auxData;

   typename SxTypeCaster<Double,B>::TRes::Type *destPtr = res.elements;
   const typename B::Type *src1Ptr = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const double *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<Double,B>::TRes::Type)(*valPtr)
               - (typename SxTypeCaster<Double,B>::TRes::Type)src1Ptr[smpIdx];
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++
               = (typename SxTypeCaster<Double,B>::TRes::Type)s
               - (typename SxTypeCaster<Double,B>::TRes::Type)*src1Ptr++;
      }

   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> - float
template<class A>
inline SXVEC<typename SxTypeCaster<A,Float>::TRes> 
operator- (const SXVEC<A> &a, const float &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Float>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Float>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const float *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,src1Ptr,valPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
            = (typename SxTypeCaster<A,Float>::TRes::Type)src1Ptr[smpIdx]  
            - (typename SxTypeCaster<A,Float>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++
               = (typename SxTypeCaster<A,Float>::TRes::Type)*src1Ptr++  
               - (typename SxTypeCaster<A,Float>::TRes::Type)s;
      }

   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> - double
template<class A>
inline SXVEC<typename SxTypeCaster<A,Double>::TRes> 
operator- (const SXVEC<A> &a, const double &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Double>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Double>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const double *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
            = (typename SxTypeCaster<A,Double>::TRes::Type)src1Ptr[smpIdx]  
            - (typename SxTypeCaster<A,Double>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++
               = (typename SxTypeCaster<A,Double>::TRes::Type)*src1Ptr++  
               - (typename SxTypeCaster<A,Double>::TRes::Type)s;
      }
         
   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> - SxComplex8
template<class A>
inline SXVEC<typename SxTypeCaster<A,Complex8>::TRes> 
operator- (const SXVEC<A> &a, const SxComplex8 &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Complex8>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Complex8>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const SxComplex8 *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
            = (typename SxTypeCaster<A,Complex8>::TRes::Type)src1Ptr[smpIdx]  
            - (typename SxTypeCaster<A,Complex8>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++
               = (typename SxTypeCaster<A,Complex8>::TRes::Type)*src1Ptr++  
               - (typename SxTypeCaster<A,Complex8>::TRes::Type)s;
      }

   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> - SxComplex16
template<class A>
inline SXVEC<typename SxTypeCaster<A,Complex16>::TRes> 
operator- (const SXVEC<A> &a, const SxComplex16 &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Complex16>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Complex16>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const SxComplex16 *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,Complex16>::TRes::Type)src1Ptr[smpIdx]  
               - (typename SxTypeCaster<A,Complex16>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++     
               = (typename SxTypeCaster<A,Complex16>::TRes::Type)*src1Ptr++  
               - (typename SxTypeCaster<A,Complex16>::TRes::Type)s;
      }
      
   VALIDATE_VECTOR (res);
   return res;
}





// --- Vector<?> = Vector<A> - Vector<B>
template<class A,class B>
inline SXVEC<typename SxTypeCaster<A,B>::TRes>
operator- (const SXVEC<A> &a, const SXVEC<B> &b) 
{
   SX_CHECK      (a.handle && b.handle);
   SX_CHECK (a.handle->size == b.handle->size,
             a.handle->size,   b.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = assign (a.handle->auxData, b.handle->auxData);
   
   typename SxTypeCaster<A,B>::TRes::Type *destPtr = res.elements;
   typename A::Type *src1Ptr = a.elements;
   typename B::Type *src2Ptr = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx = 0;
#        pragma omp parallel for shared(destPtr,src1Ptr,src2Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,B>::TRes::Type)src1Ptr[smpIdx]  
               - (typename SxTypeCaster<A,B>::TRes::Type)src2Ptr[smpIdx];
      } else
#  endif /* USE_OPENMP */
   {
      for (ssize_t i=0; i < len; i++)
         *destPtr++
            = (typename SxTypeCaster<A,B>::TRes::Type)*src1Ptr++  
            - (typename SxTypeCaster<A,B>::TRes::Type)*src2Ptr++;
   }
      
   VALIDATE_VECTOR (res);
   return res;
}



// --- Vector<?> = float * Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Float,B>::TRes> 
operator* (const Float::Type &s, const SXVEC<B> &vec)
{
   return SXVEC<typename SxTypeCaster<Float,B>::TRes> (vec) 
        * (typename SxTypeCaster<Float,B>::TRes::Type)s;
}


// --- Vector<?> = double * Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Double,B>::TRes> 
operator* (const Double::Type &s, const SXVEC<B> &vec)
{
   return SXVEC<typename SxTypeCaster<Double,B>::TRes> (vec) 
        * (typename SxTypeCaster<Double,B>::TRes::Type)s;
}

// --- Vector<?> = SxComplex8 * Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Complex8,B>::TRes> 
operator* (const Complex8::Type &s, const SXVEC<B> &vec)
{
   return SXVEC<typename SxTypeCaster<Complex8,B>::TRes> (vec) 
        * (typename SxTypeCaster<Complex8,B>::TRes::Type)s;
}

// --- Vector<?> = SxComplex16 * Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Complex16,B>::TRes> 
operator* (const Complex16::Type &s, const SXVEC<B> &vec)
{
   return SXVEC<typename SxTypeCaster<Complex16,B>::TRes> (vec) 
        * (typename SxTypeCaster<Complex16,B>::TRes::Type)s;
}



// --- Vector<?> = Vector<A> * float
template<class A>
inline SXVEC<typename SxTypeCaster<A,Float>::TRes> 
operator* (const SXVEC<A> &a, const float &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Float>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Float>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const float *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,Float>::TRes::Type)src1Ptr[smpIdx]  
               * (typename SxTypeCaster<A,Float>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++
               = (typename SxTypeCaster<A,Float>::TRes::Type)*src1Ptr++  
               * (typename SxTypeCaster<A,Float>::TRes::Type)s;
      }
      
   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> * double
template<class A>
inline SXVEC<typename SxTypeCaster<A,Double>::TRes> 
operator* (const SXVEC<A> &a, const double &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);
   SX_CHECK_NUM (s);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Double>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Double>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const double *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,Double>::TRes::Type)src1Ptr[smpIdx]  
               * (typename SxTypeCaster<A,Double>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++      
               = *src1Ptr++  * s;
//            *destPtr++      
//               = (typename SxTypeCaster<A,Double>::TRes::Type)*src1Ptr++  
//               * (typename SxTypeCaster<A,Double>::TRes::Type)s;
      }
      
   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> * SxComplex8
template<class A>
inline SXVEC<typename SxTypeCaster<A,Complex8>::TRes> 
operator* (const SXVEC<A> &a, const SxComplex8 &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Complex8>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Complex8>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const SxComplex8 *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,Complex8>::TRes::Type)src1Ptr[smpIdx]  
               * (typename SxTypeCaster<A,Complex8>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)
            *destPtr++     
               = (typename SxTypeCaster<A,Complex8>::TRes::Type)*src1Ptr++  
               * (typename SxTypeCaster<A,Complex8>::TRes::Type)s;
      }
      
   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<?> = Vector<A> * SxComplex16
template<class A>
inline SXVEC<typename SxTypeCaster<A,Complex16>::TRes> 
operator* (const SXVEC<A> &a, const SxComplex16 &s)
{
   SX_CHECK     (a.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,Complex16>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = a.handle->auxData;

   typename SxTypeCaster<A,Complex16>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const SxComplex16 *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,src1Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,Complex16>::TRes::Type)src1Ptr[smpIdx]  
               * (typename SxTypeCaster<A,Complex16>::TRes::Type)(*valPtr);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++   
               = (typename SxTypeCaster<A,Complex16>::TRes::Type)*src1Ptr++  
               * (typename SxTypeCaster<A,Complex16>::TRes::Type)s;
      }
      
   VALIDATE_VECTOR (res);
   return res;
}

// --- Matrix<?> = Matrix<A> * Vector<B> 
template<class A,class B>
inline SXVEC<typename SxTypeCaster<A,B>::TRes>
starR (const SXVEC<A> &a, const SXVEC<B> &b) 
{
   SX_CHECK      (a.handle && b.handle);
   SX_CHECK (a.handle->nCols == b.handle->size,
             a.handle->nCols,   b.handle->size);

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,B>::TRes> res (len);
   
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = assign (a.handle->auxData, b.handle->auxData);
 
   typename SXVEC<A>::Iterator aIt = a.begin ();
   typename SXVEC<B>::Iterator bIt = b.begin ();
   typename SXVEC<typename SxTypeCaster<A,B>::TRes>::Iterator resIt 
      = res.begin ();
   size_t nr = a.nRows (), 
          nc = a.nCols ();

   for (size_t ic = nc; ic; --ic)  {
      for (size_t ir = nr; ir; --ir)//fast rows!
      {  *resIt++ = *aIt++ * *bIt; }
      bIt++;
   }   
   VALIDATE_VECTOR (res);
   return res;
}

// --- Matrix<?> =  Vector<A> * Matrix<B> 
template<class A,class B>
inline SXVEC<typename SxTypeCaster<A,B>::TRes>
starL (const SXVEC<A> &a, const SXVEC<B> &b) 
{
   SX_CHECK      (a.handle && b.handle);
   SX_CHECK (a.handle->size == b.handle->nRows,
             a.handle->size,   b.handle->nRows);

   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<A,B>::TRes> res (len);
   
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = assign (a.handle->auxData, b.handle->auxData);
   
   typename SXVEC<A>::Iterator aIt;
   typename SXVEC<B>::Iterator bIt = b.begin ();
   typename SXVEC<typename SxTypeCaster<A,B>::TRes>::Iterator resIt 
      = res.begin ();
   size_t nr = b.nRows (), 
          nc = b.nCols ();

   for (size_t ic = nc; ic; --ic)  {
      aIt = a.begin ();
      for (size_t ir = nr; ir; --ir)//fast rows!
      {  *resIt++ = *aIt++ * *bIt++; }
   }   
   VALIDATE_VECTOR (res);
   return res;
}

// --- Vector<?> = Vector<A> * Vector<B>
template<class A,class B>
inline SXVEC<typename SxTypeCaster<A,B>::TRes>
operator* (const SXVEC<A> &a, const SXVEC<B> &b) 
{
   SX_CHECK      (a.handle && b.handle);
   SX_CHECK (a.handle->size > 0, a.handle->size);
   SX_CHECK (b.handle->size > 0, b.handle->size);
   if (a.handle->size > b.handle->size
       && a.handle->nCols == b.handle->size) return starR (a, b);
   if (a.handle->size < b.handle->size
       && a.handle->size == b.handle->nCols) return starL (a, b);

   SX_CHECK (a.handle->size == b.handle->size,
             a.handle->size, b.handle->size);
   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,B>::TRes> res (len);
   res.handle->nRows   = a.handle->nRows;
   res.handle->nCols   = a.handle->nCols;
   res.handle->auxData = assign (a.handle->auxData, b.handle->auxData);
   
   typename       SxTypeCaster<A,B>::TRes::Type  *destPtr = res.elements;
   const typename A::Type  *src1Ptr = a.elements;
   const typename B::Type  *src2Ptr = b.elements;
#  ifdef USE_OPENMP
#        pragma omp parallel for if (len > sxChunkSize)
         for (ssize_t smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = src1Ptr[smpIdx] * src2Ptr[smpIdx];
#else
         for (ssize_t i=0; i < len; i++)
            *destPtr++   
               = (typename SxTypeCaster<A,B>::TRes::Type)*src1Ptr++  
               * (typename SxTypeCaster<A,B>::TRes::Type)*src2Ptr++;
#  endif /* USE_OPENMP */         
      
   VALIDATE_VECTOR (res);
   return res;
}

// --- float / Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Float,B>::TRes> 
operator/ (float s, const SXVEC<B> &b)
{
   SX_CHECK     (b.handle->size);
   SX_CHECK (b.handle->size > 0, b.handle->size);

   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<Float,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = b.handle->auxData;
   
   typename SxTypeCaster<Float,B>::TRes::Type *destPtr = res.elements;
   typename                           B::Type *srcPtr  = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const float *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)  {
            SX_CHECK_DIV (srcPtr[smpIdx]);
            destPtr[smpIdx] 
               = (typename SxTypeCaster<Float,B>::TRes::Type)(*valPtr) 
               / (typename SxTypeCaster<Float,B>::TRes::Type)srcPtr[smpIdx];
         }
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)  {
            SX_CHECK_DIV (*srcPtr);
            *destPtr++ 
               = (typename SxTypeCaster<Float,B>::TRes::Type)s 
               / (typename SxTypeCaster<Float,B>::TRes::Type)*srcPtr++;
         }
      }
      
   VALIDATE_VECTOR (res);
   return res;
}

// --- double / Vector<B>
template<class B>
inline SXVEC<typename SxTypeCaster<Double,B>::TRes> 
operator/ (double s, const SXVEC<B> &b)
{
   SX_CHECK     (b.handle);
   SX_CHECK (b.handle->size > 0, b.handle->size);

   const ssize_t &len = b.handle->size;
   SXVEC<typename SxTypeCaster<Double,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = b.handle->auxData;
   
   typename SxTypeCaster<Double,B>::TRes::Type *destPtr = res.elements;
   typename                            B::Type *srcPtr  = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         const double *valPtr = &s;
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,valPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)  {
            SX_CHECK_DIV (srcPtr[smpIdx]);
            destPtr[smpIdx]
               = (typename SxTypeCaster<Double,B>::TRes::Type)(*valPtr) 
               / (typename SxTypeCaster<Double,B>::TRes::Type)srcPtr[smpIdx];
         }
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)  {
            SX_CHECK_DIV (*srcPtr);
            *destPtr++
               = (typename SxTypeCaster<Double,B>::TRes::Type)s 
               / (typename SxTypeCaster<Double,B>::TRes::Type)*srcPtr++;
         }
      }
      
   VALIDATE_VECTOR (res);
   return res;
}


// --- Vector<A> / float
template<class A>
inline SXVEC<typename SxTypeCaster<A,Float>::TRes> 
operator/ (const SXVEC<A> &a, const float &s)
{
   SX_CHECK (::fabs((Float::Type)s) > 1e-10, (Float::Type)s);
   SX_CHECK_DIV (s);
   // res = a * (1./s);
   return SXVEC<typename SxTypeCaster<A,Float>::TRes> (a) 
        * (typename SxTypeCaster<A,Float>::TRes::Type)(1. / s);
}

// --- Vector<A> / double
template<class A>
inline SXVEC<typename SxTypeCaster<A,Double>::TRes> 
operator/ (const SXVEC<A> &a, const double &s)
{
   SX_CHECK_DIV (s);
   // res = a * (1./s);
   return SXVEC<typename SxTypeCaster<A,Double>::TRes> (a) 
        * (typename SxTypeCaster<A,Double>::TRes::Type)(1. / s);
}


// --- Vector<A> / SxComplex8
template<class A>
inline SXVEC<typename SxTypeCaster<A,Complex8>::TRes> 
operator/ (const SXVEC<A> &a, const Complex8::Type &s)
{
   SX_CHECK_DIV (s);
   // res = a * (1./s);
   return SXVEC<typename SxTypeCaster<A,Complex8>::TRes> (a) 
        * (typename SxTypeCaster<A,Complex8>::TRes::Type)(SxComplex8(1.) / s);
}


// --- Vector<A> / SxComplex16
template<class A>
inline SXVEC<typename SxTypeCaster<A,Complex16>::TRes> 
operator/ (const SXVEC<A> &a, const SxComplex16 &s)
{
   SX_CHECK_DIV (s);
   // res = a * (1./s);
   return SXVEC<typename SxTypeCaster<A,Complex16>::TRes> (a) 
        * (typename SxTypeCaster<A,Complex16>::TRes::Type)(SxComplex16(1.) / s);
}



// --- Vector<A> / Vector<B>
template<class A,class B>
inline SXVEC<typename SxTypeCaster<A,B>::TRes> 
operator/ (const SXVEC<A> &a, const SXVEC<B> &b)
{
   SX_CHECK      (a.handle && b.handle);
   SX_CHECK (a.handle->size == b.handle->size,
             a.handle->size,   b.handle->size);
   SX_CHECK      (   !(a.handle->parameters & IS_TRIANGULAR)    // not yet
                  && !(b.handle->parameters & IS_TRIANGULAR));  // implemented

   const ssize_t &len = a.handle->size;
   SXVEC<typename SxTypeCaster<A,B>::TRes> res (len);
   res.handle->nRows   = b.handle->nRows;
   res.handle->nCols   = b.handle->nCols;
   res.handle->auxData = assign (a.handle->auxData, b.handle->auxData);
   
   typename SxTypeCaster<A,B>::TRes::Type *destPtr = res.elements;
   const typename A::Type *src1Ptr = a.elements;
   const typename B::Type *src2Ptr = b.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,src1Ptr,src2Ptr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)  {
            SX_CHECK_DIV (src2Ptr[smpIdx]);
            destPtr[smpIdx]
               = (typename SxTypeCaster<A,B>::TRes::Type)src1Ptr[smpIdx]  
               / (typename SxTypeCaster<A,B>::TRes::Type)src2Ptr[smpIdx];
         }
      } else
#  endif /* USE_OPENMP */
      {         
         for (ssize_t i=0; i < len; i++)  {
            SX_CHECK_DIV (*src2Ptr);
            *destPtr++ 
               = (typename SxTypeCaster<A,B>::TRes::Type)*src1Ptr++  
               / (typename SxTypeCaster<A,B>::TRes::Type)*src2Ptr++;
         }
      }
      
   VALIDATE_VECTOR (res);
   return res;
}


// ---         ? = Vector<A> . Vector<B>          scalar product
//     Matrix<?> = Vector<A> . Matrix<B>          row vector . matrix multipl.
//     Vector<?> = Matrix<A> . Vector<B>          matrix . column vector
//     Matrix<?> = Matrix<A> . Matrix<B>          matrix multiplication
/** \todo support special multiplication for trigonal matrices */
template<class A,class B>
inline SXVEC<typename SxTypeCaster<A,B>::TRes> 
operator^ (const SXVEC<A> &a, const SXVEC<B> &b)
{  
   // --- Check ranks
   SX_CHECK      (a.handle->size > 0 && a.elements);
   SX_CHECK      (b.handle->size > 0 && b.elements);
   SX_CHECK      (a.handle && b.handle);
   SX_CHECK      (   !(a.handle->parameters & IS_TRIANGULAR)    // not yet
                  && !(b.handle->parameters & IS_TRIANGULAR));  // implemented

   SXVEC<typename SxTypeCaster<A,B>::TRes> res;

// if         (a.shape->getSize()==1 && b.shape->getSize()==1)  { 
   if         (a.handle->nCols < 2 && b.handle->nCols < 2)  { 
      //catch case that b is one-element vector
      if (b.handle->nRows < 2)  {
         res = a * b(0);
      }  else  {

         SX_CHECK (a.getSize() == b.getSize(), a.getSize(), b.getSize());
      
         // vector ^ vector
         res = SXVEC<typename SxTypeCaster<A,B>::TRes> (1);
         // C++ magic: perform type cast if TRres != A, noop otherwise
         const SXVEC<typename SxTypeCaster<A,B>::TRes> &aVec = a;
         // C++ magic: perform type cast if TRres != B, noop otherwise
         const SXVEC<typename SxTypeCaster<A,B>::TRes> &bVec = b;

         *res.elements = scalarProduct (aVec.elements, 
                                        bVec.elements, 
                                        static_cast<int>(a.getSize()));
      }
   }  else if (a.handle->nCols <= 1  && b.handle->nRows == 1)    {
      
      // dyadic product res(x,y) = a(x) * b(y)
      SX_CHECK (a.handle->nRows > 0, a.handle->nRows);
      SX_CHECK (b.handle->nCols > 0, b.handle->nCols);
      
      res.resize  (a.handle->nRows * b.handle->nCols);
      res.reshape (a.handle->nRows, b.handle->nCols);
      ssize_t x,y;
      typename SxTypeCaster<A,B>::TRes::Type *destPtr = res.elements;
      typename A::Type *aPtr;
      typename B::Type *bPtr = b.elements;
      
      for (y = 0; y < b.handle->nCols; ++y, ++bPtr)  {
         aPtr = a.elements;
         for (x = 0; x < a.handle->nRows; ++x)  {
            *destPtr++ = *aPtr++ * *bPtr;
         }
      }
      
// }  else if (a.shape->getSize()==2 || b.shape->getSize()==2)  {
   }  else if (a.handle->nCols > 1   || b.handle->nCols > 1)    {

      // C++ magic: perform type cast if TRres != A, noop otherwise
      const SXVEC<typename SxTypeCaster<A,B>::TRes> &aMat = a;
      // C++ magic: perform type cast if TRres != B, noop otherwise
      const SXVEC<typename SxTypeCaster<A,B>::TRes> &bMat = b;

//    if         (a.shape->getSize()==1)     {               
      if         (a.handle->nCols < 2)     {               

         // row-vector ^ matrix
         SX_CHECK (a.getSize() == b.nRows(), a.getSize(), b.nRows());
         res = SXVEC<typename SxTypeCaster<A,B>::TRes> (1 * bMat.nCols());

         matmult (res.elements, aMat.elements, bMat.elements,
                  1, (int)aMat.getSize(), (int)bMat.nCols());
         res.reshape (1, bMat.nCols());

//    }  else if (b.shape->getSize()==1)  {
      }  else if (b.handle->nCols < 2)  {

         //     matrix ^ colVector
         SX_CHECK (a.nCols() == b.getSize(), a.nCols(), b.getSize());
         res = SXVEC<typename SxTypeCaster<A,B>::TRes> (aMat.nRows() * 1);
         matmult (res.elements, aMat.elements, bMat.elements,
                  (int)aMat.nRows(), (int)aMat.nCols(), 1);
//         ssize_t i, r, nRows = a.nRows(), nCols = a.nCols();
//         typename SxTypeCaster<A,B>::TRes::Type s, *aPtr, *bPtr;
//         for (r=0; r < nRows; ++r)  {
//            aPtr = aMat.elements + r;
//            bPtr = bMat.elements;
//            for (s=0., i=0; i < nCols; ++i, aPtr+=nCols, bPtr++)  {
//               s += *aPtr * *bPtr;
//            }
//            res(r) = s;
//         } 
         res.reshape (aMat.nRows(), 1);
         
      }  else  { 

         //     matrix ^ matrix
         SX_CHECK (a.nCols() == b.nRows(), a.nCols(), b.nRows());
         res = SXVEC<typename SxTypeCaster<A,B>::TRes> (  aMat.nRows() 
                                                        * bMat.nCols());
         matmult (res.elements, aMat.elements, bMat.elements,
                  (int)aMat.nRows(), (int)aMat.nCols(), (int)bMat.nCols());
         res.reshape (aMat.nRows(), bMat.nCols());

      }
//      //      rowVector -> SXVEC
//      if    (res.nCols() == 1 && res.nRows() > 1)  res.reshape (res.nRows());
//      else // colVector -> SXVEC
//         if (res.nRows() == 1 && res.nCols() > 1)  res.reshape (res.nCols());

   }  else  {
      printf ("operator^ not defined for ranks higher than 2!!!\n");
      SX_EXIT;
   }

   // no metadata for new matrix
   //res.handle->auxData = assign (a.handle->auxData, b.handle->auxData);
   VALIDATE_VECTOR (res);
   return res;
}




template<class T>
inline SXVEC<T> sqrt (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = sqrt (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = sqrt (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> pow (const SXVEC<T> &in, double exponent)
{
   SX_CHECK (in.handle);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = pow (srcPtr[smpIdx], exponent);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = pow (*srcPtr++, exponent);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> sin (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = sin (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = sin (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> cos (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = cos (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = cos (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> atan (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = atan (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = atan (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> sinh (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = sinh (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = sinh (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> cosh (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = cosh (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = cosh (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> exp (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = exp (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = exp (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}

template<class T>
inline SXVEC<T> log (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = log (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = log (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> derf (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = derf (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = derf (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
inline SXVEC<T> derfc (const SXVEC<T> &in)
{
   SX_CHECK (in.handle);
   SX_CHECK (in.handle->size > 0);

   const ssize_t &len = in.handle->size;
   SXVEC<T> res (len);
   res.handle->nRows   = in.handle->nRows;
   res.handle->nCols   = in.handle->nCols;
   res.handle->auxData = in.handle->auxData;

   typename T::Type *destPtr = res.elements;
   typename T::Type *srcPtr  = in.elements;
#  ifdef USE_OPENMP
      if (len > sxChunkSize)  {
         ssize_t smpIdx;
#        pragma omp parallel for shared(destPtr,srcPtr) \
         private(smpIdx) schedule(static)
         for (smpIdx=0; smpIdx < len; smpIdx++)
            destPtr[smpIdx] = derfc (srcPtr[smpIdx]);
      } else
#  endif /* USE_OPENMP */
      {
         for (ssize_t i=0; i < len; i++)
            *destPtr++ = derfc (*srcPtr++);
      }

   VALIDATE_VECTOR (res);
   return res;
}


template<class T>
void swap (SXVEC<T> &a, SXVEC<T> &b)
{
   SX_CHECK (a.getSize() == b.getSize(), a.getSize(), b.getSize());
   SX_CHECK      (   !(a.handle->parameters & IS_TRIANGULAR)    // not yet
                  && !(b.handle->parameters & IS_TRIANGULAR));  // implemented
   SXVEC<T> c (a); a = b; b = c;
}

template <class T>
SXVEC<T> mergeCols (const SXVEC<T> &x, const SXVEC<T> &y)
{
   SX_CHECK (x.getSize () > 0, x.getSize ());
   SX_CHECK (y.getSize () > 0, y.getSize ());
   SX_CHECK (x.nRows () == y.nRows (), x.nRows (), y.nRows ());
   int n = (int)x.nRows (), nx = (int)x.nCols (), ny = (int)y.nCols ();
   int nxy = nx + ny;
   SXVEC<T> res;
   res.reformat (n, nxy);
   res(SxIdx(0, n * nx - 1))       <<= x;
   res(SxIdx(n * nx, n * nxy - 1)) <<= y;
   res.handle->auxData = x.handle->auxData;
   return res;
}


#endif  /* SKIP_DOXYGEN */


template<class T>
std::ostream& operator<< (std::ostream &s, const SXVEC<T> &in)
{
   in.print (s, false);
   return s;
}



// --------------------------------------------------------------------------

template<class T>
size_t getNBytes (const SXVEC<T> &in)
{
   size_t nBytes = 0;
   if (in.handle->parameters & MANAGE_MEMORY)  {
      nBytes += in.getSize() * sizeof (typename T::Type);
   }
   nBytes += sizeof (MEMHANDLE *)          // *handle
           + sizeof (MEMHANDLE)            //  handle
           + sizeof (typename T::Type *);  // *elements
   return nBytes;
}
