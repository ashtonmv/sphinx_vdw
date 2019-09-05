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

#ifndef _SX_VECTOR3_H_
#define _SX_VECTOR3_H_

#include <stdio.h>
#include <math.h>
#include <SxComplex.h>
#include <SxList.h>
#include <SxMathLib.h>
#include <SxArray.h>
#include <iostream>
#include <SxTypeMapper.h>
#include <SxPrecision.h>
#include <SxConstants.h>
#include <SxString.h>


/**
  @ingroup Numerics
  \note Different to the vector classes, SxVector3 does not support reference
  counting. Instead, asVec3Ref may be used to make any 3 numbers behave as
  a SxVector3. For this, it is absolutely crucial that this SxVector3 remains
  memory-compatible to (T::Type)[3].
  Under no circumstances, virtual functions or additional data may be added.
  */
template<class T>
class SxVector3
{
   public:
      /** \brief The vector elements */
      typename T::Type v[3];

      inline SxVector3 ();
      inline explicit SxVector3 (const typename T::Type &xyz)  
      { 
         v[0] = xyz, v[1] = xyz, v[2] = xyz; 
      }
      inline SxVector3 (const typename T::Type &xVal, 
                        const typename T::Type &yVal, 
                        const typename T::Type &zVal)  
      { 
         v[0] = xVal, v[1] = yVal, v[2] = zVal; 
      }
      inline SxVector3 (const SxVector3<T> &in);
      template<class B>
         inline SxVector3 (const SxVector3<B> &in);
      inline explicit SxVector3 (const SxList<typename T::Type> &);
      inline explicit SxVector3 (const SxArray<typename T::Type> &);
      
      inline ~SxVector3 ();

      /** \brief Reinterpret 3 numbers in memory as a SxVector3 
        \example
        \code
SxMatrix<Double> m(3,nAtoms);
SxVector3<Double> t;
for (int ia = 0; ia < nAtoms; ia++)  {
   // m.colRef(ia) += t; // would not work
   SxVector3<Double>::toVec3Ref(m.colRef(ia).elements) += t; // works
}
       \endcode
       @return a SxVector3 reference to the memory pointed to by in.
       */
      static SxVector3<T> & toVec3Ref(typename T::Type in[3])
      {
         return reinterpret_cast<SxVector3<T> &> (*in); 
      }
      static const SxVector3<T> & toVec3Ref(const typename T::Type in[3])
      { 
         return reinterpret_cast<const SxVector3<T> &> (*in); 
      }

      inline       typename T::Type &operator() (ssize_t i);
      inline const typename T::Type &operator() (ssize_t i) const;
      inline       typename T::Type &operator() (SxAutoLoop &i);
      inline const typename T::Type &operator() (SxAutoLoop &i) const;
      /** Assign scalar value to vector */
      inline SxVector3<T>     operator= (const typename T::Type &s);
      /** Assign vector to vector */
      inline SxVector3<T> &   operator= (const SxVector3<T> &in);
      template<class B>
      inline SxVector3<T> &   operator= (const SxVector3<B> &in);
      inline void             operator+= (const SxVector3<T> &r);
      inline void             operator+= (const typename T::Type &s);
      inline void             operator-= (const typename T::Type &s);
      inline void             operator/= (const SxVector3<T> &r);
      inline void             operator*= (const typename T::Type &s);
      inline void             operator/= (const typename T::Type &s);
      /** the unary (-) operator */
      inline SxVector3<T>     operator-  () const;
      inline void             operator-= (const SxVector3<T> &r);
      inline bool             operator== (const SxVector3<T> &in) const;
      inline bool operator!= (const SxVector3<T> &in) const {
         return ! operator== (in);
      }

      inline SxVector3<T>  x (const SxVector3<T> &in) const;
      inline SxVector3<typename T::TReal> absSqr () const;
      inline void set (const typename T::Type &);
      inline typename T::Type sum () const;
      inline typename T::Type product () const;
      /** \brief Normalize the vector
        @return Reference to the normalized vector
        \note The return value allows for on-the-fly normalization
              of expressions:
        \code
n12 = (a1.x(a2)).normalize ();
        \endcode
        \note For getting a normalized vector (without changing the original)
              use
        \code
normalizedCopy = SxVector3<?> (original).normalize ();
        \endcode
        */
      inline SxVector3<T> &normalize ();
      /// Square of length of vector
      inline typename T::Real normSqr () const;
      /// Length of vector
      inline typename T::Real norm () const;
      /** minimum element */
      inline typename T::Type minval () const;
      /** maximum element */
      inline typename T::Type maxval () const;

      /// SxVector3<T> + int
      SxVector3<typename SxTypeCaster<Int, T>::TRes>
      operator+ (int);

      void print () const;

};


template<class T>
SxVector3<T>::SxVector3 ()
{
   v[0] = (typename T::Type)0;
   v[1] = (typename T::Type)0;
   v[2] = (typename T::Type)0;
}



template<class T>
   SxVector3<T>::SxVector3 (const SxVector3<T> &in)
{
   v[0] = in.v[0]; 
   v[1] = in.v[1]; 
   v[2] = in.v[2];
}

template<class T>
   template<class B>
      SxVector3<T>::SxVector3 (const SxVector3<B> &in)
{
   v[0] = (typename T::Type)in.v[0]; 
   v[1] = (typename T::Type)in.v[1]; 
   v[2] = (typename T::Type)in.v[2];
}


template<class T>
SxVector3<T>::SxVector3 (const SxList<typename T::Type> &list)
{
   SX_CHECK (list.getSize() == 3, list.getSize());
   v[0] = list(0); v[1] = list(1); v[2] = list(2);
}

template<class T>
SxVector3<T>::SxVector3 (const SxArray<typename T::Type> &array)
{
   SX_CHECK (array.getSize() == 3, array.getSize());
   v[0] = array(0); v[1] = array(1); v[2] = array(2);
}

template<class T>
SxVector3<T>::~SxVector3 ()
{
   // empty
}

template<class T>
typename T::Type &SxVector3<T>::operator() (ssize_t i)
{
   SX_CHECK (i >= 0 && i < 3, i); 
   return v[i];
}


template<class T>
const typename T::Type &SxVector3<T>::operator() (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < 3, i); 
   return v[i];
}

template<class T>
typename T::Type &SxVector3<T>::operator() (SxAutoLoop &i)
{
   SX_CHECK (i >= 0 && i < 3, i); 
   i.setLimit (3);
   return v[i];
}


template<class T>
const typename T::Type &SxVector3<T>::operator() (SxAutoLoop &i) const
{
   SX_CHECK (i >= 0 && i < 3, i); 
   i.setLimit (3);
   return v[i];
}


template<class T>
SxVector3<T> SxVector3<T>::operator= (const typename T::Type &s)
{
   v[0] = v[1] = v[2] = s;
   return *this;
}


template<class T>
SxVector3<T> &SxVector3<T>::operator= (const SxVector3<T> &in)
{
   v[0] = in.v[0];
   v[1] = in.v[1];
   v[2] = in.v[2];
   return *this;
}

template<class T>
template<class B>
SxVector3<T> &SxVector3<T>::operator= (const SxVector3<B> &in)
{
   v[0] = static_cast<typename T::Type>(in.v[0]);
   v[1] = static_cast<typename T::Type>(in.v[1]);
   v[2] = static_cast<typename T::Type>(in.v[2]);
   return *this;
}



template<class T>
void SxVector3<T>::operator+= (const SxVector3<T> &r)
{
   v[0] += r.v[0]; v[1] += r.v[1]; v[2] += r.v[2];
}

template<class T>
void SxVector3<T>::operator+= (const typename T::Type &s)
{
   v[0] += s; v[1] += s; v[2] += s;
}

template<class T>
void SxVector3<T>::operator/= (const SxVector3<T> &r)
{
   SX_CHECK_DIV (r.v[0]);
   SX_CHECK_DIV (r.v[1]);
   SX_CHECK_DIV (r.v[2]);
   v[0] /= r.v[0]; v[1] /= r.v[1]; v[2] /= r.v[2];
}

template<class T>
void SxVector3<T>::operator*= (const typename T::Type &s)
{
   v[0] *= s; v[1] *= s; v[2] *= s;
}


template<class T>
void SxVector3<T>::operator/= (const typename T::Type &s)
{
   SX_CHECK_DIV (s);
   v[0] /= s; v[1] /= s; v[2] /= s;
}


template<class T>
SxVector3<T> SxVector3<T>::operator- () const
{
   return SxVector3<T> (-v[0], -v[1], -v[2]);
}


template<class T>
void SxVector3<T>::operator-= (const SxVector3<T> &r)
{
   v[0] -= r.v[0]; v[1] -= r.v[1]; v[2] -= r.v[2];
}

template<class T>
void SxVector3<T>::operator-= (const typename T::Type &s)
{
   v[0] -= s; v[1] -= s; v[2] -= s;
}

template<class T>
bool SxVector3<T>::operator== (const SxVector3<T> &in) const
{
   for (int i=0; i<3; i++)
      if (v[i] != in.v[i] )   return false;
   return true;
}

template<>
inline bool SxVector3<Float>::operator== (const SxVector3<Float> &in) const
{
   for (int i=0; i<3; i++)
      if ( fabs(v[i] - in.v[i]) > 1e-10 )   return false;
   return true;
}

template<>
inline bool SxVector3<Double>::operator== (const SxVector3<Double> &in) const
{
   for (int i=0; i<3; i++)
      if ( fabs(v[i] - in.v[i]) > 1e-10 )   return false;
   return true;
}


template<class T>
SxVector3<T> SxVector3<T>::x (const SxVector3<T> &in) const
{
   return SxVector3<T> ( v[1]*in.v[2] - v[2]*in.v[1],
		                   v[2]*in.v[0] - v[0]*in.v[2],
		                   v[0]*in.v[1] - v[1]*in.v[0] );
}

// --- absSqr for different types
inline float absSqr (float x)       { return x*x; }
inline double absSqr (double x)     { return x*x; }
inline int absSqr (int i)           { return i*i; }
inline long int absSqr (long int i) { return i*i; }

template<class T>
inline typename T::TReal absSqr (const SxComplex<T> &z)
{ 
   return z.absSqr ();
}
// ------------------------------

template<class T>
SxVector3<typename T::TReal> SxVector3<T>::absSqr () const
{
   return SxVector3<typename T::TReal>(::absSqr (v[0]), 
                                       ::absSqr (v[1]), 
                                       ::absSqr (v[2]));
}

template<class T>
void SxVector3<T>::set (const typename T::Type &in)
{
   v[0] = in; v[1] = in; v[2] = in;
}



template<class T>
typename T::Type SxVector3<T>::sum () const
{
   return v[0] + v[1] + v[2];
}


template<class T>
typename T::Type SxVector3<T>::product () const
{
   return v[0] * v[1] * v[2];
}

template<class T>
SxVector3<T> &SxVector3<T>::normalize ()
{
   typename T::Real nrm, c;
   
   nrm = norm ();

   SX_CHECK_DIV (nrm);
   c = 1.f / nrm;
   v[0] *= c; v[1] *= c; v[2] *= c;

   SX_CHECK_NUM (v[0]);
   SX_CHECK_NUM (v[1]);
   SX_CHECK_NUM (v[2]);
   
   return *this;
}


template<class T>
void SxVector3<T>::print () const
{
   sxprintf ("[%f, %f, %f]\n", (float)v[0], (float)v[1], (float)v[2]);
}

template<class T>
typename T::Real SxVector3<T>::normSqr () const
{
   typename T::Real res = ::absSqr(v[0]) + ::absSqr(v[1]) + ::absSqr(v[2]);
   return res;
}

template<class T>
typename T::Real SxVector3<T>::norm () const
{
   return sqrt (normSqr ());
}

template <>
inline int SxVector3<Int>::norm () const
{
   // There is no sqrt for integers
   SX_EXIT;
   return -1;
}

template <class T>
typename T::Type SxVector3<T>::minval () const
{
   if (v[0] < v[1])  {  // v[0] < v[1]
      if (v[0] < v[2])  return v[0];
      else              return v[2];
   }  else  {           // v[1] <= v[0]
      if (v[1] < v[2])  return v[1];
      else              return v[2];
   }
}

// --- template specializations for complex vectors
template<>
inline SxComplex8 SxVector3<Complex8>::minval () const
{
   SX_EXIT; // complex numbers do not have minval
   return SxComplex8 ();
}

template<>
inline SxComplex16 SxVector3<Complex16>::minval () const
{
   SX_EXIT; // complex numbers do not have minval
   return SxComplex16 ();
}

template <class T>
inline typename T::Type SxVector3<T>::maxval () const
{
   if (v[0] > v[1])  {  // v[0] > v[1]
      if (v[0] > v[2])  return v[0];
      else              return v[2];
   }  else  {           // v[1] >= v[0]
      if (v[1] > v[2])  return v[1];
      else              return v[2];
   }
}

// --- template specializations for complex vectors
template<>
inline SxComplex8 SxVector3<Complex8>::maxval () const
{
   SX_EXIT; // complex numbers do not have maxval
   return SxComplex8 ();
}

template<>
inline SxComplex16 SxVector3<Complex16>::maxval () const
{
   SX_EXIT; // complex numbers do not have maxval
   return SxComplex16 ();
}


template <class T>
SxVector3<typename SxTypeCaster<Int, T>::TRes>
SxVector3<T>::operator+ (int i)
{
   return SxVector3<typename SxTypeCaster<Int, T>::TRes>
          (v[0] + i, v[1] + i, v[2] + i);
}
   
//===========================================================================
// The following functions are globally defined and will go to a
// SxOperators.h declaration file.
//===========================================================================






// --- SxVector3<?> = SxVector3<A> + SxVector3<B>
template<class A,class B>
inline SxVector3<typename SxTypeCaster<A,B>::TRes> 
operator+ (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<typename SxTypeCaster<A,B>::TRes> (
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[0] 
         + (typename SxTypeCaster<A,B>::TRes::Type)b.v[0], 
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[1] 
         + (typename SxTypeCaster<A,B>::TRes::Type)b.v[1], 
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[2] 
         + (typename SxTypeCaster<A,B>::TRes::Type)b.v[2]);
}





// --- SxVector3<?> = SxVector3<A> - SxVector3<B>
template<class A,class B>
inline SxVector3<typename SxTypeCaster<A,B>::TRes> 
operator- (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<typename SxTypeCaster<A,B>::TRes> (
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[0] 
         - (typename SxTypeCaster<A,B>::TRes::Type)b.v[0], 
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[1] 
         - (typename SxTypeCaster<A,B>::TRes::Type)b.v[1], 
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[2] 
         - (typename SxTypeCaster<A,B>::TRes::Type)b.v[2]);
}







// --- SxVector3<?> = int * SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Int,B>::TRes> 
operator* (int a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Int,B>::TRes::Type s;
   s = (typename SxTypeCaster<Int,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Int,B>::TRes> (
         s * (typename SxTypeCaster<Int,B>::TRes::Type)b.v[0], 
         s * (typename SxTypeCaster<Int,B>::TRes::Type)b.v[1], 
         s * (typename SxTypeCaster<Int,B>::TRes::Type)b.v[2]);
}

// --- SxVector3<?> = float * SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Float,B>::TRes> 
operator* (float a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Float,B>::TRes::Type s;
   s = (typename SxTypeCaster<Float,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Float,B>::TRes> (
         s * (typename SxTypeCaster<Float,B>::TRes::Type)b.v[0], 
         s * (typename SxTypeCaster<Float,B>::TRes::Type)b.v[1], 
         s * (typename SxTypeCaster<Float,B>::TRes::Type)b.v[2]);
}

// --- SxVector3<?> = double * SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Double,B>::TRes> 
operator* (double a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Double,B>::TRes::Type s;
   s = (typename SxTypeCaster<Double,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Double,B>::TRes> (
         s * (typename SxTypeCaster<Double,B>::TRes::Type)b.v[0], 
         s * (typename SxTypeCaster<Double,B>::TRes::Type)b.v[1], 
         s * (typename SxTypeCaster<Double,B>::TRes::Type)b.v[2]);
}

// --- SxVector3<?> = SxComplex8 * SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Complex8,B>::TRes> 
operator* (SxComplex8 a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Complex8,B>::TRes::Type s;
   s = (typename SxTypeCaster<Complex8,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Complex8,B>::TRes> (
         s * (typename SxTypeCaster<Complex8,B>::TRes::Type)b.v[0], 
         s * (typename SxTypeCaster<Complex8,B>::TRes::Type)b.v[1], 
         s * (typename SxTypeCaster<Complex8,B>::TRes::Type)b.v[2]);
}

// --- SxVector3<?> = SxComplex16 * SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Complex16,B>::TRes> 
operator* (SxComplex16 a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Complex16,B>::TRes::Type s;
   s = (typename SxTypeCaster<Complex16,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Complex16,B>::TRes> (
         s * (typename SxTypeCaster<Complex16,B>::TRes::Type)b.v[0], 
         s * (typename SxTypeCaster<Complex16,B>::TRes::Type)b.v[1], 
         s * (typename SxTypeCaster<Complex16,B>::TRes::Type)b.v[2]);
}


// --- SxVector3<?> = SxVector3<B> * int
template<class A>
inline SxVector3<typename SxTypeCaster<A,Int>::TRes> 
operator* (const SxVector3<A> &a, int b)
{
   typename SxTypeCaster<A,Int>::TRes::Type s;
   s = (typename SxTypeCaster<A,Int>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Int>::TRes> (
         (typename SxTypeCaster<A,Int>::TRes::Type)a.v[0] * s, 
         (typename SxTypeCaster<A,Int>::TRes::Type)a.v[1] * s, 
         (typename SxTypeCaster<A,Int>::TRes::Type)a.v[2] * s);
}


// --- SxVector3<?> = SxVector3<B> * float
template<class A>
inline SxVector3<typename SxTypeCaster<A,Float>::TRes> 
operator* (const SxVector3<A> &a, float b)
{
   typename SxTypeCaster<A,Float>::TRes::Type s;
   s = (typename SxTypeCaster<A,Float>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Float>::TRes> (
         (typename SxTypeCaster<A,Float>::TRes::Type)a.v[0] * s, 
         (typename SxTypeCaster<A,Float>::TRes::Type)a.v[1] * s, 
         (typename SxTypeCaster<A,Float>::TRes::Type)a.v[2] * s);
}


// --- SxVector3<?> = SxVector3<B> * double
template<class A>
inline SxVector3<typename SxTypeCaster<A,Double>::TRes> 
operator* (const SxVector3<A> &a, double b)
{
   typename SxTypeCaster<A,Double>::TRes::Type s;
   s = (typename SxTypeCaster<A,Double>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Double>::TRes> (
         (typename SxTypeCaster<A,Double>::TRes::Type)a.v[0] * s, 
         (typename SxTypeCaster<A,Double>::TRes::Type)a.v[1] * s, 
         (typename SxTypeCaster<A,Double>::TRes::Type)a.v[2] * s);
}


// --- SxVector3<?> = SxVector3<B> * SxComplex8
template<class A>
inline SxVector3<typename SxTypeCaster<A,Complex8>::TRes> 
operator* (const SxVector3<A> &a, SxComplex8 b)
{
   typename SxTypeCaster<A,Complex8>::TRes::Type s;
   s = (typename SxTypeCaster<A,Complex8>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Complex8>::TRes> (
         (typename SxTypeCaster<A,Complex8>::TRes::Type)a.v[0] * s, 
         (typename SxTypeCaster<A,Complex8>::TRes::Type)a.v[1] * s, 
         (typename SxTypeCaster<A,Complex8>::TRes::Type)a.v[2] * s);
}


// --- SxVector3<?> = SxVector3<B> * SxComplex16
template<class A>
inline SxVector3<typename SxTypeCaster<A,Complex16>::TRes> 
operator* (const SxVector3<A> &a, SxComplex16 b)
{
   typename SxTypeCaster<A,Complex16>::TRes::Type s;
   s = (typename SxTypeCaster<A,Complex16>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Complex16>::TRes> (
         (typename SxTypeCaster<A,Complex16>::TRes::Type)a.v[0] * s, 
         (typename SxTypeCaster<A,Complex16>::TRes::Type)a.v[1] * s, 
         (typename SxTypeCaster<A,Complex16>::TRes::Type)a.v[2] * s);
}



// --- SxVector3<?> = SxVector3<A> * SxVector3<B>
template<class A,class B>
inline SxVector3<typename SxTypeCaster<A,B>::TRes> 
operator* (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<typename SxTypeCaster<A,B>::TRes> (
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[0] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[0], 
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[1] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[1], 
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[2] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[2]);
}







// --- SxVector3<?> = int / SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Int,B>::TRes> 
operator/ (int a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Int,B>::TRes::Type s;
   s = (typename SxTypeCaster<Int,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Int,B>::TRes> (
         s / (typename SxTypeCaster<Int,B>::TRes::Type)b.v[0], 
         s / (typename SxTypeCaster<Int,B>::TRes::Type)b.v[1], 
         s / (typename SxTypeCaster<Int,B>::TRes::Type)b.v[2]);
}

// --- SxVector3<?> = float / SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Float,B>::TRes> 
operator/ (float a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Float,B>::TRes::Type s;
   s = (typename SxTypeCaster<Float,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Float,B>::TRes> (
         s / (typename SxTypeCaster<Float,B>::TRes::Type)b.v[0], 
         s / (typename SxTypeCaster<Float,B>::TRes::Type)b.v[1], 
         s / (typename SxTypeCaster<Float,B>::TRes::Type)b.v[2]);
}

// --- SxVector3<?> = double / SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Double,B>::TRes> 
operator/ (double a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Double,B>::TRes::Type s;
   s = (typename SxTypeCaster<Double,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Double,B>::TRes> (
         s / (typename SxTypeCaster<Double,B>::TRes::Type)b.v[0], 
         s / (typename SxTypeCaster<Double,B>::TRes::Type)b.v[1], 
         s / (typename SxTypeCaster<Double,B>::TRes::Type)b.v[2]);
}

// --- SxVector3<?> = SxComplex8 / SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Complex8,B>::TRes> 
operator/ (SxComplex8 a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Complex8,B>::TRes::Type s;
   s = (typename SxTypeCaster<Complex8,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Complex8,B>::TRes> (
         s / (typename SxTypeCaster<Complex8,B>::TRes::Type)b.v[0], 
         s / (typename SxTypeCaster<Complex8,B>::TRes::Type)b.v[1], 
         s / (typename SxTypeCaster<Complex8,B>::TRes::Type)b.v[2]);
}

// --- SxVector3<?> = SxComplex16 / SxVector3<B>
template<class B>
inline SxVector3<typename SxTypeCaster<Complex16,B>::TRes> 
operator/ (SxComplex16 a, const SxVector3<B> &b)
{
   typename SxTypeCaster<Complex16,B>::TRes::Type s;
   s = (typename SxTypeCaster<Complex16,B>::TRes::Type)a;
   return SxVector3<typename SxTypeCaster<Complex16,B>::TRes> (
         s / (typename SxTypeCaster<Complex16,B>::TRes::Type)b.v[0], 
         s / (typename SxTypeCaster<Complex16,B>::TRes::Type)b.v[1], 
         s / (typename SxTypeCaster<Complex16,B>::TRes::Type)b.v[2]);
}


// --- SxVector3<?> = SxVector3<B> / int
template<class A>
inline SxVector3<typename SxTypeCaster<A,Int>::TRes> 
operator/ (const SxVector3<A> &a, int b)
{
   typename SxTypeCaster<A,Int>::TRes::Type s;
   s = (typename SxTypeCaster<A,Int>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Int>::TRes> (
         (typename SxTypeCaster<A,Int>::TRes::Type)a.v[0] / s, 
         (typename SxTypeCaster<A,Int>::TRes::Type)a.v[1] / s, 
         (typename SxTypeCaster<A,Int>::TRes::Type)a.v[2] / s);
}


// --- SxVector3<?> = SxVector3<B> / float
template<class A>
inline SxVector3<typename SxTypeCaster<A,Float>::TRes> 
operator/ (const SxVector3<A> &a, float b)
{
   typename SxTypeCaster<A,Float>::TRes::Type s;
   s = (typename SxTypeCaster<A,Float>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Float>::TRes> (
         (typename SxTypeCaster<A,Float>::TRes::Type)a.v[0] / s, 
         (typename SxTypeCaster<A,Float>::TRes::Type)a.v[1] / s, 
         (typename SxTypeCaster<A,Float>::TRes::Type)a.v[2] / s);
}


// --- SxVector3<?> = SxVector3<B> / double
template<class A>
inline SxVector3<typename SxTypeCaster<A,Double>::TRes> 
operator/ (const SxVector3<A> &a, double b)
{
   typename SxTypeCaster<A,Double>::TRes::Type s;
   s = (typename SxTypeCaster<A,Double>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Double>::TRes> (
         (typename SxTypeCaster<A,Double>::TRes::Type)a.v[0] / s, 
         (typename SxTypeCaster<A,Double>::TRes::Type)a.v[1] / s, 
         (typename SxTypeCaster<A,Double>::TRes::Type)a.v[2] / s);
}


// --- SxVector3<?> = SxVector3<B> / SxComplex8
template<class A>
inline SxVector3<typename SxTypeCaster<A,Complex8>::TRes> 
operator/ (const SxVector3<A> &a, SxComplex8 b)
{
   typename SxTypeCaster<A,Complex8>::TRes::Type s;
   s = (typename SxTypeCaster<A,Complex8>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Complex8>::TRes> (
         (typename SxTypeCaster<A,Complex8>::TRes::Type)a.v[0] / s, 
         (typename SxTypeCaster<A,Complex8>::TRes::Type)a.v[1] / s, 
         (typename SxTypeCaster<A,Complex8>::TRes::Type)a.v[2] / s);
}


// --- SxVector3<?> = SxVector3<B> / SxComplex16
template<class A>
inline SxVector3<typename SxTypeCaster<A,Complex16>::TRes> 
operator/ (const SxVector3<A> &a, SxComplex16 b)
{
   typename SxTypeCaster<A,Complex16>::TRes::Type s;
   s = (typename SxTypeCaster<A,Complex16>::TRes::Type)b;
   return SxVector3<typename SxTypeCaster<A,Complex16>::TRes> (
         (typename SxTypeCaster<A,Complex16>::TRes::Type)a.v[0] / s, 
         (typename SxTypeCaster<A,Complex16>::TRes::Type)a.v[1] / s, 
         (typename SxTypeCaster<A,Complex16>::TRes::Type)a.v[2] / s);
}



// --- SxVector3<?> = SxVector3<A> / SxVector3<B>
template<class A,class B>
inline SxVector3<typename SxTypeCaster<A,B>::TRes> 
operator/ (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<typename SxTypeCaster<A,B>::TRes> (
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[0] 
         / (typename SxTypeCaster<A,B>::TRes::Type)b.v[0], 
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[1] 
         / (typename SxTypeCaster<A,B>::TRes::Type)b.v[1], 
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[2] 
         / (typename SxTypeCaster<A,B>::TRes::Type)b.v[2]);
}



// --- <?> = SxVector3<A> ^ SxVector3<B>
template<class A,class B>
typename SxTypeCaster<A,B>::Type
operator^ (const SxVector3<A> &a, const SxVector3<B> &b)
{
   typename SxTypeCaster<A,B>::TRes::Type scp;
   scp = (typename SxTypeCaster<A,B>::TRes::Type)a.v[0] 
       * (typename SxTypeCaster<A,B>::TRes::Type)b.v[0]
       + (typename SxTypeCaster<A,B>::TRes::Type)a.v[1] 
       * (typename SxTypeCaster<A,B>::TRes::Type)b.v[1]  
       + (typename SxTypeCaster<A,B>::TRes::Type)a.v[2] 
       * (typename SxTypeCaster<A,B>::TRes::Type)b.v[2];
   return scp;
}


template<class A,class B>
typename SxTypeCaster<A,B>::Type
dot (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return a ^ b;
}


template<class A,class B>
inline SxVector3<typename SxTypeCaster<A,B>::TRes> 
cross (const SxVector3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<typename SxTypeCaster<A,B>::TRes> (
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[1] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[2]
         - (typename SxTypeCaster<A,B>::TRes::Type)a.v[2] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[1],
         
           (typename SxTypeCaster<A,B>::TRes::Type)a.v[2] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[0]
         - (typename SxTypeCaster<A,B>::TRes::Type)a.v[0] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[2],

           (typename SxTypeCaster<A,B>::TRes::Type)a.v[0] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[1]
         - (typename SxTypeCaster<A,B>::TRes::Type)a.v[1] 
         * (typename SxTypeCaster<A,B>::TRes::Type)b.v[0]);
}

/// int + SxVector3<T>
template <class T>
SxVector3<typename SxTypeCaster<Int, T>::TRes>
operator+ (int i, const SxVector3<T> &vec)
{
   return SxVector3<typename SxTypeCaster<Int, T>::TRes>
          (vec.v[0] + i, vec.v[1] + i, vec.v[2] + i);
}

inline SxVector3<Int> operator%(const SxVector3<Int> &a, 
                                const SxVector3<Int> &b)
{
   SX_CHECK (b(0) != 0 && b(1) != 0 && b(2) != 0, b(0), b(1), b(2));
   return SxVector3<Int> (a(0) % b(0), a(1) % b(1), a(2) % b(2));
}

inline SxVector3<Int> operator%(const SxVector3<Int> &a, int b)
{
   SX_CHECK (b != 0);
   return SxVector3<Int> (a(0) % b, a(1) % b, a(2) % b);
}

inline SxVector3<Int> operator%(int a, const SxVector3<Int> &b)
{
   return SxVector3<Int> (a % b(0), a % b(1), a % b(2));
}

inline SxVector3<Int> & operator%=(SxVector3<Int> &a, const SxVector3<Int> &b)
{
   a(0) %= b(0); a(1) %= b(1); a(2) %= b(2);
   return a;
}

inline SxVector3<Int> & operator%=(SxVector3<Int> &a, int b)
{
   a(0) %= b; a(1) %= b; a(2) %= b;
   return a;
}

//------------------------------------------------------------------------------
inline float getAngle (const SxVector3<Float> &a, const SxVector3<Float> &b)
{
   float angle = acosf((a ^ b) / sqrtf((a^a) * (b^b)));
   if (angle < 0.f) angle += static_cast<float>(PI);
   return angle;
}

inline double getAngle (const SxVector3<Double> &a, const SxVector3<Double> &b)
{
   double angle = acos((a ^ b) / sqrt((a^a) * (b^b)));
   if (angle < 0.) angle += PI;
   return angle;
}

//------------------------------------------------------------------------------

inline SxVector3<Double> round(const SxVector3<Double> &x)
{
   return SxVector3<Double> (round(x(0)), round(x(1)), round(x(2)));
}

inline SxVector3<Double> floor(const SxVector3<Double> &x)
{
   return SxVector3<Double> (floor(x(0)), floor(x(1)), floor(x(2)));
}

inline SxVector3<Double> ceil(const SxVector3<Double> &x)
{
   return SxVector3<Double> (ceil(x(0)), ceil(x(1)), ceil(x(2)));
}


//------------------------------------------------------------------------------


template<class T>
inline std::ostream& operator<< (std::ostream &s, const SxVector3<T> &in)
{
   return s << "{" << in.v[0] << "," << in.v[1] << "," << in.v[2] << "}";
}


#endif // _SX_VECTOR3_H_

