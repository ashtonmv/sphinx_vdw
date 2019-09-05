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

#ifndef _SX_MATRIX3_H_
#define _SX_MATRIX3_H_

#include <SxVector3.h>
#include <SxComplex.h>
#include <SxError.h>
#include <SxList.h>
#include <SxMathLib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <SxString.h>

/** This class defines operators for simple 3x3 matrices. 
    Instead of the common C order of arrays here the Fortran style, 
    i.e. first dimension varies fastest, is used in order to stay consistent
    with the multidimensional SxMatrix class. The Fortran style has to be
    used since the actually matrix operations will be performed by BLAS
    function calls. BLAS was written in FORTRAN and expected array agruments
    in FORTRAN style.

    @ingroup Numerics
    @author Sixten Boeck
    @see    SxMatrix<T>
    @see    SxVector<T>
 */
template<class T>
class SxMatrix3
{
   public:
      typename T::Type m[3][3];

      inline SxMatrix3 ();
      inline SxMatrix3 (const typename T::Type &val);
      inline SxMatrix3 (const SxMatrix3<Int> &mat);
      inline SxMatrix3 (const SxMatrix3<Float> &mat);
      inline SxMatrix3 (const SxMatrix3<Double> &mat);
      inline SxMatrix3 (const SxMatrix3<Complex8> &mat);
      inline SxMatrix3 (const SxMatrix3<Complex16> &mat);
      inline SxMatrix3 (const SxVector3<T> &v1,
                        const SxVector3<T> &v2,
                        const SxVector3<T> &v3);
      inline SxMatrix3 (const typename T::Type &e00, 
                        const typename T::Type &e01, 
                        const typename T::Type &e02,
                        const typename T::Type &e10, 
                        const typename T::Type &e11, 
                        const typename T::Type &e12,
                        const typename T::Type &e20, 
                        const typename T::Type &e21, 
                        const typename T::Type &e22);
      inline explicit SxMatrix3 (const SxList<typename T::Type> &);

      inline void set (const typename T::Type &);
      inline void setCol (ssize_t i, const SxVector3<T> &colIn);

      /** Use matrix(x,y) in order to dereference an element */
      inline SxVector3<T>        operator() (ssize_t i);
      inline const SxVector3<T>  operator() (ssize_t i) const;
      /** i and j are the row and column, resp. */
      inline typename T::Type       &operator() (ssize_t i, ssize_t j);
      /** i and j are the row and column, resp. */
      inline const typename T::Type &operator() (ssize_t i, ssize_t j) const;
      /** i and j are the row and column, resp. */
      inline typename T::Type       &operator() (const SxAutoLoop &i,
                                                 const SxAutoLoop &j);
      /** i and j are the row and column, resp. */
      inline const typename T::Type &operator() (const SxAutoLoop &i,
                                                 const SxAutoLoop &j) const;

      /** Assign a scalar value to all matrix elements */
      inline SxMatrix3<T>       &operator=  (const typename T::Type &s);
      /** Assign matrix to matrix */
      inline SxMatrix3<T>       &operator=  (const SxMatrix3<T> &in);
      /** Multiply matrix by a scalar value */
      //inline SxMatrix3<T> operator*  (const T &s) const;
      inline SxMatrix3<T> operator*= (const typename T::Type &s);
      /** \brief Add a matrix elementwise */
      SxMatrix3<T> operator+= (const SxMatrix3<T> &in);
      /** \brief Subtract a matrix elementwise */
      SxMatrix3<T> operator-= (const SxMatrix3<T> &in);
      /** Divide matrix by a scalar value */
      inline SxMatrix3<T> operator/= (const typename T::Type &s);
      inline bool operator== (const SxMatrix3<T> &) const;
      inline bool operator!= (const SxMatrix3<T> &in) const;
      inline SxVector3<T>     row (ssize_t i) const;
      inline SxVector3<T>     col (ssize_t i) const;
      inline typename T::Type sum () const;
      inline typename T::Type tr  () const;
      inline SxMatrix3<typename T::TReal> absSqr () const;

      /** Evaluate inverse of matrix */
      inline SxMatrix3<T> inverse () const;
      /** Compute determinant of matrix */
      inline typename T::Type determinant () const;
      /** Return transpose of matrix */
      inline SxMatrix3<T> transpose () const;
      /** Return trace of matrix */
      inline typename T::Type trace () const;

      void print () const;
};


template<class T>
SxMatrix3<T>::SxMatrix3 ()
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (typename T::Type)0.0;
}


template<class T>
SxMatrix3<T>::SxMatrix3 (const typename T::Type &val)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = val;
}



template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<Int> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (typename T::Type)mat.m[i][j];
}
template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<Float> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (typename T::Type)mat.m[i][j];
}
template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<Double> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (typename T::Type)mat.m[i][j];
}
template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<Complex8> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (typename T::Type)mat.m[i][j];
}
template<class T>
SxMatrix3<T>::SxMatrix3 (const SxMatrix3<Complex16> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = (typename T::Type)mat.m[i][j];
}
template<>
inline SxMatrix3<Int>::SxMatrix3 (const SxMatrix3<Float> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = toInt(mat.m[i][j]);
}
template<>
inline SxMatrix3<Int>::SxMatrix3 (const SxMatrix3<Double> &mat)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = toInt(mat.m[i][j]);
}



template<class T>
SxMatrix3<T>::SxMatrix3 (const SxVector3<T> &v1,
                         const SxVector3<T> &v2,
                         const SxVector3<T> &v3)
{
    m[0][0] = v1.v[0]; m[1][0] = v2.v[0]; m[2][0] = v3.v[0];
    m[0][1] = v1.v[1]; m[1][1] = v2.v[1]; m[2][1] = v3.v[1];
    m[0][2] = v1.v[2]; m[1][2] = v2.v[2]; m[2][2] = v3.v[2];
}


template<class T>
SxMatrix3<T>::SxMatrix3 (const typename T::Type &e00, 
                         const typename T::Type &e01, 
                         const typename T::Type &e02,
                         const typename T::Type &e10, 
                         const typename T::Type &e11, 
                         const typename T::Type &e12,
                         const typename T::Type &e20, 
                         const typename T::Type &e21, 
                         const typename T::Type &e22)
{
   m[0][0] = e00; m[0][1] = e01; m[0][2] = e02;
   m[1][0] = e10; m[1][1] = e11; m[1][2] = e12;
   m[2][0] = e20; m[2][1] = e21; m[2][2] = e22;
}

template<class T>
SxMatrix3<T>::SxMatrix3 (const SxList<typename T::Type> &in)
{
   SX_CHECK (in.getSize () == 9, in.getSize ());
   typename SxList<typename T::Type>::ConstIterator it = in.begin();
   m[0][0] = *it++; m[0][1] = *it++; m[0][2] = *it++;
   m[1][0] = *it++; m[1][1] = *it++; m[1][2] = *it++;
   m[2][0] = *it++; m[2][1] = *it++; m[2][2] = *it;
}


template<class T>
void SxMatrix3<T>::set (const typename T::Type &v)
{
   m[0][0] = m[0][1] = m[0][2] =
   m[1][0] = m[1][1] = m[1][2] =
   m[2][0] = m[2][1] = m[2][2] = v;
}

template<class T>
void SxMatrix3<T>::setCol (ssize_t i, const SxVector3<T> &colIn)
{
   m[0][i] = colIn(0);
   m[1][i] = colIn(1);
   m[2][i] = colIn(2);
}


template<class T>
SxVector3<T> SxMatrix3<T>::operator() (ssize_t i)
{
   SX_CHECK (i >= 0 && i < 3, i);
   return SxVector3<T> (m[0][i], m[1][i], m[2][i]);
}

template<class T>
const SxVector3<T> SxMatrix3<T>::operator() (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < 3, i);
   return SxVector3<T> (m[0][i], m[1][i], m[2][i]);
}


template<class T>
typename T::Type &SxMatrix3<T>::operator() (ssize_t i, ssize_t j)
{
   SX_CHECK (i >= 0 && i < 3 && j >= 0 && j < 3, i, j);
   return m[i][j];
}


template<class T>
const typename T::Type& SxMatrix3<T>::operator() (ssize_t i, ssize_t j) const
{
   SX_CHECK (i >= 0 && i < 3 && j >= 0 && j < 3, i, j);
   return m[i][j];
}

template<class T>
typename T::Type &SxMatrix3<T>::operator() (const SxAutoLoop& i,
                                            const SxAutoLoop& j)
{
   i.setLimit (3);
   j.setLimit (3);
   return m[i.i][j.i];
}


template<class T>
const typename T::Type& SxMatrix3<T>::operator() (const SxAutoLoop& i,
                                                  const SxAutoLoop& j) const
{
   i.setLimit (3);
   j.setLimit (3);
   return m[i.i][j.i];
}



template<class T>
SxMatrix3<T> &SxMatrix3<T>::operator= (const typename T::Type &s)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = s;
   return *this;
}

template<class T>
SxMatrix3<T> &SxMatrix3<T>::operator= (const SxMatrix3<T> &in)
{
   if (&in == this)  return *this;

   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] = in.m[i][j];
   return *this;
}


template<class T>
SxMatrix3<T> SxMatrix3<T>::operator*= (const typename T::Type &s)
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] *= s;
   return *this;
}

template<class T>
SxMatrix3<T> SxMatrix3<T>::operator/= (const typename T::Type &s)
{
   // doesn't work properly for SxComplex<?>:
   // SX_CHECK ( fabs ((double)s) > 1e-10 , (double)s);
   int i, j;
   typename T::Type d = 1. / s;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         m[i][j] *= d;
   return *this;
}

template <class T>
SxMatrix3<T> SxMatrix3<T>::operator+= (const SxMatrix3<T> &in)
{
   int i,j;
   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
         this->m[i][j] += in.m[i][j];
   return *this;
}

template <class T>
SxMatrix3<T> SxMatrix3<T>::operator-= (const SxMatrix3<T> &in)
{
   int i,j;
   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
         this->m[i][j] -= in.m[i][j];
   return *this;
}

template<class T>
bool SxMatrix3<T>::operator!= (const SxMatrix3<T> &in) const
{
   return !operator==(in);
}

template<class T>
bool SxMatrix3<T>::operator== (const SxMatrix3<T> &in) const
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         if (m[i][j] != in.m[i][j] )   return false;
   return true;
}


template<>
inline bool SxMatrix3<Float>::operator== (const SxMatrix3<Float> &in) const
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         if ( fabs(m[i][j] - in.m[i][j]) > 1e-10 )  return false;
   return true;
}


template<>
inline bool SxMatrix3<Double>::operator== (const SxMatrix3<Double> &in) const
{
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         if ( fabs(m[i][j] - in.m[i][j]) > 1e-10 )  return false;
   return true;
}



template<class T>
SxVector3<T> SxMatrix3<T>::row (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < 3, i);
   return   SxVector3<T> (m[i][0], m[i][1], m[i][2]);
}


template<class T>
SxVector3<T> SxMatrix3<T>::col (ssize_t i) const
{
   SX_CHECK (i >= 0 && i < 3, i);
   return   SxVector3<T> (m[0][i], m[1][i], m[2][i]);
}




template<class T>
typename T::Type SxMatrix3<T>::sum () const
{
   return   m[0][0] + m[0][1] + m[0][2]
          + m[1][0] + m[1][1] + m[1][2]
          + m[2][0] + m[2][1] + m[2][2];
}

template<class T>
typename T::Type SxMatrix3<T>::tr () const
{
   return   m[0][0] + m[1][1] + m[2][2];
}


template<>
inline SxMatrix3<Int::TReal> SxMatrix3<Int>::absSqr () const
{
   SxMatrix3<Int::TReal> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j] * m[i][j];
   return res;
}

template<>
inline SxMatrix3<Float::TReal> SxMatrix3<Float>::absSqr () const
{
   SxMatrix3<Float::TReal> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j] * m[i][j];
   return res;
}

template<>
inline SxMatrix3<Double::TReal> SxMatrix3<Double>::absSqr () const
{
   SxMatrix3<Double::TReal> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j] * m[i][j];
   return res;
}

template<>
inline SxMatrix3<Complex8::TReal> SxMatrix3<Complex8>::absSqr () const
{
   SxMatrix3<Complex8::TReal> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j].absSqr ();
   return res;
}

template<>
inline SxMatrix3<Complex16::TReal> SxMatrix3<Complex16>::absSqr () const
{
   SxMatrix3<Complex16::TReal> res;
   int i, j;
   for (i=0; i < 3; i++)
      for (j=0; j < 3; j++)
         res.m[i][j] = m[i][j].absSqr ();
   return res;
}


template<class T>
typename T::Type SxMatrix3<T>::determinant() const
{
  return (m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])+
          m[1][0] * (m[2][1] * m[0][2] - m[0][1] * m[2][2])+
	       m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]));
}

template<class T>
typename T::Type SxMatrix3<T>::trace () const
{
   return (m[0][0] + m[1][1] + m[2][2]);
}


template<class T>
SxMatrix3<T> SxMatrix3<T>::transpose () const
{
   SxMatrix3<T> trans;
   int i,j;
   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
         trans.m[i][j] = m[j][i];
   
   return trans;
}


template<class T>
SxMatrix3<T> SxMatrix3<T>::inverse () const
{
  const SxMatrix3<T> &M = *this;

  SxMatrix3<T> inv (M(1).x (M(2)), 
                    M(2).x (M(0)), 
                    M(0).x (M(1)));

  SX_CHECK (fabs((double)determinant()) > 1e-10, (double)determinant());
  inv *= ((double)1.) / determinant ();
  return inv;
}

template<class T>
void SxMatrix3<T>::print () const
{
   sxprintf ("[[%f, %f, %f]\n",  (float)m[0][0], (float)m[0][1], (float)m[0][2]);
   sxprintf (" [%f, %f, %f]\n",  (float)m[1][0], (float)m[1][1], (float)m[1][2]);
   sxprintf (" [%f, %f, %f]]\n", (float)m[2][0], (float)m[2][1], (float)m[2][2]);
}



//===========================================================================
// The following functions are globally defined and will go to a
// SxOperators.h declaration file.
//===========================================================================


// SxMatrix3<?> = int * SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Int,B>::TRes> 
operator* (int a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Int,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Int,B>::TRes::Type)a
                     * (typename SxTypeCaster<Int,B>::TRes::Type)b.m[i][j];
   return res;
}


// SxMatrix3<?> = float * SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Float,B>::TRes> 
operator* (float a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Float,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Float,B>::TRes::Type)a
                     * (typename SxTypeCaster<Float,B>::TRes::Type)b.m[i][j];
   return res;
}


// SxMatrix3<?> = double * SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Double,B>::TRes> 
operator* (double a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Double,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Double,B>::TRes::Type)a
                     * (typename SxTypeCaster<Double,B>::TRes::Type)b.m[i][j];
   return res;
}


// SxMatrix3<?> = SxComplex8 * SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Complex8,B>::TRes> 
operator* (SxComplex8 a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Complex8,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Complex8,B>::TRes::Type)a
                     * (typename SxTypeCaster<Complex8,B>::TRes::Type)b.m[i][j];
   return res;
}


// SxMatrix3<?> = SxComplex16 * SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Complex16,B>::TRes> 
operator* (SxComplex16 a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Complex16,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Complex16,B>::TRes::Type)a
                     * (typename SxTypeCaster<Complex16,B>::TRes::Type)b.m[i][j];
   return res;
}






// SxMatrix3<?> = SxMatrix3<A> * int
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Int>::TRes> 
operator* (const SxMatrix3<A> &a, int b)
{
   SxMatrix3<typename SxTypeCaster<A,Int>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Int>::TRes::Type)a.m[i][j] 
                     * (typename SxTypeCaster<A,Int>::TRes::Type)b;
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> * float
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Float>::TRes> 
operator* (const SxMatrix3<A> &a, float b)
{
   SxMatrix3<typename SxTypeCaster<A,Float>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Float>::TRes::Type)a.m[i][j] 
                     * (typename SxTypeCaster<A,Float>::TRes::Type)b;
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> * double
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Double>::TRes> 
operator* (const SxMatrix3<A> &a, double b)
{
   SxMatrix3<typename SxTypeCaster<A,Double>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Double>::TRes::Type)a.m[i][j] 
                     * (typename SxTypeCaster<A,Double>::TRes::Type)b;
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> * SxComplex8
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Complex8>::TRes> 
operator* (const SxMatrix3<A> &a, SxComplex8 b)
{
   SxMatrix3<typename SxTypeCaster<A,Complex8>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Complex8>::TRes::Type)a.m[i][j] 
                     * (typename SxTypeCaster<A,Complex8>::TRes::Type)b;
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> * SxComplex16
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Complex16>::TRes> 
operator* (const SxMatrix3<A> &a, SxComplex16 b)
{
   SxMatrix3<typename SxTypeCaster<A,Complex16>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Complex16>::TRes::Type)a.m[i][j] 
                     * (typename SxTypeCaster<A,Complex16>::TRes::Type)b;
   return res;
}





// SxVector3<?> = SxMatrix3<A> ^ SxVector3<B>
template<class A,class B>
inline SxVector3<typename SxTypeCaster<A,B>::TRes> 
operator^ (const SxMatrix3<A> &a, const SxVector3<B> &b)
{
   return SxVector3<typename SxTypeCaster<A,B>::TRes> ( 
                   (typename SxTypeCaster<A,B>::TRes::Type)a.m[0][0]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[0]
                  +(typename SxTypeCaster<A,B>::TRes::Type)a.m[0][1]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[1]
                  +(typename SxTypeCaster<A,B>::TRes::Type)a.m[0][2]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[2],
                  
                   (typename SxTypeCaster<A,B>::TRes::Type)a.m[1][0]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[0]
                  +(typename SxTypeCaster<A,B>::TRes::Type)a.m[1][1]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[1]
                  +(typename SxTypeCaster<A,B>::TRes::Type)a.m[1][2]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[2],

                   (typename SxTypeCaster<A,B>::TRes::Type)a.m[2][0]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[0]
                  +(typename SxTypeCaster<A,B>::TRes::Type)a.m[2][1]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[1]
                  +(typename SxTypeCaster<A,B>::TRes::Type)a.m[2][2]
                  *(typename SxTypeCaster<A,B>::TRes::Type)b.v[2] 
          );
}

// SxVector3<?> = SxVector3<A> ^ SxMatrix3<B>
template<class A,class B>
inline SxVector3<typename SxTypeCaster<A,B>::TRes> 
operator^ (const SxVector3<A> &a, const SxMatrix3<B> &b)
{
   return SxVector3<typename SxTypeCaster<A,B>::TRes> ( 
                   (typename SxTypeCaster<A,B>::TRes::Type)b.m[0][0]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[0]
                  +(typename SxTypeCaster<A,B>::TRes::Type)b.m[1][0]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[1]
                  +(typename SxTypeCaster<A,B>::TRes::Type)b.m[2][0]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[2],
                  
                   (typename SxTypeCaster<A,B>::TRes::Type)b.m[0][1]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[0]
                  +(typename SxTypeCaster<A,B>::TRes::Type)b.m[1][1]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[1]
                  +(typename SxTypeCaster<A,B>::TRes::Type)b.m[2][1]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[2],

                   (typename SxTypeCaster<A,B>::TRes::Type)b.m[0][2]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[0]
                  +(typename SxTypeCaster<A,B>::TRes::Type)b.m[1][2]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[1]
                  +(typename SxTypeCaster<A,B>::TRes::Type)b.m[2][2]
                  *(typename SxTypeCaster<A,B>::TRes::Type)a.v[2] 
          );
}


// SxMatrix3<?> = SxMatrix3<A> ^ SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<typename SxTypeCaster<A,B>::TRes> 
operator^ (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<A,B>::TRes> res;
   int r, c;
   for (c=0; c < 3; c++)
      for (r=0; r < 3; r++)  
         res.m[r][c] = (typename SxTypeCaster<A,B>::TRes::Type)a.m[r][0] 
                     * (typename SxTypeCaster<A,B>::TRes::Type)b.m[0][c]
                     + (typename SxTypeCaster<A,B>::TRes::Type)a.m[r][1] 
                     * (typename SxTypeCaster<A,B>::TRes::Type)b.m[1][c]
                     + (typename SxTypeCaster<A,B>::TRes::Type)a.m[r][2]
                     * (typename SxTypeCaster<A,B>::TRes::Type)b.m[2][c];
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> + SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<typename SxTypeCaster<A,B>::TRes> 
operator+ (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<A,B>::TRes> res;
   int r, c;
   for (c=0; c < 3; c++)
      for (r=0; r < 3; r++)  
         res.m[r][c] = (typename SxTypeCaster<A,B>::TRes::Type)a.m[r][c] 
                     + (typename SxTypeCaster<A,B>::TRes::Type)b.m[r][c];
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> - SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<typename SxTypeCaster<A,B>::TRes> 
operator- (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<A,B>::TRes> res;
   int r, c;
   for (c=0; c < 3; c++)
      for (r=0; r < 3; r++)  
         res.m[r][c] = (typename SxTypeCaster<A,B>::TRes::Type)a.m[r][c] 
                     - (typename SxTypeCaster<A,B>::TRes::Type)b.m[r][c];
   return res;
}





// SxMatrix3<?> = SxMatrix3<A> * SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<typename SxTypeCaster<A,B>::TRes> 
operator* (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<A,B>::TRes> res;
   int r, c;
   for (c=0; c < 3; c++)
      for (r=0; r < 3; r++)  
         res.m[r][c] = (typename SxTypeCaster<A,B>::TRes::Type)a.m[r][c] 
                     * (typename SxTypeCaster<A,B>::TRes::Type)b.m[r][c];
   return res;
}






// SxMatrix3<?> = int / SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Int,B>::TRes> 
operator/ (int a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Int,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Int,B>::TRes::Type)a
                     / (typename SxTypeCaster<Int,B>::TRes::Type)b.m[i][j];
   return res;
}


// SxMatrix3<?> = float / SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Float,B>::TRes> 
operator/ (float a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Float,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Float,B>::TRes::Type)a
                     / (typename SxTypeCaster<Float,B>::TRes::Type)b.m[i][j];
   return res;
}


// SxMatrix3<?> = double / SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Double,B>::TRes> 
operator/ (double a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Double,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Double,B>::TRes::Type)a
                     / (typename SxTypeCaster<Double,B>::TRes::Type)b.m[i][j];
   return res;
}


// SxMatrix3<?> = SxComplex8 / SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Complex8,B>::TRes> 
operator/ (SxComplex8 a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Complex8,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Complex8,B>::TRes::Type)a
                     / (typename SxTypeCaster<Complex8,B>::TRes::Type)b.m[i][j];
   return res;
}


// SxMatrix3<?> = SxComplex16 / SxMatrix3<A>
template<class B>
inline SxMatrix3<typename SxTypeCaster<Complex16,B>::TRes> 
operator/ (SxComplex16 a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<Complex16,B>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<Complex16,B>::TRes::Type)a
                     / (typename SxTypeCaster<Complex16,B>::TRes::Type)b.m[i][j];
   return res;
}






// SxMatrix3<?> = SxMatrix3<A> / int
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Int>::TRes> 
operator/ (const SxMatrix3<A> &a, int b)
{
   SxMatrix3<typename SxTypeCaster<A,Int>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Int>::TRes::Type)a.m[i][j] 
                     / (typename SxTypeCaster<A,Int>::TRes::Type)b;
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> / float
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Float>::TRes> 
operator/ (const SxMatrix3<A> &a, float b)
{
   SxMatrix3<typename SxTypeCaster<A,Float>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Float>::TRes::Type)a.m[i][j] 
                     / (typename SxTypeCaster<A,Float>::TRes::Type)b;
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> / double
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Double>::TRes> 
operator/ (const SxMatrix3<A> &a, double b)
{
   SxMatrix3<typename SxTypeCaster<A,Double>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Double>::TRes::Type)a.m[i][j] 
                     / (typename SxTypeCaster<A,Double>::TRes::Type)b;
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> / SxComplex8
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Complex8>::TRes> 
operator/ (const SxMatrix3<A> &a, SxComplex8 b)
{
   SxMatrix3<typename SxTypeCaster<A,Complex8>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Complex8>::TRes::Type)a.m[i][j] 
                     / (typename SxTypeCaster<A,Complex8>::TRes::Type)b;
   return res;
}


// SxMatrix3<?> = SxMatrix3<A> / SxComplex16
template<class A>
inline SxMatrix3<typename SxTypeCaster<A,Complex16>::TRes> 
operator/ (const SxMatrix3<A> &a, SxComplex16 b)
{
   SxMatrix3<typename SxTypeCaster<A,Complex16>::TRes> res;
   int i, j;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         res.m[i][j] = (typename SxTypeCaster<A,Complex16>::TRes::Type)a.m[i][j] 
                     / (typename SxTypeCaster<A,Complex16>::TRes::Type)b;
   return res;
}



// SxMatrix3<?> = SxMatrix3<A> / SxMatrix3<B>
template<class A,class B>
inline SxMatrix3<typename SxTypeCaster<A,B>::TRes> 
operator/ (const SxMatrix3<A> &a, const SxMatrix3<B> &b)
{
   SxMatrix3<typename SxTypeCaster<A,B>::TRes> res;
   int r, c;
   for (c=0; c < 3; c++)
      for (r=0; r < 3; r++)  
         res.m[r][c] = (typename SxTypeCaster<A,B>::TRes::Type)a.m[r][c] 
                     / (typename SxTypeCaster<A,B>::TRes::Type)b.m[r][c];
   return res;
}






//------------------------------------------------------------------------------
template<class T>
std::ostream& operator<< (std::ostream &s, const SxMatrix3<T> &in)
{
   return s << "{{" << in.m[0][0] << "," << in.m[0][1] << "," << in.m[0][2] << "},"
            << " {" << in.m[1][0] << "," << in.m[1][1] << "," << in.m[1][2] << "},"
            << " {" << in.m[2][0] << "," << in.m[2][1] << "," << in.m[2][2] << "}}";
}



#endif // _SX_MATRIX3_H_
