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
#ifndef _SX_CMPLX_H_
#define _SX_CMPLX_H_

#include <cmath>

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <SxError.h>
#include <SxInfNan.h>


//using namespace std;

/** This class represents a complex number.
    If you include file SxConstants.h you also write "1. + 2.*I" instead of
    SxComplex<?> (1,2). Mind, execution of (1+2*I) is slighly slower than calling
    the constructor directly.
    @ingroup Numerics
    @author Sixten Boeck
    @version 1.0
  */
template <class T>
class SxComplex
{
   public:
      /** The real part of the complex number */
      T re;
      /** The imaginary part of the complex number */
      T im;

      // constructors
      
      inline SxComplex ()
      { 
#ifndef NDEBUG
         re = im = (T)sqrt(-1.);
#endif
      }
      inline SxComplex (const int r)                { re = (T)r; im = (T)0.0; }
      inline SxComplex (const float r)              { re = (T)r; im = (T)0.0; }
      inline SxComplex (const double r)             { re = (T)r; im = (T)0.0; }
      inline SxComplex (const T r, const T i) : re(r), im(i) { /* empty */    }
      inline SxComplex (const SxComplex<T> &c) = default;
      template <class U>
      inline SxComplex (const SxComplex<U> &c);

      /// Computes exp(I * phi)
      static inline SxComplex<T> phase(T phi);


      inline void set (const T r)            { re = r; im = (T)0.0; }
      inline void set (const T r, const T i) { re = r; im = i; }

      // 'c1 = c2'
      inline SxComplex<T> &operator= (const SxComplex<T> &c) = default;

      inline SxComplex<T> operator- ()  const {
         return SxComplex<T> (-re, -im);
      }

      inline SxComplex<T> conj () const  {
         return SxComplex<T> (re, -im);
      }

      // 'c + r'
      inline SxComplex<T> operator+ (const T &r) const {
         return SxComplex<T> ( re + r, im);
      }

      // 'c1 + c2'
      inline SxComplex<T> operator+ (const SxComplex<T> &c) const {
         return SxComplex<T> ( re + c.re, im + c.im);
      }

      // 'r + c1' globally defined

      // 'c1 += c2'
      inline void operator+= (const SxComplex<T> &c)  {
         re += c.re;
         im += c.im;
      }

      // 'c1 - r'
      inline SxComplex<T> operator- (const T &r) const {
         return SxComplex<T> ( re - r, im);
      }

      // 'c1 - c2'
      inline SxComplex<T> operator- (const SxComplex<T> &c) const {
         return SxComplex<T> ( re - c.re, im - c.im);
      }

      // 'c1 -= c2'
      inline void operator-= (const SxComplex<T> &c)  {
         re -= c.re;
         im -= c.im;
      }

      // 'r - c1' globally defined

      // 'c1 * r'
      inline SxComplex<T> operator* (const T r)  const {
         return SxComplex<T> (re * r, im * r);
      }

      // 'r * c1' globally defined

      // 'c1 * c2'
      inline SxComplex<T> operator* (const SxComplex<T> &c) const  {
         return SxComplex<T> (re * c.re  -  im * c.im,
                              re * c.im  +  im * c.re);
      }

      // 'c1 *= r'
      inline void operator*= (const T r)  {
         re *= r;
         im *= r;
      }


      // 'c1 *= c2'
      inline void operator*= (const SxComplex<T> &c)  {
         T temp = re * c.re  -  im * c.im;
         im     = re * c.im  +  im * c.re;
         re     = temp;
      }

      // 'c / r'
      inline SxComplex<T> operator/ (const T r) const {
         T t = (T)1.0 / r;
         return SxComplex<T> ( re * t, im * t);
      }

      // 'r / c' globally defined


      // 'c1 / c2'
      inline SxComplex<T> operator/ (const SxComplex<T> &c) const {
         T t = (T)1.0 / (c.re * c.re   +   c.im * c.im);
         return SxComplex<T> (t * (re * c.re  +  im * c.im),
               t * (im * c.re  -  re * c.im));
      }


      // 'c /= r'
      inline void operator/= (const T r)  {
         T t = (T)1.0 / r;
         re *= t;
         im *= t;
      }

      // 'c1 /= c2'
      inline void operator/= (const SxComplex<T> &c)  {
         T t = (T)1.0 / (c.re * c.re   +   c.im * c.im);
         T temp = t * (re * c.re  +  im * c.im);
         im          = t * (im * c.re  -  re * c.im);
         re          = temp;
      }

      // c^2
      inline T operator^ (T y) const {
         SX_CHECK ((int)y == 2, y);  // right now only y^2 is defined!
         return re*re + im*im;
      }

      inline T absSqr () const  {
         return re*re + im*im;
      }

      inline T abs  () const  {
         return ::sqrt(re*re + im*im); // '::sqrt' use that of <math.h> !
      }


      inline operator float () const {  // TODO: to be replaced with truncate
         //SX_CHECK (fabs (im) < 1e-4, im);
         return (float)re;
      }


      inline operator double () const {  // TODO: to be replaced with truncate
//         //SX_CHECK (fabs (im) < 1e-4, im);
         return (double)re;
      }


      void print () const { printf ("(%f,%f)\n", re, im); }

};

template<>
inline SxComplex<double> SxComplex<double>::phase (double phi)
{
   SxComplex<double> res;
   sincos (phi, &res.im, &res.re);
   return res;
}

template<>
inline SxComplex<float> SxComplex<float>::phase (float phi)
{
   SxComplex<float> res;
   sincosf (phi, &res.im, &res.re);
   return res;
}

template<class T>
template<class U>
SxComplex<T>::SxComplex (const SxComplex<U> &c)
{
   re = c.re;
   im = c.im;
}


template<class T>
inline SxComplex<T> operator+ (const T &r, const SxComplex<T> &c)
{
   return SxComplex<T> (r + c.re, c.im);
}

template<class T>
inline SxComplex<T> operator- (const T &r, const SxComplex<T> &c)
{
   return SxComplex<T> (r - c.re, -c.im);
}

template<class T>
inline SxComplex<T> operator* (const T &r, const SxComplex<T> &c)
{
   return SxComplex<T> (c.re * r, c.im * r);
}

template<class T>
inline SxComplex<T> operator/ (const T &r, const SxComplex<T> &c)
{
   T t = r / ( c.re * c.re + c.im * c.im );
   return SxComplex<T> (t * c.re, -t * c.im);
}


template<class T>
inline SxComplex<T> sqrt (const SxComplex<T> &c)
{
   T r = c.abs ();
   T nR, nI;
   if (r < 1e-10)
      nR = nI = r;
   else if ( c.re > 0.)  {
      nR = ::sqrt (0.5 * (r + c.re));  // '::sqrt' use that of <math.h> !
      nI = c.im / nR / 2.;
   }  else  {
      nI = ::sqrt (0.5 * (r - c.re));  // '::sqrt' use that of <math.h> !
      if (c.im < 0.)  nI = - nI;
      nR = c.im / nI / 2.;
   }
   return SxComplex<T> (nR, nI);
}


inline int               conjugate (int in)    { return in; }
inline float             conjugate (float in)  { return in; }
inline double            conjugate (double in) { return in; }
inline SxComplex<float>  conjugate (const SxComplex<float> &in)
{ return in.conj(); }
inline SxComplex<double> conjugate (const SxComplex<double> &in)
{ return in.conj(); }

inline float  toReal (float in)  { return in; }
inline double toReal (double in) { return in; }
inline float  toReal (const SxComplex<float>  &in) { return in.re; }
inline double toReal (const SxComplex<double> &in) { return in.re; }
inline float  toImag (float)  { return 0.; }
inline double toImag (double) { return 0.; }
inline float  toImag (const SxComplex<float>  &in) { return in.im; }
inline double toImag (const SxComplex<double> &in) { return in.im; }


template<class T>
std::ostream& operator<< (std::ostream &s, const SxComplex<T> &in)
{
   int width = int(s.width());
   return s << std::setw(1)     << "("
            << std::setw(width) << in.re << ","
            << std::setw(width) << in.im << ")";
}


typedef SxComplex<float>  SxComplex8;
typedef SxComplex<double> SxComplex16;

inline SxComplex16 operator*(const SxComplex8 &c8, const SxComplex16& c16)
{
   return SxComplex16 (c8.re * c16.re - c8.im * c16.im,
                       c8.re * c16.im + c8.im * c16.re);
}

inline SxComplex16 operator*(const SxComplex16& c16, const SxComplex8 &c8)
{
   return SxComplex16 (c8.re * c16.re - c8.im * c16.im,
                       c8.re * c16.im + c8.im * c16.re);
}

template<class T>
inline SxComplex<T> exp(const SxComplex<T> &in)
{
   T expre = exp(in.re);
   return SxComplex<T>(expre * cos(in.im), expre * sin(in.im));
}

inline bool isValid (const int &)
{
   return true;
}
inline bool isValid (const float &in)
{
   if (sxIsNan (in))  return false;
   if (sxIsInf (in))  return false;
   return true;
}
inline bool isValid (const double &in)
{
   if (sxIsNan (in))  return false;
   if (sxIsInf (in))  return false;
   return true;
}
inline bool isValid (const SxComplex8 &in)
{
   if (sxIsNan (in.re) || sxIsNan (in.im))  return false;
   if (sxIsInf (in.re) || sxIsInf (in.im))  return false;
   return true;
}
inline bool isValid (const SxComplex16 &in)
{
   if (sxIsNan (in.re) || sxIsNan (in.im))  return false;
   if (sxIsInf (in.re) || sxIsInf (in.im))  return false;
   return true;
}


#ifdef NDEBUG
#   define SX_CHECK_DIV(denom)                ((void)0)
#   define SX_CHECK_NUM(expr)                 ((void)0)
#else
    inline void checkDivision(int i)  {
       if (i == 0)  {
          printf ("Division by zero!\n"); SX_EXIT;
       }
    }
    inline void checkDivision(float x)  {
       if (fabs(x) < 1e-50)  {
          printf ("Division by zero! (x=%g)\n", x); SX_EXIT;
       }
    }
    inline void checkDivision(double x)  {
       if (fabs(x) < 1e-50)  {
          printf ("Division by zero! (x=%g)\n", x); SX_EXIT;
       }
    }
    inline void checkDivision(const SxComplex8 &c)  {
       if (fabs(c.re) + fabs(c.im) < 1e-50)  {
          printf ("Division by zero! (c=%g,%g)\n", c.re, c.im); SX_EXIT;
       }
    }
    inline void checkDivision(const SxComplex16 &c)  {
       if (fabs(c.re) + fabs(c.im) < 1e-50)  {
          printf ("Division by zero! (c=%g,%g)\n", c.re, c.im); SX_EXIT;
       }
    }

    inline void checkNumber (int)  { /* empty */ }
    template<class T> inline void checkNumber (const T &x)  {
       if (!isValid(x))  {
          printf ("Number is not initialized or invalid (NaN or Inf)!\n"); SX_EXIT;
       }
    }
#   define SX_CHECK_DIV(denom)        checkDivision(denom)
#   define SX_CHECK_NUM(expr)         checkNumber(expr)
#endif /* NDEBUG */



#endif // _SX_CMPLX_H_
