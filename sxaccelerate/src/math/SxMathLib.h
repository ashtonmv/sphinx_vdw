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
#ifndef _SX_MATHLIB_H_
#define _SX_MATHLIB_H_

#include <SxMath.h>
#include <SxComplex.h>
#include <SxUtil.h>
#include <SxSort.h>
#include <SxConfig.h>
#include <SxBlasLib.h>

#ifdef USE_SXMATH
   SX_EXPORT_MATH void initSPHInXMath ();
#else
   SX_EXPORT_MATH void initSPHInXMath () { }
#endif
SX_EXPORT_MATH void destroySFHIngXMath ();

//------------------------------------------------------------------------------
// Minimum, maximum, round
//------------------------------------------------------------------------------
template<class T> inline T minimum (const T a, const T b) { return a < b ? a : b; }
template<class T> inline T maximum (const T a, const T b) { return a > b ? a : b; }
template<class T> inline T sqr (const T x) { return x * x; }

//------------------------------------------------------------------------------
// in-place matrix multiplication
//------------------------------------------------------------------------------
SX_EXPORT_MATH void inPlaceRot (float *mat, const float *rotMat, 
                                int nRows, int nCols);
SX_EXPORT_MATH void inPlaceRot (double *mat, const double *rotMat, 
                                int nRows, int nCols);
SX_EXPORT_MATH void inPlaceRot (SxComplex8 *mat, const SxComplex8 *rotMat,
                                int nRows, int nCols);
SX_EXPORT_MATH void inPlaceRot (SxComplex16 *mat,  const SxComplex16 *rotMat, 
                                int nRows, int nCols);


//------------------------------------------------------------------------------
// Erf and Erfc
//------------------------------------------------------------------------------

/** Returns the error function erf(x). */
SX_EXPORT_MATH double derf (double);
/** Returns the complementary error function erfc(x). */
SX_EXPORT_MATH double derfc (double);
/** Return exp(x^2) * erfc(x) without overflow for large x */
SX_EXPORT_MATH double erfcexp (double x);
/** Returns the incomplete gamma function P(a,x). */
SX_EXPORT_MATH double gammp (double, double);
/** Returns the incomplete gamma function Q(a,x)=1-P(a,x) */
SX_EXPORT_MATH double gammq (double, double);
/** Returns the value $ln |\Gamma(x)|$ for x > 0. */
SX_EXPORT_MATH double gammln (double);
/** Returns the incomplete gamma function Q(a,x) evaluated by its continued
    fraction representation as gammcf. Also returns $ln \Gamma(a)$ as gln. */
SX_EXPORT_MATH void gcf (double *, double, double, double *);
/** Returns the incomplete gamma function P(a,x) evaluated by its series
    representation as gamser. Also returns $ln \Gamma(a)$ as gln. */
SX_EXPORT_MATH void gser (double *, double, double, double *);

#ifdef CYGWIN
   extern "C" double cbrt (double);
#endif

#ifdef WIN32
   SX_EXPORT_MATH double cbrt (double);
#endif

#ifndef HAVE_SINCOS
   SX_EXPORT_MATH void sincos (double, double *, double *);
   SX_EXPORT_MATH void sincosf (float, float *, float *);
#endif

#ifndef HAVE_ROUND
inline double round (double x) {
   return (double)( (long int)(x + (x > 0. ? +0.5 : -0.5 ) ) );
}
#endif
#ifndef HAVE_LROUND
inline long int lround (double x)
{
   return (long int) round(x);
}
#endif

// TODO: should go to some general minimization library
SX_EXPORT_MATH double lineMinimization (double y, double xT, double yT, 
                                        double dYdX, double *curvature=NULL);

#endif /* _SX_MATHLIB_H_ */
