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

// Headerguard
#ifndef _NATURALCUBICSPLINE_H_
#define _NATURALCUBICSPLINE_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxVector.h>
#include <SxArray.h>
#include <SxPtr.h>
#include <SxString.h>
#include <SxMath.h>

// -- C++-Standardheaders
#include <iostream>
#include <cmath>
using namespace std;


/** \brief Natural cubic spline interpolation

    \b SxClass = SFHIngX natural cubic spline interpolation

    \author T. Uchdorf, uchdorf@mpie.de
    \author Christoph Freysoldt, freysoldt@mpie.de
 */

class SX_EXPORT_MATH SxNaturalCubicSpline
{
   public:
      // Constructors
      // -- Constructor
      SxNaturalCubicSpline ();
      
      /// Constructor for x-y interpolation
      SxNaturalCubicSpline (const SxVector<Double> &x, const SxVector<Double> &y);
      /// Constructor for y interpolation (index serves as x)
      template <class VectorType>
      SxNaturalCubicSpline (const VectorType &y, bool estimateD2 = false);

      template <class VectorType>
      void compute (const VectorType &y, bool estimateD2 = false);

      // -- Destructor
      ~SxNaturalCubicSpline ();

      inline double getVal (double x) const
      {
         return (xVals.getSize () == 0) ? getValY (x) : getValXY (x);
      }

      template <class VectorType>
      VectorType getVals (const VectorType &x) const;

      /// Get value for i->y(i) even outside original data range
      double getValYExtra (double x) const;
      /// Get derivatives for i->y(i) even outside original data range
      double getDerivYExtra (double x) const;

  protected:

      // Methods
      void print (
            const SxVector<Double> &diag1,
            const SxVector<Double> &diag2,
            const SxVector<Double> &diag3,
            const SxVector<Double> &b);

      void compute (const SxVector<Double> &x,const SxVector<Double> &y);
      
      SxVector<Double> symmetricTridiagonalGauss (
            const SxVector<Double> & upper,
            const SxVector<Double> &mid,
            const SxVector<Double> & lower,
            const SxVector<Double> &rhs);

      // Members
      SxArray<SxVector<Double> > polyCoeff;

      SxVector<Double> xVals;

      /// Get value for i->y(i) interpolation 
      double getValY (double x) const;
      /// Get value for x(i)->y(i) interpolation 
      double getValXY (double x) const;
  public:
      const SxArray<SxVector<Double> > & getCoeff () const
      {
         return polyCoeff;
      }
};

template <class VectorType>
SxNaturalCubicSpline::SxNaturalCubicSpline (const VectorType &y,
                                            bool estimateD2)
{
   compute(y, estimateD2);
}

template <class VectorType>
void SxNaturalCubicSpline::compute (const VectorType &y, bool estimateD2)
{
   int n = (int)y.getSize () - 1;
   SX_CHECK (n > 2, n);

   // prepare internal state
   xVals.resize (0);
   polyCoeff.resize (2);

   // set up tridiagonal matrix for coefficients
   SxVector<Double> dl(n),d(n+1),du(n);
   d(0) = 1.;
   du(0) = 0.;
   // h_i := x_{i+1} - x_i
   for (int i = 1; i < n; ++i)   {
      dl(i-1) = 1.; /* h_{i-1} */
      d (i)   = 4.; /* 2(h_{i-1} + h_i) */
      du(i)   = 1.; /* h_i */
   }
   dl(n-1) = 0.;
   d(n) = 1.;

   // --- set up right-hand side
   SxVector<Double> rhs(n+1);
   if (estimateD2)
      rhs(0) = 2. * y(0) - 5. * y(1) + 4. * y(2) - y(3);
      //rhs(0) = 3. * y(0) - 9. * y(1) + 10. * y(2) - 5. * y(3) + y(4);
   else
      rhs(0) = 0.;
   for (int i = 1; i < n; ++i)  {
      // rhs(i) = 6 * ((y_{i+1}-y_i)/h_i - (y_i - y_{i-1})/h_{i-1})
      rhs(i) = 6. * (y(i+1) + y(i-1) - 2. * y(i));
   }
   if (estimateD2)
      rhs(n) = 2. * y(n) - 5. * y(n-1) + 4. * y(n-2) - y(n-3);
   else
      rhs(n) = 0.;
   // solve tdm * x = rhs
   polyCoeff(1) = symmetricTridiagonalGauss (du, d, dl, rhs);
   polyCoeff(1) /= 6.; // prefactor

   // --- copy original y-data
   SxVector<Double> &a0 = polyCoeff(0);
   a0.resize (y.getSize ());
   for (int i = 0; i < a0.getSize (); ++i)
      a0(i) = y(i);
 
}

inline double SxNaturalCubicSpline::getValY (double x) const
{
   int i = int(floor(x));
   if (i == -1 && fabs(x) < 1e-12) return polyCoeff(0)(0);
   const SxVector<Double> &data = polyCoeff(0);
   if (i == data.getSize () - 1 && fabs(x - i) < 1e-12) return data(i);
   SX_CHECK (i >= 0 && i < data.getSize () - 1, i, data.getSize ());
   double t = x -i, u = 1. - t;
   SX_CHECK (t >= -1e-12 && t <= 1.0000000000001, t);
   SX_CHECK (u >= -1e-12 && u <= 1.0000000000001, u);
   const SxVector<Double> &z = polyCoeff(1);
   // S= (z_{i+1}(x-x_i)^3 + z_i(x_{i+1}-x)^3) / (6 h_i)
   //    + (y_{i+1}/h_i - h_i*z_{i+1}/6)(x-x_i)
   //    + (y_{i}/h_i - h_i*z_{i}/6)(x_{i+1}-x)
   double s = (data(i+1) + z(i+1)*(t*t - 1.)) * t
            + (data(i) + z(i)*(u*u - 1.)) * u;
   return s;
}

inline double SxNaturalCubicSpline::getValYExtra (double x) const
{
   SX_CHECK (xVals.getSize () == 0);
   SX_CHECK (polyCoeff.getSize () == 2, polyCoeff.getSize ());
   int i = int(floor(x));
   const SxVector<Double> &data = polyCoeff(0);
   if (i < 0) i = 0;
   else if (i >= data.getSize () - 1) i = static_cast<int>(data.getSize ()) - 2;
   double t = x -i, u = 1. - t;
   const SxVector<Double> &z = polyCoeff(1);
   // S= (z_{i+1}(x-x_i)^3 + z_i(x_{i+1}-x)^3) / (6 h_i)
   //    + (y_{i+1}/h_i - h_i*z_{i+1}/6)(x-x_i)
   //    + (y_{i}/h_i - h_i*z_{i}/6)(x_{i+1}-x)
   double s = (data(i+1) + z(i+1)*(t*t - 1.)) * t
            + (data(i) + z(i)*(u*u - 1.)) * u;
   return s;
}

inline double SxNaturalCubicSpline::getDerivYExtra (double x) const
{
   SX_CHECK (xVals.getSize () == 0);
   SX_CHECK (polyCoeff.getSize () == 2, polyCoeff.getSize ());
   int i = int(floor(x));
   const SxVector<Double> &data = polyCoeff(0);
   if (i < 0) i = 0;
   else if (i >= data.getSize () - 1) i = static_cast<int>(data.getSize ()) - 2;
   double t = x -i, u = 1. - t;
   const SxVector<Double> &z = polyCoeff(1);
   // S'= (3z_{i+1}(x-x_i)^2 - 3z_i(x_{i+1}-x)^2) / (6 h_i)
   //    + (y_{i+1}/h_i - h_i*z_{i+1}/6)
   //    - (y_{i}/h_i - h_i*z_{i}/6)
   double s = (data(i+1) + z(i+1)*(3. * t*t - 1.))
            - (data(i)   + z(i)  *(3. * u*u - 1.));
   return s;
}

template <class VectorType>
VectorType SxNaturalCubicSpline::getVals (const VectorType &x) const
{
   VectorType res;
   int n = (int)x.getSize ();
   res.resize (n);
   for (int i = 0; i < n; ++i)
      res(i) = getVal (x(i));
   return res;
}


#endif /* _NATURALCUBICSPLINE_H_ */
