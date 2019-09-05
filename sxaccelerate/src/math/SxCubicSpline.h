//     
//                  S / P H I / n X
//                  http://www.sphinxlib.de
//     
//      Contact:    Sixten Boeck, boeck@mpie.de
//                  Algorithm Design and Modeling Group
//                  Computational Materials Design
//                  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//     
//      Authors:    see src/AUTHORS
//     
// ---------------------------------------------------------------------------

// Headerguard
#ifndef _CUBICSPLINE_H_
#define _CUBICSPLINE_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxVector.h>
#include <SxMatrix.h>
#include <SxArray.h>
#include <SxPtr.h>
#include <SxString.h>
#include <SxMath.h>
#include <SxError.h>
#include <SxTimer.h>

// -- C++-Standardheaders
#include <iostream>
#include <cmath>
using namespace std;


template <class V>
class SxCubicSpline
{

   typedef typename SxVecTraits<V>::MatType M;
   typedef typename SxVecTraits<V>::ScalarType S;
   typedef typename SxVecTraits<V>::ScalarType::Type s;

   public:

      //enum WorkingMode {Spline, LeastSquareFit};
      enum SplineType {Hermite, Natural, NaturalHermite, HermiteNatural};
      enum FitType {None, Normal, Extend, MirrorPlane, MirrorPoint};
      // Constructors
      // -- standart constructor
      SxCubicSpline ();

      /// Constructor for read spline
      SxCubicSpline (const V &x,
                     const V &spline);

      /// Constructor for x-y interpolation
      SxCubicSpline (const V &x, 
                     const V &y, 
                     const enum SplineType splineType,
                     const s slope1 = 0.0,
                     const s slope2 = 0.0);
      
      // Mixed Constructure due to same variable mask
      // When using Spline mode: 
      //   vec1 = x, vec2 = dx, and vec3 = y
      // When using LeastSquareFit mode:
      //   vec1 = xData, vec2 = yData, and vec3 = basis
      SxCubicSpline (const V &x, 
                     const V &dx,
                     const V &y, 
                     const enum SplineType splineType,
                     const s slope1 = 0.0,
                     const s slope2 = 0.0);

      SxCubicSpline (const V &xData, 
                     const V &yData,
                     const V &basis, 
                     const enum SplineType splineType,
                     const enum FitType fitType,
                     const s slope1 = 0.0,
                     const s slope2 = 0.0);

      // -- Destructor
      ~SxCubicSpline ();

      V getYFit () const;

      inline s getY (const s x) const;

      V getY (const V &x) const;

      inline s getdYdX (const s x) const;

      V getdYdX (const V &x) const;

      void setSpline (const V &spline);

      V getSpline () const;

      // calculate Splinepoint dependencys
      SxArray<int> getSplineDep (const enum FitType fitType);

      // set of meta data for dirac Vectors ! have to be set by hand
      int iSpecies;
      int iAtom;
      int n;
      int l;
      int m;

      inline void setXFit (const V &in) {xFit = in;};
      inline void setYFit (const V &in) {yFit = in;};
      inline void setXVals (const V &in) {xVals = in;};
      inline void setYVals (const V &in) {yVals = in;};
      inline void setSplineType (const enum SplineType in) {splineType = in;};

      M calcDRDataDRGrid ();

      
  protected:

      // Methods
      // Cubic Spline Interpolation
      void computeSpline ();
      // Least Square Fit
      void computeLeastSquareFit ();
      // Setup difference in xVals
      void setH ();
      // Setup coeff Matrix for Splinecoeff calculation
      void setT (const enum SplineType sType);
      // Setup right hand side for spline coeff calculation
      void setBeta (const enum SplineType sType);
      // get Spline Index
      inline int getSplineIdx(const typename S::Type x, const V& basis) const;
      // mirror
      void mirror ();
      // demirror
      void deMirror ();

      V symmetricTridiagonalGauss (const M &Mat, const V &rhs);

      void print (const V &diag1, const V &diag2, const V &diag3, const V &b);

      SxArray<SxList<int> > getSplineExtend (SxArray<int> &pointsPerSpline);
      SxVector<Int> getPointsPerSpline (const V& basis, const V& xData);
      V cleanZeroDepths(const V& xData, const SxVector<Int> &pointsPerSpline);
      M dady ();
      M dbdy (const M &dA, const M &dC);
      M dcdy (const enum SplineType sType);
      M dbetady (const enum SplineType sType);
      M dddy (const M &dC); 
      
      // Members
      V polyCoeff;
      V xVals;
      V yVals;
      V hVals;
      M TMat;
      V betaVals;
      V xFit;
      V yFit;
      SxArray<SxList<int> > depArray;
      enum SplineType splineType;
      enum FitType fitType;
      s leftHermite;
      s rightHermite;
};

namespace Timer {
   enum CubicSplineTimer   {
      computeSpline,
      leastSquareFit,
      solver
   };
}

SX_REGISTER_TIMERS(Timer::CubicSplineTimer)
{
   using namespace Timer;
   regTimer (computeSpline,           "spline calculation");
   regTimer (leastSquareFit,          "least square fit");
   regTimer (solver,                  "Ax = b solver");

}

#include <SxCubicSpline.hpp>

#endif /* _CUBICSPLINE_H_ */
