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
#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxMatrix.h>
#include <SxMath.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxPtr.h>
#include <SxString.h>
// -- C++-Standardheaders
#include <iostream>
#include <cmath>
using namespace std;


/** \brief ...

    \b SxPolynomialInterpolation1D = SPHInX ...

    ....

    \author Thomas Uchdorf, t.uchdorf@mpie.de */


// InterpolationMethod-enumeration
/** \brief Determines the symbol which is used for plotting a point */
enum InterpolationMethod {Lagrange, DividedDifferences};

class SX_EXPORT_MATH SxPolynomialInterpolation1D
{
   public:
      // Constructors
      // -- Constructor
      SxPolynomialInterpolation1D ();
      SxPolynomialInterpolation1D (const SxVector<Double> &x, const SxVector<Double> &y, const InterpolationMethod &method);

      // -- Destructor
      ~SxPolynomialInterpolation1D ();

      // Methods
      double getVal (const double &x) const;

   protected:
      // Methods
      void computeLagrange (const SxVector<Double> &x,const SxVector<Double> &y);

      void computeDividedDifferences (const SxVector<Double> &x,const SxVector<Double> &y);   

      double getDividedDifferencesVal (const double &x) const;

      double getLagrangeVal (const double &x) const;
      
      // Members
      InterpolationMethod computationMethod;
      
      SxVector<Double> xVals;

      SxVector<Double> coeff;
};

#endif /* _INTERPOLATION_H_ */
