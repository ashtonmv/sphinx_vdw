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
#ifndef _LINEARREGRESSION_H_
#define _LINEARREGRESSION_H_

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


/** \brief ...

    \b SxClass = SPHInX ...

    ....

    \author John Doe, johndoe@mpie.de */

class SX_EXPORT_MATH SxLinearRegression
{
   public:
      // Constructors
      // -- Constructor
      SxLinearRegression ();

      SxLinearRegression (const SxVector<Double> &x, const SxVector<Double> &y);

      // -- Destructor
      ~SxLinearRegression ();

      // Methods
      double getVal (double x);

      SxString getPolynomial ();

      double getCoeff ();
      
      double getTransl ();

      // Static methods
      static bool testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y);

      static SxString getErrStr ();

   protected:
      // Methods
      void computeLeastSquares (const SxVector<Double> &x,const SxVector<Double> &y);

      // Members
      SxVector<Double> xVals;   

      double p1, p2;
};

#endif /* _LINEARREGRESSION_H_ */
