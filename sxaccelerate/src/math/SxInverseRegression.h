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
#ifndef _INVERSEREGRESSION_H_
#define _INVERSEREGRESSION_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxVector.h>
#include <SxArray.h>
#include <SxMath.h>
#include <SxPtr.h>
#include <SxString.h>
#include <SxLinearRegression.h>

// -- C++-Standardheaders
#include <iostream>
#include <cmath>
using namespace std;


/** \brief ...

  \b SxClass = SPHInX ...

precondition:
all function values are non zero

functionality:
First one converts the inverse least square problem to a linear one:
f (x) = 1. / (a + b * x)
=> 1. / f (x) = a + b * x

Then one solves this linear least square task.

\author Thomas Uchdorf, t.uchdorf@mpie.de */
class SX_EXPORT_MATH SxInverseRegression
{

   public:

      // Constructors

      // -- Constructor
      /** \brief Creates an InverseRegression object*/
      SxInverseRegression ();

      /** \brief Creates an InverseRegression object*/      
      SxInverseRegression (const SxVector<Double> &x,
            const SxVector<Double> &y);

      // -- Destructor
      /** \brief Destroys an InverseRegression object */
      ~SxInverseRegression ();

      // Methods

      /** \brief Returns the value of the inverse regression
        function at the requested position */
      double getVal (double x);

      /** \brief Provides the formula of the inverse regression
        function as string */      
      SxString getFunc ();

      // Static methods
      static bool testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y);

      static SxString getErrStr ();


   protected:
      
      // Methods

      /** \brief Converts the inverse least square problem to a linear one
       */
      SxVector<Double> convertToLinearProblem (const SxVector<Double> &x);

      // Members

      /** \brief f (x) = 1. / (transl + coeff * x) */
      double transl, coeff;

};

#endif /* _INVERSEREGRESSION_H_ */
