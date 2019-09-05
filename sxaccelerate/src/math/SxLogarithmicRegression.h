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
#ifndef _LOGARITHMICREGRESSION_H_
#define _LOGARITHMICREGRESSION_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxVector.h>
#include <SxArray.h>
#include <SxPtr.h>
#include <SxString.h>
#include <SxLinearRegression.h>
#include <SxMath.h>

// -- C++-Standardheaders
#include <iostream>
#include <cmath>
using namespace std;


/** \brief ...

  \b SxClass = SPHInX ...

precondition:
all x-values are positive

functionality:
First one converts the logarithmic least square problem in a linear one:
f (x) = a + b * log (x)
=> f (x) = a + b * y, y := log (x)

Then one solves this linear least square task.

\author Thomas Uchdorf, t.uchdorf@mpie.de */
class SX_EXPORT_MATH SxLogarithmicRegression
{

   public:

      // Constructors

      // -- Constructor
      /** \brief Creates an LogarithmicRegression object*/
      SxLogarithmicRegression ();

      /** \brief Creates an LogarithmicRegression object*/      
      SxLogarithmicRegression (const SxVector<Double> &x,
            const SxVector<Double> &y);

      // -- Destructor
      /** \brief Destroys an LogarithmicRegression object */
      ~SxLogarithmicRegression ();

      // Methods

      /** \brief Returns the value of the logarithmic regression
        function at the requested position */
      double getVal (double x);

      /** \brief Provides the formula of the logarithmic regression
        function as string */      
      SxString getFunc ();

      // Static methods
      static bool testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y);

      static SxString getErrStr ();

   protected:
      
      // Methods

      /** \brief Converts the logarithmic least square problem to a linear one
       */
      SxVector<Double> convertToLinearProblem (const SxVector<Double> &x);

      // Members

      /** \brief f (x)= transl + coeff * log (x) */
      double transl, coeff;

};

#endif /* _LOGARITHMICREGRESSION_H_ */
