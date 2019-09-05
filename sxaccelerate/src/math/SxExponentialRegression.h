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
#ifndef _EXPONENTIALREGRESSION_H_
#define _EXPONENTIALREGRESSION_H_

// Including the header-files
// -- SPHInX-Headers
#include <SxVector.h>
#include <SxArray.h>
#include <SxMath.h>
#include <SxPtr.h>
#include <SxString.h>

// -- C++-Standardheaders
#include <iostream>
#include <cmath>
using namespace std;


/** \brief ...

  \b SxClass = SPHInX ...

precondition:
all function values are positive

functionality:
First one converts the exponential least square problem in a linear one:
f(x) = a * exp (b * x)
=> log (f (x)) = log (a) + b * x
y := log (f (x))
p[0] := log (a)
p[1] := b
=> p[0] + p[1] * x = y

Then one solves this linear least square task.

Afterwards one rechanges the problem to an exponential one:
a = exp (p[0])
b = p[1]
=> approx (x) = a * exp (b * x)

\author Thomas Uchdorf, t.uchdorf@mpie.de */
class SX_EXPORT_MATH SxExponentialRegression
{

   public:

      // Constructors

      // -- Constructor
      /** \brief Creates an ExponentialRegression object*/
      SxExponentialRegression ();

      /** \brief Creates an ExponentialRegression object*/      
      SxExponentialRegression (const SxVector<Double> &x,
            const SxVector<Double> &y);

      // -- Destructor
      /** \brief Destroys an ExponentialRegression object */
      ~SxExponentialRegression ();

      // Methods

      /** \brief Returns the value of the exponential regression
        function at the requested position */
      double getVal (double x);

      /** \brief Provides the formula of the exponential regression
        function as string */      
      SxString getFunc ();

      // Static methods
      static bool testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y);

      static SxString getErrStr ();

   protected:
      
      // Methods

      /** \brief Converts the exponential least square problem to a linear one
       */
      SxVector<Double> convertToLinearProblem (const SxVector<Double> &y);

      /** \brief Rechanges the linear problem to an exponential one */
      void convertLinSolToExpProblem ();

      /** \brief Solves the linear least square problem */
      void computeLeastSquares (const SxVector<Double> &x,const SxVector<Double> &y);

      // Members

      /** \brief f (x)= p1 * exp (p2 * x) */
      double p1, p2;

};

#endif /* _EXPONENTIALREGRESSION_H_ */
