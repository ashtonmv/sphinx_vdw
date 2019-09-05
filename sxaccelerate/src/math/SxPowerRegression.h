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
#ifndef _POWERREGRESSION_H_
#define _POWERREGRESSION_H_

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
all x- and y-values are positive

functionality:
First one converts the power least square problem to a linear one:
f (x) = a * x ^ b
<=> log (f (x)) = log (a * x ^ b)
<=> log (f (x)) = log (a) + log (x ^ b)
<=> log (f (x)) = log (a) + b * log (x)
<=> out = transl + coeff * in, transl := log (a), out := log (f (x)), in := log (x)

After the linear least square problem is solved one resubstitutes the linear
problem by a power problem
transl = log (a)
<=> a = 10 ^ transl
=> approx (x) = a * x ^ b

\author Thomas Uchdorf, t.uchdorf@mpie.de */
class SX_EXPORT_MATH SxPowerRegression
{

   public:

      // Constructors

      // -- Constructor
      /** \brief Creates an InverseRegression object*/
      SxPowerRegression ();

      /** \brief Creates an InverseRegression object*/      
      SxPowerRegression (const SxVector<Double> &x,
            const SxVector<Double> &y);

      // -- Destructor
      /** \brief Destroys an InverseRegression object */
      ~SxPowerRegression ();

      // Methods

      /** \brief Returns the value of the power regression
        function at the requested position */
      double getVal (double x);

      /** \brief Provides the formula of the power regression
        function as string */      
      SxString getFunc ();

      // Static methods
      static bool testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y);

      static SxString getErrStr ();

   protected:
      
      // Methods

      /** \brief Rechanges the linear problem to an power one */
      void convertLinSolToPowProblem ();
      
      // Members

      /** \brief f (x) = transl * x ^ coeff */
      double transl, coeff;

};

#endif /* _POWERREGRESSION_H_ */
