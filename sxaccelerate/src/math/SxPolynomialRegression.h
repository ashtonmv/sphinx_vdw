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
#ifndef _POLYNOMIAL_REGRESSION_H_
#define _POLYNOMIAL_REGRESSION_H_

// Including the header-files

// -- SPHInX-Headers
#include <SxMatrix.h>
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

    \b SxPolynomialRegression = SFHIngX ...

    ....

    \author Thomas Uchdorf, t.uchdorf@mpie.de */
class SX_EXPORT_MATH SxPolynomialRegression
{

   public:

      /** \brief Represents two different methods to obtain a regression
        polynomial respectively to solve linear equations */
      enum Solver {
         NormalEquation = 0, Householder};
      
      // Constructors

      // -- Constructor
      /** \brief Creates a PolynomialRegression object */
      SxPolynomialRegression ();

      /** \brief Creates a PolynomialRegression object - means that the
        regression polynomial is directly computed */
      SxPolynomialRegression (const SxVector<Double> &x,
            const SxVector<Double> &y, int n, Solver solver);

      // -- Destructor
      /** \brief Destroys a PolynomialRegression object */
      ~SxPolynomialRegression ();

      // Methods

      /** \brief Returns the value of the calculated polynomial at the
        requested point */
      double getVal (double x);

      /** \brief Provides a formula string representing the calculated
        polynomial */
      SxString getPolynomial ();

      // Static methods
      static bool testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y, int n);

      static SxString getErrStr ();

   protected:

      // Methods

      /** \brief computes the regression polynomial using the requested
        method for solving linear equations */
      void computeLeastSquares (const SxVector<Double> &x, 
            const SxVector<Double> &y, Solver solver);

      /** \brief Returns 1 if x is positive, 0 if x is 0 or -1 if x is negative
       */
      int sign (double x);

      /** \brief Searches for representation like Q R x = b instead of A x = b
        where Q resembles an orthogonal matrix and R an upper triangular matrix
        so that R x = Q^T b is easily solvable by backward substitution
       */
      void householder (const SxMatrix<Double> &A, const SxVector<Double> &b,
      SxMatrix<Double> &Q, SxMatrix<Double> &R, SxVector<Double> &x);

      // Members
      
      /** \brief Represents the coefficients of the computed polynomial */
      SxVector<Double> coeff;

};

#endif /* _POLYNOMIAL_REGRESSION_H_ */
