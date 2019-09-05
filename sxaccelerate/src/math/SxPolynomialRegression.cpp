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

// Including the header-file
#include <SxPolynomialRegression.h>

//------------------------------------------------------------------------------
// SxPolynomialRegression-class
//------------------------------------------------------------------------------
// Constructors
// -- Constructor
SxPolynomialRegression::SxPolynomialRegression ()
{
   
   // empty

}

SxPolynomialRegression::SxPolynomialRegression (const SxVector<Double> &x, const SxVector<Double> &y, int n, Solver solver)
   : coeff (n + 1)
{

   // Checking whether the number of x-values and y-values equals and if the 
   // order is at least two less than the number of points (otherwise it would 
   // be an interpolation respectively total nonsense)
   SX_CHECK (x.getSize () == y.getSize () && x.getSize () > 0 && n > 0 &&
             x.getSize () > n + 1, x.getSize (), y.getSize (), n);
   
   // Computing the regression polynomial
   computeLeastSquares (x, y, solver);

}

// -- Destructor
SxPolynomialRegression::~SxPolynomialRegression ()
{

   // empty

}

// Methods
void SxPolynomialRegression::computeLeastSquares (const SxVector<Double> &x, const SxVector<Double> &y, Solver solver)
{
   
   SxMatrix<Double> A (x.getSize (), coeff.getSize ());
   for (int iLin = 0; iLin < A.nRows (); ++iLin)  {

      for (int iCol = 0; iCol < A.nCols (); ++iCol)  {


         //    ( 1    x[0] x[0]^2 ... x[0]^order )
         //	A =( 1    x[1] x[1]^2 ... x[1]^order )
         //    ( ...  ...  ...    ... ...        )
         //    ( 1    x[n] x[n]^2 ... x[n]^order )
         A(iLin, iCol) = ::pow (x(iLin), iCol);

      }

   }

   SxMatrix<Double> At (A.transpose ());
 
   // A^T A
   SxMatrix<Double> lhs (At ^ A);

   // A^T y
   SxVector<Double> rhs (At ^ y);

   SxMatrix<Double> Q, R;

   // Selecting requested solution method
   if (solver == NormalEquation)  {

      // A^T A x = A^T y <=> x = (A^T A)^(-1) A^T y
      coeff = (SxVector<Double>) (lhs.inverse () ^ rhs);

   } else {

      householder (lhs, rhs, Q, R, coeff);
      
   }

}

double SxPolynomialRegression::getVal (double x)
{
  
   double ret (0.);
   
   for (int i = 0; i < coeff.getSize (); ++i)  {

      ret += coeff(i) * ::pow (x, i);
      
   }
   
   return ret;
   
}

SxString SxPolynomialRegression::getPolynomial ()
{
  
   SxString ret;
   
   for (int i = 0; i < coeff.getSize (); ++i)  {

      ret += SxString (coeff(i)) +" * x^" + SxString (i);

      if (i < coeff.getSize () - 1)  {

         ret += " + ";

      }
      
   }
   
   return ret;
   
}

int SxPolynomialRegression::sign (double x)
{

   if (x < 0)  {

      return  -1;

   } else if (x > 0)  {

      return 1;

   } else  {

      return 0;

   }

}

void SxPolynomialRegression::householder (const SxMatrix<Double> &A, const SxVector<Double> &b,
      SxMatrix<Double> &Q, SxMatrix<Double> &R, SxVector<Double> &x)
{
   // Checking whether the number of lines and columnes correspond
   SX_CHECK (A.nRows () == A.nCols (), A.nRows (), A.nCols ());
   
   SxMatrix<Double> E (A.nRows (), A.nCols ());
   SxMatrix<Double> Qi;
   SxMatrix<Double> QTrans;
   SxVector<Double> y;
   SxVector<Double> a1 (A.nRows (), 0.);
   SxVector<Double> v1 (A.nRows (), 0.);  

   // Creating an identity matrix
   //  1   0   0  ...  0   0   0
   // ...  1   0  ...  0   0   0
   // ... ... ... ... ... ... ...
   // ...  0   0  ...  0   1   0
   // ...  0   0  ...  0   0   1

   for (int iLin = 0; iLin < E.nRows (); ++iLin)  {

      for (int iCol = 0; iCol < E.nCols (); ++iCol)  {

         if (iLin == iCol)  {

            E(iLin, iCol) = 1.;

         } else  {

            E(iLin, iCol) = 0.;

         }

      }   

   }
  
   // Initializing R with A
   R.copy (A);

   // Initializing Q^T with the identity matrix
   QTrans.copy (E);

   for (int iCol = 0; iCol < A.nCols () - 2; ++iCol)  {

      // 1st step
      // Defining a1 as current column vector of R without the values which lie
      // above the diagonal
      a1 = SxVector<Double> (A.nRows (), 0.);
      
      // v1 equals a1 plus or minus the current unit vector multiplied by the
      // norm of a1
      v1 = SxVector<Double> (A.nRows (), 0.);

      // Initializing a1 and partly v1
      a1(iCol) = R(iCol, iCol);
      for (int iLin = iCol + 1; iLin < A.nRows (); ++iLin)  {

            a1(iLin) = R(iLin, iCol);                  
            v1(iLin) = a1(iLin);
      }

      // Determining the last v1 value
      v1(iCol)  = a1(iCol) + sign (a1(iCol)) * sqrt ((a1 * a1).sum ());

      // 2nd step
      // Q[i] := E - 2. / ||v1||[2]^2 * (v1 * v1^T)
      Qi = E - 2. / (v1 * v1).sum () * (v1 ^ v1.transpose ());

      // R = prod (Q[i], i = 0, ..., iCol) * A 
      R = Qi ^ R;
      
      // Q^T = prod (Q[i], i = 0, ..., iCol)
      QTrans = (Qi ^ QTrans);
      
   }

   // Computing Q (Q R x = Q Q^T A x = b)
   Q = QTrans.transpose ();
  
   // Determining rhs
   y = (QTrans ^ b);
  
   // Solving equation (Q R x = b <=> R x = Q^T b = y <=> x = R^(-1) y) via backward substitution
   x = (R.inverse () ^ y);
   
}


bool SxPolynomialRegression::testPreconditions (const SxVector<Double> &x, const SxVector<Double> &y, int n)
{

   // Checking whether the number of x-values and y-values equals and if the 
   // order is at least two less than the number of points (otherwise it would 
   // be an interpolation respectively total nonsense)
   if (x.getSize () != y.getSize () || x.getSize () <= 0 ||
         n <= 0 || x.getSize () <= n + 1)  {

      return false;

   } else   {
  
      return true;

   }
   
}

SxString SxPolynomialRegression::getErrStr ()
{

   return SxString ("The number of x-values must be positive and equal the number of y-values.\nThe order of the requested polynomial must be at least two less than the number of points\n(otherwise an interpolation would be desired).");
   
}

