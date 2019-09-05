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

#include <SxError.h>
#include <SxInterpolation1D.h>

//------------------------------------------------------------------------------
// SxInterpolation1D-class
//------------------------------------------------------------------------------

// Constructors
// -- Constructor
SxPolynomialInterpolation1D::SxPolynomialInterpolation1D ()
{
   // empty
}

SxPolynomialInterpolation1D::SxPolynomialInterpolation1D (const SxVector<Double> &x, const SxVector<Double> &y, const InterpolationMethod &method) :
   computationMethod (method), xVals (x), coeff (x.getSize ()) 
{
   SX_CHECK (x.getSize () == y.getSize ());
   if (method == Lagrange)  { 
      computeLagrange (x, y);
   } else if (method == DividedDifferences)  {
      computeDividedDifferences (x, y);   
   }
}

// -- Destructor
SxPolynomialInterpolation1D::~SxPolynomialInterpolation1D ()
{
   // empty
}

// Methods
void SxPolynomialInterpolation1D::computeLagrange (const SxVector<Double> &x,const SxVector<Double> &y)
{
   // Lagrange formula
   // p[n-1] (x) := sum (prod (x - x[j], j != i) / prod (x[i] - x[j], j != i) * y[i], i = 0 .. n)    
   // Here the coefficient of the denominator product and the adopted y-value is computed.
   // That means coeff[i] := y[i] / prod (x[i] - x[j], j != i), i = 0 .. n
   
   for (int i = 0; i < xVals.getSize (); ++i)  {
      double denominator (1.);
      for (int iLeft = 0; iLeft < i; ++iLeft)  {
         denominator *= (x(i) - x (iLeft));
      }
      for (int iRight = i + 1; iRight < xVals.getSize (); ++iRight)  {
         denominator *= (x(i) - x (iRight));
      } 
      coeff(i) = y(i) / denominator;  
   }    
}

void SxPolynomialInterpolation1D::computeDividedDifferences (const SxVector<Double> &x, const SxVector<Double> &y)
{
   /*
   // background:
   // -----------
   // x matrix m (column 0 of m equals y)
   // 0 3 
   //     \  
   // 1 2 - (2 - 3) / (1 - 0) = -1
   //     \                        \
   // 3 6 - (6 - 2) / (3 - 1) = 2  - (2 - (-1)) / (3 - 0) = 1
   // m[i][0] := y[i], i = 0 .. n
   // m[i][j] := (m[i][j - 1] - m[i-1][j-1]) / (x[i] - x[i - j])
   // p[n] (x) = sum (m[i][i] * prod (x - x[j], j = 0 .. i-1), i = 0 .. n)
   // coeff represents the diagonal values of m
   // p[n] (x) = coeff[0] + coeff[1] * (x - xVals[0]) + ... + coeff[n] * (x - xVals[0]) * ... * (x-xVals[n]) 
   */
   SxMatrix<Double> m (x.getSize (), x.getSize ());

   for (int iLin = 0; iLin < x.getSize (); ++iLin)  {
      m(iLin, 0) = y(iLin);
   }

   coeff(0) = m(0, 0);
   for (int iLin = 1; iLin < x.getSize (); ++iLin)  {
      for (int iCol = 1; iCol <= iLin; ++iCol)  {
			m(iLin, iCol) = (m(iLin, iCol - 1) - m(iLin - 1, iCol - 1)) / (x(iLin) - x(iLin - iCol));
      }
   	coeff(iLin) = m(iLin, iLin);
   }  

   //TEST
   for (int iLin = 0; iLin < x.getSize (); ++iLin)  {
      for (int iCol = 0; iCol <= iLin; ++iCol)  {
         cout << m(iLin, iCol) << " ";          
      }   
      cout << endl;      
   }
   
}

double SxPolynomialInterpolation1D::getVal (const double &x) const
{
   double ret (0.);
   if (computationMethod == Lagrange)  { 
      ret = getLagrangeVal (x);
   } else if (computationMethod == DividedDifferences)  {
      ret = getDividedDifferencesVal (x);   
   }
   return ret;   
}

double SxPolynomialInterpolation1D::getDividedDifferencesVal (const double &x) const
{
   double help (1.), ret (0.);
   for (int iLin = 0; iLin < xVals.getSize (); ++iLin) {
      if (iLin - 1 >= 0){    
         help *= (x - xVals(iLin - 1));
      }
      ret += coeff(iLin) * help;
      cout << "coeff " << coeff(iLin)<< " help: "<< help << endl;
   }
   return ret;
}

double SxPolynomialInterpolation1D::getLagrangeVal (const double &x) const
{
   double ret (0.);
   for (int i = 0; i < xVals.getSize (); ++i)  {
      double nominator (1.);
      for (int iLeft = 0; iLeft < i; ++iLeft)  {
         nominator *= (x - xVals (iLeft));
      }
      for (int iRight = i + 1; iRight < xVals.getSize (); ++iRight)  {
         nominator *= (x - xVals (iRight));
      }
      ret += coeff(i) * nominator;    
   }   
   return ret;
}

/* 
//------------------------------------------------------------------------------
// main-function
//------------------------------------------------------------------------------
int main (void)
{
   // --- initialize the SPHInX Algebra Library
   initSPHInXMath ();
   initSxUtil ();

   SxVector<Double> x(3), y(3);
   for (int i = 0; i < x.getSize (); ++i)
      x(i) = i;
   y(0) = 3;
   y(1) = 2;
   y(2) = 6;
   SxPolynomialInterpolation1D interpol (x, y, DividedDifferences);
   for (double xVal = 0.; xVal <= 2.; xVal += 0.1)
      cout << xVal << " " << interpol.getVal (xVal) << endl;   
   
   
   SxVector<Double> x(5), y(5);
   x(0) = -2.;
   x(1) = -1.;
   x(2) = 0.;
   x(3) = 1.;
   x(4) = 2.;
   y(0) = 1./5;
   y(1) = 1./2;
   y(2) = 1.;
   y(3) = 1./2;
   y(4) = 1./5;
   
   SxNaturalCubicSpline spline (x, y);
   for (int i = 0; i < x.getSize (); ++i)  {
      cout << x(i) << " " << y(i)<< endl;
   }
   cout << endl;
   
   for (double xVal = -2.; xVal <= 2.; xVal += 0.1)
      cout << xVal << " " << spline.getVal (xVal) << endl;
   return 0;
}
*/   

